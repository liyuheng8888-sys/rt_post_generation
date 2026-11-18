"""Physical regression tests: geometry, independent modal fields, and paired inputs."""
from pathlib import Path
import copy
import csv
import math
import sys
import tempfile
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import wave_packet as wp
generator=wp.import_file('packet_generator_test',wp.PROJECT/'generate_stations_source/generate_sources_stations.py')
post=wp.import_file('packet_post_test',wp.PROJECT/'postprocess_rt/postprocess_rt.py')


class PacketTests(unittest.TestCase):
    def setUp(self):self.c=wp.load_config()

    def test_fixed_geometry_and_independent_depths(self):
        rows=generator.station_rows(self.c)
        self.assertEqual(len(rows),786)
        self.assertEqual(sum(r['role']=='main' for r in rows),726)
        self.assertEqual(sum(r['role']=='main' and r['depth_role']=='holdout' for r in rows),242)
        a,b=wp.boundaries(self.c)
        self.assertTrue(all(not a<=r['burial_m']<=b for r in rows))
        self.assertEqual(len({(r['x_m'],r['y_m'],r['burial_m']) for r in rows}),786)
        bad=copy.deepcopy(self.c);bad['stations']['upper_offsets_m']=[0,20,40]
        with self.assertRaises(ValueError):generator.station_rows(bad)

    def test_integral_weighting_not_arithmetic_average(self):
        f=np.array([0.,1.,2.]);p=np.array([0.,1.,3.]);c={'R_P':np.array([0.,.2,.8])}
        value=wp.weighted_coefficients(f,p,c,['R_P'])['R_P']
        self.assertAlmostEqual(value,math.sqrt((.2**2+3*.8**2)/4))
        self.assertNotAlmostEqual(value,(.2+3*.8)/4)
        with self.assertRaises(ValueError):wp.weighted_coefficients(f,np.zeros(3),c,['R_P'])
        with self.assertRaises(ValueError):wp.weighted_coefficients([0,1,3],p,c,['R_P'])

    def test_modal_recovery_both_incident_sides_and_waves(self):
        # Independently construct X/depth polarization and plane-wave phase.
        for side,wave,angle in [('upper','P',60),('upper','SV',25),('lower','P',25),('lower','SV',30)]:
            c=copy.deepcopy(self.c);c['incident'].update(side=side,wave=wave,angle_deg=angle)
            f=np.array([30.,95.,200.]);p=math.sin(math.radians(angle))/c['model'][side]['vp' if wave=='P' else 'vs']
            for action in ('R','T'):
                layer=side if action=='R' else ('upper' if side=='lower' else 'lower')
                ds=wp.depths(c,layer);anchor=wp.boundaries(c)[0 if layer=='upper' else 1]
                direction=(-1 if side=='upper' else 1) if action=='R' else (1 if side=='upper' else -1)
                cols=[]
                for w in ('P','SV'):
                    v=c['model'][layer]['vp' if w=='P' else 'vs'];q=direction*np.conj(np.sqrt(complex(1/v**2-p*p)))
                    pol=np.array([v*p,v*q]) if w=='P' else np.array([-v*q,v*p])
                    base=anchor if abs(q.imag)<1e-12 else ds[0]
                    cols.append((np.exp(-2j*np.pi*f[:,None,None]*q*(ds[None,:,None]-base))*pol).reshape(len(f),-1))
                independent=np.stack(cols,axis=-1);known=np.array([.3+.2j,-.6+.1j])
                obs=independent@known
                G,_=wp.modal_matrix(c,f,ds,[action+'_P',action+'_SV']);inv,_=wp.modal_inverse(G)
                recovered=np.einsum('fij,fj->fi',inv,obs)
                np.testing.assert_allclose(recovered,np.broadcast_to(known,recovered.shape),rtol=1e-11,atol=1e-11)

    def test_paired_input_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            a,b=Path(tmp)/'target',Path(tmp)/'reference'
            for d in (a,b):(d/'DATA').mkdir(parents=True)
            for d in (a,b):(d/'DATA/STATIONS').write_text('S1 DB 300 2035 0 1979\n')
            (a/'DATA/FORCESOLUTION').write_text('source A')
            (b/'DATA/FORCESOLUTION').write_text('source B')
            with self.assertRaises(ValueError):post.paired_inputs(self.c,a,b)

    def test_end_to_end_packet_and_holdout_detection(self):
        c=copy.deepcopy(self.c);s=c['stations'];pr=c['processing']
        s.update(width_m=5.,step_x_m=2.5,diagnostic_y_offsets_m=[])
        pr.update(widths_m=[5.],primary_width_m=5.,padding_cycles=[6.],primary_padding_cycles=6.,tapers=[0.],primary_taper=0.,nfft=4096)
        rows=generator.station_rows(c);t=np.arange(4096)*.00025;f=np.fft.rfftfreq(len(t),t[1]-t[0])
        top,bottom=wp.boundaries(c);cx=s['center_x_m'];p=wp.horizontal_slowness(c)
        hit=wp.arrival(c,'I_P',cx,top-1e-7)
        pulse=np.exp(-.5*((t-hit)/.010)**2)*np.cos(2*np.pi*95*(t-hit))
        incident=np.fft.rfft(pulse);incident[[0,-1]]=0
        coeff=wp.theory_response(c,f)
        def field(layer,wave,direction,depth,anchor,x,amp):
            v=c['model'][layer]['vp' if wave=='P' else 'vs']
            q=direction*np.conj(np.sqrt(complex(1/v**2-p*p)))
            e=np.array([v*p,v*q]) if wave=='P' else np.array([-v*q,v*p])
            return amp[:,None]*np.exp(-2j*np.pi*f[:,None]*(p*(x-cx)+q*(depth-anchor)))*e[None,:]
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);target=base/'target';reference=base/'reference'
            text=''.join(f"{r['station']} DB {r['y_m']} {r['x_m']} 0 {r['burial_m']}\n" for r in rows)
            for folder in (target,reference):
                (folder/'DATA').mkdir(parents=True);(folder/'OUTPUT_FILES').mkdir()
                (folder/'DATA/STATIONS').write_text(text)
                (folder/'DATA/FORCESOLUTION').write_text('identical synthetic source')
                (folder/'DATA/Par_file').write_text('identical synthetic configuration')
            for r in rows:
                x,d=r['x_m'],r['burial_m'];layer=r['layer']
                if layer=='upper':
                    ref=field('upper','P',1,d,top,x,incident)
                    actual=ref+sum(field('upper',w,-1,d,top,x,incident*np.conj(coeff['R_'+w])) for w in ('P','SV'))
                else:
                    ref=field('upper','P',1,d,top,x,incident)
                    actual=sum(field('lower',w,1,d,bottom,x,incident*np.conj(coeff['T_'+w])) for w in ('P','SV'))
                # Only the withheld target row is perturbed. Fitted coefficients must
                # stay correct while its independent prediction residual detects it.
                if r['depth_role']=='holdout':actual*=1.2
                for folder,z in ((target,actual),(reference,ref)):
                    u=np.fft.irfft(z,len(t),axis=0)
                    for component,j,sign in [('FXX',0,1),('FXZ',1,-1)]:
                        np.savetxt(folder/'OUTPUT_FILES'/f"DB.{r['station']}.{component}.semd",np.column_stack((t,sign*u[:,j])),fmt='%.15e')
            result=post.analyze(c,target,reference,base/'out',plot=False)
            for row in result:self.assertLess(abs(row['relative_error_pct']),.15)
            with (base/'out/quality.csv').open(encoding='utf-8-sig') as handle:q=list(csv.DictReader(handle))
            self.assertGreater(np.median([float(r['T_holdout_residual']) for r in q]),.10)
            self.assertLess(np.median([float(r['T_fit_residual']) for r in q]),.005)


if __name__=='__main__':unittest.main()
