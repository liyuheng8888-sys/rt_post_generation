"""Shared configuration, modal geometry and spectral weighting for wave-packet RT."""
from pathlib import Path
import copy
import csv
import hashlib
import importlib.util
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT / '含软弱层的三层算例'
LAYERS = ('upper', 'middle', 'lower')


def import_file(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def solver():
    return import_file('packet_theory_core', PROJECT / '理论解/rt_coefficients_theorical_solution.py')


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def write_json(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


def write_csv(path, rows, fields=None):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fingerprint(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_config(path=None):
    c = read_json(path or ROOT/'wave_packet_config.json')
    if c.get('schema') != 'wave_packet_v1':
        raise ValueError('Expected wave_packet_v1 configuration')
    if not c['case_name'] or c['case_name'] in ('.', '..') or any(x in c['case_name'] for x in '/\\:'):
        raise ValueError('case_name must be a single directory name')
    solver().validate(c['model'], c['incident'])
    outer=c['geometry']['outer_bounds_m'];inner=c['geometry']['research_bounds_m']
    for axis in ('x','y','z'):
        lo,hi=outer[axis];a,b=inner[axis]
        if not np.isfinite([lo,hi,a,b]).all() or lo!=0 or not lo<=a<b<=hi:
            raise ValueError('Invalid outer/research bounds: '+axis)
    model=c['model']
    if not np.allclose([model['size_x_m'],model['size_y_m'],sum(model[k+'_thickness_m'] for k in LAYERS)],
                       [outer[k][1] for k in ('x','y','z')],rtol=0,atol=1e-7):
        raise ValueError('Global model dimensions/thicknesses must match outer bounds')
    a,b=boundaries(c);ref=c['incident'];rd=ref['reference_depth_m'][ref['side']]
    if not (0<ref['reference_x_m']<model['size_x_m'] and 0<rd<outer['z'][1]
            and (rd<a if ref['side']=='upper' else rd>b)):
        raise ValueError('Source reference must be inside the incident medium')
    for layer in LAYERS:
        m = c['model'][layer]
        if m['type'] != 'solid' or m.get('poroelastic_enabled', False):
            raise ValueError('Wave-packet workflow currently supports three lossless isotropic elastic solids')
    st, pr = c['stations'], c['processing']
    if st['step_x_m'] <= 0 or st['width_m'] <= 0:
        raise ValueError('Station width and spacing must be positive')
    intervals = st['width_m']/st['step_x_m']
    if not np.isclose(intervals, round(intervals)) or round(intervals) % 2:
        raise ValueError('Station width/spacing must be an even integer, including a center station')
    for key in ('upper_offsets_m', 'lower_offsets_m'):
        a = np.array(st[key], dtype=float)
        if len(a) < 2 or not np.isfinite(a).all() or np.any(a <= 0) or np.any(np.diff(a) <= 0):
            raise ValueError(key+' must contain at least two increasing positive distances')
        if not isinstance(st['fit_depth_count'], int) or not 2 <= st['fit_depth_count'] <= len(a):
            raise ValueError('fit_depth_count must be between 2 and the number of depths')
    if pr['trace_quantity'] not in ('displacement', 'velocity', 'acceleration'):
        raise ValueError('Invalid trace_quantity')
    if set(pr['component_codes']) != {'X', 'Z'} or pr['z_positive'] not in ('up', 'down'):
        raise ValueError('Only X/Z components are supported; specify the Z sign')
    for key, primary, minimum, maximum in [('widths_m','primary_width_m',0,st['width_m']),
                                          ('padding_cycles','primary_padding_cycles',0,np.inf),
                                          ('tapers','primary_taper',-1e-15,1)]:
        a = np.asarray(pr[key], float)
        if not len(a) or not np.isfinite(a).all() or np.any(a <= minimum) or np.any(a > maximum):
            raise ValueError('Invalid '+key)
        if pr[primary] not in pr[key]:
            raise ValueError(primary+' must be included in '+key)
    if not isinstance(pr['nfft'], int) or pr['nfft'] < 4 or pr['nfft'] % 2:
        raise ValueError('nfft must be an even integer >=4')
    scan=c['theory']
    if not 0<scan['spectrum_min_hz']<=scan['spectrum_max_hz'] or scan['spectrum_step_hz']<=0:
        raise ValueError('Theory frequency scan must be positive and ordered')
    if pr['minimum_samples']<4 or pr['maximum_operator_condition']<=1:
        raise ValueError('Invalid sample/conditioning limits')
    # Catch NaN in any numeric setting before coordinates or output are produced.
    json.dumps(c, allow_nan=False)
    return c


def boundaries(c):
    a = c['model']['upper_thickness_m']
    return a, a+c['model']['middle_thickness_m']


def horizontal_slowness(c):
    i = c['incident']
    return math.sin(math.radians(i['angle_deg'])) / c['model'][i['side']]['vp' if i['wave']=='P' else 'vs']


def depths(c, layer):
    a,b = boundaries(c)
    return np.array([a-d for d in c['stations']['upper_offsets_m']] if layer=='upper'
                    else [b+d for d in c['stations']['lower_offsets_m']])


def modal_info(c, phase):
    action,wave = phase.split('_')
    side = c['incident']['side']
    layer = side if action in ('I','R') else ('lower' if side=='upper' else 'upper')
    direction = (1 if side=='upper' else -1) * (-1 if action=='R' else 1)
    v = c['model'][layer]['vp' if wave=='P' else 'vs']
    q = direction*np.conj(np.sqrt(complex(1/v**2-horizontal_slowness(c)**2)))
    return layer,v,q,abs(q.imag)<1e-12 and abs(q.real)>1e-12


def output_phases(c):
    return [p for p in ('R_P','R_SV','T_P','T_SV') if modal_info(c,p)[3]]


def modal_matrix(c, f, ds, phases, anchors=None):
    """Rows are depth-major (X, depth-positive); columns are displacement modes."""
    f=np.asarray(f);ds=np.asarray(ds);p=horizontal_slowness(c)
    a,b=boundaries(c);columns=[];used=[]
    for j,phase in enumerate(phases):
        layer,v,q,propagating=modal_info(c,phase)
        anchor=(a if layer=='upper' else b) if propagating else float(ds[np.argmin(abs(ds-(a if layer=='upper' else b)))])
        if anchors is not None:anchor=anchors[j]
        used.append(anchor)
        pol=np.array([v*p,v*q]) if phase.endswith('_P') else np.array([-v*q,v*p])
        col=np.exp(-2j*np.pi*f[:,None,None]*q*(ds[None,:,None]-anchor))*pol[None,None,:]
        columns.append(col.reshape(len(f),2*len(ds)))
    return np.stack(columns,axis=-1),used


def modal_inverse(matrix, maximum_condition=1e8):
    scale=np.linalg.norm(matrix,axis=1)
    if np.any(scale < 1e-200):raise ValueError('Vanishing modal column; check depths and evanescent anchors')
    normal=matrix/scale[:,None,:]
    singular=np.linalg.svd(normal,compute_uv=False)
    condition=singular[:,0]/singular[:,-1]
    if not np.isfinite(condition).all() or np.max(condition)>maximum_condition:
        raise ValueError('Modal separation is ill-conditioned at this angle/depth configuration')
    return np.linalg.pinv(normal,rcond=1/maximum_condition)/scale[:,:,None],condition


def arrival(c,phase,x,d):
    layer,v,q,propagating=modal_info(c,phase)
    if not propagating:return None
    a,b=boundaries(c)
    if not (d<a if layer=='upper' else d>b):return None
    inc=c['incident'];side=inc['side'];p=horizontal_slowness(c)
    _,vi,qi,_=modal_info(c,'I_'+inc['wave'])
    origin=inc['peak_time_s']+p*(x-inc['reference_x_m'])
    rd=inc['reference_depth_m'][side]
    if phase.startswith('I_'):return float(origin+qi.real*(d-rd))
    near=a if side=='upper' else b
    anchor=a if layer=='upper' else b
    return float(origin+qi.real*(near-rd)+q.real*(d-anchor))


def window_bounds(c,x,padding_cycles):
    phases=['I_'+c['incident']['wave']]+output_phases(c)
    arrivals=[arrival(c,p,x,d) for p in phases for layer in ('upper','lower') for d in depths(c,layer)]
    arrivals=[t for t in arrivals if t is not None]
    pad=padding_cycles/c['incident']['frequency_hz']
    return min(arrivals)-pad,max(arrivals)+pad


def theory_response(c,frequencies):
    th=solver();f=np.asarray(frequencies,float)
    if np.any(f<0) or not np.isfinite(f).all():raise ValueError('Invalid frequencies')
    values={p:np.zeros(len(f),complex) for p in ('R_P','R_SV','T_P','T_SV')}
    for k,freq in enumerate(f):
        if freq==0:continue
        if len(f)>2000 and k%2000==0:print(f'Theory frequencies: {k}/{len(f)}',flush=True)
        coeff=th.solve(c['model'],c['incident'],frequency_hz=float(freq))['coefficients']
        for phase in values:values[phase][k]=coeff[phase]
    return values


def weighted_coefficients(f,incident_power,coefficients,phases):
    f=np.asarray(f,float);power=np.asarray(incident_power,float)
    if len(f)<2 or len(f)!=len(power) or np.any(np.diff(f)<=0) or not np.isfinite(power).all() or np.any(power<0):
        raise ValueError('Invalid reference spectrum')
    # input powers are per equally spaced FFT bin, with one-sided factors included.
    if not np.allclose(np.diff(f),np.diff(f)[0],rtol=1e-7,atol=1e-10):
        raise ValueError('Reference spectrum must use equally spaced FFT bins')
    if power.sum()<=0:raise ValueError('Reference spectrum has zero power')
    return {p:float(np.sqrt(np.sum(power*abs(coefficients[p])**2)/power.sum())) for p in phases}


def theory_cli():
    import argparse
    parser=argparse.ArgumentParser(description='Full three-layer C(f) and reference-spectrum-weighted packet coefficients')
    parser.add_argument('--config',type=Path,default=ROOT/'wave_packet_config.json')
    parser.add_argument('--reference-spectrum',type=Path,help='reference_spectrum.csv exported by packet postprocessing')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();c=load_config(args.config)
    out=args.output or ROOT/'result'/c['case_name']/'theory';out.mkdir(parents=True,exist_ok=True)
    power=None
    if args.reference_spectrum:
        with args.reference_spectrum.open(encoding='utf-8-sig') as handle:rows=list(csv.DictReader(handle))
        f=np.array([float(r['frequency_hz']) for r in rows]);power=np.array([float(r['incident_power']) for r in rows])
    else:
        s=c['theory'];f=np.arange(s['spectrum_min_hz'],s['spectrum_max_hz']+s['spectrum_step_hz']/2,s['spectrum_step_hz'])
    coeff=theory_response(c,f)
    rows=[dict(frequency_hz=float(freq),**{p+'_'+suffix:float(value) for p,z in coeff.items()
          for suffix,value in [('real',z[k].real),('imag',z[k].imag),('abs',abs(z[k]))]}) for k,freq in enumerate(f)]
    write_csv(out/'coefficients_frequency_spectrum.csv',rows)
    if power is not None:
        result=weighted_coefficients(f,power,coeff,output_phases(c))
        write_csv(out/'packet_weighted_theory.csv',[dict(phase=p,packet_coefficient=v) for p,v in result.items()])
    write_json(out/'model_parameters.json',c)
    print('Theory output:',out)
