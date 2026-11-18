#!/usr/bin/env python3
"""Generate matched three-layer/uniform inputs and fixed multi-depth packet arrays."""
from pathlib import Path
import argparse
import copy
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import wave_packet as wp
source=wp.import_file('finite_source_core',Path(__file__).with_name('finite_source_core.py'))


def station_rows(c):
    st=c['stations'];inner=c['geometry']['research_bounds_m'];outer=c['geometry']['outer_bounds_m']
    top,bottom=wp.boundaries(c);offset=outer['z'][1]-inner['z'][1]
    half=st['width_m']/2
    xs=st['center_x_m']+np.arange(-round(half/st['step_x_m']),round(half/st['step_x_m'])+1)*st['step_x_m']
    sets=[(st['plane_y_m'],xs,'main')]
    for dy in st.get('diagnostic_y_offsets_m',[]):
        if dy==0:raise ValueError('Diagnostic Y offsets must not duplicate the main plane')
        sets.append((st['plane_y_m']+dy,st['center_x_m']+np.array(st['diagnostic_x_offsets_m']),'diagnostic'))
    rows=[];seen=set()
    for y,xx,role in sets:
        for layer in ('upper','lower'):
            for depth_index,d in enumerate(wp.depths(c,layer)):
                for x in xx:
                    if not (inner['x'][0]<x<inner['x'][1] and inner['y'][0]<y<inner['y'][1]
                            and offset<d<offset+inner['z'][1]-inner['z'][0]):
                        raise ValueError('Station outside the interior research region')
                    if top<=d<=bottom:raise ValueError('No stations on interfaces or inside the weak layer')
                    key=(float(x),float(y),float(d))
                    if key in seen:raise ValueError('Duplicate station coordinates')
                    seen.add(key)
                    rows.append(dict(station=f'S{len(rows)+1:06d}',network='DB',x_m=float(x),y_m=float(y),
                                     burial_m=float(d),layer=layer,depth_index=depth_index,role=role,
                                     depth_role='fit' if depth_index<st['fit_depth_count'] else 'holdout'))
    return rows


def write_sources(path,forces,incident):
    with path.open('w',encoding='ascii') as f:
        for r in forces:
            f.write(f"FORCE {r['force_id']:06d}\ntime shift: {r['time_shift_s']:.12g}\nhdurorf0: {incident['frequency_hz']:.12g}\n"
                    f"latorUTM: {r['y_m']:.9f}\nlongorUTM: {r['x_m']:.9f}\ndepth: {r['burial_m']/1000:.12f}\n"
                    f"source time function: 1\nfactor force source: {incident['force_factor']:.12e}\n"
                    f"component dir vect source E: {r['E']:.12e}\ncomponent dir vect source N: {r['N']:.12e}\n"
                    f"component dir vect source Z_UP: {r['Z_UP']:.12e}\n")


REQUIRED = '''# Merge into your existing Par_file; this is NOT a complete solver Par_file.
# Both jobs must use the identical mesh, source, stations, dt, sampling and boundaries.
SIMULATION_TYPE = 1
ATTENUATION = .false.
ANISOTROPY = .false.
USE_FORCE_POINT_SOURCE = .true.
USE_RICKER_TIME_FUNCTION = .true.
USE_EXTERNAL_SOURCE_FILE = .false.
SUPPRESS_UTM_PROJECTION = .true.
USE_SOURCES_RECEIVERS_Z = .false.
SAVE_SEISMOGRAMS_DISPLACEMENT = .true.
SAVE_SEISMOGRAMS_VELOCITY = .true.
SAVE_SEISMOGRAMS_ACCELERATION = .true.
'''


def run(c,out,overwrite=False,plot=True):
    out=Path(out)
    if (out/'manifest.json').exists() and not overwrite:
        raise FileExistsError('Prepared case exists; use a new --output or --overwrite after reviewing inputs')
    rows=station_rows(c);m=c['model'];i=c['incident'];inner=c['geometry']['research_bounds_m']
    a=np.deg2rad(i['angle_deg']);sgn=1 if i['side']=='upper' else -1
    normal=np.array([np.sin(a),sgn*np.cos(a)])
    ref=np.array([i['reference_x_m'],i['reference_depth_m'][i['side']]])
    depth0=c['geometry']['outer_bounds_m']['z'][1]-inner['z'][1]
    research=(inner['x'][0],inner['x'][1],depth0,depth0+inner['z'][1]-inner['z'][0])
    line,centers,forces=source.finite_sources(m,i,inner,research,normal,ref)
    audit=source.source_boundary_audit(forces,m,i)
    out.mkdir(parents=True,exist_ok=True)
    uniform=copy.deepcopy(m)
    for layer in wp.LAYERS:uniform[layer]=copy.deepcopy(m[i['side']])
    station_text=''.join(f"{s['station']} {s['network']} {s['y_m']:.9f} {s['x_m']:.9f} 0.0 {s['burial_m']:.9f}\n" for s in rows)
    for name,model in [('three_layer',m),('uniform',uniform)]:
        data=out/name/'DATA';mat=out/name/'materials'
        data.mkdir(parents=True,exist_ok=True);mat.mkdir(parents=True,exist_ok=True)
        (data/'STATIONS').write_text(station_text,encoding='ascii')
        write_sources(data/'FORCESOLUTION',forces,i)
        (data/'Par_file.wave_packet_required').write_text(REQUIRED,encoding='ascii')
        source.materials(model,mat)
        wp.write_json(out/name/'model_used.json',model)
    wp.write_csv(out/'station_coordinates.csv',rows)
    wp.write_json(out/'config_used.json',c)
    starts=[dict(x_m=float(x),padding_cycles=float(pad),start_s=wp.window_bounds(c,float(x),pad)[0],
                 end_s=wp.window_bounds(c,float(x),pad)[1])
            for x in sorted({r['x_m'] for r in rows if r['role']=='main'}) for pad in c['processing']['padding_cycles']]
    wp.write_csv(out/'geometry_windows.csv',starts)
    files={name:{f:wp.fingerprint(out/name/'DATA'/f) for f in ('STATIONS','FORCESOLUTION')}
           for name in ('three_layer','uniform')}
    if files['three_layer']!=files['uniform']:raise RuntimeError('Paired input mismatch')
    wp.write_json(out/'manifest.json',dict(schema='wave_packet_prepared_v1',station_count=len(rows),
        main_station_count=sum(r['role']=='main' for r in rows),diagnostic_station_count=sum(r['role']=='diagnostic' for r in rows),
        source_count=len(forces),source_line=line,source_boundary_audit=audit,paired_input_sha256=files,
        coordinate_convention='STATIONS = station network Y X elevation burial_depth; depth positive downward in meters',
        status='Inputs prepared; not a completed SPECFEM run. Reuse identical mesh and replace only material values.',
        geometry_warning='Geometric windows do not guarantee clean packets; inspect source edges, boundaries and late layer reverberations.'))
    if plot:draw(c,rows,out)
    print(f'Prepared {out}: {len(rows)} stations ({sum(r["role"]=="main" for r in rows)} main), {len(forces)} forces')
    return rows


def draw(c,rows,out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a,b=wp.boundaries(c);fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
    ax.axhspan(a,b,color='#d99e23',label='Weak layer')
    for role,color in [('fit','#2678a4'),('holdout','#cf6737')]:
        r=[v for v in rows if v['role']=='main' and v['depth_role']==role]
        ax.scatter([v['x_m'] for v in r],[v['burial_m'] for v in r],s=7,c=color,label=role)
    half=c['processing']['primary_width_m']/2;cx=c['stations']['center_x_m']
    ax.axvspan(cx-half,cx+half,color='#2678a4',alpha=.08,label='Fixed primary segment')
    ax.invert_yaxis();ax.set(xlabel='Global X / m',ylabel='Burial depth / m',title='Wave-packet arrays: fit depths and independent holdout')
    ax.legend();ax.grid(alpha=.15);fig.savefig(out/'station_layout.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=wp.ROOT/'wave_packet_config.json')
    parser.add_argument('--output',type=Path);parser.add_argument('--overwrite',action='store_true')
    args=parser.parse_args();config=wp.load_config(args.config)
    run(config,args.output or wp.ROOT/'result'/config['case_name'],args.overwrite)
