#!/usr/bin/env python3
"""Paired-reference wave-packet RT: X/Z modal separation, displacement L2 ratios."""
from pathlib import Path
import argparse
import csv
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.signal.windows import tukey
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import wave_packet as wp


def data_dir(folder):
    for name in ('DATA','data'):
        if (folder/name/'STATIONS').is_file():return folder/name
    raise FileNotFoundError('Missing DATA/STATIONS under '+str(folder))


def paired_inputs(c,target,reference):
    td,rd=data_dir(target),data_dir(reference);checks={};warnings=[]
    for name in ('STATIONS','FORCESOLUTION','Par_file'):
        a,b=td/name,rd/name
        if name=='Par_file' and not a.exists() and not b.exists():
            warnings.append('Full Par_file unavailable: verify identical solver settings before interpreting results')
            continue
        if not a.is_file() or not b.is_file() or a.read_bytes()!=b.read_bytes():
            raise ValueError('Paired input differs or is missing: '+name)
        checks[name]=wp.fingerprint(a)
    for folder,is_reference in ((target,False),(reference,True)):
        snapshot=folder/'model_used.json'
        if snapshot.exists():
            m=wp.read_json(snapshot)
            for layer in wp.LAYERS:
                expected=c['model'][c['incident']['side'] if is_reference else layer]
                for key in ('rho','vp','vs'):
                    if not np.isclose(m[layer][key],expected[key],rtol=1e-10):
                        raise ValueError('Material metadata differs from configuration: '+str(snapshot))
        else:warnings.append('Material run metadata unavailable: '+str(folder))
    return td,checks,warnings


def select_stations(path,c):
    lookup={}
    for line in path.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):continue
        n,net,y,x,e,d=line.split();key=tuple(round(float(v),7) for v in (x,y,d))
        if key in lookup:raise ValueError('Duplicate station coordinate: '+str(key))
        lookup[key]=(net,n)
    st=c['stations'];cx=st['center_x_m'];count=round(st['width_m']/st['step_x_m'])
    xs=cx-st['width_m']/2+np.arange(count+1)*st['step_x_m']
    locations=[(float(x),st['plane_y_m'],'main') for x in xs]
    locations += [(cx+dx,st['plane_y_m']+dy,'diagnostic')
                  for dy in st.get('diagnostic_y_offsets_m',[]) for dx in st['diagnostic_x_offsets_m']]
    selected={}
    for x,y,role in locations:
        for layer in ('upper','lower'):
            for d in wp.depths(c,layer):
                key=tuple(round(float(v),7) for v in (x,y,d))
                if key not in lookup:raise ValueError('Missing configured station at '+str(key)+'; no silent depth fallback')
                selected[key]=lookup[key]
    return locations,selected


def trace_index(folder,extension,codes):
    index={}
    for path in folder.rglob('*.'+extension):
        fields=path.name.split('.')
        if len(fields)!=4 or fields[2] not in codes.values():continue
        key=tuple(fields[:3])
        if key in index:raise ValueError('Duplicate trace '+str(key))
        index[key]=path
    if not index:raise FileNotFoundError('No '+extension+' X/Z files in '+str(folder))
    return index


def read_trace(path):
    a=np.loadtxt(path,ndmin=2)
    if a.shape[0]<4 or a.shape[1]<2 or not np.isfinite(a[:,:2]).all():raise ValueError('Invalid trace '+str(path))
    t=a[:,0];dt=(t[-1]-t[0])/(len(t)-1)
    if dt<=0 or np.any(np.diff(t)<=0) or np.max(abs(np.diff(t)-dt))>.005*dt:
        raise ValueError('Invalid sample clock '+str(path))
    uniform=np.linspace(t[0],t[-1],len(t))
    return uniform,np.interp(uniform,t,a[:,1])


def read_records(c,locations,stations,target_traces,reference_traces):
    pr=c['processing'];quantity=pr['trace_quantity'];extension={'displacement':'semd','velocity':'semv','acceleration':'sema'}[quantity]
    arrays={};versions={};clock=None
    for name,folder in [('target',target_traces),('reference',reference_traces)]:
        index=trace_index(Path(folder),extension,pr['component_codes']);arrays[name]={}
        layers=('upper','lower') if name=='target' else (c['incident']['side'],)
        for layer in layers:
            records=[]
            for x,y,role in locations:
                per_depth=[]
                for d in wp.depths(c,layer):
                    net,n=stations[tuple(round(float(v),7) for v in (x,y,d))];components=[]
                    for comp in ('X','Z'):
                        key=(net,n,pr['component_codes'][comp])
                        if key not in index:raise FileNotFoundError('Missing trace '+str(key))
                        path=index[key];stat=path.stat();versions[str(path)]=(stat.st_size,stat.st_mtime_ns)
                        t,v=read_trace(path)
                        if clock is None:clock=t
                        if len(clock)!=len(t) or not np.allclose(clock,t,rtol=0,atol=1e-9):raise ValueError('Trace clocks differ')
                        if comp=='Z' and pr['z_positive']=='up':v=-v
                        components.append(v)
                    per_depth.append(np.array(components).T)
                records.append(per_depth)
            values=np.array(records)
            for _ in range({'displacement':0,'velocity':1,'acceleration':2}[quantity]):
                values=cumulative_trapezoid(values,clock,axis=2,initial=0.)
            arrays[name][layer]=values
    return clock,arrays,versions


def make_operators(c,f):
    side=c['incident']['side'];other='lower' if side=='upper' else 'upper'
    nfit=c['stations']['fit_depth_count'];ops={}
    for label,layer,phases in [('I',side,['I_'+c['incident']['wave']]),('R',side,['R_P','R_SV']),('T',other,['T_P','T_SV'])]:
        ds=wp.depths(c,layer);G,anchors=wp.modal_matrix(c,f,ds[:nfit],phases)
        inv,cond=wp.modal_inverse(G,c['processing']['maximum_operator_condition'])
        H=wp.modal_matrix(c,f,ds[nfit:],phases,anchors)[0] if len(ds)>nfit else None
        ops[label]=dict(G=G,H=H,inverse=inv,condition=cond,phases=phases,layer=layer)
    return ops


def spectrum_fit(obs,op,nfit):
    # obs shape depth, frequency, component; training depths are fixed in config.
    y=obs[:nfit].transpose(1,0,2).reshape(obs.shape[1],2*nfit)
    z=np.einsum('fij,fj->fi',op['inverse'],y)
    pred=np.einsum('fij,fj->fi',op['G'],z)
    held=None
    if op['H'] is not None:
        yh=obs[nfit:].transpose(1,0,2).reshape(obs.shape[1],-1)
        ph=np.einsum('fij,fj->fi',op['H'],z)
        held=(yh,ph)
    return z,(y,pred),held


def relative_residual(pair,weights):
    if pair is None:return None
    y,p=pair;den=float(np.sum(weights[:,None]*abs(y)**2))
    return float(np.sqrt(np.sum(weights[:,None]*abs(y-p)**2)/den)) if den>0 else None


def analyze(c,target,reference,out,target_traces=None,reference_traces=None,plot=True):
    target,reference,out=Path(target),Path(reference),Path(out)
    td,paired,warnings=paired_inputs(c,target,reference)
    locations,stations=select_stations(td/'STATIONS',c)
    t,arrays,versions=read_records(c,locations,stations,target_traces or target/'OUTPUT_FILES',reference_traces or reference/'OUTPUT_FILES')
    pr=c['processing'];nfft=pr['nfft'];dt=float(t[1]-t[0]);nfit=c['stations']['fit_depth_count']
    if nfft<len(t):raise ValueError('nfft must not truncate the time record')
    f=np.fft.rfftfreq(nfft,dt);weights=np.full(len(f),2.);weights[[0,-1]]=0.
    ops=make_operators(c,f);theory=wp.theory_response(c,f);phases=wp.output_phases(c)
    side=c['incident']['side'];other='lower' if side=='upper' else 'upper';incident_phase='I_'+c['incident']['wave']
    powers_sum={p:np.zeros(len(f)) for p in [incident_phase]+phases}
    integrals=[];quality=[];windows=[];center={}
    primary_pad=pr['primary_padding_cycles'];primary_taper=pr['primary_taper'];cx=c['stations']['center_x_m']
    for pad in pr['padding_cycles']:
        for taper in pr['tapers']:
            print(f'Packet integration: padding={pad:g}, taper={taper:g}',flush=True)
            for ix,(x,y,role) in enumerate(locations):
                if role=='diagnostic' and (pad!=primary_pad or taper!=primary_taper):continue
                start,end=wp.window_bounds(c,x,pad)
                if start<t[0] or end>t[-1]:raise ValueError('Requested wave-packet window exceeds recorded time')
                mask=(t>=start)&(t<=end)
                if mask.sum()<pr['minimum_samples']:raise ValueError('Too few samples in wave-packet window')
                win=np.zeros(len(t));win[mask]=tukey(mask.sum(),taper)
                B=np.fft.rfft(arrays['reference'][side][ix]*win[None,:,None],nfft,axis=1)
                A=np.fft.rfft(arrays['target'][side][ix]*win[None,:,None],nfft,axis=1)-B
                T=np.fft.rfft(arrays['target'][other][ix]*win[None,:,None],nfft,axis=1)
                modes={};qrow=dict(padding_cycles=pad,taper_fraction=taper,x_m=x,y_m=y,role=role)
                for label,obs in [('I',B),('R',A),('T',T)]:
                    z,fit,hold=spectrum_fit(obs,ops[label],nfit)
                    for j,phase in enumerate(ops[label]['phases']):modes[phase]=z[:,j]
                    qrow[label+'_fit_residual']=relative_residual(fit,weights)
                    qrow[label+'_holdout_residual']=relative_residual(hold,weights)
                power={p:weights*abs(modes[p])**2 for p in [incident_phase]+phases}
                ein=float(dt/nfft*power[incident_phase].sum())
                if ein<=0:raise ValueError('Zero incident wave-packet integral')
                qrow['omitted_incident_DC_fraction']=float(abs(modes[incident_phase][0])**2/(power[incident_phase].sum()+abs(modes[incident_phase][0])**2))
                raw=arrays['reference'][side][ix];ids=np.flatnonzero(mask);nedge=max(1,len(ids)//10)
                edges=np.r_[ids[:nedge],ids[-nedge:]]
                qrow['reference_edge_fraction']=float(np.sum(raw[:,edges,:]**2)/max(np.sum(raw[:,mask,:]**2),1e-300))
                before=t<start-1/c['incident']['frequency_hz']
                qrow['prearrival_rms_ratio']=float(np.sqrt(np.mean(raw[:,before,:]**2)/max(np.mean(raw[:,mask,:]**2),1e-300))) if before.any() else None
                quality.append(qrow)
                windows.append(dict(padding_cycles=pad,taper_fraction=taper,x_m=x,y_m=y,role=role,start_s=start,end_s=end,samples=int(mask.sum())))
                for phase in phases:
                    measured=float(dt/nfft*power[phase].sum())
                    predicted=float(dt/nfft*np.sum(power[incident_phase]*abs(theory[phase])**2))
                    integrals.append(dict(padding_cycles=pad,taper_fraction=taper,x_m=x,y_m=y,role=role,phase=phase,
                        incident_integral=ein,outgoing_integral=measured,predicted_integral=predicted))
                is_primary=pad==primary_pad and taper==primary_taper and role=='main' and abs(x-cx)<=pr['primary_width_m']/2+1e-8
                if is_primary:
                    for phase in powers_sum:powers_sum[phase]+=power[phase]
                if pad==primary_pad and taper==primary_taper and role=='main' and abs(x-cx)<1e-8:
                    for phase in powers_sum:
                        z=modes[phase].copy();z[[0,-1]]=0;center[phase]=np.fft.irfft(z,nfft)
                        e=dt*np.sum(center[phase]**2)
                        if not np.isclose(e,dt/nfft*power[phase].sum(),rtol=1e-10):raise RuntimeError('Parseval check failed')
    summary=[]
    for pad in pr['padding_cycles']:
        for taper in pr['tapers']:
            for width in pr['widths_m']:
                for phase in phases:
                    rows=[r for r in integrals if r['role']=='main' and r['padding_cycles']==pad and r['taper_fraction']==taper and r['phase']==phase and abs(r['x_m']-cx)<=width/2+1e-8]
                    if len(rows)<3:raise ValueError('Fewer than three positions in fixed statistical segment')
                    denominator=sum(r['incident_integral'] for r in rows)
                    value=np.sqrt(sum(r['outgoing_integral'] for r in rows)/denominator)
                    expected=np.sqrt(sum(r['predicted_integral'] for r in rows)/denominator)
                    summary.append(dict(padding_cycles=pad,taper_fraction=taper,width_m=width,phase=phase,positions=len(rows),
                        packet_coefficient=float(value),weighted_theory=float(expected),absolute_error=float(value-expected),
                        relative_error_pct=float(100*(value/expected-1)) if expected>1e-12 else None))
    primary=[r for r in summary if r['padding_cycles']==primary_pad and r['taper_fraction']==primary_taper and r['width_m']==pr['primary_width_m']]
    out.mkdir(parents=True,exist_ok=True)
    wp.write_csv(out/'RT_packet_primary.csv',primary);wp.write_csv(out/'RT_packet_sensitivity.csv',summary)
    wp.write_csv(out/'local_packet_integrals.csv',integrals);wp.write_csv(out/'quality.csv',quality);wp.write_csv(out/'station_windows.csv',windows)
    wp.write_csv(out/'reference_spectrum.csv',[dict(frequency_hz=float(freq),incident_power=float(powers_sum[incident_phase][k])) for k,freq in enumerate(f)])
    wp.write_csv(out/'frequency_response.csv',[dict(frequency_hz=float(freq),**{p+'_power':float(v[k]) for p,v in powers_sum.items()},
        **{p+'_theory_abs':float(abs(v[k])) for p,v in theory.items()}) for k,freq in enumerate(f)])
    lateral=[]
    for row in integrals:
        if row['role']!='diagnostic':continue
        central=next(r for r in integrals if r['role']=='main' and r['padding_cycles']==primary_pad and r['taper_fraction']==primary_taper and r['x_m']==row['x_m'] and r['phase']==row['phase'])
        lateral.append(dict(x_m=row['x_m'],y_m=row['y_m'],phase=row['phase'],
            incident_L2_ratio_to_center_y=float(np.sqrt(row['incident_integral']/central['incident_integral'])),
            packet_coefficient=float(np.sqrt(row['outgoing_integral']/row['incident_integral']))))
    wp.write_csv(out/'lateral_diagnostics.csv',lateral,['x_m','y_m','phase','incident_L2_ratio_to_center_y','packet_coefficient'])
    cdf=np.cumsum(powers_sum[incident_phase]);cdf/=cdf[-1]
    percentiles={f'f{q:g}_hz':float(np.interp(q/100,cdf,f)) for q in (.5,5,50,95,99.5)}
    for path,version in versions.items():
        stat=Path(path).stat()
        if (stat.st_size,stat.st_mtime_ns)!=version:raise RuntimeError('Input changed during analysis: '+path)
    wp.write_json(out/'config_used.json',c)
    wp.write_json(out/'audit.json',dict(schema='wave_packet_analysis_v1',quantity='displacement',primary=primary,
        incident_spectral_percentiles=percentiles,input_count=len(versions),input_bytes=sum(v[0] for v in versions.values()),
        paired_input_sha256=paired,warnings=warnings,fit_depth_count=nfit,
        holdout_depth_counts={k:len(wp.depths(c,k))-nfit for k in ('upper','lower')},
        maximum_operator_conditions={label:float(op['condition'].max()) for label,op in ops.items()},
        method='XZ only; incident-side target minus uniform for reflection; transmitted target unchanged; displacement L2 pooling; measured-reference spectrum weighting of full three-layer C(f).',
        scope='Given input wave packet and model; these are not single-frequency coefficients or physical energy-flux ratios.'))
    np.savez_compressed(out/'center_modal_packets.npz',t=t[0]+np.arange(nfft)*dt,**center)
    lines=['# 三层模型波包反透射验证','','主指标：位移平方积分均值之比开方；仅 X/Z，使用配对均匀参考场。',
           '理论为实测入射位移谱平方加权的完整三层响应；不等同于单频系数或物理能流比。','',
           '|波型|数值波包系数|频谱加权理论|相对偏差|','|---|---:|---:|---:|']
    for r in primary:
        error='理论近零，见绝对误差' if r['relative_error_pct'] is None else f"{r['relative_error_pct']:+.3f}%"
        lines.append(f"|{r['phase']}|{r['packet_coefficient']:.6f}|{r['weighted_theory']:.6f}|{error}|")
    lines+=['',f"入射位移平方积分 5%～95% 频率范围：{percentiles['f5_hz']:.2f}～{percentiles['f95_hz']:.2f} Hz。",
            '时间窗与台站范围由配置固定，全部敏感性组合均输出，不按理论误差选站。',
            'quality.csv 包含固定拟合深度与留出深度预测残差；没有留出深度时对应字段为空。',
            '窗口边缘贡献、前到时偏移、边界/源端部污染和侧向差异仍需核查；积分吻合不自动证明时程或逐频响应完全一致。']
    lines+=['']+warnings;(out/'report.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    if plot:draw(out,primary,t[0]+np.arange(nfft)*dt,center,wp.window_bounds(c,cx,primary_pad))
    print('Packet results:',out,flush=True)
    for row in primary:print(row,flush=True)
    return primary


def draw(out,rows,t,center,bounds):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained');x=np.arange(len(rows));w=.34
    ax.bar(x-w/2,[r['packet_coefficient'] for r in rows],w,label='Numerical packet')
    ax.bar(x+w/2,[r['weighted_theory'] for r in rows],w,label='Spectrum-weighted theory')
    ax.set(xticks=x,xticklabels=[r['phase'] for r in rows],ylabel='Displacement L2 amplitude ratio',title='Full three-layer wave-packet validation');ax.legend()
    fig.savefig(out/'packet_coefficients.png',dpi=170);plt.close(fig)
    fig,axs=plt.subplots(len(center),1,figsize=(9,2*len(center)),sharex=True,layout='constrained')
    for ax,(phase,u) in zip(np.atleast_1d(axs),center.items()):ax.plot(t,u);ax.set(ylabel=phase,xlim=bounds);ax.grid(alpha=.15)
    np.atleast_1d(axs)[-1].set_xlabel('Time / s');fig.savefig(out/'center_modal_packets.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',type=Path,default=wp.ROOT/'wave_packet_config.json')
    p.add_argument('--case-dir',type=Path);p.add_argument('--target',type=Path);p.add_argument('--reference',type=Path)
    p.add_argument('--target-traces',type=Path);p.add_argument('--reference-traces',type=Path);p.add_argument('--output',type=Path)
    a=p.parse_args();c=wp.load_config(a.config);case=a.case_dir or wp.ROOT/'result'/c['case_name']
    analyze(c,a.target or case/'three_layer',a.reference or case/'uniform',a.output or case/'packet_results',a.target_traces,a.reference_traces)
