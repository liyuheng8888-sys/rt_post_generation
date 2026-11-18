#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""三层频域求解核；命令行统一读取仓库根目录 wave_packet_config.json。"""
from pathlib import Path
import json
import csv
import argparse
import copy
import math
import numpy as np
from scipy.linalg import eig
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 以下默认值仅保留旧 solve/run Python API 兼容性。
# 波包分支的三个命令行入口统一使用 wave_packet_config.json，不修改此处。
MODEL = {
    "size_x_m": 1000.0,
    "size_y_m": 1000.0,
    "upper_thickness_m": 500.0,
    "middle_thickness_m": 1.0,
    "lower_thickness_m": 499.0,
    "upper": {"name": "上层固体", "type": "solid", "rho": 2000.0, "vp": 3000.0, "vs": 1800},
    "middle": {
        "name": "中间软弱层", "type": "solid",  # solid普通弹性 / fluid流体
        "rho": 1000.0, "vp": 1450.0, "vs": 850.0,
        "poroelastic_enabled": False,  # True: 使用下方Biot参数，忽略本层rho/vp/vs
        "biot": {
            "rho_s": 2650.0, "rho_f": 1000.0,  # 颗粒/孔隙流体密度 kg/m³
            "porosity": 0.35, "tortuosity": 2.0,
            "permeability_m2": 1.0e-10,  # 各向同性渗透率 m²
            "viscosity_Pa_s": 1.0e-3,
            "grain_bulk_Pa": 36.0e9, "fluid_bulk_Pa": 2.25e9,
            "frame_bulk_Pa": 1.0e8, "frame_shear_Pa": 2.0e7,
        },
    },
    "lower": {"name": "下层固体", "type": "solid", "rho": 2500.0, "vp": 4500.0, "vs": 2100.0},
    # 流体/多孔界面：open连续孔压与总体积流量；sealed孔隙相对流量为0。
    # 普通固体/多孔界面使用粘结、sealed。不支持部分饱和或JKD高频修正。
    "upper_middle_pore_condition": "open",
    "middle_lower_pore_condition": "sealed",
    # 材料号可改；不自动覆盖已有网格。请与materials_file逐项对应。
    "material_ids": {"upper": 1, "middle": 4, "lower": 3},
}

INCIDENT = {'side': 'lower', 'wave': 'SV', 'angle_deg': 15.0, 'frequency_hz': 95.0}

THEORY = {
    "angle_step_deg": 0.5, "max_angle_deg": 89.0,
    "frequencies_hz": [20.0, 50.0, 100.0, 150.0, 200.0],
    "spectrum_min_hz": 1.0, "spectrum_max_hz": 300.0, "spectrum_step_hz": 1.0,
}

# ==================== 用户编辑区结束 ====================

def kind(m):
    return 'poroelastic' if m.get('poroelastic_enabled',False) else m['type']


def validate(model,incident):
    for key in ('size_x_m','size_y_m','upper_thickness_m','middle_thickness_m','lower_thickness_m'):
        if not np.isfinite(model[key]) or model[key]<=0:raise ValueError(key+' must be finite and positive')
    ids=list(model['material_ids'].values())
    if len(ids)!=3 or len(set(ids))!=3 or any(type(x) is not int or x<=0 for x in ids):raise ValueError('Three distinct positive material IDs required')
    for layer in ('upper','middle','lower'):
        m=model[layer];k=kind(m)
        if k=='poroelastic':
            if layer!='middle':raise ValueError('Biot is supported for the finite middle layer only')
            b=m['biot']
            for key,val in b.items():
                if not np.isfinite(val) or val<0:raise ValueError('Invalid Biot '+key)
            for key in ('rho_s','rho_f','permeability_m2','grain_bulk_Pa','fluid_bulk_Pa','frame_bulk_Pa','frame_shear_Pa'):
                if b[key]<=0:raise ValueError('Biot '+key+' must be >0')
            if not 0<b['porosity']<1 or b['tortuosity']<1:raise ValueError('Require 0<porosity<1 and tortuosity>=1')
            alpha=1-b['frame_bulk_Pa']/b['grain_bulk_Pa']
            if not 0<alpha<=1:raise ValueError('Frame bulk modulus must be below grain bulk modulus')
            biot_parameters(m,1.0)
        elif k in ('fluid','solid'):
            if not all(np.isfinite(m[x]) for x in ('rho','vp','vs')) or m['rho']<=0 or m['vp']<=0:raise ValueError('Invalid '+layer+' density/speed')
            if k=='fluid' and m['vs']!=0:raise ValueError('Fluid vs must be zero')
            if k=='solid' and not (m['vs']>0 and m['vp']**2>4*m['vs']**2/3):raise ValueError('Solid requires positive shear and bulk moduli')
        else:raise ValueError('Unknown material type: '+k)
    if incident['side'] not in ('upper','lower'):raise ValueError('Incident side must be upper/lower')
    if incident['wave'] not in names(model[incident['side']]):raise ValueError('Unsupported incident wave in selected medium')
    if not np.isfinite(incident['angle_deg']) or not 0<=incident['angle_deg']<90:raise ValueError('Angle must be in [0,90)')
    if not np.isfinite(incident['frequency_hz']) or incident['frequency_hz']<=0:raise ValueError('Frequency must be >0')
    for key in ('upper_middle_pore_condition','middle_lower_pore_condition'):
        if model[key] not in ('open','sealed'):raise ValueError(key+' must be open/sealed')


def boundaries(model):
    a=model['upper_thickness_m'];b=a+model['middle_thickness_m']
    return a,b,b+model['lower_thickness_m']


def names(m):
    return ['P'] if kind(m)=='fluid' else ['P1','P2','SV'] if kind(m)=='poroelastic' else ['P','SV']


def speed(m,wave):
    return m['vp'] if wave=='P' else m['vs']


def biot_parameters(m,omega):
    b=m['biot'];phi=b['porosity'];rf=b['rho_f']
    rho=(1-phi)*b['rho_s']+phi*rf
    mass=b['tortuosity']*rf/phi
    alpha=1-b['frame_bulk_Pa']/b['grain_bulk_Pa']
    invM=(alpha-phi)/b['grain_bulk_Pa']+phi/b['fluid_bulk_Pa']
    if invM<=0 or rho*mass-rf**2<=0:raise ValueError('Nonpositive Biot compressibility or mass matrix')
    M=1/invM;C=alpha*M;mu=b['frame_shear_Pa']
    H=b['frame_bulk_Pa']+4*mu/3+alpha**2*M
    meff=mass+1j*b['viscosity_Pa_s']/(omega*b['permeability_m2'])
    return rho,rf,alpha,M,C,H,mu,meff


def positive_root(z):
    root=np.sqrt(complex(z))
    if root.imag<-1e-13 or (abs(root.imag)<1e-13 and root.real<0):root=-root
    return root


def modes(m,p,omega,direction):
    k=kind(m);out={}
    if k!='poroelastic':
        for wave in names(m):
            v=speed(m,wave);s=1/v;q=direction*positive_root(s*s-p*p)
            u=np.array([v*p,v*q]) if wave=='P' else np.array([-v*q,v*p])
            ux,uz=u
            if k=='fluid':tx=0;tz=m['rho']*v;pf=-tz
            else:
                mu=m['rho']*m['vs']**2;lam=m['rho']*(m['vp']**2-2*m['vs']**2)
                tx=mu*(q*ux+p*uz);tz=lam*(p*ux+q*uz)+2*mu*q*uz;pf=0
            out[wave]=dict(state=np.array([ux,uz,0,tx,tz,pf],complex),q=q,slowness=s)
        return out
    rho,rf,alpha,M,C,H,mu,meff=biot_parameters(m,omega)
    values,vectors=eig(np.array([[rho,rf],[rf,meff]],complex),np.array([[H,C],[C,M]],complex))
    order=sorted(range(2),key=lambda i:abs(values[i]))
    entries=[]
    for label,i in zip(('P1','P2'),order):
        ratio=vectors[1,i]/vectors[0,i]
        entries.append((label,values[i],ratio))
    entries.append(('SV',(rho-rf*rf/meff)/mu,-rf/meff))
    for label,s2,ratio in entries:
        s=positive_root(s2);q=direction*positive_root(s2-p*p)
        ux,uz=(p/s,q/s) if label!='SV' else (-q/s,p/s)
        wx,wz=ratio*ux,ratio*uz
        divu=p*ux+q*uz;divw=p*wx+q*wz
        pf=-M*(alpha*divu+divw)
        tx=mu*(q*ux+p*uz)
        tz=(m['biot']['frame_bulk_Pa']-2*mu/3)*divu+2*mu*q*uz-alpha*pf
        out[label]=dict(state=np.array([ux,uz,wz,tx,tz,pf],complex),q=q,slowness=s)
    return out


def interface_rows(a,b,pore):
    """Return L,R with L*state_above + R*state_below = 0."""
    ka,kb=kind(a),kind(b);L=[];R=[]
    def add(l,r):
        left=np.zeros(6);right=np.zeros(6)
        for i,v in l.items():left[i]=v
        for i,v in r.items():right[i]=v
        L.append(left);R.append(right)
    if 'poroelastic' not in (ka,kb):
        if ka==kb=='solid':
            for i in (0,1,3,4):add({i:1},{i:-1})
        else:
            for i in (1,4):add({i:1},{i:-1})
            if ka=='solid':add({3:1},{})
            if kb=='solid':add({},{3:1})
    else:
        if ka==kb:raise ValueError('Only middle layer may be poroelastic')
        # Form ordinary medium on left and porous medium on right, then reorder.
        porous_b=kb=='poroelastic';ordinary=ka if porous_b else kb
        if ordinary=='fluid':
            add({1:1},{1:-1,2:-1} if pore=='open' else {1:-1})
            add({4:1},{4:-1});add({},{3:1})
            add({4:1},{5:1}) if pore=='open' else add({},{2:1})
        else:
            if pore!='sealed':raise ValueError('Solid/poroelastic contact must be sealed (impermeable solid)')
            for i in (0,1,3,4):add({i:1},{i:-1})
            add({},{2:1})
        if not porous_b:L,R=R,L
    return np.array(L,complex),np.array(R,complex)


def flux(state):
    ux,uz,wz,tx,tz,pf=state
    return float(np.real(tx*np.conj(ux)+tz*np.conj(uz)-pf*np.conj(wz)))


def finite_layer_modes(m,p,omega,top,bottom):
    """Use a nonsingular standing-wave pair near elastic/acoustic q=0.

    With eta=q/s and V(+/-q)=E+/-eta*O, the two fields are
    (F+ + F-)/2 and (F+ - F-)/(2*eta), centered inside the layer.
    Their q=0 limits are E and O+i*omega*s*(z-center)*E.
    These are basis amplitudes, not separate up/down traveling amplitudes.
    Strongly evanescent modes retain the bounded interface-anchored basis.
    """
    entries=[]
    for direction in (1,-1):
        for wave,mode in modes(m,p,omega,direction).items():
            entries.append(dict(wave=wave,direction=direction,
                                anchor=top if direction==1 else bottom,**mode))
    if kind(m)=='poroelastic':return entries
    for wave in names(m):
        pair=[entry for entry in entries if entry['wave']==wave]
        plus,minus=pair;q=plus['q'];s=plus['slowness']
        if abs(q/s)>=1e-3 or abs((omega*q*(bottom-top)).imag)>1:continue
        even=(plus['state']+minus['state'])/2
        mu=m['rho']*m['vs']**2
        # O=s*dV/dq is evaluated analytically, including exactly q=0.
        odd=np.array([0,1,0,2*mu*p,0,0] if wave=='P'
                     else [-1,0,0,0,2*mu*p,0],complex)
        for entry,basis,state in zip(pair,('even','odd'),(even,odd)):
            entry.update(direction=0,anchor=top+(bottom-top)/2,q=q,
                         basis=basis,state=state,even_state=even,odd_state=odd)
    return entries


def basis_state(entry,omega,depth):
    """Evaluate a traveling or regularized internal basis field at depth."""
    distance=depth-entry['anchor']
    if entry.get('basis','traveling')=='traveling':
        return entry['state']*np.exp(1j*omega*entry['q']*distance)
    q=entry['q'];s=entry['slowness'];x=omega*q*distance
    even,odd=entry['even_state'],entry['odd_state']
    if entry['basis']=='even':return even*np.cos(x)+1j*(q/s)*odd*np.sin(x)
    # np.sinc(x/pi)=sin(x)/x, with its analytic value 1 at zero.
    return odd*np.cos(x)+1j*omega*s*distance*even*np.sinc(x/np.pi)


def solve(model,incident,angle_deg=None,frequency_hz=None):
    c=copy.deepcopy(incident)
    if angle_deg is not None:c['angle_deg']=float(angle_deg)
    if frequency_hz is not None:c['frequency_hz']=float(frequency_hz)
    validate(model,c)
    side=c['side'];sgn=1 if side=='upper' else -1;omega=2*np.pi*c['frequency_hz']
    p=math.sin(math.radians(c['angle_deg']))/speed(model[side],c['wave'])
    up,mid,low=(model[k] for k in ('upper','middle','lower'))
    dtop,dbottom,total=boundaries(model);h=dbottom-dtop
    unknown=[]
    # Exterior outgoing modes and both directions inside the finite layer.
    for layer,medium,directions in [('upper',up,[-1]),('middle',mid,[1,-1]),('lower',low,[1])]:
        if layer=='middle':
            unknown.extend(dict(layer=layer,**entry) for entry in finite_layer_modes(mid,p,omega,dtop,dbottom))
            continue
        for direction in directions:
            for label,entry in modes(medium,p,omega,direction).items():
                anchor=dtop if layer=='upper' or (layer=='middle' and direction==1) else dbottom
                unknown.append(dict(layer=layer,wave=label,direction=direction,anchor=anchor,**entry))
    incident_mode=modes(model[side],p,omega,sgn)[c['wave']]
    n=len(unknown);mat=[];rhs=[]
    for a,b,depth,contact in [('upper','middle',dtop,model['upper_middle_pore_condition']),('middle','lower',dbottom,model['middle_lower_pore_condition'])]:
        L,R=interface_rows(model[a],model[b],contact);block=np.zeros((len(L),n),complex);known=np.zeros(len(L),complex)
        for j,entry in enumerate(unknown):
            op=L if entry['layer']==a else R if entry['layer']==b else None
            if op is not None:block[:,j]=op@basis_state(entry,omega,depth)
        if side==a:known=L@incident_mode['state']
        if side==b:known=R@incident_mode['state']
        mat.extend(block);rhs.extend(-known)
    A=np.array(mat);b=np.array(rhs)
    if A.shape!=(n,n):raise ValueError('Interface equations do not form square system: '+str(A.shape))
    scale=np.max(abs(A),axis=1);As=A/scale[:,None];bs=b/scale
    col=np.max(abs(As),axis=0);scaled=As/col[None,:]
    solution=np.linalg.solve(scaled,bs)/col
    residual=float(np.linalg.norm(As@solution-bs)/max(np.linalg.norm(bs),1e-30))
    if not np.isfinite(solution).all() or residual>1e-8:raise ValueError('Unreliable boundary solve at angle/frequency '+str((c['angle_deg'],c['frequency_hz'],residual)))
    coeff={};energies={};incflux=abs(flux(incident_mode['state']));internal=[]
    for entry,value in zip(unknown,solution):
        layer=entry['layer'];wave=entry['wave']
        if layer=='middle':
            # Evaluate with basis_state; direction=0 denotes an even/odd pair.
            internal.append(dict(entry,coefficient=value))
            continue
        tag=('R_' if layer==side else 'T_')+wave
        coeff[tag]=value
        energies[tag]=abs(value)**2*abs(flux(entry['state']))/incflux
    er=sum(v for k,v in energies.items() if k.startswith('R_'));et=sum(v for k,v in energies.items() if k.startswith('T_'))
    # Coherent net flux in porous layer accounts for interference and dissipation.
    internal_states=[]
    for depth in (dtop,dbottom):
        state=np.zeros(6,complex)
        for entry,value in zip(unknown,solution):
            if entry['layer']=='middle':state+=value*basis_state(entry,omega,depth)
        internal_states.append(state)
    layer_loss=(flux(internal_states[0])-flux(internal_states[1]))/incflux
    return dict(coefficients=coeff,energy=energies,R_energy=er,T_energy=et,absorption=1-er-et,
                middle_flux_loss=layer_loss,residual=residual,condition_scaled=float(np.linalg.cond(scaled)),
                angle_deg=c['angle_deg'],frequency_hz=c['frequency_hz'],p=p,internal=internal)


def flatten(result):
    row={k:result[k] for k in ('angle_deg','frequency_hz','p','R_energy','T_energy','absorption','middle_flux_loss','residual','condition_scaled')}
    for name,value in result['coefficients'].items():
        for suffix,x in [('real',value.real),('imag',value.imag),('abs',abs(value))]:row[name+'_'+suffix]=float(x)
        row[name+'_energy']=result['energy'][name]
    return row


def write_table(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)


def write_weak_layer_comparison(path, rows, base_rows, phases, reference_key,
                                x_key, xlabel, title, middle_thickness):
    comparison = []
    other_key = 'frequency_hz' if x_key == 'angle_deg' else 'angle_deg'
    for actual, base in zip(rows, base_rows):
        item = {key: actual[key] for key in (x_key, other_key)}
        for phase in phases:
            for suffix in ('abs', 'real', 'imag'):
                item[f'{phase}_weak_{suffix}'] = actual[f'{phase}_{suffix}']
                item[f'{phase}_{reference_key}_{suffix}'] = base[f'{phase}_{suffix}']
            item[f'{phase}_delta_abs'] = actual[f'{phase}_abs'] - base[f'{phase}_abs']
        comparison.append(item)

    nrows = (len(phases) + 1) // 2
    fig, axes = plt.subplots(nrows, 2, figsize=(12, 3.7 * nrows), layout='constrained')
    x = [r[x_key] for r in comparison]
    for ax, phase in zip(axes.flat, phases):
        branch, wave = phase.split('_', 1)
        branch_name = '反射' if branch == 'R' else '透射'
        ax.plot(x, [r[f'{phase}_weak_abs'] for r in comparison], label=f'含 {middle_thickness:g} m 软弱层')
        ax.plot(x, [r[f'{phase}_{reference_key}_abs'] for r in comparison], ls='--', label='无软弱层')
        ax.set(title=f'{branch_name} {wave}', xlabel=xlabel, ylabel='位移幅值系数模 |C|')
        ax.grid(alpha=.2)
        ax.legend()
    for ax in axes.flat[len(phases):]:
        ax.set_visible(False)
    fig.suptitle(f'有/无软弱层的反透射对比｜{title}｜相同参考面')
    fig.savefig(path.with_suffix('.png'), dpi=170)
    plt.close(fig)
    write_table(path, comparison)


def run(model=MODEL,incident=INCIDENT,settings=THEORY,output=None):
    validate(model,incident)
    output=Path(output) if output else Path(__file__).resolve().parent/'rt_coefficients_output'/(incident['side']+'_'+incident['wave'])
    frequencies=sorted(set(settings['frequencies_hz']+[incident['frequency_hz']]))
    if not frequencies or not np.isfinite(frequencies).all() or min(frequencies)<=0:raise ValueError('All frequencies must be >0')
    if not 0<settings['angle_step_deg'] or not 0<=settings['max_angle_deg']<90:raise ValueError('Invalid angle scan')
    angles=sorted(set(np.arange(0,settings['max_angle_deg']+1e-9,settings['angle_step_deg']).tolist()+[incident['angle_deg']]))
    if not 0<settings['spectrum_min_hz']<=settings['spectrum_max_hz'] or settings['spectrum_step_hz']<=0:raise ValueError('Invalid frequency scan')
    output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    # Extend the lower halfspace through the middle layer and retain both reference planes.
    baseline=copy.deepcopy(model);baseline['middle']=copy.deepcopy(model['lower'])
    selected=solve(model,incident)
    phases=sorted(selected['coefficients'])
    case=f"{incident['side']}_{incident['wave']}"
    colors=['#246b92','#c57127','#4c8b55','#95538e']
    maxima=[]
    for frequency in frequencies:
        rows=[flatten(solve(model,incident,a,frequency)) for a in angles]
        write_table(output/f'coefficients_f{frequency:g}Hz.tsv',rows)
        fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
        for name,color in zip(selected['coefficients'],colors):axes[0].plot(angles,[r[name+'_abs'] for r in rows],label=name,color=color)
        for key,label in [('R_energy','反射能量'),('T_energy','透射能量'),('absorption','层内耗散')]:axes[1].plot(angles,[r[key] for r in rows],label=label)
        axes[0].set(ylabel='位移幅值系数模 |C|');axes[1].set(ylabel='法向能流比例')
        for ax in axes:ax.set(xlabel='入射角 / °');ax.grid(alpha=.2);ax.legend()
        thickness='/'.join(f'{model[k+"_thickness_m"]:g}' for k in ('upper','middle','lower'))
        fig.suptitle(f"{thickness} m 三层｜{incident['side']}_{incident['wave']}｜{frequency:g} Hz｜中层 {kind(model['middle'])}")
        fig.savefig(output/f'coefficients_f{frequency:g}Hz.png',dpi=170);plt.close(fig)
        base_rows=[flatten(solve(baseline,incident,a,frequency)) for a in angles]
        write_weak_layer_comparison(output/f'weak_layer_effect_angle_f{frequency:g}Hz.tsv', rows,
                                    base_rows, phases, 'no_weak', 'angle_deg', '入射角 / °',
                                    f'{case}｜{frequency:g} Hz', model['middle_thickness_m'])
        maxima.append(max(r['residual'] for r in rows))
    freq=np.arange(settings['spectrum_min_hz'],settings['spectrum_max_hz']+1e-9,settings['spectrum_step_hz'])
    rows=[flatten(solve(model,incident,frequency_hz=f)) for f in freq]
    write_table(output/'coefficients_frequency_spectrum.tsv',rows)
    fig,axes=plt.subplots(1,2,figsize=(12,4.8),layout='constrained')
    for name,color in zip(selected['coefficients'],colors):axes[0].plot(freq,[r[name+'_abs'] for r in rows],label=name,color=color)
    for key in ('R_energy','T_energy','absorption'):axes[1].plot(freq,[r[key] for r in rows],label=key)
    for ax in axes:ax.set(xlabel='频率 / Hz');ax.grid(alpha=.2);ax.legend()
    axes[0].set(ylabel='|C|');axes[1].set(ylabel='法向能流比例')
    fig.suptitle(f"三层频率响应｜{incident['side']}_{incident['wave']}｜{incident['angle_deg']:g}°")
    fig.savefig(output/'coefficients_frequency_spectrum.png',dpi=170);plt.close(fig)
    base_rows=[flatten(solve(baseline,incident,frequency_hz=f)) for f in freq]
    write_weak_layer_comparison(output/'weak_layer_effect_spectrum.tsv', rows, base_rows, phases,
                                'baseline', 'frequency_hz', '频率 / Hz',
                                f"{case}｜{incident['angle_deg']:g}°", model['middle_thickness_m'])
    snapshot=dict(model=model,incident=incident,scan=settings,
                  selected=flatten(selected),max_boundary_residual=max(maxima+[r['residual'] for r in rows]),
                  convention='exp[i omega(p x+q depth-t)], displacement coefficients at d=h_upper and d=h_upper+h_middle',
                  meaning='Effective whole-layer coherent coefficients, including all internal reverberations and mode conversions; not individual first arrivals.',
                  model_type='low-frequency Biot with constant Darcy drag' if kind(model['middle'])=='poroelastic' else 'lossless isotropic elastic/acoustic',
                  outer_geometry='Top and bottom are halfspaces for theory; physical 500/499m truncate numerical grid only; no exterior free-surface echoes included.')
    (output/'model_parameters.json').write_text(json.dumps(snapshot,ensure_ascii=False,indent=2),encoding='utf-8')
    print('Theory:',output)
    return snapshot


if __name__=='__main__':
    import sys
    sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
    from wave_packet import theory_cli
    theory_cli()
