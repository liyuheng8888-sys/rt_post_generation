"""Finite force-couple source implementation retained from the validated input generator."""
import math
import csv
import numpy as np

def kind(m):
    return 'poroelastic' if m.get('poroelastic_enabled',False) else m['type']

def names(m):
    return ['P'] if kind(m)=='fluid' else ['P','SV']

def speed(m,wave):
    return m['vp'] if wave=='P' else m['vs']

def boundaries(model):
    a=model['upper_thickness_m'];b=a+model['middle_thickness_m']
    return a,b,b+model['lower_thickness_m']

def validate(model,incident):
    for key in ('size_x_m','size_y_m','upper_thickness_m','middle_thickness_m','lower_thickness_m'):
        if not np.isfinite(model[key]) or model[key]<=0:raise ValueError(key+'必须为正数')
    ids=list(model['material_ids'].values())
    if len(ids)!=3 or len(set(ids))!=3 or any(type(x) is not int or x<=0 for x in ids):raise ValueError('三层材料号必须为不同正整数')
    for layer in ('upper','middle','lower'):
        m=model[layer];k=kind(m)
        if k=='poroelastic':
            if layer!='middle':raise ValueError('多孔开关只用于中层')
            b=m['biot']
            for key,value in b.items():
                if not np.isfinite(value) or value<0:raise ValueError('多孔参数无效：'+key)
            if not 0<b['porosity']<1 or b['tortuosity']<1:raise ValueError('孔隙率应在0~1之间，曲折度>=1')
            for key in ('rho_s','rho_f','permeability_m2','grain_bulk_Pa','fluid_bulk_Pa','frame_bulk_Pa','frame_shear_Pa'):
                if b[key]<=0:raise ValueError(key+'必须为正数')
            alpha=1-b['frame_bulk_Pa']/b['grain_bulk_Pa']
            if not 0<alpha<=1 or (alpha-b['porosity'])/b['grain_bulk_Pa']+b['porosity']/b['fluid_bulk_Pa']<=0:
                raise ValueError('多孔骨架/颗粒/流体模量不相容')
            rho=(1-b['porosity'])*b['rho_s']+b['porosity']*b['rho_f']
            if rho*b['tortuosity']*b['rho_f']/b['porosity']-b['rho_f']**2<=0:raise ValueError('多孔质量矩阵无效')
        elif k in ('solid','fluid'):
            if not all(np.isfinite(m[x]) for x in ('rho','vp','vs')) or m['rho']<=0 or m['vp']<=0:raise ValueError(layer+'密度/波速无效')
            if k=='fluid' and m['vs']!=0:raise ValueError('流体Vs必须为0')
            if k=='solid' and not (m['vs']>0 and m['vp']**2>4*m['vs']**2/3):raise ValueError('固体必须有正剪切/体积模量')
        else:raise ValueError('材料类型应为solid/fluid')
    if incident['side'] not in ('upper','lower') or incident['wave'] not in names(model[incident['side']]):raise ValueError('入射侧或波型无效')
    if not np.isfinite(incident['angle_deg']) or not 0<=incident['angle_deg']<90:raise ValueError('角度应为[0,90)')
    if not np.isfinite(incident['frequency_hz']) or incident['frequency_hz']<=0:raise ValueError('主频必须为正数')

def samples(lo,hi,step):
    if not np.isfinite([lo,hi,step]).all() or step<=0 or hi<lo:raise ValueError('Invalid sample range')
    values=np.arange(lo,hi+1e-8,step)
    if not len(values) or hi-values[-1]>1e-7:values=np.append(values,hi)
    return values

def materials(model,out):
    rows=[];poro=[]
    for layer in ('upper','middle','lower'):
        m=model[layer];k=kind(m);mid=model['material_ids'][layer]
        if k=='poroelastic':
            rows.append(f'3 {mid} 0 0 0 9999 9999 0')
            b=m['biot'];kap=b['permeability_m2']
            values=[b['rho_s'],b['rho_f'],b['porosity'],b['tortuosity'],kap,0,0,kap,0,kap,
                    b['grain_bulk_Pa'],b['fluid_bulk_Pa'],b['frame_bulk_Pa'],b['viscosity_Pa_s'],b['frame_shear_Pa']]
            poro.append(' '.join(f'{v:.16g}' for v in values))
        else:rows.append(f"{1 if k=='fluid' else 2} {mid} {m['rho']:.16g} {m['vp']:.16g} {m['vs']:.16g} 9999 9999 0")
    (out/'nummaterial_velocity_file').write_text('\n'.join(rows)+'\n',encoding='ascii')
    # No material ID in poroelastic rows; lines follow porous entries in velocity file.
    (out/'nummaterial_poroelastic_file').write_text('\n'.join(poro)+('\n' if poro else ''),encoding='ascii')

def write_csv(path,fields,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)

def finite_sources(model,incident,inner,research,normal,ref):
    """目标为外模型左右内缩100m；按顶面/底面/入射介质自动裁剪倾斜中心线。"""
    x0,x1,d0,d1=research;top,bottom,_=boundaries(model)
    margin=max(10.,incident['min_point_source_boundary_distance_m'])
    sep=incident['couple_separation_m'];pad=sep/2 if kind(model[incident['side']])=='solid' else 0.
    for key in ('source_spacing_m','source_y_spacing_m','couple_separation_m','force_factor'):
        if not np.isfinite(incident[key]) or incident[key]<=0:raise ValueError(key+'必须为正数')
    if not np.isfinite([incident['peak_time_s'],incident['uniform_time_shift_s']]).all():raise ValueError('源时刻必须有限')
    total=boundaries(model)[2]
    xmargin=incident['source_outer_x_margin_m'];top_clearance=incident['source_top_clearance_m']
    if not np.isfinite([xmargin,top_clearance]).all() or min(xmargin,top_clearance)<0:
        raise ValueError('外模型X内缩距离、顶面净距必须为有限非负数')
    # 法向和切向两对力的最大坐标偏移，确保真实力点而不仅是中心满足竖向净距。
    extent=pad*float(max(abs(normal)))
    # 外边界与层间界面分开约束。四个力点的X/埋深最大偏移均为extent。
    # 旧配置缺少新增项时也执行100m下限，避免沿用旧50m顶面净距。
    outer_clearance=incident.get('source_outer_boundary_clearance_m',100.0)
    if not np.isfinite(outer_clearance) or outer_clearance<100.0:
        raise ValueError('source_outer_boundary_clearance_m必须至少100m')
    outer_clearance=max(outer_clearance,margin)
    xclear=max(xmargin,outer_clearance);top_clearance=max(top_clearance,outer_clearance)
    xmin=xclear+extent;xmax=model['size_x_m']-xclear-extent
    dmin=top_clearance+extent;dmax=total-outer_clearance-extent
    if incident['side']=='upper':dmax=min(dmax,top-margin-extent)
    else:dmin=max(dmin,bottom+margin+extent)
    if xmin>=xmax or dmin>=dmax:raise ValueError('净距参数过大，入射介质内没有可布置震源的空间')
    ends=clipped_front(ref,normal,0.,(xmin,xmax,dmin,dmax))
    if len(ends)!=2 or np.linalg.norm(np.diff(ends,axis=0))<1e-6:
        raise ValueError('震源线与允许区域无有效交段，请检查参考点或净距参数')
    ends=sorted(ends,key=lambda point:point[0])
    if ends[0][0]>xmin+1e-6 or ends[1][0]<xmax-1e-6:
        print(f'源线已按顶面/底面/介质界面净距自动截短：X={ends[0][0]:.3f}~{ends[1][0]:.3f}m；倾角、参考点及Y范围保持不变。')
    a,b=np.array(ends);length=np.linalg.norm(b-a);tangent=np.array([-normal[1],normal[0]])
    centers=[];forces=[]
    # 保留原Y范围和间隔；若用户换模型后Y范围违反净距，明确报错而不暗中平移Y。
    source_ys=samples(inner['y'][0]+margin,inner['y'][1]-margin,incident['source_y_spacing_m'])
    if min(source_ys[0],model['size_y_m']-source_ys[-1])<outer_clearance-1e-7:
        raise ValueError('原Y范围不能满足外边界100m净距；请调整研究区域Y边界或源面Y范围，程序不会自动移动Y')
    for distance in samples(0,length,incident['source_spacing_m']):
        center=a+(b-a)*distance/length
        for y in source_ys:
            cid=len(centers)+1;centers.append(dict(center_id=cid,x_m=float(center[0]),y_m=float(y),burial_m=float(center[1])))
            pairs=[(normal,normal),(tangent,tangent)] if incident['wave']=='P' else [(normal,tangent),(tangent,normal)]
            sources=[(center,np.array([1.,0.]))] if kind(model[incident['side']])=='fluid' else [
                (center+sign*sep/2*separation,sign*force) for separation,force in pairs for sign in (-1,1)]
            for pos,force in sources:
                if min(pos[0],model['size_x_m']-pos[0],y,model['size_y_m']-y,pos[1],total-pos[1])<outer_clearance-1e-7:
                    raise ValueError('实际力点违反大模型六个外边界至少100m净距')
                if pos[1]<top_clearance-1e-7:raise ValueError('实际力点距顶面过近')
                if not (pos[1]<=top-margin+1e-7 if incident['side']=='upper' else pos[1]>=bottom+margin-1e-7):raise ValueError('实际力点距介质界面过近')
                forces.append(dict(force_id=len(forces)+1,center_id=cid,x_m=float(pos[0]),y_m=float(y),burial_m=float(pos[1]),
                                   E=float(force[0]),N=0.,Z_UP=float(-force[1]),time_shift_s=incident['uniform_time_shift_s']))
    return ends,centers,forces

def source_boundary_audit(forces,model,incident):
    """按最终写文件精度核验每个力点；全部距离均为米，深度向下。"""
    coords=np.array([[round(r['x_m'],9),round(r['y_m'],9),round(r['burial_m']/1000,12)*1000] for r in forces])
    total=boundaries(model)[2]
    gaps=np.column_stack((coords[:,0],model['size_x_m']-coords[:,0],coords[:,1],model['size_y_m']-coords[:,1],coords[:,2],total-coords[:,2]))
    labels=['xmin','xmax','ymin','ymax','top','bottom']
    required=max(100.,incident.get('source_outer_boundary_clearance_m',100.),incident['min_point_source_boundary_distance_m'])
    limits=np.array([max(required,incident['source_outer_x_margin_m'])]*2+[required]*2+[max(required,incident['source_top_clearance_m']),required])
    if not np.isfinite(gaps).all() or np.any(gaps<limits-1e-7):
        raise ValueError('写出前检查失败：实际力源点未满足外边界净距；不输出该算例')
    return dict(required_clearance_m=dict(zip(labels,limits.tolist())),
                minimum_clearance_m=dict(zip(labels,gaps.min(axis=0).tolist())),
                closest_force_id=dict(zip(labels,[forces[k]['force_id'] for k in gaps.argmin(axis=0)])),
                checked_force_count=len(forces),passed=True,
                scope='仅证明源点外边界净距合格，不代表网格Jacobian或MPI定位已通过')

def clipped_front(ref,normal,offset,bounds):
    """n·(point-ref)=offset，与X-埋深矩形的交线。"""
    x0,x1,d0,d1=bounds;points=[]
    for x in (x0,x1):
        if abs(normal[1])>1e-12:
            d=ref[1]+(offset-normal[0]*(x-ref[0]))/normal[1]
            if d0-1e-8<=d<=d1+1e-8:points.append([x,d])
    for d in (d0,d1):
        if abs(normal[0])>1e-12:
            x=ref[0]+(offset-normal[1]*(d-ref[1]))/normal[0]
            if x0-1e-8<=x<=x1+1e-8:points.append([x,d])
    unique=np.unique(np.round(points,9),axis=0) if points else []
    return np.asarray(unique[:2]).tolist()
