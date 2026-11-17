#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立三层震源/台站生成程序。只修改本文件顶部参数，然后直接运行。"""
from pathlib import Path
import json
import csv
import math
import argparse
import copy
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==================== 用户编辑区：只影响震源台站输入 ====================
# 三层材料、厚度、入射角和台站位置都在本文件内设置。
# 不读取或调用理论解程序，不在本程序计算反透射系数。
# 以下边界为网格绝对XYZ，Z向上。局部位置以小长方体左/前/顶为原点，埋深向下。
OUTER_BOUNDS_M = {"x": [0.0, 4000.0], "y": [0.0, 600.0], "z": [0.0, 4000.0]}
# 匹配 Mesh_Nested4000x600x4000_strictOuter20_inner10_weak1m_20261112 实际网格。
# 小盒随新网格沿X/Z各平移1000m；材料按原始网格：上1/中2/下3。不得修改mesh。
RESEARCH_BOUNDS_M = {"x": [1500.0, 2500.0], "y": [100.0, 500.0], "z": [1500.0, 2500.0]}
# 在参考点波峰时刻之后多少秒再画一条波峰线；None自动取到达参考界面交点的时刻。
WAVEFRONT_AFTER_REFERENCE_S = None
MODEL = {
    "size_x_m": 1000.0,
    "size_y_m": 400.0,
    "upper_thickness_m": 499.0,  # 实际网格绝对Z=2001为上界面，Z=2000为下界面
    "middle_thickness_m": 1.0,
    "lower_thickness_m": 500.0,
    "upper": {"name": "上层固体", "type": "solid", "rho": 2000.0, "vp": 3000.0, "vs": 1800.0},
    "middle": {
        "name": "中间软弱层", "type": "solid",  # solid普通弹性 / fluid流体
        "rho": 1000.0, "vp": 1450.0, "vs": 850.0,  # 当前nummaterial_velocity_file的2号材料数值
        "poroelastic_enabled": False,  # True: 使用下方Biot参数，忽略本层rho/vp/vs
        "biot": {
            "rho_s": 2500.0, "rho_f": 1000.0,  # 颗粒/孔隙流体密度 kg/m³
            "porosity": 0.35, "tortuosity": 2.0,
            "permeability_m2": 1.0e-10,  # 各向同性渗透率 m²
            "viscosity_Pa_s": 1.0e-3,
            "grain_bulk_Pa": 36.0e9, "fluid_bulk_Pa": 2.25e9,
            "frame_bulk_Pa": 1.0e8, "frame_shear_Pa": 2.0e7,
        },
    },
    "lower": {"name": "下层固体", "type": "solid", "rho": 2500.0, "vp": 4500.0, "vs": 2100.0},  # 文件中3号材料
    # 流体/多孔界面：open连续孔压与总体积流量；sealed孔隙相对流量为0。
    # 普通固体/多孔界面使用粘结、sealed。不支持部分饱和或JKD高频修正。
    "upper_middle_pore_condition": "open",
    "middle_lower_pore_condition": "sealed",
    # 材料号可改；不自动覆盖已有网格。请与materials_file逐项对应。
    "material_ids": {"upper": 1, "middle": 2, "lower": 3},
}

INCIDENT = {
    "side": "lower",                 # upper / lower，对称入口
    "wave": "P",                     # 固体P/SV，流体只能P
    "angle_deg": 30.0,                # 与界面法线夹角；水平沿+X，Y均匀
    "frequency_hz": 95.0,            # 力源Ricker主频；接收位移/加速度形状由传播决定
    "source_spacing_m": 10.0,         # 震源中心沿倾斜面的间隔，可设5m
    "source_y_spacing_m": 10.0,       # 小长方体Y方向震源间隔，可设5m
    "source_outer_boundary_clearance_m": 100.0, # 每一个实际力点距大模型六个外边界至少100m；硬下限100m
    "source_outer_x_margin_m": 100.0, # 实际力点距左右边界的净距，可进一步增大；不再仅约束源中心
    "source_top_clearance_m": 100.0,  # 实际力点距顶面净距，可进一步增大；倾角越界时自动截短X长度
    # 保持参考点和倾角不变，只延长X；Y范围/间距不变。
    # 若另一端接近介质界面或大模型底面，同样自动截断，保留下面的边界净距。
    "couple_separation_m": 1.0,       # 固体：每中心两对正负力，力点间距
    "force_factor": 1.0e15,          # 每条FORCE的幅值，不再指定解析位移幅值
    "uniform_time_shift_s": 0.0,     # 所有源写相同time shift；部分求解器会减去共同最小值
    "reference_x_m": 300.0,  # 小长方体局部X；输出全局X=1500+此值
    "reference_depth_m": {"upper": 50.0, "lower": 950.0},  # 从小盒顶面向下量；输出自动加外顶面至小盒顶面距离（当前1500m）
    "peak_time_s": 0.0,               # 记录时间轴上的源波峰时刻（用于示意/后处理），不是自动附加的时移
    "min_point_source_boundary_distance_m": 10.0,  # 强制下限10m
}

RECEIVERS = {
    "plane_y_m": 200.0,  # 小长方体局部Y，范围0~400；中心200对应全局Y=300
    "interface_hit_x_m": 500.0,  # 小长方体局部X
    "auto_reflection_line": True,    # True自动选线；False恢复固定interface_hit_x_m交点
    "auto_common_intersection": True, # 自动搜索一个公共交点；所有测线延长线必须经过它，入射线只有一条
    "auto_ray_distance_max_m": 1000.0, # 自动模式沿测线搜索的最大距离；实际仍裁在研究小盒内
    "min_clean_line_length_m": 100.0, # 各测线优先达到的连续隔离段长度
    "allow_short_clean_lines": True, # 用户允许：共同交点无法同时满足100m时，保留更短的合格段
    "min_short_clean_line_length_m": 20.0, # 短段仍至少20m；低于此值停止，不加入混波站凑长度
    "incident_min_source_distance_m": 100.0, # 入射站必须位于源面前方，且离源面至少此距离
    "source_foot_edge_clearance_m": 50.0, # 界面入射射线在有限源线上的垂足离两端至少此距离
    # 自动模式筛选每个台站：入射站也检查与反射P/SV的分离；不能只检查反射站。
    # 不满足条件的台站不写入理论测线。无可行段时明确报告，不退回混波测线。
    "reflection_isolation": "time", # time按到时分离；space严格只布置入射波束外台站
    "isolation_half_periods": 2.0,   # 入射/反射包络假定半宽=此值/主频，不是严格脉冲支撑
    "isolation_guard_s": 0.005,      # 两个理论波包之间额外留白秒数
    "incident_beam_clearance_m": 10.0,  # 反射台站到入射波束侧边的法向净距
    "candidate_step_m": 1.0,
    "min_reflected_stations": 3,      # 每种传播反射波至少有这些满足所选隔离条件的台站
    "ray_step_m": 5.0,
    "ray_distance_min_m": 20.0,
    "ray_distance_max_m": 300.0,
    "outer_boundary_margin_m": 10.0,  # 台站到研究小长方体边界的内缩距离（兼容原参数名）
    "generate_grid": False,
    "grid_spacing_m": 10.0,
    # 按研究要求，不在软弱层内部或两个界面上布站。
}

# ==================== 用户编辑区结束 ====================

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


def global_geometry(model,incident,receivers,outer,inner):
    """把研究区域局部参数转换为求解器全局XY、相对外模型顶面的埋深。"""
    for axis in ('x','y','z'):
        a,b=outer[axis];c,d=inner[axis]
        if not np.isfinite([a,b,c,d]).all() or not a<=c<d<=b:raise ValueError('研究区域必须包含在外模型内：'+axis)
        if a!=0:raise ValueError('当前全局坐标约定要求外模型XYZ最小值为0')
    if not np.allclose([model['size_x_m'],model['size_y_m'],boundaries(model)[2]],
                       [inner[k][1]-inner[k][0] for k in ('x','y','z')]):
        raise ValueError('MODEL尺寸/三层厚度之和必须等于研究小长方体尺寸')
    offset=outer['z'][1]-inner['z'][1]
    m,i,r=copy.deepcopy((model,incident,receivers))
    m['size_x_m']=outer['x'][1];m['size_y_m']=outer['y'][1]
    m['upper_thickness_m']+=offset
    m['lower_thickness_m']+=inner['z'][0]-outer['z'][0]
    i['reference_x_m']+=inner['x'][0]
    i['reference_depth_m']={k:v+offset for k,v in i['reference_depth_m'].items()}
    r['plane_y_m']+=inner['y'][0];r['interface_hit_x_m']+=inner['x'][0]
    bounds=(inner['x'][0],inner['x'][1],offset,offset+boundaries(model)[2])
    return m,i,r,bounds


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


def select_common_lines(model,incident,receivers,research,source_line,normal,ref,entries):
    """联合搜索一个公共界面交点；所有震相必须同时满足连续段/隔离要求。"""
    x0,x1,d0,d1=research;top,bottom,_=boundaries(model)
    side=incident['side'];near=top if side=='upper' else bottom
    velocity=speed(model[side],incident['wave']);slow=normal/velocity
    tangent=np.array([normal[1],-normal[0]])
    beam_lo,beam_hi=sorted(float(tangent@np.asarray(pt)) for pt in source_line)
    rm=receivers['outer_boundary_margin_m'];step=receivers['ray_step_m']
    half=receivers['isolation_half_periods']/incident['frequency_hz']
    required=2*half+receivers['isolation_guard_s']
    minimum=receivers['min_clean_line_length_m'];source_gap=receivers['incident_min_source_distance_m']
    allow_short=receivers.get('allow_short_clean_lines',False)
    short_minimum=receivers.get('min_short_clean_line_length_m',minimum)
    if not np.isfinite(short_minimum) or not 0<short_minimum<=minimum:raise ValueError('短段下限必须在(0,目标段长]内')
    edge_gap=receivers['source_foot_edge_clearance_m'];clearance=receivers['incident_beam_clearance_m']
    max_distance=receivers['auto_ray_distance_max_m'];candidate_step=receivers['candidate_step_m']
    if not np.isfinite([minimum,source_gap,edge_gap,max_distance,candidate_step,step]).all() or min(minimum,max_distance,candidate_step,step)<=0 or min(source_gap,edge_gap)<0:
        raise ValueError('自动选线长度、间隔必须为正数；源面/边缘净距必须非负')
    # 台站距离始终落在ray_step_m整数步长上，不在末尾补一个不足步长的站。
    ds=np.arange(math.ceil(receivers['ray_distance_min_m']/step)*step,max_distance+1e-8,step)
    if not len(ds):raise ValueError('自动搜索最大距离小于最小距离')
    reflect=[e for e in entries if e[0].startswith('R_')]

    def incident_present(points,margin=0.):
        u=points@tangent
        return ((points-ref)@normal>=0)&(u>=beam_lo-margin)&(u<=beam_hi+margin)

    def hit_valid(hx,edge=0.):
        points=np.column_stack((np.atleast_1d(hx),np.full(np.size(hx),near)))
        u=points@tangent
        return ((points-ref)@normal>0)&(u>=beam_lo+edge)&(u<=beam_hi-edge)

    def phase_time(points,entry):
        ph,layer,direction,depth,_=entry
        if ph.startswith('I_'):return incident['peak_time_s']+(points-ref)@slow
        return (incident['peak_time_s']+slow[1]*(near-ref[1])+slow[0]*(points[:,0]-ref[0])+
                direction[1]/speed(model[layer],ph.split('_',1)[1])*(points[:,1]-depth))

    def reflection_present(points,entry):
        _,_,direction,depth,_=entry
        distance=(points[:,1]-depth)/direction[1]
        hx=points[:,0]-distance*direction[0]
        # 不能刚越过几何波束边缘就认定竞争波不存在；按波束净距向外扩展保护带。
        return (distance>=0)&hit_valid(hx,-clearance)

    selected={};reports={};options={}
    for entry in entries:
        ph,layer,direction,depth,multiplier=entry
        options[ph]={}
        for hx in samples(x0+rm,x1-rm,candidate_step):
            if not hit_valid(hx,edge_gap)[0]:continue
            # 透射测线在远侧界面的起点作几何平移，确保延长后仍穿过近侧同一个公共点。
            # 它是相位采样线，不能解释为软弱层内部全部转换波的唯一真实路径。
            anchor_x=hx+(depth-near)*direction[0]/direction[1]
            points=np.array([anchor_x,depth])+multiplier*ds[:,None]*direction
            keep=((points[:,0]>=x0+rm)&(points[:,0]<=x1-rm)&(points[:,1]>=d0+rm)&(points[:,1]<=d1-rm))
            keep&=(points[:,1]<top-1e-8) if layer=='upper' else (points[:,1]>bottom+1e-8)
            target_time=phase_time(points,entry)
            separation=np.full(len(ds),np.inf)
            in_beam=incident_present(points,clearance)
            gap=(points-ref)@normal
            if ph.startswith('I_'):
                keep&=incident_present(points)&(gap>=source_gap-1e-8)
                competitors=reflect
            elif ph.startswith('R_'):
                ti=incident['peak_time_s']+(points-ref)@slow
                separation=np.where(in_beam,target_time-ti,np.inf)
                if receivers['reflection_isolation']=='space':keep&=~in_beam
                competitors=[e for e in reflect if e[0]!=ph]
            else:
                # 透射测线保留完整传播段；薄层内P/SV转换时延仍由后处理理论频谱处理。
                competitors=[]
            for other in competitors:
                delta=np.abs(target_time-phase_time(points,other))
                separation=np.minimum(separation,np.where(reflection_present(points,other),delta,np.inf))
            keep&=separation>=required-1e-10
            indices=np.flatnonzero(keep)
            if not len(indices):continue
            blocks=np.split(indices,np.where(np.diff(indices)>1)[0]+1)
            for block in blocks:
                span=float(ds[block[-1]]-ds[block[0]])
                if len(block)<receivers['min_reflected_stations']:continue
                min_sep=float(np.min(separation[block]));u=float(tangent@np.array([hx,near]))
                foot_edge=min(u-beam_lo,beam_hi-u)
                rank=(span,min(min_sep,1.),foot_edge,-abs(hx-receivers['interface_hit_x_m']))
                old=options[ph].get(float(hx))
                if old is None or rank>old[0]:
                    options[ph][float(hx)]=(rank,float(hx),block,points,separation,in_beam,target_time,gap)
    shared=set.intersection(*(set(by_x) for by_x in options.values()))
    def joint_rank(hx):
        ranks=[by_x[hx][0] for by_x in options.values()]
        return (min(r[0] for r in ranks),sum(r[0] for r in ranks),min(r[1] for r in ranks),
                min(r[2] for r in ranks),-abs(hx-receivers['interface_hit_x_m']))
    feasible=[hx for hx in shared if joint_rank(hx)[0]>=minimum-1e-8]
    common_x=max(feasible or shared,key=joint_rank) if shared else None
    best_span=joint_rank(common_x)[0] if common_x is not None else 0.
    shortened=not feasible and allow_short and best_span>=short_minimum-1e-8
    success=bool(feasible) or shortened
    if not success:
        print(f'公共交点约束下无可行解：要求各线连续{minimum:g}m，当前最好的最短测线只有{best_span:g}m。')
    elif shortened:
        print(f'采用用户允许的较短合格段：公共交点X={common_x:g}m，最短连续测线{best_span:g}m；源面保持不动。')
    for entry in entries:
        ph,layer,direction,depth,multiplier=entry
        if common_x is None:
            reports[ph]=dict(status='unavailable',reason='各震相没有同时可用的公共交点',feasible_segments=0)
            continue
        winner=options[ph][common_x]
        _,hx,block,points,separation,in_beam,target_time,gap=winner
        records=[]
        for k in block:
            records.append(dict(distance_m=float(ds[k]),point=points[k].tolist(),
                                phase_arrival_s=float(target_time[k]),
                                min_competing_arrival_separation_s=float(separation[k]) if np.isfinite(separation[k]) else None,
                                incident_beam_present=bool(in_beam[k]),signed_source_plane_distance_m=float(gap[k])))
        anchor_x=hx+(depth-near)*direction[0]/direction[1]
        if success:selected[ph]=dict(hit_x_m=anchor_x,common_hit_x_m=hx,records=records)
        reports[ph]=dict(status='selected' if success else 'unavailable',common_hit_x_m=hx,
                         requested_clean_length_m=minimum,short_segments_accepted=bool(shortened),
                         common_hit_depth_m=near,interface_hit_x_m=anchor_x,distance_start_m=records[0]['distance_m'],
                         distance_end_m=records[-1]['distance_m'],continuous_length_m=records[-1]['distance_m']-records[0]['distance_m'],
                         station_count=len(records),full_length_candidate_count=len(feasible),
                         feasible_segments=sum(joint_rank(h)[0]>=(short_minimum if shortened else minimum)-1e-8 for h in shared),
                         isolation_applied=not ph.startswith('T_'),required_separation_s=required if not ph.startswith('T_') else None,
                         minimum_separation_s=min((r['min_competing_arrival_separation_s'] for r in records if r['min_competing_arrival_separation_s'] is not None),default=None),
                         source_foot_edge_distance_m=winner[0][2])
        print(f'{ph}: 公共交点X={hx:g}m，沿线距所在界面{records[0]["distance_m"]:g}~{records[-1]["distance_m"]:g}m，共{len(records)}站'+('。' if success else '（仅诊断，不输出台站）。'))
    return selected,reports


def run(model=MODEL,incident=INCIDENT,receivers=RECEIVERS,out=None,outer=OUTER_BOUNDS_M,inner=RESEARCH_BOUNDS_M):
    local_inputs=copy.deepcopy(dict(model=model,incident=incident,receivers=receivers))
    validate(model,incident)
    model,incident,receivers,research=global_geometry(model,incident,receivers,outer,inner)
    x0,x1,d0,d1=research
    if kind(model['middle'])=='poroelastic':
        for layer,key in [('upper','upper_middle_pore_condition'),('lower','middle_lower_pore_condition')]:
            if kind(model[layer])=='fluid' and model[key]!='open':
                raise ValueError('SPECFEM3D 4.1.1 fluid/Biot coupling requires open pores. '
                                 'Sealed-fluid theory is supported, but numeric input needs a separate coupling patch.')
    out=Path(out) if out else Path(__file__).resolve().parent/'output'
    dtop,dbottom,total=boundaries(model);lx=model['size_x_m'];ly=model['size_y_m']
    side=incident['side'];inc=model[side];sgn=1 if side=='upper' else -1
    reference_depth=incident['reference_depth_m'][side];reference_x=incident['reference_x_m']
    margin=max(10.,incident['min_point_source_boundary_distance_m'])
    source_y=sum(inner['y'])/2
    point=[reference_x,source_y,reference_depth]
    if not (x0<reference_x<x1 and d0<reference_depth<d1):raise ValueError('源参考点必须在研究小长方体内部')
    distances=[point[0],lx-point[0],point[1],ly-point[1],point[2],total-point[2]]
    if not np.isfinite(distances).all() or min(distances)<margin-1e-9:raise ValueError(f'All point sources must be >= {margin:g}m from outer boundaries')
    if not (reference_depth<dtop if side=='upper' else reference_depth>dbottom):raise ValueError('Source-plane reference must be in the incident medium')
    if min(abs(reference_depth-dtop),abs(reference_depth-dbottom))<10:raise ValueError('Keep source reference >=10m from interfaces too')
    v=speed(inc,incident['wave']);angle=math.radians(incident['angle_deg'])
    normal=np.array([math.sin(angle),sgn*math.cos(angle)])
    slow=normal/v;ref=np.array([reference_x,reference_depth])
    inc_interface=dtop if side=='upper' else dbottom
    source_line,source_centers,forces=finite_sources(model,incident,inner,research,normal,ref)
    boundary_audit=source_boundary_audit(forces,model,incident)
    y=receivers['plane_y_m'];rm=receivers['outer_boundary_margin_m']
    if not 0<rm<min(x1-x0,inner['y'][1]-inner['y'][0],d1-d0)/2 or not inner['y'][0]+rm<=y<=inner['y'][1]-rm:raise ValueError('Receiver margins or local Y plane invalid')
    if not x0<=receivers['interface_hit_x_m']<=x1:raise ValueError('界面参考交点必须在研究区域内')
    if not 0<receivers['ray_distance_min_m']<=receivers['ray_distance_max_m']:raise ValueError('Invalid receiver distances')
    station_list=[];lookup={};groups=[]
    def add(x,d,group,number):
        if not (x0+rm-1e-8<=x<=x1-rm+1e-8 and d0+rm-1e-8<=d<=d1-rm+1e-8):return
        key=(round(float(x),7),round(float(d),7))
        if key not in lookup:
            name=f'S{len(station_list)+1:06d}';lookup[key]=name
            layer='upper' if d<dtop else 'lower' if d>dbottom else 'middle'
            station_list.append(dict(station=name,network='DB',x_m=float(x),y_m=y,burial_m=float(d),layer=layer))
        groups.append(dict(group=group,station=lookup[key],requested_number=number,x_m=float(x),y_m=y,burial_m=float(d)))
    if receivers['generate_grid']:
        grid=receivers['grid_spacing_m']
        depths=[d for d in samples(d0+rm,d1-rm,grid) if not dtop-1e-6<=d<=dbottom+1e-6]
        for d in depths:
            for x in samples(x0+rm,x1-rm,grid):add(x,d,'grid',0)
    # No receiver inside the weak layer or exactly on its interfaces.
    p=slow[0];hit_x=receivers['interface_hit_x_m'];rays=[];distance_rows=[]
    phase_entries=[('I_'+incident['wave'],side,normal,inc_interface,-1)]
    trans='lower' if side=='upper' else 'upper'
    for action,layer,dsign,depth in [('R',side,-sgn,inc_interface),('T',trans,sgn,dbottom if trans=='lower' else dtop)]:
        for wave in names(model[layer]):
            sine=p*speed(model[layer],wave)
            if sine>=1-1e-10:continue
            phase_entries.append((action+'_'+wave,layer,np.array([sine,dsign*math.sqrt(1-sine*sine)]),depth,1))
    tangent=np.array([normal[1],-normal[0]])
    beam_lo,beam_hi=sorted(float(tangent@np.asarray(pt)) for pt in source_line)
    clearance=receivers['incident_beam_clearance_m']
    if not np.isfinite(clearance) or clearance<0:raise ValueError('波束净距必须非负')
    def beam_gap(pos):
        coordinate=float(tangent@pos)
        return max(beam_lo-coordinate,coordinate-beam_hi)
    def illuminated_hit(hx):
        h=np.array([hx,inc_interface])
        return float(normal@(h-ref))>0 and beam_lo-1e-8<=tangent@h<=beam_hi+1e-8
    distances_sample=samples(receivers['ray_distance_min_m'],receivers['ray_distance_max_m'],receivers['ray_step_m'])
    def clean_points(hx,entry):
        ph,layer,direction,depth,multiplier=entry;points=[]
        for distance in distances_sample:
            pos=np.array([hx,depth])+distance*direction
            if not (x0+rm<=pos[0]<=x1-rm and d0+rm<=pos[1]<=d1-rm):continue
            if not (pos[1]<dtop if layer=='upper' else pos[1]>dbottom):continue
            if not isolated(pos,entry):continue
            points.append(pos)
        return points
    reflected=[e for e in phase_entries if e[0].startswith('R_')]
    isolation=receivers['reflection_isolation'];half=receivers['isolation_half_periods']/incident['frequency_hz']
    guard=receivers['isolation_guard_s']
    if isolation not in ('space','time') or not np.isfinite([half,guard]).all() or half<=0 or guard<0:raise ValueError('隔离方式/波包宽度无效')
    def timing(pos,entry):
        ph,layer,direction,depth,multiplier=entry
        arrival_i=incident['peak_time_s']+float(slow@(pos-ref))
        # 同介质反射相位：界面水平慢度相同，所有点直接用相位到时公式。
        arrival_r=incident['peak_time_s']+slow[1]*(inc_interface-ref[1])+p*(pos[0]-ref[0])+direction[1]/speed(model[layer],ph.split('_',1)[1])*(pos[1]-depth)
        return arrival_i,arrival_r
    def isolated(pos,entry):
        if beam_gap(pos)>clearance+1e-8:return True
        if isolation=='space':return False
        ti,tr=timing(pos,entry)
        return tr-ti>=2*half+guard-1e-9
    joint_auto=receivers['auto_reflection_line'] and receivers.get('auto_common_intersection',True)
    automatic_lines={};line_reports={}
    if joint_auto:
        automatic_lines,line_reports=select_common_lines(model,incident,receivers,research,source_line,normal,ref,phase_entries)
        if not automatic_lines:
            out.mkdir(parents=True,exist_ok=True)
            diagnosis=dict(status='no_common_intersection',local_inputs=local_inputs,line_selection=line_reports,
                           message='未生成本次STATIONS。不能用多个交点、重叠站或短于要求的段回退；该目录中若有旧结果，不属于本次生成。')
            (out/'line_selection_failed.json').write_text(json.dumps(diagnosis,ensure_ascii=False,indent=2),encoding='utf-8')
            raise ValueError('当前源面/研究区域/波包宽度下，找不到所有测线共享的合格交点。详见 '+str(out/'line_selection_failed.json'))
        receivers['interface_hit_x_by_phase_m']={ph:item['hit_x_m'] for ph,item in automatic_lines.items()}
        # 透射线穿过远侧界面的位置允许有薄层厚度引起的几何偏移，但所有延长线共点。
        receivers['common_intersection_x_depth_m']=[next(iter(automatic_lines.values()))['common_hit_x_m'],inc_interface]
        receivers['short_clean_lines_accepted']=any(r.get('short_segments_accepted',False) for r in line_reports.values())
        receivers['available_common_segment_length_m']=min(r['continuous_length_m'] for r in line_reports.values())
        receivers['isolated_phase_half_periods']={ph:receivers['isolation_half_periods'] for ph in automatic_lines if ph.startswith(('I_','R_'))}
        if 'I_'+incident['wave'] in automatic_lines:
            hit_x=automatic_lines['I_'+incident['wave']]['hit_x_m'];receivers['interface_hit_x_m']=hit_x
    candidates=[]
    if receivers['auto_reflection_line'] and not joint_auto:
        for hx in samples(x0+rm,x1-rm,receivers['candidate_step_m']):
            if not illuminated_hit(hx):continue
            counts=[len(clean_points(hx,e)) for e in reflected]
            if counts and min(counts)>=receivers['min_reflected_stations']:candidates.append(float(hx))
        if candidates:
            hit_x=candidates[(len(candidates)-1)//2]
            receivers['interface_hit_x_m']=hit_x
    counts=[len(clean_points(hit_x,e)) for e in reflected]
    reflection_available=(bool(candidates) if receivers['auto_reflection_line'] else
                          illuminated_hit(hit_x) and bool(counts) and min(counts)>=receivers['min_reflected_stations'])
    if joint_auto:
        counts=[len(automatic_lines.get(e[0],{}).get('records',[])) for e in reflected]
        reflection_available=any(counts)
    selection=dict(candidate_count=len(candidates),selected_hit_x_m=hit_x if reflection_available else None,
                   automatic=receivers['auto_reflection_line'],valid_station_counts=dict(zip([e[0] for e in reflected],counts)),
                   clearance_m=clearance,status='selected' if reflection_available else 'no_feasible_reflection_line',mode=isolation,
                   reason='按用户选择执行空间隔离或理论波包时间间隔筛选；不保证消除实际混波')
    if joint_auto:
        selection.update(strategy='common_intersection',line_selection=line_reports,
                         candidate_count=max(r['feasible_segments'] for r in line_reports.values()),
                         selected_hit_x_m=hit_x,common_intersection_x_depth_m=receivers['common_intersection_x_depth_m'],
                         interface_hit_x_by_phase_m=receivers['interface_hit_x_by_phase_m'],
                         status='selected' if len(automatic_lines)==len(phase_entries) else 'partial',
                         reason='只选一个公共交点，只有一条入射线；所有测线延长线共点。联合筛选连续合格段，核对入射/反射竞争震相的几何时间差。')
    if not reflection_available:print('当前隔离设置没有可行反射测线；未生成反射组。')
    for phase,layer,direction,depth,multiplier in phase_entries:
        phase_hit=automatic_lines[phase]['hit_x_m'] if joint_auto and phase in automatic_lines else hit_x
        anchor=np.array([phase_hit,depth]);valid=[]
        phase_distances=([r['distance_m'] for r in automatic_lines.get(phase,{}).get('records',[])] if joint_auto
                         else samples(receivers['ray_distance_min_m'],receivers['ray_distance_max_m'],receivers['ray_step_m']))
        for j,distance in enumerate(phase_distances,1):
            pos=anchor+multiplier*distance*direction
            inside=(pos[1]<dtop-1e-8 if layer=='upper' else pos[1]>dbottom+1e-8)
            if not inside:continue
            if not joint_auto and phase.startswith('R_') and (not reflection_available or not isolated(pos,(phase,layer,direction,depth,multiplier))):continue
            if phase.startswith('I_') and (float(normal@(pos-ref))<0 or beam_gap(pos)>1e-8):continue
            before=len(groups);add(*pos,'reference_'+phase,j)
            if len(groups)>before:
                valid.append(pos.tolist());st=groups[-1]['station']
                # Phase distance / speed plus p*X phase must be removed in spectral comparison.
                direction_for_phase=direction
                phase_delay=p*(pos[0]-phase_hit)+(direction_for_phase[1]/speed(model[layer],phase.split('_',1)[1]))*(pos[1]-depth)
                distance_rows.append(dict(group='reference_'+phase,station=st,distance_to_reference_interface_m=float(distance),
                                          reference_interface_depth_m=depth,reference_interface_x_m=phase_hit,
                                          harmonic_phase_delay_s=float(phase_delay),layer=layer,incident_beam_gap_m=beam_gap(pos)))
        rays.append(dict(phase=phase,layer=layer,interface_depth_m=depth,interface_hit_x_m=phase_hit,points=valid))
    if not station_list or not distance_rows:raise ValueError('研究区域内没有有效测线台站，请调整交点/距离/角度')
    window_rows=[]
    for ray in rays:
        if not ray['phase'].startswith('R_'):continue
        entry=next(e for e in phase_entries if e[0]==ray['phase'])
        for group in groups:
            if group['group']!='reference_'+ray['phase']:continue
            pos=np.array([group['x_m'],group['burial_m']]);ti,tr=timing(pos,entry)
            start,end=tr-half,tr+half
            if beam_gap(pos)<=clearance:start=max(start,ti+half+guard)
            # 也不跨过其他传播反射模式到时的中点。
            for other in reflected:
                if other[0]==entry[0]:continue
                _,to=timing(pos,other)
                if to<tr:start=max(start,(to+tr)/2)
                elif to>tr:end=min(end,(to+tr)/2)
            window_rows.append(dict(family='reference',phase=ray['phase'],station=group['station'],start_s=start,end_s=end,
                use=1,reviewed=0,incident_arrival_s=ti,reflected_arrival_s=tr,separation_s=tr-ti,
                incident_beam_gap_m=beam_gap(pos),note='几何相位到时建议窗；未考虑薄层群延迟/绕射/多次反射，必须结合实际波形核对'))
    out.mkdir(parents=True,exist_ok=True)
    failure_note=out/'line_selection_failed.json'
    if failure_note.is_file():failure_note.unlink()  # 本次成功时移除本程序上次留下的失败标记
    if joint_auto:
        audit=[dict(phase=ph,interface_hit_x_m=item['hit_x_m'],x_m=r['point'][0],burial_m=r['point'][1],
                    **{k:v for k,v in r.items() if k!='point'}) for ph,item in automatic_lines.items() for r in item['records']]
        write_csv(out/'automatic_line_separation_audit.csv',list(audit[0]) if audit else ['phase'],audit)
    write_csv(out/'reflection_window_suggestions.csv',list(window_rows[0]) if window_rows else ['family','phase','station','start_s','end_s','use','reviewed'],window_rows)
    if any(dtop-1e-8<=s['burial_m']<=dbottom+1e-8 for s in station_list):
        raise ValueError('Receiver placed in weak layer/interface: prohibited by case design')
    with (out/'STATIONS').open('w',encoding='ascii') as f:
        for s in station_list:f.write(f"{s['station']} DB {s['y_m']:.9f} {s['x_m']:.9f} 0.0 {s['burial_m']:.9f}\n")
    line_names={g['station'] for g in groups if g['group'].startswith('reference_')}
    with (out/'STATIONS_lines').open('w',encoding='ascii') as f:
        for s in station_list:
            if s['station'] in line_names:f.write(f"{s['station']} DB {s['y_m']:.9f} {s['x_m']:.9f} 0.0 {s['burial_m']:.9f}\n")
    write_csv(out/'station_coordinates.csv',list(station_list[0]),station_list)
    write_csv(out/'station_groups.csv',list(groups[0]),groups)
    write_csv(out/'ray_station_distances.csv',list(distance_rows[0]),distance_rows)
    # 真正激发波的有限面力源；位置为全局XY及相对外顶面埋深(km)。
    with (out/'FORCESOLUTION').open('w',encoding='ascii') as f:
        for row in forces:
            f.write(f"FORCE {row['force_id']:06d}\ntime shift: {row['time_shift_s']:.12g}\nhdurorf0: {incident['frequency_hz']:.12g}\n"
                    f"latorUTM: {row['y_m']:.9f}\nlongorUTM: {row['x_m']:.9f}\ndepth: {row['burial_m']/1000:.12f}\n"
                    f"source time function: 1\nfactor force source: {incident['force_factor']:.12e}\n"
                    f"component dir vect source E: {row['E']:.12e}\ncomponent dir vect source N: {row['N']:.12e}\n"
                    f"component dir vect source Z_UP: {row['Z_UP']:.12e}\n")
    write_csv(out/'sources_depth.csv',list(forces[0]),forces)
    write_csv(out/'source_centers.csv',list(source_centers[0]),source_centers)
    materials(model,out)
    requirements='''# Finite source plane extended in X; every force point >=100m from all outer boundaries.
# Use the standard, unpatched SPECFEM3D solver.
SIMULATION_TYPE = 1
GPU_MODE = .false.
USE_LDDRK = .false.
STACEY_ABSORBING_CONDITIONS = .true.
STACEY_INSTEAD_OF_FREE_SURFACE = .true.
BOTTOM_FREE_SURFACE = .false.
PML_CONDITIONS = .false.
COUPLE_WITH_INJECTION_TECHNIQUE = .false.
ATTENUATION = .false.
ANISOTROPY = .false.
USE_FORCE_POINT_SOURCE = .true.
USE_RICKER_TIME_FUNCTION = .true.
DIPOLE_SOURCE_IN_FLUID = .false.
USE_EXTERNAL_SOURCE_FILE = .false.
SUPPRESS_UTM_PROJECTION = .true.
USE_SOURCES_RECEIVERS_Z = .false.
SAVE_SEISMOGRAMS_DISPLACEMENT = .true.
SAVE_SEISMOGRAMS_VELOCITY = .true.
SAVE_SEISMOGRAMS_ACCELERATION = .true.
'''
    (out/'Par_file_required.txt').write_text(requirements,encoding='ascii')
    front_hit=next((automatic_lines[e[0]]['hit_x_m'] for e in reflected if e[0] in automatic_lines),hit_x)
    dt=float(slow@(np.array([front_hit,inc_interface])-ref)) if WAVEFRONT_AFTER_REFERENCE_S is None else float(WAVEFRONT_AFTER_REFERENCE_S)
    if not np.isfinite(dt):raise ValueError('波峰线时刻必须有限')
    # 有限源面的波峰示意：平移原始有限线段后再裁到研究盒/入射介质，不延长线段。
    def finite_front(offset):
        a,b=np.array(source_line)+offset*normal;vector=b-a
        limits=[(x0,x1),(d0,dtop if side=='upper' else d1)] if side=='upper' else [(x0,x1),(dbottom,d1)]
        tmin,tmax=0.,1.
        for axis,(lo,hi) in enumerate(limits):
            if abs(vector[axis])<1e-12:
                if not lo<=a[axis]<=hi:return []
            else:
                u,w=sorted(((lo-a[axis])/vector[axis],(hi-a[axis])/vector[axis]))
                tmin=max(tmin,u);tmax=min(tmax,w)
        return [(a+tmin*vector).tolist(),(a+tmax*vector).tolist()] if tmin<tmax else []
    fronts=[dict(label='实际有限震源面投影',time_s=incident['peak_time_s'],offset_m=0.0,points=source_line),
            dict(label='有限入射波峰线示意',time_s=incident['peak_time_s']+dt,offset_m=v*dt,points=finite_front(v*dt))]
    coordinates=dict(outer_xyz_m=outer,research_xyz_m=inner,research_x_depth_bounds_m=list(research),
                     convention='Input positions relative to research box: X from xmin, Y from ymin, depth downward from its top. Files use global X/Y and depth below outer top.',
                     local_to_global_offset_x_y_depth_m=[inner['x'][0],inner['y'][0],d0],
                     absolute_z_formula=f'Z={total:g}-burial_depth',wavefront_equation='normal dot ([X,depth]-reference) = incident_speed * (time-peak_time)')
    report=dict(schema='three_layer_v1',reflection_selection=selection,coordinate_mapping=coordinates,local_inputs=local_inputs,wavefronts=fronts,model=model,incident=incident,receivers=receivers,interface_depths_m=[dtop,dbottom],
                total_depth_m=total,source_mode='finite_plane_force_couples' if kind(inc)=='solid' else 'finite_plane_acoustic_monopoles',
                point_source_count=len(forces),active_point_source_count=len(forces),source_center_count=len(source_centers),source_boundary_audit=boundary_audit,
                source_line_x_depth_m=source_line,source_y_range_m=[min(r['y_m'] for r in forces),max(r['y_m'] for r in forces)],
                source_reference_x_y_depth_m=point,wavefront_reference_interface_x_m=front_hit,station_count=len(station_list),rays=rays,
                note='Reference lines are phase-sampling lines, not unique transmitted rays through a multiply reflecting thin layer. Compare complex spectra after depropagation, not a single pulse peak.',
                solver_status='Finite aperture point sources only; no analytic boundary injection. Forward/backward waves and edge diffraction exist; wavefront lines are geometric guides, not exact isolated one-way plane waves.')
    (out/'config_used.json').write_text(json.dumps(dict(model=model,incident=incident,receivers=receivers,local_inputs=local_inputs,coordinate_mapping=coordinates),ensure_ascii=False,indent=2),encoding='utf-8')
    (out/'geometry_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    from matplotlib.patches import Rectangle
    fig,axes=plt.subplots(1,2,figsize=(17,8.5),layout='constrained')
    for ax,zoom in zip(axes,(False,True)):
        ax.axhspan(0,total,color='#efefef',label='外部延拓/粗网格区域')
        ax.add_patch(Rectangle((x0,d0),x1-x0,d1-d0,facecolor='#eaf5fc',edgecolor='#34536a',lw=2,label='研究小长方体'))
        beam_depths=np.linspace(d0,dtop,300) if side=='upper' else np.linspace(dbottom,d1,300)
        beam_left=np.maximum(x0,(beam_lo-tangent[1]*beam_depths)/tangent[0]) if tangent[0]>0 else np.maximum(x0,(beam_hi-tangent[1]*beam_depths)/tangent[0])
        beam_right=np.minimum(x1,(beam_hi-tangent[1]*beam_depths)/tangent[0]) if tangent[0]>0 else np.minimum(x1,(beam_lo-tangent[1]*beam_depths)/tangent[0])
        ax.fill_betweenx(beam_depths,beam_left,beam_right,where=beam_right>=beam_left,color='#f5ad77',alpha=.35,label='入射波束覆盖区（需时间分离）' if isolation=='time' else '入射波束覆盖区（反射站禁入）')
        ax.axhspan(dtop,dbottom,color='#d99e23',label=f"{model['middle_thickness_m']:g}m软弱层")
        ax.axhline(dtop,color='#b68116',lw=1);ax.axhline(dbottom,color='#b68116',lw=1)
        xy=np.array([[s['x_m'],s['burial_m']] for s in station_list]);ax.scatter(xy[:,0],xy[:,1],s=.8,c='#b3c1cc',alpha=.5,label='研究区台站')
        for ray in rays:
            if ray['points']:
                xy=np.array(ray['points'])
                info=line_reports.get(ray['phase'],{})
                label=ray['phase']+'测线'
                if info.get('status')=='selected':label+=f" {info['distance_start_m']:g}–{info['distance_end_m']:g}m"
                line,=ax.plot(xy[:,0],xy[:,1],'.-',ms=2.3,lw=1.8,label=label)
                anchor=(receivers['common_intersection_x_depth_m'] if joint_auto else [ray['interface_hit_x_m'],ray['interface_depth_m']])
                ax.plot([anchor[0],xy[0,0]],[anchor[1],xy[0,1]],ls=':',lw=1,color=line.get_color())
                if not joint_auto:ax.scatter(*anchor,marker='D',s=20,color=line.get_color(),zorder=5)
                offsets={'I_P':(-22,14),'I_SV':(-22,14),'R_P':(40,12),'R_SV':(-40,35),'T_P':(-35,-32),'T_SV':(20,-48)}
                if zoom and not joint_auto:ax.annotate(ray['phase']+f"\nX={anchor[0]:g}",anchor,xytext=offsets[ray['phase']],
                                    textcoords='offset points',fontsize=7,ha='center',color=line.get_color(),
                                    arrowprops=dict(arrowstyle='-',color=line.get_color(),lw=.7))
        if joint_auto:
            common=receivers['common_intersection_x_depth_m']
            ax.scatter(*common,marker='*',s=140,color='black',zorder=7,label='唯一公共交点')
            if zoom:ax.annotate(f'公共交点\nX={common[0]:g}m，埋深={common[1]:g}m',common,xytext=(-115,-65),
                                textcoords='offset points',fontsize=9,bbox=dict(facecolor='white',alpha=.85,edgecolor='gray'),
                                arrowprops=dict(arrowstyle='->',color='black'))
        box=(x0,x1,d0,dtop) if side=='upper' else (x0,x1,dbottom,d1)
        for front,color,style in zip(fronts,('#962e91','#d54032'),('--','-')):
            points=front['points']
            if len(points)==2:
                xy=np.array(points);ax.plot(xy[:,0],xy[:,1],color=color,ls=style,lw=2.6,
                    label=front['label']+f" t={front['time_s']:.4f}s")
        arrow_start=ref+normal*70;arrow_end=arrow_start+normal*100
        ax.annotate('',xy=arrow_end,xytext=arrow_start,arrowprops=dict(arrowstyle='-|>',lw=2,color='#962e91'))
        ax.text(arrow_start[0]+15,arrow_start[1],'入射传播方向',fontsize=8,color='#962e91')
        ax.scatter(reference_x,reference_depth,marker='x',c='k',s=60,label='有限源面参考点')
        ax.set(xlim=(x0-30,x1+30) if zoom else (0,lx),ylim=(d1+30,d0-30) if zoom else (total,0),
               xlabel='全局 X / m',ylabel='相对外模型顶面的埋深 / m（向下为正）',
               title=('研究区放大' if zoom else '完整外模型与研究区')+f'；全局 Y={y:g}m')
        ax.set_aspect('equal');ax.legend(fontsize=7,loc='best',ncol=2);ax.grid(alpha=.15)
    fig.suptitle(f'三层嵌套模型：{side}侧 {incident["wave"]} 入射，角度{incident["angle_deg"]:g}°\n'
                 f'源面左端X={source_line[0][0]:g}m，右端X={source_line[1][0]:g}m；'+
                 (('所有测线延长线共用一个界面交点；实线为连续合格台站段' if joint_auto else
                   '反射台站采用时间窗分离' if isolation=='time' else '反射台站严格在入射波束外') if reflection_available else '当前隔离条件无可行反射组'),fontsize=14)
    fig.savefig(out/'stations_three_layers.png',dpi=170);plt.close(fig)
    if joint_auto:
        fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
        required=2*half+guard
        for ph,item in automatic_lines.items():
            if ph.startswith('T_'):continue
            records=item['records'];ax=axes[0 if ph.startswith('I_') else 1]
            finite=[r for r in records if r['min_competing_arrival_separation_s'] is not None]
            if finite:
                ax.plot([r['distance_m'] for r in finite],[1000*r['min_competing_arrival_separation_s'] for r in finite],
                        label=ph+f"，交点X={item['hit_x_m']:g}m")
        for ax,title in zip(axes,('入射站：与竞争反射震相的最小时间差','反射站：与入射及其他反射震相的最小时间差')):
            ax.axhline(1000*required,color='#b23a32',ls='--',label=f'要求 ≥ {1000*required:.2f} ms')
            ax.set(title=title,xlabel='沿本测线到界面的距离 / m',ylabel='几何时间差 / ms')
            ax.legend(fontsize=8);ax.grid(alpha=.2)
        fig.suptitle('自动选线隔离核查｜有限波束外不存在的竞争几何射线不计入\n'
                     '按设定波包宽度筛选；薄层群延迟、有限源绕射与三维散射仍需实际波形复核',fontsize=11)
        fig.savefig(out/'line_separation_check.png',dpi=170);plt.close(fig)
    print('Prepared:',out,'stations:',len(station_list),'active point sources:',len(forces))
    print('实际力源点距六个外边界的最小净距(m):',boundary_audit['minimum_clearance_m'])
    return report


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path)
    args=parser.parse_args();run(out=args.output)
