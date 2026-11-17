#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""独立三层反透射系数对比程序：读取理论表、台站元数据和数值波形。"""
from pathlib import Path
import argparse
import copy
import csv
import json
import math
from collections import Counter
import tempfile
import numpy as np
from scipy.signal import butter, detrend, hilbert, sosfiltfilt
from scipy.signal.windows import tukey
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ==================== 用户编辑区：只修改本文件 ====================
# 本程序不导入理论解或生成器代码；只读取它们输出的数据文件。
CASE_NAME = 'Planewave_lower_0_P_95Hz'      # 换算例时只改这里；以下三个路径自动确定
# 输入=result/CASE_NAME/prepare，波形=result/CASE_NAME/sema，输出=result/CASE_NAME/rt_result
THEORY_DIRECTORY = '../理论解/rt_coefficients_output/lower_P'  # 从prepare元数据自动选择入射侧/波型、角度、频率
RUN_MODE = 'analyze'                       # 直接运行即可；prepare只检查输入并输出XLSX计划，无需预先运行

TRACE_QUANTITY = 'acceleration'            # acceleration / velocity / displacement
FILE_EXTENSION = 'sema'                    # 对应sema / semv / semd
COMPONENT_CODES = {'X':'FXX','Y':'FXY','Z':'FXZ'}
Z_COMPONENT_POSITIVE = 'up'                # 文件Z向上为正；切换入射侧无需改动
INCLUDE_Y_IN_AMPLITUDE = True              # True：目标P/SV投影与Y合成；False：恢复纯P/SV投影
# 合成直接幅值=sqrt(s²+Y²)；合成包络=sqrt(|analytic(s)|²+|analytic(Y)|²)。缺Y时明确报错，不补零。
AMPLITUDE_SCALE = 1.0                      # 同一单位换算，不逐道归一化
DETREND = 'none'                           # none / constant / linear
BANDPASS_HZ = [70.0,120.0]                 # 与二层相同；None关闭，按实际源主频调整
FILTER_ORDER = 4

STATION_SELECTION_MODE = 'auto_segment'    # auto_segment自动最佳连续段 / manual手动距离范围
AUTO_SEGMENT_LENGTH_M = 100.0              # 优先连续100m；仅允许下方开关按新prepare明确接受的短段缩短
ALLOW_PREPARED_SHORT_SEGMENTS = True       # 仅当prepare明确接受短测线时，按其实际共同可用长度缩短；旧算例不变
AUTO_SEGMENT_MIN_VALID_FRACTION = 0.8       # 每段至少80%台站能计算；缺失/被排除台站也计入分母
# 自动模式扫描整条已有测线；下表只在manual模式生效。距离均沿测线从对应界面量起。
STATION_DISTANCE_RANGES_M = {
    'I_P':[100.0,250.0], 'R_P':[100.0,300.0], 'T_P':[100.0,300.0],
    'I_SV':[100.0,250.0], 'R_SV':[100.0,300.0], 'T_SV':[100.0,300.0],
}
PAIRING_RULE = 'all_pairs'                 # all_pairs / equal_distance / equal_travel_time
PAIRING_TOLERANCE_M = 3.0
WINDOW_HALF_PERIODS = 5.0                  # 半窗宽=本值/源主频
WINDOW_BEFORE_S = None                     # 可单独给秒数
WINDOW_AFTER_S = None
USE_GENERATOR_REFLECTION_WINDOWS = True    # 读取生成器隔离窗；新共同交点布局同时作用于入射/反射，旧布局仅反射
WINDOW_CENTER_MODE = 'theory'              # theory / waveform（理论到时附近搜索）
ARRIVAL_SEARCH_HALF_WIDTH_S = None         # None=一个主频周期
PHASE_CENTER_OFFSETS_S = {                 # 人工核对波形后可修正各相位中心，秒
    'I_P':0.0,'I_SV':0.0,'R_P':0.0,'R_SV':0.0,'T_P':0.0,'T_SV':0.0,
}
USE_THEORY_GROUP_DELAY = True              # 有同角度频谱表时，计入薄层的相位群延时
LIMIT_WINDOW_TO_TARGET_PHASE = True        # 到时中点限制拾取窗，防止拾取其他波
WINDOW_OVERRIDES_CSV = None                # family,phase,station,start_s,end_s,use,reviewed
TAPER_FRACTION = 0.10
HILBERT_PAD_FACTOR = 4
MIN_WINDOW_SAMPLES = 16
MIN_SAMPLES_PER_PERIOD = 10.0
EDGE_AMPLITUDE_LIMIT = 0.20
BEST_PAIR_METHOD = 'hilbert'               # hilbert / direct，沿用二层最佳组合定义
BEST_PAIR_ACCEPTED_ONLY = False            # False允许有混波等提示的可计算台站，但在XLSX/图中明确标记
# True只使用通过质量检查/人工复核的台站；不足150m时明确报告，不自动放宽。
# 自动段评分=段内全部入射×出射组合的系数误差RMSE；×仅指配对，系数始终是出射幅值÷入射幅值。
MIN_REFERENCE_AMPLITUDE = 1e-30
PLOT_DPI = 170
# ==================== 用户编辑区结束 ====================

ROOT=Path(__file__).resolve().parent


def settings():
    names=('CASE_NAME THEORY_DIRECTORY RUN_MODE TRACE_QUANTITY '
           'STATION_SELECTION_MODE AUTO_SEGMENT_LENGTH_M AUTO_SEGMENT_MIN_VALID_FRACTION ALLOW_PREPARED_SHORT_SEGMENTS '
           'FILE_EXTENSION COMPONENT_CODES Z_COMPONENT_POSITIVE INCLUDE_Y_IN_AMPLITUDE AMPLITUDE_SCALE DETREND BANDPASS_HZ FILTER_ORDER '
           'STATION_DISTANCE_RANGES_M PAIRING_RULE PAIRING_TOLERANCE_M WINDOW_HALF_PERIODS WINDOW_BEFORE_S '
           'WINDOW_AFTER_S USE_GENERATOR_REFLECTION_WINDOWS WINDOW_CENTER_MODE ARRIVAL_SEARCH_HALF_WIDTH_S PHASE_CENTER_OFFSETS_S '
           'USE_THEORY_GROUP_DELAY LIMIT_WINDOW_TO_TARGET_PHASE WINDOW_OVERRIDES_CSV TAPER_FRACTION '
           'HILBERT_PAD_FACTOR MIN_WINDOW_SAMPLES MIN_SAMPLES_PER_PERIOD EDGE_AMPLITUDE_LIMIT '
           'BEST_PAIR_METHOD BEST_PAIR_ACCEPTED_ONLY MIN_REFERENCE_AMPLITUDE PLOT_DPI').split()
    o={name.lower():copy.deepcopy(globals()[name]) for name in names}
    return case_paths(o)


def case_paths(o):
    name=o['case_name']
    if not name or name in ('.','..') or any(ch in name for ch in '/\\:'):
        raise ValueError('CASE_NAME只填写算例文件夹名称，不要填写路径')
    case=ROOT.parent/'result'/name
    o.update(setup_directory=str(case/'prepare'),data_directory=str(case/'sema'),output_directory=str(case/'rt_result'))
    return o


def resolve(value):
    p=Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (ROOT/p).resolve()


def read_csv(path,delimiter=','):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter=delimiter))


def write_csv(path,rows,fields=None):
    fields=fields or (list(rows[0]) if rows else ['status'])
    with Path(path).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields,extrasaction='ignore');w.writeheader();w.writerows(rows)


def medium_kind(m):
    return 'poroelastic' if m.get('poroelastic_enabled',False) else m['type']


def check_inputs(o):
    if o['station_selection_mode'] not in ('auto_segment','manual'):raise ValueError('台站选择模式应为auto_segment/manual')
    if not np.isfinite(o['auto_segment_length_m']) or o['auto_segment_length_m']<=0:raise ValueError('自动测线段长度必须为正数')
    if not 0<o['auto_segment_min_valid_fraction']<=1:raise ValueError('有效台站比例必须在(0,1]内')
    if o['run_mode'] not in ('prepare','analyze'):raise ValueError('RUN_MODE应为prepare/analyze')
    if {'acceleration':'sema','velocity':'semv','displacement':'semd'}.get(o['trace_quantity'])!=o['file_extension']:
        raise ValueError('物理量和文件扩展名不匹配')
    if o['z_component_positive'] not in ('up','down'):raise ValueError('Z_COMPONENT_POSITIVE应为up/down')
    if o['window_center_mode'] not in ('theory','waveform') or o['best_pair_method'] not in ('direct','hilbert'):
        raise ValueError('选窗/最佳组合方法无效')
    if o['pairing_rule'] not in ('all_pairs','equal_distance','equal_travel_time'):raise ValueError('配对方式无效')
    if o['detrend'] not in ('none','constant','linear'):raise ValueError('DETREND应为none/constant/linear')
    for key in ('amplitude_scale','window_half_periods','hilbert_pad_factor','min_samples_per_period','min_reference_amplitude'):
        if not np.isfinite(o[key]) or o[key]<=0:raise ValueError(key+'必须为有限正数')
    if not 0<=o['taper_fraction']<=1 or o['min_window_samples']<4:raise ValueError('渐消或最少采样点数无效')
    if o['filter_order']<1 or o['pairing_tolerance_m']<0:raise ValueError('滤波阶数或配对容差无效')
    for bounds in o['station_distance_ranges_m'].values():
        if bounds is not None and (len(bounds)!=2 or not np.isfinite(bounds).all() or not 0<=bounds[0]<=bounds[1]):raise ValueError('台站距离范围无效')


def load_context(o):
    folder=resolve(o['setup_directory'])
    c=json.loads((folder/'config_used.json').read_text(encoding='utf-8-sig'))
    rs=c.get('receivers',{})
    if (o.get('allow_prepared_short_segments',True) and o['station_selection_mode']=='auto_segment'
            and rs.get('short_clean_lines_accepted',False)):
        length=float(rs['available_common_segment_length_m'])
        if not np.isfinite(length) or length<=0:raise ValueError('prepare中的可用短段长度无效')
        if length<o['auto_segment_length_m']:
            o['requested_auto_segment_length_m']=o['auto_segment_length_m']
            o['auto_segment_length_m']=length
            print(f"本prepare明确接受短测线：后处理选段长度从{o['requested_auto_segment_length_m']:g}m改为{length:g}m。",flush=True)
    if not all(k in c for k in ('model','incident')):raise ValueError('需要本次三层生成器输出的元数据，不能用二层数据')
    model,incident=c['model'],c['incident'];side=incident['side']
    a=model['upper_thickness_m'];b=a+model['middle_thickness_m']
    inc=model[side];wave=incident['wave'];v=inc['vp'] if wave=='P' else inc['vs']
    sign=1 if side=='upper' else -1
    p=math.sin(math.radians(incident['angle_deg']))/v
    theory=resolve(o['theory_directory'])
    if not (theory/'model_parameters.json').is_file():theory=theory/(side+'_'+wave)
    snapshot=json.loads((theory/'model_parameters.json').read_text(encoding='utf-8-sig'))
    tm,ti=snapshot['model'],snapshot['incident']
    if (ti['side'],ti['wave'])!=(side,wave):raise ValueError('所选理论解入射侧/波型与数值算例不一致')
    # 理论求解器的外侧为半空间，系数在两界面定义，仅薄层厚度参与传播相位。
    # 实际到时仍使用本算例prepare里的全局界面埋深，不沿用理论表的旧坐标。
    for key in ('middle_thickness_m',):
        if not np.isclose(model[key],tm[key],rtol=1e-9,atol=1e-9):raise ValueError('理论解/数值算例厚度不一致：'+key)
    for layer in ('upper','middle','lower'):
        m,n=model[layer],tm[layer]
        if medium_kind(m)!=medium_kind(n):raise ValueError('介质类型不一致：'+layer)
        props=m['biot'] if medium_kind(m)=='poroelastic' else {k:m[k] for k in ('rho','vp','vs')}
        other=n['biot'] if medium_kind(n)=='poroelastic' else n
        for key,value in props.items():
            if not np.isclose(value,other[key],rtol=1e-9,atol=1e-12):raise ValueError('理论解/数值材料参数不一致：'+layer+'.'+key)
    if medium_kind(model['middle'])=='poroelastic':
        for key in ('upper_middle_pore_condition','middle_lower_pore_condition'):
            if model[key]!=tm[key]:raise ValueError('多孔界面条件不一致：'+key)
    frequency=float(incident['frequency_hz'])
    table=theory/f'coefficients_f{frequency:g}Hz.tsv'
    rows=read_csv(table,'\t')
    matched=[r for r in rows if np.isclose(float(r['angle_deg']),incident['angle_deg'],rtol=0,atol=1e-7)
             and np.isclose(float(r['frequency_hz']),frequency,rtol=0,atol=1e-7)]
    if len(matched)!=1:raise ValueError('理论表没有唯一的实际入射角/主频记录，请重新导出，不能使用最近角度代替')
    target=matched[0];coeff={};delays={}
    for name in target:
        if name.startswith(('R_','T_')) and name.endswith('_abs'):
            phase=name[:-4];value=float(target[name]);z=complex(float(target[phase+'_real']),float(target[phase+'_imag']))
            if not np.isfinite([z.real,z.imag,value]).all() or value<0 or not np.isclose(value,abs(z),rtol=1e-8,atol=1e-10):raise ValueError('理论复系数/模不一致：'+phase)
            coeff[phase]=value
    spectrum=theory/'coefficients_frequency_spectrum.tsv'
    if o['use_theory_group_delay'] and spectrum.is_file():
        spectrum_rows=read_csv(spectrum,'\t')
        same=[r for r in spectrum_rows if np.isclose(float(r['angle_deg']),incident['angle_deg'],atol=1e-7,rtol=0)]
        same.sort(key=lambda r:float(r['frequency_hz']))
        if len(same)>=3:
            f=np.array([float(r['frequency_hz']) for r in same])
            if np.all(np.diff(f)>0) and f[0]<=frequency<=f[-1]:
                for phase in coeff:
                    z=np.array([complex(float(r[phase+'_real']),float(r[phase+'_imag'])) for r in same])
                    derivative=np.gradient(np.unwrap(np.angle(z)),2*np.pi*f)
                    delay=float(np.interp(frequency,f,derivative))
                    if np.isfinite(delay):delays[phase]=delay
    phases={}
    for phase in ['I_'+wave]+list(coeff):
        layer=side if phase.startswith(('I_','R_')) else ('lower' if side=='upper' else 'upper')
        label=phase.split('_',1)[1];medium=model[layer];vm=medium['vp'] if label=='P' else medium['vs']
        sine=p*vm;dsign=sign if not phase.startswith('R_') else -sign
        cosine=math.sqrt(max(0,1-sine*sine));propagating=sine<1-1e-10
        direction=np.array([sine,dsign*cosine])
        phases[phase]=dict(wave=label,layer=layer,speed=vm,direction=direction,propagating=propagating,
                           theory_abs=coeff.get(phase,1.0),group_delay_s=delays.get(phase,0.0))
    group_rows=read_csv(folder/'station_groups.csv');stations={r['station']:r for r in read_csv(folder/'station_coordinates.csv')}
    group_rows=[r for r in group_rows if r['group'].startswith('reference_')]
    if not group_rows:raise ValueError('没有reference参考测线台站')
    context=dict(model=model,incident=incident,phases=phases,groups=group_rows,stations=stations,
                 geometry=json.loads((folder/'geometry_report.json').read_text(encoding='utf-8-sig')) if (folder/'geometry_report.json').exists() else {},
                 receiver_settings=c['receivers'],
                 beam_gaps={(r['group'][len('reference_'):],r['station']):float(r['incident_beam_gap_m'])
                            for r in read_csv(folder/'ray_station_distances.csv') if r.get('incident_beam_gap_m')},
                 reference_hit_x=c['receivers']['interface_hit_x_m'],
                 top=a,bottom=b,p=p,sign=sign,incident_speed=v,theory_table=str(table),theory_snapshot=str(theory/'model_parameters.json'))
    return context


def arrival(c,phase,point,o):
    info=c['phases'][phase]
    if not info['propagating']:return None
    layer='upper' if point[1]<c['top']-1e-8 else 'lower' if point[1]>c['bottom']+1e-8 else 'middle'
    if layer!=info['layer']:return None
    inc=c['incident'];refdepth=inc['reference_depth_m'][inc['side']]
    nd=c['sign']*math.cos(math.radians(inc['angle_deg']))/c['incident_speed']
    if phase.startswith('I_'):
        result=inc['peak_time_s']+c['p']*(point[0]-inc['reference_x_m'])+nd*(point[1]-refdepth)
    else:
        near=c['top'] if inc['side']=='upper' else c['bottom']
        interface=near if phase.startswith('R_') else c['bottom'] if info['layer']=='lower' else c['top']
        result=inc['peak_time_s']+nd*(near-refdepth)+c['p']*(point[0]-inc['reference_x_m'])
        result+=info['direction'][1]/info['speed']*(point[1]-interface)+info['group_delay_s']
    return float(result+o['phase_center_offsets_s'].get(phase,0.0))


def make_plan(c,o):
    overrides={}
    if o['window_overrides_csv']:
        for r in read_csv(resolve(o['window_overrides_csv'])):
            key=(r['phase'],r['station'])
            if key in overrides:raise ValueError('重复的人工窗：'+str(key))
            overrides[key]=r
    plan=[];inc=c['incident'];f0=inc['frequency_hz']
    before=o['window_before_s'] if o['window_before_s'] is not None else o['window_half_periods']/f0
    after=o['window_after_s'] if o['window_after_s'] is not None else o['window_half_periods']/f0
    if not np.isfinite([before,after]).all() or min(before,after)<=0:raise ValueError('窗宽必须为正数')
    for group in c['groups']:
        phase=group['group'][len('reference_'):]
        if phase not in c['phases']:raise ValueError('元数据包含未知相位：'+phase)
        info=c['phases'][phase]
        if not info['propagating']:continue
        name=group['station'];st=c['stations'][name]
        point=np.array([float(st['x_m']),float(st['burial_m'])])
        depth=(c['top'] if inc['side']=='upper' else c['bottom']) if phase.startswith(('I_','R_')) else c['bottom'] if info['layer']=='lower' else c['top']
        # 各线延长线共享近侧交点；透射线在远侧界面的起点有薄层几何偏移。
        x0=c['receiver_settings'].get('interface_hit_x_by_phase_m',{}).get(phase,c['reference_hit_x'])
        distance=float(np.linalg.norm(point-[x0,depth]))
        bounds=o['station_distance_ranges_m'].get(phase) if o['station_selection_mode']=='manual' else None
        if bounds is not None and not bounds[0]-1e-6<=distance<=bounds[1]+1e-6:continue
        center=arrival(c,phase,point,o)
        if center is None:continue
        others={ph:value for ph in c['phases'] if ph!=phase and (value:=arrival(c,ph,point,o)) is not None}
        coincident=[ph for ph,value in others.items() if abs(value-center)<1e-8]
        lower=max(((center+t)/2 for t in others.values() if t<center-1e-8),default=-np.inf)
        upper=min(((center+t)/2 for t in others.values() if t>center+1e-8),default=np.inf)
        target_before,target_after=before,after
        settings=c.get('receiver_settings',{})
        phase_half=settings.get('isolated_phase_half_periods',{}).get(phase)
        if o.get('use_generator_reflection_windows',True) and phase_half is not None:
            if not np.isfinite(phase_half) or phase_half<=0:raise ValueError('prepare中的震相隔离半窗宽无效')
            target_before=target_after=float(phase_half)/f0
        time_isolation=(o.get('use_generator_reflection_windows',True) and phase.startswith('R_')
                        and settings.get('reflection_isolation')=='time')
        isolation_start=-np.inf
        if time_isolation:
            target_before=target_after=float(settings['isolation_half_periods'])/f0
            gap=c.get('beam_gaps',{}).get((phase,name),-np.inf)
            if gap<=settings.get('incident_beam_clearance_m',0.0):
                ti=others.get('I_'+inc['wave'])
                if ti is not None:isolation_start=ti+target_before+float(settings['isolation_guard_s'])
            lower=max(lower,isolation_start)
        start,end=max(center-target_before,isolation_start),center+target_after
        if o['limit_window_to_target_phase']:start,end=max(start,lower),min(end,upper)
        row=dict(family='reference',phase=phase,station=name,x_m=point[0],burial_m=point[1],
                 distance_to_interface_m=distance,predicted_peak_s=center,start_s=start,end_s=end,
                 gate_start_s=lower,gate_end_s=upper,other_arrivals_s=json.dumps(others),
                 isolation_start_s=isolation_start,nominal_before_s=target_before,nominal_after_s=target_after,
                 coincident_phases=';'.join(coincident),use=1,reviewed=0,manual_window=0)
        override=overrides.pop((phase,name),{})
        for key in ('start_s','end_s'):
            if override.get(key,'')!='':row[key]=float(override[key]);row['manual_window']=1
        for key in ('use','reviewed'):
            if override.get(key,'')!='':row[key]=int(override[key])
        if row['use'] not in (0,1) or row['reviewed'] not in (0,1):raise ValueError('use/reviewed必须为0或1')
        if row['manual_window'] and o['limit_window_to_target_phase'] and (row['start_s']<lower-1e-9 or row['end_s']>upper+1e-9):
            raise ValueError('人工窗跨越相邻震相边界：'+name)
        if row['start_s']<isolation_start-1e-9:raise ValueError('人工窗进入入射波保护区：'+name)
        if row['end_s']<=row['start_s']:raise ValueError('拾取窗为空：'+name)
        plan.append(row)
    if overrides:raise ValueError('人工窗表中有不属于当前选择的台站/相位')
    if not plan:raise ValueError('没有满足测线/距离条件的台站')
    return plan


def trace_files(folder,o):
    index={}
    for path in folder.rglob('*.'+o['file_extension']):
        fields=path.name.split('.')
        if len(fields)==4:
            key=(fields[1],fields[2])
            if key in index:raise ValueError('台站分量重复：'+str(key))
            index[key]=path
    return index


def read_component(index,name,component,o):
    key=(name,o['component_codes'][component])
    if key not in index:raise FileNotFoundError('缺少台站分量：'+str(key))
    a=np.loadtxt(index[key],ndmin=2)
    if len(a)<4 or a.shape[1]<2 or not np.isfinite(a[:,:2]).all():raise ValueError('波形格式无效：'+str(key))
    t=a[:,0];dt=(t[-1]-t[0])/(len(t)-1)
    if dt<=0 or np.any(np.diff(t)<=0) or np.max(abs(np.diff(t)-dt))>.005*dt:raise ValueError('波形时间不等间隔：'+name)
    uniform=np.linspace(t[0],t[-1],len(t));x=np.interp(uniform,t,a[:,1])*o['amplitude_scale']
    if o['detrend']!='none':x=detrend(x,type=o['detrend'])
    if o['bandpass_hz'] is not None:
        lo,hi=o['bandpass_hz']
        if not 0<lo<hi<.5/dt:raise ValueError('带通频率超出有效范围')
        x=sosfiltfilt(butter(o['filter_order'],[lo,hi],fs=1/dt,btype='bandpass',output='sos'),x)
    return uniform,x


def project(index,name,info,o):
    t,x=read_component(index,name,'X',o);tz,z=read_component(index,name,'Z',o)
    if len(t)!=len(tz) or not np.allclose(t,tz,atol=1e-9,rtol=0):raise ValueError('X/Z分量时间轴不一致')
    if o['z_component_positive']=='up':z=-z
    nx,nd=info['direction'];e=np.array([nx,nd]) if info['wave']=='P' else np.array([-nd,nx])
    y=None
    if (name,o['component_codes']['Y']) in index:
        ty,y=read_component(index,name,'Y',o)
        if len(t)!=len(ty) or not np.allclose(t,ty,atol=1e-9,rtol=0):raise ValueError('Y分量时间轴不一致')
    elif o['include_y_in_amplitude']:raise FileNotFoundError('已启用Y合成但缺少Y分量：'+name)
    return t,e[0]*x+e[1]*z,-e[1]*x+e[0]*z,y


def measure(t,signal,start,end,o,y=None):
    if start<t[0]-1e-8 or end>t[-1]+1e-8:raise ValueError('window_outside_record')
    mask=(t>=start)&(t<=end);tw,x=t[mask],signal[mask]
    if len(x)<o['min_window_samples']:raise ValueError('window_too_short')
    padding=max(len(x),int(len(x)*o['hilbert_pad_factor']/2))
    env=abs(hilbert(np.pad(x*tukey(len(x),alpha=o['taper_fraction']),(padding,padding))))[padding:padding+len(x)]
    direct=abs(x);display=x
    if o['include_y_in_amplitude']:
        if y is None:raise ValueError('已启用Y合成但缺少Y分量')
        if np.shape(y)!=np.shape(signal) or not np.isfinite(y).all():raise ValueError('Y分量长度或数值无效')
        yw=y[mask]
        ey=abs(hilbert(np.pad(yw*tukey(len(x),alpha=o['taper_fraction']),(padding,padding))))[padding:padding+len(x)]
        env=np.hypot(env,ey);direct=np.hypot(x,yw);display=direct
    di,hi=int(np.argmax(direct)),int(np.argmax(env));dp,hp=float(direct[di]),float(env[hi])
    if min(dp,hp)<=o['min_reference_amplitude']:raise ValueError('zero_or_tiny_signal')
    edge=float(max(direct[0],direct[-1])/dp);n=max(2,int(np.ceil(o['taper_fraction']*len(x)/2)))
    flags=[]
    if edge>o['edge_amplitude_limit']:flags.append('window_edge_not_quiet')
    if min(di,hi)<n or max(di,hi)>=len(x)-n:flags.append('peak_at_window_edge')
    return dict(direct_peak=dp,hilbert_peak=hp,direct_peak_time_s=float(tw[di]),hilbert_peak_time_s=float(tw[hi]),
                target_direct_peak=float(max(abs(x))),amplitude_mode='目标投影+Y合成' if o['include_y_in_amplitude'] else '纯目标投影',
                edge_ratio=edge,sample_count=len(x)),flags,(tw,display,env)


def process_row(row,c,index,o):
    result={**row,'usable':0,'direct_peak':'','hilbert_peak':'','status':''}
    if not row['use']:result['status']='excluded:user_disabled';return result,None
    try:
        t,x,orth,y=project(index,row['station'],c['phases'][row['phase']],o)
        start,end=row['start_s'],row['end_s']
        if o['window_center_mode']=='waveform' and not row['manual_window']:
            width=o['arrival_search_half_width_s'] or 1/c['incident']['frequency_hz']
            mask=(t>=max(row['gate_start_s'],row['predicted_peak_s']-width))&(t<=min(row['gate_end_s'],row['predicted_peak_s']+width))
            ids=np.flatnonzero(mask)
            if len(ids)<3:raise ValueError('arrival_search_too_short')
            envelope=abs(hilbert(np.pad(x,(len(x),len(x)))))[len(x):2*len(x)]
            if o['include_y_in_amplitude']:
                ey=abs(hilbert(np.pad(y,(len(y),len(y)))))[len(y):2*len(y)]
                envelope=np.hypot(envelope,ey)
            center=float(t[ids[np.argmax(envelope[ids])]])
            before=row['predicted_peak_s']-start;after=end-row['predicted_peak_s']
            start,end=center-before,center+after
            start=max(start,row.get('isolation_start_s',-np.inf))
            if o['limit_window_to_target_phase']:start,end=max(start,row['gate_start_s']),min(end,row['gate_end_s'])
        metrics,flags,window=measure(t,x,start,end,o,y=y);result.update(metrics,start_s=start,end_s=end)
        mask=(t>=start)&(t<=end)
        reference=max(metrics['target_direct_peak'],o['min_reference_amplitude'])
        result['orthogonal_peak_ratio']=float(max(abs(orth[mask]))/reference)
        result['y_peak_ratio']=float(max(abs(y[mask]))/reference) if y is not None else ''
        if result['orthogonal_peak_ratio']>.25:flags.append('large_orthogonal_component')
        if y is not None and result['y_peak_ratio']>.1:flags.append('large_y_component')
        if row['coincident_phases']:flags.append('coincident_phases')
        default_half=o['window_half_periods']/c['incident']['frequency_hz']
        before=row.get('nominal_before_s',o['window_before_s'] or default_half)
        after=row.get('nominal_after_s',o['window_after_s'] or default_half)
        for phase,value in json.loads(row['other_arrivals_s']).items():
            if start<=value+after and end>=value-before:flags.append('possible_overlap_'+phase)
        sampling=1/((t[1]-t[0])*c['incident']['frequency_hz']);result['samples_per_period']=sampling
        hard=sampling<o['min_samples_per_period'] or bool(row['coincident_phases']) or metrics['target_direct_peak']<=o['min_reference_amplitude']
        if metrics['target_direct_peak']<=o['min_reference_amplitude']:flags.append('target_projection_tiny')
        if sampling<o['min_samples_per_period']:flags.append('insufficient_sampling')
        result['usable']=int(not hard and (not flags or row['reviewed']))
        result['status']='ok' if not flags else ('reviewed:' if row['reviewed'] and not hard else 'review:')+';'.join(flags)
        return result,(t,x,window)
    except (ValueError,OSError) as exc:
        result['status']='error:'+str(exc);return result,None


def station_pairs(rows,c,o):
    pairs=[];incident=[r for r in rows if r['phase'].startswith('I_') and r['direct_peak']!='']
    for first in incident:
        for phase,info in c['phases'].items():
            if phase.startswith('I_') or not info['propagating']:continue
            candidates=[r for r in rows if r['phase']==phase and r['direct_peak']!='']
            if o['pairing_rule']!='all_pairs':
                distance=first['distance_to_interface_m']
                target=distance if o['pairing_rule']=='equal_distance' else distance*info['speed']/c['incident_speed']
                if candidates:
                    nearest=min(candidates,key=lambda r:abs(r['distance_to_interface_m']-target))
                    candidates=[nearest] if abs(nearest['distance_to_interface_m']-target)<=o['pairing_tolerance_m'] else []
            for second in candidates:
                theory=info['theory_abs'];direct=second['direct_peak']/first['direct_peak'];hp=second['hilbert_peak']/first['hilbert_peak']
                pairs.append(dict(phase=phase,incident_station=first['station'],outgoing_station=second['station'],
                    incident_distance_m=first['distance_to_interface_m'],outgoing_distance_m=second['distance_to_interface_m'],
                    incident_direct_peak=first['direct_peak'],outgoing_direct_peak=second['direct_peak'],
                    incident_hilbert_peak=first['hilbert_peak'],outgoing_hilbert_peak=second['hilbert_peak'],
                    direct_ratio_abs=direct,hilbert_ratio_abs=hp,theory_abs=theory,
                    direct_error_pct=100*abs(direct-theory)/theory if theory>1e-12 else '',
                    hilbert_error_pct=100*abs(hp-theory)/theory if theory>1e-12 else '',accepted=int(first['usable'] and second['usable'])))
    return pairs


def selectable(row,o):
    # 缺数据、采样不足、同一时刻震相无法区分，均不能因拟合好而纳入。
    return (row.get('direct_peak','')!='' and row.get('hilbert_peak','')!=''
            and row.get('samples_per_period',0)>=o['min_samples_per_period']
            and not row.get('coincident_phases')
            and row.get('target_direct_peak',row.get('direct_peak',0))>o['min_reference_amplitude']
            and (row['usable'] or not o['best_pair_accepted_only']))


def candidate_segments(rows,phase,c,o):
    line=sorted([r for r in rows if r['phase']==phase],key=lambda r:r['distance_to_interface_m'])
    if not line:return []
    if o['station_selection_mode']=='manual':
        windows=[line]
    else:
        length=o['auto_segment_length_m'];windows=[]
        step=c['receiver_settings'].get('ray_step_m',5.0)
        for first in line:
            start=first['distance_to_interface_m'];end=start+length
            block=[r for r in line if start-1e-6<=r['distance_to_interface_m']<=end+1e-6]
            # 端点必须真实覆盖要求长度；不能用短段或有几何空洞的台站拼出150m。
            if not block or abs(block[-1]['distance_to_interface_m']-end)>1e-5:continue
            if any(b['distance_to_interface_m']-a['distance_to_interface_m']>step*1.01+1e-6 for a,b in zip(block,block[1:])):continue
            windows.append(block)
    result=[]
    for block in windows:
        valid=[r for r in block if selectable(r,o)]
        if len(valid)<2 or len(valid)/len(block)<o['auto_segment_min_valid_fraction']:continue
        result.append(dict(start=block[0]['distance_to_interface_m'],end=block[-1]['distance_to_interface_m'],
                           rows=block,valid=valid,ids={r['station'] for r in valid}))
    return result


def select_segments(rows,pairs,c,o):
    best=[];segments=[];selected_pairs={};incphase='I_'+c['incident']['wave']
    incoming=candidate_segments(rows,incphase,c,o)
    method=o['best_pair_method'];lookup={(r['phase'],r['station']):r for r in rows}
    for phase,info in c['phases'].items():
        if phase.startswith('I_'):continue
        base=dict(phase=phase,best_status='unavailable',reason='',best_method=method,theory_abs=info['theory_abs'])
        if not info['propagating']:
            base['reason']='该模态为倏逝波，没有传播测线';best.append(base);continue
        outgoing=candidate_segments(rows,phase,c,o)
        if not incoming or not outgoing:
            details=[]
            for ph in (incphase,phase):
                line=[r for r in rows if r['phase']==ph];ds=[r['distance_to_interface_m'] for r in line]
                span=max(ds)-min(ds) if ds else 0
                details.append(f"{ph}: 跨度{span:g}m，可用于筛选{sum(selectable(r,o) for r in line)}/{len(line)}站")
            base['reason']='没有满足长度/连续性/有效比例的台站段；'+'；'.join(details)
            best.append(base);continue
        available=[p for p in pairs if p['phase']==phase and selectable(lookup[incphase,p['incident_station']],o)
                   and selectable(lookup[phase,p['outgoing_station']],o)]
        # 对全体组合误差评分，不允许只拿一对的最小误差充当整段拟合质量。
        pairmap={(p['incident_station'],p['outgoing_station']):p for p in available}
        winner=None
        for first in incoming:
            for second in outgoing:
                group=[pairmap[i,j] for i in sorted(first['ids']) for j in sorted(second['ids']) if (i,j) in pairmap]
                if not group:continue
                denominator=info['theory_abs'] if info['theory_abs']>1e-12 else 1.0
                errors=np.array([(p[method+'_ratio_abs']-info['theory_abs'])/denominator for p in group])
                score=float(np.sqrt(np.mean(errors**2)))
                rank=(score,-len(group),first['start'],second['start'])
                if winner is None or rank<winner[0]:winner=(rank,first,second,group)
        if winner is None:
            base['reason']='台站段内没有满足配对规则的可计算组合';best.append(base);continue
        rank,first,second,group=winner
        item=copy.deepcopy(min(group,key=lambda p:(abs(p[method+'_ratio_abs']-info['theory_abs']),p['incident_station'],p['outgoing_station'])))
        item.update(best_status='selected',reason='',best_method=method,segment_rmse=rank[0],
                    score_unit='相对误差' if info['theory_abs']>1e-12 else '绝对系数误差（理论接近零）',
                    candidate_count=len(group),incident_segment_start_m=first['start'],incident_segment_end_m=first['end'],
                    outgoing_segment_start_m=second['start'],outgoing_segment_end_m=second['end'],
                    incident_quality=lookup[incphase,item['incident_station']]['status'],outgoing_quality=lookup[phase,item['outgoing_station']]['status'])
        for m in ('direct','hilbert'):
            values=np.array([p[m+'_ratio_abs'] for p in group])
            for label,value in [('min',values.min()),('max',values.max()),('mean',values.mean()),('median',np.median(values)),('variance',values.var())]:item[m+'_'+label]=float(value)
            for label in ('mean','median'):
                difference=abs(item[m+'_'+label]-info['theory_abs'])
                item[m+'_'+label+'_absolute_error']=difference
                item[m+'_'+label+'_error_pct']=100*difference/info['theory_abs'] if info['theory_abs']>1e-12 else ''
        best.append(item);selected_pairs[phase]=group
        for role,ph,seg in [('入射',incphase,first),('出射',phase,second)]:
            segments.append(dict(coefficient=phase,role=role,phase=ph,**seg))
    return best,segments,selected_pairs


def quality_text(status):
    replacements={'ok':'通过检查','review:':'待核对：','reviewed:':'已人工复核：','error:':'无法拾取：',
      'window_edge_not_quiet':'窗口边缘仍有明显振幅','peak_at_window_edge':'峰值靠近窗边',
      'large_orthogonal_component':'正交分量较大','large_y_component':'Y分量较大',
      'possible_overlap_':'可能混入','insufficient_sampling':'采样不足','coincident_phases':'震相同时到达',
      'window_outside_record':'窗口超出记录','window_too_short':'窗内采样点不足','zero_or_tiny_signal':'振幅过小'}
    for a,b in replacements.items():status=status.replace(a,b)
    return status


def geometry_plot(c,rows,best,segments,o,out):
    fig,axes=plt.subplots(2,2,figsize=(15,12),layout='constrained')
    colors={'I_P':'#0072B2','I_SV':'#0072B2','R_P':'#D55E00','R_SV':'#009E73','T_P':'#CC79A7','T_SV':'#6b56a5'}
    geom=c['geometry'];mapping=geom.get('coordinate_mapping',{})
    box=mapping.get('research_xyz_m')
    source=np.array(geom.get('source_line_x_depth_m',[]))
    for ax,item in zip(axes.flat,best):
        phase=item['phase']
        ax.axhline(c['top'],color='#8c681e',lw=1,label='介质界面')
        ax.axhline(c['bottom'],color='#8c681e',lw=1)
        ax.axhspan(c['top'],c['bottom'],color='#cfad54',alpha=.7)
        if source.shape==(2,2):ax.plot(source[:,0],source[:,1],color='#8c348d',ls='--',lw=2,label='震源面投影')
        for ph in c['phases']:
            line=[r for r in rows if r['phase']==ph]
            line=sorted(line,key=lambda r:r['distance_to_interface_m'])
            if line:ax.plot([r['x_m'] for r in line],[r['burial_m'] for r in line],'.-',color=colors.get(ph,'gray'),alpha=.3,ms=3,lw=.8,label=ph+'测线')
        for seg in [s for s in segments if s['coefficient']==phase]:
            line=seg['rows'];col=colors[seg['phase']]
            ax.plot([r['x_m'] for r in line],[r['burial_m'] for r in line],'-o',color=col,lw=4,ms=3,
                    label=f"选中{seg['role']}段 {seg['start']:g}–{seg['end']:g}m")
            excluded=[r for r in line if r['station'] not in seg['ids']]
            if excluded:ax.scatter([r['x_m'] for r in excluded],[r['burial_m'] for r in excluded],marker='x',s=25,color='black',zorder=4,label=seg['role']+'段内未参与计算')
            name=item['incident_station' if seg['role']=='入射' else 'outgoing_station']
            row=next(r for r in line if r['station']==name)
            ax.scatter(row['x_m'],row['burial_m'],s=150,marker='*',color=col,edgecolor='black',zorder=5)
            ax.annotate(name+f" ({row['distance_to_interface_m']:g}m)",(row['x_m'],row['burial_m']),
                        xytext=(-110,-20) if seg['role']=='入射' else (15,15),textcoords='offset points',fontsize=8,
                        arrowprops=dict(arrowstyle='-',color=col),bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
        if item['best_status']=='selected':
            score=(f"相对RMSE={100*item['segment_rmse']:.1f}%" if item['theory_abs']>1e-12 else f"绝对RMSE={item['segment_rmse']:.3g}")
            title=f"{phase} | 段内{item['candidate_count']}组合 | {item['best_method']} {score}"
        else:
            title=phase+'：未选出完整台站段'
            ax.text(.02,.02,item['reason'],transform=ax.transAxes,fontsize=8,wrap=True,bbox=dict(facecolor='white',alpha=.9))
        ax.set(title=title,xlabel='全局 X / m',ylabel='相对外顶面的埋深 / m')
        if box:
            outer=mapping['outer_xyz_m'];d0=outer['z'][1]-box['z'][1];d1=outer['z'][1]-box['z'][0]
            from matplotlib.patches import Rectangle
            ax.add_patch(Rectangle((box['x'][0],d0),box['x'][1]-box['x'][0],d1-d0,fill=False,edgecolor='#69818b',ls=':',lw=1,label='研究区'))
            ax.set_xlim(box['x'][0]-80,box['x'][1]+130);ax.set_ylim(d1+30,min(d0-30,float(source[:,1].min())-30) if source.size else d0-30)
        else:ax.invert_yaxis()
        ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.15);ax.legend(fontsize=7,ncol=2,loc='lower left')
    for ax in list(axes.flat)[len(best):]:ax.axis('off')
    fig.suptitle(f"{o['case_name']}｜最佳连续台站段与最优组合（星号）\n按理论拟合筛选；距离沿测线从界面量起",fontsize=14)
    fig.savefig(out/'台站段与最优组合.png',dpi=o['plot_dpi']);plt.close(fig)


def waveform_plots(c,best,payloads,o,out):
    selected=[b for b in best if b['best_status']=='selected']
    if not selected:return
    units={'acceleration':'m/s²','velocity':'m/s','displacement':'m'}
    for zoom,filename in [(False,'最优组合全时程.png'),(True,'最优组合震相拾取.png')]:
        fig,axes=plt.subplots(len(selected),2,figsize=(17,3.6*len(selected)),squeeze=False,layout='constrained')
        for i,item in enumerate(selected):
            for j,(phase,name,role) in enumerate([('I_'+c['incident']['wave'],item['incident_station'],'入射分母'),(item['phase'],item['outgoing_station'],'出射分子')]):
                ax=axes[i,j];row,(t,x,(tw,xx,env))=payloads[phase,name]
                ax.plot(t,x,color='#286086',lw=.8,label='滤波后极化投影')
                ax.axvspan(row['start_s'],row['end_s'],color='#9ed2d1',alpha=.35,label='实际拾取窗')
                ax.plot(tw,env,color='#d47b20',lw=1.2,label='窗内Hilbert包络')
                if o['include_y_in_amplitude']:ax.plot(tw,xx,color='#009E73',lw=.8,alpha=.7,label='目标投影+Y合成幅值')
                di,hi=np.argmax(abs(xx)),np.argmax(env)
                ax.scatter(tw[di],xx[di],marker='s',s=28,color='#0072B2',zorder=5,label='合成直接峰' if o['include_y_in_amplitude'] else '直接峰（绝对值）')
                ax.scatter(tw[hi],env[hi],marker='o',s=28,color='#D55E00',zorder=5,label='包络峰')
                for ph,value in {phase:row['predicted_peak_s'],**json.loads(row['other_arrivals_s'])}.items():
                    ax.axvline(value,color='#7b627d' if ph==phase else '#888888',lw=.6,ls='--')
                    ax.text(value,.97,ph,transform=ax.get_xaxis_transform(),fontsize=7,rotation=90,va='top',clip_on=True)
                if zoom:
                    pad=.2*(row['end_s']-row['start_s']);ax.set_xlim(max(t[0],row['start_s']-pad),min(t[-1],row['end_s']+pad))
                    visible=x[(t>=ax.get_xlim()[0])&(t<=ax.get_xlim()[1])]
                    if len(visible):
                        limit=max(max(abs(visible)),max(env))*1.22;ax.set_ylim(-limit,limit)
                else:ax.set_xlim(t[0],t[-1])
                title=f"{item['phase']}｜{role} {name}，距界面{row['distance_to_interface_m']:g}m\n直接峰={row['direct_peak']:.4g}；包络峰={row['hilbert_peak']:.4g}"
                if j==1:title+=f"\n系数：直接{item['direct_ratio_abs']:.5g} / 包络{item['hilbert_ratio_abs']:.5g} / 理论{item['theory_abs']:.5g}"
                ax.set(title=title,xlabel='记录时间 / s',ylabel=units[o['trace_quantity']])
                ax.text(.01,.02,quality_text(row['status']),transform=ax.transAxes,fontsize=7,color='#934b22',wrap=True)
                ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.15)
                # 图例移至图外，避免遮住用于拾取的峰值。
                if i==0 and j==0:handles,labels=ax.get_legend_handles_labels()
        fig.legend(handles,labels,fontsize=8,ncol=5,loc='outside lower center')
        fig.suptitle(o['case_name']+'｜'+('震相拾取放大' if zoom else '最优组合完整时程')+
                     ('｜幅值包含Y分量' if o['include_y_in_amplitude'] else '｜纯目标投影')+'（不逐道归一化）',fontsize=14)
        fig.savefig(out/filename,dpi=o['plot_dpi']);plt.close(fig)


def workbook(c,rows,best,segments,payloads,o,out,prepared=False):
    # 保持可独立运行的Python 3.10程序，沿用原有openpyxl导出依赖。
    from openpyxl import Workbook
    from openpyxl.styles import Font,PatternFill,Alignment
    wb=Workbook();wb.remove(wb.active)
    def sheet(title,headers,records):
        ws=wb.create_sheet(title);ws.append(headers)
        for values in records:
            ws.append([float(v) if isinstance(v,np.floating) else (None if isinstance(v,float) and not math.isfinite(v) else v) for v in values])
        ws.freeze_panes='A2';ws.auto_filter.ref=ws.dimensions;ws.sheet_view.showGridLines=False
        for cell in ws[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='264B63')
        ws.row_dimensions[1].height=36
        for row in ws:
            for cell in row:
                cell.alignment=Alignment(vertical='top',wrap_text=True)
                if isinstance(cell.value,(float,int)):cell.number_format='0.######'
        for col in ws.columns:
            width=max(len(str(cell.value or '')) for cell in col)
            ws.column_dimensions[col[0].column_letter].width=min(48,max(16,width*1.2))
        # Excel不会可靠地对程序写入的换行内容自动调整行高。
        for row in ws.iter_rows(min_row=2):
            lines=max(math.ceil(sum(2 if ord(ch)>255 else 1 for ch in str(cell.value or ''))/
                                ws.column_dimensions[cell.column_letter].width) for cell in row)
            ws.row_dimensions[row[0].row].height=max(30,min(180,16*lines))
        return ws
    heads=['系数','结果','入射段起点m','入射段终点m','出射段起点m','出射段终点m','段内组合数','评分方法','段内误差RMSE','误差单位',
           '最佳入射站','最佳出射站','直接系数','Hilbert系数','理论模','直接误差%','Hilbert误差%','入射拾取状态','出射拾取状态']
    keys=['phase','reason','incident_segment_start_m','incident_segment_end_m','outgoing_segment_start_m','outgoing_segment_end_m','candidate_count','best_method','segment_rmse','score_unit',
          'incident_station','outgoing_station','direct_ratio_abs','hilbert_ratio_abs','theory_abs','direct_error_pct','hilbert_error_pct','incident_quality','outgoing_quality']
    summary=[]
    for b in best:
        r=[b.get(k,'') for k in keys];r[1]='已选中（按理论拟合）' if b['best_status']=='selected' else b['reason']
        r[-2:]=[quality_text(str(v)) for v in r[-2:]];summary.append(r)
    sheet('选段与最佳组合',heads,summary)
    stats=sheet('最佳段统计与理论比较',
          ['系数','方法','组合数','理论模','均值','均值相对理论误差%','中位数','中位数相对理论误差%',
           '最小','最大','总体方差','均值绝对误差','中位数绝对误差'],
          [[b['phase'],'直接法' if m=='direct' else 'Hilbert包络法',b['candidate_count'],b['theory_abs']]+
           [b[m+'_'+k] for k in ('mean','mean_error_pct','median','median_error_pct','min','max','variance','mean_absolute_error','median_absolute_error')]
           for b in best if b['best_status']=='selected' for m in ('direct','hilbert')])
    wb.move_sheet(stats,offset=-1)  # 打开工作簿首先看到整段统计与理论比较，而不是单个最优组合。
    sheet('选中段覆盖情况',['系数','作用','震相','起点m','终点m','长度m','段内台站数','参与台站数','有效比例'],
          [[s['coefficient'],s['role'],s['phase'],s['start'],s['end'],s['end']-s['start'],len(s['rows']),len(s['valid']),len(s['valid'])/len(s['rows'])] for s in segments])
    members=[]
    for seg in segments:
        for r in seg['rows']:
            members.append([seg['coefficient'],seg['role'],seg['phase'],r['station'],r['distance_to_interface_m'],
                            r['x_m'],r['burial_m'],r['station'] in seg['ids'],quality_text(r['status'])])
    sheet('选中段台站',['系数','作用','震相','台站','沿线距界面m','全局X米','埋深米','参与计算','拾取状态'],members)
    picks=[]
    for b in best:
        if b['best_status']!='selected':continue
        for ph,name in [('I_'+c['incident']['wave'],b['incident_station']),(b['phase'],b['outgoing_station'])]:
            r=payloads[ph,name][0]
            picks.append([b['phase'],ph,name]+[r[k] for k in ('start_s','end_s','predicted_peak_s','direct_peak','direct_peak_time_s','hilbert_peak','hilbert_peak_time_s')]+[quality_text(r['status']),r['amplitude_mode'],r['y_peak_ratio']])
    sheet('最佳组合拾取',['系数','震相','台站','窗起点s','窗终点s','理论到时s','直接幅值','直接峰时s','包络幅值','包络峰时s','拾取状态','幅值口径','Y峰值与目标投影峰值比'],picks)
    notes=[['算例',o['case_name']],['状态','仅检查/计划，未处理波形' if prepared else '已处理波形'],
           ['输入元数据',o['setup_directory']],['波形目录',o['data_directory']],['理论表',c['theory_table']],
           ['入射',f"{c['incident']['side']} / {c['incident']['wave']} / {c['incident']['angle_deg']:g}度 / {c['incident']['frequency_hz']:g}Hz"],
           ['筛选方式',o['station_selection_mode']],['每段长度m',o['auto_segment_length_m']],
           ['原请求段长m',o.get('requested_auto_segment_length_m',o['auto_segment_length_m'])],
           ['短段说明','prepare明确允许较短合格测线，已按其共同可用长度选段' if 'requested_auto_segment_length_m' in o else '未按prepare缩短段长'],
           ['最低有效比例',o['auto_segment_min_valid_fraction']],['只使用通过检查台站',o['best_pair_accepted_only']],
           ['段评分','对选中两段内全部有效组合的系数误差求RMSE；先选段，再选最佳一对；理论模近零时使用绝对误差'],
           ['段统计口径','均值、中位数、总体方差均针对最佳两段内全部有效入射×出射组合的系数；不是先平均幅值再相除'],
           ['均值/中位数误差','100×|段内均值或中位数－理论模|÷理论模；不是逐组合误差的平均值。理论模接近零时百分比留空，查看绝对误差'],
           ['系数定义','出射幅值÷入射幅值；非能量比'],
           ['包含Y分量',o['include_y_in_amplitude']],
           ['幅值口径','开启Y：直接取sqrt(s²+Y²)峰值；Hilbert分别对s与Y渐消补零求解析信号，再取sqrt(|Hs|²+|HY|²)峰值。关闭Y：恢复原目标投影幅值。'],
           ['理论比较口径','理论仍为纯P/SV振幅系数；含Y数值幅值并非纯P/SV。Y很小时近似一致，Y明显时需结合Y占比判断。'],
           ['拟合含义','按理论误差筛选属于事后最优拟合，不是独立的理论验证；待核对组合不等于通过物理质量检查'],
           ['理论厚度','当前理论外侧为半空间，只核对薄层厚度；数值到时使用prepare中的实际界面埋深'],
           ['滤波',f"Butterworth {o['filter_order']}阶带通，sosfiltfilt双向零相位；频带{o['bandpass_hz']}Hz（None为关闭）"],
           ['时间窗',f"中心={o['window_center_mode']}；读取生成器隔离窗={o['use_generator_reflection_windows']}；新共同交点布局的入射/反射均按布站时波包半宽拾取"],
           ['提示','相邻震相窗口截断、滤波展宽、薄层多次波均可能影响峰值；请查看震相拾取图']]
    for phase in c['phases']:
        subset=[r for r in rows if r['phase']==phase]
        notes.append([phase+'台站统计',f"总数{len(subset)}；可测幅值{sum(r.get('direct_peak','')!='' for r in subset)}；用于筛选{sum(selectable(r,o) for r in subset)}"])
        for status,count in Counter(r.get('status','未处理') for r in subset).items():notes.append([phase+f' {count}站',quality_text(status)])
    sheet('设置与说明',['项目','内容'],notes)
    if prepared:sheet('拾取计划',['震相','台站','沿线距界面m','理论到时s','窗起点s','窗终点s'],[[r[k] for k in ('phase','station','distance_to_interface_m','predicted_peak_s','start_s','end_s')] for r in rows])
    wb.save(out/'反透射系数结果.xlsx');wb.close()


def run(o):
    check_inputs(o)
    out=resolve(o['output_directory'])
    inputs=[resolve(o[k]) for k in ('setup_directory','data_directory','theory_directory')]
    if any(out==p or out in p.parents or p in out.parents for p in inputs):raise ValueError('输出目录必须与输入目录分开')
    out.mkdir(parents=True,exist_ok=True)
    c=load_context(o);plan=make_plan(c,o);data=resolve(o['data_directory'])
    rows=[];payloads={};best=[];segments=[]
    if o['run_mode']=='analyze':
        if not data.is_dir():raise ValueError('本算例缺少sema波形文件夹：'+str(data))
        index=trace_files(data,o)
        print(f"算例 {o['case_name']}：读取{len(plan)}条测线台站记录（不处理网格台站）",flush=True)
        for n,row in enumerate(plan,1):
            result,payload=process_row(row,c,index,o);rows.append(result)
            if payload is not None:payloads[row['phase'],row['station']]=(result,payload)
            if n%50==0:print(f'已拾取 {n}/{len(plan)}',flush=True)
        pairing=copy.deepcopy(o)
        if o['station_selection_mode']=='auto_segment':pairing['pairing_rule']='all_pairs'
        pairs=station_pairs(rows,c,pairing)
        best,segments,_=select_segments(rows,pairs,c,o)
    else:rows=plan
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False})
    # 先在临时目录完成全部产物，再替换最终文件；临时过程数据不会留在输出中。
    with tempfile.TemporaryDirectory(prefix='rt_export_') as staging:
        stage=Path(staging)
        workbook(c,rows,best,segments,payloads,o,stage,prepared=o['run_mode']=='prepare')
        if o['run_mode']=='analyze':
            geometry_plot(c,rows,best,segments,o,stage);waveform_plots(c,best,payloads,o,stage)
        import shutil
        generated={p.name for p in stage.iterdir()}
        for p in stage.iterdir():shutil.copy2(p,out/p.name)
        # 只清理本程序明确拥有的旧输出，不删除用户放入的其他文件。
        owned={'window_plan.csv','对比计划说明.txt','run_report.json','amplitudes.csv','station_pairs.csv',
               'best_station_pairs.csv','反透射系数结果.png','最优台站全时程.png',
               '台站段与最优组合.png','最优组合全时程.png','最优组合震相拾取.png'}
        for name in owned-generated:
            target=out/name
            if target.is_file():target.unlink()
    for b in best:print(b['phase'],b.get('incident_segment_start_m',''),b.get('outgoing_segment_start_m',''),b['best_status'],b['reason'])
    print('结果：',out)
    return dict(status='completed',best=best,output=str(out),measured_station_phases=len(payloads))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',help='只填写result下的算例名称')
    parser.add_argument('--prepare',action='store_true',help='只检查输入，不计算幅值或选段')
    args=parser.parse_args();options=settings()
    if args.case:options['case_name']=args.case;case_paths(options)
    if args.prepare:options['run_mode']='prepare'
    run(options)
