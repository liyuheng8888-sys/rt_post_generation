#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python 3.10 震源与台站生成器。
直接修改下方“用户编辑区”，然后运行 python generate_sources_stations.py。
模型：1000 m正方体；用户输入、图件和坐标表均采用埋深，向下为正。
顶面埋深0 m，海水—海床界面埋深500 m，底部埋深1000 m。
生成有限震源面的平面波近似几何及输入文件，不执行SPECFEM求解。
"""
from __future__ import annotations

# ============================================================================
#                         用户编辑区（所有输入参数）
# ============================================================================
# 日常使用只修改本区，不需要修改下方计算代码，也不需要编辑 config.json。
# 单位：长度 m、波速 m/s、密度 kg/m³、频率 Hz、时间 s、角度 °。
# 注意：字符串保留引号；每项末尾保留逗号；Python开关写 True / False。

USER_CONFIG = {
    # ---------- 1. 入射角、波型与震源中心间距 ----------
    "incident_side": "lower",        # "lower"下层入射；"upper"上层入射
    # 两侧端点采用相同的[X,埋深]格式；只启用incident_side所选侧，另一侧保留备用。
    "upper_source_anchor_x_depth_m": [400.0, 50.0],   # 上层：0≤埋深<500-clearance
    "lower_source_anchor_x_depth_m": [400.0, 950.0], # 下层：500+clearance<埋深≤1000
    # 镜像示例：上层[400,50] ↔ 下层[400,950]；X范围均为0<X≤1000。
    "top_is_free_surface": True,      # 与实际网格边界一致；用于检查，不会修改网格或Par_file
    "bottom_is_free_surface": False,  # 底面同样可检查；通常为吸收边界，故默认False
    # 声学源不能落在所选侧的自由表面上（顶面埋深0、底面埋深1000）。
    # 以下角度、波型、间距、源型和接收器参数由上下层共用，切换侧别无需换一套键名。
    "angle_deg": 40,                # 与-x夹角大小；下层源线朝浅处延伸，上层朝深处延伸
    "incident_wave": "SV",            # 固体可选"P"或"SV"；流体只允许"P"
    "source_spacing_m": 10.0,         # 沿倾斜中心线的实际距离；可改为5.0
    "source_y_spacing_m": 10.0,       # 沿y方向的间距；整个源面加密时也改为5.0
    "source_interface_clearance_m": 50.0,  # 源中心与界面的最小竖向间隙，0<值<500

    # ---------- 2. 接收器 ----------
    # 两类台站均固定在同一个y=500 m剖面；埋深500 m处为介质分界线。
    "generate_auto_stations": True,     # 自动测线台站输出总开关：True输出；False不输出auto组
    # 选线时用波束外有效段筛选；选中后整条反射线都布站，包括波束覆盖段。
    "auto_candidate_step_m": 1.0,       # 界面交点扫描步长；越小选线越精细
    "auto_reflected_mode": "R_SV",       # 自动选线依据：R_P或R_SV（须为传播波）
    "auto_beam_clearance_m": 10.0,      # 离初始入射波束边缘的法向净距
    "auto_min_stations": 10,            # 合格反射线至少包含多少个5m间距台站
    "ray_station_naming": "distance", # distance名称含距界面距离；legacy恢复旧序号
    "generate_grid_stations":True,     # 是否生成y=500剖面密铺网格台站
    "generate_first_stations": False,    # 是否生成整个源面首次触界的first射线台站（入射、反射、透射）
    "generate_central_stations":False,  # 是否生成源面中间对应的central射线台站（入射、反射、透射）
    "generate_front_stations": False,    # 是否生成理论波阵面front台站；与first不是同一类
    "generate_interface_stations": False, # 是否生成埋深500m界面交线台站
    # 六个开关独立控制；只要central时，仅central设True，其余五项设False。
    # 关闭台站不会关闭相应理论几何计算；所有开关False时输出空台站文件，仍生成震源。
    "station_spacing_m": 10.0,        # y=500剖面的X/埋深底网格间距；可改为5.0
    "grid_station_prefix": "G",       # 第一类：平面密铺台站，例如G000001
    "line_station_prefix": "L",       # 第二类：理论线台站，例如L_first_I_P_000001
    # 先网格、再理论线、最后界面；重合点以后生成的名字覆盖，分组表同步更新。
    # 理论射线、波阵面与界面上的台站固定按5 m加密，无需另设参数。
    # STATIONS列顺序为：名称 DB y x 0 埋深_m；第5列为格式占位，不需输入高程。

    # ---------- 3. 源型、源时间函数与力幅值 ----------
    "source_kind": "auto",            # "auto"：固体用paper_couples，流体用acoustic_monopole
                                      # "paper_couples"：固体文献两对力偶，每中心4条FORCE
                                      # "plane_force"：每中心1个定向力，严格位于中心面
                                      # "acoustic_monopole"：流体每中心1个标量源，忽略力方向
    # 声学单极子要求setup/constants.h中DIPOLE_SOURCE_IN_FLUID=.false.（官方默认）。
    "couple_separation_m": 1.0,        # 正负力点间距，仅paper_couples有效
    "frequency_hz": 100.0,              # Ricker主频；source time function固定为1
    "enable_uniform_time_shift": False,  # True启用统一time shift；False时所有震源均写0
    "uniform_time_shift_s": 0.0,      # 开启后所有FORCE共用此值，单位秒，要求有限且≥0
    # 例如：开关True、数值0.3，则每条震源均写time shift: 0.3（不会给首源另设0）。
    # 官方SPECFEM要求至少一条time shift=0；全非零可能不兼容，生成时会提示。
    "force_factor": 1e15,             # 每条FORCE幅值；固体力单位N，声学源幅值需独立标定
    "force_direction": "auto",        # 自动按P/SV设置；仅plane_force可改为[1.0,-2.0,-1.0]
                                      # 向量顺序[E,N,Z_UP]；自定义方向可能激发混合波型

    # ---------- 4. 上层介质（埋深0～500 m） ----------
    "upper": {
        "name": "上层流体",            # 介质名称
        "type": "fluid",              # "fluid"流体 或 "solid"固体
        "rho": 1000.0,                 # 密度 kg/m³，必须>0
        "vp": 1500.0,                  # P波速度 m/s，必须>0
        "vs": 0.0,                     # 流体必须为0；固体须>0且vp²>4vs²/3
    },

    # ---------- 5. 下层介质（埋深500～1000 m） ----------
    "lower": {
        "name": "下层固体",            # 介质名称
        "type": "solid",              # "solid"固体 或 "fluid"流体
        "rho": 2500.0,                 # 密度 kg/m³，必须>0
        "vp": 4500.0,                  # P波速度 m/s
        "vs": 2100.0,                  # S波速度 m/s；须>0且vp²>4vs²/3
    },

    # ---------- 6. 理论波阵面台站与报告（不输出波阵面图） ----------
    "snapshot_after_first_s": 0.05,    # 首次接触界面后的传播时间，必须≥0

    # ---------- 7. SPECFEM文件格式（通常无需改动） ----------
    "force_format": "11_lines",       # "11_lines"：本次提供的格式，含源时间函数行
                                      # "legacy_10_lines"：旧格式，不含该行
    "force_depth_unit": "km",         # 源文件埋深默认km；台站埋深始终为m
                                      # 只有明确修改过求解器读入单位时才改为"m"
}

OUTPUT_DIRECTORY = "output"           # 相对路径以脚本所在目录为基准，也可填绝对路径
GENERATE_PLOTS = True                  # True：台站图及启用auto时的自动测线图；False：仅数据
PLOT_DPI = 220                         # 图片清晰度；同时输出可无限放大的PDF
STATION_ZOOM_WIDTH_M = 80.0             # 台站局部放大窗口宽度，单位m

# 固定模型约定：边长1000 m、界面埋深500 m、接收剖面y=500 m，
# 上下层震源端点均见上方参数；理论线台站间距5 m。
# 运行时不带--config：使用本区参数。
# 显式传入--config 某文件.json：使用该文件完整替代USER_CONFIG。
# --output覆盖OUTPUT_DIRECTORY；--no-plot强制关闭绘图。
# ========================== 用户编辑区结束 ==================================
#                    以下为程序实现，正常使用无需修改
# ============================================================================

import argparse
import copy
import csv
import json
import math
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 固定几何约束（与本任务模型一致，不属于日常参数入口）。
L = 1000.0                     # 正方体边长；内部几何以底部为零、向上为正
H = 500.0                      # 内部界面位置；对用户输入和输出统一转换为埋深
Y = 500.0                      # 接收器剖面y
LINE_STEP = 5.0                # 理论线上固定的加密间距


def station_switches(c):
    """旧JSON未填写开关时保留原先全部生成的行为。"""
    return {family:c.get('generate_'+family+'_stations',True)
            for family in ('grid','first','central','front','interface')} | {'auto':c.get('generate_auto_stations',False)}


def station_group_enabled(c,group):
    return station_switches(c)[group.split('_',1)[0]]


def samples(length, spacing):
    """Exact spacing, plus endpoint when remainder exists (last gap <= spacing)."""
    v = np.arange(0.0, length + 1e-8, spacing)
    if not len(v) or length - v[-1] > 1e-7:
        v = np.append(v, length)
    v[-1] = length
    return v


def segment_points(a, b, spacing):
    a, b = np.asarray(a, float), np.asarray(b, float)
    length = np.linalg.norm(b-a)
    if length < 1e-10:
        return a[None, :]
    points = a + samples(length, spacing)[:, None] * (b-a)/length
    points[-1] = b  # 精确保留端点，避免浮点舍入使边界点落到模型外。
    return points


def clip_line(origin, direction, bounds, lo=-np.inf, hi=np.inf):
    for o, d, (a, b) in zip(origin, direction, bounds):
        if abs(d) < 1e-12:
            if o < a-1e-8 or o > b+1e-8:
                return None
        else:
            q = sorted(((a-o)/d, (b-o)/d))
            lo, hi = max(lo, q[0]), min(hi, q[1])
    if hi < lo-1e-8:
        return None
    return np.asarray(origin)+lo*np.asarray(direction), np.asarray(origin)+hi*np.asarray(direction)


def speed(m, mode):
    return m['vp'] if mode == 'P' else m['vs']


def modes(m):
    return ['P', 'SV'] if m['type'] == 'solid' else ['P']


def incident_side(c):
    return c.get('incident_side', 'lower')


def source_kind(c):
    kind = c['source_kind']
    if kind == 'auto':
        return 'paper_couples' if c[incident_side(c)]['type'] == 'solid' else 'acoustic_monopole'
    return kind


def uniform_time_shift(c):
    """File value shared by every source, including all four forces of a couple."""
    enabled=c.get('enable_uniform_time_shift',False)
    value=c.get('uniform_time_shift_s',0.0)
    if not isinstance(enabled,bool):
        raise ValueError('enable_uniform_time_shift must be True or False.')
    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0:
        raise ValueError('uniform_time_shift_s must be finite and nonnegative (seconds).')
    return float(value) if enabled else 0.0


def source_anchor_internal(c):
    """Read the active side only; preserve historical anchors for old JSON files."""
    side = incident_side(c)
    key = side+'_source_anchor_x_depth_m'
    if side+'_source_anchor_xz' in c:
        raise ValueError('Old coordinate key is not accepted. Use '+key+'=[X,burial_m].')
    default = [750.0,50.0] if side == 'upper' else [750.0,1000.0]
    anchor=np.asarray(c.get(key,default),float)
    if anchor.shape!=(2,) or not np.isfinite(anchor).all():
        raise ValueError(key+' must be finite [X,burial_m].')
    return np.array([anchor[0],L-anchor[1]])


def x_depth(point):
    return [float(point[0]),float(L-point[1])]


def validate(c):
    uniform_time_shift(c)
    if c.get('ray_station_naming','distance') not in ('distance','legacy'):raise ValueError('ray_station_naming must be distance or legacy')
    for key,default in [('auto_candidate_step_m',1.),('auto_beam_clearance_m',10.)]:
        value=c.get(key,default)
        if not math.isfinite(value) or value<0 or (key=='auto_candidate_step_m' and value==0):raise ValueError(key+' invalid')
    if type(c.get('auto_min_stations',10)) is not int or c.get('auto_min_stations',10)<2:raise ValueError('auto_min_stations must be integer >=2')
    for family,enabled in station_switches(c).items():
        if not isinstance(enabled,bool):
            raise ValueError('generate_'+family+'_stations must be True or False.')
    side = incident_side(c)
    if side not in ('upper', 'lower'):
        raise ValueError('incident_side must be upper or lower.')
    anchor = source_anchor_internal(c)
    for key, default in [('top_is_free_surface',True),('bottom_is_free_surface',False)]:
        if not isinstance(c.get(key,default),bool):
            raise ValueError(key+' must be True or False.')
    # Distance from the outer boundary increases toward the interface on either side.
    inward_depth = L-anchor[1] if side == 'upper' else anchor[1]
    if not (0 < anchor[0] <= L and 0 <= inward_depth < H-c['source_interface_clearance_m']):
        bounds = '0<=burial_m<500-clearance' if side == 'upper' else '500+clearance<burial_m<=1000'
        raise ValueError(side+'_source_anchor_x_depth_m requires 0<X<=1000 and '+bounds+'.')
    free_surface = c.get('top_is_free_surface',True) if side == 'upper' else c.get('bottom_is_free_surface',False)
    if c[side]['type']=='fluid' and free_surface and inward_depth<=1e-6:
        raise ValueError('Acoustic source is on the '+side+' free surface: move its anchor inside the medium.')
    prefixes = [c.get('grid_station_prefix','G'), c.get('line_station_prefix','L')]
    if any(not isinstance(p,str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9]{0,3}',p) for p in prefixes):
        raise ValueError('Station prefixes must have 1-4 ASCII letters/digits, starting with a letter.')
    if prefixes[0].startswith(prefixes[1]) or prefixes[1].startswith(prefixes[0]):
        raise ValueError('Grid and line prefixes must differ and must not contain one another as a prefix.')
    if not 0 <= c['angle_deg'] < 90:
        raise ValueError('angle_deg must be in [0,90).')
    for key in ('source_spacing_m', 'source_y_spacing_m', 'station_spacing_m',
                'frequency_hz', 'force_factor', 'couple_separation_m'):
        if not math.isfinite(c[key]) or c[key] <= 0:
            raise ValueError(key+' must be finite and positive.')
    if not 0 < c['source_interface_clearance_m'] < H:
        raise ValueError('source_interface_clearance_m must be in (0,500).')
    if not math.isfinite(c['snapshot_after_first_s']) or c['snapshot_after_first_s'] < 0:
        raise ValueError('snapshot_after_first_s must be finite and nonnegative.')
    for key in ('upper', 'lower'):
        m = c[key]
        if m['type'] not in ('solid', 'fluid'):
            raise ValueError(key+': type must be solid or fluid.')
        if not all(math.isfinite(m[k]) for k in ('rho', 'vp', 'vs')):
            raise ValueError('Medium parameters must be finite.')
        if m['rho'] <= 0 or m['vp'] <= 0:
            raise ValueError('rho and vp must be positive.')
        if m['type'] == 'fluid' and m['vs'] != 0:
            raise ValueError('Fluid vs must be 0.')
        if m['type'] == 'solid' and not (m['vs'] > 0 and m['vp']**2 > 4*m['vs']**2/3):
            raise ValueError('Solid requires vs>0 and positive bulk modulus vp^2>4vs^2/3.')
    if c['incident_wave'] not in modes(c[side]):
        raise ValueError('Incident mode is unavailable in the selected medium (fluid supports P only).')
    if c['source_kind'] not in ('auto', 'paper_couples', 'plane_force', 'acoustic_monopole'):
        raise ValueError('Unknown source_kind.')
    fluid = c[side]['type'] == 'fluid'
    if fluid != (source_kind(c) == 'acoustic_monopole'):
        raise ValueError('Use auto or acoustic_monopole for fluid; paper_couples/plane_force for solid.')
    if c['force_format'] not in ('11_lines', 'legacy_10_lines'):
        raise ValueError('Unknown force_format.')
    if c['force_depth_unit'] not in ('km', 'm'):
        raise ValueError('force_depth_unit must be km or m.')
    if c['force_direction'] != 'auto':
        d = np.asarray(c['force_direction'], float)
        if d.shape != (3,) or not np.isfinite(d).all() or np.linalg.norm(d) == 0:
            raise ValueError('force_direction must be auto or a nonzero [E,N,Z_UP] vector.')
        if source_kind(c) != 'plane_force':
            raise ValueError('Custom direction applies only to plane_force.')


def state(m, mode, p, sign):
    """[ux,uz,txz/(i*w),tzz/(i*w)], z DOWN as in the existing theory script."""
    v = speed(m, mode)
    q = sign*np.sqrt(complex(1/v**2-p*p))
    rho, vp, vs = m['rho'], m['vp'], m['vs']
    if m['type'] == 'fluid':
        return np.array([vp*p, vp*q, 0j, rho*vp])
    mu, lam = rho*vs**2, rho*(vp**2-2*vs**2)
    if mode == 'P':
        return np.array([vp*p, vp*q, 2*mu*vp*p*q, lam/vp+2*mu*vp*q*q])
    return np.array([-vs*q, vs*p, mu*vs*(p*p-q*q), 2*mu*vs*p*q])


def coefficients(c, p):
    inc = c[incident_side(c)]
    trans = c['lower' if incident_side(c) == 'upper' else 'upper']
    sign = 1 if incident_side(c) == 'lower' else -1  # physical z UP incidence
    both_solid = inc['type'] == trans['type'] == 'solid'
    both_fluid = inc['type'] == trans['type'] == 'fluid'
    indices = [0,1,2,3] if both_solid else ([1,3] if both_fluid else [1,3,2])
    labels = [('R_'+m, inc, m, sign) for m in modes(inc)]
    labels += [('T_'+m, trans, m, -sign) for m in modes(trans)]
    columns = [(-1 if k.startswith('R') else 1)*state(m, w, p, s)[indices]
               for k, m, w, s in labels]
    a = np.column_stack(columns)
    b = state(inc, c['incident_wave'], p, -sign)[indices]
    scale = np.max(abs(a), axis=1)
    solution = np.linalg.solve(a/scale[:, None], b/scale)
    return labels, solution, float(np.linalg.norm(a@solution-b)/max(np.linalg.norm(b), 1))


def geometry(c):
    alpha = math.radians(c['angle_deg'])
    side = incident_side(c)
    sign = 1 if side == 'lower' else -1
    anchor = source_anchor_internal(c)
    inc_bounds = (0,H) if side == 'lower' else (H,L)
    trans_bounds = (H,L) if side == 'lower' else (0,H)
    tangent = np.array([-math.cos(alpha), sign*math.sin(alpha)])
    normal = np.array([math.sin(alpha), sign*math.cos(alpha)])
    end_length = anchor[0]/math.cos(alpha)
    if math.sin(alpha) > 1e-12:
        end_length = min(end_length, (abs(H-anchor[1])-c['source_interface_clearance_m'])/math.sin(alpha))
    end = np.clip(anchor+end_length*tangent, [0.0, 0.0], [L, L])
    centers = segment_points(anchor, end, c['source_spacing_m'])
    v = speed(c[side], c['incident_wave'])
    p = normal[0]/v
    # Intersect the finite bundle with the interface AND the lateral model bounds.
    hit0 = anchor[0]+(H-anchor[1])*normal[0]/normal[1]
    hit1 = end[0]+(H-end[1])*normal[0]/normal[1]
    xmin, xmax = max(0.0, hit1), min(L, hit0)
    if xmin > xmax+1e-8:
        raise ValueError('The finite incident ray bundle exits x=1000 before reaching the interface. '
                         'Reduce angle or source_interface_clearance_m.')
    first = np.array([xmin, H])
    # central优先指源面几何中点的射线，而不是裁剪后的界面覆盖区中点。
    source_middle = (anchor+end)/2
    middle_hit = source_middle+(H-source_middle[1])*normal/normal[1]
    central_is_source_midpoint = bool(0 <= middle_hit[0] <= L)
    central = middle_hit if central_is_source_midpoint else np.array([(xmin+xmax)/2, H])
    rays = []
    def add(name, ends):
        if ends is not None:
            rays.append(dict(name=name, a=np.asarray(ends[0]), b=np.asarray(ends[1])))
    for tag, hit in [('first', first), ('central', central)]:
        dist = np.dot(normal, hit-anchor)
        src = hit-dist*normal
        add(tag+'_I_'+c['incident_wave'], (src, hit))
    labels, values, residual = coefficients(c, p)
    branches = []
    for (name, medium, mode, _), value in zip(labels, values):
        vm = speed(medium, mode)
        sine = p*vm
        status = 'evanescent' if sine > 1+1e-10 else ('critical' if abs(sine-1) <= 1e-10 else 'propagating')
        branch = dict(name=name, speed=vm, sin_angle=sine, status=status,
                      coefficient_real=float(value.real), coefficient_imag=float(value.imag),
                      coefficient_abs=float(abs(value)))
        if status != 'evanescent':
            cosine = math.sqrt(max(0.0, 1-sine*sine))
            direction = np.array([sine, sign*(-1 if name.startswith('R') else 1)*cosine])
            branch['angle_deg'] = math.degrees(math.asin(min(1, sine)))
            branch['direction_xz'] = direction.tolist()
            bounds = [(0, L), inc_bounds if name.startswith('R') else trans_bounds]
            for tag, hit in [('first', first), ('central', central)]:
                add(tag+'_'+name, clip_line(hit, direction, bounds, lo=0))
        else:
            branch['decay_length_at_f0_m'] = 1/(2*math.pi*c['frequency_hz']*math.sqrt(p*p-1/vm**2))
        branches.append(branch)
    first_time = float(np.dot(normal, first-anchor)/v)
    time = first_time+c['snapshot_after_first_s']
    # Incident wavefront, clipped to the source aperture, lower medium and cube.
    add('front_I_'+c['incident_wave'], clip_line(anchor+v*time*normal, tangent,
        [(0, L), inc_bounds], lo=0, hi=end_length))
    # Outgoing front: r(h,t)=(h,H)+v_j*(t-t_hit(h))*n_j, affine in h.
    for b in branches:
        if b['status'] != 'propagating':
            continue
        d = np.array(b['direction_xz'])
        hmax = xmax if p == 0 else min(xmax, xmin+(time-first_time)/p)
        a = first+b['speed']*(time-first_time)*d
        derivative = np.array([1.0, 0.0])-b['speed']*p*d
        add('front_'+b['name'], clip_line(a, derivative,
            [(0,L), inc_bounds if b['name'].startswith('R') else trans_bounds], lo=0, hi=max(0,hmax-xmin)))
    result = dict(centers=centers, tangent=tangent, normal=normal, end=end, anchor=anchor,
                incident_side=side, incident_bounds=inc_bounds, central_is_source_midpoint=central_is_source_midpoint,
                rays=rays, branches=branches, p=p, first=first, central=central,
                first_time=first_time, snapshot_time=time, residual=residual,
                interface_footprint=[xmin, xmax], v=v)
    if c.get("generate_auto_stations",False):add_auto_rays(c,result)
    return result


def add_auto_rays(c,g):
    """Scan interface hit positions; choose the median eligible sampled ray.

    Exclusion is outside the whole incident beam's transverse strip, not just
    outside the instantaneous wavefront. It is a geometrical-optics criterion.
    """
    mode=c.get('auto_reflected_mode','R_P')
    branch=next((b for b in g['branches'] if b['name']==mode),None)
    if not mode.startswith('R_') or branch is None or branch['status']!='propagating':
        raise ValueError('Auto target must be a propagating reflected mode: '+mode)
    aperture=float(np.linalg.norm(g['end']-g['anchor']))
    margin=c.get('auto_beam_clearance_m',10.)
    minimum=c.get('auto_min_stations',10)
    candidates=[];valid=[]
    xmin,xmax=g['interface_footprint']
    def reflection_points(hit,b):
        d=np.array(b['direction_xz'])
        ends=clip_line(hit,d,[(0,L),g['incident_bounds']],lo=0)
        if ends is None:return np.empty((0,2))
        length=float(np.linalg.norm(ends[1]-hit))
        # No irregular endpoint: distance from interface is an integer multiple of 5m.
        pts=hit+np.arange(LINE_STEP,length+1e-8,LINE_STEP)[:,None]*d
        transverse=(pts-g['anchor'])@g['tangent']
        return pts[(transverse < -margin-1e-8)|(transverse > aperture+margin+1e-8)]
    for x in xmin+samples(xmax-xmin,c.get('auto_candidate_step_m',1.)):
        hit=np.array([x,H]);pts=reflection_points(hit,branch)
        row=dict(hit_x_m=float(x),station_count=len(pts),eligible=len(pts)>=minimum)
        candidates.append(row)
        if row['eligible']:valid.append((hit,pts))
    if not valid:
        raise ValueError('No auto ray outside incident beam meets auto_min_stations; change geometry/clearance/minimum. Best count='+str(max(r['station_count'] for r in candidates)))
    hit,pts=valid[(len(valid)-1)//2]
    g['auto']=hit
    dist=float(np.dot(g['normal'],hit-g['anchor']))
    g['rays'].append(dict(name='auto_I_'+c['incident_wave'],a=hit-dist*g['normal'],b=hit))
    for b in g['branches']:
        if b['status']!='propagating':continue
        if b['name'].startswith('R_'):
            ends=clip_line(hit,np.array(b['direction_xz']),[(0,L),g['incident_bounds']],lo=0)
            if ends is None:continue
            a,z=ends
        else:
            bounds=[(0,L),(0,H) if incident_side(c)=='upper' else (H,L)]
            ends=clip_line(hit,np.array(b['direction_xz']),bounds,lo=0)
            if ends is None:continue
            a,z=ends
        g['rays'].append(dict(name='auto_'+b['name'],a=a,b=z))
    g['auto_selection']=dict(target_mode=mode,scan_step_m=c.get('auto_candidate_step_m',1.),
        eligible_candidate_count=len(valid),candidate_count=len(candidates),
        selected_hit_x_m=float(hit[0]),selected_outside_beam_station_count=len(pts),
        selected_reflected_station_count=sum(len(segment_points(r['a'],r['b'],LINE_STEP)) for r in g['rays'] if r['name']=='auto_'+mode),
        output_policy="full reflected ray from interface to model boundary, including incident-beam overlap",
        beam_clearance_m=margin,selection='lower median of eligible sampled interface hits',
        limitation='Only selection uses outside-beam samples. Exported full rays include incident-beam overlap; diffraction and boundary waves are not excluded.')
    g['auto_candidates']=candidates


def force_points(c, g):
    """Yield (point3, unit_force3, center_id). Couples approximate Fig.4 in xz.

    P: M = n*n^T + t*t^T (2D isotropic expansion).
    SV: M = t*n^T + n*t^T (rotated double couple).
    Symmetric pairs are used where possible; boundary pairs are one-sided.
    """
    n, t = g['normal'], g['tangent']
    kind = source_kind(c)
    pairs = [(n,n), (t,t)] if c['incident_wave'] == 'P' else [(n,t), (t,n)]
    cid = 0
    for center in g['centers']:
        for y in samples(L, c['source_y_spacing_m']):
            cid += 1
            if kind == 'acoustic_monopole':
                # SPECFEM acoustic monopole path ignores direction (DIPOLE_SOURCE_IN_FLUID=false).
                yield np.array([center[0],y,center[1]]), np.array([1.0,0.0,0.0]), cid
                continue
            if kind == 'plane_force':
                vec2 = n if c['incident_wave'] == 'P' else t
                force = np.array([vec2[0], 0.0, vec2[1]]) if c['force_direction'] == 'auto' else np.array(c['force_direction'], float)
                yield np.array([center[0], y, center[1]]), force/np.linalg.norm(force), cid
                continue
            for separation, force2 in pairs:
                delta = c['couple_separation_m']*separation
                def inside(point):
                    z_ok = 0 <= point[1] < H if incident_side(c)=='lower' else H < point[1] <= L
                    return 0 <= point[0] <= L and z_ok
                a,b = center-delta/2, center+delta/2
                if not (inside(a) and inside(b)):
                    a,b = center, center+delta
                    if not inside(b):
                        a,b = center-delta,center
                if not (inside(a) and inside(b)):
                    raise ValueError('Couple does not fit in source medium; reduce separation.')
                for pos, sign in [(a,-1),(b,1)]:
                    yield np.array([pos[0],y,pos[1]]), sign*np.array([force2[0],0,force2[1]]),cid


def write_csv(path, header, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer=csv.writer(f); writer.writerow(header); writer.writerows(rows)


def build_stations(c,g):
    """Grid first, then ordered lines; last insertion owns a coincident position.

    Store membership by stable index, resolving names only AFTER all replacements.
    Coordinates are keyed to 1 micrometre, consistent with existing deduplication.
    """
    stations=[]; lookup={}; memberships=[]; owners=[]; replacements=[]; named_positions={}
    gp=c.get('grid_station_prefix','G'); lp=c.get('line_station_prefix','L')
    def add(x,z,group,number):
        name=f'{gp}{number:06d}' if group=='grid' else f'{lp}_{group}_{number:06d}'
        family=group.split('_',1)[0]
        if family in ('first','central','auto') and c.get('ray_station_naming','distance')=='distance':
            distance=float(np.linalg.norm(np.array([x,z])-g[family]))
            token=f'{distance:010.3f}'.replace('.','p')
            name=f"{lp}_{dict(first='F',central='C',auto='A')[family]}_{group.split('_',1)[1].replace('_','')}_D{token}"
        key=(round(float(x),6),round(float(z),6))
        if name in named_positions and named_positions[name]!=key:name+=f'n{number}'
        named_positions[name]=key
        if key not in lookup:
            lookup[key]=len(stations)
            stations.append([name,float(x),Y,float(z)])
            owners.append(group)
        else:
            idx=lookup[key]
            replacements.append([stations[idx][0],name,group,float(x),Y,float(z)])
            # Use the later sampled point, including its exact coordinates.
            stations[idx]=[name,float(x),Y,float(z)]
            owners[idx]=group
        memberships.append((group,lookup[key],name))
    count=0
    if station_group_enabled(c,'grid'):
        for z in samples(L,c['station_spacing_m']):
            for x in samples(L,c['station_spacing_m']):
                count+=1;add(x,z,'grid',count)
    for ray in g['rays']:
        if not station_group_enabled(c,ray['name']):continue
        for number,(x,z) in enumerate(segment_points(ray['a'],ray['b'],LINE_STEP),1):
            add(x,z,ray['name'],number)
    if station_group_enabled(c,'interface'):
        for number,x in enumerate(samples(L,LINE_STEP),1):
            add(x,H,'interface',number)
    # Every previous group now references the surviving name and coordinates.
    rows=[]
    for group,idx,requested in memberships:
        name,x,y,z=stations[idx]
        rows.append([group,name,x,y,z,L-z,requested])
    return stations,owners,rows,replacements


def write_output(c, g, out):
    time_shift=uniform_time_shift(c)
    if time_shift>0:
        warnings.warn(
            'All FORCE entries will use the same nonzero time shift. Official SPECFEM3D '
            'documentation requires at least one zero time shift; check solver compatibility. '
            'No individual source has been reset to zero.', UserWarning, stacklevel=2)
    out.mkdir(parents=True, exist_ok=True)
    sources = list(force_points(c,g))
    with (out/'FORCESOLUTION').open('w', encoding='ascii', newline='\n') as f:
        for i,(pos,d,cid) in enumerate(sources,1):
            x,y,z=pos
            depth=(L-z)/(1000 if c['force_depth_unit']=='km' else 1)
            frequency_label = 'hdurorf0:' if c['force_format']=='11_lines' else 'f0:'
            f.write(f'FORCE  {i:03d}\ntime shift:     {time_shift:.12g}\n{frequency_label}      {c["frequency_hz"]:.8f}\n'
                    f'latorUTM:      {y:.9f}\nlongorUTM:     {x:.9f}\ndepth:         {depth:.12f}\n')
            if c['force_format']=='11_lines':
                f.write('source time function:            1\n')
            f.write(f'factor force source:             {c["force_factor"]:.12e}\n'
                    f'component dir vect source E:     {d[0]:.12e}\n'
                    f'component dir vect source N:     {d[1]:.12e}\n'
                    f'component dir vect source Z_UP:  {d[2]:.12e}\n')
    write_csv(out/'sources_depth.csv', ['force_id','center_id','x_m','y_m','burial_m','E','N','Z_UP'],
              ([i,cid,pos[0],pos[1],L-pos[2],*d] for i,(pos,d,cid) in enumerate(sources,1)))
    stations,owners,memberships,replacements=build_stations(c,g)
    # 两类接收器必须全部在同一个剖面，禁止把理论线台站扩展为三维阵列。
    if any(station[2] != 500.0 for station in stations):
        raise ValueError('Every receiver must lie on the y=500 m plane.')
    for filename,selector in [('STATIONS',lambda group: True),
                              ('STATIONS_grid',lambda group: group=='grid'),
                              ('STATIONS_lines',lambda group: group!='grid')]:
        with (out/filename).open('w',encoding='ascii',newline='\n') as f:
            for (name,x,y,z),group in zip(stations,owners):
                if selector(group):
                    f.write(f'{name} DB {y:.9f} {x:.9f} 0.000000000 {L-z:.9f}\n')
    auto_names={row[1] for row in memberships if row[0].startswith('auto_')}
    with (out/'STATIONS_auto').open('w',encoding='ascii',newline='\n') as f:
        for name,x,y,z in stations:
            if name in auto_names:f.write(f'{name} DB {y:.9f} {x:.9f} 0.000000000 {L-z:.9f}\n')
    write_csv(out/'stations_depth.csv',['station','x_m','y_m','burial_m','station_class','owner_group'],
              ([name,x,y,L-z,'grid' if group=='grid' else 'line',group]
               for (name,x,y,z),group in zip(stations,owners)))
    write_csv(out/'station_groups.csv',['group','station','x_m','y_m','burial_m','original_name'],
              ([group,name,x,y,depth,original] for group,name,x,y,z,depth,original in memberships))
    write_csv(out/'station_name_replacements.csv',['old_name','new_name','new_group','x_m','y_m','burial_m'],
              ([old,new,group,x,y,L-z] for old,new,group,x,y,z in replacements))
    # Exact distance belongs to group membership, even when a later group renames a shared point.
    distance_rows=[]
    for group,name,x,y,z,depth,original in memberships:
        family=group.split('_',1)[0]
        if family in ('first','central','auto'):
            distance_rows.append([group,name,float(np.linalg.norm(np.array([x,z])-g[family])),original])
    write_csv(out/'ray_station_distances.csv',['group','station','distance_to_interface_m','requested_name'],distance_rows)
    write_csv(out/'auto_ray_candidates.csv',['hit_x_m','station_count','eligible'],
              ([r[k] for k in ('hit_x_m','station_count','eligible')] for r in g.get('auto_candidates',[])))
    report={k:g[k] for k in ('p','first_time','snapshot_time','residual','interface_footprint','branches')}
    report['branches']=copy.deepcopy(g['branches'])
    for branch in report['branches']:
        if 'direction_xz' in branch:
            dx,dz=branch.pop('direction_xz');branch['direction_x_depth']=[dx,-dz]
    report.update(first_contact_x_depth_m=x_depth(g['first']), central_contact_x_depth_m=x_depth(g['central']),
                  source_line_start_x_depth_m=x_depth(g['anchor']), source_line_end_x_depth_m=x_depth(g['end']),
                  coordinate_convention='All exported positions use burial depth in meters, positive downward.',
                  source_burial_range_m=[float(min(L-s[0][2] for s in sources)),float(max(L-s[0][2] for s in sources))],
                  central_is_source_midpoint=g['central_is_source_midpoint'],
                  incident_side=incident_side(c), effective_source_kind=source_kind(c),
                  uniform_time_shift_written_s=time_shift,
                  time_shift_compatibility_note=('All sources have a nonzero time shift; official SPECFEM requires at least one zero.'
                                                if time_shift>0 else 'All source time shifts are zero.'),
                  source_center_count=len(g['centers'])*len(samples(L,c['source_y_spacing_m'])),
                  force_count=len(sources),station_count=len(stations),
                  grid_station_count=owners.count('grid'), line_station_count=len(owners)-owners.count('grid'),
                  receiver_plane_y_m=Y, all_receivers_on_y500=True,
                  station_replacement_count=len(replacements),
                  station_generation=station_switches(c),
                  station_family_membership_counts={family:len({row[1] for row in memberships
                      if row[0].split('_',1)[0]==family}) for family in station_switches(c)},
                  lines=[dict(name=r['name'],a_x_depth_m=x_depth(r['a']),b_x_depth_m=x_depth(r['b'])) for r in g['rays']],
                  note='Times are phase-reference travel times; Ricker has no strict onset. '
                       'First contact is the finite plane-wave bundle limit, not first edge diffraction arrival.')
    if 'auto' in g:report.update(auto_contact_x_depth_m=x_depth(g['auto']),auto_selection=g['auto_selection'])
    (out/'theory_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    saved_config={**c,**{'generate_'+family+'_stations':enabled for family,enabled in station_switches(c).items()}}
    (out/'config_used.json').write_text(json.dumps(saved_config,ensure_ascii=False,indent=2),encoding='utf-8')
    acoustic_note = ('# Acoustic monopole: setup/constants.h must have DIPOLE_SOURCE_IN_FLUID = .false.\n'
                     '# This is a compile-time constant, NOT a Par_file setting.\n'
                     if source_kind(c)=='acoustic_monopole' else '')
    (out/'Par_file_required.txt').write_text(
        '# Merge these settings into YOUR Par_file. This is not a complete Par_file.\n'
        'SUPPRESS_UTM_PROJECTION = .true.\n'
        'USE_SOURCES_RECEIVERS_Z = .false.\n'
        'USE_FORCE_POINT_SOURCE = .true.\n'
        'USE_RICKER_TIME_FUNCTION = .true.\n'
        'USE_EXTERNAL_SOURCE_FILE = .false.\n'
        '# FORCE depth uses configured unit; STATIONS burial is always meters.\n'+acoustic_note,encoding='ascii')
    return sources,stations,report


def plot(c,g,sources,stations,out):
    """仅输出y=500的台站布置图及启用auto时的自动测线图。"""
    # 图的所有纵坐标都转换为埋深，向下为正；不改变计算几何。
    g=copy.deepcopy(g)
    for key in ('anchor','end','first','central')+ (('auto',) if 'auto' in g else ()):
        g[key]=np.asarray(x_depth(g[key]))
    for ray in g['rays']:
        ray['a']=np.asarray(x_depth(ray['a']));ray['b']=np.asarray(x_depth(ray['b']))
    stations=[[name,x,y,L-z] for name,x,y,z in stations]
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                         'axes.unicode_minus':False,'font.size':11,
                         'axes.titleweight':'bold','pdf.fonttype':42})
    colors={'I':'#c23b30','R_P':'#d58b14','R_SV':'#7755a2','T_P':'#007e99','T_SV':'#278353'}
    points=np.array([[s[1],s[3]] for s in stations],dtype=float).reshape(-1,2)
    grid=np.array([s[0].startswith(c.get('grid_station_prefix','G')) for s in stations],dtype=bool)
    gp,lp=c.get('grid_station_prefix','G'),c.get('line_station_prefix','L')
    def mode_info(key):
        action,wave=key.split('_')
        label={'I':'入射','R':'反射','T':'透射'}[action]+' '+wave+' 波'
        return label,colors['I'] if action=='I' else colors[key]
    def base(ax,source=True,labels=True):
        ax.axhspan(H,L,color='#f6f0e7',zorder=0)
        ax.axhspan(0,H,color='#eaf4fa',zorder=0)
        ax.axhline(H,color='#657888',lw=1.3,zorder=1)
        ax.set(xlim=(0,L),ylim=(L,0),xlabel='X / m',ylabel='埋深 / m（向下为正）')
        ax.set_aspect('equal');ax.set_xticks(np.arange(0,1001,200));ax.set_yticks(np.arange(0,1001,200))
        ax.grid(color='#97a6b0',alpha=.13,lw=.6)
        if labels:
            for key,z in [('upper',35),('lower',965)]:
                med=c[key]
                ax.text(18,z,f'{med["name"]}  |  Vp={med["vp"]:g}, Vs={med["vs"]:g} m/s',
                        fontsize=9,color='#51616d',va='center',bbox=dict(fc='white',ec='none',alpha=.8,pad=3))
            ax.text(980,H-13,'介质分界线：埋深500 m',ha='right',va='bottom',fontsize=9,color='#51616d',
                    bbox=dict(fc='white',ec='none',alpha=.85,pad=2))
        if source:
            ax.plot([g['anchor'][0],g['end'][0]],[g['anchor'][1],g['end'][1]],
                    color='#2c3743',lw=2.2,zorder=3)
    def save(fig,name):
        fig.savefig(out/(name+'.png'),dpi=PLOT_DPI,facecolor='white')
        fig.savefig(out/(name+'.pdf'),facecolor='white')
        plt.close(fig)
    def draw_ray(ax,r):
        key=r['name'].split('_',1)[1];label,color=mode_info(key)
        if not key.startswith('I_'):
            branch=next(b for b in g['branches'] if b['name']==key)
            label+=f' / {branch["angle_deg"]:.1f}°'
        a,b=r['a'],r['b'];d=b-a
        ax.plot([a[0],b[0]],[a[1],b[1]],color=color,lw=2.3,zorder=4)
        if np.linalg.norm(d)>1:
            ax.annotate('',xy=a+.64*d,xytext=a+.44*d,
                        arrowprops=dict(arrowstyle='-|>',color=color,lw=2.3,mutation_scale=18),zorder=5)
        return Line2D([0],[0],color=color,lw=2.5,label=label)

    if 'auto' in g:
        fig,ax=plt.subplots(figsize=(10,8))
        base(ax)
        xx,dd=np.meshgrid(np.linspace(0,L,250),np.linspace(0,L,250))
        # Convert displayed depths back to internal coordinates for beam test.
        vec=np.stack((xx, L-dd),axis=-1)-np.array([g['anchor'][0],L-g['anchor'][1]])
        transverse=vec@g['tangent'];forward=vec@g['normal']
        aperture=np.linalg.norm(g['end']-g['anchor'])
        medium=(dd<=H) if incident_side(c)=='upper' else (dd>=H)
        beam=(transverse>=0)&(transverse<=aperture)&(forward>=0)&medium
        ax.contourf(xx,dd,beam.astype(float),levels=[.5,1.5],colors=['#f4bd6a'],alpha=.35)
        handles=[]
        for r in g['rays']:
            if r['name'].startswith('auto_'):
                if r['name'].startswith('auto_R'):
                    rp=segment_points(r['a'],r['b'],LINE_STEP)
                    ax.scatter(rp[:,0],rp[:,1],s=10,c='#8b2030',zorder=8)
                handles.append(draw_ray(ax,r))
        ax.set_title('自动测线：整条反射线布站（包含橙色初始入射波束覆盖区）')
        ax.legend(handles=handles,fontsize=9)
        fig.tight_layout();save(fig,'auto_rays_y500')

    # 台站图：总览和真实接收器坐标局部放大，绝不沿y方向铺设接收器。
    if not math.isfinite(STATION_ZOOM_WIDTH_M) or not 0<STATION_ZOOM_WIDTH_M<=L:
        raise ValueError('STATION_ZOOM_WIDTH_M must be in (0,1000].')
    width=STATION_ZOOM_WIDTH_M
    zoom_reference=g['central'] if station_group_enabled(c,'central') or not station_group_enabled(c,'first') else g['first']
    center=np.clip(zoom_reference,width/2,L-width/2)
    x0,z0=center-width/2;x1,z1=center+width/2
    fig,axes=plt.subplots(1,2,figsize=(15,8.5))
    fig.subplots_adjust(left=.07,right=.98,bottom=.22,top=.85,wspace=.20)
    fig.suptitle('台站布置  |  每一个台站的 Y 坐标均为 500 m',fontsize=19,y=.97)
    enabled_names=[name for name,on in station_switches(c).items() if on]
    fig.text(.5,.92,'本次启用：'+('、'.join(enabled_names) if enabled_names else '无台站（仅生成震源与理论几何）'),ha='center',fontsize=12,color='#52616e')
    for ax in axes:
        base(ax,source=False,labels=False)
        ax.scatter(points[grid,0],points[grid,1],s=2.5,color='#9ca8b4',alpha=.7,rasterized=True,zorder=2)
        ax.scatter(points[~grid,0],points[~grid,1],s=5,color='#d65132',alpha=.9,rasterized=True,zorder=3)
    axes[0].set_title(f'A  全部台站：{len(stations):,} 个（去重后）',fontsize=13,pad=14)
    axes[0].add_patch(Rectangle((x0,z0),width,width,fill=False,ec='#212e3b',lw=1.6,zorder=5))
    axes[0].annotate('右图放大范围',xy=(x1,z1),xytext=(min(810,x1+100),min(920,z1+120)),
                     arrowprops=dict(arrowstyle='->',color='#344654'),fontsize=10,
                     bbox=dict(fc='white',ec='none',pad=3),zorder=6)
    zoom=axes[1];zoom.set(xlim=(x0,x1),ylim=(z1,z0))
    ticks=np.linspace(x0,x1,5);zticks=np.linspace(z0,z1,5)
    zoom.set_xticks(ticks);zoom.set_xticklabels([f'{x:.1f}' for x in ticks])
    zoom.set_yticks(zticks);zoom.set_yticklabels([f'{z:.1f}' for z in zticks])
    zoom.set_title(f'B  局部放大：{width:g} m × {width:g} m',fontsize=13,pad=14)
    inside=(points[:,0]>=x0)&(points[:,0]<=x1)&(points[:,1]>=z0)&(points[:,1]<=z1)
    zoom.scatter(points[inside&grid,0],points[inside&grid,1],s=25,color='#74899d',zorder=4)
    zoom.scatter(points[inside&~grid,0],points[inside&~grid,1],s=32,color='#d65132',edgecolors='white',linewidths=.5,zorder=5)
    # 仅标出一对10m网格与一对5m界面台站的间距，保持局部图易读。
    for wanted_grid,zshift,color in [(True,3.0,'#546b80'),(False,4.0,'#b64326')]:
        candidates=np.where(inside & (grid if wanted_grid else ~grid))[0]
        found=False
        target_z=center[1]-width*.25 if wanted_grid else center[1]
        for iz in sorted(set(np.round(points[candidates,1],6)),key=lambda z:abs(z-target_z)):
            ids=[i for i in candidates if abs(points[i,1]-iz)<1e-6]
            ids.sort(key=lambda i:points[i,0])
            for i,j in zip(ids[:-1],ids[1:]):
                distance=np.linalg.norm(points[j]-points[i])
                target=c['station_spacing_m'] if wanted_grid else LINE_STEP
                if abs(distance-target)<1e-5 and x0+8<points[i,0] and points[j,0]<x1-8 and z0+8<iz<z1-8:
                    a,b=points[i],points[j]
                    zoom.annotate('',xy=(b[0],iz+zshift),xytext=(a[0],iz+zshift),
                                  arrowprops=dict(arrowstyle='<->',color=color,lw=1.2),zorder=6)
                    zoom.text((a[0]+b[0])/2,iz+zshift+1.2,f'{target:g} m',
                              fontsize=9,ha='center',va='bottom',color=color,
                              bbox=dict(fc='white',ec='none',alpha=.9,pad=1),zorder=7)
                    found=True;break
            if found:break
    handles=[]
    if grid.any():handles.append(Line2D([0],[0],marker='o',ls='',color='#74899d',label=f'第一类：{gp}…  平面网格 / {c["station_spacing_m"]:g} m'))
    if (~grid).any():handles.append(Line2D([0],[0],marker='o',ls='',color='#d65132',label=f'第二类：{lp}_…  理论线加密 / 5 m'))
    handles.append(Line2D([0],[0],color='#657888',lw=1.3,label='埋深500 m：介质分界线（仍在y=500剖面内）'))
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.065),ncol=2,frameon=False,fontsize=11)
    fig.text(.5,.025,'重合位置仅保留一个台站，采用后生成的名称；STATIONS、坐标表与图件均采用埋深（m）。',ha='center',fontsize=10,color='#596774')
    save(fig,'stations_y500')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',type=Path,default=None,
                        help='可选：使用JSON文件替代脚本顶部USER_CONFIG；默认不读取JSON')
    default_output=Path(OUTPUT_DIRECTORY)
    if not default_output.is_absolute():
        default_output=Path(__file__).resolve().parent/default_output
    parser.add_argument('--output',type=Path,default=default_output,help='覆盖顶部输出目录')
    parser.add_argument('--no-plot',action='store_true',help='关闭绘图')
    args=parser.parse_args()
    c=(json.loads(args.config.read_text(encoding='utf-8-sig'))
       if args.config is not None else copy.deepcopy(USER_CONFIG))
    validate(c)
    g=geometry(c);sources,stations,report=write_output(c,g,args.output)
    if GENERATE_PLOTS and not args.no_plot: plot(c,g,sources,stations,args.output)
    print(json.dumps({k:report[k] for k in ('force_count','station_count','station_generation','station_family_membership_counts','source_burial_range_m','first_contact_x_depth_m','first_time')},ensure_ascii=False,indent=2))
    print('Output:',args.output.resolve())


if __name__=='__main__':
    if sys.version_info<(3,10): raise RuntimeError('Python 3.10 or newer required.')
    main()
