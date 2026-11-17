#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python 3.10：只用理论射线台站，从SPECFEM地震动提取反/透射系数模。
。
"""
from __future__ import annotations

# ======================== 用户编辑区（只改这里） ========================
SETUP_DIRECTORY = "../算例/Planewave_lower_20_SV/prepare"  # 本次模拟保存的生成器输出：包含config_used.json和theory_report.json
EXPECTED_INCIDENCE_ANGLE_DEG = 20
# 对比理论解：可填总目录（自动选upper_P/lower_P/lower_SV），也可直接填某个子目录。
# 例："../理论解/rt_coefficients_output/lower_SV"；绝对路径同样支持。
# None恢复原来的内置理论计算。只更换比较基准，不更改数值模拟的射线、到时或拾取窗。
# 理论表必须有本算例角度的行；不会静默取最近角度或跨临界角插值。
THEORY_DIRECTORY = "../理论解/rt_coefficients_output"
DATA_DIRECTORY = "../算例/Planewave_lower_20_SV/sema"               # .sema数据目录，可含子目录
OUTPUT_DIRECTORY = "output"                  # 相对于本脚本；同名结果重新生成
RAY_FAMILY = "auto"                          # auto自动测线；central/first；both旧两组；all读取实际存在的所有射线组
RUN_MODE = "analyze"                         # 有波形即可直接计算；prepare只是可选检查，无需先运行



NETWORK = "*"                               # 也可明确填"DB"或"SYN"
TRACE_QUANTITY = "acceleration"              # acceleration/velocity/displacement，所有相位必须同一量
FILE_EXTENSION = "sema"                     # 对应sema/semv/semd；不自动积分或混用声压
COMPONENT_CODES = {"X": "FXX", "Y": "FXY", "Z": "FXZ"}
Z_COMPONENT_POSITIVE = "up"                 # SPECFEM通常up；若数据已翻为向下，则填down
AMPLITUDE_SCALE = 1.0                        # 所有分量共同乘的单位换算系数，不做逐道归一化
TIME_COLUMN = 0                             # 从0计数，默认第1列
VALUE_COLUMN = 1                            # 默认第2列
SKIP_HEADER_LINES = 0
SAMPLE_TIME_RELATIVE_TOLERANCE = 0.005        # 兼容ASCII时间列舍入；更大不规则采样报错

# 时间窗围绕“源波峰在记录时间轴的位置 + 相位传播到时”。不要从记录第一行推算时间零点。
# 标准SPECFEM常以源中心为t=0；请结合实际源时程/波形核对，统一time shift不要重复加。
SOURCE_PEAK_TIME_S = 0.0                     # 包含实际源中心、整体时移；不是程序自动猜测值
TIME_REFERENCE_MODE = "record_time"          # record_time：上项直接是记录时间列中的源中心
# elapsed_from_start：记录从0重新计时，用SOURCE_PEAK_TIME_S-SOURCE_START_TIME_S
# stf_file：从两列源时程文件取绝对值峰位，再加STF_TO_RECORD_OFFSET_S
SOURCE_START_TIME_S = -0.012                 # 仅elapsed_from_start使用，不从主频自动猜起点
SOURCE_TIME_FUNCTION_FILE = None            # stf_file时必须提供；文件不存在会明确报错
STF_TO_RECORD_OFFSET_S = 0.0                 # 源时程时间轴到记录时间轴的显式偏移
PAIRING_RULE = "all_pairs"                 # 每个入射台站分别与每个出射台站配对；系数始终=出射幅值÷入射幅值
# 例如2个入射、3个反射台站产生6对组合；这里是组合数量相乘，不是幅值相乘。
# 也可选equal_distance/equal_travel_time，分别按等界面距离/等传播时间配对。
# 按各自测线到界面交点的沿线距离筛选，单位m，闭区间；None使用整条线。
# 例如"R_P": [100.0, 200.0]仅使用反射P线上距界面100–200m的台站。
STATION_DISTANCE_RANGES_M = {
    "I_P": [100.0, 200.0],                          # 入射P线；例如[100.0, 200.0]
    "R_P": [300.0, 600.0],                          # 反射P线
    "T_P": [200.0, 300.0],                          # 透射P线
    "T_SV": [200.0, 400.0],                         # 透射SV线
    "I_SV": [100.0, 200.0], "R_SV":[400.0, 500.0],           # 为SV入射情形保留
}
PAIRING_TOLERANCE_M = 1.6                    # 5m台站最近匹配容差；不插值/外推缺失台站
AUTO_EXPAND_SHORT_WINDOWS = True             # 短窗仅为诊断扩到足够采样点；强制标为不合格
MIN_WINDOW_SAMPLES = 16
MIN_SAMPLES_PER_PERIOD = 10.0                # 每主频周期少于此值：原始比值可算，不能通过检查
WINDOW_HALF_PERIODS = 5                   # 默认半窗宽=此值/f0；两侧可独立改
WINDOW_BEFORE_S = None                      # None用上项；显式值单位s
WINDOW_AFTER_S = None
WINDOW_CENTER_MODE = "theory"               # "theory"理论到时；"waveform"波形估计的候选到时
ARRIVAL_SEARCH_HALF_WIDTH_S = None           # waveform搜索理论到时前后多少秒；None=1个主频周期
# 台站用途由station_groups.csv决定，不根据文件名猜测，也不按哪个波最大来改变用途。
# None沿用上方全局设置；可按目标震相分别指定窗前宽、窗后宽、搜索半宽（单位秒）。
PHASE_WINDOW_SETTINGS = {
    "I_P":  {"before_s": None, "after_s": None, "search_half_width_s": None},  # 入射P台站
    "R_P":  {"before_s": None, "after_s": None, "search_half_width_s": None},  # 反射P台站
    "T_P":  {"before_s": None, "after_s": None, "search_half_width_s": None},  # 透射P台站
    "T_SV": {"before_s": None, "after_s": None, "search_half_width_s": None},  # 透射SV台站
    "I_SV": {"before_s": None, "after_s": None, "search_half_width_s": None},
    "R_SV": {"before_s": None, "after_s": None, "search_half_width_s": None},
}
LIMIT_WINDOW_TO_TARGET_PHASE = True          # 提取窗也受相邻震相到时中点限制，防止宽窗取到别的波
# 窗太窄时不跨越用途边界扩窗；同到时无法区分的波不强制给出系数。
# 手工CSV仍指定精确窗口，但若越过目标震相边界将报错，不静默修改手工窗。
# waveform：在极化投影波形的Hilbert包络上找搜索区内最高点，不是搜索整段记录。
# 为避免跳到相邻震相，搜索区还限制在目标与其他理论到时的中点之间。
# 原来的WINDOW_BEFORE_S/AFTER_S控制选定中心两侧窗宽；不改变台站配对规则。
# 手工CSV的start_s/end_s优先于此开关；估计失败退回理论中心并标记待核对。
WINDOW_OVERRIDES_CSV = None                 # 修改window_plan.csv后另存文件，在这里填路径
# override列：family,phase,station,start_s,end_s,use,reviewed
# use=0排除；reviewed=1表示已看图确认震相隔离，可接受计划中的重叠/边缘提示。
# 不会解除缺数据、零参考、窗口越界等硬错误。请不要把整张表盲目设为reviewed=1。

# 只处理选中射线组的理论台站；不读网格台站，也不以波阵面/界面混合波代替出射波。
MIN_INTERFACE_DISTANCE_M = 10.0             # 排除紧贴流固界面的接收器；可调
EXCLUDE_OUTER_BOUNDARY_STATIONS = False      # 当前无PML，不按外边界距离排除台站
MIN_OUTER_BOUNDARY_DISTANCE_M = 10.0        # 仅上项True时生效；不读取/假设PML厚度
MIN_SOURCE_DISTANCE_M = 50.0                # 入射参考接收器离源面的法向距离
OUTPUT_UNREVIEWED_RATIOS = True              # 输出待核对的原始振幅比，accepted仍为0；不会伪装成合格系数
# 若入射参考全被软提示标记，使用距界面最接近180m的非零参考计算原始比值。
# 距界面不足10m但不恰在界面上的理论台站也输出原始比值，保留too_close_interface提示。
# 缺数据、时间窗越界、零参考、用户use=0、源附近台站仍不强制计算。
REFERENCE_STATIONS = {"central": None, "first": None, "auto": None}  # 指定实际台站名；None从可用入射波中选
REFERENCE_TARGET_DISTANCE_M = 180.0          # 自动选择最靠近此“沿射线距界面距离”的合格入射台站
BASELINE_INTERVAL_S = None                  # 例如[-0.1,-0.02]，确认无震相后扣除该区间均值
DETREND = "none"                            # none/constant/linear；默认不改变原数据
BANDPASS_HZ = [80.0, 120.0]                          # 新算例默认不滤波；按实际频带设置，
FILTER_ORDER = 4
TAPER_FRACTION = 0.10                       # Hilbert前只在所选窗两端做总计10%的余弦渐消
HILBERT_PAD_FACTOR = 4                      # 双侧补零后FFT，降低周期绕回；不能修复震相重叠
EDGE_AMPLITUDE_LIMIT = 0.20                 # 窗边界幅值/峰值大于此值则标记可能截波
MIN_REFERENCE_AMPLITUDE = 1e-30             # 防止零幅值除法；量纲与输入一致
GENERATE_PLOTS = False                      # 额外过程诊断图默认关闭；最终汇总PNG始终生成
MAX_WAVEFORM_PLOTS_PER_PHASE = 3            # 每组每相位最多几张诊断图，含参考台站
PLOT_HILBERT_ERROR_COMBINATIONS = True      # 旧配置兼容项；最终达标组合已合并为唯一汇总PNG
HILBERT_ERROR_THRESHOLD_PCT = 5.0           # 严格小于此百分数；默认5%
HILBERT_ERROR_ACCEPTED_ONLY = False         # False显示所有有数值组合；图中区分合格与待核对
BEST_PAIR_METHOD = "hilbert"                # 最佳组合排序方法：hilbert包络法，或direct直接峰值法
BEST_PAIR_ACCEPTED_ONLY = True              # True只从通过质量检查的组合选；无合格组合时明确标注，不以待核对替代
BEST_PICK_CONTEXT_PERIODS = 1.0             # 旧配置兼容项，已不用于裁剪；最佳台站固定显示完整记录
# 每个射线组、每种出射波各选一对：当前上层P入射对应反射P、透射P、透射SV共三对。
# 按|数值系数模-理论系数模|/理论系数模最小选取；三对可以使用不同入射台站。
PLOT_DPI = 170
PLOT_NOTE = ""                             # 图件附注；真实数据可填算例名
# =========================== 用户编辑区结束 ============================

import argparse
import copy
import csv
import fnmatch
import json
import math
import re
import os
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.signal import butter, detrend, hilbert, sosfiltfilt
from scipy.signal.windows import tukey
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
MODEL_SIZE = 1000.0
INTERFACE_DEPTH = 500.0
GENERATED_PLOT_PATHS = []


def settings():
    """复制用户设置，便于命令行和验证程序覆盖，计算过程不修改模块全局量。"""
    names = ("SETUP_DIRECTORY DATA_DIRECTORY OUTPUT_DIRECTORY RAY_FAMILY RUN_MODE NETWORK "
             "TRACE_QUANTITY FILE_EXTENSION COMPONENT_CODES Z_COMPONENT_POSITIVE AMPLITUDE_SCALE "
             "TIME_COLUMN VALUE_COLUMN SKIP_HEADER_LINES SAMPLE_TIME_RELATIVE_TOLERANCE "
             "SOURCE_PEAK_TIME_S WINDOW_HALF_PERIODS WINDOW_BEFORE_S WINDOW_AFTER_S WINDOW_OVERRIDES_CSV "
             "MIN_INTERFACE_DISTANCE_M MIN_OUTER_BOUNDARY_DISTANCE_M MIN_SOURCE_DISTANCE_M EXCLUDE_OUTER_BOUNDARY_STATIONS OUTPUT_UNREVIEWED_RATIOS "
             "REFERENCE_STATIONS REFERENCE_TARGET_DISTANCE_M BASELINE_INTERVAL_S DETREND BANDPASS_HZ "
             "FILTER_ORDER TAPER_FRACTION HILBERT_PAD_FACTOR EDGE_AMPLITUDE_LIMIT "
             "MIN_REFERENCE_AMPLITUDE GENERATE_PLOTS MAX_WAVEFORM_PLOTS_PER_PHASE PLOT_DPI PLOT_NOTE").split()
    names += ("EXPECTED_INCIDENCE_ANGLE_DEG THEORY_DIRECTORY TIME_REFERENCE_MODE SOURCE_START_TIME_S SOURCE_TIME_FUNCTION_FILE STF_TO_RECORD_OFFSET_S "
              "STATION_DISTANCE_RANGES_M PAIRING_RULE PAIRING_TOLERANCE_M AUTO_EXPAND_SHORT_WINDOWS MIN_WINDOW_SAMPLES MIN_SAMPLES_PER_PERIOD "
              "WINDOW_CENTER_MODE ARRIVAL_SEARCH_HALF_WIDTH_S PHASE_WINDOW_SETTINGS LIMIT_WINDOW_TO_TARGET_PHASE "
              "PLOT_HILBERT_ERROR_COMBINATIONS HILBERT_ERROR_THRESHOLD_PCT HILBERT_ERROR_ACCEPTED_ONLY BEST_PAIR_METHOD BEST_PAIR_ACCEPTED_ONLY BEST_PICK_CONTEXT_PERIODS").split()
    return {name.lower(): copy.deepcopy(globals()[name]) for name in names}


def resolve(path):
    p = Path(path).expanduser()
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, fields):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def state(m, wave, p, sign):
    """位移/牵引状态，竖向向下；与生成器的理论极化基底一致。"""
    v = m["vp"] if wave == "P" else m["vs"]
    q = sign * np.sqrt(complex(1 / v**2 - p**2))
    if m["type"] == "fluid":
        return np.array([v*p, v*q, 0j, m["rho"]*v])
    mu = m["rho"]*m["vs"]**2
    lam = m["rho"]*(m["vp"]**2-2*m["vs"]**2)
    if wave == "P":
        return np.array([v*p, v*q, 2*mu*v*p*q, lam/v+2*mu*v*q*q])
    return np.array([-v*q, v*p, mu*v*(p*p-q*q), 2*mu*v*p*q])


@dataclass
class Model:
    config: dict
    report: dict
    groups: dict
    stations: dict
    phases: dict
    source: np.ndarray
    normal: np.ndarray
    tangent: np.ndarray
    aperture: float
    incident_speed: float


def load_model(folder):
    folder = Path(folder)
    c = json.loads((folder/"config_used.json").read_text(encoding="utf-8-sig"))
    r = json.loads((folder/"theory_report.json").read_text(encoding="utf-8-sig"))
    if "source_line_start_x_depth_m" not in r:
        raise ValueError("需要最新版埋深格式的theory_report.json，请用生成器重新生成。")
    side = c["incident_side"]
    inc, trans = c[side], c["lower" if side == "upper" else "upper"]
    sign = 1 if side == "upper" else -1
    angle = math.radians(c["angle_deg"])
    normal = np.array([math.sin(angle), sign*math.cos(angle)])
    tangent = np.array([-math.cos(angle), sign*math.sin(angle)])
    vi = inc["vp"] if c["incident_wave"] == "P" else inc["vs"]
    p = normal[0]/vi
    if r["incident_side"] != side or not np.isclose(p, r["p"], rtol=1e-9, atol=1e-14):
        raise ValueError("config_used与theory_report不是同一个算例。")
    labels = [("R_"+w, inc, w, -sign) for w in (["P"] if inc["type"] == "fluid" else ["P", "SV"])]
    labels += [("T_"+w, trans, w, sign) for w in (["P"] if trans["type"] == "fluid" else ["P", "SV"])]
    indices = [0,1,2,3] if inc["type"] == trans["type"] == "solid" else ([1,3] if inc["type"] == trans["type"] == "fluid" else [1,3,2])
    mat = np.column_stack([(-1 if name.startswith("R") else 1)*state(m,w,p,s)[indices] for name,m,w,s in labels])
    rhs = state(inc,c["incident_wave"],p,sign)[indices]
    scale = np.max(abs(mat),axis=1)
    coef = np.linalg.solve(mat/scale[:,None],rhs/scale)
    phases = {"I_"+c["incident_wave"]: dict(wave=c["incident_wave"], medium=inc, direction=normal,
                                               speed=vi, status="propagating", theory_abs=1.0)}
    for (name,m,w,s), value in zip(labels,coef):
        stored = next(b for b in r["branches"] if b["name"] == name)
        if not np.isclose(value,complex(stored["coefficient_real"],stored["coefficient_imag"]),rtol=1e-7,atol=1e-10):
            raise ValueError("介质参数与理论报告系数不一致，请使用同次生成的文件。")
        v = m["vp"] if w == "P" else m["vs"]
        phases[name] = dict(wave=w, medium=m, speed=v, status=stored["status"], theory_abs=float(abs(value)),
                            theory_real=float(value.real), theory_imag=float(value.imag),
                            direction=np.array([p*v,s*math.sqrt(max(0,1-(p*v)**2))]))
    stations = {}
    for line in (folder/"STATIONS_lines").read_text(encoding="utf-8-sig").splitlines():
        if not line.strip() or line.startswith("#"): continue
        name,network,y,x,elevation,depth = line.split()
        if float(y) != 500.0 or float(elevation) != 0.0:
            raise ValueError("STATIONS_lines必须全部为y=500，且第5列为0的埋深格式。")
        if name in stations: raise ValueError("理论台站名称重复："+name)
        stations[name] = dict(name=name,network=network,x_m=float(x),y_m=float(y),burial_m=float(depth))
    groups = {}
    for row in read_csv(folder/"station_groups.csv"):
        group = row["group"]
        if not group.startswith(("first_", "central_", "auto_")): continue
        name = row["station"]  # 必须使用覆盖后的名字，不能用original_name或从前缀猜归属。
        if name not in stations: raise ValueError("射线台站不在STATIONS_lines："+name)
        st = stations[name]
        if any(abs(float(row[k])-st[k]) > 2e-6 for k in ("x_m","y_m","burial_m")):
            raise ValueError("分组表与STATIONS_lines坐标不一致："+name)
        groups.setdefault(group, {})[name] = st
    a = np.asarray(r["source_line_start_x_depth_m"],float)
    b = np.asarray(r["source_line_end_x_depth_m"],float)
    if not np.allclose((b-a)/max(np.linalg.norm(b-a),1e-15),tangent,atol=1e-7):
        raise ValueError("源面几何与配置不一致。")
    distance_file=folder/'ray_station_distances.csv'
    if distance_file.exists():
        for row in read_csv(distance_file):
            group,name=row['group'],row['station']
            if group not in groups:continue
            if name not in groups[group]:raise ValueError('距离表台站不属于对应射线组：'+name)
            family=group.split('_',1)[0]
            if family+'_contact_x_depth_m' not in r:raise ValueError('距离表缺少对应界面交点：'+family)
            st=stations[name]
            actual=np.linalg.norm(np.array([st['x_m'],st['burial_m']])-r[family+'_contact_x_depth_m'])
            if abs(actual-float(row['distance_to_interface_m']))>2e-6:raise ValueError('距离表与台站坐标不一致：'+name)
    return Model(c,r,groups,stations,phases,a,normal,tangent,float(np.linalg.norm(b-a)),vi)


def select_comparison_theory(model,o):
    """Read the user's exported displacement-amplitude table, without changing geometry."""
    directory=o.get('theory_directory')
    if directory is None:
        return dict(source='internal',description='内置界面方程（生成器参数）',table_file=None,
                    parameter_file=None,selected_angle_deg=float(model.config['angle_deg']))
    folder=resolve(directory)
    if not (folder/'reflection_transmission_coefficients.txt').is_file():
        folder=folder/(model.config['incident_side']+'_'+model.config['incident_wave'])
    table=folder/'reflection_transmission_coefficients.txt'
    params=folder/'model_parameters.txt'
    if not table.is_file() or not params.is_file():
        raise ValueError('理论解目录须包含reflection_transmission_coefficients.txt和model_parameters.txt：'+str(folder))
    text=params.read_text(encoding='utf-8-sig')
    def required(pattern,label):
        match=re.search(pattern,text,re.MULTILINE)
        if match is None:raise ValueError('理论参数文件缺少或无法识别'+label+'：'+str(params))
        return match
    side=int(required(r'^INCIDENT_SIDE\s*=\s*([12])\s*$','INCIDENT_SIDE').group(1))
    wave=required(r'^INCIDENT_WAVE\s*=\s*(P|SV)\s*$','INCIDENT_WAVE').group(1)
    actual_side='upper' if side==1 else 'lower'
    c=model.config
    if (actual_side,wave)!=(c['incident_side'],c['incident_wave']):
        raise ValueError(f"所选理论解为{actual_side}_{wave}，数值算例为{c['incident_side']}_{c['incident_wave']}；请匹配理论目录与SETUP_DIRECTORY，不能用不同入射条件计算验证误差。")
    for number,key in ((1,'upper'),(2,'lower')):
        line=required(r'^MEDIA_'+str(number)+r':\s*(.+)$','MEDIA_'+str(number)).group(1)
        kind=re.search(r'\btype=(solid|fluid)\b',line)
        if kind is None or kind.group(1)!=c[key]['type']:raise ValueError('所选理论解介质类型与算例不一致：'+key)
        for prop in ('rho','vp','vs'):
            match=re.search(r'\b'+prop+r'=([^\s,]+)',line)
            value=float(match.group(1)) if match else float('nan')
            if not np.isfinite(value) or not np.isclose(value,c[key][prop],rtol=1e-9,atol=1e-9):
                raise ValueError(f"所选理论解{key}.{prop}={value}与算例{c[key][prop]}不一致。")
    lines=table.read_text(encoding='utf-8-sig').splitlines()
    header=next((line for line in lines if line.startswith('# 入射：')),None)
    if header and f'介质{side}，{wave}波' not in header:raise ValueError('理论表头与model_parameters.txt入射条件不一致。')
    reader=csv.DictReader((line for line in lines if line.strip() and not line.lstrip().startswith('#')),delimiter='\t')
    fields=reader.fieldnames or []
    expected=['angle_deg']
    for phase in model.phases:
        if phase.startswith('I_'):continue
        medium=side if phase.startswith('R_') else 3-side
        expected.extend(phase+'_M'+str(medium)+'_'+suffix for suffix in ('real','imag','abs'))
    if len(fields)!=len(set(fields)) or any(key not in fields for key in expected):
        raise ValueError('理论表缺少所需角度/模态列或列名重复；需要：'+', '.join(expected))
    rows=list(reader)
    if not rows:raise ValueError('理论系数表没有数据行。')
    try:angles=np.array([float(row['angle_deg']) for row in rows])
    except (ValueError,TypeError) as exc:raise ValueError('理论表angle_deg列存在非法数值。') from exc
    if not np.isfinite(angles).all() or np.any(angles<0) or np.any(angles>=90) or len(np.unique(angles))!=len(angles):
        raise ValueError('理论表角度须为[0,90)内的有限数且不能重复。')
    target=float(c['angle_deg']);matches=np.flatnonzero(np.isclose(angles,target,rtol=0,atol=1e-7))
    if len(matches)!=1:
        raise ValueError(f'理论表没有唯一的{target:g}°记录（范围{min(angles):g}～{max(angles):g}°）；请重新导出包含该角度的理论表，不会用最近角度替代。')
    row=rows[int(matches[0])];values={}
    for phase in model.phases:
        if phase.startswith('I_'):continue
        medium=side if phase.startswith('R_') else 3-side
        prefix=phase+'_M'+str(medium)+'_'
        try:real,imag,absolute=(float(row[prefix+k]) for k in ('real','imag','abs'))
        except (ValueError,TypeError) as exc:raise ValueError('理论系数不是有效数值：'+phase) from exc
        if not np.isfinite([real,imag,absolute]).all() or absolute<0 or not np.isclose(abs(complex(real,imag)),absolute,rtol=1e-8,atol=1e-10):
            raise ValueError('理论表的复系数与系数模不一致或非有限：'+phase)
        values[phase]=dict(theory_real=real,theory_imag=imag,theory_abs=absolute)
    # Commit only after every reference value has passed validation.
    for phase,values_for_phase in values.items():model.phases[phase].update(values_for_phase)
    return dict(source='table',description='外部理论表：'+folder.name,
                directory=str(folder),table_file=str(table),parameter_file=str(params),
                selected_angle_deg=float(angles[int(matches[0])]),angle_selection='exact',
                amplitude_definition='位移极化方向复振幅系数的模；与同物理量的出射/入射幅值比比较')


def theory_context(model,setup,o):
    """明确展示理论输入；用户指定角度只做一致性核对，不覆盖台站几何。"""
    c=model.config;expected=o.get('expected_incidence_angle_deg')
    if expected is not None:
        if not np.isfinite(expected) or not 0<=expected<90:
            raise ValueError('EXPECTED_INCIDENCE_ANGLE_DEG应为[0,90)内的角度或None')
        if not np.isclose(float(expected),c['angle_deg'],rtol=0,atol=1e-8):
            raise ValueError(f"你指定的算例角度为{expected:g}°，但SETUP_DIRECTORY的config_used.json为{c['angle_deg']:g}°。请选择本次数值模拟对应的生成器输出，不能只改理论角度。")
    comparison=select_comparison_theory(model,o)
    model.comparison_theory=comparison
    return dict(config_file=str(setup/'config_used.json'),report_file=str(setup/'theory_report.json'),
                comparison_theory=comparison,
                angle_deg=c['angle_deg'],angle_definition='入射传播方向与界面法线的夹角，单位度',
                incident_side=c['incident_side'],incident_wave=c['incident_wave'],
                upper=c['upper'],lower=c['lower'],
                coefficients={ph:{k:info.get(k,0.0) for k in ('theory_abs','theory_real','theory_imag','status')}
                              for ph,info in model.phases.items() if not ph.startswith('I_')},
                verification_limit='射线与到时由生成器元数据确定；比较基准见comparison_theory。外部表按入射侧、波型、介质参数及实际角度校验；未从波形独立反演这些参数。')


def arrival(model, phase, point):
    """有限孔径理想平面波相位到时；无覆盖返回None，非点源球面距离。"""
    point = np.asarray(point,float)
    if phase.startswith("I_"):
        length = float(np.dot(point-model.source,model.normal))
        foot = point-length*model.normal
        s = float(np.dot(foot-model.source,model.tangent))
        inc_upper = model.config["incident_side"] == "upper"
        if (point[1] > INTERFACE_DEPTH+1e-6 if inc_upper else point[1] < INTERFACE_DEPTH-1e-6): return None
        return length/model.incident_speed if length >= -1e-6 and -1e-6 <= s <= model.aperture+1e-6 else None
    info = model.phases[phase]
    if info["status"] != "propagating": return None
    nx,nd = info["direction"]
    if abs(nd) < 1e-12: return None
    length = (point[1]-INTERFACE_DEPTH)/nd
    hit_x = point[0]-length*nx
    xmin,xmax = model.report["interface_footprint"]
    if length < -1e-6 or not xmin-1e-5 <= hit_x <= xmax+1e-5: return None
    return model.report["first_time"]+model.report["p"]*(hit_x-xmin)+length/info["speed"]


def polarization(info):
    nx,nd = info["direction"]
    return np.array([nx,nd]) if info["wave"] == "P" else np.array([-nd,nx])


def validate_options(o):
    if not np.isfinite(o['best_pick_context_periods']) or o['best_pick_context_periods']<0:raise ValueError('BEST_PICK_CONTEXT_PERIODS必须非负')
    if o['best_pair_method'] not in ('hilbert','direct'):raise ValueError('BEST_PAIR_METHOD只能是hilbert或direct')
    if not isinstance(o['best_pair_accepted_only'],bool):raise ValueError('BEST_PAIR_ACCEPTED_ONLY必须为True或False')
    if not np.isfinite(o["hilbert_error_threshold_pct"]) or o["hilbert_error_threshold_pct"]<=0:raise ValueError("误差阈值必须为正数")
    if o['window_center_mode'] not in ('theory','waveform'): raise ValueError('WINDOW_CENTER_MODE只能theory或waveform')
    width=o['arrival_search_half_width_s']
    if width is not None and (not np.isfinite(width) or width<=0): raise ValueError('ARRIVAL_SEARCH_HALF_WIDTH_S必须为正数或None')
    for phase, params in o['phase_window_settings'].items():
        if phase not in ('I_P','R_P','T_P','I_SV','R_SV','T_SV'): raise ValueError('未知目标震相：'+phase)
        for key,value in params.items():
            if key not in ('before_s','after_s','search_half_width_s'): raise ValueError('未知震相窗口参数：'+key)
            if value is not None and (not np.isfinite(value) or value<=0): raise ValueError(phase+'的'+key+'必须为正数或None')
    if o['time_reference_mode'] not in ('record_time','elapsed_from_start','stf_file'): raise ValueError('未知TIME_REFERENCE_MODE')
    if o['pairing_rule'] not in ('all_pairs','equal_travel_time','equal_distance'): raise ValueError('未知PAIRING_RULE')
    for phase,bounds in o['station_distance_ranges_m'].items():
        if phase not in ('I_P','R_P','T_P','I_SV','R_SV','T_SV'):raise ValueError('未知距离筛选震相：'+phase)
        if bounds is not None and (not isinstance(bounds,(list,tuple)) or len(bounds)!=2 or
            not np.isfinite(bounds).all() or bounds[0]<0 or bounds[1]<bounds[0]):
            raise ValueError(phase+'距离范围应为None或[非负下限,不小于下限的上限]，单位m')
    if o['pairing_tolerance_m']<0 or o['min_window_samples']<4 or o['min_samples_per_period']<=0: raise ValueError('配对/采样参数无效')
    for key in ('exclude_outer_boundary_stations','output_unreviewed_ratios','limit_window_to_target_phase'):
        if not isinstance(o[key],bool): raise ValueError(key+'必须为True或False')
    if o["ray_family"] not in ("central","first","auto","both","all"): raise ValueError("RAY_FAMILY只能central/first/auto/both/all")
    if o["run_mode"] not in ("prepare","analyze"): raise ValueError("RUN_MODE只能prepare/analyze")
    expected = {"acceleration":"sema","velocity":"semv","displacement":"semd"}
    if expected.get(o["trace_quantity"]) != o["file_extension"]:
        raise ValueError("输入物理量与扩展名不一致；不支持把声压当位移或加速度相除。")
    if o["z_component_positive"] not in ("up","down"): raise ValueError("Z_COMPONENT_POSITIVE只能up/down")
    if o["detrend"] not in ("none","constant","linear"): raise ValueError("未知DETREND")
    for key in ("window_half_periods","amplitude_scale","hilbert_pad_factor","min_reference_amplitude"):
        if not np.isfinite(o[key]) or o[key] <= 0: raise ValueError(key+"必须为有限正数")
    for key in ("min_interface_distance_m","min_outer_boundary_distance_m","min_source_distance_m"):
        if not np.isfinite(o[key]) or o[key] < 0: raise ValueError(key+"必须为有限非负数")
    if not 0 <= o["taper_fraction"] < 1: raise ValueError("TAPER_FRACTION须在[0,1)")
    if not np.isfinite(o["source_peak_time_s"]): raise ValueError("源中心时间须有限")


def load_overrides(path):
    if path is None: return {}
    result = {}
    for row in read_csv(resolve(path)):
        key = (row["family"],row["phase"],row["station"])
        if key in result: raise ValueError("时间窗覆盖表有重复项："+str(key))
        result[key] = row
    return result


PLAN_FIELDS = ["family","phase","station","x_m","y_m","burial_m","distance_to_interface_m",
               "predicted_peak_s","start_s","end_s","use","reviewed","eligible","geometry_note","other_arrivals_s","manual_window",
               "role_group","target_phase","station_role","search_half_width_s","target_gate_start_s","target_gate_end_s",
               "coincident_phases","station_role_count","incident_beam_region","outside_beam_selection_margin","distance_range_selected"]


def phase_setting(o,phase,key,default):
    value=o['phase_window_settings'].get(phase,{}).get(key)
    return default if value is None else value


def phase_label(phase):
    return {'I_P':'入射P','R_P':'反射P','T_P':'透射P','I_SV':'入射SV','R_SV':'反射SV','T_SV':'透射SV'}[phase]


def target_phase_gate(theory,others):
    """时间上的用途边界；只是隔离已知震相中心，不能消除重叠波包。"""
    lo,hi=-np.inf,np.inf
    coincident=[]
    for phase,other in others.items():
        if abs(other-theory)<=1e-10: coincident.append(phase)
        elif other<theory:lo=max(lo,(theory+other)/2)
        else:hi=min(hi,(theory+other)/2)
    return lo,hi,coincident


def selected_families(model,o):
    available=[f for f in ('first','central','auto') if any(k.startswith(f+'_') for k in model.groups)]
    requested=available if o['ray_family']=='all' else (['central','first'] if o['ray_family']=='both' else [o['ray_family']])
    missing=[f for f in requested if f not in available]
    if missing or not requested:
        raise ValueError('未生成所选射线组：'+','.join(missing)+'；当前可用：'+','.join(available)+'。检查生成器开关或RAY_FAMILY。')
    for f in requested:
        if f+'_contact_x_depth_m' not in model.report:raise ValueError('理论报告缺少'+f+'界面交点，请使用同次生成的完整元数据。')
    return requested


def make_plan(model,o):
    families = selected_families(model,o)
    before = o["window_before_s"] if o["window_before_s"] is not None else o["window_half_periods"]/model.config["frequency_hz"]
    after = o["window_after_s"] if o["window_after_s"] is not None else o["window_half_periods"]/model.config["frequency_hz"]
    if not np.isfinite([before,after]).all() or min(before,after) <= 0: raise ValueError("时间窗宽必须为正")
    overrides = load_overrides(o["window_overrides_csv"])
    plan = []
    for family in families:
        if family == "central":
            middle = model.source+0.5*model.aperture*model.tangent
            hit = middle+(INTERFACE_DEPTH-middle[1])/model.normal[1]*model.normal
            if not np.allclose(hit,model.report["central_contact_x_depth_m"],atol=1e-5):
                raise ValueError("当前central组不是源面几何中点：旧版可能使用裁剪后有效孔径中点，或真正中点射线已越界。请选择first，或用更新后的生成器准备新算例；已完成模拟的元数据不能单独重生成替换。")
        hit = np.array(model.report[family+"_contact_x_depth_m"])
        for phase,info in model.phases.items():
            if info["status"] != "propagating": continue
            phase_before=phase_setting(o,phase,'before_s',before)
            phase_after=phase_setting(o,phase,'after_s',after)
            search_width=phase_setting(o,phase,'search_half_width_s',o['arrival_search_half_width_s'])
            if search_width is None:search_width=1/model.config['frequency_hz']
            for name,station in model.groups.get(family+"_"+phase,{}).items():
                point = np.array([station["x_m"],station["burial_m"]])
                distance=float(np.linalg.norm(point-hit))
                bounds=o['station_distance_ranges_m'].get(phase)
                range_selected=bounds is None or bounds[0]-1e-6<=distance<=bounds[1]+1e-6
                reasons = [] if range_selected else ['outside_selected_distance_range']
                if abs(point[1]-INTERFACE_DEPTH) < o["min_interface_distance_m"]-1e-6: reasons.append("too_close_interface")
                if o['exclude_outer_boundary_stations'] and min(*point,*(MODEL_SIZE-point)) < o["min_outer_boundary_distance_m"]-1e-6: reasons.append("too_close_outer_boundary")
                if phase.startswith("I_") and np.dot(point-model.source,model.normal) < o["min_source_distance_m"]-1e-6: reasons.append("too_close_source")
                travel = arrival(model,phase,point)
                if travel is None: reasons.append("outside_finite_aperture")
                peak = o["source_peak_time_s"]+(travel if travel is not None else 0)
                others = {ph:o["source_peak_time_s"]+t for ph in model.phases if ph != phase
                          if (t:=arrival(model,ph,point)) is not None}
                gate_lo,gate_hi,coincident=target_phase_gate(peak,others)
                transverse=float(np.dot(point-model.source,model.tangent))
                region=('覆盖区内' if arrival(model,'I_'+model.config['incident_wave'],point) is not None else '覆盖区外')
                margin=model.config.get('auto_beam_clearance_m',10.)
                outside_margin=transverse < -margin-1e-8 or transverse > model.aperture+margin+1e-8
                row = dict(distance_range_selected=int(range_selected),incident_beam_region=region,outside_beam_selection_margin=int(outside_margin),family=family,phase=phase,station=name,**{k:station[k] for k in ("x_m","y_m","burial_m")},
                           distance_to_interface_m=float(np.linalg.norm(point-hit)),predicted_peak_s=peak,
                           start_s=peak-phase_before,end_s=peak+phase_after,use=1,reviewed=0,eligible=int(not reasons),
                           geometry_note=";".join(reasons),other_arrivals_s=json.dumps(others,sort_keys=True),manual_window=0,
                           role_group=family+'_'+phase,target_phase=phase,station_role='仅拾取'+phase_label(phase),
                           search_half_width_s=search_width,target_gate_start_s=gate_lo if np.isfinite(gate_lo) else '',
                           target_gate_end_s=gate_hi if np.isfinite(gate_hi) else '',coincident_phases=';'.join(coincident))
                override = overrides.pop((family,phase,name),{})
                row['manual_window']=int(any(override.get(key,'')!='' for key in ('start_s','end_s')))
                for key in ("start_s","end_s"):
                    if override.get(key,"") != "": row[key] = float(override[key])
                for key in ("use","reviewed"):
                    if override.get(key,"") != "": row[key] = int(override[key])
                if not np.isfinite([row["start_s"],row["end_s"]]).all() or row["end_s"] <= row["start_s"]:
                    raise ValueError("时间窗无效："+name)
                if row["use"] not in (0,1) or row["reviewed"] not in (0,1): raise ValueError("use/reviewed必须0或1")
                plan.append(row)
    if overrides: raise ValueError("覆盖表含不属于当前射线组的条目："+str(next(iter(overrides))))
    if not plan: raise ValueError("当前元数据没有选中射线组；请检查station_groups.csv。")
    for row in plan:row['station_role_count']=sum(r['station']==row['station'] for r in plan)
    return plan


class TraceIndex:
    def __init__(self,folder,o):
        self.o = o
        self.files = {}
        self.count = 0
        if Path(folder).exists():
            for p in sorted(Path(folder).rglob("*."+o["file_extension"])):
                fields = p.name.split(".")
                if len(fields) != 4: continue
                network,station,component,_ = fields
                if fnmatch.fnmatchcase(network,o["network"]):
                    self.files.setdefault((station,component),[]).append(p)
                    self.count += 1

    def locate(self,station,component,optional=False):
        found = self.files.get((station,self.o["component_codes"][component]),[])
        if len(found)>1: raise ValueError("duplicate_component:"+station+":"+component+"，请指定NETWORK或缩小数据目录")
        if not found:
            if optional: return None
            raise FileNotFoundError("missing_component:"+station+":"+component)
        return found[0]


def read_trace(path,o):
    # Fortran D指数也可读；ASCII文件不会根据文件名推测坐标。
    with Path(path).open(encoding="utf-8-sig") as f:
        lines = (line.replace("D","E").replace("d","e") for i,line in enumerate(f) if i>=o["skip_header_lines"])
        a = np.loadtxt(lines,ndmin=2)
    if a.shape[0]<8 or a.shape[1]<=max(o["time_column"],o["value_column"]): raise ValueError("invalid_trace_shape:"+str(path))
    t = a[:,o["time_column"]]
    y = a[:,o["value_column"]]*o["amplitude_scale"]
    if not np.isfinite(t).all() or not np.isfinite(y).all(): raise ValueError("nonfinite_trace:"+str(path))
    dt = np.diff(t)
    step = float((t[-1]-t[0])/(len(t)-1))
    if np.any(dt<=0) or np.max(abs(dt-step))>o["sample_time_relative_tolerance"]*step: raise ValueError("nonuniform_time:"+str(path))
    # Hilbert FFT要求等间隔。仅纠正容差内的打印舍入，不补缺口、不外推。
    uniform = np.linspace(t[0],t[-1],len(t))
    y = np.interp(uniform,t,y)
    return uniform,y


def preprocess(t,y,o):
    y = np.array(y,copy=True)
    if o["baseline_interval_s"] is not None:
        a,b = o["baseline_interval_s"]
        if a<t[0] or b>t[-1] or a>=b: raise ValueError("invalid_baseline_interval")
        mask = (t>=a)&(t<=b)
        if mask.sum()<3: raise ValueError("baseline_too_short")
        y -= np.mean(y[mask])
    if o["detrend"] != "none": y = detrend(y,type=o["detrend"])
    if o["bandpass_hz"] is not None:
        low,high = o["bandpass_hz"]
        fs = 1/(t[1]-t[0])
        if not 0<low<high<fs/2: raise ValueError("bandpass_outside_nyquist")
        y = sosfiltfilt(butter(o["filter_order"],[low,high],btype="bandpass",fs=fs,output="sos"),y)
    return y


def project_trace(index,station,info,o):
    tx,x = read_trace(index.locate(station,"X"),o)
    tz,z = read_trace(index.locate(station,"Z"),o)
    if tx.shape!=tz.shape or not np.allclose(tx,tz,rtol=0,atol=(tx[1]-tx[0])*1e-4): raise ValueError("component_time_mismatch:"+station)
    x,z = preprocess(tx,x,o),preprocess(tz,z,o)
    if o["z_component_positive"] == "up": z = -z
    e = polarization(info)
    scalar = e[0]*x+e[1]*z
    orthogonal = -e[1]*x+e[0]*z
    ypath = index.locate(station,"Y",optional=True)
    y = None
    if ypath:
        ty,y = read_trace(ypath,o)
        if tx.shape!=ty.shape or not np.allclose(tx,ty,rtol=0,atol=(tx[1]-tx[0])*1e-4): raise ValueError("component_time_mismatch:"+station)
        y = preprocess(ty,y,o)
    return tx,scalar,orthogonal,y


def measure(t,signal,start,end,o):
    dt = t[1]-t[0]
    if start<t[0]-dt/2 or end>t[-1]+dt/2: raise ValueError("window_outside_record")
    mask = (t>=start)&(t<=end)
    tw,x = t[mask],signal[mask]
    if len(x)<o['min_window_samples']: raise ValueError(f"window_too_short:points={len(x)},required={o['min_window_samples']},dt={dt:.9g}")
    weights = tukey(len(x),alpha=o["taper_fraction"])
    padding = max(len(x),int(len(x)*o["hilbert_pad_factor"]/2))
    env = abs(hilbert(np.pad(x*weights,(padding,padding))))[padding:padding+len(x)]
    di,hi = int(np.argmax(abs(x))),int(np.argmax(env))
    dp,hp = float(abs(x[di])),float(env[hi])
    flags=[]
    if dp<=o["min_reference_amplitude"] or hp<=o["min_reference_amplitude"]: flags.append("zero_or_tiny_signal")
    edge = float(max(abs(x[0]),abs(x[-1]))/max(dp,o["min_reference_amplitude"]))
    if edge>o["edge_amplitude_limit"]: flags.append("window_edge_not_quiet")
    edge_n=max(2,int(np.ceil(o["taper_fraction"]*len(x)/2)))
    if min(di,hi)<edge_n or max(di,hi)>=len(x)-edge_n: flags.append("peak_at_window_edge")
    return dict(direct_peak=dp,hilbert_peak=hp,direct_peak_time_s=float(tw[di]),hilbert_peak_time_s=float(tw[hi]),
                edge_ratio=edge,sample_count=len(x)),flags,(tw,x,env)


def should_measure(row,o):
    if not row['use'] or not row.get('distance_range_selected',1): return False
    # 恰在界面上的混合记录不能当作分离的出射波；近界面记录可供原始比值诊断。
    if abs(row['burial_m']-INTERFACE_DEPTH)<=1e-6: return False
    if row['eligible']: return True
    return (o['output_unreviewed_ratios'] and row['geometry_note']=='too_close_interface')


def has_measured_amplitude(row,o):
    if row['status'].startswith(('error:','excluded:')): return False
    return all(isinstance(row.get(key),(int,float)) and np.isfinite(row[key]) and row[key]>o['min_reference_amplitude']
               for key in ('direct_peak','hilbert_peak'))


def explain_status(status):
    labels={'outside_selected_distance_range':'不在用户指定的沿线距界面范围内','too_close_source':'入射台站距源面小于MIN_SOURCE_DISTANCE_M',
            'too_close_interface':'台站距界面过近或恰在界面',
            'target_phase_window_too_short':'目标波允许窗口的采样点不足，不能跨到其他波扩窗',
            'window_outside_record':'时间窗超出已记录的时间范围',
            'missing_component':'找不到所需分量文件',
            'insufficient_sampling':'每个主频周期的采样点不足',
            'auto_expanded_short_window':'原始窗口太短，已作诊断性扩窗',
            'window_edge_not_quiet':'窗边仍有明显振幅，可能截断波包',
            'peak_at_window_edge':'峰值靠近窗边',
            'window_limited_to_target_phase':'为避开其他震相而缩窗',
            'large_orthogonal_component':'正交分量较大，可能含其他波',
            'large_y_component':'Y分量较大',
            'user_disabled':'用户关闭此台站',
            'zero_or_tiny_signal':'振幅为零或过小'}
    found=[text for flag,text in labels.items() if flag in status]
    if 'possible_overlap_' in status:found.append('目标波可能与其他波包重叠')
    return '；'.join(found) or status


def choose_window_center(row,t,signal,f0,o):
    """以理论时间约束波包包络峰搜索；不把其他震相理论到时当作目标波到时。"""
    theory=float(row['predicted_peak_s'])
    start,end=float(row['start_s']),float(row['end_s'])
    metadata=dict(window_center_mode_used='theory',selected_center_s=theory,estimated_peak_s='',center_offset_s=0.,
                  arrival_pick_status='not_requested',arrival_search_start_s='',arrival_search_end_s='')
    if row.get('manual_window',0):
        center=(start+end)/2
        metadata.update(window_center_mode_used='manual',selected_center_s=center,center_offset_s=center-theory,arrival_pick_status='manual_csv_priority')
        return start,end,metadata,[]
    if o['window_center_mode']=='theory':return start,end,metadata,[]
    width=row.get('search_half_width_s') or (o['arrival_search_half_width_s'] if o['arrival_search_half_width_s'] is not None else 1/f0)
    lo,hi=theory-width,theory+width
    for other in json.loads(row['other_arrivals_s']).values():
        if other<theory-1e-10:lo=max(lo,(theory+other)/2)
        elif other>theory+1e-10:hi=min(hi,(theory+other)/2)
    lo,hi=max(lo,float(t[0])),min(hi,float(t[-1]))
    metadata.update(arrival_search_start_s=lo,arrival_search_end_s=hi)
    indices=np.flatnonzero((t>=lo)&(t<=hi))
    failure=None
    if theory<t[0] or theory>t[-1]:failure='theory_outside_record'
    elif len(indices)<3:failure='search_too_few_samples'
    elif max(abs(signal[indices]))<=o['min_reference_amplitude']:failure='search_zero_signal'
    if failure:
        metadata.update(window_center_mode_used='theory_fallback',arrival_pick_status=failure)
        return start,end,metadata,['arrival_pick_failed_'+failure]
    # 使用整条预处理后的有符号记录的解析信号，避免把短搜索窗截断造成假包络峰。
    padding=len(signal)
    envelope=abs(hilbert(np.pad(signal,(padding,padding))))[padding:padding+len(signal)]
    local=int(np.argmax(envelope[indices]));picked=float(t[indices[local]])
    flags=[]
    if local in (0,len(indices)-1):flags.append('arrival_pick_at_search_edge')
    metadata.update(window_center_mode_used='waveform',estimated_peak_s=picked,selected_center_s=picked,
                    center_offset_s=picked-theory,arrival_pick_status='search_edge' if flags else 'candidate_envelope_peak')
    return start+(picked-theory),end+(picked-theory),metadata,flags


def constrain_target_window(row,start,end,o):
    """限制幅值提取本身，而不仅是中心搜索；不把别的震相峰算到本台站用途上。"""
    if not o['limit_window_to_target_phase']:return start,end,[]
    lo,hi,coincident=target_phase_gate(float(row['predicted_peak_s']),json.loads(row['other_arrivals_s']))
    if coincident:raise ValueError('target_phase_not_separable:coincident_with_'+','.join(coincident))
    clipped_start,clipped_end=max(start,lo),min(end,hi)
    if clipped_end<=clipped_start:raise ValueError('target_phase_window_empty')
    clipped=clipped_start>start+1e-10 or clipped_end<end-1e-10
    if clipped and row.get('manual_window',0):raise ValueError('manual_window_crosses_target_phase_gate')
    return clipped_start,clipped_end,['window_limited_to_target_phase'] if clipped else []


def process_row(row,model,index,o):
    result = dict(row,status="",usable=0,direct_peak="",hilbert_peak="")
    if not should_measure(row,o):
        result["status"] = "excluded:"+(row["geometry_note"] or "user_disabled")
        return result,None
    try:
        if row['station'] not in model.groups.get(row['family']+'_'+row['phase'],{}):
            raise ValueError('station_not_assigned_to_target_phase')
        t,s,orth,y = project_trace(index,row["station"],model.phases[row["phase"]],o)
        dt=t[1]-t[0]
        start,end,center_info,pick_flags=choose_window_center(row,t,s,model.config['frequency_hz'],o)
        result.update(center_info)
        result.update(requested_start_s=start,requested_end_s=end,sample_dt_s=float(dt),
                      samples_per_period=float(1/(dt*model.config['frequency_hz'])))
        start,end,gate_flags=constrain_target_window(row,start,end,o)
        auto_flags=list(pick_flags)+gate_flags
        if ((t>=start)&(t<=end)).sum()<o['min_window_samples'] and o['auto_expand_short_windows']:
            center=float(center_info['selected_center_s'])
            half=o['min_window_samples']*dt/2+dt*1e-4
            if not row.get('manual_window',0):
                expanded_start,expanded_end,_=constrain_target_window(row,center-half,center+half,o)
                start,end=min(start,expanded_start),max(end,expanded_end)
                auto_flags.append('auto_expanded_short_window')
        result.update(start_s=start,end_s=end)
        if o['limit_window_to_target_phase'] and ((t>=start)&(t<=end)).sum()<o['min_window_samples']:
            raise ValueError('target_phase_window_too_short:cannot_expand_across_other_phase')
        metrics,flags,window = measure(t,s,start,end,o)
        flags.extend(auto_flags)
        if result['samples_per_period']<o['min_samples_per_period']: flags.append('insufficient_sampling')
        if not row['eligible']: flags.extend(row['geometry_note'].split(';'))
        result.update(metrics)
        mask = (t>=start)&(t<=end)
        result["orthogonal_peak_ratio"] = float(max(abs(orth[mask]))/max(metrics["direct_peak"],o["min_reference_amplitude"]))
        result["y_peak_ratio"] = "" if y is None else float(max(abs(y[mask]))/max(metrics["direct_peak"],o["min_reference_amplitude"]))
        if result["orthogonal_peak_ratio"]>.25: flags.append("large_orthogonal_component")
        if y is not None and result["y_peak_ratio"]>.1: flags.append("large_y_component")
        for phase,other in json.loads(row["other_arrivals_s"]).items():
            default_half = o["window_half_periods"]/model.config["frequency_hz"]
            before=phase_setting(o,phase,'before_s',o['window_before_s'] or default_half)
            after=phase_setting(o,phase,'after_s',o['window_after_s'] or default_half)
            if start <= other+after and end >= other-before: flags.append("possible_overlap_"+phase)
        hard = any(f in flags for f in ('zero_or_tiny_signal','insufficient_sampling','auto_expanded_short_window'))
        result["usable"] = int(row['eligible'] and not hard and (not flags or row["reviewed"]))
        result["status"] = ("reviewed:" if flags and row["reviewed"] and not hard else "review:")+";".join(flags) if flags else "ok"
        return result,(t,s,window)
    except (ValueError,FileNotFoundError,OSError) as exc:
        result["status"] = "error:"+str(exc)
        return result,None


AMPLITUDE_FIELDS = PLAN_FIELDS+["status","usable","direct_peak","hilbert_peak","direct_peak_time_s",
                              "hilbert_peak_time_s","edge_ratio","sample_count","orthogonal_peak_ratio","y_peak_ratio",
                              "requested_start_s","requested_end_s","sample_dt_s","samples_per_period"]
AMPLITUDE_FIELDS += ['window_center_mode_used','selected_center_s','estimated_peak_s','center_offset_s',
                     'arrival_pick_status','arrival_search_start_s','arrival_search_end_s']
COEFFICIENT_FIELDS = ["family","phase","station","reference_station","x_m","burial_m","distance_to_interface_m",
                      "direct_ratio_abs","hilbert_ratio_abs","theory_abs","direct_relative_error_pct","hilbert_relative_error_pct",
                      "accepted","result_quality","status","reference_status"]
SUMMARY_FIELDS = ["family","phase","wave_status","reference_station","candidate_count","measured_count","accepted_count",
                  "direct_median_abs","direct_q25","direct_q75","hilbert_median_abs","hilbert_q25","hilbert_q75",
                  "theory_abs","direct_relative_error_pct","hilbert_relative_error_pct",
                  "raw_ratio_count","raw_direct_median_abs","raw_hilbert_median_abs","reference_quality","status"]


def relative_error(value,theory):
    return "" if theory < 1e-12 or value == "" else 100*(value-theory)/theory


def coefficients_from_rows(rows,model,o):
    pairs,summaries,references = [],[],{}
    families = selected_families(model,o)
    for family in families:
        incident = [r for r in rows if r["family"]==family and r["phase"].startswith("I_")]
        specified = o["reference_stations"].get(family)
        if specified:
            chosen = [r for r in incident if r["station"]==specified]
            if not chosen: raise ValueError("指定参考台站不属于选中入射射线："+specified)
            reference = chosen[0] if chosen[0]["usable"] else None
            if reference is None and o['output_unreviewed_ratios'] and chosen[0]['eligible'] and chosen[0]['use'] and has_measured_amplitude(chosen[0],o):
                reference=chosen[0]
        else:
            available = [r for r in incident if r["usable"]]
            reference = min(available,key=lambda r:abs(r["distance_to_interface_m"]-o["reference_target_distance_m"])) if available else None
            if reference is None and o['output_unreviewed_ratios']:
                measured=[r for r in incident if r['eligible'] and r['use'] and has_measured_amplitude(r,o)]
                if measured: reference=min(measured,key=lambda r:abs(r['distance_to_interface_m']-o['reference_target_distance_m']))
        references[family] = reference
        for phase,info in model.phases.items():
            if phase.startswith("I_"): continue
            subset = [r for r in rows if r["family"]==family and r["phase"]==phase and should_measure(r,o)]
            phase_pairs = []
            for row in subset:
                measured = has_measured_amplitude(row,o)
                d = row["direct_peak"]/reference["direct_peak"] if measured and reference else ""
                h = row["hilbert_peak"]/reference["hilbert_peak"] if measured and reference else ""
                pair = {k:row[k] for k in ("family","phase","station","x_m","burial_m","distance_to_interface_m")}
                pair.update(reference_station=reference["station"] if reference else (specified or ""),direct_ratio_abs=d,hilbert_ratio_abs=h,
                            theory_abs=info["theory_abs"],direct_relative_error_pct=relative_error(d,info["theory_abs"]),
                            hilbert_relative_error_pct=relative_error(h,info["theory_abs"]),accepted=int(bool(reference) and reference['usable'] and row["usable"]),
                            result_quality=('prepared_only' if row['status']=='prepared_only' else 'accepted' if reference and reference['usable'] and row['usable'] else ('unreviewed_raw_ratio' if d!='' else 'unavailable')),
                            status=row["status"],reference_status=reference["status"] if reference else "no_usable_incident_reference")
                phase_pairs.append(pair)
            pairs.extend(phase_pairs)
            accepted = [r for r in phase_pairs if r["accepted"]]
            raw = [r for r in phase_pairs if r['direct_ratio_abs']!='' and r['hilbert_ratio_abs']!='']
            summary = dict(family=family,phase=phase,wave_status=info["status"],reference_station=reference["station"] if reference else "",
                           candidate_count=len(subset),measured_count=sum(r["direct_peak"]!="" for r in subset),accepted_count=len(accepted),
                           theory_abs=info["theory_abs"],raw_ratio_count=len(raw),
                           reference_quality=('accepted' if reference and reference['usable'] else ('unreviewed' if reference else 'unavailable')),
                           status="ok" if accepted else ("unreviewed_raw_ratios_only" if raw else "no_accepted_measurement"))
            if info["status"] != "propagating": summary["status"] = "not_applicable_"+info["status"]
            elif not reference: summary["status"] = "no_usable_incident_reference"
            elif not subset: summary["status"] = "no_eligible_station"
            for kind in ("direct","hilbert"):
                if raw: summary['raw_'+kind+'_median_abs']=float(np.median([r[kind+'_ratio_abs'] for r in raw]))
                if accepted:
                    q = np.quantile([r[kind+"_ratio_abs"] for r in accepted],[.25,.5,.75])
                    summary.update({kind+"_q25":float(q[0]),kind+"_median_abs":float(q[1]),kind+"_q75":float(q[2]),
                                    kind+"_relative_error_pct":relative_error(float(q[1]),info["theory_abs"])})
            summaries.append(summary)
    return pairs,summaries,references


def configure_plots():
    plt.rcParams.update({"font.sans-serif":["Microsoft YaHei","SimHei","DejaVu Sans"],"axes.unicode_minus":False,
                         "pdf.fonttype":42,"font.size":10})


def save_plot(fig,path,o):
    for extension in ('.png','.pdf'):
        target=path.with_suffix(extension)
        options=dict(bbox_inches='tight',facecolor='white')
        if extension=='.png': options['dpi']=o['plot_dpi']
        try:
            fig.savefig(target,**options)
        except OSError as exc:
            if exc.errno not in (13,22): raise
            # Windows图片/PDF查看器可能占用旧文件；不能让绘图失败丢失整次结果。
            target=path.with_name(path.name+'_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')).with_suffix(extension)
            fig.savefig(target,**options)
            print('原图被占用或不能覆盖，本次图另存：',target)
        GENERATED_PLOT_PATHS.append(str(target))
    plt.close(fig)


def plot_plan(plan,model,out,o):
    configure_plots()
    families = sorted(set(r["family"] for r in plan))
    fig,axes = plt.subplots(1,len(families),figsize=(7*len(families),6.7),squeeze=False,constrained_layout=True)
    colors={"I_P":"#b34232","I_SV":"#b34232","R_P":"#ba7d19","R_SV":"#8a57a8","T_P":"#007d99","T_SV":"#218457"}
    for ax,family in zip(axes[0],families):
        ax.axhspan(0,500,color="#eaf4fa");ax.axhspan(500,1000,color="#f6f0e7")
        ax.axhline(500,color="#788793",lw=1)
        end = model.source+model.aperture*model.tangent
        ax.plot([model.source[0],end[0]],[model.source[1],end[1]],color="#323d49",lw=2,label="源面投影")
        for phase in model.phases:
            selected = [r for r in plan if r["family"]==family and r["phase"]==phase and should_measure(r,o)]
            if selected: ax.scatter([r["x_m"] for r in selected],[r["burial_m"] for r in selected],s=10,color=colors[phase],label=phase)
        ax.set(xlim=(0,1000),ylim=(1000,0),xlabel="X / m",ylabel="埋深 / m（向下为正）",
               title=("源面中间对应射线" if family=="central" else "自动完整测线（含波束覆盖段）" if family=="auto" else "整个源面首次触界射线")+" | y=500 m")
        ax.set_aspect("equal");ax.legend(loc="upper center",bbox_to_anchor=(.5,-.13),ncol=3,frameon=False)
    fig.suptitle("后处理候选台站：仅理论射线，不使用网格台站"+("\n"+o["plot_note"] if o["plot_note"] else ""),fontsize=15)
    save_plot(fig,out/"selected_stations",o)


def plot_waveform(row,payload,out,o):
    t,signal,(tw,x,env) = payload
    fig,axes = plt.subplots(2,1,figsize=(10,6.6),constrained_layout=True)
    fig.suptitle(f'{row["family"]} / {row["station"]}\n本条处理用途：仅拾取{phase_label(row["phase"])}（{row["phase"]}）'+("\n"+o["plot_note"] if o["plot_note"] else ""),fontsize=12)
    axes[0].plot(t,signal,color="#3d5a72",lw=.85)
    axes[0].axvspan(row["start_s"],row["end_s"],color="#f3d08d",alpha=.5,label="提取时间窗")
    axes[0].axvline(row["predicted_peak_s"],color="#b04535",ls="--",label=phase_label(row['phase'])+" 目标理论到时")
    if row.get('window_center_mode_used') in ('waveform','manual'):
        axes[0].axvline(row['selected_center_s'],color='#218457',ls='-.',label='实际选窗中心')
    if row.get('estimated_peak_s','')!='':
        axes[0].axvspan(row['arrival_search_start_s'],row['arrival_search_end_s'],color='#218457',alpha=.10,label='候选到时搜索区')
    for phase,value in json.loads(row["other_arrivals_s"]).items():
        if t[0]<=value<=t[-1]: axes[0].axvline(value,color="#999",ls=":",label=phase+" 干扰参考（不拾取）")
    axes[0].legend(fontsize=8,ncol=3,loc="upper right")
    axes[1].plot(tw,x,color="#3d5a72",label="投影后的有符号波形")
    axes[1].plot(tw,env,color="#cf6a24",lw=1.5,label="Hilbert包络")
    axes[1].axvline(row['predicted_peak_s'],color='#b04535',ls='--',lw=.8,label='理论到时')
    if row.get('selected_center_s','')!='':axes[1].axvline(row['selected_center_s'],color='#218457',ls='-.',lw=.8,label='选窗中心')
    axes[1].axvline(row["direct_peak_time_s"],color="#3d5a72",ls=":")
    axes[1].axvline(row["hilbert_peak_time_s"],color="#cf6a24",ls=":")
    axes[1].set_title(f'直接峰值={row["direct_peak"]:.5g}；包络峰值={row["hilbert_peak"]:.5g}',fontsize=10)
    axes[1].legend(fontsize=9)
    unit={"acceleration":"m/s²","velocity":"m/s","displacement":"m"}[o["trace_quantity"]]
    for ax in axes: ax.set(xlabel="记录时间 / s",ylabel="极化投影幅值 / "+unit);ax.grid(alpha=.15)
    fig.text(.01,-.035,row["status"],fontsize=8,color="#8b4932",wrap=True)
    save_plot(fig,out/(row["family"]+"_"+row["phase"]+"_"+row["station"]),o)


def plot_coefficients(pairs,summaries,out,o):
    configure_plots()
    for family in sorted(set(s["family"] for s in summaries)):
        ss=[s for s in summaries if s["family"]==family]
        fig,axes=plt.subplots(len(ss),1,figsize=(10,3.3*len(ss)),squeeze=False,constrained_layout=True)
        for ax,s in zip(axes[:,0],ss):
            if s['wave_status'] != 'propagating':
                ax.axis('off')
                meaning='倏逝波' if s['wave_status']=='evanescent' else '临界掠射波'
                ax.text(.5,.5,f'{s["phase"]}：{meaning}\n不按远场传播波峰值比计算\n理论复系数模 = {s["theory_abs"]:.6g}',
                        ha='center',va='center',transform=ax.transAxes,color='#8d523b',fontsize=13)
                continue
            valid=[r for r in pairs if r["family"]==family and r["phase"]==s["phase"] and r["accepted"]]
            raw=[r for r in pairs if r['family']==family and r['phase']==s['phase'] and not r['accepted'] and r['direct_ratio_abs']!='']
            ax.axhline(s["theory_abs"],color="#424e59",ls="--",label="理论系数模")
            if valid:
                ax.scatter([r["distance_to_interface_m"] for r in valid],[r["direct_ratio_abs"] for r in valid],s=18,color="#237a99",label="直接峰值比")
                ax.scatter([r["distance_to_interface_m"] for r in valid],[r["hilbert_ratio_abs"] for r in valid],s=20,marker="^",color="#d27b22",label="Hilbert包络峰值比")
            if raw:
                ax.scatter([r['distance_to_interface_m'] for r in raw],[r['direct_ratio_abs'] for r in raw],s=20,facecolors='none',edgecolors='#237a99',label='原始直接比（待核对）')
                ax.scatter([r['distance_to_interface_m'] for r in raw],[r['hilbert_ratio_abs'] for r in raw],s=23,marker='^',facecolors='none',edgecolors='#d27b22',label='原始包络比（待核对）')
            if not valid:
                ax.text(.5,.92,"仅原始振幅比，尚无合格系数" if raw else "无可计算数据："+s['status'],ha='center',va='top',transform=ax.transAxes,color='#9b493a',bbox=dict(fc='white',ec='none',alpha=.9))
            ax.set(title=s["phase"],xlabel="沿射线距界面交点 / m",ylabel="反/透射系数模")
            ax.legend(loc="best",fontsize=8);ax.grid(alpha=.15)
        fig.suptitle(f"{family} | 实心：合格结果；空心：原始振幅比（待核对）"+("\n"+o["plot_note"] if o["plot_note"] else ""),fontsize=14)
        save_plot(fig,out/("coefficients_"+family),o)


def resolve_source_clock(o):
    mode=o['time_reference_mode']
    source_start=None
    if mode=='record_time': peak=float(o['source_peak_time_s'])
    elif mode=='elapsed_from_start':
        source_start=float(o['source_start_time_s'])
        peak=float(o['source_peak_time_s'])-source_start
    else:
        if not o['source_time_function_file']: raise ValueError('stf_file模式必须填写SOURCE_TIME_FUNCTION_FILE')
        a=np.loadtxt(resolve(o['source_time_function_file']),ndmin=2)
        if a.shape[1]!=2 or len(a)<3 or not np.isfinite(a).all() or np.any(np.diff(a[:,0])<=0) or not np.max(abs(a[:,1]))>0:
            raise ValueError('源时程必须是两列递增时间、有限非零幅值')
        source_start=float(a[0,0]);peak=float(a[np.argmax(abs(a[:,1])),0])+o['stf_to_record_offset_s']
    if not np.isfinite(peak):raise ValueError('源中心时间不是有限数')
    return dict(mode=mode,source_start_in_original_clock_s=source_start,effective_source_peak_record_s=peak,
                formula='arrival_record_time = effective_source_peak_record_s + travel_time; do not subtract negative record start twice')


def export_pairs(rows,model,o,out):
    """每个入射台站独立作分母；支持全组合及旧几何匹配，不按结果接近理论程度筛选。"""
    families=sorted(set(r['family'] for r in rows))
    rows=[r for r in rows if r.get('distance_range_selected',1) and r.get('use',1)]
    pairs=[]
    outputs=[ph for ph in model.phases if not ph.startswith('I_')]
    for incident in (r for r in rows if r['phase'].startswith('I_')):
        family=incident['family'];di=incident['distance_to_interface_m']
        for phase in outputs:
            info=model.phases[phase]
            target=di*info['speed']/model.incident_speed if o['pairing_rule']=='equal_travel_time' else di
            candidates=[r for r in rows if r['family']==family and r['phase']==phase]
            status='matched';match=None
            if info['status']!='propagating':status='not_applicable_'+info['status']
            elif not candidates: status='no_station_on_ray'
            else:
                lo=min(r['distance_to_interface_m'] for r in candidates);hi=max(r['distance_to_interface_m'] for r in candidates)
                if target<lo-1e-6 or target>hi+1e-6:status='target_outside_ray'
                else:
                    nearest=min(candidates,key=lambda r:abs(r['distance_to_interface_m']-target))
                    if abs(nearest['distance_to_interface_m']-target)>o['pairing_tolerance_m']:status='no_match_within_tolerance'
                    else:match=nearest
            if o['pairing_rule']=='all_pairs':
                matches=candidates if info['status']=='propagating' else []
                status='all_pairs' if matches else ('no_station_in_selected_range' if info['status']=='propagating' else 'not_applicable_'+info['status'])
            else:matches=[match]
            for match in matches:
                values=bool(match) and has_measured_amplitude(incident,o) and has_measured_amplitude(match,o)
                accepted=int(bool(values) and incident['usable'] and match['usable'])
                allowed=values and (o['output_unreviewed_ratios'] or accepted)
                d=match['direct_peak']/incident['direct_peak'] if allowed else ''
                h=match['hilbert_peak']/incident['hilbert_peak'] if allowed else ''
                quality='accepted' if accepted else ('unreviewed_raw_ratio' if allowed else 'unavailable')
                reasons=[]
                if not match:reasons.append('未找到几何配对台站：'+status)
                if not has_measured_amplitude(incident,o):reasons.append('入射分母不可用：'+explain_status(incident['status']))
                if match and not has_measured_amplitude(match,o):reasons.append('出射分子不可用：'+explain_status(match['status']))
                if values and not accepted:
                    if not incident['usable']:reasons.append('入射波待核对：'+explain_status(incident['status']))
                    if not match['usable']:reasons.append('出射波待核对：'+explain_status(match['status']))
                    if not o['output_unreviewed_ratios']:reasons.append('配置禁止输出待核对比值')
                if incident['status']=='prepared_only':
                    quality='prepared_only';reasons=['仅生成台站与时间窗计划，尚未读取波形计算系数']
                quality_zh={'prepared_only':'仅准备，未计算','accepted':'通过检查','unreviewed_raw_ratio':'已有原始比值，待核对','unavailable':'无法计算比值'}[quality]
                p=dict(family=family,phase=phase,incident_station=incident['station'],outgoing_station=match['station'] if match else '',
                       pairing_rule=o['pairing_rule'],incident_distance_m=di,target_distance_m=target if o['pairing_rule']!='all_pairs' else '',
                       actual_distance_m=match['distance_to_interface_m'] if match else '',
                       distance_mismatch_m=match['distance_to_interface_m']-target if match and o['pairing_rule']!='all_pairs' else '',
                       incident_time_to_interface_s=di/model.incident_speed,
                       outgoing_time_from_interface_s=match['distance_to_interface_m']/info['speed'] if match else '',
                       incident_x_m=incident['x_m'],incident_burial_m=incident['burial_m'],
                       outgoing_x_m=match['x_m'] if match else '',outgoing_burial_m=match['burial_m'] if match else '',
                       incident_direct_peak=incident['direct_peak'],incident_hilbert_peak=incident['hilbert_peak'],
                       outgoing_direct_peak=match['direct_peak'] if match else '',outgoing_hilbert_peak=match['hilbert_peak'] if match else '',
                       direct_ratio_abs=d,hilbert_ratio_abs=h,theory_abs=info['theory_abs'],accepted=accepted,
                       result_quality=quality,result_quality_zh=quality_zh,quality_reason='；'.join(reasons) or '通过当前程序检查',
                       match_status=status,incident_status=incident['status'],outgoing_status=match['status'] if match else '')
                pairs.append(p)
    write_csv(out/'station_pairs.csv',pairs,list(pairs[0]) if pairs else ['family','phase','incident_station'])
    stats=[]
    for family in families:
        for phase in outputs:
            all_pairs=[p for p in pairs if p['family']==family and p['phase']==phase]
            for scope in ('raw','accepted'):
                valid=[p for p in all_pairs if p['direct_ratio_abs']!='' and (scope=='raw' or p['accepted'])]
                for method in ('direct','hilbert'):
                    v=np.array([p[method+'_ratio_abs'] for p in valid],float)
                    s=dict(family=family,phase=phase,scope=scope,method=method,pair_count=len(v),
                           unique_outgoing_station_count=len(set(p['outgoing_station'] for p in valid)),
                           unique_incident_station_count=len(set(p['incident_station'] for p in valid)),
                           min='',max='',median='',mean='',variance_population='',variance_sample='',theory_abs=model.phases[phase]['theory_abs'])
                    if len(v):s.update(min=float(v.min()),max=float(v.max()),median=float(np.median(v)),mean=float(v.mean()),
                                       variance_population=float(v.var(ddof=0)),variance_sample=float(v.var(ddof=1)) if len(v)>1 else '')
                    stats.append(s)
    write_csv(out/'paired_coefficient_statistics.csv',stats,list(stats[0]) if stats else ['family'])
    names={'family':'射线组','phase':'波型','scope':'统计范围','method':'方法','pair_count':'配对数量',
           'unique_incident_station_count':'不同入射台站数','unique_outgoing_station_count':'不同出射台站数','min':'最小值','max':'最大值','median':'中位数','mean':'平均值',
           'variance_population':'总体方差_N','variance_sample':'样本方差_N减1','theory_abs':'理论模'}
    translated=[{names[k]:({'raw':'原始比值_含待核对','accepted':'仅合格','direct':'直接','hilbert':'Hilbert'}.get(v,v) if isinstance(v,str) else v) for k,v in s.items()} for s in stats]
    write_csv(out/'配对系数统计_最大最小中位平均方差.csv',translated,list(names.values()))
    return dict(total_pairs=len(pairs),measured_pairs=sum(p['direct_ratio_abs']!='' for p in pairs),accepted_pairs=sum(p['accepted'] for p in pairs))


def validate_case_metadata(model,plan,data,out):
    """若算例保留DATA输入，核对真实台站坐标和源主频，拒绝同名异坐标的混用。"""
    folders=[Path(data).parent/'data',Path(data).parent/'DATA']
    folder=next((p for p in folders if p.is_dir()),folders[0])
    errors=[]
    report={'status':'not_verified','input_folder':str(folder),'checks':{},'errors':errors}
    station_file=folder/'STATIONS'
    if station_file.exists():
        actual={}
        for line in station_file.read_text(encoding='utf-8-sig').splitlines():
            parts=line.split('#',1)[0].split()
            if not parts:continue
            if len(parts)!=6:raise ValueError('实际STATIONS格式应为6列：'+str(station_file))
            key=(parts[0],parts[1])
            if key in actual:raise ValueError('实际STATIONS台站/网络重复：'+str(key))
            actual[key]=np.array([float(v.replace('D','e').replace('d','e')) for v in parts[2:]])
        requested={r['station'] for r in plan if r['use']}
        missing=[];mismatch=[]
        for name in sorted(requested):
            st=model.stations[name];key=(name,st['network'])
            if key not in actual:missing.append(name);continue
            expected=[st['y_m'],st['x_m'],0.,st['burial_m']]
            if not np.allclose(actual[key],expected,rtol=0,atol=2e-6):mismatch.append(name)
        report['checks']['stations']={'requested':len(requested),'missing_count':len(missing),'coordinate_mismatch_count':len(mismatch),
                                      'missing_examples':missing[:5],'mismatch_examples':mismatch[:5]}
        if missing or mismatch:errors.append(f'台站元数据与实际STATIONS不一致：缺少{len(missing)}个，同名坐标不同{len(mismatch)}个。')
    source_file=folder/'FORCESOLUTION'
    if source_file.exists():
        content=source_file.read_text(encoding='utf-8-sig')
        values=re.findall(r'^\s*(?:hdurorf0|f0)\s*:+\s*([^\s]+)',content,flags=re.MULTILINE|re.IGNORECASE)
        if values:
            frequencies=np.array([float(v.replace('D','e').replace('d','e')) for v in values])
            report['checks']['source_frequency']={'metadata_hz':model.config['frequency_hz'],
                'actual_min_hz':float(frequencies.min()),'actual_max_hz':float(frequencies.max())}
            if not np.allclose(frequencies,model.config['frequency_hz'],rtol=1e-7,atol=1e-9):
                errors.append(f"源主频不一致：元数据{model.config['frequency_hz']:g}Hz，实际FORCESOLUTION为{frequencies.min():g}–{frequencies.max():g}Hz。")
    report['status']='mismatch' if errors else ('checked_available_inputs' if report['checks'] else 'not_verified')
    (out/'case_consistency.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    if errors:raise ValueError('算例配置混用，已停止计算。'+' '.join(errors)+' 请将SETUP_DIRECTORY指向本次模拟对应的生成器输出，详见case_consistency.json。')
    if not report['checks']:print('未找到可核对的实际STATIONS/FORCESOLUTION；请自行确认元数据与波形属于同次模拟。')
    return report


def input_diagnostics(model,index,o):
    result={'source_frequency_hz':model.config['frequency_hz'],'warnings':[]}
    paths=next((paths for (name,component),paths in index.files.items() if name in model.stations and component==o['component_codes']['X']),None)
    if paths:
        t,_=read_trace(paths[0],o);dt=float(t[1]-t[0]);f0=model.config['frequency_hz']
        result.update(example_file=str(paths[0]),samples=len(t),start_s=float(t[0]),end_s=float(t[-1]),dt_s=dt,
                      sampling_frequency_hz=1/dt,nyquist_hz=1/(2*dt),samples_per_period=1/(dt*f0))
        if result['samples_per_period']<o['min_samples_per_period']:
            result['warnings'].append(f'主频{f0:g}Hz，记录间隔{dt:g}s，每周期仅{1/(dt*f0):.3f}点。扩窗/插值不能恢复未保存的峰值，需更密输出。')
        width=sum(o[k] if o[k] is not None else o['window_half_periods']/f0 for k in ('window_before_s','window_after_s'))
        if width/dt<o['min_window_samples']:result['warnings'].append(f'请求窗宽{width:g}s，约{width/dt:.2f}个采样间隔，不足{ o["min_window_samples"]}点。自动扩窗仅输出诊断比值。')
    else:result['warnings'].append('未找到理论台站X分量文件，请核对目录、NETWORK和台站名。')
    par=resolve(o['data_directory']).parent/'data'/'Par_file'
    if par.exists():
        entries={}
        for line in par.read_text(encoding='utf-8-sig').splitlines():
            line=line.split('#')[0]
            if '=' in line:
                k,v=line.split('=',1);entries[k.strip()]=v.strip()
        result['solver_parameters']={k:entries.get(k) for k in ('DT','NSTEP','NTSTEP_BETWEEN_OUTPUT_SAMPLE','ATTENUATION')}
        if result.get('samples') and all(entries.get(k) for k in ('DT','NSTEP','NTSTEP_BETWEEN_OUTPUT_SAMPLE')):
            total=int(entries['NSTEP']);stride=int(entries['NTSTEP_BETWEEN_OUTPUT_SAMPLE'])
            expected=(total+stride-1)//stride
            result['expected_samples_from_par_file']=expected
            if result['samples']<expected:result['warnings'].append(f'样本道仅{result["samples"]}点，Par_file预期约{expected}点；请核对是否模拟未完成或数据未复制完整。')
        if entries.get('ATTENUATION','').lower()=='.true.':result['warnings'].append('Par_file开启ATTENUATION；实际是否衰减需结合材料Q值，不能直接假定为论文无衰减介质。')
    return result


def plot_pair_results(out,o):
    records=read_csv(out/'station_pairs.csv')
    for family in sorted(set(r['family'] for r in records)):
        phases=sorted(set(r['phase'] for r in records if r['family']==family))
        fig,axes=plt.subplots(len(phases),1,figsize=(10,3.1*len(phases)),squeeze=False,constrained_layout=True)
        for ax,phase in zip(axes[:,0],phases):
            rr=[r for r in records if r['family']==family and r['phase']==phase]
            valid=[r for r in rr if r['direct_ratio_abs']!='']
            if rr:ax.axhline(float(rr[0]['theory_abs']),color='#555',ls='--',label='理论模')
            for method,color,marker,label in [('direct','#237a99','o','直接'),('hilbert','#d27b22','^','Hilbert')]:
                for accepted in ('0','1'):
                    subset=[r for r in valid if r['accepted']==accepted]
                    if subset:ax.scatter([float(r['incident_distance_m']) for r in subset],[float(r[method+'_ratio_abs']) for r in subset],
                                         s=22,marker=marker,edgecolors=color,facecolors=color if accepted=='1' else 'none',
                                         label=label+('（合格）' if accepted=='1' else '（待核对）'))
            ax.set(title=phase,xlabel='入射台站沿射线距界面 / m',ylabel='配对振幅比')
            ax.grid(alpha=.15);ax.legend(fontsize=8)
            if not valid:ax.text(.5,.5,'无可计算配对，原因见station_pairs.csv',ha='center',transform=ax.transAxes)
        fig.suptitle(f'{family}：每个入射台站使用自己的振幅作分母\n配对规则：{o["pairing_rule"]}；空心点尚未通过检查',fontsize=13)
        save_plot(fig,out/('paired_coefficients_'+family),o)


def hilbert_error_matches(records,threshold,accepted_only=False):
    selected=[]
    for r in records:
        if accepted_only and str(r['accepted'])!='1':continue
        if r.get('hilbert_ratio_abs','')=='':continue
        value,theory=float(r['hilbert_ratio_abs']),float(r['theory_abs'])
        if not np.isfinite([value,theory]).all() or theory<=1e-12:continue
        error=100*abs(value-theory)/theory
        if error<threshold and not np.isclose(error,threshold,rtol=0,atol=1e-10):
            selected.append(dict(r,hilbert_error_pct=error))
    return selected


def best_station_pairs(records,families,phases,o):
    """在已计算组合中逐类取最小误差，不改变波形、窗口、距离筛选或系数。"""
    results=[];method=o['best_pair_method']
    for family in families:
        for phase in phases:
            if phase.startswith('I_'):continue
            eligible=[]
            for row in records:
                if row['family']!=family or row['phase']!=phase:continue
                if o['best_pair_accepted_only'] and str(row['accepted'])!='1':continue
                if row.get(method+'_ratio_abs','')=='':continue
                value,theory=float(row[method+'_ratio_abs']),float(row['theory_abs'])
                if not np.isfinite([value,theory]).all() or theory<=1e-12:continue
                eligible.append(dict(row,best_method=method,best_error_pct=100*abs(value-theory)/theory))
            if eligible:
                winner=min(eligible,key=lambda r:(r['best_error_pct'],r['incident_station'],r['outgoing_station']))
                winner.update(best_status='selected',best_candidate_count=len(eligible))
            else:
                winner=dict(family=family,phase=phase,best_method=method,best_status='无可用合格组合或理论模接近零' if o['best_pair_accepted_only'] else '无可计算组合或理论模接近零',best_candidate_count=0)
            results.append(winner)
    return results


def best_pick_details(best,rows,model,index,o):
    """用实际已选窗口重现最佳组合拾取，校验与计算值一致；不重新搜索峰窗。"""
    details=[];payloads={}
    lookup={(r['family'],r['phase'],r['station']):r for r in rows}
    for pair in best:
        if pair['best_status']!='selected':continue
        key=(pair['family'],pair['phase']);payloads[key]=[]
        for role,phase,name in [('入射分母','I_'+model.config['incident_wave'],pair['incident_station']),
                                ('出射分子',pair['phase'],pair['outgoing_station'])]:
            row=lookup[(pair['family'],phase,name)]
            t,signal,_,_=project_trace(index,name,model.phases[phase],o)
            metrics,_,window=measure(t,signal,row['start_s'],row['end_s'],o)
            for field in ('direct_peak','hilbert_peak','direct_peak_time_s','hilbert_peak_time_s'):
                if not np.isclose(metrics[field],row[field],rtol=1e-10,atol=1e-12):
                    raise ValueError('最佳台站拾取重现与原计算不一致：'+name+' '+field)
            detail=dict(family=pair['family'],coefficient_phase=pair['phase'],role=role,station=name,target_phase=phase,
                        distance_to_interface_m=row['distance_to_interface_m'],predicted_peak_s=row['predicted_peak_s'],
                        selected_center_s=row['selected_center_s'],start_s=row['start_s'],end_s=row['end_s'],
                        **metrics,status=row['status'],usable=row['usable'],
                        envelope_to_direct_peak=metrics['hilbert_peak']/metrics['direct_peak'],
                        direct_ratio_abs=float(pair['direct_ratio_abs']),hilbert_ratio_abs=float(pair['hilbert_ratio_abs']),
                        theory_abs=float(pair['theory_abs']),
                        bandpass_hz=str(o['bandpass_hz']),quantity=o['trace_quantity'],
                        extraction='直接=max|窗内有符号投影波形|；包络=窗内波形经Tukey渐消、双侧补零、Hilbert解析信号取模后的最大值；系数=出射峰值/入射峰值')
            details.append(detail)
            # 展示完整记录；峰值仍严格来自原来的实际拾取窗。
            payloads[key].append((detail,t,signal,window))
    return details,payloads


def plot_best_pick(ax,payload,o):
    detail,t,signal,(tw,x,env)=payload
    direct_index=int(np.argmax(abs(x)));hilbert_index=int(np.argmax(env))
    ax.axvspan(detail['start_s'],detail['end_s'],color='#d7edf3',alpha=.6,label='实际拾取窗')
    ax.plot(t,signal,color='#226b9c',lw=.9,label='全时程投影波形（按当前预处理设置）')
    ax.plot(tw,env,color='#c75b23',lw=1.5,label='拾取窗内Hilbert包络（含渐消）')
    ax.scatter(tw[direct_index],x[direct_index],s=65,marker='s',color='#163d77',zorder=6,label='直接峰：取绝对值')
    ax.scatter(tw[hilbert_index],env[hilbert_index],s=65,marker='o',color='#d72d58',zorder=7,label='包络峰')
    ax.vlines(tw[direct_index],0,x[direct_index],color='#163d77',lw=1.5,zorder=5)
    ax.vlines(tw[hilbert_index],0,env[hilbert_index],color='#d72d58',lw=1.5,zorder=5)
    ax.axhline(0,color='#999',lw=.6,zorder=0)
    ax.axvline(detail['predicted_peak_s'],color='#32804b',ls='--',lw=1,label='理论到时')
    for edge in (detail['start_s'],detail['end_s']):ax.axvline(edge,color='#6c929d',lw=.8)
    units={'acceleration':'m/s²','velocity':'m/s','displacement':'m'}
    unit=units[o['trace_quantity']]
    ax.set(xlabel='记录时间 / s（全时程）',ylabel='幅值 / '+unit,xlim=(t[0],t[-1]))
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)
    ax.set_title(f"{detail['role']}：{detail['station']}（{detail['target_phase']}）\n"
                 f"直接峰 |A|={detail['direct_peak']:.6g} {unit}，t={detail['direct_peak_time_s']:.6f} s\n"
                 f"包络峰 A={detail['hilbert_peak']:.6g} {unit}，t={detail['hilbert_peak_time_s']:.6f} s",fontsize=10)
    ax.legend(fontsize=7,loc='lower left',ncol=2);ax.grid(alpha=.15)
    ax.margins(y=.25)


def plot_best_full_records(root,model,o,best,payloads):
    """宽幅全时程图，每对台站上下排列；保留原始幅值比例，不逐道归一化。"""
    configure_plots()
    heights=[8 if payloads.get((p['family'],p['phase'])) else 1 for p in best]
    fig=plt.figure(figsize=(18,max(4,sum(heights))),layout='constrained')
    grid=fig.add_gridspec(max(1,len(best)),1,height_ratios=heights or [1])
    filt='未带通滤波' if o['bandpass_hz'] is None else f"带通 {o['bandpass_hz']} Hz"
    fig.suptitle(f"最优台站全时程与幅值拾取｜{model.config['incident_side']}层 {model.config['incident_wave']} 入射 {model.config['angle_deg']:g}°\n"
                 +filt+'；'+getattr(model,'comparison_theory',{}).get('description','内置界面方程')
                 +'；显示完整记录，蓝色阴影为实际拾取窗，峰值仅在窗内提取；不逐道归一化。',fontsize=15)
    for i,pair in enumerate(best):
        panel=fig.add_subfigure(grid[i,0])
        payload=payloads.get((pair['family'],pair['phase']),[])
        heading=f"{pair['family']} / {phase_label(pair['phase'])}"
        if not payload:
            ax=panel.subplots();ax.axis('off')
            info=model.phases[pair['phase']]
            reason='倏逝模态，无远场传播测线拾取' if info['status']=='evanescent' else pair['best_status']
            ax.text(.5,.5,heading+'：'+reason,ha='center',va='center',transform=ax.transAxes)
            continue
        method='Hilbert' if pair['best_method']=='hilbert' else '直接峰值'
        quality='通过检查' if str(pair['accepted'])=='1' else '待核对'
        panel.suptitle(heading+f"｜{method}误差最小组合（{quality}）｜直接幅值比={float(pair['direct_ratio_abs']):.8g}；"
                       f"包络幅值比={float(pair['hilbert_ratio_abs']):.8g}；理论模={float(pair['theory_abs']):.8g}",fontsize=12)
        axes=panel.subplots(2,1)
        for ax,item in zip(axes,payload):
            plot_best_pick(ax,item,o)
            ax.legend(fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.19),ncol=3)
            detail=item[0]
            ax.text(.99,.97,f"实际窗：[{detail['start_s']:.6f}, {detail['end_s']:.6f}] s",
                    ha='right',va='top',transform=ax.transAxes,fontsize=9,
                    bbox=dict(facecolor='white',edgecolor='#ddd',alpha=.9))
    target=root/'最优台站全时程.png'
    fig.savefig(target,dpi=o['plot_dpi'],bbox_inches='tight',facecolor='white')
    plt.close(fig);GENERATED_PLOT_PATHS.append(str(target))
    return target


def plot_hilbert_error_combinations(out,model,o,final_path=None,pick_payloads=None):
    """Each point at (incident distance, outgoing distance) identifies one pair.
    The geometry panel retains the source plane projection, interface and rays.
    Selection by agreement is descriptive and is not independent validation.
    """
    configure_plots()
    records=read_csv(out/'station_pairs.csv')
    best={(r['family'],r['phase']):r for r in best_station_pairs(records,selected_families(model,o),model.phases,o)}
    threshold=o['hilbert_error_threshold_pct']
    folder=out/'hilbert_error_combinations'
    if final_path is None:folder.mkdir(exist_ok=True)
    panel_count=len(selected_families(model,o))*sum(not ph.startswith('I_') for ph in model.phases)
    if final_path is not None:
        combined=plt.figure(figsize=(25,8.4*panel_count),layout='constrained')
        layout=combined.add_gridspec(panel_count,3,width_ratios=[1,1.1,1.3])
        case_label=f"{model.config['incident_side']}层 {model.config['incident_wave']} 入射；入射角 {model.config['angle_deg']:g}°（相对界面法线）"
        reference=getattr(model,'comparison_theory',{}).get('description','内置界面方程')
        combined.suptitle('反透射系数：Hilbert误差达标组合'+('（仅准备，尚未计算）' if o['run_mode']=='prepare' else '')+'\n'+case_label+'；比较基准：'+reference,fontsize=18)
    panel=0
    for family in selected_families(model,o):
        for phase,info in model.phases.items():
            if phase.startswith('I_'):continue
            subset=[r for r in records if r['family']==family and r['phase']==phase]
            hits=hilbert_error_matches(subset,threshold,o['hilbert_error_accepted_only'])
            accepted=sum(str(r['accepted'])=='1' for r in hits)
            if final_path is None:
                fig,(geo,pairs_ax)=plt.subplots(1,2,figsize=(16,8),gridspec_kw={'width_ratios':[1,1.15]},layout='constrained')
            else:
                fig=combined;geo=fig.add_subplot(layout[panel,0]);pairs_ax=fig.add_subplot(layout[panel,1])
                wave_layout=layout[panel,2].subgridspec(2,1)
                wave_axes=[fig.add_subplot(wave_layout[i,0]) for i in range(2)]
                payload=(pick_payloads or {}).get((family,phase),[])
                if payload:
                    for ax,item in zip(wave_axes,payload):plot_best_pick(ax,item,o)
                else:
                    for ax in wave_axes:
                        ax.axis('off')
                    reason='透射P为倏逝波，无远场传播测线拾取' if info['status']=='evanescent' else '没有可展示的最佳组合拾取（未计算或无可用组合）'
                    wave_axes[0].text(.5,.5,reason,ha='center',va='center',transform=wave_axes[0].transAxes,wrap=True)
                panel+=1
            heading=(f'{family} · {phase_label(phase)} | Hilbert系数相对理论模误差 < {threshold:g}%\n'
                         f'{len(hits)}组满足数值条件：{accepted}组通过波形检查，{len(hits)-accepted}组待核对'+('\n'+o['plot_note'] if o['plot_note'] else ''))
            if final_path is None:fig.suptitle(heading,fontsize=16)
            geo.axhspan(0,500,color='#eaf4fa');geo.axhspan(500,1000,color='#f6f0e7')
            geo.axhline(500,color='#596d7a',lw=1.6,label='介质界面（埋深500 m）')
            a=model.source;b=a+model.aperture*model.tangent
            geo.plot([a[0],b[0]],[a[1],b[1]],color='#263442',lw=4,label='震源面在y=500剖面的投影')
            colors={'I_P':'#bb4a37','I_SV':'#bb4a37','R_P':'#cb8b19','R_SV':'#8855a0','T_P':'#21859b','T_SV':'#32834c'}
            for ray in model.report['lines']:
                if not ray['name'].startswith(family+'_'):continue
                ph=ray['name'].split('_',1)[1];aa,bb=ray['a_x_depth_m'],ray['b_x_depth_m']
                geo.plot([aa[0],bb[0]],[aa[1],bb[1]],color=colors.get(ph,'#888'),lw=1.6,label=phase_label(ph)+'测线')
            if hits:
                ins={r['incident_station']:r for r in hits};outs={r['outgoing_station']:r for r in hits}
                geo.scatter([float(r['incident_x_m']) for r in ins.values()],[float(r['incident_burial_m']) for r in ins.values()],s=28,marker='o',facecolors='none',edgecolors='#b33d2e',label='入选组合的入射台站',zorder=5)
                geo.scatter([float(r['outgoing_x_m']) for r in outs.values()],[float(r['outgoing_burial_m']) for r in outs.values()],s=30,marker='x',color='#14637c',label='入选组合的出射台站',zorder=6)
                norm=matplotlib.colors.Normalize(0,threshold)
                for flag,marker,label in [('1','o','通过波形检查'),('0','x','待核对原始比值')]:
                    group=[r for r in hits if str(r['accepted'])==flag]
                    if group:pairs_ax.scatter([float(r['incident_distance_m']) for r in group],[float(r['actual_distance_m']) for r in group],c=[r['hilbert_error_pct'] for r in group],cmap='viridis',norm=norm,s=30,marker=marker,label=label)
                fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm,cmap='viridis'),ax=pairs_ax,label='Hilbert系数相对误差 / %',shrink=.8)
                pairs_ax.legend(loc='best')
            else:
                message=('理论系数接近零，相对误差不适用' if info['theory_abs']<=1e-12 else
                         '该波型不是传播波，无对应射线组合' if info['status']!='propagating' else
                         '未找到满足条件的组合（或尚无可计算的波形）')
                pairs_ax.text(.5,.5,message,ha='center',va='center',transform=pairs_ax.transAxes,wrap=True)
            winner=best[(family,phase)]
            best_title='最佳组合：'+winner['best_status']
            if winner['best_status']=='selected':
                for prefix,tag,color in [('incident','I*','#db3372'),('outgoing','O*','#e8ab00')]:
                    xy=(float(winner[prefix+'_x_m']),float(winner[prefix+'_burial_m']))
                    geo.scatter(*xy,s=200,marker='*',color=color,edgecolors='#222',linewidths=.8,zorder=9,label='最佳组合：'+('入射台站' if prefix=='incident' else '出射台站'))
                    geo.annotate(tag,xy,xytext=(-25,-15) if prefix=='incident' else (12,12),textcoords='offset points',fontsize=10,fontweight='bold',arrowprops={'arrowstyle':'-','color':'#333'},zorder=10)
                pairs_ax.scatter(float(winner['incident_distance_m']),float(winner['actual_distance_m']),s=230,marker='*',color='#ef4277',edgecolors='#222',zorder=10,label='最小误差组合（星号）')
                method_label='Hilbert' if winner['best_method']=='hilbert' else '直接峰值'
                quality='通过检查' if str(winner['accepted'])=='1' else '待核对原始比值'
                best_title=(f"{method_label}最小误差组合（{quality}；候选{winner['best_candidate_count']}对）\n"
                            f"I*: {winner['incident_station']}（距界面{float(winner['incident_distance_m']):g} m）\n"
                            f"O*: {winner['outgoing_station']}（距界面{float(winner['actual_distance_m']):g} m）\n"
                            f"数值={float(winner[winner['best_method']+'_ratio_abs']):.8g}，理论={float(winner['theory_abs']):.8g}，误差={winner['best_error_pct']:.6g}%")
                pairs_ax.legend(loc='best',fontsize=8)
            geo.set(xlim=(0,1000),ylim=(1000,0),xlabel='X / m',ylabel='埋深 / m（向下为正）',title='空间位置：y=500 m剖面')
            if final_path is not None:geo.set_title(heading+f'\n理论系数模={info["theory_abs"]:.6g}',fontsize=11)
            geo.set_aspect('equal');geo.legend(loc='best',fontsize=8)
            pairs_ax.set(xlabel='入射台站沿测线到界面的距离 / m',ylabel=phase_label(phase)+'台站沿测线到界面的距离 / m',
                         title=best_title)
            pairs_ax.title.set_fontsize(10)
            for ax in (geo,pairs_ax):ax.grid(alpha=.15)
            filt='未启用带通滤波' if o['bandpass_hz'] is None else f"Butterworth带通 {o['bandpass_hz']} Hz，N={o['filter_order']}，前后向滤波"
            if final_path is None:fig.supxlabel(f'误差 = |Hilbert系数模 − 理论系数模| / 理论系数模 × 100%；理论模={info["theory_abs"]:.6g}\n'
                          +filt+'。按理论符合度筛选不等于独立验证通过。',fontsize=10)
            if final_path is None:save_plot(fig,folder/(family+'_'+phase+'_hilbert_error'),o)
    if final_path is not None:
        combined.supxlabel(f'普通点：Hilbert相对误差严格小于{o["hilbert_error_threshold_pct"]:g}%；星号：所选方法在候选组合中误差最小的一对（不受5%门槛限制）。\n'+filt+'；最小误差组合是事后筛选结果，不代替全部组合统计。',fontsize=10)
        staged=out/'最终结果图_写入中.png'
        combined.savefig(staged,dpi=o['plot_dpi'],bbox_inches='tight',facecolor='white')
        plt.close(combined);os.replace(staged,final_path)
        GENERATED_PLOT_PATHS.append(str(final_path))


def prepare_output_folder(root):
    root=Path(root).resolve();root.mkdir(parents=True,exist_ok=True)
    process=root/'过程文件';process.mkdir(exist_ok=True)
    if process.is_symlink():raise ValueError('过程文件目录不能为符号链接')
    pending=[]
    history=process/('旧版输出_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'))
    for item in list(root.iterdir()):
        if item.name in ('过程文件','反透射系数结果.xlsx','反透射系数结果.png','最优台站全时程.png'):continue
        if item.is_symlink():continue
        generated=(item.is_file() and (item.suffix.lower() in ('.csv','.json','.txt','.png','.pdf') or item.name.startswith('STATIONS'))) or item.name in ('waveform_checks','hilbert_error_combinations')
        if generated:
            target=process/item.name
            if item.resolve().parent!=root or root not in target.resolve().parents:raise ValueError('输出归档路径越界')
            if target.exists():
                history.mkdir(exist_ok=True)
                target.rename(history/item.name)
            try:item.rename(target)
            except PermissionError:
                if not item.is_file():raise
                shutil.copy2(item,target);pending.append(str(item))
    # 归档已停用的旧表和旧图，防止把历史结果误当作本次结果；不删除数据。
    legacy_names={'coefficients_by_station.csv','coefficients_summary.csv','计算结果_先看这个.csv',
                  '全部组合系数_每个入射乘每个出射.csv','台站配对系数_每个入射台站一行.csv',
                  'waveform_checks','hilbert_error_combinations','plotted_waveforms.csv'}
    for item in list(process.iterdir()):
        old_plot=item.suffix.lower() in ('.png','.pdf') and item.name.startswith(('coefficients_','paired_coefficients_','selected_stations'))
        if item.name not in legacy_names and not old_plot:continue
        if item.is_symlink() or item.resolve().parent!=process.resolve():raise ValueError('旧版输出归档路径越界')
        history.mkdir(exist_ok=True)
        target=history/item.name
        if root not in target.resolve().parents:raise ValueError('旧版输出归档目标越界')
        item.rename(target)
    (process/'待归档占用文件.json').write_text(json.dumps(pending,ensure_ascii=False,indent=2),encoding='utf-8')
    if pending:print('以下旧过程文件被占用，已复制进过程文件目录；关闭后重新运行即可归档：',pending)
    return process


def export_result_workbook(out,root,o,report):
    # Runtime stays Python 3.10; workbook is a typed snapshot of this run's numerical results.
    try:
        from openpyxl import Workbook
        from openpyxl.styles import Font, PatternFill, Alignment
        from openpyxl.utils import get_column_letter
    except ImportError as exc:
        raise RuntimeError('导出Excel需要openpyxl，请运行 python -m pip install openpyxl') from exc
    wb=Workbook();ws=wb.active;ws.title='结果概览'
    overview=[['反透射系数计算结果','数值 / 说明'],['运行模式',o['run_mode']],['射线组',o['ray_family']],
              ['配对方式',o['pairing_rule']],['组合总数',report['pairing']['total_pairs']],
              ['可计算组合数',report['pairing']['measured_pairs']],['通过检查组合数',report['pairing']['accepted_pairs']],
              ['Hilbert误差筛选阈值（%）',o['hilbert_error_threshold_pct']],
              ['滤波', '未启用' if o['bandpass_hz'] is None else 'Butterworth '+str(o['bandpass_hz'])+' Hz'],
              ['结果说明','直接/Hilbert两列均为系数模；空白表示不可用，不是零。'],
              ['统计说明','原始比值含待核对组合；仅合格统计另列。方差是组合的描述统计。'],
              ['数据性质','仅生成计划，未计算实际系数' if o['run_mode']=='prepare' else '本次输入波形的计算结果']]
    for row in overview:ws.append(row)
    context=report['theory_context']
    ws.append(['理论入射角（度，与界面法线夹角）',context['angle_deg']])
    ws.append(['理论入射侧 / 波型',context['incident_side']+' / '+context['incident_wave']])
    ws.append(['理论参数文件',context['config_file']])
    comparison=context['comparison_theory']
    ws.append(['对比理论来源',comparison['description']])
    ws.append(['外部理论系数表',comparison['table_file'] or '未使用外部表'])
    ws.append(['外部理论参数文件',comparison['parameter_file'] or '未使用外部表'])
    ws.append(['最佳台站组合','见“最佳台站组合”页；每种出射波独立选一对，可使用不同入射台站。'])
    ws.append(['最佳组合筛选方法',o['best_pair_method']+'；'+('仅通过检查的组合' if o['best_pair_accepted_only'] else '全部有数值的组合（含待核对）')])
    sheet=wb.create_sheet('最佳台站组合')
    best_columns={'family':'射线组','phase':'出射波型','best_method':'排序方法','best_candidate_count':'候选组合数',
                  'incident_station':'入射台站_I','outgoing_station':'出射台站_O',
                  'incident_distance_m':'入射沿线距界面_m','actual_distance_m':'出射沿线距界面_m',
                  'direct_ratio_abs':'直接峰值系数模','hilbert_ratio_abs':'Hilbert系数模','theory_abs':'理论系数模',
                  'best_error_pct':'所选方法相对误差_%','result_quality_zh':'质量状态','quality_reason':'质量原因','best_status':'筛选状态'}
    sheet.append(list(best_columns.values()))
    numeric={'best_candidate_count','incident_distance_m','actual_distance_m','direct_ratio_abs','hilbert_ratio_abs','theory_abs','best_error_pct'}
    for row in report['best_station_pairs']:
        sheet.append([None if row.get(k,'')=='' else float(row[k]) if k in numeric else row[k] for k in best_columns])
    sheet=wb.create_sheet('最佳组合幅值拾取')
    pick_columns={'family':'射线组','coefficient_phase':'对应系数波型','role':'分子或分母','station':'台站名','target_phase':'仅拾取的目标震相',
                  'distance_to_interface_m':'沿测线距界面_m','predicted_peak_s':'理论到时_s','selected_center_s':'选窗中心_s',
                  'start_s':'实际窗起点_s','end_s':'实际窗终点_s','sample_count':'窗内采样点数',
                  'direct_peak':'直接峰值_绝对值','direct_peak_time_s':'直接峰值时刻_s',
                  'hilbert_peak':'Hilbert包络峰值','hilbert_peak_time_s':'包络峰值时刻_s',
                  'envelope_to_direct_peak':'包络峰值除以直接峰值','direct_ratio_abs':'该组合直接系数','hilbert_ratio_abs':'该组合Hilbert系数',
                  'theory_abs':'理论系数模','usable':'通过提取检查','status':'详细检查状态',
                  'bandpass_hz':'带通频带_Hz','quantity':'输入物理量','extraction':'实际拾取算法'}
    sheet.append(list(pick_columns.values()))
    for row in report['best_pick_details']:sheet.append([row.get(k) for k in pick_columns])
    if not report['best_pick_details']:sheet.append(['未计算或没有可用最佳组合，未生成幅值拾取记录'])
    ws.append(['最佳组合幅值拾取','见同名工作表、最优台站全时程.png和汇总图右栏：完整记录中标出分母/分子的实际窗口、峰值、峰时；不做逐道归一化。'])
    sheet=wb.create_sheet('理论解与介质');sheet.append(['项目','数值 / 说明','单位 / 定义'])
    sheet.append(['比较基准',comparison['description'],''])
    sheet.append(['外部理论系数表',comparison['table_file'],'实际用于误差及最优组合排序'])
    sheet.append(['外部理论参数文件',comparison['parameter_file'],'核对入射侧、波型和介质参数'])
    sheet.append(['入射角',context['angle_deg'],context['angle_definition']])
    sheet.append(['入射侧',context['incident_side'],'upper上层；lower下层'])
    sheet.append(['入射波型',context['incident_wave'],'P或SV'])
    for side,label in [('upper','上层'),('lower','下层')]:
        for key,unit in [('type','fluid流体；solid固体'),('vp','m/s'),('vs','m/s'),('rho','kg/m³')]:
            sheet.append([label+' '+key,context[side][key],unit])
    for phase,values in context['coefficients'].items():
        for key,label in [('theory_abs','系数模（与数值幅值比比较）'),('theory_real','复系数实部'),('theory_imag','复系数虚部'),('status','传播状态')]:
            sheet.append([phase+' '+label,values[key],'无量纲' if key!='status' else ''])
    sheet.append(['参数文件',context['config_file'],'本次数值模拟对应的生成器输出'])
    sheet.append(['理论核对文件',context['report_file'],'与重新计算的系数一致性核对'])
    sheet.append(['核对范围',context['verification_limit'],''])
    mappings=[('系数统计','配对系数统计_最大最小中位平均方差.csv'),('全部台站组合','station_pairs.csv'),
              ('台站与时间窗','台站用途与拾取结果.csv'),('距离筛选','distance_selection_summary.csv')]
    for title,filename in mappings:
        source=out/filename
        if not source.exists():continue
        with source.open(encoding='utf-8-sig',newline='') as f:
            reader=csv.reader(f);headers=next(reader,[]);sheet=wb.create_sheet(title);sheet.append(headers)
            part=1
            for row in reader:
                if sheet.max_row>=1048576:
                    part+=1;sheet=wb.create_sheet(title+'_'+str(part));sheet.append(headers)
                typed=[]
                for value in row:
                    if value=='':typed.append(None);continue
                    try:
                        numeric=float(value)
                        typed.append(numeric if math.isfinite(numeric) else value)
                    except ValueError:typed.append(value)
                sheet.append(typed)
    sheet=wb.create_sheet('运行参数');sheet.append(['参数','本次值'])
    for key,value in o.items():sheet.append([key,json.dumps(value,ensure_ascii=False) if isinstance(value,(dict,list)) else str(value)])
    for sheet in wb:
        sheet.freeze_panes='A2';sheet.auto_filter.ref=sheet.dimensions
        sheet.sheet_view.showGridLines=False
        sheet.row_dimensions[1].height=34
        for cell in sheet[1]:
            cell.fill=PatternFill('solid',fgColor='24465B');cell.font=Font(name='微软雅黑',bold=True,color='FFFFFF')
            cell.alignment=Alignment(wrap_text=True,vertical='center')
        for col in range(1,sheet.max_column+1):
            text=[str(sheet.cell(r,col).value or '') for r in range(1,min(sheet.max_row,80)+1)]
            width=min(54,max(16,max((len(v.encode('utf-8'))*.65 for v in text),default=16)))
            sheet.column_dimensions[get_column_letter(col)].width=width
        for row in sheet.iter_rows(min_row=2):
            for cell in row:
                cell.font=Font(name='微软雅黑',size=10)
                if cell.data_type=='f':cell.data_type='s' # preserve source text, not accidental formulas
                if isinstance(cell.value,(int,float)):
                    cell.number_format='0.000000E+00' if cell.value and (abs(cell.value)<1e-6 or abs(cell.value)>=1e9) else '0.##########'
                else:cell.alignment=Alignment(vertical='top',wrap_text=False)
        sheet.sheet_properties.pageSetUpPr.fitToPage=True
    ws.column_dimensions['A'].width=32;ws.column_dimensions['B'].width=100
    staged=out/'最终工作簿_写入中.xlsx';wb.save(staged);wb.close()
    target=root/'反透射系数结果.xlsx';os.replace(staged,target)
    return target


def run(o):
    o=copy.deepcopy(o)
    GENERATED_PLOT_PATHS.clear()
    validate_options(o)
    setup,data,out = (resolve(o[k]) for k in ("setup_directory","data_directory","output_directory"))
    if o['run_mode']=='analyze' and not data.is_dir():raise ValueError('波形目录不存在，请将DATA_DIRECTORY改为本次新台站模拟的sema目录：'+str(data))
    if out in (setup,data): raise ValueError("输出目录不能与元数据/原始波形目录相同。")
    model=load_model(setup)
    context=theory_context(model,setup,o)
    root=out
    out=prepare_output_folder(root)
    (out/'theory_context.json').write_text(json.dumps(context,ensure_ascii=False,indent=2),encoding='utf-8')
    print(f"理论解：{context['incident_side']}层 {context['incident_wave']} 入射，入射角={context['angle_deg']:g}°（相对界面法线）；参数来源：{context['config_file']}")
    print('比较基准：',context['comparison_theory']['table_file'] or context['comparison_theory']['description'])
    clock=resolve_source_clock(o)
    o['source_peak_time_s']=clock['effective_source_peak_record_s']
    (out/'time_reference.json').write_text(json.dumps(clock,ensure_ascii=False,indent=2),encoding='utf-8')
    plan=make_plan(model,o)
    consistency=validate_case_metadata(model,plan,data,out) if o['run_mode']=='analyze' else {'status':'prepare_only'}
    write_csv(out/"window_plan.csv",plan,PLAN_FIELDS)
    write_csv(out/'station_roles.csv',plan,['family','station','role_group','target_phase','station_role','station_role_count',
              'x_m','y_m','burial_m','predicted_peak_s','search_half_width_s','target_gate_start_s','target_gate_end_s','coincident_phases','incident_beam_region','outside_beam_selection_margin'])
    range_summary=[]
    for family in selected_families(model,o):
        for phase in model.phases:
            subset=[r for r in plan if r['family']==family and r['phase']==phase]
            selected_rows=[r for r in subset if r['distance_range_selected'] and r['use']]
            bounds=o['station_distance_ranges_m'].get(phase)
            range_summary.append(dict(family=family,phase=phase,min_distance_m=bounds[0] if bounds else '',
                max_distance_m=bounds[1] if bounds else '',total_stations=len(subset),
                in_range_stations=len(selected_rows),eligible_for_measurement=sum(should_measure(r,o) for r in selected_rows)))
    write_csv(out/'distance_selection_summary.csv',range_summary,list(range_summary[0]))
    selected={r["station"]:model.stations[r["station"]] for r in plan if should_measure(r,o)}
    with (out/"STATIONS_selected").open("w",encoding="ascii") as f:
        for name,s in selected.items(): f.write(f'{name} {s["network"]} 500.000000000 {s["x_m"]:.9f} 0.000000000 {s["burial_m"]:.9f}\n')
    (out/"processing_config.json").write_text(json.dumps(o,ensure_ascii=False,indent=2),encoding="utf-8")
    if o["generate_plots"]: plot_plan(plan,model,out,o)
    index=TraceIndex(data,o)
    diagnostics=input_diagnostics(model,index,o) if o['run_mode']=='analyze' else {'source_frequency_hz':model.config['frequency_hz'],'warnings':[],'status':'prepare_only'}
    (out/'input_diagnostics.json').write_text(json.dumps(diagnostics,ensure_ascii=False,indent=2),encoding='utf-8')
    for warning in diagnostics['warnings']:print('数据诊断：',warning)
    rows=[]
    for i,row in enumerate(plan):
        result,_=process_row(row,model,index,o) if o["run_mode"]=="analyze" else (dict(row,status="prepared_only",usable=0,direct_peak="",hilbert_peak=""),None)
        rows.append(result)
    write_csv(out/"amplitudes.csv",rows,AMPLITUDE_FIELDS)
    write_csv(out/'window_centers.csv',rows,['family','phase','station','station_role','target_phase','predicted_peak_s','window_center_mode_used',
              'estimated_peak_s','selected_center_s','center_offset_s','arrival_pick_status','arrival_search_start_s',
              'arrival_search_end_s','target_gate_start_s','target_gate_end_s','requested_start_s','requested_end_s','start_s','end_s','status'])
    role_names={'family':'射线组','station':'台站名','station_role':'本条记录的用途','role_group':'用途来源分组',
                'distance_range_selected':'是否在指定距离范围内','incident_beam_region':'初始入射波束覆盖情况','outside_beam_selection_margin':'满足波束外选线净距',
                'station_role_count':'该台站承担的处理用途数量','target_phase':'唯一目标震相','predicted_peak_s':'目标理论到时_s',
                'selected_center_s':'实际选窗中心_s','start_s':'实际窗起点_s','end_s':'实际窗终点_s',
                'direct_peak_time_s':'目标直接峰时_s','hilbert_peak_time_s':'目标包络峰时_s',
                'direct_peak':'目标直接幅值','hilbert_peak':'目标包络幅值','usable':'是否通过检查','status':'处理状态'}
    write_csv(out/'台站用途与拾取结果.csv',[{name:r.get(key,'') for key,name in role_names.items()} for r in rows],list(role_names.values()))
    pair_report=export_pairs(rows,model,o,out)
    pairs,summaries,refs=coefficients_from_rows(rows,model,o)
    # 单参考仅供可选波形诊断选择参考台站，不再输出旧版系数表或重复图件。
    if o["generate_plots"]:
        wave_dir=out/"waveform_checks";wave_dir.mkdir(exist_ok=True)
        # 单次运行的诊断图清单，旧图不会被误列为当前结果。
        plot_rows=[]
        original_plan={(r['family'],r['phase'],r['station']):r for r in plan}
        for family in refs:
            for phase in model.phases:
                available=[r for r in rows if r["family"]==family and r["phase"]==phase and r["direct_peak"]!=""]
                available.sort(key=lambda r:(not r["usable"],abs(r["distance_to_interface_m"]-o["reference_target_distance_m"])))
                if phase.startswith("I_") and refs[family]:
                    available=[refs[family]]+[r for r in available if r["station"]!=refs[family]["station"]]
                for row in available[:o["max_waveform_plots_per_phase"]]:
                    _,payload=process_row(original_plan[(row['family'],row['phase'],row['station'])],model,index,o)
                    if payload is not None: plot_waveform(row,payload,wave_dir,o);plot_rows.append(row)
        write_csv(out/"plotted_waveforms.csv",plot_rows,["family","phase","station","status"])
    counts={}
    for row in rows:
        status=row["status"].split(":")[0];counts[status]=counts.get(status,0)+1
    report=dict(run_mode=o['run_mode'],input_data=str(data),setup=str(setup),recognized_component_files=index.count,selected_unique_stations=len(selected),
                pairing=pair_report,time_reference=clock,input_diagnostics=diagnostics,theory_context=context,
                case_consistency=consistency,
                target_phase_window_limit=o['limit_window_to_target_phase'],
                role_assignment='station_groups.csv; each family/phase/station is a separate target; other arrivals are interference references only',
                phase_station_rows=len(plan),status_counts=counts,quantity=o["trace_quantity"],families=o["ray_family"],
                accepted_coefficient_rows=pair_report["accepted_pairs"],
                raw_ratio_rows=pair_report["measured_pairs"],
                exclude_outer_boundary_stations=o['exclude_outer_boundary_stations'],
                plot_files=list(GENERATED_PLOT_PATHS),
                reference_quality={f:('accepted' if r and r['usable'] else 'unreviewed' if r else 'unavailable') for f,r in refs.items()},
                reference_stations={f:r["station"] if r else None for f,r in refs.items()},
                interpretation="Both measured ratios are coefficient magnitudes, not signed/complex coefficients or energy ratios.",
                limitation="Hilbert removes constant phase-rotation peak bias; it does not separate overlapping modes or repair finite-aperture/boundary contamination.",
                timing="Use the trace time column. SOURCE_PEAK_TIME_S is explicit; no automatic time shift is added.")
    (out/"run_report.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    lines=["本次后处理结果",f"输入物理量：{o['trace_quantity']}；射线选择：{o['ray_family']}",
           f"候选理论台站：{len(selected)}；识别的分量文件：{index.count}",
           "直接方法=max|极化投影波形|之比；Hilbert方法=max|解析信号|之比。两列均为系数模。",""]
    lines += [f"台站组合：{pair_report['total_pairs']}；可计算：{pair_report['measured_pairs']}；通过检查：{pair_report['accepted_pairs']}",
              "最终结果在上一级：反透射系数结果.xlsx、反透射系数结果.png。",
              "最优台站全时程.png完整展示最佳组合的入射/出射波形，标出窗内直接峰与包络峰的幅值、峰时和拾取窗。",
              "Excel的全部台站组合页逐对列出直接/Hilbert系数及质量原因；系数统计页给出两种方法的最大、最小、均值、中位数和方差。",
              "当前目录只保存计算过程，旧版输出子目录是历史归档，不是本次结果。",
              "amplitudes.csv保存各台站目标震相的幅值与状态；window_plan.csv保存提取窗。",
              "空白数值表示不可用，绝不是0。合格统计只使用accepted=1；原始统计含待核对比值。",
              "按理论符合度筛选的图示不能代替波形质量检查；Hilbert不能消除混波。",
              *diagnostics['warnings']]
    (out/"结果说明.txt").write_text("\n".join(lines),encoding="utf-8")
    report['best_station_pairs']=best_station_pairs(read_csv(out/'station_pairs.csv'),selected_families(model,o),model.phases,o)
    report['best_pick_details'],pick_payloads=best_pick_details(report['best_station_pairs'],rows,model,index,o)
    workbook=export_result_workbook(out,root,o,report)
    plot_hilbert_error_combinations(out,model,o,final_path=root/'反透射系数结果.png',pick_payloads=pick_payloads)
    report['best_full_record_image']=str(plot_best_full_records(root,model,o,report['best_station_pairs'],pick_payloads))
    report['final_workbook']=str(workbook);report['final_image']=str(root/'反透射系数结果.png')
    report['process_directory']=str(out)
    report['plot_files']=list(GENERATED_PLOT_PATHS)
    (out/'run_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    return report,summaries


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--family",choices=["central","first","auto","both","all"])
    parser.add_argument('--center',choices=['theory','waveform'],help='覆盖顶部WINDOW_CENTER_MODE')
    parser.add_argument("--setup");parser.add_argument("--data");parser.add_argument("--output")
    theory_args=parser.add_mutually_exclusive_group()
    theory_args.add_argument('--theory-directory',help='对比理论总目录或upper_P/lower_P/lower_SV子目录')
    theory_args.add_argument('--internal-theory',action='store_true',help='使用原内置理论计算，不读取外部理论表')
    parser.add_argument("--windows");parser.add_argument("--prepare",action="store_true")
    parser.add_argument("--no-plots",action="store_true",help="关闭额外过程诊断图；最终汇总PNG仍生成")
    parser.add_argument("--analyze",action="store_true",help="读取真实波形并计算；需配置同次模拟的数据目录")
    args=parser.parse_args();o=settings()
    for arg,key in [("family","ray_family"),("setup","setup_directory"),("data","data_directory"),("output","output_directory"),("windows","window_overrides_csv")]:
        value=getattr(args,arg)
        if value is not None: o[key]=value
    if args.prepare and args.analyze:parser.error("--prepare与--analyze不能同时指定")
    if args.prepare:o["run_mode"]="prepare"
    if args.analyze:o["run_mode"]="analyze"
    if args.center:o['window_center_mode']=args.center
    if args.theory_directory is not None:o['theory_directory']=args.theory_directory
    if args.internal_theory:o['theory_directory']=None
    if args.no_plots:o["generate_plots"]=False
    try:
        report,_=run(o)
    except Exception as exc:
        out=resolve(o['output_directory'])/'过程文件';out.mkdir(parents=True,exist_ok=True)
        error={'status':'failed','error_type':type(exc).__name__,'message':str(exc),'time':datetime.now().isoformat(),
               'note':'本次未完成；目录中已有结果可能来自之前运行，不应当作本次成功结果。'}
        (out/'last_run_status.json').write_text(json.dumps(error,ensure_ascii=False,indent=2),encoding='utf-8')
        print('本次失败，原因已保存：',out/'last_run_status.json')
        raise
    (resolve(o['output_directory'])/'过程文件'/'last_run_status.json').write_text(json.dumps({'status':'completed','run_mode':o['run_mode'],'time':datetime.now().isoformat()},indent=2),encoding='utf-8')
    print(json.dumps(report,ensure_ascii=False,indent=2))
    print("输出目录：",resolve(o["output_directory"]))


if __name__ == "__main__":
    if sys.version_info < (3,10): raise RuntimeError("需要Python 3.10或更高版本")
    main()
