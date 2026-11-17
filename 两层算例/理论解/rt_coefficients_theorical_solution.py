#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
两介质平面波反射/透射系数计算器（Python 3.10 兼容版）
- 两侧介质类型均可独立设置为 solid（固体）或 fluid（流体）
- 可选择从介质1（上侧）或介质2（下侧）入射
- 固体入射可选 P 或 SV；流体仅允许 P 波
- 只计算当前选定的这一对介质，不默认海水，也不自动比较多种海床
- 输出一张系数模曲线图、系数 TXT 表、临界角 TXT 表和参数快照

振幅约定：各模式的标量位移振幅，P/SV 极化方向按脚本中的状态向量定义。
图和表输出系数模 |C|；TXT 同时保留复系数实部、虚部和模。
"""

from pathlib import Path
import math
import sys

if sys.version_info < (3, 10):
    raise RuntimeError("本脚本需要 Python 3.10 或更高版本。")

import numpy as np
import matplotlib
matplotlib.use("Agg")  # 适用于 Linux 服务器/无图形界面环境
import matplotlib.pyplot as plt

# ====================== 用户可编辑区 ======================
# type 可以分别设为 "solid" 或 "fluid"；两侧可任意组合。
# 固体需要 rho、vp、vs；流体只需要 rho、vp，vs 固定为 0。
MEDIA_1 = {
    "name": "上层流体",
    "type": "fluid",
    "rho": 1000.0,       # kg/m^3
    "vp": 1500.0,        # P波速度 m/s
    "vs": 0,        # S波速度 m/s
}
MEDIA_2 = {
    "name": "下层固体",
    "type": "solid",
    "rho": 2500.0,
    "vp": 4500.0,
    "vs": 2100.0,
}

INCIDENT_SIDE = 2       # 1=从介质1入射；2=从介质2入射
INCIDENT_WAVE = "P"     # 固体可选 "P" 或 "SV"；流体只能选 "P"

ANGLE_STEP_DEG = 0.5
MAX_ANGLE_DEG = 89.5    # 相对界面法线的入射角；不取90°端点
OUTPUT_DIR = Path(__file__).resolve().parent / "rt_coefficients_output"
# =========================================================


def validate_media(m, idx):
    kind = str(m.get("type", "")).lower()
    if kind not in ("solid", "fluid"):
        raise ValueError(f"MEDIA_{idx}['type'] 必须是 'solid' 或 'fluid'")
    m["type"] = kind
    for key in ("rho", "vp"):
        if key not in m or float(m[key]) <= 0:
            raise ValueError(f"MEDIA_{idx} 的 {key} 必须是正数")
    if kind == "solid":
        if "vs" not in m or float(m["vs"]) <= 0:
            raise ValueError(f"固体 MEDIA_{idx} 必须提供大于0的 vs")
        if float(m["vp"]) <= float(m["vs"]):
            raise ValueError(f"固体 MEDIA_{idx} 通常应满足 vp > vs")
    else:
        m["vs"] = 0.0


def modes(m):
    return ["P", "SV"] if m["type"].lower() == "solid" else ["P"]


def wave_speed(m, mode):
    if mode == "P":
        return float(m["vp"])
    if mode == "SV" and m["type"].lower() == "solid":
        return float(m["vs"])
    raise ValueError(f"{m['name']} 不支持 {mode} 波")


def q_vertical(v, p, direction):
    """
    垂向慢度 q。时间因子采用 exp[i*w*(p*x+q*z-t)]，z向下为正。
    主平方根使倏逝波的根具有正虚部；direction 控制向上/向下的衰减分支。
    """
    q0 = np.sqrt(complex(1.0 / (v * v) - p * p, 0.0))
    return direction * q0


def full_state(m, mode, p, direction):
    """
    返回全状态列 [ux, uz, sigma_xz/(i*w), sigma_zz/(i*w)]。
    模态标量振幅按位移极化方向归一化。
    """
    vp = float(m["vp"])
    rho = float(m["rho"])
    q = q_vertical(wave_speed(m, mode), p, direction)

    if m["type"].lower() == "fluid":
        # 流体仅有P波；剪切应力为零。
        return np.array([vp*p, vp*q, 0j, rho*vp], dtype=complex)

    vs = float(m["vs"])
    lam = rho * (vp*vp - 2.0*vs*vs)
    mu = rho * vs*vs

    if mode == "P":
        ux = vp*p
        uz = vp*q
        txz = 2.0*mu*vp*p*q
        tzz = lam/vp + 2.0*mu*vp*q*q
        return np.array([ux, uz, txz, tzz], dtype=complex)

    if mode == "SV":
        # SV极化取 e_SV=(-vs*q, vs*p)，与波矢 (p,q) 正交。
        ux = -vs*q
        uz = vs*p
        txz = mu*vs*(p*p - q*q)
        tzz = 2.0*mu*vs*p*q
        return np.array([ux, uz, txz, tzz], dtype=complex)

    raise ValueError(f"未知波型：{mode}")


def interface_state(m1, m2, m, mode, p, direction):
    """根据界面类型选择需要连续的状态分量。"""
    state = full_state(m, mode, p, direction)
    both_solid = m1["type"] == "solid" and m2["type"] == "solid"
    both_fluid = m1["type"] == "fluid" and m2["type"] == "fluid"
    if both_solid:
        # ux、uz、剪切应力、法向应力连续
        return state[[0, 1, 2, 3]]
    if both_fluid:
        # 两侧均为流体：法向位移和法向应力连续
        return state[[1, 3]]
    # 流固界面：法向位移、法向应力连续；固体侧剪切应力为零。
    return state[[1, 3, 2]]


def solve_linear(A, b):
    """行缩放后求解；临界角附近若矩阵退化则使用最小二乘兜底。"""
    A = np.asarray(A, dtype=complex)
    b = np.asarray(b, dtype=complex)
    scales = np.max(np.abs(A), axis=1)
    scales[scales == 0] = 1.0
    As = A / scales[:, None]
    bs = b / scales
    try:
        return np.linalg.solve(As, bs)
    except np.linalg.LinAlgError:
        return np.linalg.lstsq(As, bs, rcond=None)[0]


def solve_one_angle(m1, m2, incident_side, incident_wave, angle_deg):
    """求指定入射角下，所有反射和透射模态的复振幅系数。"""
    media = {1: m1, 2: m2}
    inc_m = media[incident_side]
    other_side = 2 if incident_side == 1 else 1
    trans_m = media[other_side]

    if incident_wave not in modes(inc_m):
        raise ValueError(f"介质{incident_side}为{inc_m['type']}，不支持 {incident_wave} 入射")

    v_inc = wave_speed(inc_m, incident_wave)
    theta = math.radians(angle_deg)
    p = math.sin(theta) / v_inc

    # z向下为正：从介质1入射向下(+1)，从介质2入射向上(-1)。
    incident_direction = +1 if incident_side == 1 else -1
    reflected_direction = -incident_direction
    transmitted_direction = incident_direction

    inc_state = interface_state(
        m1, m2, inc_m, incident_wave, p, incident_direction
    )

    ref_modes = modes(inc_m)
    trans_modes = modes(trans_m)
    ref_states = [
        interface_state(m1, m2, inc_m, mode, p, reflected_direction)
        for mode in ref_modes
    ]
    trans_states = [
        interface_state(m1, m2, trans_m, mode, p, transmitted_direction)
        for mode in trans_modes
    ]

    if incident_side == 1:
        # state_1 - state_2 = 0:
        # B_ref*x_ref - B_trans*x_trans = -state_inc
        A = np.column_stack(ref_states + [-s for s in trans_states])
        rhs = -inc_state
    else:
        # state_1 - state_2 = 0:
        # -B_ref*x_ref + B_trans*x_trans = state_inc
        A = np.column_stack([-s for s in ref_states] + trans_states)
        rhs = inc_state

    coeffs = solve_linear(A, rhs)
    result = {}
    for i, mode in enumerate(ref_modes):
        result[f"R_{mode}_M{incident_side}"] = coeffs[i]
    for j, mode in enumerate(trans_modes):
        result[f"T_{mode}_M{other_side}"] = coeffs[len(ref_modes) + j]
    return result


def get_critical_angles(m1, m2, incident_side, incident_wave):
    """对当前入射情况，逐一列出所有反射/透射模式的临界角或无临界角。"""
    media = {1: m1, 2: m2}
    inc_m = media[incident_side]
    other_side = 2 if incident_side == 1 else 1
    trans_m = media[other_side]
    v_inc = wave_speed(inc_m, incident_wave)
    rows = []

    for mode in modes(inc_m):
        v_out = wave_speed(inc_m, mode)
        label = f"反射 {mode}（介质{incident_side}）"
        if v_out > v_inc + 1e-10:
            angle = math.degrees(math.asin(v_inc / v_out))
            rows.append((label, v_inc, v_out, angle, "存在"))
        else:
            rows.append((label, v_inc, v_out, None, "无有限临界角（出射波速度不大于入射波速度）"))

    for mode in modes(trans_m):
        v_out = wave_speed(trans_m, mode)
        label = f"透射 {mode}（介质{other_side}）"
        if v_out > v_inc + 1e-10:
            angle = math.degrees(math.asin(v_inc / v_out))
            rows.append((label, v_inc, v_out, angle, "存在"))
        else:
            rows.append((label, v_inc, v_out, None, "无有限临界角（出射波速度不大于入射波速度）"))
    return rows


def mode_label(key):
    kind, mode, side = key.split("_")
    action = "反射" if kind == "R" else "透射"
    return f"{action}{mode}（介质{side}）"


def main():
    validate_media(MEDIA_1, 1)
    validate_media(MEDIA_2, 2)
    if INCIDENT_SIDE not in (1, 2):
        raise ValueError("INCIDENT_SIDE 只能设为 1 或 2")
    if ANGLE_STEP_DEG <= 0 or not (0 < MAX_ANGLE_DEG < 90):
        raise ValueError("ANGLE_STEP_DEG 必须大于0，MAX_ANGLE_DEG 必须在0到90度之间")
    if INCIDENT_WAVE.upper() not in ("P", "SV"):
        raise ValueError("INCIDENT_WAVE 只能设为 'P' 或 'SV'")
    inc_m = MEDIA_1 if INCIDENT_SIDE == 1 else MEDIA_2
    incident_wave = INCIDENT_WAVE.upper()
    if incident_wave not in modes(inc_m):
        raise ValueError(f"介质{INCIDENT_SIDE}是{inc_m['type']}，只能选择 P 波入射")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    angles = np.round(np.arange(0, MAX_ANGLE_DEG + 1e-9, ANGLE_STEP_DEG), 8)
    rows = []
    for angle in angles:
        row = {"angle_deg": float(angle)}
        row.update(solve_one_angle(MEDIA_1, MEDIA_2, INCIDENT_SIDE, incident_wave, float(angle)))
        rows.append(row)

    # 当前这对介质只生成一套系数结果图。
    coeff_keys = [k for k in rows[0].keys() if k != "angle_deg"]
    plt.rcParams["font.sans-serif"] = [
        "Noto Sans CJK SC", "Noto Sans CJK JP", "SimHei", "DejaVu Sans"
    ]
    plt.rcParams["axes.unicode_minus"] = False
    fig, ax = plt.subplots(figsize=(10.5, 6.5))
    for key in coeff_keys:
        y = np.array([abs(r[key]) for r in rows], dtype=float)
        y[~np.isfinite(y)] = np.nan
        ax.plot(angles, y, linewidth=1.7, label=mode_label(key))

    crit_rows = get_critical_angles(MEDIA_1, MEDIA_2, INCIDENT_SIDE, incident_wave)
    for label, _, _, critical, status in crit_rows:
        if critical is not None and critical <= MAX_ANGLE_DEG:
            ax.axvline(critical, linestyle=":", linewidth=1.2, alpha=0.8)
            ax.text(critical, 0.98, f"{critical:.2f}°",
                    transform=ax.get_xaxis_transform(), rotation=90,
                    va="top", ha="right", fontsize=8)

    ax.set_xlabel("入射角（°，相对界面法线）")
    ax.set_ylabel("反射/透射系数模 |C|")
    ax.set_title(
        f"{MEDIA_1['name']}（介质1）—{MEDIA_2['name']}（介质2）\n"
        f"介质{INCIDENT_SIDE}的{incident_wave}波入射"
    )
    ax.set_xlim(0, MAX_ANGLE_DEG)
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    plot_path = OUTPUT_DIR / "reflection_transmission_coefficients.png"
    fig.savefig(plot_path, dpi=220, bbox_inches="tight")
    plt.close(fig)

    # 系数 TXT：每个模态列出实部、虚部和模，便于识别超临界后的复系数。
    table_path = OUTPUT_DIR / "reflection_transmission_coefficients.txt"
    columns = []
    for key in coeff_keys:
        columns.extend([f"{key}_real", f"{key}_imag", f"{key}_abs"])
    lines = [
        "# 两介质平面波反射/透射系数表（TSV制表符分隔）",
        f"# 介质1：{MEDIA_1['name']} type={MEDIA_1['type']} rho={MEDIA_1['rho']} kg/m^3 Vp={MEDIA_1['vp']} m/s Vs={MEDIA_1.get('vs', 0)} m/s",
        f"# 介质2：{MEDIA_2['name']} type={MEDIA_2['type']} rho={MEDIA_2['rho']} kg/m^3 Vp={MEDIA_2['vp']} m/s Vs={MEDIA_2.get('vs', 0)} m/s",
        f"# 入射：介质{INCIDENT_SIDE}，{incident_wave}波；入射角范围0~{MAX_ANGLE_DEG}度，步长{ANGLE_STEP_DEG}度",
        "# 系数为位移极化方向上的复振幅系数；_abs为系数模。角度相对界面法线。",
        "\t".join(["angle_deg"] + columns),
    ]
    for row in rows:
        vals = [f"{row['angle_deg']:.5f}"]
        for key in coeff_keys:
            value = complex(row[key])
            vals.extend([f"{value.real:.12g}", f"{value.imag:.12g}", f"{abs(value):.12g}"])
        lines.append("\t".join(vals))
    table_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    # 临界角 TXT：明确列出存在/不存在的各个出射模式临界角。
    critical_path = OUTPUT_DIR / "critical_angles.txt"
    crit_lines = [
        "临界角计算结果",
        "临界角相对界面法线定义；满足 sin(theta_c)=V_incident/V_outgoing。",
        "仅当出射模式速度大于入射模式速度时，0~90度内存在有限临界角。",
        "",
        f"介质1：{MEDIA_1['name']}，类型={MEDIA_1['type']}，rho={MEDIA_1['rho']}，Vp={MEDIA_1['vp']}，Vs={MEDIA_1.get('vs', 0)}",
        f"介质2：{MEDIA_2['name']}，类型={MEDIA_2['type']}，rho={MEDIA_2['rho']}，Vp={MEDIA_2['vp']}，Vs={MEDIA_2.get('vs', 0)}",
        f"入射：介质{INCIDENT_SIDE}，{incident_wave}波",
        "",
        "出射模式\\t入射波速度(m/s)\\t出射波速度(m/s)\\t临界角(deg)\\t状态",
    ]
    for label, vi, vo, critical, status in crit_rows:
        crit_text = f"{critical:.8f}" if critical is not None else "—"
        crit_lines.append(f"{label}\t{vi:.8g}\\t{vo:.8g}\\t{crit_text}\\t{status}")
    crit_lines += [
        "",
        "说明：超过某出射模式的临界角后，该模式的垂向慢度为虚数，对应倏逝波；系数一般为复数。",
        "系数图显示模 |C|，系数表保留实部、虚部和模。",
    ]
    critical_path.write_text("\n".join(crit_lines) + "\n", encoding="utf-8")

    params_path = OUTPUT_DIR / "model_parameters.txt"
    params = [
        "本次计算使用的介质参数",
        f"MEDIA_1: name={MEDIA_1['name']}, type={MEDIA_1['type']}, rho={MEDIA_1['rho']} kg/m^3, vp={MEDIA_1['vp']} m/s, vs={MEDIA_1.get('vs', 0)} m/s",
        f"MEDIA_2: name={MEDIA_2['name']}, type={MEDIA_2['type']}, rho={MEDIA_2['rho']} kg/m^3, vp={MEDIA_2['vp']} m/s, vs={MEDIA_2.get('vs', 0)} m/s",
        f"INCIDENT_SIDE={INCIDENT_SIDE}",
        f"INCIDENT_WAVE={incident_wave}",
        f"ANGLE_STEP_DEG={ANGLE_STEP_DEG}",
        f"MAX_ANGLE_DEG={MAX_ANGLE_DEG}",
        "",
        "编辑脚本顶部的 MEDIA_1、MEDIA_2、INCIDENT_SIDE、INCIDENT_WAVE 和角度参数即可重新计算。",
        "type='solid' 时需要填写 rho、vp、vs；type='fluid' 时填写 rho、vp，vs自动设为0。",
    ]
    params_path.write_text("\n".join(params) + "\n", encoding="utf-8")

    print("计算完成。本次只输出当前选定的这一对介质的一套结果：")
    for path in (plot_path, table_path, critical_path, params_path):
        print(path.resolve())
    print("\n临界角：")
    for label, _, _, critical, status in crit_rows:
        print(f"  {label}: {critical:.4f}°" if critical is not None else f"  {label}: {status}")


if __name__ == "__main__":
    main()
