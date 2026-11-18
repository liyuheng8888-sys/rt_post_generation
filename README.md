# 三层软弱夹层波包反透射验证

本分支把理论解、震源台站生成和后处理统一为波包验证流程。三个入口共用根目录的 `wave_packet_config.json`，不再分别修改材料、入射角和主频。原有两层程序保持不变，`main` 的 `1.0.0` 标签保留。

当前支持无材料耗散、各向同性、三个弹性固体层，支持上/下侧 P 或 SV 入射。主指标为**位移平方积分均值之比开方**，仅使用 X/Z，不使用 Y、Hilbert 峰值或单频峰值。中间层内的相干多次反透射通过完整三层复系数保留。该指标是位移波包等效幅值比，不是物理能流系数，也不代表逐频相位或结构峰值已验证。

## 三个入口

|用途|程序|
|---|---|
|完整三层频率响应及实测频谱加权理论|`含软弱层的三层算例/理论解/rt_coefficients_theorical_solution.py`|
|配对源输入及多深度台站|`含软弱层的三层算例/generate_stations_source/generate_sources_stations.py`|
|波包分解、积分比与核查|`含软弱层的三层算例/postprocess_rt/postprocess_rt.py`|

需要 Python 3.10 或更新版本。在仓库根目录安装依赖：

```shell
python -m pip install -r requirements.txt
```

## 1 修改统一配置并生成配对输入

```shell
python "含软弱层的三层算例/generate_stations_source/generate_sources_stations.py"
```

默认生成到 `result/WavePacket_upper_60_P_95Hz/`。所有位置参数采用全局坐标，埋深从外模型顶面向下量。当前上下界面埋深为 1999、2000 m。

|阵列|距界面距离|实际埋深|用途|
|---|---|---|---|
|上侧|20、40 m|1979、1959 m|入射侧为上侧时提取入射及反射；反之提取透射|
|上侧|70 m|1929 m|独立预测核查，不参与主拟合|
|下侧|20、40 m|2020、2040 m|入射侧为上侧时提取透射；反之提取入射及反射|
|下侧|70 m|2070 m|独立预测核查，不参与主拟合|

主阵列 Y=300 m，X=1885～2185 m，间距 2.5 m，共 726 站。固定主统计段为 X=1985～2085 m、41 个位置；100、200、300 m 全部范围作为敏感性对照。另在 Y=250、350 m，各选 5 个 X 位置与同样六个深度，增加 60 个侧向诊断站，总计 **786 站**。诊断站也仅读取 X/Z。

输出包括 `station_coordinates.csv`、布站图、几何窗口表、配置快照以及：

```text
three_layer/DATA/STATIONS
three_layer/DATA/FORCESOLUTION
three_layer/DATA/Par_file.wave_packet_required
three_layer/materials/nummaterial_velocity_file
uniform/DATA/STATIONS
uniform/DATA/FORCESOLUTION
uniform/DATA/Par_file.wave_packet_required
uniform/materials/nummaterial_velocity_file
```

两组的 STATIONS、FORCESOLUTION 逐字节相同；均匀模型所有材料号都使用入射侧介质的参数。材料号、网格拓扑和分区不变。`Par_file.wave_packet_required` 是需合入既有 Par_file 的参数片段，**不是完整的超算输入目录**。把材料文件放入你的实际网格读取路径，再重新生成相应数据库；不要把它留在无效的 prepare 快照目录。Q=9999 按无衰减处理，求解器同时关闭 ATTENUATION。

复用相同网格、震源、时间步、输出采样、边界与台站，仅改变材料值。程序不会修改网格，也不会提交超算作业。默认要求保存原生位移 `.semd`；同时建议保存 `.semv`、`.sema` 供后续工程指标分析。输出完整波包，包含层内多次波的尾部，但避免把外边界或源端部干扰纳入目标窗。

默认台站布局针对当前波包与几何；改变入射角、波速或有效频带后，应重新核查空间采样、模式分解条件数、有限源覆盖和几何窗口。20/40/70 m 不是适用于所有模型的固定物理常数。第三排预测和侧向诊断用来检验假设，不通过理论误差挑选台站。

## 2 计算理论频率响应

```shell
python "含软弱层的三层算例/理论解/rt_coefficients_theorical_solution.py"
```

生成 `result/<case_name>/theory/coefficients_frequency_spectrum.csv`，默认扫描 1～300 Hz。这个扫描只用于检查频率响应，后处理会在记录的实际 FFT 频点直接求解，不通过粗频率表插值。若没有实际入射谱，程序只输出 C(f)，不会把一个假定源谱冒充实测波包权重。

## 3 返回超算波形后计算波包系数

把两组对应波形分别放到 `three_layer/OUTPUT_FILES/` 和 `uniform/OUTPUT_FILES/`，文件名采用 SPECFEM 的 `网络.台站.分量.semd`，默认分量码为 FXX、FXZ。两组实际完整 Par_file 放在各自 DATA 下，以便比较；缺失时会在审计文件明确提示。实际材料参数仍需与已运行模型核对。

```shell
python "含软弱层的三层算例/postprocess_rt/postprocess_rt.py"
```

可用 `--config`、`--case-dir`、`--target`、`--reference`、`--target-traces`、`--reference-traces`、`--output` 指定其他路径。后处理不依赖历史结果表、旧后处理脚本或旧目录。缺少配置要求的台站、分量、时段，或配对输入不一致时会报错，不会自动选取更接近理论的台站。波形时间轴必须一致，Z 正方向由统一配置指定。

默认流程为：入射侧三层场减均匀参考场 → 固定深度逐频分解反射 P/SV；出射侧直接分解三层场 → 保留倏逝项作为分解中的附加模态，但不报告其远场传播系数。均匀场的入射模态是分母。逐频传播相位归算至相应界面后，对各传播模态的位移波包积分。

设 E 是每个位置的位移平方积分，则数值结果为 `sqrt(sum(E_out)/sum(E_in))`。理论结果为 `sqrt(sum_s,sum_f |U_in|² |C(f)|² / sum_s,sum_f |U_in|²)`。权重由测得的入射谱决定，无按理论误差调权。默认采用主频周期数控制几何窗前后留白，但主指标不是该主频的单频系数。全部 5 种窗宽、2 种窗形和3种统计范围都输出。

主要结果在 `result/<case_name>/packet_results/`：

|文件|内容|
|---|---|
|`RT_packet_primary.csv`|固定主窗口与主台站段的波包结果|
|`RT_packet_sensitivity.csv`|所有预设窗口、窗形和统计范围|
|`quality.csv`|拟合残差、独立深度预测残差、积分初始状态和窗边缘贡献|
|`lateral_diagnostics.csv`|侧向诊断点的入射整体幅值差与波包系数|
|`reference_spectrum.csv`|主统计段实测入射谱平方，供理论程序独立加权|
|`frequency_response.csv`|数值模态谱平方与完整三层理论频率响应|
|`center_modal_packets.npz`、PNG|中心位置模态波包，保留物理幅值，不逐道归一化|
|`audit.json`、`config_used.json`、`report.md`|审计、实际参数与结果说明|

再次运行理论入口时，可独立复算同一入射频谱对应的理论波包系数：

```shell
python "含软弱层的三层算例/理论解/rt_coefficients_theorical_solution.py" --reference-spectrum "result/WavePacket_upper_60_P_95Hz/packet_results/reference_spectrum.csv"
```

如历史记录只有加速度，将配置的 `trace_quantity` 改为 `acceleration`，程序会先按初始静止条件在完整记录上连续积分两次，再截窗；不作任意去趋势或带通。初始漂移或静态分量明显时不能仅凭小的波包比误差认定积分可靠，应优先使用原生位移输出。

## 验证和历史数据复算

```shell
python -m unittest discover -s tests -v
python -m unittest discover -s "含软弱层的三层算例/理论解" -p test_critical_angles.py -v
python tests/run_existing_data_regression.py --data-root "原三层算例/result"
```

`examples/existing_95Hz_acceleration.json` 仅用于现有两深度加速度数据的回归；它明确关闭不存在的70 m留出排及侧向点，不冒充新增阵列结果。默认新配置始终要求六排主阵列。原始波形、网格、计算输出与临时文件不上传 Git。

方法参考：Peng, H., Sripanich, Y., Vasconcelos, I., & Trampert, J. (2021). *Benchmarking wave equation solvers using interface conditions: the case of porous media*. GJI, 224(1), 355–376. https://doi.org/10.1093/gji/ggaa468 。本项目采用其位移平方积分思想，双界面完整理论、配对参考与多深度分解为本项目的适配实现。
