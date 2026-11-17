# Python 3.10 震源与台站生成器

适用模型：1000 m正方体。用户输入、坐标表和图件统一采用埋深d，向下为正；顶面d=0 m、海水—海床界面d=500 m、底部d=1000 m。
本脚本将文件夹内张文月等（2021）文献图4的二维点源叠加思路沿 y 方向扩展为有限震源面。
生成输入、理论几何与图件；没有运行 SPECFEM3D，也没有据此证明数值误差达到文献的 0.5%。

## 运行

在本目录打开终端：

```powershell
py -3.10 -m pip install numpy matplotlib
py -3.10 generate_sources_stations.py
```

也可以用已经指向 Python 3.10 的 `python`。不需要 MATLAB、SciPy 或联网运行。
**默认直接修改Python脚本最上方的“用户编辑区”**：`USER_CONFIG`包含全部物理与文件格式参数，`OUTPUT_DIRECTORY`设置输出目录，`GENERATE_PLOTS`控制绘图，`PLOT_DPI`设置图片清晰度，`STATION_ZOOM_WIDTH_M`设置台站局部放大窗口宽度。各项已标注中文含义、单位与可选值。直接运行时不读取config.json。

如需使用独立JSON配置，可显式指定（此时JSON完整替代脚本顶部USER_CONFIG）：

```powershell
py -3.10 generate_sources_stations.py --config config.json --output output_my_case
```

修改参数后重新运行。相同输出目录内的同名生成文件会更新，原始文献与“理论解”不受影响。

## 参数入口：脚本顶部 USER_CONFIG

### 台站生成开关

在顶部“2. 接收器”中，以下开关独立控制是否生成台站，Python使用`True`/`False`：

| 开关 | 控制的台站 |
|---|---|
| generate_grid_stations | y=500剖面密铺网格 |
| generate_first_stations | first组入射、反射和透射射线 |
| generate_central_stations | central组入射、反射和透射射线 |
| generate_front_stations | front理论波阵面（不是first） |
| generate_interface_stations | 埋深500m的界面交线 |

例如，只生成central理论射线台站：

```python
"generate_grid_stations": False,
"generate_first_stations": False,
"generate_central_stations": True,
"generate_front_stations": False,
"generate_interface_stations": False,
```

默认全部True，保留原来的行为。旧JSON不含这些字段时也默认全部True；JSON配置中布尔值写小写true/false。
开关在台站采样与重合命名之前生效；关闭的组不会出现在STATIONS、坐标表、分组表或重命名记录中。保留下来的组仍遵循后生成名称覆盖规则，因此关闭某组可能改变重合位置的最终名字。应使用同一次生成的台站与分组文件配套进行新模拟和后处理，不要替换已有算例的元数据。

STATIONS_grid和STATIONS_lines仍会每次写出；相应类别关闭且无台站时文件为空，避免同一输出目录残留上次的内容。全部关闭时STATIONS也为空，但仍生成震源和理论几何，绘图支持零台站。
台站布置图只画实际生成的台站；首触、中央射线和理论波阵面仍写入报告，但不再生成geometry_section和wavefronts_y500图。
theory_report.json中的station_generation记录开关，station_family_membership_counts按归属统计每组的去重台站数；跨组共享台站会分别计数，因此各组归属数之和可能大于station_count。
仅生成central时，后处理选RAY_FAMILY="central"；仅生成first时选"first"；两组均生成时可选"both"。

以下参数均在脚本顶部，可直接修改；config.json仅保留为可选配置示例。命令行`--output`优先于顶部输出目录，`--no-plot`可强制关闭绘图。

| 参数 | 含义 |
|---|---|
| incident_side | lower从下层向右上入射；upper从上层向右下入射 |
| upper_source_anchor_x_depth_m | 上层源中心线端点[X,埋深]，单位m；脚本当前[400,50] |
| lower_source_anchor_x_depth_m | 下层源中心线端点[X,埋深]，单位m；脚本当前[750,1000]，可直接调整 |
| top_is_free_surface | 顶面是否为自由表面，默认True；上层入射时检查声学源 |
| bottom_is_free_surface | 底面是否为自由表面，默认False；下层入射时检查声学源；两项均不修改实际网格或求解器边界 |
| angle_deg | 投影与−x方向的夹角大小α，范围0≤α<90°；下层源线向浅处延伸，上层源线向深处延伸 |
| incident_wave | 入射介质为固体时可选P/SV；为流体时只允许P |
| source_spacing_m | 沿倾斜震源中心线的真实弧长间距，默认10，可改5 |
| source_y_spacing_m | 沿 y 的震源中心间距，默认10；全平面改5时这项也改5 |
| station_spacing_m | y=500 平面 X、埋深方向的底网格间距，默认10，可改5 |
| grid_station_prefix | 第一类网格台站前缀，默认G，例如G000001 |
| line_station_prefix | 第二类理论线台站前缀，默认L，例如L_first_I_P_000001 |
| source_interface_clearance_m | 震源中心线距离界面的最小竖向间隙，默认50 m，防止震源面穿越两种介质 |
| source_kind | auto：按入射介质自动选择；固体用paper_couples两对力偶，流体用acoustic_monopole标量源；固体还可选plane_force定向力 |
| couple_separation_m | 每对力偶中正负力间距，默认1 m；这是力偶尺寸，不是中心间距 |
| frequency_hz | Ricker 主频，默认5 Hz |
| enable_uniform_time_shift | True开启统一时移；False关闭，全部写0，默认False |
| uniform_time_shift_s | 开启后写入所有FORCE的同一time shift，单位s，必须有限且非负，默认0.0 |
| force_factor | 每条FORCE幅值，默认1e15；固体力单位N，声学标量源幅值需独立标定，不能与固体力幅值直接比较 |
| force_direction | 默认auto；仅plane_force可设置[E,N,Z_UP]自定义非零向量 |
| force_format | 11_lines 对应用户本次模板；legacy_10_lines 对应旧文件夹中无 source time function 行的格式 |
| force_depth_unit | 默认km，遵循官方埋深模式约定；m仅供明确修改过读取单位的本地版本 |
| snapshot_after_first_s | 有限孔径首触后多少秒计算front台站及报告中的理论波阵面，默认0.05 s |
| upper/lower | 上下介质的name、type、rho（kg/m³）、vp与vs（m/s） |

上下层均可设fluid或solid，也均可作为入射侧。fluid必须vs=0且只能P入射；solid须满足正密度、正剪切模量与正体积模量。
将`incident_side`改为`"upper"`，保留`source_kind="auto"`即可生成上层入射算例。源的位置、法向、反/透射介质、波速、系数、波阵面与图标题均随之切换。
上层端点直接填埋深，例如`upper_source_anchor_x_depth_m=[750,50]`表示X=750 m、埋深50 m。范围为`0<X≤1000`、`0≤埋深<500-clearance`。若上层为流体且顶面为自由表面，埋深必须大于0，否则生成器报错，防止出现声学源被置零的警告。原来的`upper_source_anchor_xz`参数不再接受，避免误读旧配置。

流体采用每中心一个声学单极子，输出仍为FORCESOLUTION。要求求解器`setup/constants.h`中的编译常量`DIPOLE_SOURCE_IN_FLUID=.false.`（所核对的官方默认值）。该常量不是Par_file参数。默认声学路径忽略方向向量，不能使用只改变方向正负的固体力偶配置；因此流体源型显式设置成paper_couples/plane_force会报错，提示使用auto/acoustic_monopole。
依据：[官方声学源初始化代码](https://github.com/SPECFEM/specfem3d/blob/master/src/specfem3D/setup_sources_receivers.f90)、[官方编译常量](https://github.com/SPECFEM/specfem3d/blob/master/setup/constants.h.in)。这是有限平面声源阵列的推广，并非声学域中逐项照搬文献力偶。

理论线上的台站固定每5 m一个，并补齐端点（最后一段可能小于5 m）；与底网格重合时去重，分组表保留所有归属。

## 上下层对称调整

只需用 `incident_side` 选择入射侧，调整对应的 `upper_source_anchor_x_depth_m` 或 `lower_source_anchor_x_depth_m`。两项都填 `[X,埋深]`（m），不使用向上为正的内部坐标；只校验和使用选中侧的端点。

- 上层端点：`0<X≤1000`，`0≤埋深<500-source_interface_clearance_m`。
- 下层端点：`0<X≤1000`，`500+source_interface_clearance_m<埋深≤1000`。
- 镜像位置示例：上层 `[400,50]` 对应下层 `[400,950]`。为保留既有算例，脚本没有自动改动下层原有的 `[750,1000]` 默认位置。
- `angle_deg`、`incident_wave`、源间距、界面间隙、源型、幅值、时移和全部台站参数由两侧共用。角度取正值，源线均向负X延伸并靠近界面；改变入射侧时，无需手动翻转角度或自动力方向。
- `top_is_free_surface` 对应上层，`bottom_is_free_surface` 对应下层。如果入射侧是流体且外侧边界为自由表面，端点必须放在介质内部。固体可选P/SV，流体只能P，这是介质差别，不是参数入口限制。
- `auto_reflected_mode` 两侧都可用：流体反射选 `R_P`；固体反射可选传播状态的 `R_P` 或 `R_SV`。自动选线仍要求存在满足净距和最少台站数的有效射线，否则会报错，不能保证任意几何都可选线。
- 旧JSON缺少下层端点或底面开关时，仍使用 `[750,1000]` 和 `False`；其余旧配置用法保持不变。

绘图只保留 `stations_y500.png/pdf` 和启用自动测线时的 `auto_rays_y500.png/pdf`。不再生成 `geometry_section`、`wavefronts_y500`（PNG/PDF均不生成）。旧输出目录中的历史图不会自动删除。`generate_front_stations` 与 `snapshot_after_first_s` 仍控制front台站和报告数据，不依赖绘图。

## 两类台站、重名覆盖与分组

**所有台站都在y=500 m这个唯一的X–埋深剖面内，包括网格台站、理论射线台站、波阵面台站和界面线台站。** 两类名称只区分采样方式，不代表位于两个空间平面。
这里的“界面线台站”是y=500剖面与埋深500 m介质界面的交线，即(x,500,500)，并不是在整个水平界面铺设二维台站。
STATIONS第三列是y，所有记录均为500。震源面仍按原要求沿y展开；震源与接收器是不同对象。

第一类先生成，覆盖y=500整个剖面，默认X、埋深间距10 m，由`station_spacing_m`修改，名称为`G000001`、`G000002`等。
第二类后生成，每条理论线按5 m采样，名称包含线的标识及在线上的序号，例如：

- `L_first_I_P_000001`：首触点对应的入射P波射线。
- `L_central_R_SV_000001`：源面几何中点对应的反射SV波射线。
- `L_front_T_P_000001`：指定时刻透射P波的波阵面。
- `L_interface_000001`：界面线。

可在顶部修改G/L前缀（1～4位英文字母或数字，首位为字母），两者不能相同或互为开头。
去重按(X,埋深)坐标四舍五入到小数点后6位，即1微米精度；y固定500。重合点只输出一条记录，使用**后生成的名称和坐标**。
生成顺序为：网格 → 首触入射线 → 中央入射线 → 各反/透射模态的首触、中央射线 → 入射波阵面 → 各出射波阵面 → 界面线。没有实射线的倏逝模态不参与射线生成。
因此第二类覆盖第一类，第二类内部也后者覆盖前者；界面线最后生成，重合时使用`L_interface_…`。

`STATIONS`是最终去重全集；`STATIONS_grid`与`STATIONS_lines`按最终名称归属拆分，二者无重复且合并正好等于全集。
网格中的重合位置仍有接收器，只是归入第二类并使用第二类名称；序号可能存在跳号。
`station_groups.csv`的station列始终是最终有效名字，original_name列保留原请求名字，一点可属于多条线；不能把original_name当成求解器实际台站名。
`station_name_replacements.csv`记录覆盖过程；该表的中间new_name还可能被后续覆盖，最终名字以STATIONS或station_groups.csv的station列为准。

## 坐标：震源、台站、坐标表与图件均用埋深

统一时间偏移可在脚本顶部设置，例如：

```python
"enable_uniform_time_shift": True,
"uniform_time_shift_s": 0.3,
```

此时每条FORCE的time shift都为0.3 s，包括每组力偶的全部四个力点；关闭开关时全部写0。不会给首源或任意单独震源另设0，源间相对时差始终为0。
官方SPECFEM3D说明要求至少一条time shift为0，因此全非零文件可能被拒绝或无法实现预期整体延时；生成器保留用户指定值并发出兼容性提示，实际行为需核对所用版本。依据：[官方FORCESOLUTION说明](https://specfem3d.readthedocs.io/en/latest/05_running_the_solver/)。
此参数仅控制文件字段，不修改求解器或其内部时间基准。理论图与报告的first_time/snapshot_time仍为相对源面发射的传播参考时间，未计入该字段；报告另存uniform_time_shift_written_s记录实际写入值。几何和台站位置不受共同时间偏移影响。

`STATIONS` 每行：

```text
台站名 DB y x 0.0 埋深_m
```

例如X=750 m、Y=500 m、埋深50 m的台站写为`名称 DB 500 750 0 50`。第5列固定为0，是格式占位；第6列直接填写相对模型顶面的埋深。海床界面埋深500 m，模型底部埋深1000 m。

源和台站共用同一个坐标开关，所以必须设置：

```text
SUPPRESS_UTM_PROJECTION = .true.
USE_SOURCES_RECEIVERS_Z = .false.
```

`FORCESOLUTION` 中 `latorUTM=y`、`longorUTM=x`，默认 `depth=埋深_m/1000`，单位km。例如埋深50 m写为`depth: 0.050000000000`。
源与台站的埋深基准均为模型顶面，不是海床。图件与CSV一律使用m，FORCESOLUTION遵循源文件单位km。
`component dir vect source Z_UP`是力方向分量，向上为正，不是位置字段；不得因位置改为埋深而随意翻转力方向。
用户模板的 `hdurorf0::` 改为标准单冒号 `hdurorf0:`。科学计数法e格式可被Fortran读取。
必须让源和台站共同使用埋深开关；`Par_file_required.txt` 是需要合并的参数片段，不是完整Par_file。

坐标与深度依据：[官方坐标约定](https://danielpeter.github.io/specfem3d/A_reference_frame)、[官方求解器输入说明](https://specfem3d.readthedocs.io/en/latest/05_running_the_solver/)。后者源字段列表有早期x/y文字混淆，本脚本采用坐标附录及台站示例一致的 longitude=x、latitude=y 约定。
旧数据FORCESOLUTION为10行格式，本次默认按用户新给出的11行格式生成；应与实际安装版本匹配。

## 几何公式和输出波型

以下公式全部采用(X,d)，d为埋深。设α为源面投影与−X方向夹角：

```text
下层入射：A=(lower_source_anchor_x_depth_m)，σ=-1
上层入射：A=(用户X,用户埋深)，σ=+1
源面切向 t=(-cosα, σ*sinα)
入射方向 n=(sinα, σ*cosα)
r_source(s)=A+s*t
s_max=min(A_x/cosα, (abs(500-A_d)-clearance)/sinα)
p=sinα/v_incident
x_hit=x_source+(500-d_source)*n_x/n_d
t_hit=(500-d_source)/(v_incident*n_d)
sinθ_j=p*v_j
反射方向=(p*v_j, -σ*sqrt(1-(p*v_j)^2))
透射方向=(p*v_j,  σ*sqrt(1-(p*v_j)^2))
```

α=0时跳过s_max的第二个限制。源面沿y从0延伸到1000；接收器全部固定y=500。
只在x_hit落在[0,1000]的有限孔径中选最早接触点；若整束波先离开侧边界则报错。
相对于水平界面法线的入射角等于α。

透射波在这里就是折射波，并非另一个独立波型。流体不生成SV；固固界面包含反射P/SV和透射P/SV。
平面各向同性P-SV问题不激发SH，所以不画SH射线。

若`p*v_j>1`，对应倏逝波，无实数体波射线；报告给出主频处的衰减长度
`1 / [2π f0 sqrt(p²−1/v_j²)]`，界面也按5 m布站。此时不能把界面附近测到的幅值当成远场传播波系数。
临界条件`p*v_j=1`会标记critical并给出掠射界面方向。

输出同时包含首触点与源面几何中点的入射/反射/透射射线，以及指定时刻的波阵面，全部5 m加密。
如果源面中点的入射射线在触界前已经离开模型，生成器仅保留有效孔径中央作为绘图参考，并设置`central_is_source_midpoint=false`；后处理会拒绝将它冒充源面中点。此时仍可使用first组。
波阵面由每个界面点h的到时推算：
`r_j(h,t)=(h,500)+v_j*(t-t_hit(h))*n_j`，并裁剪到介质、模型范围和已到达的有效孔径。
报告的到时是震源中心面出发的**相位传播参考时间**，未加SPECFEM的Ricker内部时间偏移。
Ricker并无严格有限起点；有限源边缘还会产生更早绕射，因此这里的“首触”特指理想平面波射线束的首触。

反/透射系数由位移与应力连续条件求解，与原“理论解”脚本相同约定：系数求解采用向下为正的竖向分量，绘图也采用向下为正的埋深。
流固界面连续法向位移、法向应力，并令固体剪应力为零；固固界面连续两分量位移与牵引。
复系数为各模态标量**位移振幅**比，不是能量比或压力比；符号依赖所定义的极化基底。

## 文献源型的三维推广

默认auto在固体中选择paper_couples，在每个中心配置两对正负力，P采用 `M=n⊗n+t⊗t`，SV采用 `M=t⊗n+n⊗t`；流体中自动选择前述单极子。
这是图4的二维膨胀源/双力偶在入射平面内的旋转与沿y方向复制，不是三维各向同性爆炸源（后者还需要yy分量）。
内部用对称力偶；中心在外边界时用保持同一力矩的单侧力偶，保证所有实际FORCE位于选定的入射介质内。下层中心线底端始终是(X,埋深)=(750,1000)，上层使用其专用端点；真实力点可偏离中心面约一个力偶间距。
使用plane_force时所有实际力点严格位于同一平面，自动P力方向为n、SV为t。
用户示例[1,-2,-1]会引入非零y极化及混合波型，若选择该向量，报告仍只表示配置的理想入射P/SV模式，不能解释为实际波场只有该模式。

各源同时起振。有限面及离散点源只能近似平面波，也会向反方向辐射；靠近边缘、震源或边界处不能直接当作纯平面波。
5 Hz时默认下层P波长约900 m，已接近模型边长；因此是否有足够干净的时空窗必须通过实际波形验证。
上下层端点均可配置；下层为兼容旧算例仍默认[750,1000]，可调整到介质内部，例如[400,950]。源面沿y仍覆盖0～1000 m。顶面和底面检查不会代替实际网格与PML检查。
改变源中心间距会改变单位面积的源数，改变力偶间距会改变力矩F*d；跨算例比较绝对幅值须做源强归一化，或使用同配置的均匀介质参考运行。

按文献第3节，后续应采用吸收边界、合适网格与时间步，按到时分离震相，再比较极化投影后的位移Hilbert包络峰值比。
若比较ux/uz分量需消除极化投影因子；流体压力与固体位移不能直接作振幅比。文献第8页明确排除了快衰波的这类平面传播波验证。

## 文件清单

- FORCESOLUTION、STATIONS：输入文件；STATIONS_grid和STATIONS_lines为去重后的两类台站子集。
- sources_depth.csv、stations_depth.csv：X、Y、burial_m（埋深m）核对表，不输出位置高程。
- station_groups.csv：同一台站在底网格、首触射线、中央射线、波阵面、界面的归属。
- station_name_replacements.csv：台站重合时的名称覆盖记录。
- theory_report.json：首触位置/到时、波速、角度、复系数、传播/倏逝状态、理论线端点。
- config_used.json：此次运行参数快照。
- stations_y500.png/pdf：y=500全部台站分布与局部放大；两类点使用不同颜色，并标出实际间距。
- Par_file_required.txt：需合并的开关。

所有文件均位于output子目录。没有更改旧算例，也没有生成或替换网格和材料文件。
绘图全部采用y=500的二维视角，纵轴是埋深，0在上、1000在下。旧坐标表保存在previous_coordinate_format子目录，仅供历史参考；本次使用根目录的*_depth.csv。旧图在previous_figures中仅供历史参考。

theory_report.json的位置字段使用`*_x_depth_m`，波的方向字段使用`direction_x_depth`。其中源线起止点和source_burial_range_m均为埋深m。直接运行使用脚本顶部的上层入射示例；可选config.json仍为下层入射示例。


## 自动选取初始入射波束之外的反射测线

顶部generate_auto_stations=True启用auto组；原first、central及其他开关独立。
只需要auto组时，关闭grid、first、central、front、interface五个开关。
auto_reflected_mode选择用于选线的反射模式（默认R_P）。
auto_candidate_step_m是界面交点扫描间隔（默认1m），auto_beam_clearance_m是波束边缘净距（默认10m），auto_min_stations是合格反射线最少台站数（默认10）。

判据使用整个有限源面的初始入射波束，不使用某一时刻的波阵面。设源线起点A、沿源线单位向量t、长度a，对反射线上点r计算s=(r-A)·t。仅保留s < -净距或s > a+净距的点。这个充分条件保守地排除了初始波束的整个横向条带。反射方向、转换波方向继续由Snell定律和已有弹性界面方程确定。

在源面可照射的界面区间，以auto_candidate_step_m扫描候选交点；各反射射线从界面起按5m整数距离采样，仅保留模型内且满足判据的点。台站数达到门槛才是合格候选。按交点X排序，选合格候选列表的中位项；偶数项取较小X的中位项。这是离散扫描近似，减小步长可提高选线分辨率，不声称精确枚举连续无限射线。选点对应连续有限源面的几何射线，不要求恰好穿过某一个离散FORCE位置。

选中后生成配套auto_I、auto_R、auto_T；反射台站从界面交点一直布置到模型边界，包含初始入射波束覆盖区。波束外判据只决定选哪条线，不裁剪最终输出。入射线用于后处理分母，透射线用于透射系数。所有传播反射分支都输出完整线段；主选反射模式的波束外部分保证达到auto_min_stations。不存在满足条件的线（例如正入射时反射线与入射束重合）时明确报错，不能伪造隔离台站。

输出auto_ray_candidates.csv保存所有扫描候选和合格标记，theory_report.json的auto_selection记录所选交点、门槛、扫描精度和数量；auto_rays_y500.png/pdf中橙色表示初始波束，整条反射线均显示台站，包括橙色覆盖区内的部分。示例在output_auto_example，未覆盖已有模拟元数据。

几何波束外不等于有限频率波场为零：有限孔径绕射、离散点源辐射、自由面和边界反射仍可能到达接收器。该模式减少理想直接入射波干扰，不能保证只存在某一反射波或保证P/S完全分离。有限孔径衍射参考：https://pmc.ncbi.nlm.nih.gov/articles/PMC10981538/

## 射线台站名中的距离

ray_station_naming="distance"为默认，"legacy"保留原序号方式。
例如L_A_RP_D000300p000：L为理论线前缀，A为auto（F为first，C为central），RP为反射P，D后表示沿该射线到界面交点的距离300.000m；p是小数点。不是垂直埋深，也不是从保留线段的第一个台站起算。名称保留毫米，精确距离写在ray_station_distances.csv中；舍入造成同名时增加n序号区分。
重合点仍以后生成名称覆盖；同一物理位置属于多个测线时，最终名字只代表最后所属线，必须用距离表中的group和精确距离查该点在其他测线中的含义。

新台站需要重新进行SPECFEM模拟，旧记录不能仅靠改名复用。后处理RAY_FAMILY="auto"读取auto组；SETUP_DIRECTORY必须指向本次新生成的目录。both继续表示first和central，不自动包含auto。


### 自动测线输出开关与完整反射线

在脚本顶部“2. 接收器”设置`"generate_auto_stations": True`输出auto台站；改成False不生成auto组。通过`--config`运行时，须在指定JSON里设置同名字段（JSON使用true/false）；脚本顶部配置此时不生效。示例config.json已补上此字段。
输出STATIONS_auto只含自动组台站（共享点使用去重后的名字），总STATIONS及STATIONS_lines也包含它们；关闭开关时STATIONS_auto重新写为空文件，分组表中没有auto组。
自动选线的波束外判据保持不变，最终反射台站按距界面0、5、10……米生成，并补齐模型边界端点（最后一段可小于5m）。报告selected_outside_beam_station_count为筛选时波束外有效点数，selected_reflected_station_count为主反射测线完整输出的点数，二者不应混淆。
波束覆盖段用于保留完整测线，在后处理中仍可能出现入射/反射重叠标记；位于界面上的点也会生成，但后处理仍按原规则判断是否可用。
