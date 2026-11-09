"""Build Chinese notes from the inspected PPT and saved aggregate results only."""
from pathlib import Path
import hashlib
import json
import re
import sys
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "outputs/_directory_audit/ppt_record_20261105"
PPT = ROOT / "11.4 fatty liver .pptx"
INVENTORY = json.loads((AUDIT / "ppt_source_inventory.json").read_text(encoding="utf-8"))
PROTECTED = {PPT: hashlib.md5(PPT.read_bytes()).hexdigest()}
SOURCES = []


def adopted(prefix):
    module = next(p for p in (ROOT / "outputs").iterdir() if p.name.startswith(prefix))
    run = (module / "LATEST_RUN.txt").read_text(encoding="utf-8-sig").strip().splitlines()[0]
    return module / "current" / run


def read(prefix, name):
    p = adopted(prefix) / "tables" / name
    PROTECTED[p] = hashlib.md5(p.read_bytes()).hexdigest()
    SOURCES.append(p)
    return pd.read_csv(p)


def f(value, digits=3):
    return "—" if pd.isna(value) else f"{float(value):.{digits}f}"


def q(value):
    return "—" if pd.isna(value) else f"{float(value):.3g}"


def table(headers, rows):
    rows = [[str(x).replace("|", "/").replace("\n", " ") for x in row] for row in rows]
    return "\n".join(["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"] + ["| " + " | ".join(row) + " |" for row in rows])


forest = read("25_", "Figure25_displayed_overall_SGBs_EN.csv").sort_values("display_order")
grade = read("33_", "B_Displayed_ordinal_SGBs_EN.csv").head(6)
grade_comparison = read("33_", "B_Fresh_vs_previous_results_EN.csv")
candidates = read("50_", "All_selected_candidates_EN.csv").sort_values("display_order")
host = read("50_", "Host_display_EN.csv")
host_summary = read("50_", "Host_summary_EN.csv")
host_filter = read("50_", "Host_filter_EN.csv")
paths = read("44_", "Figure44_Adult_Model4_all81_selection_order_EN.csv").sort_values("display_rank")
path_summary = read("44_", "Pathway_redesign_result_summary_EN.csv")
phylum = read("40_", "Significant_SGB_phylum_counts_and_percentages_EN.csv")
counts = read("24_", "00_Cohort_counts_EN.csv")
genus_summary = read("38_", "Analysis_sample_and_discovery_summary_EN.csv")
discovery = read("45_", "Analysis_populations_EN.csv")
ldm_summary = read("45_", "Model_and_global_summary_EN.csv")
assert len(forest) == 16 and len(grade) == 6 and len(candidates) == 10 and len(paths) == 81
assert int((candidates.P_BH < .05).sum()) == 2
assert grade_comparison.samples_filters_status_and_significance_identical.all()
assert len(INVENTORY["slides"]) == 7 and [r["slide"] for r in INVENTORY["slides"]] == list(range(8, 15))
assert PROTECTED[PPT] == INVENTORY["ppt_md5"]

forest_table = table(["SGB", "调整后丰度比", "95% CI", "BH q"], [
    [r.SGB_id, f(r.abundance_ratio), f"{f(r.ratio_CI_low)}–{f(r.ratio_CI_high)}", q(r.P_BH)] for r in forest.itertuples()
])
grade_table = table(["SGB", "等级每增加1级的调整后丰度比", "95% CI", "BH q"], [
    [r.SGB_id, f(r.abundance_ratio), f"{f(r.ratio_CI_low)}–{f(r.ratio_CI_high)}", q(r.P_BH)] for r in grade.itertuples()
])
candidate_table = table(["SGB", "脂肪变组方向", "LDM频率尺度差值（百分点）", "完整595项检验的BH q", "层级"], [
    [r.SGB_id, "较低" if r.beta_percentage_points < 0 else "较高", f(r.beta_percentage_points, 4), f(r.P_BH, 5), "q<0.05" if r.P_BH < .05 else "0.05≤q<0.10，探索性"] for r in candidates.itertuples()
])
phylum_table = table(["树枝显示分组", "SGB数", "占379个SGB的比例"], [
    [r.branch_colour_group, r.SGBs, f"{f(r.percentage_of_significant_SGBs, 2)}%"] for r in phylum.itertuples()
])
traits = [("BMI18", "BMI", "kg/m²"), ("WAIST18", "腰围", "cm"), ("syst18", "收缩压", "mmHg"),
          ("GLUK18", "空腹血糖", "mmol/L"), ("HbA1c18", "HbA1c", "mmol/mol"), ("TRIGLY18", "甘油三酯", "mmol/L"),
          ("HDLKOL18", "HDL胆固醇", "mmol/L"), ("ALAT18", "ALT", "U/L"), ("ASAT18", "AST", "U/L")]
zh = {code: label for code, label, unit in traits}
trait_table = table(["代谢指标", "原始单位", "每项相关分析实际N"], [
    [label, unit, int(host_summary.loc[(host_summary["mode"] == "Source-style") & (host_summary.phenotype == code), "N_max"].iloc[0])] for code, label, unit in traits
])
links = host.loc[(host["mode"] == "Source-style") & (host.P_BH < .05)].sort_values("P_BH")
assert len(links) == 13
link_table = table(["SGB", "指标", "偏Spearman ρ", "近似95% CI", "BH q", "N"], [
    [r.feature_id.replace("t__", ""), zh[r.phenotype], f(r.rho), f"{f(r.CI_low)}–{f(r.CI_high)}", q(r.P_BH), int(r.N)] for r in links.itertuples()
])
conditional = host.loc[(host["mode"] == "Steatosis-conditional") & (host.P_BH < .10)].sort_values("P_BH")
assert len(conditional) == 2
conditional_table = table(["SGB", "指标", "额外调整脂肪变后的ρ", "近似95% CI", "BH q"], [
    [r.feature_id.replace("t__", ""), zh[r.phenotype], f(r.rho), f"{f(r.CI_low)}–{f(r.CI_high)}", q(r.P_BH)] for r in conditional.itertuples()
])
path_table = table(["图中顺序", "通路ID", "图中主要分类", "成人Model 4丰度比", "BH q"], [
    [int(r.display_rank), r.feature_id, r.display_major, f(r.abundance_ratio), q(r.P_BH)] for r in paths.itertuples()
])

doc = f"""# PPT中Fig. S1起各图的方法与基本数据记录

记录日期：2026-11-05。

来源：`{PPT.as_posix()}`。本记录按当前PPT的实际顺序，覆盖第8–14页的7张图片：Fig. S1、Fig. S2、Fig. S3、Fig. S4、Fig. 5、Fig. 6，以及末页尚未编号的功能通路图。前面的队列流程、Table 1/2和Fig. 2–4仅用于交代基本研究背景，不另作逐图展开。

本记录由PPT文字、嵌入图片及对应的已保存结果表核对整理。PPT未修改，未重新拟合模型、重新置换或重新计算BH校正。以下比值、区间及q值均来自已有结果，表中仅调整显示精度。统计记录遵循用户指定的Althouse等（2021）AHA统计报告指南，具体模型、筛选和阈值属于本项目的实施选择。

## 一、研究对象、分析人群和共同方法

### 1. 基本研究情况

研究使用Young Finns Study项目中的临床、超声肝脏检查与粪便宏基因组数据。这里的疾病状态是“超声检出的肝脂肪变”，不是根据代谢条件另外构建的疾病分类。G0、G1、G2为从年长到年轻的三代标签。基线家系要求至少观察到两代成员，不要求每个家系三代全部齐全。

主要家系标签（primary family）用于处理聚类和家庭内比较。标签本身并不确认具体的一度亲属关系或共同居住情况。后面的“家系内”均按该标签定义。

### 2. 各分析使用的人数

| 分析人群 | 人数 | 脂肪变 | 无脂肪变 | 家系数 | 主要用途 |
| --- | --- | --- | --- | --- | --- |
| 两代以上临床队列 | 5,660 | 1,054 | 4,606 | 1,693 | 前部临床描述背景 |
| 两代以上微生物组基线队列 | 3,217 | 645 | 2,572 | 1,104 | Fig. S2检出率与平均相对丰度等描述轨道 |
| 全龄Model 3完全案例 | 3,193 | 640 | 2,553 | 1,104 | Fig. S1–S4的主分析；末页功能图a及c |
| 成人Model 3/4共同完全案例 | 2,040 | 447 | 1,593 | 1,030 | 成人敏感性分析；末页功能图b的成人列 |
| 家系内成人发现人群 | 731 | 320 | 411 | 289 | Fig. 5；候选菌来自该人群的Model 4 |
| 每家一名成人的宿主相关人群 | 最多1,104 | 本记录未另行统计分组病例数 | — | 1,104 | Fig. 6；各代谢指标使用各自完整测量人数 |

人数差异来自分析目的、年龄限制、协变量完整性、家系内是否同时存在两种脂肪变状态等条件，不能把各图的N互换使用。3,217人到全龄Model 3的3,193人减少24人，病例从645人减少到640人；本组模型采用完全案例，未使用多重插补。

全龄Model 3人群的代际构成为：G0 1,121人（303例），G1 1,019人（291例），G2 1,053人（46例）；包括421名未成年人。成人共同样本为G0 713人、G1 834人、G2 493人。

### 3. 模型名称与协变量

本范围主要涉及Model 3和Model 4。代码中的B1对应展示名称Model 3，B2对应Model 4；“成人Model 3”与“全龄Model 3”使用同一调整层级，但分析人群不同。

- Model 3：代际、各代内年龄的自然样条（通常3个自由度）、性别、建库/制备批次、测序池、采样季节的正弦和余弦项、采便前3个月抗生素使用。
- Model 4：在Model 3基础上加入吸烟、饮酒、饮食评分、体力活动、对数转换后的总能量摄入和连续BMI。成人Model 3与Model 4使用相同完整样本。
- 年龄、技术和生活方式变量按既定模型纳入；常量或完全线性依赖项按设计矩阵处理，不按单因素P值筛选。

MaAsLin2疾病关联模型采用主要家系随机截距；Fig. 5的LDM则使用家系固定指标和家系内置换。两者处理家系结构的方式不同。

### 4. 相对丰度、丰度比和多重比较

微生物组基线保留3,949个非零SGB。SGB为物种水平基因组分箱，不等同于菌株。原始丰度在完整微生物谱分母上转换为相对丰度；筛选或只展示部分菌后不重新归一化。

MaAsLin2使用官方1.24.1版的LOG线性混合模型，预先相对化后设normalization=NONE。对每个特征，在相应模型的分析样本中，将零值替换为该特征最小正值的一半，正值不变，然后取log2。模型可概括为：

`log2(零值替换后的相对丰度) = 截距 + β×脂肪变状态 + 协变量项 + 家系随机截距 + 误差`。

调整后丰度比为`2^β`，95% CI由log2尺度区间取2的幂得到。比值小于1表示脂肪变组较低，大于1表示较高。它描述零值替换后的相对丰度尺度，不能直接解释为绝对菌量、疾病OR或风险比。

BH q值来自相应模型的完整合格检验集合。只展示16个菌、10个候选菌或81条通路时，不在展示子集内重新计算BH。q<0.05用于常规发现；Fig. 5扩展到q<0.10的部分为探索性候选。未达到阈值表示证据未达到所定标准，不能据此断言没有关联。

## 二、Fig. S1：总体人群的属/基因组属水平脂肪变关联

PPT位置：第8页。原始图：模块39的`Figure39_Main_All_genus_volcano`。

研究问题：总体人群中，哪些属或基因组属组的相对丰度与超声脂肪变有关？

基本数据：全龄Model 3的3,193人，其中640例；共1,104个家系。由3,949个基线SGB中具有明确单一属归属的3,572个SGB聚合成2,295个属/基因组属组；377个归属模糊或未分类SGB不被强行分配到某一属。每个组要求在本模型样本中的检出率≥10%、方差非零，最终检验364组。

计算方法：同一受试者中属于同一明确属组的SGB相对丰度相加，保留原完整谱分母。逐属运行Model 3的MaAsLin2 LOG混合模型，家系随机截距；对完整364项检验作BH校正。

图形含义：每个点是一组属/基因组属，横轴为调整后丰度比（对数坐标），纵轴为−log10（原始P）。蓝色为脂肪变组较低且BH q<0.05，橙色为较高且BH q<0.05，灰色为未达到该阈值。点不是单个受试者。

基本结果：364组中265组达到q<0.05，包括250组较低、15组较高；其余99组未达到阈值。265个发现中61组具有命名属标签，204组为GGB标签。横向虚线把BH界限换算到原始P尺度，原图采用`0.05×265/364≈0.03640`，不是普通P=0.05线。

记录时应明确：本图是属/基因组属层级，不能与后面379个SGB的发现数相加或当成相同分类单位。

## 三、Fig. S2：脂肪变相关SGB在既有参考树上的分布

PPT位置：第9页。原始图：模块40的`Figure40_SGB_circular_tree_short_lines_20261103`。

研究问题：与脂肪变有关的SGB分布在哪些分类分支，它们的检出率、平均丰度和关联方向如何？

基本数据：全龄Model 3在587个合格SGB中发现379个BH q<0.05的SGB，其中343个在脂肪变组较低、36个较高。关联模型用3,193人（640例，1,104个家系）；检出率和平均相对丰度用基线3,217人。587个检验来自“任一基线代际检出率≥10%”的692个候选并集，再按实际模型样本的总体检出率≥10%等条件筛选。

计算方法：疾病效应与Fig. S3一样来自既有Model 3 MaAsLin2结果，BH保留完整587项检验。树使用原微生物组RDS中既有参考树，裁剪到全部379个发现SGB；保留其参考相对枝长。没有重新推断一棵研究人群的系统发育树。现有资料未记录该参考树的精确数据库版本和建树算法。

从内向外的环及外侧柱：

1. 检出率：`检出相对丰度>0的人数 / 3,217 ×100%`，蓝色。
2. 平均相对丰度：对3,217人的相对丰度求均值，包含未检出的零值；橙色强度用于显示。
3. 调整后脂肪变丰度比：`2^β`，蓝色小于1、白色等于1、红色大于1；颜色在log2尺度以1为中心。
4. 外侧柱：`−log10（原始P）`，不是BH q、基因型数量或菌株数。

树枝颜色的分组数量如下，Archaea和Other phyla为显示用宽分类组：

{phylum_table}

46个SGB具有唯一明确的物种名，显示斜体物种名及SGB编号；其余333个只保留SGB编号，其中87个存在冲突物种注释，246个为未命名或不完整注释。不根据名称相似性强行补物种名。

本图展示参考树上的关联分布，不提供菌株传播、传播方向或宿主遗传力的证据。

## 四、Fig. S3：总体人群中16个代表SGB的调整后丰度比

PPT位置：第10页。原始图：模块25的`Figure25_Overall_abundance_ratios`。

研究问题：展示总体Model 3中部分证据较强SGB的效应大小和不确定性。

基本数据：3,193人、640例脂肪变、1,104个家系。完整发现仍为587项检验、379个q<0.05；图中仅按BH q、原始P和编号的既定顺序选择16个有效发现展示，其中12个方向较低、4个较高。

计算方法：复用全龄Model 3的官方MaAsLin2 LOG混合模型结果及家系随机截距。点为`2^β`，横线为95% CI；参考线为丰度比1。16个菌的显示不定义一个新的多重检验集合，也不是独立验证集。

图中16个SGB的基本数值：

{forest_table}

图注中模型统一记为Model 3。年龄应按实际拟合写作“各代内年龄样条”，季节实际由周期正弦/余弦项表示。比值反映零值替换后的相对丰度，不表示绝对菌量的倍数变化。

## 五、Fig. S4：不同超声脂肪变等级下6个SGB的丰度分布

PPT位置：第11页。嵌入图片匹配模块33历史运行中统一模型标签前的`Figure33_Grade_abundance_profiles`；图片内部仍写“Ordinal B1”。本记录统一称“有序Model 3”。11月1日完整重算的保存核对表确认样本、筛选、模型状态、效应、P和BH q与该旧分析一致。

研究问题：SGB丰度是否与观察到的超声脂肪变等级有关？

基本数据与编码：全龄Model 3共3,193人，无脂肪变2,553人、轻度602人、明显38人。原超声编码1/2/3分别映射为No/Mild/Clear，并记为0/1/2；二分类脂肪变定义与“等级>0”一致。成人共同样本另有2,040人，无脂肪变1,593人、轻度421人、明显26人。

计算方法：以0/1/2等级分数作为MaAsLin2 LOG混合模型中的数值预测变量，调整Model 3协变量并设家系随机截距。即把共同模型中的二分类脂肪变项替换为等级分数；`2^β`表示等级每增加1级的相对丰度比。此编码对相邻等级在log2丰度上的变化作共同斜率假设，不能视为实测肝脂含量的等距增加。

筛选条件为：相对丰度至少0.0002（即0.02%）的检测在分析样本中达到至少10%，再检查变异等条件。全龄模型共348个合格SGB，完整348项有序检验作BH校正，226个达到q<0.05（211个负向、15个正向），246个达到q<0.10。不能把本图的348项检验与Fig. S2/S3的587项检验混为同一筛选集合。

图中只展示按有序Model 3的BH q、原始P和编号排序的前6个SGB。箱线图为原始相对丰度百分比的描述分布，箱体是中位数及四分位范围，灰点是受试者；q来自调整后的模型，不是对图中原始箱体重新做的组间检验。

{grade_table}

项目另保存三水平分类暴露模型（轻度与无脂肪变、明显与无脂肪变）；本页标注的q来自有序Model 3，不应替换为某一分类对比的q。明显组只有38人，解释这一级分布时须保留其样本量信息。这里分析横断面的观察等级，不表示随访中的病情进展。

## 六、Fig. 5：家系内脂肪变相关候选SGB的丰度分布

PPT位置：第12页。原始图：模块50的`Figure50_Family_abundance_profiles_restored_v3`，候选来源于模块45的家系内重新发现。

研究问题：在同一家系同时存在有、无脂肪变成人的情况下，哪些SGB仍有家系内关联，候选菌在两种状态下的家庭均值如何分布？

基本数据：成人Model 4完整样本731人，其中320例脂肪变、411例无脂肪变，来自289个家系。G0/G1/G2分别314/246/171人。每个保留家系在完整分析样本中均有两种脂肪变状态；家庭大小不要求为2，所以不能把731人表述为固定数量的一对一病例对照配对。

发现方法：从全部SGB重新筛选，不继承总体分析的379个阳性名单。按本分析人群总体检出率≥10%且携带者至少5人等条件，成人主模型检验595个SGB。使用官方LDM 6.0.1，加入家系指标与Model 4个体协变量，以频率尺度和反正弦平方根（AST）尺度的omnibus检验推断；置换仅在家系内进行，最多10,000次，达到100次超越后可顺序停止。P为Monte Carlo/顺序置换估计，最小加一P分辨率为1/10,001。常规BH基于全部595个已报告omnibus P，不能改用软件另一种q值替代。

常规q<0.05有2个SGB。查看q<0.05结果后，项目将候选扩展到q<0.10用于后续刻画，最终10个，其中8个属于0.05≤q<0.10的探索性层级；这不是预先固定的10个发现，也不是独立验证。

描述图计算方法：对每个SGB，分别在同一家系的无脂肪变成员和脂肪变成员中计算平均相对丰度百分比。每个家系贡献两个状态均值，灰线连接同一家系的两个均值；灰点代表家系状态均值。箱体、中位线和1.5倍IQR须线总结这些家系均值，蓝色为无脂肪变，橙色为脂肪变。伪对数坐标用于同时容纳零值和较大值，零值被保留。

这张图完整使用731人、289家系的原发现矩阵。用于绘图的家庭均值不替代LDM中的个体观测；图上的q是个体层面的调整后LDM结果。本次恢复箱线图后没有额外运行配对检验或bootstrap。

{candidate_table}

上表频率尺度差值以百分点表示，是调整后的描述性方向/幅度信息；omnibus P联合两个尺度，不对应一个唯一的效应单位。不能把这些差值与Fig. S3的LOG丰度比直接比较大小。10个候选在频率尺度上9个较低、1个较高；只有Yanshouia hominis（SGB53461）具有该候选表中的明确物种名，其余保留SGB标签及未解决的原注释。

## 七、Fig. 6：候选SGB与9项代谢指标的偏Spearman相关

PPT位置：第13页。原始图：模块50的`Figure50_Candidate_host_phenotypes_v2`。

研究问题：上述10个候选SGB与肥胖、血糖、血脂及肝酶等测量指标是否相关？加入脂肪变状态后，相关性证据如何变化？

基本数据：从3,217人的微生物组基线中，先保留年龄、性别、BMI完整且BMI有效的成人2,792人。按家系及样本编号固定排序后，用种子20261104在每个家系随机选择1名成人，得到1,104人/1,104家系。选择在查看菌丰度、特定代谢指标缺失或关联结果之前完成。它与发现分析来自重叠队列，不是外部或独立验证样本。

每个代谢指标按实际有效测量作完全案例分析，未插补；实际人数如下：

{trait_table}

计算方法：菌相对丰度包含零值，先对菌丰度、代谢指标及控制变量取秩，重复值取平均秩。分别将菌秩和指标秩对控制变量秩回归，两个残差之间的Pearson相关系数即偏Spearman ρ。原结果中的P采用双侧t近似，区间采用近似Fisher-z方法。

- a图：调整年龄、性别、BMI；分析BMI自身时只调整年龄和性别。
- b图：在a图对应调整项基础上，额外加入超声脂肪变状态。
- 背景筛选：在固定的1,104人样本中，检出率严格>20%的374个SGB参与检验。每种调整方式的BH集合是374×9=3,366项，不是只对展示的10×9=90项校正。不能估计的预定检验在校正内部按P=1处理，展示统计保持缺失。
- 显示：蓝色为负相关、红色为正相关；星号表示q<0.05，加号表示0.05≤q<0.10，叉号表示未检验或不可估计。两个面板使用同一相关系数色标。

SGB53461（Yanshouia hominis）在该样本的检出率约19.66%，未达到严格>20%的门槛，所以9项相关均标记未检验。每个面板显示90格，其中81格可估计。

a图候选菌中有13条q<0.05关联；q<0.10累计27条。13条常规阈值关联如下，区间为原保存结果中的近似95% CI：

{link_table}

b图额外调整脂肪变后，没有候选关联达到q<0.05，2条达到0.05≤q<0.10：

{conditional_table}

这些对照说明结果对额外条件调整的敏感性。不能把阈值变化直接解释为关联不存在、脂肪变的中介作用或因果路径。偏相关反映横断面关联。

## 八、末页未编号图片：功能通路关联、分类与来源菌

PPT位置：第14页。该页目前只有图片、没有图题或图注。原始图精确匹配模块44的`Figure44_Pathway_function_and_sources_compact`。本记录暂称“功能通路与来源菌图”，不擅自设定正式Fig.编号。

研究问题：宏基因组推断的功能通路丰度与脂肪变有何关联？成人充分调整后保留哪些通路？这些通路的已观测HUMAnN分层信号由哪些分类来源贡献？

### a. 全部功能通路的总体关联火山图

基本数据：全龄Model 3共3,193人、640例脂肪变。原始供应数据含620条未分层community通路。先用全部620条供应通路固定每个人的分母，再根据≥10%非零检出率和非零方差等筛选，得到420条合格通路；不能先挑通路再重新归一化。

逐通路采用MaAsLin2 LOG混合模型、主要家系随机截距及Model 3调整集，BH覆盖完整420项。420项中246条q<0.05：87条在脂肪变组较低，159条较高。图a展示全部有效检验，横轴为调整后丰度比的对数坐标，纵轴为−log10（原始P）；蓝/橙/灰分别表示较低的BH发现、较高的BH发现、未达到阈值。横向虚线为BH发现中的最大原始P，约0.02667，不是普通P=0.05阈值。

### b. 成人Model 4的81条通路及其多模型效应热图

成人Model 3与Model 4使用相同2,040人（447例，1,030家系），分别检验413条合格通路。成人Model 3有174条q<0.05，成人Model 4有81条q<0.05。图b完整展示这81条，47条在脂肪变组较低、34条较高；其中79条也属于全龄Model 3的246条发现。没有另设“前10条”或“前30条”展示限制。

六列依次为全龄Model 3、成人Model 3、成人Model 4、G0、G1、G2。最后三列来自同一个全龄Model 3中的代际特异脂肪变系数，仍是全龄分析，不是将成人发现人群再分别重拟合三代。

热图颜色表示丰度比，在log2尺度以比值1为白色中心，蓝色较低、橙色较高。星号来自对应完整效应检验集合的原BH q<0.05。总体模型各自校正完整420或413项；代际特异效应在完整420×3=1,260项集合内校正。某一代出现星号不等同于代际交互作用显著。项目另保存420项联合两自由度代际差异检验，94项q<0.05，但不能把该计数作为图b星号的含义。

81行用有来源记录的历史MetaCyc谱系组织，按主要类别、子类和更深类别分组；最深同类内根据原六列log2系数进行完整连接聚类，不按行标准化、不把缺失填0。b/c保持完全相同的行顺序。用于显示的主类为生物合成60条、降解10条、能量代谢6条、未匹配5条。该分类和排序是描述性组织，不是通路富集检验。

### c. 无脂肪变与脂肪变组的通路来源比例

本面板使用原全龄3,193人的实际HUMAnN物种分层通路量，按其原注释聚合到属。该人群不同于用于筛选81条通路的成人Model 4人群。

对通路p和疾病组g，属t的来源比例为：

`来源比例(g,p,t) = 100 × 该组该属分层通路量的均值 / 该组该通路全部分类来源分层量均值之和`。

分母包括所有其余命名属和未分类来源；不使用community总通路量作为分母，因为HUMAnN的community量不要求等于物种分层量之和。图中比例是未调整的描述量，没有给来源组间比较附加P值。

按跨全部81条通路合并的来源份额选择前8个命名属，不按脂肪变组间差异挑选：Bacteroides、Faecalibacterium、Escherichia、Blautia、Alistipes、Methanobrevibacter、Parabacteroides、Prevotella。其余保留为Other named genera和Unclassified；两组各10个来源列。蓝色强度表示来源百分比。

功能结果来自DNA宏基因组通路丰度，描述功能潜力，不能直接写成基因表达、代谢物浓度或实测通量。来源图也不把属级信号强行归到某个SGB或菌株。历史MetaCyc分类存在5条未匹配通路，保留Unmapped块；数据库标签本身不确认样本中生物来源或宿主组织来源。

## 九、当前PPT的标签与版本记录

1. Fig. S1是属/基因组属分析；Fig. S2–S4和Fig. 5/6的行是SGB，分类层级不同。
2. Fig. S2页面的参考说明仍出现“GEE method”和“Genus model”，与这张SGB树的脂肪变关联模型口径不一致。该树的关联估计实际来自MaAsLin2 Model 3；本记录按结果和方法文件说明，未修改该页。
3. Fig. S3当前文字已不再写B1，但也未写出模型名称和样本量。本记录补记为全龄Model 3、3,193人、640例。
4. Fig. S4嵌入图片仍写“Ordinal B1”，是统一模型显示标签之前的图。B1应对应Model 3；已保存的新旧结果核对表确认分析数值一致。
5. Fig. 5的10个候选不等于10个q<0.05发现。2个达到q<0.05，另外8个是查看结果后采用q<0.10扩展的探索性候选。
6. Fig. 6每家选择1人，最多1,104人；各指标的N为1,095–1,104。叉号尤其包括SGB53461未通过>20%检出率门槛，不表示关联为零。
7. 末页暂未编号、无图注，应在后续正式排版时另行确定图号；本次只整理文本。

## 十、可直接用于个人记录的简要流程

先在两代以上微生物组队列中形成全龄Model 3的3,193人分析集，进行属/基因组属关联筛查（Fig. S1）和SGB关联筛查，并以参考树及16菌森林图呈现结果（Fig. S2/S3）。随后用原超声三级编码分析严重程度关联（Fig. S4）。另在有、无脂肪变成人并存的家系中进行不依赖前述阳性名单的全SGB家系内LDM发现，得到2个q<0.05发现，扩展为10个q<0.10探索性候选，展示家庭状态均值（Fig. 5）。再在每家1名成人的固定样本中，进行这10个候选与9项宿主指标的偏Spearman刻画，BH覆盖全部背景SGB与9指标（Fig. 6）。功能通路作为单独的数据层，使用完整供应通路分母、全量模型和成人Model 4的81条发现，展示多模型效应及实际分类来源（末页）。

这些都是同一研究队列中的关联分析或后续刻画；总体关联、家系内关联、宿主相关和功能来源的目标与效应尺度不同。保留各自样本量、家系处理、完整检验集合、效应量及区间，避免串用分母或把同队列结果写成独立验证。

## 附录A：末页b/c全部81条通路（按实际显示顺序）

下表比值及q来自成人Model 4；完整六列效应、95% CI、原始P、q和来源矩阵可在附录B的模块44结果表中查阅。

{path_table}

## 附录B：核对来源与复现位置

主要统计报告指南：`E:/YFS/workspace/STATISTICAL_ANALYSIS_GUIDE.md`。原文：Althouse AD et al. Circulation. 2021;144:e70–e91. DOI:10.1161/CIRCULATIONAHA.121.055393。

PPT源文件MD5：`{PROTECTED[PPT]}`。已提取的页码、图注、图片哈希及匹配路径保存于：`{(AUDIT / 'ppt_source_inventory.json').as_posix()}`。

直接用于本记录表格和人数核对的已保存结果如下：
"""

for p in dict.fromkeys(SOURCES):
    doc += f"\n- [{p.name}]({p.as_posix()})"

doc += "\n\n其他关键复现资料：\n"
extra = [
    adopted("25_") / "reports/Methods_and_references_EN.md",
    adopted("38_") / "reports/Methods_and_references_EN.md",
    adopted("39_") / "tables/Figure39_All_tested_genus_volcano_source_EN.csv",
    adopted("40_") / "reports/Figure40_legend_EN.md",
    adopted("40_") / "tables/Figure40_Three_ring_source_EN.csv",
    adopted("33_") / "reports/Fresh_reanalysis_summary_zh.md",
    adopted("33_") / "tables/B_All_grade_associations_EN.csv",
    adopted("45_") / "reports/Methods_results_and_limits_EN.md",
    adopted("45_") / "tables/All_tested_features_EN.csv",
    adopted("50_") / "reports/Methods_results_and_limits_EN.md",
    adopted("50_") / "tables/Family_state_descriptive_summary_EN.csv",
    adopted("50_") / "tables/Host_tests_EN.csv",
    adopted("44_") / "reports/Methods_results_and_figure_legend_EN.md",
    adopted("44_") / "tables/Figure44b_Adult_Model4_all81_effects_EN.csv",
    adopted("44_") / "tables/Figure44c_observed_source_matrix_EN.csv",
    ROOT / "R/reanalysis_20261024.R",
    ROOT / "R/family_candidate_followup.R",
]
for p in extra:
    assert p.exists(), str(p)
    doc += f"\n- [{p.name}]({p.as_posix()})"
doc += "\n\n记录构建脚本：`" + Path(__file__).as_posix() + "`。本脚本只读取已保存结果并组织文字，不运行统计模型。\n"

# Record the resolved remaining image versions; the PPT itself remains untouched.
INVENTORY["slides"][3]["images"][0]["pixel_match"] = str(ROOT / "outputs/33_Steatosis_severity/history/ultrasound_grade_Mei_style_v1_20261028_20261028_164128/audit/model_labels_20261030/before/Figure33_Grade_abundance_profiles.tiff")
INVENTORY["slides"][6]["images"][0]["matches"] = [str(adopted("44_") / "figures/Figure44_Pathway_function_and_sources_compact.tiff")]
(AUDIT / "ppt_source_inventory.json").write_text(json.dumps(INVENTORY, ensure_ascii=False, indent=2), encoding="utf-8")

destination = ROOT / "PPT_FigS1_onward_methods_data_20261105.md"
destination.write_text(doc, encoding="utf-8")
plain = re.sub(r"^#{1,6}\s+", "", doc, flags=re.M)
plain = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r"\1（\2）", plain).replace("`", "")
destination.with_suffix(".txt").write_text(plain, encoding="utf-8-sig")
for p, md5 in PROTECTED.items():
    assert hashlib.md5(p.read_bytes()).hexdigest() == md5, f"Protected source changed: {p}"
validation = {"status": "PASS", "slides": list(range(8, 15)), "figures": 7,
              "ppt_unchanged": True, "source_tables_unchanged": True,
              "models_refitted": False, "source_fingerprints": {str(p): h for p, h in PROTECTED.items()},
              "outputs": [str(destination), str(destination.with_suffix('.txt'))]}
(AUDIT / "text_record_validation.json").write_text(json.dumps(validation, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps({"status": "PASS", "characters": len(doc), "lines": len(doc.splitlines()),
                  "outputs": validation["outputs"]}, ensure_ascii=False, indent=2))
