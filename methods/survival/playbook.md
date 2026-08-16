# ICI biomarker survival playbook / ICI 生物标志物生存分析手册

**Scope / 范围.** Time-to-event analysis of a prespecified transcriptomic marker
(here: **TACSTD2** / TROP2 and **CLDN4**) in immune-checkpoint inhibitor (ICI)
cohorts. English and Chinese are kept in lockstep. Every numeric result in
§13 was computed from public GEO files by `demo/02_analysis.py`. Nothing else
in this document is a TACSTD2 or CLDN4 effect size.

**范围。** 预先指定的转录组标志物（本文：**TACSTD2**/TROP2 与 **CLDN4**）在免疫检查点
抑制剂（ICI）队列中的时间–事件分析。中英文同步。第 13 节所有数字均由
`demo/02_analysis.py` 从公开 GEO 文件当场计算。本文其余部分不是 TACSTD2 或 CLDN4
的效应量。

Outputs live only under `methods/survival/`.

---

## 1. Write the estimand before touching the matrix / 先写清研究问题

A survival analysis that starts with a Kaplan–Meier plot of a data-driven cut
is already biased. Freeze, in writing:

1. **Population:** histology, line of therapy, drug/regimen, biopsy timing
   (pretreatment vs on-treatment).
2. **Endpoint:** OS, PFS, or TTR, with the exact event and censoring rules
   (see §2). Do not switch after seeing *p*.
3. **Marker:** gene identifier, transformation (`log2(TPM+1)`), and whether
   the primary model is continuous. Dichotomisation is a sensitivity analysis.
4. **Comparators:** a clinical-only model (PD-L1, TMB, line, histology) and,
   if space allows, one known-direction immune gene (this demo uses CD8A).
5. **Validation:** a locked external cohort whose cut, transformation, and
   direction are not re-tuned.

生存分析如果从“最优切点”的 Kaplan–Meier 图开始，就已经有偏。先书面冻结：

1. **目标人群：** 组织学类型、治疗线次、药物/方案、活检时点（治疗前 vs 治疗中）。
2. **终点：** OS、PFS 或 TTR，以及精确的事件与删失规则（见第 2 节）。看到 *p* 之后不得更换。
3. **标志物：** 基因标识、变换（`log2(TPM+1)`）、以及主模型是否为连续变量。二分只作敏感性分析。
4. **对照：** 仅临床模型（PD-L1、TMB、线次、组织学），以及（如有余力）一个方向已知的免疫基因
   （本示例用 CD8A）。
5. **验证：** 锁定的外部队列，切点、变换和方向均不再重调。

A single-arm ICI series cannot tell a **predictive** marker (treatment
interaction) from a **prognostic** marker (outcome regardless of drug). That
distinction needs a randomised comparison such as OAK/POPLAR (§10).

单臂 ICI 序列无法区分**预测**标志物（治疗交互）和**预后**标志物（与药物无关的结局）。
该区分需要 OAK/POPLAR 这类随机对照（第 10 节）。

---

## 2. Time-to-event definitions / 时间–事件定义

| Endpoint | Clock starts | Event | Typical censoring | ICI-specific traps |
|---|---|---|---|---|
| **OS** | ICI start (or randomisation) | Death from any cause | Alive at last contact | Cleanest endpoint. Needs long follow-up. Post-progression therapy dilutes a PFS signal but is part of the OS estimand. |
| **PFS** | ICI start (or randomisation) | RECIST progression **or** death | Alive and progression-free at last scan | Scan-interval bias; pseudoprogression; iRECIST vs RECIST; clinical progression without imaging. |
| **TTR / TTP** | ICI start | Response (TTR) or progression (TTP) | No response / no progression; death is **not** the event | Death is a competing risk. 1−KM that censors death overestimates incidence (§8). |

| 终点 | 计时起点 | 事件 | 典型删失 | ICI 特有陷阱 |
|---|---|---|---|---|
| **OS** | ICI 开始（或随机化） | 任何原因死亡 | 末次随访仍存活 | 最干净的终点。需要足够随访。进展后治疗会稀释 PFS 信号，但这属于 OS 的估计目标。 |
| **PFS** | ICI 开始（或随机化） | RECIST 进展 **或** 死亡 | 末次影像时仍无进展且存活 | 扫描间隔偏倚；假性进展；iRECIST vs RECIST；无影像的临床进展。 |
| **TTR / TTP** | ICI 开始 | 缓解（TTR）或进展（TTP） | 未缓解 / 未进展；死亡**不是**事件 | 死亡是竞争风险。把死亡当删失的 1−KM 会高估发生率（第 8 节）。 |

**Unit of time.** GEO depositors mix days and months. GSE135222 stores
`pfs.time` in **days**; GSE190265 stores `time_PFS` in **months**. Convert
explicitly (`days / (365.25/12)`) and write the unit on every axis. Do not
meta-analyse a hazard ratio whose time unit was never checked.

**时间单位。** GEO 提交者混用日与月。GSE135222 的 `pfs.time` 是**天**；GSE190265 的
`time_PFS` 是**月**。必须显式换算（`天数 / (365.25/12)`），并在每条轴上写明单位。
未核对时间单位的风险比不得做荟萃分析。

**Event indicator.** In both demo series, `1` = progression (or death if that
is how the depositor coded PFS) and `0` = censored. Confirm this against the
series matrix rather than assuming Kaplan–Meier convention. A swapped
indicator silently inverts every HR.

**事件指示。** 两个示例队列中 `1` = 进展（或提交者把死亡并入 PFS），`0` = 删失。
对照 series matrix 核实，不要假定。指示符对调会无声地反转所有 HR。

**What these two GEO series do not give you.** No OS. No cause of death. No
scan dates. No iRECIST. No PD-L1. No TMB. No treatment line in GSE190265.
You cannot reconstruct TTR.

**这两个 GEO 序列不能提供的。** 无 OS、无死亡原因、无扫描日期、无 iRECIST、无 PD-L1、
无 TMB、GSE190265 无治疗线次。无法重建 TTR。

---

## 3. Cox PH is the primary model / 主模型是 Cox 比例风险

Fit the marker as a **continuous** covariate, scaled to 1 SD in the
**training** cohort:

```text
h(t | z) = h0(t) exp(β z),    z = (x − mean_train) / sd_train
```

Report:

- HR per 1 SD with a 95% CI and a 1-df partial likelihood-ratio *p*;
- events per variable (EPV);
- Harrell *C* **and** Uno *C* (IPCW) — Harrell *C* is biased under heavy
  censoring;
- a Schoenfeld PH test (§4);
- the same model on a locked external cohort, using the **training** mean/SD.

将标志物作为**连续**协变量拟合，并按**训练**队列的 1 个标准差缩放：

```text
h(t | z) = h0(t) exp(β z),    z = (x − mean_train) / sd_train
```

报告：

- 每 1 SD 的 HR、95% CI，以及 1 自由度偏似然比 *p*；
- 每个变量的事件数（EPV）；
- Harrell *C* **和** Uno *C*（IPCW）——重度删失下 Harrell *C* 有偏；
- Schoenfeld 比例风险检验（第 4 节）；
- 在锁定外部队列上用**训练**均值/标准差重复同一模型。

**EPV rule of thumb.** Peduzzi / Vittinghoff: keep EPV ≥ 10 before adding
covariates. GSE135222 has 21 events → univariable only, or one extra
covariate at most. Age + sex already drops EPV to 7. The templates refuse a
fit below a configurable EPV floor.

**EPV 经验规则。** Peduzzi / Vittinghoff：加入协变量前保持 EPV ≥ 10。GSE135222 有
21 个事件 → 只做单变量，或最多再加一个协变量。年龄 + 性别已使 EPV 降到 7。
模板在 EPV 低于可配置门槛时拒绝拟合。

**Do not** screen 20 000 genes, pick the best Cox *p*, and then report that
*p*. That is the same crime as an optimal cut (§5), at a larger scale.

**禁止**筛 20 000 个基因、挑最小的 Cox *p* 再报告该 *p*。这与“最优切点”（第 5 节）
是同一类错误，只是规模更大。

Templates: `templates/python/cox_ph.py`, `templates/R/cox_ph.R`.

---

## 4. Proportional-hazards assumptions / 比例风险假定

Cox reports **one** HR. That number is a weighted average of time-varying
effects. ICI curves often separate late (or cross). Check, in this order:

1. **Schoenfeld residuals** vs time (`cox.zph` / `proportional_hazard_test`),
   with more than one time transform (KM, rank, identity). A non-significant
   test in *n* = 27 does **not** prove PH — it is underpowered.
2. **Piecewise HR** split at a prespecified time (median event time, or 3 / 6
   months), not at a time chosen to maximise the contrast.
3. **Log-log survival plots** by a prespecified split (median, not optimal).
4. If PH is implausible, report **RMST** to a prespecified τ inside the
   follow-up support, or a time-varying coefficient. Do not quote a single HR
   as if it were constant.

Cox 只报告**一个** HR。该数字是时变效应的加权平均。ICI 曲线常晚期分开（或交叉）。
按此顺序检查：

1. **Schoenfeld 残差**对时间（`cox.zph` / `proportional_hazard_test`），并用多种
   时间变换（KM、秩、恒等）。*n* = 27 时检验不显著**不能**证明 PH——只是效能不足。
2. **分段 HR**，在预先指定的时间切开（事件时间的中位数，或 3 / 6 个月），
   而不是选让对比最大的时间。
3. 按预先指定的切分（中位数，不是最优切点）画 **log-log 生存图**。
4. 若 PH 不可信，报告随访支持内预先指定 τ 的 **RMST**，或时变系数。
   不要把单个 HR 当作恒定值引用。

---

## 5. Median cut vs “optimal” cut / 中位数切点 vs “最优”切点

A median (or other prespecified quantile) split is legitimate as a
**secondary** display. A cut chosen to maximise the log-rank statistic is a
multiple-testing procedure. The naive *p* at that cut is not a type-I error.

中位数（或其他预先指定分位数）切分可作为**次要**展示。为最大化 log-rank 而选的切点
是多重检验。该切点上的朴素 *p* 不是 I 类错误。

Required companion numbers whenever an “optimal” cut is shown:

1. How many candidate cuts were searched (and the minimum arm fraction).
2. A **permutation** (or `maxstat` / Lausen–Schumacher) *p* for the maximised
   statistic.
3. A **bootstrap** of the selected percentile — if the cut wanders across
   half the distribution, it is not a threshold.
4. The same cut, frozen, on an external cohort. Re-optimising the cut in
   validation is not validation.

只要展示“最优”切点，就必须同时给出：

1. 搜索了多少候选切点（以及每臂最小比例）。
2. 对最大化统计量的**置换**（或 `maxstat` / Lausen–Schumacher）*p*。
3. 所选分位数的**自助法**——若切点在分布的一半范围内游荡，它就不是阈值。
4. 将该切点冻结后用于外部队列。在验证集上重新优化切点不是验证。

Templates: `templates/python/cutpoint.py`, `templates/R/cutpoint.R`.

---

## 6. Time-dependent ROC / 时依 ROC

A single AUC that treats “progressed vs not” as a binary label ignores
censoring and time. Use a **cumulative/dynamic** AUC with IPCW
(`timeROC`, `sksurv.metrics.cumulative_dynamic_auc`) at prespecified times
that still have events **and** patients at risk.

把“进展 vs 未进展”当二分类的单个 AUC 忽略了删失和时间。应在仍有事件**且**仍有
处于风险患者的预先指定时点，使用带 IPCW 的**累积/动态** AUC
（`timeROC`，`sksurv.metrics.cumulative_dynamic_auc`）。

**Orientation.** These functions expect a **risk** score (higher = earlier
event). TACSTD2 / CLDN4 may be protective, null, or harmful. Always pass the
Cox linear predictor, not the raw TPM. Passing a protective marker as a risk
score reports `1 − AUC` and will be misread as “worse than chance”.

**方向。** 这些函数需要**风险**分数（越高 = 事件越早）。TACSTD2 / CLDN4 可能保护、
无效或有害。始终传入 Cox 线性预测变量，而不是原始 TPM。把保护性标志物当风险分数
会报告 `1 − AUC`，并被误读为“差于随机”。

**Support.** IPCW is undefined after the last event. GSE135222’s last PFS
event is before 12 months, so a 12-month AUC cannot be computed — report the
drop, do not impute.

**支持区间。** 最后一个事件之后 IPCW 无定义。GSE135222 的末次 PFS 事件在 12 个月
之前，因此无法计算 12 个月 AUC——报告该缺失，不要填补。

Apparent AUC on the same *n* = 27 used to fit the Cox model is optimistic.
Nest the calculation inside cross-validation or apply a locked model
externally.

在用于拟合 Cox 的同一 *n* = 27 上得到的表观 AUC 偏乐观。应嵌套在交叉验证内，
或把锁定模型用于外部队列。

Templates: `templates/python/time_dependent_auc.py`, `templates/R/time_dependent_auc.R`.

---

## 7. DCB ≥ 6 months vs continuous PFS / 持久临床获益 vs 连续 PFS

Durable clinical benefit (DCB) is often coded as PFS ≥ 6 months (sometimes
≥ 12). It is a **coarsened** version of the same clock, not an independent
endpoint.

持久临床获益（DCB）常被编码为 PFS ≥ 6 个月（有时 ≥ 12）。它是同一时钟的**粗化**，
不是独立终点。

Problems:

1. **Guarantee-time bias.** Patients censored before 6 months cannot be
   labelled. Dropping them conditions on future follow-up. Report the
   unclassifiable count every time. In GSE135222 every censored patient was
   followed past 6 months, so the trap is empty — that fact still has to be
   shown, not assumed.
2. **Power.** Binary DCB discards event times. A Cox model on all patients
   uses every event.
3. **Landmark analysis** is the honest way to ask “given that a patient
   reached 6 months event-free, does the marker still associate with later
   PFS?”. Restart the clock at the landmark. Do not put DCB on the left-hand
   side of a logistic regression and call it survival.

问题：

1. **保证时间偏倚。** 6 个月前被删失的患者无法标记。丢掉他们等于以未来随访为条件。
   每次都要报告不可分类人数。GSE135222 中每例删失患者随访都超过 6 个月，陷阱是空的
   ——这一事实仍须展示，不能假定。
2. **效能。** 二分类 DCB 丢掉了事件时间。对全部患者的 Cox 使用每一个事件。
3. **路标分析**才是诚实的问法：“若患者在 6 个月时仍无事件，标志物是否仍与此后 PFS
   相关？”在路标处重新计时。不要把 DCB 放在 logistic 回归左侧并称之为生存分析。

DCB remains useful as a **secondary, prespecified** binary display, and as a
bridge to papers that only published DCB. It is not a substitute for PFS.

DCB 仍可用作**次要、预先指定**的二分类展示，以及与只发表了 DCB 的论文对接。
它不能代替 PFS。

---

## 8. Competing risks / 竞争风险

Use a competing-risks model only when a competing event is **coded**. Typical
ICI cases: death without documented progression (for TTP), treatment switch,
or an irAE that stops the drug when irAE-free survival is the endpoint.

仅在竞争事件**被编码**时才使用竞争风险模型。ICI 中的典型情形：无记录进展的死亡
（对 TTP）、换药，或当终点是无 irAE 生存时因 irAE 停药。

- **Aalen–Johansen** CIF for incidence.
- **Fine–Gray** (`cmprsk::crr`) for a covariate effect on the CIF. The
  subdistribution hazard is not a cause-specific hazard; say so.
- **Cause-specific Cox** if the scientific question is the instantaneous
  rate among those still at risk.

- **Aalen–Johansen** CIF 用于发生率。
- **Fine–Gray**（`cmprsk::crr`）用于协变量对 CIF 的效应。子分布风险不是原因特异
  风险；必须写明。
- 若科学问题是仍处于风险者的瞬时速率，用**原因特异 Cox**。

`1 − KM` that censors the competing event **overestimates** incidence. The
demo shows this on **simulated** data with a known truth (GEO files here have
no cause field). Templates refuse to run when only `{0, 1}` are present.

把竞争事件当删失的 `1 − KM` **高估**发生率。示例在**已知真相的模拟数据**上展示
这一点（此处 GEO 文件没有原因字段）。当只有 `{0, 1}` 时，模板拒绝运行。

Templates: `templates/python/competing_risks.py`, `templates/R/competing_risks.R`.

---

## 9. Small-*n* pitfalls / 小样本陷阱

GSE135222 (*n* = 27, 21 PFS events) and GSE190265 (*n* = 43, 35 PFS events)
are the largest **open GEO** NSCLC ICI series that deposit both a full
RNA-seq matrix **and** a usable PFS clock **and** TACSTD2/CLDN4. They are
still too small to claim a biomarker.

GSE135222（*n* = 27，21 个 PFS 事件）和 GSE190265（*n* = 43，35 个 PFS 事件）
是同时具备完整 RNA-seq 矩阵、可用 PFS 时钟以及 TACSTD2/CLDN4 的最大**公开 GEO**
NSCLC ICI 序列。它们仍然太小，不足以声称生物标志物。

| Trap | Why it bites here |
|---|---|
| Schoenfeld events for 80% power, two-sided α = 0.05, equal split | **65** events to detect HR = 0.50; **196** events to detect HR = 0.67 (Schoenfeld formula, computed in the demo). 21 and 35 events cannot rule in or out a modest effect. |
| Apparent *C* / AUC | In-sample concordance on 21 events is a coin-flip with a wide interval. |
| Multivariable Cox | Age + sex already violates EPV ≥ 10 on GSE135222. |
| Optimal cut | 13-ish candidate cuts, each a test. Naive *p* < 0.05 is expected under the null often enough to fool a figure legend. |
| Transfer of a z-score | Gene-wise mean/SD do not match across libraries. A frozen cut can put 10% or 90% of the next cohort in the “high” arm. |
| Multiple genes × multiple endpoints × multiple cuts | TACSTD2, CLDN4, their mean, DCB, PFS, median, optimal: that is already a family. Report the family, or correct it. |
| Purity / histology | Bulk TACSTD2 and CLDN4 track epithelial content. A “predictive” HR may be a LUSC-vs-LUAD or tumour-fraction HR. Neither GEO series has purity or a reliable histology field for adjustment. |
| No control arm | Any association is prognostic until a randomised comparison says otherwise. |

| 陷阱 | 为何在这里会咬人 |
|---|---|
| 80% 效能、双侧 α = 0.05、均等分组所需的 Schoenfeld 事件数 | 检出 HR = 0.50 需 **65** 个事件；检出 HR = 0.67 需 **196** 个事件（Schoenfeld 公式，示例中计算）。21 和 35 个事件既不能证实也不能排除中等效应。 |
| 表观 *C* / AUC | 21 个事件上的样本内一致性接近抛硬币，区间很宽。 |
| 多变量 Cox | 在 GSE135222 上加年龄 + 性别已违反 EPV ≥ 10。 |
| 最优切点 | 大约 13 个候选切点，每个都是一次检验。零假设下朴素 *p* < 0.05 的频率足以骗过图注。 |
| z 分数迁移 | 基因均值/标准差在文库间不匹配。冻结切点可能把下一队列的 10% 或 90% 放进“高”组。 |
| 多基因 × 多终点 × 多切点 | TACSTD2、CLDN4、它们的均值、DCB、PFS、中位数、最优切点：这已经是一个家族。要么报告家族，要么校正。 |
| 纯度 / 组织学 | bulk TACSTD2 与 CLDN4 追随上皮含量。“预测”HR 可能只是 LUSC-vs-LUAD 或肿瘤比例的 HR。两个 GEO 序列都没有纯度或可靠组织学字段可供校正。 |
| 无对照臂 | 在随机比较给出相反证据之前，任何关联都只是预后关联。 |

A null result on these series is **not** evidence of no effect. A positive
result on these series is **not** evidence of a usable test.

这些序列上的阴性结果**不是**无效应的证据。这些序列上的阳性结果**也不是**可用检测的证据。

---

## 10. Why OAK / POPLAR matter, and why they are not the demo / 为何 OAK/POPLAR 重要，以及为何不是本示例

|  | GSE135222 / GSE190265 | POPLAR + OAK |
|---|---|---|
| Design | Single-arm ICI | Randomised atezolizumab vs docetaxel (POPLAR NCT01903993, *n* = 287, Fehrenbacher *Lancet* 2016; OAK NCT02008227, *n* = 850, Rittmeyer *Lancet* 2017) |
| Primary endpoint | PFS only, investigator / depositor coded | OS (powered), PFS secondary |
| Events | 21 and 35 PFS events | Hundreds of OS events — enough for EPV, PH diagnostics, time-dependent ROC, and a treatment × marker interaction |
| Marker measurability | Open RNA-seq TPM including TACSTD2 and CLDN4 | RNA exists (EGA / Roche access, not GEO). This playbook does **not** analyse it and does **not** invent TACSTD2/CLDN4 HRs |
| What you can claim | Prognostic association with PFS, at best, with intervals that include the null | Predictive vs prognostic, if the interaction is prespecified and the expression data are actually in hand |

|  | GSE135222 / GSE190265 | POPLAR + OAK |
|---|---|---|
| 设计 | 单臂 ICI | 随机 atezolizumab vs 多西他赛（POPLAR NCT01903993，*n* = 287，Fehrenbacher *Lancet* 2016；OAK NCT02008227，*n* = 850，Rittmeyer *Lancet* 2017） |
| 主要终点 | 仅 PFS，研究者/提交者编码 | OS（按效能设计），PFS 为次要终点 |
| 事件数 | 21 与 35 个 PFS 事件 | 数百个 OS 事件——足够做 EPV、PH 诊断、时依 ROC，以及治疗 × 标志物交互 |
| 标志物可测性 | 开放 RNA-seq TPM，含 TACSTD2 与 CLDN4 | RNA 存在（EGA / Roche 权限，不在 GEO）。本手册**不**分析它，也**不**编造 TACSTD2/CLDN4 的 HR |
| 能声称什么 | 充其量是与 PFS 的预后关联，且区间常包含无效值 | 若交互预先指定且表达数据确实在手，可区分预测与预后 |

Trial-level OS results (atezolizumab vs docetaxel) are in those *Lancet*
papers. They are **not** TACSTD2 or CLDN4 results. Do not paste a published
OS HR next to a GEO forest plot as if it were the same estimand.

试验水平的 OS 结果（atezolizumab vs 多西他赛）见上述 *Lancet* 论文。它们**不是**
TACSTD2 或 CLDN4 的结果。不要把已发表的 OS HR 贴到 GEO 森林图旁边，假装是同一估计目标。

The correct use of this playbook on OAK/POPLAR, **once expression is
accessible**, is: lock the transformation and direction on a training slice
or on POPLAR, test on OAK, fit `treatment + z + treatment:z` for OS, and
report the interaction. Until that file is in the working directory, the
honest sentence is “not analysed here”.

一旦表达数据可及，本手册在 OAK/POPLAR 上的正确用法是：在训练切片或 POPLAR 上锁定
变换与方向，在 OAK 上测试，对 OS 拟合 `treatment + z + treatment:z`，并报告交互。
在该文件进入工作目录之前，诚实的句子是“此处未分析”。

---

## 11. Public OS / PFS cohorts checked for this playbook / 本手册核对过的公开 OS/PFS 队列

Checked on 2026-08-16 against GEO series matrices and supplementary files.
“Usable for TACSTD2/CLDN4 survival” means: (i) a time and an event column,
(ii) a gene-level matrix that actually contains both genes, (iii) ICI-treated
patients. No numbers below are taken from memory.

于 2026-08-16 对照 GEO series matrix 与补充文件核对。“可用于 TACSTD2/CLDN4 生存”
指：(i) 有时间和事件列，(ii) 基因矩阵确实包含这两个基因，(iii) 患者接受了 ICI。
下列数字均非凭记忆填写。

| Accession | ICI setting | Clock on GEO | TACSTD2 / CLDN4 | Verdict |
|---|---|---|---|---|
| **GSE135222** | NSCLC, anti–PD-1/PD-L1, *n* = 27 | PFS days + event; **no OS** | Both in the RSEM TPM matrix | **Demo discovery** |
| **GSE190265** | NSCLC, anti–PD-1, *n* = 43 | PFS months + event; **no OS** | Both in the TPM matrix | **Demo validation** |
| GSE136961 | NSCLC, anti–PD-1, *n* = 21 | PFS **and** OS deposited (18 PFS events, 14 OS events; parsed from `characteristics: survival`) | Series matrix has **0** expression rows; submitter lost raw data for 11/21 samples | Clock exists; genes are not measurable here |
| GSE93157 | Mixed MEL / NSCLC / HNSCC, nivo/pembro, *n* = 65 | PFS months + event (45 events); **no OS** | NanoString PanCancer 730 Immune: **neither gene present** | Clock exists; wrong panel |
| GSE126044 | NSCLC, anti–PD-1, *n* = 16 | Response only | RNA-Access; not used | No time-to-event |
| GSE166449 | Advanced lung, pre-treatment | Response assembly; no PFS/OS in the series matrix | Not used | No time-to-event on GEO |
| GSE207422 | Neoadjuvant anti–PD-1 + chemo NSCLC | Pathologic / sampling fields; not a PFS/OS clock | scRNA, not a bulk survival table | Wrong estimand |
| POPLAR / OAK | Randomised 2L NSCLC, atezolizumab vs docetaxel | OS + PFS in the trials | RNA not on GEO | Right design; data not open in this repo |
| TCGA-LUAD / LUSC | Mostly untreated primary resection | OS / PFI | Both genes present | **Not ICI.** Prognostic ≠ predictive |
| IMvigor210 | Urothelial, atezolizumab | OS public via `IMvigor210CoreBiologies` | Both genes present in that package | Wrong disease for a lung claim |

| 登录号 | ICI 设置 | GEO 上的时钟 | TACSTD2 / CLDN4 | 结论 |
|---|---|---|---|---|
| **GSE135222** | NSCLC，抗 PD-1/PD-L1，*n* = 27 | PFS 天 + 事件；**无 OS** | RSEM TPM 矩阵中两者都有 | **示例发现集** |
| **GSE190265** | NSCLC，抗 PD-1，*n* = 43 | PFS 月 + 事件；**无 OS** | TPM 矩阵中两者都有 | **示例验证集** |
| GSE136961 | NSCLC，抗 PD-1，*n* = 21 | 存有 PFS **和** OS（18 个 PFS 事件，14 个 OS 事件；从 `characteristics: survival` 解析） | series matrix **0** 行表达；提交者丢失 11/21 例原始数据 | 时钟在；基因不可测 |
| GSE93157 | 混合 MEL / NSCLC / HNSCC，nivo/pembro，*n* = 65 | PFS 月 + 事件（45 个事件）；**无 OS** | NanoString PanCancer 730 Immune：**两基因都不在** | 时钟在；面板不对 |
| GSE126044 | NSCLC，抗 PD-1，*n* = 16 | 仅缓解 | RNA-Access；未用 | 无时间–事件 |
| GSE166449 | 晚期肺，治疗前 | 缓解汇编；series matrix 无 PFS/OS | 未用 | GEO 上无时间–事件 |
| GSE207422 | 新辅助抗 PD-1 + 化疗 NSCLC | 病理/取样字段；不是 PFS/OS 时钟 | scRNA，不是 bulk 生存表 | 估计目标不对 |
| POPLAR / OAK | 随机 2L NSCLC，atezolizumab vs 多西他赛 | 试验中有 OS + PFS | RNA 不在 GEO | 设计正确；本仓库未开放数据 |
| TCGA-LUAD / LUSC | 多为未治疗的原发切除 | OS / PFI | 两基因都有 | **不是 ICI。** 预后 ≠ 预测 |
| IMvigor210 | 尿路上皮，atezolizumab | OS 可通过 `IMvigor210CoreBiologies` 公开 | 该包中两基因都有 | 病种不对，不能当肺的结论 |

If a later agent finds another open lung ICI series with OS **and** both
genes, add it to this table with a computed *n* / event count. Do not add it
from a paper abstract.

若后续工作找到另一个同时具备 OS **和** 两基因的开放肺 ICI 序列，用计算得到的
*n* / 事件数加入此表。不要从论文摘要抄入。

---

## 12. Recommended reporting order / 建议报告顺序

1. Cohort flow: *n*, events, reverse-KM follow-up, time unit, event coding.
2. Continuous Cox, HR per SD, partial LR *p*, EPV, Schoenfeld, RMST.
3. Prespecified median-split KM as a picture, not as the primary *p*.
4. If an optimal cut is shown: naive *p*, permutation *p*, bootstrap of the
   percentile, frozen-cut external HR.
5. Time-dependent AUC at times that still have support, using the Cox LP.
6. DCB ≥ 6 months as a secondary binary, with the unclassifiable count and
   (if anyone remains) a landmark Cox.
7. Locked external cohort. Same transformation. No re-cut.
8. A one-sentence limit: single-arm, small events, no PD-L1/TMB/purity,
   not OAK/POPLAR.

1. 队列流程：*n*、事件数、反向 KM 随访、时间单位、事件编码。
2. 连续 Cox、每 SD 的 HR、偏 LR *p*、EPV、Schoenfeld、RMST。
3. 预先指定的中位数切分 KM 作为图，而不是作为主要 *p*。
4. 若展示最优切点：朴素 *p*、置换 *p*、分位数自助法、冻结切点的外部 HR。
5. 在仍有支持的时点用 Cox 线性预测变量做时依 AUC。
6. DCB ≥ 6 个月作为次要二分类，附不可分类人数，以及（若还有人）路标 Cox。
7. 锁定外部队列。同一变换。不再切。
8. 一句限制：单臂、事件少、无 PD-L1/TMB/纯度、不是 OAK/POPLAR。

---

## 13. Worked example (computed, not invented) / 实算示例（计算所得，非编造）

Source: `demo/results/demo_results.md` and `demo/results/demo_results.json`,
produced by `demo/02_analysis.py` from the GEO files listed in
`demo/data/provenance.json`. Marker = `log2(TPM+1)`, then z-scored. HR is
per 1 SD. Full printout is in that markdown file; this section keeps only
the numbers that teach a methods point.

来源：`demo/results/demo_results.md` 与 `demo/results/demo_results.json`，由
`demo/02_analysis.py` 根据 `demo/data/provenance.json` 所列 GEO 文件生成。
标志物 = `log2(TPM+1)` 再 z 化。HR 为每 1 SD。完整打印见该 markdown；本节只保留
能说明方法问题的数字。

### 13.1 Cohort clocks / 队列时钟

| | GSE135222 | GSE190265 |
|---|---|---|
| *n* / PFS events / censored | 27 / 21 / 6 | 43 / 35 / 8 |
| median PFS, months (95% CI) | 1.94 (1.22–5.52) | 3.50 (1.90–5.60) |
| reverse-KM follow-up, months | 10.64 | 19.50 |
| PFS at 6 months (n at risk) | 25.9% (7) | 32.6% (14) |
| OS deposited? | no | no |
| DCB ≥ 6 mo / NDB / unclassifiable | 7 / 20 / 0 | 14 / 29 / 0 |

Schoenfeld events needed for 80% power, two-sided α = 0.05, equal split:
**65** (HR = 0.50), **196** (HR = 0.67). Neither cohort reaches the smaller
target.

80% 效能、双侧 α = 0.05、均等分组所需 Schoenfeld 事件数：**65**（HR = 0.50），
**196**（HR = 0.67）。两个队列都达不到较小的那个目标。

### 13.2 Continuous Cox — primary / 连续 Cox——主分析

| Marker | Cohort | HR / SD (95% CI) | Wald *p* | partial LR *p* | Harrell *C* | Uno *C* |
|---|---|---|---|---|---|---|
| TACSTD2 | GSE135222 | 1.07 (0.68–1.69) | 0.78 | 0.77 | 0.562 | 0.567 |
| TACSTD2 | GSE190265 (discovery-scaled) | 1.17 (0.82–1.68) | 0.39 | 0.38 | 0.511 | 0.520 |
| CLDN4 | GSE135222 | 1.12 (0.72–1.75) | 0.61 | 0.60 | 0.550 | 0.552 |
| CLDN4 | GSE190265 (discovery-scaled) | 1.10 (0.79–1.51) | 0.58 | 0.58 | 0.499 | 0.496 |
| CD8A (comparator) | GSE135222 | 0.73 (0.49–1.08) | 0.12 | 0.12 | 0.628 | 0.632 |
| CD8A (comparator) | GSE190265 (discovery-scaled) | 0.75 (0.56–1.02) | 0.066 | 0.068 | 0.574 | 0.567 |

Every TACSTD2 / CLDN4 interval includes 1. That is a **null on a
hopelessly small clock**, not a proof of no effect. Age + sex adjustment on
GSE135222 drops EPV to 7.0; the templates would refuse that model at the
default floor of 10.

所有 TACSTD2 / CLDN4 区间都包含 1。这是**在过小时钟上的阴性**，不是无效应的证明。
在 GSE135222 上校正年龄 + 性别使 EPV 降到 7.0；模板在默认门槛 10 下会拒绝该模型。

### 13.3 Why the “optimal” cut lies / “最优”切点如何说谎

On GSE135222, 13 candidate cuts, minimum 25% per arm:

| Marker | Median-split HR (*p*) | Max-selected HR (naive *p*) | Permutation *p* (1000) | Bootstrap cut percentiles |
|---|---|---|---|---|
| TACSTD2 | 1.83 (0.17) | 2.18 (0.081) | **0.34** | 15–85 |
| CLDN4 | 1.05 (0.91) | 1.48 (0.39) | **0.90** | 4–85 |
| TACSTD2+CLDN4 mean-z (exploratory) | 2.36 (0.048) | 3.23 (0.0066) | **0.068** | 19–78 |
| CD8A (comparator) | 0.34 (0.020) | 0.16 (0.0012) | 0.016 | 7–81 |

The exploratory mean-z median split is *p* = 0.048 and the max-selected
naive *p* is 0.0066. After permuting the score, *p* = 0.068. That is the
entire point of §5. Do not put the 0.0066 in a figure legend.

探索性 mean-z 的中位数切分 *p* = 0.048，最大选择切点的朴素 *p* = 0.0066。
置换分数后 *p* = 0.068。这就是第 5 节的全部要点。不要把 0.0066 写进图注。

### 13.4 Frozen cuts do not travel / 冻结切点走不出去

Discovery-scaled z on GSE190265:

| Marker | Validation mean (SD) | % above discovery “optimal” cut | Frozen-cut HR | Own-median HR |
|---|---|---|---|---|
| TACSTD2 | −0.58 (0.99) | 9% (4 / 39) | 1.78 (0.61–5.18) | 1.08 (0.55–2.09) |
| CLDN4 | −0.86 (0.52) | **0%** (0 / 43) | not estimable | 1.02 (0.52–2.01) |

A cut that puts nobody in the high arm of the next library is not a
biomarker threshold. It is a batch difference.

把下一文库的高组放成空集的切点不是生物标志物阈值，只是批次差异。

### 13.5 Time-dependent AUC and orientation / 时依 AUC 与方向

GSE135222, Cox linear predictor, IPCW. 12-month AUC dropped (last event at
8.4 months).

| Marker | AUC(3 mo) | AUC(6 mo) | If raw marker is passed as risk |
|---|---|---|---|
| TACSTD2 | 0.559 | 0.571 | same (log-HR = +0.066) |
| CLDN4 | 0.612 | 0.443 | same (log-HR = +0.116) |
| CD8A | **0.635** | **0.736** | **0.365 / 0.264** (log-HR = −0.313) |

CD8A is the orientation lesson: a protective marker fed in as a “risk”
score reports `1 − AUC`. TACSTD2 / CLDN4 happen to have a positive (null)
log-HR, so the trap is invisible unless you also run a gene that goes the
other way.

CD8A 是方向课：把保护性标志物当“风险”分数会报告 `1 − AUC`。TACSTD2 / CLDN4
的 log-HR 恰好为正（且接近无效），除非再跑一个方向相反的基因，否则看不出这个陷阱。

### 13.6 DCB ≥ 6 months / 持久获益

GSE135222 binary DCB (7 vs 20, 0 unclassifiable): TACSTD2 apparent AUC =
0.429 (Mann–Whitney *p* = 0.61); CLDN4 AUC = 0.557 (*p* = 0.69). Landmark
Cox at 6 months is not estimable (7 patients, 1 later event). GSE190265
binary DCB (14 vs 29, 0 unclassifiable) is equally null for both genes
(AUC 0.46 and 0.49). These binary AUCs agree, to the reported precision,
with the independent DCB-only analysis of GSE135222 in this repository
(PR 87). They are the same clock, coarsened.

GSE135222 二分类 DCB（7 vs 20，0 例不可分类）：TACSTD2 表观 AUC = 0.429
（Mann–Whitney *p* = 0.61）；CLDN4 AUC = 0.557（*p* = 0.69）。6 个月路标 Cox
不可估计（7 例，1 个后续事件）。GSE190265 二分类 DCB（14 vs 29，0 例不可分类）
对两基因同样为阴性（AUC 0.46 与 0.49）。这些二分类 AUC 在报告精度内与本仓库
对 GSE135222 的独立 DCB 分析（PR 87）一致。它们是同一时钟的粗化。

### 13.7 Competing risks (simulation only) / 竞争风险（仅模拟）

On *n* = 400 simulated patients (185 progressions, 139 competing deaths,
76 censored), 12-month progression risk is 45.3% by Aalen–Johansen vs
60.3% by 1−KM that censors death: +15.0 percentage points. GEO files here
have no cause field.

在 *n* = 400 的模拟患者中（185 进展，139 竞争死亡，76 删失），12 个月进展风险
Aalen–Johansen 为 45.3%，把死亡当删失的 1−KM 为 60.3%：高估 15.0 个百分点。
此处 GEO 文件没有原因字段。

Figures: `demo/results/fig_km_tacstd2_cldn4_median.png`,
`demo/results/fig_cutpoint_profile.png`.

---

## 14. Templates / 模板

| Task | Python (lifelines / sksurv) | R (survival / survminer) |
|---|---|---|
| Continuous Cox + PH + EPV guard | `templates/python/cox_ph.py` | `templates/R/cox_ph.R` |
| Median vs max-selected cut + permutation | `templates/python/cutpoint.py` | `templates/R/cutpoint.R` |
| Time-dependent AUC (Cox LP) | `templates/python/time_dependent_auc.py` | `templates/R/time_dependent_auc.R` |
| Competing risks (refuses if no competing event) | `templates/python/competing_risks.py` | `templates/R/competing_risks.R` |

CSV columns expected by the templates: `time`, `event` (0/1), `marker`.
Competing-risks templates expect `cause` (0 = censored, 1 = interest, 2+ =
competing).

模板期望的 CSV 列：`time`、`event`（0/1）、`marker`。竞争风险模板期望 `cause`
（0 = 删失，1 = 目标事件，2+ = 竞争）。

```bash
# demo (writes the numbers quoted in §13)
cd methods/survival/demo
python3 01_fetch_geo.py
python3 02_analysis.py
```

---

## 15. What this playbook does not claim / 本手册不声称

- That TACSTD2 or CLDN4 is, or is not, an ICI biomarker.
- Any OS result for these genes (no open GEO lung ICI series checked here
  deposits OS **and** a usable expression matrix for both genes).
- Any OAK/POPLAR gene-level HR.
- That a permutation-corrected cut is clinically actionable.
- That DCB and PFS are independent confirmations of each other.

- TACSTD2 或 CLDN4 是或不是 ICI 生物标志物。
- 这些基因的任何 OS 结果（此处核对过的开放 GEO 肺 ICI 序列中，没有同时存有 OS
  **和** 两基因可用表达矩阵的）。
- 任何 OAK/POPLAR 基因水平 HR。
- 经置换校正的切点具有临床可操作性。
- DCB 与 PFS 可以互相独立确证。
