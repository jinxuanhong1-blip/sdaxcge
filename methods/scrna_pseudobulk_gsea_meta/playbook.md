# Methods: malignant-cell patient pseudobulk GSEA (TACSTD2-high vs low)

中文与英文并列。**Additive only.** User A8 (TCGA-LUAD/LUSC keratin / tight-junction up, Hallmark EMT down) is **taken as given** and is **not** re-run. This extra asks the same prerank question in public **lung tumor scRNA**, after collapsing **malignant cells to one row per patient**.

---

## 范围 / Scope

| In | Out |
|---|---|
| Public GEO lung **tumor** scRNA with a processed whole-transcriptome matrix | TCGA A8 (not re-audited) |
| Patient-level **malignant** (or author epithelial-in-tumor) pseudobulk | Mixed-cell bulk; cell-level GSEA (pseudoreplication) |
| TACSTD2 median-split prerank GSEA for keratin / TJ / EMT / IFN / MHC-I | FASTQ / SRA / dbGaP |
| Meta NES across cohorts that pass the n gate | A pooled “new A8 NES” quoted as if it were TCGA |

Named series (locked before NES): **GSE207422, GSE241934, GSE253013, GSE131907, GSE291670, GSE205335**.

---

## 锁定设计 / Locked design (before NES)

| 项目 Item | 选择 Choice | 诚实限制 Limitation |
|---|---|---|
| Unit | **Patient** (sample if 1:1 with patient) | Multiple tumor sites from one patient are summed, not averaged as independent. |
| Compartment | Author `Malignant` if deposited; else author tumor `Epi`; else marker malignant-like | Marker gates are not CopyKAT. GSE207422 and GSE291670 have **no** public malignant table. |
| Min cells | ≥ **30** malignant/Epi cells per patient | Patients below the gate are dropped, not imputed. |
| Transform | UMI sum → `log2(CPM+1)` | Not log1p CP10k (that is a cell-level scale). |
| Primary split | TACSTD2 **median** high vs low, Welch *t* | Quartile only if **both arms ≥ 6**. |
| Complementary | Spearman ρ of every gene vs continuous TACSTD2 | Uses every gated patient; not a two-group *t*. |
| GSEA | Same prerank engine as A8 (`gsea_core.py`): weighted KS *p*=1, 1000 gene-set permutations, seed=42 | **Gene-set** permutation, not sample permutation. |
| Primary sets | A8 freeze 12 (keratin/TJ/EMT) **plus** Hallmark IFN-α, Hallmark IFN-γ, `CUSTOM_MHC_I_ANTIGEN_PRESENTATION` (15) | BH-FDR **within these 15** per contrast. Hallmark EMT decides the EMT arm. |
| MHC-I set | Locked custom classical MHC-I / processing panel (no MHC-II) | Not an MSigDB term. |
| IFN / MHC-I | Reported; **undirected** | They do **not** enter the A8-style `supportive` verdict. |
| Meta | Among cohorts with **both arms ≥ 5**: median NES, and Stouffer of **NES signs** \(z=\sum\mathrm{sign}(\mathrm{NES}_i)/\sqrt{k}\) | Sign Stouffer ignores NES magnitude. Do not quote a TCGA NES. n is hypothesis-generating. |
| GSE253013 | Inventory only | Public RDS is 9.3 GB; 16 GB RAM cannot load it. Prior public extract: 9 treatment-naïve LUAD, ~6–7 tumors with ≥30 malignant-like cells — **below** the meta n gate even if loaded. |

Positive NES = enriched in TACSTD2-high.

`TACSTD2` is not a member of the 15 primary sets (checked at runtime).

---

## 恶性细胞规则 / Malignant-cell rules

| Accession | Malignant definition | Patient key | Tumor-only rule |
|---|---|---|---|
| GSE207422 | Marker: epithelial lineage max-score **and** normal-lung score ≤ 75th percentile of epithelial cells. Lineage markers as in the public A3 extract. | GEO `Patient` (1 sample / patient) | All deposited scRNA samples |
| GSE241934 | Author `major.cell.type == Epi` on resected tumors (IIT + real-world matrices) | `sampleID` | No normal-lung matrix deposited |
| GSE253013 | Not GSEA'd (see above) | MRC00x | Tumor + adjacent exist; no ICI labels |
| GSE131907 | Author `Cell_subtype == Malignant cells` | `Sample` | Drop `nLung` / `nLN` |
| GSE291670 | Same marker malignant-like rule as GSE207422 | GEO sample (MPR-1… / Non-MPR-1…) | 6 tumors only |
| GSE205335 | Author `lineage.sub == Malignant cells` | SOFT `patient` | Drop Normal Lung / Normal LN / Normal Brain; sum remaining sites |

---

## 功效与措辞 / Power and wording

- GSEA is **run** if both median arms ≥ 3 (so a number exists).
- A cohort is **underpowered** if either arm < 5. It is shown with NES/FDR/n and **excluded from meta**.
- Verdict labels match A8 (Hallmark EMT only): `supportive`, `keratin_TJ_up_Hallmark_EMT_opposite`, `keratin_TJ_up_Hallmark_EMT_null`, `partial`, `mixed`, `null`, `contradicts_TJ_KRT`, `n_too_small`.
- Write the observed **NES / FDR / n_high / n_low**. Do not write “A8 conserved in scRNA” unless a gated cohort is `supportive`.
- Meta median NES is a **descriptive** concordance statistic, not a new TCGA result.

---

## 复现 / Reproduce

```bash
pip install -r requirements.txt
bash scripts/scrna_pseudobulk_gsea_meta/download.sh /tmp/scrna_pb_gsea_data
python3 scripts/scrna_pseudobulk_gsea_meta/pseudobulk.py
python3 scripts/scrna_pseudobulk_gsea_meta/analyze.py
```

Outputs: `results/scrna_pseudobulk_gsea_meta/`.

---

## 中文摘要

本手册只补一层：公开肺癌肿瘤 scRNA，把恶性细胞按患者加总，再按 TACSTD2 中位数做与 A8 相同引擎的 prerank GSEA（角蛋白 / 紧密连接 / EMT / IFN / MHC-I）。**不重跑 TCGA A8。** 主对比是患者级中位数，不是细胞级。Meta 只用两臂都 ≥5 的队列，报告中位 NES 和 NES 符号的 Stouffer，并写明 n。GSE253013 的 9.3 GB RDS 在 16 GB 内存下无法载入，且患者数过关后仍不够 meta，只做清单。
