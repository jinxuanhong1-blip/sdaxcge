# FINDING — CLDN4-only IFN/MHC/TJ regulons (GSE148071 AUCell proxy)

**Method actually run: AUCell + public TF–target priors. Not full pySCENIC.** pySCENIC is not installed here. cisTarget motif rankings were not downloaded. Regulons are the union of TRRUST v2 + DoRothEA + CollecTRI targets, scored with an Aibar-style recovery curve on `log1p(CP10k)`. Companion gene-set AUCell (ISG / MHC-I APM / TJ) is pathway activity, not a TF regulon.

**Verdict (honest n=38 patients, not 42):** in putative-malignant epithelium, CLDN4-high cells (within-patient Q4) have **higher TJ regulon AUCell** than CLDN4-low (Q1) — **37/38**, p=1.5×10⁻¹¹ — and **modestly higher IFN regulon AUCell** — **30/38**, p=2.1×10⁻⁶. **MHC master-TF regulons do not move** (CIITA/RFX union **17/38**, p=0.62). NLRC5 has **zero** public prior edges. This is **not** a TJ-up / IFN-down tradeoff (ΔTJ vs ΔIFN ρ=+0.13, p=0.43).

CLDN4 is the only splitter. TACSTD2 is a companion (37/38 track the split) and is **never a gate**. GSE207422 and GSE131907 were not used.

---

## English

### Honest n

| Item | n | Note |
| --- | ---: | --- |
| Patients deposited | **42** | Wu et al., *Nat Commun* 2021, PMID 33953163; stage III/IV NSCLC biopsies |
| GEO processed matrices | 42 | `GSE148071_RAW.tar` / `GSM*_P*_exp.txt.gz` |
| Cells after UMI≥200 | 89,887 | all lineages |
| Epithelial (putative malignant) | 51,215 | marker-argmax; **not** CopyKAT |
| Mal-like | 42,519 | epithelial AND normal-lung score ≤ sample Q75 |
| Occupancy-eligible | **41** | ≥20 mal-like and ≥8 cells per CLDN4 tail |
| Dropped for occupancy | **1** | **P37** (12 mal-like; CLDN4 all zero) |
| Dropped for no CLDN4 split | **3** | **P4, P11, P42** (Q1=Q2=Q3=0) |
| **Primary n** | **38** | occupancy **and** Q3>Q1; **unit of inference** |
| Primary mal-like cells | 42,168 | P3+P17+P41+P1+P25 = 53.7% of cells; test is still patient-level |

Do **not** cite n=42. Cell-level Spearman (n=42,519) is exploratory (pseudoreplication).

### What was scored

| Program | TFs / set | Prior targets used | Held out |
| --- | --- | ---: | --- |
| IFN union | STAT1, STAT2, IRF1, IRF3, IRF7, IRF9 | 996 | CLDN4 |
| MHC union | CIITA, RFX5, RFXANK, RFXAP | 48 | CLDN4 |
| MHC+IRF1 | MHC union + IRF1 | 192 | CLDN4 |
| TJ union | GRHL2, ELF3, KLF4, OVOL1 (GRHL1/OVOL2 <5 targets) | 167 | CLDN4 |
| NLRC5 | MHC-I master TF | **0** (no TRRUST/DoRothEA/CollecTRI edges) | — |
| Control TF | HIF1A | 505 | CLDN4 |
| Gene-set TJ | 12 TJ/polarity genes | — | **CLDN4** |
| Gene-set IFN | 40-gene type-I ISG core | — | — |
| Gene-set MHC | 21-gene MHC-I APM | — | — |
| Control set | 10-gene OXPHOS | — | — |

Public priors list CLDN4 as a GRHL2 and STAT2 target; those edges are dropped before AUCell.

### Primary table — same cells, paired patients (n=38)

Mal-like epithelium. High = CLDN4 Q4; low = CLDN4 Q1. Mean is the mean of patient-level means. Wilcoxon is two-sided exact. Full file: `results/primary_table.tsv`.

| Feature | n | high vs low | Δ median | high>low | p | Layer |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| **prior_TJ_union** | **38** | 0.079 vs 0.068 | +0.0087 | **37/38** | **1.5e-11** | AUCell TF-target |
| prior_ELF3 | 38 | 0.066 vs 0.047 | +0.014 | 35/38 | 5.5e-09 | AUCell TF-target |
| prior_KLF4 | 38 | 0.083 vs 0.076 | +0.0054 | 36/38 | 1.0e-10 | AUCell TF-target |
| prior_OVOL1 | 38 | 0.080 vs 0.055 | +0.017 | 32/38 | 3.1e-08 | AUCell TF-target |
| prior_GRHL2 | 37 | 0.033 vs 0.023 | +0.0076 | 27/37 | 3.4e-05 | AUCell TF-target |
| gs_tj_no_cldn4 | 38 | 0.042 vs 0.024 | +0.014 | 33/38 | 9.2e-09 | AUCell gene set |
| mean_tj_no_cldn4 | 38 | 0.288 vs 0.193 | +0.069 | **38/38** | 7.3e-12 | mean log1p(CP10k) |
| **prior_IFN_union** | **38** | 0.032 vs 0.030 | +0.0013 | **30/38** | **2.1e-06** | AUCell TF-target |
| prior_STAT1 | 38 | 0.034 vs 0.032 | +0.0018 | 35/38 | 1.1e-07 | AUCell TF-target |
| prior_IRF1 | 38 | 0.043 vs 0.039 | +0.0037 | 33/38 | 5.0e-06 | AUCell TF-target |
| prior_IRF7 | 38 | 0.058 vs 0.048 | +0.0082 | 31/38 | 2.1e-06 | AUCell TF-target |
| gs_ifn_isg | 38 | 0.062 vs 0.058 | +0.0021 | 25/38 | 0.028 | AUCell gene set |
| **prior_MHC_union** | **38** | 0.072 vs 0.073 | −0.0008 | **17/38** | **0.62** | AUCell TF-target |
| prior_CIITA | 38 | 0.097 vs 0.098 | −0.0001 | 18/38 | 0.80 | AUCell TF-target |
| prior_RFX5 | 38 | 0.050 vs 0.056 | −0.0014 | 16/38 | 0.34 | AUCell TF-target |
| prior_MHC_IRF1_union | 38 | 0.047 vs 0.044 | +0.0022 | 28/38 | 0.012 | AUCell (IFN-adjacent) |
| gs_mhc1_apm | 38 | 0.176 vs 0.173 | +0.0007 | 20/38 | 0.70 | AUCell gene set |
| mean_mhc1_apm | 38 | 0.773 vs 0.728 | +0.036 | 27/38 | 0.010 | mean log1p(CP10k) |
| gs_oxphos (control) | 38 | 0.165 vs 0.158 | +0.0008 | 21/38 | **0.19** | AUCell gene set |
| prior_HIF1A (control TF) | 38 | 0.053 vs 0.049 | +0.0031 | 36/38 | 8.3e-08 | AUCell TF-target |
| companion TACSTD2 | 38 | 1.64 vs 0.86 | +0.71 | 37/38 | 1.5e-11 | RNA, not a gate |
| library size (UMI) | 38 | 6586 vs 4997 | +1408 | 35/38 | 2.3e-08 | depth caveat |

### Sample-level Spearman (still n=38)

Patient-mean CLDN4 vs patient-mean AUCell (not the paired tails):

| Contrast | ρ | p |
| --- | ---: | ---: |
| CLDN4 vs prior_TJ_union | +0.50 | 0.0014 |
| CLDN4 vs prior_IFN_union | +0.45 | 0.0052 |
| CLDN4 vs prior_MHC_union | +0.18 | 0.29 |
| CLDN4 vs gs_tj_no_cldn4 | +0.57 | 2.2e-04 |
| CLDN4 vs gs_ifn_isg | +0.11 | 0.51 |
| CLDN4 vs gs_mhc1_apm | +0.02 | 0.89 |
| CLDN4 vs gs_oxphos | −0.38 | 0.021 |
| Δ prior_TJ vs Δ prior_IFN | +0.13 | 0.43 |

Paired tails detect a within-tumor IFN lift that patient-mean ISG AUCell does not (ρ=+0.11, p=0.51). MHC stays null either way.

### What is not supported

- **n=42** as the claim n. Four biopsies cannot be split on CLDN4.
- **Full pySCENIC / motif-supported regulons.** This is an AUCell / TF-target **proxy**.
- **MHC-I master-TF activity higher in CLDN4-high.** CIITA/RFX AUCell and MHC-I APM AUCell are null. NLRC5 was not scorable. The MHC+IRF1 union is an IFN-adjacent set.
- **IFN down in CLDN4-high.** Direction is up, not down. Deltas of TJ and IFN are not a tradeoff.
- **Author malignant IDs.** Epithelium is marker-argmax; mal-like drops the top normal-lung quartile only.
- **Dual-high (TACSTD2 AND CLDN4).** Not run. TACSTD2 tracks the CLDN4 split and is reported only as a companion.
- **ICI response / histology.** GEO has age and sex only.
- **Depth-free specificity.** CLDN4-high cells have more UMI (35/38). Rank-based OXPHOS AUCell is NS (p=0.19); mean OXPHOS and HIF1A AUCell are not. Do not over-read small IFN AUCell deltas as a large program shift.

### Files

- `results/primary_table.tsv` — compact primary table
- `results/n_table.tsv` / `eligibility.tsv` / `per_patient.tsv` / `paired_high_low.tsv` / `tests.tsv`
- `results/regulon_membership.tsv` / `resources/tf_targets_public.tsv`
- `results/figures/fig_paired_regulons.png` / `fig_paired_genesets.png` / `fig_honest_n.png` / `fig_delta_scatter.png`
- `results/summary.json`

### Reproduce

```bash
python3 methods/gse148071_scenic_cldn4/scripts/download.py
python3 methods/gse148071_scenic_cldn4/scripts/download_priors.py
python3 methods/gse148071_scenic_cldn4/scripts/analyze.py
```

---

## 中文

**方法：AUCell + 公开 TF–靶基因先验，不是完整 pySCENIC。** 本环境未安装 pySCENIC，也未下载 cisTarget。Regulon 是 TRRUST v2 / DoRothEA / CollecTRI 的并集，用 Aibar 回收曲线打分。

**结论（诚实 n=38，不是 42）：** GSE148071 假定恶性上皮里，同一患者内 CLDN4 高（Q4）比低（Q1）的 **TJ regulon AUCell 更高（37/38，p=1.5×10⁻¹¹）**，**IFN regulon 也略高（30/38，p=2.1×10⁻⁶）**。**MHC 主控 TF（CIITA/RFX）不升高（17/38，p=0.62）**。NLRC5 在公开先验里 **0 条边**。这不是 TJ 升 / IFN 降的此消彼长（Δ ρ=+0.13，p=0.43）。

只按 CLDN4 分组，不用 TACSTD2 做门控。不用 GSE207422 / GSE131907。存档 42 例里，P37 细胞不够，P4/P11/P42 的 CLDN4 四分位全是 0，不能当高低对比。上皮是 marker-argmax，不是 CopyKAT。CLDN4 高细胞 UMI 更深（35/38）；OXPHOS 的 AUCell 对照不显著（p=0.19），不要把 IFN 的小 Δ 写成大程序改变。
