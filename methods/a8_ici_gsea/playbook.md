# Methods: keratin / tight-junction / EMT GSEA in public ICI lung bulk

中文与英文并列。This playbook is **additive** to the TCGA-only A8 GSEA slide. It does not re-run TCGA-LUAD/LUSC.

---

## 范围 / Scope

Ask whether **pre-treatment TACSTD2-high** ICI-treated lung tumors enrich the same locked keratin / tight-junction programs and deplete Hallmark EMT that user A8 reported in TCGA.

- **In scope:** open GEO whole-transcriptome **tumor bulk** from human lung cancers treated with PD-1 / PD-L1 (± chemo or SBRT), with `TACSTD2` on the deposited matrix.
- **Out of scope:** TCGA (already A8); FASTQ/SRA; dbGaP / OAK / POPLAR; immune-only panels that lack `TACSTD2`; post-treatment-only matrices used as a primary split; GSE207422 (already GSEA'd in the `hunt_tj_gsea` slice).
- **Not a predictive claim.** A program that tracks TACSTD2 inside a treated cohort is not a TROP2 × ICI interaction.

---

## 锁定设计 / Locked design (before NES)

| 项目 Item | 选择 Choice | 诚实限制 Limitation |
|---|---|---|
| 主对比 Primary | TACSTD2 **median** high vs low, Welch *t* | Quartiles on n≈16 are 4 vs 4. Median is the small-n lock. |
| 敏感性 Sensitivity | Quartile only if **both arms ≥ 6** | Still discards the middle half. |
| 连续互补 Continuous | Spearman ρ of every gene vs TACSTD2, then prerank | Uses every sample. |
| 疗效对比 Response | Responder vs NR Welch *t* only if **both arms ≥ 5** | Deposited labels only. No invented RECIST. |
| GSEA | Preranked weighted KS, *p*=1, 1000 gene-set permutations, seed=42 | Same engine as A8. **Gene-set** permutation, not sample permutation. |
| FDR | BH within the **12 primary sets** per contrast | Not nested Broad FDR. |
| EMT 判定 | **Hallmark EMT** only | Do not quote GOBP EMT as Hallmark. |
| 基因集 Sets | `data/genesets/a8_sets.json` primary 12 | Identical to the TCGA A8 freeze. `TACSTD2` is not in those sets. |

Positive NES = enriched in TACSTD2-high (or in responders).

Do **not** pool NES across cohorts. Do **not** quote a TCGA NES from this run.

---

## 队列与变换 / Cohorts and transforms

Processed public matrices only. Files >2 GB and FASTQ are refused.

| Accession | Transform | Response label used here |
|---|---|---|
| GSE126044 | counts → log2(CPM+1) | GEO responder / non-responder |
| GSE135222 | TPM → log2(TPM+1) if not already logged | GEO-derived DCB = PFS ≥ 180 days |
| GSE166449 | TPM → log2(TPM+1) if not already logged | GEO responder / non-responder |
| GSE253564 | pre-treatment FPKM → log2(FPKM+1) | Only if MPR is **deposited** on GEO; otherwise TACSTD2 split only |
| GSE190265 France3 | TPM → log2(TPM+1) if not already logged | Deposited DCB / PFS in the France3 sample-info table |
| GSE283829 | supplementary raw counts → log2(CPM+1) | GEO `disease stage` as CR vs PD when that field is RECIST-like |

GSE190266 (France4) is listed and skipped if `TACSTD2` is absent from the deposited TPM file.

If a matrix is already on a log2 scale (max < 25 and median < 12), do not log again.

---

## 功效与措辞 / Power and wording

n = 16–32 is **hypothesis-generating**. A modest true effect is not ruled out by FDR>0.05. Write the observed **NES / FDR / n_high / n_low**.

Verdict labels (same as A8):

- `supportive`: Hallmark EMT down **and** a primary TJ set up **and** a primary keratin/barrier set up (FDR<0.05).
- `keratin_TJ_up_Hallmark_EMT_opposite`: TJ+KRT up, Hallmark EMT significantly up.
- `keratin_TJ_up_Hallmark_EMT_null`: TJ+KRT up, Hallmark EMT not significant.
- `partial` / `mixed` / `null` / `contradicts_TJ_KRT`: as in the A8 report.
- `n_too_small` / `no_deposited_response_label`: contrast not run.

Do not write “conserved in ICI tumors” unless the median-split Hallmark EMT, a TJ set, and a keratin set all pass in the **same** leftover cohort.

---

## 复现 / Reproduce

```bash
pip install -r requirements.txt
bash scripts/a8_ici_gsea/download.sh /tmp/a8_ici_gsea_data
python3 scripts/a8_ici_gsea/analyze.py
```

Outputs: `results/a8_ici_gsea/REPORT.md`, `tables/`, `figures/`.

---

## 中文摘要

本手册只补一层：**公开、治疗前、ICI 处理过的肺癌 bulk**，用与 A8 相同的 12 个基因集做 prerank GSEA。不重跑 TCGA。主对比是 TACSTD2 中位数高低，不是四分位（除非两臂都 ≥6）。疗效对比仅在两臂都 ≥5 且标签来自 GEO 时才跑。报告真实的 NES、BH-FDR 和 n。Hallmark EMT 说了算，GOBP EMT 不算 Hallmark。禁止把各队列 NES 加权成一个“ICI 总 NES”。
