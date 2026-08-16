# B5 HNSC: public ICI CLDN4 versus ORR

## Result

**Not estimable from eligible public data.**

The prespecified question was whether pretreatment tumor **CLDN4** expression is
associated with RECIST objective response to immune-checkpoint inhibition in
HNSCC, defining ORR as CR/PR versus SD/PD. No screened public cohort has both:

1. a CLDN4 expression measurement; and
2. a patient-level RECIST response label joinable to that measurement.

Therefore, **zero primary association tests were run**. There is no RECIST
effect estimate, p-value, AUC, optimized cutpoint, or plot. This is a
data-availability result—not evidence for or against CLDN4 as a response
biomarker.

The key near-match remains **NIVACTOR/GSE212549**: its 80-sample Clariom D
matrix measures CLDN4 (official Entrez 1364 mapping: `TC0700007993.hg.1`;
observed normalized range 3.47–7.53). The paper reports 79 evaluable
expression-profiled cases (PR=12, SD=14, PD=53), but neither GEO nor the
article supplement provides the sample-to-response key. Assigning labels from
aggregate counts or figures would be fabrication.

The NIVACTOR testing series **GSE212550** also measures CLDN4 (same probe;
n=20; range 2.82–4.78) but exposes only long-term versus short-term survival
(LTS=12, STS=8), not RECIST.

Conversely, **CLB-IHN/GSE159067** exposes all 102 RECIST labels
(CR=5, PR=6, SD=27, PD=64), but its 2,559-transcript HTG panel does not contain
CLDN4. GSE93157 has the same failure in its five-case HNSCC subset.

See `candidate_audit.tsv` for the seven-cohort screen and
`analysis_status.json` for the machine-readable conclusion.

## Leftover: GSE179730 (not RECIST ORR)

Liu / Knochelmann et al. 2021 ([PMID 34755131](https://pubmed.ncbi.nlm.nih.gov/34755131/),
[PMC8561238](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8561238/),
[GSE179730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179730),
NCT03021993) is a neoadjuvant nivolumab oral-cavity SCC leftover. It was added
because it is the only screened public matrix that simultaneously measures
CLDN4 and has a locked patient-level outcome table.

Locked facts, independently read from supplement Table S2 / S3
(`liu_table_s2.tsv`):

- Trial n=12. Outcome column: Responder=4, Stable=3, Progressor=5.
  Table S1’s header (4/4/4) does not match Table S2; Table S2 is used.
- Pretreatment RNA-seq exists for 11 tumors (`*.Pre` columns). Pt06/HN006 is a
  Responder without pretreatment RNA (Table S3), so leftover n=3 vs 8.
- GEO phenotype is only Pre / Post / Recurrent. Labels are not in GEO.
- Deposited `GSE179730_RNAseq-combinedCPM.txt.gz` is linear CPM (column sums
  ≈1e6) despite GEO saying log2 CPM. Values used: `log2(CPM+1)`.
- CLDN4 is nonzero in 6/11 pretreatment tumors.
- Radiographic scan1-versus-scan2 change never reaches RECIST PR (−30%) in the
  11 RNA cases (0 radiographic CR/PR). Radiographic RECIST ORR is therefore
  not testable. Table S2 “Responder” is a pathologic/hybrid size-change class
  (final pathologic size versus pretreatment imaging), not radiographic RECIST
  ORR. Pt06 is labeled Responder at −29.5% pathologic change.

Leftover CLDN4-only exact permutation Mann–Whitney (no cutpoint, no TACSTD2):

| Contrast | What it is | n | AUC | rank-biserial | exact p |
|---|---|---|---|---|---|
| Responder vs Stable+Progressor | leftover pathologic-ORR analog; **not RECIST** | 3 vs 8 | 0.25 | −0.50 | 0.248 |
| Responder+Stable vs Progressor | clinical benefit / DCR analog; **not ORR** | 6 vs 5 | 0.133 | −0.73 | 0.045 |

The leftover pathologic-ORR analog is **inconclusive**. The DCR contrast is
reported only so a p=0.045 clinical-benefit result cannot be mistaken for a
hidden ORR finding; it is not the B5 primary claim and is not treated as
validation. Expression is sparse and n=11 is exploratory.

## Guardrails

- Primary endpoint frozen as RECIST ORR: CR/PR versus SD/PD.
- CLDN4 frozen as HGNC CLDN4 / Entrez 1364 / ENSG00000189143.
- No substitution of CLDN3, TACSTD2, another epithelial marker, disease
  control, survival, or pathologic response **for the primary claim**.
- No NIVACTOR response labels inferred from plots, aggregate counts,
  signatures, or publication totals.
- No median or outcome-optimized expression split.
- GSE301741 is not treated as ORR evidence: it uses neoadjuvant pathologic
  tumor regression and does not expose response in GEO sample metadata.
- GSE190575 is additionally a pembrolizumab–afatinib combination cohort and
  cannot isolate an ICI-monotherapy association.
- GSE212550 is survival-selected (LTS vs STS) and is not an ORR cohort.
- GSE179730 leftover tests stay labeled leftover and do not change the
  primary “not estimable” RECIST result.

## Reproduce

Only Python 3 standard-library modules, `curl`, `tar`, and the public processed
files are required.

```bash
bash results/w200/B5_HNSC/download.sh
```

`download.sh` retrieves processed GEO matrices (not FASTQ), obtains the official
Bioconductor Clariom D annotation database, and runs `analyze.py`. Downloaded
inputs are intentionally ignored by Git. `checksums.sha256` records the exact
files used for this result.

Generated files:

- `candidate_audit.tsv` — cohort-level eligibility audit
- `liu_table_s2.tsv` — locked Liu Table S2 (not generated; source table)
- `leftover_per_sample.tsv` — GSE179730 pretreatment CLDN4 and locked outcomes
- `leftover_statistics.tsv` — leftover CLDN4 tests, labeled not RECIST
- `analysis_status.json` — final status and directly verified facts
- `checksums.sha256` — SHA-256 input manifest
- `figures/README.md` — why no association figure exists

## Public sources

- [GSE159067](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159067):
  CLB-IHN targeted expression and patient-level RECIST response.
- [GSE93157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE93157):
  NanoString pan-cancer anti-PD-1 cohort with five HNSCC cases.
- [GSE190575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190575):
  ALPHA pembrolizumab–afatinib cohort.
- [GSE212549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212549):
  NIVACTOR whole-transcriptome training cohort.
- [GSE212550](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212550):
  NIVACTOR testing cohort (LTS vs STS only).
- [Serafini et al., 2024](https://doi.org/10.1136/jitc-2023-007823):
  NIVACTOR clinical outcome totals and supplement.
- [GSE179730](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179730):
  neoadjuvant nivolumab oral-cavity RNA-seq leftover.
- [Liu et al., 2021](https://doi.org/10.1016/j.xcrm.2021.100411):
  Table S2/S3 outcomes and RNA inventory.
- [GSE301741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741):
  neoadjuvant pembrolizumab single-cell cohort.
- [Mints et al., 2026](https://doi.org/10.1016/j.xcrm.2026.102715):
  pathologic-response definition and CLDN4 epithelial-state context.
- [Bioconductor Clariom D annotation](https://bioconductor.org/packages/clariomdhumantranscriptcluster.db/):
  official probe-to-Entrez mapping.

Data and source pages were accessed on 2026-08-16.

## 中文摘要

结论：**现有符合条件的公共数据无法估计 CLDN4 与 HNSCC ICI 客观缓解率
（RECIST ORR）的关联。** GSE159067 和 GSE93157 有患者级 RECIST 疗效，但检测
面板不含 CLDN4；GSE212549 检测了 CLDN4，却未公开可与表达矩阵对应的患者级疗效
标签；GSE212550 同样有 CLDN4，但只有生存长短分组。因此未进行主分析统计检验，
也未用替代基因、替代终点、推断标签或结局优化阈值来制造主结果。

GSE179730 是新辅助纳武利尤单抗的遗留队列：有 CLDN4 和 Table S2 患者级结局，
但该结局是病理/混合肿瘤缩小，不是影像 RECIST ORR（RNA 可及病例中影像 PR=0）。
遗留的病理 ORR 类似对比（3 vs 8）无结论（AUC=0.25，精确 p=0.248）。临床获益
对比（p=0.045）是疾病控制率，不是 ORR，不能当作 B5 主结论。该记录只表示公共
RECIST 数据不可识别，并不表示 CLDN4 有效或无效。
