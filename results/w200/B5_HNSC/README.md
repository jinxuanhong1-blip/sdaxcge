# B5 HNSC: public ICI CLDN4 versus ORR

## Result

**Not estimable from eligible public data.**

The prespecified question was whether pretreatment tumor **CLDN4** expression is
associated with RECIST objective response to immune-checkpoint inhibition in
HNSCC, defining ORR as CR/PR versus SD/PD. No screened public cohort has both:

1. a CLDN4 expression measurement; and
2. a patient-level RECIST response label joinable to that measurement.

Therefore, zero association tests were run. There is no effect estimate, p-value,
AUC, optimized cutpoint, or plot. This is a data-availability result—not evidence
for or against CLDN4 as a response biomarker.

The key near-match is **NIVACTOR/GSE212549**: its 80-sample Clariom D matrix
measures CLDN4 (official Entrez 1364 mapping:
`TC0700007993.hg.1`; observed normalized range 3.47–7.53). The paper reports 79
evaluable expression-profiled cases (PR=12, SD=14, PD=53), but neither GEO nor
the article supplement provides the sample-to-response key. Assigning labels
from aggregate counts or figures would be fabrication.

Conversely, **CLB-IHN/GSE159067** exposes all 102 RECIST labels
(CR=5, PR=6, SD=27, PD=64), but its 2,559-transcript HTG panel does not contain
CLDN4. GSE93157 has the same failure in its five-case HNSCC subset.

See `candidate_audit.tsv` for the complete five-cohort screen and
`analysis_status.json` for the machine-readable conclusion.

## Guardrails

- Primary endpoint frozen as RECIST ORR: CR/PR versus SD/PD.
- CLDN4 frozen as HGNC CLDN4 / Entrez 1364 / ENSG00000189143.
- No substitution of CLDN3, another epithelial marker, disease control,
  survival, or pathologic response.
- No response labels inferred from plots, aggregate counts, signatures, or
  publication claims.
- No median or outcome-optimized expression split.
- GSE301741 is not treated as ORR evidence: it uses neoadjuvant pathologic tumor
  regression and does not expose response in GEO sample metadata.
- GSE190575 is additionally a pembrolizumab–afatinib combination cohort and
  cannot isolate an ICI-monotherapy association.

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
- [Serafini et al., 2024](https://doi.org/10.1136/jitc-2023-007823):
  NIVACTOR clinical outcome totals and supplement.
- [GSE301741](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE301741):
  neoadjuvant pembrolizumab single-cell cohort.
- [Mints et al., 2026](https://doi.org/10.1016/j.xcrm.2026.102715):
  pathologic-response definition and CLDN4 epithelial-state context.
- [Bioconductor Clariom D annotation](https://bioconductor.org/packages/clariomdhumantranscriptcluster.db/):
  official probe-to-Entrez mapping.

Data and source pages were accessed on 2026-08-16.

## 中文摘要

结论：**现有符合条件的公共数据无法估计 CLDN4 与 HNSCC ICI 客观缓解率
（ORR）的关联。** GSE159067 和 GSE93157 有患者级 RECIST 疗效，但检测面板
不含 CLDN4；GSE212549 检测了 CLDN4，却未公开可与表达矩阵对应的患者级疗效
标签。因此未进行统计检验，也未用替代基因、替代终点、推断标签或结局优化阈值
来制造结果。该结论只表示公共数据不可识别，并不表示 CLDN4 有效或无效。
