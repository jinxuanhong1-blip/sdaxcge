# align_cldn4_ifn

Public test of whether CLDN4 (and TACSTD2) loss opens IFN / MHC-I / APM.

**Start here:** [WRITEUP.md](WRITEUP.md)

| path | what |
|---|---|
| `WRITEUP.md` | honest result |
| `tables/verdict.tsv` | one row per contrast |
| `tables/core6_per_gene.tsv` | IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A |
| `tables/geneset_results.tsv` | set-level competitive + permutation |
| `tables/ifn_apm_per_gene.tsv` | full ISG + APM gene matrix |
| `figures/` | QC, CORE6 heatmap, IFNα bars |
| `scripts/` | download + analysis |

Reproduce: `python3 scripts/05_analyze.py && python3 scripts/06_report.py` after `bash scripts/02_download.sh`.
