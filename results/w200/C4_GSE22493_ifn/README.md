# C4 / GSE22493 — IFN / MHC-I / APM

Honest verdict: **does not support** CLDN4 siRNA opening IFN / MHC-I / APM.
Full bilingual writeup: `notes/w200/C4_GSE22493_ifn/WRITEUP.md`.

| File | Contents |
|---|---|
| `priority_genes.tsv` | IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A (deposited log2 KD/control) |
| `apm_genes.tsv` | 16-gene MHC-I / APM panel |
| `ifn_genes.tsv` | 73-gene IFN_IMMUNE list |
| `geneset_stats.tsv` | Mann-Whitney vs background |
| `cldn4_diagnostic.tsv` | CLDN4 probe 17169: deposited vs ScanArray |
| `sensitivity_deposited_vs_scanarray.tsv` | Priority genes, two quantitations |
| `probe_level.tsv` | Every priority/APM/CLDN4 probe × array |
| `sample_table.tsv` | Channel assignment |
| `key_stats.json` | Machine-readable summary |
| `fig1_priority_apm_log2.png` | Per-array points + gene medians |

```bash
python3 scripts/w200/C4_GSE22493_ifn/download_data.py
python3 scripts/w200/C4_GSE22493_ifn/run_analysis.py
```
