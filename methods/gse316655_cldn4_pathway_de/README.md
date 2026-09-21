# GSE316655 CLDN4-related NHEJ / STING / IFN / APM

Public 10x MTX/TSV only (Cell Ranger `GRCh38_and_mm10-2020-A`). Human CD45+ cells from SK-MEL-5 tumors and blood in NSG-SGM3 mice humanized with cord-blood CD34+ cells. Liu et al., *Sci Immunol* 2026, PMID 41931598.

```bash
python3 download.py --outdir data/raw
python3 analyze.py --raw data/raw --outdir .
```

Raw matrices stay in `data/` and are not committed. Read `FINDING.md` for the counts. Tables and figures are the result.
