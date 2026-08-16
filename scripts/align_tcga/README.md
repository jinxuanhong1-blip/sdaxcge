# USER-ALIGN scripts

Reproduce the TACSTD2 vs immune / cytotoxic / exhaustion Spearman and
purity-partial analysis. All inputs are public/open.

```bash
# 1. download matrices (Xena HiSeqV2, Aran CPE, GEO durvalumab FPKMs)
bash scripts/align_tcga/download_data.sh

# 2. optional: confirm the ESTIMATE Python port matches the official sample
python3 scripts/align_tcga/estimate_score.py --selftest
#    (needs results/align_tcga/data/sample_input.txt from the ESTIMATE
#     R-package tarball; not committed)

# 3. TCGA LUAD+LUSC (ESTIMATE / ABSOLUTE / CPE partials)
python3 scripts/align_tcga/run_tcga.py

# 4. OncoSG LUAD 2020 (cBioPortal public API)
python3 scripts/align_tcga/run_oncosg.py

# 5. GSE253564 + GSE248378 durvalumab NSCLC
python3 scripts/align_tcga/run_durva.py
```

Outputs land in `results/align_tcga/`. See `results/align_tcga/WRITEUP.md`.
