# Methods

- Date: 2026-09-21
- Paper: Yamamoto / Webb / Bitler, Mol Cancer Ther 2022, DOI 10.1158/1535-7163.MCT-21-0827, PMC8988515
- PRIDE Archive v2 `search/projects?keyword=` for CLDN4, claudin-4, claudin 4, shCLDN4, Bitler, claudin, OVCAR3
- EBI Search domain `pride`, phrase queries `"CLDN4"`, `"claudin-4"`, `"claudin 4"`, `shCLDN4`, `Bitler`
- OmicsDI `ws/dataset/search` restricted to `source:pride` and `source:massive`
- Supplements: Figshare collection 10.1158/1535-7163.c.6543310, files 39984958 (S1), 39984961 (S2), 39984964 (S3), 39984967 (S4)
- Table S1 header row is Excel row 2. q<0.05 count was checked against the paper’s 1,582
- `kd_logic_score`: −1 manuscript reports a significant protein drop after shCLDN4; 0 not named; no +1 in this set. Numeric RPPA folds are null
- `tcga_sign`: sign of Table S1 log2(CLDN4-high/CLDN4-low) only when q<0.05, else 0. This column is RNA, not RPPA
- Full xlsx files are downloaded under `results/nhej_sting_rppa/cache/` at runtime and are not committed
- Script: `scripts/nhej_sting_rppa/score_cldn4_kd_rppa.py`
