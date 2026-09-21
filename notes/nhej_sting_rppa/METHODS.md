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
- Deepen script: `scripts/nhej_sting_rppa/deepen_table_s1.py`
- Gene sets, MSigDB symbols vendored under `scripts/nhej_sting_rppa/genesets/`: Reactome NHEJ (NHEJ), Reactome STING-mediated induction of host immune responses, Hallmark interferon alpha, Hallmark interferon gamma
- NHEJ core drops symbols starting with H2, H3, or H4. Those histones are tested on a separate row
- Predicted NHEJ-down under CLDN4-low: q<0.05 and log2(CLDN4-high/CLDN4-low)>0. Predicted STING/IFN-up: q<0.05 and log2<0
- Fisher exact is two-sided on the q<0.05 slice: set-in-predicted-direction versus set-in-other-direction, against the rest of the 1,582 q<0.05 genes
- Bitler supplement scrape, 2026-09-21: Europe PMC XML for 15 Bitler CLDN4/RPPA papers; open-access supplementary zips text-searched for RPPA, PRKDC, 53BP1, XRCC1. Closed NIHMS packages (DUSP PMC9357222, VDX-111 PMC12481673, Jordan CCR PMC7923250) returned HTML/403 and were not parsed as tables. No fold from another experiment was written into the CLDN4 knockdown score
