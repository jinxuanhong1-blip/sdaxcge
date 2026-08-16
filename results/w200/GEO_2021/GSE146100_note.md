# GSE146100 single-cell side audit

[GSE146100](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146100) contains three resected lung-adenocarcinoma nodules from one woman after three pembrolizumab cycles. W2 responded; W1 and W3 did not.

In its normalized matrix (11,612 cells), detectable TACSTD2/CLDN4 co-expression was 2.17% in W1 (88/4,062), 4.68% in W2 (174/3,716), and 14.81% in W3 (568/3,834). Requiring at least two additional epithelial markers (EPCAM, KRT7/8/18/19, MUC1) changed these to 2.14%, 4.60%, and 14.08%. W3 therefore contains a clear epithelial-like TACSTD2+/CLDN4+ population. This is not a consistent response association: nonresponding W1 has the lowest frequency, while responding W2 is intermediate.

These cells cannot be called malignant from this matrix alone. TACSTD2 and CLDN4 are epithelial, not tumor-specific; normal epithelium, ambient RNA, and doublets remain alternatives. Indeed, 184 co-positive cells also detect PTPRC. GEO supplies no cell annotations, raw counts, or pretreatment sample here. The lesions are not independent patients, so cell-level significance tests would be pseudoreplication. This is a descriptive side result, not patient-level biomarker evidence.

Reproduce with `python3 analyze.py GSE146100_NormData.txt.gz`; counts are in `marker_summary.csv`. Study: [Zhang et al., 2021](https://doi.org/10.1136/jitc-2020-002312).
