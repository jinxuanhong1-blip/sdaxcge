"""Pre-specified gene lists for the scRNA TLS / B vs TACSTD2 / CLDN4 meta.

TLS12 is the Coppola / 12-chemokine signature (CCL2/3/4/5/8/18/19/21,
CXCL9/10/11/13). Scoring recipe is stated in playbook.md: mean of per-gene
within-cohort z-scores of patient-mean log1p(CP10K). This is the GEP18 analog,
not ssGSEA.
"""

from __future__ import annotations

TARGETS = ["TACSTD2", "CLDN4"]

TLS12 = [
    "CCL2",
    "CCL3",
    "CCL4",
    "CCL5",
    "CCL8",
    "CCL18",
    "CCL19",
    "CCL21",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CXCL13",
]

TLS_STRUCT = ["MS4A1", "CD79A", "CD79B", "CR2", "CXCR5", "SELL", "LAMP3", "CCR7"]

LINEAGE = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCR1"],
    "B": ["CD79A", "CD79B", "MS4A1", "CD19"],
    "plasma": ["JCHAIN", "MZB1", "SDC1"],
    "myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "C1QA"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}

NORMAL_LUNG = [
    "SFTPA1",
    "SFTPA2",
    "SFTPC",
    "SFTPB",
    "AGER",
    "SCGB1A1",
    "SCGB3A1",
    "SCGB3A2",
    "TPPP3",
    "FOXJ1",
    "CAPS",
]

QC = ["PTPRC", "MKI67"]

PANEL = sorted(
    set(TARGETS)
    | set(TLS12)
    | set(TLS_STRUCT)
    | {g for genes in LINEAGE.values() for g in genes}
    | set(NORMAL_LUNG)
    | set(QC)
)
