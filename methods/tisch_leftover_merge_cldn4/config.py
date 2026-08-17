"""Leftover TISCH2 NSCLC objects (not GSE131907 / 148071 / 205335 / 207422).

CLDN4-only. Public TISCH2 h5 + CellMetainfo only. Skip files >2 GB.
Skip a series when both CLDN4 (epithelial/malignant) and T/NK are absent.
"""

from __future__ import annotations

TISCH_BASE = "https://tisch.compbio.cn/static/data"
MAX_H5_BYTES = 2 * 1024 * 1024 * 1024  # 2 GB

# Already used in the 131907/148071/205335/207422 set — not leftover.
EXCLUDED_SET = {
    "GSE131907",
    "GSE148071",
    "GSE205335",
    "GSE207422",
}

MALIGNANT = {"Malignant"}
EPITHELIAL_NONMALIGNANT = {
    "Epithelial",
    "Alveolar",
    "AT1",
    "AT2",
    "Basal",
    "Ciliated",
    "Club",
}
TNK_LINEAGES = {
    "CD8T",
    "CD8Tex",
    "CD4Tconv",
    "Treg",
    "Tprolif",
    "TMKI67",
    "NK",
    "Tcell",
    "NKT",
    "ILC",
}

GENES = ["CLDN4", "EPCAM", "PTPRC", "CD3D", "CD8A", "NKG7"]
MIN_EPI_CELLS = 20
MIN_TNK_CELLS = 20
MIN_N_SPEARMAN = 5
MIN_N_MERGE = 8
MIN_N_Q4 = 6

TUMOR_TOKENS = ("tumor", "tlung", "metast", "mbrain", "mln", "pleural")
EXCLUDE_TISSUE_TOKENS = ("normal", "nat", "nln", "pbmc", "blood", "broncho")

# User-requested leftover hunt. download_h5 is decided after metainfo + size.
DATASETS = {
    "NSCLC_GSE117570": {
        "gse": "GSE117570",
        "pmid": "31033233",
        "paper": "Song et al. 2019, Nat Commun (10x NSCLC; tumor + NAT)",
        "ici_gallery": "None",
        "expected": "malignant + T/NK; few eligible patients",
    },
    "NSCLC_GSE139555": {
        "gse": "GSE139555",
        "pmid": "32103181",
        "paper": "Wu/Zhang T-cell atlas (immune-sorted)",
        "ici_gallery": "None",
        "expected": "T/immune-sorted; 0 epithelial/malignant",
    },
    "NSCLC_GSE146100": {
        "gse": "GSE146100",
        "pmid": "33820821",
        "paper": "Yang et al. 2021, 1 patient / 3 LUAD nodules after pembrolizumab",
        "ici_gallery": "Immunotherapy",
        "expected": "epithelial (no Malignant call); 1 patient",
    },
    "NSCLC_GSE149655": {
        "gse": "GSE149655",
        "pmid": "32891189",
        "paper": "2-patient 10x (alveolar/club; no Malignant call)",
        "ici_gallery": "None",
        "expected": "epithelial_like; n_patients=2",
    },
    "NSCLC_GSE143423": {
        "gse": "GSE143423",
        "pmid": None,
        "paper": "Brain-metastasis 10x (3 patients)",
        "ici_gallery": "None",
        "expected": "malignant + few T/NK",
    },
    "NSCLC_GSE127471": {
        "gse": "GSE127471",
        "pmid": "31061481",
        "paper": "PBMC-only 10x (1 patient)",
        "ici_gallery": "None",
        "expected": "PBMC only; 0 epithelial/malignant",
    },
}
