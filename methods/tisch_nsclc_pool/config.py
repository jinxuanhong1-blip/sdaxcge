"""TISCH2 NSCLC catalog for the additive TACSTD2/CLDN4 vs T/NK pool.

Public objects only. TISCH2 h5 + CellMetainfo are the only inputs.
Treatment tags are from the TISCH2 gallery (2026-08-16), not invented GEO recodes.
"""

from __future__ import annotations

TISCH_BASE = "https://tisch.compbio.cn/static/data"

# TISCH2 major-lineage vocabulary (strip trailing spaces at load time)
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

GENES = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC"]
TARGET_GENES = ["TACSTD2", "CLDN4"]

MIN_EPI_CELLS = 20
MIN_TNK_CELLS = 20
MIN_N_FOR_SPEARMAN = 5
MIN_N_FOR_FISHER_Z = 6

# Tumor-like tissue tokens (lowercase substring match). NAT/normal/PBMC/nLN excluded.
TUMOR_TOKENS = (
    "tumor",
    "tlung",
    "metast",
    "mbrain",
    "mln",
    "pleural",
)
EXCLUDE_TISSUE_TOKENS = (
    "normal",
    "nat",
    "nln",
    "pbmc",
    "blood",
    "broncho",
)

# Gallery rows from https://tisch.compbio.cn/gallery/?cancer=NSCLC (17 datasets).
# ici_gallery: TISCH2 "Treatment" column. CellMetainfo ICI columns are audited at runtime.
DATASETS = {
    "NSCLC_GSE151537": {
        "pmid": "33691136",
        "paper": "Chiou et al. 2021, Nat Commun (TCR / sorted T cells)",
        "ici_gallery": "Immunotherapy",
        "priority": "user",
        "download_h5": True,
        "skip_reason_expected": "T-sorted only (no epithelial/malignant cells)",
    },
    "NSCLC_GSE146100": {
        "pmid": "33820821",
        "paper": "Yang et al. 2021, multiple primary LUAD after pembrolizumab (1 patient, 3 nodules)",
        "ici_gallery": "Immunotherapy",
        "priority": "user",
        "download_h5": True,
    },
    "NSCLC_GSE117570": {
        "pmid": "31033233",
        "paper": "Song et al. 2019, Nat Commun (Lambrechts-adjacent 10x NSCLC)",
        "ici_gallery": "None",
        "priority": "user",
        "download_h5": True,
    },
    "NSCLC_GSE127465": {
        "pmid": "30979687",
        "paper": "Zilionis et al. 2019, Immunity",
        "ici_gallery": "None",
        "priority": "user",
        "download_h5": True,
    },
    "NSCLC_GSE131907": {
        "pmid": "32385277",
        "paper": "Kim et al. 2020, Nat Commun (LUAD atlas)",
        "ici_gallery": "None",
        "priority": "user",
        "download_h5": True,
    },
    "NSCLC_GSE148071": {
        "pmid": "33953163",
        "paper": "Wu et al. 2021, Nat Commun (advanced NSCLC)",
        "ici_gallery": "None",
        "priority": "user",
        "download_h5": True,
    },
    "NSCLC_EMTAB6149": {
        "pmid": "29988129",
        "paper": "Lambrechts et al. 2018, Nat Med",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": True,
        "skip_reason_expected": "CellMetainfo has no Patient/Sample column",
    },
    "NSCLC_GSE127471": {
        "pmid": "31061481",
        "paper": "PBMC-only 10x (1 patient)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": False,
        "skip_reason_expected": "PBMC only; 0 epithelial/malignant cells",
    },
    "NSCLC_GSE139555": {
        "pmid": "32103181",
        "paper": "Wu/Zhang T-cell atlas (immune-sorted)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": False,
        "skip_reason_expected": "T/immune-sorted; 0 epithelial/malignant cells",
    },
    "NSCLC_GSE143423": {
        "pmid": None,
        "paper": "Brain-metastasis 10x (3 patients)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": True,
    },
    "NSCLC_GSE99254": {
        "pmid": "29942094",
        "paper": "Guo et al. 2018, Nat Med (T-sorted Smart-seq2)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": False,
        "skip_reason_expected": "T-sorted; 0 epithelial/malignant cells",
    },
    "NSCLC_GSE149655": {
        "pmid": "32891189",
        "paper": "2-patient 10x (alveolar-rich; no Malignant call)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": True,
    },
    "NSCLC_GSE150660": {
        "pmid": "32675368",
        "paper": "Metastasis 10x (2 patients)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": True,
    },
    "NSCLC_GSE153935": {
        "pmid": None,
        "paper": "12-patient 10x (epithelial + alveolar present)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": True,
    },
    "NSCLC_GSE162498": {
        "pmid": "33514641",
        "paper": "11-patient 10x (epithelial present; T-rich)",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": True,
    },
    "NSCLC_GSE176021_aPD1": {
        "pmid": "34290408",
        "paper": "Caushi et al. 2021, Nature (MANA T cells, anti-PD-1)",
        "ici_gallery": "Immunotherapy",
        "priority": "other_nsclc",
        "download_h5": False,
        "skip_reason_expected": "T-sorted (~817k cells); 0 epithelial/malignant; CellMetainfo has no Response column",
    },
    "NSCLC_GSE179373": {
        "pmid": "35584630",
        "paper": "T-only; Treatment=Systemic therapy in CellMetainfo",
        "ici_gallery": "None",
        "priority": "other_nsclc",
        "download_h5": False,
        "skip_reason_expected": "T-sorted; 0 epithelial/malignant cells",
    },
}
