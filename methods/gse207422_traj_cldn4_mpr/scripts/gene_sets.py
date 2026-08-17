"""Locked gene sets for the GSE207422 CLDN4-only PAGA / DPT slice.

CLDN4 is the readout. It is NOT in the barrier/keratin score (no circularity).
TACSTD2 is a companion only — never a dual-high gate.
A3-malignant-like uses the same normal-lung zero-UMI panel as the given
A3 / dual-high slices. Author CopyKAT barcodes are not public.
"""

from __future__ import annotations

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T": ["CD3D", "CD3E", "CD2"],
    "NK": ["NKG7", "GNLY", "FGFBP2"],
    "B": ["CD79A", "MS4A1"],
    "plasma": ["IGHG1", "MZB1"],
    "myeloid": ["LYZ", "CD68", "CD14"],
    "neutrophil": ["CSF3R"],
    "fibroblast": ["COL1A1", "DCN"],
    "endothelial": ["VWF", "PECAM1"],
    "mast": ["KIT"],
}

# A3 given malignant-like: epithelial AND zero UMI on this panel.
A3_NORMAL_LUNG = ("SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3")

STATES = {
    "AT2": ("SFTPC", "SFTPB", "SFTPA1", "NAPSA", "LAMP3", "ABCA3"),
    "AT1": ("AGER", "PDPN", "CAV1"),
    "club": ("SCGB1A1", "SCGB3A2", "SCGB3A1"),
    "basal": ("KRT5", "KRT15", "TP63", "NGFR"),
    "ciliated": ("FOXJ1", "TPPP3", "PIFO"),
    # CLDN4 excluded on purpose.
    "barrier_keratin": ("KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR"),
    "malignant_like": ("CEACAM5", "CEACAM6", "MKI67"),
}

FOCAL = ("CLDN4",)
COMPARATOR = ("TACSTD2",)
CONTROLS = ("SFTPC", "KRT8", "SCGB1A1", "KRT5", "EPCAM", "PTPRC")
QC_NEG = ("PTPRC", "PECAM1", "COL1A1")

# Hu et al. 2023: 3 pre-biopsy TN; 12 post-surgery (MPR n=4 incl. pCR P06; NMPR n=8).
PAPER_GROUP = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",  # pCR
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}

MARKER_GENES = sorted(
    {
        *FOCAL,
        *COMPARATOR,
        *CONTROLS,
        *QC_NEG,
        *A3_NORMAL_LUNG,
        *(g for genes in LINEAGES.values() for g in genes),
        *(g for genes in STATES.values() for g in genes),
        "KRT17",
        "CDH1",
        "SFTPA2",
        "CD3G",
        "TRAC",
        "CD4",
        "CD8A",
        "CD8B",
        "KLRD1",
        "KLRF1",
        "NCR1",
        "NCAM1",
        "CSF3R",
    }
)
