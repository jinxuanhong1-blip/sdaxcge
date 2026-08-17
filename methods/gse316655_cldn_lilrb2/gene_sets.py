"""Compact human gene sets for GSE316655 (GRCh38 names after prefix strip).

Sets are short, a priori, and documented. Missing genes are dropped at runtime
and recorded — they are not imputed.
"""

from __future__ import annotations

# Tight-junction claudins (protein-coding CLDN1–25 + CLDN34). CLDND* excluded.
CLAUDINS = [
    "CLDN1",
    "CLDN2",
    "CLDN3",
    "CLDN4",
    "CLDN5",
    "CLDN6",
    "CLDN7",
    "CLDN8",
    "CLDN9",
    "CLDN10",
    "CLDN11",
    "CLDN12",
    "CLDN14",
    "CLDN15",
    "CLDN16",
    "CLDN17",
    "CLDN18",
    "CLDN19",
    "CLDN20",
    "CLDN22",
    "CLDN23",
    "CLDN24",
    "CLDN25",
    "CLDN34",
]

# Paper axis + common barrier partners (not treated as LILRB ligands).
BARRIER = ["TACSTD2", "EPCAM", "CDH1", "F11R", "OCLN", "TJP1"]

LILR = [
    "LILRA1",
    "LILRA2",
    "LILRA4",
    "LILRA5",
    "LILRA6",
    "LILRB1",
    "LILRB2",
    "LILRB3",
    "LILRB4",
    "LILRB5",
]

LINEAGE = {
    "myeloid": ["LYZ", "CD14", "CD68", "CSF1R", "FCGR3A", "S100A8", "S100A9", "CST3"],
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "FGFBP2"],
    "B": ["CD79A", "MS4A1", "CD19", "IGHM"],
    "pDC": ["IL3RA", "CLEC4C", "LILRA4"],
    "melanoma": ["MLANA", "PMEL", "TYR", "MITF", "DCT"],
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
}

# Compact NF-κB / inflammatory myeloid program (not full Hallmark).
NFKB = [
    "NFKB1",
    "NFKB2",
    "RELA",
    "RELB",
    "REL",
    "NFKBIA",
    "NFKBIZ",
    "TNF",
    "IL1B",
    "CXCL8",
    "CXCL2",
    "CCL3",
    "CCL4",
    "TNFAIP3",
    "ICAM1",
    "SOD2",
    "BIRC3",
    "PTGS2",
]

# Compact JAK/STAT (STAT1 + STAT3 axes).
STAT = [
    "STAT1",
    "STAT2",
    "STAT3",
    "STAT4",
    "STAT5A",
    "STAT5B",
    "STAT6",
    "JAK1",
    "JAK2",
    "JAK3",
    "SOCS1",
    "SOCS3",
    "IRF1",
    "IRF9",
    "IL6ST",
    "IL6",
]

# MHC class I antigen-presentation cassette.
MHCI = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "NLRC5"]

MHCII = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74"]

CYTOTOX = ["GZMA", "GZMB", "GZMH", "GZMK", "PRF1", "GNLY", "NKG7", "IFNG", "FASLG"]

EXHAUST = ["PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX", "TOX2"]

SUPPRESS_MYELOID = ["ARG1", "IDO1", "IL10", "TGFB1", "CD274", "SIRPA", "VSIR"]

# Cheap ligand–receptor pairs (curated; not a full CellPhoneDB run).
# ligand, receptor, note
LR_PAIRS = [
    ("CLDN4", "LILRB2", "paper CLDN–LILRB2"),
    ("CLDN4", "LILRB5", "paper CLDN–LILRB5"),
    ("CLDN1", "LILRB2", "other claudin–LILRB2"),
    ("CLDN3", "LILRB2", "other claudin–LILRB2"),
    ("CLDN7", "LILRB2", "other claudin–LILRB2"),
    ("CLDN18", "LILRB2", "paper CLDN18.2–LILRB2"),
    ("CLDN18", "LILRB5", "paper CLDN18–LILRB5"),
    ("HLA-A", "LILRB1", "classical MHC-I–LILRB1"),
    ("HLA-B", "LILRB1", "classical MHC-I–LILRB1"),
    ("HLA-C", "LILRB1", "classical MHC-I–LILRB1"),
    ("HLA-A", "LILRB2", "classical MHC-I–LILRB2"),
    ("HLA-B", "LILRB2", "classical MHC-I–LILRB2"),
    ("HLA-C", "LILRB2", "classical MHC-I–LILRB2"),
    ("HLA-G", "LILRB1", "nonclassical MHC-I–LILRB1"),
    ("HLA-G", "LILRB2", "nonclassical MHC-I–LILRB2"),
    ("HLA-E", "LILRB1", "HLA-E–LILRB1"),
    ("HLA-E", "LILRB2", "HLA-E–LILRB2"),
    ("HLA-F", "LILRB1", "HLA-F–LILRB1"),
    ("HLA-F", "LILRB2", "HLA-F–LILRB2"),
    ("ANGPTL2", "LILRB2", "known LILRB2 ligand"),
    ("CD274", "PDCD1", "PD-L1–PD-1"),
    ("CD47", "SIRPA", "don't-eat-me"),
    ("HLA-E", "KLRC1", "NKG2A"),
    ("LGALS9", "HAVCR2", "galectin-9–TIM3"),
    ("NECTIN2", "TIGIT", "nectin–TIGIT"),
    ("PVR", "TIGIT", "PVR–TIGIT"),
]

PANEL_EXTRA = [
    "PTPRC",
    "CD45",
    "CD4",
    "CD8A",
    "CD8B",
    "FOXP3",
    "IL2RA",
    "NCAM1",
    "ITGAM",
    "ITGAX",
    "MRC1",
    "CD163",
    "C1QA",
    "C1QB",
    "MSR1",
    "FCGR1A",
    "HLA-G",
    "HLA-E",
    "HLA-F",
    "ANGPTL2",
    "CD47",
    "SIRPA",
    "KLRC1",
    "LGALS9",
    "NECTIN2",
    "PVR",
    "CD80",
    "CD86",
    "CXCL9",
    "CXCL10",
    "CCL5",
    "IFIT1",
    "ISG15",
    "MX1",
]


def all_panel_genes() -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    chunks = [
        CLAUDINS,
        BARRIER,
        LILR,
        NFKB,
        STAT,
        MHCI,
        MHCII,
        CYTOTOX,
        EXHAUST,
        SUPPRESS_MYELOID,
        PANEL_EXTRA,
    ]
    for block in chunks:
        for g in block:
            if g not in seen:
                seen.add(g)
                out.append(g)
    for genes in LINEAGE.values():
        for g in genes:
            if g not in seen:
                seen.add(g)
                out.append(g)
    for lig, rec, _ in LR_PAIRS:
        for g in (lig, rec):
            if g not in seen:
                seen.add(g)
                out.append(g)
    return out
