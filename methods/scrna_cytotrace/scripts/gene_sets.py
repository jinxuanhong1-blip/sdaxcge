"""Gene sets for CytoTRACE-like stemness / cycling / keratin scoring.

Tirosh S/G2M lists are the Seurat/Scanpy defaults (Tirosh et al. Science 2016).
Keratin / alveolar sets are canonical lung-epithelial differentiation markers.
CytoTRACE2 (Kang et al. 2024) is NOT run; see METHODS.md.
"""

# Tirosh et al. Science 2016 — Seurat cc.genes$s.genes (MLF1IP = CENPU)
S_GENES = [
    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2",
    "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", "CENPU", "HELLS",
    "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2",
    "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6",
    "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1", "CLSPN", "POLA1",
    "CHAF1B", "BRIP1", "E2F8",
]

# Tirosh G2M — Seurat cc.genes$g2m.genes (FAM64A = PIMREG; HN1 = JPT1)
G2M_GENES = [
    "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80",
    "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A",
    "PIMREG", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11",
    "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", "JPT1",
    "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5",
    "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", "ANLN",
    "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5",
    "CENPA",
]

# Simple-epithelium keratins (LUAD-typical)
KERATIN_SIMPLE = ["KRT7", "KRT8", "KRT18", "KRT19"]

# Basal / squamous keratins (more differentiated squamous / basal program)
KERATIN_BASAL = ["KRT5", "KRT6A", "KRT6B", "KRT14", "KRT17"]

KERATIN_ALL = KERATIN_SIMPLE + KERATIN_BASAL

# Differentiated alveolar / club lung epithelium
ALVEOLAR_DIFF = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "NAPSA", "AGER",
    "SCGB1A1", "SCGB3A2",
]

# Sparse cancer-stemness markers (secondary; often dropout-heavy)
STEM_MARKERS = [
    "SOX2", "MYC", "EZH2", "BMI1", "PROM1", "ALDH1A1", "POU5F1", "NANOG",
    "CD44", "ALDH1A3",
]

# Lineage for GSE207422 malignant-like call (no author CopyKAT on GEO)
LINEAGE = {
    "Epithelial": [
        "EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "KRT17",
        "ELF3", "CDH1", "MUC1",
    ],
    "T/NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1", "IGHM"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "FCN1", "ITGAX"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5", "RAMP2"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN"],
}

NORMAL_LUNG = [
    "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SFTPD", "AGER", "NAPSA",
    "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS",
]

FOCAL = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC", "MKI67"]


def all_panel_genes() -> list[str]:
    genes: list[str] = []
    seen: set[str] = set()
    for block in (
        FOCAL,
        S_GENES,
        G2M_GENES,
        KERATIN_ALL,
        ALVEOLAR_DIFF,
        STEM_MARKERS,
        NORMAL_LUNG,
        *[LINEAGE[k] for k in LINEAGE],
    ):
        for g in block:
            if g not in seen:
                seen.add(g)
                genes.append(g)
    return genes
