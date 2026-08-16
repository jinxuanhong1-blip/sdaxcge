"""Public treatment-naive LUAD/NSCLC scRNA joint-object settings.

This is extra atlas n (epithelial TACSTD2/CLDN4 vs T/NK and B/TLS-like
fractions). It is not an ICI / MPR / response test.
"""

from __future__ import annotations

from pathlib import Path

ROOT = Path("/workspace")
DATA = Path("/tmp/scrna_harmony_naive")
EXTRACTED = DATA / "extracted"
RESULTS = ROOT / "results" / "scrna_harmony_naive"

# Tumor sites in Kim et al. (exclude nLung / nLN).
GSE131907_TUMOR_ORIGINS = ("tLung", "tL/B", "mLN", "PE", "mBrain")

# Leader/GSE154826 is CD45+ CITE-seq; epithelium is the authors'
# epi_endo_fibro_doublet leak (two digest patients only). Not a whole-tumor
# epithelial atlas — omitted from the joint object.
EXCLUDED = {
    "GSE154826": (
        "Leader et al. Cancer Cell 2021. CD45-bead CITE-seq. Epithelium is "
        "restricted to the epi_endo_fibro_doublet gate; Trop-2 ADT is not on "
        "the panel. Omitted because a joint epithelial score would be leak, "
        "not a tumor-epithelial compartment."
    )
}

TARGETS = ["TACSTD2", "CLDN4"]

LINEAGES = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "T": ["CD3D", "CD3E", "CD2", "CD8A", "CD4"],
    "NK": ["NKG7", "GNLY", "FGFBP2", "KLRD1", "KLRF1"],
    "B": ["CD79A", "MS4A1", "CD19", "CD79B"],
    "plasma": ["MZB1", "JCHAIN", "XBP1"],
    "myeloid": ["LYZ", "CD68", "CD14", "FCGR3A"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN"],
    "endothelial": ["VWF", "PECAM1", "RAMP2"],
}

NORMAL_LUNG = ["SFTPA2", "SFTPC", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3", "FOXJ1"]
TLS_CHEMOKINE = ["CXCL13", "CCL19", "CCL21", "CXCL9", "CXCL10"]
EXTRA = [
    "PTPRC",
    "MKI67",
    "STMN1",
    "NKX2-1",
    "NAPSA",
    "KRT5",
    "ACTA2",
    "FAP",
    "PDPN",
    "NCAM1",
    "PRF1",
    "GZMB",
    "CR2",
    "CXCR5",
    "LTB",
    "SFTPA1",
]

GENE_PANEL = sorted(
    set(TARGETS)
    | {g for genes in LINEAGES.values() for g in genes}
    | set(NORMAL_LUNG)
    | set(TLS_CHEMOKINE)
    | set(EXTRA)
)

MIN_EPI = 20
MIN_TNK = 20
MIN_B = 10
HARMONY_MAX_PER_DONOR = 2500
HARMONY_NPC = 20
RANDOM_SEED = 0

URLS = {
    "GSE131907_annot": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ),
    "GSE131907_umi": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/"
        "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ),
    "GSE148071_tar": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/suppl/"
        "GSE148071_RAW.tar"
    ),
    "GSE148071_series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148071/matrix/"
        "GSE148071_series_matrix.txt.gz"
    ),
    "GSE253013_rds": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/suppl/"
        "GSE253013_all_luad_garnett_temp.rds.gz"
    ),
    "GSE253013_series": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253013/matrix/"
        "GSE253013_series_matrix.txt.gz"
    ),
    "GSE127465_meta": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
        "GSE127465_human_cell_metadata_54773x25.tsv.gz"
    ),
    "GSE127465_genes": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
        "GSE127465_gene_names_human_41861.tsv.gz"
    ),
    "GSE127465_mtx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE127nnn/GSE127465/suppl/"
        "GSE127465_human_counts_normalized_54773x41861.mtx.gz"
    ),
}

FILES = {
    "GSE131907_annot": DATA / "GSE131907_Lung_Cancer_cell_annotation.txt.gz",
    "GSE131907_umi": DATA / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
    "GSE148071_tar": DATA / "GSE148071_RAW.tar",
    "GSE148071_series": DATA / "GSE148071_series_matrix.txt.gz",
    "GSE253013_rds": DATA / "GSE253013_all_luad_garnett_temp.rds.gz",
    "GSE253013_series": DATA / "GSE253013_series_matrix.txt.gz",
    "GSE127465_meta": DATA / "GSE127465_human_cell_metadata_54773x25.tsv.gz",
    "GSE127465_genes": DATA / "GSE127465_gene_names_human_41861.tsv.gz",
    "GSE127465_mtx": DATA / "GSE127465_human_counts_normalized_54773x41861.mtx.gz",
}
