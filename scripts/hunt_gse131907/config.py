"""GSE131907 hunt configuration.

This atlas (Kim et al., Nat Commun 2020) is dissociated 10x scRNA-seq.
There is no spatial coordinate. "Neighborhood" and "TLS" are sample-level
co-occurrence / chemokine proxies, not histology or Visium adjacency.
"""

from pathlib import Path

DATA_DIR = Path("/tmp/gse131907")
RESULTS_DIR = Path("/workspace/results/hunt_gse131907")
GENESET_JSON = Path("/workspace/data/genesets/primary_sets.json")

ANNOT_FILE = DATA_DIR / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
LOG2TPM_FILE = DATA_DIR / "GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt.gz"
UMI_FILE = DATA_DIR / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
SERIES_FILE = DATA_DIR / "GSE131907_series_matrix.txt.gz"

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907"
DOWNLOADS = {
    ANNOT_FILE.name: f"{GEO_BASE}/suppl/{ANNOT_FILE.name}",
    LOG2TPM_FILE.name: f"{GEO_BASE}/suppl/{LOG2TPM_FILE.name}",
    UMI_FILE.name: f"{GEO_BASE}/suppl/{UMI_FILE.name}",
    SERIES_FILE.name: f"{GEO_BASE}/matrix/{SERIES_FILE.name}",
}

# Author labels
CD8_SUBTYPES = ["Cytotoxic CD8+ T", "Exhausted CD8+ T", "Naive CD8+ T", "CD8 low T"]
GC_B_SUBTYPES = ["GC B cells in the DZ", "GC B cells in the LZ"]
TUMOR_EPI_SUBTYPES = ["Malignant cells", "tS1", "tS2", "tS3"]
TUMOR_ORIGINS = ["tLung", "tL/B", "mBrain", "mLN", "PE"]
IMMUNE_TYPES = [
    "T lymphocytes",
    "NK cells",
    "B lymphocytes",
    "Myeloid cells",
    "MAST cells",
]

MIN_EPI = 20
MIN_IMMUNE = 20
RANDOM_SEED = 0
GSEA_PERMS = 1000

# Lineage / TLS genes streamed from the matrix (all cells).
LINEAGE = ["EPCAM", "KRT8", "KRT18", "KRT19", "PTPRC", "CD3D", "CD3E"]
CYTO = ["CD8A", "CD8B", "NKG7", "GZMB", "NCAM1", "PRF1"]
TLS_CHEMOKINE = [
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18",
    "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13",
]
TLS_STRUCT = ["CR2", "CXCR5", "MS4A1", "CD79A", "CD79B", "SELL", "LAMP3", "CCR7"]
FOCAL = ["TACSTD2", "CLDN4"]

SELECTED = sorted(set(FOCAL + LINEAGE + CYTO + TLS_CHEMOKINE + TLS_STRUCT))
