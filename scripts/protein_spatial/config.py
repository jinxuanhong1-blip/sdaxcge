"""Shared paths and gene identifiers for TROP2/CLDN4 protein+spatial analyses."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = Path("/workspace/data/protein_spatial")
RESULTS = ROOT / "results" / "protein_spatial"
NOTES = ROOT / "notes"
FIGDIR = RESULTS / "figures"

# Official identifiers (do not invent)
TACSTD2_ENSG = "ENSG00000184292"
CLDN4_ENSG = "ENSG00000189143"
TACSTD2_SYMBOL = "TACSTD2"
CLDN4_SYMBOL = "CLDN4"
TACSTD2_UNIPROT = "P09758"
CLDN4_UNIPROT = "O14493"
# Mouse (MGI / UniProt reviewed)
TACSTD2_MOUSE = "Tacstd2"
CLDN4_MOUSE = "Cldn4"
TACSTD2_MOUSE_UNIPROT = "Q8BGV3"
CLDN4_MOUSE_UNIPROT = "O35114"

CPTAC_BASE = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized"
)
CPTAC_COHORTS = ("LUAD", "LSCC")
CPTAC_PROTEIN_SUFFIX = (
    "_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
)

# RNA-derived immune contexture columns from the same freeze phenotype table
IMMUNE_SCORE_COLS = [
    "ESTIMATE_ImmuneScore",
    "ESTIMATE_StromalScore",
    "xCell_immune_score",
    "xCell_T_cell_CD8+",
    "xCell_T_cell_CD4+_(non-regulatory)",
    "xCell_T_cell_regulatory_(Tregs)",
    "xCell_Macrophage_M1",
    "xCell_Macrophage_M2",
    "xCell_NK_cell",
    "CIBERSORT_T_cell_CD8+",
    "CIBERSORT_T_cell_CD4+_memory_activated",
    "CIBERSORT_T_cell_regulatory_(Tregs)",
    "CIBERSORT_Macrophage_M1",
    "CIBERSORT_Macrophage_M2",
    "CIBERSORT_NK_cell_activated",
    "CIBERSORT_B_cell_plasma",
]

# Protein-intrinsic immune markers (Ensembl gene IDs, version stripped)
PROTEIN_IMMUNE_MARKERS = {
    "CD8A": "ENSG00000153563",
    "CD4": "ENSG00000010610",
    "CD3E": "ENSG00000198851",
    "PTPRC": "ENSG00000081237",
    "CD274": "ENSG00000120217",
    "GZMB": "ENSG00000100453",
    "PRF1": "ENSG00000180644",
    "CD68": "ENSG00000129226",
    "MS4A1": "ENSG00000156738",
    "NCAM1": "ENSG00000149294",
}

VISIUM_IMMUNE_GENES = [
    "CD8A",
    "CD4",
    "CD3E",
    "PTPRC",
    "CD68",
    "CD274",
    "GZMB",
    "PRF1",
    "MS4A1",
    "FOXP3",
    "NCAM1",
]

MAX_BYTES = 2 * 1024 ** 3  # skip files >2GB
