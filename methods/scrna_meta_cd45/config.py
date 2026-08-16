"""Paths and public accessions for the CD45+ / immune-library extra table.

This analysis is immune-compartment only. It is not malignant-cell TACSTD2/CLDN4.
User A3 (epithelial / CopyKAT GSE207422) is taken as given and is not re-run.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
RESULTS = ROOT / "results" / "extra" / "scrna_meta_cd45"
METHODS = Path(__file__).resolve().parent

GSE243013_DIR = DATA / "gse243013"
GSE229353_DIR = DATA / "gse229353"
GSE154826_DIR = DATA / "gse154826"

GSE243013_META = GSE243013_DIR / "GSE243013_NSCLC_immune_scRNA_metadata.csv.gz"
GSE243013_GENES = GSE243013_DIR / "GSE243013_genes.csv.gz"
GSE243013_EXTRACT = GSE243013_DIR / "extracted_tacstd2_cldn4_cells.FROM_PR142.tsv"
GSE243013_MTX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE243nnn/GSE243013/suppl/"
    "GSE243013_NSCLC_immune_scRNA_counts.mtx.gz"
)
GSE243013_MTX_BYTES = 7_123_039_063

GSE229353_RAW = GSE229353_DIR / "raw"
GSE229353_SUPP_TXT = DATA / "gse229353" / "s41698_supp_table1.txt"

# Paper Supplementary Table 1 (Hui et al., npj Precis Oncol 2023; public PDF).
# RVT% and MPR are not on the GEO series matrix.
GSE229353_TABLE_S1 = {
    "P01": {"age": 68, "sex": "Male", "therapy": "NAC", "author_n_cells": 4057,
            "pathology": "SCC", "rvt_pct": 40, "pathologic_response": "non-MPR"},
    "P02": {"age": 69, "sex": "Male", "therapy": "NAPC", "author_n_cells": 3045,
            "pathology": "SCC", "rvt_pct": 90, "pathologic_response": "non-MPR"},
    "P03": {"age": 51, "sex": "Female", "therapy": "NAPC", "author_n_cells": 2182,
            "pathology": "SCC", "rvt_pct": 0, "pathologic_response": "MPR"},
    "P04": {"age": 66, "sex": "Female", "therapy": "NAPC", "author_n_cells": 2888,
            "pathology": "AD", "rvt_pct": 80, "pathologic_response": "non-MPR"},
    "P05": {"age": 67, "sex": "Male", "therapy": "NAPC", "author_n_cells": 4474,
            "pathology": "SCC", "rvt_pct": 0, "pathologic_response": "MPR"},
    "P06": {"age": 62, "sex": "Male", "therapy": "NAPC", "author_n_cells": 7579,
            "pathology": "SCC", "rvt_pct": 0, "pathologic_response": "MPR"},
    "P07": {"age": 46, "sex": "Male", "therapy": "NAPC", "author_n_cells": 2636,
            "pathology": "AD", "rvt_pct": 60, "pathologic_response": "non-MPR"},
}

GENES = ("TACSTD2", "CLDN4")
