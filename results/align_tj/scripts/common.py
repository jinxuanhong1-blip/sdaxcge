"""Shared helpers for the USER-ALIGN TROP2 tight-junction analysis.

All outputs live under results/align_tj/.
"""
import os
import gzip
import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(BASE, "data")
TABLES = os.path.join(BASE, "tables")
FIGURES = os.path.join(BASE, "figures")
for d in (DATA, TABLES, FIGURES):
    os.makedirs(d, exist_ok=True)

# The user's private hypothesis gene sets (what we are trying to align public data to)
USER_TJ_GENES = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]
USER_TFS = ["ELF3", "GRHL1", "KLF4", "TFAP2A"]
TROP2 = "TACSTD2"

# Themes the user expects enriched in TROP2-high tumors
THEME_KEYWORDS = {
    "keratinization": ["keratin"],
    "skin_barrier": ["skin", "epiderm", "cornif"],
    "tight_junction": ["tight junction", "apical junction", "cell-cell junction",
                        "bicellular tight junction", "cell junction"],
    "emt": ["epithelial mesenchymal transition", "epithelial-mesenchymal"],
}


def read_xena_matrix(path):
    """Read a UCSC Xena gene x sample matrix (gene symbols in first column)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index.name = "gene"
    # collapse duplicate gene rows by mean if any
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    return df


def tumor_samples(cols):
    """Keep TCGA primary tumor samples (barcode sample-type code 01)."""
    keep = []
    for c in cols:
        parts = c.split("-")
        if len(parts) >= 4:
            code = parts[3][:2]
            if code == "01":  # primary solid tumor
                keep.append(c)
    return keep
