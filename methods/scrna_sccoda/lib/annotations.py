"""Lineage markers, collapses, and TACSTD2 / TLS definitions.

These maps are the combinatorial annotation axis. Every collapse is
declared here so the grid cannot silently invent a cell-type name.
"""

from __future__ import annotations

# Hu et al. Genome Medicine 2023 canonical markers (GSE207422 paper).
HU_LINEAGES: dict[str, list[str]] = {
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

NORMAL_LUNG = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]

MARKER_PANEL = sorted(
    {
        "TACSTD2",
        "CLDN4",
        "CXCL13",
        "CCL19",
        "CCL21",
        "MS4A1",
        "CD79A",
        *NORMAL_LUNG,
        *(g for genes in HU_LINEAGES.values() for g in genes),
    }
)

# Author major labels in GSE241934 (NEOTIDE).
AUTHOR_MAJOR_CANON = {
    "T": "T",
    "NK": "NK",
    "B": "B",
    "Myeloid": "myeloid",
    "Epi": "epithelial",
    "Fibro": "fibroblast",
    "Endo": "endothelial",
    "Mast": "mast",
}

# Fine labels that are TLS-proximal in dissociated scRNA (B + plasma + Tfh).
TLS_FINE_PREFIXES = (
    "B_",
    "B_naive",
    "Plasma",
    "CD4T_Tfh_",
)

# Primary compositional endpoints. Names must match a collapse output.
PRIMARY_COMPARTMENTS = ("TNK", "TLS")

# How each annotation scheme collapses raw labels into testable compartments.
# Keys are raw labels (after that scheme's naming); values are collapsed names.
# Unlisted labels stay as themselves (author_major) or map to Other.

COLLAPSE = {
    "author_major": {
        # identity plus TNK / TLS rollups added as extra schemes
        "T": "T",
        "NK": "NK",
        "B": "B",
        "myeloid": "myeloid",
        "epithelial": "epithelial",
        "fibroblast": "fibroblast",
        "endothelial": "endothelial",
        "mast": "mast",
    },
    "author_collapsed": {
        "T": "TNK",
        "NK": "TNK",
        "B": "TLS",
        "myeloid": "myeloid",
        "epithelial": "epithelial",
        "fibroblast": "stromal",
        "endothelial": "stromal",
        "mast": "mast",
    },
    "hu_markers": {
        "T": "T",
        "NK": "NK",
        "B": "B",
        "plasma": "plasma",
        "epithelial": "epithelial",
        "myeloid": "myeloid",
        "neutrophil": "neutrophil",
        "fibroblast": "fibroblast",
        "endothelial": "endothelial",
        "mast": "mast",
        "other": "other",
    },
    "hu_collapsed": {
        "T": "TNK",
        "NK": "TNK",
        "B": "TLS",
        "plasma": "TLS",
        "epithelial": "epithelial",
        "myeloid": "myeloid",
        "neutrophil": "myeloid",
        "fibroblast": "stromal",
        "endothelial": "stromal",
        "mast": "mast",
        "other": "other",
    },
    "drmref": {
        "CD8+ T cells": "T",
        "CD4+ T cells": "T",
        "NK cells": "NK",
        "B cells": "B",
        "Plasma cells": "plasma",
        "Malignant cells": "epithelial",
        "Mono/Macro": "myeloid",
        "Neutrophils": "neutrophil",
        "pDCs": "myeloid",
        "Mast cells": "mast",
        "Fibroblasts": "fibroblast",
        "Endothelial cells": "endothelial",
    },
    "drmref_collapsed": {
        "CD8+ T cells": "TNK",
        "CD4+ T cells": "TNK",
        "NK cells": "TNK",
        "B cells": "TLS",
        "Plasma cells": "TLS",
        "Malignant cells": "epithelial",
        "Mono/Macro": "myeloid",
        "Neutrophils": "myeloid",
        "pDCs": "myeloid",
        "Mast cells": "mast",
        "Fibroblasts": "stromal",
        "Endothelial cells": "stromal",
    },
    "marker_coarse": {
        "TNK": "TNK",
        "TLS": "TLS",
        "epithelial": "epithelial",
        "other": "other",
    },
    "author_garnett_limited": {
        "TNK": "TNK",
        "epithelial": "epithelial",
        "other": "other",
    },
}


def collapse_label(scheme: str, raw: str) -> str:
    mapping = COLLAPSE[scheme]
    return mapping.get(raw, "other")


def is_tls_fine(label: str) -> bool:
    if label in {"Plasma", "B_naive"} or label.startswith("B_"):
        return True
    if label.startswith("CD4T_Tfh_"):
        return True
    return False


def is_tnk_fine(label: str) -> bool:
    if label.startswith(("CD4T_", "CD8T_", "NK_", "T_prolifer", "ILC")):
        return True
    if label in {"T", "NK"}:
        return True
    return False
