"""Documented gene sets for B4 extras / GSE126044 index scoring.

The 7-gene module is CLDN1, CLDN4, CLDN7, F11R, TJP1, TJP2, OCLN.
The 5-gene module is CLDN1, CLDN4, CLDN7, F11R, PARD3.
Related NSCLC barrier lists and public TJ pathway sets are retained as comparators.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# User-requested modules (item 6)
# ---------------------------------------------------------------------------
CLDN4_ALONE = ["CLDN4"]

# Explicit wave-2 request: CLDN1/4/7/F11R/PARD3 mean-z
CLDN147_F11R_PARD3 = ["CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3"]

# Compact basal/squamous keratin module used as the keratin contrast
# (not from Cho 2020). These are the lung-relevant KRTs in
# REACTOME_KERATINIZATION / GO keratinization.
KRT_BASAL_SQUAMOUS = [
    "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT14", "KRT15", "KRT16", "KRT17", "KRT19",
]

# ---------------------------------------------------------------------------
# Author / paper / database TJ lists (item 1)
# ---------------------------------------------------------------------------
# Chae et al. 2018 Sci Rep 8:1023. Table 2 Cellular Barrier Molecule genes,
# tight-junction class: CLDN1/5/7, JAM1/2 (F11R/JAM2), TJP1/2.
# Related NSCLC immune-infiltration paper (not ICI-response, not Cho 2020).
CHAE2018_CBM_TJ = ["CLDN1", "CLDN5", "CLDN7", "F11R", "JAM2", "TJP1", "TJP2"]

# User PPT / claim-page 7-gene module named on prior B3/B4 claim pages.
# Not published by Cho 2020. Included because that is the list the PPT
# analysis is believed to have used (PR #90 reported MW p=0.019 for this
# z-mean). Recomputed independently here; not tuned.
USER_PPT_7GENE = ["CLDN1", "CLDN4", "CLDN7", "F11R", "TJP1", "TJP2", "OCLN"]

# Repo 15-gene structural TJ (PR #69).
TJ_15 = [
    "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    "OCLN",
    "TJP1", "TJP2", "TJP3",
    "F11R", "JAM2", "JAM3",
    "MARVELD2", "MARVELD3",
    "CGN", "CGNL1",
]

# Frozen MSigDB REACTOME_TIGHT_JUNCTION_INTERACTIONS (R-HSA-420029), 30 genes.
# Official symbols; GRCh37 aliases applied at match time (PALS1->MPP5, PATJ->INADL).
REACTOME_TJ = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8",
    "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN16",
    "CLDN17", "CLDN18", "CLDN19", "CLDN20", "CLDN22", "CLDN23",
    "CRB3", "F11R", "PALS1", "PARD3", "PARD6A", "PARD6B", "PARD6G",
    "PATJ", "PRKCI",
]

# ---------------------------------------------------------------------------
# Additional TJ comparators (pathway / prior-PR lists)
# ---------------------------------------------------------------------------
CUSTOM_TJ_CORE = [
    "AMOT", "AMOTL1", "CGN", "CGNL1", "CLDN1", "CLDN3", "CLDN4", "CLDN7",
    "CLDN8", "CRB3", "F11R", "INADL", "JAM2", "JAM3", "MAGI1", "MAGI3",
    "MARVELD2", "MPDZ", "OCLN", "PARD3", "PARD6A", "PARD6B", "PATJ",
    "TJP1", "TJP2", "TJP3",
]

CORE_TJ_40 = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8",
    "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15", "CLDN16",
    "CLDN17", "CLDN18", "CLDN19", "CLDN20", "CLDN23",
    "OCLN", "MARVELD2", "MARVELD3",
    "TJP1", "TJP2", "TJP3", "SYMPK", "CGN", "CGNL1",
    "F11R", "JAM2", "JAM3",
    "PARD3", "PARD6A", "PARD6B", "PRKCZ", "PRKCI",
    "INADL", "PATJ", "MPP5", "CRB3", "AMOT",
]

# Positive-control immune genes (labels should separate R vs NR).
CD8A_ALONE = ["CD8A"]
CD8 = ["CD8A", "CD8B"]

# GRCh37 / older-symbol aliases present in GSE126044 FeatureCounts (GENCODE 19).
ALIASES = {
    "PALS1": ["MPP5"],
    "PATJ": ["INADL"],
    "FYB1": ["FYB"],
    "F11R": ["JAM1", "JAMA"],
    "HLA-F": ["HLAF"],
    "HLA-G": ["HLAG"],
}


def parse_kegg_genes(text: str) -> list[str]:
    import re

    genes: list[str] = []
    in_gene = False
    for line in text.splitlines():
        if line.startswith("GENE"):
            in_gene = True
            rest = line[4:].strip()
        elif in_gene and (line.startswith(" ") or line.startswith("\t")):
            rest = line.strip()
        elif in_gene:
            break
        else:
            continue
        m = re.match(r"(\d+)\s+([A-Za-z0-9-]+)", rest)
        if m:
            genes.append(m.group(2))
    seen: set[str] = set()
    out: list[str] = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            out.append(g)
    return out
