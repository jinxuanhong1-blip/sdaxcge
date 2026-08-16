"""Gene panel and shared configuration for analysis A11.

A11: Do the CD47, Galectin, Nectin and TGF-beta immune-evasion axes differ
between TACSTD2 (TROP2)-high and TACSTD2-low lung tumours?

Gene symbols are HGNC-approved; legacy aliases are listed so that older
datasets (which predate the NECTIN/CCN renaming) still map correctly.
"""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "w200" / "A11"

TARGET = "TACSTD2"

# ---------------------------------------------------------------------------
# The four axes named in the question. Each axis holds the ligands/receptors
# that are actually part of that signalling module.
# ---------------------------------------------------------------------------
AXES = {
    "CD47": ["CD47", "SIRPA", "THBS1"],
    "Galectin": ["LGALS1", "LGALS3", "LGALS9", "LGALS8", "LGALS3BP",
                 "HAVCR2", "CD44", "LAG3"],
    "Nectin": ["PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4",
               "TIGIT", "CD226", "CD96"],
    "TGFB": ["TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3",
             "LTBP1", "SERPINE1", "SMAD7", "SKIL", "JUNB"],
}

# Ligands only: the tumour-cell-side molecules that an ADC/immunotherapy
# combination would actually have to engage. Kept separate from receptors,
# which are mostly expressed by immune cells and therefore track TME
# composition rather than tumour biology.
AXIS_LIGANDS = {
    "CD47": ["CD47"],
    "Galectin": ["LGALS1", "LGALS3", "LGALS9", "LGALS3BP"],
    "Nectin": ["PVR", "NECTIN2", "NECTIN4"],
    "TGFB": ["TGFB1", "TGFB2", "TGFB3"],
}

# Pan-fibroblast TGF-beta response signature (Mariathasan et al. 2018,
# Nature 554:544). Read-out of pathway *activity* rather than ligand level.
F_TBRS = ["ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2",
          "CTPS1", "RFLNB", "FSTL3", "HSPB1", "IGFBP3", "JUNB", "LTBP2",
          "MATN2", "MYL9", "MYLK", "NID2", "NOTCH3", "PALLD", "PDLIM7",
          "PMP22", "PPP1R12A", "CAVIN3", "RHOA", "TAGLN", "TGFB1I1",
          "TGFBR2", "TNC", "TPM1", "VIM"]

# Compartment markers - used to test whether any TACSTD2 association is
# simply tumour-content (purity) confounding.
COMPARTMENT = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"],
    "immune": ["PTPRC"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB"],
    "endothelial": ["PECAM1", "VWF", "CLDN5"],
    "myeloid": ["CD68", "CD163", "ITGAM", "LYZ"],
    "Tcell": ["CD3E", "CD3D", "CD2", "CD8A", "FOXP3"],
    "NK": ["NKG7", "KLRD1", "GNLY"],
    "B_plasma": ["MS4A1", "CD79A", "MZB1"],
}

# Contextual checkpoints and lineage markers (reported, not primary).
CONTEXT = ["CD274", "PDCD1LG2", "PDCD1", "CTLA4", "IDO1", "MKI67",
           "NAPSA", "SFTPC", "SFTPB", "TP63", "KRT5", "CEACAM5", "MUC1",
           "ERBB2", "ERBB3", "EGFR"]

# Negative-control genes: no expected biological link to TACSTD2. Used to
# calibrate how much apparent signal the analysis generates by chance /
# by composition alone.
NEG_CONTROL = ["ACTB", "GAPDH", "RPL13A", "TBP", "PPIA", "B2M", "UBC",
               "SDHA", "HPRT1", "PGK1"]

# Legacy symbol -> current symbol, for older annotation sets.
ALIASES = {
    "PVRL1": "NECTIN1", "PVRL2": "NECTIN2", "PVRL3": "NECTIN3",
    "PVRL4": "NECTIN4", "CTGF": "CCN2", "FAM101B": "RFLNB",
    "PRKCDBP": "CAVIN3", "TROP2": "TACSTD2", "CD112": "NECTIN2",
    "CD155": "PVR", "TACSTD1": "EPCAM", "MLLT4": "AFDN",
}


def all_genes():
    """Every gene symbol used anywhere in A11, de-duplicated."""
    g = {TARGET}
    for v in AXES.values():
        g.update(v)
    for v in COMPARTMENT.values():
        g.update(v)
    g.update(F_TBRS)
    g.update(CONTEXT)
    g.update(NEG_CONTROL)
    return sorted(g)


def gene_to_axis():
    """Map each axis gene to its axis name (for annotating result tables)."""
    m = {}
    for axis, genes in AXES.items():
        for g in genes:
            m[g] = axis
    return m


# Primary panel = the genes the question is actually about. FDR is computed
# over this set, so that context/control genes do not dilute the correction.
def primary_panel():
    g = []
    for axis in ["CD47", "Galectin", "Nectin", "TGFB"]:
        g.extend(AXES[axis])
    return sorted(set(g))
