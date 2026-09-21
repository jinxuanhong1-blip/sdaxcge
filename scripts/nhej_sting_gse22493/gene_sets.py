#!/usr/bin/env python3
"""Pre-specified panels for the GSE22493 NHEJ / STING / IFN slice.

Symbols are current HGNC symbols (NCBI gene_info
Symbol_from_nomenclature_authority). Historical names that still appear
on GPL10555 are remapped in run_analysis.py; they are not alternate
panel entries.

NHEJ is classical non-homologous end joining (Ku, DNA-PKcs, XRCC4-LIG4
complex, end processing, gap-fill). It does not include MRN (shared with
homologous recombination), alt-EJ (PARP1, LIG3, XRCC1), or 53BP1-Shieldin
pathway choice.

STING is the cGAS-STING signaling axis, including two negative regulators
(TREX1, ENPP1) and two DNA sensors that feed the same axis (DDX41, IFI16).
The IFN transcriptional program is a separate panel, so a STING shift is
not just an ISG shift. IFNB1 and CXCL10 stay in IFN.

IFN_IMMUNE is the same 73-gene list used in the earlier GSE22493
IFN/MHC-I slice, with one symbol update: DDX58 -> RIGI.

BREAK_MARKERS (H2AX, TP53BP1) are an adjunct. They are not NHEJ enzymes
and are not part of the three-panel score.
"""

from __future__ import annotations

# (symbol, role)
NHEJ = [
    ("XRCC6", "Ku70; DNA end binding"),
    ("XRCC5", "Ku80; DNA end binding"),
    ("PRKDC", "DNA-PKcs"),
    ("XRCC4", "XRCC4-LIG4 complex"),
    ("LIG4", "XRCC4-LIG4 complex"),
    ("NHEJ1", "XLF; XRCC4-LIG4 complex"),
    ("PAXX", "XRCC4-LIG4 complex"),
    ("DCLRE1C", "Artemis; end processing"),
    ("APLF", "end processing"),
    ("PNKP", "end processing"),
    ("APTX", "end processing"),
    ("POLL", "gap-fill polymerase"),
    ("POLM", "gap-fill polymerase"),
]

# (symbol, role). Negative regulators are labelled so a positive log-ratio
# is not read as pathway activation.
STING = [
    ("CGAS", "DNA sensor"),
    ("STING1", "adaptor"),
    ("TBK1", "kinase"),
    ("IKBKE", "kinase"),
    ("IRF3", "transcription factor"),
    ("DDX41", "DNA sensor"),
    ("IFI16", "DNA sensor"),
    ("TREX1", "negative regulator; cytosolic DNA"),
    ("ENPP1", "negative regulator; cGAMP"),
]

# Same list as scripts/w200/C4_GSE22493_ifn/gene_sets.py IFN_IMMUNE,
# except DDX58 is entered as current symbol RIGI.
IFN = [
    "ISG15",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "RSAD2",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFIT5",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "IFI6",
    "IFI27",
    "IFI35",
    "IFI44",
    "IFI44L",
    "IFI16",
    "RIGI",
    "IFIH1",
    "DDX60",
    "XAF1",
    "HERC5",
    "USP18",
    "CMPK2",
    "BST2",
    "GBP1",
    "GBP2",
    "GBP4",
    "GBP5",
    "EIF2AK2",
    "ZBP1",
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF7",
    "IRF9",
    "JAK2",
    "SOCS1",
    "SOCS3",
    "B2M",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "NLRC5",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "HLA-E",
    "HLA-F",
    "TAPBP",
    "CXCL10",
    "CXCL11",
    "CXCL9",
    "CCL5",
    "CCL2",
    "IL6",
    "TNF",
    "NFKB1",
    "RELA",
    "IL15",
    "IL32",
    "IFNB1",
    "IFNG",
    "IFNAR1",
    "IFNAR2",
    "IFNGR1",
    "IFNGR2",
    "IFNL1",
]

# Proximal cGAS-STING transcriptional outputs that live in the IFN panel.
STING_OUTPUTS = ["IFNB1", "CXCL10"]

BREAK_MARKERS = [
    ("H2AX", "DNA-break histone; was H2AFX"),
    ("TP53BP1", "53BP1; pathway choice, not a ligase"),
]

PERTURBATION = ["CLDN4"]

# Platform tokens that must resolve this way. If NCBI drops a synonym,
# the run stops instead of silently calling the gene absent.
REQUIRED_REMAPS = {
    "G22P1": "XRCC6",
    "C6orf150": "CGAS",
    "C2orf13": "APLF",
    "DDX58": "RIGI",
    "G1P2": "ISG15",
    "H2AFX": "H2AX",
    "G1P3": "IFI6",
    "PRKR": "EIF2AK2",
    "cig5": "RSAD2",
    "C1orf29": "IFI44L",
    "NOD27": "NLRC5",
    "NK4": "IL32",
    "HSXIAPAF1": "XAF1",
}
