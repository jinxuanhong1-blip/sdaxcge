"""Frozen, documented gene sets for B4 TJ recompute.

Nothing here was chosen to reproduce p=0.019. Two public TJ collections
are used as written by their source databases.

Primary documented TJ score
---------------------------
Reactome R-HSA-420029 "Tight junction interactions"
MSigDB: REACTOME_TIGHT_JUNCTION_INTERACTIONS (M17967)
https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/REACTOME_TIGHT_JUNCTION_INTERACTIONS
https://reactome.org/content/detail/R-HSA-420029

This is the Reactome *tight-junction interactions* pathway (claudins +
CRB3/PALS1/PATJ + Par3/Par6/aPKC + JAM-A/F11R). It is not the parent
pathway R-HSA-421270 (cell-cell junction organization), which also
contains adherens junctions and desmosomes.

Secondary documented TJ score
-----------------------------
KEGG pathway hsa04530 "Tight junction"
https://www.kegg.jp/pathway/hsa04530
MSigDB analog: KEGG_TIGHT_JUNCTION

KEGG's set is broader than structural TJ proteins (actin/ARP2/3, myosins,
tubulins, kinases). It is reported because it is the standard named
"tight junction" pathway, not because it is a clean barrier signature.

Aliases
-------
GRCh37 / GENCODE 19 (GSE126044) uses older symbols for two Reactome genes:
  PALS1 -> MPP5
  PATJ  -> INADL
"""

from __future__ import annotations

REACTOME_TJ_ID = "R-HSA-420029"
REACTOME_TJ_MSIGDB = "REACTOME_TIGHT_JUNCTION_INTERACTIONS"
REACTOME_TJ_MSIGDB_URL = (
    "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/"
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS"
)
REACTOME_TJ_REACTOME_URL = "https://reactome.org/content/detail/R-HSA-420029"

# Official MSigDB geneSymbols for REACTOME_TIGHT_JUNCTION_INTERACTIONS
# (downloaded 2026-08-16 from MSigDB JSON endpoint; 30 unique symbols).
REACTOME_TJ_GENES = [
    "CLDN1",
    "CLDN2",
    "CLDN3",
    "CLDN4",
    "CLDN5",
    "CLDN6",
    "CLDN7",
    "CLDN8",
    "CLDN9",
    "CLDN10",
    "CLDN11",
    "CLDN12",
    "CLDN14",
    "CLDN15",
    "CLDN16",
    "CLDN17",
    "CLDN18",
    "CLDN19",
    "CLDN20",
    "CLDN22",
    "CLDN23",
    "CRB3",
    "F11R",
    "PALS1",
    "PARD3",
    "PARD6A",
    "PARD6B",
    "PARD6G",
    "PATJ",
    "PRKCI",
]

KEGG_TJ_ID = "hsa04530"
KEGG_TJ_URL = "https://www.kegg.jp/pathway/hsa04530"
KEGG_REST_URL = "https://rest.kegg.jp/get/hsa04530"

# Symbol used in a matrix -> official / set symbol
SYMBOL_ALIASES = {
    "MPP5": "PALS1",
    "INADL": "PATJ",
    "JAM1": "F11R",
    "SLC9A3R1": "NHERF1",
}

# Official symbol -> older symbols that may appear in GRCh37 matrices
SYMBOL_FALLBACKS = {
    "PALS1": ["MPP5"],
    "PATJ": ["INADL"],
    "F11R": ["JAM1"],
    "NHERF1": ["SLC9A3R1"],
}

CONTROLS = {
    "CLDN4": ["CLDN4"],
    "CD8A": ["CD8A"],
    "CD8": ["CD8A", "CD8B"],
    "TACSTD2": ["TACSTD2"],
}

# GRCh37 / GENCODE 19 Ensembl IDs. mygene.info current IDs sometimes point at
# post-hg19 replacements (e.g. CLDN7 ENSG00000288292) that are absent from
# GSE135222 (hg19). These are the IDs used on that matrix.
GRCH37_ENSEMBL = {
    "CLDN1": "ENSG00000163347",
    "CLDN2": "ENSG00000165376",
    "CLDN3": "ENSG00000165215",
    "CLDN4": "ENSG00000189143",
    "CLDN5": "ENSG00000184113",
    "CLDN6": "ENSG00000184697",
    "CLDN7": "ENSG00000181885",
    "CLDN8": "ENSG00000156284",
    "CLDN9": "ENSG00000213937",
    "CLDN10": "ENSG00000134873",
    "CLDN11": "ENSG00000013297",
    "CLDN12": "ENSG00000157224",
    "CLDN14": "ENSG00000159261",
    "CLDN15": "ENSG00000106404",
    "CLDN16": "ENSG00000113946",
    "CLDN17": "ENSG00000156282",
    "CLDN18": "ENSG00000066405",
    "CLDN19": "ENSG00000164007",
    "CLDN20": "ENSG00000171217",
    "CLDN22": "ENSG00000177300",
    "CLDN23": "ENSG00000198642",
    "CRB3": "ENSG00000130545",
    "F11R": "ENSG00000158769",
    "PALS1": "ENSG00000072415",
    "MPP5": "ENSG00000072415",
    "PARD3": "ENSG00000148498",
    "PARD6A": "ENSG00000102981",
    "PARD6B": "ENSG00000124171",
    "PARD6G": "ENSG00000178184",
    "PATJ": "ENSG00000132849",
    "INADL": "ENSG00000132849",
    "PRKCI": "ENSG00000163558",
    "CD8A": "ENSG00000153563",
    "CD8B": "ENSG00000172116",
    "TACSTD2": "ENSG00000184292",
    "CTTN": "ENSG00000085795",
    "SCRIB": "ENSG00000180900",
    "ACTG1": "ENSG00000184009",
    "ACTN4": "ENSG00000130402",
}

# Extra historical / alternate Ensembl IDs seen on hg19 matrices.
GRCH37_ENSEMBL_ALTS = {
    "CLDN7": ["ENSG00000181885", "ENSG00000288292"],
    "CLDN23": ["ENSG00000253958", "ENSG00000198642"],
    "CTTN": ["ENSG00000085795", "ENSG00000143376"],
    "ACTG1": ["ENSG00000184009", "ENSG00000291420"],
    "ACTN4": ["ENSG00000130402", "ENSG00000282844"],
    "SCRIB": ["ENSG00000180900", "ENSG00000274287"],
}
