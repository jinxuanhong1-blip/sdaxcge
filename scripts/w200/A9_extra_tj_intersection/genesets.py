"""Pre-specified tight-junction catalog for the A9 extra-cohort supplement.

This is NOT the User A9 slide list (CLDN1/4/7, F11R, PARD3). That five-gene
intersection is taken as given and is never used as an inclusion filter here.

Primary universe = published TJ sets already used in the sibling A8/A9 GSEA
GMT (Enrichr GO_Cellular_Component_2021 + Reactome_2022), plus one documented
canonical JAM (JAM2) that is in KEGG Tight junction but missing from those two
sets.

  GOCC_BICELLULAR_TIGHT_JUNCTION  GO:0005923
  REACTOME_TIGHT_JUNCTION_INTERACTIONS  R-HSA-420029
  CANONICAL_ADDITIONS  JAM2 only

KEGG Tight junction is intentionally NOT used as the universe: it is dominated
by actin, myosin, MAPK, and tubulin genes that are not junction strand proteins.
Those genes are not tested and are not counted in the overlap.

TACSTD2 is the splitter, not a tested TJ gene.
"""

from __future__ import annotations

# User A9 slide intersection — context annotation only. Not the test set.
A9_SLIDE_GENES = ("CLDN1", "CLDN4", "CLDN7", "F11R", "PARD3")

TARGET = "TACSTD2"

# GO_Cellular_Component_2021: bicellular tight junction (GO:0005923)
# Source: data/genesets/claim_A8A9_genesets.gmt on cursor/claim-a8a9-trop2-tj-gsea-2be8
GOCC_BICELLULAR_TIGHT_JUNCTION = (
    "AMOT",
    "AMOTL1",
    "AMOTL2",
    "ANK3",
    "AOC1",
    "APC",
    "BVES",
    "CDH5",
    "CGN",
    "CLDN1",
    "CLDN10",
    "CLDN11",
    "CLDN12",
    "CLDN14",
    "CLDN15",
    "CLDN16",
    "CLDN17",
    "CLDN18",
    "CLDN19",
    "CLDN2",
    "CLDN20",
    "CLDN22",
    "CLDN23",
    "CLDN24",
    "CLDN25",
    "CLDN3",
    "CLDN34",
    "CLDN4",
    "CLDN5",
    "CLDN6",
    "CLDN7",
    "CLDN8",
    "CLDN9",
    "CLMP",
    "CXADR",
    "CYTH1",
    "CYTH3",
    "DDX58",
    "DLG1",
    "ECT2",
    "EPCAM",
    "EPPK1",
    "F11R",
    "FRMD4A",
    "FRMD4B",
    "FRMPD2",
    "IGSF5",
    "JAM3",
    "JAML",
    "MAGI2",
    "MARVELD2",
    "MARVELD3",
    "MICALL2",
    "MPP7",
    "MTDH",
    "OCLN",
    "PARD3",
    "PARD6A",
    "PARD6B",
    "PATJ",
    "POF1B",
    "RAB13",
    "RAP2B",
    "RAP2C",
    "RAPGEF2",
    "SH3BP1",
    "SHROOM2",
    "STRN",
    "TBCD",
    "TGFBR1",
    "TJP1",
    "TJP2",
    "TJP3",
    "UBN1",
    "VAPA",
    "WNK3",
    "WNK4",
    "YBX3",
)

# Reactome_2022: Tight Junction Interactions R-HSA-420029
REACTOME_TIGHT_JUNCTION_INTERACTIONS = (
    "CLDN1",
    "CLDN10",
    "CLDN11",
    "CLDN12",
    "CLDN14",
    "CLDN15",
    "CLDN16",
    "CLDN17",
    "CLDN18",
    "CLDN19",
    "CLDN2",
    "CLDN20",
    "CLDN22",
    "CLDN23",
    "CLDN3",
    "CLDN4",
    "CLDN5",
    "CLDN6",
    "CLDN7",
    "CLDN8",
    "CLDN9",
    "CRB3",
    "F11R",
    "PALS1",
    "PARD3",
    "PARD6A",
    "PARD6B",
    "PARD6G",
    "PATJ",
    "PRKCI",
)

# Canonical TJ strand / JAM gene in KEGG_TIGHT_JUNCTION but not in the two sets above.
CANONICAL_ADDITIONS = ("JAM2",)

# Catalog symbol -> current HGNC / cBioPortal symbol when the catalog name is retired.
# Rows are renamed back to the catalog symbol after fetch so the universe stays fixed.
CURRENT_SYMBOL = {
    "BVES": "POPDC1",  # HGNC: previous symbol BVES
    "DDX58": "RIGI",  # HGNC: previous symbol DDX58
}

TJ_UNIVERSE = tuple(
    sorted(
        set(GOCC_BICELLULAR_TIGHT_JUNCTION)
        | set(REACTOME_TIGHT_JUNCTION_INTERACTIONS)
        | set(CANONICAL_ADDITIONS)
    )
)

GENESET_PROVENANCE = {
    "gocc_bicellular_tight_junction": {
        "id": "GO:0005923",
        "source": "Enrichr GO_Cellular_Component_2021 (same GMT as claim A8/A9)",
        "n": len(GOCC_BICELLULAR_TIGHT_JUNCTION),
    },
    "reactome_tight_junction_interactions": {
        "id": "R-HSA-420029",
        "source": "Enrichr Reactome_2022 (same GMT as claim A8/A9)",
        "n": len(REACTOME_TIGHT_JUNCTION_INTERACTIONS),
    },
    "canonical_additions": {
        "genes": list(CANONICAL_ADDITIONS),
        "reason": "JAM2 is a canonical JAM in KEGG Tight junction; GOCC/Reactome lists F11R/JAM3/JAML but not JAM2",
    },
    "excluded": {
        "kegg_tight_junction_nonjunctional": (
            "KEGG Tight junction actin/MAPK/tubulin/myosin members were not added. "
            "They are pathway genes, not a TJ-gene catalog."
        ),
        "a9_slide_list_is_not_the_universe": (
            "CLDN1/4/7, F11R, PARD3 are annotated when present but do not define "
            "who is tested or who counts as a recurrent TJ gene."
        ),
    },
    "n_universe": len(TJ_UNIVERSE),
    "target_excluded": TARGET,
}


def membership(gene: str) -> dict:
    return {
        "in_gocc_go0005923": gene in GOCC_BICELLULAR_TIGHT_JUNCTION,
        "in_reactome_r_hsa_420029": gene in REACTOME_TIGHT_JUNCTION_INTERACTIONS,
        "in_canonical_additions": gene in CANONICAL_ADDITIONS,
        "a9_slide_gene": gene in A9_SLIDE_GENES,
    }
