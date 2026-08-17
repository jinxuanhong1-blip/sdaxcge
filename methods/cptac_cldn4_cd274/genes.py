"""Locked protein list for CPTAC CLDN4 vs CD274 / IFN.

Not a search. CD274 is the primary endpoint. IFN ligands, receptors,
signaling, and a short ISG cassette are scored only if the TMT row exists
and has ≥8 non-missing tumors. MHC-I (HLA-A/B/C, B2M) is an extra IFN-axis
cassette, same four genes as the DepMap CLDN4–IFN page (PR #304).
"""

from __future__ import annotations

# symbol -> preferred Ensembl gene id(s). Alternates seen in GENCODE / freeze.
ENSEMBL: dict[str, list[str]] = {
    "CLDN4": ["ENSG00000189143"],
    "CD274": ["ENSG00000120217"],
    # IFN ligands (often missing from TMT)
    "IFNG": ["ENSG00000111537"],
    "IFNA1": ["ENSG00000197919"],
    "IFNB1": ["ENSG00000171855"],
    # IFN receptors
    "IFNAR1": ["ENSG00000142166"],
    "IFNAR2": ["ENSG00000159110"],
    "IFNGR1": ["ENSG00000027697"],
    "IFNGR2": ["ENSG00000159128"],
    # IFN signaling
    "STAT1": ["ENSG00000115415"],
    "STAT2": ["ENSG00000170581"],
    "IRF1": ["ENSG00000125347"],
    "IRF9": ["ENSG00000213928"],
    "JAK1": ["ENSG00000162434"],
    "JAK2": ["ENSG00000096968"],
    # Canonical ISGs
    "ISG15": ["ENSG00000187608"],
    "MX1": ["ENSG00000157601"],
    "OAS1": ["ENSG00000089127"],
    "IFI35": ["ENSG00000068079"],
    # MHC-I extra (IFN axis; not IFN ligands). IDs confirmed unique in freeze v1.2.
    "HLA-A": ["ENSG00000206503"],
    "HLA-B": ["ENSG00000234745"],
    "HLA-C": ["ENSG00000204525"],
    "B2M": ["ENSG00000166710"],
}

PRIMARY = ["CD274"]

IFN_LIGANDS = ["IFNG", "IFNA1", "IFNB1"]
IFN_RECEPTORS = ["IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2"]
IFN_SIGNALING = ["STAT1", "STAT2", "IRF1", "IRF9", "JAK1", "JAK2"]
IFN_ISG = ["ISG15", "MX1", "OAS1", "IFI35"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]

# Mean-z score members if ≥4 present with ≥8 observations.
IFN_SCORE_MEMBERS = IFN_SIGNALING + IFN_ISG
MHC1_SCORE_MEMBERS = MHC1

CLASS = {
    **{g: "primary" for g in PRIMARY},
    **{g: "ifn_ligand" for g in IFN_LIGANDS},
    **{g: "ifn_receptor" for g in IFN_RECEPTORS},
    **{g: "ifn_signaling" for g in IFN_SIGNALING},
    **{g: "ifn_isg" for g in IFN_ISG},
    **{g: "mhc1_extra" for g in MHC1},
}

ALL_ENDPOINTS = PRIMARY + IFN_LIGANDS + IFN_RECEPTORS + IFN_SIGNALING + IFN_ISG + MHC1
