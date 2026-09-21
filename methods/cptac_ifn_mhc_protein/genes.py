"""Locked protein list for CPTAC LUAD/LSCC: CLDN4 and TACSTD2 vs MHC-I/IFN/CD8A.

Not a search. Predictors and endpoints are fixed before the correlations
are read. IFN ligands and B2M are on the list so absence is reported.
Antigen-processing genes are an extra MHC-I cassette, not a second primary.
"""

from __future__ import annotations

# symbol -> Ensembl gene id. Version suffixes in the freeze are matched by prefix.
ENSEMBL: dict[str, list[str]] = {
    "CLDN4": ["ENSG00000189143"],
    "TACSTD2": ["ENSG00000184292"],
    "CD8A": ["ENSG00000153563"],
    "CD8B": ["ENSG00000172116"],
    "IFNG": ["ENSG00000111537"],
    "IFNA1": ["ENSG00000197919"],
    "IFNB1": ["ENSG00000171855"],
    "IFNAR1": ["ENSG00000142166"],
    "IFNAR2": ["ENSG00000159110"],
    "IFNGR1": ["ENSG00000027697"],
    "IFNGR2": ["ENSG00000159128"],
    "STAT1": ["ENSG00000115415"],
    "STAT2": ["ENSG00000170581"],
    "IRF1": ["ENSG00000125347"],
    "IRF9": ["ENSG00000213928"],
    "JAK1": ["ENSG00000162434"],
    "JAK2": ["ENSG00000096968"],
    "ISG15": ["ENSG00000187608"],
    "MX1": ["ENSG00000157601"],
    "OAS1": ["ENSG00000089127"],
    "IFI35": ["ENSG00000068079"],
    "HLA-A": ["ENSG00000206503"],
    "HLA-B": ["ENSG00000234745"],
    "HLA-C": ["ENSG00000204525"],
    "B2M": ["ENSG00000166710"],
    "TAP1": ["ENSG00000168394"],
    "TAP2": ["ENSG00000204267"],
    "TAPBP": ["ENSG00000231925"],
    "PSMB8": ["ENSG00000204264"],
    "PSMB9": ["ENSG00000240065"],
    "NLRC5": ["ENSG00000140853"],
}

PREDICTORS = ["CLDN4", "TACSTD2"]

IFN_LIGANDS = ["IFNG", "IFNA1", "IFNB1"]
IFN_RECEPTORS = ["IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2"]
IFN_SIGNALING = ["STAT1", "STAT2", "IRF1", "IRF9", "JAK1", "JAK2"]
IFN_ISG = ["ISG15", "MX1", "OAS1", "IFI35"]
MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
# Extra antigen-processing proteins. Not folded into the primary MHC-I score.
APM = ["TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"]
CD8 = ["CD8A", "CD8B"]

IFN_SCORE_MEMBERS = IFN_SIGNALING + IFN_ISG
MHC1_SCORE_MEMBERS = ["HLA-A", "HLA-B", "HLA-C"]  # B2M is absent in this freeze

CLASS = {
    **{g: "predictor" for g in PREDICTORS},
    **{g: "cd8" for g in CD8},
    **{g: "ifn_ligand" for g in IFN_LIGANDS},
    **{g: "ifn_receptor" for g in IFN_RECEPTORS},
    **{g: "ifn_signaling" for g in IFN_SIGNALING},
    **{g: "ifn_isg" for g in IFN_ISG},
    **{g: "mhc1" for g in MHC1},
    **{g: "apm_extra" for g in APM},
}

# Gene-level endpoints scored against each predictor.
ENDPOINTS = CD8 + IFN_LIGANDS + IFN_RECEPTORS + IFN_SIGNALING + IFN_ISG + MHC1 + APM

# Primary support family: these scores, both predictors, both cohorts.
PRIMARY_SCORES = ["MHC1_protein", "IFN_core_protein", "CD8A_protein"]

# Decomposition of the locked list. Not new genes. Used only to see whether
# the IFN-core null is ISG dilution. Minimum quantified members for a score.
SENSITIVITY_PANELS: dict[str, tuple[list[str], int]] = {
    "IFN_signaling_protein": (IFN_SIGNALING, 4),
    "IFN_isg_protein": (IFN_ISG, 3),
    "IFN_receptor_protein": (["IFNAR1", "IFNGR1"], 2),
    "APM_protein": (["TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5"], 4),
}
