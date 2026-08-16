"""Mouse gene sets for the C4 analog (GSE334497 4T1 Trop2 KO).

Claim C4 (private): CLDN4 KD opens IFN / MHC-I / APM.
This analog tests Trop2 KO (upstream of the claudin barrier) in 4T1 tumors.

Symbols are mouse. Sets are intentionally small and pre-specified.
"""

# Perturbation / barrier genes (not a test set; QC + analog check)
QC_GENES = [
    "Tacstd2",  # CRISPR target; must drop
    "Cldn4",    # C4 analog hinge: does Trop2 KO lower Cldn4?
    "Cldn7",    # paper's stated functional mediator (not CLDN4)
    "Cldn1",
    "Ocln",
    "Tjp1",
    "F11r",
    "Epcam",
]

# User-requested focal genes
FOCAL = ["Cldn4", "Cxcl9"]

# Private-claim CORE6 (human IFI27 OAS2 IFIT1 MX1 ISG15 HLA-A) → mouse
USER_CORE6 = [
    "Ifi27",
    "Ifi27l2a",
    "Ifi27l2b",
    "Oas2",
    "Ifit1",
    "Mx1",
    "Isg15",
    "H2-K1",
    "H2-D1",
]

# Compact type-I ISG panel (mouse)
ISG_CORE = [
    "Isg15", "Ifi27", "Ifi27l2a", "Ifi27l2b", "Ifi44",
    "Ifit1", "Ifit2", "Ifit3",
    "Ifitm1", "Ifitm2", "Ifitm3",
    "Mx1", "Mx2",
    "Oas1a", "Oas1b", "Oas1g", "Oas2", "Oas3", "Oasl1", "Oasl2",
    "Rsad2", "Usp18", "Bst2", "Xaf1",
    "Stat1", "Stat2", "Irf7", "Irf9",
    "Ddx58", "Ifih1", "Samd9l", "Herc6",
    "Epsti1", "Cmpk2", "Parp9", "Dtx3l", "Ly6e",
    "Sp100", "Sp110", "Plscr1",
]

# MHC-I antigen-processing / presentation (mouse)
MHC1_APM = [
    "H2-K1", "H2-D1", "H2-Q4", "H2-Q6", "H2-Q7", "H2-T23",
    "B2m", "Tap1", "Tap2", "Tapbp",
    "Psmb8", "Psmb9", "Psmb10", "Psme1", "Psme2",
    "Nlrc5", "Erap1", "Calr", "Pdia3", "Canx", "Irf1",
]

# IFN-inducible CXC chemokines (includes user Cxcl9)
CXCL_IFN = ["Cxcl9", "Cxcl10", "Cxcl11"]

# T-cell / cytotoxicity — infiltration confound for bulk tumor RNA
T_CYT = [
    "Cd3e", "Cd3d", "Cd8a", "Cd4", "Cd2",
    "Gzmb", "Gzma", "Prf1", "Ifng", "Nkg7", "Klrd1",
]

# Negative-control sets (should not systematically follow the IFN claim)
CTRL_OXPHOS = [
    "Cox4i1", "Cox5a", "Cox6a1", "Cox7a2", "Ndufa1", "Ndufa4",
    "Ndufb2", "Ndufs2", "Sdha", "Sdhb", "Uqcrc1", "Uqcrc2",
    "Atp5a1", "Atp5b", "Atp5c1", "Cycs", "Vdac1", "Idh2",
    "Mdh2", "Fh1",
]

CTRL_MYC = [
    "Myc", "Ncl", "Npm1", "Nolc1", "Nop56", "Fbl",
    "Rpl5", "Rpl8", "Rps2", "Rps6", "Eif4a1", "Eif4e",
    "Cad", "Odc1", "Ldha", "Hk2", "Tyms", "Rrm1",
    "Cct2", "Hsp90ab1",
]

SETS = {
    "USER_CORE6": USER_CORE6,
    "ISG_CORE": ISG_CORE,
    "MHC1_APM": MHC1_APM,
    "CXCL_IFN": CXCL_IFN,
    "T_CYT": T_CYT,
    "CTRL_OXPHOS": CTRL_OXPHOS,
    "CTRL_MYC": CTRL_MYC,
}

# Claim-facing order for reporting
REPORT_ORDER = [
    "USER_CORE6",
    "ISG_CORE",
    "MHC1_APM",
    "CXCL_IFN",
    "T_CYT",
    "CTRL_OXPHOS",
    "CTRL_MYC",
]
