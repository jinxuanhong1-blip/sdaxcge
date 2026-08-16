"""Curated marker gene sets for the CLDN4-loss analysis.

Symbols given in HUMAN (UPPERCASE). For mouse we title-case on the fly
(mouse orthologs use Title case, e.g. CLDN4 -> Cldn4). A few mouse-only
aliases are added explicitly where the ortholog symbol differs.
"""

# The primary hypothesis gene.
TARGET = ["TACSTD2"]  # Trop-2

# Tight junction, adherens junction, desmosome, and related junction genes.
JUNCTION = [
    # claudins (CLDN4 itself is the perturbed gene; keep as internal control)
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7",
    "CLDN8", "CLDN9", "CLDN10", "CLDN11", "CLDN12", "CLDN14", "CLDN15",
    "CLDN18", "CLDN23",
    # tight-junction scaffold / integral
    "TJP1", "TJP2", "TJP3", "OCLN", "MARVELD2", "MARVELD3",
    "F11R", "JAM2", "JAM3", "CGN", "CGNL1", "CLDND1",
    # adherens junction
    "CDH1", "CDH2", "CTNNB1", "CTNNA1", "CTNND1", "JUP", "VCL",
    # desmosome
    "DSP", "DSG2", "DSG3", "DSC2", "DSC3", "PKP2", "PKP3", "PERP",
    # epithelial polarity / junction-associated
    "CRB3", "PARD3", "EPCAM",
]

# Interferon / antiviral / antigen-presentation immune genes (classic ISG core
# plus type I/II IFN machinery and inflammatory chemokines).
IFN_IMMUNE = [
    # type I ISG core
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2",
    "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3",
    "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L", "IFI16",
    "DDX58", "IFIH1", "DDX60", "XAF1", "HERC5", "USP18", "CMPK2",
    "BST2", "GBP1", "GBP2", "GBP4", "GBP5", "EIF2AK2", "ZBP1",
    # signal transducers / regulators
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "JAK2", "SOCS1", "SOCS3",
    # antigen presentation / MHC
    "B2M", "TAP1", "TAP2", "PSMB8", "PSMB9", "NLRC5", "HLA-A", "HLA-B",
    "HLA-C", "HLA-E", "HLA-F", "TAPBP",
    # inflammatory chemokines/cytokines (IFN-inducible)
    "CXCL10", "CXCL11", "CXCL9", "CCL5", "CCL2", "IL6", "TNF",
    "NFKB1", "RELA", "IL15", "IL32",
    # ligands
    "IFNB1", "IFNG", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNL1",
]

# Mouse-specific symbol fixes where the ortholog differs from title-cased human.
MOUSE_ALIASES = {
    "HLA-A": ["H2-K1", "H2-D1"],
    "HLA-B": ["H2-K1"],
    "HLA-C": ["H2-D1"],
    "HLA-E": ["H2-T23"],
    "HLA-F": [],
    "B2M": ["B2m"],
    "IFI44L": [],   # no clear mouse ortholog
    "IFI27": ["Ifi27", "Ifi27l2a", "Ifi27l2b"],
    "MX1": ["Mx1"],
    "MX2": ["Mx2"],
    "IFIT1": ["Ifit1", "Ifit1bl1"],
    "OASL": ["Oasl1", "Oasl2"],
    "OAS1": ["Oas1a", "Oas1b", "Oas1g"],
    "GBP1": ["Gbp2", "Gbp3"],
    "IL32": [],  # no mouse ortholog
    "IFNL1": [],
}


def human_symbols():
    return {"TACSTD2": TARGET, "JUNCTION": JUNCTION, "IFN_IMMUNE": IFN_IMMUNE}


def _to_mouse(sym):
    """Return list of candidate mouse symbols for a human symbol."""
    if sym in MOUSE_ALIASES:
        aliases = MOUSE_ALIASES[sym]
        base = sym[0].upper() + sym[1:].lower() if "-" not in sym else None
        out = list(aliases)
        if base and base not in out:
            out.append(base)
        return out
    # default: Title-case (first letter upper, rest lower)
    return [sym[0].upper() + sym[1:].lower()]


def mouse_symbols():
    out = {}
    for name, lst in human_symbols().items():
        s = set()
        for h in lst:
            for m in _to_mouse(h):
                if m:
                    s.add(m)
        out[name] = sorted(s)
    return out


if __name__ == "__main__":
    import json
    print("HUMAN:", json.dumps({k: len(v) for k, v in human_symbols().items()}))
    print("MOUSE junction sample:", mouse_symbols()["JUNCTION"][:10])
    print("MOUSE ifn sample:", mouse_symbols()["IFN_IMMUNE"][:10])
