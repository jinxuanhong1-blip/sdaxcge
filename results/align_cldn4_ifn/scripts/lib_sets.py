#!/usr/bin/env python3
"""Gene-set definitions and human<->mouse symbol handling.

The user's primary readout is a six-gene IFN/MHC-I panel (IFI27, OAS2, IFIT1,
MX1, ISG15, HLA-A). We test that panel verbatim plus wider IFN / antigen-
processing-and-presentation (APM) sets, and two IFN-unrelated hallmark sets
that act as negative controls for a global "everything moves" artefact.
"""
import os

DATA = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "data"))

# ---------------------------------------------------------------- user panel
USER_CORE6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# Type-I interferon-stimulated genes (canonical, literature-standard panel).
ISG_CORE = [
    "ISG15", "IFI6", "IFI27", "IFI44", "IFI44L", "IFIT1", "IFIT2", "IFIT3",
    "IFIT5", "IFITM1", "IFITM2", "IFITM3", "MX1", "MX2", "OAS1", "OAS2",
    "OAS3", "OASL", "RSAD2", "USP18", "BST2", "XAF1", "STAT1", "STAT2",
    "IRF7", "IRF9", "DDX58", "IFIH1", "SAMD9", "SAMD9L", "HERC5", "HERC6",
    "EPSTI1", "CMPK2", "PARP9", "DTX3L", "LY6E", "SP100", "SP110", "PLSCR1",
]

# MHC class I antigen processing & presentation machinery.
MHC1_APM = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "NLRC5",
    "ERAP1", "ERAP2", "CALR", "PDIA3", "CANX", "IRF1",
]

# Mouse counterparts (MGI orthologs; MHC-I is not 1:1 so classical H2 loci are
# used in place of HLA-A/B/C and H2-T23/H2-Q in place of HLA-E/F).
MOUSE_OVERRIDE = {
    "IFI27": ["Ifi27", "Ifi27l2a", "Ifi27l2b"],
    "HLA-A": ["H2-K1", "H2-D1"],
    "HLA-B": ["H2-K1", "H2-D1"],
    "HLA-C": ["H2-L", "H2-D1"],
    "HLA-E": ["H2-T23"],
    "HLA-F": ["H2-Q7", "H2-Q6"],
    "DDX58": ["Ddx58", "Rigi"],
    "ERAP1": ["Erap1"],
    "ERAP2": ["Lnpep"],
    "H2-L": ["H2-D1"],
}


def to_mouse(symbols):
    """Map human symbols to mouse symbols (override table, else Title-case)."""
    out = []
    for s in symbols:
        if s in MOUSE_OVERRIDE:
            out.extend(MOUSE_OVERRIDE[s])
        else:
            out.append(s[0] + s[1:].lower() if s.isupper() else s.capitalize())
    seen, uniq = set(), []
    for s in out:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return uniq


def read_gmt(path):
    sets = {}
    with open(path) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) > 2:
                sets[p[0]] = [g for g in p[2:] if g]
    return sets


def gene_sets(species="human"):
    """Return {set_name: [symbols]} for the requested species."""
    hs = read_gmt(os.path.join(DATA, "h.all.v2023.2.Hs.symbols.gmt"))
    re_ = read_gmt(os.path.join(DATA, "c2.cp.reactome.v2023.2.Hs.symbols.gmt"))
    mm = read_gmt(os.path.join(DATA, "mh.all.v2023.2.Mm.symbols.gmt"))

    reactome_pick = {
        "REACTOME_IFN_ALPHA_BETA":
            "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
        "REACTOME_MHC1_PEPTIDE_LOADING":
            "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
    }

    sets = {
        "USER_CORE6": USER_CORE6,
        "ISG_CORE": ISG_CORE,
        "MHC1_APM": MHC1_APM,
        "HALLMARK_IFN_ALPHA": hs["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        "HALLMARK_IFN_GAMMA": hs["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        # negative / specificity controls
        "CTRL_HALLMARK_MYC_V1": hs["HALLMARK_MYC_TARGETS_V1"],
        "CTRL_HALLMARK_OXPHOS": hs["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
    }
    for short, full in reactome_pick.items():
        if full in re_:
            sets[short] = re_[full]

    if species == "mouse":
        out = {k: to_mouse(v) for k, v in sets.items()}
        # prefer MSigDB's own mouse hallmark symbols where available
        for k, full in [("HALLMARK_IFN_ALPHA", "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
                        ("HALLMARK_IFN_GAMMA", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
                        ("CTRL_HALLMARK_MYC_V1", "HALLMARK_MYC_TARGETS_V1"),
                        ("CTRL_HALLMARK_OXPHOS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")]:
            if full in mm:
                out[k] = mm[full]
        return out
    return sets


# Ordered for reporting: IFN/APM sets first, controls last.
SET_ORDER = ["USER_CORE6", "ISG_CORE", "MHC1_APM", "HALLMARK_IFN_ALPHA",
             "HALLMARK_IFN_GAMMA", "REACTOME_IFN_ALPHA_BETA",
             "REACTOME_MHC1_PEPTIDE_LOADING",
             "CTRL_HALLMARK_MYC_V1", "CTRL_HALLMARK_OXPHOS"]


if __name__ == "__main__":
    for sp in ("human", "mouse"):
        gs = gene_sets(sp)
        print(f"--- {sp}")
        for k in SET_ORDER:
            if k in gs:
                print(f"  {k:32s} n={len(gs[k]):4d}  {gs[k][:6]}")
