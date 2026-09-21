"""Mouse gene modules for the public KP/KL wave beyond GSE137244.

NHEJ_core is the canonical non-homologous end-joining machinery
(Reactome R-HSA-5693571 core members; mouse symbols). MRN-complex
genes are not included.

STING_core is the cGAS–STING pathway machinery, not the downstream
interferon program. STING_ISG is reported separately so Sting1
silencing is not averaged away inside an ISG score.

Aliases cover older GEO annotations (Tmem173, Mb21d1).
"""

from __future__ import annotations

# Canonical symbol -> aliases accepted in a matrix (first name is canonical).
NHEJ_CORE = {
    "Xrcc5": ["Xrcc5", "Ku80"],
    "Xrcc6": ["Xrcc6", "Ku70"],
    "Prkdc": ["Prkdc", "DNA-PKcs", "Prkdc"],
    "Dclre1c": ["Dclre1c", "Artemis"],
    "Lig4": ["Lig4"],
    "Xrcc4": ["Xrcc4"],
    "Nhej1": ["Nhej1", "Xlf"],
    "Poll": ["Poll"],
    "Polm": ["Polm"],
    "Pnkp": ["Pnkp"],
    "Aptx": ["Aptx"],
    "Paxx": ["Paxx"],
    "Aplf": ["Aplf"],
}

STING_CORE = {
    "Cgas": ["Cgas", "Mb21d1"],
    "Sting1": ["Sting1", "Tmem173", "Sting"],
    "Tbk1": ["Tbk1"],
    "Irf3": ["Irf3"],
    "Ikbke": ["Ikbke", "Ikke"],
}

STING_ISG = {
    "Ifnb1": ["Ifnb1", "Ifnb"],
    "Cxcl10": ["Cxcl10"],
    "Isg15": ["Isg15"],
    "Ifit1": ["Ifit1"],
    "Rsad2": ["Rsad2"],
    "Mx1": ["Mx1"],
    "Stat1": ["Stat1"],
    "Irf7": ["Irf7"],
}

CLDN4 = {"Cldn4": ["Cldn4"]}

QC = {
    "Stk11": ["Stk11", "Lkb1"],
    "Actb": ["Actb"],
    "Gapdh": ["Gapdh"],
}

MODULES = {
    "NHEJ": NHEJ_CORE,
    "STING_core": STING_CORE,
    "STING_ISG": STING_ISG,
    "Cldn4": CLDN4,
    "QC": QC,
}


def all_aliases() -> set[str]:
    out: set[str] = set()
    for mod in MODULES.values():
        for aliases in mod.values():
            out.update(aliases)
    return out
