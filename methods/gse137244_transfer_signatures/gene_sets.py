"""Frozen mouse gene lists for scoring GSE137244 signatures on other data.

These lists are the transferable objects. They are not re-fit on a query
single-cell matrix. Score a query by the mean of log-normalized expression
of the genes that are present; leave the absent genes out of that mean.
"""

from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "gene_sets"

# Epithelial tight-junction core. Same 18 structural genes as the GSE137244
# deep note. Not chosen to hit a target delta.
TJ_EPITHELIAL = [
    "Cldn1",
    "Cldn3",
    "Cldn4",
    "Cldn7",
    "Ocln",
    "Marveld2",
    "Marveld3",
    "Tjp1",
    "Tjp2",
    "Tjp3",
    "F11r",
    "Jam2",
    "Jam3",
    "Cgn",
    "Cgnl1",
    "Crb3",
    "Ildr1",
    "Lsr",
]

# Shorter TJ list used for the TISMO TJ score. Tacstd2 is not inside it.
TJ_TISMO = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]

# Published epithelial leading edge of the TJ core on this series
# (Welch prerank, prior GSE137244 deep note). Vendored so a later scRNA
# score does not re-cut the edge.
CLDN4_TJ_EDGE = [
    "Cldn4",
    "Cldn1",
    "Cldn7",
    "Cgnl1",
    "Marveld2",
    "Marveld3",
    "Tjp1",
    "Tjp2",
    "Ildr1",
    "Cldn3",
    "Ocln",
]

# Enzymatic c-NHEJ. PAXX is included when the matrix has that symbol.
# 53BP1 (Trp53bp1) is the pathway factor used with this core previously.
NHEJ_CORE = [
    "Xrcc6",
    "Xrcc5",
    "Prkdc",
    "Xrcc4",
    "Lig4",
    "Nhej1",
    "Dclre1c",
    "Paxx",
    "Poll",
    "Polm",
    "Trp53bp1",
]

# DNA-end sensing and 53BP1-pathway genes. Separate from the enzymatic core
# so the sweep can see whether a wider NHEJ list changes the sign.
NHEJ_EXTENDED = [
    "Trp53bp1",
    "Rnf8",
    "Rnf168",
    "Rif1",
    "Mad2l2",
    "Atm",
    "H2afx",
    "Mre11a",
    "Rad50",
    "Nbn",
    "Parp1",
    "Pnkp",
    "Aplf",
    "Aptx",
]

# cGAS–STING machinery. This matrix uses Mb21d1 and Tmem173.
# ISGs are not in this list.
STING_CORE = ["Mb21d1", "Tmem173", "Tbk1", "Ikbke", "Irf3"]

# Compact epithelial ISG list used on prior public KL vs KP contrasts.
IFN_COMPACT = [
    "Stat1",
    "Stat2",
    "Irf1",
    "Irf7",
    "Irf9",
    "Isg15",
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Mx1",
    "Oasl2",
    "Rsad2",
    "Ifih1",
    "Ddx58",
    "Ifnb1",
]

# OAS antiviral enzymes. Included in the sweep because this family can move
# opposite the compact ISG list. It is not the primary IFN signature.
IFN_OAS = ["Oas1a", "Oas1b", "Oas1g", "Oas2", "Oas3", "Oasl1", "Oasl2"]

# MHC-I antigen-presentation machinery for a tumor-cell score.
# Immunoproteasome subunits are included because that is the Deng limb.
# CD4/CD8/KLRC genes from the broad KEGG list are not in this set.
APM_MHCI = [
    "B2m",
    "H2-K1",
    "H2-D1",
    "H2-L",
    "Tap1",
    "Tap2",
    "Tapbp",
    "Psmb8",
    "Psmb9",
    "Psmb10",
    "Erap1",
    "Nlrc5",
    "Calr",
    "Canx",
    "Pdia3",
    "Psme1",
    "Psme2",
]

# Explicit ortholog renames used only when the requested symbol is absent
# and the mouse symbol is present. No other ortholog guessing.
ALIASES = {
    "wars1": "Wars",
    "rigi": "Ddx58",
    "tmt1b": "Mettl7b",
    "cgas": "Mb21d1",
    "sting1": "Tmem173",
}


def read_symbol_file(name: str) -> list[str]:
    lines = (SOURCE / name).read_text().splitlines()
    # First line is a source note, not a gene.
    return [ln.strip() for ln in lines[1:] if ln.strip()]


def match_symbols(requested: list[str], universe: list[str]) -> tuple[list[str], list[str]]:
    lookup = {g.lower(): g for g in universe}
    used: list[str] = []
    missing: list[str] = []
    seen: set[str] = set()
    for raw in requested:
        key = raw.lower()
        hit = lookup.get(key)
        if hit is None and key in ALIASES:
            hit = lookup.get(ALIASES[key].lower())
        if hit is None or hit in seen:
            if hit is None:
                missing.append(raw)
            continue
        used.append(hit)
        seen.add(hit)
    return used, missing
