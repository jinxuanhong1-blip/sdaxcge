"""Pre-specified mouse gene sets for GSE179502 (Cldn4-only).

A8 Hallmark / MHC / TJ families are human. They are mapped to mouse
symbols that exist on the GSE179502 mm10 features table. Cldn4 is
held out of TJ. Tacstd2 is audit-only and is never a gate.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Paper Fig. 4c ATII markers + standard AT2 identity (Murray/Winslow 2022).
AT2 = [
    "Sftpc",
    "Sftpb",
    "Sftpa1",
    "Lamp3",
    "Slc34a2",
    "Napsa",
    "Lyz2",
    "Abca3",
    "Cxcl15",
    "Hc",
    "Lpcat1",
    "Etv5",
    "Cebpa",
    "Nkx2-1",
]

# Compact leftover IFN / MHC / TJ cores (same names as other mouse slices).
IFN_CORE = [
    "Ifit1",
    "Ifit2",
    "Ifit3",
    "Isg15",
    "Mx1",
    "Stat1",
    "Irf7",
    "Cxcl10",
    "Bst2",
    "Rsad2",
]
MHC_CORE = [
    "H2-K1",
    "H2-D1",
    "H2-Q4",
    "B2m",
    "Tap1",
    "Tap2",
    "Psmb8",
    "Psmb9",
    "Nlrc5",
]
TJ_CORE = [
    "Cldn3",
    "Cldn7",
    "Cldn18",
    "Ocln",
    "Tjp1",
    "Tjp2",
    "F11r",
    "Marveld2",
    "Marveld3",
    "Cdh1",
    "Cgn",
    "Pard3",
]

HLA_TO_H2 = {
    "HLA-A": ["H2-K1"],
    "HLA-B": ["H2-D1"],
    "HLA-C": ["H2-Q4"],
    "HLA-E": ["H2-T23"],
    "HLA-F": ["H2-Q10"],
    "HLA-G": ["H2-Q7"],
    "HLA-DMA": ["H2-DMa"],
    "HLA-DQA1": ["H2-Aa"],
    "HLA-DRB1": ["H2-Eb1"],
}

# Human leftover aliases that do not title-case onto mm10.
ALIASES = {
    "OAS1": ["Oas1a", "Oas1g"],
    "OASL": ["Oasl1", "Oasl2"],
    "OAS2": ["Oas2"],
    "OAS3": ["Oas3"],
    "MX1": ["Mx1"],
    "MX2": ["Mx2"],
    "IFI27": ["Ifi27", "Ifi27l2a"],
    "IFI44L": ["Ifi44"],
    "WARS1": ["Wars"],
    "C1R": ["C1ra"],
    "C1S": ["C1s1"],
    "FCGR1A": ["Fcgr1"],
    "SAMD9": ["Samd9l"],
    "SAMD9L": ["Samd9l"],
    "MARCHF1": ["Marchf1", "March1"],
    "TENT5A": ["Tent5a", "Fam46a"],
    "PALS1": ["Mpp5"],
    "PATJ": ["Inadl"],
    "PNP": ["Pnp"],
    "MT2A": ["Mt2"],
}


def human_to_mouse_candidates(symbol: str) -> list[str]:
    s = symbol.strip()
    if not s:
        return []
    if s in HLA_TO_H2:
        return list(HLA_TO_H2[s])
    if s in ALIASES:
        return list(ALIASES[s])
    if s.startswith("MIR"):
        return []
    if "-" in s and s.upper() == s:
        head, tail = s.split("-", 1)
        return [head.capitalize() + "-" + tail]
    if s.upper() == s:
        return [s.capitalize()]
    return [s[0].upper() + s[1:]]


def map_human_set(human_genes: list[str], present: set[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for g in human_genes:
        for m in human_to_mouse_candidates(g):
            if m in present and m not in seen:
                seen.add(m)
                out.append(m)
    return out


def load_a8_mouse(present: set[str]) -> dict[str, list[str]]:
    a8 = json.loads((HERE / "data" / "a8_families.json").read_text())
    sets = a8["sets"]
    ifn_h = sorted(
        set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
    )
    mhc_h = list(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"])
    tj_h = sorted(
        set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    )
    for g in a8.get("focal_genes", []):
        if g not in sets.get("KRT_EPITHELIAL", []):
            tj_h.append(g)
    ifn = map_human_set(ifn_h, present)
    mhc = map_human_set(mhc_h, present)
    tj = map_human_set(tj_h, present)
    for drop in ("Cldn4", "Stk11"):
        tj = [g for g in tj if g != drop]
    return {"IFN_A8": ifn, "MHC_A8": mhc, "TJ_A8": tj}


def present_subset(genes: list[str], present: set[str]) -> list[str]:
    return [g for g in genes if g in present]
