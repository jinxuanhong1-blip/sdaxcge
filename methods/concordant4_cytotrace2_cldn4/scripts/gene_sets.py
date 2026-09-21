"""Locked scores. CLDN4 is the readout and is not inside barrier or stemness.

Stemness sets are MSigDB c2.cgp v2023.2.Hs. EPCAM and the barrier/keratin
genes are removed so stemness is not the epithelial gate and not the barrier score.
"""
from __future__ import annotations

from pathlib import Path

FOCAL = ("CLDN4",)
BARRIER_KERATIN = ("KRT8", "KRT18", "KRT19", "KRT7", "CDKN1A", "PLAUR")
# Dropped from every stemness set before scoring.
STEMNESS_DROP = set(BARRIER_KERATIN) | {"CLDN4", "TACSTD2", "EPCAM"}

_DIR = Path(__file__).resolve().parents[1] / "data" / "genesets"


def _load(name: str) -> tuple[str, ...]:
    path = _DIR / f"{name}.txt"
    genes = []
    for line in path.read_text().splitlines():
        g = line.strip().upper()
        if g and g not in STEMNESS_DROP and g not in genes:
            genes.append(g)
    return tuple(genes)


def stemness_sets() -> dict[str, tuple[str, ...]]:
    return {
        "benporath_es1": _load("BENPORATH_ES_1"),
        "wong_esc": _load("WONG_EMBRYONIC_STEM_CELL_CORE"),
        "benporath_core9": _load("BENPORATH_ES_CORE_NINE"),
    }
