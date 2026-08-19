"""Public gene lists for the TISMO / GEMM cell-line analysis.

Tight-junction *core* is the user-specified 9-gene program.
GO / Reactome sets are Enrichr libraries (human symbols); mapped to mouse
by case-insensitive match against each matrix.
"""

from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# User-specified core (mouse symbols). Tacstd2 is the ranking gene, not in the score.
TJ_CORE = [
    "Cldn3",
    "Cldn4",
    "Cldn6",
    "Cldn7",
    "Cdh1",
    "F11r",
    "Ocln",
    "Tjp1",
    "Nectin2",
]

FOCUS_GENES = ["Tacstd2", "Cldn4", *TJ_CORE]


def load_human_sets(path: Path | None = None) -> dict[str, dict]:
    p = path or (HERE / "go_sets_human.json")
    return json.loads(p.read_text())


def to_mouse(human_genes: list[str], universe: list[str]) -> list[str]:
    """Map human symbols onto a mouse (or mixed-case) expression universe."""
    lookup = {g.lower(): g for g in universe}
    out = []
    seen = set()
    for h in human_genes:
        key = str(h).lower()
        if key in lookup and lookup[key] not in seen:
            out.append(lookup[key])
            seen.add(lookup[key])
    return out


def gene_sets_for_universe(universe: list[str]) -> dict[str, list[str]]:
    human = load_human_sets()
    sets = {
        "TJ_CORE": [g for g in TJ_CORE if g in set(universe) or g.lower() in {u.lower() for u in universe}],
    }
    # remap core through the same lookup so Cldn4 vs CLDN4 works
    sets["TJ_CORE"] = to_mouse(TJ_CORE, universe)
    for name, rec in human.items():
        sets[name] = to_mouse(rec["genes"], universe)
    return sets
