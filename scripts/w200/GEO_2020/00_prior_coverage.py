#!/usr/bin/env python3
"""Inventory GSE accessions already touched by prior waves of this repo.

Reads the pre-extracted /tmp/inv/branch_gse.tsv (git-grep of every remote
branch) plus path-level mentions. Writes a coverage table so "leftover"
is defined against evidence, not memory.
"""
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / "results" / "w200" / "GEO_2020"
OUT.mkdir(parents=True, exist_ok=True)

BRANCH_TSV = Path("/tmp/inv/branch_gse.tsv")
PATH_TSV = Path("/tmp/inv/path_gse.tsv")

# Waves that already claimed GEO lung ICI / TACSTD2 / CLDN4 work.
PRIOR_WAVES = {
    "cursor/fable-geo-2019-2021-c71e": "geo_2019_2021",
    "cursor/fable-geo-2022-2023-99b1": "geo_2022_2023",
    "cursor/fable-geo-2026-lung-ici-1c1b": "geo_2026",
    "cursor/geo-remaining-sweep-b8cf": "geo_remaining_sweep",
    "cursor/fable-ici-bulk-tacstd2-cldn4-108c": "ici_bulk",
    "cursor/ici-lung-response-1ce3": "ici_lung_response",
    "cursor/fable-ici-bulk-tacstd2-cldn4-108c": "ici_bulk",
}


def load_pairs(path):
    pairs = []
    if not path.exists():
        return pairs
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        branch, acc = line.split("\t", 1)
        branch = branch.replace("origin/", "")
        pairs.append((branch, acc.strip()))
    return pairs


def main():
    mentioned = defaultdict(set)
    for branch, acc in load_pairs(BRANCH_TSV):
        mentioned[acc].add(branch)
    in_path = defaultdict(set)
    for branch, acc in load_pairs(PATH_TSV):
        in_path[acc].add(branch)

    rows = []
    for acc in sorted(set(mentioned) | set(in_path)):
        waves = sorted(mentioned.get(acc, set()) | in_path.get(acc, set()))
        prior = [PRIOR_WAVES[w] for w in waves if w in PRIOR_WAVES]
        rows.append({
            "accession": acc,
            "n_branches": len(waves),
            "in_dedicated_path": "Y" if acc in in_path else "N",
            "prior_geo_waves": ";".join(sorted(set(prior))),
            "branches": ";".join(waves),
        })

    out = OUT / "prior_gse_coverage.tsv"
    cols = ["accession", "n_branches", "in_dedicated_path", "prior_geo_waves", "branches"]
    with out.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print(f"Wrote {len(rows)} accessions to {out}")
    print(f"  with dedicated files: {sum(1 for r in rows if r['in_dedicated_path']=='Y')}")
    print(f"  mentioned in a prior GEO wave: {sum(1 for r in rows if r['prior_geo_waves'])}")


if __name__ == "__main__":
    main()
