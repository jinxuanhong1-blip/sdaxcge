#!/usr/bin/env python3
"""Filter/rank search results for open, processed matrix files < 2GB.

Reads results/fable_zenodo/search_results.json, scores each record for relevance
to lung ICI / TACSTD2 / CLDN4 and for presence of downloadable processed matrices,
and writes a ranked shortlist to results/fable_zenodo/shortlist.json plus a
markdown table to notes/fable_zenodo/shortlist.md.
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "fable_zenodo"
NOTES = ROOT / "notes" / "fable_zenodo"

MAX_BYTES = 2 * 1024 ** 3  # 2 GB

REL_TERMS = {
    "lung": 3, "nsclc": 3, "adenocarcinoma": 2, "pulmonary": 2,
    "immune checkpoint": 3, "immunotherapy": 3, "ici": 2, "anti-pd": 3,
    "pd-1": 2, "pd-l1": 2, "pembrolizumab": 2, "nivolumab": 2, "atezolizumab": 2,
    "tacstd2": 6, "trop2": 6, "cldn4": 6, "claudin": 3,
    "single cell": 2, "single-cell": 2, "scrna": 2, "sc-rna": 2,
}

MATRIX_EXT = (
    ".h5ad", ".rds", ".h5", ".loom", ".mtx", ".mtx.gz", ".csv", ".csv.gz",
    ".tsv", ".tsv.gz", ".txt.gz", ".xlsx", ".h5seurat", ".zip", ".tar.gz",
    ".feather", ".parquet",
)

OPEN_ACCESS = {"open", "public", "embargoed_open", "cc-by", "cc-by-4.0", "cc0-1.0"}


def is_open(rec):
    ar = (rec.get("access_right") or "").lower()
    if ar in ("restricted", "closed", "private"):
        return False
    if ar in ("open", "public"):
        return True
    # figshare/osf usually public if listed; default treat unknown-non-restricted as maybe
    return ar not in ("restricted", "closed", "private")


def matrix_files(rec):
    out = []
    for f in rec.get("files") or []:
        key = (f.get("key") or "").lower()
        size = f.get("size")
        if not key:
            continue
        if key.endswith(MATRIX_EXT):
            out.append(f)
    return out


def rel_score(rec):
    # Score on the record's own title only; the query text is a fuzzy-match
    # attribution and produces false positives if included.
    text = (rec.get("title") or "").lower()
    s = 0
    hits = []
    for term, w in REL_TERMS.items():
        if term in text:
            s += w
            hits.append(term)
    return s, hits


def main():
    records = json.loads((RESULTS / "search_results.json").read_text())
    scored = []
    for rec in records:
        s, hits = rel_score(rec)
        if s < 3:
            continue
        mfiles = matrix_files(rec)
        downloadable = [
            f for f in mfiles
            if f.get("size") is not None and f["size"] <= MAX_BYTES
        ]
        rec2 = dict(rec)
        rec2["rel_score"] = s
        rec2["rel_hits"] = hits
        rec2["open"] = is_open(rec)
        rec2["matrix_files"] = mfiles
        rec2["downloadable_matrix_files"] = downloadable
        rec2["has_downloadable_matrix"] = bool(downloadable)
        scored.append(rec2)

    # Rank: has core gene term > has downloadable matrix > open > rel_score
    def core(rec):
        return any(t in rec["rel_hits"] for t in ("tacstd2", "trop2", "cldn4"))

    scored.sort(key=lambda r: (
        core(r), r["has_downloadable_matrix"], r["open"], r["rel_score"]
    ), reverse=True)

    (RESULTS / "shortlist.json").write_text(json.dumps(scored, indent=2, ensure_ascii=False))

    lines = ["# Shortlist of candidate records", "",
             "| # | source | rel | core | open | matrix? | title | doi |",
             "|---|--------|-----|------|------|---------|-------|-----|"]
    for i, r in enumerate(scored[:60], 1):
        title = (r.get("title") or "")[:70].replace("|", "/")
        doi = r.get("doi") or ""
        lines.append(
            f"| {i} | {r['source']} | {r['rel_score']} | "
            f"{'Y' if core(r) else ''} | {'Y' if r['open'] else 'N'} | "
            f"{'Y' if r['has_downloadable_matrix'] else ''} | {title} | {doi} |"
        )
    (NOTES / "shortlist.md").write_text("\n".join(lines) + "\n")

    print(f"Scored (rel>=3): {len(scored)}")
    print(f"With downloadable matrix: {sum(1 for r in scored if r['has_downloadable_matrix'])}")
    print(f"Core gene term (TACSTD2/TROP2/CLDN4): {sum(1 for r in scored if core(r))}")
    print("Top 15:")
    for i, r in enumerate(scored[:15], 1):
        print(f" {i:2d}. [{r['source']}] rel={r['rel_score']} core={core(r)} open={r['open']} "
              f"matrix={r['has_downloadable_matrix']} :: {(r.get('title') or '')[:75]}")


if __name__ == "__main__":
    main()
