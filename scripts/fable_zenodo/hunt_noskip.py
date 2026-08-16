#!/usr/bin/env python3
"""Hunt remaining open processed matrices that could contain TACSTD2/CLDN4.

No size cap. Restricted records are still skipped. Writes a candidate table
to results/noskip/zenodo/hunt_candidates.json.
"""
import json
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "noskip" / "zenodo"
OUT.mkdir(parents=True, exist_ok=True)

UA = "fable-noskip/1.0 (research)"
PROCESSED_HINTS = (
    ".h5ad", ".h5", ".h5seurat", ".loom", ".mtx", ".mtx.gz",
    ".csv", ".csv.gz", ".tsv", ".tsv.gz", ".txt.gz",
    ".rds", ".RDS", ".parquet", ".feather",
)
ARCHIVE_HINTS = (".zip", ".tar.gz", ".tar.xz", ".tgz")
LUNG_HINTS = (
    "lung", "nsclc", "luad", "lusc", "pulmonary", "adenocarcinoma",
    "ici", "immunotherapy", "checkpoint", "anti-pd", "pd-1", "pd-l1",
    "tacstd2", "trop2", "cldn4", "claudin",
)


def human(n):
    if n is None:
        return "?"
    x = float(n)
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if x < 1024:
            return f"{x:.1f}{u}"
        x /= 1024
    return f"{x:.1f}PB"


def title_ok(title):
    t = (title or "").lower()
    return any(k in t for k in LUNG_HINTS)


def file_kind(name):
    n = (name or "").lower()
    if n.endswith(PROCESSED_HINTS):
        return "processed_table"
    if n.endswith(ARCHIVE_HINTS):
        return "archive"
    return "other"


def main():
    recs = json.loads((ROOT / "results" / "fable_zenodo" / "search_results.json").read_text())
    rows = []
    for r in recs:
        if str(r.get("access_right", "")).lower() in ("restricted", "closed", "private", "embargoed"):
            continue
        if not title_ok(r.get("title")):
            continue
        for f in r.get("files") or []:
            kind = file_kind(f.get("key"))
            if kind == "other":
                continue
            rows.append({
                "source": r["source"],
                "id": r["id"],
                "doi": r.get("doi"),
                "title": r.get("title"),
                "license": r.get("license"),
                "access_right": r.get("access_right"),
                "file": f.get("key"),
                "size": f.get("size"),
                "size_h": human(f.get("size")),
                "kind": kind,
                "link": f.get("link"),
            })
    rows.sort(key=lambda x: -(x["size"] or 0))
    (OUT / "hunt_candidates.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False))
    print(f"Candidates: {len(rows)}")
    for row in rows[:40]:
        print(f"  {row['size_h']:>9} [{row['source']} {row['id']}] {row['kind']:16} {row['file'][:50]}")


if __name__ == "__main__":
    main()
