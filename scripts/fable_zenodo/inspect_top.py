#!/usr/bin/env python3
"""Print detailed metadata for top shortlist candidates that have downloadable matrices."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "fable_zenodo"

scored = json.loads((RESULTS / "shortlist.json").read_text())


def human(n):
    if n is None:
        return "?"
    for u in ["B", "KB", "MB", "GB"]:
        if n < 1024:
            return f"{n:.1f}{u}"
        n /= 1024
    return f"{n:.1f}TB"


cands = [r for r in scored if r.get("has_downloadable_matrix")]
cands.sort(key=lambda r: r["rel_score"], reverse=True)
for r in cands[:25]:
    print("=" * 90)
    print(f"[{r['source']}] rel={r['rel_score']} open={r['open']} access={r.get('access_right')}")
    print(f"TITLE: {r.get('title')}")
    print(f"DOI:   {r.get('doi')}  | {r.get('doi_url')}")
    print(f"LICENSE: {r.get('license')}  DATE: {r.get('publication_date')}")
    print(f"HTML:  {r.get('html')}")
    print(f"hits: {r.get('rel_hits')}")
    for f in r.get("downloadable_matrix_files", [])[:12]:
        print(f"   - {f.get('key')}  [{human(f.get('size'))}]")
    if len(r.get("downloadable_matrix_files", [])) > 12:
        print(f"   ... (+{len(r['downloadable_matrix_files'])-12} more)")
