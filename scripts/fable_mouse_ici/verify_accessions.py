#!/usr/bin/env python3
"""Verify GEO / ArrayExpress accessions for the mouse-lung ICI RNA slice.

For each GEO series this queries the public GEO ``acc.cgi`` endpoint
(``targ=self&form=text``) which returns a MINiML-flavoured plain-text record
containing the series title, organism(s), summary/overall-design, and the list
of supplementary files. For the ArrayExpress/BioStudies accession we hit the
BioStudies JSON API.

Output: notes/fable_mouse_ici/accession_verification.json (machine readable) and
a short human-readable log to stdout.
"""
import json
import re
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

GEO_ACCS = [
    "GSE239485", "GSE133604", "GSE222158", "GSE129297", "GSE241978",
    "GSE330658", "GSE197260", "GSE297630", "GSE297632",
]
AE_ACCS = ["E-MTAB-13704"]

NOTES = Path(__file__).resolve().parents[2] / "notes" / "fable_mouse_ici"
NOTES.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0 (fable-mouse-ici verifier)"}


def fetch(url, tries=4):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 ** i)
    raise last


def parse_geo_text(txt):
    """Parse the GEO ``form=text`` record into a dict of lists keyed by field."""
    fields = {}
    for line in txt.splitlines():
        m = re.match(r"^!(\w+?)_(\w+) = (.*)$", line)
        if not m:
            continue
        key = f"{m.group(2)}"
        fields.setdefault(key, []).append(m.group(3))
    return fields


def verify_geo(acc):
    url = (
        "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
        f"acc={acc}&targ=self&form=text&view=quick"
    )
    rec = {"accession": acc, "source": "GEO", "url": url}
    try:
        txt = fetch(url)
    except Exception as e:  # noqa: BLE001
        rec["error"] = f"fetch failed: {e}"
        return rec
    if "Could not find" in txt or "no items found" in txt.lower() or not txt.strip():
        rec["exists"] = False
        rec["raw_head"] = txt[:300]
        return rec
    f = parse_geo_text(txt)
    rec["exists"] = True
    rec["title"] = " | ".join(f.get("title", []))
    rec["organism"] = f.get("sample_organism", []) or f.get("platform_organism", [])
    rec["summary"] = " ".join(f.get("summary", []))
    rec["overall_design"] = " ".join(f.get("overall_design", []))
    rec["type"] = f.get("type", [])
    supp = f.get("supplementary_file", [])
    rec["supplementary_files"] = supp
    rec["n_samples"] = len(f.get("sample_id", []))
    return rec


def verify_ae(acc):
    """BioStudies/ArrayExpress: fetch study info + file list."""
    rec = {"accession": acc, "source": "ArrayExpress/BioStudies"}
    info_url = f"https://www.ebi.ac.uk/biostudies/api/v1/studies/{acc}/info"
    try:
        info = json.loads(fetch(info_url))
        rec["exists"] = True
        rec["title"] = info.get("title")
        rec["n_files"] = info.get("files")
        rec["file_size_bytes"] = info.get("fileSize")
        rec["release_date"] = info.get("released")
    except Exception as e:  # noqa: BLE001
        rec["exists"] = False
        rec["error"] = str(e)
    return rec


def main():
    results = []
    for acc in GEO_ACCS:
        print(f"[verify] {acc} ...", file=sys.stderr)
        rec = verify_geo(acc)
        results.append(rec)
        time.sleep(0.5)
    for acc in AE_ACCS:
        print(f"[verify] {acc} ...", file=sys.stderr)
        results.append(verify_ae(acc))
        time.sleep(0.5)

    out = NOTES / "accession_verification.json"
    out.write_text(json.dumps(results, indent=2))
    print(f"wrote {out}")

    for r in results:
        print("=" * 70)
        print(r.get("accession"), "| exists:", r.get("exists"))
        print("  title:", (r.get("title") or "")[:140])
        print("  organism:", r.get("organism"))
        print("  type:", r.get("type"))
        nsupp = len(r.get("supplementary_files", []) or [])
        print("  n_samples:", r.get("n_samples"), "| n_supp_files:", nsupp)
        for s in (r.get("supplementary_files") or []):
            print("     supp:", s)


if __name__ == "__main__":
    main()
