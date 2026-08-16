#!/usr/bin/env python3
"""Fetch GEO series metadata (esummary) for every 2020-search UID."""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT_DIR = ROOT / "results" / "w200" / "GEO_2020"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def esummary(uids):
    params = {"db": "gds", "id": ",".join(uids), "retmode": "json"}
    url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(6):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                return json.load(r)
        except Exception as e:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  esummary retry {attempt} after {wait}s ({e})", file=sys.stderr)
            time.sleep(wait)
    return {}


def main():
    data = json.loads((OUT_DIR / "search_uids.json").read_text())
    uids = data["uids"]
    records = []
    for i in range(0, len(uids), 20):
        batch = uids[i:i + 20]
        res = esummary(batch)
        result = res.get("result", {})
        for uid in result.get("uids", []):
            r = result[uid]
            records.append({
                "uid": uid,
                "accession": r.get("accession"),
                "title": r.get("title"),
                "summary": r.get("summary"),
                "taxon": r.get("taxon"),
                "n_samples": r.get("n_samples"),
                "gdsType": r.get("gdsType"),
                "gpl": r.get("gpl"),
                "pdat": r.get("pdat"),
                "suppFile": r.get("suppFile"),
                "ftplink": r.get("ftplink"),
                "pubmedids": r.get("pubmedids"),
                "seriestitle": r.get("seriestitle"),
                "entrytype": r.get("entrytype"),
            })
        print(f"fetched {min(i + 20, len(uids))}/{len(uids)}")
        time.sleep(0.35)

    records.sort(key=lambda x: x.get("accession") or "")
    (OUT_DIR / "candidates_metadata.json").write_text(json.dumps(records, indent=2))
    print(f"Wrote {len(records)} records to candidates_metadata.json")


if __name__ == "__main__":
    main()
