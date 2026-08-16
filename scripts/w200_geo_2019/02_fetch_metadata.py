#!/usr/bin/env python3
"""Fetch esummary + SOFT-header metadata for every UID found by 01_search_geo.py.

esummary gives title/summary/taxon/platform/n_samples/PubMed ids/supp files.
The SOFT "brief" header additionally gives the true GEO submission date, the
public-release date, the overall design and the per-series supplementary file
list, which is what the triage in 03 actually needs.
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "w200" / "GEO_2019" / "search"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def _get(url, timeout=120, tries=5):
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            wait = 2 ** attempt
            print(f"  retry {attempt} in {wait}s ({exc})", file=sys.stderr)
            time.sleep(wait)
    return b""


def esummary(uids):
    out = {}
    for i in range(0, len(uids), 100):
        chunk = uids[i:i + 100]
        url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"})
        data = json.loads(_get(url) or b"{}")
        res = data.get("result", {})
        for uid in res.get("uids", []):
            out[uid] = res[uid]
        time.sleep(0.4)
    return out


def soft_header(gse):
    """Series-level SOFT header (targ=self, brief) -> dict of !Series_* fields."""
    url = f"{GEO_ACC}?" + urllib.parse.urlencode(
        {"acc": gse, "targ": "self", "form": "text", "view": "brief"})
    txt = (_get(url) or b"").decode("utf-8", "replace")
    fields = {}
    for line in txt.splitlines():
        m = re.match(r"^!(\w+)\s*=\s*(.*)$", line.strip())
        if m:
            fields.setdefault(m.group(1), []).append(m.group(2).strip())
    return fields


def main():
    search = json.loads((OUT_DIR / "search_uids.json").read_text())
    uids = search["uids"]
    print(f"esummary for {len(uids)} UIDs ...")
    summ = esummary(uids)

    records = []
    for uid, s in summ.items():
        gse = "GSE" + str(s.get("gse", "")).strip()
        if not s.get("gse"):
            continue
        print(f"  SOFT {gse}")
        soft = soft_header(gse)
        records.append({
            "uid": uid,
            "accession": gse,
            "title": s.get("title", ""),
            "summary": s.get("summary", ""),
            "gdsType": s.get("gdstype", ""),
            "taxon": s.get("taxon", ""),
            "n_samples": s.get("n_samples", 0),
            "platform": s.get("gpl", ""),
            "pubmed": ";".join(str(p) for p in s.get("pubmedids", []) or []),
            "pdat": s.get("pdat", ""),
            "ftplink": s.get("ftplink", ""),
            "suppfile": s.get("suppfile", ""),
            "soft_submission_date": "; ".join(soft.get("Series_submission_date", [])),
            "soft_status": "; ".join(soft.get("Series_status", [])),
            "soft_last_update": "; ".join(soft.get("Series_last_update_date", [])),
            "soft_type": "; ".join(soft.get("Series_type", [])),
            "soft_overall_design": " ".join(soft.get("Series_overall_design", [])),
            "soft_supplementary_file": "; ".join(soft.get("Series_supplementary_file", [])),
            "soft_platform_id": "; ".join(soft.get("Series_platform_id", [])),
            "soft_sample_count": len(soft.get("Series_sample_id", [])),
        })
        time.sleep(0.35)

    records.sort(key=lambda r: r["accession"])
    (OUT_DIR / "candidates_metadata.json").write_text(json.dumps(records, indent=2))
    print(f"\nWrote {OUT_DIR / 'candidates_metadata.json'} ({len(records)} series)")


if __name__ == "__main__":
    main()
