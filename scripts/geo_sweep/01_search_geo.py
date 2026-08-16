#!/usr/bin/env python3
"""
Step 1 of the GEO remaining-series sweep.

Searches NCBI GEO DataSets (db=gds) for OPEN human/mouse *lung* Series (GSE)
that textually mention any of the target terms:
    TACSTD2, TROP2, CLDN4, claudin-4, ICI

Terms are quoted so we match literal text mentions (title/summary/sample
characteristics) rather than platform gene-annotation presence. The classic
reference set is excluded. Results are written to a raw metadata JSON that the
later catalog / download / analysis steps consume.

Only the NCBI E-utilities HTTP API is used (no external python deps beyond stdlib).
Outputs go ONLY under notes/geo_sweep/.
"""
import json
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NOTES = Path("notes/geo_sweep")
NOTES.mkdir(parents=True, exist_ok=True)

# The 5 target terms exactly as given in the task. Quoted to force literal-text
# matching and avoid the gene-symbol platform-annotation false positives that
# an unquoted TACSTD2 search produces (~338 platform hits vs 4 real mentions).
TERMS = ['"TACSTD2"', '"TROP2"', '"CLDN4"', '"claudin-4"', '"ICI"']

ORG_FILTER = '("Homo sapiens"[Organism] OR "Mus musculus"[Organism])'
BASE_FILTER = f'lung AND {ORG_FILTER} AND gse[ETYP]'

# Classic reference set to exclude.
CLASSIC = {
    "GSE126044", "GSE135222", "GSE136961", "GSE166449",
    "GSE93157", "GSE207422", "GSE205335",
}


def http_get(url, tries=5):
    last = None
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed after {tries} tries: {url}\n{last}")


def esearch(term, filt):
    q = f"{term} AND {filt}"
    url = (
        f"{EUTILS}/esearch.fcgi?db=gds&retmax=1000&term="
        + urllib.parse.quote(q)
    )
    data = http_get(url)
    root = ET.fromstring(data)
    ids = [e.text for e in root.findall(".//IdList/Id")]
    count = root.findtext("Count")
    return ids, int(count or 0)


def esummary(uids):
    """Fetch DocSummaries for a batch of gds UIDs, return list of dicts."""
    out = []
    for start in range(0, len(uids), 200):
        chunk = uids[start:start + 200]
        url = (
            f"{EUTILS}/esummary.fcgi?db=gds&version=2.0&id="
            + ",".join(chunk)
        )
        data = http_get(url)
        root = ET.fromstring(data)
        for ds in root.findall(".//DocumentSummary"):
            rec = {}
            for child in ds:
                if list(child):
                    # nested (e.g. Samples); capture count only
                    rec[child.tag] = [
                        {gc.tag: gc.text for gc in item} for item in child
                    ]
                else:
                    rec[child.tag] = child.text
            rec["_uid"] = ds.get("uid")
            out.append(rec)
        time.sleep(0.4)
    return out


def main():
    term_to_uids = {}
    all_uids = set()
    for term in TERMS:
        ids, count = esearch(term, BASE_FILTER)
        term_to_uids[term] = ids
        all_uids.update(ids)
        print(f"{term:16s} -> {count} hits ({len(ids)} ids returned)")
        time.sleep(0.4)

    all_uids = sorted(all_uids)
    print(f"\nUnique UIDs across all terms: {len(all_uids)}")

    summ = esummary(all_uids)
    by_uid = {r["_uid"]: r for r in summ}

    # attach matched terms per record
    for term, ids in term_to_uids.items():
        clean = term.strip('"')
        for uid in ids:
            if uid in by_uid:
                by_uid[uid].setdefault("_matched_terms", []).append(clean)

    # normalise / filter to GSE series, drop classic set
    records = []
    excluded_classic = []
    for uid, r in by_uid.items():
        acc = r.get("Accession") or ""
        if not acc.startswith("GSE"):
            continue
        if acc in CLASSIC:
            excluded_classic.append(acc)
            continue
        records.append(r)

    records.sort(key=lambda r: r.get("Accession", ""))
    print(f"GSE records (excl. classic): {len(records)}")
    print(f"Classic excluded encountered: {sorted(set(excluded_classic))}")

    out = {
        "terms": TERMS,
        "base_filter": BASE_FILTER,
        "classic_excluded_set": sorted(CLASSIC),
        "classic_encountered": sorted(set(excluded_classic)),
        "n_records": len(records),
        "records": records,
    }
    outfile = NOTES / "search_raw.json"
    outfile.write_text(json.dumps(out, indent=2))
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
