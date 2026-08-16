#!/usr/bin/env python3
"""Live NCBI GEO E-utilities search for NEW 2025–2026 human lung IO scRNA-seq.

Every accession is taken from the live esearch/esummary JSON. Nothing is invented.
"""
import json
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = "results/hunt_geo2026_sc/geo_esearch_raw.json"

IO = (
    '("immune checkpoint"[All Fields] OR immunotherapy[All Fields] OR '
    'anti-PD-1[All Fields] OR "anti PD-1"[All Fields] OR anti-PD-L1[All Fields] OR '
    'pembrolizumab[All Fields] OR nivolumab[All Fields] OR atezolizumab[All Fields] OR '
    'durvalumab[All Fields] OR sintilimab[All Fields] OR tislelizumab[All Fields] OR '
    'camrelizumab[All Fields] OR toripalimab[All Fields] OR ipilimumab[All Fields] OR '
    'cemiplimab[All Fields] OR sugemalimab[All Fields] OR avelumab[All Fields] OR '
    '"PD-1 blockade"[All Fields] OR "checkpoint blockade"[All Fields] OR ICI[All Fields] OR '
    'ICB[All Fields] OR "immune checkpoint inhibitor"[All Fields])'
)
LUNG = (
    '(lung[All Fields] OR NSCLC[All Fields] OR LUAD[All Fields] OR LUSC[All Fields] OR '
    '"non-small cell"[All Fields] OR "small cell lung"[All Fields] OR SCLC[All Fields])'
)
SC = (
    '("single cell"[All Fields] OR single-cell[All Fields] OR scRNA-seq[All Fields] OR '
    'scRNAseq[All Fields] OR snRNA-seq[All Fields] OR "single nucleus"[All Fields] OR '
    '10x[All Fields] OR CITE-seq[All Fields])'
)
DATE = '("2025/01/01"[PDAT] : "2026/12/31"[PDAT])'
QUERY = f"{LUNG} AND {IO} AND {SC} AND {DATE} AND gse[ETYP] AND Homo sapiens[ORGN]"


def fetch(url, retries=5):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "hunt_geo2026_sc/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode()
        except Exception:
            if i == retries - 1:
                raise
            time.sleep(2 ** (i + 1))


def main():
    url = (f"{EUTILS}/esearch.fcgi?db=gds&term={urllib.parse.quote(QUERY)}"
           f"&retmax=0&usehistory=y&retmode=json")
    js = json.loads(fetch(url))["esearchresult"]
    count = int(js["count"])
    webenv, qk = js["webenv"], js["querykey"]
    print(f"esearch hit count: {count}")
    print(f"query: {QUERY}")

    docs = {}
    retstep = 200
    for start in range(0, count, retstep):
        url = (f"{EUTILS}/esummary.fcgi?db=gds&query_key={qk}&WebEnv={webenv}"
               f"&retstart={start}&retmax={retstep}&retmode=json")
        res = json.loads(fetch(url))["result"]
        for uid in res.get("uids", []):
            docs[uid] = res[uid]
        print(f"fetched {min(start + retstep, count)}/{count}")
        time.sleep(0.35)

    with open(OUT, "w") as f:
        json.dump({"query": QUERY, "count": count, "docs": docs}, f, indent=1)
    print(f"wrote {len(docs)} docsums to {OUT}")


if __name__ == "__main__":
    main()
