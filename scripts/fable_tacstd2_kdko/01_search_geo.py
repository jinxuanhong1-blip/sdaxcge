#!/usr/bin/env python3
"""Exhaustively search NCBI GEO for TACSTD2/TROP2 loss-of-function datasets.

Casts a wide net across several query variants against the GEO DataSets (gds)
database, pulls esummary records, dedupes to GSE series, and writes a candidate
table. No filtering on "KD/KO" is done here beyond the query terms; triage of
which series are true perturbation experiments is done manually downstream.

Outputs:
  results/fable_tacstd2_kdko/geo_search_raw.json   (all esummary docs, deduped)
  results/fable_tacstd2_kdko/geo_search_candidates.tsv
"""
import json
import time
import sys
import urllib.parse
import urllib.request
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OUT = Path(__file__).resolve().parents[2] / "results" / "fable_tacstd2_kdko"
OUT.mkdir(parents=True, exist_ok=True)

# Wide net of query variants. Field-restricted to Title/Description so we catch
# datasets *about* TACSTD2/TROP2 perturbation rather than any dataset that merely
# measures the gene.
QUERIES = [
    "(TACSTD2[Title] OR TROP2[Title]) AND (knockdown OR knockout OR shRNA OR siRNA OR sgRNA OR CRISPR OR silencing OR depletion OR knock-out OR knock-down OR deletion OR loss)",
    "(TACSTD2[Description] OR TROP2[Description]) AND (knockdown OR knockout OR shRNA OR siRNA OR sgRNA OR CRISPR OR silencing OR depletion OR knock-out OR knock-down OR deletion)",
    "TACSTD2[Title] AND (knockdown OR knockout OR shRNA OR siRNA OR sgRNA OR CRISPR OR silencing OR depletion OR overexpression OR KO OR KD)",
    "TROP2[Title] AND (knockdown OR knockout OR shRNA OR siRNA OR sgRNA OR CRISPR OR silencing OR depletion OR overexpression OR KO OR KD)",
    "(TACSTD2[Description] OR TROP2[Description]) AND (knockdown OR knockout OR shRNA OR siRNA OR sgRNA OR CRISPR OR silencing OR depletion OR overexpression)",
    "(Trop2[Title] OR Tacstd2[Title]) AND (knockout OR deficient OR null OR mutant OR CRISPR)",  # mouse
]


def eget(endpoint, params):
    url = f"{EUTILS}/{endpoint}?" + urllib.parse.urlencode(params)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return r.read().decode()
        except Exception as e:  # noqa
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)
    return None


def esearch_uids(term):
    txt = eget("esearch.fcgi", {"db": "gds", "term": term, "retmax": 500, "retmode": "json"})
    data = json.loads(txt)
    return data["esearchresult"].get("idlist", [])


def esummary(uids):
    if not uids:
        return {}
    txt = eget("esummary.fcgi", {"db": "gds", "id": ",".join(uids), "retmode": "json"})
    data = json.loads(txt)
    return data.get("result", {})


def main():
    all_uids = set()
    per_query = {}
    for q in QUERIES:
        uids = esearch_uids(q)
        per_query[q] = uids
        all_uids.update(uids)
        print(f"[search] {len(uids):4d} uids :: {q}", file=sys.stderr)
        time.sleep(0.4)

    print(f"[search] total unique uids: {len(all_uids)}", file=sys.stderr)

    # Fetch summaries in batches
    uids = sorted(all_uids, key=int)
    docs = {}
    for i in range(0, len(uids), 100):
        batch = uids[i:i + 100]
        res = esummary(batch)
        for k, v in res.items():
            if k == "uids":
                continue
            docs[k] = v
        time.sleep(0.4)

    (OUT / "geo_search_raw.json").write_text(json.dumps(docs, indent=2))

    # Build candidate table
    rows = []
    for uid, d in docs.items():
        acc = d.get("accession", "")
        rows.append({
            "uid": uid,
            "accession": acc,
            "entrytype": d.get("entrytype", ""),
            "gdstype": d.get("gdstype", ""),
            "title": d.get("title", "").replace("\t", " ").replace("\n", " "),
            "taxon": d.get("taxon", ""),
            "n_samples": d.get("n_samples", ""),
            "gpl": d.get("gpl", ""),
            "suppfile": d.get("suppfile", "").replace("\t", " "),
            "summary": d.get("summary", "").replace("\t", " ").replace("\n", " ")[:600],
        })
    # Prefer GSE series
    rows.sort(key=lambda r: (r["entrytype"] != "GSE", r["accession"]))

    tsv = OUT / "geo_search_candidates.tsv"
    cols = ["accession", "entrytype", "taxon", "n_samples", "gdstype", "gpl", "suppfile", "title"]
    with tsv.open("w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"[out] {len(rows)} records -> {tsv}", file=sys.stderr)
    n_gse = sum(1 for r in rows if r["entrytype"] == "GSE")
    print(f"[out] GSE series: {n_gse}", file=sys.stderr)


if __name__ == "__main__":
    main()
