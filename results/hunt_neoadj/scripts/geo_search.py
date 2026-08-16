#!/usr/bin/env python3
"""Search GEO (GDS) for public lung neoadjuvant IO scRNA/bulk 2023-2026.

Uses NCBI E-utilities only. Does not invent accessions: every GSE is returned
by esearch/esummary. Writes raw JSON + a curated candidate table.
"""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/hunt_neoadj/metadata")
OUT.mkdir(parents=True, exist_ok=True)

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Date window requested by the hunt. Seeds older than 2023 are still recorded
# separately; this search is the 2023-2026 sweep.
QUERIES = {
    "scrna_neoadj_lung_io_2023_2026": (
        '('
        'lung[Title] OR NSCLC[Title] OR "non-small cell lung"[Title] OR '
        '"non small cell lung"[Title] OR LUAD[Title] OR LUSC[Title]'
        ') AND ('
        'neoadjuvant OR "neo-adjuvant" OR "preoperative immunotherapy" OR '
        '"induction immunotherapy" OR "pre-operative"'
        ') AND ('
        'immunotherapy OR "PD-1" OR "PD-L1" OR PD1 OR PDL1 OR pembrolizumab OR '
        'nivolumab OR atezolizumab OR camrelizumab OR sintilimab OR toripalimab OR '
        'tislelizumab OR durvalumab OR "immune checkpoint" OR "checkpoint blockade" OR ICI'
        ') AND ('
        '"single-cell" OR "single cell" OR scRNA OR scRNAseq OR "scRNA-seq" OR 10x OR '
        '"single-cell RNA"'
        ') AND 2023:2026[PDAT]'
    ),
    "bulk_or_scrna_neoadj_lung_io_2023_2026": (
        '('
        'lung[Title] OR NSCLC[Title] OR "non-small cell lung"[Title] OR '
        '"non small cell lung"[Title] OR LUAD[Title] OR LUSC[Title]'
        ') AND ('
        'neoadjuvant OR "neo-adjuvant" OR "induction immunotherapy"'
        ') AND ('
        'immunotherapy OR "PD-1" OR "PD-L1" OR PD1 OR pembrolizumab OR nivolumab OR '
        'atezolizumab OR camrelizumab OR sintilimab OR toripalimab OR tislelizumab OR '
        'durvalumab OR "immune checkpoint" OR ICI'
        ') AND ('
        'RNA-seq OR RNAseq OR transcriptome OR "expression profiling" OR "single-cell" OR scRNA'
        ') AND 2023:2026[PDAT]'
    ),
    "scrna_lung_io_response_2023_2026": (
        '('
        'lung[Title] OR NSCLC[Title] OR "non-small cell lung"[Title]'
        ') AND ('
        'immunotherapy OR "PD-1" OR pembrolizumab OR nivolumab OR "checkpoint"'
        ') AND ('
        '"single-cell" OR scRNA OR "scRNA-seq"'
        ') AND ('
        'MPR OR pCR OR "pathologic response" OR responder OR neoadjuvant'
        ') AND 2023:2026[PDAT]'
    ),
}

# Also fetch the four seed accessions via accession search (existence check).
SEEDS = ["GSE207422", "GSE243013", "GSE146100", "GSE229353"]


def get(url: str) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "hunt-neoadj/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read()


def esearch(term: str, retmax: int = 200) -> list[str]:
    q = urllib.parse.urlencode(
        {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"}
    )
    data = json.loads(get(f"{ESEARCH}?{q}"))
    return data["esearchresult"].get("idlist", [])


def esummary(ids: list[str]) -> dict:
    if not ids:
        return {}
    out = {}
    for i in range(0, len(ids), 40):
        chunk = ids[i : i + 40]
        q = urllib.parse.urlencode(
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"}
        )
        data = json.loads(get(f"{ESUMMARY}?{q}"))
        out.update(data.get("result", {}))
        time.sleep(0.35)
    return out


def main():
    all_ids = set()
    query_hits = {}
    for name, term in QUERIES.items():
        ids = esearch(term)
        query_hits[name] = {"term": term, "n": len(ids), "ids": ids}
        all_ids.update(ids)
        print(f"{name}: {len(ids)} hits", flush=True)
        time.sleep(0.35)

    seed_hits = {}
    for acc in SEEDS:
        ids = esearch(f"{acc}[ACCN]")
        seed_hits[acc] = ids
        all_ids.update(ids)
        print(f"seed {acc}: ids={ids}", flush=True)
        time.sleep(0.35)

    summaries = esummary(sorted(all_ids))
    records = []
    for uid, rec in summaries.items():
        if uid == "uids" or not isinstance(rec, dict):
            continue
        acc = rec.get("accession") or ""
        records.append(
            {
                "uid": uid,
                "accession": acc,
                "title": rec.get("title"),
                "summary": rec.get("summary"),
                "gdstype": rec.get("gdstype"),
                "n_samples": rec.get("n_samples"),
                "taxon": rec.get("taxon"),
                "pdat": rec.get("pdat"),
                "ptech": rec.get("ptech"),
                "gpl": rec.get("gpl"),
                "gse": rec.get("gse"),
                "pubmedids": rec.get("pubmedids"),
                "ftp": rec.get("ftplink"),
            }
        )
    records.sort(key=lambda r: (r.get("pdat") or "", r.get("accession") or ""))

    payload = {
        "queries": query_hits,
        "seed_hits": seed_hits,
        "n_unique_uids": len(all_ids),
        "records": records,
    }
    (OUT / "geo_search_raw.json").write_text(json.dumps(payload, indent=2))

    # Compact TSV of GSE-level records only
    gses = [r for r in records if str(r.get("accession", "")).startswith("GSE")]
    lines = ["accession\tpdat\tn_samples\tgdstype\ttitle"]
    for r in gses:
        title = (r.get("title") or "").replace("\t", " ").replace("\n", " ")
        lines.append(
            f"{r['accession']}\t{r.get('pdat')}\t{r.get('n_samples')}\t{r.get('gdstype')}\t{title}"
        )
    (OUT / "geo_search_gse.tsv").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(gses)} GSE records")


if __name__ == "__main__":
    main()
