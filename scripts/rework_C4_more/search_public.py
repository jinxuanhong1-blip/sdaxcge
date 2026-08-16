#!/usr/bin/env python3
"""Live hunt for public CLDN4/Cldn4 KD/KO transcriptomes.

Writes JSON dumps under results/rework/C4_more/. Does not invent accessions.
Known-excluded GSE IDs are listed only so they can be filtered, not reused.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/rework/C4_more")
EXCLUDED = {"GSE50927", "GSE207704", "GSE22493", "GSE245459"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch(url: str, retries: int = 4) -> str:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "C4-more-hunt/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8")
        except Exception as e:  # noqa: BLE001
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"fetch failed: {url}: {last}")


def esearch(db: str, term: str, retmax: int = 200) -> dict:
    q = urllib.parse.quote(term)
    url = f"{EUTILS}/esearch.fcgi?db={db}&term={q}&retmax={retmax}&retmode=json"
    return json.loads(fetch(url))


def esummary(db: str, ids: list[str]) -> dict:
    out: dict = {"result": {"uids": []}}
    if not ids:
        return out
    for i in range(0, len(ids), 12):
        batch = ",".join(ids[i : i + 12])
        url = f"{EUTILS}/esummary.fcgi?db={db}&id={batch}&retmode=json"
        chunk = json.loads(fetch(url))
        out["result"]["uids"].extend(chunk.get("result", {}).get("uids", []))
        for uid in chunk.get("result", {}).get("uids", []):
            out["result"][uid] = chunk["result"][uid]
        time.sleep(0.35)
    return out


def slim_gds(rec: dict) -> dict:
    samples = rec.get("samples") or []
    return {
        "accession": rec.get("accession"),
        "title": rec.get("title"),
        "summary": rec.get("summary"),
        "gdstype": rec.get("gdstype"),
        "taxon": rec.get("taxon"),
        "n_samples": rec.get("n_samples"),
        "pdat": rec.get("pdat"),
        "pubmedids": rec.get("pubmedids"),
        "bioproject": rec.get("bioproject"),
        "sample_titles": [s.get("title") for s in samples],
        "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={rec.get('accession')}",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    queries = {
        "CLDN4_gse": "CLDN4[All Fields] AND gse[Entry Type]",
        "Cldn4_gse": "Cldn4[All Fields] AND gse[Entry Type]",
        "claudin4_quoted_gse": '"claudin 4"[All Fields] AND gse[Entry Type]',
        "claudin-4_gse": "claudin-4[All Fields] AND gse[Entry Type]",
        "CLDN4_kdko_gse": (
            "CLDN4[All Fields] AND (knockdown OR knockout OR siRNA OR shRNA "
            "OR CRISPR OR silencing) AND gse[Entry Type]"
        ),
        "Cldn4_kdko_gse": (
            "Cldn4[All Fields] AND (knockdown OR knockout OR KO OR siRNA OR "
            "shRNA OR CRISPR) AND gse[Entry Type]"
        ),
        "Cldn-null_gse": "Cldn-null AND gse[Entry Type]",
        "GSE22421": "GSE22421[All Fields]",
        "GSE274940": "GSE274940[All Fields]",
        "GSE245459": "GSE245459[All Fields]",
        "GSE334497": "GSE334497[All Fields]",
    }
    search_dump = {}
    all_uids: set[str] = set()
    for name, term in queries.items():
        rec = esearch("gds", term)
        search_dump[name] = {
            "term": term,
            "count": rec.get("esearchresult", {}).get("count"),
            "idlist": rec.get("esearchresult", {}).get("idlist", []),
            "querytranslation": rec.get("esearchresult", {}).get("querytranslation"),
        }
        for uid in rec.get("esearchresult", {}).get("idlist", []):
            if str(uid).startswith("200"):  # GSE series UIDs
                all_uids.add(str(uid))
        time.sleep(0.35)

    summaries = esummary("gds", sorted(all_uids))
    series = []
    for uid in summaries.get("result", {}).get("uids", []):
        rec = summaries["result"][uid]
        if rec.get("entrytype") != "GSE":
            continue
        series.append(slim_gds(rec))

    pubmed_terms = {
        "Cldn4_KO_RNAseq": "Cldn4 knockout RNA-seq",
        "CLDN4_CRISPR_RNAseq": "CLDN4 CRISPR RNA-seq",
        "claudin4_silencing_transcriptome": "claudin-4 silencing transcriptome",
        "CLDN4_KD_GEO": (
            'CLDN4 AND (knockdown OR knockout) AND (RNA-seq OR microarray) '
            'AND ("Gene Expression Omnibus" OR GSE OR ArrayExpress)'
        ),
    }
    pubmed = {}
    pmids: list[str] = []
    for name, term in pubmed_terms.items():
        rec = esearch("pubmed", term, retmax=30)
        pubmed[name] = {
            "term": term,
            "count": rec.get("esearchresult", {}).get("count"),
            "idlist": rec.get("esearchresult", {}).get("idlist", []),
        }
        pmids.extend(rec.get("esearchresult", {}).get("idlist", []))
        time.sleep(0.35)
    pm_sum = esummary("pubmed", sorted(set(pmids)))
    papers = []
    for uid in pm_sum.get("result", {}).get("uids", []):
        rec = pm_sum["result"][uid]
        papers.append(
            {
                "pmid": uid,
                "title": rec.get("title"),
                "source": rec.get("source"),
                "pubdate": rec.get("pubdate"),
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/",
            }
        )

    sra = {}
    for name, term in {
        "CLDN4_KO_RNAseq": "CLDN4 knockout RNA-seq",
        "Cldn4_KO_transcriptome": "Cldn4 knockout transcriptome",
        "CLDN4_KD_RNAseq": "CLDN4 knockdown RNA-seq",
        "CLDN4_shRNA_RNAseq": "CLDN4 shRNA RNA-seq",
    }.items():
        rec = esearch("sra", term, retmax=20)
        sra[name] = {
            "term": term,
            "count": rec.get("esearchresult", {}).get("count"),
            "idlist": rec.get("esearchresult", {}).get("idlist", []),
        }
        time.sleep(0.35)

    payload = {
        "audit_date": "2026-08-16",
        "excluded_known": sorted(EXCLUDED),
        "geo_queries": search_dump,
        "geo_series": series,
        "pubmed_queries": pubmed,
        "pubmed_hits": papers,
        "sra_queries": sra,
        "verdict": "none",
    }
    (OUT / "live_search.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUT / 'live_search.json'} n_series={len(series)} n_papers={len(papers)}")


if __name__ == "__main__":
    main()
