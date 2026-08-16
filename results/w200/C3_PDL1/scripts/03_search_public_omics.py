#!/usr/bin/env python3
"""Live search of public omics for C3 (CLDN4-loss or TROP2-ADC vs CD274).

Sources (all queried live, 2026-08-16):
  - NCBI GEO via esearch/esummary/elink
  - ArrayExpress / BioStudies
  - OmicsDI
  - SRA text search (presence only)

This script does not download matrices. It writes an inventory the
analyst then screens. Known C4 accessions are looked up explicitly so
they cannot be missed.

Outputs
-------
omics/search_log.md
omics/geo_hits.tsv
omics/biostudies_hits.tsv
omics/omicsdi_hits.tsv
omics/known_accessions.tsv
logs/omics_search_raw/*.json
"""

from __future__ import annotations

import json
import os
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RAW = os.path.join(ROOT, "logs", "omics_search_raw")
OMICS = os.path.join(ROOT, "omics")
UA = "C3-PDL1-catalog/1.0 (research; public omics inventory)"

KNOWN = [
    # from C4 work in this repo — must be re-read for CD274, not IFN
    ("GSE50927", "CLDN4_KD", "mouse lung Cldn4 KD; C4 analog that opened IFN"),
    ("GSE207704", "CLDN4_KD", "cancer-line CLDN4 perturbation; C4 analog that closed IFN"),
    ("GSE22493", "CLDN4_KD", "SKOV-3 CLDN4 siRNA vs overexpression; C4 analog, IFN not opened"),
    ("GSE22421", "CLDN4_LIGAND", "SKOV-3 C-CPE (CLDN3/4 binder), not genetic KD"),
    ("GSE274940", "CLDN_NULL", "EpH4 complete Cldn-null, not CLDN4-specific"),
    ("GSE245459", "TROP2_KD", "SKOV-3 TACSTD2/TROP2 shRNA, not CLDN4, not ADC"),
    ("GSE334497", "TROP2_KO", "4T1 Trop2 KO RNA-seq; C4 analog, IFN/APM not opened"),
]


def get(url: str, retries: int = 4) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    last = None
    for i in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90) as fh:
                return fh.read()
        except Exception as exc:
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"GET failed {url}: {last}")


def get_json(url: str):
    return json.loads(get(url).decode("utf-8", "replace"))


def geo_esearch(term: str, retmax: int = 80) -> list[str]:
    q = urllib.parse.quote(term)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=gds&term={q}&retmax={retmax}&retmode=json"
    )
    payload = get_json(url)
    time.sleep(0.35)
    return payload.get("esearchresult", {}).get("idlist", [])


def geo_esummary(ids: list[str]) -> list[dict]:
    if not ids:
        return []
    out = []
    for i in range(0, len(ids), 20):
        chunk = ",".join(ids[i:i + 20])
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=gds&id={chunk}&retmode=json"
        )
        payload = get_json(url)
        time.sleep(0.35)
        res = payload.get("result", {})
        for uid in res.get("uids", []):
            rec = res.get(uid, {})
            out.append({
                "uid": uid,
                "accession": rec.get("accession", ""),
                "title": rec.get("title", ""),
                "summary": rec.get("summary", ""),
                "gpl": rec.get("gpl", ""),
                "n_samples": rec.get("n_samples", ""),
                "gdstype": rec.get("gdstype", ""),
                "taxon": rec.get("taxon", ""),
                "pdat": rec.get("pdat", ""),
                "ftplink": rec.get("ftplink", ""),
            })
    return out


def biostudies(query: str, page_size: int = 50) -> list[dict]:
    q = urllib.parse.quote(query)
    url = (
        "https://www.ebi.ac.uk/biostudies/api/v1/search"
        f"?query={q}&pageSize={page_size}"
    )
    try:
        payload = get_json(url)
    except Exception as exc:
        return [{"error": str(exc), "query": query}]
    hits = []
    for rec in payload.get("hits", []) or payload.get("embedded", {}).get("hits", []) or []:
        hits.append({
            "accession": rec.get("accession") or rec.get("id", ""),
            "title": rec.get("title", ""),
            "type": rec.get("type", ""),
            "release": rec.get("releaseTime") or rec.get("release", ""),
            "files": rec.get("files", ""),
            "query": query,
        })
    return hits


def omicsdi(query: str, size: int = 40) -> list[dict]:
    q = urllib.parse.quote(query)
    url = f"https://www.omicsdi.org/ws/dataset/search?query={q}&size={size}"
    try:
        payload = get_json(url)
    except Exception as exc:
        return [{"error": str(exc), "query": query}]
    hits = []
    for rec in payload.get("datasets", []) or []:
        hits.append({
            "id": rec.get("id", ""),
            "source": rec.get("source", ""),
            "title": rec.get("title", ""),
            "description": (rec.get("description") or "")[:400],
            "query": query,
        })
    return hits


def sra_count(term: str) -> int:
    q = urllib.parse.quote(term)
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=sra&term={q}&retmax=0&retmode=json"
    )
    payload = get_json(url)
    time.sleep(0.35)
    return int(payload.get("esearchresult", {}).get("count", 0))


GEO_QUERIES = [
    ("G_CLDN4_KD", 'GSE[Entry Type] AND (CLDN4 OR "claudin 4" OR "claudin-4" OR Cldn4) AND (knockdown OR knockout OR siRNA OR shRNA OR CRISPR OR silencing)'),
    ("G_CLDN4_PDL1", 'GSE[Entry Type] AND (CLDN4 OR "claudin-4" OR Cldn4) AND (PD-L1 OR PDL1 OR CD274)'),
    ("G_CLDN4_RNA", 'GSE[Entry Type] AND (CLDN4[Title] OR "claudin-4"[Title] OR Cldn4[Title]) AND (Expression profiling by array[DataSet Type] OR Expression profiling by high throughput sequencing[DataSet Type])'),
    ("G_TROP2_ADC", 'GSE[Entry Type] AND (sacituzumab OR govitecan OR datopotamab OR "Dato-DXd" OR "DS-1062" OR SKB264 OR "MK-2870" OR "TROP2 ADC" OR "anti-TROP2")'),
    ("G_TROP2_PDL1", 'GSE[Entry Type] AND (TROP2 OR TACSTD2 OR "TROP-2") AND (PD-L1 OR PDL1 OR CD274)'),
    ("G_TROP2_KD", 'GSE[Entry Type] AND (TACSTD2 OR TROP2 OR "TROP-2") AND (knockdown OR knockout OR siRNA OR shRNA OR CRISPR)'),
    ("G_SN38_PDL1", 'GSE[Entry Type] AND (SN-38 OR SN38 OR irinotecan) AND (PD-L1 OR PDL1 OR CD274)'),
    ("G_CPE_CLDN", 'GSE[Entry Type] AND ("C-CPE" OR "clostridium perfringens enterotoxin") AND (CLDN4 OR claudin)'),
]


def write_tsv(path: str, rows: list[dict], cols: list[str]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ")[:500] for c in cols) + "\n")


def main() -> None:
    os.makedirs(RAW, exist_ok=True)
    os.makedirs(OMICS, exist_ok=True)
    log = []

    # GEO
    geo_all = {}
    for qid, term in GEO_QUERIES:
        sys.stderr.write(f"[GEO {qid}] {term}\n")
        ids = geo_esearch(term)
        recs = geo_esummary(ids)
        with open(os.path.join(RAW, f"geo_{qid}.json"), "w") as fh:
            json.dump({"qid": qid, "term": term, "n_ids": len(ids), "recs": recs}, fh, indent=1)
        log.append(f"- GEO `{qid}`: {len(ids)} uids / {len(recs)} summaries — `{term}`")
        for r in recs:
            key = r.get("accession") or r.get("uid")
            if key not in geo_all:
                geo_all[key] = dict(r)
                geo_all[key]["queries"] = []
            geo_all[key]["queries"].append(qid)
        sys.stderr.write(f"  {len(ids)} ids\n")

    # explicit known accessions
    known_rows = []
    for acc, kind, note in KNOWN:
        ids = geo_esearch(f"{acc}[Accession]")
        recs = geo_esummary(ids)
        with open(os.path.join(RAW, f"known_{acc}.json"), "w") as fh:
            json.dump({"accession": acc, "recs": recs}, fh, indent=1)
        title = recs[0]["title"] if recs else "NOT_FOUND"
        summary = recs[0]["summary"] if recs else ""
        known_rows.append({
            "accession": acc, "kind": kind, "note": note,
            "n_geo": len(recs), "title": title, "summary": summary[:400],
        })
        log.append(f"- known `{acc}` ({kind}): n={len(recs)} title={title[:80]}")

    write_tsv(
        os.path.join(OMICS, "geo_hits.tsv"),
        [{**v, "queries": "|".join(v.get("queries", []))} for v in geo_all.values()],
        ["accession", "uid", "queries", "taxon", "n_samples", "gdstype", "pdat", "gpl", "title", "summary"],
    )
    write_tsv(
        os.path.join(OMICS, "known_accessions.tsv"),
        known_rows,
        ["accession", "kind", "n_geo", "note", "title", "summary"],
    )

    # BioStudies
    bs_queries = [
        "CLDN4 AND (knockdown OR knockout OR siRNA OR CRISPR)",
        "sacituzumab OR datopotamab OR SKB264",
        "TROP2 AND (PD-L1 OR CD274)",
        "E-MTAB AND CLDN4",
    ]
    bs_rows = []
    for q in bs_queries:
        sys.stderr.write(f"[BioStudies] {q}\n")
        hits = biostudies(q)
        with open(os.path.join(RAW, f"biostudies_{abs(hash(q))}.json"), "w") as fh:
            json.dump({"query": q, "hits": hits}, fh, indent=1)
        log.append(f"- BioStudies `{q}`: {len(hits)} hits")
        for h in hits:
            h = dict(h)
            h["query"] = q
            bs_rows.append(h)
        time.sleep(0.3)
    write_tsv(
        os.path.join(OMICS, "biostudies_hits.tsv"),
        bs_rows,
        ["accession", "query", "type", "release", "title"],
    )

    # OmicsDI
    od_queries = [
        "CLDN4 knockdown",
        "CLDN4 knockout RNA-seq",
        "sacituzumab govitecan RNA-seq",
        "datopotamab deruxtecan",
        "TROP2 ADC PD-L1",
    ]
    od_rows = []
    for q in od_queries:
        sys.stderr.write(f"[OmicsDI] {q}\n")
        hits = omicsdi(q)
        with open(os.path.join(RAW, f"omicsdi_{abs(hash(q))}.json"), "w") as fh:
            json.dump({"query": q, "hits": hits}, fh, indent=1)
        log.append(f"- OmicsDI `{q}`: {len(hits)} hits")
        od_rows.extend(hits)
        time.sleep(0.3)
    write_tsv(
        os.path.join(OMICS, "omicsdi_hits.tsv"),
        od_rows,
        ["id", "source", "query", "title", "description"],
    )

    # SRA presence
    sra_terms = [
        "CLDN4 AND (knockout OR knockdown) AND RNA-seq",
        "Cldn4 AND (knockout OR knockdown) AND RNA-seq",
        "sacituzumab AND RNA-seq",
        "datopotamab AND RNA-seq",
        "SKB264 OR MK-2870",
    ]
    sra_rows = []
    for t in sra_terms:
        n = sra_count(t)
        sra_rows.append({"term": t, "count": n})
        log.append(f"- SRA `{t}`: {n}")
    write_tsv(os.path.join(OMICS, "sra_counts.tsv"), sra_rows, ["term", "count"])

    with open(os.path.join(OMICS, "search_log.md"), "w") as fh:
        fh.write("# Live public-omics search log (C3)\n\n")
        fh.write("Queried live. Nothing here is treated as supporting C3.\n\n")
        fh.write("\n".join(log) + "\n")
        fh.write(f"\nUnique GEO accessions after union: {len(geo_all)}\n")

    sys.stderr.write(f"unique GEO {len(geo_all)}\n")


if __name__ == "__main__":
    main()
