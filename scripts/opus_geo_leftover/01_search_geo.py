#!/usr/bin/env python3
"""Systematic GEO search for human lung-cancer immune-checkpoint-inhibitor bulk expression series.

Queries NCBI Entrez (gds) with several complementary term sets, merges the hits, and writes a
candidate table. Deliberately excludes the seven accessions already covered by sibling slices.
"""
from __future__ import annotations

import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parents[2] / "results" / "opus_geo_leftover"
OUT.mkdir(parents=True, exist_ok=True)

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TOOL = "opus_geo_leftover"
EMAIL = "cloud-agent@example.org"

EXCLUDED = {
    "GSE126044",
    "GSE135222",
    "GSE136961",
    "GSE166449",
    "GSE93157",
    "GSE207422",
    "GSE205335",
}

# Each query is intentionally over-broad; filtering happens downstream.
DISEASE = (
    '("lung"[All Fields] OR "NSCLC"[All Fields] OR "non-small cell lung"[All Fields] OR '
    '"adenocarcinoma of lung"[All Fields] OR "squamous cell lung"[All Fields] OR '
    '"SCLC"[All Fields] OR "small cell lung"[All Fields] OR "pulmonary"[All Fields])'
)
THERAPY = (
    '("immunotherapy"[All Fields] OR "immune checkpoint"[All Fields] OR "checkpoint '
    'inhibitor"[All Fields] OR "PD-1"[All Fields] OR "PD1"[All Fields] OR "PD-L1"[All Fields] OR '
    '"PDL1"[All Fields] OR "CTLA-4"[All Fields] OR "CTLA4"[All Fields] OR '
    '"nivolumab"[All Fields] OR "pembrolizumab"[All Fields] OR "atezolizumab"[All Fields] OR '
    '"durvalumab"[All Fields] OR "avelumab"[All Fields] OR "ipilimumab"[All Fields] OR '
    '"cemiplimab"[All Fields] OR "camrelizumab"[All Fields] OR "sintilimab"[All Fields] OR '
    '"tislelizumab"[All Fields] OR "toripalimab"[All Fields] OR "anti-PD-1"[All Fields] OR '
    '"anti-PD1"[All Fields] OR "ICI"[All Fields] OR "ICB"[All Fields])'
)
BASE_FILTER = (
    '"Homo sapiens"[Organism] AND "gse"[Entry Type]'
)

QUERIES = {
    "lung_ici_broad": f"{DISEASE} AND {THERAPY} AND {BASE_FILTER}",
    "lung_ici_response": (
        f"{DISEASE} AND {THERAPY} AND "
        '("response"[All Fields] OR "responder"[All Fields] OR "non-responder"[All Fields] OR '
        '"resistance"[All Fields] OR "efficacy"[All Fields] OR "outcome"[All Fields] OR '
        '"survival"[All Fields]) AND '
        f"{BASE_FILTER}"
    ),
    "lung_ici_rnaseq": (
        f"{DISEASE} AND {THERAPY} AND "
        '("Expression profiling by high throughput sequencing"[DataSet Type] OR '
        '"Expression profiling by array"[DataSet Type]) AND '
        f"{BASE_FILTER}"
    ),
    "nsclc_neoadjuvant": (
        f"{DISEASE} AND "
        '("neoadjuvant"[All Fields] OR "adjuvant"[All Fields] OR "perioperative"[All Fields]) AND '
        f"{THERAPY} AND {BASE_FILTER}"
    ),
}


def eget(endpoint: str, params: dict) -> str:
    params = {**params, "tool": TOOL, "email": EMAIL}
    url = f"{EUTILS}/{endpoint}?{urllib.parse.urlencode(params)}"
    last = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=120) as fh:
                return fh.read().decode("utf-8", "replace")
        except Exception as exc:  # noqa: BLE001 - transient NCBI failures are expected
            last = exc
            time.sleep(2 ** attempt)
    raise RuntimeError(f"E-utilities failed for {endpoint}: {last}")


def esearch(term: str, retmax: int = 1000) -> list[str]:
    raw = eget("esearch.fcgi", {"db": "gds", "term": term, "retmax": retmax, "retmode": "json"})
    return json.loads(raw)["esearchresult"].get("idlist", [])


def esummary(uids: list[str]) -> dict:
    out: dict = {}
    for i in range(0, len(uids), 200):
        chunk = uids[i : i + 200]
        raw = eget(
            "esummary.fcgi", {"db": "gds", "id": ",".join(chunk), "retmode": "json"}
        )
        payload = json.loads(raw).get("result", {})
        for uid in payload.get("uids", []):
            out[uid] = payload[uid]
        time.sleep(0.4)
    return out


def main() -> int:
    hits: dict[str, set[str]] = {}
    for name, term in QUERIES.items():
        uids = esearch(term)
        print(f"[search] {name}: {len(uids)} uids", file=sys.stderr)
        for uid in uids:
            hits.setdefault(uid, set()).add(name)
        time.sleep(0.5)

    summaries = esummary(sorted(hits))
    rows = []
    for uid, meta in summaries.items():
        acc = meta.get("accession", "")
        if not acc.startswith("GSE"):
            continue
        if acc in EXCLUDED:
            continue
        rows.append(
            {
                "accession": acc,
                "uid": uid,
                "title": (meta.get("title") or "").strip(),
                "summary": re.sub(r"\s+", " ", meta.get("summary") or "").strip(),
                "gdsType": meta.get("gdstype", ""),
                "n_samples": meta.get("n_samples"),
                "pdat": meta.get("pdat", ""),
                "taxon": meta.get("taxon", ""),
                "pubmed": ";".join(str(x) for x in (meta.get("pubmedids") or [])),
                "supp_file": meta.get("suppfile", ""),
                "ftplink": meta.get("ftplink", ""),
                "sample_titles": " | ".join(
                    (s.get("title") or "") for s in (meta.get("samples") or [])[:40]
                ),
                "queries": ";".join(sorted(hits[uid])),
            }
        )

    rows.sort(key=lambda r: r["accession"])
    (OUT / "geo_search_candidates.json").write_text(json.dumps(rows, indent=2))
    print(f"[search] wrote {len(rows)} GSE candidates", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
