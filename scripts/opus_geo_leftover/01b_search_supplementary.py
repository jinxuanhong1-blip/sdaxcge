#!/usr/bin/env python3
"""Second-pass GEO search: trial-oriented and outcome-oriented phrasings the first pass can miss.

Some lung ICI series describe themselves only via a trial name, a "biomarker" framing, or an
outcome word without ever using the words "lung" and "immunotherapy" together in the summary.
Hits from here are merged into the same candidate pool used by the triage step.
"""
from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"

spec = importlib.util.spec_from_file_location(
    "s1", Path(__file__).with_name("01_search_geo.py")
)
s1 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s1)

EXTRA_QUERIES = {
    "trials": (
        '("CheckMate"[All Fields] OR "KEYNOTE"[All Fields] OR "IMpower"[All Fields] OR '
        '"OAK trial"[All Fields] OR "POPLAR"[All Fields] OR "PACIFIC"[All Fields] OR '
        '"NADIM"[All Fields] OR "CheckMate 816"[All Fields] OR "LCMC3"[All Fields]) AND '
        '"Homo sapiens"[Organism] AND "gse"[Entry Type]'
    ),
    "biomarker_framing": (
        '("lung"[All Fields] OR "NSCLC"[All Fields]) AND '
        '("predictive biomarker"[All Fields] OR "predict response"[All Fields] OR '
        '"biomarker of response"[All Fields] OR "transcriptomic signature"[All Fields] OR '
        '"gene expression signature"[All Fields]) AND '
        '("PD-1"[All Fields] OR "PD-L1"[All Fields] OR "immunotherapy"[All Fields] OR '
        '"checkpoint"[All Fields]) AND "Homo sapiens"[Organism] AND "gse"[Entry Type]'
    ),
    "tumour_biopsy_ici": (
        '("pretreatment"[All Fields] OR "pre-treatment"[All Fields] OR "baseline biopsy"[All Fields] '
        'OR "on-treatment"[All Fields] OR "tumor biopsy"[All Fields] OR "tumour biopsy"[All Fields]) '
        'AND ("PD-1"[All Fields] OR "PD-L1"[All Fields] OR "immune checkpoint"[All Fields]) AND '
        '("lung"[All Fields] OR "NSCLC"[All Fields]) AND "Homo sapiens"[Organism] AND '
        '"gse"[Entry Type]'
    ),
    "pancancer_ici": (
        '("pan-cancer"[All Fields] OR "solid tumors"[All Fields] OR "solid tumours"[All Fields] OR '
        '"multiple tumor types"[All Fields] OR "advanced cancer"[All Fields]) AND '
        '("immune checkpoint"[All Fields] OR "anti-PD-1"[All Fields] OR "anti-PD-L1"[All Fields]) '
        'AND ("response"[All Fields] OR "responder"[All Fields]) AND '
        '"Homo sapiens"[Organism] AND "gse"[Entry Type]'
    ),
}


def main() -> int:
    existing = {r["accession"]: r for r in json.loads((OUT / "geo_search_candidates.json").read_text())}
    hits: dict[str, set[str]] = {}
    for name, term in EXTRA_QUERIES.items():
        uids = s1.esearch(term)
        print(f"[search2] {name}: {len(uids)} uids")
        for uid in uids:
            hits.setdefault(uid, set()).add(name)
        time.sleep(0.5)

    summaries = s1.esummary(sorted(hits))
    added = 0
    for uid, meta in summaries.items():
        acc = meta.get("accession", "")
        if not acc.startswith("GSE") or acc in s1.EXCLUDED:
            continue
        if acc in existing:
            existing[acc]["queries"] += ";" + ";".join(sorted(hits[uid]))
            continue
        existing[acc] = {
            "accession": acc,
            "uid": uid,
            "title": (meta.get("title") or "").strip(),
            "summary": " ".join((meta.get("summary") or "").split()),
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
        added += 1

    rows = sorted(existing.values(), key=lambda r: r["accession"])
    (OUT / "geo_search_candidates.json").write_text(json.dumps(rows, indent=2))
    print(f"[search2] added {added} new GSEs; pool now {len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
