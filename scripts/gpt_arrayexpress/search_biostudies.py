#!/usr/bin/env python3
"""Reproduce the lung ICI / TACSTD2 / CLDN4 ArrayExpress search."""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

API = "https://www.ebi.ac.uk/biostudies/api/v1"
PAGE_SIZE = 100
GEO_RE = re.compile(r"\b(?:GSE|E-GEOD-)\d+\b", re.I)

LUNG_TERMS = ['lung', '"lung cancer"', '"non-small cell lung"', "NSCLC", "LLC"]
ICI_TERMS = [
    '"anti-PD-1"', "anti-PD1", '"anti-PD-L1"', "anti-PDL1",
    "aPD1", "aPDL1", '"PD-1"', "PD1", '"PD-L1"', "PDL1", "ICI",
    "checkpoint", '"checkpoint inhibitor"', '"immune checkpoint inhibitors"',
    '"checkpoint blockade"', '"immune checkpoint blockade"', "immunotherapy", "nivolumab",
    "pembrolizumab", "atezolizumab", "durvalumab", "cemiplimab",
    "ipilimumab",
]
TARGET_TERMS = ["TACSTD2", "TROP2", '"TROP-2"', "CLDN4", '"claudin 4"', '"claudin-4"']

# These decisions are intentionally explicit and auditable. They are applied only
# after API discovery; the script never manufactures an accession.
DECISIONS = {
    "E-MTAB-13704": ("include", "in-situ lung GEMM anti-PD-L1 combination experiment; processed counts"),
    "E-MTAB-15883": ("include", "direct lung-cancer-model anti-PD-1 experiment; processed data"),
    "E-MTAB-9451": ("include_context", "NSCLC immune profiling motivated by checkpoint response; no ICI-treated samples"),
    "E-MTAB-10633": ("metadata_only", "whole-lung anti-PD-L1/TGF-beta-trap experiment; no processed files"),
    "E-MTAB-8867": ("metadata_only", "three NSCLC patients in ICI-myocarditis cohort; no processed files"),
    "E-MTAB-13708": ("exclude", "lung cancer immunotherapy, but vector immunotherapy rather than ICI"),
    "E-MTAB-10027": ("exclude", "breast-cancer model; lung appears only as metastatic site"),
    "E-MTAB-13770": ("exclude", "melanoma cohort; lung is incidental text"),
    "E-MTAB-3218": ("exclude", "renal-cell carcinoma cohort"),
    "E-MTAB-3732": ("exclude", "generic expression atlas; no lung ICI experiment"),
    "E-MTAB-16855": ("exclude", "lung epithelial antiviral study; not cancer ICI"),
}


def get_json(url: str, retries: int = 4) -> dict:
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "gpt-arrayexpress/1.0"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def search(query: str) -> list[dict]:
    hits, page = [], 1
    while True:
        params = urllib.parse.urlencode({"query": query, "pageSize": PAGE_SIZE, "page": page})
        payload = get_json(f"{API}/arrayexpress/search?{params}")
        batch = payload.get("hits", [])
        hits.extend(batch)
        if not batch or len(hits) >= int(payload.get("totalHits", len(hits))):
            return hits
        page += 1


def walk(node, ancestors=()):
    if isinstance(node, dict):
        yield node, ancestors
        label = str(node.get("type") or node.get("accno") or "")
        for value in node.values():
            yield from walk(value, ancestors + (label,))
    elif isinstance(node, list):
        for value in node:
            yield from walk(value, ancestors)


def flatten_values(node) -> str:
    values = []
    if isinstance(node, dict):
        for value in node.values():
            values.append(flatten_values(value))
    elif isinstance(node, list):
        for value in node:
            values.append(flatten_values(value))
    elif node is not None:
        values.append(str(node))
    return " ".join(values)


def extract_files(study: dict) -> list[dict]:
    files = []
    for node, ancestors in walk(study):
        if node.get("type") != "file" or "path" not in node:
            continue
        attrs = {a.get("name"): a.get("value") for a in node.get("attributes", [])}
        context = " / ".join(ancestors)
        description = str(attrs.get("Description", ""))
        processed = "processed" in (context + " " + description).lower()
        files.append({
            "path": node["path"],
            "size": int(node.get("size") or 0),
            "description": description,
            "processed": processed,
            "download_url": f"https://www.ebi.ac.uk/biostudies/files/{study['accno']}/{urllib.parse.quote(node['path'])}",
        })
    return files


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=Path("results/gpt_arrayexpress/search_manifest.json"))
    args = parser.parse_args()

    queries = []
    for lung in LUNG_TERMS:
        queries.extend(f"{term} AND {lung}" for term in ICI_TERMS)
        queries.extend(f"{term} AND {lung}" for term in TARGET_TERMS)
    # Exact target-only searches guard against records whose lung context is
    # tokenized unexpectedly; full metadata is inspected below.
    queries.extend(TARGET_TERMS)
    queries = list(dict.fromkeys(queries))

    found: dict[str, dict] = {}
    query_hits: dict[str, list[str]] = {}
    hit_by_accession = defaultdict(list)
    for query in queries:
        hits = search(query)
        query_hits[query] = [hit["accession"] for hit in hits]
        for hit in hits:
            accession = hit["accession"]
            found[accession] = hit
            hit_by_accession[accession].append(query)

    studies = []
    for accession in sorted(found):
        study = get_json(f"{API}/studies/{urllib.parse.quote(accession)}")
        text = flatten_values(study)
        files = extract_files(study)
        geo_matches = sorted(set(GEO_RE.findall(text)))
        decision, rationale = DECISIONS.get(accession, ("exclude", "query hit failed manual disease/intervention relevance review"))
        if accession.upper().startswith("E-GEOD-") or geo_matches:
            decision, rationale = "exclude_geo_mirror", "GEO-mirrored series"
        studies.append({
            "accession": accession,
            "title": found[accession].get("title"),
            "queries": hit_by_accession[accession],
            "decision": decision,
            "rationale": rationale,
            "geo_matches": geo_matches,
            "study_url": f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}",
            "processed_files": [f for f in files if f["processed"]],
            "all_files": files,
        })

    payload = {
        "searched_at_utc": datetime.now(timezone.utc).isoformat(),
        "api": f"{API}/arrayexpress/search",
        "page_size": PAGE_SIZE,
        "queries": query_hits,
        "study_count": len(studies),
        "studies": studies,
        "notes": [
            "ArrayExpress is hosted in BioStudies; the collection endpoint was searched with complete pagination.",
            "E-GEOD accessions and any record containing a GSE/E-GEOD identifier were excluded as GEO mirrors.",
            "The repository contained no pre-existing GEO accession list at search time.",
            "The <2GB rule is applied per open processed file by download_biostudies.py.",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote {args.out}: {len(studies)} unique studies from {len(queries)} queries")


if __name__ == "__main__":
    main()
