#!/usr/bin/env python3
"""ArrayExpress leftover hunt: lung ICI / TACSTD2 / CLDN4 not already covered.

The first-wave hunt (PR #3, scripts/gpt_arrayexpress) already reviewed a
hard-coded accession set. This script re-queries BioStudies ArrayExpress,
excludes that set plus GEO mirrors, and writes an honest leftover catalog.
E-MTAB-13530 is skipped by request (Visium NSCLC; parallel leftover slice).
"""

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

# First-wave ArrayExpress hunt decisions (PR #3). Do not re-claim these.
COVERED = {
    "E-MTAB-13704": ("include", "GEMM in-situ lung aPD-L1 combinations; processed counts"),
    "E-MTAB-15883": ("include", "lung-cancer-model anti-PD-1; processed data"),
    "E-MTAB-9451": ("include_context", "NSCLC ICOS/Treg profiling motivated by checkpoint response; no ICI-treated samples"),
    "E-MTAB-10633": ("metadata_only", "whole-lung aPD-L1/TGF-beta-trap; no processed files"),
    "E-MTAB-8867": ("metadata_only", "ICI-myocarditis cohort with 3 NSCLC patients; no processed files"),
    "E-MTAB-13708": ("exclude", "murine lung cancer vector immunotherapy, not ICI"),
    "E-MTAB-10027": ("exclude", "breast-cancer model; lung is metastatic site"),
    "E-MTAB-13770": ("exclude", "melanoma ICI; lung is incidental text"),
    "E-MTAB-3218": ("exclude", "RCC nivolumab biopsies"),
    "E-MTAB-3732": ("exclude", "generic expression atlas"),
    "E-MTAB-16855": ("exclude", "lung epithelial antiviral study; not cancer ICI"),
}

# User: skip E-MTAB-13530 if already hunted. Companion Visium NSCLC atlas;
# not ICI-labeled; parallel B6 Visium leftover slice owns it.
SKIP = {
    "E-MTAB-13530": "skip_requested: Visium NSCLC atlas, not ICI-labeled; already assigned / hunted elsewhere",
}

# Manual leftover decisions applied after API discovery. Never invent accessions.
LEFTOVER_DECISIONS = {
    "E-MTAB-15784": (
        "leftover_analyze",
        "NSCLC patient-derived organoid + autologous TIL co-culture scRNA; processed MTX <2GB; TIL therapy, not ICI",
    ),
    "E-MTAB-15782": (
        "leftover_duplicate",
        "same study as E-MTAB-15784; processed matrix is identical in size to LCP89_CO (co-culture only)",
    ),
    "E-MTAB-12508": (
        "leftover_metadata",
        "KP lung + FLT3L/aCD40 DC therapy; paper title says refractory to checkpoint inhibitors; no processed files",
    ),
    "E-MTAB-13710": (
        "leftover_metadata",
        "human NSCLC FRC/TLS companion to excluded E-MTAB-13708; no ICI treatment; no processed files",
    ),
    "E-MTAB-16381": (
        "leftover_metadata",
        "LLC-OVA tumor-infiltrating NK cells + mifepristone; not ICI; no processed files",
    ),
    "E-MTAB-13526": (
        "leftover_too_large",
        "companion scRNA atlas to skipped E-MTAB-13530; not ICI; annotated h5ads 45–58 GB exceed <2GB rule",
    ),
    "E-MTAB-14799": (
        "exclude",
        "TILT-123 oncolytic adenovirus PBMCs from mixed solid tumors; not lung ICI",
    ),
    "E-MTAB-11377": ("exclude", "Trop2+ mouse intestine; not lung"),
    "E-MTAB-11382": ("exclude", "Trop2+ mouse intestinal adenoma; not lung"),
    "E-MTAB-11466": ("exclude", "TROP2+ human CRC; not lung"),
    "E-MTAB-16433": ("exclude", "TROP2 ADC in CRC PDOX; not lung"),
    "E-MTAB-16843": ("exclude", "sacituzumab govitecan CRC organoids; not lung"),
    "E-MTAB-16849": ("exclude", "sacituzumab govitecan CRC liver mets; not lung"),
    "E-MTAB-16835": ("exclude", "murine CRC; not lung"),
    "E-MTAB-16836": ("exclude", "CRC organoid chemo; not lung"),
    "E-MTAB-16583": ("exclude", "CRC primary/liver met scRNA; not lung"),
    "E-MTAB-16585": ("exclude", "CRC PDOX scRNA; not lung"),
    "E-MTAB-14030": ("exclude", "allergy immunotherapy (AIT) lung CD4 T cells; not cancer ICI"),
    "E-MTAB-6044": ("exclude", "IL-22 immunotherapy in bacterial lung infection; not cancer ICI"),
    "E-MTAB-8294": ("exclude", "AXL vs NK/CTL cytotoxicity in lung lines; raw arrays only; not ICI"),
    "E-MTAB-3266": ("exclude", "CD8/CD103 TIL prognostic arrays in lung; raw only; not ICI"),
    "E-MTAB-3448": ("exclude", "early-stage NSCLC bulk array; treatment-naive; not ICI leftover"),
    "E-MTAB-3665": ("exclude", "all-stage NSCLC bulk array; treatment-naive; not ICI leftover"),
    "E-MTAB-6043": ("exclude", "NSCLC microarray meta-dataset; treatment-naive; not ICI leftover"),
    "E-MTAB-10745": ("exclude", "PSME4 / immunosuppressive TME; not a lung ICI cohort"),
    "E-MTAB-11734": ("exclude", "TLR agonist tumor digests; not lung ICI"),
    "E-MTAB-10399": ("exclude", "combined LUAD/SCLC RNAseq; not ICI"),
    "E-MTAB-13968": ("exclude", "lung PDX scRNA with minimal TME; not ICI"),
    "E-MTAB-13962": ("exclude", "PC9 + NIH-3T3 co-culture; not ICI"),
}

QUERIES = [
    'lung AND immunotherapy',
    'lung AND "checkpoint"',
    'lung AND "checkpoint inhibitor"',
    'lung AND "checkpoint blockade"',
    'lung AND "immune checkpoint"',
    'NSCLC AND immunotherapy',
    'NSCLC AND "PD-1"',
    'NSCLC AND "PD-L1"',
    'lung AND nivolumab',
    'lung AND pembrolizumab',
    'lung AND atezolizumab',
    'lung AND durvalumab',
    'lung AND "anti-PD-1"',
    'lung AND "anti-PD-L1"',
    'lung AND PD1',
    'lung AND PDL1',
    'lung AND TACSTD2',
    'lung AND TROP2',
    'lung AND CLDN4',
    '"claudin-4" AND lung',
    '"claudin 4"',
    'TACSTD2',
    'TROP2',
    'CLDN4',
    '"small cell lung" AND immunotherapy',
    'LLC AND immunotherapy',
    'E-MTAB-13526',
    'E-MTAB-13530',
    'E-MTAB-15784',
    'E-MTAB-15782',
    'E-MTAB-12508',
    'E-MTAB-13710',
]


def get_json(url: str, retries: int = 4) -> dict:
    request = urllib.request.Request(
        url, headers={"Accept": "application/json", "User-Agent": "ae-leftover/1.0"}
    )
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                return json.load(response)
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep(2**attempt)
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
    if isinstance(node, dict):
        return " ".join(flatten_values(v) for v in node.values())
    if isinstance(node, list):
        return " ".join(flatten_values(v) for v in node)
    return "" if node is None else str(node)


def extract_files(study: dict) -> list[dict]:
    files = []
    for node, ancestors in walk(study):
        if node.get("type") != "file" or "path" not in node:
            continue
        attrs = {a.get("name"): a.get("value") for a in node.get("attributes", [])}
        context = " / ".join(ancestors)
        description = str(attrs.get("Description", ""))
        processed = "processed" in (context + " " + description).lower()
        files.append(
            {
                "path": node["path"],
                "size": int(node.get("size") or 0),
                "description": description,
                "processed": processed,
                "download_url": (
                    f"https://www.ebi.ac.uk/biostudies/files/"
                    f"{study['accno']}/{urllib.parse.quote(node['path'])}"
                ),
            }
        )
    return files


def decide(accession: str, text: str, geo_matches: list[str]) -> tuple[str, str]:
    if accession in SKIP:
        return "skip_requested", SKIP[accession]
    if accession in COVERED:
        prior, why = COVERED[accession]
        return f"already_covered_{prior}", why
    if accession.upper().startswith("E-GEOD-") or geo_matches:
        return "exclude_geo_mirror", "GEO-mirrored series; not an ArrayExpress leftover"
    if accession in LEFTOVER_DECISIONS:
        return LEFTOVER_DECISIONS[accession]
    return (
        "exclude",
        "query hit failed leftover lung ICI / TACSTD2 / CLDN4 relevance review",
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/w200/AE_leftover/search_manifest.json"),
    )
    args = parser.parse_args()

    found: dict[str, dict] = {}
    query_hits: dict[str, list[str]] = {}
    hit_by_accession: dict[str, list[str]] = defaultdict(list)
    for query in QUERIES:
        hits = search(query)
        query_hits[query] = [h["accession"] for h in hits]
        for hit in hits:
            acc = hit["accession"]
            found[acc] = hit
            hit_by_accession[acc].append(query)
        print(f"{len(hits):4d}  {query}", flush=True)

    studies = []
    for accession in sorted(found):
        study = get_json(f"{API}/studies/{urllib.parse.quote(accession)}")
        text = flatten_values(study)
        files = extract_files(study)
        geo_matches = sorted(set(GEO_RE.findall(text)))
        decision, rationale = decide(accession, text, geo_matches)
        studies.append(
            {
                "accession": accession,
                "title": found[accession].get("title"),
                "queries": hit_by_accession[accession],
                "decision": decision,
                "rationale": rationale,
                "geo_matches": geo_matches,
                "study_url": f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}",
                "n_files": len(files),
                "n_processed": sum(1 for f in files if f["processed"]),
                "processed_bytes": sum(f["size"] for f in files if f["processed"]),
                "processed_files": [
                    {k: f[k] for k in ("path", "size", "description", "download_url")}
                    for f in files
                    if f["processed"]
                ],
                "metadata_files": [
                    {k: f[k] for k in ("path", "size", "description", "download_url")}
                    for f in files
                    if f["path"].endswith((".idf.txt", ".sdrf.txt"))
                ],
            }
        )

    counts = defaultdict(int)
    for s in studies:
        counts[s["decision"]] += 1

    payload = {
        "searched_at_utc": datetime.now(timezone.utc).isoformat(),
        "api": f"{API}/arrayexpress/search",
        "queries": query_hits,
        "study_count": len(studies),
        "decision_counts": dict(counts),
        "covered_first_wave": COVERED,
        "skipped": SKIP,
        "honest_bottom_line": (
            "No additional ArrayExpress lung ICI cohort with processed expression "
            "and ICI response labels was found beyond the first-wave set. "
            "The only leftover with open processed matrices under 2 GB that is "
            "lung + immune is E-MTAB-15784 (TIL–organoid co-culture, not ICI)."
        ),
        "studies": studies,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote {args.out}: {len(studies)} studies; counts={dict(counts)}")


if __name__ == "__main__":
    main()
