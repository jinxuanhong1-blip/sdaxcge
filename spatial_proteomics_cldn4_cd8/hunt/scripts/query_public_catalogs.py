#!/usr/bin/env python3
"""Re-query public catalogs for CLDN4 protein + CD8 multiplex imaging.

Does not download EGA/dbGaP. Writes JSON under hunt/api_snapshots/.
"""
from __future__ import annotations

import json
import ssl
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parents[1] / "api_snapshots"
OUT.mkdir(parents=True, exist_ok=True)
CTX = ssl.create_default_context()


def get_json(url: str, timeout: int = 60) -> dict:
    req = urllib.request.Request(url, headers={"User-Agent": "CLDN4-CD8-spatial-hunt/1.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=CTX) as r:
        return json.loads(r.read().decode())


def geo_esearch(term: str) -> dict:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "term": term, "retmax": 50, "retmode": "json"}
    )
    return get_json(url)


def geo_esummary(ids: list[str]) -> dict:
    if not ids:
        return {}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "id": ",".join(ids), "retmode": "json"}
    )
    return get_json(url)


def zenodo_search(q: str, size: int = 15) -> dict:
    url = "https://zenodo.org/api/records?" + urllib.parse.urlencode({"q": q, "size": size})
    return get_json(url)


def epmc_search(q: str, page_size: int = 15) -> dict:
    url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode(
        {"query": q, "format": "json", "pageSize": page_size, "resultType": "lite"}
    )
    return get_json(url)


def figshare_search(q: str) -> list:
    url = "https://api.figshare.com/v2/articles/search"
    payload = json.dumps({"search_for": q, "item_type": 3, "page_size": 20}).encode()
    req = urllib.request.Request(
        url,
        data=payload,
        headers={"Content-Type": "application/json", "User-Agent": "CLDN4-CD8-spatial-hunt/1.0"},
    )
    with urllib.request.urlopen(req, timeout=60, context=CTX) as r:
        return json.loads(r.read().decode())


def main() -> None:
    geo_terms = [
        'CLDN4[All Fields] AND (IMC OR "imaging mass cytometry" OR CODEX OR PhenoCycler OR MIBI OR CosMx OR Xenium) AND (gse[Entry Type])',
        '"Claudin-4"[All Fields] AND (IMC OR CODEX OR PhenoCycler OR MIBI OR CosMx) AND (gse[Entry Type])',
    ]
    geo = {}
    ids: set[str] = set()
    for t in geo_terms:
        try:
            geo[t] = geo_esearch(t)
            ids.update(geo[t].get("esearchresult", {}).get("idlist", []))
        except Exception as e:
            geo[t] = {"error": str(e)}
    (OUT / "geo_esearch.json").write_text(json.dumps(geo, indent=2))
    try:
        (OUT / "geo_esummary.json").write_text(json.dumps(geo_esummary(sorted(ids)), indent=2))
    except Exception as e:
        (OUT / "geo_esummary.json").write_text(json.dumps({"error": str(e)}, indent=2))

    zenodo_qs = [
        'metadata.title:"CLDN4"',
        'metadata.description:"Claudin-4" AND metadata.description:CODEX',
        'metadata.description:"Claudin-4" AND (IMC OR "imaging mass")',
        '"Claudin-4" AND "antibody panel"',
    ]
    zen = {}
    for q in zenodo_qs:
        try:
            zen[q] = zenodo_search(q)
        except Exception as e:
            zen[q] = {"error": str(e)}
    (OUT / "zenodo_quoted.json").write_text(json.dumps(zen, indent=2))

    epmc_qs = [
        '"Claudin-4" AND PhenoCycler',
        '"Claudin 4" AND "imaging mass cytometry"',
        '"Claudin-4" AND "multiplex immunofluorescence" AND CD8',
    ]
    epmc = {}
    for q in epmc_qs:
        try:
            epmc[q] = epmc_search(q)
        except Exception as e:
            epmc[q] = {"error": str(e)}
    (OUT / "europepmc.json").write_text(json.dumps(epmc, indent=2))

    try:
        fig = figshare_search("CLDN4 CODEX OR IMC OR multiplex Claudin-4")
        (OUT / "figshare.json").write_text(json.dumps(fig, indent=2))
    except Exception as e:
        (OUT / "figshare.json").write_text(json.dumps({"error": str(e)}, indent=2))

    print(f"Wrote snapshots to {OUT}")
    print(f"GEO unique GDS ids: {len(ids)}")


if __name__ == "__main__":
    main()
