#!/usr/bin/env python3
"""Search NCBI GEO for 2021 human lung ICI series. Real accessions only."""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path(__file__).resolve().parent
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

LUNG = (
    "lung",
    "NSCLC",
    "non-small cell lung",
    "lung adenocarcinoma",
    "LUAD",
    "lung squamous",
    "LUSC",
    "small cell lung",
    "SCLC",
    "pulmonary carcinoma",
)
ICI = (
    "immune checkpoint",
    "checkpoint inhibitor",
    "checkpoint blockade",
    "immunotherapy",
    "PD-1",
    "PD1",
    "PD-L1",
    "PDL1",
    "CTLA-4",
    "CTLA4",
    "nivolumab",
    "pembrolizumab",
    "atezolizumab",
    "durvalumab",
    "avelumab",
    "ipilimumab",
    "cemiplimab",
)
COMMON = '"Homo sapiens"[Organism] AND "gse"[Entry Type]'
DATE = '("2021/01/01"[PDAT] : "2021/12/31"[PDAT])'


def or_fields(terms: tuple[str, ...]) -> str:
    return "(" + " OR ".join(f'"{term}"[All Fields]' for term in terms) + ")"


def esearch(term: str) -> dict:
    params = {"db": "gds", "term": term, "retmax": "2000", "retmode": "json"}
    url = f"{EUTILS}/esearch.fcgi?" + urllib.parse.urlencode(params)
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                return json.load(response).get("esearchresult", {})
        except Exception as error:  # noqa: BLE001
            time.sleep(2**attempt)
            last = error
    raise RuntimeError(last)


def esummary(uids: list[str]) -> list[dict]:
    records = []
    for start in range(0, len(uids), 20):
        batch = uids[start : start + 20]
        params = {"db": "gds", "id": ",".join(batch), "retmode": "json"}
        url = f"{EUTILS}/esummary.fcgi?" + urllib.parse.urlencode(params)
        for attempt in range(5):
            try:
                with urllib.request.urlopen(url, timeout=90) as response:
                    payload = json.load(response)
                break
            except Exception as error:  # noqa: BLE001
                time.sleep(2**attempt)
                last = error
        else:
            raise RuntimeError(last)
        result = payload.get("result", {})
        for uid in result.get("uids", []):
            rec = result[uid]
            records.append(
                {
                    "uid": uid,
                    "accession": rec.get("accession"),
                    "title": rec.get("title"),
                    "summary": rec.get("summary"),
                    "taxon": rec.get("taxon"),
                    "n_samples": rec.get("n_samples"),
                    "gdsType": rec.get("gdsType"),
                    "pdat": rec.get("pdat"),
                    "suppFile": rec.get("suppFile"),
                    "pubmedids": rec.get("pubmedids"),
                    "ftplink": rec.get("ftplink"),
                }
            )
        print(f"fetched {min(start + 20, len(uids))}/{len(uids)}")
        time.sleep(0.4)
    records.sort(key=lambda item: item.get("accession") or "")
    return records


def main() -> None:
    term = f"{or_fields(LUNG)} AND {or_fields(ICI)} AND {COMMON} AND {DATE}"
    result = esearch(term)
    uids = set(result.get("idlist", []))
    print(f"combined count={result.get('count')} fetched={len(uids)}")
    for drug in (
        "nivolumab",
        "pembrolizumab",
        "atezolizumab",
        "durvalumab",
        "ipilimumab",
        "avelumab",
        "cemiplimab",
    ):
        extra = esearch(f'{or_fields(LUNG)} AND "{drug}"[All Fields] AND {COMMON} AND {DATE}')
        new = set(extra.get("idlist", []))
        print(f"  {drug}: count={extra.get('count')} new={len(new - uids)}")
        uids |= new
        time.sleep(0.4)

    records = esummary(sorted(uids))
    (OUT / "search_uids_2021.json").write_text(
        json.dumps({"query": term, "uids": sorted(uids)}, indent=2)
    )
    (OUT / "candidates_2021.json").write_text(json.dumps(records, indent=2))
    print(f"Wrote {len(records)} 2021 records")


if __name__ == "__main__":
    main()
