#!/usr/bin/env python3
"""Search GEO for human lung ICI blood/PBMC/plasma series. Do not drop blood."""
from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

OUT = Path("results/noskip/blood_all")
OUT.mkdir(parents=True, exist_ok=True)

QUERIES = {
    "lung_ici_blood": (
        '(NSCLC[Title] OR "non-small cell lung"[Title] OR "lung cancer"[Title] '
        'OR "lung adenocarcinoma"[Title] OR "small cell lung"[Title] '
        'OR LUAD[Title] OR LUSC[Title]) '
        'AND (PD-1 OR PD-L1 OR CTLA-4 OR nivolumab OR pembrolizumab OR atezolizumab '
        'OR durvalumab OR avelumab OR ipilimumab OR cemiplimab OR tislelizumab '
        'OR camrelizumab OR sintilimab OR toripalimab OR "immune checkpoint" '
        'OR immunotherapy OR ICI OR ICB) '
        'AND (blood OR PBMC OR PBMCs OR plasma OR serum OR "whole blood" '
        'OR PAXgene OR platelet OR TEP OR "peripheral blood" OR circulating) '
        'AND "Homo sapiens"[Organism] AND gse[Entry Type]'
    ),
    "lung_ici_pbmc": (
        '(NSCLC OR "lung cancer" OR "non-small cell lung") '
        'AND (nivolumab OR pembrolizumab OR atezolizumab OR durvalumab '
        'OR "anti-PD-1" OR "anti-PD-L1" OR "PD-1 blockade") '
        'AND (PBMC OR PBMCs OR "peripheral blood mononuclear") '
        'AND "Homo sapiens"[Organism] AND gse[Entry Type]'
    ),
    "named_accessions": (
        "(GSE285888[ACCN] OR GSE305086[ACCN] OR GSE111414[ACCN] OR GSE152590[ACCN] "
        "OR GSE216297[ACCN] OR GSE235048[ACCN] OR GSE266219[ACCN] OR GSE225620[ACCN] "
        "OR GSE295969[ACCN] OR GSE306542[ACCN] OR GSE310370[ACCN] OR GSE207715[ACCN] "
        "OR GSE237087[ACCN] OR GSE202417[ACCN] OR GSE242780[ACCN] OR GSE190905[ACCN] "
        "OR GSE247754[ACCN] OR GSE136961[ACCN] OR GSE161537[ACCN] OR GSE274588[ACCN] "
        "OR GSE315510[ACCN] OR GSE278119[ACCN] OR GSE244416[ACCN] OR GSE93157[ACCN] "
        "OR GSE115821[ACCN] OR GSE145996[ACCN] OR GSE179351[ACCN] OR GSE183996[ACCN] "
        "OR GSE186143[ACCN] OR GSE201322[ACCN] OR GSE99254[ACCN] OR GSE139327[ACCN] "
        "OR GSE165383[ACCN] OR GSE295601[ACCN] OR GSE308745[ACCN]) AND gse[Entry Type]"
    ),
}

ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"


def get_json(url: str, retries: int = 5) -> dict:
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "noskip-blood-ici/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"failed {url}: {last}")


def esearch(term: str) -> list[str]:
    q = urllib.parse.urlencode(
        {"db": "gds", "term": term, "retmax": 400, "retmode": "json", "sort": "relevance"}
    )
    data = get_json(f"{ESEARCH}?{q}")
    return data.get("esearchresult", {}).get("idlist", [])


def esummary(ids: list[str]) -> dict:
    out = {}
    for i in range(0, len(ids), 40):
        chunk = ids[i : i + 40]
        q = urllib.parse.urlencode(
            {"db": "gds", "id": ",".join(chunk), "retmode": "json"}
        )
        data = get_json(f"{ESUMMARY}?{q}")
        out.update(data.get("result", {}))
        time.sleep(0.35)
    return out


def main() -> None:
    uid_to_queries: dict[str, list[str]] = {}
    for name, term in QUERIES.items():
        print(f"search {name} ...", flush=True)
        ids = esearch(term)
        print(f"  {len(ids)} uids", flush=True)
        for uid in ids:
            uid_to_queries.setdefault(uid, []).append(name)
        time.sleep(0.35)

    uids = sorted(uid_to_queries)
    print(f"unique uids: {len(uids)}", flush=True)
    summaries = esummary(uids)

    rows = []
    for uid in uids:
        rec = summaries.get(uid) or {}
        if not isinstance(rec, dict):
            continue
        acc = rec.get("accession") or ""
        if not acc.startswith("GSE"):
            # GDS records sometimes wrap GSE
            ext = rec.get("extrelations") or []
            for rel in ext:
                tacc = (rel or {}).get("targetacc") or ""
                if tacc.startswith("GSE"):
                    acc = tacc
                    break
        rows.append(
            {
                "uid": uid,
                "accession": acc,
                "title": rec.get("title") or "",
                "summary": rec.get("summary") or "",
                "taxon": rec.get("taxon") or "",
                "n_samples": rec.get("n_samples"),
                "gdstype": rec.get("gdstype") or "",
                "pdat": rec.get("pdat") or "",
                "ftp": rec.get("ftplink") or "",
                "suppfile": rec.get("suppfile") or "",
                "queries": ";".join(uid_to_queries.get(uid, [])),
            }
        )

    (OUT / "geo_search_raw.json").write_text(json.dumps({"queries": QUERIES, "rows": rows}, indent=2))
    # TSV
    keys = ["accession", "title", "taxon", "n_samples", "gdstype", "pdat", "suppfile", "ftp", "queries", "summary"]
    lines = ["\t".join(keys)]
    seen = set()
    for r in rows:
        acc = r["accession"]
        if not acc or acc in seen:
            continue
        seen.add(acc)
        lines.append("\t".join(str(r.get(k, "")).replace("\t", " ").replace("\n", " ") for k in keys))
    (OUT / "geo_search_candidates.tsv").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(seen)} unique GSE to {OUT / 'geo_search_candidates.tsv'}")


if __name__ == "__main__":
    main()
