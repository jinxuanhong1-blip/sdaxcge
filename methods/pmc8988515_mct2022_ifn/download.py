#!/usr/bin/env python3
"""Download the expression tables scored for PMC8988515 and refresh the repository hunt."""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
UA = "Mozilla/5.0 (compatible; pmc8988515-ifn/1.0)"

FILES = {
    "table_s1.xlsx": "https://ndownloader.figshare.com/files/39984958",
    "table_s2.xlsx": "https://ndownloader.figshare.com/files/39984961",
    "table_s3.xlsx": "https://ndownloader.figshare.com/files/39984964",
    "table_s4.xlsx": "https://ndownloader.figshare.com/files/39984967",
    "tcga_ov_hiseqv2.gz": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.OV.sampleMap/HiSeqV2.gz",
    "coscia_moesm1558.xlsx": "https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fncomms12645/MediaObjects/41467_2016_BFncomms12645_MOESM1558_ESM.xlsx",
}


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"keep {dest.name} ({dest.stat().st_size} bytes)")
        return
    print(f"get {url}")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    print(f"  wrote {dest.name} {dest.stat().st_size}")


def get_json(url: str, data: bytes | None = None, headers: dict | None = None):
    h = {"User-Agent": UA, "Accept": "application/json"}
    if headers:
        h.update(headers)
    req = urllib.request.Request(url, data=data, headers=h)
    with urllib.request.urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode())


def geo_count(term: str) -> int:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(
        {"db": "gds", "retmax": "0", "retmode": "json", "term": term}
    )
    d = get_json(url)
    return int(d["esearchresult"]["count"])


def biostudies_hits(query: str) -> tuple[int, str]:
    url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(
        {"query": query, "pageSize": "5"}
    )
    d = get_json(url)
    titles = "; ".join(
        f"{h.get('accession')} {(h.get('title') or '')[:80]}" for h in (d.get("hits") or [])[:3]
    )
    return int(d.get("totalHits") or 0), titles


def pride_n(keyword: str) -> int:
    url = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?" + urllib.parse.urlencode(
        {"keyword": keyword, "pageSize": "5"}
    )
    d = get_json(url)
    if isinstance(d, list):
        return len(d)
    return int(d.get("totalElements") or d.get("page", {}).get("totalElements") or 0)


def node_n(keyword: str) -> int:
    body = json.dumps(
        {
            "queryWord": keyword,
            "advQueryWord": "",
            "sortKey": "modifiedDate",
            "sortType": "desc",
            "pageNum": 1,
            "pageSize": 5,
        }
    ).encode()
    d = get_json(
        "https://www.biosino.org/node/api/app/browse/search",
        data=body,
        headers={"Content-Type": "application/json"},
    )
    page = (d.get("data") or {}).get("pageInfo") or {}
    return int(page.get("totalElements") or 0)


def write_inventory() -> None:
    rows = []

    def add(**kwargs):
        rows.append(kwargs)

    geo_queries = [
        ("Bitler[Author] AND (CLDN4 OR claudin-4)", "author + CLDN4"),
        ("Yamamoto[Author] AND (CLDN4 OR claudin-4) AND ovarian", "author + CLDN4 + ovarian"),
        ('PMC8988515 OR MCT-21-0827 OR "Loss of Claudin-4 Reduces DNA Damage"', "paper title / PMC"),
    ]
    for term, role in geo_queries:
        n = geo_count(term)
        add(
            repository="GEO",
            query=term,
            n_hits=n,
            accession="",
            role=role,
            download_status="no series",
            scored="no",
            note="esearch db=gds",
        )

    ae_queries = [
        (
            "collection:arrayexpress AND (Bitler OR Yamamoto) AND (CLDN4 OR claudin)",
            "author + CLDN4 in ArrayExpress",
        ),
        (
            "accession:E-MTAB-8107 OR accession:E-MTAB-6149 OR accession:E-MTAB-6653",
            "Qian 2020 blueprint cited as ref 17",
        ),
    ]
    for term, role in ae_queries:
        n, titles = biostudies_hits(term)
        add(
            repository="ArrayExpress",
            query=term,
            n_hits=n,
            accession="E-MTAB-8107;E-MTAB-6149;E-MTAB-6653" if "E-MTAB-8107" in term else "",
            role=role,
            download_status="listed via BioStudies API; ftp.ebi.ac.uk transfer failed in this environment"
            if "E-MTAB-8107" in term
            else "no ArrayExpress study",
            scored="no",
            note=titles,
        )

    for kw in ["Bitler claudin", "Yamamoto OVCAR3 CLDN4", "PMC8988515"]:
        n = pride_n(kw)
        add(
            repository="PRIDE",
            query=kw,
            n_hits=n,
            accession="",
            role="author / paper keyword",
            download_status="no project",
            scored="no",
            note="pride ws archive v2 search/projects",
        )
    add(
        repository="PRIDE",
        query="PXD003668 (cited Coscia 2016, ref 22; not an author deposit)",
        n_hits=1,
        accession="PXD003668",
        role="cited ovarian cell-line proteome",
        download_status="PRIDE files are RAW; scored the Nature processed LFQ table (MOESM1558)",
        scored="yes",
        note="30 columns of log2 MaxLFQ; OvCa subset drops ME180, ME180C13, IOSE397, IOSE7576",
    )

    for kw in ["CLDN4", "claudin-4", "Bitler", "PMC8988515", "MCT-21-0827", "OVCAR3 CLDN4"]:
        n = node_n(kw)
        add(
            repository="NODE",
            query=kw,
            n_hits=n,
            accession="",
            role="public browse full text",
            download_status="no record",
            scored="no",
            note="POST /node/api/app/browse/search; lung cancer control query returns hits",
        )

    add(
        repository="figshare (AACR supplement, not GEO/AE/PRIDE/NODE)",
        query="10.1158/1535-7163.22522216",
        n_hits=1,
        accession="mct-21-0827 table S1",
        role="author transcriptome table of TCGA OV CLDN4 high vs low",
        download_status="downloaded",
        scored="yes",
        note="full gene table with log2(high/low), p, q. S2 BRCA status, S3 drug viability, S4 ex vivo histology. No RPPA matrix.",
    )
    add(
        repository="UCSC Xena (TCGA OV HiSeqV2; not GEO/AE/PRIDE/NODE)",
        query="TCGA.OV.sampleMap/HiSeqV2",
        n_hits=1,
        accession="TCGA-OV Firehose RSEM log2(norm_count+1)",
        role="public matrix behind the paper Fig 1 median split",
        download_status="downloaded",
        scored="yes",
        note="paper text Low n=154 High n=153 on cBioPortal Firehose Legacy 2020-12-20; this file is the Xena freeze",
    )

    TABLES.mkdir(parents=True, exist_ok=True)
    cols = ["repository", "query", "n_hits", "accession", "role", "download_status", "scored", "note"]
    lines = ["\t".join(cols)]
    for row in rows:
        lines.append("\t".join(str(row.get(c, "")).replace("\t", " ") for c in cols))
    out = TABLES / "inventory.tsv"
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out} ({len(rows)} rows)")


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for name, url in FILES.items():
        fetch(url, DATA / name)
    write_inventory()


if __name__ == "__main__":
    main()
