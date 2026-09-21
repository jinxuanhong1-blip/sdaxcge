#!/usr/bin/env python3
"""Inventory NODE and CNCB (GSA, GSA-Human, OMIX, BioProject) for
CLDN4 knockdown datasets in lung, ovarian, or breast cancer.

Writes query counts and hit tables. Does not download sequence files.
Controlled-access FASTQ is recorded and left in place.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "results" / "tables"
TABLES.mkdir(parents=True, exist_ok=True)

UA = {"User-Agent": "Mozilla/5.0", "Accept": "application/json"}


def get_json(url: str, data: dict | None = None, timeout: int = 60) -> dict:
    body = None
    headers = dict(UA)
    if data is not None:
        body = json.dumps(data).encode()
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode())


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    lines = ["\t".join(columns)]
    for row in rows:
        vals = []
        for col in columns:
            text = str(row.get(col, "") if row.get(col) is not None else "")
            text = text.replace("\t", " ").replace("\n", " ").replace("\r", " ")
            vals.append(text)
        lines.append("\t".join(vals))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def node_search(query: str, page_size: int = 20) -> tuple[int, list[dict]]:
    url = "https://www.biosino.org/node/api/app/browse/search"
    got: list[dict] = []
    total = 0
    page = 1
    while True:
        payload = {
            "queryWord": query,
            "advQueryWord": "",
            "sortKey": "modifiedDate",
            "sortType": "desc",
            "pageNum": page,
            "pageSize": page_size,
        }
        data = get_json(url, payload)
        info = data["data"]["pageInfo"]
        total = int(info.get("totalElements") or 0)
        batch = info.get("content") or []
        got.extend(batch)
        if not batch or len(got) >= total:
            break
        page += 1
        if page > 20:
            break
    return total, got


def cncb_search(db: str, query: str, page_size: int = 50, max_pages: int = 8) -> tuple[int, list[dict]]:
    base = "https://ngdc.cncb.ac.cn/search/api/specific"
    got: list[dict] = []
    total = 0
    for start in range(0, page_size * max_pages, page_size):
        url = base + "?" + urllib.parse.urlencode(
            {"db": db, "q": query, "start": start, "length": page_size}
        )
        data = get_json(url)
        block = data["result"]["data"]
        total = int(block.get("recordsTotal") or 0)
        batch = block.get("data") or []
        if not batch:
            break
        got.extend(batch)
        if len(got) >= total:
            break
    return total, got


def mentions_cldn4(blob: str) -> bool:
    text = blob.lower()
    return any(token in text for token in ("cldn4", "claudin-4", "claudin 4", "claudin4"))


def main() -> None:
    node_queries = [
        "CLDN4",
        "claudin-4",
        "claudin",
        "CLDN",
        "knockdown",
        "shRNA",
        "siRNA",
        "PRJCA023797",
        "HRA006761",
        "CRAD",
        "ELFN1",
    ]
    node_rows = []
    node_hits = []
    node_cldn = []
    for query in node_queries:
        total, items = node_search(query, page_size=30)
        cldn4_n = 0
        for item in items:
            blob = json.dumps(item, ensure_ascii=False)
            if mentions_cldn4(blob):
                cldn4_n += 1
                node_hits.append(
                    {
                        "query": query,
                        "type": item.get("type"),
                        "id": item.get("id"),
                        "name": item.get("name"),
                        "description": (item.get("description") or "")[:400],
                        "organism": item.get("organism"),
                        "tissue": item.get("tissue"),
                    }
                )
            if query == "CLDN":
                node_cldn.append(
                    {
                        "type": item.get("type"),
                        "id": item.get("id"),
                        "name": item.get("name"),
                        "description": (item.get("description") or "")[:240],
                        "mentions_cldn4": "yes" if mentions_cldn4(blob) else "no",
                    }
                )
        node_rows.append(
            {
                "database": "NODE",
                "query": query,
                "total": total,
                "returned": len(items),
                "mentions_cldn4": cldn4_n,
            }
        )
        print(f"NODE {query}: total={total} cldn4={cldn4_n}")

    cncb_plan = [
        ("gsa", "CLDN4"),
        ("gsa", "claudin-4"),
        ("gsa", "Cldn4"),
        ("hra", "CLDN4"),
        ("hra", "claudin-4"),
        ("hra", "claudin"),
        ("hra", "Cldn4"),
        ("hra", "PRJCA023797"),
        ("hra", "knockdown"),
        ("hra", "shRNA"),
        ("hra", "siRNA"),
        ("omix", "CLDN4"),
        ("omix", "claudin-4"),
        ("omix", "claudin"),
        ("omix", "CLDN"),
        ("omix", "knockdown"),
        ("bioproject", "CLDN4"),
        ("bioproject", "claudin-4"),
        ("bioproject", "CLDN4 silencing"),
        ("bioproject", "CLDN4 siRNA"),
    ]
    cncb_rows = []
    cncb_hits = []
    omix_cldn = []
    hra_kd_cldn4 = []
    for db, query in cncb_plan:
        # knockdown/shRNA/siRNA catalogs are large; keep one page of 50 plus a full
        # scan only when the query is the HRA knockdown family, capped at 200.
        max_pages = 4 if query in {"knockdown", "shRNA", "siRNA"} else 6
        page_size = 50
        total, items = cncb_search(db, query, page_size=page_size, max_pages=max_pages)
        cldn4_n = 0
        prjca_n = 0
        for item in items:
            blob = json.dumps(item, ensure_ascii=False)
            is_cldn4 = mentions_cldn4(blob)
            is_cn = item.get("id", "").startswith(("PRJCA", "CRA", "HRA", "OMIX", "OEP")) or (
                "PRJCA" in blob or "HRA" in blob[:80]
            )
            if is_cldn4:
                cldn4_n += 1
            if str(item.get("id", "")).startswith("PRJCA"):
                prjca_n += 1
            keep = is_cldn4 or (db == "omix" and query in {"claudin", "CLDN"}) or (
                db == "hra" and query == "PRJCA023797"
            )
            if keep and not (db == "omix" and query in {"claudin", "CLDN"}):
                cncb_hits.append(
                    {
                        "database": db,
                        "query": query,
                        "id": item.get("id"),
                        "title": item.get("title"),
                        "url": item.get("url"),
                        "description": (item.get("description") or "")[:400],
                        "mentions_cldn4": "yes" if is_cldn4 else "no",
                    }
                )
            if db == "omix" and query in {"claudin", "CLDN"}:
                omix_cldn.append(
                    {
                        "query": query,
                        "id": item.get("id"),
                        "title": item.get("title"),
                        "url": item.get("url"),
                        "description": (item.get("description") or "")[:300],
                        "mentions_cldn4": "yes" if is_cldn4 else "no",
                    }
                )
            if db == "hra" and query in {"knockdown", "shRNA", "siRNA"} and is_cldn4:
                hra_kd_cldn4.append(
                    {
                        "query": query,
                        "id": item.get("id"),
                        "title": item.get("title"),
                        "description": (item.get("description") or "")[:300],
                    }
                )
        cncb_rows.append(
            {
                "database": db,
                "query": query,
                "total": total,
                "returned": len(items),
                "mentions_cldn4": cldn4_n,
                "prjca_ids": prjca_n,
            }
        )
        print(f"CNCB {db} {query}: total={total} returned={len(items)} cldn4={cldn4_n}")

    # Public run catalog for the only Chinese CLDN4-titled lung study.
    run_rows = []
    page = 1
    while page <= 8:
        url = "https://ngdc.cncb.ac.cn/gsa-human/ajaxb/runinstudy?" + urllib.parse.urlencode(
            {"accession": "HRA006761", "pageNo": page, "pageSize": 10}
        )
        data = get_json(url)
        views = data.get("runViews") or []
        if not views:
            break
        for row in views:
            run_rows.append(
                {
                    "study": "HRA006761",
                    "bioproject": "PRJCA023797",
                    "access": "Controlled",
                    "run": row.get("runAcc"),
                    "file": row.get("runFileName"),
                    "file_size_gb": row.get("runFileSize"),
                    "experiment": row.get("expName"),
                    "experiment_title": row.get("expTitle"),
                    "strategy": row.get("strategyName"),
                    "platform": row.get("platform"),
                    "sample": row.get("sampleAcc"),
                    "biosample": row.get("biosampleAcc"),
                }
            )
        page += 1

    papers = [
        {
            "pmid": "26058359",
            "year": "2015",
            "organ": "breast",
            "perturbation": "CLDN4 shRNA and overexpression in MCF-7",
            "cldn4_knockdown": "yes",
            "affiliation": "Affiliated Hospital of Guangdong Medical College, Zhanjiang, China",
            "omics_accession": "",
            "in_node_or_gsa": "no",
            "note": "Functional assays and xenografts. No GEO, HRA, OMIX, or OEP in the PubMed record.",
        },
        {
            "pmid": "25871476",
            "year": "2015",
            "organ": "breast",
            "perturbation": "CLDN4 monoclonal antibody blockade of channel formation",
            "cldn4_knockdown": "no",
            "affiliation": "Academy of Military Medical Sciences, Beijing, China",
            "omics_accession": "",
            "in_node_or_gsa": "no",
            "note": "Vasculogenic mimicry in MDA-MB-231. Antibody blockade, not a deposited knockdown transcriptome.",
        },
        {
            "pmid": "30808546",
            "year": "2019",
            "organ": "breast",
            "perturbation": "PAK4 knockdown; CLDN4 restoration",
            "cldn4_knockdown": "no",
            "affiliation": "China Medical University, Shenyang, China",
            "omics_accession": "",
            "in_node_or_gsa": "no",
            "note": "CLDN4 is the downstream gene, not the silenced gene.",
        },
        {
            "pmid": "33006362",
            "year": "2020",
            "organ": "lung",
            "perturbation": "CRAD knockdown in A549 and H1299; microarray readout",
            "cldn4_knockdown": "no",
            "affiliation": "Jining Medical University, Shandong, China",
            "omics_accession": "",
            "in_node_or_gsa": "no",
            "note": "Abstract says Claudin 4 rose after CRAD silencing. GEO has no CRAD-knockdown lung series. NODE query CRAD is recorded in node_queries.tsv.",
        },
        {
            "pmid": "32779991",
            "year": "2020",
            "organ": "ovarian",
            "perturbation": "ELFN1-AS1 knockdown; CLDN4 is the miR-497-3p target",
            "cldn4_knockdown": "no",
            "affiliation": "Jiangxi Maternal and Child Health Hospital, Nanchang, China",
            "omics_accession": "",
            "in_node_or_gsa": "no",
            "note": "PMC8291874. Rescue experiments around the axis. No omics accession.",
        },
        {
            "pmid": "38629624",
            "year": "2024",
            "organ": "lung",
            "perturbation": "none; baseline MPE scRNA-seq, CLDN4 used as a tumor-cell marker",
            "cldn4_knockdown": "no",
            "affiliation": "Shanghai Pulmonary Hospital, Tongji University, Shanghai, China",
            "omics_accession": "HRA006761 / PRJCA023797",
            "in_node_or_gsa": "GSA-Human, controlled",
            "note": "Not a knockdown. FASTQ not downloaded.",
        },
    ]

    non_chinese_open = [
        {
            "database": "GEO (also indexed by CNCB GSA/BioProject as INSDC)",
            "accession": "GSE22493 / PRJNA128653",
            "organ": "ovarian",
            "design": "SKOV-3 CLDN4 siRNA lentivirus vs CLDN4-high control, microarray, n=3",
            "country": "USA",
            "institution": "Brigham and Women's Hospital, Harvard Medical School",
            "open": "yes",
            "downloaded_here": "no",
            "reason": "Not a Chinese NODE/GSA deposit. Series matrix stays at GEO.",
        },
        {
            "database": "GEO / SRA (GSA search returns SRX16067609 and SRX16067610)",
            "accession": "GSE207704 / PRJNA856719",
            "organ": "breast",
            "design": "T47D and MCF7 parental vs CLDN4 knockout RNA-seq",
            "country": "Japan",
            "institution": "Fukushima Medical University",
            "open": "yes",
            "downloaded_here": "no",
            "reason": "Not a Chinese paper. PMID 37059993.",
        },
        {
            "database": "GEO / SRA (GSA search returns SRX352051, SRX352053, SRX352054)",
            "accession": "GSE50927 / PRJNA219385",
            "organ": "mouse lung",
            "design": "WT vs Cldn4 knockout lungs after ventilator-induced lung injury",
            "country": "USA",
            "institution": "University of Southern California",
            "open": "yes",
            "downloaded_here": "no",
            "reason": "Not cancer and not a Chinese deposit. PMID 25106430.",
        },
    ]

    write_tsv(
        TABLES / "node_queries.tsv",
        node_rows,
        ["database", "query", "total", "returned", "mentions_cldn4"],
    )
    write_tsv(
        TABLES / "node_cldn4_hits.tsv",
        node_hits,
        ["query", "type", "id", "name", "description", "organism", "tissue"],
    )
    write_tsv(
        TABLES / "node_cldn_records.tsv",
        node_cldn,
        ["type", "id", "name", "description", "mentions_cldn4"],
    )
    write_tsv(
        TABLES / "cncb_queries.tsv",
        cncb_rows,
        ["database", "query", "total", "returned", "mentions_cldn4", "prjca_ids"],
    )
    write_tsv(
        TABLES / "cncb_cldn4_hits.tsv",
        cncb_hits,
        ["database", "query", "id", "title", "url", "description", "mentions_cldn4"],
    )
    write_tsv(
        TABLES / "omix_cldn.tsv",
        omix_cldn,
        ["query", "id", "title", "url", "description", "mentions_cldn4"],
    )
    write_tsv(
        TABLES / "hra_knockdown_cldn4.tsv",
        hra_kd_cldn4,
        ["query", "id", "title", "description"],
    )
    write_tsv(
        TABLES / "hra006761_files.tsv",
        run_rows,
        [
            "study",
            "bioproject",
            "access",
            "run",
            "file",
            "file_size_gb",
            "experiment",
            "experiment_title",
            "strategy",
            "platform",
            "sample",
            "biosample",
        ],
    )
    write_tsv(
        TABLES / "chinese_papers.tsv",
        papers,
        [
            "pmid",
            "year",
            "organ",
            "perturbation",
            "cldn4_knockdown",
            "affiliation",
            "omics_accession",
            "in_node_or_gsa",
            "note",
        ],
    )
    write_tsv(
        TABLES / "open_non_chinese_kd.tsv",
        non_chinese_open,
        [
            "database",
            "accession",
            "organ",
            "design",
            "country",
            "institution",
            "open",
            "downloaded_here",
            "reason",
        ],
    )

    unique_runs = sorted({row["run"] for row in run_rows})
    size_gb = 0.0
    for row in run_rows:
        try:
            size_gb += float(row["file_size_gb"] or 0)
        except ValueError:
            pass
    gate = {
        "question": "CLDN4 knockdown lung, ovarian, or breast datasets in NODE or GSA-family China archives",
        "node_cldn4_records": sum(row["mentions_cldn4"] for row in node_rows if row["query"] == "CLDN4"),
        "hra_cldn4_records": next(row["total"] for row in cncb_rows if row["database"] == "hra" and row["query"] == "CLDN4"),
        "omix_cldn4_records": next(row["total"] for row in cncb_rows if row["database"] == "omix" and row["query"] == "CLDN4"),
        "gsa_cldn4_records_insdc_mirror": next(
            row["total"] for row in cncb_rows if row["database"] == "gsa" and row["query"] == "CLDN4"
        ),
        "chinese_cldn4_knockdown_accessions": 0,
        "files_downloaded": 0,
        "nearest_chinese_lung_study": {
            "accession": "HRA006761",
            "bioproject": "PRJCA023797",
            "pmid": "38629624",
            "access": "Controlled",
            "design": "16 baseline MPE scRNA-seq samples, recurrent vs non-recurrent NSCLC, not a CLDN4 knockdown",
            "public_run_rows": len(run_rows),
            "unique_runs": len(unique_runs),
            "listed_fastq_gb": round(size_gb, 2),
            "downloaded": False,
        },
        "hra_knockdown_records_mentioning_cldn4": len(hra_kd_cldn4),
    }
    (ROOT / "results" / "gate.json").write_text(json.dumps(gate, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
