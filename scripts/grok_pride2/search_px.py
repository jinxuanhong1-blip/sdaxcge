#!/usr/bin/env python3
"""Exhaustive PRIDE / MassIVE / ProteomeXchange / OmicsDI search.

Scope: lung + immunotherapy / PD-1 / TROP2 / CLDN4 proteomics.
Excluded (already known): PXD042091, PXD059688, PXD019061.
Protein-table discovery only — never download raw MS in this script.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

UA = "grok-pride2/1.0 (research; protein-tables-only; no-raw-ms)"
EXCLUDED = {"PXD042091", "PXD059688", "PXD019061"}
ROOT = Path(__file__).resolve().parents[2]
NOTES = ROOT / "notes" / "grok_pride2"
RESULTS = ROOT / "results" / "grok_pride2"
NOTES.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)

LUNG_RE = re.compile(
    r"\b(lung|nsclc|luad|lusc|sclc|pulmonary|bronch|pleural|"
    r"adenocarcinoma of the lung|lung adenocarcinoma|lung squamous|"
    r"non[- ]small[- ]cell lung|small[- ]cell lung)\b",
    re.I,
)
IO_RE = re.compile(
    r"\b(immunotherap(?:y|ies)|immuno[- ]oncolog(?:y|ies)|checkpoint|"
    r"anti[- ]?pd[- ]?[l1]|pd[- ]?1|pd[- ]?l1|pd[- ]?l2|pdcd1|cd274|"
    r"pembrolizumab|nivolumab|atezolizumab|durvalumab|cemiplimab|"
    r"ipilimumab|tremelimumab|camrelizumab|sintilimab|tislelizumab|"
    r"toripalimab|sugemalimab|keynote|checkmate|impower|pacific|"
    r"ici\b|icb\b|immune checkpoint)\b",
    re.I,
)
TROP2_RE = re.compile(r"\b(trop[- ]?2|tacstd2|ega[mp]1|tumor[- ]associated calcium)\b", re.I)
CLDN4_RE = re.compile(r"\b(cldn4|claudin[- ]?4)\b", re.I)
PROTEOME_SOURCES = {
    "pride",
    "massive",
    "iprox",
    "jpost",
    "panorama",
    "peptideatlas",
    "passel",
    "proteomexchange",
}

LITERATURE_SEEDS = [
    # Core ICI / PD-1 lung proteomes beyond the excluded three
    "PXD031474",  # iProX DIA plasma first-line anti-PD-1 NSCLC
    "PXD039141",  # iProX serum glycoproteins anti-PD-1/PD-L1 LUAD/LUSC
    "PXD058967",  # iProX multi-omics unresectable NSCLC immunotherapy
    "PXD066804",  # jPOST LDN myeloid resistance first-line anti-PD1/PD-L1
    "PXD019573",  # PSME4/PA200 NSCLC + Durvalumab association
    "PXD028364",  # PSME4 companion
    "PXD037365",  # PSME4 companion
    "PXD020191",  # Lehtiö NSCLC proteogenomics / immune evasion
    "PXD020548",  # Lehtiö DIA validation cohort
    "MSV000085049",  # MassIVE global companion of excluded PXD019061
    "PXD044086",
    "PXD047616",
    "PXD054390",
    "PXD027259",
    "PXD010154",  # CPTAC LUAD (TROP2/CLDN4 protein context)
    "PXD013455",
    "PXD016000",
    "PXD030014",
    "IPX0003936000",
    "IPX0004653000",
    "JPST003908",
]

PRIDE_KEYWORDS = [
    "immunotherapy",
    "pembrolizumab",
    "nivolumab",
    "atezolizumab",
    "durvalumab",
    "ipilimumab",
    "camrelizumab",
    "sintilimab",
    "tislelizumab",
    "PD-1",
    "PD-L1",
    "PD1",
    "PDL1",
    "checkpoint",
    "TROP2",
    "TACSTD2",
    "CLDN4",
    "claudin-4",
    "NSCLC immunotherapy",
    "lung immunotherapy",
    "lung PD-1",
    "lung PD-L1",
]

OMICSDI_QUERIES = [
    "(source:pride OR source:massive OR source:iprox OR source:jpost OR source:panorama) AND (lung OR NSCLC OR LUAD OR LUSC) AND (immunotherapy OR pembrolizumab OR nivolumab OR atezolizumab OR durvalumab)",
    "(source:pride OR source:massive OR source:iprox OR source:jpost) AND (lung OR NSCLC) AND (PD-1 OR PD-L1 OR PD1 OR checkpoint)",
    "(source:pride OR source:massive OR source:iprox OR source:jpost) AND (TROP2 OR TACSTD2 OR CLDN4) AND (lung OR NSCLC OR cancer)",
    "id:PXD* AND NSCLC AND (immunotherapy OR PD-1 OR PD-L1)",
    "id:MSV* AND (lung OR NSCLC) AND (immunotherapy OR PD-1 OR PD-L1 OR TROP2 OR CLDN4)",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def http_json(url: str, timeout: int = 25, retries: int = 2):
    last = None
    for i in range(retries):
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                raw = resp.read()
                if not raw:
                    return None
                return json.loads(raw.decode("utf-8", errors="replace"))
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(0.8 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def http_headers(url: str, timeout: int = 45):
    req = urllib.request.Request(url, headers={"User-Agent": UA}, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return {k.lower(): v for k, v in resp.headers.items()}
    except Exception:
        return {}


def norm_acc(value: str | None) -> str | None:
    if not value:
        return None
    value = value.strip().upper()
    m = re.search(r"\b(PXD\d{5,}|MSV\d{9}|PASS\d+|JPST\d+|IPX\d+)\b", value)
    return m.group(1) if m else None


def text_blob(obj) -> str:
    parts: list[str] = []

    def walk(x):
        if isinstance(x, dict):
            for v in x.values():
                walk(v)
        elif isinstance(x, list):
            for v in x:
                walk(v)
        elif isinstance(x, str):
            parts.append(x)

    walk(obj)
    return "\n".join(parts)


def score_record(title: str, description: str, keywords: str, extra: str = "") -> dict:
    blob = " ".join([title or "", description or "", keywords or "", extra or ""])
    lung = bool(LUNG_RE.search(blob))
    io = bool(IO_RE.search(blob))
    trop2 = bool(TROP2_RE.search(blob))
    cldn4 = bool(CLDN4_RE.search(blob))
    if lung and io:
        tier = "core_lung_io"
    elif lung and (trop2 or cldn4):
        tier = "lung_trop2_cldn4"
    elif io and (trop2 or cldn4):
        tier = "io_trop2_cldn4_not_lung"
    elif lung:
        tier = "lung_adjacent"
    elif io:
        tier = "io_not_lung"
    else:
        tier = "unrelated"
    return {
        "lung": lung,
        "io": io,
        "trop2": trop2,
        "cldn4": cldn4,
        "tier": tier,
        "keep": tier in {"core_lung_io", "lung_trop2_cldn4"},
    }


def pride_search_keyword(keyword: str, page_size: int = 100, max_pages: int = 8) -> list[dict]:
    hits = []
    page = 0
    while page < max_pages:
        qs = urllib.parse.urlencode(
            {"keyword": keyword, "pageSize": page_size, "page": page}
        )
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?{qs}"
        data = http_json(url) or []
        if not isinstance(data, list) or not data:
            break
        hits.extend(data)
        if len(data) < page_size:
            break
        page += 1
        time.sleep(0.05)
    return hits


def omicsdi_search(query: str, size: int = 100, max_records: int = 300) -> list[dict]:
    hits = []
    start = 0
    while start < max_records:
        qs = urllib.parse.urlencode({"query": query, "size": size, "start": start})
        url = f"https://www.omicsdi.org/ws/dataset/search?{qs}"
        data = http_json(url) or {}
        datasets = data.get("datasets") or []
        if not datasets:
            break
        hits.extend(datasets)
        total = int(data.get("count") or 0)
        start += size
        if start >= total or len(datasets) < size:
            break
        time.sleep(0.2)
    return hits


def px_dataset(accession: str) -> dict | None:
    qs = urllib.parse.urlencode(
        {"ID": accession, "outputMode": "JSON", "test": "no"}
    )
    url = f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?{qs}"
    try:
        return http_json(url, timeout=45)
    except Exception:
        return None


def pride_project(accession: str) -> dict | None:
    try:
        return http_json(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}")
    except Exception:
        return None


def pride_files(accession: str) -> list[dict]:
    try:
        data = http_json(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files")
        return data if isinstance(data, list) else []
    except Exception:
        return []


def extract_title_desc(rec: dict) -> tuple[str, str, str]:
    title = rec.get("title") or rec.get("name") or ""
    desc = (
        rec.get("projectDescription")
        or rec.get("description")
        or rec.get("summary")
        or rec.get("abstract")
        or ""
    )
    if isinstance(desc, list):
        desc = " ".join(str(x) for x in desc)
    kws = rec.get("keywords") or rec.get("projectTags") or rec.get("keywordsList") or []
    if isinstance(kws, list):
        kws = ";".join(str(x) for x in kws)
    return str(title), str(desc), str(kws)


def main() -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    seen: dict[str, dict] = {}
    query_log = []

    def add_hit(acc: str, source: str, query: str, rec: dict | None = None):
        acc = norm_acc(acc) or acc
        if not acc:
            return
        item = seen.setdefault(
            acc,
            {
                "accession": acc,
                "sources": set(),
                "queries": set(),
                "raw_hits": [],
                "excluded": acc in EXCLUDED,
            },
        )
        item["sources"].add(source)
        item["queries"].add(query)
        if rec:
            item["raw_hits"].append({"source": source, "record": rec})

    log(f"[pride] keyword searches: {len(PRIDE_KEYWORDS)}")
    for kw in PRIDE_KEYWORDS:
        try:
            hits = pride_search_keyword(kw)
        except Exception as exc:
            query_log.append({"api": "pride", "query": kw, "error": str(exc), "n": 0})
            log(f"  PRIDE fail {kw}: {exc}")
            continue
        query_log.append({"api": "pride", "query": kw, "n": len(hits)})
        log(f"  PRIDE {kw!r}: {len(hits)}")
        for rec in hits:
            acc = rec.get("accession") or rec.get("projectAccession")
            if acc:
                add_hit(acc, "pride", f"keyword:{kw}", rec)

    log("[omicsdi] queries")
    for q in OMICSDI_QUERIES:
        try:
            hits = omicsdi_search(q)
        except Exception as exc:
            query_log.append({"api": "omicsdi", "query": q, "error": str(exc), "n": 0})
            log(f"  OmicsDI fail: {exc}")
            continue
        query_log.append({"api": "omicsdi", "query": q, "n": len(hits)})
        log(f"  OmicsDI {q[:70]}...: {len(hits)}")
        for rec in hits:
            source = (rec.get("source") or "").lower()
            if source and source not in PROTEOME_SOURCES and "pride" not in source and "massive" not in source:
                # keep only proteomics repositories from OmicsDI
                if source not in {"iprox", "jpost", "panorama", "peptideatlas", "passel"}:
                    continue
            acc = rec.get("id") or rec.get("accession")
            add_hit(acc, f"omicsdi:{source or 'unknown'}", q, rec)

    log("[seeds] literature / companion accessions")
    for acc in LITERATURE_SEEDS:
        add_hit(acc, "literature_seed", "seed", {"accession": acc})

    log(f"[prefilter] {len(seen)} unique accessions")
    prelim = []
    for acc, item in seen.items():
        rec0 = item["raw_hits"][0]["record"] if item["raw_hits"] else {}
        title, desc, kws = extract_title_desc(rec0)
        extra = text_blob(rec0)
        score = score_record(title, desc, kws, extra)
        # Enrich seeds, excluded knowns, keep hits, and ICI-drug titles
        # even if the snippet omitted the word "lung".
        drugish = bool(
            re.search(
                r"pembrolizumab|nivolumab|atezolizumab|durvalumab|camrelizumab|"
                r"sintilimab|tislelizumab|ipilimumab|anti[- ]?pd",
                " ".join([title, desc, kws]),
                re.I,
            )
        )
        enrich = (
            item["excluded"]
            or acc in {norm_acc(x) or x for x in LITERATURE_SEEDS}
            or score["keep"]
            or drugish
            or any(s.startswith("literature") for s in item["sources"])
        )
        prelim.append((acc, item, score, title, desc, kws, enrich))

    to_enrich = [p for p in prelim if p[6]]
    log(f"[enrich] {len(to_enrich)} accessions (skip {len(prelim) - len(to_enrich)} unrelated)")
    rows = []
    for acc, item, pre_score, title0, desc0, kws0, enrich in sorted(prelim, key=lambda x: x[0]):
        title, desc, kws = title0, desc0, kws0
        host = ""
        organisms = []
        px = None
        pride = None
        if enrich and re.match(r"^(PXD|MSV|PASS|JPST|IPX)", acc):
            try:
                px = px_dataset(acc)
            except Exception as exc:
                log(f"  PX fail {acc}: {exc}")
            time.sleep(0.04)
            if acc.startswith("PXD"):
                try:
                    pride = pride_project(acc)
                except Exception as exc:
                    log(f"  PRIDE project fail {acc}: {exc}")
                time.sleep(0.04)
        rec0 = item["raw_hits"][0]["record"] if item["raw_hits"] else {}
        if pride or px:
            title, desc, kws = extract_title_desc(pride or px or rec0)
        if px:
            title = title or px.get("title") or title0
            desc = desc or px.get("description") or px.get("summary") or desc0
            host = px.get("hostingRepository") or px.get("repositoryName") or ""
        if pride:
            orgs = pride.get("organisms") or []
            organisms = [
                (o.get("name") if isinstance(o, dict) else str(o)) for o in orgs
            ]
            host = host or "PRIDE"
        extra = text_blob(pride or px or rec0)
        score = score_record(title, desc, kws, extra)
        row = {
            "accession": acc,
            "excluded_known": item["excluded"],
            "tier": score["tier"],
            "keep_for_download": bool(score["keep"] and not item["excluded"]),
            "lung": score["lung"],
            "io": score["io"],
            "trop2": score["trop2"],
            "cldn4": score["cldn4"],
            "title": title,
            "description": (desc or "")[:2000],
            "keywords": kws,
            "host": host,
            "organisms": ";".join(organisms),
            "sources": ";".join(sorted(item["sources"])),
            "queries": " | ".join(sorted(item["queries"]))[:500],
            "pride_url": f"https://www.ebi.ac.uk/pride/archive/projects/{acc}"
            if acc.startswith("PXD")
            else "",
            "px_url": f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}",
        }
        rows.append(row)

    # File inventory for keep + excluded (document why excluded)
    file_rows = []
    targets = [r for r in rows if r["keep_for_download"] or r["excluded_known"]]
    log(f"[files] listing for {len(targets)} keep/excluded projects")
    for r in targets:
        acc = r["accession"]
        if not acc.startswith("PXD"):
            continue
        files = pride_files(acc)
        time.sleep(0.1)
        proteinish = 0
        rawish = 0
        for f in files:
            name = f.get("fileName") or ""
            cat = ((f.get("fileCategory") or {}).get("value") or "").upper()
            size = int(f.get("fileSizeBytes") or 0)
            locs = f.get("publicFileLocations") or []
            http = ""
            for loc in locs:
                val = loc.get("value") or ""
                if val.startswith("http"):
                    http = val
                elif val.startswith("ftp://ftp.pride.ebi.ac.uk"):
                    http = "https" + val[3:]
            is_raw = cat == "RAW" or re.search(
                r"\.(raw|wiff|wiff\.scan|mzml|mzxml|d|raw\.gz)$", name, re.I
            )
            is_protein = bool(
                re.search(
                    r"(protein|proteingroup|proteinGroups|pg_matrix|diann-output|"
                    r"quant|maxquant|pd_proteins|proteome.?discover|mzTab|"
                    r"search\.zip|txt\.zip|tables|supplement)",
                    name,
                    re.I,
                )
            ) and not is_raw
            if is_raw:
                rawish += 1
            if is_protein or (cat in {"SEARCH", "RESULT", "OTHER", "PEAK"} and not is_raw):
                proteinish += 1
            file_rows.append(
                {
                    "accession": acc,
                    "fileName": name,
                    "category": cat,
                    "size_bytes": size,
                    "is_raw": is_raw,
                    "protein_table_candidate": is_protein and not is_raw,
                    "download_http": http,
                    "tier": r["tier"],
                    "keep_for_download": r["keep_for_download"],
                }
            )
        r["n_files"] = len(files)
        r["n_raw_files"] = rawish
        r["n_proteinish_files"] = proteinish

    # Write outputs
    catalog_path = RESULTS / "px_catalog.csv"
    fields = [
        "accession",
        "excluded_known",
        "tier",
        "keep_for_download",
        "lung",
        "io",
        "trop2",
        "cldn4",
        "title",
        "host",
        "organisms",
        "n_files",
        "n_raw_files",
        "n_proteinish_files",
        "sources",
        "keywords",
        "description",
        "pride_url",
        "px_url",
        "queries",
    ]
    with catalog_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["tier"], x["accession"])):
            w.writerow(r)

    keep_path = RESULTS / "px_keep.csv"
    keep_rows = [r for r in rows if r["keep_for_download"]]
    with keep_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(sorted(keep_rows, key=lambda x: x["accession"]))

    files_path = RESULTS / "px_file_inventory.csv"
    if file_rows:
        with files_path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(file_rows[0].keys()))
            w.writeheader()
            w.writerows(file_rows)

    summary = {
        "generated_utc": stamp,
        "n_unique_accessions": len(rows),
        "n_excluded_known": sum(1 for r in rows if r["excluded_known"]),
        "tier_counts": dict(sorted(defaultdict(int, {k: 0 for k in []}).items())),
        "keep_accessions": [r["accession"] for r in keep_rows],
        "excluded_known": sorted(EXCLUDED),
        "query_log": query_log,
    }
    tier_counts = defaultdict(int)
    for r in rows:
        tier_counts[r["tier"]] += 1
    summary["tier_counts"] = dict(sorted(tier_counts.items()))
    (RESULTS / "px_search_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    (NOTES / "search_log.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    log(json.dumps(summary["tier_counts"], indent=2))
    log("keep: " + ", ".join(summary["keep_accessions"]))
    log("wrote " + str(catalog_path))


if __name__ == "__main__":
    main()
