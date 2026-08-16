#!/usr/bin/env python3
"""Curate high-confidence lung IO datasets and download protein tables only."""

from __future__ import annotations

import csv
import json
import time
import urllib.request
from pathlib import Path

UA = "grok-pride2/1.0 (research; protein-tables-only)"
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "grok_pride2"
TABLES = RESULTS / "protein_tables"
TABLES.mkdir(parents=True, exist_ok=True)

# High-confidence lung + ICI / PD-1 / immunotherapy MS datasets beyond excluded 3.
# confidence: A = patient ICI-treated proteome; B = lung IO mechanism / immunopeptidome;
# C = companion / cell-line / adjacent immune-evasion atlas.
CURATED = [
    {
        "accession": "PXD078247",
        "confidence": "A",
        "why": "Serum proteome of NSCLC patients treated with ICB monotherapy",
        "host": "PRIDE/iProX?",
    },
    {
        "accession": "PXD062630",
        "confidence": "A",
        "why": "Urine EV proteome predicts anti-PD1 PFS in advanced NSCLC (n=33)",
        "host": "PRIDE",
    },
    {
        "accession": "PXD031474",
        "confidence": "A",
        "why": "DIA plasma biomarkers, first-line anti-PD-1 NSCLC (iProX)",
        "host": "iProX",
    },
    {
        "accession": "PXD039141",
        "confidence": "A",
        "why": "Serum glycoproteins before/during anti-PD-1/PD-L1 in metastatic LUAD/LUSC (iProX)",
        "host": "iProX",
    },
    {
        "accession": "PXD058967",
        "confidence": "A",
        "why": "Plasma multi-omics (metabolome+proteome) unresectable NSCLC immunotherapy (iProX)",
        "host": "iProX",
    },
    {
        "accession": "PXD066804",
        "confidence": "A",
        "why": "LDN myeloid proteome, first-line anti-PD1/PD-L1 resistance in NSCLC (jPOST)",
        "host": "jPOST",
    },
    {
        "accession": "PXD040761",
        "confidence": "B",
        "why": "Durvalumab + chemo proteome in LUAD cell lines",
        "host": "PRIDE",
    },
    {
        "accession": "PXD019573",
        "confidence": "B",
        "why": "PSME4/PA200 NSCLC proteome; ratio associated with durvalumab response",
        "host": "PRIDE",
    },
    {
        "accession": "PXD028364",
        "confidence": "C",
        "why": "A549 immunopeptidome companion of PSME4/PA200 series",
        "host": "PRIDE",
    },
    {
        "accession": "PXD037365",
        "confidence": "C",
        "why": "A549 immunopeptidome companion of PSME4/PA200 series",
        "host": "PRIDE",
    },
    {
        "accession": "PXD065735",
        "confidence": "A",
        "why": "SCLC ERBB2 immune evasion / ICI resistance (PRIDE); MassIVE MSV000098579",
        "host": "PRIDE",
    },
    {
        "accession": "MSV000098579",
        "confidence": "A",
        "why": "MassIVE mirror of PXD065735 SCLC ICI resistance",
        "host": "MassIVE",
    },
    {
        "accession": "PXD058303",
        "confidence": "B",
        "why": "SCLC immunopeptidome; ATAD2 as immunotherapy target",
        "host": "PRIDE/iProX?",
    },
    {
        "accession": "PXD081026",
        "confidence": "B",
        "why": "NSCLC MUC1-MAPK-Elk-1-PD-L1 axis proteome / nanoparticle study",
        "host": "PRIDE/iProX?",
    },
    {
        "accession": "PXD044086",
        "confidence": "B",
        "why": "Human lung cancer membrane TMT to find new immunotherapy targets (PD-1/PD-L1 context)",
        "host": "PRIDE",
    },
    {
        "accession": "PXD034772",
        "confidence": "B",
        "why": "Lung cancer immunopeptidome vs T-cell infiltration / ICI relevance",
        "host": "PRIDE",
    },
    {
        "accession": "PXD043057",
        "confidence": "B",
        "why": "MHC-I immunopeptidome of human NSCLC in ICB-failure context",
        "host": "PRIDE",
    },
    {
        "accession": "PXD022949",
        "confidence": "B",
        "why": "HLA-I immunopeptidome EGFR-mutant LUAD (+ melanoma)",
        "host": "PRIDE",
    },
    {
        "accession": "PXD027766",
        "confidence": "B",
        "why": "sHLA peptidome of pleural effusions; lung/tumor antigens",
        "host": "PRIDE",
    },
    {
        "accession": "PXD042991",
        "confidence": "B",
        "why": "Multi-target molecule restores immune response in lung cancer (mouse PISA proteome)",
        "host": "PRIDE",
    },
    {
        "accession": "PXD049438",
        "confidence": "B",
        "why": "Time-serial proteome of lung cancer TIL cell therapy",
        "host": "PRIDE",
    },
    {
        "accession": "PXD020191",
        "confidence": "C",
        "why": "Lehtiö NSCLC proteogenomics; immune-evasion subtypes (not ICI-treated)",
        "host": "PRIDE",
    },
    {
        "accession": "PXD020548",
        "confidence": "C",
        "why": "Lehtiö DIA validation NSCLC tissue cohort",
        "host": "PRIDE",
    },
    {
        "accession": "MSV000085049",
        "confidence": "C",
        "why": "MassIVE global companion of excluded PXD019061 (not a new study)",
        "host": "MassIVE",
    },
]

# Explicit protein-table URLs (HTTP, small enough). No raw MS.
DIRECT = [
    {
        "accession": "PXD062630",
        "name": "011321-Koranyi-NSCLC-urine-EVs-quant_reanalysis_with_bacteria_viruses_fungi.mzTab",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/04/PXD062630/011321-Koranyi-NSCLC-urine-EVs-quant_reanalysis_with_bacteria_viruses_fungi.mzTab",
    },
    {
        "accession": "PXD042991",
        "name": "20221213_nU5_Exp480_1_SApcher_PISA_2022_mouse_proteins.xlsx",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD042991/20221213_nU5_Exp480_1_SApcher_PISA_2022_mouse_proteins.xlsx",
    },
    {
        "accession": "PXD065735",
        "name": "20250507_130654_4351_PG_Report.tsv",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/10/PXD065735/20250507_130654_4351_PG_Report.tsv",
    },
    {
        "accession": "PXD019573",
        "name": "DataDepositionTable.xlsx",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/03/PXD019573/DataDepositionTable.xlsx",
    },
    {
        "accession": "PXD027766",
        "name": "Samples_raw_files_and_description.txt",
        "url": "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/05/PXD027766/Samples_raw_files_and_description.txt",
    },
]


def http_json(url: str, timeout: int = 30):
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read()
        return json.loads(raw.decode("utf-8", errors="replace")) if raw else None


def download(url: str, dest: Path, max_bytes: int = 80 * 1024 * 1024) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return {"ok": True, "reason": "exists", "bytes": dest.stat().st_size}
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    written = 0
    try:
        with urllib.request.urlopen(req, timeout=180) as resp, tmp.open("wb") as fh:
            cl = resp.headers.get("Content-Length")
            if cl and int(cl) > max_bytes:
                return {"ok": False, "reason": f"too large {cl}", "bytes": 0}
            while True:
                chunk = resp.read(256 * 1024)
                if not chunk:
                    break
                written += len(chunk)
                if written > max_bytes:
                    tmp.unlink(missing_ok=True)
                    return {"ok": False, "reason": "stream cap", "bytes": written}
                fh.write(chunk)
        tmp.replace(dest)
        return {"ok": True, "reason": "downloaded", "bytes": written}
    except Exception as exc:  # noqa: BLE001
        tmp.unlink(missing_ok=True)
        return {"ok": False, "reason": str(exc), "bytes": written}


def ftp_from_px(acc: str) -> str | None:
    try:
        data = http_json(
            f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}&outputMode=JSON&test=no"
        )
    except Exception:
        return None
    if not data:
        return None
    for key in ("fullDatasetLinkList", "datasetLink", "datasetFiles"):
        val = data.get(key)
        if isinstance(val, list):
            for item in val:
                if isinstance(item, dict):
                    u = item.get("value") or item.get("url") or ""
                    if "ftp.pride" in u or "ftp.ebi" in u:
                        return u.replace("ftp://", "https://")
                elif isinstance(item, str) and "ftp" in item:
                    return item.replace("ftp://", "https://")
        elif isinstance(val, str) and "ftp" in val:
            return val.replace("ftp://", "https://")
    return None


def pride_files(acc: str) -> list[dict]:
    try:
        data = http_json(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{acc}/files")
        return data if isinstance(data, list) else []
    except Exception:
        return []


def main() -> None:
    catalog = {r["accession"]: r for r in csv.DictReader((RESULTS / "px_catalog.csv").open())}
    rows = []
    for item in CURATED:
        acc = item["accession"]
        meta = catalog.get(acc, {})
        rows.append(
            {
                **item,
                "title": meta.get("title") or "",
                "description": (meta.get("description") or "")[:1500],
                "excluded_known": acc in {"PXD042091", "PXD059688", "PXD019061"},
                "px_url": f"https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID={acc}",
            }
        )
    out = RESULTS / "curated_lung_io.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # Fetch file lists for curated PXDs missing from inventory
    extra_files = []
    for item in CURATED:
        acc = item["accession"]
        if not acc.startswith("PXD"):
            continue
        files = pride_files(acc)
        time.sleep(0.05)
        print(f"{acc}: {len(files)} pride files", flush=True)
        for f in files:
            name = f.get("fileName") or ""
            cat = ((f.get("fileCategory") or {}).get("value") or "").upper()
            size = int(f.get("fileSizeBytes") or 0)
            http = ""
            for loc in f.get("publicFileLocations") or []:
                val = loc.get("value") or ""
                if val.startswith("http"):
                    http = val
                elif val.startswith("ftp://ftp.pride.ebi.ac.uk"):
                    http = "https" + val[3:]
            extra_files.append(
                {
                    "accession": acc,
                    "fileName": name,
                    "category": cat,
                    "size_bytes": size,
                    "download_http": http,
                }
            )
    with (RESULTS / "curated_file_inventory.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["accession", "fileName", "category", "size_bytes", "download_http"]
        )
        w.writeheader()
        w.writerows(extra_files)

    log = []
    for spec in DIRECT:
        dest = TABLES / spec["accession"] / spec["name"]
        print(f"GET {spec['accession']} {spec['name']}", flush=True)
        result = download(spec["url"], dest)
        log.append({**spec, **result, "local": str(dest.relative_to(ROOT)) if dest.exists() else ""})
        time.sleep(0.1)

    # Additional small protein-like files from freshly fetched lists
    protein_name = (
        "protein",
        "pg_report",
        "mztab",
        "proteins.xlsx",
        "proteins.txt",
        "proteingroups",
        "sdrf",
    )
    skip_ext = (".raw", ".wiff", ".mzml", ".mzxml", ".mgf", ".dta", ".apl", ".dat", ".msf", ".sne")
    for f in extra_files:
        name = f["fileName"]
        low = name.lower()
        if any(low.endswith(ext) for ext in skip_ext):
            continue
        if f["size_bytes"] and f["size_bytes"] > 80 * 1024 * 1024:
            continue
        if not any(tok in low for tok in protein_name):
            continue
        if not f["download_http"]:
            continue
        dest = TABLES / f["accession"] / name
        if dest.exists():
            continue
        print(f"GET extra {f['accession']} {name}", flush=True)
        result = download(f["download_http"], dest)
        log.append(
            {
                "accession": f["accession"],
                "name": name,
                "url": f["download_http"],
                **result,
                "local": str(dest.relative_to(ROOT)) if dest.exists() else "",
            }
        )
        time.sleep(0.1)

    (RESULTS / "curated_download_log.json").write_text(json.dumps(log, indent=2), encoding="utf-8")
    print("curated rows", len(rows), "downloads", len(log))


if __name__ == "__main__":
    main()
