#!/usr/bin/env python3
"""Download public TROP2 / CLDN4 lung IHC-mIF scoring tables and supplements.

Outputs only under results/grok_ihc/. Network failures are logged, not raised.
"""
from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "results" / "grok_ihc" / "raw"
EXTRACTED = ROOT / "results" / "grok_ihc" / "extracted"
INVENTORY = ROOT / "results" / "grok_ihc" / "inventory"
NOTES = ROOT / "notes" / "grok_ihc"

UA = (
    "Mozilla/5.0 (compatible; grok-ihc-slice/1.0; research; "
    "+https://github.com/) Python-urllib"
)
TIMEOUT = 90


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(url: str, dest: Path | None = None, retries: int = 3) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True) if dest else None
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
            with urlopen(req, timeout=TIMEOUT) as resp:
                data = resp.read()
                info = {
                    "url": url,
                    "final_url": getattr(resp, "url", url),
                    "status": getattr(resp, "status", 200),
                    "content_type": resp.headers.get("Content-Type", ""),
                    "bytes": len(data),
                    "sha256": sha256_bytes(data),
                    "retrieved_at": utc_now(),
                    "ok": True,
                    "error": "",
                    "path": "",
                }
                if dest:
                    dest.write_bytes(data)
                    info["path"] = str(dest.relative_to(ROOT))
                return info | {"_data": data}
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last_err = exc
            time.sleep(min(2**attempt, 8))
    return {
        "url": url,
        "final_url": url,
        "status": getattr(last_err, "code", None),
        "content_type": "",
        "bytes": 0,
        "sha256": "",
        "retrieved_at": utc_now(),
        "ok": False,
        "error": f"{type(last_err).__name__}: {last_err}",
        "path": "",
        "_data": b"",
    }


def json_get(url: str) -> tuple[dict | list | None, dict]:
    rec = fetch(url)
    if not rec["ok"]:
        return None, rec
    try:
        return json.loads(rec["_data"].decode("utf-8", errors="replace")), rec
    except json.JSONDecodeError as exc:
        rec["ok"] = False
        rec["error"] = f"JSONDecodeError: {exc}"
        return None, rec


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def safe_name(name: str) -> str:
    name = re.sub(r"[^\w.\-+=]+", "_", name)
    return name[:180] or "unnamed"


def download_named(url: str, dest_dir: Path, filename: str) -> dict:
    dest = dest_dir / safe_name(filename)
    rec = fetch(url, dest)
    rec.pop("_data", None)
    rec["filename"] = filename
    return rec


def figshare_article(article_id: int, dest_dir: Path) -> list[dict]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    recs = []
    payload, meta = json_get(f"https://api.figshare.com/v2/articles/{article_id}")
    write_json(dest_dir / f"figshare_{article_id}_meta.json", payload or {"error": meta})
    recs.append({**{k: v for k, v in meta.items() if k != "_data"}, "filename": f"figshare_{article_id}_meta.json"})
    if not payload:
        return recs
    for f in payload.get("files", []):
        name = f.get("name") or f"file_{f.get('id')}"
        url = f.get("download_url") or f"https://ndownloader.figshare.com/files/{f.get('id')}"
        recs.append(download_named(url, dest_dir, name))
    return recs


def europepmc_fulltext_links(pmcid: str) -> list[str]:
    url = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
    rec = fetch(url)
    if not rec["ok"]:
        return []
    text = rec["_data"].decode("utf-8", errors="replace")
    hrefs = re.findall(r'xlink:href="([^"]+)"', text)
    hrefs += re.findall(r'<supplementary-material[^>]+href="([^"]+)"', text)
    return hrefs


def pmc_oa_package(pmcid: str, dest_dir: Path) -> list[dict]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    recs = []
    rec = fetch(f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmcid}")
    xml = rec["_data"].decode("utf-8", errors="replace") if rec["ok"] else rec.get("error", "")
    (dest_dir / f"{pmcid}_oa.xml").write_text(xml, encoding="utf-8")
    recs.append({**{k: v for k, v in rec.items() if k != "_data"}, "filename": f"{pmcid}_oa.xml"})
    links = re.findall(r'href="(ftp://[^"]+|https://[^"]+)"', xml)
    for link in links:
        http = link.replace("ftp://ftp.ncbi.nlm.nih.gov/", "https://ftp.ncbi.nlm.nih.gov/")
        name = http.rstrip("/").split("/")[-1]
        recs.append(download_named(http, dest_dir, name))
    return recs


def extract_zip(zip_path: Path, out_dir: Path) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    names = []
    if not zip_path.exists() or zip_path.stat().st_size < 20:
        return names
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for info in zf.infolist():
                if info.is_dir():
                    continue
                target = out_dir / Path(info.filename).name
                target.write_bytes(zf.read(info))
                names.append(str(target.relative_to(ROOT)))
    except zipfile.BadZipFile:
        pass
    return names


def main() -> int:
    RAW.mkdir(parents=True, exist_ok=True)
    EXTRACTED.mkdir(parents=True, exist_ok=True)
    INVENTORY.mkdir(parents=True, exist_ok=True)
    NOTES.mkdir(parents=True, exist_ok=True)

    all_recs: list[dict] = []

    # --- Bessede 2024 CCR figshare supplements (clinical tables; mIF images NOT public)
    bessede_ids = {
        25232284: "S1_OAK_POPLAR_characteristics",
        25232281: "S2_OAK_POPLAR_multivariate",
        25232278: "S3_BIP_mIHF_characteristics",
        25232275: "S4_BIP_plasma_proteomics_characteristics",
    }
    for aid, tag in bessede_ids.items():
        dest = RAW / "bessede_ccr2024" / tag
        all_recs.extend(figshare_article(aid, dest))

    # Bessede collection landing / data availability pages
    for url, name in [
        (
            "https://aacrjournals.org/clincancerres/article/30/4/779/734210/TROP2-Is-Associated-with-Primary-Resistance-to",
            "bessede_ccr2024_landing.html",
        ),
        (
            "https://aacrjournals.org/clincancerres/article-pdf/30/4/779/3359782/779.pdf",
            "bessede_ccr2024_article.pdf",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "bessede_ccr2024", name))

    # --- Hashimoto 2025 Sci Rep TROP2 IHC + Nivo-Ipi
    for url, name in [
        (
            "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-025-19362-3/MediaObjects/41598_2025_19362_MOESM1_ESM.docx",
            "hashimoto2025_sci_rep_MOESM1.docx",
        ),
        (
            "https://www.nature.com/articles/s41598-025-19362-3.pdf",
            "hashimoto2025_sci_rep_article.pdf",
        ),
        (
            "https://www.nature.com/articles/s41598-025-19362-3",
            "hashimoto2025_sci_rep_landing.html",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "hashimoto2025_scirep", name))

    # --- Dum 2022 Pathobiology TROP2 TMA 18,563 tumors (includes lung)
    all_recs.extend(pmc_oa_package("PMC9393818", RAW / "dum2022_trop2_tma"))
    for url, name in [
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9393818/pdf/pat-0089-0245.pdf",
            "dum2022_pathobiology_article.pdf",
        ),
        (
            "https://karger.com/pat/article-pdf/89/4/245/3668885/000522206.pdf",
            "dum2022_karger_article.pdf",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "dum2022_trop2_tma", name))

    # --- Omori 2021 J Cancer Res Clin Oncol TROP2 IHC lung ± ICI mention
    all_recs.extend(pmc_oa_package("PMC11801068", RAW / "omori2021_trop2_lung"))
    for url, name in [
        (
            "https://static-content.springer.com/esm/art%3A10.1007%2Fs00432-021-03784-3/MediaObjects/432_2021_3784_MOESM1_ESM.docx",
            "omori2021_MOESM1.docx",
        ),
        (
            "https://static-content.springer.com/esm/art%3A10.1007%2Fs00432-021-03784-3/MediaObjects/432_2021_3784_MOESM2_ESM.pdf",
            "omori2021_MOESM2.pdf",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "omori2021_trop2_lung", name))

    # --- Inamura 2017 Oncotarget TROP2 IHC lung subtypes
    for url, name in [
        (
            "https://www.oncotarget.com/article/15647/pdf/",
            "inamura2017_oncotarget.pdf",
        ),
        (
            "https://www.oncotarget.com/index.php?journal=oncotarget&page=article&op=downloadSuppFile&path%5B%5D=15647&path%5B%5D=22680",
            "inamura2017_supplement.pdf",
        ),
        (
            "https://www.oncotarget.com/article/15647/",
            "inamura2017_landing.html",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "inamura2017_oncotarget", name))

    # --- Pak 2012 WJSO EpCAM/TROP2 NSCLC (open)
    for url, name in [
        (
            "https://wjso.biomedcentral.com/counter/pdf/10.1186/1477-7819-10-53.pdf",
            "pak2012_wjso.pdf",
        ),
        (
            "https://wjso.biomedcentral.com/articles/10.1186/1477-7819-10-53",
            "pak2012_landing.html",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "pak2012_wjso", name))

    # --- PLOS ONE Kuo 2025 Trop-2 NSCLC (IHC H-scores exist; IPD request-only)
    for url, name in [
        (
            "https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0321555&type=printable",
            "kuo2025_plosone.pdf",
        ),
        (
            "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0321555",
            "kuo2025_landing.html",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s001",
            "kuo2025_s001",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s002",
            "kuo2025_s002",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s003",
            "kuo2025_s003",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s004",
            "kuo2025_s004",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s005",
            "kuo2025_s005",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s006",
            "kuo2025_s006",
        ),
        (
            "https://doi.org/10.1371/journal.pone.0321555.s007",
            "kuo2025_s007",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "kuo2025_plosone", name))

    # --- Frontiers review compilation table (TROP2 lung literature)
    for url, name in [
        (
            "https://www.frontiersin.org/articles/10.3389/fonc.2025.1638054/full",
            "frontiers2025_trop2_review.html",
        ),
        (
            "https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2025.1638054/full",
            "frontiers2025_trop2_review_alt.html",
        ),
        (
            "https://www.frontiersin.org/api/v2/articles/1638054/file/DataSheet1.XLSX",
            "frontiers2025_DataSheet1.xlsx",
        ),
        (
            "https://www.frontiersin.org/articles/10.3389/fonc.2025.1638054/full#supplementary-material",
            "frontiers2025_supp_anchor.html",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "frontiers2025_trop2_review", name))

    # --- POR 2023 CLDN H-scores in rare lung cancers (ACC/MEC)
    for url, name in [
        (
            "https://www.por-journal.com/journals/pathology-and-oncology-research/articles/10.3389/pore.2023.1611328/pdf",
            "gyulai2023_por_cldn_acc_mec.pdf",
        ),
        (
            "https://www.por-journal.com/journals/pathology-and-oncology-research/articles/10.3389/pore.2023.1611328/full",
            "gyulai2023_por_landing.html",
        ),
        (
            "https://www.frontiersin.org/articles/10.3389/pore.2023.1611328/full",
            "gyulai2023_frontiers_landing.html",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "gyulai2023_cldn_rare_lung", name))

    # --- POR 2016 CLDN1-7 NSCLC IHC scores (includes CLDN4)
    for url, name in [
        (
            "https://www.por-journal.com/articles/10.1007/s12253-016-0115-0/pdf",
            "moldvay2016_por_cldn_nsclc.pdf",
        ),
        (
            "https://link.springer.com/content/pdf/10.1007/s12253-016-0115-0.pdf",
            "moldvay2016_springer.pdf",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "moldvay2016_cldn_nsclc", name))

    # --- Jung 2014 KJTCVS CLDN4 lung ADC (open)
    for url, name in [
        (
            "https://www.jchestsurg.org/upload/pdf/kjtcs-47-262.pdf",
            "jung2014_kjtcs_cldn4.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4157471/pdf/kjtcs-47-262.pdf",
            "jung2014_pmc.pdf",
        ),
    ]:
        all_recs.append(download_named(url, RAW / "jung2014_cldn4_luad", name))

    # --- Human Protein Atlas TACSTD2 / CLDN4 IHC summaries
    hpa_dir = RAW / "hpa"
    for url, name in [
        ("https://www.proteinatlas.org/ENSG00000184292-TACSTD2.tsv", "HPA_TACSTD2.tsv"),
        ("https://www.proteinatlas.org/ENSG00000184292-TACSTD2.json", "HPA_TACSTD2.json"),
        ("https://www.proteinatlas.org/ENSG00000189143-CLDN4.tsv", "HPA_CLDN4.tsv"),
        ("https://www.proteinatlas.org/ENSG00000189143-CLDN4.json", "HPA_CLDN4.json"),
        ("https://www.proteinatlas.org/download/pathology.tsv.zip", "HPA_pathology.tsv.zip"),
        (
            "https://www.proteinatlas.org/search/gene_name:TACSTD2%20OR%20gene_name:CLDN4?format=tsv",
            "HPA_search_TACSTD2_CLDN4.tsv",
        ),
        (
            "https://www.proteinatlas.org/ENSG00000184292-TACSTD2/pathology/lung+cancer",
            "HPA_TACSTD2_lung_cancer.html",
        ),
        (
            "https://www.proteinatlas.org/ENSG00000189143-CLDN4/pathology/lung+cancer",
            "HPA_CLDN4_lung_cancer.html",
        ),
    ]:
        all_recs.append(download_named(url, hpa_dir, name))

    zip_path = hpa_dir / "HPA_pathology.tsv.zip"
    if zip_path.exists():
        extract_zip(zip_path, EXTRACTED / "hpa")

    # --- Zenodo ADC TROP2 H-score record (expected restricted)
    for url, name in [
        ("https://zenodo.org/api/records/18543127", "zenodo_18543127.json"),
        ("https://zenodo.org/api/records/18494664", "zenodo_18494664.json"),
        ("https://zenodo.org/records/18543127/files/cohort_data.csv?download=1", "zenodo_18543127_cohort_data.csv"),
        ("https://zenodo.org/records/18494664/files/cohort_data.csv?download=1", "zenodo_18494664_cohort_data.csv"),
    ]:
        all_recs.append(download_named(url, RAW / "zenodo_adc_trop2", name))

    # --- Figshare / EuropePMC discovery dumps
    searches = {
        "figshare_trop2_ihc": "https://api.figshare.com/v2/articles/search",
    }
    # figshare search is POST; skip POST and use GET search
    payload, meta = json_get(
        "https://api.figshare.com/v2/articles?search_for=TROP2%20IHC%20NSCLC&page_size=20"
    )
    write_json(INVENTORY / "figshare_search_trop2_ihc_nsclc.json", payload or {"error": meta})
    all_recs.append({**{k: v for k, v in meta.items() if k != "_data"}, "filename": "figshare_search_trop2_ihc_nsclc.json"})

    payload, meta = json_get(
        "https://api.figshare.com/v2/articles?search_for=CLDN4%20immunohistochemistry%20lung&page_size=20"
    )
    write_json(INVENTORY / "figshare_search_cldn4_lung.json", payload or {"error": meta})
    all_recs.append({**{k: v for k, v in meta.items() if k != "_data"}, "filename": "figshare_search_cldn4_lung.json"})

    payload, meta = json_get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=TROP2%20IHC%20NSCLC%20supplementary&resultType=lite&format=json&pageSize=25"
    )
    write_json(INVENTORY / "europepmc_trop2_ihc_nsclc.json", payload or {"error": meta})
    all_recs.append({**{k: v for k, v in meta.items() if k != "_data"}, "filename": "europepmc_trop2_ihc_nsclc.json"})

    payload, meta = json_get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=CLDN4%20IHC%20lung%20cancer%20H-score&resultType=lite&format=json&pageSize=25"
    )
    write_json(INVENTORY / "europepmc_cldn4_ihc_lung.json", payload or {"error": meta})
    all_recs.append({**{k: v for k, v in meta.items() if k != "_data"}, "filename": "europepmc_cldn4_ihc_lung.json"})

    # Clean recs for CSV
    rows = []
    for rec in all_recs:
        rec = {k: v for k, v in rec.items() if k != "_data"}
        rows.append(rec)

    inv_json = INVENTORY / "download_log.json"
    write_json(inv_json, {"generated_at": utc_now(), "n": len(rows), "records": rows})

    inv_csv = INVENTORY / "download_log.csv"
    fields = [
        "ok",
        "status",
        "bytes",
        "filename",
        "path",
        "url",
        "final_url",
        "content_type",
        "sha256",
        "error",
        "retrieved_at",
    ]
    with inv_csv.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)

    ok_n = sum(1 for r in rows if r.get("ok") and r.get("bytes", 0) > 200)
    fail_n = sum(1 for r in rows if not r.get("ok"))
    print(f"downloads logged: {len(rows)}  ok>200B: {ok_n}  fail: {fail_n}")
    print(f"inventory: {inv_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
