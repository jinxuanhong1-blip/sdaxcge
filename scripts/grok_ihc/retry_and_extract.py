#!/usr/bin/env python3
"""Retry failed open downloads and extract scoring / clinicopathologic tables."""
from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
from datetime import datetime, timezone
from html.parser import HTMLParser
from io import BytesIO
from pathlib import Path
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "results" / "grok_ihc" / "raw"
EXT = ROOT / "results" / "grok_ihc" / "extracted"
INV = ROOT / "results" / "grok_ihc" / "inventory"
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/126.0.0.0 Safari/537.36"


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def fetch(url: str, dest: Path, retries: int = 3) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for i in range(retries):
        try:
            req = Request(url, headers={"User-Agent": UA, "Accept": "*/*"})
            with urlopen(req, timeout=90) as resp:
                data = resp.read()
                dest.write_bytes(data)
                return {
                    "ok": True,
                    "url": url,
                    "final_url": getattr(resp, "url", url),
                    "status": getattr(resp, "status", 200),
                    "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "path": str(dest.relative_to(ROOT)),
                    "content_type": resp.headers.get("Content-Type", ""),
                    "error": "",
                    "retrieved_at": utc_now(),
                }
        except (HTTPError, URLError, TimeoutError, OSError) as exc:
            last = exc
    return {
        "ok": False,
        "url": url,
        "path": str(dest.relative_to(ROOT)),
        "bytes": 0,
        "error": f"{type(last).__name__}: {last}",
        "retrieved_at": utc_now(),
    }


def looks_like_html(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 20:
        return False
    head = path.read_bytes()[:200].lstrip().lower()
    return head.startswith(b"<!doctype") or head.startswith(b"<html")


def docx_tables(path: Path) -> list[list[list[str]]]:
    tables = []
    with zipfile.ZipFile(path) as zf:
        xml = zf.read("word/document.xml")
    root = ET.fromstring(xml)
    ns = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
    for tbl in root.findall(".//w:tbl", ns):
        rows = []
        for tr in tbl.findall("w:tr", ns):
            cells = []
            for tc in tr.findall("w:tc", ns):
                texts = [t.text or "" for t in tc.findall(".//w:t", ns)]
                cells.append("".join(texts).strip())
            if any(cells):
                rows.append(cells)
        if rows:
            tables.append(rows)
    return tables


def write_tables_csv(tables: list[list[list[str]]], dest_dir: Path, prefix: str) -> list[str]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    out = []
    for i, tbl in enumerate(tables, 1):
        p = dest_dir / f"{prefix}_table{i:02d}.csv"
        with p.open("w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerows(tbl)
        out.append(str(p.relative_to(ROOT)))
    return out


def html_tables(path: Path) -> list[list[list[str]]]:
    class T(HTMLParser):
        def __init__(self):
            super().__init__()
            self.tables = []
            self.cur = None
            self.row = None
            self.cell = None
            self.buf = []

        def handle_starttag(self, tag, attrs):
            if tag == "table":
                self.cur = []
            elif tag == "tr" and self.cur is not None:
                self.row = []
            elif tag in {"td", "th"} and self.row is not None:
                self.cell = True
                self.buf = []

        def handle_endtag(self, tag):
            if tag in {"td", "th"} and self.cell:
                self.row.append(re.sub(r"\s+", " ", "".join(self.buf)).strip())
                self.cell = False
            elif tag == "tr" and self.row is not None:
                if any(self.row):
                    self.cur.append(self.row)
                self.row = None
            elif tag == "table" and self.cur is not None:
                if self.cur:
                    self.tables.append(self.cur)
                self.cur = None

        def handle_data(self, data):
            if self.cell:
                self.buf.append(data)

    p = T()
    p.feed(path.read_text(encoding="utf-8", errors="replace"))
    return p.tables


def extract_hpa_cancer_ihc(html_path: Path) -> list[dict]:
    text = html_path.read_text(encoding="utf-8", errors="replace")
    # HPA cancer pages embed counts like High/Medium/Low/Not detected
    rows = []
    # Try JSON blobs
    for m in re.finditer(r"cancer.*?lung", text, re.I):
        pass
    # Antibody staining summary sentences
    for m in re.finditer(
        r"(lung cancer|Lung cancer)[^.]{0,400}",
        text,
        re.I,
    ):
        rows.append({"snippet": re.sub(r"\s+", " ", m.group(0))[:500]})
    # table-like High Medium Low Not detected
    for m in re.finditer(
        r"(High|high)[^0-9]{0,40}(\d+)[^0-9]{0,40}(Medium|medium)[^0-9]{0,40}(\d+)[^0-9]{0,40}(Low|low)[^0-9]{0,40}(\d+)[^0-9]{0,40}(Not detected|not detected)[^0-9]{0,40}(\d+)",
        text,
    ):
        rows.append(
            {
                "high": m.group(2),
                "medium": m.group(4),
                "low": m.group(6),
                "not_detected": m.group(8),
                "context": re.sub(r"\s+", " ", m.group(0))[:300],
            }
        )
    return rows


def main() -> None:
    EXT.mkdir(parents=True, exist_ok=True)
    INV.mkdir(parents=True, exist_ok=True)
    retries = []

    extra = [
        (
            "https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC9393818&blobtype=pdf",
            RAW / "dum2022_trop2_tma" / "dum2022_europepmc.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9393818/",
            RAW / "dum2022_trop2_tma" / "dum2022_pmc.html",
        ),
        (
            "https://europepmc.org/articles/PMC9393818?pdf=render",
            RAW / "dum2022_trop2_tma" / "dum2022_europepmc_render.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410292/pdf/oncotarget-08-28725.pdf",
            RAW / "inamura2017_oncotarget" / "inamura2017_pmc.pdf",
        ),
        (
            "https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC5410292&blobtype=pdf",
            RAW / "inamura2017_oncotarget" / "inamura2017_europepmc.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5410292/",
            RAW / "inamura2017_oncotarget" / "inamura2017_pmc.html",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4157471/",
            RAW / "jung2014_cldn4_luad" / "jung2014_pmc.html",
        ),
        (
            "https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC4157471&blobtype=pdf",
            RAW / "jung2014_cldn4_luad" / "jung2014_europepmc.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12411510/bin/DataSheet1.xlsx",
            RAW / "frontiers2025_trop2_review" / "DataSheet1.xlsx",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC12411510/",
            RAW / "frontiers2025_trop2_review" / "frontiers2025_pmc.html",
        ),
        (
            "https://v23.proteinatlas.org/download/pathology.tsv.zip",
            RAW / "hpa" / "HPA_v23_pathology.tsv.zip",
        ),
        (
            "https://v19.proteinatlas.org/download/pathology.tsv.zip",
            RAW / "hpa" / "HPA_v19_pathology.tsv.zip",
        ),
        (
            "https://www.proteinatlas.org/ENSG00000184292-TACSTD2.xml",
            RAW / "hpa" / "HPA_TACSTD2.xml",
        ),
        (
            "https://www.proteinatlas.org/ENSG00000189143-CLDN4.xml",
            RAW / "hpa" / "HPA_CLDN4.xml",
        ),
        (
            "https://v23.proteinatlas.org/ENSG00000184292-TACSTD2.xml",
            RAW / "hpa" / "HPA_v23_TACSTD2.xml",
        ),
        (
            "https://v23.proteinatlas.org/ENSG00000189143-CLDN4.xml",
            RAW / "hpa" / "HPA_v23_CLDN4.xml",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7463615/pdf/cancers-12-03328.pdf",
            RAW / "mito2020_cancers" / "mito2020_cancers.pdf",
        ),
        (
            "https://www.mdpi.com/2072-6694/12/11/3328/pdf",
            RAW / "mito2020_cancers" / "mito2020_mdpi.pdf",
        ),
        (
            "https://www.mdpi.com/2072-6694/12/11/3328",
            RAW / "mito2020_cancers" / "mito2020_landing.html",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11801068/",
            RAW / "omori2021_trop2_lung" / "omori2021_pmc.html",
        ),
        (
            "https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC11801068&blobtype=pdf",
            RAW / "omori2021_trop2_lung" / "omori2021_europepmc.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11999141/",
            RAW / "kuo2025_plosone" / "kuo2025_pmc.html",
        ),
        (
            "https://pdfs.semanticscholar.org/b6e2/48fd2fb95f89ceaceabc278cf70b8e3753d5.pdf",
            RAW / "bessede_ccr2024" / "bessede_ccr2024_semantic.pdf",
        ),
        (
            "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10870750/",
            RAW / "bessede_ccr2024" / "bessede_pmc.html",
        ),
        (
            "https://europepmc.org/backend/ptpmcrender.fcgi?accid=PMC10870750&blobtype=pdf",
            RAW / "bessede_ccr2024" / "bessede_europepmc.pdf",
        ),
    ]
    for url, dest in extra:
        rec = fetch(url, dest)
        retries.append(rec)
        print(f"{'OK' if rec['ok'] else 'FAIL'} {rec.get('bytes', 0):8d} {dest.name} {rec.get('error','')}")

    # unzip HPA pathology
    for zp in [
        RAW / "hpa" / "HPA_v23_pathology.tsv.zip",
        RAW / "hpa" / "HPA_v19_pathology.tsv.zip",
    ]:
        if zp.exists() and not looks_like_html(zp) and zp.stat().st_size > 1000:
            out = EXT / "hpa"
            out.mkdir(parents=True, exist_ok=True)
            try:
                with zipfile.ZipFile(zp) as zf:
                    zf.extractall(out)
            except zipfile.BadZipFile:
                pass

    extracted = []

    # Bessede docx
    for docx in (RAW / "bessede_ccr2024").rglob("*.docx"):
        tables = docx_tables(docx)
        paths = write_tables_csv(tables, EXT / "bessede_ccr2024", docx.stem)
        extracted.append({"source": str(docx.relative_to(ROOT)), "n_tables": len(tables), "csvs": paths})

    # Hashimoto / Omori docx
    for docx in list((RAW / "hashimoto2025_scirep").glob("*.docx")) + list(
        (RAW / "omori2021_trop2_lung").glob("*.docx")
    ):
        tables = docx_tables(docx)
        paths = write_tables_csv(tables, EXT / docx.parent.name, docx.stem)
        extracted.append({"source": str(docx.relative_to(ROOT)), "n_tables": len(tables), "csvs": paths})

    # HTML tables
    html_targets = [
        RAW / "hashimoto2025_scirep" / "hashimoto2025_sci_rep_landing.html",
        RAW / "pak2012_wjso" / "pak2012_landing.html",
        RAW / "gyulai2023_cldn_rare_lung" / "gyulai2023_por_landing.html",
        RAW / "kuo2025_plosone" / "kuo2025_landing.html",
        RAW / "kuo2025_plosone" / "kuo2025_pmc.html",
        RAW / "dum2022_trop2_tma" / "dum2022_pmc.html",
        RAW / "inamura2017_oncotarget" / "inamura2017_pmc.html",
        RAW / "jung2014_cldn4_luad" / "jung2014_pmc.html",
        RAW / "mito2020_cancers" / "mito2020_landing.html",
        RAW / "omori2021_trop2_lung" / "omori2021_pmc.html",
        RAW / "bessede_ccr2024" / "bessede_pmc.html",
        RAW / "frontiers2025_trop2_review" / "frontiers2025_pmc.html",
    ]
    for hp in html_targets:
        if not hp.exists() or looks_like_html(hp) is False and hp.stat().st_size < 500:
            continue
        if not hp.exists():
            continue
        tables = html_tables(hp)
        # keep reasonably sized tables
        tables = [t for t in tables if len(t) >= 2 and max(len(r) for r in t) >= 2]
        paths = write_tables_csv(tables, EXT / hp.parent.name / "html_tables", hp.stem)
        extracted.append({"source": str(hp.relative_to(ROOT)), "n_tables": len(tables), "csvs": paths})

    # HPA XML cancer IHC
    hpa_rows = []
    for xmlp, gene in [
        (RAW / "hpa" / "HPA_v23_TACSTD2.xml", "TACSTD2"),
        (RAW / "hpa" / "HPA_v23_CLDN4.xml", "CLDN4"),
        (RAW / "hpa" / "HPA_TACSTD2.xml", "TACSTD2"),
        (RAW / "hpa" / "HPA_CLDN4.xml", "CLDN4"),
    ]:
        if not xmlp.exists() or looks_like_html(xmlp):
            continue
        try:
            root = ET.parse(xmlp).getroot()
        except ET.ParseError:
            continue
        for tissue in root.findall(".//tissue"):
            name = tissue.attrib.get("organ") or tissue.attrib.get("name") or ""
            if "lung" not in name.lower() and tissue.attrib.get("ontologyTerms", "").find("lung") < 0:
                # still collect cancer tissues
                pass
            level = tissue.attrib.get("level") or tissue.attrib.get("staining")
            if level or "cancer" in (tissue.attrib.get("category", "") + name).lower():
                hpa_rows.append({"gene": gene, "tissue": name, **tissue.attrib})
        for cancer in root.findall(".//cancer"):
            hpa_rows.append({"gene": gene, "tag": "cancer", **cancer.attrib, "text": (cancer.text or "")[:200]})
        for patient in root.findall(".//patient"):
            organ = patient.attrib.get("sex", "")
            # patient cancer IHC
            hpa_rows.append({"gene": gene, "tag": "patient", **patient.attrib})

    if hpa_rows:
        p = EXT / "hpa" / "hpa_xml_cancer_ihc.csv"
        p.parent.mkdir(parents=True, exist_ok=True)
        keys = sorted({k for r in hpa_rows for k in r})
        with p.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(hpa_rows)
        extracted.append({"source": "hpa_xml", "n_rows": len(hpa_rows), "csvs": [str(p.relative_to(ROOT))]})

    # Filter HPA pathology.tsv for TACSTD2 / CLDN4 lung
    for tsv in (EXT / "hpa").glob("*.tsv"):
        rows = []
        with tsv.open(encoding="utf-8", errors="replace") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            for row in r:
                g = (row.get("Gene name") or row.get("Gene") or "").upper()
                c = (row.get("Cancer") or row.get("cancer") or "").lower()
                if g in {"TACSTD2", "CLDN4", "TROP2"} and ("lung" in c or c == ""):
                    rows.append(row)
                elif g in {"TACSTD2", "CLDN4"}:
                    rows.append(row)
        if rows:
            outp = EXT / "hpa" / f"{tsv.stem}_TACSTD2_CLDN4.csv"
            with outp.open("w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
                w.writeheader()
                w.writerows(rows)
            extracted.append({"source": str(tsv.relative_to(ROOT)), "n_rows": len(rows), "csvs": [str(outp.relative_to(ROOT))]})

    # HPA HTML snippets
    for hp, gene in [
        (RAW / "hpa" / "HPA_TACSTD2_lung_cancer.html", "TACSTD2"),
        (RAW / "hpa" / "HPA_CLDN4_lung_cancer.html", "CLDN4"),
    ]:
        if hp.exists():
            snippets = extract_hpa_cancer_ihc(hp)
            outp = EXT / "hpa" / f"{gene}_lung_html_snippets.json"
            outp.parent.mkdir(parents=True, exist_ok=True)
            outp.write_text(json.dumps(snippets[:50], indent=2, ensure_ascii=False), encoding="utf-8")

    # Copy HPA search TSV
    src = RAW / "hpa" / "HPA_search_TACSTD2_CLDN4.tsv"
    if src.exists() and src.stat().st_size > 100:
        dest = EXT / "hpa" / "HPA_search_TACSTD2_CLDN4.tsv"
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(src.read_bytes())

    summary = {
        "generated_at": utc_now(),
        "retries": retries,
        "extracted": extracted,
    }
    (INV / "extract_log.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"extracted sources: {len(extracted)}")


if __name__ == "__main__":
    main()
