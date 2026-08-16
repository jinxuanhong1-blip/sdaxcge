#!/usr/bin/env python3
"""Hunt NGDC/GSA for truly OPEN processed Chinese lung ICI data.

Scope (honest):
  - OMIX Open-access records (processed tables; HTTPS, no DAC).
  - GSA-Human records flagged Open (usually raw FASTQ, not processed).
  - Do NOT request or download HRA Controlled / DAC data.
  - Question: is TACSTD2 and/or CLDN4 measurable against MPR (or any ICI
    response label) in an open processed lung cohort?

Outputs land in results/hunt_gsa_open/.
"""

from __future__ import annotations

import csv
import html
import json
import os
import re
import sys
import zipfile
from collections import Counter
from typing import Iterable

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from ngdc_client import download, fetch, fetch_many, setup_logging  # noqa: E402

OUT = os.path.join(ROOT, "results", "hunt_gsa_open")
CAT = os.path.join(ROOT, "data", "catalogs")
RAW = os.path.join(ROOT, "data", "raw")

LUNG_RE = re.compile(
    r"\b(lung|nsclc|sclc|luad|lusc|pulmonary|bronchogenic)\b", re.I
)
ICI_RE = re.compile(
    r"(immunotherap|immune checkpoint|checkpoint inhibit|anti-pd|"
    r"\bpd-?1\b|\bpd-?l1\b|ctla-?4|pembrolizumab|nivolumab|ipilimumab|"
    r"camrelizumab|sintilimab|tislelizumab|toripalimab|durvalumab|"
    r"atezolizumab|adebrelimab|sugemalimab|penpulimab|serplulimab|"
    r"cemiplimab|envafolimab|cadonilimab|ak112|ivonescimab|"
    r"chemoimmuno|chemo-immuno|immuno-chemo|icb\b|ici\b)",
    re.I,
)
NEOADJ_RE = re.compile(r"(neoadjuvant|perioperative|resectable)", re.I)
MPR_RE = re.compile(
    r"\b(mpr|pcr|mpr/pcr|pathologic(?:al)? (?:complete )?response|"
    r"major patholog|residual viable tumor|%rvt|ypT)\b",
    re.I,
)
GENE_TITLE_RE = re.compile(
    r"(TACSTD2|TROP-?2|CLDN4|claudin[- ]?4|GA733-1|EGP-?1)", re.I
)
# processed-looking OMIX data types / file titles
PROCESSED_RE = re.compile(
    r"(expression|rna-?seq|transcript|nanostring|count|tpm|fpkm|rpkm|"
    r"matrix|proteom|metabolom|biomarker|clinical|phenotype|response|"
    r"xlsx|xls|csv|tsv|txt)",
    re.I,
)
RAWSEQ_RE = re.compile(r"\.(fq|fastq|bam|cram|sra)(\.gz)?$", re.I)

GENE_ALIASES = {
    "TACSTD2": [
        "TACSTD2",
        "TROP2",
        "TROP-2",
        "GA733-1",
        "GA7331",
        "EGP-1",
        "EGP1",
        "M1S1",
        "ENSG00000184292",
    ],
    "CLDN4": [
        "CLDN4",
        "CLAUDIN4",
        "CLAUDIN-4",
        "CLAUDIN 4",
        "CPE-R",
        "CPER",
        "ENSG00000189143",
    ],
}

MAX_DOWNLOAD_BYTES = 40 * 1024 * 1024  # 40 MB — processed tables, not WES/FASTQ


def read_tsv(path: str) -> list[dict]:
    with open(path, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: str, rows: list[dict], fields: list[str] | None = None) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if not rows and not fields:
        open(path, "w").close()
        return
    fields = fields or list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def hay(row: dict, extra: str = "") -> str:
    return " ".join(
        [
            row.get("title") or "",
            row.get("organism") or "",
            row.get("bioproject") or "",
            extra,
        ]
    )


def flags(text: str) -> dict:
    return {
        "lung": bool(LUNG_RE.search(text)),
        "ici": bool(ICI_RE.search(text)),
        "neoadjuvant": bool(NEOADJ_RE.search(text)),
        "mpr_keyword": bool(MPR_RE.search(text)),
        "gene_keyword": bool(GENE_TITLE_RE.search(text)),
    }


def strip_html(raw: str) -> str:
    t = re.sub(r"<script.*?</script>", " ", raw, flags=re.S)
    t = re.sub(r"<style.*?</style>", " ", t, flags=re.S)
    t = html.unescape(re.sub(r"<[^>]+>", " ", t))
    t = re.sub(r"[ \t]+", " ", t)
    t = re.sub(r"\n\s*\n+", "\n", t)
    return t


def parse_omix_page(raw: str) -> dict:
    """Pull summary fields + file table from an OMIX release page."""
    text = strip_html(raw)
    def grab(label: str) -> str:
        m = re.search(rf"{label}\s+(.*?)(?:\n|$)", text)
        # labels sit as their own tokens; value is the following non-empty run
        m = re.search(
            rf"\b{label}\s+(.+?)(?=\s+(?:Title|Description|Organism|Data Type|"
            rf"Data Accessibility|BioProject|Release Date|Submitter|"
            rf"Organization|Submission Date|Files? & Download|File ID|Paper Title)\b)",
            text,
            re.S,
        )
        return (m.group(1).strip() if m else "")

    # file download links
    files = []
    for url in set(re.findall(r"https://download\.cncb\.ac\.cn/OMIX/[^\"'\s]+", raw)):
        files.append(
            {
                "url": url.rstrip("/"),
                "name": url.rstrip("/").split("/")[-1],
            }
        )
    accs = sorted(set(re.findall(r"(?:OMIX|HRA|CRA|GSE|PRJCA)\d+", raw)))
    pubs = []
    for m in re.finditer(
        r"<tr>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>\s*<td[^>]*>(.*?)</td>",
        raw,
        re.S,
    ):
        a, b, c = (_text_cell(m.group(i)) for i in (1, 2, 3))
        if a and a not in ("Paper Title", "File ID", "OMIX ID") and len(a) > 8:
            pubs.append({"title": a, "journal": b, "year": c})

    return {
        "title": grab("Title"),
        "description": grab("Description"),
        "organism": grab("Organism"),
        "data_type": grab("Data Type"),
        "accessibility": grab("Data Accessibility"),
        "bioproject": grab("BioProject"),
        "release_date": grab("Release Date"),
        "submitter": grab("Submitter"),
        "organization": grab("Organization"),
        "files": files,
        "related_accessions": accs,
        "publications": pubs[:8],
        "text": text[:20000],
    }


def _text_cell(s: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", " ", s)).strip()


def parse_hra_page(raw: str) -> dict:
    text = strip_html(raw)
    def grab(label: str) -> str:
        m = re.search(rf"{re.escape(label)}\s*:\s*(.+)", text)
        return m.group(1).strip() if m else ""

    files = []
    for url in set(
        re.findall(r"https://download\.cncb\.ac\.cn/gsa-human/HRA\d+", raw)
    ):
        files.append({"url": url, "name": url.rstrip("/").split("/")[-1]})
    acc = ""
    m = re.search(r"\b(HRA\d+)\b", text)
    if m:
        acc = m.group(1)
    return {
        "accession": acc,
        "title": grab("Title"),
        "description": grab("Description")[:4000],
        "accessibility": grab("Data Accessibility") or (
            "Open access" if re.search(r"Open access", text) else ""
        ),
        "bioproject": grab("BioProject"),
        "study_type": grab("Study type"),
        "release_date": grab("Release date"),
        "organization": grab("Organization"),
        "files": files,
        "text": text[:15000],
    }


def size_from_omix_html(raw: str, file_id: str) -> str:
    # "2.9 MB" near the file id
    m = re.search(
        rf"{re.escape(file_id)}.*?(?:File Size|)\s*([0-9.]+)\s*(KB|MB|GB|B)",
        strip_html(raw),
        re.S | re.I,
    )
    return m.group(0)[-20:] if m else ""


def parse_size_token(raw: str, filename: str) -> tuple[float, str]:
    """Return (bytes_estimate, display) from nearby '2.9 MB' text."""
    t = strip_html(raw)
    stem = filename.rsplit(".", 1)[0]
    m = re.search(
        rf"{re.escape(stem)}.*?([0-9.]+)\s*(KB|MB|GB|B)\b", t, re.S | re.I
    )
    if not m:
        return 0.0, ""
    n = float(m.group(1))
    unit = m.group(2).upper()
    mult = {"B": 1, "KB": 1024, "MB": 1024**2, "GB": 1024**3}[unit]
    return n * mult, f"{n} {unit}"


def shortlist_omix(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        if r["access"] != "Open":
            continue
        f = flags(hay(r))
        # keep anything lung, ICI, gene-named, or MPR-named
        if not (f["lung"] or f["ici"] or f["gene_keyword"] or f["mpr_keyword"]):
            continue
        rec = dict(r)
        rec.update(f)
        rec["why"] = "+".join(k for k, v in f.items() if v)
        out.append(rec)
    return out


def shortlist_hra(rows: list[dict]) -> list[dict]:
    out = []
    for r in rows:
        f = flags(hay(r))
        interesting = f["lung"] or f["ici"] or f["gene_keyword"] or f["mpr_keyword"]
        if not interesting:
            continue
        rec = dict(r)
        rec.update(f)
        rec["why"] = "+".join(k for k, v in f.items() if v)
        out.append(rec)
    return out


def scan_text_for_genes(text: str) -> dict:
    hits = {}
    up = text.upper()
    for gene, aliases in GENE_ALIASES.items():
        found = [a for a in aliases if a.upper() in up]
        hits[gene] = found
    return hits


def iter_table_files(path: str):
    """Yield (source_name, header_or_first_col_chunk) from xlsx/xls/csv/tsv/txt/zip."""
    lower = path.lower()
    if zipfile.is_zipfile(path) and not lower.endswith((".xlsx", ".xls")):
        try:
            with zipfile.ZipFile(path) as zf:
                for info in zf.infolist():
                    if info.is_dir() or info.file_size > MAX_DOWNLOAD_BYTES:
                        continue
                    name = info.filename
                    if not re.search(r"\.(csv|tsv|txt|xlsx|xls)$", name, re.I):
                        continue
                    dest = path + ".__unz__" + os.path.basename(name)
                    if not os.path.exists(dest):
                        with zf.open(info) as src, open(dest, "wb") as dst:
                            dst.write(src.read())
                    yield from iter_table_files(dest)
        except zipfile.BadZipFile:
            return
        return

    if lower.endswith((".csv", ".tsv", ".txt")):
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                yield path, fh.read(2_000_000)
        except OSError:
            return
        return

    if lower.endswith((".xlsx", ".xls")):
        try:
            import openpyxl

            wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
            for ws in wb.worksheets:
                buf = []
                for i, row in enumerate(ws.iter_rows(values_only=True)):
                    buf.append("\t".join("" if c is None else str(c) for c in row))
                    if i >= 4000:
                        break
                yield f"{path}::{ws.title}", "\n".join(buf)
            wb.close()
        except Exception:
            # xls (old) — try xlrd
            try:
                import xlrd

                book = xlrd.open_workbook(path)
                for sh in book.sheets():
                    buf = []
                    for i in range(min(sh.nrows, 4000)):
                        buf.append(
                            "\t".join(str(sh.cell_value(i, j)) for j in range(sh.ncols))
                        )
                    yield f"{path}::{sh.name}", "\n".join(buf)
            except Exception:
                return


def scan_file_for_targets(path: str) -> dict:
    gene_hits = {g: [] for g in GENE_ALIASES}
    pheno_hits = []
    n_tables = 0
    for src, chunk in iter_table_files(path):
        n_tables += 1
        gh = scan_text_for_genes(chunk)
        for g, aliases in gh.items():
            for a in aliases:
                gene_hits[g].append(f"{os.path.basename(src)}:{a}")
        if MPR_RE.search(chunk) or re.search(
            r"\b(response|responder|PFS|OS|RECIST|DCR|ORR|durable)\b", chunk, re.I
        ):
            # keep a short evidence snippet
            m = MPR_RE.search(chunk) or re.search(
                r".{0,40}\b(response|responder|RECIST|PFS|ORR).{0,40}", chunk, re.I
            )
            if m:
                pheno_hits.append(m.group(0).replace("\n", " ")[:160])
    return {
        "n_tables": n_tables,
        "TACSTD2": sorted(set(gene_hits["TACSTD2"])),
        "CLDN4": sorted(set(gene_hits["CLDN4"])),
        "phenotype_snippets": pheno_hits[:8],
    }


def list_https_index(url: str) -> list[dict]:
    r = fetch(url if url.endswith("/") else url + "/")
    if not r["ok"]:
        return []
    out = []
    for name in re.findall(r'href="([^"?/][^"]*)"', r["text"]):
        if name in ("../",) or name.startswith("?"):
            continue
        out.append({"name": name.rstrip("/"), "url": url.rstrip("/") + "/" + name})
    return out


def main() -> None:
    setup_logging()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(RAW, exist_ok=True)

    omix = read_tsv(os.path.join(CAT, "omix_release.tsv"))
    hra = read_tsv(os.path.join(CAT, "hra_finished.tsv"))

    # --- catalog-level honesty stats ---
    stats = {
        "omix_total": len(omix),
        "omix_open": sum(1 for r in omix if r["access"] == "Open"),
        "omix_controlled": sum(1 for r in omix if r["access"] == "Controlled"),
        "hra_total": len(hra),
        "hra_open": sum(1 for r in hra if r["access"] == "Open"),
        "hra_controlled": sum(1 for r in hra if r["access"] == "Controlled"),
    }
    omix_lung = [r for r in omix if LUNG_RE.search(r["title"] or "")]
    omix_lung_ici = [r for r in omix_lung if ICI_RE.search(r["title"] or "")]
    stats["omix_lung_title"] = len(omix_lung)
    stats["omix_lung_open"] = sum(1 for r in omix_lung if r["access"] == "Open")
    stats["omix_lung_ici_title"] = len(omix_lung_ici)
    stats["omix_lung_ici_open"] = sum(1 for r in omix_lung_ici if r["access"] == "Open")
    hra_lung_ici = [
        r
        for r in hra
        if LUNG_RE.search(r["title"] or "") and ICI_RE.search(r["title"] or "")
    ]
    stats["hra_lung_ici_title"] = len(hra_lung_ici)
    stats["hra_lung_ici_open"] = sum(1 for r in hra_lung_ici if r["access"] == "Open")

    sl_omix = shortlist_omix(omix)
    sl_hra = shortlist_hra(hra)
    write_tsv(os.path.join(OUT, "omix_open_shortlist.tsv"), sl_omix)
    write_tsv(os.path.join(OUT, "hra_all_lung_ici_gene.tsv"), sl_hra)

    # Fetch ALL open shortlist OMIX pages (title-level filter is leaky).
    print(f"Fetching {len(sl_omix)} open OMIX detail pages...")
    omix_urls = [r["url"] for r in sl_omix]
    omix_pages = fetch_many(omix_urls, workers=6, progress_every=50)

    # Fetch HRA pages only for OPEN interesting records + all lung+ICI (to document control).
    hra_fetch = [
        r
        for r in sl_hra
        if r["access"] == "Open" or (r.get("lung") and r.get("ici"))
    ]
    print(f"Fetching {len(hra_fetch)} HRA detail pages (open interesting + all lung+ICI)...")
    hra_pages = fetch_many([r["url"] for r in hra_fetch], workers=4, progress_every=20)

    # Parse OMIX pages; upgrade flags from description.
    omix_details = []
    download_jobs = []
    for r in sl_omix:
        page = omix_pages.get(r["url"]) or {}
        parsed = parse_omix_page(page.get("text") or "") if page.get("ok") else {}
        blob = " ".join(
            [
                r.get("title") or "",
                parsed.get("description") or "",
                parsed.get("data_type") or "",
                parsed.get("title") or "",
            ]
        )
        f = flags(blob)
        rec = {
            **r,
            "page_ok": bool(page.get("ok")),
            "description": (parsed.get("description") or "")[:1500],
            "data_type": parsed.get("data_type") or "",
            "accessibility_page": parsed.get("accessibility") or "",
            "organization": parsed.get("organization") or "",
            "n_files": len(parsed.get("files") or []),
            "file_names": ";".join(x["name"] for x in parsed.get("files") or []),
            "related_accessions": ";".join(parsed.get("related_accessions") or []),
            "publications": " | ".join(
                p.get("title", "") for p in (parsed.get("publications") or [])
            )[:500],
        }
        rec.update({f"desc_{k}": v for k, v in f.items()})
        rec["priority"] = int(
            f["lung"]
            and (f["ici"] or f["neoadjuvant"] or f["mpr_keyword"] or f["gene_keyword"])
        )
        omix_details.append(rec)

        if not (f["lung"] or f["gene_keyword"]):
            continue
        if not (
            f["ici"]
            or f["neoadjuvant"]
            or f["mpr_keyword"]
            or f["gene_keyword"]
            or PROCESSED_RE.search(blob)
        ):
            # lung-only open records: still note, don't download everything
            if not (f["ici"] or f["gene_keyword"] or f["mpr_keyword"]):
                continue

        for fi in parsed.get("files") or []:
            nbytes, disp = parse_size_token(page.get("text") or "", fi["name"])
            if RAWSEQ_RE.search(fi["name"]):
                continue
            if nbytes and nbytes > MAX_DOWNLOAD_BYTES:
                rec.setdefault("skipped_large", "")
                rec["skipped_large"] = (rec.get("skipped_large") or "") + f"{fi['name']}({disp});"
                continue
            dest = os.path.join(RAW, "omix", r["accession"], fi["name"])
            download_jobs.append(
                {
                    "accession": r["accession"],
                    "url": fi["url"],
                    "dest": dest,
                    "size_disp": disp,
                    "priority": rec["priority"],
                    "why": rec.get("why"),
                }
            )

    write_tsv(os.path.join(OUT, "omix_open_details.tsv"), omix_details)

    # HRA details
    hra_details = []
    for r in hra_fetch:
        page = hra_pages.get(r["url"]) or {}
        parsed = parse_hra_page(page.get("text") or "") if page.get("ok") else {}
        listing = []
        if r["access"] == "Open" and (r.get("lung") or r.get("ici") or r.get("gene_keyword")):
            listing = list_https_index(
                f"https://download.cncb.ac.cn/gsa-human/{r['accession']}"
            )
        file_kinds = Counter()
        for item in listing:
            n = item["name"]
            if n.startswith("HRR") or n.startswith("HRA"):
                file_kinds["run_or_nested"] += 1
            elif n.endswith(".fq.gz") or n.endswith(".fastq.gz"):
                file_kinds["fastq"] += 1
            elif n in ("md5sum.txt",):
                file_kinds["md5"] += 1
            else:
                file_kinds["other"] += 1
        # peek first HRR child if present
        example = ""
        for item in listing:
            if item["name"].startswith("HRR"):
                child = list_https_index(item["url"])
                example = ",".join(c["name"] for c in child[:6])
                for c in child:
                    if c["name"].endswith(".fq.gz") or c["name"].endswith(".fastq.gz"):
                        file_kinds["fastq"] += 1
                break
        rec = {
            **r,
            "page_ok": bool(page.get("ok")),
            "description": (parsed.get("description") or "")[:1200],
            "accessibility_page": parsed.get("accessibility") or "",
            "study_type": parsed.get("study_type") or "",
            "release_date_page": parsed.get("release_date") or "",
            "n_listing": len(listing),
            "listing_kinds": json.dumps(file_kinds),
            "example_run_files": example,
            "processed_on_hra": int(
                file_kinds.get("other", 0) > 0 and file_kinds.get("fastq", 0) == 0
            ),
        }
        hra_details.append(rec)
    write_tsv(os.path.join(OUT, "hra_details.tsv"), hra_details)

    # Download open processed OMIX files (priority first).
    download_jobs.sort(key=lambda x: -x["priority"])
    # de-dup
    seen = set()
    uniq = []
    for j in download_jobs:
        if j["url"] in seen:
            continue
        seen.add(j["url"])
        uniq.append(j)
    download_jobs = uniq
    print(f"Downloading {len(download_jobs)} open OMIX files (cap {MAX_DOWNLOAD_BYTES} B)...")
    dl_log = []
    for j in download_jobs:
        ok = download(j["url"], j["dest"])
        size = os.path.getsize(j["dest"]) if ok and os.path.exists(j["dest"]) else 0
        # refuse if we accidentally pulled something huge
        if size > MAX_DOWNLOAD_BYTES * 1.5:
            os.remove(j["dest"])
            ok = False
            size = 0
            note = "deleted_oversize"
        else:
            note = "ok" if ok else "fail"
        dl_log.append({**j, "ok": int(ok), "bytes": size, "note": note})
        print(f"  {j['accession']} {j['url'].split('/')[-1]} {note} {size}")
    write_tsv(os.path.join(OUT, "downloads.tsv"), dl_log)

    # Scan downloaded tables
    scans = []
    for j in dl_log:
        if not j["ok"]:
            continue
        sc = scan_file_for_targets(j["dest"])
        scans.append(
            {
                "accession": j["accession"],
                "file": os.path.basename(j["dest"]),
                "bytes": j["bytes"],
                "n_tables": sc["n_tables"],
                "TACSTD2": ";".join(sc["TACSTD2"]),
                "CLDN4": ";".join(sc["CLDN4"]),
                "phenotype_snippets": " || ".join(sc["phenotype_snippets"]),
            }
        )
    write_tsv(os.path.join(OUT, "gene_phenotype_scans.tsv"), scans)

    # BioProject pages for high-priority accessions
    prjs = sorted(
        {
            r["bioproject"]
            for r in omix_details
            if r.get("desc_lung") and (r.get("desc_ici") or r.get("desc_gene_keyword"))
        }
        | {
            r["bioproject"]
            for r in hra_details
            if r["access"] == "Open" and (r.get("lung") or r.get("gene_keyword"))
        }
    )
    prjs = [p for p in prjs if p.startswith("PRJCA")]
    print(f"Fetching {len(prjs)} BioProject pages...")
    bp_rows = []
    for p in prjs:
        url = f"https://ngdc.cncb.ac.cn/bioproject/browse/{p}"
        page = fetch(url)
        text = strip_html(page.get("text") or "") if page.get("ok") else ""
        accs = sorted(set(re.findall(r"(?:OMIX|HRA|CRA|GSE)\d+", text)))
        f = flags(text)
        bp_rows.append(
            {
                "bioproject": p,
                "url": url,
                "page_ok": int(bool(page.get("ok"))),
                "linked_accessions": ";".join(accs),
                **{f"bp_{k}": v for k, v in f.items()},
                "snippet": text[text.find("描述") : text.find("描述") + 400]
                if "描述" in text
                else text[:400],
            }
        )
    write_tsv(os.path.join(OUT, "bioprojects.tsv"), bp_rows)

    # Verdict table: can we actually test TACSTD2/CLDN4 vs MPR?
    verdicts = []
    for r in omix_details:
        if not (r.get("desc_lung") or r.get("lung")):
            continue
        if not (
            r.get("desc_ici")
            or r.get("ici")
            or r.get("desc_gene_keyword")
            or r.get("gene_keyword")
            or r.get("desc_mpr_keyword")
        ):
            continue
        file_scans = [s for s in scans if s["accession"] == r["accession"]]
        has_t = any(s["TACSTD2"] for s in file_scans)
        has_c = any(s["CLDN4"] for s in file_scans)
        has_mpr = bool(r.get("desc_mpr_keyword") or r.get("mpr_keyword")) or any(
            MPR_RE.search(s.get("phenotype_snippets") or "") for s in file_scans
        )
        has_resp = any(s.get("phenotype_snippets") for s in file_scans)
        verdicts.append(
            {
                "source": "OMIX",
                "accession": r["accession"],
                "bioproject": r["bioproject"],
                "access": "Open",
                "title": r["title"][:180],
                "data_type": r.get("data_type") or "",
                "n_files_downloaded": len(file_scans),
                "has_TACSTD2": int(has_t),
                "has_CLDN4": int(has_c),
                "has_MPR_label": int(has_mpr),
                "has_any_response_label": int(has_resp),
                "usable_for_TACSTD2_vs_MPR": int(has_t and has_mpr),
                "usable_for_CLDN4_vs_MPR": int(has_c and has_mpr),
                "note": (
                    "open processed"
                    if file_scans
                    else "open page but no processed file downloaded"
                ),
            }
        )
    for r in hra_details:
        if r["access"] != "Open":
            continue
        if not (r.get("lung") or r.get("gene_keyword") or r.get("ici")):
            continue
        kinds = json.loads(r.get("listing_kinds") or "{}")
        is_fastq = kinds.get("fastq", 0) > 0
        verdicts.append(
            {
                "source": "HRA",
                "accession": r["accession"],
                "bioproject": r["bioproject"],
                "access": "Open",
                "title": r["title"][:180],
                "data_type": r.get("study_type") or "",
                "n_files_downloaded": 0,
                "has_TACSTD2": 0,
                "has_CLDN4": 0,
                "has_MPR_label": int(bool(r.get("mpr_keyword"))),
                "has_any_response_label": int(bool(r.get("ici"))),
                "usable_for_TACSTD2_vs_MPR": 0,
                "usable_for_CLDN4_vs_MPR": 0,
                "note": (
                    "OPEN but RAW FASTQ only — not processed"
                    if is_fastq
                    else f"OPEN listing={r.get('listing_kinds')}"
                ),
            }
        )
    write_tsv(os.path.join(OUT, "verdicts.tsv"), verdicts)

    stats["omix_files_downloaded"] = sum(1 for j in dl_log if j["ok"])
    stats["omix_files_failed"] = sum(1 for j in dl_log if not j["ok"])
    stats["any_TACSTD2_in_open_processed"] = int(any(s["TACSTD2"] for s in scans))
    stats["any_CLDN4_in_open_processed"] = int(any(s["CLDN4"] for s in scans))
    stats["any_usable_TACSTD2_vs_MPR"] = int(
        any(v["usable_for_TACSTD2_vs_MPR"] for v in verdicts)
    )
    stats["any_usable_CLDN4_vs_MPR"] = int(
        any(v["usable_for_CLDN4_vs_MPR"] for v in verdicts)
    )
    with open(os.path.join(OUT, "stats.json"), "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2)
    print(json.dumps(stats, indent=2))
    print("Done. Results in", OUT)


if __name__ == "__main__":
    main()
