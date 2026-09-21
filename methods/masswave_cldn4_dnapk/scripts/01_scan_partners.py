#!/usr/bin/env python3
"""Scan open claudin-4 protein tables for DNA-PKcs, Ku, and STING.

Inputs are local files downloaded by 02_download.py. No RAW files.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(parents=True, exist_ok=True)

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
# Gene symbols and protein-name phrases. "Ku" alone is not used: it matches
# Kunitz (SPINT2) in the canine FASTA headers.
TARGET_RE = re.compile(
    r"\bPRKDC\b|\bXRCC4\b|\bXRCC5\b|\bXRCC6\b|\bSTING1\b|\bTMEM173\b|"
    r"\bMB21D1\b|\bCGAS\b|"
    r"DNA-dependent protein kinase|DNA dependent protein kinase|\bDNA-PKcs\b|"
    r"\bKu70\b|\bKu80\b|\bKu86\b|"
    r"X-ray repair cross|"
    r"stimulator of interferon genes",
    re.I,
)
TARGETS = ["PRKDC", "XRCC5", "XRCC6", "XRCC4", "STING1", "TMEM173", "CGAS", "MB21D1"]


def sha1(path: Path) -> str:
    h = hashlib.sha1()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_shared_strings(z: zipfile.ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in z.namelist():
        return []
    root = ET.fromstring(z.read("xl/sharedStrings.xml"))
    out = []
    for si in root.findall(f"{NS}si"):
        out.append("".join((t.text or "") for t in si.iter(f"{NS}t")))
    return out


def sheet_rows(z: zipfile.ZipFile, sheet: str, strings: list[str]) -> list[dict]:
    root = ET.fromstring(z.read(sheet))
    rows = []
    for row in root.iter(f"{NS}row"):
        vals = {}
        for cell in row.findall(f"{NS}c"):
            ref = cell.get("r") or ""
            match = re.match(r"([A-Z]+)", ref)
            if not match:
                continue
            col = match.group(1)
            node = cell.find(f"{NS}v")
            if node is None or node.text is None:
                val = ""
            elif cell.get("t") == "s":
                val = strings[int(node.text)]
            else:
                val = node.text
            vals[col] = val
        if vals:
            rows.append(vals)
    return rows


def scan_pride_coip(path: Path) -> dict:
    hits = []
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        fasta = row.get("Fasta headers") or ""
        blob = " ".join(
            [
                row.get("Gene names") or "",
                row.get("Protein names") or "",
                fasta,
                row.get("Majority protein IDs") or "",
            ]
        )
        found = sorted({m.group(0) for m in TARGET_RE.finditer(blob)})
        if not found:
            continue
        genes = sorted(set(re.findall(r"GN=([A-Za-z0-9_-]+)", fasta)))
        hits.append(
            {
                "source": "PXD031094 proteinGroups_Cldn4.txt",
                "genes": ";".join(genes),
                "matched": ";".join(found),
                "unique_peptides": row.get("Unique peptides") or "",
                "lfq_cldn4": ";".join(row.get(f"LFQ intensity Cldn4_{i}") or "" for i in range(1, 5)),
                "lfq_gfp": ";".join(row.get(f"LFQ intensity GFPctrl_{i}") or "" for i in range(1, 5)),
                "protein": fasta.split(";")[0][:180],
            }
        )
    # PARP1 is recorded as a DNA-repair neighbor that is not a requested target.
    parp = []
    for row in rows:
        fasta = row.get("Fasta headers") or ""
        if re.search(r"GN=PARP1\b", fasta):
            parp.append(
                {
                    "gene": "PARP1",
                    "unique_peptides": row.get("Unique peptides") or "",
                    "lfq_cldn4": [row.get(f"LFQ intensity Cldn4_{i}") or "0" for i in range(1, 5)],
                    "lfq_gfp": [row.get(f"LFQ intensity GFPctrl_{i}") or "0" for i in range(1, 5)],
                }
            )
    controls = []
    for symbol in ("CLDN4", "OAS1"):
        for row in rows:
            fasta = row.get("Fasta headers") or ""
            if re.search(rf"GN={symbol}\b", fasta):
                controls.append(
                    {
                        "gene": symbol,
                        "ids": row.get("Majority protein IDs") or "",
                        "unique_peptides": row.get("Unique peptides") or "",
                        "lfq_cldn4": [row.get(f"LFQ intensity Cldn4_{i}") or "0" for i in range(1, 5)],
                        "lfq_gfp": [row.get(f"LFQ intensity GFPctrl_{i}") or "0" for i in range(1, 5)],
                    }
                )
    return {
        "n_protein_groups": len(rows),
        "target_hits": hits,
        "parp1_rows": parp,
        "controls": controls,
    }


def scan_plos(path: Path) -> dict:
    z = zipfile.ZipFile(path)
    strings = load_shared_strings(z)
    rows = sheet_rows(z, "xl/worksheets/sheet1.xml", strings)
    # Row 2 (index 1) is the header. BL-Cldn4 names are column M.
    names = []
    target_hits = []
    for vals in rows[2:]:
        name = (vals.get("M") or "").strip()
        if not name:
            continue
        names.append(
            {
                "rank": vals.get("K") or "",
                "accession": vals.get("L") or "",
                "protein": name,
                "av_n_psm_opn": vals.get("N") or "",
            }
        )
        if TARGET_RE.search(name):
            target_hits.append(name)
    out = RESULTS / "plos2015_bl_cldn4_proteins.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rank", "accession", "protein", "av_n_psm_opn"])
        writer.writeheader()
        writer.writerows(names)
    return {
        "n_bl_cldn4_proteins": len(names),
        "target_hits": target_hits,
        "table": str(out.relative_to(ROOT)),
    }


def docx_text(path: Path) -> str:
    z = zipfile.ZipFile(path)
    xml = z.read("word/document.xml").decode("utf-8", errors="replace")
    text = re.sub(r"</w:p>", "\n", xml)
    text = re.sub(r"<[^>]+>", "", text)
    for src, dst in (("&amp;", "&"), ("&lt;", "<"), ("&gt;", ">")):
        text = text.replace(src, dst)
    return text


def scan_docx(path: Path) -> dict:
    text = docx_text(path)
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    # Header then repeating quartets: gene, control, cldn4, log2fc.
    start = 0
    for i, ln in enumerate(lines):
        if ln == "Alternate ID":
            start = i + 4
            break
    records = []
    i = start
    while i + 3 < len(lines):
        gene, ctrl, bait, fc = lines[i : i + 4]
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", gene):
            break
        if not re.fullmatch(r"-?\d+(\.\d+)?", ctrl):
            break
        records.append(
            {
                "gene": gene,
                "control_sample_1": ctrl,
                "cldn4_sample_2": bait,
                "log2fc_cldn4_over_control": fc,
            }
        )
        i += 4
    target_hits = [r for r in records if TARGET_RE.search(r["gene"])]
    above = []
    for rec in records:
        try:
            if float(rec["log2fc_cldn4_over_control"]) > 1:
                above.append(rec["gene"])
        except ValueError:
            pass
    out = RESULTS / "crc24_bioid_suppst1.tsv"
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["gene", "control_sample_1", "cldn4_sample_2", "log2fc_cldn4_over_control"],
        )
        writer.writeheader()
        writer.writerows(records)
    return {
        "n_rows": len(records),
        "n_log2fc_gt_1": len(above),
        "target_hits": [r["gene"] for r in target_hits],
        "genes": [r["gene"] for r in records],
        "table": str(out.relative_to(ROOT)),
    }


def main() -> None:
    coip_path = DATA / "proteinGroups_Cldn4.txt"
    plos_path = DATA / "plos2015_s2_all_proteins.xlsx"
    docx_path = DATA / "crc24_suppst1_bioid.docx"
    coip = scan_pride_coip(coip_path)
    plos = scan_plos(plos_path)
    s3_path = DATA / "plos2015_s3_enriched.xlsx"
    s3_strings = load_shared_strings(zipfile.ZipFile(s3_path))
    s3_hits = sorted({s.strip() for s in s3_strings if TARGET_RE.search(s)})
    docx = scan_docx(docx_path)
    presence = []
    for source, hits in [
        ("PXD031094_Cldn4_CoIP", coip["target_hits"]),
        ("PLOS2015_BL_Cldn4", plos["target_hits"]),
        ("PLOS2015_S3_strings", s3_hits),
        ("CRC2024_BioID_S1", docx["target_hits"]),
    ]:
        presence.append(
            {
                "source": source,
                "n_target_hits": len(hits),
                "hits": ";".join(str(h) for h in hits),
            }
        )
    with (RESULTS / "target_presence.tsv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source", "n_target_hits", "hits"])
        writer.writeheader()
        writer.writerows(presence)
    with (RESULTS / "pxd031094_target_hits.tsv").open("w", newline="") as fh:
        fields = ["source", "genes", "matched", "unique_peptides", "lfq_cldn4", "lfq_gfp", "protein"]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(coip["target_hits"])
    summary = {
        "targets": TARGETS,
        "pxd031094": {
            "file": coip_path.name,
            "bytes": coip_path.stat().st_size,
            "sha1": sha1(coip_path),
            "pride_api_sha1": "4dbde6c0d827105412ef1f34cbb969b8dd8acb36",
            "sha1_matches_pride_api": sha1(coip_path) == "4dbde6c0d827105412ef1f34cbb969b8dd8acb36",
            "n_protein_groups": coip["n_protein_groups"],
            "n_target_hits": len(coip["target_hits"]),
            "parp1_rows": coip["parp1_rows"],
            "controls": coip["controls"],
        },
        "plos2015_s2": {
            "file": plos_path.name,
            "bytes": plos_path.stat().st_size,
            "sha1": sha1(plos_path),
            "n_bl_cldn4_proteins": plos["n_bl_cldn4_proteins"],
            "n_target_hits": len(plos["target_hits"]),
        },
        "plos2015_s3": {
            "file": s3_path.name,
            "bytes": s3_path.stat().st_size,
            "sha1": sha1(s3_path),
            "n_shared_strings": len(s3_strings),
            "n_target_hits": len(s3_hits),
            "target_hits": s3_hits,
        },
        "crc2024_bioid": {
            "file": docx_path.name,
            "bytes": docx_path.stat().st_size,
            "sha1": sha1(docx_path),
            "figshare_md5": "70bcbba8b14311102973f11eda6cc41d",
            "n_rows": docx["n_rows"],
            "n_log2fc_gt_1": docx["n_log2fc_gt_1"],
            "n_target_hits": len(docx["target_hits"]),
        },
    }
    (RESULTS / "scan_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
