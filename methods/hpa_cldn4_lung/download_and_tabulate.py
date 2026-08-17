#!/usr/bin/env python3
"""Download official Human Protein Atlas CLDN4 (ENSG00000189143) public tables.

Produces ordinal IHC High/Medium/Low/Not-detected counts for lung-cancer TMA
patients, split by SNOMED-labeled adenocarcinoma (LUAD) vs squamous carcinoma
(LUSC) when those labels exist. Also extracts consensus tissue RNA nTPM.

Semi-quantitative ordinal scores only. Does not invent H-scores, pair cores,
or add ICI / survival labels that HPA does not publish.
"""

from __future__ import annotations

import csv
import hashlib
import json
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

GENE = "CLDN4"
ENSEMBL = "ENSG00000189143"
HPA = "https://www.proteinatlas.org"
HPA_VERSION = "25.1"

URLS = {
    "CLDN4.xml": f"{HPA}/{ENSEMBL}.xml",
    "CLDN4.tsv": f"{HPA}/{ENSEMBL}.tsv",
    "CLDN4.json": f"{HPA}/{ENSEMBL}.json",
    "cancer_data.tsv.zip": f"{HPA}/download/tsv/cancer_data.tsv.zip",
    "rna_tissue_consensus.tsv.zip": f"{HPA}/download/tsv/rna_tissue_consensus.tsv.zip",
}

# HPA lung-cancer TMA uses SNOMED histology strings, not TCGA LUAD/LUSC codes.
HISTOLOGY_MAP = {
    "adenocarcinoma, nos": "LUAD",
    "squamous cell carcinoma, nos": "LUSC",
}

STAIN_ORDER = ("High", "Medium", "Low", "Not detected")


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path, retries: int = 4) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    delay = 4
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url, headers={"User-Agent": "hpa-cldn4-lung/1.0 (public HPA extract)"}
            )
            with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    out.write(chunk)
            return
        except Exception as exc:  # noqa: BLE001
            last_err = exc
            if attempt == retries:
                raise
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(last_err)


def local(tag: str) -> str:
    return tag.split("}")[-1]


def xml_text(el: ET.Element | None) -> str:
    return (el.text or "").strip() if el is not None else ""


def xml_first(el: ET.Element, name: str) -> ET.Element | None:
    for child in el:
        if local(child.tag) == name:
            return child
    return None


def xml_all(el: ET.Element, name: str) -> list[ET.Element]:
    return [child for child in el if local(child.tag) == name]


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def unzip_first(zip_path: Path, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        if not names:
            raise RuntimeError(f"empty zip: {zip_path}")
        zf.extract(names[0], dest_dir)
        return dest_dir / names[0]


def snomed_from_patient(patient: ET.Element) -> list[str]:
    desc: list[str] = []
    seen: set[str] = set()
    for snomed in patient.iter():
        if local(snomed.tag) != "snomed":
            continue
        d = snomed.attrib.get("tissueDescription") or xml_text(snomed)
        code = snomed.attrib.get("snomedCode") or ""
        item = f"{d}|{code}" if code else d
        if item and item not in seen:
            seen.add(item)
            desc.append(item)
    return desc


def histology_label(snomed: list[str]) -> tuple[str, str]:
    """Return (raw SNOMED histology, LUAD/LUSC/other bucket)."""
    raw_parts: list[str] = []
    for item in snomed:
        name = item.split("|", 1)[0]
        low = name.lower()
        if low in {"lung", "bronchus"} or name.startswith("Lung|") or "T-28000" in item:
            continue
        if low.startswith("normal tissue"):
            continue
        raw_parts.append(name)
    raw = "; ".join(dict.fromkeys(raw_parts)) if raw_parts else "unspecified"
    bucket = "other_labeled"
    for part in raw_parts:
        mapped = HISTOLOGY_MAP.get(part.lower())
        if mapped:
            bucket = mapped
            break
    if raw == "unspecified":
        bucket = "unspecified"
    return raw, bucket


def parse_levels(patient: ET.Element) -> dict[str, str]:
    out: dict[str, str] = {}
    for level in xml_all(patient, "level"):
        out[level.attrib.get("type", "level")] = xml_text(level)
    return out


def image_urls_from_patient(patient: ET.Element) -> list[str]:
    urls: list[str] = []
    for el in patient.iter():
        if local(el.tag) == "imageUrl":
            u = xml_text(el)
            if u:
                urls.append(u)
    return urls


def extract_cancer_patients(xml_path: Path) -> tuple[str, list[dict]]:
    root = ET.parse(xml_path).getroot()
    entry = next(el for el in root if local(el.tag) == "entry")
    hpa_version = entry.attrib.get("version", "")
    entry_url = entry.attrib.get("url", "")
    patients: list[dict] = []
    for antibody in [el for el in entry if local(el.tag) == "antibody"]:
        ab_id = antibody.attrib.get("id", "")
        for te in antibody:
            if local(te.tag) != "tissueExpression":
                continue
            if te.attrib.get("assayType") != "cancer":
                continue
            for data in xml_all(te, "data"):
                tissue_el = xml_first(data, "tissue")
                tissue = xml_text(tissue_el)
                if tissue.strip().lower() != "lung cancer":
                    continue
                for patient in xml_all(data, "patient"):
                    levels = parse_levels(patient)
                    snomed = snomed_from_patient(patient)
                    raw_hist, bucket = histology_label(snomed)
                    urls = image_urls_from_patient(patient)
                    patients.append(
                        {
                            "gene": GENE,
                            "ensembl": ENSEMBL,
                            "hpa_entry_version": hpa_version,
                            "antibody_id": ab_id,
                            "tissue": tissue,
                            "patient_id": xml_text(xml_first(patient, "patientId")),
                            "sex": xml_text(xml_first(patient, "sex")),
                            "age": xml_text(xml_first(patient, "age")),
                            "staining": levels.get("staining", ""),
                            "intensity": levels.get("intensity", ""),
                            "quantity": xml_text(xml_first(patient, "quantity")),
                            "location": xml_text(xml_first(patient, "location")),
                            "snomed": "; ".join(snomed),
                            "histology_from_snomed": raw_hist,
                            "histology_bucket": bucket,
                            "n_public_images": str(len(urls)),
                            "image_urls": " | ".join(urls),
                            "source_url": entry_url or f"{HPA}/{ENSEMBL}-{GENE}",
                            "score_type": "ordinal_IHC",
                            "note": "HPA published ordinal IHC only. Not an H-score.",
                        }
                    )
    return hpa_version, patients


def extract_consensus_rna_xml(xml_path: Path) -> list[dict]:
    root = ET.parse(xml_path).getroot()
    entry = next(el for el in root if local(el.tag) == "entry")
    rows: list[dict] = []
    for rna in [el for el in entry if local(el.tag) == "rnaExpression"]:
        assay = rna.attrib.get("assayType", "")
        source = rna.attrib.get("source", "")
        tech = rna.attrib.get("technology", "")
        if assay not in {"consensusTissue", "tissue"}:
            continue
        for data in xml_all(rna, "data"):
            tissue_el = xml_first(data, "tissue")
            tissue = xml_text(tissue_el)
            level_el = xml_first(data, "level")
            if level_el is None:
                continue
            unit = level_el.attrib.get("unitRNA", "")
            value = level_el.attrib.get("expRNA", "")
            if unit.lower() != "ntpm" and unit != "nTPM":
                # keep rows that still look like nTPM consensus
                if assay != "consensusTissue":
                    continue
            rows.append(
                {
                    "gene": GENE,
                    "ensembl": ENSEMBL,
                    "assay_type": assay,
                    "source": source,
                    "technology": tech,
                    "tissue": tissue,
                    "unit": unit or "nTPM",
                    "nTPM": value,
                    "note": "Companion RNA from official HPA gene XML. Not IHC.",
                }
            )
    return rows


def count_stains(rows: list[dict]) -> dict[str, int]:
    c = Counter((r.get("staining") or "").strip() for r in rows)
    return {k: int(c.get(k, 0)) for k in STAIN_ORDER}


def knowledge_lung_row(cancer_tsv: Path) -> dict:
    with cancer_tsv.open(newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = (row.get("Gene") or "").strip()
            name = (row.get("Gene name") or "").strip()
            cancer = (row.get("Cancer") or "").strip().lower()
            if cancer == "lung cancer" and (gene == ENSEMBL or name == GENE):
                return row
    raise RuntimeError("CLDN4 lung cancer row missing from cancer_data.tsv")


def consensus_from_tsv(tsv_path: Path) -> list[dict]:
    rows: list[dict] = []
    with tsv_path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = (row.get("Gene") or "").strip()
            name = (row.get("Gene name") or "").strip()
            if gene != ENSEMBL and name != GENE:
                continue
            tissue = (row.get("Tissue") or "").strip()
            ntpm = (row.get("nTPM") or row.get("NX") or "").strip()
            rows.append(
                {
                    "gene": GENE,
                    "ensembl": ENSEMBL,
                    "assay_type": "consensusTissue",
                    "source": "HPA+GTEx consensus",
                    "technology": "RNAseq",
                    "tissue": tissue,
                    "unit": "nTPM",
                    "nTPM": ntpm,
                    "note": "Official rna_tissue_consensus.tsv. Not IHC.",
                }
            )
    return rows


def main() -> int:
    here = Path(__file__).resolve().parent
    cache = here / ".cache"
    tables = here / "tables"
    cache.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    retrieved = utc_now()
    manifest: list[dict] = []
    for name, url in URLS.items():
        dest = cache / name
        if not dest.exists() or dest.stat().st_size == 0:
            download(url, dest)
        manifest.append(
            {
                "file": name,
                "url": url,
                "bytes": str(dest.stat().st_size),
                "sha256": sha256_file(dest),
                "retrieved_utc": retrieved,
                "hpa_version_expected": HPA_VERSION,
            }
        )

    xml_path = cache / "CLDN4.xml"
    hpa_version, patients = extract_cancer_patients(xml_path)
    if not patients:
        raise RuntimeError("no CLDN4 lung-cancer patients in gene XML")

    cancer_zip = cache / "cancer_data.tsv.zip"
    cancer_tsv = unzip_first(cancer_zip, cache / "unzip_cancer")
    know = knowledge_lung_row(cancer_tsv)
    know_row = {
        "gene": GENE,
        "ensembl": ENSEMBL,
        "cancer": know.get("Cancer", "lung cancer"),
        "High": know.get("High", ""),
        "Medium": know.get("Medium", ""),
        "Low": know.get("Low", ""),
        "Not detected": know.get("Not detected", ""),
        "n": str(
            sum(
                int(know.get(k) or 0)
                for k in ("High", "Medium", "Low", "Not detected")
            )
        ),
        "hpa_version": HPA_VERSION,
        "source_file": "cancer_data.tsv",
        "source_url": URLS["cancer_data.tsv.zip"],
        "score_type": "knowledge_based_ordinal_IHC_counts",
        "note": "Official HPA knowledge table. Not split by LUAD/LUSC. Not an H-score.",
    }

    rna_zip = cache / "rna_tissue_consensus.tsv.zip"
    rna_tsv = unzip_first(rna_zip, cache / "unzip_rna")
    consensus = consensus_from_tsv(rna_tsv)
    if not consensus:
        # fallback: gene XML consensusTissue rows
        consensus = [
            r
            for r in extract_consensus_rna_xml(xml_path)
            if r["assay_type"] == "consensusTissue"
        ]
    if not consensus:
        raise RuntimeError("no CLDN4 consensus nTPM rows")

    def ntpm_float(row: dict) -> float:
        return float(row["nTPM"])

    consensus_sorted = sorted(consensus, key=ntpm_float, reverse=True)
    for i, row in enumerate(consensus_sorted, start=1):
        row["rank_desc"] = str(i)
        row["n_tissues"] = str(len(consensus_sorted))
    lung_rows = [r for r in consensus_sorted if r["tissue"].strip().lower() == "lung"]
    if not lung_rows:
        raise RuntimeError("lung nTPM missing from consensus table")
    lung_ntpm = ntpm_float(lung_rows[0])
    other = [r for r in consensus_sorted if r["tissue"].strip().lower() != "lung"]
    other_vals = [ntpm_float(r) for r in other]
    other_vals_sorted = sorted(other_vals)
    mid = len(other_vals_sorted) // 2
    if len(other_vals_sorted) % 2:
        other_median = other_vals_sorted[mid]
    else:
        other_median = (other_vals_sorted[mid - 1] + other_vals_sorted[mid]) / 2
    other_mean = sum(other_vals) / len(other_vals)
    other_max_row = max(other, key=ntpm_float)
    lung_rank = int(lung_rows[0]["rank_desc"])

    # Count tables
    antibodies = sorted({p["antibody_id"] for p in patients})
    count_rows: list[dict] = []

    def add_count(scope: str, subset: list[dict], extra: str = "") -> None:
        stains = count_stains(subset)
        n = len(subset)
        count_rows.append(
            {
                "gene": GENE,
                "ensembl": ENSEMBL,
                "antibody_id": ",".join(antibodies),
                "scope": scope,
                "histology_bucket": extra or scope,
                "High": str(stains["High"]),
                "Medium": str(stains["Medium"]),
                "Low": str(stains["Low"]),
                "Not_detected": str(stains["Not detected"]),
                "n_patients": str(n),
                "pct_High": f"{100 * stains['High'] / n:.1f}" if n else "",
                "pct_Medium": f"{100 * stains['Medium'] / n:.1f}" if n else "",
                "pct_Low": f"{100 * stains['Low'] / n:.1f}" if n else "",
                "pct_Not_detected": f"{100 * stains['Not detected'] / n:.1f}" if n else "",
                "score_type": "ordinal_IHC_counts",
                "note": (
                    "Semi-quantitative HPA staining categories only. "
                    "No H-score was calculated from intensity x quantity."
                ),
            }
        )

    add_count("all_lung_cancer_TMA", patients, "all_labeled_plus_other")
    for bucket, label in (
        ("LUAD", "LUAD_SNOMED_adenocarcinoma_NOS"),
        ("LUSC", "LUSC_SNOMED_squamous_cell_carcinoma_NOS"),
        ("other_labeled", "other_SNOMED_not_LUAD_LUSC"),
        ("unspecified", "unspecified"),
    ):
        add_count(label, [p for p in patients if p["histology_bucket"] == bucket], bucket)

    write_tsv(
        tables / "01_download_manifest.tsv",
        manifest,
        ["file", "url", "bytes", "sha256", "retrieved_utc", "hpa_version_expected"],
    )
    write_tsv(
        tables / "02_knowledge_lung_cancer_counts.tsv",
        [know_row],
        [
            "gene",
            "ensembl",
            "cancer",
            "High",
            "Medium",
            "Low",
            "Not detected",
            "n",
            "hpa_version",
            "source_file",
            "source_url",
            "score_type",
            "note",
        ],
    )
    write_tsv(
        tables / "03_ihc_patients.tsv",
        patients,
        [
            "gene",
            "ensembl",
            "hpa_entry_version",
            "antibody_id",
            "tissue",
            "patient_id",
            "sex",
            "age",
            "staining",
            "intensity",
            "quantity",
            "location",
            "snomed",
            "histology_from_snomed",
            "histology_bucket",
            "n_public_images",
            "image_urls",
            "source_url",
            "score_type",
            "note",
        ],
    )
    write_tsv(
        tables / "04_ihc_counts_luad_lusc.tsv",
        count_rows,
        [
            "gene",
            "ensembl",
            "antibody_id",
            "scope",
            "histology_bucket",
            "High",
            "Medium",
            "Low",
            "Not_detected",
            "n_patients",
            "pct_High",
            "pct_Medium",
            "pct_Low",
            "pct_Not_detected",
            "score_type",
            "note",
        ],
    )
    write_tsv(
        tables / "05_ntpm_tissue_consensus.tsv",
        consensus_sorted,
        [
            "gene",
            "ensembl",
            "assay_type",
            "source",
            "technology",
            "tissue",
            "unit",
            "nTPM",
            "rank_desc",
            "n_tissues",
            "note",
        ],
    )
    write_tsv(
        tables / "06_ntpm_lung_vs_other.tsv",
        [
            {
                "gene": GENE,
                "ensembl": ENSEMBL,
                "unit": "nTPM",
                "lung_nTPM": f"{lung_ntpm:g}",
                "lung_rank_desc": str(lung_rank),
                "n_tissues": str(len(consensus_sorted)),
                "other_tissues_n": str(len(other)),
                "other_median_nTPM": f"{other_median:g}",
                "other_mean_nTPM": f"{other_mean:.4g}",
                "other_max_tissue": other_max_row["tissue"],
                "other_max_nTPM": other_max_row["nTPM"],
                "lung_over_other_median": f"{(lung_ntpm / other_median):.3g}"
                if other_median
                else "",
                "source": "rna_tissue_consensus.tsv",
                "source_url": URLS["rna_tissue_consensus.tsv.zip"],
                "note": "Consensus nTPM is RNA, not IHC. HPA consensus = max(HPA, GTEx) per tissue.",
            }
        ],
        [
            "gene",
            "ensembl",
            "unit",
            "lung_nTPM",
            "lung_rank_desc",
            "n_tissues",
            "other_tissues_n",
            "other_median_nTPM",
            "other_mean_nTPM",
            "other_max_tissue",
            "other_max_nTPM",
            "lung_over_other_median",
            "source",
            "source_url",
            "note",
        ],
    )

    summary = {
        "gene": GENE,
        "ensembl": ENSEMBL,
        "hpa_version": HPA_VERSION,
        "xml_entry_version": hpa_version,
        "retrieved_utc": retrieved,
        "n_lung_cancer_patients": len(patients),
        "antibodies": antibodies,
        "knowledge_counts": {
            k: know_row[k] for k in ("High", "Medium", "Low", "Not detected", "n")
        },
        "ihc_counts": count_rows,
        "lung_nTPM": lung_ntpm,
        "lung_rank_desc": lung_rank,
        "n_tissues": len(consensus_sorted),
        "other_median_nTPM": other_median,
        "other_max_tissue": other_max_row["tissue"],
        "other_max_nTPM": float(other_max_row["nTPM"]),
        "pages": {
            "gene": f"{HPA}/{ENSEMBL}-{GENE}",
            "lung_cancer": f"{HPA}/{ENSEMBL}-{GENE}/pathology/lung+cancer",
            "tissue": f"{HPA}/{ENSEMBL}-{GENE}/tissue",
        },
    }
    (tables / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
