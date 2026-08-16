#!/usr/bin/env python3
"""Download official Human Protein Atlas tables and extract lung TACSTD2/CLDN4 IHC/protein.

Honest public slices only. No invented H-scores, ICI labels, or paired cores.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from datetime import datetime, timezone
from pathlib import Path

GENES = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
}
GENE_IDS = set(GENES.values())
GENE_NAMES = set(GENES)

HPA_CURRENT = "https://www.proteinatlas.org"
HPA_V23 = "https://v23.proteinatlas.org"
HPA_VERSION_CURRENT = "25.1"
HPA_VERSION_V23 = "23.0"

RESPIRATORY = {"lung", "bronchus", "nasopharynx"}
LUNG_CANCER = {"lung cancer", "lung ac", "lung sqcc"}

URLS = {
    "normal_ihc_data.tsv.zip": f"{HPA_CURRENT}/download/tsv/normal_ihc_data.tsv.zip",
    "cancer_data.tsv.zip": f"{HPA_CURRENT}/download/tsv/cancer_data.tsv.zip",
    "cancer_prognostic_data.tsv.zip": f"{HPA_CURRENT}/download/tsv/cancer_prognostic_data.tsv.zip",
    "cancer_cptac.tsv.zip": f"{HPA_CURRENT}/download/tsv/cancer_cptac.tsv.zip",
    "ms_tissue.tsv.zip": f"{HPA_CURRENT}/download/tsv/ms_tissue.tsv.zip",
    "ms_tissue_sample_data.tsv.zip": f"{HPA_CURRENT}/download/tsv/ms_tissue_sample_data.tsv.zip",
    "v23_normal_tissue.tsv.zip": f"{HPA_V23}/download/normal_tissue.tsv.zip",
    "v23_pathology.tsv.zip": f"{HPA_V23}/download/pathology.tsv.zip",
    "TACSTD2.tsv": f"{HPA_CURRENT}/ENSG00000184292.tsv",
    "CLDN4.tsv": f"{HPA_CURRENT}/ENSG00000189143.tsv",
    "TACSTD2.json": f"{HPA_CURRENT}/ENSG00000184292.json",
    "CLDN4.json": f"{HPA_CURRENT}/ENSG00000189143.json",
    "TACSTD2.xml": f"{HPA_CURRENT}/ENSG00000184292.xml",
    "CLDN4.xml": f"{HPA_CURRENT}/ENSG00000189143.xml",
}


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
    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "w200-hpa-lung/1.0"})
            with urllib.request.urlopen(req, timeout=120) as resp, dest.open("wb") as out:
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


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def unzip_member(zip_path: Path, dest_dir: Path) -> Path:
    dest_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        if not names:
            raise RuntimeError(f"empty zip: {zip_path}")
        zf.extract(names[0], dest_dir)
        return dest_dir / names[0]


def gene_match(row: dict) -> bool:
    gene = row.get("Gene") or row.get("Ensembl") or ""
    name = row.get("Gene name") or row.get("Gene") or ""
    return gene in GENE_IDS or name in GENE_NAMES


def tissue_is_respiratory(value: str) -> bool:
    return (value or "").strip().lower() in RESPIRATORY


def cancer_is_lung(value: str) -> bool:
    v = (value or "").strip().lower()
    return v in LUNG_CANCER or v.startswith("lung ")


def snomed_from_patient(patient: ET.Element) -> list[str]:
    desc = []
    for snomed in patient.iter():
        if local(snomed.tag) == "snomed":
            d = snomed.attrib.get("tissueDescription") or xml_text(snomed)
            code = snomed.attrib.get("snomedCode") or ""
            if d or code:
                desc.append(f"{d}|{code}" if code else d)
    # unique preserve order
    seen = set()
    out = []
    for item in desc:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def image_urls_from_patient(patient: ET.Element) -> list[str]:
    urls = []
    for el in patient.iter():
        if local(el.tag) == "imageUrl":
            u = xml_text(el)
            if u:
                urls.append(u)
    return urls


def parse_levels(patient: ET.Element) -> dict[str, str]:
    out = {}
    for level in xml_all(patient, "level"):
        out[level.attrib.get("type", "level")] = xml_text(level)
    return out


def extract_xml_ihc(xml_path: Path, gene_name: str, ensembl: str) -> tuple[list[dict], list[dict], list[dict]]:
    root = ET.parse(xml_path).getroot()
    entry = next(el for el in root if local(el.tag) == "entry")
    hpa_version = entry.attrib.get("version", "")
    entry_url = entry.attrib.get("url", "")

    summaries = []
    for el in entry.iter():
        if local(el.tag) == "summary":
            summaries.append(
                {
                    "gene_name": gene_name,
                    "ensembl": ensembl,
                    "hpa_entry_version": hpa_version,
                    "summary_type": el.attrib.get("type", "subcellular"),
                    "summary_text": xml_text(el),
                    "source_url": entry_url,
                }
            )

    normal_ab = []
    cancer_patients = []
    for antibody in [el for el in entry if local(el.tag) == "antibody"]:
        ab_id = antibody.attrib.get("id", "")
        ab_release = antibody.attrib.get("releaseVersion", "")
        ab_date = antibody.attrib.get("releaseDate", "")
        for te in antibody:
            if local(te.tag) != "tissueExpression":
                continue
            assay = te.attrib.get("assayType", "")
            for data in xml_all(te, "data"):
                tissue_el = xml_first(data, "tissue")
                tissue = xml_text(tissue_el)
                organ = tissue_el.attrib.get("organ", "") if tissue_el is not None else ""
                if assay == "tissue" and not tissue_is_respiratory(tissue):
                    continue
                if assay == "cancer" and not cancer_is_lung(tissue):
                    continue
                for tc in xml_all(data, "tissueCell"):
                    level_el = xml_first(tc, "level")
                    normal_ab.append(
                        {
                            "gene_name": gene_name,
                            "ensembl": ensembl,
                            "hpa_entry_version": hpa_version,
                            "antibody_id": ab_id,
                            "antibody_release_version": ab_release,
                            "antibody_release_date": ab_date,
                            "assay_type": assay,
                            "tissue": tissue,
                            "organ": organ,
                            "cell_type": xml_text(xml_first(tc, "cellType")),
                            "staining": xml_text(level_el),
                            "staining_count": level_el.attrib.get("count", "") if level_el is not None else "",
                            "quantity": xml_text(xml_first(tc, "quantity")),
                            "location": xml_text(xml_first(tc, "location")),
                            "source_url": entry_url,
                        }
                    )
                if assay != "cancer":
                    continue
                for patient in xml_all(data, "patient"):
                    levels = parse_levels(patient)
                    snomed = snomed_from_patient(patient)
                    histologies = [s.split("|", 1)[0] for s in snomed if s.startswith(("Adenocarcinoma", "Squamous", "Carcinoma", "Small cell", "Large cell", "Bronchiolo", "Carcinoid", "Mesotheli"))]
                    if not histologies:
                        histologies = [s.split("|", 1)[0] for s in snomed if not s.startswith("Lung|") and "T-28000" not in s]
                    urls = image_urls_from_patient(patient)
                    cancer_patients.append(
                        {
                            "gene_name": gene_name,
                            "ensembl": ensembl,
                            "hpa_entry_version": hpa_version,
                            "antibody_id": ab_id,
                            "antibody_release_version": ab_release,
                            "antibody_release_date": ab_date,
                            "tissue": tissue,
                            "organ": organ,
                            "patient_id": xml_text(xml_first(patient, "patientId")),
                            "sex": xml_text(xml_first(patient, "sex")),
                            "age": xml_text(xml_first(patient, "age")),
                            "staining": levels.get("staining", ""),
                            "intensity": levels.get("intensity", ""),
                            "quantity": xml_text(xml_first(patient, "quantity")),
                            "location": xml_text(xml_first(patient, "location")),
                            "snomed": "; ".join(snomed),
                            "histology_from_snomed": "; ".join(dict.fromkeys(histologies)),
                            "n_public_images": str(len(urls)),
                            "image_urls": " | ".join(urls),
                            "source_url": entry_url,
                            "note": "Ordinal IHC only. Not an H-score. No ICI outcome. Patients are not paired across genes or antibodies.",
                        }
                    )
    return summaries, normal_ab, cancer_patients


def extract_rna_lung(xml_path: Path, gene_name: str, ensembl: str) -> list[dict]:
    root = ET.parse(xml_path).getroot()
    entry = next(el for el in root if local(el.tag) == "entry")
    rows = []
    for rna in [el for el in entry if local(el.tag) == "rnaExpression"]:
        assay = rna.attrib.get("assayType", "")
        source = rna.attrib.get("source", "")
        tech = rna.attrib.get("technology", "")
        for data in xml_all(rna, "data"):
            tissue_el = xml_first(data, "tissue")
            tissue = xml_text(tissue_el)
            if tissue.lower() != "lung":
                continue
            level_el = xml_first(data, "level")
            rows.append(
                {
                    "gene_name": gene_name,
                    "ensembl": ensembl,
                    "modality": "RNA",
                    "source": source,
                    "technology": tech,
                    "assay_type": assay,
                    "tissue": tissue,
                    "unit": level_el.attrib.get("unitRNA", "") if level_el is not None else "",
                    "value": level_el.attrib.get("expRNA", "") if level_el is not None else "",
                    "sample_id": "",
                    "sex": "",
                    "age": "",
                    "note": "Companion RNA from the same HPA gene XML. Not IHC/protein.",
                }
            )
            for sample in xml_all(data, "RNASample"):
                rows.append(
                    {
                        "gene_name": gene_name,
                        "ensembl": ensembl,
                        "modality": "RNA",
                        "source": source,
                        "technology": tech,
                        "assay_type": assay,
                        "tissue": tissue,
                        "unit": sample.attrib.get("unitRNA", ""),
                        "value": sample.attrib.get("expRNA", ""),
                        "sample_id": sample.attrib.get("sampleId", ""),
                        "sex": sample.attrib.get("sex", ""),
                        "age": sample.attrib.get("age", ""),
                        "note": "Companion RNA from the same HPA gene XML. Not IHC/protein.",
                    }
                )
    return rows


def main() -> int:
    repo = Path(__file__).resolve().parents[2]
    out = repo / "results" / "w200" / "HPA_lung"
    raw = out / "raw"
    tables = out / "tables"
    cache = out / ".cache"
    raw.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)

    retrieved = utc_now()
    manifest = []
    for name, url in URLS.items():
        dest = cache / name
        print(f"download {url}", flush=True)
        download(url, dest)
        rec = {
            "file": name,
            "url": url,
            "bytes": str(dest.stat().st_size),
            "sha256": sha256_file(dest),
            "retrieved_utc": retrieved,
        }
        manifest.append(rec)
        print(f"  {rec['bytes']} bytes sha256={rec['sha256'][:12]}...", flush=True)

    # copy compact official gene files into raw/
    for name in ("TACSTD2.tsv", "CLDN4.tsv", "TACSTD2.json", "CLDN4.json"):
        (raw / name).write_bytes((cache / name).read_bytes())

    unzipped = cache / "unzipped"
    paths = {}
    for zip_name, key in [
        ("normal_ihc_data.tsv.zip", "normal_ihc"),
        ("cancer_data.tsv.zip", "cancer"),
        ("cancer_prognostic_data.tsv.zip", "prognostic"),
        ("cancer_cptac.tsv.zip", "cptac"),
        ("ms_tissue.tsv.zip", "ms_tissue"),
        ("ms_tissue_sample_data.tsv.zip", "ms_sample"),
        ("v23_normal_tissue.tsv.zip", "v23_normal"),
        ("v23_pathology.tsv.zip", "v23_pathology"),
    ]:
        paths[key] = unzip_member(cache / zip_name, unzipped)

    # 01 gene identity from official gene TSV
    identity_rows = []
    keep_cols = [
        "Gene",
        "Gene synonym",
        "Ensembl",
        "Gene description",
        "Uniprot",
        "Evidence",
        "Antibody",
        "Reliability (IH)",
        "Reliability (IF)",
        "Subcellular location",
        "Subcellular main location",
        "Protein tissue specificity",
        "Protein tissue distribution",
        "Protein tissue specific Intensity",
        "Protein cell type specificity",
        "Protein cell type distribution",
        "RNA tissue specificity",
        "RNA tissue distribution",
    ]
    for gene_name in ("TACSTD2", "CLDN4"):
        row = read_tsv(cache / f"{gene_name}.tsv")[0]
        identity_rows.append(
            {
                **{c: row.get(c, "") for c in keep_cols},
                "hpa_version": HPA_VERSION_CURRENT,
                "source_url": URLS[f"{gene_name}.tsv"],
                "gene_page": f"{HPA_CURRENT}/{row.get('Ensembl')}-{gene_name}",
                "lung_tissue_page": f"{HPA_CURRENT}/{row.get('Ensembl')}-{gene_name}/tissue/lung",
                "lung_cancer_page": f"{HPA_CURRENT}/{row.get('Ensembl')}-{gene_name}/cancer/lung+cancer",
            }
        )
    write_tsv(
        tables / "01_gene_identity.tsv",
        identity_rows,
        keep_cols
        + ["hpa_version", "source_url", "gene_page", "lung_tissue_page", "lung_cancer_page"],
    )

    # 02 official v25 normal IHC, respiratory tissues only
    normal_rows = []
    for row in read_tsv(paths["normal_ihc"]):
        if gene_match(row) and tissue_is_respiratory(row.get("Tissue", "")):
            normal_rows.append(
                {
                    **row,
                    "hpa_version": HPA_VERSION_CURRENT,
                    "source_file": "normal_ihc_data.tsv",
                    "source_url": URLS["normal_ihc_data.tsv.zip"],
                    "score_type": "knowledge_based_ordinal_IHC",
                    "note": "Official HPA cell-type IHC annotation. Not an H-score.",
                }
            )
    write_tsv(
        tables / "02_normal_ihc_knowledge_respiratory.tsv",
        normal_rows,
        [
            "Gene",
            "Gene name",
            "Tissue",
            "IHC tissue name",
            "Cell type",
            "Level",
            "Reliability",
            "hpa_version",
            "source_file",
            "source_url",
            "score_type",
            "note",
        ],
    )

    # 03 official v25 cancer IHC lung
    cancer_rows = []
    for row in read_tsv(paths["cancer"]):
        if gene_match(row) and cancer_is_lung(row.get("Cancer", "")):
            high = int(row.get("High") or 0)
            med = int(row.get("Medium") or 0)
            low = int(row.get("Low") or 0)
            nd = int(row.get("Not detected") or 0)
            n = high + med + low + nd
            cancer_rows.append(
                {
                    **row,
                    "n_patients_in_counts": str(n),
                    "hpa_version": HPA_VERSION_CURRENT,
                    "source_file": "cancer_data.tsv",
                    "source_url": URLS["cancer_data.tsv.zip"],
                    "score_type": "knowledge_based_ordinal_IHC_counts",
                    "note": "HPA TMA ordinal counts (High/Medium/Low/Not detected). Not H-scores. No ICI labels. Not paired TACSTD2+CLDN4.",
                }
            )
    write_tsv(
        tables / "03_cancer_ihc_knowledge_lung.tsv",
        cancer_rows,
        [
            "Gene",
            "Gene name",
            "Cancer",
            "High",
            "Medium",
            "Low",
            "Not detected",
            "n_patients_in_counts",
            "hpa_version",
            "source_file",
            "source_url",
            "score_type",
            "note",
        ],
    )

    # 04 v23 vs v25 lung IHC comparison
    v23_normal = [
        r
        for r in read_tsv(paths["v23_normal"])
        if gene_match(r) and tissue_is_respiratory(r.get("Tissue", ""))
    ]
    v23_path = [
        r for r in read_tsv(paths["v23_pathology"]) if gene_match(r) and cancer_is_lung(r.get("Cancer", ""))
    ]
    compare = []
    for r in normal_rows:
        key = (r["Gene name"], r["Tissue"].lower(), r["Cell type"].lower())
        match = next(
            (
                x
                for x in v23_normal
                if (x["Gene name"], x["Tissue"].lower(), x["Cell type"].lower()) == key
            ),
            None,
        )
        compare.append(
            {
                "table": "normal_IHC",
                "gene_name": r["Gene name"],
                "tissue_or_cancer": r["Tissue"],
                "cell_type": r["Cell type"],
                "v25_level_or_counts": r["Level"],
                "v23_level_or_counts": match["Level"] if match else "ABSENT_IN_V23",
                "same": "yes" if match and match["Level"] == r["Level"] else "no",
            }
        )
    for r in cancer_rows:
        match = next((x for x in v23_path if x["Gene name"] == r["Gene name"]), None)
        v25c = f"H{r['High']}/M{r['Medium']}/L{r['Low']}/ND{r['Not detected']}"
        v23c = (
            f"H{match['High']}/M{match['Medium']}/L{match['Low']}/ND{match['Not detected']}"
            if match
            else "ABSENT_IN_V23"
        )
        compare.append(
            {
                "table": "cancer_IHC",
                "gene_name": r["Gene name"],
                "tissue_or_cancer": r["Cancer"],
                "cell_type": "tumor cells (knowledge-based counts)",
                "v25_level_or_counts": v25c,
                "v23_level_or_counts": v23c,
                "same": "yes" if v25c == v23c else "no",
            }
        )
    write_tsv(
        tables / "04_v23_vs_v25_lung_ihc.tsv",
        compare,
        [
            "table",
            "gene_name",
            "tissue_or_cancer",
            "cell_type",
            "v25_level_or_counts",
            "v23_level_or_counts",
            "same",
        ],
    )

    # XML antibody-level and patient-level
    summaries = []
    normal_ab = []
    cancer_patients = []
    rna_rows = []
    for gene_name, ensembl in GENES.items():
        s, n, c = extract_xml_ihc(cache / f"{gene_name}.xml", gene_name, ensembl)
        summaries.extend(s)
        normal_ab.extend(n)
        cancer_patients.extend(c)
        rna_rows.extend(extract_rna_lung(cache / f"{gene_name}.xml", gene_name, ensembl))

    write_tsv(
        tables / "05_hpa_xml_summaries.tsv",
        summaries,
        [
            "gene_name",
            "ensembl",
            "hpa_entry_version",
            "summary_type",
            "summary_text",
            "source_url",
        ],
    )
    write_tsv(
        tables / "06_antibody_normal_ihc_respiratory.tsv",
        [r for r in normal_ab if r["assay_type"] == "tissue"],
        [
            "gene_name",
            "ensembl",
            "hpa_entry_version",
            "antibody_id",
            "antibody_release_version",
            "antibody_release_date",
            "tissue",
            "organ",
            "cell_type",
            "staining",
            "quantity",
            "location",
            "source_url",
        ],
    )
    write_tsv(
        tables / "07_antibody_cancer_ihc_patients.tsv",
        cancer_patients,
        [
            "gene_name",
            "ensembl",
            "hpa_entry_version",
            "antibody_id",
            "antibody_release_version",
            "antibody_release_date",
            "tissue",
            "patient_id",
            "sex",
            "age",
            "staining",
            "intensity",
            "quantity",
            "location",
            "histology_from_snomed",
            "snomed",
            "n_public_images",
            "image_urls",
            "source_url",
            "note",
        ],
    )

    # 08 MS tissue protein
    ms_rows = []
    for row in read_tsv(paths["ms_tissue"]):
        if gene_match(row) and row.get("Tissue", "").lower() == "lung":
            ms_rows.append(
                {
                    **row,
                    "hpa_version": HPA_VERSION_CURRENT,
                    "source_file": "ms_tissue.tsv",
                    "source_url": URLS["ms_tissue.tsv.zip"],
                    "modality": "mass_spectrometry_bulk_tissue",
                    "note": "HPA bulk tissue MS intensity. Not IHC. Not an H-score.",
                }
            )
    write_tsv(
        tables / "08_ms_tissue_lung.tsv",
        ms_rows,
        [
            "Gene",
            "Gene name",
            "Tissue",
            "Intensity",
            "hpa_version",
            "source_file",
            "source_url",
            "modality",
            "note",
        ],
    )
    ms_rep = []
    for row in read_tsv(paths["ms_sample"]):
        if gene_match(row) and row.get("Tissue", "").lower() == "lung":
            ms_rep.append(
                {
                    **row,
                    "hpa_version": HPA_VERSION_CURRENT,
                    "source_file": "ms_tissue_sample_data.tsv",
                    "source_url": URLS["ms_tissue_sample_data.tsv.zip"],
                    "modality": "mass_spectrometry_bulk_tissue_replicate",
                    "note": "HPA bulk tissue MS replicate intensity. Not IHC.",
                }
            )
    write_tsv(
        tables / "09_ms_tissue_lung_replicates.tsv",
        ms_rep,
        [
            "Gene",
            "Gene name",
            "sample_name",
            "Tissue",
            "replicate_nr",
            "Intensity",
            "hpa_version",
            "source_file",
            "source_url",
            "modality",
            "note",
        ],
    )

    # 09 CPTAC via HPA: report present and explicitly absent lung rows
    cptac_all = [r for r in read_tsv(paths["cptac"]) if gene_match(r)]
    cptac_lung = [r for r in cptac_all if cancer_is_lung(r.get("Cancer", ""))]
    cptac_out = []
    for gene_name, ensembl in GENES.items():
        hits = [r for r in cptac_lung if r.get("Gene name") == gene_name]
        if hits:
            for r in hits:
                cptac_out.append(
                    {
                        **r,
                        "present_in_hpa_cptac_lung": "yes",
                        "hpa_version": HPA_VERSION_CURRENT,
                        "source_file": "cancer_cptac.tsv",
                        "source_url": URLS["cancer_cptac.tsv.zip"],
                        "note": "HPA-hosted CPTAC tumor-vs-normal MS logFC. Not IHC. Treatment-naive CPTAC, not ICI.",
                    }
                )
        else:
            other = sorted({r["Cancer"] for r in cptac_all if r.get("Gene name") == gene_name})
            cptac_out.append(
                {
                    "Cancer": "Lung AC; Lung SQCC",
                    "Gene": ensembl,
                    "Gene name": gene_name,
                    "p-value adjusted": "",
                    "logFC": "",
                    "present_in_hpa_cptac_lung": "no",
                    "hpa_version": HPA_VERSION_CURRENT,
                    "source_file": "cancer_cptac.tsv",
                    "source_url": URLS["cancer_cptac.tsv.zip"],
                    "note": "No lung CPTAC row in the official HPA cancer_cptac.tsv. Other cancers present: "
                    + (", ".join(other) if other else "none"),
                }
            )
    write_tsv(
        tables / "10_cptac_hpa_lung.tsv",
        cptac_out,
        [
            "Cancer",
            "Gene",
            "Gene name",
            "p-value adjusted",
            "logFC",
            "present_in_hpa_cptac_lung",
            "hpa_version",
            "source_file",
            "source_url",
            "note",
        ],
    )

    # prognostic RNA (not IHC) lung only, labeled
    prog = []
    for row in read_tsv(paths["prognostic"]):
        if gene_match(row) and "lung" in (row.get("Cancer") or "").lower():
            prog.append(
                {
                    **row,
                    "hpa_version": HPA_VERSION_CURRENT,
                    "source_file": "cancer_prognostic_data.tsv",
                    "source_url": URLS["cancer_prognostic_data.tsv.zip"],
                    "modality": "TCGA_RNA_survival_association",
                    "note": "HPA RNA-survival labels. Not IHC/protein. Not ICI outcome. Every lung row in this official file is unprognostic.",
                }
            )
    write_tsv(
        tables / "11_rna_prognostic_lung_companion.tsv",
        prog,
        list(prog[0].keys()) if prog else ["note"],
    )

    # transparent recount of XML patients vs knowledge-based cancer_data.tsv
    from collections import Counter

    antibody_counts = []
    by_ab = {}
    for row in cancer_patients:
        by_ab.setdefault((row["gene_name"], row["antibody_id"]), []).append(row)
    knowledge = {(r["Gene name"], r["Cancer"]): r for r in cancer_rows}
    for (gene_name, ab_id), rs in by_ab.items():
        stain = Counter(r["staining"] for r in rs)
        hist = Counter(r["histology_from_snomed"] or "unspecified" for r in rs)
        k = knowledge.get((gene_name, "lung cancer"))
        k_counts = (
            f"H{k['High']}/M{k['Medium']}/L{k['Low']}/ND{k['Not detected']}" if k else ""
        )
        a_counts = (
            f"H{stain.get('High', 0)}/M{stain.get('Medium', 0)}/"
            f"L{stain.get('Low', 0)}/ND{stain.get('Not detected', 0)}"
        )
        antibody_counts.append(
            {
                "gene_name": gene_name,
                "antibody_id": ab_id,
                "n_patients": str(len(rs)),
                "high": str(stain.get("High", 0)),
                "medium": str(stain.get("Medium", 0)),
                "low": str(stain.get("Low", 0)),
                "not_detected": str(stain.get("Not detected", 0)),
                "xml_counts": a_counts,
                "knowledge_cancer_data_counts": k_counts,
                "matches_knowledge_table": "yes" if a_counts == k_counts else "no",
                "histology_from_snomed": "; ".join(f"{h}={n}" for h, n in hist.most_common()),
                "note": "Recount of official gene XML patient staining. Knowledge table is a single collective score, not the mean of all antibodies.",
            }
        )
    write_tsv(
        tables / "14_cancer_ihc_counts_by_antibody.tsv",
        antibody_counts,
        [
            "gene_name",
            "antibody_id",
            "n_patients",
            "high",
            "medium",
            "low",
            "not_detected",
            "xml_counts",
            "knowledge_cancer_data_counts",
            "matches_knowledge_table",
            "histology_from_snomed",
            "note",
        ],
    )

    write_tsv(
        tables / "12_rna_lung_companion.tsv",
        rna_rows,
        [
            "gene_name",
            "ensembl",
            "modality",
            "source",
            "technology",
            "assay_type",
            "tissue",
            "unit",
            "value",
            "sample_id",
            "sex",
            "age",
            "note",
        ],
    )

    # gaps
    gaps = [
        {
            "requested": "Paired TACSTD2 + CLDN4 IHC on the same lung cores",
            "status": "not_in_HPA_tables",
            "detail": "HPA publishes each gene separately. Patient IDs do not overlap between TACSTD2 antibodies and CLDN4 CAB002610 in the lung-cancer XML extract.",
        },
        {
            "requested": "Numeric H-score or mIF co-localization",
            "status": "not_in_HPA_tables",
            "detail": "Public HPA IHC is ordinal (not detected / low / medium / high) plus intensity/quantity/location. No H-score column.",
        },
        {
            "requested": "ICI / PD-1 / PD-L1 outcome labels",
            "status": "not_in_HPA_tables",
            "detail": "HPA normal and cancer IHC TMAs have no immunotherapy annotation.",
        },
        {
            "requested": "CLDN4 CPTAC lung tumor-vs-normal in HPA cancer_cptac.tsv",
            "status": "absent",
            "detail": "Official file has TACSTD2 Lung AC and Lung SQCC rows only. CLDN4 appears in other cancers, not lung.",
        },
        {
            "requested": "Current-site standalone normal_tissue.tsv / pathology.tsv",
            "status": "retired_on_v25_site",
            "detail": "v25 uses normal_ihc_data.tsv and cancer_data.tsv. v23 copies remain at v23.proteinatlas.org and match these lung IHC calls.",
        },
    ]
    write_tsv(tables / "13_not_available.tsv", gaps, ["requested", "status", "detail"])

    write_tsv(
        out / "download_manifest.tsv",
        manifest,
        ["file", "url", "bytes", "sha256", "retrieved_utc"],
    )

    # compact source catalog
    catalog = [
        {
            "source": "HPA v25.1 normal_ihc_data.tsv",
            "url": URLS["normal_ihc_data.tsv.zip"],
            "what_was_kept": "TACSTD2/CLDN4 rows for lung, bronchus, nasopharynx",
            "output": "tables/02_normal_ihc_knowledge_respiratory.tsv",
        },
        {
            "source": "HPA v25.1 cancer_data.tsv",
            "url": URLS["cancer_data.tsv.zip"],
            "what_was_kept": "TACSTD2/CLDN4 lung cancer ordinal counts",
            "output": "tables/03_cancer_ihc_knowledge_lung.tsv",
        },
        {
            "source": "HPA gene XML ENSG00000184292 / ENSG00000189143",
            "url": f"{HPA_CURRENT}/ENSG00000184292.xml ; {HPA_CURRENT}/ENSG00000189143.xml",
            "what_was_kept": "Per-antibody respiratory IHC, lung-cancer patient ordinal scores + public image URLs, antibody vs knowledge recount",
            "output": "tables/06_antibody_normal_ihc_respiratory.tsv ; tables/07_antibody_cancer_ihc_patients.tsv ; tables/14_cancer_ihc_counts_by_antibody.tsv",
        },
        {
            "source": "HPA v25.1 ms_tissue.tsv / ms_tissue_sample_data.tsv",
            "url": URLS["ms_tissue.tsv.zip"],
            "what_was_kept": "TACSTD2/CLDN4 lung MS intensities",
            "output": "tables/08_ms_tissue_lung.tsv ; tables/09_ms_tissue_lung_replicates.tsv",
        },
        {
            "source": "HPA v25.1 cancer_cptac.tsv",
            "url": URLS["cancer_cptac.tsv.zip"],
            "what_was_kept": "Lung rows if present; explicit absence for CLDN4 lung",
            "output": "tables/10_cptac_hpa_lung.tsv",
        },
        {
            "source": "HPA v23 normal_tissue.tsv / pathology.tsv",
            "url": URLS["v23_normal_tissue.tsv.zip"],
            "what_was_kept": "Lung IHC comparison only; full archives not committed",
            "output": "tables/04_v23_vs_v25_lung_ihc.tsv",
        },
    ]
    write_tsv(out / "source_catalog.tsv", catalog, ["source", "url", "what_was_kept", "output"])

    # do not keep large caches in git; leave a note
    (cache / "README.txt").write_text(
        "Local download cache for scripts/w200_hpa_lung/download_and_extract.py. Not committed.\n",
        encoding="utf-8",
    )

    print(f"wrote tables under {tables}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
