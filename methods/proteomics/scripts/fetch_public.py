#!/usr/bin/env python3
"""Download public LUAD/LSCC gene-abundance / RNA tables.

Sources
-------
1. LinkedOmics CPTAC-LUAD / CPTAC-LSCC CCT matrices (gene symbols; TMT log2-ratio).
2. GDC open STAR gene-counts via HTTPS or the public S3 bucket
   ``s3://gdc-cptac-phs001287-2-open/{file_id}/`` (Ensembl ENSG00000184292 / ENSG00000189143).
3. Optional local CPTAC pan-cancer gene-abundance files (Ensembl rows), if the user
   already mirrored ``LUAD_*gene_abundance*.txt`` / ``LSCC_*gene_abundance*.txt``.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

LINKEDOMICS = {
    "LUAD_protein_tumor": "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct",
    "LUAD_protein_normal": "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_proteome_ratio_NArm_NORMAL.cct",
    "LUAD_rna_tumor": "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_rnaseq_uq_rpkm_log2_NArm_TUMOR.cct",
    "LUAD_clinical": "https://linkedomics.org/data_download/CPTAC-LUAD/HS_CPTAC_LUAD_cli.tsi",
    "LSCC_protein_tumor": "https://linkedomics.org/data_download/CPTAC-LSCC/HS_CPTAC_LSCC_2020_proteome_ratio_NArm_TUMOR.cct",
    "LSCC_protein_normal": "https://linkedomics.org/data_download/CPTAC-LSCC/HS_CPTAC_LSCC_2020_proteome_ratio_NArm_NORMAL.cct",
    "LSCC_rna_tumor": "https://linkedomics.org/data_download/CPTAC-LSCC/HS_CPTAC_LSCC_2020_rnaseq_uq_fpkm_log2_NArm_TUMOR.cct",
    "LSCC_clinical": "https://linkedomics.org/data_download/CPTAC-LSCC/HS_CPTAC_LSCC_2020_clinical_phenotypes_TUMOR.tsi",
    "LSCC_molecular": "https://linkedomics.org/data_download/CPTAC-LSCC/HS_CPTAC_LSCC_2020_molecular_phenotypes_TUMOR.tsi",
}

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"
S3_BUCKET = "s3://gdc-cptac-phs001287-2-open"
S3_HTTP = "https://gdc-cptac-phs001287-2-open.s3.amazonaws.com"

# CPTAC pan-cancer freeze (Ensembl rows). Often CDS/dbGaP-gated; scripts still parse the format.
PANCANCER_KEYS = {
    "LUAD": "LUAD/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LSCC": "LSCC/LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
}


def _download(url: str, dest: Path, timeout: int = 120) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "methods-proteomics/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r, open(tmp, "wb") as f:
        shutil.copyfileobj(r, f)
    tmp.replace(dest)
    return dest


def fetch_linkedomics(cache: Path) -> dict[str, Path]:
    out = {}
    for key, url in LINKEDOMICS.items():
        dest = cache / "linkedomics" / Path(url).name
        print(f"[linkedomics] {key} -> {dest.name}", flush=True)
        out[key] = _download(url, dest)
    return out


def gdc_query_lung_star(histology: str, size: int = 8) -> list[dict]:
    diagnosis = {
        "LUAD": ["Adenocarcinoma, NOS"],
        "LSCC": ["Squamous cell carcinoma, NOS"],
        "lung": ["Adenocarcinoma, NOS", "Squamous cell carcinoma, NOS"],
    }[histology]
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": ["CPTAC-3"]}},
            {"op": "in", "content": {"field": "cases.primary_site", "value": ["Bronchus and lung"]}},
            {"op": "in", "content": {"field": "cases.diagnoses.primary_diagnosis", "value": diagnosis}},
            {"op": "in", "content": {"field": "files.analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "files.access", "value": ["open"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": ["Primary Tumor"]}},
        ],
    }
    body = {
        "filters": filters,
        "fields": "file_id,file_name,file_size,cases.submitter_id,cases.diagnoses.primary_diagnosis,cases.samples.sample_type",
        "format": "JSON",
        "size": str(size),
    }
    req = urllib.request.Request(
        GDC_FILES,
        data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(req, timeout=60) as r:
        payload = json.load(r)
    hits = []
    for h in payload.get("data", {}).get("hits", []):
        case = (h.get("cases") or [{}])[0]
        hits.append(
            {
                "file_id": h["file_id"],
                "file_name": h.get("file_name"),
                "file_size": h.get("file_size"),
                "submitter_id": case.get("submitter_id"),
                "diagnosis": (case.get("diagnoses") or [{}])[0].get("primary_diagnosis"),
                "s3_uri": f"{S3_BUCKET}/{h['file_id']}/{h.get('file_name')}",
                "s3_http": f"{S3_HTTP}/{h['file_id']}/{h.get('file_name')}",
                "gdc_http": f"{GDC_DATA}/{h['file_id']}",
            }
        )
    return hits


def fetch_star_file(hit: dict, cache: Path) -> Path:
    dest = cache / "star" / f"{hit['submitter_id']}_{hit['file_id'][:8]}.star.tsv"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Prefer unsigned S3, then GDC HTTPS
    aws = shutil.which("aws")
    if aws:
        cmd = [aws, "s3", "cp", "--no-sign-request", hit["s3_uri"], str(dest)]
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=120)
            if dest.exists() and dest.stat().st_size > 0:
                return dest
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pass
    try:
        return _download(hit["s3_http"], dest)
    except Exception:
        return _download(hit["gdc_http"], dest)


def write_manifest(paths: dict, star_hits: list[dict], dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "linkedomics": {k: str(v) for k, v in paths.items()},
        "gdc_star": star_hits,
        "pancancer_s3_note": (
            "CPTAC pan-cancer gene-abundance TSVs use Ensembl row IDs "
            "(ENSG00000184292 / ENSG00000189143). The historical bucket "
            "cptac-pancancer-data is often access-gated; drop local copies "
            "matching *gene_abundance*.txt into --pancancer-dir."
        ),
        "pancancer_keys": PANCANCER_KEYS,
    }
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache-dir", type=Path, default=Path("methods/proteomics/demo/cache"))
    ap.add_argument("--star-per-histology", type=int, default=2)
    ap.add_argument("--skip-star", action="store_true")
    ap.add_argument("--manifest", type=Path, default=None)
    args = ap.parse_args(argv)

    cache = args.cache_dir
    paths = fetch_linkedomics(cache)
    star_hits: list[dict] = []
    if not args.skip_star:
        for histo in ("LUAD", "LSCC"):
            hits = gdc_query_lung_star(histo, size=args.star_per_histology)
            print(f"[gdc] {histo} STAR hits: {len(hits)}", flush=True)
            for h in hits:
                p = fetch_star_file(h, cache)
                h["local_path"] = str(p)
                star_hits.append(h)
                print(f"  {h['submitter_id']} {h['file_id']} -> {p}", flush=True)
    man = args.manifest or (cache / "fetch_manifest.json")
    write_manifest(paths, star_hits, man)
    print(f"manifest -> {man}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
