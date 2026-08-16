#!/usr/bin/env python3
"""Fetch TCGA-COAD RNA-seq (STAR - Counts) from the NCI GDC open-access API.

This builds a reproducible gene x sample TPM matrix plus a sample manifest
(sample type: tumor vs normal). Nothing here is synthetic: every value comes
from a GDC Gene Expression Quantification file (GENCODE v36, workflow
"STAR - Counts").

Outputs (under --outdir, default: data/):
  - coad_sample_manifest.tsv : one row per RNA-seq file with case/sample ids
                               and sample_type.
  - coad_gene_annotation.tsv : gene_id, gene_name, gene_type (GENCODE v36).
  - coad_star_tpm.parquet    : genes (rows) x aliquot barcodes (cols), TPM.
"""
from __future__ import annotations

import argparse
import io
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests

GDC_FILES = "https://api.gdc.cancer.gov/files"
GDC_DATA = "https://api.gdc.cancer.gov/data"

FILTERS = {
    "op": "and",
    "content": [
        {"op": "in", "content": {"field": "cases.project.project_id", "value": ["TCGA-COAD"]}},
        {"op": "in", "content": {"field": "data_category", "value": ["Transcriptome Profiling"]}},
        {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
        {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        {"op": "in", "content": {"field": "access", "value": ["open"]}},
    ],
}


def get_manifest() -> pd.DataFrame:
    """Return one row per RNA-seq file with sample metadata."""
    fields = [
        "file_id",
        "file_name",
        "cases.submitter_id",
        "cases.samples.submitter_id",
        "cases.samples.sample_type",
    ]
    rows: list[dict] = []
    size = 500
    frm = 0
    while True:
        params = {
            "filters": json.dumps(FILTERS),
            "format": "JSON",
            "fields": ",".join(fields),
            "size": str(size),
            "from": str(frm),
        }
        r = requests.get(GDC_FILES, params=params, timeout=120)
        r.raise_for_status()
        data = r.json()["data"]
        hits = data["hits"]
        for h in hits:
            case = (h.get("cases") or [{}])[0]
            samples = case.get("samples") or [{}]
            sample = samples[0]
            rows.append(
                {
                    "file_id": h["file_id"],
                    "file_name": h.get("file_name"),
                    "case_barcode": case.get("submitter_id"),
                    "sample_barcode": sample.get("submitter_id"),
                    "sample_type": sample.get("sample_type"),
                }
            )
        total = data["pagination"]["total"]
        frm += len(hits)
        if frm >= total or not hits:
            break
    df = pd.DataFrame(rows).drop_duplicates("file_id").reset_index(drop=True)
    return df


def _download_one(file_id: str, retries: int = 4) -> tuple[str, pd.DataFrame]:
    delay = 4
    last = None
    for _ in range(retries):
        try:
            r = requests.get(f"{GDC_DATA}/{file_id}", timeout=180)
            r.raise_for_status()
            df = pd.read_csv(
                io.BytesIO(r.content),
                sep="\t",
                comment="#",
                header=0,
            )
            # Drop STAR summary rows (N_unmapped, ...) which have no gene_id "ENSG".
            df = df[df["gene_id"].astype(str).str.startswith("ENSG")].copy()
            return file_id, df
        except Exception as exc:  # noqa: BLE001 - retry any transient failure
            last = exc
            time.sleep(delay)
            delay *= 2
    raise RuntimeError(f"failed to download {file_id}: {last}")


def build_matrix(manifest: pd.DataFrame, workers: int = 8) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Download every file and assemble a genes x sample TPM matrix."""
    tpm_cols: dict[str, pd.Series] = {}
    annotation: pd.DataFrame | None = None
    id_to_barcode = dict(zip(manifest["file_id"], manifest["sample_barcode"]))

    done = 0
    total = len(manifest)
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futures = {ex.submit(_download_one, fid): fid for fid in manifest["file_id"]}
        for fut in as_completed(futures):
            fid, df = fut.result()
            if annotation is None:
                annotation = df[["gene_id", "gene_name", "gene_type"]].reset_index(drop=True)
            s = pd.Series(df["tpm_unstranded"].values, index=df["gene_id"].values)
            barcode = id_to_barcode[fid]
            # If two files map to the same aliquot barcode, keep file_id-suffixed uniqueness.
            col = barcode if barcode not in tpm_cols else f"{barcode}|{fid[:8]}"
            tpm_cols[col] = s
            done += 1
            if done % 25 == 0 or done == total:
                print(f"  downloaded {done}/{total}", file=sys.stderr, flush=True)

    matrix = pd.DataFrame(tpm_cols)
    matrix = matrix.loc[annotation["gene_id"].values]
    return matrix, annotation


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="data")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=0, help="debug: limit number of files")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("Querying GDC file manifest for TCGA-COAD STAR - Counts ...", file=sys.stderr)
    manifest = get_manifest()
    if args.limit:
        manifest = manifest.head(args.limit)
    print(f"  {len(manifest)} files", file=sys.stderr)
    manifest.to_csv(outdir / "coad_sample_manifest.tsv", sep="\t", index=False)

    print("Downloading expression files ...", file=sys.stderr)
    matrix, annotation = build_matrix(manifest, workers=args.workers)
    annotation.to_csv(outdir / "coad_gene_annotation.tsv", sep="\t", index=False)
    matrix.to_parquet(outdir / "coad_star_tpm.parquet")

    print(
        f"Wrote matrix: {matrix.shape[0]} genes x {matrix.shape[1]} samples", file=sys.stderr
    )


if __name__ == "__main__":
    main()
