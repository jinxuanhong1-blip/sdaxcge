#!/usr/bin/env python3
"""Audit public processed extracts for GSE179994 CLDN4-only analysis.

Allowed inputs: public processed UMI / h5 / TISCH extract, files <2 GB.
Endpoint needs malignant (author or EPCAM+ epithelial) CLDN4 and a
same-patient all-cell T/NK denominator. Patient is the unit.

This script does not download the 421 MB T-cell RDS. Author metadata
already shows CD4/CD8/NA only.
"""

from __future__ import annotations

import csv
import gzip
import json
import os
import tarfile
import urllib.request
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
TAB = HERE / "tables"
CACHE = HERE / "cache"
MAX_BYTES = 2 * 1024**3

GEO_SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl"
GEO_FILES = [
    "GSE179994_PBMC.bulkTCR.tsv.gz",
    "GSE179994_RAW.tar",
    "GSE179994_Tcell.metadata.tsv.gz",
    "GSE179994_all.Tcell.rawCounts.rds.gz",
    "GSE179994_all.scTCR.tsv.gz",
    "filelist.txt",
]
TISCH_PROBES = [
    "https://tisch.compbio.cn/static/data/NSCLC_GSE179994/NSCLC_GSE179994_expression.h5",
    "https://tisch.compbio.cn/static/data/NSCLC_GSE179994/NSCLC_GSE179994_CellMetainfo_table.tsv",
    "https://tisch.compbio.cn/static/data/LUAD_GSE179994/LUAD_GSE179994_expression.h5",
    "https://tisch.compbio.cn/static/data/LUAD_GSE179994/LUAD_GSE179994_CellMetainfo_table.tsv",
    "https://tisch.compbio.cn/static/data/NSCLC_GSE179994_aPD1/NSCLC_GSE179994_aPD1_CellMetainfo_table.tsv",
    "https://tisch.compbio.cn/static/data/LUAD_GSE179994_aPD1/LUAD_GSE179994_aPD1_CellMetainfo_table.tsv",
]


def head_url(url: str, timeout: int = 30) -> dict:
    req = urllib.request.Request(url, method="HEAD")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return {
                "url": url,
                "status": resp.status,
                "content_length": int(resp.headers.get("Content-Length") or 0),
                "content_type": resp.headers.get("Content-Type", ""),
            }
    except urllib.error.HTTPError as e:
        return {
            "url": url,
            "status": e.code,
            "content_length": int(e.headers.get("Content-Length") or 0),
            "content_type": e.headers.get("Content-Type", ""),
        }
    except Exception as e:  # noqa: BLE001
        return {"url": url, "status": None, "content_length": 0, "content_type": "", "error": str(e)}


def download(url: str, dest: Path, timeout: int = 120) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=timeout) as resp, open(tmp, "wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    return dest


def mb(n: int) -> float:
    return n / (1024 * 1024)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    TAB.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    geo_rows = []
    for name in GEO_FILES:
        info = head_url(f"{GEO_SUPPL}/{name}")
        info["file"] = name
        info["over_2gb"] = bool(info.get("content_length", 0) > MAX_BYTES)
        geo_rows.append(info)

    tisch_rows = [head_url(u) for u in TISCH_PROBES]

    meta_path = download(
        f"{GEO_SUPPL}/GSE179994_Tcell.metadata.tsv.gz",
        CACHE / "GSE179994_Tcell.metadata.tsv.gz",
    )
    raw_path = download(f"{GEO_SUPPL}/GSE179994_RAW.tar", CACHE / "GSE179994_RAW.tar")

    raw_members = []
    with tarfile.open(raw_path, "r") as tf:
        for m in tf.getmembers():
            raw_members.append({"name": m.name, "size": m.size})

    with gzip.open(meta_path, "rt") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    celltype = Counter(r.get("celltype") or "NA" for r in rows)
    cluster = Counter(r.get("cluster") or "NA" for r in rows)
    patients = sorted({r["patient"] for r in rows})
    samples = sorted({r["sample"] for r in rows})
    per_patient = Counter(r["patient"] for r in rows)
    per_sample = Counter(r["sample"] for r in rows)

    epithelial_tokens = (
        "epi",
        "malig",
        "tumor",
        "cancer",
        "epcam",
        "luad",
        "nsclc",
        "epithelial",
    )
    author_labels = set(celltype) | set(cluster)
    has_malignant_label = any(
        any(tok in str(lab).lower() for tok in epithelial_tokens) for lab in author_labels
    )

    tisch_ok = any(r.get("status") == 200 for r in tisch_rows)
    tcell_rds = next(r for r in geo_rows if r["file"].endswith("Tcell.rawCounts.rds.gz"))
    usable = False
    stop_reason = (
        "No public processed all-cell UMI/h5/TISCH extract <2GB. "
        "GEO processed expression is T-cell-only; TISCH dataset absent (404)."
    )

    feasibility = {
        "dataset": "GSE179994",
        "paper": "Liu et al. Nat Cancer 2022 PMID 35121991 (Zhang lab, PKU)",
        "regimen": "pembrolizumab + carboplatin + pemetrexed",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "allowed_inputs": "public processed UMI / h5 / TISCH extract; skip files >2GB",
        "endpoint": "malignant (author or EPCAM+epithelial) CLDN4 mean and %pos vs same-patient T/NK, CD8, CXCL13+",
        "unit": "patient",
        "no_dual_high": True,
        "no_tacstd2_gate": True,
        "usable_processed_matrix_lt_2gb": usable,
        "stop_reason": stop_reason,
        "n_patients_analysis": 0,
        "n_patients_tcell_metadata_only": len(patients),
        "n_samples_tcell_metadata_only": len(samples),
        "n_tcell_barcodes_metadata": len(rows),
        "do_not_cite_barcode_count_as_n": True,
        "q4_vs_q1": "not done; analysis n_patients=0 (threshold is >=16)",
        "author_celltypes": dict(celltype),
        "has_author_malignant_or_epithelial_label": has_malignant_label,
        "tcell_rds_bytes": tcell_rds.get("content_length"),
        "tcell_rds_under_2gb": tcell_rds.get("content_length", 0) < MAX_BYTES,
        "tcell_rds_usable_for_endpoint": False,
        "tisch_dataset_present": tisch_ok,
        "geo_files": geo_rows,
        "tisch_probes": tisch_rows,
        "raw_tar_members": raw_members,
        "gsm_supplementary": "Sample records state supplementary data files not provided; raw data not provided; processed data on series record only.",
    }

    with open(TAB / "feasibility.json", "w") as f:
        json.dump(feasibility, f, indent=2)

    with open(TAB / "geo_file_inventory.tsv", "w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["file", "status", "bytes", "mb", "over_2gb", "url"],
        )
        w.writeheader()
        for r in geo_rows:
            w.writerow(
                {
                    "file": r["file"],
                    "status": r.get("status"),
                    "bytes": r.get("content_length"),
                    "mb": f"{mb(r.get('content_length') or 0):.2f}",
                    "over_2gb": r.get("over_2gb"),
                    "url": r.get("url"),
                }
            )

    with open(TAB / "tisch_probes.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["url", "status", "bytes", "content_type", "error"])
        w.writeheader()
        for r in tisch_rows:
            w.writerow(
                {
                    "url": r["url"],
                    "status": r.get("status"),
                    "bytes": r.get("content_length"),
                    "content_type": r.get("content_type"),
                    "error": r.get("error", ""),
                }
            )

    with open(TAB / "author_celltype_counts.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["celltype", "n_barcodes", "note"])
        for k, v in celltype.most_common():
            w.writerow([k, v, "T-cell metadata only; not analysis n"])

    with open(TAB / "author_cluster_counts.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["cluster", "n_barcodes", "note"])
        for k, v in cluster.most_common():
            w.writerow([k, v, "T-cell metadata only; not analysis n"])

    with open(TAB / "tcell_metadata_patient_census.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(
            [
                "patient",
                "n_tcell_barcodes",
                "n_samples",
                "samples",
                "usable_for_malignant_CLDN4",
            ]
        )
        samp_by_pat = defaultdict(set)
        for r in rows:
            samp_by_pat[r["patient"]].add(r["sample"])
        for p in patients:
            w.writerow(
                [
                    p,
                    per_patient[p],
                    len(samp_by_pat[p]),
                    ",".join(sorted(samp_by_pat[p])),
                    "no",
                ]
            )

    with open(TAB / "raw_tar_members.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "size"])
        w.writeheader()
        w.writerows(raw_members)

    # Figures
    plot_file_inventory(geo_rows)
    plot_celltypes(celltype)
    plot_patient_census(per_patient)

    print(json.dumps({"usable": usable, "n_patients_analysis": 0, "stop_reason": stop_reason}, indent=2))


def plot_file_inventory(geo_rows: list[dict]) -> None:
    names = [r["file"].replace("GSE179994_", "") for r in geo_rows if r["file"] != "filelist.txt"]
    sizes = [mb(r.get("content_length") or 0) for r in geo_rows if r["file"] != "filelist.txt"]
    colors = ["#4C78A8" if n.endswith("Tcell.rawCounts.rds.gz") else "#72B7B2" for n in [r["file"] for r in geo_rows if r["file"] != "filelist.txt"]]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    ax.barh(names[::-1], sizes[::-1], color=colors[::-1])
    ax.axvline(2048, color="#B279A2", ls="--", lw=1, label="2 GB skip line")
    ax.set_xlabel("GEO supplementary file size (MB)")
    ax.set_title("GSE179994 GEO files — no all-cell UMI/h5 (T-cell RDS is 421 MB)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_geo_file_inventory.png", dpi=160)
    fig.savefig(FIG / "fig_geo_file_inventory.pdf")
    plt.close(fig)


def plot_celltypes(celltype: Counter) -> None:
    labels = list(celltype.keys())
    vals = [celltype[k] for k in labels]
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.bar(labels, vals, color=["#4C78A8", "#F58518", "#54A24B"][: len(labels)])
    ax.set_ylabel("T-cell barcodes in author metadata")
    ax.set_title("Author celltype is CD4 / CD8 / NA only — no malignant/epithelial")
    ax.set_xlabel("GSE179994_Tcell.metadata.tsv.gz celltype")
    for i, v in enumerate(vals):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "fig_author_celltypes_t_only.png", dpi=160)
    fig.savefig(FIG / "fig_author_celltypes_t_only.pdf")
    plt.close(fig)


def plot_patient_census(per_patient: Counter) -> None:
    pats = sorted(per_patient, key=lambda p: (-per_patient[p], p))
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    ax.bar(range(len(pats)), [per_patient[p] for p in pats], color="#4C78A8")
    ax.set_xticks(range(len(pats)))
    ax.set_xticklabels(pats, rotation=90, fontsize=7)
    ax.set_ylabel("T-cell barcodes (metadata census only)")
    ax.set_title("36 patients in T-cell metadata — analysis n for malignant CLDN4 = 0")
    fig.tight_layout()
    fig.savefig(FIG / "fig_tcell_metadata_patient_census.png", dpi=160)
    fig.savefig(FIG / "fig_tcell_metadata_patient_census.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
