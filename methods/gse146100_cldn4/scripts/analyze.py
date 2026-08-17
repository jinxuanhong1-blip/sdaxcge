#!/usr/bin/env python3
"""CLDN4-only: epithelial (malignant proxy) vs T/NK on GSE146100.

Honest n
--------
GEO / TISCH / Zhang et al. (JITC 2021, PMID 33820821) record **1 patient**
and **3 nodules** (W1, W2, W3). This is not a 3-patient cohort. n=3 is
descriptive only. Spearman is not computed (need n>=5). Cell-level tests
are not used as evidence.

CLDN4-only
----------
TACSTD2 is extracted only as a present/absent sanity gene. It is never a
gate. Dual-high (TACSTD2+CLDN4) is not run.

Malignant
---------
TISCH Celltype (malignancy) has no Malignant label on this object
(Immune vs Stromal only). Primary epithelial set = TISCH major-lineage
Epithelial. Sensitivity = marker rule (any of EPCAM/KRT8/KRT18/KRT19 > 0
and PTPRC == 0) on the TISCH cells and on the GEO NormData matrix.
"""

from __future__ import annotations

import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

TISCH_H5 = DATA / "NSCLC_GSE146100_expression.h5"
TISCH_META = DATA / "NSCLC_GSE146100_CellMetainfo_table.tsv"
GEO = DATA / "GSE146100_NormData.txt.gz"

SAMPLE_META = {
    "W1": {
        "gsm": "GSM4365352",
        "response": "NR",
        "genotype": "EGFR L858R",
        "title": "Non-responded W1",
    },
    "W2": {
        "gsm": "GSM4365353",
        "response": "R",
        "genotype": "KRAS G12C",
        "title": "Responded W2",
    },
    "W3": {
        "gsm": "GSM4365354",
        "response": "NR",
        "genotype": "EGFR L858R/R77H",
        "title": "Non-responded W3",
    },
}

TNK_LINEAGES = ("CD4Tconv", "CD8T", "NK", "Treg")
EPI_MARKERS = ("EPCAM", "KRT8", "KRT18", "KRT19")
WANTED = (
    "CLDN4",
    "TACSTD2",
    "EPCAM",
    "PTPRC",
    "CD3D",
    "CD3E",
    "CD8A",
    "NKG7",
    "KRT8",
    "KRT18",
    "KRT19",
)
MIN_N_SPEARMAN = 5


def _as_str(arr) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode())
        else:
            out.append(str(x))
    return np.array(out, dtype=object)


def load_tisch_expr() -> tuple[pd.DataFrame, dict]:
    with h5py.File(TISCH_H5, "r") as f:
        g = f["matrix"]
        barcodes = _as_str(g["barcodes"][:])
        genes = _as_str(g["features/name"][:])
        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = [int(x) for x in g["shape"][:]]
    n_cells = len(barcodes)
    n_genes = len(genes)
    csr = sparse.csr_matrix((data, indices, indptr), shape=(n_cells, n_genes))
    present = [g for g in WANTED if g in set(genes)]
    missing = [g for g in WANTED if g not in set(genes)]
    table = {}
    for gene in present:
        i = int(np.where(genes == gene)[0][0])
        table[gene] = np.asarray(csr[:, i].todense()).ravel()
    expr = pd.DataFrame(table, index=barcodes)
    audit = {
        "path": str(TISCH_H5),
        "bytes": TISCH_H5.stat().st_size,
        "declared_shape": shape,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "genes_present": present,
        "genes_absent": missing,
        "normalization": "TISCH2 / MAESTRO log2(TPM/10+1); not re-normalized",
    }
    return expr, audit


def load_geo_expr() -> tuple[pd.DataFrame, dict]:
    wanted = set(WANTED)
    rows: dict[str, np.ndarray] = {}
    with gzip.open(GEO, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cells = header[1:]
        for line in handle:
            gene = line.split("\t", 1)[0]
            if gene in wanted:
                rows[gene] = np.array(line.rstrip("\n").split("\t")[1:], dtype=np.float32)
                if len(rows) == len(wanted):
                    break
    expr = pd.DataFrame(rows, index=cells)
    expr["sample"] = [c.split("_")[0] for c in expr.index]
    audit = {
        "path": str(GEO),
        "bytes": GEO.stat().st_size,
        "n_cells": int(expr.shape[0]),
        "n_genes_extracted": int(expr.shape[1] - 1),
        "genes_present": [g for g in WANTED if g in expr.columns],
        "genes_absent": [g for g in WANTED if g not in expr.columns],
        "normalization": "author Seurat log-norm (GSE146100_NormData.txt.gz)",
    }
    return expr, audit


def marker_epithelial(frame: pd.DataFrame) -> pd.Series:
    epi_any = np.zeros(len(frame), dtype=bool)
    for gene in EPI_MARKERS:
        if gene in frame.columns:
            epi_any |= frame[gene].to_numpy() > 0
    ptprc = frame["PTPRC"].to_numpy() if "PTPRC" in frame.columns else 0
    return pd.Series(epi_any & (ptprc == 0), index=frame.index)


def _stats(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    if len(values) == 0:
        return {
            "n_cells": 0,
            "mean": float("nan"),
            "median": float("nan"),
            "pct_pos": float("nan"),
        }
    return {
        "n_cells": int(len(values)),
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "pct_pos": float(100.0 * np.mean(values > 0)),
    }


def per_sample_tisch(joined: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample in ("W1", "W2", "W3"):
        sub = joined[joined["Sample"] == sample]
        epi = sub[sub["Celltype (major-lineage)"] == "Epithelial"]
        tnk = sub[sub["Celltype (major-lineage)"].isin(TNK_LINEAGES)]
        malig = sub[sub["Celltype (malignancy)"] == "Malignant"]
        marker = marker_epithelial(sub)
        meta = SAMPLE_META[sample]
        epi_s = _stats(epi["CLDN4"].to_numpy())
        marker_s = _stats(sub.loc[marker, "CLDN4"].to_numpy())
        tnk_s = _stats(tnk["CLDN4"].to_numpy())
        rows.append(
            {
                "sample": sample,
                "gsm": meta["gsm"],
                "response": meta["response"],
                "genotype": meta["genotype"],
                "patient": "Patient 1",
                "n_cells": int(len(sub)),
                "n_malignant_tisch": int(len(malig)),
                "n_epithelial": epi_s["n_cells"],
                "n_marker_epithelial": marker_s["n_cells"],
                "n_tnk": int(len(tnk)),
                "frac_tnk": float(len(tnk) / len(sub)) if len(sub) else float("nan"),
                "n_CD4Tconv": int((sub["Celltype (major-lineage)"] == "CD4Tconv").sum()),
                "n_CD8T": int((sub["Celltype (major-lineage)"] == "CD8T").sum()),
                "n_NK": int((sub["Celltype (major-lineage)"] == "NK").sum()),
                "n_Treg": int((sub["Celltype (major-lineage)"] == "Treg").sum()),
                "cldn4_epi_mean": epi_s["mean"],
                "cldn4_epi_median": epi_s["median"],
                "cldn4_epi_pctpos": epi_s["pct_pos"],
                "cldn4_marker_mean": marker_s["mean"],
                "cldn4_marker_median": marker_s["median"],
                "cldn4_marker_pctpos": marker_s["pct_pos"],
                "cldn4_tnk_mean": tnk_s["mean"],
                "cldn4_tnk_pctpos": tnk_s["pct_pos"],
                "epcam_epi_mean": float(epi["EPCAM"].mean()) if len(epi) else float("nan"),
                "ptprc_epi_mean": float(epi["PTPRC"].mean()) if len(epi) else float("nan"),
                "tacstd2_epi_mean": float(epi["TACSTD2"].mean()) if len(epi) else float("nan"),
                "tacstd2_used_as_gate": False,
            }
        )
    return pd.DataFrame(rows)


def per_sample_geo(geo: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sample in ("W1", "W2", "W3"):
        sub = geo[geo["sample"] == sample]
        marker = marker_epithelial(sub)
        tnk_m = (sub[["CD3D", "CD8A", "NKG7"]].max(axis=1) > 0) & (sub["EPCAM"] == 0)
        meta = SAMPLE_META[sample]
        marker_s = _stats(sub.loc[marker, "CLDN4"].to_numpy())
        rows.append(
            {
                "sample": sample,
                "gsm": meta["gsm"],
                "response": meta["response"],
                "genotype": meta["genotype"],
                "patient": "Patient 1",
                "n_cells": int(len(sub)),
                "n_marker_epithelial": marker_s["n_cells"],
                "n_tnk_marker": int(tnk_m.sum()),
                "frac_tnk_marker": float(tnk_m.mean()) if len(sub) else float("nan"),
                "cldn4_marker_mean": marker_s["mean"],
                "cldn4_marker_median": marker_s["median"],
                "cldn4_marker_pctpos": marker_s["pct_pos"],
            }
        )
    return pd.DataFrame(rows)


def lineage_cldn4(joined: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (sample, lin), sub in joined.groupby(
        ["Sample", "Celltype (major-lineage)"], observed=True
    ):
        s = _stats(sub["CLDN4"].to_numpy())
        rows.append(
            {
                "sample": sample,
                "lineage": lin,
                "n_cells": s["n_cells"],
                "cldn4_mean": s["mean"],
                "cldn4_median": s["median"],
                "cldn4_pctpos": s["pct_pos"],
            }
        )
    return pd.DataFrame(rows).sort_values(["sample", "lineage"])


def make_figures(tisch_tab: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    order = ["W1", "W2", "W3"]
    tab = tisch_tab.set_index("sample").loc[order]
    colors = {"W1": "#4C78A8", "W2": "#54A24B", "W3": "#E45756"}
    col = [colors[s] for s in order]
    labels = [f"{s}\n{tab.loc[s, 'response']}" for s in order]

    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.4))
    ax[0].bar(labels, tab["cldn4_epi_mean"], color=col)
    ax[0].set_ylabel("Mean CLDN4 (TISCH log2)")
    ax[0].set_title("Epithelial CLDN4")
    ax[1].bar(labels, tab["cldn4_epi_pctpos"], color=col)
    ax[1].set_ylabel("% CLDN4+")
    ax[1].set_title("Epithelial % CLDN4+")
    ax[2].bar(labels, 100 * tab["frac_tnk"], color=col)
    ax[2].set_ylabel("% T/NK of all cells")
    ax[2].set_title("T/NK fraction")
    fig.suptitle("GSE146100 · 1 patient / 3 nodules · descriptive only", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_bars_cldn4_tnk.png", dpi=150, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_bars_cldn4_tnk.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.6, 4.2))
    for s in order:
        ax.scatter(
            tab.loc[s, "cldn4_epi_mean"],
            100 * tab.loc[s, "frac_tnk"],
            s=90,
            color=colors[s],
            zorder=3,
        )
        ax.annotate(
            f"{s} ({tab.loc[s, 'response']})",
            (tab.loc[s, "cldn4_epi_mean"], 100 * tab.loc[s, "frac_tnk"]),
            textcoords="offset points",
            xytext=(6, 6),
            fontsize=9,
        )
    ax.set_xlabel("TISCH epithelial mean CLDN4")
    ax.set_ylabel("T/NK fraction (% of cells)")
    ax.set_title("n=3 nodules, 1 patient — no Spearman")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_scatter_cldn4_vs_tnk.png", dpi=150, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_scatter_cldn4_vs_tnk.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    for path in (TISCH_H5, TISCH_META, GEO):
        if not path.exists():
            raise SystemExit(f"missing {path}; run scripts/download.py")

    expr, h5_audit = load_tisch_expr()
    meta = pd.read_csv(TISCH_META, sep="\t")
    if not meta["Cell"].isin(expr.index).all():
        raise SystemExit("TISCH Cell ids do not match h5 barcodes")
    joined = meta.set_index("Cell").join(expr, how="inner")
    geo, geo_audit = load_geo_expr()

    tisch_tab = per_sample_tisch(joined)
    geo_tab = per_sample_geo(geo)
    lin_tab = lineage_cldn4(joined)
    tisch_tab.to_csv(TABLES / "per_nodule_tisch.tsv", sep="\t", index=False)
    geo_tab.to_csv(TABLES / "per_nodule_geo_marker.tsv", sep="\t", index=False)
    lin_tab.to_csv(TABLES / "cldn4_by_lineage.tsv", sep="\t", index=False)

    lineage_counts = (
        joined.groupby(["Sample", "Celltype (major-lineage)"])
        .size()
        .rename("n_cells")
        .reset_index()
    )
    lineage_counts.to_csv(TABLES / "lineage_counts.tsv", sep="\t", index=False)

    sample_meta = pd.DataFrame(
        [
            {
                "sample": s,
                **m,
                "patient": "Patient 1",
                "n_patients_in_series": 1,
                "n_nodules_in_series": 3,
            }
            for s, m in SAMPLE_META.items()
        ]
    )
    sample_meta.to_csv(TABLES / "sample_metadata.tsv", sep="\t", index=False)

    keep = [
        "Sample",
        "Patient",
        "Tissue",
        "Cluster",
        "Celltype (malignancy)",
        "Celltype (major-lineage)",
        "Celltype (minor-lineage)",
        "CLDN4",
        "EPCAM",
        "PTPRC",
        "TACSTD2",
    ]
    joined[keep].to_csv(TABLES / "cell_annotation_tisch.tsv", sep="\t")

    n_units = int(tisch_tab.shape[0])
    spearman_note = (
        f"Spearman not computed: n_nodules={n_units} < {MIN_N_SPEARMAN}. "
        "n_patients=1. Descriptive only."
    )
    contrast = pd.DataFrame(
        [
            {
                "contrast": "epithelial_CLDN4_mean_vs_fracTNK",
                "unit": "nodule",
                "n_units": n_units,
                "n_patients": 1,
                "rho": "",
                "p": "",
                "note": spearman_note,
            },
            {
                "contrast": "epithelial_CLDN4_pctpos_vs_fracTNK",
                "unit": "nodule",
                "n_units": n_units,
                "n_patients": 1,
                "rho": "",
                "p": "",
                "note": spearman_note,
            },
        ]
    )
    contrast.to_csv(TABLES / "vs_tnk.tsv", sep="\t", index=False)

    one_row = {
        "dataset": "GSE146100",
        "pmid": 33820821,
        "n_patients": 1,
        "n_nodules": 3,
        "framed_as_n3_patients": True,
        "honest_n": "1 patient / 3 nodules (not a 3-patient cohort)",
        "matrix": "TISCH2 h5 + CellMetainfo; GEO NormData sensitivity",
        "malignant_label": "none (TISCH Immune/Stromal only)",
        "epithelial_proxy_n": int(tisch_tab["n_epithelial"].sum()),
        "tnk_definition": "TISCH major-lineage CD4Tconv+CD8T+NK+Treg",
        "dual_high_run": False,
        "tacstd2_used_as_gate": False,
        "spearman_computed": False,
        "W1_cldn4_epi_mean": float(tisch_tab.loc[tisch_tab["sample"] == "W1", "cldn4_epi_mean"].iloc[0]),
        "W2_cldn4_epi_mean": float(tisch_tab.loc[tisch_tab["sample"] == "W2", "cldn4_epi_mean"].iloc[0]),
        "W3_cldn4_epi_mean": float(tisch_tab.loc[tisch_tab["sample"] == "W3", "cldn4_epi_mean"].iloc[0]),
        "W1_frac_tnk": float(tisch_tab.loc[tisch_tab["sample"] == "W1", "frac_tnk"].iloc[0]),
        "W2_frac_tnk": float(tisch_tab.loc[tisch_tab["sample"] == "W2", "frac_tnk"].iloc[0]),
        "W3_frac_tnk": float(tisch_tab.loc[tisch_tab["sample"] == "W3", "frac_tnk"].iloc[0]),
    }
    pd.DataFrame([one_row]).to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    make_figures(tisch_tab)

    sanity = {
        "n_patients": 1,
        "n_nodules": 3,
        "patient_ids": sorted(joined["Patient"].astype(str).unique().tolist()),
        "sample_ids": ["W1", "W2", "W3"],
        "tisch_n_cells": int(joined.shape[0]),
        "geo_n_cells": int(geo.shape[0]),
        "tisch_n_genes": h5_audit["n_genes"],
        "tisch_malignant_n": int((joined["Celltype (malignancy)"] == "Malignant").sum()),
        "tisch_epithelial_n": int((joined["Celltype (major-lineage)"] == "Epithelial").sum()),
        "tisch_tnk_n": int(joined["Celltype (major-lineage)"].isin(TNK_LINEAGES).sum()),
        "malignancy_labels": joined["Celltype (malignancy)"].value_counts().to_dict(),
        "major_lineage_counts": joined["Celltype (major-lineage)"].value_counts().to_dict(),
        "cldn4_present": "CLDN4" in joined.columns,
        "tacstd2_present": "TACSTD2" in joined.columns,
        "tacstd2_used_as_gate": False,
        "dual_high_run": False,
        "spearman_computed": False,
        "response_confounded_with_genotype": True,
        "pretreatment_samples": False,
        "author_celltype_labels": "TISCH2 CellMetainfo (not author GEO labels)",
        "cnv_available": False,
    }
    (ROOT / "sanity.json").write_text(json.dumps(sanity, indent=2) + "\n")

    summary = {
        "dataset": "GSE146100",
        "pmid": 33820821,
        "task": "ADDITIVE CLDN4-only: epithelial (malignant proxy) vs T/NK",
        "honest_verdict": (
            "n_patients=1, n_nodules=3. This is not a 3-patient cohort. "
            "TISCH has no Malignant call; epithelial is the proxy. "
            "CLDN4 is epithelial-restricted. The three nodule points are "
            "W3-high CLDN4 / low T/NK, W1 mid CLDN4 / high T/NK, W2 low-mid "
            "CLDN4 / mid T/NK. Spearman is not computed. Dual-high was not run. "
            "No patient-level or response-level claim."
        ),
        "n_patients": 1,
        "n_nodules": 3,
        "framed_as_n3_patients": (
            "Task text said n=3 patients. Public record is 1 patient / 3 nodules."
        ),
        "malignant_definition": (
            "TISCH Epithelial (no Malignant label). Marker-rule sensitivity "
            "on TISCH and GEO NormData."
        ),
        "tnk_definition": "TISCH CD4Tconv + CD8T + NK + Treg",
        "dual_high_run": False,
        "spearman_computed": False,
        "per_nodule_tisch": tisch_tab.to_dict(orient="records"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
    }
    (ROOT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (ROOT / "audit.json").write_text(
        json.dumps(
            {
                "tisch_h5": h5_audit,
                "geo_normdata": geo_audit,
                "skipped": [
                    "SRA/FASTQ",
                    "inferCNV / CopyKAT",
                    "Spearman / MWU as inferential tests",
                    "dual-high TACSTD2+CLDN4 gate",
                    "TACSTD2 as a CLDN4 gate",
                ],
            },
            indent=2,
        )
        + "\n"
    )

    manifest = pd.DataFrame(
        [
            {
                "file": str(GEO),
                "bytes": GEO.stat().st_size,
                "used": True,
                "note": "GEO author Seurat log-norm; marker-rule sensitivity; not committed",
            },
            {
                "file": str(TISCH_H5),
                "bytes": TISCH_H5.stat().st_size,
                "used": True,
                "note": "TISCH2 expression.h5 <2GB; primary CLDN4 values; not committed",
            },
            {
                "file": str(TISCH_META),
                "bytes": TISCH_META.stat().st_size,
                "used": True,
                "note": "TISCH2 CellMetainfo; 1 patient / 3 samples; no Malignant label",
            },
        ]
    )
    manifest.to_csv(ROOT / "file_manifest.tsv", sep="\t", index=False)
    print("Done.")
    print((ROOT / "summary.json").read_text())


if __name__ == "__main__":
    main()
