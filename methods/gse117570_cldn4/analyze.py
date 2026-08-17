#!/usr/bin/env python3
"""ADDITIVE CLDN4-only: GSE117570 malignant CLDN4 vs same-patient T/NK.

Public TISCH2 NSCLC_GSE117570 expression.h5 + CellMetainfo (Song et al. 2019,
Cancer Medicine, PMID 31033233). Patient is the unit. No dual-high gate.
No TACSTD2 claim. Spearman is not computed unless n>=5 eligible patients.

Downloads stay under $GSE117570_CLDN4_DATA (default /tmp/gse117570_cldn4)
and are not committed.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE117570_CLDN4_DATA", "/tmp/gse117570_cldn4"))

TISCH_BASE = "https://tisch.compbio.cn/static/data/NSCLC_GSE117570"
H5_URL = f"{TISCH_BASE}/NSCLC_GSE117570_expression.h5"
META_URL = f"{TISCH_BASE}/NSCLC_GSE117570_CellMetainfo_table.tsv"
GEO_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117570/"
    "matrix/GSE117570_series_matrix.txt.gz"
)

TNK_LINEAGES = {"CD8T", "CD4Tconv", "NK", "CD8Tex", "Treg", "Tprolif", "TMKI67", "Tcell", "NKT", "ILC"}
MALIGNANT = {"Malignant"}
MIN_MAL = 20
MIN_TNK = 20
MIN_TNK_PAIRED = 5
MIN_N_SPEARMAN = 5

# CLDN4-only. Companion genes are audit / restriction controls, not dual-high.
GENES = [
    "CLDN4",
    "EPCAM",
    "KRT19",
    "KRT7",
    "NAPSA",
    "TACSTD2",
    "PTPRC",
    "CD3D",
    "CD8A",
    "NKG7",
    "NKX2-1",
    "TP63",
]

CLINICAL = {
    "P1": {
        "histology": "LUAD",
        "tnm": "pT1b pN1",
        "sex": "Male",
        "smoking": "former (quit 1993)",
        "gsm_tumor": "GSM3304007",
    },
    "P2": {
        "histology": "LUSC",
        "tnm": "pT2a pN0",
        "sex": "Female",
        "smoking": "former (quit Oct 2017)",
        "gsm_tumor": "GSM3304009",
    },
    "P3": {
        "histology": "LUAD",
        "tnm": "pT2a pNx",
        "sex": "Female",
        "smoking": "former",
        "gsm_tumor": "GSM3304011",
    },
    "P4": {
        "histology": "LUAD",
        "tnm": "pT1a pN0",
        "sex": "Female",
        "smoking": "unknown",
        "gsm_tumor": "GSM3304013",
    },
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def curl_download(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return {"url": url, "path": str(dest), "bytes": dest.stat().st_size, "skipped": True, "ok": True}
    cmd = ["curl", "-fL", "--retry", "3", "--retry-delay", "4", "-o", str(dest), url]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        return {"url": url, "path": str(dest), "ok": False, "error": str(e)}
    if dest.exists() and dest.stat().st_size > 0:
        return {"url": url, "path": str(dest), "bytes": dest.stat().st_size, "skipped": False, "ok": True}
    return {"url": url, "path": str(dest), "ok": False, "error": "empty file"}


def _as_str(arr) -> np.ndarray:
    return np.array([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr], dtype=object)


def load_tisch_h5(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, dict]:
    audit: dict = {"path": str(path), "bytes": path.stat().st_size}
    with h5py.File(path, "r") as f:
        if "matrix" in f and all(k in f["matrix"] for k in ("data", "indices", "indptr")):
            grp = f["matrix"]
        else:
            raise ValueError(f"unexpected h5 keys: {list(f.keys())}")
        genes = _as_str(grp["features"]["name"][:]) if "features" in grp else None
        barcodes = _as_str(grp["barcodes"][:])
        shape = tuple(int(x) for x in grp["shape"][:])
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
    n_genes, n_cells = int(shape[0]), int(shape[1])
    csc = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))
    audit.update({"n_genes": n_genes, "n_cells": n_cells, "shape": [n_genes, n_cells]})
    present = [g for g in wanted if g in set(genes)]
    audit["genes_present"] = present
    audit["genes_absent"] = [g for g in wanted if g not in set(genes)]
    gene_index = {g: int(np.where(genes == g)[0][0]) for g in present}
    table = {"barcode": barcodes}
    for g, i in gene_index.items():
        table[g] = np.asarray(csc.getrow(i).todense()).ravel()
    return pd.DataFrame(table).set_index("barcode"), audit


def parse_geo_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as f:
        lines = f.read().splitlines()
    titles = None
    rows = []
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [x.strip('"') for x in line.split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession") and titles:
            acc = [x.strip('"') for x in line.split("\t")[1:]]
            rows.append(("gsm", acc))
        elif line.startswith("!Sample_characteristics_ch1") and titles:
            vals = [x.strip('"') for x in line.split("\t")[1:]]
            key = vals[0].split(":", 1)[0].strip() if vals and ":" in vals[0] else "characteristic"
            cleaned = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
            rows.append((key, cleaned))
    if not titles:
        return pd.DataFrame()
    out = {"sample": titles}
    for k, v in rows:
        out[k] = v
    return pd.DataFrame(out)


def mean_pct(v: np.ndarray) -> tuple[float, float]:
    if len(v) == 0:
        return float("nan"), float("nan")
    return float(np.mean(v)), float(np.mean(v > 0) * 100.0)


def fmt_p(p: float) -> str:
    if p != p:
        return "—"
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3f}"


def fmt_rho(x: float) -> str:
    if x != x:
        return "—"
    return f"{x:+.3f}"


def download_all() -> dict:
    DATA.mkdir(parents=True, exist_ok=True)
    files = {
        "h5": (H5_URL, DATA / "NSCLC_GSE117570_expression.h5"),
        "meta": (META_URL, DATA / "NSCLC_GSE117570_CellMetainfo_table.tsv"),
        "geo": (GEO_MATRIX_URL, DATA / "GSE117570_series_matrix.txt.gz"),
    }
    prov = {"generated_at": datetime.now(timezone.utc).isoformat(), "files": {}}
    for key, (url, dest) in files.items():
        print(f">> {key} {url}", flush=True)
        rec = curl_download(url, dest)
        if rec.get("ok") and dest.exists():
            rec["sha256"] = sha256_file(dest)
        prov["files"][key] = rec
        print(f"   ok={rec.get('ok')} bytes={rec.get('bytes')}", flush=True)
    return prov


def analyze() -> dict:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    prov = download_all()
    h5_path = DATA / "NSCLC_GSE117570_expression.h5"
    meta_path = DATA / "NSCLC_GSE117570_CellMetainfo_table.tsv"
    geo_path = DATA / "GSE117570_series_matrix.txt.gz"
    if not h5_path.exists() or not meta_path.exists():
        raise SystemExit("TISCH files missing")

    meta = pd.read_csv(meta_path, sep="\t")
    cell_col = meta.columns[0]
    meta = meta.set_index(cell_col)
    meta.index = meta.index.astype(str)
    lin_col = "Celltype (major-lineage)"
    meta["_lineage"] = meta[lin_col].astype(str).str.strip()
    meta["_is_mal"] = meta["_lineage"].isin(MALIGNANT)
    meta["_is_tnk"] = meta["_lineage"].isin(TNK_LINEAGES)
    meta["_is_cd8"] = meta["_lineage"].eq("CD8T")
    meta["_is_cd4"] = meta["_lineage"].eq("CD4Tconv")
    meta["_is_nk"] = meta["_lineage"].eq("NK")

    expr, h5_audit = load_tisch_h5(h5_path, GENES)
    common = meta.index.intersection(expr.index)
    meta = meta.loc[common]
    expr = expr.loc[common]
    for g in h5_audit["genes_present"]:
        meta[g] = expr[g].to_numpy()

    geo = parse_geo_matrix(geo_path) if geo_path.exists() else pd.DataFrame()

    # Detectability: any CLDN4 > 0 in the patient (all tissues).
    det_rows = []
    for pat, sub in meta.groupby("Patient", sort=True):
        det_rows.append(
            {
                "patient": pat,
                "n_cells": int(len(sub)),
                "CLDN4_n_pos": int((sub["CLDN4"] > 0).sum()) if "CLDN4" in sub else 0,
                "CLDN4_max": float(sub["CLDN4"].max()) if "CLDN4" in sub else float("nan"),
                "EPCAM_n_pos": int((sub["EPCAM"] > 0).sum()) if "EPCAM" in sub else 0,
                "TACSTD2_n_pos": int((sub["TACSTD2"] > 0).sum()) if "TACSTD2" in sub else 0,
                "cldn4_detectable": bool((sub["CLDN4"] > 0).any()) if "CLDN4" in sub else False,
            }
        )
    detect = pd.DataFrame(det_rows)

    # Per-sample composition + scores
    sample_rows = []
    for (pat, source, sample), sub in meta.groupby(["Patient", "Source", "Sample"], sort=True):
        mal = sub.loc[sub["_is_mal"]]
        tnk = sub.loc[sub["_is_tnk"]]
        cldn4_mal_mean, cldn4_mal_pct = mean_pct(mal["CLDN4"].to_numpy() if "CLDN4" in mal else np.array([]))
        cldn4_tnk_mean, cldn4_tnk_pct = mean_pct(tnk["CLDN4"].to_numpy() if "CLDN4" in tnk and len(tnk) else np.array([]))
        epcam_mal_mean, epcam_mal_pct = mean_pct(mal["EPCAM"].to_numpy() if "EPCAM" in mal else np.array([]))
        krt19_mal_mean, krt19_mal_pct = mean_pct(mal["KRT19"].to_numpy() if "KRT19" in mal else np.array([]))
        napsa_mal_mean, napsa_mal_pct = mean_pct(mal["NAPSA"].to_numpy() if "NAPSA" in mal else np.array([]))
        sample_rows.append(
            {
                "patient": pat,
                "source": source,
                "sample": sample,
                "n_cells": int(len(sub)),
                "n_malignant": int(sub["_is_mal"].sum()),
                "n_tnk": int(sub["_is_tnk"].sum()),
                "n_cd8": int(sub["_is_cd8"].sum()),
                "n_cd4": int(sub["_is_cd4"].sum()),
                "n_nk": int(sub["_is_nk"].sum()),
                "frac_tnk": float(sub["_is_tnk"].mean()) if len(sub) else float("nan"),
                "frac_cd8": float(sub["_is_cd8"].mean()) if len(sub) else float("nan"),
                "frac_nk": float(sub["_is_nk"].mean()) if len(sub) else float("nan"),
                "CLDN4_mal_mean": cldn4_mal_mean,
                "CLDN4_mal_pctpos": cldn4_mal_pct,
                "CLDN4_tnk_mean": cldn4_tnk_mean,
                "CLDN4_tnk_pctpos": cldn4_tnk_pct,
                "EPCAM_mal_mean": epcam_mal_mean,
                "EPCAM_mal_pctpos": epcam_mal_pct,
                "KRT19_mal_mean": krt19_mal_mean,
                "KRT19_mal_pctpos": krt19_mal_pct,
                "NAPSA_mal_mean": napsa_mal_mean,
                "NAPSA_mal_pctpos": napsa_mal_pct,
            }
        )
    samples = pd.DataFrame(sample_rows)

    # Patient-level primary = Tumor only
    tumor = samples.loc[samples["source"] == "Tumor"].copy()
    tumor = tumor.merge(detect[["patient", "cldn4_detectable", "CLDN4_n_pos", "CLDN4_max"]], on="patient", how="left")
    for k, col in [
        ("histology", "histology"),
        ("tnm", "tnm"),
        ("sex", "sex"),
        ("smoking", "smoking"),
        ("gsm_tumor", "gsm_tumor"),
    ]:
        tumor[col] = tumor["patient"].map(lambda p: CLINICAL[p][k])
    tumor["eligible_spearman"] = (
        (tumor["n_malignant"] >= MIN_MAL)
        & (tumor["n_tnk"] >= MIN_TNK)
        & tumor["cldn4_detectable"]
    )
    tumor["eligible_paired"] = (
        (tumor["n_malignant"] >= MIN_MAL)
        & (tumor["n_tnk"] >= MIN_TNK_PAIRED)
        & tumor["cldn4_detectable"]
    )
    tumor["mal_gt_tnk"] = tumor["CLDN4_mal_mean"] > tumor["CLDN4_tnk_mean"]
    tumor["note"] = ""
    tumor.loc[~tumor["cldn4_detectable"], "note"] = "CLDN4 empty in all patient cells (TISCH h5)"
    tumor.loc[tumor["cldn4_detectable"] & (tumor["n_tnk"] < MIN_TNK), "note"] = (
        f"T/NK n<{MIN_TNK} in tumor; fraction underpowered"
    )
    tumor.loc[tumor["eligible_spearman"], "note"] = "eligible for Spearman (n still <5 globally)"

    # Lineage restriction (all tissues, and tumor-only)
    lin_rows = []
    for source, sub0 in [("all", meta), ("Tumor", meta.loc[meta["Source"] == "Tumor"])]:
        for lin, sub in sub0.groupby("_lineage"):
            mean, pct = mean_pct(sub["CLDN4"].to_numpy())
            lin_rows.append(
                {
                    "source": source,
                    "lineage": lin,
                    "n_cells": int(len(sub)),
                    "CLDN4_mean": mean,
                    "CLDN4_pctpos": pct,
                }
            )
    lineage = pd.DataFrame(lin_rows)

    # Stats
    det_tumor = tumor.loc[tumor["cldn4_detectable"]].copy()
    elig = tumor.loc[tumor["eligible_spearman"]].copy()
    paired = tumor.loc[tumor["eligible_paired"]].copy()
    paired_any_tnk = tumor.loc[tumor["cldn4_detectable"] & (tumor["n_tnk"] >= 1) & (tumor["n_malignant"] >= MIN_MAL)].copy()

    stats_rows = []

    def add_stat(name, n, rho=np.nan, p=np.nan, note="", extra=None):
        rec = {"contrast": name, "n": int(n), "rho": rho, "p": p, "ok": bool(n >= MIN_N_SPEARMAN and rho == rho), "note": note}
        if extra:
            rec.update(extra)
        stats_rows.append(rec)

    if len(elig) >= MIN_N_SPEARMAN:
        rho, p = stats.spearmanr(elig["CLDN4_mal_mean"], elig["frac_tnk"])
        add_stat("CLDN4_mal_mean_vs_fracTNK_eligible", len(elig), float(rho), float(p))
    else:
        add_stat(
            "CLDN4_mal_mean_vs_fracTNK_eligible",
            len(elig),
            note=f"n={len(elig)}<{MIN_N_SPEARMAN}; Spearman not computed",
        )

    if len(det_tumor) >= MIN_N_SPEARMAN:
        rho, p = stats.spearmanr(det_tumor["CLDN4_mal_mean"], det_tumor["frac_tnk"])
        add_stat("CLDN4_mal_mean_vs_fracTNK_detectable", len(det_tumor), float(rho), float(p), "includes T/NK-thin patients")
    else:
        add_stat(
            "CLDN4_mal_mean_vs_fracTNK_detectable",
            len(det_tumor),
            note=f"n={len(det_tumor)}<{MIN_N_SPEARMAN}; Spearman not computed (P3/P4 T/NK-thin)",
        )

    # Paired Wilcoxon: malignant CLDN4 vs T/NK CLDN4, same patient
    for label, df in (
        ("paired_mal_vs_tnk_CLDN4_tnk>=5", paired),
        ("paired_mal_vs_tnk_CLDN4_tnk>=1", paired_any_tnk),
    ):
        if len(df) >= 3:
            try:
                w = stats.wilcoxon(df["CLDN4_mal_mean"], df["CLDN4_tnk_mean"], alternative="greater", zero_method="wilcox")
                p_w = float(w.pvalue)
            except ValueError:
                p_w = float("nan")
        else:
            p_w = float("nan")
        add_stat(
            label,
            len(df),
            p=p_w,
            note="Wilcoxon signed-rank malignant > T/NK; n<6 is descriptive",
            extra={
                "n_mal_gt_tnk": int(df["mal_gt_tnk"].sum()) if len(df) else 0,
                "median_mal": float(df["CLDN4_mal_mean"].median()) if len(df) else float("nan"),
                "median_tnk": float(df["CLDN4_tnk_mean"].median()) if len(df) else float("nan"),
            },
        )

    stats_df = pd.DataFrame(stats_rows)

    # Lineage composition for figures
    tumor_cells = meta.loc[meta["Source"] == "Tumor"].copy()

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "dataset": "GSE117570",
        "tisch_id": "NSCLC_GSE117570",
        "pmid": "31033233",
        "paper": "Song et al. 2019, Cancer Medicine",
        "unit": "patient",
        "gene": "CLDN4",
        "dual_high": False,
        "n_patients_catalog": 4,
        "n_patients_cldn4_detectable": int(detect["cldn4_detectable"].sum()),
        "n_patients_tumor_mal>=20": int((tumor["n_malignant"] >= MIN_MAL).sum()),
        "n_patients_eligible_spearman": int(tumor["eligible_spearman"].sum()),
        "n_patients_eligible_paired_tnk>=5": int(tumor["eligible_paired"].sum()),
        "n_patients_paired_tnk>=1_detectable": int(len(paired_any_tnk)),
        "n_cells_aligned": int(len(meta)),
        "n_tumor_cells": int((meta["Source"] == "Tumor").sum()),
        "n_malignant_tumor": int(((meta["Source"] == "Tumor") & meta["_is_mal"]).sum()),
        "n_tnk_tumor": int(((meta["Source"] == "Tumor") & meta["_is_tnk"]).sum()),
        "p1_cldn4_empty": bool(not detect.loc[detect.patient == "P1", "cldn4_detectable"].iloc[0]),
        "spearman_computed": bool(stats_df.loc[stats_df.contrast == "CLDN4_mal_mean_vs_fracTNK_eligible", "ok"].iloc[0]),
        "h5": h5_audit,
        "provenance": {k: {kk: vv for kk, vv in rec.items() if kk != "sha256"} | {"sha256": rec.get("sha256")} for k, rec in prov["files"].items()},
        "geo_traits_n": int(len(geo)),
        "ici_labels": False,
        "treatment": "treatment-naive",
        "statistics": stats_df.to_dict(orient="records"),
        "patients": tumor.to_dict(orient="records"),
    }

    detect.to_csv(TABLES / "detectability.tsv", sep="\t", index=False)
    samples.to_csv(TABLES / "per_sample.tsv", sep="\t", index=False)
    tumor.to_csv(TABLES / "per_patient_tumor.tsv", sep="\t", index=False)
    lineage.to_csv(TABLES / "cldn4_by_lineage.tsv", sep="\t", index=False)
    stats_df.to_csv(TABLES / "stats.tsv", sep="\t", index=False)
    if len(geo):
        geo.to_csv(TABLES / "geo_sample_traits.tsv", sep="\t", index=False)
    one = pd.DataFrame(
        [
            {
                "dataset": "GSE117570 TISCH NSCLC",
                "unit": "patient",
                "n_catalog": 4,
                "n_CLDN4_detectable": int(detect["cldn4_detectable"].sum()),
                "n_eligible_spearman": int(tumor["eligible_spearman"].sum()),
                "CLDN4_vs_fracTNK_rho": "",
                "CLDN4_vs_fracTNK_p": "",
                "note": "Spearman not computed (eligible n=1; detectable n=3). P1 CLDN4 empty in TISCH h5.",
                "paired_mal_gt_tnk": f"{int(paired_any_tnk['mal_gt_tnk'].sum())}/{len(paired_any_tnk)}",
                "dual_high": False,
                "ici_labels": False,
            }
        ]
    )
    one.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (TABLES / "download_provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    make_figures(tumor, samples, lineage, tumor_cells, detect)
    write_finding(tumor, samples, lineage, detect, stats_df, summary, h5_audit, prov)
    print("wrote", HERE / "FINDING.md")
    return summary


def _save(fig, name: str) -> None:
    fig.tight_layout()
    fig.savefig(FIGURES / f"{name}.png", dpi=160)
    fig.savefig(FIGURES / f"{name}.pdf")
    plt.close(fig)


def make_figures(tumor: pd.DataFrame, samples: pd.DataFrame, lineage: pd.DataFrame, tumor_cells: pd.DataFrame, detect: pd.DataFrame) -> None:
    patients = ["P1", "P2", "P3", "P4"]
    colors = {"P1": "#7f7f7f", "P2": "#1f77b4", "P3": "#d62728", "P4": "#2ca02c"}
    t = tumor.set_index("patient").loc[patients]

    # fig1 paired malignant vs T/NK CLDN4
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.array([0, 1])
    for p in patients:
        y = [t.loc[p, "CLDN4_mal_mean"], t.loc[p, "CLDN4_tnk_mean"]]
        ls = "--" if not t.loc[p, "cldn4_detectable"] else "-"
        ax.plot(x, y, ls, color=colors[p], marker="o", label=f"{p} nTNK={int(t.loc[p,'n_tnk'])}")
    ax.set_xticks([0, 1], ["Malignant", "T/NK"])
    ax.set_ylabel("CLDN4 mean (TISCH log2(TPM/10+1))")
    ax.set_title("Same-patient CLDN4: malignant vs T/NK (tumor)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(-0.2, 1.2)
    _save(fig, "fig1_paired_malig_vs_tnk")

    # fig2 scatter vs fraction
    fig, ax = plt.subplots(figsize=(6.2, 4.4))
    for p in patients:
        m = "x" if not t.loc[p, "cldn4_detectable"] else ("o" if t.loc[p, "n_tnk"] >= MIN_TNK else "s")
        ax.scatter(t.loc[p, "CLDN4_mal_mean"], t.loc[p, "frac_tnk"], s=90, c=colors[p], marker=m, zorder=3)
        ax.annotate(p, (t.loc[p, "CLDN4_mal_mean"], t.loc[p, "frac_tnk"]), textcoords="offset points", xytext=(6, 4), fontsize=9)
    ax.set_xlabel("Malignant CLDN4 mean")
    ax.set_ylabel("Tumor T/NK fraction")
    ax.set_title("Patient-level malignant CLDN4 vs T/NK fraction")
    ax.scatter([], [], marker="o", c="k", label="CLDN4+ and T/NK n≥20")
    ax.scatter([], [], marker="s", c="k", label="CLDN4+ and T/NK n<20")
    ax.scatter([], [], marker="x", c="k", label="CLDN4 empty (P1)")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, "fig2_malig_cldn4_vs_frac_tnk")

    # fig3 lineage violin (tumor, drop P1 for CLDN4 biology; show all in caption via detect fig)
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    order = ["Malignant", "Epithelial", "CD4Tconv", "CD8T", "NK", "Mono/Macro", "DC", "B", "Plasma", "Endothelial"]
    data, labels = [], []
    use = tumor_cells.loc[tumor_cells["Patient"] != "P1"] if "CLDN4" in tumor_cells else tumor_cells
    for lin in order:
        v = use.loc[use["_lineage"] == lin, "CLDN4"].to_numpy()
        if len(v) == 0:
            continue
        data.append(v)
        labels.append(f"{lin}\nn={len(v)}")
    ax.violinplot(data, showmeans=True, showextrema=False)
    ax.set_xticks(range(1, len(labels) + 1), labels, fontsize=7)
    ax.set_ylabel("CLDN4 (TISCH log2(TPM/10+1))")
    ax.set_title("Tumor CLDN4 by TISCH lineage (P2–P4; P1 gene-empty)")
    _save(fig, "fig3_cldn4_by_lineage")

    # fig4 composition
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    lin_order = ["Malignant", "Epithelial", "CD4Tconv", "CD8T", "NK", "Mono/Macro", "DC", "B", "Plasma", "Endothelial"]
    pal = plt.cm.tab20(np.linspace(0, 1, len(lin_order)))
    bottom = np.zeros(4)
    ct = pd.crosstab(tumor_cells["Patient"], tumor_cells["_lineage"], normalize="index").reindex(patients).fillna(0)
    for i, lin in enumerate(lin_order):
        vals = ct[lin].to_numpy() if lin in ct.columns else np.zeros(4)
        ax.bar(patients, vals, bottom=bottom, color=pal[i], label=lin, width=0.7)
        bottom = bottom + vals
    ax.set_ylabel("Fraction of tumor cells")
    ax.set_title("Tumor composition (TISCH major-lineage)")
    ax.legend(frameon=False, fontsize=7, bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.subplots_adjust(right=0.72)
    _save(fig, "fig4_tumor_composition")

    # fig5 %pos mal vs tnk
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x = np.arange(4)
    w = 0.35
    ax.bar(x - w / 2, t["CLDN4_mal_pctpos"], w, color="#1f77b4", label="Malignant %pos")
    ax.bar(x + w / 2, t["CLDN4_tnk_pctpos"].fillna(0), w, color="#ff7f0e", label="T/NK %pos")
    ax.set_xticks(x, patients)
    ax.set_ylabel("CLDN4 percent positive")
    ax.set_title("CLDN4 % positive, same-patient compartments")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, "fig5_pctpos_malig_vs_tnk")

    # fig6 tumor vs NAT T/NK
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    nat = samples.loc[samples["source"] == "NAT"].set_index("patient")
    for p in patients:
        ax.plot([0, 1], [nat.loc[p, "frac_tnk"], t.loc[p, "frac_tnk"]], "-o", color=colors[p], label=p)
    ax.set_xticks([0, 1], ["NAT", "Tumor"])
    ax.set_ylabel("T/NK fraction")
    ax.set_title("Same-patient T/NK fraction: NAT vs tumor")
    ax.legend(frameon=False, fontsize=8)
    ax.set_xlim(-0.2, 1.2)
    _save(fig, "fig6_nat_vs_tumor_tnk")

    # fig7 detectability
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    vals = detect.set_index("patient").loc[patients, "CLDN4_n_pos"].to_numpy()
    ax.bar(patients, vals, color=[colors[p] for p in patients])
    ax.set_ylabel("Cells with CLDN4 > 0")
    ax.set_title("CLDN4 detectability in TISCH h5 (all tissues)")
    for i, p in enumerate(patients):
        ax.text(i, vals[i] + max(vals) * 0.02, f"n={int(detect.set_index('patient').loc[p,'n_cells'])}", ha="center", fontsize=8)
    _save(fig, "fig7_cldn4_detectability")

    # fig8 dual bars mean CLDN4 + T/NK n
    fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
    ax2 = ax1.twinx()
    ax1.bar(np.arange(4) - 0.18, t["CLDN4_mal_mean"], 0.36, color="#1f77b4", label="Mal. CLDN4 mean")
    ax2.bar(np.arange(4) + 0.18, t["n_tnk"], 0.36, color="#ff7f0e", label="Tumor T/NK n")
    ax1.set_xticks(np.arange(4), patients)
    ax1.set_ylabel("Malignant CLDN4 mean")
    ax2.set_ylabel("Tumor T/NK cell count")
    ax1.set_title("Malignant CLDN4 mean and T/NK count")
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8)
    _save(fig, "fig8_cldn4_mean_and_tnk_n")


def write_finding(
    tumor: pd.DataFrame,
    samples: pd.DataFrame,
    lineage: pd.DataFrame,
    detect: pd.DataFrame,
    stats_df: pd.DataFrame,
    summary: dict,
    h5_audit: dict,
    prov: dict,
) -> None:
    t = tumor.set_index("patient")
    det = detect.set_index("patient")
    lin_all = lineage.loc[lineage["source"] == "all"].set_index("lineage")
    lin_tum = lineage.loc[lineage["source"] == "Tumor"].set_index("lineage")

    def row(p: str) -> str:
        r = t.loc[p]
        return (
            f"| {p} | {r.histology} | {r.tnm} | {r.sex} | {int(r.n_cells)} | {int(r.n_malignant)} | "
            f"{r.CLDN4_mal_mean:.3f} | {r.CLDN4_mal_pctpos:.1f} | {int(r.n_tnk)} | {r.frac_tnk:.3f} | "
            f"{int(r.n_cd8)}/{int(r.n_cd4)}/{int(r.n_nk)} | {r.CLDN4_tnk_mean:.3f} | {r.CLDN4_tnk_pctpos:.1f} | {r.note} |"
        )

    n_det = int(detect["cldn4_detectable"].sum())
    n_elig = int(tumor["eligible_spearman"].sum())
    n_pair5 = int(tumor["eligible_paired"].sum())
    pair1 = tumor.loc[tumor["cldn4_detectable"] & (tumor["n_tnk"] >= 1)]
    n_pair1 = int(len(pair1))
    n_gt = int(pair1["mal_gt_tnk"].sum())
    p_pair1 = stats_df.loc[stats_df.contrast == "paired_mal_vs_tnk_CLDN4_tnk>=1", "p"]
    p_pair1_s = fmt_p(float(p_pair1.iloc[0])) if len(p_pair1) and pd.notna(p_pair1.iloc[0]) else "n=3, descriptive"

    p1_n = int(det.loc["P1", "n_cells"])
    p1_pos = int(det.loc["P1", "CLDN4_n_pos"])
    h5_bytes = prov["files"]["h5"].get("bytes", 0)
    meta_bytes = prov["files"]["meta"].get("bytes", 0)

    mal_mean = lin_tum.loc["Malignant", "CLDN4_mean"] if "Malignant" in lin_tum.index else float("nan")
    mal_pct = lin_tum.loc["Malignant", "CLDN4_pctpos"] if "Malignant" in lin_tum.index else float("nan")
    tnk_lin = ["CD4Tconv", "CD8T", "NK"]
    tnk_n = int(sum(lin_tum.loc[x, "n_cells"] for x in tnk_lin if x in lin_tum.index))
    tnk_mean = float(
        np.average(
            [lin_tum.loc[x, "CLDN4_mean"] for x in tnk_lin if x in lin_tum.index],
            weights=[lin_tum.loc[x, "n_cells"] for x in tnk_lin if x in lin_tum.index],
        )
    ) if tnk_n else float("nan")

    md = f"""# GSE117570 — malignant CLDN4 vs same-patient T/NK (CLDN4-only)

**Additive only.** TISCH NSCLC catalog already listed this object in `methods/tisch_nsclc_pool/` (sample-level, TACSTD2+CLDN4, Spearman skipped at n=2). This folder does **not** rewrite that pool. It re-scores **CLDN4 only**, with **patient as the unit**, on the public TISCH extract.

Song et al., *Cancer Medicine* 2019 ([PMID 31033233](https://pubmed.ncbi.nlm.nih.gov/31033233/); GEO [GSE117570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117570); TISCH `NSCLC_GSE117570`). Four treatment-naïve early NSCLC cases, tumor + adjacent normal, 10x 3′. No ICI / RECIST / MPR column on GEO or in TISCH CellMetainfo. **No dual-high.** TACSTD2 is an audit gene only (P1 library is alive).

Matrix: TISCH2 `expression.h5` ({h5_bytes/1e6:.1f} MB) + `CellMetainfo_table.tsv` ({meta_bytes/1e6:.2f} MB). GEO series matrix is metadata (0 expression rows). GEO `GSE117570_RAW.tar` UMI tables exist (~19 MB) and were **not** required; TISCH already carries lineage + patient.

## Verdict

CLDN4 is epithelial/malignant-restricted in the three patients where the gene is present. Same-patient T/NK cells are near floor. That is a **compartment** contrast, not an ICI or dual-high claim.

A patient-level Spearman of malignant CLDN4 vs T/NK **fraction** is **empty**. Honest n does not reach 5.

- Catalog patients: **4** (P1–P4).
- Patients with any CLDN4 > 0 in the TISCH h5: **{n_det}** (P2, P3, P4).
- **P1 is gene-empty for CLDN4** (0/{p1_n} cells, tumor+NAT). EPCAM is also 0/{p1_n}. KRT19 / KRT7 / NAPSA / TACSTD2 are present, so this is not a dead library. P1 cannot be used as a CLDN4-low biological point.
- Tumor malignant ≥20: all 4. Tumor T/NK ≥20 **and** CLDN4 detectable: **{n_elig}** (P2 only).
- Spearman (need n≥5): **not computed**.
- Paired malignant > T/NK CLDN4, detectable patients with ≥1 T/NK: **{n_gt}/{n_pair1}** (P2, P3, P4). P3 has 3 T/NK cells; P4 has 8.

## Honest n

| item | n | source |
|---|---:|---|
| GEO / TISCH patients | **4** | P1–P4; 8 samples (tumor+NAT) |
| TISCH cells aligned | **{summary['n_cells_aligned']}** | CellMetainfo ∩ h5 barcodes |
| Tumor cells | {summary['n_tumor_cells']} | TISCH `Source=Tumor` |
| Tumor malignant | {summary['n_malignant_tumor']} | major-lineage `Malignant` |
| Tumor T/NK | {summary['n_tnk_tumor']} | CD4Tconv + CD8T + NK (only T/NK labels present) |
| Patients with CLDN4 > 0 anywhere | **{n_det}** | P2 / P3 / P4 |
| P1 CLDN4-positive cells | **{p1_pos}** | empty row in this extract |
| Eligible for Spearman (mal≥20, T/NK≥20, CLDN4 detectable) | **{n_elig}** | P2 only |
| Eligible paired T/NK≥5 + CLDN4 detectable | **{n_pair5}** | P2 (52 T/NK) and P4 (8 T/NK) |
| Spearman computed | **0** | n&lt;5 |
| ICI / RECIST / MPR labels | **0** | treatment-naïve; none on GEO or TISCH |

Do not write this as n=4 CLDN4-vs-infiltration. Do not recycle P1 as CLDN4-low / T/NK-high: P1 T/NK fraction is high (0.316) but CLDN4 was not measured.

## One-row table

| dataset | unit | n catalog | n CLDN4 detectable | n Spearman-eligible | CLDN4 vs frac T/NK | paired mal&gt;T/NK | dual-high | ICI |
|---|---|---:|---:|---:|---|---|---|---|
| GSE117570 TISCH NSCLC | patient | 4 | **{n_det}** | **{n_elig}** | empty (n&lt;5) | {n_gt}/{n_pair1} (P3 T/NK=3) | no | no |

## Per-patient tumor (unit = patient)

TISCH values are MAESTRO `log2(TPM/10+1)`. T/NK fraction = n_TNK / n_tumor_cells.

| Patient | Histology | TNM | Sex | n tumor | n mal | CLDN4 mal mean | CLDN4 mal %pos | n T/NK | frac T/NK | CD8/CD4/NK | CLDN4 T/NK mean | CLDN4 T/NK %pos | note |
|---|---|---|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|
{row("P1")}
{row("P2")}
{row("P3")}
{row("P4")}

P2 is the only LUSC. P3 and P4 tumors are T-poor after dissociation (3 and 8 T/NK cells). Paper text already notes CD8/NK are low in these tumors versus NAT.

## Paired compartment (malignant CLDN4 vs same-patient T/NK CLDN4)

Primary question as asked: malignant CLDN4 versus T/NK **in the same patient**, not a dual-high gate.

| floor | n | # mal &gt; T/NK | median mal | median T/NK | Wilcoxon p (greater) |
|---|---:|---:|---:|---:|---|
| CLDN4 detectable, T/NK ≥1 | {n_pair1} | {n_gt} | {pair1['CLDN4_mal_mean'].median():.3f} | {pair1['CLDN4_tnk_mean'].median():.3f} | {p_pair1_s} (n=3; descriptive) |
| CLDN4 detectable, T/NK ≥5 | {n_pair5} | {int(t.loc[t['eligible_paired'],'mal_gt_tnk'].sum())} | {t.loc[t['eligible_paired'],'CLDN4_mal_mean'].median():.3f} | {t.loc[t['eligible_paired'],'CLDN4_tnk_mean'].median():.3f} | n=2, not tested |

P2: 1.440 vs 0.041 (52 T/NK). P3: 2.226 vs 0 (3 T/NK). P4: 1.989 vs 0.144 (8 T/NK). P1 both 0 because the gene is missing.

Tumor-lineage restriction (all four patients pooled, so P1 zeros dilute malignant %pos): malignant CLDN4 mean {mal_mean:.3f} ({mal_pct:.1f}% pos, n={int(lin_tum.loc['Malignant','n_cells']) if 'Malignant' in lin_tum.index else 0}); T/NK lineages sit at the floor (CD4Tconv / CD8T / NK each &lt;1% pos in the full object). Extra figure `fig3` drops P1 so the restriction is visible.

## vs T/NK fraction (secondary; empty)

| contrast | n | ρ | p | note |
|---|---:|---|---|---|
| Malignant CLDN4 mean vs frac T/NK, Spearman-eligible | {n_elig} | — | — | n=1 (P2) |
| Same, all CLDN4-detectable tumors | {n_det} | — | — | n=3&lt;5; P3/P4 T/NK-thin |

Descriptive ranks among P2/P3/P4: P2 has the lowest malignant CLDN4 (1.44) and the highest T/NK fraction (0.040); P3/P4 are CLDN4-high and almost T-empty. That is three points. It is not a Spearman.

## Extra figures

| file | what |
|---|---|
| `figures/fig1_paired_malig_vs_tnk.png` | same-patient malignant vs T/NK CLDN4 means |
| `figures/fig2_malig_cldn4_vs_frac_tnk.png` | patient scatter; P1 marked gene-empty |
| `figures/fig3_cldn4_by_lineage.png` | tumor violin, P2–P4 |
| `figures/fig4_tumor_composition.png` | TISCH lineage stack |
| `figures/fig5_pctpos_malig_vs_tnk.png` | % positive, both compartments |
| `figures/fig6_nat_vs_tumor_tnk.png` | paired NAT vs tumor T/NK fraction |
| `figures/fig7_cldn4_detectability.png` | cells with CLDN4&gt;0 (P1 = 0) |
| `figures/fig8_cldn4_mean_and_tnk_n.png` | malignant mean and raw T/NK counts |

## Methods (short)

1. Inputs: TISCH2 `NSCLC_GSE117570_expression.h5` and `CellMetainfo_table.tsv` only. GEO series matrix for diagnosis / TNM / sex / smoking. No FASTQ, no dual-high, no TACSTD2 test.
2. Values are TISCH MAESTRO `log2(TPM/10+1)`. Positive = value &gt; 0.
3. Malignant = TISCH major-lineage `Malignant`. T/NK = `CD4Tconv`, `CD8T`, `NK` (the only T/NK labels in this object).
4. Unit = **patient**. Primary tissue = `Source=Tumor`. NAT is extra (composition / T/NK drop), not a CLDN4 score.
5. Detectable = at least one cell in that patient with CLDN4 &gt; 0. P1 fails.
6. Spearman only if ≥5 patients with ≥20 malignant and ≥20 T/NK and CLDN4 detectable. That filter leaves **P2**.
7. Paired compartment = patient mean CLDN4 in malignant vs T/NK. Wilcoxon is reported only as a label; n=3 is not a powered test.
8. TISCH called 1,066 “Malignant” cells in P2 NAT. Those cells are not used for the tumor score. They are a TISCH over-call risk, not a second tumor.

## How to read this

- **Public processed object exists** (TISCH h5 &lt; 2 GB). This is not an empty accession.
- **Honest n for CLDN4 vs T/NK fraction is 1 eligible / 3 detectable / 4 catalog.** Spearman is empty.
- **Honest paired compartment is 3/3** among patients where CLDN4 exists, with P3/P4 T/NK counts of 3 and 8.
- **P1 cannot support a CLDN4-low story.** The gene is absent from the extract (and so is EPCAM).
- **No dual-high.** No ICI endpoint. No claim that CLDN4 tracks infiltration in this series.

This slice adds a patient-level CLDN4-only table and extra figures. It does not add a testable correlation.

## Files

- `analyze.py` — download TISCH + GEO matrix, score, write this page
- `tables/per_patient_tumor.tsv`, `per_sample.tsv`, `detectability.tsv`, `cldn4_by_lineage.tsv`, `stats.tsv`, `one_row.tsv`, `geo_sample_traits.tsv`, `summary.json`
- `figures/fig1_*.png` … `fig8_*.png` (and PDF)

```bash
python3 -m pip install -r methods/gse117570_cldn4/requirements.txt
python3 methods/gse117570_cldn4/analyze.py
```
"""
    (HERE / "FINDING.md").write_text(md)


if __name__ == "__main__":
    analyze()
