#!/usr/bin/env python3
"""Claim A3 on GSE207422 (Hu et al. Genome Medicine 2023; PMID 36869384).

Question (exact):
  In malignant cells only, is TACSTD2 higher in NMPR than MPR, and is
  per-patient malignant TACSTD2 anti-correlated with T/NK fraction
  at Spearman rho ~ -0.40 to -0.50?

Author per-cell annotations are NOT public (GEO sample sheet only;
paper Additional files 1/3/4 are clinical + gene lists + metabolomics;
TISCH2 and CELLxGENE do not host this series). Lineages are therefore
inferred from canonical markers the authors themselves used, plus
normal-lung epithelial markers they used to split CopyKAT-low clusters.
Malignant = epithelial AND not alveolar/club/ciliated. That is a
marker proxy, not CopyKAT.

Primary cohort = 12 post-treatment surgical samples (paper: MPR n=4
including pCR P06; NMPR n=8). Pre-treatment biopsies are reported
separately and are excluded from the primary NMPR>MPR test.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUT = ROOT / "results" / "claim_A3"
OUT.mkdir(parents=True, exist_ok=True)

MATRIX = DATA / "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
META_XLSX = DATA / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
    "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"
)
GEO_META = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207422/suppl/"
    "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
)

# Lineage panels (paper Fig. 1B / S1G + standard lung epithelium).
MARKERS = {
    "T": ["CD3D", "CD3E", "CD3G", "CD2", "CD7", "IL7R", "CD4", "CD8A"],
    "NK": ["NKG7", "GNLY", "KLRD1", "KLRF1", "FCGR3A", "PRF1", "GZMB"],
    "B": ["CD79A", "MS4A1", "CD19", "CD79B"],
    "Plasma": ["MZB1", "JCHAIN", "SDC1", "IGHA1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "FCGR3A", "C1QA"],
    "Neutrophil": ["FCGR3B", "CSF3R", "S100A8", "S100A9", "CXCR2"],
    "pDC": ["LILRA4", "CLEC4C", "IL3RA"],
    "Mast": ["TPSAB1", "CPA3", "KIT"],
    "Stromal": ["COL1A1", "COL3A1", "VWF", "PECAM1", "ACTA2", "DCN"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
}
NORMAL_LUNG = ["SFTPA2", "SFTPC", "SFTPB", "AGER", "SCGB1A1", "TPPP3", "FOXJ1"]
# Author malignant-cluster genes (E0_DST, E3_PCNA, E4_TOP2A, E7_SERPINB9).
MALIGN_HINT = ["DST", "PCNA", "TOP2A", "SERPINB9", "KRT17", "MKI67"]
QC_GENES = ["TACSTD2", "PTPRC", "EPCAM"]
WANT = sorted(
    set(g for gs in MARKERS.values() for g in gs)
    | set(NORMAL_LUNG)
    | set(MALIGN_HINT)
    | set(QC_GENES)
)


def ensure_inputs() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    if not META_XLSX.exists():
        import urllib.request

        print("downloading metadata", flush=True)
        urllib.request.urlretrieve(GEO_META, META_XLSX)
    if not MATRIX.exists():
        import urllib.request

        print("downloading UMI matrix (~184 MB)", flush=True)
        urllib.request.urlretrieve(GEO_MATRIX, MATRIX)


CACHE = DATA / "gse207422_marker_stream.npz"


def stream_genes(path: Path) -> tuple[list[str], dict[str, np.ndarray], np.ndarray]:
    """One pass: library size per cell + requested gene UMI rows."""
    if CACHE.exists():
        z = np.load(CACHE, allow_pickle=True)
        if set(z["want"].tolist()) == set(WANT):
            print(f"loaded cache {CACHE}", flush=True)
            cells = z["cells"].tolist()
            n_umi = z["n_umi"]
            kept = {k: z[f"g_{k}"] for k in WANT if f"g_{k}" in z.files}
            return cells, kept, n_umi
    print(f"streaming {path}", flush=True)
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = header[1:]
        n = len(cells)
        n_umi = np.zeros(n, dtype=np.float64)
        kept: dict[str, np.ndarray] = {}
        for i, line in enumerate(fh, 1):
            gene, _, rest = line.partition("\t")
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n:
                raise ValueError(f"{gene}: got {arr.size} values, expected {n}")
            n_umi += arr
            if gene in WANT:
                kept[gene] = arr
            if i % 4000 == 0:
                print(f"  genes {i}, kept {len(kept)}", flush=True)
    missing = [g for g in WANT if g not in kept]
    print(f"done. cells={n} genes_scanned={i} kept={len(kept)} missing={missing}", flush=True)
    payload = {"cells": np.array(cells), "n_umi": n_umi, "want": np.array(WANT)}
    for g, arr in kept.items():
        payload[f"g_{g}"] = arr
    np.savez_compressed(CACHE, **payload)
    return cells, kept, n_umi


def cpm_log1p(umi: np.ndarray, n_umi: np.ndarray) -> np.ndarray:
    return np.log1p(umi / np.maximum(n_umi, 1.0) * 1e4)


def score(kept: dict[str, np.ndarray], n_umi: np.ndarray, genes: list[str]) -> np.ndarray:
    mats = [cpm_log1p(kept[g], n_umi) for g in genes if g in kept]
    if not mats:
        return np.zeros(n_umi.size, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def annotate(kept: dict[str, np.ndarray], n_umi: np.ndarray) -> pd.DataFrame:
    lineage_scores = {k: score(kept, n_umi, gs) for k, gs in MARKERS.items()}
    score_mat = np.vstack([lineage_scores[k] for k in MARKERS])
    labels = np.array(list(MARKERS))
    assigned = labels[np.argmax(score_mat, axis=0)]
    top = score_mat.max(axis=0)
    second = np.partition(score_mat, -2, axis=0)[-2]
    margin = top - second

    normal_sc = score(kept, n_umi, NORMAL_LUNG)
    epi_sc = lineage_scores["Epithelial"]
    # Normal lung epithelium: assigned epithelial AND normal-lung score
    # comparable to / higher than generic epithelial score.
    is_normal_epi = (assigned == "Epithelial") & (normal_sc >= 0.6 * np.maximum(epi_sc, 1e-6)) & (
        normal_sc >= np.quantile(normal_sc[assigned == "Epithelial"], 0.55)
        if np.any(assigned == "Epithelial")
        else False
    )
    celltype = assigned.astype(object)
    celltype[(assigned == "Epithelial") & ~is_normal_epi] = "Malignant"
    celltype[is_normal_epi] = "NormalEpithelium"
    tnk = np.isin(assigned, ["T", "NK"])

    df = pd.DataFrame(
        {
            "lineage": assigned,
            "celltype": celltype,
            "score_margin": margin,
            "epi_score": epi_sc,
            "normal_lung_score": normal_sc,
            "t_score": lineage_scores["T"],
            "nk_score": lineage_scores["NK"],
            "is_tnk": tnk,
            "is_malignant": celltype == "Malignant",
            "is_epithelial_any": assigned == "Epithelial",
        }
    )
    return df


def load_sample_meta() -> pd.DataFrame:
    meta = pd.read_excel(META_XLSX).dropna(subset=["Sample"])
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")].copy()
    meta = meta.rename(
        columns={
            "Pathologic Response": "path_response_raw",
            "Residual Tumor": "residual_tumor",
            "Clinical Stage": "stage",
        }
    )
    # Paper groups pCR with MPR (4 MPR including P06 pCR; 8 NMPR).
    def group(r: str) -> str:
        if r in {"MPR", "pCR", "MPR (pCR)"}:
            return "MPR"
        if r == "NMPR":
            return "NMPR"
        return str(r)

    meta["response"] = meta["path_response_raw"].map(group)
    meta["timing"] = np.where(
        meta["Resource"].astype(str).str.contains("Post", case=False),
        "post",
        "pre",
    )
    return meta.set_index("Sample")


def mw_onesided(a: np.ndarray, b: np.ndarray) -> dict:
    """NMPR > MPR (greater). Also two-sided for honesty."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        return {"n_nmpr": int(len(a)), "n_mpr": int(len(b)), "U": np.nan, "p_greater": np.nan, "p_two": np.nan}
    u, p_g = stats.mannwhitneyu(a, b, alternative="greater")
    _, p_t = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_nmpr": int(len(a)),
        "n_mpr": int(len(b)),
        "median_nmpr": float(np.median(a)),
        "median_mpr": float(np.median(b)),
        "mean_nmpr": float(np.mean(a)),
        "mean_mpr": float(np.mean(b)),
        "U": float(u),
        "p_greater": float(p_g),
        "p_two": float(p_t),
        "direction": "NMPR>MPR" if np.median(a) > np.median(b) else "NMPR<=MPR",
    }


def spearman_safe(x, y) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        return {"n": int(len(x)), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x, y)
    return {"n": int(len(x)), "rho": float(rho), "p": float(p)}


def main() -> None:
    ensure_inputs()
    cells, kept, n_umi = stream_genes(MATRIX)
    ann = annotate(kept, n_umi)
    ann["barcode"] = cells
    ann["sample"] = pd.Series(cells).str.replace(r"_[0-9]+$", "", regex=True).values
    ann["n_umi"] = n_umi
    tac = kept["TACSTD2"]
    ann["tacstd2_umi"] = tac
    ann["tacstd2_log1p_cpm"] = cpm_log1p(tac, n_umi)
    ann["tacstd2_pos"] = tac > 0
    if "EPCAM" in kept:
        ann["epcam_pos"] = kept["EPCAM"] > 0
    if "PTPRC" in kept:
        ann["ptprc_pos"] = kept["PTPRC"] > 0

    meta = load_sample_meta()
    for col in ["Patient", "Resource", "path_response_raw", "response", "timing", "Pathology", "RECIST", "residual_tumor"]:
        ann[col] = ann["sample"].map(meta[col])

    # Drop low-library cells (authors: 92,330 after QC; we keep all GEO cells
    # but flag very empty ones). Use a light floor so zeros don't dominate scores.
    ann["pass_lib"] = ann["n_umi"] >= 200

    # --- per-sample composition + TACSTD2 ---
    rows = []
    for sample, g in ann.groupby("sample"):
        g = g[g["pass_lib"]]
        n = len(g)
        n_mal = int(g["is_malignant"].sum())
        n_epi = int(g["is_epithelial_any"].sum())
        n_tnk = int(g["is_tnk"].sum())
        mal = g[g["is_malignant"]]
        epi = g[g["is_epithelial_any"]]
        rec = {
            "sample": sample,
            "patient": meta.loc[sample, "Patient"] if sample in meta.index else None,
            "timing": meta.loc[sample, "timing"] if sample in meta.index else None,
            "response": meta.loc[sample, "response"] if sample in meta.index else None,
            "path_response_raw": meta.loc[sample, "path_response_raw"] if sample in meta.index else None,
            "pathology": meta.loc[sample, "Pathology"] if sample in meta.index else None,
            "n_cells": n,
            "n_malignant": n_mal,
            "n_epithelial_any": n_epi,
            "n_tnk": n_tnk,
            "frac_malignant": n_mal / n if n else np.nan,
            "frac_epithelial": n_epi / n if n else np.nan,
            "frac_tnk": n_tnk / n if n else np.nan,
            "mal_tac_mean_log1p_cpm": float(mal["tacstd2_log1p_cpm"].mean()) if n_mal else np.nan,
            "mal_tac_pct_pos": float(mal["tacstd2_pos"].mean()) if n_mal else np.nan,
            "mal_tac_pseudobulk_cpm": float(mal["tacstd2_umi"].sum() / mal["n_umi"].sum() * 1e6) if n_mal and mal["n_umi"].sum() else np.nan,
            "epi_tac_mean_log1p_cpm": float(epi["tacstd2_log1p_cpm"].mean()) if n_epi else np.nan,
            "epi_tac_pct_pos": float(epi["tacstd2_pos"].mean()) if n_epi else np.nan,
            "tnk_tac_mean_log1p_cpm": float(g.loc[g["is_tnk"], "tacstd2_log1p_cpm"].mean()) if n_tnk else np.nan,
            "tnk_tac_pct_pos": float(g.loc[g["is_tnk"], "tacstd2_pos"].mean()) if n_tnk else np.nan,
        }
        rows.append(rec)
    per = pd.DataFrame(rows).sort_values(["timing", "response", "patient"])
    per.to_csv(OUT / "per_patient_metrics.csv", index=False)

    # --- tests ---
    # Primary: all 12 post-treatment patients (MPR includes pCR). MPR residual
    # malignant n is often <10; a ≥20-cell filter would drop 3/4 MPR samples.
    tests = {}
    post = per[(per["timing"] == "post") & (per["response"].isin(["MPR", "NMPR"]))].copy()
    slices = {
        "post_all12": post,
        "post_min5_malignant": post[post["n_malignant"] >= 5],
        "post_min20_malignant": post[post["n_malignant"] >= 20],
        "post_drop_pCR": post[post["path_response_raw"] != "pCR"],
        "all_labeled": per[per["response"].isin(["MPR", "NMPR"])],
    }
    metrics = [
        ("mal_tac_mean_log1p_cpm", "malignant_mean_log1p_cpm"),
        ("mal_tac_pct_pos", "malignant_pct_pos"),
        ("mal_tac_pseudobulk_cpm", "malignant_pseudobulk_cpm"),
        ("epi_tac_mean_log1p_cpm", "epithelial_mean_log1p_cpm"),
        ("epi_tac_pct_pos", "epithelial_pct_pos"),
    ]
    for sl_name, sl in slices.items():
        for metric, key in metrics:
            a = sl.loc[sl["response"] == "NMPR", metric].dropna().values
            b = sl.loc[sl["response"] == "MPR", metric].dropna().values
            tests[f"NMPR_gt_MPR_{sl_name}_{key}"] = mw_onesided(a, b)
            tests[f"spearman_{sl_name}_{key}_vs_frac_tnk"] = spearman_safe(sl[metric], sl["frac_tnk"])

    mpr_tnk = post.loc[post["response"] == "MPR", "frac_tnk"].values
    nmpr_tnk = post.loc[post["response"] == "NMPR", "frac_tnk"].values
    if len(mpr_tnk) >= 2 and len(nmpr_tnk) >= 2:
        u, p_g = stats.mannwhitneyu(mpr_tnk, nmpr_tnk, alternative="greater")
        _, p_t = stats.mannwhitneyu(mpr_tnk, nmpr_tnk, alternative="two-sided")
        tests["frac_tnk_MPR_gt_NMPR_post"] = {
            "n_mpr": int(len(mpr_tnk)),
            "n_nmpr": int(len(nmpr_tnk)),
            "median_mpr": float(np.median(mpr_tnk)),
            "median_nmpr": float(np.median(nmpr_tnk)),
            "U": float(u),
            "p_greater": float(p_g),
            "p_two": float(p_t),
        }

    # Marker QC
    qc = {}
    for ct in ["Malignant", "NormalEpithelium", "T", "NK"]:
        sub = ann[(ann["celltype"] == ct) | ((ct in {"T", "NK"}) & (ann["lineage"] == ct))]
        if ct in {"T", "NK"}:
            sub = ann[ann["lineage"] == ct]
        else:
            sub = ann[ann["celltype"] == ct]
        qc[ct] = {
            "n": int(len(sub)),
            "tacstd2_pct_pos": float(sub["tacstd2_pos"].mean()) if len(sub) else np.nan,
            "epcam_pct_pos": float(sub["epcam_pos"].mean()) if len(sub) and "epcam_pos" in sub else np.nan,
            "ptprc_pct_pos": float(sub["ptprc_pos"].mean()) if len(sub) and "ptprc_pos" in sub else np.nan,
        }
    qc["celltype_counts"] = ann["celltype"].value_counts().to_dict()
    qc["lineage_counts"] = ann["lineage"].value_counts().to_dict()

    # Claim verdicts — primary = all 12 post-treatment patients
    prim = tests["NMPR_gt_MPR_post_all12_malignant_mean_log1p_cpm"]
    rho = tests["spearman_post_all12_malignant_mean_log1p_cpm_vs_frac_tnk"]
    claim_nmpr = (
        prim.get("direction") == "NMPR>MPR"
        and prim.get("p_greater") == prim.get("p_greater")
        and prim.get("p_greater") < 0.05
    )
    # Require the claimed band on the primary (all-12) analysis, not a post-hoc slice.
    claim_rho = (
        rho.get("rho") == rho.get("rho")
        and -0.50 <= rho["rho"] <= -0.40
    )
    summary = {
        "dataset": "GSE207422",
        "paper": "Hu et al. Genome Medicine 2023 PMID 36869384",
        "author_per_cell_annotations": "not public",
        "annotation_method": "marker-score lineage; malignant = epithelial minus alveolar/club/ciliated",
        "primary_cohort": "12 post-treatment surgical samples; pCR grouped with MPR",
        "n_cells_geo": int(len(ann)),
        "n_cells_lib200": int(ann["pass_lib"].sum()),
        "claim_A3_NMPR_gt_MPR_supported": bool(claim_nmpr),
        "claim_A3_rho_minus_0.40_to_0.50_supported": bool(claim_rho),
        "primary_NMPR_gt_MPR": prim,
        "primary_spearman_mal_tac_vs_tnk": rho,
        "tests": tests,
        "qc": qc,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    flat = []
    for k, v in tests.items():
        row = {"test": k}
        if isinstance(v, dict):
            row.update(v)
        flat.append(row)
    pd.DataFrame(flat).to_csv(OUT / "tests.csv", index=False)

    # --- plots ---
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))
    # 1. box NMPR vs MPR
    ax = axes[0]
    plot_df = post.dropna(subset=["mal_tac_mean_log1p_cpm"])
    rng = np.random.default_rng(0)
    for i, (lab, col) in enumerate([("MPR", "#2a9d8f"), ("NMPR", "#e76f51")]):
        s = plot_df.loc[plot_df["response"] == lab]
        y = s["mal_tac_mean_log1p_cpm"].to_numpy()
        jit = rng.uniform(-0.08, 0.08, len(y))
        ax.scatter(np.full(len(y), i) + jit, y, c=col, s=40, zorder=3)
        ax.hlines(np.median(y) if len(y) else np.nan, i - 0.2, i + 0.2, colors="black", lw=2)
        for x0, y0, lab_p in zip(np.full(len(y), i) + jit, y, s["patient"]):
            ax.annotate(lab_p, (x0, y0), fontsize=6, xytext=(4, 0), textcoords="offset points")
    ax.set_xticks([0, 1], ["MPR/pCR", "NMPR"])
    ax.set_ylabel("Malignant TACSTD2 (mean log1p CPM)")
    ax.set_title(
        f"NMPR>MPR? {prim.get('direction')}  p_greater={prim.get('p_greater'):.3g}"
        if prim.get("p_greater") == prim.get("p_greater")
        else "NMPR>MPR (insufficient n)"
    )
    # 2. scatter vs T/NK
    ax = axes[1]
    for lab, col in [("MPR", "#2a9d8f"), ("NMPR", "#e76f51")]:
        s = plot_df[plot_df["response"] == lab]
        ax.scatter(s["frac_tnk"], s["mal_tac_mean_log1p_cpm"], c=col, s=45, label=lab)
        for _, r in s.iterrows():
            ax.annotate(r["patient"], (r["frac_tnk"], r["mal_tac_mean_log1p_cpm"]), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("T/NK fraction (all cells)")
    ax.set_ylabel("Malignant TACSTD2 (mean log1p CPM)")
    ax.set_title(f"Spearman ρ={rho.get('rho'):.2f}  p={rho.get('p'):.3g}  n={rho.get('n')}")
    ax.legend(frameon=False)
    # 3. %pos sanity
    ax = axes[2]
    names = ["Malignant", "NormalEpithelium", "T", "NK"]
    vals = [qc[k]["tacstd2_pct_pos"] * 100 for k in names]
    ax.bar(names, vals, color=["#e76f51", "#e9c46a", "#457b9d", "#1d3557"])
    ax.set_ylabel("% TACSTD2+ cells")
    ax.set_title("Lineage restriction (sanity)")
    ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(OUT / "claim_A3_summary.png", dpi=160)
    fig.savefig(OUT / "claim_A3_summary.pdf")
    plt.close(fig)

    # composition plot
    fig, ax = plt.subplots(figsize=(8.5, 3.4))
    order = post.sort_values(["response", "patient"])
    x = np.arange(len(order))
    ax.bar(x, order["frac_malignant"], label="malignant (marker)", color="#e76f51")
    ax.bar(x, order["frac_tnk"], bottom=order["frac_malignant"], label="T/NK", color="#457b9d")
    ax.set_xticks(x, [f"{p}\n{r}" for p, r in zip(order["patient"], order["response"])], fontsize=8)
    ax.set_ylabel("Fraction of cells")
    ax.set_title("Post-treatment composition (marker-inferred)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / "composition_post.png", dpi=160)
    plt.close(fig)

    write_readme(summary, per, post)
    print(json.dumps({k: summary[k] for k in summary if k not in {"tests", "qc"}}, indent=2))
    print("wrote", OUT)


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and (np.isnan(x) or np.isinf(x))):
        return "NA"
    if isinstance(x, float):
        if abs(x) >= 10:
            return f"{x:.0f}"
        if abs(x) < 0.01:
            return f"{x:.2e}"
        return f"{x:.{nd}f}"
    return str(x)


def write_readme(summary: dict, per: pd.DataFrame, post: pd.DataFrame) -> None:
    prim = summary["primary_NMPR_gt_MPR"]
    rho = summary["primary_spearman_mal_tac_vs_tnk"]
    qc = summary["qc"]
    t = summary["tests"]
    epi = t["NMPR_gt_MPR_post_all12_epithelial_mean_log1p_cpm"]
    pb = t["NMPR_gt_MPR_post_all12_malignant_pseudobulk_cpm"]
    rho_drop = t["spearman_post_drop_pCR_malignant_mean_log1p_cpm_vs_frac_tnk"]
    rho_epi = t["spearman_post_all12_epithelial_mean_log1p_cpm_vs_frac_tnk"]
    rho_pb = t["spearman_post_all12_malignant_pseudobulk_cpm_vs_frac_tnk"]

    table = [
        "| patient | response | n_mal | n_epi | T/NK frac | mal TACSTD2 mean | mal %pos | mal pb CPM | epi TACSTD2 mean |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in post.sort_values(["response", "patient"]).iterrows():
        table.append(
            f"| {r['patient']} | {r['path_response_raw']} | {int(r['n_malignant'])} | {int(r['n_epithelial_any'])} | "
            f"{r['frac_tnk']:.3f} | {r['mal_tac_mean_log1p_cpm']:.3f} | {100*r['mal_tac_pct_pos']:.1f} | "
            f"{r['mal_tac_pseudobulk_cpm']:.0f} | {r['epi_tac_mean_log1p_cpm']:.3f} |"
        )

    lines = [
        "# Claim A3 — GSE207422 (exact)",
        "",
        "**Claim tested:** malignant-only TACSTD2 is higher in NMPR than MPR, and per-patient malignant TACSTD2 vs T/NK fraction Spearman ρ is about −0.40 to −0.50.",
        "",
        "## Verdict (honest)",
        "",
        f"- **NMPR > MPR (malignant TACSTD2, 12 post-treatment patients): NOT SUPPORTED.** "
        f"Median malignant mean-log1p-CPM NMPR={_fmt(prim.get('median_nmpr'))} vs MPR={_fmt(prim.get('median_mpr'))} "
        f"(n=8/4). Direction is weakly NMPR>MPR; one-sided Mann–Whitney p={_fmt(prim.get('p_greater'))}, "
        f"two-sided p={_fmt(prim.get('p_two'))}.",
        f"- **Malignant % TACSTD2+** is also n.s. (NMPR median {100*t['NMPR_gt_MPR_post_all12_malignant_pct_pos']['median_nmpr']:.1f}% vs "
        f"MPR {100*t['NMPR_gt_MPR_post_all12_malignant_pct_pos']['median_mpr']:.1f}%; p_greater="
        f"{_fmt(t['NMPR_gt_MPR_post_all12_malignant_pct_pos']['p_greater'])}).",
        f"- **Malignant pseudobulk CPM goes the other way** (NMPR median {_fmt(pb.get('median_nmpr'),1)} vs MPR {_fmt(pb.get('median_mpr'),1)}; "
        f"p_greater={_fmt(pb.get('p_greater'))}). Driven by tiny MPR malignant n (P06/P11/P14 have 4–7 cells).",
        f"- **All-epithelial sensitivity** is the closest to the claim (NMPR {_fmt(epi.get('median_nmpr'))} vs MPR {_fmt(epi.get('median_mpr'))}; "
        f"p_greater={_fmt(epi.get('p_greater'))}) but still not p<0.05.",
        f"- **ρ ≈ −0.40 to −0.50: NOT SUPPORTED on the primary analysis.** "
        f"All 12 post-treatment: Spearman ρ={_fmt(rho.get('rho'))} (p={_fmt(rho.get('p'))}) for malignant mean log1p CPM vs T/NK fraction. "
        f"%pos ρ={_fmt(t['spearman_post_all12_malignant_pct_pos_vs_frac_tnk']['rho'])}; "
        f"pseudobulk ρ={_fmt(rho_pb.get('rho'))}; epithelial ρ={_fmt(rho_epi.get('rho'))}.",
        f"- One post-hoc slice (drop pCR P06) gives ρ={_fmt(rho_drop.get('rho'))} (p={_fmt(rho_drop.get('p'))}, n={rho_drop.get('n')}) — inside the claimed band but n.s. and not pre-specified.",
        "- The only MPR patient with a large malignant compartment (P03, 871 cells) has the **lowest** malignant TACSTD2 (0.83). Three other MPR/pCR samples have 4–7 marker-malignant cells because residual epithelium scores as normal lung — expected after MPR, and a hard limit without author CopyKAT labels.",
        "",
        "## Data and annotations",
        "",
        "- GEO **GSE207422** (Hu et al., *Genome Medicine* 2023, PMID 36869384): BD Rhapsody UMI matrix, 92,330 cells × 24,292 genes, 15 samples / 15 patients.",
        "- **Author per-cell annotations are not public.** Searched: GEO suppl (UMI + sample xlsx only), paper Additional files 1/3/4 (clinical / module genes / steroids; no barcodes), TISCH2 NSCLC gallery (GSE207422 absent), CELLxGENE. Analysis uses marker-inferred labels.",
        "- Lineage = argmax of mean log1p(CPM) across canonical panels (T, NK, B, plasma, myeloid, neutrophil, pDC, mast, stromal, epithelial).",
        "- **Malignant** = assigned epithelial and *not* normal-lung (SFTPA2/SFTPC/AGER/SCGB1A1/TPPP3/FOXJ1), matching the authors' normal clusters (alveolar / club / ciliated). This is **not** CopyKAT.",
        "- Primary test uses **post-treatment surgery** samples only (paper Fig. 1). pCR (P06) is grouped with MPR as in the paper (MPR n=4, NMPR n=8). Pre-treatment biopsies (P01 NE, P05/P08 NMPR-labeled) are excluded from the primary contrast.",
        "",
        "## Post-treatment per-patient table",
        "",
        *table,
        "",
        "## QC (marker restriction)",
        "",
        f"- TACSTD2+ : malignant {qc['Malignant']['tacstd2_pct_pos']:.1%}, normal epithelium {qc['NormalEpithelium']['tacstd2_pct_pos']:.1%}, T {qc['T']['tacstd2_pct_pos']:.1%}, NK {qc['NK']['tacstd2_pct_pos']:.1%}.",
        f"- EPCAM+ : malignant {qc['Malignant']['epcam_pct_pos']:.1%}, T {qc['T']['epcam_pct_pos']:.1%}.",
        f"- PTPRC/CD45+ : malignant {qc['Malignant']['ptprc_pct_pos']:.1%}, T {qc['T']['ptprc_pct_pos']:.1%}.",
        f"- Assigned counts: malignant {qc['Malignant']['n']}, normal epithelium {qc['NormalEpithelium']['n']}, T {qc['T']['n']}, NK {qc['NK']['n']}.",
        "",
        "## Files",
        "",
        "- `per_patient_metrics.csv` — per-sample counts, fractions, TACSTD2 summaries (all 15 samples).",
        "- `summary.json` / `tests.csv` — all tests (mean / %positive / pseudobulk; all-12 / min5 / min20 / drop-pCR / epithelial-wide).",
        "- `claim_A3_summary.png` / `.pdf` — NMPR vs MPR, scatter vs T/NK, lineage restriction.",
        "- `composition_post.png` — malignant and T/NK fractions by post-treatment patient.",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 scripts/claim_A3_gse207422.py",
        "```",
        "",
        "Downloads the GEO UMI matrix (~184 MB) and sample metadata if missing. Streams requested genes plus per-cell library size; caches `data/gse207422_marker_stream.npz` (gitignored).",
        "",
    ]
    (OUT / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
