#!/usr/bin/env python3
"""Merged GSE131907+GSE205335 NicheNet-style ligand activity (CLDN4-only).

Sender = CLDN4-high vs low malignant (CLDN4 only; no dual-high).
Receiver = same-patient T/NK. Patient is the unit.
GSE207422 is not re-run (PR #334 n=12 signed tests were NS).

Prior = published NicheNet-v2 ligand–target matrix (Zenodo 7074291),
scored in Python. R / nichenetr is not required.
"""
from __future__ import annotations

import gc
import gzip
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import average_precision_score, roc_auc_score

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
CACHE = Path("/tmp/merge_131907_205335_nichenet_cldn4")
OUT = ROOT / "results"
FIG = ROOT / "figures"
OUT.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "scripts"))
from gene_sets import (  # noqa: E402
    CYTOTOXICITY,
    EXHAUSTION,
    EXTRA_TNK,
    GSE131907_MALIGNANT_SUBTYPES,
    GSE131907_TNK_TYPES,
    GSE131907_TUMOR_ORIGINS,
    IFN,
)
from lib_stats import (  # noqa: E402
    assign_quartiles,
    q4_vs_q1,
    random_effects_dl,
    rank_biserial_pool,
    spearman,
)

DETECT_FRAC = 0.10
MIN_CELLS = 20
EMPIRICAL_P = 0.15
TOP_N = 15
PRIMARY_SETS = ("a_priori_ifn", "a_priori_cytotoxicity")
KEEP_META = {
    "barcode",
    "Index",
    "Sample",
    "sample",
    "orig.ident",
    "patient",
    "patient_id",
    "cohort",
    "origin",
    "histology",
    "recist",
    "compartment",
    "cell_type",
    "cell_subtype",
    "Cell_type",
    "Cell_subtype",
    "Sample_Origin",
    "Sample_Origin_geo",
    "tissue_origin_abbrevation",
    "nCount_RNA",
    "CLDN4",
    "geo_accession",
    "tumor_stage",
    "is_tumor",
    "lineage.sub",
    "lineage.total",
    "tissue",
    "gsm",
    "cancer_subtype",
    "cldn4_high",
    "dual_high_companion",
}
CONCAT_META = [
    "cohort",
    "patient",
    "Sample",
    "origin",
    "histology",
    "recist",
    "compartment",
    "nCount_RNA",
    "CLDN4",
    "cldn4_high",
    "dual_high_companion",
]


def mem_report(tag: str) -> None:
    rss = "?"
    try:
        with open("/proc/self/status") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    rss = line.split()[1]
                    break
    except OSError:
        pass
    print(f"[mem] {tag} VmRSS={rss} kB", flush=True)


def used_log_genes(panel_genes: list[str], lr: pd.DataFrame) -> list[str]:
    """Genes that need log1p(CP10k): CLDN4/TACSTD2, programs, and LR pairs in the panel."""
    want = {"CLDN4", "TACSTD2"}
    want.update(IFN)
    want.update(CYTOTOXICITY)
    want.update(EXHAUSTION)
    want.update(EXTRA_TNK)
    want.update(lr["from"].astype(str))
    want.update(lr["to"].astype(str))
    return [g for g in panel_genes if g in want]


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3g}"


def fmt_num(value: float, digits: int = 3) -> str:
    if not np.isfinite(value):
        return "NA"
    return f"{value:+.{digits}f}"


def log1p_cp10k(umi: np.ndarray, ncount: np.ndarray) -> np.ndarray:
    ncount = np.maximum(ncount.astype(float), 1.0)
    return np.log1p(1e4 * umi.astype(float) / ncount).astype(np.float32)


def ligand_activity(lt: pd.DataFrame, geneset: set[str], background: list[str], ligands: list[str]) -> pd.DataFrame:
    genes = [g for g in background if g in lt.index]
    if len(genes) < 20 or not ligands:
        return pd.DataFrame()
    y = np.array([1 if g in geneset else 0 for g in genes], dtype=int)
    if y.sum() < 2 or y.sum() > len(y) - 2:
        return pd.DataFrame()
    use = [L for L in ligands if L in lt.columns]
    if not use:
        return pd.DataFrame()
    X = lt.loc[genes, use].to_numpy(dtype=np.float32)
    rows = []
    for j, L in enumerate(use):
        s = X[:, j]
        if not np.isfinite(s).all() or np.allclose(s, s[0]):
            continue
        pear = stats.pearsonr(s.astype(np.float64), y)
        try:
            auroc = float(roc_auc_score(y, s))
        except ValueError:
            auroc = np.nan
        try:
            aupr = float(average_precision_score(y, s))
        except ValueError:
            aupr = np.nan
        rows.append(
            {
                "ligand": L,
                "pearson": float(pear.statistic),
                "pearson_p": float(pear.pvalue),
                "auroc": auroc,
                "aupr": aupr,
                "n_background": len(genes),
                "n_geneset": int(y.sum()),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    return out.sort_values("pearson", ascending=False).reset_index(drop=True)


def expressed_genes(sub: pd.DataFrame, genes: list[str], frac: float = DETECT_FRAC) -> list[str]:
    """Detection on raw UMI or log1p(CP10k); both are 0 iff UMI is 0. Skip NaN cells."""
    keep = []
    if len(sub) == 0:
        return keep
    for g in genes:
        if g in sub.columns:
            col = g
        elif f"log_{g}" in sub.columns:
            col = f"log_{g}"
        else:
            continue
        vals = sub[col]
        if getattr(vals, "notna", None) is not None:
            mask = vals.notna()
            if not bool(mask.any()):
                continue
            vals = vals[mask]
        if len(vals) == 0:
            continue
        if float((vals > 0).mean()) >= frac:
            keep.append(g)
    return keep


def mean_score(df: pd.DataFrame, genes: list[str], prefix: str = "log_") -> pd.Series:
    cols = [f"{prefix}{g}" for g in genes if f"{prefix}{g}" in df.columns]
    if not cols:
        return pd.Series(np.nan, index=df.index)
    return df[cols].mean(axis=1)


def parse_gse131907_series(path: Path) -> pd.DataFrame:
    fields: dict[str, list[str]] = {}
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if not line.startswith("!Sample_"):
                continue
            key, _, rest = line.rstrip("\n").partition("\t")
            key = key.removeprefix("!Sample_")
            values = [v.strip('"') for v in rest.split("\t")]
            if key == "characteristics_ch1":
                prefix, _, _ = values[0].partition(": ")
                key = prefix.strip().lower().replace(" ", "_")
                values = [v.split(": ", 1)[-1] if ": " in v else v for v in values]
            if key not in fields:
                fields[key] = values
    n = len(fields["title"])
    return pd.DataFrame({k: v for k, v in fields.items() if len(v) == n})


def parse_gse205335_soft(path: Path) -> pd.DataFrame:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    current["description"] = descriptions[0] if descriptions else ""
                    current["title"] = titles[0] if titles else ""
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
                descriptions, titles = [], []
            elif current is not None and line.startswith("!Sample_title = "):
                titles.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_description = "):
                descriptions.append(line.split(" = ", 1)[1])
            elif current is not None and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, item = value.split(": ", 1)
                    current[key] = item
        if current is not None:
            current["description"] = descriptions[0] if descriptions else ""
            current["title"] = titles[0] if titles else ""
            records.append(current)
    metadata = pd.DataFrame(records)
    if "platform" in metadata.columns:
        read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
        metadata["orig.ident"] = (
            metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
        )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def add_log_cols(cells: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    """log1p(CP10k) as float32 for the requested genes only."""
    ncount = cells["nCount_RNA"].to_numpy()
    logs = {}
    for g in genes:
        if g in cells.columns:
            logs[f"log_{g}"] = log1p_cp10k(cells[g].to_numpy(), ncount)
    if logs:
        cells = pd.concat([cells, pd.DataFrame(logs, index=cells.index)], axis=1)
    return cells


def slim_cells(cells: pd.DataFrame, used_genes: list[str], panel_genes: list[str]) -> pd.DataFrame:
    """Log only used genes, then drop raw UMIs except CLDN4 (needed for %pos)."""
    present = [g for g in used_genes if g in cells.columns]
    ncount = cells["nCount_RNA"].to_numpy()
    logs = {f"log_{g}": log1p_cp10k(cells[g].to_numpy(), ncount) for g in present}
    meta = [c for c in cells.columns if c in KEEP_META]
    out = cells.loc[:, meta].copy()
    if logs:
        out = pd.concat([out, pd.DataFrame(logs, index=out.index)], axis=1)
    dropped_raw = [g for g in panel_genes if g != "CLDN4" and g in cells.columns]
    print(
        f"slim cells={len(out)} cols={out.shape[1]} logs={len(logs)} dropped_raw={len(dropped_raw)}",
        flush=True,
    )
    return out


def concat_cohorts(c131: pd.DataFrame, c205: pd.DataFrame) -> pd.DataFrame:
    """Align metadata + log columns; missing logs stay float32 NaN (no unused raw UMIs)."""
    logs = sorted({c for df in (c131, c205) for c in df.columns if c.startswith("log_")})
    frames = []
    for df in (c131, c205):
        cols = [c for c in CONCAT_META if c in df.columns]
        meta = df.loc[:, cols].copy()
        block = np.empty((len(df), len(logs)), dtype=np.float32)
        for j, col in enumerate(logs):
            if col in df.columns:
                block[:, j] = np.asarray(df[col].to_numpy(), dtype=np.float32)
            else:
                block[:, j] = np.nan
        frames.append(pd.concat([meta, pd.DataFrame(block, index=meta.index, columns=logs)], axis=1))
    return pd.concat(frames, ignore_index=True)


def load_gse131907(panel: pd.DataFrame) -> pd.DataFrame:
    ann_path = CACHE / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    ann = pd.read_csv(ann_path, sep="\t")
    print(f"GSE131907 ann columns={list(ann.columns)} n={len(ann)}", flush=True)
    barcode_col = None
    for cand in ("Index", "barcode", "cell_barcode", "Cell", "cell"):
        if cand in ann.columns:
            barcode_col = cand
            break
    if barcode_col is None:
        if ann.columns[0] == "Unnamed: 0" or ann.index.name:
            ann = ann.reset_index()
        barcode_col = ann.columns[0]
    ann = ann.rename(columns={barcode_col: "barcode"})
    ann["barcode"] = ann["barcode"].astype(str)
    panel = panel.copy()
    panel.index = panel.index.astype(str)
    cells = ann.merge(panel, left_on="barcode", right_index=True, how="inner")
    print(f"GSE131907 joined {len(cells)} / ann {len(ann)} / panel {len(panel)}", flush=True)

    sm = parse_gse131907_series(CACHE / "GSE131907_series_matrix.txt.gz")
    title_col = "title" if "title" in sm.columns else sm.columns[0]
    pid_col = None
    for cand in ("patient_id", "patient", "patientid"):
        if cand in sm.columns:
            pid_col = cand
            break
    if pid_col is None:
        for c in sm.columns:
            if "patient" in c.lower():
                pid_col = c
                break
    sample_map = sm.rename(columns={title_col: "Sample"})
    if pid_col:
        sample_map = sample_map.rename(columns={pid_col: "patient_id"})
    else:
        sample_map["patient_id"] = sample_map["Sample"]
    keep = [c for c in ["Sample", "patient_id", "geo_accession", "tumor_stage"] if c in sample_map.columns]
    sample_map = sample_map[keep].drop_duplicates("Sample")
    sample_map.to_csv(OUT / "gse131907_sample_metadata.tsv", sep="\t", index=False)

    if "Sample" not in cells.columns:
        for cand in ("sample", "orig.ident"):
            if cand in cells.columns:
                cells = cells.rename(columns={cand: "Sample"})
                break
    cells = cells.merge(sample_map, on="Sample", how="left")
    if cells["patient_id"].isna().any():
        cells["patient_id"] = cells["patient_id"].fillna(cells["Sample"])
    origin_col = "Sample_Origin" if "Sample_Origin" in cells.columns else None
    if origin_col is None:
        for cand in ("Sample_Origin_geo", "tissue_origin_abbrevation"):
            if cand in cells.columns:
                origin_col = cand
                break
    if origin_col is None:
        raise RuntimeError(f"no origin column in GSE131907: {list(cells.columns)[:30]}")
    cells["origin"] = cells[origin_col].astype(str)
    cells["cell_type"] = cells["Cell_type"].astype(str) if "Cell_type" in cells.columns else ""
    cells["cell_subtype"] = cells["Cell_subtype"].astype(str) if "Cell_subtype" in cells.columns else ""
    cells["cohort"] = "GSE131907"
    cells["patient"] = cells["patient_id"].astype(str)
    cells["is_tumor"] = cells["origin"].isin(GSE131907_TUMOR_ORIGINS)
    cells["compartment"] = np.select(
        [
            cells["is_tumor"] & cells["cell_subtype"].isin(GSE131907_MALIGNANT_SUBTYPES),
            cells["cell_type"].isin(GSE131907_TNK_TYPES),
        ],
        ["Malignant", "T/NK"],
        default="Other",
    )
    cells["histology"] = "LUAD"
    cells["recist"] = ""
    return cells


def load_gse205335(panel: pd.DataFrame) -> pd.DataFrame:
    ident = pd.read_csv(CACHE / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    print(f"GSE205335 ident columns={list(ident.columns)} n={len(ident)}", flush=True)
    if "barcode" not in ident.columns:
        ident = ident.rename(columns={ident.columns[0]: "barcode"})
    ident["barcode"] = ident["barcode"].astype(str)
    panel = panel.copy()
    panel.index = panel.index.astype(str)
    cells = ident.merge(panel, left_on="barcode", right_index=True, how="inner")
    print(f"GSE205335 joined {len(cells)} / ident {len(ident)} / panel {len(panel)}", flush=True)
    meta = parse_gse205335_soft(CACHE / "GSE205335_family.soft.gz")
    meta.to_csv(OUT / "gse205335_gsm_sample_metadata.csv", index=False)
    merge_cols = [c for c in ["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype"] if c in meta.columns]
    cells = cells.merge(meta[merge_cols], on="orig.ident", how="left")
    if cells["patient"].isna().any():
        raise RuntimeError("GSE205335 patient metadata failed to join")
    cells["cohort"] = "GSE205335"
    cells["patient"] = cells["patient"].astype(str)
    cells["Sample"] = cells["orig.ident"].astype(str)
    cells["origin"] = cells["tissue"].astype(str) if "tissue" in cells.columns else ""
    cells["histology"] = cells["cancer_subtype"].astype(str) if "cancer_subtype" in cells.columns else ""
    cells["recist"] = cells["recist"].astype(str) if "recist" in cells.columns else ""
    cells["compartment"] = np.select(
        [
            cells["lineage.sub"].eq("Malignant cells"),
            cells["lineage.total"].eq("T/NK cells"),
        ],
        ["Malignant", "T/NK"],
        default="Other",
    )
    return cells


def patient_table(cells: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (cohort, pid), g in cells.groupby(["cohort", "patient"], observed=True):
        mg = g[g["compartment"] == "Malignant"]
        tg = g[g["compartment"] == "T/NK"]
        rec = {
            "cohort": cohort,
            "patient": pid,
            "histology": str(g["histology"].iloc[0]) if "histology" in g else "",
            "recist": str(g["recist"].iloc[0]) if "recist" in g else "",
            "n_samples": int(g["Sample"].nunique()) if "Sample" in g else 1,
            "n_cells": int(len(g)),
            "n_malignant": int(len(mg)),
            "n_tnk": int(len(tg)),
            "tnk_frac": float((g["compartment"] == "T/NK").mean()),
        }
        if len(mg):
            rec["mal_CLDN4_mean"] = float(mg["log_CLDN4"].mean())
            rec["mal_CLDN4_pct_pos"] = float((mg["CLDN4"] > 0).mean())
            rec["n_cldn4_high"] = int(mg["cldn4_high"].sum()) if "cldn4_high" in mg else 0
            rec["n_cldn4_low"] = int((~mg["cldn4_high"]).sum()) if "cldn4_high" in mg else 0
            rec["pct_cldn4_high"] = float(mg["cldn4_high"].mean()) if "cldn4_high" in mg else np.nan
            if "log_TACSTD2" in mg.columns:
                rec["mal_TACSTD2_companion"] = float(mg["log_TACSTD2"].mean())
        if len(tg):
            rec["tnk_cyto"] = float(mean_score(tg, CYTOTOXICITY).mean())
            rec["tnk_ifn"] = float(mean_score(tg, IFN).mean())
            rec["tnk_exh"] = float(mean_score(tg, EXHAUSTION).mean())
            rec["tnk_ifn_minus_cyto"] = rec["tnk_ifn"] - rec["tnk_cyto"]
        rows.append(rec)
    return pd.DataFrame(rows)


def plot_scatter(df: pd.DataFrame, x, y, xlabel, ylabel, title, path: Path) -> None:
    colors = {"GSE131907": "#2c7bb6", "GSE205335": "#d7191c"}
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for cohort, sub in df.groupby("cohort"):
        ax.scatter(sub[x], sub[y], c=colors.get(cohort, "0.4"), s=42, label=f"{cohort} n={len(sub)}", zorder=3)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_q4_boxes(patients: pd.DataFrame, metric: str, ylabel: str, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.8), sharey=True)
    for ax, cohort in zip(axes, ["GSE131907", "GSE205335"]):
        sub = patients[patients["cohort"] == cohort].copy()
        if sub.empty or "q_pct" not in sub.columns:
            ax.axis("off")
            continue
        q1 = sub.loc[sub["q_pct"] == "Q1", metric].dropna().to_numpy()
        q4 = sub.loc[sub["q_pct"] == "Q4", metric].dropna().to_numpy()
        bp = ax.boxplot([q1, q4], tick_labels=[f"Q1\nn={len(q1)}", f"Q4\nn={len(q4)}"], patch_artist=True, widths=0.55)
        for patch, color in zip(bp["boxes"], ["#6a8aaa", "#b2182b"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.55)
        rng = np.random.default_rng(0)
        for i, vals in enumerate((q1, q4), start=1):
            ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c="black", s=16, zorder=3)
        ax.set_title(cohort, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].set_ylabel(ylabel)
    fig.suptitle(f"CLDN4 %pos Q4 vs Q1 — {ylabel}", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_ligand_bars(top: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.4))
    for ax, gs, title in [
        (axes[0], "a_priori_ifn", "Prior activity vs IFN set"),
        (axes[1], "a_priori_cytotoxicity", "Prior activity vs cytotoxicity set"),
    ]:
        sub = top[(top["setting"] == "merged_cldn4_high_to_TNK") & (top["geneset"] == gs)].head(12)
        if sub.empty:
            ax.axis("off")
            continue
        ax.barh(sub["ligand"][::-1], sub["pearson"][::-1], color="#4C72B0")
        ax.set_xlabel("Pearson (prior vs gene-set membership)")
        ax.set_title(title, fontsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Merged GSE131907+GSE205335 · CLDN4-high malignant → T/NK", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def plot_signed_forest(tbl: pd.DataFrame, title: str, path: Path) -> None:
    show = tbl.head(12).copy()
    if show.empty:
        return
    fig, ax = plt.subplots(figsize=(7.2, 0.42 * len(show) + 1.6))
    y = np.arange(len(show))
    colors = np.where(show["pooled_rho"] >= 0, "#b2182b", "#2166ac")
    ax.barh(y, show["pooled_rho"], color=colors, alpha=0.85)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.ligand}  p={fmt_p(r.p)}" for r in show.itertuples()], fontsize=8)
    ax.axvline(0, color="0.3", lw=0.8)
    ax.set_xlabel("Pooled Spearman ρ (ligand in CLDN4-high mal. vs T/NK program)")
    ax.set_title(title, fontsize=9)
    ax.invert_yaxis()
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def write_finding(summary: dict, patients: pd.DataFrame, primary: pd.DataFrame, signed: pd.DataFrame) -> None:
    n131 = int((patients["cohort"] == "GSE131907").sum())
    n205 = int((patients["cohort"] == "GSE205335").sum())
    n_samp_131 = int(summary["n"]["gse131907_eligible_samples"])
    n_pat_geo_131 = int(summary["n"]["gse131907_unique_patients"])

    def row_set(gs: str) -> pd.DataFrame:
        if primary.empty:
            return primary
        return primary[(primary["setting"] == "merged_cldn4_high_to_TNK") & (primary["geneset"] == gs)]

    ifn = row_set("a_priori_ifn").head(8)
    cyto = row_set("a_priori_cytotoxicity").head(8)

    def md_table(df: pd.DataFrame) -> list[str]:
        if df.empty:
            return ["_(empty)_", ""]
        lines = ["| Rank | Ligand | Pearson | AUROC | pearson p |", "| --- | --- | ---: | ---: | ---: |"]
        for i, r in enumerate(df.itertuples(), start=1):
            lines.append(
                f"| {i} | **{r.ligand}** | {r.pearson:.3f} | {r.auroc:.3f} | {fmt_p(r.pearson_p)} |"
            )
        lines.append("")
        return lines

    tests = summary.get("patient_tests", {})
    lines = [
        "# FINDING — merged GSE131907+GSE205335 NicheNet-style CLDN4-only ligand activity",
        "",
        "ADDITIVE. **CLDN4-only.** No dual-high TACSTD2×CLDN4. **GSE207422 is not re-run**",
        "(PR #334: n=12 signed tests were NS). Sender = CLDN4-high vs low **malignant**.",
        "Receiver = **same-patient T/NK**. Patient is the unit. GSE131907 GEO extract is",
        "sample-level; cells from eligible samples are collapsed to unique `patient_id`.",
        "",
        "Prior = published **NicheNet-v2** ligand–target matrix (Browaeys et al.;",
        "Zenodo 10.5281/zenodo.7074291, `ligand_target_matrix_nsga2r_final.rds`).",
        "Converted and scored in **Python**. R / `nichenetr` was **not** installed and was **not** run.",
        "",
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| GSE131907 eligible samples | {n_samp_131} | tumor-origin, ≥{MIN_CELLS} author-malig, ≥{MIN_CELLS} T/NK |",
        f"| GSE131907 unique patients | **{n_pat_geo_131}** | **unit** after collapse |",
        f"| GSE205335 patients | **{n205}** | locked ≥{MIN_CELLS} author-malig + T/NK (PR #320) |",
        f"| Merged patients | **{n131 + n205}** | two cohorts; not one mixed bag of cells |",
        f"| Malignant cells (eligible) | {summary['n']['malignant_cells']:,} | author labels |",
        f"| CLDN4-high malignant | {summary['n']['cldn4_high_malignant']:,} | log1p(CP10k) ≥ cohort malignant median |",
        f"| Dual-high companion (not used) | {summary['n']['dual_high_companion_not_used']:,} | TACSTD2 and CLDN4 ≥ median |",
        f"| T/NK cells (eligible) | {summary['n']['tnk_cells']:,} | same-patient receivers |",
        f"| Potential ligands (merged) | **{summary['n']['potential_ligands']}** | ≥10% CLDN4-high + T/NK receptor in v2 LR |",
        f"| Background genes | {summary['n']['background_genes']} | T/NK-expressed ∩ prior targets |",
        f"| IFN / cytotoxicity genes in prior | {summary['n']['ifn_in_prior']} / {summary['n']['cyto_in_prior']} | a priori lists |",
        "",
        "Cells are counts, not replicates. p-values are descriptive.",
        "",
        "## Signed patient-level T/NK programs (not GSE207422)",
        "",
        "CLDN4 score = malignant **%pos** (same primary cut as PR #320). Programs = mean",
        "log1p(CP10k) of a priori IFN / cytotoxicity lists in same-patient T/NK.",
        "",
        "| test | n | result |",
        "| --- | ---: | --- |",
    ]
    for key, label in [
        ("spearman_pct_ifn_merged", "Spearman CLDN4 %pos vs T/NK IFN (DL merge)"),
        ("spearman_pct_cyto_merged", "Spearman CLDN4 %pos vs T/NK cytotoxicity (DL merge)"),
        ("spearman_pct_tnkfrac_merged", "Spearman CLDN4 %pos vs T/NK fraction (DL merge)"),
        ("q4q1_ifn_merged", "Q4 vs Q1 T/NK IFN (pooled r, tails only)"),
        ("q4q1_cyto_merged", "Q4 vs Q1 T/NK cytotoxicity (pooled r, tails only)"),
    ]:
        t = tests.get(key, {})
        if not t:
            continue
        if "pooled_rho" in t:
            lines.append(
                f"| {label} | {t.get('n_patients_total', 'NA')} | "
                f"ρ={fmt_num(t['pooled_rho'])} p={fmt_p(t.get('p', np.nan))} I²={t.get('I2', float('nan')):.0f}% |"
            )
        elif "r_rb" in t:
            lines.append(
                f"| {label} | {t.get('n_compared', 'NA')} | "
                f"r={fmt_num(t['r_rb'])} p={fmt_p(t.get('p', np.nan))} |"
            )
    lines += [
        "",
        "Within-cohort rows are in `results/patient_level_tests.tsv`.",
        "GSE205335 Q4 mixes SCLC with ADC (PR #362); that mix is reported, not hidden.",
        "",
        "## Ligand activity (unsigned NicheNet-v2 prior recovery)",
        "",
        "Primary sender = pooled CLDN4-high malignant from eligible patients in both",
        "cohorts. Receiver background = T/NK-expressed genes present in the prior.",
        "Rank = Pearson of the ligand’s prior target scores vs gene-set membership.",
        "These ranks recover **prior structure** (ISG / MHC ligands sit next to IFN",
        "genes). They are **not** by themselves a signed CLDN4→IFN inductive axis.",
        "",
        f"**A priori IFN ({summary['n']['ifn_in_prior']} genes in prior). Top 8 of {summary['n']['potential_ligands']}:**",
        "",
    ]
    lines += md_table(ifn)
    lines += [
        f"**A priori cytotoxicity ({summary['n']['cyto_in_prior']} genes). Top 8:**",
        "",
    ]
    lines += md_table(cyto)

    if not signed.empty:
        lines += [
            "## Signed ligand vs same-patient T/NK program (patient-level)",
            "",
            "Per-patient mean ligand log1p(CP10k) in **CLDN4-high malignant** cells vs",
            "same-patient T/NK IFN or cytotoxicity. Within-cohort Spearman, then",
            "DerSimonian–Laird on Fisher-z. This is the signed, patient-unit test.",
            "",
        ]
        for vs, label in [("tnk_ifn", "T/NK IFN"), ("tnk_cyto", "T/NK cytotoxicity")]:
            sub = signed[signed["vs"] == vs].sort_values("p").head(8)
            lines += [
                f"**vs {label} (lowest p, unadjusted):**",
                "",
                "| Ligand | k | N | pooled ρ | p | I² |",
                "| --- | ---: | ---: | ---: | ---: | ---: |",
            ]
            for r in sub.itertuples():
                lines.append(
                    f"| {r.ligand} | {int(r.k)} | {int(r.n_patients_total)} | "
                    f"{fmt_num(r.pooled_rho)} | {fmt_p(r.p)} | {r.I2:.0f}% |"
                )
            lines.append("")

    lines += [
        "## What this is not",
        "",
        "- Not GSE207422 NicheNet (PR #286 dual-high; PR #334 CLDN4-only n=12 NS).",
        "- Not a full-transcriptome `nichenetr` R run. Background is T/NK-expressed genes ∩ prior.",
        "- Not TACSTD2∩CLDN4 dual-high senders. Companion dual-high count is reported only.",
        "- Not CellChat / LIANA (those are other PRs on the single cohorts).",
        "- Cells are not n. GSE131907 sample n and unique-patient n are both stated.",
        "",
        "## Files",
        "",
        "- `results/ligand_activity_table.tsv` — primary table (prior Pearson/AUROC + setting)",
        "- `results/ligand_activity_all.tsv` — all settings × gene sets",
        "- `results/ligand_vs_tnk_program_merged.tsv` — signed patient-level ligand tests",
        "- `results/patient_scores.tsv`, `n_table.tsv`, `patient_level_tests.tsv`",
        "- Extra figures: `figures/`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/merge_131907_205335_nichenet_cldn4/scripts/00_download.py",
        "python3 methods/merge_131907_205335_nichenet_cldn4/scripts/01_convert_prior.py",
        "python3 methods/merge_131907_205335_nichenet_cldn4/scripts/02_extract.py",
        "python3 methods/merge_131907_205335_nichenet_cldn4/scripts/03_analyze.py",
        "```",
        "",
    ]
    (ROOT / "FINDING.md").write_text("\n".join(lines) + "\n")
    print("wrote FINDING.md", flush=True)


def _mark_cldn4(cells: pd.DataFrame, name: str) -> pd.DataFrame:
    mal = cells["compartment"] == "Malignant"
    if "log_CLDN4" not in cells.columns:
        raise RuntimeError(f"CLDN4 missing in {name}")
    med = float(cells.loc[mal, "log_CLDN4"].median()) if mal.any() else np.nan
    cells["cldn4_high"] = mal & (cells["log_CLDN4"] >= med)
    if "log_TACSTD2" in cells.columns:
        med_t = float(cells.loc[mal, "log_TACSTD2"].median())
        cells["dual_high_companion"] = mal & (cells["log_CLDN4"] >= med) & (cells["log_TACSTD2"] >= med_t)
    else:
        cells["dual_high_companion"] = False
    print(f"{name} malig={int(mal.sum())} CLDN4-high={int(cells['cldn4_high'].sum())} med={med:.3f}", flush=True)
    return cells


def load_and_slim(cohort: str, panel_path: Path, lr: pd.DataFrame) -> pd.DataFrame:
    mem_report(f"before {cohort} panel")
    panel = pd.read_parquet(panel_path)
    panel_genes = [c for c in panel.columns if c != "nCount_RNA"]
    print(f"panel {cohort}={panel.shape} genes={len(panel_genes)}", flush=True)
    if cohort == "GSE131907":
        cells = load_gse131907(panel)
    elif cohort == "GSE205335":
        cells = load_gse205335(panel)
    else:
        raise ValueError(cohort)
    del panel
    gc.collect()
    mem_report(f"after {cohort} join")
    used = used_log_genes(panel_genes, lr)
    print(f"{cohort} used_log_genes={len(used)} / panel={len(panel_genes)}", flush=True)
    cells = slim_cells(cells, used, panel_genes)
    gc.collect()
    mem_report(f"after {cohort} slim")
    return cells


def main() -> int:
    lr = pd.read_csv(DATA / "lr_network.tsv", sep="\t")
    mem_report("start")

    c131 = load_and_slim("GSE131907", CACHE / "gse131907_panel.parquet", lr)
    c131 = _mark_cldn4(c131, "GSE131907")
    gc.collect()
    c205 = load_and_slim("GSE205335", CACHE / "gse205335_panel.parquet", lr)
    c205 = _mark_cldn4(c205, "GSE205335")
    gc.collect()

    # Eligible patients: ≥20 malignant and ≥20 T/NK after pooling samples
    def eligible_ids(cells: pd.DataFrame) -> set[str]:
        keep = set()
        for pid, g in cells.groupby("patient"):
            n_m = int((g["compartment"] == "Malignant").sum())
            n_t = int((g["compartment"] == "T/NK").sum())
            if n_m >= MIN_CELLS and n_t >= MIN_CELLS:
                keep.add(str(pid))
        return keep

    e131 = eligible_ids(c131)
    e205 = eligible_ids(c205)
    print(f"eligible patients GSE131907={len(e131)} GSE205335={len(e205)}", flush=True)

    # Sample-level n for GSE131907 (honest: GEO is sample-level)
    samp_rows = []
    for s, g in c131.groupby("Sample"):
        samp_rows.append(
            {
                "Sample": s,
                "patient": g["patient"].iloc[0],
                "origin": g["origin"].iloc[0],
                "n_malignant": int((g["compartment"] == "Malignant").sum()),
                "n_tnk": int((g["compartment"] == "T/NK").sum()),
            }
        )
    samp = pd.DataFrame(samp_rows)
    samp["eligible_sample"] = (samp["n_malignant"] >= MIN_CELLS) & (samp["n_tnk"] >= MIN_CELLS)
    samp.to_csv(OUT / "gse131907_sample_n.tsv", sep="\t", index=False)
    n_samp_elig = int(samp["eligible_sample"].sum())
    n_pat_from_samp = int(samp.loc[samp["eligible_sample"], "patient"].nunique())

    c131 = c131[c131["patient"].isin(e131)].copy()
    c205 = c205[c205["patient"].isin(e205)].copy()
    gc.collect()

    def patient_n(df: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for (cohort, pid), g in df.groupby(["cohort", "patient"], observed=True):
            rows.append(
                {
                    "cohort": cohort,
                    "patient": pid,
                    "n_samples": int(g["Sample"].nunique()) if "Sample" in g.columns else 1,
                    "n_cells": int(len(g)),
                    "tnk_frac": float((g["compartment"] == "T/NK").mean()),
                }
            )
        return pd.DataFrame(rows)

    # T/NK fraction uses all cells (including Other). Do not concat unused UMIs.
    keep_n = pd.concat([patient_n(c131), patient_n(c205)], ignore_index=True)
    c131 = c131[c131["compartment"].isin(["Malignant", "T/NK"])].copy()
    c205 = c205[c205["compartment"].isin(["Malignant", "T/NK"])].copy()
    gc.collect()
    cells = concat_cohorts(c131, c205)
    del c131, c205
    gc.collect()
    mem_report("after eligible concat")
    print(f"eligible Malignant+T/NK cells {len(cells)}", flush=True)

    # Recompute CLDN4-high on eligible malignant only (same rule, eligible slice)
    meds = {}
    dual_n = 0
    for cohort, g in cells.groupby("cohort"):
        mal = g["compartment"] == "Malignant"
        med = float(g.loc[mal, "log_CLDN4"].median())
        meds[cohort] = med
        idx = g.index
        cells.loc[idx, "cldn4_high"] = mal.to_numpy() & (g["log_CLDN4"].to_numpy() >= med)
        if "log_TACSTD2" in g.columns:
            med_t = float(g.loc[mal, "log_TACSTD2"].median())
            cells.loc[idx, "dual_high_companion"] = (
                mal.to_numpy()
                & (g["log_CLDN4"].to_numpy() >= med)
                & (g["log_TACSTD2"].to_numpy() >= med_t)
            )
        dual_n += int(cells.loc[idx, "dual_high_companion"].sum())
        print(f"eligible {cohort} med_CLDN4={med:.3f} high={int(cells.loc[idx, 'cldn4_high'].sum())}", flush=True)

    # Refresh patient malignant/TNK program scores on the eligible Malignant+T/NK slice
    # (n_cells / tnk_frac already include Other from the pre-drop table).
    prog = patient_table(cells)
    drop_prog = [c for c in ("n_cells", "n_samples", "tnk_frac") if c in prog.columns]
    patients = prog.drop(columns=drop_prog).merge(keep_n, on=["cohort", "patient"], how="left")
    for cohort, sub_idx in patients.groupby("cohort").groups.items():
        patients.loc[sub_idx, "q_pct"] = assign_quartiles(patients.loc[sub_idx, "mal_CLDN4_pct_pos"])
        patients.loc[sub_idx, "q_mean"] = assign_quartiles(patients.loc[sub_idx, "mal_CLDN4_mean"])
    patients = patients.sort_values(["cohort", "patient"])
    patients.to_csv(OUT / "patient_scores.tsv", sep="\t", index=False)

    # Signed program tests
    tests = []
    cohort_effects = {"ifn": [], "cyto": [], "tnkfrac": [], "q_ifn": [], "q_cyto": []}
    for cohort, sub in patients.groupby("cohort"):
        for x, y, tag in [
            ("mal_CLDN4_pct_pos", "tnk_ifn", "ifn"),
            ("mal_CLDN4_pct_pos", "tnk_cyto", "cyto"),
            ("mal_CLDN4_pct_pos", "tnk_frac", "tnkfrac"),
            ("mal_CLDN4_mean", "tnk_ifn", "ifn_mean"),
            ("mal_CLDN4_mean", "tnk_cyto", "cyto_mean"),
        ]:
            rho, p, n = spearman(sub[x], sub[y])
            tests.append({"cohort": cohort, "test": f"spearman_{x}__{y}", "n": n, "rho": rho, "p": p})
            if tag in cohort_effects:
                cohort_effects[tag].append({"cohort": cohort, "rho": rho, "n": n, "p": p})
        for metric, tag in [("tnk_ifn", "q_ifn"), ("tnk_cyto", "q_cyto"), ("tnk_frac", "q_tnkfrac")]:
            q = q4_vs_q1(sub["mal_CLDN4_pct_pos"], sub[metric])
            if q:
                tests.append({"cohort": cohort, "test": f"q4q1_pct__{metric}", **q})
                if tag in cohort_effects:
                    cohort_effects[tag].append(q)

    merged_tests = {}
    for tag, key in [("ifn", "spearman_pct_ifn_merged"), ("cyto", "spearman_pct_cyto_merged"), ("tnkfrac", "spearman_pct_tnkfrac_merged")]:
        xs = cohort_effects[tag]
        if xs:
            pooled = random_effects_dl([x["rho"] for x in xs], [x["n"] for x in xs])
            merged_tests[key] = pooled
            tests.append({"cohort": "MERGED", "test": key, **pooled})
    for tag, key in [("q_ifn", "q4q1_ifn_merged"), ("q_cyto", "q4q1_cyto_merged")]:
        xs = cohort_effects[tag]
        pooled = rank_biserial_pool(xs)
        if pooled.get("k"):
            merged_tests[key] = pooled
            tests.append({"cohort": "MERGED", "test": key, **pooled})
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(OUT / "patient_level_tests.tsv", sep="\t", index=False)

    mem_report("before prior")
    lt = pd.read_parquet(CACHE / "prior_ligand_target.parquet")
    if lt.shape[0] < lt.shape[1]:
        print(f"warning: lt shape {lt.shape} (expected genes x ligands)", flush=True)
    lt.index = lt.index.astype(str)
    lt.columns = lt.columns.astype(str)
    print(f"prior full {lt.shape}", flush=True)
    mem_report("after prior load")

    ligands_all = sorted(set(lr["from"].astype(str)) & set(lt.columns))
    receptors_all = sorted(set(lr["to"].astype(str)))
    high = cells.loc[cells["cldn4_high"]]
    low = cells.loc[(cells["compartment"] == "Malignant") & (~cells["cldn4_high"])]
    tnk = cells.loc[cells["compartment"] == "T/NK"]

    expr_lig_high = expressed_genes(high, ligands_all)
    expr_rec_tnk = expressed_genes(tnk, receptors_all)
    potential = sorted(set(lr.loc[lr["from"].isin(expr_lig_high) & lr["to"].isin(expr_rec_tnk), "from"]))
    print(
        f"expressed ligands CLDN4-high={len(expr_lig_high)} T/NK receptors={len(expr_rec_tnk)} potential={len(potential)}",
        flush=True,
    )

    pot_by_cohort = {}
    for cohort in ("GSE131907", "GSE205335"):
        h = high[high["cohort"] == cohort]
        t = tnk[tnk["cohort"] == cohort]
        pot_by_cohort[cohort] = sorted(
            set(
                lr.loc[
                    lr["from"].isin(expressed_genes(h, ligands_all))
                    & lr["to"].isin(expressed_genes(t, receptors_all)),
                    "from",
                ]
            )
        )
        print(f"potential {cohort}={len(pot_by_cohort[cohort])}", flush=True)

    log_genes = [c[4:] for c in tnk.columns if c.startswith("log_")]
    tnk_expr = expressed_genes(tnk, [g for g in log_genes if g in lt.index])
    background = sorted(set(tnk_expr) | set(CYTOTOXICITY) | set(IFN) | set(EXHAUSTION) | set(EXTRA_TNK))
    background = [g for g in background if g in lt.index]
    print(f"background genes {len(background)}", flush=True)

    # Subset the 33354 × 1226 prior to background × candidate ligands before scoring.
    cand_ligs = sorted(
        set(potential) | set(pot_by_cohort["GSE131907"]) | set(pot_by_cohort["GSE205335"])
    )
    use_ligs = [L for L in cand_ligs if L in lt.columns]
    use_genes = [g for g in background if g in lt.index]
    lt = lt.loc[use_genes, use_ligs].astype(np.float32)
    gc.collect()
    mem_report("after prior subset")
    print(f"prior subset {lt.shape} (background x candidate ligands)", flush=True)

    # Empirical T/NK gene association with CLDN4 %pos (patient means, DL merge)
    tnk_genes = sorted(
        set(CYTOTOXICITY + IFN + EXHAUSTION + EXTRA_TNK + tnk_expr)
        & {c.replace("log_", "") for c in tnk.columns if c.startswith("log_")}
    )
    gene_pat = []
    for (cohort, pid), g in tnk.groupby(["cohort", "patient"]):
        rec = {"cohort": cohort, "patient": pid}
        for gene in tnk_genes:
            col = f"log_{gene}"
            if col in g.columns:
                rec[gene] = float(g[col].mean())
        gene_pat.append(rec)
    gene_pat_df = pd.DataFrame(gene_pat).merge(
        patients[["cohort", "patient", "mal_CLDN4_pct_pos", "q_pct"]],
        on=["cohort", "patient"],
        how="inner",
    )
    emp_rows = []
    for gene in tnk_genes:
        rhos, ns, ps = [], [], []
        for cohort, sub in gene_pat_df.groupby("cohort"):
            if gene not in sub.columns:
                continue
            rho, p, n = spearman(sub["mal_CLDN4_pct_pos"], sub[gene])
            rhos.append(rho)
            ns.append(n)
            ps.append(p)
        if len(rhos) < 1:
            continue
        pooled = random_effects_dl(rhos, ns) if len(rhos) >= 1 else {}
        emp_rows.append(
            {
                "gene": gene,
                "in_ifn": gene in IFN,
                "in_cyto": gene in CYTOTOXICITY,
                "in_exh": gene in EXHAUSTION,
                "k": pooled.get("k", len(rhos)),
                "n": pooled.get("n_patients_total", int(np.nansum(ns))),
                "pooled_rho": pooled.get("pooled_rho", rhos[0] if rhos else np.nan),
                "p": pooled.get("p", ps[0] if ps else np.nan),
                "I2": pooled.get("I2", 0.0),
            }
        )
    emp = pd.DataFrame(emp_rows).sort_values("p")
    emp.to_csv(OUT / "tnk_gene_vs_cldn4_pct.tsv", sep="\t", index=False)
    emp_up = set(emp.loc[(emp["p"] < EMPIRICAL_P) & (emp["pooled_rho"] > 0), "gene"])
    emp_down = set(emp.loc[(emp["p"] < EMPIRICAL_P) & (emp["pooled_rho"] < 0), "gene"])

    sets = {
        "a_priori_ifn": set(IFN) & set(lt.index),
        "a_priori_cytotoxicity": set(CYTOTOXICITY) & set(lt.index),
        "a_priori_exhaustion_sensitivity": set(EXHAUSTION) & set(lt.index),
        "empirical_CLDN4_TNK_up": emp_up & set(lt.index),
        "empirical_CLDN4_TNK_down": emp_down & set(lt.index),
    }

    activities = []
    settings = [
        ("merged_cldn4_high_to_TNK", potential),
        ("GSE131907_cldn4_high_to_TNK", pot_by_cohort["GSE131907"]),
        ("GSE205335_cldn4_high_to_TNK", pot_by_cohort["GSE205335"]),
    ]
    for set_name, gs in sets.items():
        for setting, ligs in settings:
            act = ligand_activity(lt, gs, background, ligs)
            if act.empty:
                continue
            act.insert(0, "setting", setting)
            act.insert(1, "geneset", set_name)
            activities.append(act)
    act_all = pd.concat(activities, ignore_index=True) if activities else pd.DataFrame()
    if not act_all.empty:
        act_all.to_csv(OUT / "ligand_activity_all.tsv", sep="\t", index=False)

    tops = []
    if not act_all.empty:
        for (setting, gs), sub in act_all.groupby(["setting", "geneset"]):
            chunk = sub.head(TOP_N).copy()
            chunk["rank"] = np.arange(1, len(chunk) + 1)
            tops.append(chunk)
    top = pd.concat(tops, ignore_index=True) if tops else pd.DataFrame()
    if not top.empty:
        top.to_csv(OUT / "top_ligands.tsv", sep="\t", index=False)

    primary = pd.DataFrame()
    if not act_all.empty:
        primary = act_all[
            (act_all["setting"] == "merged_cldn4_high_to_TNK")
            & (act_all["geneset"].isin(PRIMARY_SETS + ("empirical_CLDN4_TNK_up", "empirical_CLDN4_TNK_down")))
        ].copy()
        primary["rank_in_set"] = primary.groupby("geneset")["pearson"].rank(ascending=False, method="first").astype(int)
        primary = primary.sort_values(["geneset", "rank_in_set"])
        primary.to_csv(OUT / "ligand_activity_table.tsv", sep="\t", index=False)
        primary[primary["geneset"].isin(PRIMARY_SETS)].to_csv(OUT / "top_ligands_primary.tsv", sep="\t", index=False)

    # Patient-level ligand means in CLDN4-high cells
    focus = sorted(set(potential) | set(pot_by_cohort["GSE131907"]) | set(pot_by_cohort["GSE205335"]))
    lig_pat = []
    for (cohort, pid), g in high.groupby(["cohort", "patient"]):
        rec = {"cohort": cohort, "patient": pid, "n_cldn4_high": int(len(g))}
        for L in focus:
            col = f"log_{L}"
            if col in g.columns:
                rec[L] = float(g[col].mean())
        lig_pat.append(rec)
    lig_pat_df = pd.DataFrame(lig_pat)
    lig_pat_df.to_csv(OUT / "cldn4_high_ligand_patient_means.tsv", sep="\t", index=False)
    merged_lig = lig_pat_df.merge(
        patients[["cohort", "patient", "tnk_cyto", "tnk_ifn", "tnk_exh", "tnk_frac", "mal_CLDN4_pct_pos"]],
        on=["cohort", "patient"],
    )

    signed_rows = []
    for L in focus:
        if L not in merged_lig.columns:
            continue
        for y in ("tnk_ifn", "tnk_cyto", "tnk_exh", "tnk_frac"):
            rhos, ns, ps = [], [], []
            for cohort, sub in merged_lig.groupby("cohort"):
                rho, p, n = spearman(sub[L], sub[y])
                if n >= 4 and np.isfinite(rho):
                    rhos.append(rho)
                    ns.append(n)
                    ps.append(p)
            if not rhos:
                continue
            pooled = random_effects_dl(rhos, ns)
            signed_rows.append({"ligand": L, "vs": y, **pooled})
    signed = pd.DataFrame(signed_rows)
    if not signed.empty:
        signed = signed.sort_values(["vs", "p"])
        signed.to_csv(OUT / "ligand_vs_tnk_program_merged.tsv", sep="\t", index=False)

    # High vs low ligand delta per patient (sender contrast)
    delta_rows = []
    for (cohort, pid), gmal in cells[cells["compartment"] == "Malignant"].groupby(["cohort", "patient"]):
        h = gmal[gmal["cldn4_high"]]
        l = gmal[~gmal["cldn4_high"]]
        if len(h) < 5 or len(l) < 5:
            continue
        rec = {"cohort": cohort, "patient": pid, "n_high": int(len(h)), "n_low": int(len(l))}
        for L in focus:
            col = f"log_{L}"
            if col in gmal.columns:
                rec[L] = float(h[col].mean() - l[col].mean())
        delta_rows.append(rec)
    delta_df = pd.DataFrame(delta_rows)
    if not delta_df.empty:
        delta_df.to_csv(OUT / "cldn4_high_minus_low_ligand_patient.tsv", sep="\t", index=False)
        hl = []
        for L in focus:
            if L not in delta_df.columns:
                continue
            for cohort, sub in delta_df.groupby("cohort"):
                v = sub[L].dropna()
                if len(v) < 4:
                    continue
                w = stats.wilcoxon(v.to_numpy(), alternative="two-sided", zero_method="wilcox")
                hl.append(
                    {
                        "ligand": L,
                        "cohort": cohort,
                        "n": int(len(v)),
                        "median_delta": float(v.median()),
                        "mean_delta": float(v.mean()),
                        "p": float(w.pvalue),
                    }
                )
        if hl:
            pd.DataFrame(hl).sort_values("p").to_csv(OUT / "cldn4_high_vs_low_ligand_wilcoxon.tsv", sep="\t", index=False)

    n_tbl = pd.DataFrame(
        [
            {"item": "gse131907_eligible_samples", "n": n_samp_elig, "note": f"≥{MIN_CELLS} malig + T/NK; sample-level GEO"},
            {"item": "gse131907_unique_patients", "n": int(len(e131)), "note": "unit after collapse; from eligible samples"},
            {"item": "gse131907_patients_from_eligible_samples", "n": n_pat_from_samp, "note": "should match unique patients"},
            {"item": "gse205335_patients", "n": int(len(e205)), "note": f"≥{MIN_CELLS} author-malig + T/NK"},
            {"item": "merged_patients", "n": int(len(patients)), "note": "GSE131907 patients + GSE205335 patients"},
            {"item": "malignant_cells", "n": int((cells["compartment"] == "Malignant").sum()), "note": "eligible patients only"},
            {"item": "cldn4_high_malignant", "n": int(cells["cldn4_high"].sum()), "note": "CLDN4-only; cohort malignant median"},
            {"item": "dual_high_companion_not_used", "n": int(dual_n), "note": "TACSTD2∩CLDN4; not used as sender"},
            {"item": "tnk_cells", "n": int((cells["compartment"] == "T/NK").sum()), "note": "same-patient receivers"},
            {"item": "potential_ligands", "n": int(len(potential)), "note": "merged CLDN4-high + T/NK receptor"},
            {"item": "potential_ligands_GSE131907", "n": int(len(pot_by_cohort["GSE131907"])), "note": ""},
            {"item": "potential_ligands_GSE205335", "n": int(len(pot_by_cohort["GSE205335"])), "note": ""},
            {"item": "background_genes", "n": int(len(background)), "note": "T/NK-expressed ∩ prior"},
            {"item": "ifn_in_prior", "n": int(len(sets["a_priori_ifn"])), "note": ""},
            {"item": "cyto_in_prior", "n": int(len(sets["a_priori_cytotoxicity"])), "note": ""},
            {"item": "empirical_TNK_up", "n": int(len(emp_up)), "note": f"merged Spearman p<{EMPIRICAL_P}, ρ>0"},
            {"item": "empirical_TNK_down", "n": int(len(emp_down)), "note": f"merged Spearman p<{EMPIRICAL_P}, ρ<0"},
        ]
    )
    n_tbl.to_csv(OUT / "n_table.tsv", sep="\t", index=False)

    summary = {
        "datasets": ["GSE131907", "GSE205335"],
        "skipped": {
            "GSE207422": "PR #334 already scored; n=12 signed tests NS; not re-run",
            "nichenetr_R": "not installed; NicheNet-v2 ligand-target prior scored in Python",
            "dual_high": "TACSTD2 not used to call senders",
        },
        "prior": "NicheNet-v2 ligand_target_matrix_nsga2r_final.rds (Zenodo 7074291)",
        "sender_definition": "CLDN4-only malignant high (log1p CP10k >= cohort malignant median)",
        "unit": "patient (GSE131907 collapsed from eligible samples)",
        "n": {r["item"]: int(r["n"]) for r in n_tbl.to_dict(orient="records")},
        "cutoffs": meds,
        "patient_tests": merged_tests,
        "empirical_sets": {"up": sorted(emp_up), "down": sorted(emp_down), "p_cutoff": EMPIRICAL_P},
        "caveats": [
            "Unsigned prior ranks recover prior structure, not a signed inductive axis",
            "GSE205335 Q4 includes SCLC mixed with ADC",
            "GSE131907 has no ICI/MPR labels",
            "Patient-level ligand tests are unadjusted across ligands",
            "Cells are descriptive counts",
        ],
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    plot_scatter(
        patients,
        "mal_CLDN4_pct_pos",
        "tnk_ifn",
        "Malignant CLDN4 %pos",
        "T/NK IFN mean log1p(CP10k)",
        "Merged · same-patient T/NK IFN vs malignant CLDN4",
        FIG / "fig1_cldn4_vs_tnk_ifn.png",
    )
    plot_scatter(
        patients,
        "mal_CLDN4_pct_pos",
        "tnk_cyto",
        "Malignant CLDN4 %pos",
        "T/NK cytotoxicity mean log1p(CP10k)",
        "Merged · same-patient T/NK cytotoxicity vs malignant CLDN4",
        FIG / "fig2_cldn4_vs_tnk_cyto.png",
    )
    plot_scatter(
        patients,
        "mal_CLDN4_pct_pos",
        "tnk_frac",
        "Malignant CLDN4 %pos",
        "T/NK fraction",
        "Merged · same-patient T/NK fraction vs malignant CLDN4",
        FIG / "fig3_cldn4_vs_tnk_frac.png",
    )
    plot_q4_boxes(patients, "tnk_ifn", "T/NK IFN", FIG / "fig4_q4q1_tnk_ifn.png")
    plot_q4_boxes(patients, "tnk_cyto", "T/NK cytotoxicity", FIG / "fig5_q4q1_tnk_cyto.png")
    if not top.empty:
        plot_ligand_bars(top, FIG / "fig6_ligand_activity_top.png")
    if not signed.empty:
        plot_signed_forest(
            signed[signed["vs"] == "tnk_ifn"].sort_values("p"),
            "Signed: CLDN4-high ligand vs same-patient T/NK IFN",
            FIG / "fig7_signed_ligand_vs_ifn.png",
        )
        plot_signed_forest(
            signed[signed["vs"] == "tnk_cyto"].sort_values("p"),
            "Signed: CLDN4-high ligand vs same-patient T/NK cytotoxicity",
            FIG / "fig8_signed_ligand_vs_cyto.png",
        )

    write_finding(summary, patients, primary, signed if not signed.empty else pd.DataFrame())
    print("wrote results to", OUT, flush=True)
    print(json.dumps({k: summary[k] for k in ("n", "patient_tests")}, indent=2, default=str), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
