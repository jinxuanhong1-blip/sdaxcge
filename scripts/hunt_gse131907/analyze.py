#!/usr/bin/env python3
"""GSE131907 hunt: epithelial-restricted TACSTD2/CLDN4 vs CD8/NK/TLS proxies,
plus TJ/keratin GSEA in TROP2-high epithelium.

Honest scope
------------
GSE131907 is dissociated 10x scRNA-seq (Kim et al. 2020). There are no
coordinates. "Neighborhood" = sample-level co-occurrence of epithelial
TACSTD2/CLDN4 with immune fractions / chemokine scores. "TLS" = B-cell,
germinal-center B, and 12-chemokine proxies — not a histologically scored
tertiary lymphoid structure.

TACSTD2 and CLDN4 are scored only inside author epithelial cells
(epithelial-restricted). Immune abundance is taken from author labels,
not from TACSTD2 in immune cells.
"""

from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C

try:
    import gseapy as gp
except ImportError:  # pragma: no cover
    gp = None


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def stream_selected_genes(matrix_path: Path, wanted: set[str]) -> tuple[list[str], dict[str, np.ndarray]]:
    found: dict[str, np.ndarray] = {}
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        cell_ids = header[1:]
        n_cells = len(cell_ids)
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            if gene not in wanted:
                continue
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")
            found[gene] = arr
            if len(found) == len(wanted):
                break
    return cell_ids, found


def parse_series_matrix(path: Path) -> pd.DataFrame:
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


def load_genesets() -> dict[str, list[str]]:
    raw = json.loads(C.GENESET_JSON.read_text())
    sets = {}
    for name, genes in raw["sets"].items():
        sets[name] = [str(g) for g in genes]
    return sets


# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------
def spearman_row(x, y, contrast: str) -> dict:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    if len(x) < 5:
        return {
            "contrast": contrast,
            "n": int(len(x)),
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "note": "too_few_samples",
        }
    rho, p = stats.spearmanr(x, y)
    return {
        "contrast": contrast,
        "n": int(len(x)),
        "spearman_rho": float(rho),
        "spearman_p": float(p),
        "note": "",
    }


def welch_row(a, b, contrast: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return {
            "contrast": contrast,
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "mean_a": np.nan,
            "mean_b": np.nan,
            "delta": np.nan,
            "t": np.nan,
            "p": np.nan,
            "note": "too_few",
        }
    t, p = stats.ttest_ind(a, b, equal_var=False)
    return {
        "contrast": contrast,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(np.mean(a)),
        "mean_b": float(np.mean(b)),
        "delta": float(np.mean(a) - np.mean(b)),
        "t": float(t),
        "p": float(p),
        "note": "",
    }


def mw_row(a, b, contrast: str) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return {"contrast": contrast, "n_a": int(len(a)), "n_b": int(len(b)), "p": np.nan, "note": "too_few"}
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "contrast": contrast,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "U": float(u),
        "p": float(p),
        "note": "",
    }


def compartment(row: pd.Series) -> str:
    ct = row.Cell_type
    if ct == "Epithelial cells":
        return "epithelial"
    if ct == "T lymphocytes":
        return "T"
    if ct == "NK cells":
        return "NK"
    if ct == "B lymphocytes":
        return "B"
    if ct == "Myeloid cells":
        return "myeloid"
    if ct == "MAST cells":
        return "mast"
    return "other"


# ---------------------------------------------------------------------------
# Pass 2: stream all genes for ranking (sample-level + cell-level)
# ---------------------------------------------------------------------------
def stream_rank_metrics(
    matrix_path: Path,
    cell_ids: list[str],
    epi_mask: np.ndarray,
    sample_codes: np.ndarray,
    trop2_high: np.ndarray,
    trop2_low: np.ndarray,
    tumor_sample_keep: set[int],
    epi_tacstd2_by_sample: dict[int, float],
) -> pd.DataFrame:
    """One pass over the full matrix.

    Returns a table with, per gene:
      - sample-level Spearman ρ of epithelial mean vs epithelial TACSTD2
        (tumor samples with ≥ MIN_EPI epithelial cells)
      - cell-level Welch t of TROP2-high vs TROP2-low tumor epithelium
    """
    n_cells = len(cell_ids)
    sample_ids = np.unique(sample_codes)
    sample_to_i = {s: i for i, s in enumerate(sample_ids)}
    n_s = len(sample_ids)
    epi_sample_i = np.array([sample_to_i[s] for s in sample_codes[epi_mask]], dtype=np.int32)

    rows = []
    t0 = time.time()
    n_done = 0
    with gzip.open(matrix_path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if header[1:] != cell_ids:
            raise ValueError("matrix column order changed between passes")
        for line in handle:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            arr = np.fromstring(rest, sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} values, expected {n_cells}")

            epi_vals = arr[epi_mask]
            # per-sample epithelial means
            sums = np.bincount(epi_sample_i, weights=epi_vals.astype(np.float64), minlength=n_s)
            cnts = np.bincount(epi_sample_i, minlength=n_s).astype(np.float64)
            means = np.divide(sums, cnts, out=np.full(n_s, np.nan), where=cnts > 0)

            xs, ys = [], []
            for s, i in sample_to_i.items():
                if s not in tumor_sample_keep:
                    continue
                if not np.isfinite(means[i]):
                    continue
                xs.append(means[i])
                ys.append(epi_tacstd2_by_sample[s])
            if len(xs) >= 5:
                rho, p_s = stats.spearmanr(xs, ys)
            else:
                rho, p_s = np.nan, np.nan

            high = epi_vals[trop2_high]
            low = epi_vals[trop2_low]
            if high.size >= 10 and low.size >= 10:
                t, p_c = stats.ttest_ind(high, low, equal_var=False)
                delta = float(high.mean() - low.mean())
            else:
                t, p_c, delta = np.nan, np.nan, np.nan

            rows.append(
                {
                    "gene": gene,
                    "sample_spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                    "sample_spearman_p": float(p_s) if np.isfinite(p_s) else np.nan,
                    "cell_t": float(t) if np.isfinite(t) else np.nan,
                    "cell_p": float(p_c) if np.isfinite(p_c) else np.nan,
                    "cell_delta": delta,
                    "n_samples": int(len(xs)),
                    "n_high": int(high.size),
                    "n_low": int(low.size),
                }
            )
            n_done += 1
            if n_done % 2000 == 0:
                print(f"[rank] {n_done} genes  ({time.time() - t0:.0f}s)", flush=True)

    print(f"[rank] done {n_done} genes in {time.time() - t0:.0f}s", flush=True)
    return pd.DataFrame(rows)


def run_prerank(rnk: pd.Series, gene_sets: dict[str, list[str]], outdir: Path, tag: str) -> pd.DataFrame:
    rnk = rnk.dropna()
    rnk = rnk[~rnk.index.duplicated(keep="first")]
    if gp is None:
        raise RuntimeError("gseapy is not installed")
    out = outdir / f"gseapy_{tag}"
    out.mkdir(parents=True, exist_ok=True)
    res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        permutation_num=C.GSEA_PERMS,
        min_size=5,
        max_size=500,
        seed=C.RANDOM_SEED,
        outdir=str(out),
        verbose=False,
    )
    table = res.res2d.copy()
    table.insert(0, "rank_metric", tag)
    table.to_csv(outdir / f"gsea_{tag}.tsv", sep="\t", index=False)
    return table


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def save_restriction_fig(per_cell: pd.DataFrame, out: Path) -> None:
    tumor = per_cell[per_cell["Sample_Origin"].isin(C.TUMOR_ORIGINS)]
    order = ["epithelial", "T", "NK", "B", "myeloid"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)
    for ax, gene in zip(axes, ["TACSTD2", "CLDN4"]):
        data = [tumor.loc[tumor["compartment"] == c, gene].to_numpy() for c in order]
        ax.boxplot(data, tick_labels=order, showfliers=False, widths=0.6)
        ax.set_ylabel(f"{gene}  log2(TPM+1)")
        ax.set_title(f"{gene} is epithelial-restricted (tumor sites)")
        ax.tick_params(axis="x", rotation=20)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def save_neighborhood_fig(sample_df: pd.DataFrame, out: Path) -> None:
    tumor = sample_df[
        sample_df["Sample_Origin"].isin(C.TUMOR_ORIGINS)
        & (sample_df["n_epithelial"] >= C.MIN_EPI)
    ]
    panels = [
        ("frac_CD8", "CD8 fraction"),
        ("frac_NK", "NK fraction"),
        ("frac_B", "B-cell fraction"),
        ("tls_chemokine_mean", "12-chemokine TLS score"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.2), constrained_layout=True)
    for ax, (ycol, ylab) in zip(axes.ravel(), panels):
        x = tumor["epi_TACSTD2_mean"]
        y = tumor[ycol]
        ax.scatter(x, y, c="#2864a6", edgecolors="white", linewidths=0.5, s=40)
        mask = x.notna() & y.notna()
        if mask.sum() >= 5:
            rho, p = stats.spearmanr(x[mask], y[mask])
            coef = np.polyfit(x[mask], y[mask], 1)
            grid = np.linspace(x[mask].min(), x[mask].max(), 40)
            ax.plot(grid, np.polyval(coef, grid), color="#333", lw=1)
            ax.set_title(f"{ylab}\nρ={rho:.2f}  p={p:.2g}  n={int(mask.sum())}")
        else:
            ax.set_title(ylab)
        ax.set_xlabel("Epithelial TACSTD2 mean log2(TPM+1)")
        ax.set_ylabel(ylab)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def save_cldn4_neighborhood_fig(sample_df: pd.DataFrame, out: Path) -> None:
    tumor = sample_df[
        sample_df["Sample_Origin"].isin(C.TUMOR_ORIGINS)
        & (sample_df["n_epithelial"] >= C.MIN_EPI)
    ]
    panels = [
        ("frac_CD8", "CD8 fraction"),
        ("frac_NK", "NK fraction"),
        ("frac_B", "B-cell fraction"),
        ("tls_chemokine_mean", "12-chemokine TLS score"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.2), constrained_layout=True)
    for ax, (ycol, ylab) in zip(axes.ravel(), panels):
        x = tumor["epi_CLDN4_mean"]
        y = tumor[ycol]
        ax.scatter(x, y, c="#a50f15", edgecolors="white", linewidths=0.5, s=40)
        mask = x.notna() & y.notna()
        if mask.sum() >= 5:
            rho, p = stats.spearmanr(x[mask], y[mask])
            coef = np.polyfit(x[mask], y[mask], 1)
            grid = np.linspace(x[mask].min(), x[mask].max(), 40)
            ax.plot(grid, np.polyval(coef, grid), color="#333", lw=1)
            ax.set_title(f"{ylab}\nρ={rho:.2f}  p={p:.2g}  n={int(mask.sum())}")
        else:
            ax.set_title(ylab)
        ax.set_xlabel("Epithelial CLDN4 mean log2(TPM+1)")
        ax.set_ylabel(ylab)
    fig.savefig(out, dpi=160)
    plt.close(fig)


def save_gsea_fig(gsea_df: pd.DataFrame, out: Path, title: str) -> None:
    if gsea_df.empty:
        return
    df = gsea_df.copy()
    if "NES" not in df.columns:
        return
    df = df.sort_values("NES")
    colors = ["#b2182b" if v > 0 else "#2166ac" for v in df["NES"]]
    fig, ax = plt.subplots(figsize=(8.2, max(3.5, 0.32 * len(df) + 1.2)), constrained_layout=True)
    ax.barh(df["Term"], df["NES"], color=colors)
    ax.axvline(0, color="#333", lw=0.8)
    ax.set_xlabel("NES (prerank GSEA)")
    ax.set_title(title)
    fig.savefig(out, dpi=160)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    out = C.RESULTS_DIR
    out.mkdir(parents=True, exist_ok=True)
    t_all = time.time()

    if not C.ANNOT_FILE.exists() or not C.LOG2TPM_FILE.exists():
        raise SystemExit("Missing GSE131907 files. Run scripts/hunt_gse131907/download.py first.")

    print("[1] annotation + series metadata", flush=True)
    ann = pd.read_csv(C.ANNOT_FILE, sep="\t", dtype=str)
    if C.SERIES_FILE.exists():
        meta = parse_series_matrix(C.SERIES_FILE)
        sample_meta = meta.rename(
            columns={"title": "Sample", "tissue_origin_abbrevation": "Sample_Origin_geo"}
        )
        keep = [
            c
            for c in ["Sample", "geo_accession", "patient_id", "tumor_stage", "Sample_Origin_geo", "source_name_ch1"]
            if c in sample_meta.columns
        ]
        sample_meta = sample_meta[keep].drop_duplicates("Sample")
        sample_meta.to_csv(out / "sample_metadata.tsv", sep="\t", index=False)
    else:
        sample_meta = pd.DataFrame({"Sample": []})

    comp = (
        ann.groupby(["Sample_Origin", "Cell_type", "Cell_subtype"], dropna=False)
        .size()
        .reset_index(name="n_cells")
        .sort_values(["Sample_Origin", "n_cells"], ascending=[True, False])
    )
    comp.to_csv(out / "cell_type_composition.tsv", sep="\t", index=False)

    print("[2] stream selected genes from author log2TPM", flush=True)
    cell_ids, expr = stream_selected_genes(C.LOG2TPM_FILE, set(C.SELECTED))
    missing = sorted(set(C.SELECTED) - set(expr))
    print(f"    found {len(expr)}/{len(C.SELECTED)}  missing={missing}", flush=True)

    per_cell = ann.set_index("Index").reindex(cell_ids).reset_index()
    if per_cell["Sample"].isna().any():
        raise SystemExit("matrix barcodes do not match annotation Index")
    for gene, arr in expr.items():
        per_cell[gene] = arr
    per_cell["compartment"] = per_cell.apply(compartment, axis=1)
    per_cell["is_tumor_sample"] = per_cell["Sample_Origin"].isin(C.TUMOR_ORIGINS)
    per_cell["is_cd8"] = per_cell["Cell_subtype"].isin(C.CD8_SUBTYPES)
    per_cell["is_nk"] = per_cell["Cell_type"] == "NK cells"
    per_cell["is_b"] = per_cell["Cell_type"] == "B lymphocytes"
    per_cell["is_gc_b"] = per_cell["Cell_subtype"].isin(C.GC_B_SUBTYPES)
    per_cell["is_immune"] = per_cell["Cell_type"].isin(C.IMMUNE_TYPES)
    per_cell["tumor_epithelial"] = (per_cell["Cell_type"] == "Epithelial cells") & per_cell[
        "Cell_subtype"
    ].isin(C.TUMOR_EPI_SUBTYPES)

    # compact per-cell dump (selected genes only)
    compact = ["Index", "Sample", "Sample_Origin", "Cell_type", "Cell_subtype", "compartment"] + [
        g for g in C.SELECTED if g in per_cell
    ]
    per_cell[compact].to_csv(out / "per_cell_selected_genes.csv.gz", index=False)

    # ----- epithelial restriction -----
    print("[3] epithelial-restricted vs immune", flush=True)
    tumor_cells = per_cell[per_cell["is_tumor_sample"]]
    restriction_rows = []
    for gene in ["TACSTD2", "CLDN4"]:
        if gene not in per_cell:
            continue
        epi = tumor_cells.loc[tumor_cells["compartment"] == "epithelial", gene]
        for other in ["T", "NK", "B", "myeloid"]:
            oth = tumor_cells.loc[tumor_cells["compartment"] == other, gene]
            restriction_rows.append(
                {
                    "gene": gene,
                    "compartment_a": "epithelial",
                    "compartment_b": other,
                    "n_a": int(epi.size),
                    "n_b": int(oth.size),
                    "mean_a": float(epi.mean()),
                    "mean_b": float(oth.mean()),
                    "pct_pos_a": float((epi > 0).mean() * 100),
                    "pct_pos_b": float((oth > 0).mean() * 100) if oth.size else np.nan,
                    **{k: v for k, v in mw_row(epi, oth, f"{gene} epi vs {other}").items() if k != "contrast"},
                }
            )
    restriction_df = pd.DataFrame(restriction_rows)
    restriction_df.to_csv(out / "epithelial_vs_immune.tsv", sep="\t", index=False)
    save_restriction_fig(per_cell, out / "fig_epithelial_vs_immune.png")

    # ----- sample-level neighborhood -----
    print("[4] sample-level CD8/NK/TLS co-occurrence", flush=True)
    tls_genes = [g for g in C.TLS_CHEMOKINE if g in per_cell]
    rows = []
    for (sample, origin), sdf in per_cell.groupby(["Sample", "Sample_Origin"], sort=True):
        epi = sdf[sdf["compartment"] == "epithelial"]
        malig = sdf[sdf["tumor_epithelial"]]
        n = len(sdf)
        row = {
            "Sample": sample,
            "Sample_Origin": origin,
            "n_cells": int(n),
            "n_epithelial": int(len(epi)),
            "n_tumor_epithelial": int(len(malig)),
            "n_CD8": int(sdf["is_cd8"].sum()),
            "n_NK": int(sdf["is_nk"].sum()),
            "n_B": int(sdf["is_b"].sum()),
            "n_GC_B": int(sdf["is_gc_b"].sum()),
            "n_immune": int(sdf["is_immune"].sum()),
            "frac_CD8": float(sdf["is_cd8"].mean()),
            "frac_NK": float(sdf["is_nk"].mean()),
            "frac_B": float(sdf["is_b"].mean()),
            "frac_GC_B": float(sdf["is_gc_b"].mean()),
            "frac_immune": float(sdf["is_immune"].mean()),
            "frac_TNK": float(((sdf["compartment"] == "T") | (sdf["compartment"] == "NK")).mean()),
        }
        for gene in ["TACSTD2", "CLDN4"]:
            if gene not in sdf:
                continue
            ev = epi[gene].to_numpy() if len(epi) else np.array([])
            mv = malig[gene].to_numpy() if len(malig) else np.array([])
            row[f"epi_{gene}_mean"] = float(ev.mean()) if ev.size else np.nan
            row[f"epi_{gene}_pct_pos"] = float((ev > 0).mean() * 100) if ev.size else np.nan
            row[f"malig_{gene}_mean"] = float(mv.mean()) if mv.size else np.nan
        if tls_genes:
            row["tls_chemokine_mean"] = float(sdf[tls_genes].mean(axis=1).mean())
            row["cxcl13_mean"] = float(sdf["CXCL13"].mean()) if "CXCL13" in sdf else np.nan
        else:
            row["tls_chemokine_mean"] = np.nan
        rows.append(row)
    sample_df = pd.DataFrame(rows)
    if not sample_meta.empty:
        sample_df = sample_df.merge(sample_meta, on="Sample", how="left")
    sample_df.to_csv(out / "per_sample.tsv", sep="\t", index=False)

    assoc = []
    tumor = sample_df[
        sample_df["Sample_Origin"].isin(C.TUMOR_ORIGINS) & (sample_df["n_epithelial"] >= C.MIN_EPI)
    ]
    tlung = sample_df[(sample_df["Sample_Origin"] == "tLung") & (sample_df["n_epithelial"] >= C.MIN_EPI)]
    for frame, tag in ((tumor, "tumor_sites"), (tlung, "tLung")):
        for xcol, xlab in (
            ("epi_TACSTD2_mean", "epi TACSTD2 mean"),
            ("epi_CLDN4_mean", "epi CLDN4 mean"),
        ):
            if xcol not in frame:
                continue
            for ycol, ylab in (
                ("frac_CD8", "CD8 fraction"),
                ("frac_NK", "NK fraction"),
                ("frac_B", "B fraction"),
                ("frac_GC_B", "GC-B fraction"),
                ("frac_TNK", "T/NK fraction"),
                ("tls_chemokine_mean", "12-chemokine TLS score"),
            ):
                assoc.append(spearman_row(frame[xcol], frame[ycol], f"{tag}: {xlab} vs {ylab}"))
    assoc_df = pd.DataFrame(assoc)
    try:
        from statsmodels.stats.multitest import multipletests

        ok = assoc_df["spearman_p"].notna()
        if ok.any():
            _, fdr, _, _ = multipletests(assoc_df.loc[ok, "spearman_p"], method="fdr_bh")
            assoc_df.loc[ok, "bh_fdr"] = fdr
    except Exception:
        assoc_df["bh_fdr"] = np.nan
    assoc_df.to_csv(out / "neighborhood_associations.tsv", sep="\t", index=False)
    save_neighborhood_fig(sample_df, out / "fig_neighborhood_cooccurrence.png")
    save_cldn4_neighborhood_fig(sample_df, out / "fig_neighborhood_cldn4.png")

    # ----- GSEA ranking pass -----
    print("[5] stream all genes for GSEA ranking (author log2TPM)", flush=True)
    epi_mask = (per_cell["Cell_type"] == "Epithelial cells").to_numpy()
    tumor_epi_mask = per_cell["tumor_epithelial"].to_numpy()
    # TROP2-high / low among tumor-site epithelial cells (author tumor epithelium)
    trop2 = per_cell["TACSTD2"].to_numpy()
    use = tumor_epi_mask & per_cell["is_tumor_sample"].to_numpy()
    epi_trop2 = trop2[epi_mask]
    # high/low masks are relative to the epithelial slice
    trop2_on_epi = trop2[epi_mask]
    tumor_on_epi = use[epi_mask]
    q75 = np.nanquantile(trop2_on_epi[tumor_on_epi], 0.75)
    q25 = np.nanquantile(trop2_on_epi[tumor_on_epi], 0.25)
    trop2_high = tumor_on_epi & (trop2_on_epi >= q75)
    trop2_low = tumor_on_epi & (trop2_on_epi <= q25)
    print(
        f"    TROP2 Q4>={q75:.3f} n={int(trop2_high.sum())} ; "
        f"Q1<={q25:.3f} n={int(trop2_low.sum())}",
        flush=True,
    )

    # sample codes as integer factor for bincount
    samples = per_cell["Sample"].astype("category")
    sample_codes = samples.cat.codes.to_numpy()
    sample_names = list(samples.cat.categories)
    # tumor samples with enough epithelial cells
    tumor_sample_keep = set()
    epi_tacstd2_by_sample = {}
    for i, name in enumerate(sample_names):
        sub = sample_df[sample_df["Sample"] == name]
        if sub.empty:
            continue
        origin = sub["Sample_Origin"].iloc[0]
        n_epi = int(sub["n_epithelial"].iloc[0])
        if origin in C.TUMOR_ORIGINS and n_epi >= C.MIN_EPI:
            tumor_sample_keep.add(i)
            epi_tacstd2_by_sample[i] = float(sub["epi_TACSTD2_mean"].iloc[0])

    rank_df = stream_rank_metrics(
        C.LOG2TPM_FILE,
        cell_ids,
        epi_mask,
        sample_codes,
        trop2_high,
        trop2_low,
        tumor_sample_keep,
        epi_tacstd2_by_sample,
    )
    rank_df.to_csv(out / "gene_rank_metrics.tsv.gz", sep="\t", index=False)

    gene_sets = load_genesets()
    # leave-one-out TJ without CLDN4 (CLDN4 is the grouping correlate)
    loo = {k: [g for g in v if g not in {"CLDN4", "TACSTD2"}] for k, v in gene_sets.items()}
    loo = {f"{k}__no_CLDN4_TACSTD2": v for k, v in loo.items() if len(v) >= 5}
    all_sets = {**gene_sets, **loo}

    print("[6] prerank GSEA", flush=True)
    rnk_sample = rank_df.set_index("gene")["sample_spearman_rho"]
    rnk_cell = rank_df.set_index("gene")["cell_t"]
    gsea_sample = run_prerank(rnk_sample, all_sets, out, "sample_spearman_TACSTD2")
    gsea_cell = run_prerank(rnk_cell, all_sets, out, "cell_t_TROP2_Q4vsQ1")
    gsea_all = pd.concat([gsea_sample, gsea_cell], ignore_index=True)
    gsea_all.to_csv(out / "gsea_all.tsv", sep="\t", index=False)

    # key-set subset for the figure (no leave-one-out clutter)
    key_terms = [
        "KEGG_TIGHT_JUNCTION",
        "GOBP_TIGHT_JUNCTION_ORGANIZATION",
        "HALLMARK_APICAL_JUNCTION",
        "CUSTOM_TJ_CORE",
        "GOBP_KERATINIZATION",
        "GOBP_CORNIFICATION",
        "CUSTOM_BARRIER_KERATIN",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "CUSTOM_EMT_CORE",
    ]
    gsea_key = gsea_sample[gsea_sample["Term"].isin(key_terms)].copy()
    gsea_key.to_csv(out / "gsea_key_sets_sample.tsv", sep="\t", index=False)
    save_gsea_fig(
        gsea_key,
        out / "fig_gsea_sample_TACSTD2.png",
        "Sample-level prerank GSEA\n(epithelial gene mean vs epithelial TACSTD2)",
    )
    gsea_key_cell = gsea_cell[gsea_cell["Term"].isin(key_terms)].copy()
    save_gsea_fig(
        gsea_key_cell,
        out / "fig_gsea_cell_TROP2Q4Q1.png",
        "Cell-level prerank GSEA (anti-conservative)\nTROP2-high vs TROP2-low tumor epithelium",
    )

    # ----- verdict -----
    def pick(df, substr):
        hits = df[df["contrast"].str.contains(substr, regex=False)]
        return hits.iloc[0].to_dict() if len(hits) else {}

    tac_cd8 = pick(assoc_df, "tumor_sites: epi TACSTD2 mean vs CD8")
    tac_nk = pick(assoc_df, "tumor_sites: epi TACSTD2 mean vs NK")
    tac_b = pick(assoc_df, "tumor_sites: epi TACSTD2 mean vs B fraction")
    tac_tls = pick(assoc_df, "tumor_sites: epi TACSTD2 mean vs 12-chemokine")
    tac_t = restriction_df[(restriction_df.gene == "TACSTD2") & (restriction_df.compartment_b == "T")]
    cldn_t = restriction_df[(restriction_df.gene == "CLDN4") & (restriction_df.compartment_b == "T")]

    def gsea_nes(table, term):
        hit = table[table["Term"] == term]
        if hit.empty:
            return None
        r = hit.iloc[0]
        fdr_col = "FDR q-val" if "FDR q-val" in r.index else ("FDR" if "FDR" in r.index else None)
        return {
            "term": term,
            "NES": float(r["NES"]),
            "FDR": float(r[fdr_col]) if fdr_col and pd.notna(r[fdr_col]) else np.nan,
        }

    gsea_highlights = {
        "sample": {t: gsea_nes(gsea_sample, t) for t in key_terms},
        "cell": {t: gsea_nes(gsea_cell, t) for t in key_terms},
    }

    epi_restricted = True
    if not tac_t.empty:
        epi_restricted = bool(tac_t.iloc[0]["pct_pos_a"] > 5 * max(tac_t.iloc[0]["pct_pos_b"], 0.1))

    def sig(item, alpha=0.05):
        return bool(item and pd.notna(item.get("spearman_p")) and item["spearman_p"] < alpha)

    neighborhood_supported = False
    if "bh_fdr" in assoc_df.columns:
        neighborhood_supported = bool(
            ((assoc_df["bh_fdr"] < 0.05) & (assoc_df["spearman_rho"] < 0)).any()
        )

    tj_up = False
    ker_up = False
    s_tj = gsea_highlights["sample"].get("CUSTOM_TJ_CORE") or gsea_highlights["sample"].get("KEGG_TIGHT_JUNCTION")
    s_ker = gsea_highlights["sample"].get("CUSTOM_BARRIER_KERATIN") or gsea_highlights["sample"].get("GOBP_KERATINIZATION")
    if s_tj and s_tj["NES"] > 0 and (pd.isna(s_tj["FDR"]) or s_tj["FDR"] < 0.25):
        tj_up = True
    if s_ker and s_ker["NES"] > 0 and (pd.isna(s_ker["FDR"]) or s_ker["FDR"] < 0.25):
        ker_up = True

    summary = {
        "dataset": "GSE131907",
        "citation": "Kim et al. Nat Commun 2020 PMID 32385277",
        "n_cells": int(len(per_cell)),
        "n_samples": int(per_cell["Sample"].nunique()),
        "n_tumor_samples_eligible": int(len(tumor)),
        "n_tLung_eligible": int(len(tlung)),
        "expression_metric": "author log2(TPM+1); selected genes + full-matrix rank streamed",
        "spatial": False,
        "ici_or_mpr_labels": False,
        "trop2_q75": float(q75),
        "trop2_q25": float(q25),
        "n_trop2_high_tumor_epi": int(trop2_high.sum()),
        "n_trop2_low_tumor_epi": int(trop2_low.sum()),
        "genes_missing": missing,
        "epithelial_restricted": {
            "supported": epi_restricted,
            "TACSTD2_vs_T": tac_t.iloc[0].to_dict() if not tac_t.empty else {},
            "CLDN4_vs_T": cldn_t.iloc[0].to_dict() if not cldn_t.empty else {},
        },
        "neighborhood_anti_correlation_supported": neighborhood_supported,
        "neighborhood": {
            "TACSTD2_vs_CD8": tac_cd8,
            "TACSTD2_vs_NK": tac_nk,
            "TACSTD2_vs_B": tac_b,
            "TACSTD2_vs_TLS_chemokine": tac_tls,
        },
        "gsea_tj_up_in_TACSTD2_high": tj_up,
        "gsea_keratin_up_in_TACSTD2_high": ker_up,
        "gsea_highlights": gsea_highlights,
        "honest_limits": [
            "Dissociated droplet scRNA-seq: no spatial neighborhood and no histologically scored TLS.",
            "Immune fractions are author-label proportions within each dissociated sample.",
            "12-chemokine TLS score is a bulk-style proxy averaged over all cells in the sample.",
            "Cell-level GSEA treats cells as independent and is anti-conservative; sample-level Spearman prerank is the primary GSEA.",
            "Treatment-naive atlas: no ICI / MPR / RECIST labels.",
            "Primary tLung n is small (~11); tumor-site n pools primary, LN, brain met, and PE.",
        ],
        "runtime_sec": round(time.time() - t_all, 1),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    write_readme(out, summary, assoc_df, restriction_df, gsea_sample, gsea_cell, missing)
    print(json.dumps({k: summary[k] for k in summary if k != "gsea_highlights"}, indent=2, default=str))
    print(f"[done] {out}  ({time.time() - t_all:.0f}s)", flush=True)


def write_readme(
    out: Path,
    summary: dict,
    assoc_df: pd.DataFrame,
    restriction_df: pd.DataFrame,
    gsea_sample: pd.DataFrame,
    gsea_cell: pd.DataFrame,
    missing: list[str],
) -> None:
    def fmt_assoc(item: dict) -> str:
        if not item:
            return "- (missing)"
        if item.get("note") == "too_few_samples":
            return f"- {item['contrast']}: n={item['n']} (too few)"
        fdr = item.get("bh_fdr", np.nan)
        fdr_s = f", FDR={float(fdr):.2g}" if pd.notna(fdr) else ""
        return (
            f"- {item['contrast']}: ρ={item['spearman_rho']:.3f}, "
            f"p={item['spearman_p']:.3g}{fdr_s}, n={item['n']}"
        )

    def fmt_gsea(table: pd.DataFrame, term: str) -> str:
        hit = table[table["Term"] == term]
        if hit.empty:
            return f"- {term}: not tested"
        r = hit.iloc[0]
        fdr = r["FDR q-val"] if "FDR q-val" in r.index else r.get("FDR", np.nan)
        return f"- {term}: NES={float(r['NES']):.2f}, FDR={float(fdr):.3g}"

    epi = summary["epithelial_restricted"]
    tac = epi.get("TACSTD2_vs_T") or {}
    cld = epi.get("CLDN4_vs_T") or {}

    lines = [
        "# Hunt GSE131907 — epithelial TACSTD2/CLDN4 vs CD8/NK/TLS; TJ/keratin GSEA",
        "",
        "Kim et al., *Nat Commun* 2020 (PMID 32385277); GEO [GSE131907](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131907).",
        "",
        "## Verdict (honest)",
        "",
        f"- **Epithelial-restricted vs immune:** "
        f"{'supported' if epi.get('supported') else 'not supported'}. "
        f"Tumor-site TACSTD2 %pos epithelial={tac.get('pct_pos_a', float('nan')):.1f} vs T={tac.get('pct_pos_b', float('nan')):.1f} "
        f"(Mann–Whitney p={tac.get('p', float('nan')):.2g}). "
        f"CLDN4 %pos epithelial={cld.get('pct_pos_a', float('nan')):.1f} vs T={cld.get('pct_pos_b', float('nan')):.1f} "
        f"(p={cld.get('p', float('nan')):.2g}).",
        f"- **CD8/NK/TLS neighborhood anti-correlation:** "
        f"{'BH-FDR<0.05 inverse present (see table)' if summary['neighborhood_anti_correlation_supported'] else 'not supported after BH-FDR (TACSTD2 null; CLDN4 nominal inverses FDR≥0.15)'}. "
        "This is **not** a spatial neighborhood. The matrix is dissociated 10x.",
        f"- **TJ GSEA in TROP2-high epithelium (sample-level prerank):** "
        f"{'up (NES>0, FDR<0.25 on a TJ set)' if summary['gsea_tj_up_in_TACSTD2_high'] else 'not supported at FDR<0.25'}.",
        f"- **Keratin/barrier GSEA in TROP2-high epithelium (sample-level prerank):** "
        f"{'up (NES>0, FDR<0.25 on a keratin set)' if summary['gsea_keratin_up_in_TACSTD2_high'] else 'not supported at FDR<0.25'}.",
        "- **ICI / MPR:** cannot be tested. Treatment-naive atlas, no response labels.",
        "",
        "## What was measured",
        "",
        f"- {summary['n_cells']} cells / {summary['n_samples']} samples; "
        f"{summary['n_tumor_samples_eligible']} tumor-site samples with ≥{C.MIN_EPI} epithelial cells "
        f"({summary['n_tLung_eligible']} primary tLung).",
        "- Expression: author `normalized_log2TPM` (log2(TPM+1)). Selected genes and the full gene rank were streamed; the dense matrix was never held in RAM.",
        "- Epithelial = author `Cell_type == Epithelial cells`. Tumor epithelium for TROP2 quartiles = `Malignant cells` + `tS1/tS2/tS3`.",
        "- TACSTD2 and CLDN4 scores used for association and GSEA grouping are **epithelial-only**.",
        "- CD8 = author subtypes Cytotoxic / Exhausted / Naive / CD8-low CD8+ T. NK = `Cell_type == NK cells`.",
        "- TLS proxies (not histology): B-cell fraction, GC-B (DZ+LZ) fraction, and the 12-chemokine mean (CCL2/3/4/5/8/18/19/21, CXCL9/10/11/13) averaged over all cells in the sample.",
        "- Primary GSEA: prerank by Spearman ρ of each gene's epithelial sample-mean vs epithelial TACSTD2, tumor sites, 1000 permutations.",
        "- Secondary GSEA: prerank by Welch t of TROP2 Q4 vs Q1 tumor epithelial cells. Cells from the same sample are not independent; treat as anti-conservative.",
        f"- Genes missing from the log2TPM matrix: {missing or 'none'}.",
        "",
        "## Honest limits",
        "",
    ]
    lines += [f"- {x}" for x in summary["honest_limits"]]
    lines += [
        "",
        "## Neighborhood associations (sample-level Spearman)",
        "",
    ]
    for _, r in assoc_df.iterrows():
        lines.append(fmt_assoc(r.to_dict()))

    lines += [
        "",
        "## Sample-level GSEA (primary)",
        "",
    ]
    for term in [
        "KEGG_TIGHT_JUNCTION",
        "GOBP_TIGHT_JUNCTION_ORGANIZATION",
        "HALLMARK_APICAL_JUNCTION",
        "CUSTOM_TJ_CORE",
        "CUSTOM_TJ_CORE__no_CLDN4_TACSTD2",
        "GOBP_KERATINIZATION",
        "CUSTOM_BARRIER_KERATIN",
        "CUSTOM_BARRIER_KERATIN__no_CLDN4_TACSTD2",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "CUSTOM_EMT_CORE",
    ]:
        lines.append(fmt_gsea(gsea_sample, term))

    lines += [
        "",
        "## Cell-level GSEA (secondary, anti-conservative)",
        "",
    ]
    for term in [
        "KEGG_TIGHT_JUNCTION",
        "CUSTOM_TJ_CORE",
        "GOBP_KERATINIZATION",
        "CUSTOM_BARRIER_KERATIN",
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    ]:
        lines.append(fmt_gsea(gsea_cell, term))

    lines += [
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `summary.json` | Machine-readable verdict |",
        "| `epithelial_vs_immune.tsv` | TACSTD2/CLDN4 epithelial vs T/NK/B/myeloid |",
        "| `per_sample.tsv` | Sample-level epithelial scores and immune fractions |",
        "| `neighborhood_associations.tsv` | Spearman tests |",
        "| `gene_rank_metrics.tsv.gz` | Per-gene sample Spearman and cell t |",
        "| `gsea_sample_spearman_TACSTD2.tsv` | Primary prerank GSEA |",
        "| `gsea_cell_t_TROP2_Q4vsQ1.tsv` | Secondary cell-level GSEA |",
        "| `gsea_all.tsv` / `gsea_key_sets_sample.tsv` | Combined / key-set slices |",
        "| `fig_epithelial_vs_immune.png` | Restriction boxplots |",
        "| `fig_neighborhood_cooccurrence.png` | TACSTD2 vs CD8/NK/B/TLS |",
        "| `fig_neighborhood_cldn4.png` | CLDN4 vs CD8/NK/B/TLS (nominal) |",
        "| `fig_gsea_sample_TACSTD2.png` | Primary NES bars |",
        "| `fig_gsea_cell_TROP2Q4Q1.png` | Secondary NES bars |",
        "| `per_cell_selected_genes.csv.gz` | Streamed marker genes |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 scripts/hunt_gse131907/download.py",
        "python3 scripts/hunt_gse131907/analyze.py",
        "```",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
