"""Depth control for the CLDN4-high vs low contrast.

CLDN4 UMI>0 cells in this object are much deeper than CLDN4 UMI=0 cells.
The caliper match is a sensitivity added after that gap was measured, before
the matched test statistics were used as a biological claim.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


def match_pairs(meta: pd.DataFrame, caliper: float, min_pairs: int) -> tuple[np.ndarray, pd.DataFrame]:
    """Greedy 1:1 match on log1p(nUMI) within sample. Deterministic."""
    keep = np.zeros(len(meta), dtype=bool)
    rows = []
    for sample, sub in meta.groupby("sample", sort=True):
        high = sub[sub["cldn4_umi"] > 0].copy()
        low = sub[sub["cldn4_umi"] <= 0].copy()
        high["_lu"] = np.log1p(high["n_umi"].to_numpy())
        low["_lu"] = np.log1p(low["n_umi"].to_numpy())
        high = high.sort_values(["_lu", "barcode"])
        low = low.sort_values(["barcode"])
        used: set = set()
        pairs = []
        for ix, rec in high.iterrows():
            cand = low.loc[~low.index.isin(used)]
            if cand.empty:
                break
            delta = (cand["_lu"] - rec["_lu"]).abs()
            j = delta.idxmin()
            d = float(delta.loc[j])
            if d <= caliper:
                used.add(j)
                pairs.append((ix, j, d, float(rec["_lu"]), float(cand.loc[j, "_lu"])))
        if len(pairs) < min_pairs:
            continue
        for ix, j, d, lu_h, lu_l in pairs:
            keep[int(ix)] = True
            keep[int(j)] = True
            rows.append(
                {
                    "sample": sample,
                    "origin": sub["origin"].iloc[0],
                    "barcode_high": meta.at[ix, "barcode"],
                    "barcode_low": meta.at[j, "barcode"],
                    "abs_d_log1p_umi": d,
                    "log1p_umi_high": lu_h,
                    "log1p_umi_low": lu_l,
                    "n_umi_high": float(meta.at[ix, "n_umi"]),
                    "n_umi_low": float(meta.at[j, "n_umi"]),
                }
            )
    return keep, pd.DataFrame(rows)


def depth_summary(meta: pd.DataFrame, paired_means_fn) -> pd.DataFrame:
    high = meta["cldn4_umi"].to_numpy() > 0
    sample = meta["sample"].to_numpy()
    rows = []
    for col, label in (("n_umi", "n_umi"), ("n_genes", "n_genes")):
        # report log1p for UMI, raw for n_genes
        values = np.log1p(meta[col].to_numpy()) if col == "n_umi" else meta[col].to_numpy().astype(float)
        test = paired_means_fn(values, high, sample)
        rows.append(
            {
                "measure": "log1p_n_umi" if col == "n_umi" else "n_genes",
                "n_samples": test["n"],
                "median_delta_high_minus_low": test["median"],
                "n_samples_delta_pos": test["n_pos"],
                "p_wilcoxon": test["p"],
            }
        )
    # arm medians
    arm_rows = []
    for s, sub in meta.groupby("sample", sort=True):
        h = sub[sub["cldn4_umi"] > 0]
        lo = sub[sub["cldn4_umi"] <= 0]
        if len(h) < 5 or len(lo) < 5:
            continue
        arm_rows.append(
            {
                "sample": s,
                "origin": sub["origin"].iloc[0],
                "n_high": int(len(h)),
                "n_low": int(len(lo)),
                "median_umi_high": float(h["n_umi"].median()),
                "median_umi_low": float(lo["n_umi"].median()),
                "median_genes_high": float(h["n_genes"].median()),
                "median_genes_low": float(lo["n_genes"].median()),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(arm_rows)


def _program_table(scores, n_used, high, sample, paired_means_fn, bh_fn, program_order, controls):
    rows = []
    for name in program_order + controls:
        test = paired_means_fn(scores[name], high, sample)
        rows.append(
            {
                "family": "control" if name in controls else "program",
                "program": name,
                "n_genes_in_matrix": n_used[name],
                "n_samples": test["n"],
                "median_delta_high_minus_low": test["median"],
                "mean_delta": test["mean"],
                "n_samples_delta_pos": test["n_pos"],
                "n_samples_delta_neg": test["n_neg"],
                "p_wilcoxon": test["p"],
            }
        )
    df = pd.DataFrame(rows)
    mask = df["family"].eq("program")
    df["q_bh"] = np.nan
    df.loc[mask, "q_bh"] = bh_fn(df.loc[mask, "p_wilcoxon"].tolist())
    return df


def run_depth_match(
    meta: pd.DataFrame,
    counts: np.ndarray,
    lognorm: np.ndarray,
    symbols: list[str],
    ensembl: list[str],
    programs: dict,
    raw_embeddings: dict[str, np.ndarray],
    *,
    tables: Path,
    figs: Path,
    caliper: float,
    min_pairs: int,
    helpers: dict,
):
    paired_means_fn = helpers["paired_means"]
    bh_fn = helpers["bh"]
    loso = helpers["loso_projection"]
    mean_rho = helpers["mean_gene_rho"]
    diffusion_fn = helpers["diffusion_map"]
    orient_fn = helpers["orient_components"]
    plot_deltas_fn = helpers["plot_deltas"]
    plot_paired_fn = helpers["plot_paired"]
    plot_genes_fn = helpers["plot_genes"]
    write = helpers["write_tsv"]
    program_order = helpers["PROGRAM_ORDER"]
    controls = helpers["CONTROLS"]
    seed = helpers["SEED"]

    summary, arms = depth_summary(meta, paired_means_fn)
    write(summary, tables / "depth_gap_tests.tsv")
    write(arms, tables / "depth_gap_by_sample.tsv")

    keep, pairs = match_pairs(meta, caliper, min_pairs)
    write(pairs, tables / "depth_caliper_pairs.tsv")
    print(
        f"depth caliper {caliper}: {len(pairs)} pairs, "
        f"{pairs['sample'].nunique() if len(pairs) else 0} samples, "
        f"{int(keep.sum())} cells",
        flush=True,
    )
    if pairs.empty:
        return

    sub_meta = meta.loc[keep].reset_index(drop=True)
    sub_counts = counts[keep]
    sub_log = lognorm[keep]
    high = sub_meta["cldn4_umi"].to_numpy() > 0
    sample = sub_meta["sample"].to_numpy()
    origin_of = sub_meta.groupby("sample")["origin"].first().to_dict()
    scores, n_used = helpers["program_scores"](sub_log, symbols, programs)
    prog = _program_table(scores, n_used, high, sample, paired_means_fn, bh_fn, program_order, controls)
    write(prog, tables / "program_paired_tests_matched.tsv")
    mask = prog["family"].eq("program")
    plot_deltas_fn(
        prog.loc[mask, "program"].tolist(),
        prog.loc[mask, "median_delta_high_minus_low"].tolist(),
        prog.loc[mask, "q_bh"].tolist(),
        figs / "fig_program_paired_matched.png",
        f"Depth-matched program scores (caliper {caliper} on log1p nUMI)",
        "Median paired difference of mean log1p(CP10k)",
    )

    # PCA + diffusion refit on the matched cells, CLDN4 held out.
    cldn_col = symbols.index("CLDN4") if "CLDN4" in symbols else None
    detected = (sub_counts > 0).sum(axis=0)
    hvg_ok = detected >= 10
    if cldn_col is not None:
        hvg_ok[cldn_col] = False
    var = sub_log[:, hvg_ok].var(axis=0)
    n_hvg = min(2000, int(hvg_ok.sum()))
    hvg_local = np.flatnonzero(hvg_ok)[np.argsort(var)[-n_hvg:]]
    x = sub_log[:, hvg_local].astype(np.float64)
    x = (x - x.mean(0)) / np.clip(x.std(0), 1e-6, None)
    n_pcs = min(30, x.shape[0] - 1, x.shape[1])
    pca = PCA(n_components=n_pcs, svd_solver="full", random_state=seed)
    x_pca = pca.fit_transform(x)
    n_show = min(10, n_pcs)
    evals, diff = diffusion_fn(x_pca, n_comps=n_show + 1, n_neighbors=15)
    x_pca_o = orient_fn(x_pca[:, :n_show], high, sample)
    n_dc = min(n_show, diff.shape[1])
    diff_o = orient_fn(diff[:, :n_dc], high, sample)
    geom_rows = []
    for kind, arr in (("PC", x_pca_o), ("DC", diff_o)):
        tests = [paired_means_fn(arr[:, j], high, sample) for j in range(arr.shape[1])]
        qs = bh_fn([t["p"] for t in tests])
        for j, (t, q) in enumerate(zip(tests, qs)):
            geom_rows.append(
                {
                    "space": kind,
                    "component": f"{kind}{j + 1}",
                    "n_samples": t["n"],
                    "median_delta_high_minus_low": t["median"],
                    "n_samples_delta_pos": t["n_pos"],
                    "p_wilcoxon": t["p"],
                    "q_bh": q,
                    "variance_explained": float(pca.explained_variance_ratio_[j]) if kind == "PC" else np.nan,
                }
            )
    geom = pd.DataFrame(geom_rows)
    geom["diffusion_eigenvalue_trivial"] = float(evals[0]) if len(evals) else np.nan
    geom["n_cells"] = int(keep.sum())
    geom["n_pairs"] = int(len(pairs))
    write(geom, tables / "pca_diffusion_paired_tests_matched.tsv")
    plot_deltas_fn(
        geom["component"].tolist(),
        geom["median_delta_high_minus_low"].tolist(),
        geom["q_bh"].tolist(),
        figs / "fig_pca_diffusion_paired_matched.png",
        "Depth-matched PCA / diffusion (CLDN4 held out)",
        "Median paired difference after orienting median ≥ 0",
    )

    embed_rows = []
    holdout_proj = None
    for name, emb in raw_embeddings.items():
        if emb is None or emb.shape[0] != len(meta):
            continue
        sub = emb[keep]
        proj, _ = loso(sub, high, sample)
        test = paired_means_fn(proj, high, sample)
        embed_rows.append(
            {
                "embedding": name,
                "n_samples": test["n"],
                "n_pairs": int(len(pairs)),
                "median_delta_high_minus_low": test["median"],
                "mean_delta": test["mean"],
                "n_samples_delta_pos": test["n_pos"],
                "n_samples_delta_neg": test["n_neg"],
                "p_wilcoxon": test["p"],
            }
        )
        if name == "gf_v1_holdout":
            holdout_proj = proj
            plot_paired_fn(
                test,
                origin_of,
                figs / "fig_geneformer_v1_loso_paired_matched.png",
                f"Geneformer V1 holdout, depth-matched  n={test['n']}",
                "Leave-one-sample-out projection",
            )
    embed_df = pd.DataFrame(embed_rows)
    if len(embed_df):
        primary = embed_df["embedding"].isin(["gf_v1_holdout", "scgpt_holdout"])
        embed_df["q_bh_holdouts"] = np.nan
        if primary.any():
            embed_df.loc[primary, "q_bh_holdouts"] = bh_fn(embed_df.loc[primary, "p_wilcoxon"].tolist())
        write(embed_df, tables / "embedding_loso_tests_matched.tsv")

    if holdout_proj is not None:
        axis_rows = []
        for name in program_order + controls:
            rhos = []
            for s in pd.unique(sample):
                m = (sample == s) & np.isfinite(holdout_proj) & np.isfinite(scores[name])
                if int((m & high).sum()) < 5 or int((m & ~high).sum()) < 5:
                    continue
                if np.std(holdout_proj[m]) < 1e-8 or np.std(scores[name][m]) < 1e-8:
                    continue
                from scipy.stats import rankdata

                rho = np.corrcoef(rankdata(holdout_proj[m]), rankdata(scores[name][m]))[0, 1]
                rhos.append(float(rho))
            rhos_a = np.asarray(rhos, dtype=float)
            p = np.nan
            if rhos_a.size >= 6 and np.any(rhos_a != 0):
                import warnings
                from scipy.stats import wilcoxon

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    p = float(wilcoxon(rhos_a, alternative="two-sided", zero_method="wilcox").pvalue)
            elif rhos_a.size >= 6:
                p = 1.0
            axis_rows.append(
                {
                    "program": name,
                    "family": "control" if name in controls else "program",
                    "n_samples": int(rhos_a.size),
                    "median_spearman_vs_gf_holdout_loso": float(np.median(rhos_a)) if rhos_a.size else np.nan,
                    "p_wilcoxon": p,
                }
            )
        axis = pd.DataFrame(axis_rows)
        axis["q_bh"] = np.nan
        msk = axis["family"].eq("program")
        axis.loc[msk, "q_bh"] = bh_fn(axis.loc[msk, "p_wilcoxon"].tolist())
        write(axis, tables / "program_vs_geneformer_axis_matched.tsv")

        detect = (sub_counts > 0).mean(axis=0)
        rho_all = mean_rho(holdout_proj, sub_log, sample, high)
        gene_tbl = pd.DataFrame(
            {
                "symbol": symbols,
                "ensembl": ensembl,
                "mean_within_sample_spearman": rho_all,
                "detection_frac": detect,
            }
        )
        gene_tbl = gene_tbl[detect >= 0.05].sort_values("mean_within_sample_spearman", ascending=False)
        write(gene_tbl.head(40), tables / "nearest_genes_expression_axis_top_matched.tsv")
        write(
            gene_tbl.tail(40).sort_values("mean_within_sample_spearman"),
            tables / "nearest_genes_expression_axis_bottom_matched.tsv",
        )
        plot_genes_fn(
            gene_tbl.head(20),
            figs / "fig_nearest_genes_axis_matched.png",
            "Genes tracking the depth-matched Geneformer axis",
            "mean_within_sample_spearman",
        )
