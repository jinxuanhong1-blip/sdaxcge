#!/usr/bin/env python3
"""GSE243013 public TCR + MPR + residual TACSTD2/CLDN4 (immune-only).

Unit = sample/patient. Combinatorial: expanded vs non-expanded.
Immune-compartment TACSTD2/CLDN4 is residual/ambient detection, not a tumor score.
E2 (tumor-epithelial TACSTD2 vs expanded CXCL13+) is not estimable.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from stats_util import mw_row, paired_mw, spearman_row, write_json, write_tsv

# Author names that literally contain CXCL13. There is no CD8T_Tex_CXCL13
# in the public TCR table (verified).
CXCL13_CLUSTERS = ("CD4T_Tfh_CXCL13", "CD4T_Th1-like_CXCL13")
CD8_TEX_CLUSTERS = ("CD8T_Tex_HAVCR2", "CD8T_terminal_Tex_LAYN")
CD8_SUBSTR = "CD8T"


def bin_response(x: str, rule: str) -> str:
    sl = str(x).strip().lower().replace(" ", "")
    if sl in {"nan", "none", "", "unknowm", "unknown"}:
        return "unknown"
    is_pcr = sl == "pcr"
    is_mpr = sl == "mpr"
    is_non = sl in {"non-mpr", "nonmpr", "nmpr"}
    if rule == "pcr_or_mpr_vs_nonmpr":
        if is_pcr or is_mpr:
            return "MPR-any"
        if is_non:
            return "non-MPR"
    elif rule == "pcr_vs_nonpcr":
        if is_pcr:
            return "pCR"
        if is_mpr or is_non:
            return "non-pCR"
    elif rule == "mpr_only_vs_nonmpr":
        if is_pcr:
            return "drop"
        if is_mpr:
            return "MPR-only"
        if is_non:
            return "non-MPR"
    return "unknown"


def sample_meta(meta: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "sampleID",
        "cancer_type",
        "pathological_response",
        "pathological_response_rate",
        "radiological_response",
        "anti-PD1_therapy",
        "chemotherapy",
        "gender",
        "age",
    ]
    have = [c for c in cols if c in meta.columns]
    s = meta[have].drop_duplicates("sampleID")
    s["path_bin"] = s["pathological_response"].map(lambda x: bin_response(x, "pcr_or_mpr_vs_nonmpr"))
    s["path_pcr"] = s["pathological_response"].map(lambda x: bin_response(x, "pcr_vs_nonpcr"))
    s["path_mpr_only"] = s["pathological_response"].map(lambda x: bin_response(x, "mpr_only_vs_nonmpr"))
    return s


def endpoints(tcr: pd.DataFrame, min_tcr: int) -> pd.DataFrame:
    df = tcr.copy()
    df["clone_id"] = df["TRA_cdr3"].fillna("") + "|" + df["TRB_cdr3"].fillna("")
    df = df[df["clone_id"] != "|"].copy()
    df["is_cxcl13"] = df["sub_cell_type"].isin(CXCL13_CLUSTERS)
    df["is_cd8"] = df["sub_cell_type"].astype(str).str.startswith(CD8_SUBSTR)
    df["is_cd8_tex"] = df["sub_cell_type"].isin(CD8_TEX_CLUSTERS)
    df["clone_size"] = df.groupby(["sampleID", "clone_id"], observed=True)["cellID"].transform("size")
    df["expanded"] = df["clone_size"] >= 2
    df["expanded_ge3"] = df["clone_size"] >= 3
    df["author_expanded"] = df["expansion"].astype(str).str.lower().eq("expanded")

    rows = []
    for sid, g in df.groupby("sampleID", observed=True):
        n_tcr = int(len(g))
        n_cd8 = int(g["is_cd8"].sum())
        n_cx = int(g["is_cxcl13"].sum())
        n_tex = int(g["is_cd8_tex"].sum())
        n_exp = int(g["expanded"].sum())
        n_non = n_tcr - n_exp
        n_exp_cx = int((g["is_cxcl13"] & g["expanded"]).sum())
        n_non_cx = int((g["is_cxcl13"] & ~g["expanded"]).sum())
        n_exp_tex = int((g["is_cd8_tex"] & g["expanded"]).sum())
        n_non_tex = int((g["is_cd8_tex"] & ~g["expanded"]).sum())
        n_exp_cx_clones = int(g.loc[g["is_cxcl13"] & g["expanded"], "clone_id"].nunique())
        n_exp_tex_clones = int(g.loc[g["is_cd8_tex"] & g["expanded"], "clone_id"].nunique())
        n_clones = int(g["clone_id"].nunique())
        n_exp_clones = int(g.loc[g["expanded"], "clone_id"].nunique())
        # combinatorial rates (cell-level within patient, then tested at patient level)
        frac_cx_exp = (n_exp_cx / n_exp) if n_exp else np.nan
        frac_cx_non = (n_non_cx / n_non) if n_non else np.nan
        frac_tex_exp = (n_exp_tex / n_exp) if n_exp else np.nan
        frac_tex_non = (n_non_tex / n_non) if n_non else np.nan
        # log odds ratio CXCL13 | expanded vs not (Haldane-Anscombe)
        a, b = n_exp_cx + 0.5, (n_exp - n_exp_cx) + 0.5
        c, d = n_non_cx + 0.5, (n_non - n_non_cx) + 0.5
        lor_cx = float(np.log((a / b) / (c / d))) if n_exp and n_non else np.nan
        rows.append(
            {
                "sampleID": sid,
                "n_tcr": n_tcr,
                "n_cd8": n_cd8,
                "n_clones": n_clones,
                "n_exp_clones": n_exp_clones,
                "frac_exp_cells": n_exp / n_tcr if n_tcr else np.nan,
                "frac_exp_cells_ge3": float(g["expanded_ge3"].mean()) if n_tcr else np.nan,
                "frac_author_exp": float(g["author_expanded"].mean()) if n_tcr else np.nan,
                "n_cxcl13": n_cx,
                "frac_cxcl13": n_cx / n_tcr if n_tcr else np.nan,
                "n_exp_cxcl13_cells": n_exp_cx,
                "n_exp_cxcl13_clones": n_exp_cx_clones,
                "exp_cxcl13_per_k_cd8": (n_exp_cx_clones / n_cd8 * 1000.0) if n_cd8 else np.nan,
                "exp_cxcl13_per_k_tcr": (n_exp_cx_clones / n_tcr * 1000.0) if n_tcr else np.nan,
                "frac_cxcl13_that_are_exp": (n_exp_cx / n_cx) if n_cx else np.nan,
                "frac_cxcl13_among_expanded": frac_cx_exp,
                "frac_cxcl13_among_nonexpanded": frac_cx_non,
                "lor_cxcl13_exp_vs_non": lor_cx,
                "n_cd8_tex": n_tex,
                "frac_cd8_tex": n_tex / n_tcr if n_tcr else np.nan,
                "n_exp_tex_cells": n_exp_tex,
                "n_exp_tex_clones": n_exp_tex_clones,
                "exp_tex_per_k_cd8": (n_exp_tex_clones / n_cd8 * 1000.0) if n_cd8 else np.nan,
                "frac_tex_among_expanded": frac_tex_exp,
                "frac_tex_among_nonexpanded": frac_tex_non,
                "pass_depth": n_tcr >= min_tcr,
            }
        )
    return pd.DataFrame(rows)


def run_contrasts(end: pd.DataFrame, metrics: list[str]) -> pd.DataFrame:
    rows = []
    work = end[end["pass_depth"]].copy()
    for metric in metrics:
        a = work.loc[work["path_bin"] == "MPR-any", metric]
        b = work.loc[work["path_bin"] == "non-MPR", metric]
        rows.append(mw_row(a, b, metric, "MPR-any_vs_non-MPR"))
        for hist in ("LUAD", "LUSC"):
            h = work[work["cancer_type"] == hist]
            rows.append(
                mw_row(
                    h.loc[h["path_bin"] == "MPR-any", metric],
                    h.loc[h["path_bin"] == "non-MPR", metric],
                    metric,
                    f"{hist}_MPR-any_vs_non-MPR",
                )
            )
        rows.append(
            mw_row(
                work.loc[work["path_pcr"] == "pCR", metric],
                work.loc[work["path_pcr"] == "non-pCR", metric],
                metric,
                "pCR_vs_non-pCR",
            )
        )
        rows.append(
            mw_row(
                work.loc[work["path_mpr_only"] == "MPR-only", metric],
                work.loc[work["path_mpr_only"] == "non-MPR", metric],
                metric,
                "MPR-only_vs_non-MPR",
            )
        )
    return pd.DataFrame(rows)


def residual_corrs(end: pd.DataFrame, residual_cols: list[str], tcr_metrics: list[str]) -> pd.DataFrame:
    work = end[end["pass_depth"]].copy()
    rows = []
    for rx in residual_cols:
        if rx not in work.columns:
            continue
        for ty in tcr_metrics:
            rows.append(spearman_row(work[rx], work[ty], rx, ty))
            for hist in ("LUAD", "LUSC"):
                h = work[work["cancer_type"] == hist]
                r = spearman_row(h[rx], h[ty], rx, ty)
                r["stratum"] = hist
                rows.append(r)
    out = pd.DataFrame(rows)
    if "stratum" not in out.columns:
        out["stratum"] = "all"
    out["stratum"] = out["stratum"].fillna("all")
    return out


def make_figures(end: pd.DataFrame, outdir: Path) -> None:
    work = end[end["pass_depth"] & end["path_bin"].isin(["MPR-any", "non-MPR"])].copy()
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6))
    specs = [
        ("exp_cxcl13_per_k_cd8", "Expanded CXCL13+ clones / 1k CD8"),
        ("frac_cxcl13_among_expanded", "CXCL13+ fraction among expanded cells"),
        ("frac_exp_cells", "Expanded-cell fraction"),
    ]
    order = ["MPR-any", "non-MPR"]
    for ax, (col, title) in zip(axes, specs):
        data = [work.loc[work["path_bin"] == g, col].dropna().to_numpy() for g in order]
        ax.boxplot(data, tick_labels=order, widths=0.55)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel(col, fontsize=8)
        ax.tick_params(labelsize=8)
    fig.suptitle("GSE243013 patient-level TCR endpoints (depth-pass)", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "fig_gse243013_mpr_boxplots.png", dpi=140)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(4.4, 4.0))
    ok = work.dropna(subset=["TACSTD2_frac_pos", "exp_cxcl13_per_k_cd8"])
    colors = {"MPR-any": "#2a6f97", "non-MPR": "#c44536"}
    for g, sub in ok.groupby("path_bin"):
        ax.scatter(
            sub["TACSTD2_frac_pos"],
            sub["exp_cxcl13_per_k_cd8"],
            s=14,
            alpha=0.7,
            c=colors.get(g, "gray"),
            label=f"{g} n={len(sub)}",
        )
    ax.set_xlabel("Residual TACSTD2 fraction positive (immune MTX)")
    ax.set_ylabel("Expanded CXCL13+ clones / 1k CD8")
    ax.legend(fontsize=8)
    ax.set_title("Not a tumor-epithelial score", fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir / "fig_gse243013_residual_tacstd2_vs_cxcl13.png", dpi=140)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--tcr", type=Path, required=True)
    p.add_argument("--meta", type=Path, required=True)
    p.add_argument("--residual", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--min-tcr", type=int, default=50)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    tcr = pd.read_csv(args.tcr, low_memory=False)
    meta = pd.read_csv(args.meta, low_memory=False)
    residual = pd.read_csv(args.residual, sep="\t")

    sm = sample_meta(meta)
    end = endpoints(tcr, args.min_tcr)
    end = end.merge(sm, on="sampleID", how="left")
    keep_res = [
        "sampleID",
        "n_cells",
        "TACSTD2_frac_pos",
        "TACSTD2_mean_cpm",
        "CLDN4_frac_pos",
        "CLDN4_mean_cpm",
    ]
    end = end.merge(residual[keep_res], on="sampleID", how="left")
    write_tsv(end, args.out / "patient_endpoints.tsv")

    metrics = [
        "exp_cxcl13_per_k_cd8",
        "exp_cxcl13_per_k_tcr",
        "frac_cxcl13",
        "frac_cxcl13_that_are_exp",
        "frac_cxcl13_among_expanded",
        "frac_cxcl13_among_nonexpanded",
        "lor_cxcl13_exp_vs_non",
        "frac_exp_cells",
        "frac_exp_cells_ge3",
        "exp_tex_per_k_cd8",
        "frac_cd8_tex",
        "frac_tex_among_expanded",
        "n_exp_cxcl13_clones",
        "n_exp_tex_clones",
    ]
    stats = run_contrasts(end, metrics)
    write_tsv(stats, args.out / "stats_vs_mpr.tsv")

    combo = []
    work = end[end["pass_depth"]].copy()
    combo.append(
        paired_mw(
            work["frac_cxcl13_among_expanded"],
            work["frac_cxcl13_among_nonexpanded"],
            "frac_CXCL13_expanded_vs_nonexpanded",
        )
    )
    combo.append(
        paired_mw(
            work["frac_tex_among_expanded"],
            work["frac_tex_among_nonexpanded"],
            "frac_CD8Tex_expanded_vs_nonexpanded",
        )
    )
    write_tsv(pd.DataFrame(combo), args.out / "stats_combinatorial_paired.tsv")

    residual_metrics = [
        "exp_cxcl13_per_k_cd8",
        "frac_cxcl13",
        "frac_cxcl13_among_expanded",
        "frac_exp_cells",
        "exp_tex_per_k_cd8",
        "lor_cxcl13_exp_vs_non",
    ]
    cor = residual_corrs(
        end,
        ["TACSTD2_frac_pos", "CLDN4_frac_pos", "TACSTD2_mean_cpm", "CLDN4_mean_cpm"],
        residual_metrics,
    )
    write_tsv(cor, args.out / "stats_vs_residual_tacstd2_cldn4.tsv")

    # cluster census (descriptive, not a test)
    census = (
        tcr.groupby("sub_cell_type")
        .agg(n_cells=("cellID", "size"), n_samples=("sampleID", "nunique"))
        .reset_index()
    )
    write_tsv(census, args.out / "cluster_census.tsv")

    n_tcr_samples = int(end.shape[0])
    n_pass = int(end["pass_depth"].sum())
    n_mpr = int(((end["pass_depth"]) & (end["path_bin"] == "MPR-any")).sum())
    n_non = int(((end["pass_depth"]) & (end["path_bin"] == "non-MPR")).sum())
    spec = {
        "dataset": "GSE243013",
        "pmid": 40147443,
        "compartment": "CD45+ immune only; no malignant epithelium",
        "unit": "sample/patient",
        "n_tcr_cells": int(len(tcr)),
        "n_tcr_samples": n_tcr_samples,
        "n_pass_depth": n_pass,
        "min_tcr_cells": args.min_tcr,
        "n_MPR_any_pass": n_mpr,
        "n_nonMPR_pass": n_non,
        "cxcl13_clusters": list(CXCL13_CLUSTERS),
        "cd8_tex_clusters": list(CD8_TEX_CLUSTERS),
        "note_no_CD8T_Tex_CXCL13": True,
        "clone_key": "TRA_cdr3|TRB_cdr3 amino-acid, within sample",
        "expansion_primary": "clone_size >= 2",
        "e2_tumor_tacstd2": "not estimable: no paired tumor-epithelial score",
        "residual_tacstd2_source": str(args.residual),
        "residual_is_not_tumor_score": True,
        "forbidden": "Wilcoxon on cells labeled MPR; using residual TACSTD2 as E2",
    }
    write_json(spec, args.out / "run_info.json")
    (args.out / "e2_status.txt").write_text(
        "not estimable: GSE243013 is immune-only. Residual TACSTD2/CLDN4 "
        "detection in the CD45+ MTX is not a tumor-epithelial score.\n"
    )
    make_figures(end, args.out)
    print(
        f"GSE243013 wrote {args.out}  samples={n_tcr_samples} "
        f"pass={n_pass} MPR-any={n_mpr} non-MPR={n_non}"
    )


if __name__ == "__main__":
    main()
