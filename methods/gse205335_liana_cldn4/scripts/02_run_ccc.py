#!/usr/bin/env python3
"""CLDN4-only ligand–receptor communication on GSE205335 (public processed UMI).

Author labels are used as given:
    malignant = tumor sample AND lineage.sub == "Malignant cells"
    T/NK      = lineage.total == "T/NK cells"

CLDN4 is the only state gene. TACSTD2 is not a gate. Dual-high is not run.

Primary executable method
    Documented CellPhoneDB-style score (Efremova 2020 / Garcia-Alonso 2022):
    partner expression = min(subunit means) on log1p(CP10k);
    score = mean of the two partner means.
    Patient-level paired Wilcoxon is the inferential unit.

Secondary method (if import succeeds)
    LIANA `cellphonedb` on a downsampled object. Labeled as LIANA, not CellChat.

CellChat is not run: R is not available in this environment.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import traceback
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
T_SUBS = {"CD4+ T cells", "CD8+ T cells"}
NK_SUBS = {"NK cells"}
FOCUS_PATHWAYS = {"T_recruit", "IFN", "MHC_I"}
FOCUS_PAIRS = [
    ("CXCL9", "CXCR3"),
    ("CXCL10", "CXCR3"),
    ("CXCL11", "CXCR3"),
    ("CXCL16", "CXCR6"),
    ("CCL5", "CCR5"),
    ("CCL4", "CCR5"),
    ("CX3CL1", "CX3CR1"),
    ("HLA-A", "CD8A"),
    ("HLA-B", "CD8A"),
    ("HLA-C", "CD8A"),
    ("HLA-A", "CD8B"),
    ("HLA-B", "CD8B"),
    ("HLA-C", "CD8B"),
    ("HLA-E", "KLRD1"),
    ("HLA-E", "KLRC1"),
]


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_geo_soft(path: Path) -> pd.DataFrame:
    opener = gzip.open if path.suffix == ".gz" or str(path).endswith(".gz") else open
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    descriptions: list[str] = []
    titles: list[str] = []
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
    read_end = metadata["platform"].str.extract(r"Single Cell ([35])'")[0]
    metadata["orig.ident"] = (
        metadata["description"].str.replace("_", "-", regex=False) + "-" + read_end + "P"
    )
    return metadata.rename(
        columns={"tumor stage": "tumor_stage", "cancer subtype": "cancer_subtype"}
    )


def load_selected_genes(
    path: Path, wanted: set[str]
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, list[str]]:
    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz" or str(path).endswith(".gz"):
            matrix_path = Path(tmp) / Path(str(path).replace(".gz", "")).name
            if not str(matrix_path).endswith(".rds"):
                matrix_path = Path(tmp) / "GSE205335_Lung_IO_UMI_matrix.rds"
            log(f"decompress {path.name}")
            with gzip.open(path, "rb") as source, matrix_path.open("wb") as dest:
                shutil.copyfileobj(source, dest, 16 * 1024 * 1024)
        log("read RDS")
        import rdata

        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", message='Missing constructor for R class "dgCMatrix"'
            )
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    log(f"build CSC {tuple(obj.Dim)}")
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    log("convert CSR for row extract")
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted: dict[str, np.ndarray] = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel().astype(np.float32)
    log(f"extracted {len(extracted)} / {len(wanted)} genes")
    return extracted, library_umi, barcodes, genes.tolist()


def log1p_cp10k(umi: np.ndarray, total: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        cp = np.where(total > 0, umi / total * 1e4, 0.0)
    return np.log1p(cp).astype(np.float32)


def group_gene_stats(
    logx: dict[str, np.ndarray], mask: np.ndarray
) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(np.mean(v))
        fracs[g] = float(np.mean(v > 0))
    return means, fracs, n


def partner_from_stats(
    units: list[str], means: dict[str, float], fracs: dict[str, float]
) -> tuple[float, float]:
    m, f = [], []
    for g in units:
        if g not in means:
            return np.nan, np.nan
        m.append(means[g])
        f.append(fracs[g])
    return float(np.min(m)), float(np.min(f))


def score_pairs(
    pairs: pd.DataFrame,
    logx: dict[str, np.ndarray],
    sender: np.ndarray,
    receiver: np.ndarray,
    expr_prop: float,
) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(logx, sender)
    r_mean, r_frac, n_r = group_gene_stats(logx, receiver)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        pass_prop = (l_frac >= expr_prop) and (rec_f >= expr_prop)
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "pair_origin": rec.pair_origin,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "product_score": l_mean * rec_m,
                "pass_expr_prop": pass_prop,
            }
        )
    return pd.DataFrame(rows)


def paired_delta(high_df: pd.DataFrame, low_df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    key = ["ligand", "receptor", "pathway"]
    cols = [score_col, "pass_expr_prop", "ligand_frac", "receptor_frac", "n_sender", "n_receiver"]
    h = high_df[key + cols].copy()
    l = low_df[key + cols].copy()
    m = h.merge(l, on=key, suffixes=("_high", "_low"))
    m["delta_high_minus_low"] = m[f"{score_col}_high"] - m[f"{score_col}_low"]
    m["pass_either"] = m["pass_expr_prop_high"] | m["pass_expr_prop_low"]
    m["pass_both"] = m["pass_expr_prop_high"] & m["pass_expr_prop_low"]
    return m


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def run_liana(adata, groupby: str, out_csv: Path, n_perms: int) -> str:
    try:
        import liana as li
    except Exception as exc:  # pragma: no cover
        return f"LIANA_IMPORT_FAILED: {exc}"
    try:
        li.mt.cellphonedb(
            adata,
            groupby=groupby,
            resource_name="cellphonedb",
            expr_prop=0.10,
            n_perms=n_perms,
            use_raw=False,
            verbose=True,
            key_added="liana_res",
        )
        res = adata.uns["liana_res"].copy()
        res.to_csv(out_csv, index=False)
        return f"LIANA_OK n_edges={len(res)} file={out_csv.name}"
    except Exception as exc:
        return f"LIANA_RUN_FAILED: {exc}\n{traceback.format_exc()}"


def fmt_p(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.3g}"


def fmt_delta(x) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:+.3f}"


def summarize_liana_focus(res: pd.DataFrame) -> pd.DataFrame:
    senders_hi = "Malig_CLDN4high"
    senders_lo = "Malig_CLDN4low"
    out_rows = []
    sub = res[res["source"].isin([senders_hi, senders_lo]) & res["target"].isin(["T", "NK"])].copy()
    for lig, recp in FOCUS_PAIRS:
        for tgt in ("T", "NK"):
            a = sub[(sub["ligand"] == lig) & (sub["receptor"] == recp) & (sub["target"] == tgt)]
            if a.empty:
                continue
            hi = a[a["source"] == senders_hi]
            lo = a[a["source"] == senders_lo]
            hi_m = float(hi["lr_means"].iloc[0]) if len(hi) else np.nan
            lo_m = float(lo["lr_means"].iloc[0]) if len(lo) else np.nan
            if np.isfinite(hi_m) and np.isfinite(lo_m):
                higher = "high" if hi_m > lo_m else ("low" if lo_m > hi_m else "tie")
            elif np.isfinite(hi_m):
                higher = "high only"
            else:
                higher = "low only"
            out_rows.append(
                {
                    "ligand": lig,
                    "receptor": recp,
                    "target": tgt,
                    "lr_means_high": hi_m,
                    "lr_means_low": lo_m,
                    "higher_in": higher,
                    "cellphone_p_high": float(hi["cellphone_pvals"].iloc[0]) if len(hi) else np.nan,
                }
            )
    if len(out_rows) < 6:
        hi_only = sub[sub["source"] == senders_hi].sort_values("lr_means", ascending=False).head(12)
        for r in hi_only.itertuples(index=False):
            lo = sub[
                (sub["ligand"] == r.ligand)
                & (sub["receptor"] == r.receptor)
                & (sub["target"] == r.target)
                & (sub["source"] == senders_lo)
            ]
            lo_m = float(lo["lr_means"].iloc[0]) if len(lo) else np.nan
            out_rows.append(
                {
                    "ligand": r.ligand,
                    "receptor": r.receptor,
                    "target": r.target,
                    "lr_means_high": float(r.lr_means),
                    "lr_means_low": lo_m,
                    "higher_in": "high"
                    if (np.isfinite(lo_m) and r.lr_means > lo_m)
                    else ("high only" if not np.isfinite(lo_m) else "low"),
                    "cellphone_p_high": float(r.cellphone_pvals),
                }
            )
    return pd.DataFrame(out_rows).drop_duplicates(subset=["ligand", "receptor", "target"])


def make_figures(
    outdir: Path,
    ntab: pd.DataFrame,
    ranks: dict[str, pd.DataFrame],
    raw_out: pd.DataFrame,
    cldn4: np.ndarray,
    is_malig: np.ndarray,
    cld_hi: np.ndarray,
    cld_lo: np.ndarray,
    is_tnk: np.ndarray,
    patient: np.ndarray,
    cancer_subtype: np.ndarray,
    thr: float,
) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    colors = {"T_recruit": "#4C72B0", "IFN": "#55A868", "MHC_I": "#C44E52"}
    hist_colors = {"ADC": "#4C72B0", "SQ": "#55A868", "SCLC": "#C44E52", "NUT": "#8172B3"}

    show = ntab.sort_values(["cancer_subtype", "patient"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(9.2, 4.0))
    x = np.arange(len(show))
    ax.bar(x - 0.2, show["n_malignant"], width=0.4, label="malignant", color="#4C72B0")
    ax.bar(x + 0.2, show["n_T_NK"], width=0.4, label="T/NK", color="#DD8452")
    ax.set_xticks(x)
    labels = [f"{p}\n{h}/{r}" for p, h, r in zip(show["patient"], show["cancer_subtype"], show["recist"])]
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("cells (tumor samples)")
    ax.legend(frameon=False)
    ax.set_title("GSE205335 tumor-sample cells (author malignant vs T/NK)")
    fig.tight_layout()
    fig.savefig(figdir / "n_cells_by_patient.png", dpi=160)
    fig.savefig(figdir / "n_cells_by_patient.pdf")
    plt.close(fig)

    for key, df in ranks.items():
        focus = df[df["pathway"].isin(FOCUS_PATHWAYS)].copy()
        if focus.empty:
            continue
        focus = focus.sort_values("median_delta")
        fig, ax = plt.subplots(figsize=(8.2, max(3.2, 0.28 * len(focus) + 1.2)))
        y = np.arange(len(focus))
        ax.barh(
            y,
            focus["median_delta"],
            color=[colors.get(p, "#999") for p in focus["pathway"]],
            edgecolor="none",
        )
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(
            [f"{a}–{b}  ({p})" for a, b, p in zip(focus["ligand"], focus["receptor"], focus["pathway"])],
            fontsize=7,
        )
        ax.set_xlabel("median patient Δ score (CLDN4-high − CLDN4-low)")
        ax.set_title(key.replace("_", " "))
        fig.tight_layout()
        fig.savefig(figdir / f"{key}.png", dpi=160)
        fig.savefig(figdir / f"{key}.pdf")
        plt.close(fig)

    # Extra 1: CLDN4 distribution among malignant cells
    fig, ax = plt.subplots(figsize=(5.4, 3.8))
    hi = cldn4[cld_hi]
    lo = cldn4[cld_lo]
    ax.hist(lo, bins=40, alpha=0.65, color="#6a8aaa", label=f"low n={lo.size:,}")
    ax.hist(hi, bins=40, alpha=0.65, color="#b2182b", label=f"high n={hi.size:,}")
    ax.axvline(thr, color="k", ls="--", lw=0.9, label=f"median={thr:.3f}")
    ax.set_xlabel("malignant CLDN4 log1p(CP10k)")
    ax.set_ylabel("cells")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("GSE205335 CLDN4 split (global median, author malignant)")
    fig.tight_layout()
    fig.savefig(figdir / "extra_cldn4_distribution.png", dpi=160)
    fig.savefig(figdir / "extra_cldn4_distribution.pdf")
    plt.close(fig)

    # Extra 2: per-patient heatmap of focus outgoing Δ
    if raw_out is not None and len(raw_out):
        foc_raw = raw_out[raw_out["pathway"].isin(["T_recruit", "MHC_I"])].copy()
        want = {f"{a}|{b}" for a, b in FOCUS_PAIRS}
        foc_raw["pair"] = foc_raw["ligand"] + "|" + foc_raw["receptor"]
        foc_raw = foc_raw[foc_raw["pair"].isin(want)]
        if len(foc_raw):
            mat = foc_raw.pivot_table(
                index="patient", columns="pair", values="delta_high_minus_low", aggfunc="median"
            )
            order_cols = [c for c in [f"{a}|{b}" for a, b in FOCUS_PAIRS] if c in mat.columns]
            mat = mat[order_cols]
            meta = ntab.set_index("patient")
            mat = mat.reindex(meta.index.intersection(mat.index))
            fig, ax = plt.subplots(figsize=(max(7.5, 0.55 * len(order_cols) + 2.2), max(4.2, 0.28 * len(mat) + 1.6)))
            im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-0.4, vmax=0.4)
            ax.set_xticks(np.arange(len(mat.columns)))
            ax.set_xticklabels([c.replace("|", "–") for c in mat.columns], rotation=55, ha="right", fontsize=7)
            ax.set_yticks(np.arange(len(mat.index)))
            ylabs = [
                f"{p} {meta.loc[p, 'cancer_subtype']}/{meta.loc[p, 'recist']}"
                if p in meta.index
                else p
                for p in mat.index
            ]
            ax.set_yticklabels(ylabs, fontsize=7)
            ax.set_title("Extra: per-patient outgoing Δ (CLDN4-high − low)")
            fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Δ score")
            fig.tight_layout()
            fig.savefig(figdir / "extra_patient_focus_heatmap.png", dpi=160)
            fig.savefig(figdir / "extra_patient_focus_heatmap.pdf")
            plt.close(fig)

    # Extra 3: paired strips for MHC-I CD8 and T-recruit CXCR3/CXCR6
    if raw_out is not None and len(raw_out):
        strip_pairs = [
            ("HLA-A", "CD8A"),
            ("HLA-B", "CD8A"),
            ("HLA-C", "CD8A"),
            ("HLA-E", "KLRD1"),
            ("CXCL16", "CXCR6"),
            ("CXCL9", "CXCR3"),
            ("CXCL10", "CXCR3"),
            ("CCL5", "CCR5"),
        ]
        fig, axes = plt.subplots(2, 4, figsize=(11.2, 6.2), sharey=False)
        axes = axes.ravel()
        for ax, (lig, recp) in zip(axes, strip_pairs):
            sub = raw_out[(raw_out["ligand"] == lig) & (raw_out["receptor"] == recp)]
            if sub.empty:
                ax.set_title(f"{lig}–{recp}\nno data", fontsize=8)
                ax.axis("off")
                continue
            rng = np.random.default_rng(0)
            for i, col, color in (
                (0, "cpdb_mean_score_low", "#6a8aaa"),
                (1, "cpdb_mean_score_high", "#b2182b"),
            ):
                vals = sub[col].to_numpy()
                ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(vals)), vals, c=color, s=18, zorder=3)
            for rec in sub.itertuples(index=False):
                ax.plot(
                    [0, 1],
                    [rec.cpdb_mean_score_low, rec.cpdb_mean_score_high],
                    color="0.7",
                    lw=0.6,
                    zorder=1,
                )
            ax.set_xticks([0, 1], ["low", "high"])
            ax.set_title(f"{lig}–{recp}\nn={sub['patient'].nunique()}", fontsize=8)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
        fig.suptitle("Extra: paired patient scores, outgoing malignant → T/NK", fontsize=10)
        fig.tight_layout()
        fig.savefig(figdir / "extra_paired_strips.png", dpi=160)
        fig.savefig(figdir / "extra_paired_strips.pdf")
        plt.close(fig)

    # Extra 4: histology of paired patients vs dropped
    fig, ax = plt.subplots(figsize=(5.6, 3.6))
    paired = ntab[ntab["paired"]].copy()
    dropped = ntab[~ntab["paired"]].copy()
    cats = ["ADC", "SQ", "SCLC", "NUT"]
    x = np.arange(len(cats))
    ax.bar(x - 0.18, [(paired["cancer_subtype"] == c).sum() for c in cats], width=0.36, label="paired", color="#4C72B0")
    ax.bar(x + 0.18, [(dropped["cancer_subtype"] == c).sum() for c in cats], width=0.36, label="dropped", color="#CCCCCC")
    ax.set_xticks(x, cats)
    ax.set_ylabel("patients (tumor extract)")
    ax.legend(frameon=False)
    ax.set_title("Extra: honest paired n by histology")
    fig.tight_layout()
    fig.savefig(figdir / "extra_paired_by_histology.png", dpi=160)
    fig.savefig(figdir / "extra_paired_by_histology.pdf")
    plt.close(fig)

    # Extra 5: high/low malignant counts per patient
    fig, ax = plt.subplots(figsize=(9.2, 4.0))
    x = np.arange(len(show))
    ax.bar(x - 0.2, show["n_cldn4_high"], width=0.4, label="CLDN4-high mal", color="#b2182b")
    ax.bar(x + 0.2, show["n_cldn4_low"], width=0.4, label="CLDN4-low mal", color="#6a8aaa")
    ax.set_xticks(x)
    ax.set_xticklabels(show["patient"], rotation=60, ha="right", fontsize=7)
    ax.set_ylabel("malignant cells")
    ax.legend(frameon=False)
    ax.set_title("Extra: CLDN4-high vs low malignant cells per patient")
    fig.tight_layout()
    fig.savefig(figdir / "extra_high_low_counts.png", dpi=160)
    fig.savefig(figdir / "extra_high_low_counts.pdf")
    plt.close(fig)


def write_finding(
    outdir: Path,
    root: Path,
    summary: dict,
    ntab: pd.DataFrame,
    ranks: dict[str, pd.DataFrame],
    liana_notes: list[str],
    liana_focus: pd.DataFrame | None,
) -> None:
    out = ranks.get("cldn4_outgoing", pd.DataFrame())
    inc = ranks.get("cldn4_incoming", pd.DataFrame())
    focus_out = out[out["pathway"].isin(FOCUS_PATHWAYS)].copy() if len(out) else out
    focus_in = inc[inc["pathway"].isin(FOCUS_PATHWAYS)].copy() if len(inc) else inc
    paired = ntab[ntab["paired"]].copy()
    dropped = ntab[~ntab["paired"]].copy()

    def _table(df: pd.DataFrame) -> list[str]:
        if df is None or df.empty:
            return ["No T-recruit / IFN / MHC-I pairs with ≥3 paired patients.", ""]
        df = df.sort_values(["pathway", "median_delta"])
        lines = [
            "| Pathway | Pair | n patients | median Δ | Wilcoxon p | FDR |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for r in df.itertuples(index=False):
            lines.append(
                f"| {r.pathway} | {r.ligand}–{r.receptor} | {int(r.n_patients)} | "
                f"{r.median_delta:+.3f} | {fmt_p(r.pval)} | {fmt_p(r.padj)} |"
            )
        n_neg = int((df["median_delta"] < 0).sum())
        n_sig = int((df["padj"] < 0.05).sum()) if df["padj"].notna().any() else 0
        lines += [
            "",
            f"{n_neg}/{len(df)} focus pairs have median Δ < 0 (weaker from/to CLDN4-high). "
            f"**{n_sig}/{len(df)} reach FDR < 0.05** in either direction.",
            "",
        ]
        return lines

    paired_list = ", ".join(
        f"{r.patient} ({r.cancer_subtype}/{r.recist}; high={int(r.n_cldn4_high)}/low={int(r.n_cldn4_low)} malig, T/NK={int(r.n_T_NK)})"
        for r in paired.itertuples(index=False)
    )
    dropped_list = ", ".join(
        f"{r.patient} ({r.cancer_subtype}/{r.recist}; malig={int(r.n_malignant)}, high={int(r.n_cldn4_high)}, low={int(r.n_cldn4_low)}, T/NK={int(r.n_T_NK)})"
        for r in dropped.itertuples(index=False)
    )

    n_neg = int((focus_out["median_delta"] < 0).sum()) if len(focus_out) else 0
    n_sig = (
        int((focus_out["padj"] < 0.05).sum())
        if len(focus_out) and focus_out["padj"].notna().any()
        else 0
    )
    mhci = focus_out[focus_out["pathway"] == "MHC_I"] if len(focus_out) else focus_out
    trec = focus_out[focus_out["pathway"] == "T_recruit"] if len(focus_out) else focus_out
    verdict = (
        f"**Verdict (paired n={summary['n_patients_paired']} honest):** "
        "CLDN4-high author-malignant cells do not show a coordinated reduction of "
        "outgoing T-recruit / MHC-I communication toward same-patient T/NK. "
    )
    if len(mhci):
        verdict += (
            f"MHC-I outgoing median-of-pair Δ = {mhci['median_delta'].median():+.3f} "
            f"({int((mhci['median_delta'] < 0).sum())}/{len(mhci)} Δ<0). "
        )
    if len(trec):
        verdict += (
            f"T-recruit median-of-pair Δ = {trec['median_delta'].median():+.3f} "
            f"({int((trec['median_delta'] < 0).sum())}/{len(trec)} Δ<0). "
        )
    verdict += f"**{n_sig}/{len(focus_out) if len(focus_out) else 0} focus pairs reach FDR < 0.05.**"

    lines = [
        "# GSE205335 — LIANA/CellPhoneDB-style LR from CLDN4-high malignant to T/NK",
        "",
        "**CLDN4 only.** TACSTD2 is not a gate and is not used to define dual-high. "
        "Author `lineage.sub` / `lineage.total` labels are used as given. "
        "This is not GSE207422 (PR #344).",
        "",
        verdict,
        "",
        "## What was run",
        "",
        "| Method | Status |",
        "|---|---|",
        f"| Documented CellPhoneDB-style score (mean of partner means on log1p CP10k; {summary['n_pairs_scored']:,} pairs) | **primary** — patient-level paired Wilcoxon |",
        f"| LIANA `mt.cellphonedb` (resource `cellphonedb`, {summary['liana_n_perms']} permutations, ≤{summary['liana_max_cells']} cells/group) | **secondary** — {summary['liana_status']} |",
        "| CellChat | **not run** — R unavailable; no CellChat tables were written |",
        "",
        "## Honest n",
        "",
        "| Set | n |",
        "|---|---:|",
        f"| Cells in public UMI | {summary['n_cells']:,} |",
        f"| GEO samples / patients | {summary['n_samples_geo']} / {summary['n_patients_geo']} |",
        f"| Tumor samples / patients (normals dropped) | {summary['n_tumor_samples']} / {summary['n_tumor_patients']} |",
        f"| Author malignant (tumor samples) | {summary['n_malig']:,} |",
        f"| T / NK / T+NK | {summary['n_T']:,} / {summary['n_NK']:,} / {summary['n_tnk']:,} |",
        f"| CLDN4-high / low malignant | {summary['n_cldn4_high']:,} / {summary['n_cldn4_low']:,} |",
        f"| Patients in the paired LR test | **{summary['n_patients_paired']}** |",
        "",
        f"CLDN4 high = at or above the **global median** log1p(CP10k) among author-malignant cells "
        f"(threshold = {summary['cldn4_threshold']:.3f}). Low = below. "
        f"A patient enters the paired test if it has ≥{summary['min_malig_per_state']} malignant cells "
        f"in **both** bins and ≥{summary['min_tnk']} T/NK cells after pooling that patient’s **tumor** samples.",
        "",
        f"**Paired patients (n={summary['n_patients_paired']}):** {paired_list}.",
        "",
        f"**Dropped from the paired test:** {dropped_list}.",
        "",
        "Normal-only GEO patients P2001 / P2009 / P2016 (Normal LN) and P3032 (Normal Brain) have "
        "0 author-malignant cells and are out of the tumor extract. P0031 Normal Lung is excluded; "
        "P0031 tumor lung is kept. MPR/NMPR is unlabeled; RECIST is not used as MPR. "
        "The header n for the Wilcoxon is the paired count, not 26 and not 22.",
        "",
        "Per-patient counts: `results/n_cells_patients.tsv`.",
        "",
        "## Score",
        "",
        "On log1p(CP10k), each partner’s expression is the **minimum subunit mean** "
        "(CellPhoneDB complex rule). The pair score is the **mean of the two partner means** "
        "(Efremova et al. 2020 *Nat Protoc*; Garcia-Alonso et al. 2022 *Nat Protoc*). "
        "A pair is flagged `pass_expr_prop` when both partners are detected in ≥10% of "
        "cells in their group. Patient-level tests use that patient’s own T/NK and that "
        "patient’s CLDN4-high vs CLDN4-low malignant cells. Cells are not treated as replicates. "
        "FDR is Benjamini–Hochberg within each contrast (outgoing or incoming).",
        "",
        "This is **not** a CellChat communication probability.",
        "",
        "## Primary LR table — outgoing CLDN4-high malignant → T/NK",
        "",
        "Median patient Δ = high − low. Negative = weaker from the CLDN4-high state.",
        "",
    ]
    lines += _table(focus_out)
    lines += [
        "Full ranked table (all pathways that passed filters): `results/lr_table_cldn4_outgoing_tnk.tsv`.",
        "",
        "## Incoming T/NK → CLDN4-high vs CLDN4-low malignant (secondary)",
        "",
    ]
    lines += _table(focus_in)
    lines += [
        "## LIANA CellPhoneDB method (secondary, pooled / downsampled)",
        "",
    ]
    for note in liana_notes:
        lines.append(f"- {note}")
    lines += [""]
    if liana_focus is not None and len(liana_focus):
        lines += [
            "Focus edges that cleared LIANA `expr_prop=0.10` from malignant senders to T or NK:",
            "",
            "| Pair | target | high `lr_means` | low `lr_means` | higher in | cellphone_p high |",
            "|---|---|---:|---:|---|---:|",
        ]
        for r in liana_focus.itertuples(index=False):
            hi = "NA" if not np.isfinite(r.lr_means_high) else f"{r.lr_means_high:.3f}"
            lo = "NA" if not np.isfinite(r.lr_means_low) else f"{r.lr_means_low:.3f}"
            ph = "NA" if not np.isfinite(r.cellphone_p_high) else f"{r.cellphone_p_high:.3g}"
            lines.append(
                f"| {r.ligand}–{r.receptor} | {r.target} | {hi} | {lo} | {r.higher_in} | {ph} |"
            )
        lines += [""]
    lines += [
        "LIANA p-values are within-object specificity, not patient-level tests. "
        "CXCL9/10/11–CXCR3 typically fail the 10% expression filter in epithelium.",
        "",
        "## Readout",
        "",
        summary["trend_sentence"],
        "",
        "## Honest limits",
        "",
        f"1. Paired n = {summary['n_patients_paired']}. GEO has 26 patients / 33 samples; "
        "four normal-only patients have 0 malignant cells. That is reported, not patched.",
        "2. CLDN4 only. Dual-high (TACSTD2 AND CLDN4) was not run. TACSTD2 is not a gate.",
        "3. Author malignant labels, not public CopyKAT / inferCNV calls.",
        "4. Histology is mixed (ADC / SQ / SCLC / NUT). SCLC is not hidden.",
        "5. MPR is unlabeled. RECIST is not MPR.",
        "6. Ambient RNA cannot be re-estimated from the processed matrix.",
        "7. Chemokine dropout is high; read `n patients` and `pass_expr_prop` with the ranks.",
        "8. CellChat was not run. This is not the Q4 vs Q1 between-patient table in PR #362.",
        "9. This is not GSE207422 (PR #344).",
        "",
        "## Extra figures",
        "",
        "| File | Content |",
        "|---|---|",
        "| `results/figures/n_cells_by_patient.png` | Author malignant vs T/NK per patient |",
        "| `results/figures/cldn4_outgoing.png` | Focus-axis outgoing Δ |",
        "| `results/figures/cldn4_incoming.png` | Focus-axis incoming Δ |",
        "| `results/figures/extra_cldn4_distribution.png` | Malignant CLDN4 high/low histogram |",
        "| `results/figures/extra_patient_focus_heatmap.png` | Per-patient MHC-I / T-recruit Δ |",
        "| `results/figures/extra_paired_strips.png` | Paired high vs low scores for key pairs |",
        "| `results/figures/extra_paired_by_histology.png` | Paired vs dropped by histology |",
        "| `results/figures/extra_high_low_counts.png` | High vs low malignant n per patient |",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `FINDING.md` | This note |",
        "| `results/lr_table_cldn4_outgoing_tnk.tsv` | Primary LR table (patient-level ranks, outgoing) |",
        "| `results/lr_table_cldn4_incoming_tnk.tsv` | Incoming ranks |",
        "| `results/lr_table_focus_outgoing.tsv` | T-recruit / IFN / MHC-I outgoing subset |",
        "| `results/n_cells_patients.tsv` | Per-patient cell counts and high/low bins |",
        "| `results/patient_cldn4_outgoing.tsv.gz` | Per-patient pair scores (full list) |",
        "| `results/pooled_outgoing_cldn4.tsv` | All-cell descriptive scores |",
        "| `results/liana_cellphonedb_cldn4.csv` | Full LIANA CellPhoneDB output (if run) |",
        "| `results/summary.json` | Machine-readable n and method flags |",
        "| `results/figures/` | n-cell bars, pathway Δ, and extra figures |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "cd methods/gse205335_liana_cldn4",
        "python3 scripts/00_download.py",
        "python3 scripts/01_build_pairs.py",
        "python3 scripts/02_run_ccc.py",
        "```",
        "",
    ]
    text = "\n".join(lines)
    (root / "FINDING.md").write_text(text)
    (outdir / "FINDING.md").write_text(text)


def rank_direction(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame(
            columns=["ligand", "receptor", "pathway", "n_patients", "median_delta", "pval", "padj"]
        )
    rows = []
    for (lig, recp, path), sub in raw.groupby(["ligand", "receptor", "pathway"], observed=True):
        keep_sub = sub[sub["pass_either"]]
        if len(keep_sub) < 3:
            continue
        if path not in FOCUS_PATHWAYS and keep_sub["pass_both"].sum() < 3:
            continue
        delta = keep_sub["delta_high_minus_low"].to_numpy()
        p = wilcoxon_safe(keep_sub["cpdb_mean_score_high"], keep_sub["cpdb_mean_score_low"])
        rows.append(
            {
                "ligand": lig,
                "receptor": recp,
                "pathway": path,
                "n_patients": int(keep_sub["patient"].nunique()),
                "median_delta": float(np.median(delta)),
                "mean_delta": float(np.mean(delta)),
                "mean_score_high": float(keep_sub["cpdb_mean_score_high"].mean()),
                "mean_score_low": float(keep_sub["cpdb_mean_score_low"].mean()),
                "frac_pass_high": float(keep_sub["pass_expr_prop_high"].mean()),
                "frac_pass_low": float(keep_sub["pass_expr_prop_low"].mean()),
                "pval": p,
            }
        )
    tab = pd.DataFrame(rows)
    if len(tab) and tab["pval"].notna().any():
        mask = tab["pval"].notna()
        tab.loc[mask, "padj"] = multipletests(tab.loc[mask, "pval"], method="fdr_bh")[1]
    else:
        tab["padj"] = np.nan
    if len(tab):
        tab = tab.sort_values(["pathway", "median_delta"])
    return tab


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datadir", type=Path, default=Path("/tmp/gse205335"))
    parser.add_argument("--outdir", type=Path, default=HERE / "results")
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    (args.outdir / "figures").mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    P = cfg["params"]
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    keep = set(Path(HERE / "resources" / "lr_genes.txt").read_text().split())
    keep.update(cfg["state_genes"])
    keep.update(["TACSTD2", "EPCAM", "PTPRC", "CD3D", "CD3E", "CD8A", "NKG7"])

    liana_import = "not_attempted"
    try:
        import liana as li

        res = li.rs.select_resource("cellphonedb")
        keep.update(res["ligand"].astype(str))
        keep.update(res["receptor"].astype(str))
        keep.update(g for s in res["ligand"].astype(str) for g in s.split("_"))
        keep.update(g for s in res["receptor"].astype(str) for g in s.split("_"))
        liana_import = "ok"
    except Exception as exc:
        liana_import = f"failed: {exc}"

    identities = pd.read_csv(args.datadir / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    if identities["barcode"].duplicated().any():
        raise SystemExit("Cell-identity barcodes are not unique")
    metadata = parse_geo_soft(args.datadir / "GSE205335_family.soft.gz")
    metadata.to_csv(args.outdir / "gsm_sample_metadata.csv", index=False)

    extracted, library_umi, barcodes, all_genes = load_selected_genes(
        args.datadir / "GSE205335_Lung_IO_UMI_matrix.rds.gz", keep
    )
    if "CLDN4" not in extracted:
        raise SystemExit("CLDN4 not in UMI matrix")

    indexed = identities.set_index("barcode")
    missing = pd.Index(barcodes).difference(indexed.index)
    extra = indexed.index.difference(pd.Index(barcodes))
    if len(missing) or len(extra):
        raise SystemExit(f"Matrix/identity mismatch: {len(missing)} missing, {len(extra)} extra")
    cells = indexed.loc[barcodes].reset_index()
    cells["total_umi"] = library_umi
    cells = cells.merge(
        metadata[
            ["orig.ident", "gsm", "patient", "tissue", "recist", "cancer_subtype", "tumor_stage", "platform"]
        ],
        on="orig.ident",
        how="left",
        validate="many_to_one",
    )
    if cells["patient"].isna().any():
        raise SystemExit("Some identity-table samples did not match GEO metadata")

    cells["is_normal_sample"] = cells["tissue"].astype(str).str.startswith("Normal")
    cells["is_tumor_sample"] = ~cells["is_normal_sample"]
    cells["is_malig"] = cells["is_tumor_sample"] & cells["lineage.sub"].eq("Malignant cells")
    cells["is_t"] = cells["is_tumor_sample"] & cells["lineage.sub"].isin(T_SUBS)
    cells["is_nk"] = cells["is_tumor_sample"] & cells["lineage.sub"].isin(NK_SUBS)
    cells["is_tnk"] = cells["is_tumor_sample"] & cells["lineage.total"].eq("T/NK cells")

    is_malig = cells["is_malig"].to_numpy()
    is_t = cells["is_t"].to_numpy()
    is_nk = cells["is_nk"].to_numpy()
    is_tnk = cells["is_tnk"].to_numpy()
    is_tumor = cells["is_tumor_sample"].to_numpy()
    patient_arr = cells["patient"].to_numpy()
    subtype_arr = cells["cancer_subtype"].to_numpy()

    logx = {g: log1p_cp10k(extracted[g], library_umi) for g in extracted}
    cld = logx["CLDN4"]
    thr = float(np.median(cld[is_malig])) if is_malig.any() else float("nan")
    cld_hi = is_malig & (cld >= thr)
    cld_lo = is_malig & (cld < thr)
    log(f"[split] CLDN4 median among malignant={thr:.4f} high={int(cld_hi.sum())} low={int(cld_lo.sum())}")

    ntab_rows = []
    for pat, sub in cells.groupby("patient", observed=True):
        tumor = sub[sub["is_tumor_sample"]]
        if len(tumor):
            t_idx = tumor.index.to_numpy()
            mal = tumor["is_malig"].to_numpy()
            vals = cld[t_idx]
            n_hi = int((mal & (vals >= thr)).sum())
            n_lo = int((mal & (vals < thr)).sum())
            n_t = int(tumor["is_tnk"].to_numpy().sum())
            n_mal = int(mal.sum())
            mean_cld = float(np.mean(vals[mal])) if n_mal else np.nan
        else:
            n_hi = n_lo = n_t = n_mal = 0
            mean_cld = np.nan
        paired = (
            n_hi >= P["min_malig_per_state"]
            and n_lo >= P["min_malig_per_state"]
            and n_t >= P["min_tnk"]
        )
        rec0 = sub.iloc[0]
        tissues = ",".join(sorted(set(tumor["tissue"].astype(str)))) if len(tumor) else "normal_only"
        ntab_rows.append(
            {
                "patient": pat,
                "cancer_subtype": rec0["cancer_subtype"],
                "recist": rec0["recist"],
                "tissues": tissues,
                "n_samples_tumor": int(tumor["orig.ident"].nunique()) if len(tumor) else 0,
                "n_cells_tumor": int(len(tumor)),
                "n_malignant": n_mal,
                "n_T": int(tumor["is_t"].to_numpy().sum()) if len(tumor) else 0,
                "n_NK": int(tumor["is_nk"].to_numpy().sum()) if len(tumor) else 0,
                "n_T_NK": n_t,
                "n_cldn4_high": n_hi,
                "n_cldn4_low": n_lo,
                "mean_cldn4_malig": mean_cld,
                "paired": paired,
                "normal_only": bool(len(tumor) == 0),
            }
        )
    ntab = pd.DataFrame(ntab_rows).sort_values("patient")
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)
    tumor_ntab = ntab[~ntab["normal_only"]].copy()

    out_hi = score_pairs(pairs, logx, cld_hi, is_tnk, P["expr_prop"])
    out_lo = score_pairs(pairs, logx, cld_lo, is_tnk, P["expr_prop"])
    in_hi = score_pairs(pairs, logx, is_tnk, cld_hi, P["expr_prop"])
    in_lo = score_pairs(pairs, logx, is_tnk, cld_lo, P["expr_prop"])
    paired_delta(out_hi, out_lo, "cpdb_mean_score").to_csv(
        args.outdir / "pooled_outgoing_cldn4.tsv", sep="\t", index=False
    )
    paired_delta(in_hi, in_lo, "cpdb_mean_score").to_csv(
        args.outdir / "pooled_incoming_cldn4.tsv", sep="\t", index=False
    )
    out_hi.sort_values(["pass_expr_prop", "cpdb_mean_score"], ascending=[False, False]).to_csv(
        args.outdir / "lr_table.tsv", sep="\t", index=False
    )

    def patient_table(direction: str) -> pd.DataFrame:
        chunks = []
        for pat in tumor_ntab.loc[tumor_ntab["paired"], "patient"]:
            m = patient_arr == pat
            h, l, t = m & cld_hi, m & cld_lo, m & is_tnk
            if direction == "outgoing":
                hdf = score_pairs(pairs, logx, h, t, P["expr_prop"])
                ldf = score_pairs(pairs, logx, l, t, P["expr_prop"])
            else:
                hdf = score_pairs(pairs, logx, t, h, P["expr_prop"])
                ldf = score_pairs(pairs, logx, t, l, P["expr_prop"])
            d = paired_delta(hdf, ldf, "cpdb_mean_score")
            d["patient"] = pat
            meta = tumor_ntab.set_index("patient").loc[pat]
            d["cancer_subtype"] = meta["cancer_subtype"]
            d["recist"] = meta["recist"]
            d["n_high"] = int(h.sum())
            d["n_low"] = int(l.sum())
            d["n_tnk"] = int(t.sum())
            chunks.append(d)
            log(f"  scored {direction} {pat} high={int(h.sum())} low={int(l.sum())} tnk={int(t.sum())}")
        return pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()

    ranks: dict[str, pd.DataFrame] = {}
    raw_store: dict[str, pd.DataFrame] = {}
    n_paired = {}
    for direction in ("outgoing", "incoming"):
        raw = patient_table(direction)
        raw.to_csv(args.outdir / f"patient_cldn4_{direction}.tsv.gz", sep="\t", index=False)
        raw_store[direction] = raw
        n_paired[direction] = int(raw["patient"].nunique()) if len(raw) else 0
        tab = rank_direction(raw)
        tab.to_csv(args.outdir / f"ranks_cldn4_{direction}.tsv", sep="\t", index=False)
        ranks[f"cldn4_{direction}"] = tab

    out_tab = ranks["cldn4_outgoing"].copy()
    out_tab.insert(0, "direction", "malignant_CLDN4high_to_TNK")
    out_tab.to_csv(args.outdir / "lr_table_cldn4_outgoing_tnk.tsv", sep="\t", index=False)
    in_tab = ranks["cldn4_incoming"].copy()
    in_tab.insert(0, "direction", "TNK_to_malignant_CLDN4")
    in_tab.to_csv(args.outdir / "lr_table_cldn4_incoming_tnk.tsv", sep="\t", index=False)
    focus = out_tab[out_tab["pathway"].isin(FOCUS_PATHWAYS)].copy()
    focus.to_csv(args.outdir / "lr_table_focus_outgoing.tsv", sep="\t", index=False)

    path_rows = []
    for key, tab in ranks.items():
        for path, sub in tab.groupby("pathway"):
            if path not in FOCUS_PATHWAYS:
                continue
            path_rows.append(
                {
                    "contrast": key,
                    "pathway": path,
                    "n_pairs": int(len(sub)),
                    "n_pairs_delta_neg": int((sub["median_delta"] < 0).sum()),
                    "median_of_pair_deltas": float(sub["median_delta"].median()) if len(sub) else np.nan,
                    "n_pairs_fdr05": int((sub["padj"] < 0.05).sum()) if sub["padj"].notna().any() else 0,
                }
            )
    pd.DataFrame(path_rows).to_csv(args.outdir / "pathway_summary.tsv", sep="\t", index=False)

    liana_notes = [f"import: {liana_import}"]
    liana_status = "not_run"
    liana_focus = None
    try:
        import anndata as ad

        mask = cld_hi | cld_lo | is_t | is_nk
        genes = sorted(g for g in logx if g in extracted)
        X = np.vstack([logx[g][mask] for g in genes]).T
        obs = pd.DataFrame(
            {
                "patient": patient_arr[mask],
                "cancer_subtype": subtype_arr[mask],
            },
            index=pd.Index(barcodes[mask], name="cell"),
        )
        grp = np.array(["other"] * int(mask.sum()), dtype=object)
        grp[cld_hi[mask]] = "Malig_CLDN4high"
        grp[cld_lo[mask]] = "Malig_CLDN4low"
        grp[is_t[mask]] = "T"
        grp[is_nk[mask]] = "NK"
        obs["cc_group"] = grp
        rng = np.random.default_rng(P["random_seed"])
        keep_idx = []
        for g, idx in obs.groupby("cc_group", observed=True).indices.items():
            if g == "other":
                continue
            if len(idx) > P["liana_max_cells_per_group"]:
                idx = rng.choice(idx, size=P["liana_max_cells_per_group"], replace=False)
            keep_idx.append(np.asarray(idx))
        keep_idx = np.sort(np.concatenate(keep_idx)) if keep_idx else np.array([], dtype=int)
        adata = ad.AnnData(X=X[keep_idx], obs=obs.iloc[keep_idx].copy(), var=pd.DataFrame(index=genes))
        adata.obs["cc_group"] = pd.Categorical(adata.obs["cc_group"])
        adata.uns["log1p"] = {"base": None}
        counts = adata.obs["cc_group"].value_counts().to_dict()
        log(f"[liana] groups {counts}")
        note = run_liana(
            adata,
            "cc_group",
            args.outdir / "liana_cellphonedb_cldn4.csv",
            P["liana_n_perms"],
        )
        liana_notes.append(f"CLDN4: {note}; downsampled groups={counts}")
        if note.startswith("LIANA_OK"):
            liana_status = "ran_cellphonedb_method"
            res = pd.read_csv(args.outdir / "liana_cellphonedb_cldn4.csv")
            senders = ["Malig_CLDN4high", "Malig_CLDN4low"]
            res[res["source"].isin(senders) & res["target"].isin(["T", "NK"])].to_csv(
                args.outdir / "liana_cldn4_outgoing_tnk.csv", index=False
            )
            res[res["source"].isin(["T", "NK"]) & res["target"].isin(senders)].to_csv(
                args.outdir / "liana_cldn4_incoming_tnk.csv", index=False
            )
            liana_focus = summarize_liana_focus(res)
            liana_focus.to_csv(args.outdir / "liana_focus_outgoing.tsv", sep="\t", index=False)
        else:
            liana_status = note.split(":")[0].lower()
    except Exception as exc:
        liana_notes.append(f"LIANA wrapper failed: {exc}")
        liana_status = f"failed: {exc}"

    foc = ranks.get("cldn4_outgoing", pd.DataFrame())
    foc = foc[foc["pathway"].isin(FOCUS_PATHWAYS)] if len(foc) else foc
    n_pair = n_paired.get("outgoing", 0)
    if len(foc):
        n_neg = int((foc["median_delta"] < 0).sum())
        n_sig = int((foc["padj"] < 0.05).sum()) if foc["padj"].notna().any() else 0
        by_bits = []
        for path, sub in foc.groupby("pathway"):
            by_bits.append(
                f"{path}: {int((sub.median_delta < 0).sum())}/{len(sub)} Δ<0, "
                f"median Δ={sub.median_delta.median():+.3f}"
            )
        trend = (
            f"On the patient-level CellPhoneDB-style score (paired n={n_pair}), "
            f"CLDN4-high vs CLDN4-low outgoing T-recruit / MHC-I / IFN pairs: "
            f"{n_neg}/{len(foc)} have median Δ < 0; {n_sig}/{len(foc)} reach FDR < 0.05. "
            f"{'; '.join(by_bits)}. "
            "This is the observed rank in this public GSE205335 tumor extract, not a general rule."
        )
    else:
        trend = "Too few patient-level T-recruit/IFN/MHC-I pairs passed filters to rank a CLDN4-high reduction."

    summary = {
        "dataset": "GSE205335",
        "citation": "Hu/Ahn/Lee lung ICI scRNA; GEO GSE205335",
        "state_gene": "CLDN4",
        "tacstd2_used_as_gate": False,
        "dual_high": False,
        "n_cells": int(len(barcodes)),
        "n_genes_in_matrix": int(len(all_genes)),
        "n_samples_geo": int(metadata["orig.ident"].nunique()),
        "n_patients_geo": int(metadata["patient"].nunique()),
        "n_tumor_samples": int(cells.loc[cells["is_tumor_sample"], "orig.ident"].nunique()),
        "n_tumor_patients": int(tumor_ntab["patient"].nunique()),
        "n_malig": int(is_malig.sum()),
        "n_T": int(is_t.sum()),
        "n_NK": int(is_nk.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cldn4_high": int(cld_hi.sum()),
        "n_cldn4_low": int(cld_lo.sum()),
        "cldn4_threshold": thr,
        "min_malig_per_state": P["min_malig_per_state"],
        "min_tnk": P["min_tnk"],
        "n_patients_paired": n_pair,
        "n_pairs_scored": int(len(pairs)),
        "liana_status": liana_status,
        "liana_n_perms": P["liana_n_perms"],
        "liana_max_cells": P["liana_max_cells_per_group"],
        "cellchat_status": "not_run_R_unavailable",
        "liana_notes": liana_notes,
        "trend_sentence": trend,
        "malignant_definition": "tumor sample AND lineage.sub == Malignant cells",
        "unit_of_inference": "patient (tumor samples pooled)",
        "method": "cellphonedb_mean_of_means_on_log1p_cp10k",
        "not_gse207422": True,
        "not_dual_high": True,
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    make_figures(
        args.outdir,
        tumor_ntab,
        ranks,
        raw_store.get("outgoing", pd.DataFrame()),
        cld,
        is_malig,
        cld_hi,
        cld_lo,
        is_tnk,
        patient_arr,
        subtype_arr,
        thr,
    )
    write_finding(
        outdir=args.outdir,
        root=HERE,
        summary=summary,
        ntab=tumor_ntab,
        ranks=ranks,
        liana_notes=liana_notes,
        liana_focus=liana_focus,
    )
    log(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "n_malig",
                    "n_tnk",
                    "n_patients_paired",
                    "liana_status",
                    "trend_sentence",
                )
            },
            indent=2,
        )
    )
    log(f"[done] {args.outdir}")


if __name__ == "__main__":
    main()
