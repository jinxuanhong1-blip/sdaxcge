#!/usr/bin/env python3
"""GSE180581: DNA-PKcs vs Ku70 vs Ku80 depletion, prerank GSEA.

Matched contrast inside each monoallelic HEK293T background:
  DNA-PKcs: DNA-PKcs+/- siDNA-PKcs vs DNA-PKcs+/- siControl
  Ku70:     Ku70+/- siKu70 vs Ku70+/- siControl
  Ku80:     Ku80+/- siKu80 vs Ku80+/- siControl

Positive NES = the set is enriched at the knockdown end of a Welch t
rank (log2(median-of-ratios + 1), KD minus matched siControl).

Engine: scripts/gse180581_dnapk_subunit_gsea/gsea_core.py
weighted KS p=1, 1000 gene-set permutations, seed=42.
BH-FDR is within the scored panel of that run (primary or sensitivity).

GEO sample titles for the DNA-PKcs arm say "KuDNA-PKcs". Characteristics
and the growth protocol identify those six libraries as DNA-PKcs+/-
treated with siControl or siDNA-PKcs. This script follows that.
"""
from __future__ import annotations

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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank  # noqa: E402

DATA = ROOT / "methods" / "gse180581_dnapk_subunit_gsea" / "data"
OUT = ROOT / "methods" / "gse180581_dnapk_subunit_gsea"
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

SUBUNIT_GENES = ("XRCC6", "XRCC5", "PRKDC")
DROP_FOR_SENSITIVITY = set(SUBUNIT_GENES)
FAMILY_ORDER = ["DNA-sensing", "STING", "IFN", "APM"]
CONTRASTS = [
    {
        "id": "DNA-PKcs",
        "kd": ["A19", "A20", "A21"],
        "ctrl": ["A16", "A17", "A18"],
        "target": "PRKDC",
    },
    {
        "id": "Ku70",
        "kd": ["A07", "A08", "A09"],
        "ctrl": ["A04", "A05", "A06"],
        "target": "XRCC6",
    },
    {
        "id": "Ku80",
        "kd": ["A13", "A14", "A15"],
        "ctrl": ["A10", "A11", "A12"],
        "target": "XRCC5",
    },
]
FOCUS = {
    "target": ["XRCC6", "XRCC5", "PRKDC"],
    "DNA-sensing": [
        "CGAS",
        "STING1",
        "IFI16",
        "ZBP1",
        "AIM2",
        "DDX41",
        "DHX9",
        "DHX36",
        "MRE11",
        "TREX1",
        "POLR3A",
        "LRRFIP1",
        "HMGB1",
    ],
    "IFN": [
        "IFNB1",
        "IFNA1",
        "IFNG",
        "IFNL1",
        "IFNL2",
        "IFNL3",
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "STAT1",
        "STAT2",
        "IRF1",
        "IRF3",
        "IRF7",
        "CXCL10",
        "CXCL9",
    ],
    "APM": [
        "HLA-A",
        "HLA-B",
        "HLA-C",
        "HLA-E",
        "B2M",
        "TAP1",
        "TAP2",
        "TAPBP",
        "PSMB8",
        "PSMB9",
        "PSMB10",
        "NLRC5",
        "ERAP1",
        "CALR",
    ],
}
MIN_COUNT_IN_UNIVERSE = 10
MIN_SAMPLES_FOR_UNIVERSE = 3


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def deseq_size_factors(counts: pd.DataFrame) -> pd.Series:
    mat = counts.to_numpy(dtype=np.float64)
    logged = np.full(mat.shape, np.nan, dtype=np.float64)
    positive = mat > 0
    logged[positive] = np.log(mat[positive])
    valid = positive.any(axis=1)
    log_geo = np.full(mat.shape[0], np.nan)
    log_geo[valid] = np.nanmean(logged[valid], axis=1)
    ok = np.isfinite(log_geo)
    if int(ok.sum()) < 100:
        raise SystemExit("too few genes for median-of-ratios size factors")
    ratios = mat[ok] / np.exp(log_geo[ok])[:, None]
    sf = np.median(ratios, axis=0)
    sf = sf / np.exp(np.mean(np.log(sf)))
    if np.any(~np.isfinite(sf)) or np.any(sf <= 0):
        raise SystemExit(f"bad size factors: {sf}")
    return pd.Series(sf, index=counts.columns, name="size_factor")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    counts = pd.read_csv(DATA / "GSE180581_counts_ensembl.csv.gz")
    samples = pd.read_csv(DATA / "samples.tsv", sep="\t")
    mapping = pd.read_csv(DATA / "ensembl_to_symbol.tsv", sep="\t")
    payload = json.loads((DATA / "genesets.json").read_text())
    sample_ids = samples["sample_id"].tolist()
    if list(counts.columns[1:]) != sample_ids:
        raise SystemExit(
            f"count columns {list(counts.columns[1:])} != samples {sample_ids}"
        )
    counts = counts.set_index("ensembl")[sample_ids].astype(np.float64)
    sym = mapping.drop_duplicates("ensembl").set_index("ensembl")["symbol"]
    counts = counts.join(sym, how="left")
    return counts, samples, payload


def symbol_matrix(counts: pd.DataFrame, norm: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Sum raw and normalized counts onto NCBI symbols."""
    mapped = counts["symbol"].notna()
    raw = counts.loc[mapped].drop(columns="symbol")
    raw = raw.groupby(counts.loc[mapped, "symbol"]).sum()
    nrm = norm.loc[mapped.index[mapped]]
    nrm = nrm.groupby(counts.loc[mapped, "symbol"]).sum()
    nrm = nrm.reindex(raw.index)
    return raw, nrm


def welch_table(log_mat: pd.DataFrame, raw: pd.DataFrame, kd: list[str], ctrl: list[str]) -> pd.DataFrame:
    eh = log_mat[kd].to_numpy(dtype=np.float64)
    el = log_mat[ctrl].to_numpy(dtype=np.float64)
    rh = raw[kd].to_numpy(dtype=np.float64)
    rl = raw[ctrl].to_numpy(dtype=np.float64)
    ok = np.isfinite(eh).all(axis=1) & np.isfinite(el).all(axis=1)
    vh = np.var(eh, axis=1, ddof=1)
    vl = np.var(el, axis=1, ddof=1)
    ok &= (vh + vl) > 1e-8
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    se = np.sqrt(vh / eh.shape[1] + vl / el.shape[1])
    tstat = np.full(log_mat.shape[0], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat[ok] = (mh[ok] - ml[ok]) / se[ok]
    welch = stats.ttest_ind(eh, el, axis=1, equal_var=False, nan_policy="omit")
    out = pd.DataFrame(
        {
            "log2FC": mh - ml,
            "t": tstat,
            "welch_p": welch.pvalue,
            "mean_log2_kd": mh,
            "mean_log2_ctrl": ml,
            "mean_count_kd": rh.mean(axis=1),
            "mean_count_ctrl": rl.mean(axis=1),
        },
        index=log_mat.index,
    )
    out.index.name = "symbol"
    return out


def rank_from_table(tab: pd.DataFrame) -> pd.Series:
    s = tab["t"].replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
    return s


def run_gsea(rank: pd.Series, gene_sets: dict[str, list[str]]) -> pd.DataFrame:
    g = gsea_prerank(rank, gene_sets, nperm=NPERM, seed=SEED)
    if g.empty:
        raise SystemExit("no gene set passed size filters")
    g["fdr"] = bh_fdr(g["nom_p"])
    return g


def sets_for(payload: dict, drop_subunits: bool) -> dict[str, list[str]]:
    out = {}
    for sid, meta in payload["sets"].items():
        genes = list(meta["genes"])
        if drop_subunits:
            genes = [g for g in genes if g not in DROP_FOR_SENSITIVITY]
        out[sid] = genes
    return out


def family_summary(gsea: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    if "family" in gsea.columns:
        df = gsea
    else:
        df = gsea.merge(meta[["family"]], left_on="term", right_index=True, how="left")
    rows = []
    for contrast, sub in df.groupby("contrast", sort=False):
        for family in FAMILY_ORDER:
            s = sub[sub["family"] == family]
            if s.empty:
                continue
            rows.append(
                {
                    "contrast": contrast,
                    "family": family,
                    "n_sets": int(len(s)),
                    "median_nes": float(s["nes"].median()),
                    "min_nes": float(s["nes"].min()),
                    "max_nes": float(s["nes"].max()),
                    "n_fdr05_up": int(((s["fdr"] < 0.05) & (s["nes"] > 0)).sum()),
                    "n_fdr05_down": int(((s["fdr"] < 0.05) & (s["nes"] < 0)).sum()),
                    "n_nom05_up": int(((s["nom_p"] < 0.05) & (s["nes"] > 0)).sum()),
                    "n_nom05_down": int(((s["nom_p"] < 0.05) & (s["nes"] < 0)).sum()),
                    "n_nes_positive": int((s["nes"] > 0).sum()),
                }
            )
    return pd.DataFrame(rows)


def plot_nes_heatmap(wide: pd.DataFrame, path: Path) -> None:
    contrast_ids = [c["id"] for c in CONTRASTS]
    nes_cols = [f"nes_{c}" for c in contrast_ids]
    fdr_cols = [f"fdr_{c}" for c in contrast_ids]
    wide = wide.dropna(subset=nes_cols, how="all").reset_index(drop=True)
    mat = wide[nes_cols].to_numpy(dtype=float)
    fdr = wide[fdr_cols].to_numpy(dtype=float)
    fig_h = max(6.5, 0.32 * len(wide) + 1.4)
    fig, ax = plt.subplots(figsize=(7.4, fig_h))
    vmax = np.nanmax(np.abs(mat))
    vmax = max(1.5, float(np.ceil(vmax * 10) / 10))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(contrast_ids)))
    ax.set_xticklabels(contrast_ids, fontsize=11)
    ax.set_yticks(range(len(wide)))
    ax.set_yticklabels(wide["label"].tolist(), fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if not np.isfinite(val):
                continue
            star = "*" if fdr[i, j] < 0.05 else ""
            color = "white" if abs(val) > 0.62 * vmax else "black"
            ax.text(j, i, f"{val:+.2f}{star}", ha="center", va="center", fontsize=7, color=color)
    # family brackets via left ticks already labeled; add a thin color strip
    fam_colors = {
        "DNA-sensing": "#4C78A8",
        "STING": "#F58518",
        "IFN": "#E45756",
        "APM": "#54A24B",
    }
    for i, fam in enumerate(wide["family"]):
        ax.add_patch(
            plt.Rectangle(
                (-0.72, i - 0.5),
                0.12,
                1,
                color=fam_colors[fam],
                clip_on=False,
                linewidth=0,
            )
        )
    ax.set_xlim(-0.5, len(contrast_ids) - 0.5)
    ax.set_title(
        "GSE180581 prerank NES\nKD minus matched siControl  (* BH-FDR < 0.05 within panel)",
        fontsize=11,
    )
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("NES (positive = up in KD)")
    handles = [
        plt.Line2D([0], [0], marker="s", color="none", markerfacecolor=fam_colors[f], markersize=8, label=f)
        for f in FAMILY_ORDER
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.18, 1.0), frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_qc(log_mat: pd.DataFrame, samples: pd.DataFrame, focus: pd.DataFrame, path: Path) -> None:
    # Top variable genes among those with a real dynamic range. The full
    # rank universe includes sparse genes whose variance is sampling noise.
    spread = log_mat.max(axis=1) - log_mat.min(axis=1)
    eligible = log_mat.loc[spread >= 0.5]
    var = eligible.var(axis=1).sort_values(ascending=False)
    top = eligible.loc[var.index[:2000]]
    x = top.to_numpy().T
    x = x - x.mean(axis=0, keepdims=True)
    _u, s, vt = np.linalg.svd(x, full_matrices=False)
    pcs = x @ vt[:2].T
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))
    group_color = {
        "parental_siControl": "#9E9E9E",
        "Ku70_siControl": "#F4A3A0",
        "Ku70_siKu70": "#E64B35",
        "Ku80_siControl": "#8FD0C4",
        "Ku80_siKu80": "#00A087",
        "DNAPKcs_siControl": "#A9B6D6",
        "DNAPKcs_siDNAPKcs": "#3C5488",
    }
    ax = axes[0]
    for _, row in samples.iterrows():
        i = list(log_mat.columns).index(row["sample_id"])
        ax.scatter(
            pcs[i, 0],
            pcs[i, 1],
            s=46,
            color=group_color[row["group"]],
            label=row["group"],
            zorder=3,
        )
        ax.text(pcs[i, 0], pcs[i, 1], row["sample_id"], fontsize=6, ha="left", va="bottom")
    handles, labels = ax.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax.legend(uniq.values(), uniq.keys(), fontsize=6, frameon=False, loc="best")
    var_exp = (s[:2] ** 2) / (s ** 2).sum()
    ax.set_xlabel(f"PC1 ({100 * var_exp[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({100 * var_exp[1]:.1f}%)")
    ax.set_title("Top 2000 variable genes with range ≥ 0.5 log2")
    ax.axhline(0, color="#dddddd", lw=0.6)
    ax.axvline(0, color="#dddddd", lw=0.6)

    ax = axes[1]
    sub = focus[focus["class"] == "target"].copy()
    genes = list(SUBUNIT_GENES)
    xloc = np.arange(len(genes))
    width = 0.24
    colors = {"DNA-PKcs": "#3C5488", "Ku70": "#E64B35", "Ku80": "#00A087"}
    for i, cid in enumerate([c["id"] for c in CONTRASTS]):
        vals = [
            float(sub[(sub["contrast"] == cid) & (sub["symbol"] == g)]["log2FC"].iloc[0])
            for g in genes
        ]
        ax.bar(xloc + (i - 1) * width, vals, width=width, color=colors[cid], label=cid)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(xloc)
    ax.set_xticklabels(["Ku70\nXRCC6", "Ku80\nXRCC5", "DNA-PKcs\nPRKDC"])
    ax.set_ylabel("log2FC (KD − matched siControl)")
    ax.set_title("On-target depletion")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_focus_heatmap(focus: pd.DataFrame, path: Path) -> None:
    order = []
    for cls, genes in FOCUS.items():
        if cls == "target":
            continue
        order.extend(genes)
    contrast_ids = [c["id"] for c in CONTRASTS]
    mat = np.full((len(order), len(contrast_ids)), np.nan)
    expressed = np.zeros_like(mat, dtype=bool)
    for i, gene in enumerate(order):
        for j, cid in enumerate(contrast_ids):
            hit = focus[(focus["symbol"] == gene) & (focus["contrast"] == cid)]
            if hit.empty:
                continue
            mat[i, j] = float(hit["log2FC"].iloc[0])
            expressed[i, j] = float(hit["mean_count_kd"].iloc[0]) >= 1 or float(hit["mean_count_ctrl"].iloc[0]) >= 1
    fig, ax = plt.subplots(figsize=(6.2, 0.28 * len(order) + 1.6))
    show = mat.copy()
    show[~expressed] = np.nan
    vmax = np.nanmax(np.abs(show)) if np.isfinite(show).any() else 1
    vmax = max(1.0, float(np.ceil(vmax * 10) / 10))
    im = ax.imshow(np.ma.masked_invalid(show), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(contrast_ids)))
    ax.set_xticklabels(contrast_ids)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7)
    for i in range(show.shape[0]):
        for j in range(show.shape[1]):
            if not expressed[i, j] or not np.isfinite(mat[i, j]):
                ax.text(j, i, "·", ha="center", va="center", color="#888888", fontsize=8)
            else:
                val = mat[i, j]
                color = "white" if abs(val) > 0.62 * vmax else "black"
                ax.text(j, i, f"{val:+.2f}", ha="center", va="center", fontsize=6, color=color)
    ax.set_title("Focus genes, log2FC\n· mean raw count < 1 in both arms")
    cbar = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("log2FC (KD − matched siControl)")
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_sensitivity(wide: pd.DataFrame, path: Path) -> None:
    """Primary NES vs NES after removing XRCC5, XRCC6, and PRKDC."""
    terms = [
        ("REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES", "Reactome STING"),
        ("CUSTOM_CYTOSOLIC_DNA_SENSORS", "Custom DNA sensors"),
        ("REACTOME_CYTOSOLIC_SENSORS_OF_PATHOGEN_ASSOCIATED_DNA", "Reactome DNA sensors"),
    ]
    colors = {"DNA-PKcs": "#3C5488", "Ku70": "#E64B35", "Ku80": "#00A087"}
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 4.0), sharey=True)
    x = np.arange(len(terms))
    width = 0.36
    for ax, cid in zip(axes, [c["id"] for c in CONTRASTS]):
        primary, dropped = [], []
        for term, _label in terms:
            row = wide.loc[wide["term"] == term].iloc[0]
            primary.append(float(row[f"nes_{cid}"]))
            dropped.append(float(row[f"nes_drop_subunits_{cid}"]))
        ax.bar(x - width / 2, primary, width, color=colors[cid], label="With subunit genes")
        ax.bar(
            x + width / 2,
            dropped,
            width,
            color=colors[cid],
            alpha=0.35,
            edgecolor=colors[cid],
            linewidth=0.8,
            label="XRCC5/XRCC6/PRKDC removed",
        )
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels([lab for _t, lab in terms], rotation=20, ha="right", fontsize=8)
        ax.set_title(cid)
        ax.set_ylabel("NES" if cid == "DNA-PKcs" else "")
    axes[0].legend(frameon=False, fontsize=7, loc="lower left")
    fig.suptitle("DNA-sensing / STING NES is carried by the depleted subunit gene", fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    fig.savefig(path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    log("load")
    counts, samples, payload = load_inputs()
    sample_ids = samples["sample_id"].tolist()
    sf = deseq_size_factors(counts[sample_ids])
    sf.to_csv(TAB / "size_factors.tsv", sep="\t", header=True)
    norm = counts[sample_ids].div(sf, axis=1)
    raw_all, norm_all = symbol_matrix(counts, norm)
    keep = (raw_all >= MIN_COUNT_IN_UNIVERSE).sum(axis=1) >= MIN_SAMPLES_FOR_UNIVERSE
    raw_sym = raw_all.loc[keep]
    norm_sym = norm_all.loc[keep]
    log_mat = np.log2(norm_sym + 1.0)
    log_all = np.log2(norm_all + 1.0)
    log(
        f"symbols in rank universe {log_mat.shape[0]} "
        f"(mapped rows {int(counts['symbol'].notna().sum())}, unmapped ENSG {int(counts['symbol'].isna().sum())})"
    )
    meta = pd.DataFrame(payload["sets"]).T
    meta.index.name = "term"
    gene_sets = sets_for(payload, drop_subunits=False)
    gene_sets_sens = sets_for(payload, drop_subunits=True)

    gsea_rows = []
    sens_rows = []
    de_rows = []
    focus_rows = []
    ranks = {}
    for contrast in CONTRASTS:
        cid = contrast["id"]
        log(f"contrast {cid}")
        tab = welch_table(log_mat, raw_sym, contrast["kd"], contrast["ctrl"])
        rank = rank_from_table(tab)
        ranks[cid] = rank
        tab = tab.loc[rank.index].copy()
        tab.insert(0, "symbol", tab.index)
        tab.insert(0, "contrast", cid)
        de_rows.append(tab.reset_index(drop=True))
        g = run_gsea(rank, gene_sets)
        g.insert(0, "contrast", cid)
        g.insert(1, "run", "primary")
        gsea_rows.append(g)
        gs = run_gsea(rank, gene_sets_sens)
        gs.insert(0, "contrast", cid)
        gs.insert(1, "run", "drop_XRCC5_XRCC6_PRKDC")
        sens_rows.append(gs)
        for cls, genes in FOCUS.items():
            for gene in genes:
                if gene not in raw_all.index:
                    focus_rows.append(
                        {
                            "contrast": cid,
                            "class": cls,
                            "symbol": gene,
                            "present": False,
                            "in_rank_universe": False,
                        }
                    )
                    continue
                mean_kd = float(raw_all.loc[gene, contrast["kd"]].mean())
                mean_ctrl = float(raw_all.loc[gene, contrast["ctrl"]].mean())
                if gene not in tab["symbol"].values:
                    focus_rows.append(
                        {
                            "contrast": cid,
                            "class": cls,
                            "symbol": gene,
                            "present": True,
                            "in_rank_universe": False,
                            "log2FC": float(
                                log_all.loc[gene, contrast["kd"]].mean()
                                - log_all.loc[gene, contrast["ctrl"]].mean()
                            ),
                            "mean_count_kd": mean_kd,
                            "mean_count_ctrl": mean_ctrl,
                        }
                    )
                    continue
                rec = tab.loc[tab["symbol"] == gene].iloc[0].to_dict()
                rec["class"] = cls
                rec["present"] = True
                rec["in_rank_universe"] = True
                focus_rows.append(rec)
        on = tab.loc[tab["symbol"] == contrast["target"]].iloc[0]
        log(
            f"  {contrast['target']} log2FC={on['log2FC']:+.3f} t={on['t']:+.2f} "
            f"p={on['welch_p']:.3g} rank n={len(rank)}"
        )

    gsea = pd.concat(gsea_rows, ignore_index=True)
    sens = pd.concat(sens_rows, ignore_index=True)
    de = pd.concat(de_rows, ignore_index=True)
    focus = pd.DataFrame(focus_rows)
    gsea = gsea.merge(meta[["family", "label", "library"]], left_on="term", right_index=True, how="left")
    sens = sens.merge(meta[["family", "label"]], left_on="term", right_index=True, how="left")

    fam = family_summary(gsea, meta)
    fam_s = family_summary(sens, meta)
    fam_s.insert(1, "run", "drop_XRCC5_XRCC6_PRKDC")
    fam.insert(1, "run", "primary")

    # wide NES comparison
    pieces = []
    for cid in [c["id"] for c in CONTRASTS]:
        sub = gsea[gsea["contrast"] == cid][
            ["term", "nes", "fdr", "nom_p", "n_set_in_rank", "mean_stat"]
        ].rename(
            columns={
                "nes": f"nes_{cid}",
                "fdr": f"fdr_{cid}",
                "nom_p": f"nom_p_{cid}",
                "n_set_in_rank": f"n_{cid}",
                "mean_stat": f"mean_t_{cid}",
            }
        )
        pieces.append(sub.set_index("term"))
    wide = pieces[0].join(pieces[1]).join(pieces[2])
    sens_nes = sens.pivot(index="term", columns="contrast", values="nes")
    sens_nes = sens_nes.rename(columns={c: f"nes_drop_subunits_{c}" for c in sens_nes.columns})
    wide = wide.join(sens_nes)
    wide = meta[["family", "label", "library"]].join(wide)
    wide = wide.reset_index().rename(columns={"index": "term"})
    wide["_fam"] = wide["family"].map({f: i for i, f in enumerate(FAMILY_ORDER)})
    wide = wide.sort_values(["_fam", "label"]).drop(columns="_fam")

    # rank correlations
    ids = [c["id"] for c in CONTRASTS]
    corr_rows = []
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            joined = pd.concat([ranks[a].rename(a), ranks[b].rename(b)], axis=1, join="inner")
            rho, p = stats.spearmanr(joined[a], joined[b])
            nes_a = gsea.loc[gsea["contrast"] == a, ["term", "nes"]].set_index("term")["nes"]
            nes_b = gsea.loc[gsea["contrast"] == b, ["term", "nes"]].set_index("term")["nes"]
            both = pd.concat([nes_a.rename("a"), nes_b.rename("b")], axis=1, join="inner")
            rho_n, p_n = stats.spearmanr(both["a"], both["b"])
            corr_rows.append(
                {
                    "a": a,
                    "b": b,
                    "spearman_welch_t": float(rho),
                    "spearman_welch_t_p": float(p),
                    "n_genes": int(len(joined)),
                    "spearman_panel_nes": float(rho_n),
                    "spearman_panel_nes_p": float(p_n),
                    "n_sets": int(len(both)),
                }
            )
    corr = pd.DataFrame(corr_rows)

    breadth = []
    for cid, sub in de.groupby("contrast"):
        breadth.append(
            {
                "contrast": cid,
                "n_rank": int(len(sub)),
                "n_abs_log2fc_gt1": int((sub["log2FC"].abs() > 1).sum()),
                "n_abs_log2fc_gt1_welch_p_lt_0.05": int(
                    ((sub["log2FC"].abs() > 1) & (sub["welch_p"] < 0.05)).sum()
                ),
                "n_welch_p_lt_0.05": int((sub["welch_p"] < 0.05).sum()),
            }
        )
    breadth = pd.DataFrame(breadth)

    # Sample-sample Pearson on genes with mean raw count >= 20.
    expressed = raw_sym.index[raw_sym.mean(axis=1) >= 20]
    log_expr = log_mat.loc[log_mat.index.intersection(expressed)]
    sample_corr = pd.DataFrame(
        np.corrcoef(log_expr.to_numpy().T),
        index=log_expr.columns,
        columns=log_expr.columns,
    )
    rep_rows = []
    for group, sub in samples.groupby("group"):
        cols = sub["sample_id"].tolist()
        c = sample_corr.loc[cols, cols].to_numpy()
        iu = np.triu_indices(len(cols), k=1)
        rep_rows.append(
            {
                "group": group,
                "n": len(cols),
                "n_genes_mean_count_ge_20": int(log_expr.shape[0]),
                "min_pearson": float(c[iu].min()) if len(iu[0]) else np.nan,
                "median_pearson": float(np.median(c[iu])) if len(iu[0]) else np.nan,
            }
        )
    reps = pd.DataFrame(rep_rows)
    sample_corr.to_csv(TAB / "sample_pearson.tsv", sep="\t")

    universe = set(log_mat.index)
    cov_rows = []
    for sid, meta_row in meta.iterrows():
        members = [g for g in gene_sets[sid] if g]
        in_u = sorted(set(members) & universe)
        cov_rows.append(
            {
                "term": sid,
                "family": meta_row["family"],
                "label": meta_row["label"],
                "n_set": len(set(members)),
                "n_in_rank_universe": len(in_u),
                "scored": len(in_u) >= 8,
            }
        )
    coverage = pd.DataFrame(cov_rows)

    coverage.to_csv(TAB / "geneset_coverage.tsv", sep="\t", index=False)
    pd.Series(raw_sym.sum(), name="library_size").to_csv(TAB / "library_sizes.tsv", sep="\t", header=True)
    gsea.to_csv(TAB / "gsea_primary.tsv", sep="\t", index=False)
    sens.to_csv(TAB / "gsea_drop_subunits.tsv", sep="\t", index=False)
    wide.to_csv(TAB / "gsea_subunit_comparison.tsv", sep="\t", index=False)
    fam.to_csv(TAB / "gsea_family_summary.tsv", sep="\t", index=False)
    pd.concat([fam, fam_s], ignore_index=True).to_csv(
        TAB / "gsea_family_summary_with_sensitivity.tsv", sep="\t", index=False
    )
    focus.to_csv(TAB / "focus_genes.tsv", sep="\t", index=False)
    corr.to_csv(TAB / "subunit_rank_correlation.tsv", sep="\t", index=False)
    breadth.to_csv(TAB / "de_breadth.tsv", sep="\t", index=False)
    reps.to_csv(TAB / "replicate_correlation.tsv", sep="\t", index=False)
    for cid, sub in de.groupby("contrast"):
        safe = cid.replace("-", "")
        sub.to_csv(TAB / f"de_{safe}.tsv.gz", sep="\t", index=False, compression="gzip")

    plot_nes_heatmap(wide, FIG / "fig_nes_heatmap.png")
    plot_qc(log_mat, samples, focus, FIG / "fig_qc_pca_targets.png")
    plot_focus_heatmap(focus, FIG / "fig_focus_log2fc.png")
    plot_sensitivity(wide, FIG / "fig_subunit_gene_sensitivity.png")

    summary = {
        "nperm": NPERM,
        "seed": SEED,
        "n_symbols_universe": int(log_mat.shape[0]),
        "n_unmapped_ensembl": int(counts["symbol"].isna().sum()),
        "contrasts": CONTRASTS,
        "family_summary": fam.to_dict(orient="records"),
        "correlation": corr.to_dict(orient="records"),
        "breadth": breadth.to_dict(orient="records"),
        "replicates": reps.to_dict(orient="records"),
        "size_factors": sf.round(4).to_dict(),
        "headline": wide[
            [
                "term",
                "family",
                "label",
                "nes_DNA-PKcs",
                "fdr_DNA-PKcs",
                "nes_Ku70",
                "fdr_Ku70",
                "nes_Ku80",
                "fdr_Ku80",
            ]
        ].to_dict(orient="records"),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2))
    log("wrote tables and figures")
    print(fam.to_string(index=False))
    print(corr.to_string(index=False))
    print(breadth.to_string(index=False))
    print(reps.to_string(index=False))


if __name__ == "__main__":
    main()
