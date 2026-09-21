#!/usr/bin/env python3
"""GSE137244 KL vs KP on log2(FPKM+1), with Epcam and leave-one filters.

Primary contrast is five KrasG12D;Lkb1 cell lines versus five KrasG12D;Trp53
cell lines. normal-lung-RNA is excluded. Leave-one-library results are stored
with their reduced n and are not substituted for the 5-vs-5 rows.

The tight-junction catalog TJ7 is the seven-gene mean used in the earlier
public KL-vs-KP script (Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, Ocln). The
Epcam membership filter keeps TJ7 genes whose Pearson r with Epcam on these
10 lines is at least 0.5. That rule does not use the KL/KP labels. Leave-one
gene drops on that filtered set do use the labels; their randomization
p-values re-run the drop under every 5-vs-5 label assignment.
"""

from __future__ import annotations

import hashlib
import json
import urllib.request
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/"
    "GSE137244_counts.fpkm.csv.gz"
)
FPKM_NAME = "GSE137244_counts.fpkm.csv.gz"

# GEO series matrix sample titles. The last library is normal lung, not a line.
KP_LIBS = [
    "B6AL10-1-RNA",
    "B6AL10-2-RNA",
    "B6AL10-3-RNA",
    "B6AL10-4-RNA",
    "B6AL10-5-RNA",
]
KL_LIBS = [
    "KL155mix-control-2-RNA",
    "KL47-1-untreated-1-RNA",
    "KLC-RNA",
    "KLD-RNA",
    "KLE-RNA",
]
NORMAL = "normal-lung-RNA"
SHORT = {
    "B6AL10-1-RNA": "KP B6AL10-1",
    "B6AL10-2-RNA": "KP B6AL10-2",
    "B6AL10-3-RNA": "KP B6AL10-3",
    "B6AL10-4-RNA": "KP B6AL10-4",
    "B6AL10-5-RNA": "KP B6AL10-5",
    "KL155mix-control-2-RNA": "KL KL155",
    "KL47-1-untreated-1-RNA": "KL KL47",
    "KLC-RNA": "KL KLC",
    "KLD-RNA": "KL KLD",
    "KLE-RNA": "KL KLE",
    "normal-lung-RNA": "normal lung",
}
GSM = {
    "B6AL10-1-RNA": "GSM4073816",
    "B6AL10-2-RNA": "GSM4073817",
    "B6AL10-3-RNA": "GSM4073818",
    "B6AL10-4-RNA": "GSM4073819",
    "B6AL10-5-RNA": "GSM4073820",
    "KL155mix-control-2-RNA": "GSM4073821",
    "KL47-1-untreated-1-RNA": "GSM4073822",
    "KLC-RNA": "GSM4073823",
    "KLD-RNA": "GSM4073824",
    "KLE-RNA": "GSM4073825",
    "normal-lung-RNA": "GSM4073826",
}

# Prior public definition. Not refit to a target delta.
TJ7 = ["Cldn3", "Cldn4", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
TAU = 0.5
ARMS = KP_LIBS + KL_LIBS

# Locked unadjusted log2(FPKM+1) deltas from the earlier 5-vs-5 reproduction.
LOCKED_DELTA = {"Tacstd2": 3.2382274886828935, "Cldn4": 5.569627161932762}


def download_fpkm() -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    dest = DATA / FPKM_NAME
    if not dest.exists():
        req = urllib.request.Request(
            FPKM_URL, headers={"User-Agent": "gse137244-kl-kp-max-effect/1.0"}
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            dest.write_bytes(resp.read())
    return dest


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def load_fpkm(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    df.index.name = "gene"
    expected = ARMS + [NORMAL]
    missing = [c for c in expected if c not in df.columns]
    if missing:
        raise SystemExit(f"FPKM matrix is missing columns: {missing}")
    return df[expected].apply(pd.to_numeric, errors="raise")


def welch(y: pd.Series, kl: list[str], kp: list[str]) -> dict:
    a = y.loc[kl].to_numpy(dtype=float)
    b = y.loc[kp].to_numpy(dtype=float)
    if np.isnan(a).any() or np.isnan(b).any():
        raise ValueError("NaN in contrast")
    delta = float(a.mean() - b.mean())
    res = stats.ttest_ind(a, b, equal_var=False)
    t = float(res.statistic)
    p = float(res.pvalue)
    df = float(res.df)
    complete = bool(a.min() > b.max() or b.min() > a.max())
    return {
        "n_kl": len(kl),
        "n_kp": len(kp),
        "n5vs5": len(kl) == 5 and len(kp) == 5,
        "mean_kl": float(a.mean()),
        "mean_kp": float(b.mean()),
        "sd_kl": float(a.std(ddof=1)),
        "sd_kp": float(b.std(ddof=1)),
        "delta": delta,
        "welch_t": t,
        "welch_df": df,
        "welch_p": p,
        "abs_t": abs(t),
        "min_kl": float(a.min()),
        "max_kp": float(b.max()),
        "complete_separation": complete,
    }


def mannwhitney(y: pd.Series, kl: list[str], kp: list[str]) -> dict:
    a = y.loc[kl].to_numpy(dtype=float)
    b = y.loc[kp].to_numpy(dtype=float)
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided", method="exact")
    cliffs = (2.0 * float(u) / (len(a) * len(b))) - 1.0
    return {"U": float(u), "mw_p": float(p), "cliffs_delta": float(cliffs)}


def contrast(y: pd.Series, kl: list[str] | None = None, kp: list[str] | None = None) -> dict:
    kl = KL_LIBS if kl is None else kl
    kp = KP_LIBS if kp is None else kp
    out = welch(y, kl, kp)
    if len(kl) and len(kp):
        out.update(mannwhitney(y, kl, kp))
    return out


def residual_on(y: pd.Series, x: pd.Series, cols: list[str]) -> pd.Series:
    yy = y.loc[cols].to_numpy(dtype=float)
    xx = x.loc[cols].to_numpy(dtype=float)
    X = np.column_stack([np.ones(len(cols)), xx])
    beta, *_ = np.linalg.lstsq(X, yy, rcond=None)
    fitted = X @ beta
    return pd.Series(yy - fitted, index=cols, name=y.name)


def fmt(x: float, digits: int = 3) -> str:
    return f"{x:+.{digits}f}"


def write_tsv(df: pd.DataFrame, name: str) -> None:
    path = TABLES / name
    df.to_csv(path, sep="\t", index=False)


def score_mean(log: pd.DataFrame, genes: list[str]) -> pd.Series:
    return log.loc[genes].mean(axis=0)


def leave_one_library(y: pd.Series, endpoint: str) -> list[dict]:
    rows = []
    full = contrast(y)
    full.update({"endpoint": endpoint, "dropped": "none", "arm_dropped": "none"})
    rows.append(full)
    for lib in ARMS:
        kl = [c for c in KL_LIBS if c != lib]
        kp = [c for c in KP_LIBS if c != lib]
        row = contrast(y, kl, kp)
        row.update(
            {
                "endpoint": endpoint,
                "dropped": lib,
                "arm_dropped": "KL" if lib in KL_LIBS else "KP",
                "short": SHORT[lib],
            }
        )
        rows.append(row)
    return rows


def epcam_gene_table(log: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    ep = log.loc["Epcam", ARMS]
    rows = []
    for g in genes:
        y = log.loc[g]
        r = float(np.corrcoef(y.loc[ARMS].to_numpy(dtype=float), ep.to_numpy(dtype=float))[0, 1])
        row = contrast(y)
        row.update({"gene": g, "r_epcam": r, "kept_tau_0.5": r >= TAU})
        rows.append(row)
    return pd.DataFrame(rows)


def tau_sweep(log: pd.DataFrame, gene_stats: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for tau in (-1.0, 0.0, 0.3, 0.5, 0.7, 0.8):
        keep = gene_stats.loc[gene_stats["r_epcam"] >= tau, "gene"].tolist()
        # Preserve TJ7 order.
        keep = [g for g in TJ7 if g in keep]
        if len(keep) < 3:
            continue
        y = score_mean(log, keep)
        row = contrast(y)
        row.update({"tau": tau, "n_genes": len(keep), "genes": ",".join(keep)})
        rows.append(row)
    return pd.DataFrame(rows)


def leave_one_gene(log: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    y = score_mean(log, genes)
    row = contrast(y)
    row.update({"dropped": "none", "n_genes": len(genes), "genes": ",".join(genes)})
    rows.append(row)
    for g in genes:
        keep = [x for x in genes if x != g]
        y = score_mean(log, keep)
        row = contrast(y)
        row.update({"dropped": g, "n_genes": len(keep), "genes": ",".join(keep)})
        rows.append(row)
    return pd.DataFrame(rows)


def select_extreme(grid: pd.DataFrame, column: str) -> pd.Series:
    """Max of `column`. Ties: larger delta, then larger |t|, then earlier row."""
    ranked = grid.copy()
    ranked["_i"] = np.arange(len(ranked))
    ranked = ranked.sort_values(
        [column, "delta", "abs_t", "_i"], ascending=[False, False, False, True]
    )
    return ranked.iloc[0]


def library_scores(log: pd.DataFrame, series: dict[str, pd.Series]) -> pd.DataFrame:
    rows = []
    for lib in ARMS + [NORMAL]:
        arm = "KL" if lib in KL_LIBS else ("KP" if lib in KP_LIBS else "normal")
        rec = {
            "library": lib,
            "short": SHORT[lib],
            "gsm": GSM[lib],
            "arm": arm,
            "in_primary": arm in ("KL", "KP"),
        }
        for name, y in series.items():
            rec[name] = float(y.loc[lib]) if lib in y.index else float("nan")
        rows.append(rec)
    return pd.DataFrame(rows)


def epcam_sample_gate(log: pd.DataFrame, endpoints: dict[str, pd.Series]) -> pd.DataFrame:
    """Gates that drop a line change n. They are recorded and not used as results."""
    ep_fpkm = (2.0 ** log.loc["Epcam", ARMS]) - 1.0
    # Thresholds at each library's Epcam FPKM, plus just below the minimum.
    thresholds = sorted(set(np.round(ep_fpkm.to_numpy(dtype=float), 6)))
    rows = []
    for thr in thresholds:
        keep = [lib for lib in ARMS if float(ep_fpkm.loc[lib]) + 1e-9 >= thr]
        kl = [c for c in keep if c in KL_LIBS]
        kp = [c for c in keep if c in KP_LIBS]
        if len(kl) < 2 or len(kp) < 2:
            continue
        for name, y in endpoints.items():
            row = contrast(y, kl, kp)
            row.update(
                {
                    "endpoint": name,
                    "epcam_fpkm_min": thr,
                    "libraries_kept": ",".join(keep),
                }
            )
            rows.append(row)
    return pd.DataFrame(rows)


def permutation_leaveone(
    log: pd.DataFrame, genes: list[str]
) -> pd.DataFrame:
    """Exact 5-vs-5 label randomization.

    The gene list is fixed (label-free Epcam filter). For each assignment of
    five libraries to KL, re-select the leave-one-gene drop that maximizes
    |Welch t| and the drop that maximizes delta. 'none' is an allowed drop.
    """
    libs = list(ARMS)
    gene_vectors = {g: log.loc[g, libs].to_numpy(dtype=float) for g in genes}
    full = np.vstack([gene_vectors[g] for g in genes]).mean(axis=0)
    drops = {"none": full}
    for g in genes:
        keep = [x for x in genes if x != g]
        drops[g] = np.vstack([gene_vectors[x] for x in keep]).mean(axis=0)

    def _t_stat(vec: np.ndarray, mask: np.ndarray) -> tuple[float, float]:
        a = vec[mask]
        b = vec[~mask]
        delta = float(a.mean() - b.mean())
        t = float(stats.ttest_ind(a, b, equal_var=False).statistic)
        if not np.isfinite(t):
            t = 0.0 if delta == 0.0 else float(np.sign(delta) * np.inf)
        return delta, t

    rows = []
    for comb in combinations(range(len(libs)), 5):
        kl_idx = np.array(comb, dtype=int)
        mask = np.zeros(len(libs), dtype=bool)
        mask[kl_idx] = True
        cands = []
        for name, vec in drops.items():
            delta, t = _t_stat(vec, mask)
            cands.append({"name": name, "delta": delta, "t": t, "abs_t": abs(t)})
        fixed = next(c for c in cands if c["name"] == "none")
        # On a tie, the alphabetically last drop name is kept. Observed maxima are unique.
        win_t = max(cands, key=lambda c: (c["abs_t"], c["delta"], c["name"]))
        win_d = max(cands, key=lambda c: (c["delta"], c["abs_t"], c["name"]))
        rows.append(
            {
                "kl_libraries": ",".join(libs[i] for i in kl_idx),
                "is_observed": set(libs[i] for i in kl_idx) == set(KL_LIBS),
                "abs_t_tj_epcam": fixed["abs_t"],
                "loo_maxt_abs_t": win_t["abs_t"],
                "loo_maxt_delta": win_t["delta"],
                "loo_maxt_dropped": win_t["name"],
                "loo_maxt_t": win_t["t"],
                "loo_maxd_delta": win_d["delta"],
                "loo_maxd_abs_t": win_d["abs_t"],
                "loo_maxd_dropped": win_d["name"],
                "loo_maxd_t": win_d["t"],
            }
        )
    return pd.DataFrame(rows)


def perm_p(null: pd.Series, observed: float) -> dict:
    n = int(len(null))
    n_ge = int(np.sum(null.to_numpy(dtype=float) >= observed - 1e-12))
    return {"n_perm": n, "n_ge": n_ge, "p": n_ge / n}


def plot_libraries(scores: pd.DataFrame, path_stem: Path) -> None:
    endpoints = [
        ("Epcam", "Epcam"),
        ("Cldn4", "Cldn4"),
        ("Tacstd2", "Tacstd2"),
        ("TJ7", "TJ7"),
        ("TJ_EPCAM", "TJ Epcam r≥0.5"),
    ]
    fig, ax = plt.subplots(figsize=(8.4, 4.6))
    rng = np.random.default_rng(0)
    for i, (key, label) in enumerate(endpoints):
        for arm, color, dx in (
            ("KP", "#4C78A8", -0.14),
            ("KL", "#E45756", 0.14),
        ):
            sub = scores.loc[scores["arm"] == arm]
            jitter = rng.uniform(-0.04, 0.04, len(sub))
            ax.scatter(
                np.full(len(sub), i) + dx + jitter,
                sub[key],
                s=36,
                color=color,
                zorder=3,
                label=arm if i == 0 else None,
            )
            ax.hlines(sub[key].mean(), i + dx - 0.08, i + dx + 0.08, color=color, lw=1.6, zorder=2)
        normal = scores.loc[scores["arm"] == "normal", key]
        ax.scatter(
            [i],
            normal,
            marker="D",
            s=28,
            color="#7F7F7F",
            zorder=4,
            label="normal lung (not in test)" if i == 0 else None,
        )
    ax.set_xticks(range(len(endpoints)))
    ax.set_xticklabels([lab for _, lab in endpoints])
    ax.set_ylabel("log2(FPKM+1)")
    ax.set_title("GSE137244 cell lines, n = 5 KL vs 5 KP")
    ax.legend(frameon=False, loc="upper left")
    ax.set_xlim(-0.6, len(endpoints) - 0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path_stem.with_suffix(".png"), dpi=160)
    fig.savefig(path_stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_grid(rows: pd.DataFrame, path_stem: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.2), sharey=True)
    y = np.arange(len(rows))
    colors = ["#E45756" if flag else "#B0B0B0" for flag in rows["n5vs5"]]
    axes[0].barh(y, rows["delta"], color=colors)
    axes[1].barh(y, rows["welch_t"], color=colors)
    axes[0].axvline(0, color="black", lw=0.6)
    axes[1].axvline(0, color="black", lw=0.6)
    axes[0].set_xlabel("KL − KP  Δ  log2(FPKM+1)")
    axes[1].set_xlabel("Welch t")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(rows["label"])
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("Red bars keep n = 5 vs 5. Gray bars drop one library.", fontsize=11)
    fig.tight_layout()
    fig.savefig(path_stem.with_suffix(".png"), dpi=160)
    fig.savefig(path_stem.with_suffix(".pdf"))
    plt.close(fig)


def finding_md(ctx: dict) -> str:
    c4 = ctx["cldn4_raw"]
    ta = ctx["tacstd2_raw"]
    ep = ctx["epcam_raw"]
    tj7 = ctx["tj7"]
    tje = ctx["tj_epcam"]
    maxt = ctx["loo_maxt"]
    maxd = ctx["loo_maxd"]
    p_fixed = ctx["p_fixed"]
    p_maxt = ctx["p_maxt"]
    p_maxd = ctx["p_maxd"]
    c4_rel = ctx["cldn4_rel"]
    ta_rel = ctx["tacstd2_rel"]
    c4_res = ctx["cldn4_res"]
    ta_res = ctx["tacstd2_res"]
    drop3_c4 = ctx["drop3_cldn4"]
    drop3_ta = ctx["drop3_tacstd2"]
    return f"""# GSE137244 KL vs KP, n = 5 vs 5

Deng et al. cell-line RNA-seq (GEO GSE137244). Five KrasG12D;Lkb1 libraries (KL155, KL47, KLC, KLD, KLE) versus five KrasG12D;Trp53 libraries (B6AL10-1 through B6AL10-5). `normal-lung-RNA` (GSM4073826) is normal lung, not a tumor cell line, and is excluded. Scale is log2(FPKM+1). Every primary row uses all ten cell lines.

The locked unadjusted contrasts reproduce on this file: Tacstd2 Δ = {fmt(ta['delta'])}, Cldn4 Δ = {fmt(c4['delta'])}. Both completely separate. Exact two-sided Mann–Whitney p = 2/252 = 0.00794.

## Cldn4 and Tacstd2

Among adjustments that keep n = 5 vs 5, the unadjusted values are the joint maximum of Δ and of Welch |t|.

| endpoint | filter | Δ | Welch t | Welch p | MW p |
|---|---|---:|---:|---:|---:|
| Cldn4 | none | {fmt(c4['delta'])} | {fmt(c4['welch_t'])} | {c4['welch_p']:.3g} | {c4['mw_p']:.5f} |
| Cldn4 | minus log2(Epcam+1) | {fmt(c4_rel['delta'])} | {fmt(c4_rel['welch_t'])} | {c4_rel['welch_p']:.3g} | {c4_rel['mw_p']:.5f} |
| Cldn4 | residual on Epcam | {fmt(c4_res['delta'])} | {fmt(c4_res['welch_t'])} | {c4_res['welch_p']:.3g} | {c4_res['mw_p']:.5f} |
| Tacstd2 | none | {fmt(ta['delta'])} | {fmt(ta['welch_t'])} | {ta['welch_p']:.3g} | {ta['mw_p']:.5f} |
| Tacstd2 | minus log2(Epcam+1) | {fmt(ta_rel['delta'])} | {fmt(ta_rel['welch_t'])} | {ta_rel['welch_p']:.3g} | {ta_rel['mw_p']:.5f} |
| Tacstd2 | residual on Epcam | {fmt(ta_res['delta'])} | {fmt(ta_res['welch_t'])} | {ta_res['welch_p']:.3g} | {ta_res['mw_p']:.5f} |

Epcam itself is higher in KL (Δ = {fmt(ep['delta'])}, Welch t = {fmt(ep['welch_t'])}, p = {ep['welch_p']:.3g}), so subtracting it or residualizing on it removes part of the KL-versus-KP difference. KP B6AL10-3 has Epcam FPKM {ctx['epcam_fpkm_b6al10_3']:.2f}. Any Epcam sample gate above that value drops B6AL10-3 and leaves n = 5 vs 4. The gate that still keeps five and five does not remove any line. Those gates are in `tables/epcam_sample_gate.tsv` and are not the result.

Leave-one library keeps KL > KP for Cldn4, Tacstd2, TJ7, and the Epcam-filtered TJ mean in every drop. The largest |t| is always the drop of B6AL10-3, which is 5 vs 4: Cldn4 Δ = {fmt(drop3_c4['delta'])}, t = {fmt(drop3_c4['welch_t'])}; Tacstd2 Δ = {fmt(drop3_ta['delta'])}, t = {fmt(drop3_ta['welch_t'])}. That drop raises |t| and lowers Δ. It is not a 5-vs-5 result.

## Tight junction

TJ7 is the mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, and Ocln. On this matrix Δ = {fmt(tj7['delta'])}, Welch t = {fmt(tj7['welch_t'])}, exact MW p = {tj7['mw_p']:.5f}. The handoff value TJ +3.03 is a different average and is not recomputed here.

Epcam membership filter, τ = 0.5, on those seven genes: keep a gene if its Pearson r with Epcam across the ten lines is at least 0.5. The rule does not use the KL/KP labels. It removes Cdh1 (r = {ctx['r_cdh1']:+.3f}, Δ = {fmt(ctx['d_cdh1'])}) and F11r (r = {ctx['r_f11r']:+.3f}, Δ = {fmt(ctx['d_f11r'])}). The remaining genes are {ctx['tj_epcam_genes']}.

That five-gene mean (TJ_EPCAM) is Δ = {fmt(tje['delta'])}, Welch t = {fmt(tje['welch_t'])}, complete separation, exact MW p = {tje['mw_p']:.5f}. On the τ sweep of TJ7, τ = 0.3 selects the same genes, and no other τ in {{-1, 0, 0.3, 0.5, 0.7, 0.8}} has a larger Δ or a larger |t| while keeping at least three genes. Because the correlations are estimated on these same ten profiles, the gene list is internal to this matrix. Given the list, the label-randomization p for |Welch t| is {p_fixed['n_ge']}/{p_fixed['n_perm']} = {p_fixed['p']:.5f}.

Leave-one gene on TJ_EPCAM, with the choice made on these labels:

- Maximum |t|: drop {maxt['dropped']}. Genes {maxt['genes']}. Δ = {fmt(maxt['delta'])}, Welch t = {fmt(maxt['welch_t'])}. Still n = 5 vs 5 and completely separated. The parametric Welch p ({maxt['welch_p']:.3g}) is not adjusted for the search. The exact randomization p that re-picks the drop under every 5-vs-5 labeling is {p_maxt['n_ge']}/{p_maxt['n_perm']} = {p_maxt['p']:.5f}.
- Maximum Δ: drop {maxd['dropped']}. Genes {maxd['genes']}. Δ = {fmt(maxd['delta'])}, Welch t = {fmt(maxd['welch_t'])}. Randomization p for that maximized Δ is {p_maxd['n_ge']}/{p_maxd['n_perm']} = {p_maxd['p']:.5f}.

Both leave-one rows are above TJ7 and above TJ_EPCAM on the metric they optimize, and each also improves the other metric relative to TJ_EPCAM. They are the two Pareto rows of a six-candidate grid (the filtered mean, plus one drop for each of its five genes). They are not a second external cohort.

## What this does not change

The locked cell-line statement stays Tacstd2 Δ = +3.24, Cldn4 Δ = +5.57, n = 5 vs 5, Mann–Whitney p = 0.00794. Epcam adjustment does not raise those two deltas or their Welch |t|. The TJ increase is a filter on a pre-listed seven-gene mean, not a new experiment. No library was deleted to change n. No value was imputed. Private 8-KL matrices and TISMO LLC were not used.
"""


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    path = download_fpkm()
    digest = sha256(path)
    fpkm = load_fpkm(path)
    log = np.log2(fpkm + 1.0)

    for gene, locked in LOCKED_DELTA.items():
        got = contrast(log.loc[gene])["delta"]
        if abs(got - locked) > 1e-9:
            raise SystemExit(f"{gene} delta {got} != locked {locked}")

    ep = log.loc["Epcam"]
    cldn4 = log.loc["Cldn4"]
    tac = log.loc["Tacstd2"]
    cldn4_rel = cldn4 - ep
    tac_rel = tac - ep
    cldn4_res = residual_on(cldn4, ep, ARMS)
    tac_res = residual_on(tac, ep, ARMS)
    # Residuals were fit on the 10 lines only; put NaN on normal lung for the score table.
    cldn4_res = cldn4_res.reindex(ARMS + [NORMAL])
    tac_res = tac_res.reindex(ARMS + [NORMAL])

    gene_stats = epcam_gene_table(log, TJ7 + ["Epcam", "Tacstd2"])
    tj_genes_only = gene_stats[gene_stats["gene"].isin(TJ7)].copy()
    kept = [g for g in TJ7 if float(tj_genes_only.set_index("gene").loc[g, "r_epcam"]) >= TAU]
    if len(kept) < 3:
        raise SystemExit(f"Epcam filter kept too few genes: {kept}")

    tj7 = score_mean(log, TJ7)
    tj_epcam = score_mean(log, kept)
    loo_gene = leave_one_gene(log, kept)
    maxt = select_extreme(loo_gene, "abs_t")
    maxd = select_extreme(loo_gene, "delta")

    single_rows = []
    for endpoint, y, filt in (
        ("Cldn4", cldn4, "none"),
        ("Cldn4", cldn4_rel, "minus_log2_Epcam"),
        ("Cldn4", cldn4_res, "residual_on_Epcam"),
        ("Tacstd2", tac, "none"),
        ("Tacstd2", tac_rel, "minus_log2_Epcam"),
        ("Tacstd2", tac_res, "residual_on_Epcam"),
        ("Epcam", ep, "none"),
        ("TJ7", tj7, "none"),
        ("TJ_EPCAM", tj_epcam, f"TJ7_r_epcam_ge_{TAU}"),
    ):
        row = contrast(y.loc[ARMS] if y.isna().any() else y)
        # contrast indexes the series; residuals lack normal lung but contain ARMS.
        if set(ARMS).issubset(y.index):
            row = contrast(y)
        row.update({"endpoint": endpoint, "filter": filt, "uses_labels_to_choose": False})
        single_rows.append(row)
    singles = pd.DataFrame(single_rows)

    loo_lib_rows = []
    for name, y in (
        ("Cldn4", cldn4),
        ("Tacstd2", tac),
        ("TJ7", tj7),
        ("TJ_EPCAM", tj_epcam),
        ("Cldn4_minus_Epcam", cldn4_rel),
        ("Tacstd2_minus_Epcam", tac_rel),
    ):
        loo_lib_rows.extend(leave_one_library(y, name))
    loo_lib = pd.DataFrame(loo_lib_rows)

    gate_endpoints = {"Cldn4": cldn4, "Tacstd2": tac, "TJ7": tj7, "TJ_EPCAM": tj_epcam}
    gates = epcam_sample_gate(log, gate_endpoints)
    sweep = tau_sweep(log, tj_genes_only)

    print("Running exact 5-vs-5 label randomization...")
    perm = permutation_leaveone(log, kept)
    obs = perm.loc[perm["is_observed"]].iloc[0]
    p_fixed = perm_p(perm["abs_t_tj_epcam"], float(obs["abs_t_tj_epcam"]))
    p_maxt = perm_p(perm["loo_maxt_abs_t"], float(obs["loo_maxt_abs_t"]))
    p_maxd = perm_p(perm["loo_maxd_delta"], float(obs["loo_maxd_delta"]))

    scores = library_scores(
        log,
        {
            "Epcam": ep,
            "Cldn4": cldn4,
            "Tacstd2": tac,
            "Cldn4_minus_Epcam": cldn4_rel,
            "Tacstd2_minus_Epcam": tac_rel,
            "TJ7": tj7,
            "TJ_EPCAM": tj_epcam,
        },
    )

    # Display grid: 5-vs-5 filters plus the two leave-one-gene maxima and the
    # one leave-one-library row that maximizes |t| (disclosed as n=5 vs 4).
    display = []

    def add_display(label: str, row: dict, n5: bool | None = None) -> None:
        display.append(
            {
                "label": label,
                "delta": row["delta"],
                "welch_t": row["welch_t"],
                "n5vs5": row["n5vs5"] if n5 is None else n5,
                "n_kl": row["n_kl"],
                "n_kp": row["n_kp"],
            }
        )

    by = singles.set_index(["endpoint", "filter"])
    add_display("Cldn4 unadjusted", by.loc[("Cldn4", "none")])
    add_display("Cldn4 − Epcam", by.loc[("Cldn4", "minus_log2_Epcam")])
    add_display("Cldn4 residual on Epcam", by.loc[("Cldn4", "residual_on_Epcam")])
    add_display("Tacstd2 unadjusted", by.loc[("Tacstd2", "none")])
    add_display("Tacstd2 − Epcam", by.loc[("Tacstd2", "minus_log2_Epcam")])
    add_display("Tacstd2 residual on Epcam", by.loc[("Tacstd2", "residual_on_Epcam")])
    add_display("TJ7", by.loc[("TJ7", "none")])
    add_display("TJ Epcam r≥0.5", by.loc[("TJ_EPCAM", f"TJ7_r_epcam_ge_{TAU}")])
    add_display(f"TJ Epcam, drop {maxt['dropped']} (max |t|)", maxt)
    add_display(f"TJ Epcam, drop {maxd['dropped']} (max Δ)", maxd)
    drop3 = loo_lib[(loo_lib["endpoint"] == "Cldn4") & (loo_lib["dropped"] == "B6AL10-3-RNA")].iloc[0]
    drop3_ta = loo_lib[(loo_lib["endpoint"] == "Tacstd2") & (loo_lib["dropped"] == "B6AL10-3-RNA")].iloc[0]
    add_display("Cldn4 drop B6AL10-3 (not 5vs5)", drop3)
    add_display("Tacstd2 drop B6AL10-3 (not 5vs5)", drop3_ta)
    display_df = pd.DataFrame(display)

    samples = pd.DataFrame(
        [
            {
                "library": lib,
                "short": SHORT[lib],
                "gsm": GSM[lib],
                "arm": "KL" if lib in KL_LIBS else ("KP" if lib in KP_LIBS else "excluded_normal_lung"),
                "in_primary_n5vs5": lib in ARMS,
            }
            for lib in ARMS + [NORMAL]
        ]
    )

    gstat = tj_genes_only.set_index("gene")
    ctx = {
        "cldn4_raw": by.loc[("Cldn4", "none")].to_dict(),
        "tacstd2_raw": by.loc[("Tacstd2", "none")].to_dict(),
        "epcam_raw": by.loc[("Epcam", "none")].to_dict(),
        "cldn4_rel": by.loc[("Cldn4", "minus_log2_Epcam")].to_dict(),
        "tacstd2_rel": by.loc[("Tacstd2", "minus_log2_Epcam")].to_dict(),
        "cldn4_res": by.loc[("Cldn4", "residual_on_Epcam")].to_dict(),
        "tacstd2_res": by.loc[("Tacstd2", "residual_on_Epcam")].to_dict(),
        "tj7": by.loc[("TJ7", "none")].to_dict(),
        "tj_epcam": by.loc[("TJ_EPCAM", f"TJ7_r_epcam_ge_{TAU}")].to_dict(),
        "loo_maxt": maxt.to_dict(),
        "loo_maxd": maxd.to_dict(),
        "p_fixed": p_fixed,
        "p_maxt": p_maxt,
        "p_maxd": p_maxd,
        "drop3_cldn4": drop3.to_dict(),
        "drop3_tacstd2": drop3_ta.to_dict(),
        "epcam_fpkm_b6al10_3": float(fpkm.loc["Epcam", "B6AL10-3-RNA"]),
        "r_cdh1": float(gstat.loc["Cdh1", "r_epcam"]),
        "r_f11r": float(gstat.loc["F11r", "r_epcam"]),
        "d_cdh1": float(gstat.loc["Cdh1", "delta"]),
        "d_f11r": float(gstat.loc["F11r", "delta"]),
        "tj_epcam_genes": ", ".join(kept),
    }

    # Joint-max checks the prose depends on.
    if not (
        ctx["cldn4_raw"]["abs_t"] > ctx["cldn4_rel"]["abs_t"]
        and ctx["cldn4_raw"]["abs_t"] > ctx["cldn4_res"]["abs_t"]
        and ctx["cldn4_raw"]["delta"] > ctx["cldn4_rel"]["delta"]
        and ctx["cldn4_raw"]["delta"] > ctx["cldn4_res"]["delta"]
    ):
        raise SystemExit("Cldn4 unadjusted is not the joint maximum; update the finding rule")
    if not (
        ctx["tacstd2_raw"]["abs_t"] > ctx["tacstd2_rel"]["abs_t"]
        and ctx["tacstd2_raw"]["abs_t"] > ctx["tacstd2_res"]["abs_t"]
        and ctx["tacstd2_raw"]["delta"] > ctx["tacstd2_rel"]["delta"]
        and ctx["tacstd2_raw"]["delta"] > ctx["tacstd2_res"]["delta"]
    ):
        raise SystemExit("Tacstd2 unadjusted is not the joint maximum; update the finding rule")
    if ctx["loo_maxt"]["dropped"] == "none" or ctx["loo_maxd"]["dropped"] == "none":
        raise SystemExit("Leave-one did not beat the filtered mean; update the finding")
    tje_delta = float(ctx["tj_epcam"]["delta"])
    tje_t = float(ctx["tj_epcam"]["abs_t"])
    if not (
        float(ctx["loo_maxt"]["delta"]) > tje_delta
        and float(ctx["loo_maxt"]["abs_t"]) > tje_t
        and float(ctx["loo_maxd"]["delta"]) > tje_delta
        and float(ctx["loo_maxd"]["abs_t"]) > tje_t
    ):
        raise SystemExit("A leave-one row does not improve both metrics; update the finding")
    if float(ctx["drop3_cldn4"]["delta"]) >= float(ctx["cldn4_raw"]["delta"]):
        raise SystemExit("Dropping B6AL10-3 did not lower the Cldn4 delta; update the finding")
    if float(ctx["drop3_tacstd2"]["delta"]) >= float(ctx["tacstd2_raw"]["delta"]):
        raise SystemExit("Dropping B6AL10-3 did not lower the Tacstd2 delta; update the finding")
    best_sweep = sweep.sort_values(["abs_t", "delta"], ascending=False).iloc[0]
    if abs(float(best_sweep["tau"]) - TAU) > 1e-9 and set(str(best_sweep["genes"]).split(",")) != set(kept):
        raise SystemExit("tau=0.5 is not the sweep maximum; update the finding")

    write_tsv(samples, "samples.tsv")
    write_tsv(scores, "library_scores.tsv")
    write_tsv(singles, "endpoint_filters.tsv")
    write_tsv(tj_genes_only, "tj7_gene_epcam.tsv")
    write_tsv(sweep, "tj7_tau_sweep.tsv")
    write_tsv(loo_gene, "tj_epcam_leaveone_gene.tsv")
    write_tsv(loo_lib, "leaveone_library.tsv")
    write_tsv(gates, "epcam_sample_gate.tsv")
    write_tsv(perm, "permutation_labelings.tsv")
    write_tsv(display_df, "display_grid.tsv")

    summary = {
        "accession": "GSE137244",
        "file": FPKM_NAME,
        "sha256": digest,
        "scale": "log2(FPKM+1)",
        "n_kl": 5,
        "n_kp": 5,
        "excluded": NORMAL,
        "tj7_genes": TJ7,
        "tau": TAU,
        "tj_epcam_genes": kept,
        "locked_delta_reproduced": LOCKED_DELTA,
        "contrasts": {
            "Cldn4": ctx["cldn4_raw"],
            "Tacstd2": ctx["tacstd2_raw"],
            "Epcam": ctx["epcam_raw"],
            "TJ7": ctx["tj7"],
            "TJ_EPCAM": ctx["tj_epcam"],
            "TJ_EPCAM_loo_max_abs_t": ctx["loo_maxt"],
            "TJ_EPCAM_loo_max_delta": ctx["loo_maxd"],
        },
        "randomization": {
            "abs_t_tj_epcam": p_fixed,
            "loo_max_abs_t": p_maxt,
            "loo_max_delta": p_maxd,
            "observed_loo_maxt_dropped": obs["loo_maxt_dropped"],
            "observed_loo_maxd_dropped": obs["loo_maxd_dropped"],
        },
    }
    # Series are not JSON-safe if any numpy types remain inside the contrast dicts.
    def convert(obj):
        if isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [convert(v) for v in obj]
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.bool_,)):
            return bool(obj)
        return obj

    summary = convert(summary)
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (ROOT / "FINDING.md").write_text(finding_md(convert(ctx)))

    plot_libraries(scores, FIGURES / "library_scores")
    plot_grid(display_df, FIGURES / "filter_grid")

    print(json.dumps(summary["randomization"], indent=2))
    print("TJ_EPCAM genes", kept)
    print("max |t| drop", maxt["dropped"], float(maxt["delta"]), float(maxt["welch_t"]))
    print("max delta drop", maxd["dropped"], float(maxd["delta"]), float(maxd["welch_t"]))
    print("sha256", digest)
    print("wrote", TABLES, FIGURES)


if __name__ == "__main__":
    main()
