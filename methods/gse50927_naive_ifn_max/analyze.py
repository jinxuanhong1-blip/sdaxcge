#!/usr/bin/env python3
"""GSE50927 naive lung: max-effect IFN/ISG panel on the deposited edgeR column.

Contrast is Cldn4 KO no VILI (GSM1232581) minus WT no VILI (GSM1232580).
The public processed table has one logFC per gene. n = 1 vs 1.
Design text says duplicates; those duplicates are not columns in this table.
This script does not promote n.

Panels are locked lists already used in this repository, plus the MSigDB
mouse Hallmark IFN-α / IFN-γ sets (and their intersection and union).
No gene is added or dropped because of its logFC.

Selection, applied only to the unfiltered detected-gene rows:
  * max joint = maximum of (mean logFC × fraction of genes with logFC > 0)
  * max joint module = same, restricted to panels with >= 20 detected genes
  * max concordance = maximum fraction up, tie-broken by higher mean logFC

A logCPM > 0 row is a sensitivity on the same gene list. It is not allowed
to define a new panel.
"""

from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
SETDIR = HERE / "genesets"
TAB = HERE / "tables"
FIG = HERE / "figures"
DATA = Path(os.environ.get("GSE50927_NAIVE_IFN_DATA", "/tmp/gse50927_naive_ifn"))
URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/"
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
)
MODULE_MIN_N = 20
UNIT = "genes on one KO-vs-WT edgeR column; not mice"

plt.rcParams.update({
    "font.size": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def log(msg: str) -> None:
    print(msg, flush=True)


def fetch_table() -> Path:
    dest = DATA / "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    DATA.mkdir(parents=True, exist_ok=True)
    log(f"download {URL}")
    urllib.request.urlretrieve(URL, dest)
    return dest


def load_edger(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={"Marker.Symbol": "symbol"}).dropna(subset=["symbol"])
    df["symbol"] = df["symbol"].astype(str)
    for c in ("logFC", "logCPM", "PValue", "FDR"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    # Same rule as the earlier GSE50927 scripts: keep the higher-abundance row.
    df = df.sort_values("logCPM", ascending=False).drop_duplicates("symbol")
    return df.set_index("symbol", drop=False)


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def default_mouse(human: str) -> list[str]:
    if "-" in human:
        return [human[0] + human[1:].lower()]
    return [human.capitalize()]


def first_present(candidates: list[str], index: pd.Index, used: set[str]) -> str | None:
    for sym in candidates:
        if sym in index and sym not in used:
            return sym
    return None


def map_human(genes: list[str], overrides: dict[str, list[str]], index: pd.Index) -> list[str]:
    """First listed ortholog that is on the table. One mouse gene, once."""
    out: list[str] = []
    used: set[str] = set()
    for human in genes:
        cands = overrides[human] if human in overrides else default_mouse(human)
        hit = first_present(cands, index, used)
        if hit is None:
            continue
        used.add(hit)
        out.append(hit)
    return out


def pairwise_concordance(n_up: int, n_down: int, n_zero: int) -> float:
    n = n_up + n_down + n_zero
    if n < 2:
        return float("nan")
    # Same-sign pairs. Zeros agree with nobody, including other zeros:
    # a zero is not a directional call.
    agree = n_up * (n_up - 1) / 2 + n_down * (n_down - 1) / 2
    total = n * (n - 1) / 2
    return agree / total


def score_panel(df: pd.DataFrame, genes: list[str], min_logcpm: float | None) -> dict:
    """min_logcpm is an exclusive floor: keep genes with logCPM > floor."""
    listed = list(dict.fromkeys(genes))
    detected, missing, low = [], [], []
    for g in listed:
        if g not in df.index:
            missing.append(g)
            continue
        if min_logcpm is not None and not (float(df.at[g, "logCPM"]) > min_logcpm):
            low.append(g)
            continue
        detected.append(g)
    if not detected:
        raise RuntimeError("panel resolved to zero genes")
    sub = df.loc[detected, ["logFC", "logCPM", "PValue", "FDR"]].copy()
    fc = sub["logFC"].astype(float)
    n = int(len(fc))
    n_up = int((fc > 0).sum())
    n_down = int((fc < 0).sum())
    n_zero = n - n_up - n_down
    bg = df.loc[~df.index.isin(detected), "logFC"].astype(float).dropna()
    nz = fc[fc != 0]
    wilcoxon_p = float(stats.wilcoxon(nz, alternative="greater").pvalue) if len(nz) else float("nan")
    sign_p = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    mwu_p = float(stats.mannwhitneyu(fc, bg, alternative="greater").pvalue)
    mean = float(fc.mean())
    conc = n_up / n
    return {
        "n_listed": len(listed),
        "n_detected": n,
        "n_missing": len(missing),
        "n_low_cpm_excluded": len(low),
        "mean_logFC": mean,
        "median_logFC": float(fc.median()),
        "n_up": n_up,
        "n_down": n_down,
        "n_zero": n_zero,
        "concordance_frac_up": conc,
        "pairwise_concordance": pairwise_concordance(n_up, n_down, n_zero),
        "joint_mean_x_concordance": mean * conc,
        "sign_test_p_up_genes": sign_p,
        "wilcoxon_p_up_genes": wilcoxon_p,
        "mwu_p_panel_gt_bg_genes": mwu_p,
        "genes_up": ",".join(fc[fc > 0].sort_values(ascending=False).index),
        "genes_down": ",".join(fc[fc < 0].sort_values().index),
        "genes_low_cpm_excluded": ",".join(low),
        "genes_missing": ",".join(missing),
        "detected_genes": detected,
        "fc": fc,
        "sub": sub,
    }


def gene_table(panel: str, filt: str, scored: dict) -> pd.DataFrame:
    sub = scored["sub"].copy()
    sub.insert(0, "panel", panel)
    sub.insert(1, "filter", filt)
    sub.insert(2, "symbol", sub.index)
    sub["direction"] = np.where(sub["logFC"] > 0, "UP", np.where(sub["logFC"] < 0, "DOWN", "ZERO"))
    sub["low_count"] = sub["logCPM"] <= 0
    return sub.reset_index(drop=True)


def pick(rows: list[dict], key: str, eligible) -> dict:
    pool = [r for r in rows if eligible(r)]
    if not pool:
        raise RuntimeError(f"no panel eligible for {key}")
    return max(pool, key=lambda r: (r[key], r["mean_logFC"], r["concordance_frac_up"], -r["n_detected"]))


def pareto(rows: list[dict]) -> list[dict]:
    front = []
    for r in rows:
        dominated = False
        for o in rows:
            if o is r:
                continue
            ge_mean = o["mean_logFC"] >= r["mean_logFC"]
            ge_conc = o["concordance_frac_up"] >= r["concordance_frac_up"]
            strict = o["mean_logFC"] > r["mean_logFC"] or o["concordance_frac_up"] > r["concordance_frac_up"]
            if ge_mean and ge_conc and strict:
                dominated = True
                break
        if not dominated:
            front.append(r)
    front.sort(key=lambda r: -r["mean_logFC"])
    return front


def draw_lollipop(ax, scored: dict, title: str) -> None:
    fc = scored["fc"].sort_values()
    colors, markers = [], []
    y = np.arange(len(fc))
    for sym, val in fc.items():
        low = float(scored["sub"].at[sym, "logCPM"]) <= 0
        colors.append("#2166ac" if val < 0 else "#b2182b")
        markers.append("o" if not low else "D")
    ax.axvline(0, color="#333333", lw=0.6)
    ax.axvline(scored["mean_logFC"], color="#b2182b", lw=1.0, ls="--",
               label=f"mean {scored['mean_logFC']:+.3f}")
    for i, (sym, val) in enumerate(fc.items()):
        ax.plot([0, val], [i, i], color=colors[i], lw=1.4, solid_capstyle="round")
        ax.scatter([val], [i], s=28 if markers[i] == "o" else 36, c=colors[i],
                   marker=markers[i], zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(list(fc.index), fontsize=8)
    ax.set_xlabel("author edgeR logFC (Cldn4 KO − WT)")
    ax.set_title(title, loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="lower right")


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    df = load_edger(fetch_table())
    locked = json.loads((SETDIR / "locked_panels.json").read_text())
    gmt = read_gmt(SETDIR / "mh.all.v2023.2.Mm.symbols.gmt")
    ifna = gmt["HALLMARK_INTERFERON_ALPHA_RESPONSE"]
    ifng = gmt["HALLMARK_INTERFERON_GAMMA_RESPONSE"]
    fable_ext = locked["fable_isg_human"] + locked["fable_mhc_apm_human"]

    mouse_panels: dict[str, tuple[str, list[str]]] = {
        "hallmark_ifn_alpha": ("hallmark", ifna),
        "hallmark_ifn_gamma": ("hallmark", ifng),
        "hallmark_ifn_intersection": ("hallmark", sorted(set(ifna) & set(ifng))),
        "hallmark_ifn_union": ("hallmark", sorted(set(ifna) | set(ifng))),
        "compact26": ("compact", locked["compact26_mouse"]),
        "fable_isg": ("fable", map_human(locked["fable_isg_human"], locked["fable_mouse_first_hit"], df.index)),
        "fable_extended": ("fable", map_human(fable_ext, locked["fable_mouse_first_hit"], df.index)),
        "fable_core8": ("fable", map_human(locked["fable_core8_human"], locked["fable_mouse_first_hit"], df.index)),
        "ifn_alpha_type1": ("kdko", map_human(locked["ifn_alpha_type1_human"], locked["kdko_mouse_exceptions"], df.index)),
        "ifn_gamma_type2": ("kdko", map_human(locked["ifn_gamma_type2_human"], locked["kdko_mouse_exceptions"], df.index)),
        "ayers_ifng6": ("ayers", map_human(locked["ayers_ifng6_human"], locked["kdko_mouse_exceptions"], df.index)),
    }

    filters = (("detected", None), ("logCPM_gt_0", 0.0), ("logCPM_gt_1", 1.0))
    rows: list[dict] = []
    scored_cache: dict[tuple[str, str], dict] = {}
    for name, (family, genes) in mouse_panels.items():
        for filt, floor in filters:
            scored = score_panel(df, genes, floor)
            scored_cache[(name, filt)] = scored
            rec = {
                "panel": name,
                "family": family,
                "filter": filt,
                "n_listed_resolved": len(list(dict.fromkeys(genes))),
                "n_detected": scored["n_detected"],
                "n_missing_from_table": scored["n_missing"] if floor is None else scored["n_missing"],
                "n_low_cpm_excluded": scored["n_low_cpm_excluded"],
                "mean_logFC_KO_minus_WT": scored["mean_logFC"],
                "median_logFC_KO_minus_WT": scored["median_logFC"],
                "n_up": scored["n_up"],
                "n_down": scored["n_down"],
                "concordance_frac_up": scored["concordance_frac_up"],
                "pairwise_concordance": scored["pairwise_concordance"],
                "joint_mean_x_concordance": scored["joint_mean_x_concordance"],
                "sign_test_p_up_genes": scored["sign_test_p_up_genes"],
                "wilcoxon_p_up_genes": scored["wilcoxon_p_up_genes"],
                "mwu_p_panel_gt_background_genes": scored["mwu_p_panel_gt_bg_genes"],
                "inference_unit": UNIT,
                "genes_down": scored["genes_down"],
                "genes_low_cpm_excluded": scored["genes_low_cpm_excluded"],
                "genes_missing": scored["genes_missing"],
                "mean_logFC": scored["mean_logFC"],
            }
            rows.append(rec)

    detected_rows = [r for r in rows if r["filter"] == "detected"]
    winners = {
        "max_joint": pick(detected_rows, "joint_mean_x_concordance", lambda r: True),
        "max_joint_module_n_ge_20": pick(
            detected_rows, "joint_mean_x_concordance", lambda r: r["n_detected"] >= MODULE_MIN_N
        ),
        "max_concordance": pick(detected_rows, "concordance_frac_up", lambda r: True),
        "max_mean": pick(detected_rows, "mean_logFC", lambda r: True),
        # Concordance gate is the gene-level sign test, not a cutoff placed
        # between two fractions. A 5/6 panel can lead the mean and still fail it.
        "max_mean_sign_concordant": pick(
            detected_rows, "mean_logFC", lambda r: r["sign_test_p_up_genes"] < 0.05
        ),
    }
    front = pareto(detected_rows)

    # Locked biological checks. These fail the run if the wrong table is loaded.
    cldn4 = df.loc["Cldn4"]
    compact = scored_cache[("compact26", "detected")]
    ayers = scored_cache[("ayers_ifng6", "detected")]
    assert abs(float(cldn4["logFC"]) - (-6.061296179)) < 1e-6
    assert compact["n_detected"] == 26 and compact["n_up"] == 23
    assert abs(compact["mean_logFC"] - 1.0237861997692308) < 1e-9
    assert ayers["n_detected"] == 6 and ayers["n_up"] == 5
    assert winners["max_joint"]["panel"] == "ayers_ifng6"
    assert winners["max_joint_module_n_ge_20"]["panel"] == "compact26"
    assert winners["max_concordance"]["panel"] == "fable_isg"
    assert winners["max_mean_sign_concordant"]["panel"] == "compact26"

    score_df = pd.DataFrame(rows).drop(columns=["mean_logFC"])
    score_df.to_csv(TAB / "panel_scores.tsv", sep="\t", index=False)

    membership_rows = []
    for name, (_family, genes) in mouse_panels.items():
        for i, g in enumerate(list(dict.fromkeys(genes))):
            on = g in df.index
            membership_rows.append({
                "panel": name,
                "order": i,
                "symbol": g,
                "on_table": on,
                "logFC": float(df.at[g, "logFC"]) if on else np.nan,
                "logCPM": float(df.at[g, "logCPM"]) if on else np.nan,
                "author_FDR": float(df.at[g, "FDR"]) if on else np.nan,
                "author_PValue": float(df.at[g, "PValue"]) if on else np.nan,
            })
    pd.DataFrame(membership_rows).to_csv(TAB / "panel_membership.tsv", sep="\t", index=False)

    headline_names = []
    for key in ("max_joint", "max_joint_module_n_ge_20", "max_concordance"):
        name = winners[key]["panel"]
        if name not in headline_names:
            headline_names.append(name)
    gene_frames = []
    for name in headline_names:
        for filt in ("detected", "logCPM_gt_0"):
            gene_frames.append(gene_table(name, filt, scored_cache[(name, filt)]))
    pd.concat(gene_frames, ignore_index=True).to_csv(TAB / "headline_genes.tsv", sep="\t", index=False)

    def brief(r: dict) -> dict:
        return {
            "panel": r["panel"],
            "n_detected": r["n_detected"],
            "n_up": r["n_up"],
            "n_down": r["n_down"],
            "mean_logFC_KO_minus_WT": r["mean_logFC_KO_minus_WT"],
            "concordance_frac_up": r["concordance_frac_up"],
            "joint_mean_x_concordance": r["joint_mean_x_concordance"],
            "genes_down": r["genes_down"],
        }

    one = []
    for label, r in winners.items():
        one.append({"role": label, **brief(r), "n_mice": "1 vs 1", "inference_unit": UNIT})
    for r in front:
        one.append({"role": "pareto_detected", **brief(r), "n_mice": "1 vs 1", "inference_unit": UNIT})
    pd.DataFrame(one).drop_duplicates().to_csv(TAB / "one_row.tsv", sep="\t", index=False)

    summary = {
        "series": "GSE50927",
        "contrast": "naive whole lung, Cldn4 KO minus WT",
        "gsm_ko": "GSM1232581",
        "gsm_wt": "GSM1232580",
        "n": "1 vs 1",
        "n_note": "One deposited edgeR column. Design-text duplicates are not in this table and are not used.",
        "logFC_sign": "author edgeR logFC is KO minus WT; Cldn4 = -6.061296179 confirms the sign",
        "cldn4_logFC": float(cldn4["logFC"]),
        "cldn4_logCPM": float(cldn4["logCPM"]),
        "cldn4_author_FDR": float(cldn4["FDR"]),
        "selection_rule": (
            "Locked panels only. No gene added or removed for its logFC. "
            "max_joint maximizes mean_logFC * (n_up/n). "
            "max_joint_module applies that rule to panels with at least 20 detected genes. "
            "max_concordance maximizes n_up/n, tie-broken by mean logFC. "
            "max_mean_sign_concordant maximizes mean logFC among panels whose "
            "gene-level sign test (fraction up > 1/2) has p < 0.05."
        ),
        "inference_unit": UNIT,
        "winners": {k: brief(v) for k, v in winners.items()},
        "pareto_detected": [brief(r) for r in front],
        "compact26_logCPM_gt_0": brief(next(r for r in rows if r["panel"] == "compact26" and r["filter"] == "logCPM_gt_0")),
        "ayers_ifng6_logCPM_gt_0": brief(next(r for r in rows if r["panel"] == "ayers_ifng6" and r["filter"] == "logCPM_gt_0")),
    }
    (TAB / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Figure 1 — the two non-dominated headlines that answer the request.
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 7.2), gridspec_kw={"width_ratios": [1.15, 2.35]})
    a = scored_cache[("ayers_ifng6", "detected")]
    c = scored_cache[("compact26", "detected")]
    draw_lollipop(
        axes[0], a,
        "Ayers IFN-γ 6 — largest mean\n"
        "mean {m:+.3f} · {u}/{n} up · gene sign p = 0.11".format(
            m=a["mean_logFC"], u=a["n_up"], n=a["n_detected"]),
    )
    draw_lollipop(
        axes[1], c,
        "IFN compact 26 — largest mean among sign-concordant panels\n"
        "mean {m:+.3f} · {u}/{n} up · gene sign p = 4.4×10$^{{-5}}$".format(
            m=c["mean_logFC"], u=c["n_up"], n=c["n_detected"]),
    )
    fig.suptitle(
        "GSE50927 naive lung  ·  Cldn4 KO − WT  ·  n = 1 vs 1\n"
        "Diamonds are logCPM ≤ 0. Gene concordance is not a mouse-level test.",
        fontsize=11, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(FIG / "fig1_max_effect_lollipop.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "fig1_max_effect_lollipop.pdf", bbox_inches="tight")
    plt.close(fig)

    # Figure 2 — ranked mean, with concordance written on the bar.
    # Pareto panels stay red. Everything a higher-mean panel also beats on
    # concordance is gray, so the tradeoff is readable without stacked labels.
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    order = sorted(detected_rows, key=lambda r: r["mean_logFC_KO_minus_WT"])
    front_names = {r["panel"] for r in front}
    y = np.arange(len(order))
    ax.barh(
        y,
        [r["mean_logFC_KO_minus_WT"] for r in order],
        color=["#b2182b" if r["panel"] in front_names else "#d0d0d0" for r in order],
        height=0.72,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([r["panel"] for r in order], fontsize=8)
    for i, r in enumerate(order):
        ax.text(
            r["mean_logFC_KO_minus_WT"] + 0.03, i,
            "{u}/{n} up   joint {j:.2f}".format(
                u=r["n_up"], n=r["n_detected"], j=r["joint_mean_x_concordance"]),
            va="center", ha="left", fontsize=7.5, color="#222222",
        )
    ax.set_xlim(0, 2.35)
    ax.set_xlabel("mean logFC (Cldn4 KO − WT)")
    ax.set_title("Locked IFN/ISG panels · naive lung · n = 1 vs 1\nRed = Pareto front of mean logFC and fraction up")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_panel_pareto.png", dpi=160)
    fig.savefig(FIG / "fig2_panel_pareto.pdf")
    plt.close(fig)

    log("winners")
    for k, r in winners.items():
        log(f"  {k}: {r['panel']} mean={r['mean_logFC_KO_minus_WT']:+.4f} "
            f"up={r['n_up']}/{r['n_detected']} joint={r['joint_mean_x_concordance']:.4f}")
    log("pareto: " + ", ".join(r["panel"] for r in front))
    log(f"wrote {TAB} and {FIG}")


if __name__ == "__main__":
    main()
