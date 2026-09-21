#!/usr/bin/env python3
"""GSE274940: EpH4 wild type vs Cldn-null bulk RNA-seq.

The series is the RNA-seq source data for Kashihara et al., Science Advances
2025 (PMID 41171911, doi:10.1126/sciadv.adx7431). The genotype is a multiplex
CRISPR knockout of the eight claudins EpH4 expresses (Cldn3, Cldn4, Cldn7,
Cldn8, Cldn9, Cldn12, Cldn23, Cldn25). It is not a CLDN4-only knockout.

Question this script answers from the deposited expected-count matrix:
  1. Is Cldn4 mRNA actually lower in the Cldn-null libraries?
  2. What happens to IFN and tight-junction genes in the same contrast?

Statistics match the CLDN4/TACSTD2 KD-KO wave: log2(median-of-ratios
normalized counts + 1), Welch t-test across libraries, Benjamini-Hochberg
across genes, plus a sample-level gene-set score and a competitive
Mann-Whitney on t statistics. n = 3 vs 3. With that n the two-sided
Mann-Whitney p-value cannot fall below 0.1, so Welch is the primary test.
"""
from __future__ import annotations

import hashlib
import os
import urllib.request

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data", "gse274940")
RES = os.path.join(ROOT, "results", "gse274940")
COUNT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274940/suppl/"
    "GSE274940_raw_counts.csv.gz"
)
COUNT_NAME = "GSE274940_raw_counts.csv.gz"

SAMPLES = ["WT1", "WT2", "WT3", "KO1", "KO2", "KO3"]
WT = ["WT1", "WT2", "WT3"]
KO = ["KO1", "KO2", "KO3"]
# GEO sample titles EpH4 WT_1..3 and EpH4 Cldn-null_1..3, in accession order.
GSM = {
    "WT1": "GSM8462258",
    "WT2": "GSM8462259",
    "WT3": "GSM8462260",
    "KO1": "GSM8462261",
    "KO2": "GSM8462262",
    "KO3": "GSM8462263",
}
# Endogenous EpH4 claudins cut by the multiplex CRISPR (paper Fig. 1C).
KO_TARGETS = ["Cldn3", "Cldn4", "Cldn7", "Cldn8", "Cldn9", "Cldn12", "Cldn23", "Cldn25"]
EXPR_MIN = 10.0  # max(group mean raw count) required to enter the Welch screen


def bh_fdr(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.sum() == 0:
        return out
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order]
    n = len(ranked)
    adj = ranked * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(n)
    res[order] = adj
    out[ok] = res
    return out


def download_counts() -> str:
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, COUNT_NAME)
    if not os.path.exists(path):
        urllib.request.urlretrieve(COUNT_URL, path)
    return path


def load_counts(path: str) -> tuple[pd.DataFrame, pd.Series]:
    raw = pd.read_csv(path)
    missing = [c for c in ["ALIAS", *SAMPLES] if c not in raw.columns]
    if missing:
        raise SystemExit(f"count matrix missing columns: {missing}")
    raw = raw.dropna(subset=["ALIAS"]).copy()
    raw["ALIAS"] = raw["ALIAS"].astype(str)
    n_ids = raw.groupby("ALIAS")["ENSEMBL"].nunique()
    counts = raw.groupby("ALIAS")[SAMPLES].sum().astype(float)
    if counts.index.duplicated().any():
        raise SystemExit("ALIAS collapse left duplicate symbols")
    return counts, n_ids


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq2 median-of-ratios size factors (Anders & Huber 2010)."""
    sub = counts.loc[(counts > 0).all(axis=1)]
    log_geo = np.log(sub.to_numpy()).mean(axis=1, keepdims=True)
    log_sf = np.median(np.log(sub.to_numpy()) - log_geo, axis=0)
    return pd.Series(np.exp(log_sf), index=counts.columns, name="size_factor")


def welch_ci(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float, float]:
    """Mean(a)-mean(b), Welch df, and 95% CI for that difference."""
    na, nb = len(a), len(b)
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    diff = float(np.mean(a) - np.mean(b))
    se2 = va / na + vb / nb
    if se2 <= 0:
        return diff, np.nan, diff, diff
    df = se2 ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    tcrit = float(stats.t.ppf(0.975, df))
    se = np.sqrt(se2)
    return diff, float(df), diff - tcrit * se, diff + tcrit * se


def gene_tests(logm: pd.DataFrame, tested: pd.Index) -> pd.DataFrame:
    a = logm.loc[tested, KO].to_numpy(dtype=float)
    b = logm.loc[tested, WT].to_numpy(dtype=float)
    t_stat, p = stats.ttest_ind(a, b, axis=1, equal_var=False, nan_policy="omit")
    rows = []
    for i, gene in enumerate(tested):
        diff, df, lo, hi = welch_ci(a[i], b[i])
        ai, bi = a[i], b[i]
        if np.ptp(np.r_[ai, bi]) > 0:
            mw = float(stats.mannwhitneyu(ai, bi, alternative="two-sided").pvalue)
        else:
            mw = np.nan
        sdp = np.sqrt((np.var(ai, ddof=1) + np.var(bi, ddof=1)) / 2)
        rows.append({
            "gene": gene,
            "log2FC": diff,
            "ci95_low": lo,
            "ci95_high": hi,
            "df_welch": df,
            "cohens_d": float(diff / sdp) if sdp > 0 else np.nan,
            "t_stat": float(t_stat[i]),
            "p_welch": float(p[i]),
            "p_mannwhitney": mw,
        })
    out = pd.DataFrame(rows).set_index("gene")
    out["q_welch_BH"] = bh_fdr(out["p_welch"].to_numpy())
    return out


def load_sets() -> dict[str, list[str]]:
    path = os.path.join(HERE, "gene_sets_mouse.tsv")
    sets = {}
    for row in pd.read_csv(path, sep="\t").itertuples(index=False):
        genes = [g.strip() for g in str(row.mouse_symbols).split(";") if g.strip()]
        sets[row.set_name] = genes
    return sets


def set_tests(logm: pd.DataFrame, gstats: pd.DataFrame,
              sets: dict[str, list[str]]) -> pd.DataFrame:
    all_t = gstats["t_stat"]
    rows = []
    for name, genes in sets.items():
        present = [g for g in genes if g in logm.index]
        tested = [g for g in present if g in gstats.index]
        rec: dict[str, object] = {
            "set_name": name,
            "n_genes_listed": len(genes),
            "n_genes_in_matrix": len(present),
            "n_genes_tested": len(tested),
            "genes_missing": ";".join(g for g in genes if g not in logm.index),
            "n_WT": 3,
            "n_KO": 3,
        }
        if len(tested) >= 3:
            sub = logm.loc[tested]
            sd = sub.std(axis=1, ddof=1).replace(0, np.nan)
            z = sub.sub(sub.mean(axis=1), axis=0).div(sd, axis=0)
            score = z.mean(axis=0, skipna=True)
            sa, sb = score[KO].to_numpy(), score[WT].to_numpy()
            tt, pp = stats.ttest_ind(sa, sb, equal_var=False)
            sdp = np.sqrt((np.var(sa, ddof=1) + np.var(sb, ddof=1)) / 2)
            rec.update({
                "sample_score_KO": float(np.mean(sa)),
                "sample_score_WT": float(np.mean(sb)),
                "sample_delta_KO_minus_WT": float(np.mean(sa) - np.mean(sb)),
                "sample_t": float(tt),
                "sample_p": float(pp),
                "sample_cohens_d": float((np.mean(sa) - np.mean(sb)) / sdp) if sdp else np.nan,
            })
        else:
            rec.update({k: np.nan for k in (
                "sample_score_KO", "sample_score_WT", "sample_delta_KO_minus_WT",
                "sample_t", "sample_p", "sample_cohens_d")})
        t_in = all_t.reindex(tested).dropna()
        t_out = all_t.drop(index=t_in.index, errors="ignore").dropna()
        if len(t_in) >= 3 and len(t_out) >= 20:
            u = stats.mannwhitneyu(t_in, t_out, alternative="two-sided")
            rec.update({
                "comp_median_t_in_set": float(np.median(t_in)),
                "comp_median_t_background": float(np.median(t_out)),
                "comp_mannwhitney_p": float(u.pvalue),
                "comp_median_log2FC": float(gstats.loc[t_in.index, "log2FC"].median()),
                "n_log2FC_positive": int((gstats.loc[t_in.index, "log2FC"] > 0).sum()),
                "n_log2FC_negative": int((gstats.loc[t_in.index, "log2FC"] < 0).sum()),
                "n_welch_p_lt_0.05": int((gstats.loc[t_in.index, "p_welch"] < 0.05).sum()),
            })
        else:
            rec.update({k: np.nan for k in (
                "comp_median_t_in_set", "comp_median_t_background",
                "comp_mannwhitney_p", "comp_median_log2FC",
                "n_log2FC_positive", "n_log2FC_negative", "n_welch_p_lt_0.05")})
        rows.append(rec)
    df = pd.DataFrame(rows)
    df["sample_q_BH"] = bh_fdr(df["sample_p"].to_numpy())
    df["comp_q_BH"] = bh_fdr(df["comp_mannwhitney_p"].to_numpy())
    return df


def mRNA_call(mean_wt: float, mean_ko: float, pct: float, p: float) -> str:
    if mean_wt < 1 and mean_ko < 1:
        return "undetectable_both"
    if mean_wt < 1 and mean_ko >= 10:
        return "induced_from_undetectable"
    if mean_wt >= 10 and pct >= 80:
        return "not_depleted"
    if mean_wt >= 10 and pct < 50 and p < 0.05:
        return "depleted"
    if mean_wt >= 10 and pct < 80 and p < 0.05:
        return "reduced"
    if mean_wt >= 10 and pct < 80:
        return "lower_point_estimate_not_significant"
    return "see_numbers"


def validity_table(counts: pd.DataFrame, cpm: pd.DataFrame,
                   gstats: pd.DataFrame) -> pd.DataFrame:
    symbols = sorted(
        {g for g in counts.index if str(g).startswith("Cldn") and str(g)[4:].isdigit()},
        key=lambda g: int(g[4:]),
    )
    rows = []
    for gene in symbols:
        mean_wt = float(counts.loc[gene, WT].mean())
        mean_ko = float(counts.loc[gene, KO].mean())
        cpm_wt = float(cpm.loc[gene, WT].mean())
        cpm_ko = float(cpm.loc[gene, KO].mean())
        pct = float(100.0 * cpm_ko / cpm_wt) if cpm_wt > 0 else np.nan
        if gene in gstats.index:
            st = gstats.loc[gene]
            log2fc, p, q = float(st.log2FC), float(st.p_welch), float(st.q_welch_BH)
            lo, hi = float(st.ci95_low), float(st.ci95_high)
        else:
            log2fc = p = q = lo = hi = np.nan
        rows.append({
            "gene": gene,
            "crispr_target": gene in KO_TARGETS,
            "mean_raw_WT": mean_wt,
            "mean_raw_KO": mean_ko,
            "mean_CPM_WT": cpm_wt,
            "mean_CPM_KO": cpm_ko,
            "pct_CPM_remaining": pct,
            "log2FC_KO_minus_WT": log2fc,
            "ci95_low": lo,
            "ci95_high": hi,
            "p_welch": p,
            "q_welch_BH": q,
            "mRNA_call": mRNA_call(mean_wt, mean_ko, pct, p),
            "n_WT": 3,
            "n_KO": 3,
        })
    return pd.DataFrame(rows)


def write_tsv(df: pd.DataFrame, name: str, index: bool = False) -> None:
    path = os.path.join(RES, name)
    df.to_csv(path, sep="\t", index=index)


def make_figures(valid: pd.DataFrame, gstats: pd.DataFrame, sets: dict[str, list[str]]) -> None:
    targets = valid[valid["gene"].isin(KO_TARGETS)].copy()
    # Cldn25 is undetectable; keep it in the table but not as a percent bar.
    plot_df = targets[targets["mean_raw_WT"] >= 1].sort_values("pct_CPM_remaining")

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    colors = ["#E45756" if g == "Cldn4" else "#4C78A8" for g in plot_df["gene"]]
    y = np.arange(len(plot_df))
    ax.barh(y, plot_df["pct_CPM_remaining"], color=colors, height=0.72)
    ax.axvline(100, color="#333333", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["gene"])
    ax.set_xlabel("Cldn-null CPM as % of wild-type CPM")
    ax.set_title("GSE274940 EpH4: targeted claudin mRNA\nCldn4 is not depleted (red)")
    for yi, pct, call in zip(y, plot_df["pct_CPM_remaining"], plot_df["mRNA_call"]):
        ax.text(pct + 1.5, yi, f"{pct:.0f}%  {call.replace('_', ' ')}", va="center", fontsize=8)
    ax.set_xlim(0, max(plot_df["pct_CPM_remaining"]) * 1.45)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.01, 0.01,
        "Cldn25 is also a CRISPR target; expected counts are 0 in all 6 libraries, so it has no bar.",
        fontsize=8, color="#333333",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(os.path.join(RES, "fig_cldn_mRNA_remaining.png"), dpi=160)
    plt.close(fig)

    # Forest of IFN-I and TJ-scaffold genes that were tested.
    focus = []
    for set_name, color in (
        ("TJ_SCAFFOLD_NO_KO_CLDN", "#54A24B"),
        ("IFN_ALPHA_TYPE1", "#F58518"),
        ("ANTIGEN_PRESENTATION_MHC1", "#B279A2"),
    ):
        for gene in sets[set_name]:
            if gene in gstats.index:
                focus.append((gene, set_name, color))
    # one row per gene; IFN and MHC overlap on Stat1/Tap/Psmb — keep first set
    seen = set()
    rows = []
    for gene, set_name, color in focus:
        if gene in seen:
            continue
        seen.add(gene)
        st = gstats.loc[gene]
        rows.append((gene, set_name, color, st.log2FC, st.ci95_low, st.ci95_high))
    rows.sort(key=lambda r: r[3])
    fig, ax = plt.subplots(figsize=(8.2, 9.2))
    y = np.arange(len(rows))
    for i, (gene, set_name, color, fc, lo, hi) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=color, lw=1.4)
        ax.plot(fc, i, "o", color=color, ms=4.5)
    ax.axvline(0, color="#333333", lw=1, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=7)
    ax.set_xlabel("log2FC, Cldn-null minus WT  (95% Welch CI)")
    ax.set_title("IFN-I (orange), MHC-I/APM (purple),\nTJ scaffolds excluding KO claudins (green)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "fig_ifn_tj_forest.png"), dpi=160)
    plt.close(fig)


def main() -> None:
    os.makedirs(RES, exist_ok=True)
    path = download_counts()
    digest = hashlib.sha256(open(path, "rb").read()).hexdigest()
    write_tsv(pd.DataFrame([{
        "file": COUNT_NAME,
        "url": COUNT_URL,
        "sha256": digest,
        "bytes": os.path.getsize(path),
    }]), "downloads.tsv")

    counts, n_ids = load_counts(path)
    if "Cldn4" not in counts.index:
        raise SystemExit("Cldn4 is absent from the count matrix")
    sf = median_ratio_size_factors(counts)
    norm = counts.div(sf, axis=1)
    logm = np.log2(norm + 1.0)
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6

    mean_wt = counts[WT].mean(axis=1)
    mean_ko = counts[KO].mean(axis=1)
    tested = counts.index[(np.maximum(mean_wt, mean_ko) >= EXPR_MIN)]
    gstats = gene_tests(logm, tested)
    gstats.insert(0, "n_ensembl_ids", n_ids.reindex(gstats.index).astype(int))
    gstats.insert(1, "mean_raw_WT", mean_wt.reindex(gstats.index))
    gstats.insert(2, "mean_raw_KO", mean_ko.reindex(gstats.index))
    gstats.insert(3, "mean_CPM_WT", cpm[WT].mean(axis=1).reindex(gstats.index))
    gstats.insert(4, "mean_CPM_KO", cpm[KO].mean(axis=1).reindex(gstats.index))
    gstats["n_WT"] = 3
    gstats["n_KO"] = 3
    gstats["contrast"] = "Cldn-null minus WT"
    gstats["log2_scale"] = "log2(median-of-ratios normalized expected count + 1)"

    sets = load_sets()
    setdf = set_tests(logm, gstats, sets)
    valid = validity_table(counts, cpm, gstats)

    focus_genes = []
    for genes in sets.values():
        focus_genes.extend(genes)
    focus_genes = list(dict.fromkeys(KO_TARGETS + ["Cldn1", "Tacstd2"] + focus_genes))
    focus = gstats.reindex([g for g in focus_genes if g in gstats.index]).copy()
    focus.insert(0, "in_crispr_targets", focus.index.isin(KO_TARGETS))

    meta = pd.DataFrame([{
        "sample": s,
        "gsm": GSM[s],
        "genotype": "WT" if s in WT else "Cldn-null",
        "library_expected_counts": float(lib[s]),
        "size_factor": float(sf[s]),
        "day": 7,
        "cell_line": "EpH4",
        "organism": "Mus musculus",
    } for s in SAMPLES])

    write_tsv(meta, "sample_metadata.tsv")
    per_sample_rows = []
    for gene in KO_TARGETS + ["Cldn1"]:
        rec = {"gene": gene, "crispr_target": gene in KO_TARGETS}
        for s in SAMPLES:
            rec[f"raw_{s}"] = float(counts.loc[gene, s])
            rec[f"CPM_{s}"] = float(cpm.loc[gene, s])
        per_sample_rows.append(rec)
    write_tsv(pd.DataFrame(per_sample_rows), "cldn_target_per_sample.tsv")
    write_tsv(valid, "cldn_mRNA_validity.tsv")
    write_tsv(gstats.sort_values("p_welch"), "gene_stats_CldnNull_vs_WT.tsv", index=True)
    write_tsv(focus, "focus_genes.tsv", index=True)
    write_tsv(setdf, "set_stats.tsv")
    make_figures(valid, gstats, sets)

    c4 = valid.set_index("gene").loc["Cldn4"]
    print(f"sha256 {digest}")
    print("size factors", sf.round(3).to_dict())
    print("Cldn4 call", c4["mRNA_call"],
          f"pct {c4['pct_CPM_remaining']:.2f}",
          f"log2FC {c4['log2FC_KO_minus_WT']:+.3f}",
          f"CI {c4['ci95_low']:+.3f} {c4['ci95_high']:+.3f}",
          f"p {c4['p_welch']:.4g} q {c4['q_welch_BH']:.4g}")
    print(valid.loc[valid["crispr_target"],
                    ["gene", "mean_raw_WT", "mean_raw_KO", "pct_CPM_remaining",
                     "log2FC_KO_minus_WT", "p_welch", "q_welch_BH", "mRNA_call"]]
          .to_string(index=False))
    print(setdf[["set_name", "n_genes_tested", "sample_delta_KO_minus_WT",
                 "sample_p", "sample_q_BH", "comp_median_log2FC",
                 "comp_mannwhitney_p", "comp_q_BH",
                 "n_log2FC_positive", "n_log2FC_negative", "n_welch_p_lt_0.05"]]
          .to_string(index=False))
    # A few single genes the write-up quotes.
    for g in ["Cldn1", "Cldn3", "Cldn4", "Cldn23", "Tacstd2", "Ocln", "Tjp1",
              "F11r", "Cdh1", "Cd274", "B2m", "H2-K1", "Stat1", "Isg15",
              "Ifit1", "Mx1", "Oasl2", "Cxcl10", "Psmb9"]:
        if g not in gstats.index:
            print(f"{g}: not tested")
            continue
        r = gstats.loc[g]
        print(f"{g:10} log2FC {r.log2FC:+7.3f}  CI {r.ci95_low:+7.3f},{r.ci95_high:+7.3f}"
              f"  p {r.p_welch:.3g}  q {r.q_welch_BH:.3g}"
              f"  CPM {r.mean_CPM_WT:.2f}->{r.mean_CPM_KO:.2f}")


if __name__ == "__main__":
    main()
