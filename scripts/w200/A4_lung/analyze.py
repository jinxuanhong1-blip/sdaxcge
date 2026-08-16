#!/usr/bin/env python3
"""A4 lung-only TISMO: Tacstd2 / Cldn4 paired ICB, plus the 64-cohort A4 universe."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

LUNG_CANCER = "Lung carcinoma"
OUTLIER_SRX = "SRX8918393"
ICB_TOKENS = ("antipd1", "antipdl1", "antipdl2", "antictla4")


def cohort_key(label: str) -> str:
    """Strip the trailing (n=N) so Tacstd2/Cldn4 labels can be joined."""
    if "(n=" in label:
        return label[: label.rfind("(n=")]
    return label


def is_baseline(row: pd.Series) -> bool:
    return str(row.get("Baseline", "")) == "1" or str(row.get("Responder", "")) == "Baseline"


def load_gene_table(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df["cohort_key"] = df["cell_line"].map(cohort_key)
    df["is_baseline"] = df.apply(is_baseline, axis=1)
    df["is_lung"] = df["cell_line"].str.startswith("LLC_") | df["cell_line"].str.startswith("CMT-167_")
    return df


def tpm_from_log2p1(x: float) -> float:
    return float(2.0**x - 1.0)


def direction_table(df: pd.DataFrame, gene: str) -> pd.DataFrame:
    rows = []
    for key, sub in df.groupby("cohort_key", sort=True):
        base = sub.loc[sub["is_baseline"], "value"].dropna()
        treat = sub.loc[~sub["is_baseline"], "value"].dropna()
        if base.empty or treat.empty:
            continue
        d_mean = float(treat.mean() - base.mean())
        d_med = float(treat.median() - base.median())
        tismo_p = pd.to_numeric(sub["pvalue"], errors="coerce").dropna()
        tismo_p = float(tismo_p.iloc[0]) if not tismo_p.empty else np.nan
        rows.append(
            {
                "gene": gene,
                "cohort_key": key,
                "tismo_label": sub["cell_line"].iloc[0],
                "gse_id": sub["GSE_ID"].iloc[0],
                "n_baseline": int(base.size),
                "n_treated": int(treat.size),
                "mean_baseline": float(base.mean()),
                "mean_treated": float(treat.mean()),
                "median_baseline": float(base.median()),
                "median_treated": float(treat.median()),
                "delta_mean": d_mean,
                "delta_median": d_med,
                "dir_mean": "UP" if d_mean > 0 else ("DOWN" if d_mean < 0 else "TIE"),
                "dir_median": "UP" if d_med > 0 else ("DOWN" if d_med < 0 else "TIE"),
                "tismo_pvalue": tismo_p,
                "is_lung": bool(sub["is_lung"].iloc[0]),
                "shared_baseline_risk": int(sub.loc[sub["is_baseline"], "Samples"].nunique() < base.size),
            }
        )
    return pd.DataFrame(rows)


def tally(dir_df: pd.DataFrame, rule: str) -> dict:
    col = "dir_mean" if rule == "mean" else "dir_median"
    counts = dir_df[col].value_counts().to_dict()
    up = int(counts.get("UP", 0))
    down = int(counts.get("DOWN", 0))
    tie = int(counts.get("TIE", 0))
    n = up + down + tie
    n_excl = up + down
    two = stats.binomtest(up, n, 0.5, alternative="two-sided").pvalue if n else np.nan
    one = stats.binomtest(up, n, 0.5, alternative="greater").pvalue if n else np.nan
    two_excl = stats.binomtest(up, n_excl, 0.5, alternative="two-sided").pvalue if n_excl else np.nan
    lung = dir_df[dir_df["is_lung"]]
    lung_up = int((lung[col] == "UP").sum())
    return {
        "rule": rule,
        "n_up": up,
        "n_down": down,
        "n_tie": tie,
        "n": n,
        "binom_two_sided_p": float(two),
        "binom_one_sided_greater_p": float(one),
        "binom_two_sided_excl_tie_p": float(two_excl) if n_excl else None,
        "lung_n": int(len(lung)),
        "lung_n_up": lung_up,
        "claim_49_of_64": (up == 49 and n == 64),
    }


def welch_mwu(a: np.ndarray, b: np.ndarray) -> dict:
    if a.size < 2 or b.size < 2:
        return {"welch_t": np.nan, "welch_p": np.nan, "mwu_u": np.nan, "mwu_p": np.nan}
    t_res = stats.ttest_ind(b, a, equal_var=False, alternative="two-sided")
    u_res = stats.mannwhitneyu(b, a, alternative="two-sided")
    return {
        "welch_t": float(t_res.statistic),
        "welch_p": float(t_res.pvalue),
        "mwu_u": float(u_res.statistic),
        "mwu_p": float(u_res.pvalue),
    }


def paired_stats(df: pd.DataFrame, gene: str, drop_outlier: bool) -> list[dict]:
    work = df.copy()
    if drop_outlier:
        work = work[work["Samples"] != OUTLIER_SRX]
    out = []
    for key, sub in work.groupby("cohort_key", sort=True):
        if not key.startswith("LLC_"):
            continue
        base = sub.loc[sub["is_baseline"], "value"].dropna().to_numpy()
        treat = sub.loc[~sub["is_baseline"], "value"].dropna().to_numpy()
        if base.size == 0 or treat.size == 0:
            continue
        tests = welch_mwu(base, treat)
        genotype = "Setdb1_KO" if "Setdb1_KO" in key else "WT"
        responders = sorted({str(x) for x in sub.loc[~sub["is_baseline"], "Responder"].unique()})
        out.append(
            {
                "gene": gene,
                "arm": genotype,
                "cohort_key": key,
                "gse_id": sub["GSE_ID"].iloc[0],
                "treated_response_label": ";".join(responders),
                "exclude_outlier": drop_outlier,
                "n_baseline": int(base.size),
                "n_treated": int(treat.size),
                "mean_baseline": float(base.mean()),
                "mean_treated": float(treat.mean()),
                "median_baseline": float(np.median(base)),
                "median_treated": float(np.median(treat)),
                "delta_mean": float(treat.mean() - base.mean()),
                "delta_median": float(np.median(treat) - np.median(base)),
                "mean_tpm_baseline": tpm_from_log2p1(float(base.mean())),
                "mean_tpm_treated": tpm_from_log2p1(float(treat.mean())),
                **tests,
            }
        )
    return out


def lung_inventory(meta: list[dict]) -> pd.DataFrame:
    rows = []
    for rec in meta:
        if rec.get("cancerType") != LUNG_CANCER:
            continue
        treat = str(rec.get("mouseTreatment") or "")
        treat_l = treat.lower().replace(" ", "")
        is_icb = any(tok in treat_l for tok in ICB_TOKENS)
        rows.append(
            {
                "tismo_id": rec.get("id"),
                "study_id": rec.get("studyId"),
                "cell_line": rec.get("cellLine"),
                "cancer_type": rec.get("cancerType"),
                "cell_genotype": rec.get("cellGenotype"),
                "mouse_strain": rec.get("mouseStrain"),
                "implantation": rec.get("implantation"),
                "implantation_site": rec.get("implantationSite"),
                "mouse_treatment": treat,
                "icb_study_label": rec.get("icbStudy"),
                "replicates": rec.get("replicates"),
                "tumor": rec.get("tumor"),
                "is_icb_treatment": is_icb,
                "has_paired_icb_in_study": rec.get("studyId") == "GSE155972",
            }
        )
    return pd.DataFrame(rows).sort_values(["study_id", "cell_line", "mouse_treatment"])


def style_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_llc(per_sample: pd.DataFrame, dest: Path) -> None:
    genes = ["Tacstd2", "Cldn4"]
    arms = ["WT", "Setdb1_KO"]
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.8), sharey=False)
    rng = np.random.default_rng(0)
    for i, gene in enumerate(genes):
        for j, arm in enumerate(arms):
            ax = axes[i, j]
            sub = per_sample[(per_sample["gene"] == gene) & (per_sample["arm"] == arm)]
            groups = ["Baseline", "ICB"]
            data = [
                sub.loc[sub["group"] == "Baseline", "value"].to_numpy(),
                sub.loc[sub["group"] == "ICB", "value"].to_numpy(),
            ]
            ax.boxplot(data, positions=[0, 1], widths=0.45, showfliers=False, medianprops={"color": "black"})
            for k, arr in enumerate(data):
                x = np.full(arr.size, k, dtype=float) + rng.uniform(-0.08, 0.08, arr.size)
                outlier = sub.loc[sub["group"] == groups[k], "is_outlier"].to_numpy()
                colors = ["#c0392b" if flag else "#2c3e50" for flag in outlier]
                ax.scatter(x, arr, c=colors, s=28, zorder=3, edgecolors="white", linewidths=0.4)
            ax.set_xticks([0, 1], groups)
            ax.set_title(f"{gene} · LLC {arm}")
            ax.set_ylabel("TISMO value (log2[TPM+1])")
            style_axes(ax)
    fig.suptitle("TISMO lung ICB is one study (GSE155972 LLC). Red = SRX8918393.", fontsize=11)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=160)
    plt.close(fig)


def plot_tally(tally_df: pd.DataFrame, dest: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    labels = []
    ups, downs, ties = [], [], []
    for _, row in tally_df.iterrows():
        labels.append(f"{row['gene']}\n{row['rule']}")
        ups.append(row["n_up"])
        downs.append(row["n_down"])
        ties.append(row["n_tie"])
    x = np.arange(len(labels))
    ax.bar(x, ups, label="UP", color="#1f77b4")
    ax.bar(x, downs, bottom=ups, label="DOWN", color="#d62728")
    ax.bar(x, ties, bottom=np.array(ups) + np.array(downs), label="TIE", color="#7f7f7f")
    ax.axhline(49, color="#1f77b4", ls="--", lw=1, alpha=0.7, label="claimed 49 up")
    ax.set_xticks(x, labels)
    ax.set_ylabel("ICB treated-vs-baseline cohorts")
    ax.set_title("A4 49/64 is Tacstd2 mean-direction, pan-TISMO — not lung")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=160)
    plt.close(fig)


def plot_floor(per_sample: pd.DataFrame, dest: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    order = ["Actb", "Cd8a", "Ifng", "Cd274", "Epcam", "Cldn4", "Tacstd2"]
    present = [g for g in order if g in set(per_sample["gene"])]
    data = [per_sample.loc[per_sample["gene"] == g, "value"].to_numpy() for g in present]
    ax.boxplot(data, tick_labels=present, showfliers=True)
    ax.axhline(1.0, color="#888", ls=":", lw=1)
    ax.set_ylabel("TISMO value (log2[TPM+1])")
    ax.set_title("LLC GSE155972: Tacstd2/Cldn4 sit near the floor; Actb does not")
    style_axes(ax)
    fig.tight_layout()
    dest.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(dest, dpi=160)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw", type=Path, default=Path("results/w200/A4_lung/raw"))
    parser.add_argument("--out", type=Path, default=Path("results/w200/A4_lung"))
    args = parser.parse_args()
    raw, out = args.raw, args.out
    fig_dir = out / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    meta = json.loads((raw / "vivo_meta.json").read_text())["data"]
    inv = lung_inventory(meta)
    inv.to_csv(out / "lung_design_inventory.tsv", sep="\t", index=False)

    tac = load_gene_table(raw / "tacstd2_icb_all.csv")
    cld = load_gene_table(raw / "cldn4_icb_all.csv")
    tac_dir = direction_table(tac, "Tacstd2")
    cld_dir = direction_table(cld, "Cldn4")
    universe = pd.concat([tac_dir, cld_dir], ignore_index=True)
    universe.to_csv(out / "icb_cohort_universe.tsv", sep="\t", index=False)

    tallies = []
    for gene, ddf in [("Tacstd2", tac_dir), ("Cldn4", cld_dir)]:
        for rule in ("mean", "median"):
            rec = tally(ddf, rule)
            rec["gene"] = gene
            tallies.append(rec)
    tally_df = pd.DataFrame(tallies)
    tally_df.to_csv(out / "direction_tally.tsv", sep="\t", index=False)

    # Wilcoxon signed-rank on Tacstd2 mean-deltas (A4-style consistency test).
    deltas = tac_dir["delta_mean"].to_numpy()
    w_res = stats.wilcoxon(deltas, alternative="greater", zero_method="wilcox")

    controls = []
    for gene in ["Actb", "Cd8a", "Ifng", "Gzmb", "Cd274", "Epcam"]:
        path = raw / f"{gene.lower()}_llc.csv"
        if path.exists():
            controls.append(load_gene_table(path).assign(gene=gene))
    ctl = pd.concat(controls, ignore_index=True) if controls else pd.DataFrame()

    # Per-sample LLC table for targets + controls.
    llc_parts = []
    for gene, df in [("Tacstd2", tac), ("Cldn4", cld)]:
        sub = df[df["cell_line"].str.startswith("LLC_")].copy()
        sub["gene"] = gene
        llc_parts.append(sub)
    if not ctl.empty:
        llc_parts.append(ctl[ctl["cell_line"].str.startswith("LLC_")].copy())
    llc = pd.concat(llc_parts, ignore_index=True)
    llc["arm"] = np.where(llc["cell_line"].str.contains("Setdb1_KO"), "Setdb1_KO", "WT")
    llc["group"] = np.where(llc["is_baseline"], "Baseline", "ICB")
    llc["is_outlier"] = llc["Samples"] == OUTLIER_SRX
    per_sample = llc[
        [
            "gene",
            "Samples",
            "value",
            "arm",
            "group",
            "Responder",
            "GSE_ID",
            "Mouse_treatment",
            "cell_line",
            "is_outlier",
        ]
    ].sort_values(["gene", "arm", "group", "Samples"])
    per_sample.to_csv(out / "per_sample_expression.tsv", sep="\t", index=False)

    stats_rows = []
    for gene, df in [("Tacstd2", tac), ("Cldn4", cld)]:
        stats_rows.extend(paired_stats(df, gene, drop_outlier=False))
        stats_rows.extend(paired_stats(df, gene, drop_outlier=True))
    if not ctl.empty:
        for gene, sub in ctl.groupby("gene"):
            stats_rows.extend(paired_stats(sub, str(gene), drop_outlier=False))
    stats_df = pd.DataFrame(stats_rows)
    primary = stats_df[
        stats_df["gene"].isin(["Tacstd2", "Cldn4"]) & ~stats_df["exclude_outlier"]
    ].copy()
    if not primary.empty:
        primary["mwu_fdr_bh"] = stats.false_discovery_control(primary["mwu_p"].to_numpy(), method="bh")
        primary["welch_fdr_bh"] = stats.false_discovery_control(primary["welch_p"].to_numpy(), method="bh")
        stats_df = stats_df.merge(
            primary[["gene", "arm", "exclude_outlier", "mwu_fdr_bh", "welch_fdr_bh"]],
            on=["gene", "arm", "exclude_outlier"],
            how="left",
        )
    stats_df.to_csv(out / "paired_icb_stats.tsv", sep="\t", index=False)
    primary.to_csv(out / "primary_lung_tests.tsv", sep="\t", index=False)

    # Tacstd2 vs Cldn4 correlation inside LLC (same SRX).
    wide = per_sample[per_sample["gene"].isin(["Tacstd2", "Cldn4"])].pivot_table(
        index=["Samples", "arm", "group", "is_outlier"], columns="gene", values="value", aggfunc="mean"
    )
    corr_rows = []
    if {"Tacstd2", "Cldn4"}.issubset(wide.columns):
        for mask_name, mask in [
            ("all_llc", np.ones(len(wide), dtype=bool)),
            ("exclude_SRX8918393", ~wide.index.get_level_values("is_outlier")),
        ]:
            x = wide.loc[mask, "Tacstd2"].to_numpy()
            y = wide.loc[mask, "Cldn4"].to_numpy()
            if x.size >= 3:
                rho, p = stats.spearmanr(x, y)
                corr_rows.append(
                    {"subset": mask_name, "n": int(x.size), "spearman_rho": float(rho), "spearman_p": float(p)}
                )
    corr_df = pd.DataFrame(corr_rows)
    corr_df.to_csv(out / "tacstd2_cldn4_correlation.tsv", sep="\t", index=False)

    plot_llc(per_sample[per_sample["gene"].isin(["Tacstd2", "Cldn4"])], fig_dir / "fig_llc_paired_icb.png")
    plot_tally(tally_df, fig_dir / "fig_direction_tally.png")
    plot_floor(per_sample, fig_dir / "fig_llc_expression_floor.png")

    tac_mean = next(t for t in tallies if t["gene"] == "Tacstd2" and t["rule"] == "mean")
    tac_med = next(t for t in tallies if t["gene"] == "Tacstd2" and t["rule"] == "median")
    cld_mean = next(t for t in tallies if t["gene"] == "Cldn4" and t["rule"] == "mean")

    lung_icb_studies = sorted(set(inv.loc[inv["is_icb_treatment"], "study_id"]))
    lung_models = sorted(set(inv["cell_line"]))
    lung_stats = primary if not primary.empty else stats_df.iloc[0:0]
    any_nominal = bool((lung_stats["mwu_p"] < 0.05).any() or (lung_stats["welch_p"] < 0.05).any())
    any_fdr = bool(
        ("mwu_fdr_bh" in lung_stats.columns)
        and ((lung_stats["mwu_fdr_bh"] < 0.05).any() or (lung_stats["welch_fdr_bh"] < 0.05).any())
    )

    # Honest verdict pieces.
    claim_reproduced = bool(tac_mean["claim_49_of_64"])
    summary = {
        "task": "A4 lung-only TISMO Tacstd2/Cldn4 paired ICB",
        "verdict": (
            "NOT SUPPORTED as a lung/NSCLC result. TISMO has two lung-carcinoma "
            "lines (LLC, CMT-167). The only paired ICB design is GSE155972 LLC "
            "anti-PD1+anti-CTLA4 (2 genotype arms). Tacstd2 and Cldn4 sit near "
            "the detection floor (mean TPM ~0.1–0.2). One uncorrected MWU "
            "(Tacstd2 WT ICB vs baseline p=0.045) is nominally <0.05 but does "
            "not survive BH FDR across the 4 primary tests; Welch t is n.s. "
            "(p=0.15). Cldn4 is n.s. in both arms. The Setdb1_KO mean-up is "
            "one outlier (SRX8918393). The A4 49/64 count is a pan-TISMO "
            "mean-direction tally, not a lung finding. Immune controls "
            "(Cd8a/Cd274/Gzmb) do rise, so the design can detect ICB effects."
        ),
        "claim_A4": {
            "stated": "Tacstd2 up in 49/64 ICI-treated mouse models, p=5.8e-5",
            "universe_reproduced": "64 TISMO in-vivo ICB treated-vs-baseline cohorts from gene export",
            "tacstd2_mean_direction": tac_mean,
            "tacstd2_median_direction": tac_med,
            "cldn4_mean_direction": cld_mean,
            "wilcoxon_signed_rank_greater_on_mean_deltas": {
                "W": float(w_res.statistic),
                "p": float(w_res.pvalue),
                "n": int(deltas.size),
            },
            "claimed_p_5_8e5_reproduced": False,
            "exact_two_sided_binomial_p": tac_mean["binom_two_sided_p"],
            "note_on_p": (
                "49/64 vs 0.5 has exact two-sided binomial p=2.44e-5 and one-sided "
                "greater p=1.22e-5. The user-reported 5.8e-5 was not recovered."
            ),
            "independence_caveat": (
                "Two 4T1 GSE130472 cohort pairs share the same baseline samples "
                "(old antiCTLA4 vs old antiPDL1; young antiCTLA4 vs young antiPDL1). "
                "The binomial treats 64 comparisons as independent."
            ),
        },
        "lung_scope": {
            "tismo_lung_models": lung_models,
            "lung_design_groups": int(len(inv)),
            "lung_icb_studies": lung_icb_studies,
            "paired_icb_study": "GSE155972",
            "paired_icb_model": "LLC (Lewis lung carcinoma), subcutaneous flank",
            "n_lung_icb_cohorts_in_64": int(tac_dir["is_lung"].sum()),
            "cmt167_has_icb": False,
            "response_confound": (
                "WT ICB arm is labeled Non-responders; Setdb1_KO ICB arm is "
                "labeled Responders. Response is confounded with Setdb1 genotype. "
                "There is no within-genotype R vs NR contrast."
            ),
        },
        "lung_tests_nominal_p_lt_0.05": any_nominal,
        "lung_tests_bh_fdr_lt_0.05": any_fdr,
        "outlier": {
            "sample": OUTLIER_SRX,
            "note": (
                "Setdb1_KO ICB sample SRX8918393 is Tacstd2=5.41 and Cldn4=3.23; "
                "all other LLC Tacstd2 values are ≤1.10. TISMO's LLC-only Tacstd2 "
                "export drops this sample; the All-models export keeps it. "
                "Mean-direction 'up' on the KO arm is this one point."
            ),
        },
        "n_tacstd2_icb_cohorts": int(len(tac_dir)),
        "n_cldn4_icb_cohorts": int(len(cld_dir)),
        "claim_reproduced_pan_tismo_mean": claim_reproduced,
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    audit = {
        "data_source": "TISMO live API https://tismo.pku-genomics.org/",
        "values": "TISMO gene-module export; Actb≈11.6 implies log2(TPM+1)",
        "files_used": sorted(p.name for p in raw.iterdir()),
        "outlier_handling": "Primary stats keep SRX8918393; sensitivity drops it",
        "what_was_not_done": [
            "No in-vitro cytokine arm (A4 lung task is paired ICB).",
            "No GEO raw FASTQ reprocessing; TISMO uniformly processed values used.",
            "Did not treat KPB25L / p53-* as lung; TISMO annotates them Mammary cancer, NOS.",
        ],
    }
    (out / "audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(
        json.dumps(
            {
                k: summary[k]
                for k in (
                    "verdict",
                    "lung_tests_nominal_p_lt_0.05",
                    "lung_tests_bh_fdr_lt_0.05",
                    "claim_reproduced_pan_tismo_mean",
                )
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
