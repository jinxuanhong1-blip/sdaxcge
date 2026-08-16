#!/usr/bin/env python3
"""CPTAC LSCC/LUSC TACSTD2 + CLDN4 protein, RNA, and immune analysis.

Treatment-naive surgical cohort (Satpathy et al. Cell 2021). No ICI labels.
Open S3 freeze v1.2. Real pairwise-complete stats only.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ensembl_map import PHENOTYPE_IMMUNE, SIGNATURES, TARGETS
from lib_stats import (
    bh_fdr,
    mannwhitney_two_sided,
    match_ensembl_row,
    median_split,
    safe_log10_p,
    signature_score,
    spearman_pair,
    strip_ensembl_version,
    wilcoxon_signed_rank,
)

DATA = ROOT / "data" / "grok_cptac_lscc"
TABLES = ROOT / "results" / "grok_cptac_lscc" / "tables"
FIGS = ROOT / "results" / "grok_cptac_lscc" / "figures"
NOTES = ROOT / "notes" / "grok_cptac_lscc"

PROTEIN_T = DATA / "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"
PROTEIN_N = DATA / "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Normal.txt"
RNA_T = DATA / "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
RNA_N = DATA / "LSCC_RNAseq_gene_RSEM_coding_UQ_1500_log2_Normal.txt"
PHENO = DATA / "LSCC_phenotype.txt"
SURV = DATA / "LSCC_survival.txt"
META = DATA / "LSCC_meta.txt"

sns.set_theme(style="whitegrid", context="talk")


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df.set_index(df.columns[0])
    df.index = df.index.astype(str)
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def extract_targets(mat: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    rows = {}
    info = {}
    for symbol, ens in TARGETS.items():
        raw = match_ensembl_row(mat.index, ens)
        if raw is None:
            rows[symbol] = pd.Series(np.nan, index=mat.columns, name=symbol)
            info[symbol] = {"ensembl": ens, "row": None, "n": 0, "n_na": int(mat.shape[1])}
            continue
        s = mat.loc[raw]
        if isinstance(s, pd.DataFrame):
            s = s.iloc[0]
        s = pd.to_numeric(s, errors="coerce")
        s.name = symbol
        rows[symbol] = s
        info[symbol] = {
            "ensembl": ens,
            "row": raw,
            "n": int(s.notna().sum()),
            "n_na": int(s.isna().sum()),
            "min": float(s.min()) if s.notna().any() else np.nan,
            "median": float(s.median()) if s.notna().any() else np.nan,
            "max": float(s.max()) if s.notna().any() else np.nan,
        }
    return pd.DataFrame(rows), info


def reindex_ensembl(mat: pd.DataFrame) -> pd.DataFrame:
    out = mat.copy()
    out.index = strip_ensembl_version(out.index)
    out = out[~out.index.duplicated(keep="first")]
    return out


def load_meta() -> pd.DataFrame:
    meta = pd.read_csv(META, sep="\t")
    meta = meta[meta["case_id"] != "data_type"].copy()
    meta = meta.set_index("case_id")
    return meta


def load_pheno() -> pd.DataFrame:
    ph = pd.read_csv(PHENO, sep="\t")
    ph = ph.set_index(ph.columns[0])
    ph.index.name = "case_id"
    keep = [c for c in PHENOTYPE_IMMUNE if c in ph.columns]
    missing = [c for c in PHENOTYPE_IMMUNE if c not in ph.columns]
    if missing:
        print("WARN missing phenotype columns:", missing, flush=True)
    return ph[keep].apply(pd.to_numeric, errors="coerce")


def load_surv() -> pd.DataFrame:
    s = pd.read_csv(SURV, sep="\t").set_index("case_id")
    return s.apply(pd.to_numeric, errors="coerce")


def corr_table(features: pd.DataFrame, scores: pd.DataFrame, feature_kind: str, score_kind: str) -> pd.DataFrame:
    rows = []
    for feat in features.columns:
        for sc in scores.columns:
            res = spearman_pair(features[feat], scores[sc])
            rows.append(
                {
                    "feature": feat,
                    "feature_kind": feature_kind,
                    "score": sc,
                    "score_kind": score_kind,
                    "n": res["n"],
                    "rho": res["rho"],
                    "p": res["p"],
                }
            )
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["p"])
    return out.sort_values(["fdr", "p"], na_position="last")


def tumor_vs_normal(tumor: pd.Series, normal: pd.Series, gene: str, layer: str) -> dict:
    paired = wilcoxon_signed_rank(tumor, normal)
    unpaired = mannwhitney_two_sided(tumor, normal)
    return {
        "gene": gene,
        "layer": layer,
        "n_tumor": int(pd.to_numeric(tumor, errors="coerce").notna().sum()),
        "n_normal": int(pd.to_numeric(normal, errors="coerce").notna().sum()),
        "median_tumor": float(pd.to_numeric(tumor, errors="coerce").median()),
        "median_normal": float(pd.to_numeric(normal, errors="coerce").median()),
        "paired_n": paired["n"],
        "paired_W": paired["W"],
        "paired_p": paired["p"],
        "paired_median_tumor_minus_normal": paired["median_delta"],
        "unpaired_n_tumor": unpaired["n_a"],
        "unpaired_n_normal": unpaired["n_b"],
        "unpaired_U": unpaired["U"],
        "unpaired_p": unpaired["p"],
        "unpaired_rank_biserial": unpaired["rank_biserial"],
    }


def high_low_immune(feature: pd.Series, scores: pd.DataFrame, gene: str, layer: str) -> pd.DataFrame:
    grp = median_split(feature)
    rows = []
    for sc in scores.columns:
        high = scores.loc[grp == "high", sc]
        low = scores.loc[grp == "low", sc]
        res = mannwhitney_two_sided(high, low)
        rows.append(
            {
                "gene": gene,
                "layer": layer,
                "score": sc,
                "n_high": res["n_a"],
                "n_low": res["n_b"],
                "U": res["U"],
                "p": res["p"],
                "rank_biserial_high_vs_low": res["rank_biserial"],
                "median_high": float(pd.to_numeric(high, errors="coerce").median()) if high.notna().any() else np.nan,
                "median_low": float(pd.to_numeric(low, errors="coerce").median()) if low.notna().any() else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    out["fdr"] = bh_fdr(out["p"])
    return out


def survival_tests(feature: pd.Series, surv: pd.DataFrame, gene: str, layer: str) -> list[dict]:
    rows = []
    for time_col, event_col, endpoint in [
        ("OS_days", "OS_event", "OS"),
        ("PFS_days", "PFS_event", "PFS"),
    ]:
        d = pd.concat(
            [
                feature.rename("value"),
                surv[time_col].rename("time"),
                surv[event_col].rename("event"),
            ],
            axis=1,
        ).dropna()
        d = d[(d["time"] > 0)]
        n = int(d.shape[0])
        n_event = int(d["event"].sum())
        rec = {
            "gene": gene,
            "layer": layer,
            "endpoint": endpoint,
            "n": n,
            "n_event": n_event,
            "logrank_p": np.nan,
            "n_high": np.nan,
            "n_low": np.nan,
            "events_high": np.nan,
            "events_low": np.nan,
            "cox_hr_per_unit": np.nan,
            "cox_p": np.nan,
            "cox_hr_ci_low": np.nan,
            "cox_hr_ci_high": np.nan,
        }
        if n < 10 or n_event < 5:
            rows.append(rec)
            continue
        grp = median_split(d["value"])
        high = d[grp == "high"]
        low = d[grp == "low"]
        rec["n_high"] = int(high.shape[0])
        rec["n_low"] = int(low.shape[0])
        rec["events_high"] = int(high["event"].sum())
        rec["events_low"] = int(low["event"].sum())
        if high.shape[0] >= 5 and low.shape[0] >= 5:
            lr = logrank_test(high["time"], low["time"], high["event"], low["event"])
            rec["logrank_p"] = float(lr.p_value)
        try:
            cph = CoxPHFitter()
            cph.fit(d[["value", "time", "event"]], duration_col="time", event_col="event")
            rec["cox_hr_per_unit"] = float(cph.hazard_ratios_["value"])
            rec["cox_p"] = float(cph.summary.loc["value", "p"])
            rec["cox_hr_ci_low"] = float(np.exp(cph.confidence_intervals_.loc["value"].iloc[0]))
            rec["cox_hr_ci_high"] = float(np.exp(cph.confidence_intervals_.loc["value"].iloc[1]))
        except Exception as exc:
            rec["cox_error"] = str(exc)
        rows.append(rec)
    return rows


def savefig(fig: plt.Figure, name: str) -> None:
    path = FIGS / name
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("wrote", path, flush=True)


def plot_tumor_normal(tumor: pd.Series, normal: pd.Series, title: str, ylabel: str, fname: str) -> None:
    t = pd.to_numeric(tumor, errors="coerce").dropna()
    n = pd.to_numeric(normal, errors="coerce").dropna()
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    sns.boxplot(
        data=pd.DataFrame(
            {
                "value": pd.concat([t, n], ignore_index=True),
                "group": ["Tumor"] * len(t) + ["NAT"] * len(n),
            }
        ),
        x="group",
        y="value",
        hue="group",
        palette={"Tumor": "#c0392b", "NAT": "#2980b9"},
        legend=False,
        ax=ax,
        width=0.55,
    )
    sns.stripplot(
        data=pd.DataFrame(
            {
                "value": pd.concat([t, n], ignore_index=True),
                "group": ["Tumor"] * len(t) + ["NAT"] * len(n),
            }
        ),
        x="group",
        y="value",
        color="black",
        size=3,
        alpha=0.35,
        ax=ax,
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(ylabel)
    savefig(fig, fname)


def plot_scatter(x: pd.Series, y: pd.Series, xlabel: str, ylabel: str, title: str, fname: str) -> None:
    d = pd.concat([x.rename("x"), y.rename("y")], axis=1).dropna()
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ax.scatter(d["x"], d["y"], s=28, alpha=0.75, c="#1f4e79", edgecolors="none")
    if d.shape[0] >= 4:
        res = spearman_pair(d["x"], d["y"])
        ax.set_title(f"{title}\nSpearman ρ={res['rho']:.3f} p={res['p']:.2e} n={res['n']}")
    else:
        ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    savefig(fig, fname)


def plot_heatmap(corr: pd.DataFrame, fname: str, title: str) -> None:
    sub = corr.copy()
    mat = sub.pivot(index="score", columns="feature", values="rho")
    fig, ax = plt.subplots(figsize=(7.2, max(6.5, 0.28 * mat.shape[0] + 2)))
    sns.heatmap(
        mat,
        cmap="RdBu_r",
        center=0,
        vmin=-0.6,
        vmax=0.6,
        ax=ax,
        cbar_kws={"label": "Spearman ρ"},
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("")
    savefig(fig, fname)


def plot_volcano(corr: pd.DataFrame, fname: str, title: str) -> None:
    d = corr.dropna(subset=["rho", "p"]).copy()
    d["mlogp"] = d["p"].map(safe_log10_p)
    fig, ax = plt.subplots(figsize=(7.5, 5.6))
    colors = {"TACSTD2": "#c0392b", "CLDN4": "#2980b9"}
    for feat, g in d.groupby("feature"):
        ax.scatter(g["rho"], g["mlogp"], s=26, alpha=0.8, label=feat, c=colors.get(feat, "gray"))
    ax.axhline(-np.log10(0.05), ls="--", lw=1, color="gray")
    ax.axvline(0, ls=":", lw=1, color="gray")
    ax.set_xlabel("Spearman ρ")
    ax.set_ylabel("−log10 p")
    ax.set_title(title)
    ax.legend(frameon=False)
    # label top hits
    top = d.sort_values("p").head(8)
    for _, r in top.iterrows():
        ax.annotate(
            f"{r['feature'][:3]}:{r['score']}",
            (r["rho"], r["mlogp"]),
            fontsize=7,
            alpha=0.85,
        )
    savefig(fig, fname)


def plot_km(feature: pd.Series, surv: pd.DataFrame, time_col: str, event_col: str, title: str, fname: str) -> None:
    d = pd.concat(
        [feature.rename("value"), surv[time_col].rename("time"), surv[event_col].rename("event")],
        axis=1,
    ).dropna()
    d = d[d["time"] > 0]
    grp = median_split(d["value"])
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    km = KaplanMeierFitter()
    for label, color in [("low", "#2980b9"), ("high", "#c0392b")]:
        sub = d[grp == label]
        if sub.empty:
            continue
        km.fit(sub["time"], sub["event"], label=f"{label} (n={len(sub)}, ev={int(sub['event'].sum())})")
        km.plot(ax=ax, ci_show=True, color=color)
    high = d[grp == "high"]
    low = d[grp == "low"]
    if len(high) >= 5 and len(low) >= 5:
        lr = logrank_test(high["time"], low["time"], high["event"], low["event"])
        ax.set_title(f"{title}\nlog-rank p={lr.p_value:.3g}")
    else:
        ax.set_title(title)
    ax.set_xlabel("Days")
    ax.set_ylabel("Survival probability")
    savefig(fig, fname)


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    NOTES.mkdir(parents=True, exist_ok=True)

    protein_t = load_matrix(PROTEIN_T)
    protein_n = load_matrix(PROTEIN_N)
    rna_t = load_matrix(RNA_T)
    rna_n = load_matrix(RNA_N)
    pheno = load_pheno()
    surv = load_surv()
    meta = load_meta()

    prot_t_feat, prot_t_info = extract_targets(protein_t)
    prot_n_feat, prot_n_info = extract_targets(protein_n)
    rna_t_feat, rna_t_info = extract_targets(rna_t)
    rna_n_feat, rna_n_info = extract_targets(rna_n)

    # Align all tumor tables on protein tumor case IDs (108 LSCC tumors).
    cases = prot_t_feat.index.intersection(rna_t_feat.index).intersection(pheno.index)
    prot_t_feat = prot_t_feat.loc[cases]
    rna_t_feat = rna_t_feat.loc[cases]
    pheno = pheno.loc[pheno.index.intersection(cases)]
    surv = surv.loc[surv.index.intersection(cases)]
    meta = meta.loc[meta.index.intersection(cases)]

    missingness = []
    for layer, info, n_samples in [
        ("protein_tumor", prot_t_info, protein_t.shape[1]),
        ("protein_normal", prot_n_info, protein_n.shape[1]),
        ("rna_tumor", rna_t_info, rna_t.shape[1]),
        ("rna_normal", rna_n_info, rna_n.shape[1]),
    ]:
        for gene, rec in info.items():
            missingness.append({"layer": layer, "gene": gene, "n_samples_in_matrix": n_samples, **rec})
    miss_df = pd.DataFrame(missingness)
    miss_df.to_csv(TABLES / "missingness.tsv", sep="\t", index=False)

    # Tumor vs NAT
    tvn_rows = []
    for gene in TARGETS:
        tvn_rows.append(tumor_vs_normal(prot_t_feat[gene], prot_n_feat[gene], gene, "protein"))
        tvn_rows.append(tumor_vs_normal(rna_t_feat[gene], rna_n_feat[gene], gene, "rna"))
    tvn = pd.DataFrame(tvn_rows)
    tvn["paired_fdr"] = bh_fdr(tvn["paired_p"])
    tvn["unpaired_fdr"] = bh_fdr(tvn["unpaired_p"])
    tvn.to_csv(TABLES / "tumor_vs_normal.tsv", sep="\t", index=False)

    # Protein-RNA and TACSTD2-CLDN4 coexpression
    co_rows = []
    for gene in TARGETS:
        res = spearman_pair(prot_t_feat[gene], rna_t_feat[gene])
        co_rows.append({"comparison": f"{gene}_protein_vs_rna", "n": res["n"], "rho": res["rho"], "p": res["p"]})
    for layer, feat in [("protein", prot_t_feat), ("rna", rna_t_feat)]:
        res = spearman_pair(feat["TACSTD2"], feat["CLDN4"])
        co_rows.append({"comparison": f"TACSTD2_vs_CLDN4_{layer}", "n": res["n"], "rho": res["rho"], "p": res["p"]})
    co = pd.DataFrame(co_rows)
    co["fdr"] = bh_fdr(co["p"])
    co.to_csv(TABLES / "coexpression_protein_rna.tsv", sep="\t", index=False)

    # RNA signature scores on tumor RNA (Ensembl, version stripped)
    rna_ens = reindex_ensembl(rna_t)
    sig_scores = {}
    sig_cov = {}
    for name, mapping in SIGNATURES.items():
        ids = list(mapping.values())
        score, cov = signature_score(rna_ens, ids)
        sig_scores[name] = score
        sig_cov[name] = cov
    sig_df = pd.DataFrame(sig_scores).loc[cases]
    sig_df.to_csv(TABLES / "rna_signature_scores.tsv", sep="\t")
    pd.DataFrame([{"signature": k, "coverage": v} for k, v in sig_cov.items()]).to_csv(
        TABLES / "rna_signature_coverage.tsv", sep="\t", index=False
    )

    # Also try protein-level signatures where genes are quantified
    prot_ens = reindex_ensembl(protein_t)
    prot_sig = {}
    prot_sig_cov = {}
    for name, mapping in SIGNATURES.items():
        score, cov = signature_score(prot_ens, list(mapping.values()))
        prot_sig[name] = score
        prot_sig_cov[name] = cov
    prot_sig_df = pd.DataFrame(prot_sig).loc[prot_sig[next(iter(prot_sig))].index.intersection(cases)]
    prot_sig_df.to_csv(TABLES / "protein_signature_scores.tsv", sep="\t")
    pd.DataFrame([{"signature": k, "coverage": v} for k, v in prot_sig_cov.items()]).to_csv(
        TABLES / "protein_signature_coverage.tsv", sep="\t", index=False
    )

    # Correlations vs freeze deconvolution / HALLMARK / PROGENy
    corr_prot_pheno = corr_table(prot_t_feat, pheno, "protein", "phenotype_immune")
    corr_rna_pheno = corr_table(rna_t_feat, pheno, "rna", "phenotype_immune")
    corr_prot_sig = corr_table(prot_t_feat, sig_df, "protein", "rna_signature")
    corr_rna_sig = corr_table(rna_t_feat, sig_df, "rna", "rna_signature")
    corr_prot_psig = corr_table(prot_t_feat, prot_sig_df, "protein", "protein_signature")
    corr_prot_pheno.to_csv(TABLES / "spearman_protein_vs_phenotype.tsv", sep="\t", index=False)
    corr_rna_pheno.to_csv(TABLES / "spearman_rna_vs_phenotype.tsv", sep="\t", index=False)
    corr_prot_sig.to_csv(TABLES / "spearman_protein_vs_rna_signatures.tsv", sep="\t", index=False)
    corr_rna_sig.to_csv(TABLES / "spearman_rna_vs_rna_signatures.tsv", sep="\t", index=False)
    corr_prot_psig.to_csv(TABLES / "spearman_protein_vs_protein_signatures.tsv", sep="\t", index=False)

    # Median-split MWU on key immune scores
    key_scores = pheno[
        [
            c
            for c in [
                "CIBERSORT_T_cell_CD8+",
                "ESTIMATE_ImmuneScore",
                "xCell_T_cell_CD8+",
                "xCell_immune_score",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                "HALLMARK_INFLAMMATORY_RESPONSE",
                "PROGENy_JAK-STAT",
                "xCell_T_cell_regulatory_(Tregs)",
                "CIBERSORT_Macrophage_M2",
                "HALLMARK_TGF_BETA_SIGNALING",
            ]
            if c in pheno.columns
        ]
    ]
    hl_parts = []
    for gene in TARGETS:
        hl_parts.append(high_low_immune(prot_t_feat[gene], key_scores, gene, "protein"))
        hl_parts.append(high_low_immune(rna_t_feat[gene], key_scores, gene, "rna"))
    hl = pd.concat(hl_parts, ignore_index=True)
    hl.to_csv(TABLES / "median_split_mwu_immune.tsv", sep="\t", index=False)

    # Survival (treatment-naive OS/PFS — not ICI outcome)
    surv_rows = []
    for gene in TARGETS:
        surv_rows.extend(survival_tests(prot_t_feat[gene], surv, gene, "protein"))
        surv_rows.extend(survival_tests(rna_t_feat[gene], surv, gene, "rna"))
    surv_df = pd.DataFrame(surv_rows)
    surv_df["logrank_fdr"] = bh_fdr(surv_df["logrank_p"])
    surv_df["cox_fdr"] = bh_fdr(surv_df["cox_p"])
    surv_df.to_csv(TABLES / "survival_os_pfs.tsv", sep="\t", index=False)

    # Per-sample feature table
    sample = pd.concat(
        [
            prot_t_feat.add_prefix("protein_"),
            rna_t_feat.add_prefix("rna_"),
            pheno[
                [
                    c
                    for c in [
                        "CIBERSORT_T_cell_CD8+",
                        "ESTIMATE_ImmuneScore",
                        "ESTIMATE_StromalScore",
                        "xCell_T_cell_CD8+",
                        "xCell_immune_score",
                        "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                        "TMB",
                    ]
                    if c in pheno.columns
                ]
            ],
            surv,
            meta[["Age", "Sex", "Stage", "Tobacco_smoking_history", "TP53_mutation", "CDKN2A_mutation"]],
        ],
        axis=1,
    )
    sample.to_csv(TABLES / "sample_level_features.tsv", sep="\t")

    # Figures
    plot_tumor_normal(
        prot_t_feat["TACSTD2"],
        prot_n_feat["TACSTD2"],
        "TACSTD2 protein (LSCC vs NAT)",
        "log2 reference-intensity protein",
        "tacstd2_protein_tumor_vs_nat.png",
    )
    plot_tumor_normal(
        prot_t_feat["CLDN4"],
        prot_n_feat["CLDN4"],
        "CLDN4 protein (LSCC vs NAT)",
        "log2 reference-intensity protein",
        "cldn4_protein_tumor_vs_nat.png",
    )
    plot_tumor_normal(
        rna_t_feat["TACSTD2"],
        rna_n_feat["TACSTD2"],
        "TACSTD2 RNA (LSCC vs NAT)",
        "log2 RSEM UQ (coding)",
        "tacstd2_rna_tumor_vs_nat.png",
    )
    plot_tumor_normal(
        rna_t_feat["CLDN4"],
        rna_n_feat["CLDN4"],
        "CLDN4 RNA (LSCC vs NAT)",
        "log2 RSEM UQ (coding)",
        "cldn4_rna_tumor_vs_nat.png",
    )
    plot_scatter(
        rna_t_feat["TACSTD2"],
        prot_t_feat["TACSTD2"],
        "TACSTD2 RNA (log2 RSEM UQ)",
        "TACSTD2 protein (log2 intensity)",
        "TACSTD2 protein vs RNA",
        "tacstd2_protein_vs_rna.png",
    )
    plot_scatter(
        rna_t_feat["CLDN4"],
        prot_t_feat["CLDN4"],
        "CLDN4 RNA (log2 RSEM UQ)",
        "CLDN4 protein (log2 intensity)",
        "CLDN4 protein vs RNA",
        "cldn4_protein_vs_rna.png",
    )
    plot_scatter(
        prot_t_feat["TACSTD2"],
        prot_t_feat["CLDN4"],
        "TACSTD2 protein",
        "CLDN4 protein",
        "TACSTD2 vs CLDN4 protein",
        "tacstd2_vs_cldn4_protein.png",
    )
    plot_scatter(
        rna_t_feat["TACSTD2"],
        rna_t_feat["CLDN4"],
        "TACSTD2 RNA",
        "CLDN4 RNA",
        "TACSTD2 vs CLDN4 RNA",
        "tacstd2_vs_cldn4_rna.png",
    )
    plot_scatter(
        prot_t_feat["TACSTD2"],
        pheno["ESTIMATE_ImmuneScore"],
        "TACSTD2 protein",
        "ESTIMATE ImmuneScore",
        "TACSTD2 protein vs ESTIMATE ImmuneScore",
        "tacstd2_protein_vs_estimate_immune.png",
    )
    plot_scatter(
        prot_t_feat["TACSTD2"],
        pheno["CIBERSORT_T_cell_CD8+"],
        "TACSTD2 protein",
        "CIBERSORT CD8 T cells",
        "TACSTD2 protein vs CIBERSORT CD8",
        "tacstd2_protein_vs_cibersort_cd8.png",
    )
    plot_scatter(
        prot_t_feat["CLDN4"],
        pheno["ESTIMATE_ImmuneScore"],
        "CLDN4 protein",
        "ESTIMATE ImmuneScore",
        "CLDN4 protein vs ESTIMATE ImmuneScore",
        "cldn4_protein_vs_estimate_immune.png",
    )
    plot_scatter(
        rna_t_feat["TACSTD2"],
        sig_df["GEP_Tcell_inflamed"],
        "TACSTD2 RNA",
        "GEP T-cell inflamed (RNA)",
        "TACSTD2 RNA vs GEP",
        "tacstd2_rna_vs_gep.png",
    )
    plot_heatmap(corr_prot_pheno, "heatmap_protein_vs_phenotype.png", "Protein vs freeze immune scores (ρ)")
    plot_heatmap(corr_rna_pheno, "heatmap_rna_vs_phenotype.png", "RNA vs freeze immune scores (ρ)")
    plot_volcano(corr_prot_pheno, "volcano_protein_vs_phenotype.png", "Protein vs immune phenotype")
    plot_volcano(corr_rna_pheno, "volcano_rna_vs_phenotype.png", "RNA vs immune phenotype")
    plot_km(
        prot_t_feat["TACSTD2"],
        surv,
        "OS_days",
        "OS_event",
        "OS by TACSTD2 protein median (treatment-naive, not ICI)",
        "km_os_tacstd2_protein.png",
    )
    plot_km(
        prot_t_feat["CLDN4"],
        surv,
        "OS_days",
        "OS_event",
        "OS by CLDN4 protein median (treatment-naive, not ICI)",
        "km_os_cldn4_protein.png",
    )
    plot_km(
        rna_t_feat["TACSTD2"],
        surv,
        "OS_days",
        "OS_event",
        "OS by TACSTD2 RNA median (treatment-naive, not ICI)",
        "km_os_tacstd2_rna.png",
    )

    # Cohort summary JSON for the writeup
    summary = {
        "cohort": "CPTAC LSCC (LUSC synonym)",
        "citation": "Satpathy et al. Cell 2021 PMID 34358469",
        "treatment": "treatment-naive surgical resections; no prior chemo/RT; no ICI labels",
        "n_protein_tumor": int(protein_t.shape[1]),
        "n_protein_normal": int(protein_n.shape[1]),
        "n_rna_tumor": int(rna_t.shape[1]),
        "n_rna_normal": int(rna_n.shape[1]),
        "n_aligned_tumor_cases": int(len(cases)),
        "n_protein_genes": int(protein_t.shape[0]),
        "n_rna_genes": int(rna_t.shape[0]),
        "targets": {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"},
        "missingness": miss_df.to_dict(orient="records"),
        "signature_coverage_rna": sig_cov,
        "signature_coverage_protein": prot_sig_cov,
        "n_os_events": int(surv["OS_event"].sum()) if "OS_event" in surv else None,
        "n_pfs_events": int(surv["PFS_event"].sum()) if "PFS_event" in surv else None,
        "meta_sex": meta["Sex"].value_counts(dropna=False).to_dict() if "Sex" in meta else {},
        "meta_stage": meta["Stage"].value_counts(dropna=False).to_dict() if "Stage" in meta else {},
    }
    (TABLES / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    print(json.dumps(summary, indent=2, default=str))
    print("done", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
