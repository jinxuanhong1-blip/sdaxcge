#!/usr/bin/env python3
"""GSE135222: TACSTD2 (TROP2) and CLDN4 vs durable clinical benefit (DCB).

Cohort: 27 advanced NSCLC patients treated with anti-PD-1/PD-L1
(Jung et al., Nat Commun 2019 PMID 31537801; RNA-seq subset of the SMC
cohort also reported in Kim et al., Clin Epigenetics 2020 PMID 32762727).
Expression is RSEM TPM (hg19) from the GEO supplementary file.

Primary DCB labels are the author-assigned benefit column in Kim 2020
Supplementary Table 1 (Y/N), matched by patient ID. Jung/Kim define DCB as
RECIST v1.1 PR or SD lasting >6 months (not a pure PFS-time cutoff).
PFS ≥ 183 days disagrees with the authors on one RNA-seq patient
(NSCLC1708: PFS 168 d, event, author DCB=Y) and is kept only as sensitivity.

Primary tests: two-sided Mann-Whitney U and AUC on log2(TPM+1) for TACSTD2
and CLDN4 vs author DCB; BH-FDR across those two genes. Survival (Cox PH
on continuous expression; log-rank on median split) does not use the DCB
cutoff. CD8A/CXCL9/GZMB/CD274/IFNG are positive-control ICI genes, not
part of the primary family.

Usage: python3 scripts/gse135222_tacstd2_cldn4_dcb.py
Outputs to results/w200/GSE135222/.
"""

from __future__ import annotations

import gzip
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from scipy import stats

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "GSE135222"
OUT_DIR = ROOT / "results" / "w200" / "GSE135222"

EXP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/"
    "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/"
    "GSE135222_series_matrix.txt.gz"
)
KIM_XLSX_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1186%2Fs13148-020-00907-4/MediaObjects/"
    "13148_2020_907_MOESM6_ESM.xlsx"
)

# Author-assigned DCB from Kim 2020 Table 1 (RNA-seq patients only).
# Used if the xlsx download fails; also used to verify a successful download.
AUTHOR_DCB = {
    1554: 0, 378: 1, 1352: 1, 1809: 0, 1327: 1, 947: 1, 1510: 1,
    1358: 1, 1528: 1, 1203: 0, 990: 0, 1155: 0, 573: 0, 825: 0,
    1508: 0, 1066: 0, 1412: 0, 1164: 0, 1017: 0, 1425: 0, 1145: 0,
    1708: 1, 1104: 0, 1619: 0, 1873: 0, 1079: 0, 1401: 0,
}

PRIMARY = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}
CONTROLS = {
    "CD8A": "ENSG00000153563",
    "CXCL9": "ENSG00000138755",
    "GZMB": "ENSG00000100453",
    "CD274": "ENSG00000120217",
    "IFNG": "ENSG00000111537",
}
ALL_GENES = {**PRIMARY, **CONTROLS}

DCB_CUTOFF_DAYS = 183
N_BOOT = 2000
N_SIM = 2000
SEED = 20260816


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def load_clinical() -> pd.DataFrame:
    path = fetch(MATRIX_URL, DATA_DIR / "GSE135222_series_matrix.txt.gz")
    lines = gzip.open(path, "rt").read().splitlines()

    def row(key):
        for ln in lines:
            if ln.startswith(key):
                return [f.strip('"') for f in ln.split("\t")[1:]]
        raise KeyError(key)

    titles = row("!Sample_title")
    gsm = row("!Sample_geo_accession")
    char_rows = [ln for ln in lines if ln.startswith("!Sample_characteristics_ch1")]
    meta = {}
    for ln in char_rows:
        fields = [f.strip('"') for f in ln.split("\t")[1:]]
        key = fields[0].split(":")[0].strip()
        meta[key] = [f.split(":", 1)[1].strip() for f in fields]

    df = pd.DataFrame(
        {
            "sample": [t.replace("NSCLC ", "NSCLC") for t in titles],
            "patient_id": [int(t.replace("NSCLC ", "")) for t in titles],
            "gsm": gsm,
            "sex": meta["gender"],
            "age": [int(x) for x in meta["age"]],
            "pfs_event": [int(x) for x in meta["progression-free survival (pfs)"]],
            "pfs_days": [float(x) for x in meta["pfs.time"]],
        }
    )
    early = df[(df.pfs_event == 0) & (df.pfs_days < DCB_CUTOFF_DAYS)]
    assert early.empty, "Patient censored before 183 d: PFS-derived DCB ambiguous"
    df["DCB_pfs183"] = (df.pfs_days >= DCB_CUTOFF_DAYS).astype(int)
    return df


def load_author_dcb() -> pd.Series:
    """Kim 2020 Table 1 benefit (Y=1/N=0), keyed by integer patient ID."""
    hardcoded = pd.Series(AUTHOR_DCB, name="DCB_author")
    try:
        path = fetch(KIM_XLSX_URL, DATA_DIR / "kim2020_13148_2020_907_MOESM6.xlsx")
        raw = pd.read_excel(path, sheet_name="Table1", header=None)
        # Row 2 is the header in 0-based indexing of the printed table.
        tab = raw.iloc[3:].copy()
        tab.columns = ["patient_id", "PFS", "PFS_event", "benefit"]
        tab = tab.dropna(subset=["patient_id"])
        tab["patient_id"] = tab["patient_id"].astype(int)
        tab["DCB_author"] = tab["benefit"].map({"Y": 1, "N": 0})
        series = tab.set_index("patient_id")["DCB_author"]
        overlap = series.index.intersection(hardcoded.index)
        assert (series.loc[overlap] == hardcoded.loc[overlap]).all(), (
            "Downloaded Kim 2020 Table 1 does not match hardcoded author DCB"
        )
        return series
    except Exception as exc:
        print(f"WARNING: using hardcoded author DCB ({exc})")
        return hardcoded


def load_expression() -> pd.DataFrame:
    path = fetch(EXP_URL, DATA_DIR / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")
    exp = pd.read_csv(path, sep="\t", index_col=0)
    exp.index = [re.sub(r"\..*$", "", g) for g in exp.index]
    out = {}
    for name, ensg in ALL_GENES.items():
        sub = exp.loc[[ensg]]
        assert len(sub) == 1, f"{name}/{ensg}: expected 1 row, got {len(sub)}"
        out[f"{name}_TPM"] = sub.iloc[0]
    return pd.DataFrame(out)


def auc_mannwhitney(values: np.ndarray, labels: np.ndarray) -> float:
    pos, neg = values[labels == 1], values[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return np.nan
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(u / (len(pos) * len(neg)))


def bootstrap_auc_ci(values, labels, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    pos_idx, neg_idx = np.where(labels == 1)[0], np.where(labels == 0)[0]
    if len(pos_idx) < 2 or len(neg_idx) < 2:
        return np.array([np.nan, np.nan])
    aucs = []
    for _ in range(n_boot):
        idx = np.concatenate(
            [
                rng.choice(pos_idx, len(pos_idx), replace=True),
                rng.choice(neg_idx, len(neg_idx), replace=True),
            ]
        )
        aucs.append(auc_mannwhitney(values[idx], labels[idx]))
    return np.percentile(aucs, [2.5, 97.5])


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def mwu_auc_row(feature, values, labels, role):
    pos, neg = values[labels == 1], values[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return {
            "feature": feature,
            "role": role,
            "n_DCB": int(labels.sum()),
            "n_NDB": int((1 - labels).sum()),
            "median_DCB": np.nan,
            "median_NDB": np.nan,
            "delta_median": np.nan,
            "mannwhitney_U": np.nan,
            "p_two_sided": np.nan,
            "AUC_DCB": np.nan,
            "AUC_95CI_lo": np.nan,
            "AUC_95CI_hi": np.nan,
        }
    mwu = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = auc_mannwhitney(values, labels)
    lo, hi = bootstrap_auc_ci(values, labels)
    return {
        "feature": feature,
        "role": role,
        "n_DCB": int(labels.sum()),
        "n_NDB": int((1 - labels).sum()),
        "median_DCB": float(np.median(pos)),
        "median_NDB": float(np.median(neg)),
        "delta_median": float(np.median(pos) - np.median(neg)),
        "mannwhitney_U": float(mwu.statistic),
        "p_two_sided": float(mwu.pvalue),
        "AUC_DCB": auc,
        "AUC_95CI_lo": float(lo),
        "AUC_95CI_hi": float(hi),
    }


def cox_and_logrank(df: pd.DataFrame, col: str, feature: str) -> dict:
    work = df[["pfs_days", "pfs_event", col]].dropna().copy()
    work["z"] = stats.zscore(work[col], ddof=1)
    cph = CoxPHFitter()
    cph.fit(work[["pfs_days", "pfs_event", "z"]], duration_col="pfs_days", event_col="pfs_event")
    s = cph.summary.loc["z"]
    high = work[col] >= work[col].median()
    lr = logrank_test(
        work.loc[high, "pfs_days"],
        work.loc[~high, "pfs_days"],
        event_observed_A=work.loc[high, "pfs_event"],
        event_observed_B=work.loc[~high, "pfs_event"],
    )
    return {
        "feature": feature,
        "n": len(work),
        "n_events": int(work.pfs_event.sum()),
        "cox_HR_per_SD": float(s["exp(coef)"]),
        "cox_HR_95CI_lo": float(s["exp(coef) lower 95%"]),
        "cox_HR_95CI_hi": float(s["exp(coef) upper 95%"]),
        "cox_p": float(s["p"]),
        "logrank_p_median_split": float(lr.p_value),
        "n_high": int(high.sum()),
        "n_low": int((~high).sum()),
    }


def power_mwu(n_dcb, n_ndb, true_auc, n_sim=N_SIM, seed=SEED) -> float:
    """Approximate power of two-sided MWU when both groups are normal.

    For equal-variance normals, AUC = Phi(delta / sqrt(2)).
    """
    if true_auc <= 0.5:
        return np.nan
    from scipy.special import ndtri

    delta = ndtri(true_auc) * np.sqrt(2.0)
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_sim):
        pos = rng.normal(delta, 1.0, n_dcb)
        neg = rng.normal(0.0, 1.0, n_ndb)
        p = stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue
        hits += p < 0.05
    return hits / n_sim


def roc_points(values, labels):
    order = np.argsort(-values)
    y = labels[order]
    tpr = np.concatenate([[0.0], np.cumsum(y) / y.sum(), [1.0]])
    fpr = np.concatenate([[0.0], np.cumsum(1 - y) / (1 - y).sum(), [1.0]])
    return fpr, tpr


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    clin = load_clinical()
    author = load_author_dcb()
    expr = load_expression()
    df = clin.merge(expr, left_on="sample", right_index=True, how="inner")
    assert len(df) == 27, f"expected 27 samples, got {len(df)}"
    df["DCB_author"] = df["patient_id"].map(author)
    assert df["DCB_author"].notna().all(), "missing author DCB for an RNA-seq sample"
    df["DCB_author"] = df["DCB_author"].astype(int)
    df["DCB"] = df["DCB_author"]  # primary label

    for name in ALL_GENES:
        df[f"{name}_log2TPM1"] = np.log2(df[f"{name}_TPM"] + 1)

    z_pri = {n: stats.zscore(df[f"{n}_log2TPM1"], ddof=1) for n in PRIMARY}
    df["combined_z"] = np.mean(list(z_pri.values()), axis=0)
    z_imm = {n: stats.zscore(df[f"{n}_log2TPM1"], ddof=1) for n in ("CD8A", "CXCL9", "GZMB")}
    df["immune_z"] = np.mean(list(z_imm.values()), axis=0)

    df = df.sort_values("pfs_days", ascending=False)
    df.to_csv(OUT_DIR / "sample_data.csv", index=False)

    # --- Author vs PFS-derived DCB ---
    disc = df.loc[
        df.DCB_author != df.DCB_pfs183,
        ["sample", "gsm", "pfs_days", "pfs_event", "DCB_author", "DCB_pfs183"],
    ]
    disc.to_csv(OUT_DIR / "author_vs_pfs183_discrepancy.csv", index=False)

    labels = df["DCB"].to_numpy()
    features = (
        [(n, f"{n}_log2TPM1", "primary") for n in PRIMARY]
        + [("combined_z (exploratory)", "combined_z", "exploratory")]
        + [(n, f"{n}_log2TPM1", "positive_control") for n in CONTROLS]
        + [("immune_z CD8A/CXCL9/GZMB (control)", "immune_z", "positive_control")]
    )

    rows = [mwu_auc_row(lab, df[col].to_numpy(), labels, role) for lab, col, role in features]
    res = pd.DataFrame(rows)
    pri = res["role"] == "primary"
    res.loc[pri, "q_BH_2genes"] = bh_fdr(res.loc[pri, "p_two_sided"].to_numpy())
    res.to_csv(OUT_DIR / "stats.csv", index=False)

    # --- Sensitivity: DCB cutoffs on PFS time ---
    sens = []
    for days, name in [(90, "PFS>=90d"), (168, "PFS>=168d (24wk)"), (183, "PFS>=183d"), (365, "PFS>=365d")]:
        lab = (df.pfs_days >= days).astype(int).to_numpy()
        for feat, col, role in features:
            if role == "positive_control":
                continue
            row = mwu_auc_row(feat, df[col].to_numpy(), lab, role)
            row["dcb_definition"] = name
            sens.append(row)
    # author labels as a row-block for comparison
    for feat, col, role in features:
        if role == "positive_control":
            continue
        row = mwu_auc_row(feat, df[col].to_numpy(), labels, role)
        row["dcb_definition"] = "author_Kim2020"
        sens.append(row)
    pd.DataFrame(sens).to_csv(OUT_DIR / "sensitivity_dcb_cutoffs.csv", index=False)

    # --- Survival ---
    surv_feats = [(n, f"{n}_log2TPM1") for n in list(PRIMARY) + list(CONTROLS)] + [
        ("combined_z", "combined_z"),
        ("immune_z", "immune_z"),
    ]
    surv = pd.DataFrame([cox_and_logrank(df, col, lab) for lab, col in surv_feats])
    surv.to_csv(OUT_DIR / "cox_km.csv", index=False)

    # --- Spearman vs PFS time (ignores censoring) ---
    sp_rows = []
    for lab, col, role in features:
        r, pv = stats.spearmanr(df[col], df["pfs_days"])
        sp_rows.append({"feature": lab, "role": role, "spearman_r_vs_pfs_days": r, "p": pv})
    r_tc, p_tc = stats.spearmanr(df["TACSTD2_log2TPM1"], df["CLDN4_log2TPM1"])
    sp_rows.append(
        {
            "feature": "TACSTD2 vs CLDN4",
            "role": "coexpression",
            "spearman_r_vs_pfs_days": r_tc,
            "p": p_tc,
        }
    )
    pd.DataFrame(sp_rows).to_csv(OUT_DIR / "spearman_pfs.csv", index=False)

    # --- Leave-one-out influence on primary genes ---
    loo = []
    for drop in df["sample"]:
        sub = df[df["sample"] != drop]
        lab = sub["DCB"].to_numpy()
        for name in PRIMARY:
            v = sub[f"{name}_log2TPM1"].to_numpy()
            row = mwu_auc_row(name, v, lab, "primary")
            row["left_out"] = drop
            loo.append(row)
    loo_df = pd.DataFrame(loo)
    loo_df.to_csv(OUT_DIR / "leave_one_out.csv", index=False)

    # --- Power ---
    n_dcb, n_ndb = int(labels.sum()), int((1 - labels).sum())
    pow_rows = []
    for auc in (0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90):
        pow_rows.append(
            {
                "n_DCB": n_dcb,
                "n_NDB": n_ndb,
                "true_AUC": auc,
                "power_mwu_alpha0.05": power_mwu(n_dcb, n_ndb, auc),
            }
        )
    pd.DataFrame(pow_rows).to_csv(OUT_DIR / "power.csv", index=False)

    # ---- Plots ----
    rng = np.random.default_rng(SEED)

    def box_panel(ax, values, labels_, title, ylabel):
        groups = [values[labels_ == g] for g in (0, 1)]
        ax.boxplot(
            groups,
            tick_labels=[f"NDB\n(n={n_ndb})", f"DCB\n(n={n_dcb})"],
            showfliers=False,
        )
        for i, g in enumerate(groups):
            ax.scatter(
                np.full(len(g), i + 1) + rng.uniform(-0.08, 0.08, len(g)),
                g,
                s=25,
                alpha=0.8,
                color="#1f77b4" if i == 0 else "#d62728",
                zorder=3,
            )
        ax.set_title(title)
        ax.set_ylabel(ylabel)

    plot_specs = [
        ("TACSTD2", "TACSTD2_log2TPM1", "log2(TPM+1)"),
        ("CLDN4", "CLDN4_log2TPM1", "log2(TPM+1)"),
        ("combined_z", "combined_z", "mean z-score"),
        ("immune_z (control)", "immune_z", "mean z-score"),
    ]
    fig, axes = plt.subplots(1, 4, figsize=(14, 4))
    for ax, (lab, col, ylab) in zip(axes, plot_specs):
        row = res[res.feature.str.startswith(lab.split()[0])].iloc[0]
        box_panel(
            ax,
            df[col].to_numpy(),
            labels,
            f"{lab}\nMWU p={row.p_two_sided:.3f}, AUC={row.AUC_DCB:.2f}",
            ylab,
        )
    fig.suptitle(
        "GSE135222 (NSCLC, anti-PD-1/PD-L1): expression vs author DCB (Kim 2020)"
    )
    fig.tight_layout()
    fig.savefig(OUT_DIR / "boxplots_dcb.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 5, figsize=(15, 3.6))
    for ax, name in zip(axes, CONTROLS):
        row = res[res.feature == name].iloc[0]
        box_panel(
            ax,
            df[f"{name}_log2TPM1"].to_numpy(),
            labels,
            f"{name}\nMWU p={row.p_two_sided:.3f}, AUC={row.AUC_DCB:.2f}",
            "log2(TPM+1)",
        )
    fig.suptitle("Positive-control ICI genes vs author DCB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "boxplots_controls.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    for lab, col, role in features:
        if role == "exploratory":
            continue
        fpr, tpr = roc_points(df[col].to_numpy(), labels)
        row = res[res.feature == lab].iloc[0]
        ls = "-" if role == "primary" else "--"
        ax.plot(fpr, tpr, ls=ls, label=f"{lab} (AUC={row.AUC_DCB:.2f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC: author DCB as positive class")
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "roc_dcb.png", dpi=200)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    km_cols = [
        ("TACSTD2", "TACSTD2_log2TPM1"),
        ("CLDN4", "CLDN4_log2TPM1"),
        ("immune_z", "immune_z"),
    ]
    for ax, (lab, col) in zip(axes, km_cols):
        high = df[col] >= df[col].median()
        km = KaplanMeierFitter()
        km.fit(df.loc[high, "pfs_days"], df.loc[high, "pfs_event"], label=f"{lab} high")
        km.plot_survival_function(ax=ax, ci_show=False, color="#d62728")
        km.fit(df.loc[~high, "pfs_days"], df.loc[~high, "pfs_event"], label=f"{lab} low")
        km.plot_survival_function(ax=ax, ci_show=False, color="#1f77b4")
        prow = surv[surv.feature == lab].iloc[0]
        ax.set_title(
            f"{lab} median split\nlog-rank p={prow.logrank_p_median_split:.3f}; "
            f"Cox HR/SD={prow.cox_HR_per_SD:.2f} p={prow.cox_p:.3f}"
        )
        ax.set_xlabel("PFS (days)")
        ax.set_ylabel("PFS probability")
    fig.suptitle("GSE135222 Kaplan–Meier (median split) and Cox (per 1 SD)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "km_curves.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(
        df["TACSTD2_log2TPM1"],
        df["CLDN4_log2TPM1"],
        c=np.where(labels == 1, "#d62728", "#1f77b4"),
        s=40,
        alpha=0.85,
    )
    for _, r in df.iterrows():
        if r.sample in {"NSCLC947", "NSCLC1708", "NSCLC990", "NSCLC1401"}:
            ax.annotate(r.sample.replace("NSCLC", ""), (r.TACSTD2_log2TPM1, r.CLDN4_log2TPM1), fontsize=7)
    ax.set_xlabel("TACSTD2 log2(TPM+1)")
    ax.set_ylabel("CLDN4 log2(TPM+1)")
    ax.set_title(f"TACSTD2 vs CLDN4 (ρ={r_tc:.2f}, p={p_tc:.3f})\nred=author DCB, blue=NDB")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "scatter_tacstd2_cldn4.png", dpi=200)
    plt.close(fig)

    print("=== Primary (author DCB) ===")
    print(res.to_string(index=False))
    print("\n=== Discrepancy vs PFS>=183 ===")
    print(disc.to_string(index=False) if len(disc) else "(none)")
    print("\n=== Cox / log-rank ===")
    print(surv.to_string(index=False))
    print("\n=== Power ===")
    print(pd.DataFrame(pow_rows).to_string(index=False))
    print(f"\nTACSTD2 vs CLDN4 Spearman r={r_tc:.3f} p={p_tc:.4f}")


if __name__ == "__main__":
    main()
