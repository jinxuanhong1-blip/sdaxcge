#!/usr/bin/env python3
"""GSE135222: TACSTD2 (TROP2) and CLDN4 expression vs durable clinical benefit (DCB).

Cohort: 27 advanced NSCLC patients treated with anti-PD-1/PD-L1
(Jung et al., Nat Commun 2019; GEO GSE135222). Expression is RSEM TPM (hg19)
from the GEO supplementary file.

DCB definition: GEO provides PFS event (1 = progression/event, 0 = censored)
and pfs.time in days, but no explicit DCB label. We use the standard
definition DCB = PFS >= 183 days (6 months), NDB = PFS < 183 days.
This is unambiguous in this cohort because every censored patient has
follow-up >= 205 days (i.e., no patient is censored before 183 days).

Statistics: two-sided Mann-Whitney U on log2(TPM+1); AUC (DCB as positive
class) with 2000-sample stratified bootstrap 95% CI; BH-FDR across the two
genes. An exploratory combined score (mean of per-gene z-scores) is reported
separately and was not part of the two-gene primary comparison.

Usage: python3 scripts/gse135222_tacstd2_cldn4_dcb.py
Outputs to results/w200/GSE135222/.
"""

import gzip
import io
import re
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
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

GENES = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
}
DCB_CUTOFF_DAYS = 183  # 6 months
N_BOOT = 2000
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
            "gsm": gsm,
            "sex": meta["gender"],
            "age": [int(x) for x in meta["age"]],
            "pfs_event": [int(x) for x in meta["progression-free survival (pfs)"]],
            "pfs_days": [float(x) for x in meta["pfs.time"]],
        }
    )
    # Sanity check that DCB is unambiguous: no censoring before the cutoff.
    censored_early = df[(df.pfs_event == 0) & (df.pfs_days < DCB_CUTOFF_DAYS)]
    assert censored_early.empty, "Patient censored before 183 d: DCB ambiguous"
    df["DCB"] = (df.pfs_days >= DCB_CUTOFF_DAYS).astype(int)
    return df


def load_expression() -> pd.DataFrame:
    path = fetch(EXP_URL, DATA_DIR / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")
    exp = pd.read_csv(path, sep="\t", index_col=0)
    exp.index = [re.sub(r"\..*$", "", g) for g in exp.index]
    out = {}
    for name, ensg in GENES.items():
        sub = exp.loc[[ensg]]
        assert len(sub) == 1, f"{name}/{ensg}: expected 1 row, got {len(sub)}"
        out[f"{name}_TPM"] = sub.iloc[0]
    return pd.DataFrame(out)


def auc_mannwhitney(values: np.ndarray, labels: np.ndarray) -> float:
    """AUC for DCB (labels==1) as positive class, via U statistic."""
    pos, neg = values[labels == 1], values[labels == 0]
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return u / (len(pos) * len(neg))


def bootstrap_auc_ci(values, labels, n_boot=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    pos_idx, neg_idx = np.where(labels == 1)[0], np.where(labels == 0)[0]
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


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    clin = load_clinical()
    expr = load_expression()
    df = clin.merge(expr, left_on="sample", right_index=True, how="inner")
    assert len(df) == 27, f"expected 27 samples, got {len(df)}"

    for name in GENES:
        df[f"{name}_log2TPM1"] = np.log2(df[f"{name}_TPM"] + 1)

    # Exploratory combined score: mean of per-gene z-scores of log2(TPM+1).
    z = {n: stats.zscore(df[f"{n}_log2TPM1"], ddof=1) for n in GENES}
    df["combined_z"] = np.mean(list(z.values()), axis=0)

    df = df.sort_values("pfs_days", ascending=False)
    df.to_csv(OUT_DIR / "sample_data.csv", index=False)

    labels = df["DCB"].to_numpy()
    n_dcb, n_ndb = int(labels.sum()), int((1 - labels).sum())

    rows = []
    features = [(n, f"{n}_log2TPM1") for n in GENES] + [
        ("combined_z (exploratory)", "combined_z")
    ]
    for label, col in features:
        v = df[col].to_numpy()
        mwu = stats.mannwhitneyu(v[labels == 1], v[labels == 0], alternative="two-sided")
        auc = auc_mannwhitney(v, labels)
        lo, hi = bootstrap_auc_ci(v, labels)
        rows.append(
            {
                "feature": label,
                "n_DCB": n_dcb,
                "n_NDB": n_ndb,
                "median_DCB": np.median(v[labels == 1]),
                "median_NDB": np.median(v[labels == 0]),
                "delta_median": np.median(v[labels == 1]) - np.median(v[labels == 0]),
                "mannwhitney_U": mwu.statistic,
                "p_two_sided": mwu.pvalue,
                "AUC_DCB": auc,
                "AUC_95CI_lo": lo,
                "AUC_95CI_hi": hi,
            }
        )
    res = pd.DataFrame(rows)
    # BH-FDR across the two primary genes only.
    p = res.loc[:1, "p_two_sided"].to_numpy()
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * 2 / np.arange(1, 3))[::-1])[::-1]
    res.loc[:1, "q_BH_2genes"] = pd.Series(np.clip(q, 0, 1), index=order).sort_index().values
    res.to_csv(OUT_DIR / "stats.csv", index=False)

    # Spearman with PFS time (secondary, ignores censoring; interpret with care).
    sp_rows = []
    for label, col in features:
        r, pv = stats.spearmanr(df[col], df["pfs_days"])
        sp_rows.append({"feature": label, "spearman_r_vs_pfs_days": r, "p": pv})
    pd.DataFrame(sp_rows).to_csv(OUT_DIR / "spearman_pfs.csv", index=False)

    # ---- Plots ----
    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    for ax, (label, col) in zip(axes, features):
        groups = [df.loc[df.DCB == g, col].to_numpy() for g in (0, 1)]
        ax.boxplot(
            groups, tick_labels=[f"NDB\n(n={n_ndb})", f"DCB\n(n={n_dcb})"], showfliers=False
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
        row = res[res.feature == label].iloc[0]
        ax.set_title(f"{label}\nMWU p={row.p_two_sided:.3f}, AUC={row.AUC_DCB:.2f}")
        ax.set_ylabel("log2(TPM+1)" if col != "combined_z" else "mean z-score")
    fig.suptitle("GSE135222 (NSCLC, anti-PD-1/PD-L1): expression vs DCB (PFS \u2265 6 mo)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "boxplots_dcb.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 5))
    for label, col in features:
        v = df[col].to_numpy()
        thr = np.unique(v)
        tpr = [np.mean(v[labels == 1] >= t) for t in thr[::-1]]
        fpr = [np.mean(v[labels == 0] >= t) for t in thr[::-1]]
        row = res[res.feature == label].iloc[0]
        ax.plot([0] + fpr + [1], [0] + tpr + [1], label=f"{label} (AUC={row.AUC_DCB:.2f})")
    ax.plot([0, 1], [0, 1], "k--", lw=0.8)
    ax.set_xlabel("False positive rate")
    ax.set_ylabel("True positive rate")
    ax.set_title("ROC: DCB as positive class")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "roc_dcb.png", dpi=200)
    plt.close(fig)

    print(res.to_string(index=False))
    print()
    print(pd.DataFrame(sp_rows).to_string(index=False))


if __name__ == "__main__":
    main()
