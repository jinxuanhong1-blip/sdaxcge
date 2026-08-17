#!/usr/bin/env python3
"""GSE135222 leftover: CLDN4 vs CD8A, CD274, and ICI response.

CLDN4 only. Patient is the unit. Public GEO RNA-seq TPM + Kim 2020 author DCB.

Cohort: 27 advanced NSCLC patients treated with anti-PD-1/PD-L1
(Jung et al., Nat Commun 2019 PMID 31537801; RNA-seq slice of the SMC
cohort also in Kim et al., Clin Epigenetics 2020 PMID 32762727).

Usage: python3 methods/gse135222_cldn4/analyze.py
"""

from __future__ import annotations

import gzip
import json
import re
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"
TABLES = HERE / "tables"
PROC = HERE / "processed"
FIGS = HERE / "figures"

EXP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/"
    "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"
)
MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/"
    "GSE135222_series_matrix.txt.gz"
)

# Author-assigned DCB from Kim 2020 Table 1 (RNA-seq patients only).
# Y=1 / N=0, keyed by integer patient ID in the GEO title (NSCLC####).
AUTHOR_DCB = {
    1554: 0,
    378: 1,
    1352: 1,
    1809: 0,
    1327: 1,
    947: 1,
    1510: 1,
    1358: 1,
    1528: 1,
    1203: 0,
    990: 0,
    1155: 0,
    573: 0,
    825: 0,
    1508: 0,
    1066: 0,
    1412: 0,
    1164: 0,
    1017: 0,
    1425: 0,
    1145: 0,
    1708: 1,
    1104: 0,
    1619: 0,
    1873: 0,
    1079: 0,
    1401: 0,
}

GENES = {
    "CLDN4": "ENSG00000189143",
    "CD8A": "ENSG00000153563",
    "CD274": "ENSG00000120217",
}

DCB_CUTOFF_DAYS = 183
SEED = 20260817


def fetch(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.exists():
        urllib.request.urlretrieve(url, dest)
    return dest


def load_clinical() -> pd.DataFrame:
    path = fetch(MATRIX_URL, CACHE / "GSE135222_series_matrix.txt.gz")
    lines = gzip.open(path, "rt").read().splitlines()

    def row(key: str) -> list[str]:
        for ln in lines:
            if ln.startswith(key):
                return [f.strip('"') for f in ln.split("\t")[1:]]
        raise KeyError(key)

    titles = row("!Sample_title")
    gsm = row("!Sample_geo_accession")
    char_rows = [ln for ln in lines if ln.startswith("!Sample_characteristics_ch1")]
    meta: dict[str, list[str]] = {}
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
    if not early.empty:
        raise AssertionError("Patient censored before 183 d: PFS-derived DCB ambiguous")
    df["DCB_pfs183"] = (df.pfs_days >= DCB_CUTOFF_DAYS).astype(int)
    df["DCB_author"] = df["patient_id"].map(AUTHOR_DCB)
    if df["DCB_author"].isna().any():
        missing = df.loc[df["DCB_author"].isna(), "sample"].tolist()
        raise AssertionError(f"missing author DCB for {missing}")
    df["DCB_author"] = df["DCB_author"].astype(int)
    return df


def load_expression() -> pd.DataFrame:
    path = fetch(EXP_URL, CACHE / "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")
    exp = pd.read_csv(path, sep="\t", index_col=0)
    exp.index = [re.sub(r"\..*$", "", g) for g in exp.index]
    out = {}
    for name, ensg in GENES.items():
        if ensg not in exp.index:
            raise KeyError(f"{name}/{ensg} not in GSE135222 matrix")
        sub = exp.loc[[ensg]]
        if len(sub) != 1:
            raise AssertionError(f"{name}/{ensg}: expected 1 row, got {len(sub)}")
        out[f"{name}_TPM"] = sub.iloc[0]
    return pd.DataFrame(out)


def cliff_delta(pos: np.ndarray, neg: np.ndarray) -> float:
    if len(pos) == 0 or len(neg) == 0:
        return float("nan")
    u = stats.mannwhitneyu(pos, neg, alternative="two-sided").statistic
    return float(2.0 * u / (len(pos) * len(neg)) - 1.0)


def spearman_row(x: pd.Series, y: pd.Series, contrast: str) -> dict:
    mask = x.notna() & y.notna()
    n = int(mask.sum())
    if n < 3:
        return {
            "contrast": contrast,
            "n": n,
            "metric": "Spearman_rho",
            "effect": np.nan,
            "p": np.nan,
            "note": "n<3",
        }
    r, p = stats.spearmanr(x[mask], y[mask])
    return {
        "contrast": contrast,
        "n": n,
        "metric": "Spearman_rho",
        "effect": float(r),
        "p": float(p),
        "note": "log2(TPM+1)",
    }


def response_row(values: pd.Series, labels: pd.Series, contrast: str, endpoint: str) -> dict:
    work = pd.DataFrame({"v": values, "y": labels}).dropna()
    pos = work.loc[work.y == 1, "v"].to_numpy()
    neg = work.loc[work.y == 0, "v"].to_numpy()
    n_pos, n_neg = int(len(pos)), int(len(neg))
    n = n_pos + n_neg
    if n_pos == 0 or n_neg == 0:
        return {
            "contrast": contrast,
            "n": n,
            "n_pos": n_pos,
            "n_neg": n_neg,
            "metric": "Cliff_delta_pos_minus_neg",
            "effect": np.nan,
            "p": np.nan,
            "AUC": np.nan,
            "median_pos": np.nan,
            "median_neg": np.nan,
            "endpoint": endpoint,
            "note": "one class empty",
        }
    mwu = stats.mannwhitneyu(pos, neg, alternative="two-sided")
    auc = float(mwu.statistic / (n_pos * n_neg))
    return {
        "contrast": contrast,
        "n": n,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "metric": "Cliff_delta_pos_minus_neg",
        "effect": cliff_delta(pos, neg),
        "p": float(mwu.pvalue),
        "AUC": auc,
        "median_pos": float(np.median(pos)),
        "median_neg": float(np.median(neg)),
        "endpoint": endpoint,
        "note": "log2(TPM+1); pos=DCB, neg=NDB",
    }


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    PROC.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    clin = load_clinical()
    expr = load_expression()
    df = clin.merge(expr, left_on="sample", right_index=True, how="inner")
    if len(df) != 27:
        raise AssertionError(f"expected 27 RNA-seq patients, got {len(df)}")

    present = {name: f"{name}_TPM" in df.columns for name in GENES}
    for name in GENES:
        df[f"{name}_log2TPM1"] = np.log2(df[f"{name}_TPM"] + 1.0)

    patient_cols = [
        "sample",
        "patient_id",
        "gsm",
        "sex",
        "age",
        "pfs_event",
        "pfs_days",
        "DCB_author",
        "DCB_pfs183",
        "CLDN4_TPM",
        "CD8A_TPM",
        "CD274_TPM",
        "CLDN4_log2TPM1",
        "CD8A_log2TPM1",
        "CD274_log2TPM1",
    ]
    patient = df[patient_cols].sort_values("patient_id")
    patient.to_csv(PROC / "GSE135222_patient.tsv", sep="\t", index=False)

    tests = []
    tests.append(spearman_row(df["CLDN4_log2TPM1"], df["CD8A_log2TPM1"], "CLDN4 vs CD8A"))
    tests.append(spearman_row(df["CLDN4_log2TPM1"], df["CD274_log2TPM1"], "CLDN4 vs CD274"))
    tests.append(
        response_row(
            df["CLDN4_log2TPM1"],
            df["DCB_author"],
            "CLDN4 vs author DCB",
            "Kim2020_benefit_Y_N",
        )
    )
    tests.append(
        response_row(
            df["CLDN4_log2TPM1"],
            df["DCB_pfs183"],
            "CLDN4 vs PFS>=183d (sensitivity)",
            "GEO_pfs.time_ge_183",
        )
    )
    tests_df = pd.DataFrame(tests)
    tests_df.to_csv(TABLES / "leftover_tests.tsv", sep="\t", index=False)

    n_dcb = int((df.DCB_author == 1).sum())
    n_ndb = int((df.DCB_author == 0).sum())
    n_dcb_pfs = int((df.DCB_pfs183 == 1).sum())
    n_ndb_pfs = int((df.DCB_pfs183 == 0).sum())
    n_disagree = int((df.DCB_author != df.DCB_pfs183).sum())

    inventory = pd.DataFrame(
        [
            {
                "cohort": "GSE135222",
                "tissue": "advanced NSCLC tumor RNA-seq (RSEM TPM, hg19)",
                "treatment": "anti-PD-1/PD-L1",
                "n_patients": int(len(df)),
                "CLDN4_present": True,
                "CD8A_present": present["CD8A"],
                "CD274_present": present["CD274"],
                "response_present": True,
                "response_source": "Kim 2020 Table 1 benefit Y/N (not a GEO characteristic)",
                "n_DCB_author": n_dcb,
                "n_NDB_author": n_ndb,
                "n_DCB_pfs183": n_dcb_pfs,
                "n_NDB_pfs183": n_ndb_pfs,
                "n_label_disagree_author_vs_pfs183": n_disagree,
                "unit": "patient",
                "gene_scope": "CLDN4 only",
            }
        ]
    )
    inventory.to_csv(TABLES / "inventory.tsv", sep="\t", index=False)

    primary = tests_df[tests_df.contrast.isin(["CLDN4 vs CD8A", "CLDN4 vs CD274", "CLDN4 vs author DCB"])].copy()
    primary.to_csv(TABLES / "primary.tsv", sep="\t", index=False)

    summary = {
        "cohort": "GSE135222",
        "n": int(len(df)),
        "n_DCB_author": n_dcb,
        "n_NDB_author": n_ndb,
        "CLDN4_vs_CD8A_rho": float(tests_df.loc[tests_df.contrast == "CLDN4 vs CD8A", "effect"].iloc[0]),
        "CLDN4_vs_CD8A_p": float(tests_df.loc[tests_df.contrast == "CLDN4 vs CD8A", "p"].iloc[0]),
        "CLDN4_vs_CD274_rho": float(tests_df.loc[tests_df.contrast == "CLDN4 vs CD274", "effect"].iloc[0]),
        "CLDN4_vs_CD274_p": float(tests_df.loc[tests_df.contrast == "CLDN4 vs CD274", "p"].iloc[0]),
        "CLDN4_vs_authorDCB_cliff": float(tests_df.loc[tests_df.contrast == "CLDN4 vs author DCB", "effect"].iloc[0]),
        "CLDN4_vs_authorDCB_p": float(tests_df.loc[tests_df.contrast == "CLDN4 vs author DCB", "p"].iloc[0]),
        "CLDN4_vs_authorDCB_AUC": float(tests_df.loc[tests_df.contrast == "CLDN4 vs author DCB", "AUC"].iloc[0]),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    rng = np.random.default_rng(SEED)
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.8))

    ax = axes[0]
    ax.scatter(df["CLDN4_log2TPM1"], df["CD8A_log2TPM1"], s=28, c="#1f77b4", alpha=0.85)
    r = summary["CLDN4_vs_CD8A_rho"]
    p = summary["CLDN4_vs_CD8A_p"]
    ax.set_xlabel("CLDN4 log2(TPM+1)")
    ax.set_ylabel("CD8A log2(TPM+1)")
    ax.set_title(f"CLDN4 vs CD8A\nρ={r:.2f}, p={p:.2f}, n=27")

    ax = axes[1]
    ax.scatter(df["CLDN4_log2TPM1"], df["CD274_log2TPM1"], s=28, c="#1f77b4", alpha=0.85)
    r = summary["CLDN4_vs_CD274_rho"]
    p = summary["CLDN4_vs_CD274_p"]
    ax.set_xlabel("CLDN4 log2(TPM+1)")
    ax.set_ylabel("CD274 log2(TPM+1)")
    ax.set_title(f"CLDN4 vs CD274\nρ={r:.2f}, p={p:.2f}, n=27")

    ax = axes[2]
    groups = [
        df.loc[df.DCB_author == 0, "CLDN4_log2TPM1"].to_numpy(),
        df.loc[df.DCB_author == 1, "CLDN4_log2TPM1"].to_numpy(),
    ]
    ax.boxplot(groups, tick_labels=[f"NDB\n(n={n_ndb})", f"DCB\n(n={n_dcb})"], showfliers=False)
    for i, g in enumerate(groups):
        ax.scatter(
            np.full(len(g), i + 1) + rng.uniform(-0.08, 0.08, len(g)),
            g,
            s=22,
            alpha=0.85,
            color="#1f77b4" if i == 0 else "#d62728",
            zorder=3,
        )
    ax.set_ylabel("CLDN4 log2(TPM+1)")
    ax.set_title(
        f"CLDN4 vs author DCB\nδ={summary['CLDN4_vs_authorDCB_cliff']:.2f}, "
        f"p={summary['CLDN4_vs_authorDCB_p']:.2f}"
    )

    fig.suptitle("GSE135222 leftover — CLDN4 only (NSCLC anti-PD-1/PD-L1)")
    fig.tight_layout()
    fig.savefig(FIGS / "CLDN4_vs_CD8A_CD274_DCB.png", dpi=200)
    plt.close(fig)

    print(tests_df.to_string(index=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
