"""ICI response axis in public NSCLC anti-PD-1/PD-L1 RNA-seq cohorts.

Cohorts (both pre-treatment tumour RNA-seq under PD-1/PD-L1 blockade):
  * GSE126044 (Cho et al. 2020) n=16, binary responder / non-responder.
  * GSE135222 (Jung et al. 2019, Nat Commun) n=27, progression-free survival.

Neither cohort carries EGFR/ALK genotype, so these establish the *response* axis for
TACSTD2/CLDN4 in the broader NSCLC-ICI setting; the EGFR/ALK-mutant specificity is
carried by the TCGA analysis.  This is stated explicitly in the write-up.

Outputs:
  tables/ici_gse126044_target_response.csv
  tables/ici_gse135222_target_pfs.csv
  tables/ici_target_expression_<gse>.csv (per-sample, small)
  figures/ici_gse126044_response_box.png
  figures/ici_gse135222_km_<target>.png
"""
from __future__ import annotations

import re
import gzip
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

import common as C

TARGETS = list(C.TARGETS.keys())
IMMUNE = list(C.IMMUNE_GENES.keys())


# ---------------------------------------------------------------------------
# GEO series-matrix metadata parser
# ---------------------------------------------------------------------------
def parse_series_matrix(path: str) -> pd.DataFrame:
    titles, accs, chars = [], [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                accs = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                chars.append([x.strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    meta = pd.DataFrame({"title": titles, "gsm": accs})
    for row in chars:
        # each row is "key: value" repeated across samples; use first key
        key = row[0].split(":")[0].strip().lower().replace(" ", "_").replace(".", "")
        meta[key] = [v.split(":", 1)[1].strip() if ":" in v else v for v in row]
    return meta


def log2_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def sig_table(expr_samples_by_genes: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=expr_samples_by_genes.index)
    for name, genes in C.SIGNATURES.items():
        out[f"sig_{name}"] = C.signature_score(expr_samples_by_genes, genes)
    return out


# ---------------------------------------------------------------------------
# GSE126044 : responder vs non-responder
# ---------------------------------------------------------------------------
def run_gse126044() -> pd.DataFrame:
    counts = pd.read_csv(f"{C.DATA}/GSE126044_counts.txt.gz", sep="\t", index_col=0)
    logcpm = log2_cpm(counts)                      # genes x samples
    expr = logcpm.T                                # samples x genes

    meta = parse_series_matrix(f"{C.DATA}/GSE126044_series_matrix.txt.gz")
    meta["key"] = meta["title"].str.replace("RNA-seq_", "", regex=False)
    resp_map = dict(zip(meta["key"], meta["patient_response"]))
    expr = expr.loc[[s for s in expr.index if s in resp_map]]
    responder = expr.index.map(lambda s: 1 if resp_map[s] == "responder" else 0)

    sigs = sig_table(expr[[g for g in IMMUNE if g in expr.columns]])
    per = pd.concat([expr[TARGETS], sigs], axis=1)
    per["response"] = responder.values
    per.to_csv(f"{C.TABLES}/ici_target_expression_gse126044.csv")

    rows = []
    feats = TARGETS + [f"sig_{s}" for s in C.SIGNATURES]
    for f in feats:
        r = per.loc[per["response"] == 1, f].dropna()
        nr = per.loc[per["response"] == 0, f].dropna()
        u, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        auc = u / (len(r) * len(nr))               # AUC for responder>non-responder
        rows.append({"cohort": "GSE126044", "feature": f,
                     "n_resp": len(r), "n_nonresp": len(nr),
                     "median_resp": r.median(), "median_nonresp": nr.median(),
                     "mwu_U": u, "auc_resp_high": auc, "p": p})
    res = pd.DataFrame(rows)
    res["q_BH"] = C.benjamini_hochberg(res["p"].values)
    res.to_csv(f"{C.TABLES}/ici_gse126044_target_response.csv", index=False)

    # figure
    fig, axes = plt.subplots(1, len(TARGETS), figsize=(3.0 * len(TARGETS), 3.4))
    for ax, tgt in zip(np.atleast_1d(axes), TARGETS):
        data = [per.loc[per["response"] == 0, tgt].dropna().values,
                per.loc[per["response"] == 1, tgt].dropna().values]
        bp = ax.boxplot(data, tick_labels=["non-resp", "resp"], showfliers=False,
                        patch_artist=True, widths=0.6)
        for patch, col in zip(bp["boxes"], ["#d95f02", "#1b9e77"]):
            patch.set_facecolor(col); patch.set_alpha(0.6)
        for i, d in enumerate(data):
            x = np.random.normal(i + 1, 0.06, size=len(d))
            ax.scatter(x, d, s=14, color="k", alpha=0.5, zorder=3)
        pval = res.loc[res.feature == tgt, "p"].values[0]
        ax.set_title(f"{tgt}\nMWU p={pval:.3f}", fontsize=10)
        ax.set_ylabel("log2(CPM+1)", fontsize=8)
        ax.tick_params(labelsize=8)
    fig.suptitle("GSE126044: anti-PD-1 NSCLC (n=16)", fontsize=10)
    fig.tight_layout()
    fig.savefig(f"{C.FIGURES}/ici_gse126044_response_box.png", dpi=150)
    plt.close(fig)
    return res


# ---------------------------------------------------------------------------
# GSE135222 : progression-free survival
# ---------------------------------------------------------------------------
def run_gse135222() -> pd.DataFrame:
    expr_raw = pd.read_csv(f"{C.DATA}/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
                           sep="\t", index_col=0)
    # rows Ensembl (versioned), values TPM/FPKM; map to symbols for our genes
    ens2sym = {v: k for k, v in {**C.TARGETS, **C.IMMUNE_GENES}.items()}
    idx_nover = expr_raw.index.to_series().str.replace(r"\.\d+$", "", regex=True)
    keep = idx_nover.isin(ens2sym.keys())
    sub = expr_raw[keep.values].copy()
    sub.index = idx_nover[keep.values].map(ens2sym).values
    expr = np.log2(sub.T.astype(float) + 1.0)      # samples x symbols

    meta = parse_series_matrix(f"{C.DATA}/GSE135222_series_matrix.txt.gz")
    meta["key"] = meta["title"].str.replace(" ", "", regex=False)
    # find pfs event + time columns robustly
    pfs_col = [c for c in meta.columns if c.startswith("progression")][0]
    time_col = [c for c in meta.columns if "pfstime" in c][0]
    meta["pfs_event"] = meta[pfs_col].astype(int)
    meta["pfs_time"] = meta[time_col].astype(float)
    mmap = meta.set_index("key")[["pfs_event", "pfs_time"]]

    common_ids = [s for s in expr.index if s in mmap.index]
    expr = expr.loc[common_ids]
    surv = mmap.loc[common_ids]

    sigs = sig_table(expr[[g for g in IMMUNE if g in expr.columns]])
    per = pd.concat([expr[TARGETS], sigs, surv], axis=1)
    per.to_csv(f"{C.TABLES}/ici_target_expression_gse135222.csv")

    rows = []
    kmf = KaplanMeierFitter()
    for tgt in TARGETS + [f"sig_{s}" for s in C.SIGNATURES]:
        x = per[tgt].astype(float)
        # Spearman vs PFS time
        rho, prho = stats.spearmanr(x, per["pfs_time"])
        # median split KM + logrank
        med = x.median()
        high = x > med
        lr = logrank_test(per.loc[high, "pfs_time"], per.loc[~high, "pfs_time"],
                          per.loc[high, "pfs_event"], per.loc[~high, "pfs_event"])
        rows.append({"cohort": "GSE135222", "feature": tgt, "n": len(per),
                     "spearman_rho_vs_pfstime": rho, "spearman_p": prho,
                     "logrank_p_medsplit": lr.p_value,
                     "median_pfs_high": per.loc[high, "pfs_time"].median(),
                     "median_pfs_low": per.loc[~high, "pfs_time"].median()})
        if tgt in TARGETS:
            fig, ax = plt.subplots(figsize=(4, 3.4))
            for grp, lab, col in [(high, f"{tgt} high", "#d95f02"),
                                  (~high, f"{tgt} low", "#1b9e77")]:
                kmf.fit(per.loc[grp, "pfs_time"], per.loc[grp, "pfs_event"], label=lab)
                kmf.plot_survival_function(ax=ax, ci_show=False, color=col)
            ax.set_title(f"GSE135222 PFS by {tgt}\nlogrank p={lr.p_value:.3f}", fontsize=10)
            ax.set_xlabel("PFS time (days)", fontsize=9)
            ax.set_ylabel("PFS probability", fontsize=9)
            ax.tick_params(labelsize=8)
            fig.tight_layout()
            fig.savefig(f"{C.FIGURES}/ici_gse135222_km_{tgt}.png", dpi=150)
            plt.close(fig)
    res = pd.DataFrame(rows)
    res.to_csv(f"{C.TABLES}/ici_gse135222_target_pfs.csv", index=False)
    return res


def main() -> None:
    np.random.seed(C.RANDOM_SEED)
    r1 = run_gse126044()
    r2 = run_gse135222()
    print("=== GSE126044 responder vs non-responder ===")
    print(r1[r1.feature.isin(TARGETS + ["sig_CD8_effector"])]
          [["feature", "median_resp", "median_nonresp", "auc_resp_high", "p", "q_BH"]]
          .to_string(index=False))
    print("\n=== GSE135222 PFS ===")
    print(r2[r2.feature.isin(TARGETS + ["sig_CD8_effector"])]
          [["feature", "spearman_rho_vs_pfstime", "spearman_p", "logrank_p_medsplit"]]
          .to_string(index=False))


if __name__ == "__main__":
    main()
