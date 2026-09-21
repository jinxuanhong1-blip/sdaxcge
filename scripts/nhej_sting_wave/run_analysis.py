#!/usr/bin/env python3
"""CLDN4-low tumors vs mutation burden, c-NHEJ, STING/IFN/APM, and CD8.

Two cohorts, two layers:

1. MCT reproduction (Yamamoto et al., Mol Cancer Ther 2022, Fig. 4F).
   cBioPortal Firehose Legacy RNA-seq V2 RSEM and MUTATION_COUNT.
   OV also has CPTAC CLDN4 protein on the same study. The published
   comparison is mutation-count bins 2–69 vs 70–1899. LUAD is the same
   contrast on luad_tcga (the lung extension). This layer does not use GSVA.

2. Logic-wave scores on Xena STAR log2(TPM+1), primary tumor, one sample
   per patient. GSVA is the Bioconductor GSVA Gaussian-kernel path
   (Hänzelmann et al. 2013; bandwidth = sd/4, logit of the kernel CDF,
   column ranks with ties.method='last', tau=1, maxDiff=True,
   absRanking=False). LUAD correlations are also reported after a
   keratin partial Spearman (KRT18+KRT19, the locked adjustment; the
   KRT8+KRT18+KRT19 triad is a sensitivity).

Raw matrices are not written. Run download_data.py first.
"""

import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA_DIR = os.environ.get("NHEJ_STING_DATA", "/tmp/nhej_sting_data")
OUT_DIR = os.environ.get(
    "NHEJ_STING_OUT",
    os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
        "results",
        "nhej_sting_wave",
    ),
)

# Pre-specified classical NHEJ core (Ku/DNA-PKcs/Lig4/XLF/PAXX/Artemis).
# XRCC1 is not in this set: it is BER/SSBR, reported only as a single-gene
# MCT effector. TP53BP1 is the MCT foci effector and is scored separately.
CNHEJ_CORE = [
    "XRCC6",
    "XRCC5",
    "PRKDC",
    "LIG4",
    "XRCC4",
    "NHEJ1",
    "PAXX",
    "DCLRE1C",
]
# Proximal cGAS–STING machinery, deliberately without the IFN/ISG output
# so this score is not a copy of the Hallmark IFN sets.
STING_PROXIMAL = ["MB21D1", "STING1", "TBK1", "IKBKE", "IRF3", "IRF7"]
# MHC-I antigen-presentation machinery. Overlap with Hallmark IFN is expected
# (B2M, TAP1, PSMB8/9) and is recorded, not hidden.
APM_GENES = [
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "B2M",
    "TAP1",
    "TAP2",
    "TAPBP",
    "PSMB8",
    "PSMB9",
    "PSMB10",
    "PSME1",
    "PSME2",
    "ERAP1",
    "ERAP2",
    "NLRC5",
    "CALR",
    "CANX",
    "PDIA3",
]
# HGNC updates that post-date GENCODE v36, plus a few older aliases.
ALIASES = {
    "WARS1": "WARS",
    "MARCHF1": "MARCH1",
    "TENT5A": "FAM46A",
    "STING1": "TMEM173",
    "PAXX": "C9orf142",
    "MB21D1": "CGAS",
    "METTL7B": "METTL7B",
}

# Primary Spearman family (BH within cohort, and separately within the
# LUAD keratin-adjusted family). Sensitivities are reported but not in this FDR.
PRIMARY_ENDPOINTS = [
    "log10_mutation_count",
    "gsva_cNHEJ_core",
    "gsva_STING_proximal",
    "gsva_IFN_alpha",
    "gsva_IFN_gamma",
    "gsva_APM",
    "CD8A",
    "MCP_CD8_T",
]

# Logic-wave sign of Spearman(CLDN4, endpoint). Negative means CLDN4-low
# tumors have the higher endpoint.
EXPECTED_SIGN = {
    "log10_mutation_count": -1,
    "gsva_cNHEJ_core": +1,
    "gsva_STING_proximal": -1,
    "gsva_IFN_alpha": -1,
    "gsva_IFN_gamma": -1,
    "gsva_APM": -1,
    "CD8A": -1,
    "MCP_CD8_T": -1,
    "gsva_KEGG_NHEJ": +1,
    "gsva_cNHEJ_plus_TP53BP1": +1,
    "TP53BP1": +1,
    "XRCC1": +1,
}

ENDPOINT_LABEL = {
    "log10_mutation_count": "log10 mutation count",
    "gsva_cNHEJ_core": "c-NHEJ core GSVA",
    "gsva_STING_proximal": "STING proximal GSVA",
    "gsva_IFN_alpha": "IFN-α GSVA",
    "gsva_IFN_gamma": "IFN-γ GSVA",
    "gsva_APM": "APM GSVA",
    "CD8A": "CD8A",
    "MCP_CD8_T": "MCP-counter CD8 T",
    "gsva_KEGG_NHEJ": "KEGG NHEJ GSVA",
    "gsva_cNHEJ_plus_TP53BP1": "c-NHEJ+TP53BP1 GSVA",
    "TP53BP1": "TP53BP1",
    "XRCC1": "XRCC1",
}


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    order = np.argsort(p[ok])
    ranked = p[ok][order]
    n = len(ranked)
    adj = ranked * n / (np.arange(n) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out = np.empty(n)
    out[order] = adj
    q[ok] = out
    return q


def spearman(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    n = int(mask.sum())
    if n < 5:
        return n, np.nan, np.nan
    r, p = stats.spearmanr(x[mask], y[mask])
    return n, float(r), float(p)


def fisher_ci(r, n, k_cov=0):
    """Fisher z 95% CI. k_cov extra covariates for a partial correlation."""
    if not np.isfinite(r) or n - k_cov - 3 <= 0 or abs(r) >= 1:
        return np.nan, np.nan
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n - k_cov - 3)
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def partial_spearman(x, y, cov):
    """Spearman partial correlation: Pearson of rank-residuals.

    cov is an (n, k) array. Rows with any non-finite value are dropped.
    The p-value uses df = n - 2 - k.
    """
    cov = np.asarray(cov, dtype=float)
    if cov.ndim == 1:
        cov = cov[:, None]
    mask = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(cov), axis=1)
    n = int(mask.sum())
    k = cov.shape[1]
    if n < k + 5:
        return n, np.nan, np.nan
    rx = stats.rankdata(x[mask])
    ry = stats.rankdata(y[mask])
    rc = np.column_stack([stats.rankdata(cov[mask, j]) for j in range(k)])

    def resid(v):
        A = np.column_stack([np.ones(n), rc])
        beta, *_ = np.linalg.lstsq(A, v, rcond=None)
        return v - A @ beta

    xr, yr = resid(rx), resid(ry)
    if np.std(xr) == 0 or np.std(yr) == 0:
        return n, np.nan, np.nan
    r = float(np.corrcoef(xr, yr)[0, 1])
    df = n - 2 - k
    if df <= 0 or abs(r) >= 1:
        return n, r, np.nan
    tstat = r * np.sqrt(df / (1.0 - r * r))
    p = float(2 * stats.t.sf(abs(tstat), df))
    return n, r, p


def col_ranks_last(Z):
    """Column ranks, 1 = lowest. Ties: earliest row index gets the highest rank.

    This is matrixStats colRanks(ties.method='last'), which GSVA uses so that
    ranks match a stable order() tie break.
    """
    p, n = Z.shape
    idx = np.arange(p)
    ranks = np.empty((p, n), dtype=np.int32)
    for j in range(n):
        order = np.lexsort((-idx, Z[:, j]))
        ranks[order, j] = np.arange(1, p + 1)
    return ranks


def gaussian_logit_ecdf(X):
    """GSVA row norm: logit of a Gaussian-kernel CDF, bandwidth sd/4.

    Matches kernel_estimation.c row_d (SIGMA_FACTOR=4) up to the precomputed
    normal-CDF grid, which is replaced here by erf.
    """
    n_genes, n = X.shape
    sd = X.std(axis=1, ddof=1)
    h = np.maximum(sd / 4.0, 1e-3)
    Z = np.empty((n_genes, n), dtype=np.float64)
    from scipy.special import erf

    inv_sqrt2 = 1.0 / np.sqrt(2.0)
    batch = 128
    n_batch = int(np.ceil(n_genes / batch))
    for b, s in enumerate(range(0, n_genes, batch)):
        xb = X[s : s + batch]
        hb = h[s : s + batch]
        diff = (xb[:, :, None] - xb[:, None, :]) / hb[:, None, None]
        left = (0.5 * (1.0 + erf(diff * inv_sqrt2))).mean(axis=2)
        left = np.clip(left, 1e-12, 1.0 - 1e-12)
        Z[s : s + batch] = np.log(left / (1.0 - left))
        if b % 25 == 0:
            print(f"    kernel batch {b + 1}/{n_batch}", flush=True)
    return Z


def gsva_scores(ranks, gene_idx):
    """Modified Kuiper ES: max(walk) + min(walk). tau=1, absRanking=False."""
    p, n = ranks.shape
    k = len(gene_idx)
    if k < 3 or k >= p:
        return np.full(n, np.nan)
    scores = np.empty(n, dtype=np.float64)
    gene_idx = np.asarray(gene_idx, dtype=np.int32)
    for j in range(n):
        r = ranks[:, j].astype(np.float64)
        dec = (p - r + 1).astype(np.int32)
        sym = np.abs(p / 2.0 - r)
        step_in = np.zeros(p, dtype=np.float64)
        step_out = np.ones(p, dtype=np.float64)
        pos = dec[gene_idx] - 1
        step_in[pos] = sym[gene_idx]
        step_out[pos] = 0.0
        cin = np.cumsum(step_in)
        cout = np.cumsum(step_out)
        if cin[-1] <= 0 or cout[-1] <= 0:
            scores[j] = np.nan
            continue
        wlk = cin / cin[-1] - cout / cout[-1]
        scores[j] = float(wlk.max() + wlk.min())
    return scores


def _self_test():
    z = np.array([[5.0], [5.0], [1.0]])
    r = col_ranks_last(z)
    if list(r[:, 0]) != [3, 2, 1]:
        raise RuntimeError(f"col_ranks_last failed: {r[:, 0]}")
    # Gene set = the three highest genes. The walk should finish positive.
    ranks = np.array([[4], [3], [2], [1]], dtype=np.int32)
    s = gsva_scores(ranks, np.array([0, 1, 2]))
    if not np.isfinite(s[0]) or s[0] <= 0:
        raise RuntimeError(f"gsva walk failed: {s}")
    # Partial Spearman of y=x with an independent covariate stays near 1.
    rng = np.random.default_rng(0)
    x = rng.normal(size=200)
    n_, r_, p_ = partial_spearman(x, x + rng.normal(scale=0.01, size=200), rng.normal(size=200))
    if n_ != 200 or r_ < 0.95:
        raise RuntimeError(f"partial spearman failed: {n_} {r_}")


def read_gmt_line(path, name):
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == name:
                return [g for g in parts[2:] if g]
    raise KeyError(name)


def resolve(symbols, available):
    """Map each requested symbol onto a matrix symbol. Return used, missing."""
    used, missing = [], []
    for g in symbols:
        if g in available:
            used.append(g)
        elif ALIASES.get(g) in available:
            used.append(ALIASES[g])
        else:
            missing.append(g)
    # unique, preserve order
    seen = set()
    uniq = []
    for g in used:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq, missing


def load_star(cohort, symbols_keep):
    """log2(TPM+1), genes x primary-tumor samples, one sample per patient."""
    probemap = pd.read_csv(
        os.path.join(DATA_DIR, "gencode.v36.probemap"),
        sep="\t",
        usecols=["id", "gene"],
    )
    id2sym = dict(zip(probemap["id"], probemap["gene"]))
    path = os.path.join(DATA_DIR, f"TCGA-{cohort}.star_tpm.tsv.gz")
    header = pd.read_csv(path, sep="\t", nrows=0)
    by_patient = {}
    n_tumor_rows = 0
    for s in header.columns[1:]:
        parts = s.split("-")
        if len(parts) < 4:
            continue
        code = parts[3][:2]
        if code != "01":
            continue
        n_tumor_rows += 1
        patient = "-".join(parts[:3])
        # Lexicographically first aliquot, so the choice does not depend on
        # column order in the Xena file.
        if patient not in by_patient or s < by_patient[patient]:
            by_patient[patient] = s
    tumor_cols = [by_patient[p] for p in sorted(by_patient)]
    expr = pd.read_csv(path, sep="\t", index_col=0, usecols=[header.columns[0], *tumor_cols])
    expr.index = expr.index.map(lambda i: id2sym.get(i, i))
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[expr.index.isin(symbols_keep)]
    expr = (
        expr.assign(_mean=expr.mean(axis=1))
        .sort_values("_mean", ascending=False)
        .drop(columns="_mean")
    )
    expr = expr[~expr.index.duplicated(keep="first")]
    meta = {
        "n_primary_tumor_columns": n_tumor_rows,
        "n_patients": int(expr.shape[1]),
    }
    return expr, meta


def load_mcp_cd8():
    probemap = pd.read_csv(
        os.path.join(DATA_DIR, "gencode.v36.probemap"),
        sep="\t",
        usecols=["id", "gene"],
    )
    ens2sym = {i.split(".")[0]: g for i, g in zip(probemap["id"], probemap["gene"])}
    sig = pd.read_csv(os.path.join(DATA_DIR, "mcpcounter_genes.txt"), sep="\t")
    sig.columns = [c.strip('"') for c in sig.columns]
    genes = []
    for _, row in sig.iterrows():
        if row["Cell population"] != "CD8 T cells":
            continue
        genes.append(ens2sym.get(str(row["ENSEMBL ID"]), str(row["HUGO symbols"]).strip('"')))
    return sorted(set(genes))


def load_cbioportal(name):
    df = pd.read_csv(os.path.join(DATA_DIR, name), sep="\t")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df[df["sample_id"].str.contains(r"-01$", regex=True)].copy()
    df = df.drop_duplicates("patient_id", keep="first")
    return df


def load_purity():
    path = os.path.join(DATA_DIR, "tcga_absolute_purity.txt")
    df = pd.read_csv(path, sep="\t", low_memory=False)
    # Prefer the 15-character sample barcode in `array` (TCGA-XX-XXXX-01).
    sample = df["array"].astype(str) if "array" in df.columns else df["sample"].astype(str)
    out = pd.DataFrame(
        {
            "sample": sample,
            "purity": pd.to_numeric(df["purity"], errors="coerce"),
            "call": df["call status"].astype(str) if "call status" in df.columns else "called",
        }
    )
    out = out.dropna(subset=["purity"])
    out = out[out["sample"].str.slice(13, 15) == "01"]
    out["patient_id"] = out["sample"].str.slice(0, 12)
    out = out.sort_values(["patient_id", "call"])
    # 'called' sorts before other statuses only if we rank it. Keep called rows first.
    out["_called"] = out["call"].eq("called")
    out = out.sort_values(["patient_id", "_called"], ascending=[True, False])
    out = out.drop_duplicates("patient_id", keep="first")
    return out.set_index("patient_id")["purity"]


def mutation_frame(mut_file, rsem_file):
    mut = load_cbioportal(mut_file)
    rsem = load_cbioportal(rsem_file)
    m = mut.merge(rsem, on="patient_id", suffixes=("_mut", "_rsem"))
    m = m.rename(columns={"value_mut": "mutation_count", "value_rsem": "cldn4_rsem"})
    m["cldn4_log2"] = np.log2(m["cldn4_rsem"].clip(lower=0) + 1.0)
    m["log10_mutation_count"] = np.nan
    pos = m["mutation_count"] >= 1
    m.loc[pos, "log10_mutation_count"] = np.log10(m.loc[pos, "mutation_count"])
    return m


def _mw(a, b):
    if len(a) < 5 or len(b) < 5:
        return np.nan, np.nan, np.nan
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    t, pt = stats.ttest_ind(a, b, equal_var=False)
    return float(u), float(p), float(pt)


def mct_tests(df, cohort, mut_source, expr_source):
    """Published 2–69 vs 70–1899 bins, plus a within-cohort median split.

    The 70 cut is the OV median of the PanCan count in this extract. It is
    not a median split in LUAD, where mutation counts sit higher.
    """
    sub = df[df["mutation_count"].between(2, 1899)].copy()
    sub["bin"] = np.where(sub["mutation_count"] <= 69, "low_2_69", "high_70_1899")
    low = sub.loc[sub["bin"] == "low_2_69", "cldn4_log2"].to_numpy()
    high = sub.loc[sub["bin"] == "high_70_1899", "cldn4_log2"].to_numpy()
    u, p_mw, p_t = _mw(low, high)
    n_sp, r_sp, p_sp = spearman(sub["cldn4_log2"].to_numpy(), sub["log10_mutation_count"].to_numpy())
    n_all, r_all, p_all = spearman(df["cldn4_log2"].to_numpy(), df["log10_mutation_count"].to_numpy())
    ge1 = df[df["mutation_count"] >= 1]
    med = float(ge1["mutation_count"].median()) if len(ge1) else np.nan
    le = ge1.loc[ge1["mutation_count"] <= med, "cldn4_log2"].to_numpy()
    gt = ge1.loc[ge1["mutation_count"] > med, "cldn4_log2"].to_numpy()
    _, p_med, _ = _mw(le, gt)
    row = {
        "cohort": cohort,
        "mutation_source": mut_source,
        "expression_source": expr_source,
        "n_joined": int(len(df)),
        "n_mutation_0": int((df["mutation_count"] == 0).sum()),
        "n_in_paper_range_2_1899": int(len(sub)),
        "n_low_2_69": int((sub["bin"] == "low_2_69").sum()) if len(sub) else 0,
        "n_high_70_1899": int((sub["bin"] == "high_70_1899").sum()) if len(sub) else 0,
        "median_cldn4_log2_low": float(np.median(low)) if len(low) else np.nan,
        "median_cldn4_log2_high": float(np.median(high)) if len(high) else np.nan,
        "median_cldn4_rsem_low": float(sub.loc[sub["bin"] == "low_2_69", "cldn4_rsem"].median())
        if (sub["bin"] == "low_2_69").any()
        else np.nan,
        "median_cldn4_rsem_high": float(sub.loc[sub["bin"] == "high_70_1899", "cldn4_rsem"].median())
        if (sub["bin"] == "high_70_1899").any()
        else np.nan,
        "mannwhitney_p_paper_bins": p_mw,
        "welch_p_log2_paper_bins": p_t,
        "high_bin_has_lower_cldn4": bool(np.median(high) < np.median(low))
        if len(low) and len(high)
        else False,
        "spearman_n_count_ge1": n_all,
        "spearman_rho_count_ge1": r_all,
        "spearman_p_count_ge1": p_all,
        "spearman_rho_paper_range": r_sp,
        "spearman_p_paper_range": p_sp,
        "within_cohort_median_mutation_count": med,
        "n_at_or_below_cohort_median": int(len(le)),
        "n_above_cohort_median": int(len(gt)),
        "median_cldn4_log2_at_or_below_cohort_median": float(np.median(le)) if len(le) else np.nan,
        "median_cldn4_log2_above_cohort_median": float(np.median(gt)) if len(gt) else np.nan,
        "mannwhitney_p_cohort_median_split": p_med,
        "mutation_min": float(df["mutation_count"].min()),
        "mutation_max": float(df["mutation_count"].max()),
    }
    return row, sub


def protein_tests(mut_file, mut_source):
    mut = load_cbioportal(mut_file).rename(columns={"value": "mutation_count"})
    prot = load_cbioportal("ov_cldn4_protein.tsv").rename(columns={"value": "cldn4_protein"})
    m = mut.merge(prot[["patient_id", "cldn4_protein"]], on="patient_id", how="inner")
    m = m.dropna(subset=["cldn4_protein", "mutation_count"])
    m["log10_mutation_count"] = np.nan
    pos = m["mutation_count"] >= 1
    m.loc[pos, "log10_mutation_count"] = np.log10(m.loc[pos, "mutation_count"])
    n, r, p = spearman(m["cldn4_protein"].to_numpy(), m["log10_mutation_count"].to_numpy())
    sub = m[m["mutation_count"].between(2, 1899)].copy()
    sub["bin"] = np.where(sub["mutation_count"] <= 69, "low_2_69", "high_70_1899")
    low = sub.loc[sub["bin"] == "low_2_69", "cldn4_protein"].to_numpy()
    high = sub.loc[sub["bin"] == "high_70_1899", "cldn4_protein"].to_numpy()
    if len(low) >= 5 and len(high) >= 5:
        _, p_mw = stats.mannwhitneyu(low, high, alternative="two-sided")
    else:
        p_mw = np.nan
    return {
        "cohort": "OV",
        "layer": "CPTAC_protein",
        "mutation_source": mut_source,
        "n_joined": int(len(m)),
        "n_in_paper_range_2_1899": int(len(sub)),
        "n_low_2_69": int((sub["bin"] == "low_2_69").sum()) if len(sub) else 0,
        "n_high_70_1899": int((sub["bin"] == "high_70_1899").sum()) if len(sub) else 0,
        "median_protein_low": float(np.median(low)) if len(low) else np.nan,
        "median_protein_high": float(np.median(high)) if len(high) else np.nan,
        "mannwhitney_p_paper_bins": float(p_mw) if np.isfinite(p_mw) else np.nan,
        "spearman_n_count_ge1": n,
        "spearman_rho_count_ge1": r,
        "spearman_p_count_ge1": p,
        "mutation_min": float(m["mutation_count"].min()),
        "mutation_max": float(m["mutation_count"].max()),
        "mutation_median": float(m["mutation_count"].median()),
    }, m


def detect_filter(expr, frac=0.10):
    detected = (expr > 0).mean(axis=1)
    keep = expr.loc[(detected >= frac) & (expr.std(axis=1) > 0)]
    return keep


def run_gsva(expr):
    """expr: genes x samples, already filtered. Returns score DataFrame."""
    genes = expr.index.to_list()
    X = expr.to_numpy(dtype=np.float64)
    print(f"  GSVA Gaussian kernel on {X.shape[0]} genes x {X.shape[1]} samples")
    Z = gaussian_logit_ecdf(X)
    print("  column ranks")
    ranks = col_ranks_last(Z)
    available = set(genes)
    index = {g: i for i, g in enumerate(genes)}

    ifn_a = read_gmt_line(os.path.join(DATA_DIR, "msigdb_hallmark_2020.txt"), "Interferon Alpha Response")
    ifn_g = read_gmt_line(os.path.join(DATA_DIR, "msigdb_hallmark_2020.txt"), "Interferon Gamma Response")
    kegg = read_gmt_line(os.path.join(DATA_DIR, "kegg_2021_human.txt"), "Non-homologous end-joining")

    sets = {
        "cNHEJ_core": CNHEJ_CORE,
        "cNHEJ_plus_TP53BP1": CNHEJ_CORE + ["TP53BP1"],
        "KEGG_NHEJ": kegg,
        "STING_proximal": STING_PROXIMAL,
        "IFN_alpha": ifn_a,
        "IFN_gamma": ifn_g,
        "APM": APM_GENES,
    }
    coverage = []
    scores = {"patient_id": [s[:12] for s in expr.columns], "sample_id": list(expr.columns)}
    for name, symbols in sets.items():
        used, missing = resolve(symbols, available)
        idx = np.array([index[g] for g in used], dtype=np.int32)
        print(f"  {name}: {len(used)}/{len(symbols)} genes in matrix; missing {missing}")
        scores[f"gsva_{name}"] = gsva_scores(ranks, idx)
        coverage.append(
            {
                "geneset": name,
                "n_requested": len(symbols),
                "n_used": len(used),
                "used": ",".join(used),
                "missing": ",".join(missing),
            }
        )
    return pd.DataFrame(scores), pd.DataFrame(coverage), available


def quartile_contrast(cldn4, y):
    mask = np.isfinite(cldn4) & np.isfinite(y)
    n = int(mask.sum())
    if n < 20:
        return {"n": n}
    r = stats.rankdata(cldn4[mask])
    # Q1 = lowest CLDN4, Q4 = highest. Rank quartiles avoid qcut edge ties.
    q = np.ceil(r / n * 4).astype(int)
    q = np.clip(q, 1, 4)
    yv = y[mask]
    low, high = yv[q == 1], yv[q == 4]
    u, p = stats.mannwhitneyu(low, high, alternative="two-sided")
    return {
        "n": n,
        "n_Q1": int((q == 1).sum()),
        "n_Q4": int((q == 4).sum()),
        "median_Q1_CLDN4low": float(np.median(low)),
        "median_Q4_CLDN4high": float(np.median(high)),
        "delta_Q1_minus_Q4": float(np.median(low) - np.median(high)),
        "mannwhitney_p": float(p),
    }


def correlate_block(df, cldn4, endpoints, cohort, adjustments):
    """adjustments: list of (tag, cov_matrix or None)."""
    rows = []
    for ep in endpoints:
        y = df[ep].to_numpy(dtype=float)
        x = cldn4
        n, r, p = spearman(x, y)
        lo, hi = fisher_ci(r, n, 0)
        expected = EXPECTED_SIGN.get(ep, np.nan)
        row = {
            "cohort": cohort,
            "endpoint": ep,
            "family": "primary" if ep in PRIMARY_ENDPOINTS else "sensitivity",
            "n": n,
            "rho": r,
            "p": p,
            "ci95_low": lo,
            "ci95_high": hi,
            "expected_sign": expected,
            "sign_matches_logic": bool(np.sign(r) == expected) if np.isfinite(r) else False,
        }
        for tag, cov in adjustments:
            if cov is None:
                continue
            nn, rr, pp = partial_spearman(x, y, cov)
            lo2, hi2 = fisher_ci(rr, nn, cov.shape[1] if cov.ndim == 2 else 1)
            row[f"n_{tag}"] = nn
            row[f"rho_{tag}"] = rr
            row[f"p_{tag}"] = pp
            row[f"ci95_low_{tag}"] = lo2
            row[f"ci95_high_{tag}"] = hi2
        qrow = quartile_contrast(x, y)
        for k, v in qrow.items():
            row[f"q_{k}"] = v
        rows.append(row)
    out = pd.DataFrame(rows)
    for cohort_name, sub in out.groupby("cohort"):
        pass
    # FDR within this cohort call (one cohort per call).
    prim = out["family"] == "primary"
    out.loc[prim, "q"] = bh(out.loc[prim, "p"].to_numpy())
    out.loc[~prim, "q"] = np.nan
    if "p_krt18_19" in out.columns:
        out.loc[prim, "q_krt18_19"] = bh(out.loc[prim, "p_krt18_19"].to_numpy())
    if "p_krt_triad" in out.columns:
        out.loc[prim, "q_krt_triad"] = bh(out.loc[prim, "p_krt_triad"].to_numpy())
    if "p_purity" in out.columns:
        out.loc[prim, "q_purity"] = bh(out.loc[prim, "p_purity"].to_numpy())
    if "p_mut" in out.columns:
        # mutation adjustment is only defined for non-mutation endpoints
        m = prim & out["endpoint"].ne("log10_mutation_count")
        out.loc[m, "q_mut"] = bh(out.loc[m, "p_mut"].to_numpy())
    return out


def star_patient_table(cohort, expr, scores, mut, purity):
    """One row per STAR primary tumor, with mutation count and purity joined."""
    df = scores.copy()
    df["CLDN4"] = expr.loc["CLDN4"].to_numpy() if "CLDN4" in expr.index else np.nan
    for g in ["KRT8", "KRT18", "KRT19", "CD8A", "CD8B", "TP53BP1", "XRCC1"]:
        df[g] = expr.loc[g].to_numpy() if g in expr.index else np.nan
    mcp_genes = [g for g in load_mcp_cd8() if g in expr.index]
    if mcp_genes:
        df["MCP_CD8_T"] = expr.loc[mcp_genes].mean(axis=0).to_numpy()
    else:
        df["MCP_CD8_T"] = np.nan
    df["mcp_cd8_genes"] = ",".join(mcp_genes)
    mut_s = mut.set_index("patient_id")["mutation_count"]
    df["mutation_count"] = df["patient_id"].map(mut_s)
    df["log10_mutation_count"] = np.nan
    pos = df["mutation_count"] >= 1
    df.loc[pos, "log10_mutation_count"] = np.log10(df.loc[pos, "mutation_count"].to_numpy())
    df["purity"] = df["patient_id"].map(purity)
    df["cohort"] = cohort
    return df, mcp_genes


def plot_mct(ov_sub, lu_sub, ov_prot, path):
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6))

    def box(ax, sub, ylab, title):
        data = [
            sub.loc[sub["bin"] == "low_2_69", ylab].dropna(),
            sub.loc[sub["bin"] == "high_70_1899", ylab].dropna(),
        ]
        bp = ax.boxplot(data, tick_labels=["2–69", "70–1899"], patch_artist=True, widths=0.6)
        for patch, color in zip(bp["boxes"], ["#4C78A8", "#E45756"]):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        ax.set_xlabel("Mutation count (Firehose)")
        ax.set_title(title)
        ns = [len(d) for d in data]
        ax.set_xticks([1, 2])
        ax.set_xticklabels([f"2–69\nn={ns[0]}", f"70–1899\nn={ns[1]}"])

    box(
        axes[0],
        ov_sub,
        "cldn4_log2",
        "OV mRNA  (paper bins)",
    )
    axes[0].set_ylabel("CLDN4 log2(RSEM+1)")
    box(axes[1], lu_sub, "cldn4_log2", "LUAD mRNA  (same bins)")
    axes[1].set_ylabel("CLDN4 log2(RSEM+1)")

    if ov_prot is not None and len(ov_prot):
        axes[2].scatter(
            ov_prot["mutation_count"],
            ov_prot["cldn4_protein"],
            s=12,
            c="#4C78A8",
            alpha=0.75,
            linewidths=0,
        )
        axes[2].set_xscale("log")
        axes[2].set_xlabel("Mutation count")
        axes[2].set_ylabel("CLDN4 protein (CPTAC)")
        axes[2].set_title(f"OV protein  n={len(ov_prot)}")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_forest(corr, path):
    eps = PRIMARY_ENDPOINTS
    fig, ax = plt.subplots(figsize=(8.2, 5.4))
    y = np.arange(len(eps))
    series = [
        ("OV raw", "OV", "rho", "ci95_low", "ci95_high", "#4C78A8", -0.22),
        ("LUAD raw", "LUAD", "rho", "ci95_low", "ci95_high", "#F58518", 0.0),
        ("LUAD KRT18/19", "LUAD", "rho_krt18_19", "ci95_low_krt18_19", "ci95_high_krt18_19", "#54A24B", 0.22),
    ]
    for label, cohort, rho, lo, hi, color, dy in series:
        sub = corr[(corr["cohort"] == cohort) & (corr["endpoint"].isin(eps))].set_index("endpoint")
        xs, xerr_lo, xerr_hi, ys = [], [], [], []
        for i, ep in enumerate(eps):
            if ep not in sub.index or not np.isfinite(sub.loc[ep, rho]):
                continue
            r = sub.loc[ep, rho]
            xs.append(r)
            xerr_lo.append(r - sub.loc[ep, lo])
            xerr_hi.append(sub.loc[ep, hi] - r)
            ys.append(i + dy)
        ax.errorbar(
            xs,
            ys,
            xerr=[xerr_lo, xerr_hi],
            fmt="o",
            color=color,
            label=label,
            ms=5,
            lw=1,
        )
    ax.axvline(0, color="0.4", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([ENDPOINT_LABEL[e] for e in eps])
    ax.set_xlabel("Spearman rho with CLDN4  (95% Fisher CI)")
    ax.set_xlim(-0.55, 0.55)
    ax.legend(frameon=False, loc="lower right")
    ax.set_title("CLDN4 vs NHEJ–STING endpoints")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_quartiles(tables, path):
    """Q1 vs Q4 medians for the logic-wave scores, OV and LUAD."""
    eps = [
        "log10_mutation_count",
        "gsva_cNHEJ_core",
        "gsva_STING_proximal",
        "gsva_IFN_gamma",
        "gsva_APM",
        "CD8A",
    ]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.2), sharey=False)
    for ax, cohort in zip(axes, ["OV", "LUAD"]):
        df = tables[cohort]
        x = np.arange(len(eps))
        q1, q4 = [], []
        for ep in eps:
            c = df["CLDN4"].to_numpy()
            y = df[ep].to_numpy(dtype=float)
            mask = np.isfinite(c) & np.isfinite(y)
            r = stats.rankdata(c[mask])
            n = int(mask.sum())
            q = np.clip(np.ceil(r / n * 4).astype(int), 1, 4)
            yv = y[mask]
            q1.append(np.median(yv[q == 1]))
            q4.append(np.median(yv[q == 4]))
        # z within the two medians so different units share an axis: show delta only
        # Better: grouped bars of the two medians after within-endpoint z is misleading.
        # Show Q1-Q4 delta of the rank of y, which is unit-free.
        deltas = []
        for ep in eps:
            c = df["CLDN4"].to_numpy()
            y = df[ep].to_numpy(dtype=float)
            mask = np.isfinite(c) & np.isfinite(y)
            r = stats.rankdata(c[mask])
            n = mask.sum()
            q = np.clip(np.ceil(r / n * 4).astype(int), 1, 4)
            yv = stats.rankdata(y[mask]) / n
            deltas.append(np.median(yv[q == 1]) - np.median(yv[q == 4]))
        colors = ["#E45756" if d > 0 else "#4C78A8" for d in deltas]
        ax.barh(np.arange(len(eps)), deltas, color=colors)
        ax.axvline(0, color="0.3", lw=0.8)
        ax.set_yticks(np.arange(len(eps)))
        ax.set_yticklabels([ENDPOINT_LABEL[e] for e in eps])
        ax.set_xlabel("Median rank(endpoint): CLDN4 Q1 − Q4")
        ax.set_title(f"{cohort}  (positive = higher in CLDN4-low)")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main():
    _self_test()
    os.makedirs(OUT_DIR, exist_ok=True)
    purity = load_purity()
    print(f"ABSOLUTE purity patients: {len(purity)}")

    specs = [
        ("OV", "PanCan_MUTATION_COUNT", "Firehose_RSEM", "ov_mutation_count_pancan.tsv", "ov_cldn4_rsem.tsv"),
        ("OV", "PanCan_MUTATION_COUNT", "PanCan_RSEM", "ov_mutation_count_pancan.tsv", "ov_cldn4_rsem_pancan.tsv"),
        ("OV", "Firehose_MUTATION_COUNT", "Firehose_RSEM", "ov_mutation_count_firehose.tsv", "ov_cldn4_rsem.tsv"),
        ("LUAD", "PanCan_MUTATION_COUNT", "Firehose_RSEM", "luad_mutation_count_pancan.tsv", "luad_cldn4_rsem.tsv"),
        ("LUAD", "PanCan_MUTATION_COUNT", "PanCan_RSEM", "luad_mutation_count_pancan.tsv", "luad_cldn4_rsem_pancan.tsv"),
        ("LUAD", "Firehose_MUTATION_COUNT", "Firehose_RSEM", "luad_mutation_count_firehose.tsv", "luad_cldn4_rsem.tsv"),
    ]
    mct_rows = []
    subs = {}
    for cohort, mut_source, expr_source, mut_file, rsem_file in specs:
        frame = mutation_frame(mut_file, rsem_file)
        row, sub = mct_tests(frame, cohort, mut_source, expr_source)
        mct_rows.append(row)
        subs[(cohort, mut_source, expr_source)] = sub
        print(
            f"MCT {cohort} {mut_source} x {expr_source}: "
            f"n={row['n_joined']} bins {row['n_low_2_69']}/{row['n_high_70_1899']} "
            f"MW={row['mannwhitney_p_paper_bins']:.3g} "
            f"rho={row['spearman_rho_count_ge1']:.3f} p={row['spearman_p_count_ge1']:.3g}"
        )
    prot_rows = []
    prot_frames = {}
    for source, mut_file in [
        ("Firehose_MUTATION_COUNT", "ov_mutation_count_firehose.tsv"),
        ("PanCan_MUTATION_COUNT", "ov_mutation_count_pancan.tsv"),
    ]:
        prot_row, prot_df = protein_tests(mut_file, source)
        prot_rows.append(prot_row)
        prot_frames[source] = prot_df
        print(
            f"MCT protein x {source}: n={prot_row['n_joined']} "
            f"rho={prot_row['spearman_rho_count_ge1']:.3f} p={prot_row['spearman_p_count_ge1']:.3g}"
        )
    pd.DataFrame(mct_rows).to_csv(os.path.join(OUT_DIR, "mct_mutation_burden.csv"), index=False)
    pd.DataFrame(prot_rows).to_csv(os.path.join(OUT_DIR, "mct_ov_protein.csv"), index=False)
    # Figure uses the count whose range matches Fig. 4F (PanCan, max ~1893)
    # joined to Firehose RSEM, and the protein join whose n is nearest 82.
    plot_mct(
        subs[("OV", "PanCan_MUTATION_COUNT", "Firehose_RSEM")],
        subs[("LUAD", "PanCan_MUTATION_COUNT", "Firehose_RSEM")],
        prot_frames["Firehose_MUTATION_COUNT"],
        os.path.join(OUT_DIR, "fig1_mct_mutation_bins.png"),
    )
    star_mut = {
        "OV": load_cbioportal("ov_mutation_count_pancan.tsv").rename(columns={"value": "mutation_count"}),
        "LUAD": load_cbioportal("luad_mutation_count_pancan.tsv").rename(columns={"value": "mutation_count"}),
    }

    # Genes required on top of the filtered transcriptome.
    extra = set(CNHEJ_CORE + STING_PROXIMAL + APM_GENES + [
        "CLDN4", "KRT8", "KRT18", "KRT19", "CD8A", "CD8B", "TP53BP1", "XRCC1", "PAXX", "TMEM173", "C9orf142",
    ])
    mcp_genes = load_mcp_cd8()
    extra.update(mcp_genes)

    corr_parts = []
    coverage_parts = []
    tables = {}
    summaries = {}
    for cohort in ["OV", "LUAD"]:
        print(f"\n=== STAR {cohort} ===")
        # Load a broad set: all genes that pass the detection filter. We need
        # the background, so load every symbolled gene by passing None-like
        # via a two-step read. load_star filters to symbols_keep; use all
        # probemap symbols.
        probemap = pd.read_csv(
            os.path.join(DATA_DIR, "gencode.v36.probemap"), sep="\t", usecols=["gene"]
        )
        expr_all, meta = load_star(cohort, set(probemap["gene"].dropna()))
        n_before = expr_all.shape[0]
        expr = detect_filter(expr_all, 0.10)
        # Keep curated single genes and the small sets even if a member is
        # sparse. Hallmark IFN genes stay under the detection filter.
        force = [
            g
            for g in (set(CNHEJ_CORE + STING_PROXIMAL + APM_GENES + list(extra)) )
            if g in expr_all.index and g not in expr.index
        ]
        if force:
            expr = pd.concat([expr, expr_all.loc[force]])
        print(f"  genes {n_before} -> {expr.shape[0]} after detection filter; patients {meta['n_patients']}")
        if "CLDN4" not in expr.index:
            raise RuntimeError(f"CLDN4 absent after filter in {cohort}")
        scores, coverage, _available = run_gsva(expr)
        coverage["cohort"] = cohort
        coverage_parts.append(coverage)
        table, mcp_used = star_patient_table(cohort, expr, scores, star_mut[cohort], purity)
        tables[cohort] = table
        cldn4 = table["CLDN4"].to_numpy(dtype=float)
        endpoints = [c for c in [
            "log10_mutation_count",
            "gsva_cNHEJ_core",
            "gsva_cNHEJ_plus_TP53BP1",
            "gsva_KEGG_NHEJ",
            "gsva_STING_proximal",
            "gsva_IFN_alpha",
            "gsva_IFN_gamma",
            "gsva_APM",
            "CD8A",
            "MCP_CD8_T",
            "TP53BP1",
            "XRCC1",
        ] if c in table.columns]
        adjustments = []
        if cohort == "LUAD":
            k18 = table[["KRT18", "KRT19"]].to_numpy(dtype=float)
            triad = table[["KRT8", "KRT18", "KRT19"]].to_numpy(dtype=float)
            adjustments.append(("krt18_19", k18))
            adjustments.append(("krt_triad", triad))
        pur = table["purity"].to_numpy(dtype=float)
        adjustments.append(("purity", pur[:, None]))
        # Mutation-count adjustment of the immune/NHEJ scores (not of itself).
        mutv = table["log10_mutation_count"].to_numpy(dtype=float)
        adjustments.append(("mut", mutv[:, None]))
        if cohort == "LUAD":
            both = np.column_stack([table[["KRT18", "KRT19"]].to_numpy(dtype=float), mutv])
            adjustments.append(("krt18_19_mut", both))
        corr = correlate_block(table, cldn4, endpoints, cohort, adjustments)
        corr_parts.append(corr)
        # Keratin collinearity, reported so the adjustment is interpretable.
        k_rows = {}
        for g in ["KRT8", "KRT18", "KRT19"]:
            n, r, p = spearman(cldn4, table[g].to_numpy(dtype=float))
            k_rows[g] = {"n": n, "rho": r, "p": p}
        summaries[cohort] = {
            "star": meta,
            "n_genes_gsva": int(expr.shape[0]),
            "mcp_cd8_genes": mcp_used,
            "cldn4_vs_keratin": k_rows,
            "n_with_mutation": int(np.isfinite(table["log10_mutation_count"]).sum()),
            "n_with_purity": int(np.isfinite(table["purity"]).sum()),
        }
        table.drop(columns=["mcp_cd8_genes"], errors="ignore").to_csv(
            os.path.join(OUT_DIR, f"sample_scores_{cohort}.csv"), index=False
        )

    corr_all = pd.concat(corr_parts, ignore_index=True)
    corr_all.to_csv(os.path.join(OUT_DIR, "endpoint_correlations.csv"), index=False)
    pd.concat(coverage_parts, ignore_index=True).to_csv(
        os.path.join(OUT_DIR, "geneset_coverage.csv"), index=False
    )
    plot_forest(corr_all, os.path.join(OUT_DIR, "fig2_spearman_forest.png"))
    plot_quartiles(tables, os.path.join(OUT_DIR, "fig3_cldn4low_minus_high.png"))

    # Quartile table (primary endpoints).
    q_rows = []
    for _, row in corr_all.iterrows():
        if row["family"] != "primary":
            continue
        q_rows.append({k: row[k] for k in row.index if k.startswith("q_") or k in ("cohort", "endpoint")})
    pd.DataFrame(q_rows).to_csv(os.path.join(OUT_DIR, "quartile_contrasts.csv"), index=False)

    with open(os.path.join(OUT_DIR, "cohort_summary.json"), "w") as fh:
        json.dump(
            {
                "mct_mrna": mct_rows,
                "mct_protein": prot_rows,
                "star": summaries,
                "primary_endpoints": PRIMARY_ENDPOINTS,
                "expected_sign": EXPECTED_SIGN,
                "gsva": {
                    "kcdf": "Gaussian",
                    "bandwidth": "sd/4",
                    "row_transform": "logit(kernel CDF)",
                    "ties": "matrixStats last",
                    "tau": 1,
                    "maxDiff": True,
                    "absRanking": False,
                    "detection_filter": "log2(TPM+1)>0 in >=10% of primary tumors",
                },
            },
            fh,
            indent=2,
            default=float,
        )
    print(f"\nWrote {OUT_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
