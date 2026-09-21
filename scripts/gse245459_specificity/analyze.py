#!/usr/bin/env python3
"""GSE245459 SKOV3 shTACSTD2: is the IFN/APM shift specific?

Prior slice (PR #97) reported that untreated shTACSTD2 lowers CLDN4 and
closes IFN/APM, with a junction dip it called short of a scaffold collapse.
This script does not redo that claim. It asks whether the IFN/APM movement
is a focused program or the same thing as a global / epithelial crash.

FPKM is compositional (columns are relative abundance). An absolute
transcriptional shutdown is not identifiable from FPKM sums. A global crash
in this matrix would be a broad relative shift: most expressed genes down,
housekeeping down, and the lost mass absorbed by a small set of survivors.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "gse245459_specificity"
RAW = OUT / "raw" / "GSE245459_fpkm.anno.txt.gz"
TABLES = OUT / "tables"
FIGS = OUT / "figures"

URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245459/"
    "suppl/GSE245459_fpkm.anno.txt.gz"
)

SAMPLE_COLS = [
    "shNC1",
    "shNC2",
    "shNC3",
    "shNCDDP1",
    "shNCDDP2",
    "shNCDDP3",
    "sh1",
    "sh2",
    "sh3",
    "shDDP1",
    "shDDP2",
    "shDDP3",
]

GROUPS = {
    "shNC": ["shNC1", "shNC2", "shNC3"],
    "shTACSTD2": ["sh1", "sh2", "sh3"],
    "shNC_DDP": ["shNCDDP1", "shNCDDP2", "shNCDDP3"],
    "shTACSTD2_DDP": ["shDDP1", "shDDP2", "shDDP3"],
}

CONTRASTS = [
    ("KD_noDDP", "shTACSTD2", "shNC"),
    ("KD_DDP", "shTACSTD2_DDP", "shNC_DDP"),
]

PSEUDO = 1.0
MIN_MEAN_FPKM = 1.0
MIN_GROUP_FPKM = 0.3
N_PERM = 5000
SEED = 245459

# Same panels as the prior C4 analog, so medians can be checked against it.
C4_IFN_MHCI = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]
APM = [
    "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5",
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G",
    "CALR", "CANX", "ERAP1", "ERAP2", "PDIA3",
]
IFN = [
    "ISG15", "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL", "RSAD2",
    "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFITM1", "IFITM2", "IFITM3",
    "IFI6", "IFI27", "IFI35", "IFI44", "IFI44L", "IFI16",
    "DDX58", "IFIH1", "DDX60", "XAF1", "HERC5", "USP18", "CMPK2", "BST2",
    "GBP1", "GBP2", "GBP4", "GBP5", "EIF2AK2", "ZBP1",
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "JAK1", "JAK2",
    "SOCS1", "SOCS3", "CXCL10", "CXCL11", "CXCL9", "CCL5",
    "IFNB1", "IFNG", "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "IFNL1",
]
JUNCTION = [
    "CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7",
    "CLDN8", "CLDN9", "CLDN10", "CLDN12", "CLDN18",
    "TJP1", "TJP2", "TJP3", "OCLN", "MARVELD2", "MARVELD3",
    "F11R", "JAM2", "JAM3", "CGN", "CDH1", "CTNNB1", "CTNNA1", "EPCAM", "CRB3",
]

# Specificity controls. The knockdown gene TACSTD2 is not used as evidence
# that an epithelial program collapsed.
EPITHELIAL = [
    "KRT7", "KRT8", "KRT18", "KRT19", "KRT80", "KRT5", "KRT6A", "KRT6B",
    "KRT14", "KRT17", "KRT23", "EPCAM", "CDH1", "CDH3", "CLDN1", "CLDN3",
    "CLDN4", "CLDN7", "OCLN", "TJP1", "F11R", "MUC1", "MUC16", "ELF3",
    "GRHL2", "ESRP1", "ESRP2", "SPINT1", "SPINT2", "CD24", "S100A14", "CRB3",
]
KERATIN = [
    "KRT7", "KRT8", "KRT18", "KRT19", "KRT80", "KRT5", "KRT6A", "KRT6B",
    "KRT14", "KRT17", "KRT23",
]
HOUSEKEEPING = [
    "ACTB", "GAPDH", "PPIA", "RPLP0", "RPL13A", "RPS18", "RPS27", "HPRT1",
    "TBP", "UBC", "YWHAZ", "SDHA", "PGK1", "NONO", "HMBS", "GUSB",
]
RIBOSOME = [
    "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8", "RPL10", "RPL11",
    "RPL13", "RPL13A", "RPL18", "RPL19", "RPL23", "RPL27", "RPL30",
    "RPS2", "RPS3", "RPS6", "RPS8", "RPS14", "RPS18", "RPS19", "RPS27",
]
CELL_CYCLE = [
    "MKI67", "PCNA", "TOP2A", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
    "CDK1", "CCNB1", "CCNB2", "CCNA2", "AURKA", "AURKB", "BIRC5", "E2F1",
    "FOXM1", "CDC20", "PLK1",
]
MESENCHYMAL = [
    "VIM", "CDH2", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2", "TWIST1",
    "SERPINE1", "TGFB1", "ACTA2", "SPARC", "COL1A1", "COL1A2",
]
APOPTOSIS = [
    "BAX", "BAK1", "BBC3", "PMAIP1", "BCL2", "BCL2L1", "CASP3", "CASP7",
    "CASP8", "CASP9", "GADD45A", "DDIT3", "ATF4", "FAS", "TNFRSF10B",
]

SETS = {
    "C4_IFN_MHCI": C4_IFN_MHCI,
    "C4_IFN_MHCI_noHLA": [g for g in C4_IFN_MHCI if g != "HLA-A"],
    "APM": APM,
    "APM_noHLA": [g for g in APM if g != "HLA-A"],
    "IFN": IFN,
    "JUNCTION": JUNCTION,
    "EPITHELIAL": EPITHELIAL,
    "EPITHELIAL_noTarget": [g for g in EPITHELIAL if g not in ("TACSTD2", "CLDN4")],
    "KERATIN": KERATIN,
    "HOUSEKEEPING": HOUSEKEEPING,
    "RIBOSOME": RIBOSOME,
    "CELL_CYCLE": CELL_CYCLE,
    "MESENCHYMAL": MESENCHYMAL,
    "APOPTOSIS": APOPTOSIS,
}


def ensure_raw() -> None:
    RAW.parent.mkdir(parents=True, exist_ok=True)
    if RAW.exists() and RAW.stat().st_size > 1_000_000:
        return
    import urllib.request

    print(f"GET {URL}")
    urllib.request.urlretrieve(URL, RAW)


def load_fpkm() -> pd.DataFrame:
    df = pd.read_csv(RAW, sep="\t", compression="gzip", low_memory=False)
    for c in SAMPLE_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["GeneName"] = df["GeneName"].astype(str).str.strip()
    df = df[df["GeneName"].notna() & (df["GeneName"] != "") & (df["GeneName"] != "nan")]
    if "Biotype" in df.columns:
        coding = df[df["Biotype"].astype(str) == "protein_coding"].copy()
        rest = df[df["Biotype"].astype(str) != "protein_coding"].copy()
        df = pd.concat([coding, rest], ignore_index=True)
    df["mean_fpkm"] = df[SAMPLE_COLS].mean(axis=1)
    df = df.sort_values("mean_fpkm", ascending=False)
    df = df.drop_duplicates("GeneName", keep="first")
    return df.set_index("GeneName")


def welch_contrast(logx: pd.DataFrame, treat: list[str], ctrl: list[str]) -> pd.DataFrame:
    a = logx[treat].to_numpy(dtype=float)
    b = logx[ctrl].to_numpy(dtype=float)
    n_a = np.isfinite(a).sum(axis=1)
    n_b = np.isfinite(b).sum(axis=1)
    mean_a = np.nanmean(a, axis=1)
    mean_b = np.nanmean(b, axis=1)
    var_a = np.nanvar(a, axis=1, ddof=1)
    var_b = np.nanvar(b, axis=1, ddof=1)
    log2fc = mean_a - mean_b
    se2 = var_a / np.maximum(n_a, 1) + var_b / np.maximum(n_b, 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        tvals = log2fc / np.sqrt(se2)
        df = se2**2 / (
            (var_a / np.maximum(n_a, 1)) ** 2 / np.maximum(n_a - 1, 1)
            + (var_b / np.maximum(n_b, 1)) ** 2 / np.maximum(n_b - 1, 1)
        )
        pvals = 2 * stats.t.sf(np.abs(tvals), df)
    bad = (n_a < 2) | (n_b < 2) | ~np.isfinite(se2) | (se2 <= 0)
    tvals = tvals.copy()
    pvals = pvals.copy()
    tvals[bad] = np.nan
    pvals[bad] = np.nan
    out = pd.DataFrame(
        {
            "gene": logx.index,
            "mean_treat_log2": mean_a,
            "mean_ctrl_log2": mean_b,
            "log2FC": log2fc,
            "t": tvals,
            "p": pvals,
        }
    )
    mask = out["p"].notna()
    q = np.full(len(out), np.nan)
    if mask.sum() > 0:
        q[mask.to_numpy()] = multipletests(out.loc[mask, "p"], method="fdr_bh")[1]
    out["q"] = q
    return out.set_index("gene")


def perm_median(obs_genes: list[str], pool: pd.Series, n_perm: int, rng: np.random.Generator) -> dict:
    """Two-sided permutation of the set median against same-sized draws."""
    present = [g for g in obs_genes if g in pool.index]
    if len(present) < 3:
        return {"n": len(present), "median": None, "perm_p": None, "null_median": None}
    obs = float(pool.loc[present].median())
    universe = pool.drop(index=present, errors="ignore").to_numpy()
    k = len(present)
    null = np.empty(n_perm)
    for i in range(n_perm):
        draw = rng.choice(universe, size=k, replace=False)
        null[i] = np.median(draw)
    center = float(np.median(null))
    p = (np.sum(np.abs(null - center) >= np.abs(obs - center)) + 1) / (n_perm + 1)
    tail = (np.sum(null <= obs) + 1) / (n_perm + 1)
    return {
        "n": k,
        "median": obs,
        "null_median": center,
        "perm_p": float(p),
        "lower_tail_p": float(tail),
        "null_p05": float(np.quantile(null, 0.05)),
        "null_p95": float(np.quantile(null, 0.95)),
    }


def matched_perm(
    obs_genes: list[str],
    de: pd.DataFrame,
    expressed: pd.Series,
    n_perm: int,
    rng: np.random.Generator,
) -> dict:
    """Expression-decile matched permutation of the set median."""
    bg = de.loc[expressed].dropna(subset=["log2FC"]).copy()
    bg["decile"] = pd.qcut(bg["mean_expr"], 10, labels=False, duplicates="drop")
    present = [g for g in obs_genes if g in bg.index]
    if len(present) < 3:
        return {"n": len(present), "matched_perm_p": None}
    set_df = bg.loc[present]
    obs = float(set_df["log2FC"].median())
    counts = set_df["decile"].value_counts().to_dict()
    pools = {}
    for d, k in counts.items():
        pool = bg.loc[(bg["decile"] == d) & (~bg.index.isin(present)), "log2FC"].to_numpy()
        if len(pool) < k:
            return {"n": len(present), "matched_perm_p": None, "median": obs}
        pools[d] = (pool, k)
    null = np.empty(n_perm)
    for i in range(n_perm):
        draws = [rng.choice(pool, size=k, replace=False) for pool, k in pools.values()]
        null[i] = np.median(np.concatenate(draws))
    center = float(np.median(null))
    p = (np.sum(np.abs(null - center) >= np.abs(obs - center)) + 1) / (n_perm + 1)
    tail = (np.sum(null <= obs) + 1) / (n_perm + 1)
    return {
        "n": len(present),
        "median": obs,
        "matched_null_median": center,
        "matched_perm_p": float(p),
        "matched_lower_tail_p": float(tail),
    }


def mw_vs_bg(set_fc: pd.Series, bg_fc: pd.Series) -> dict:
    if len(set_fc) < 3 or len(bg_fc) < 20:
        return {"mw_p": None, "delta_median": None}
    u, p = stats.mannwhitneyu(set_fc, bg_fc, alternative="two-sided")
    return {
        "mw_U": float(u),
        "mw_p": float(p),
        "delta_median": float(set_fc.median() - bg_fc.median()),
    }


def _sign_counts(fc: pd.Series) -> dict:
    fc = fc.dropna()
    n_up = int((fc > 0).sum())
    n_down = int((fc < 0).sum())
    binom_p = None
    if len(fc):
        binom_p = float(stats.binomtest(n_up, n=int(len(fc)), p=0.5, alternative="two-sided").pvalue)
    return {
        "n": int(len(fc)),
        "median": float(fc.median()) if len(fc) else None,
        "mean": float(fc.mean()) if len(fc) else None,
        "n_up": n_up,
        "n_down": n_down,
        "binom_p": binom_p,
    }


def program_row(name: str, genes: list[str], de: pd.DataFrame, eligible_idx: pd.Index, rng: np.random.Generator) -> dict:
    """Specificity on genes with mean FPKM >= 0.3 in either contrasted group.

    `median_log2FC_all` uses every symbol present in the matrix, including
    genes below that floor, so it can be checked against the prior slice.
    Tests use the eligible subset. The null is other eligible genes.
    """
    all_fc = de.reindex(genes)["log2FC"].dropna()
    all_stats = _sign_counts(all_fc)
    pool = de.loc[eligible_idx, "log2FC"].dropna()
    present = [g for g in genes if g in pool.index]
    set_fc = pool.loc[present]
    bg_fc = pool.drop(index=present, errors="ignore")
    stats_e = _sign_counts(set_fc)
    mw = mw_vs_bg(set_fc, bg_fc)
    perm = perm_median(present, pool, N_PERM, rng)
    matched = matched_perm(genes, de, eligible_idx, N_PERM, rng)
    return {
        "set": name,
        "n_in_set": len(genes),
        "n_all": all_stats["n"],
        "median_log2FC_all": all_stats["median"],
        "n_up_all": all_stats["n_up"],
        "n_down_all": all_stats["n_down"],
        "binom_p_all": all_stats["binom_p"],
        "n_eligible": stats_e["n"],
        "median_log2FC": stats_e["median"],
        "mean_log2FC": stats_e["mean"],
        "n_up": stats_e["n_up"],
        "n_down": stats_e["n_down"],
        "binom_p": stats_e["binom_p"],
        "bg_median_log2FC": float(bg_fc.median()) if len(bg_fc) else None,
        **mw,
        "perm_p": perm.get("perm_p"),
        "lower_tail_p": perm.get("lower_tail_p"),
        "null_p05": perm.get("null_p05"),
        "null_p95": perm.get("null_p95"),
        "matched_perm_p": matched.get("matched_perm_p"),
        "matched_lower_tail_p": matched.get("matched_lower_tail_p"),
        "matched_null_median": matched.get("matched_null_median"),
    }


def direct_mw(name_a: str, genes_a: list[str], name_b: str, genes_b: list[str], de: pd.DataFrame, eligible_idx: pd.Index) -> dict:
    pool = de.loc[eligible_idx, "log2FC"].dropna()
    a = pool.reindex([g for g in genes_a if g in pool.index]).dropna()
    b = pool.reindex([g for g in genes_b if g in pool.index]).dropna()
    # Drop shared genes from the reference so the test is not partly circular.
    b = b.drop(index=a.index.intersection(b.index), errors="ignore")
    if len(a) < 3 or len(b) < 3:
        return {
            "contrast_sets": f"{name_a}_vs_{name_b}",
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "median_a": float(a.median()) if len(a) else None,
            "median_b": float(b.median()) if len(b) else None,
            "mw_p": None,
        }
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "contrast_sets": f"{name_a}_vs_{name_b}",
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(a.median()),
        "median_b": float(b.median()),
        "delta_median_a_minus_b": float(a.median() - b.median()),
        "mw_U": float(u),
        "mw_p": float(p),
    }


def global_row(contrast: str, de: pd.DataFrame, eligible_idx: pd.Index) -> dict:
    fc = de.loc[eligible_idx, "log2FC"].dropna()
    q = de.loc[fc.index, "q"]
    return {
        "contrast": contrast,
        "n_expressed": int(len(fc)),
        "median_log2FC": float(fc.median()),
        "mean_log2FC": float(fc.mean()),
        "q25": float(fc.quantile(0.25)),
        "q75": float(fc.quantile(0.75)),
        "frac_lt_0": float((fc < 0).mean()),
        "frac_lt_m0_5": float((fc < -0.5).mean()),
        "frac_lt_m1": float((fc < -1).mean()),
        "frac_gt_0_5": float((fc > 0.5).mean()),
        "frac_gt_1": float((fc > 1).mean()),
        "n_q05_down": int(((q < 0.05) & (fc < 0)).sum()),
        "n_q05_up": int(((q < 0.05) & (fc > 0)).sum()),
        "spearman_expr_vs_lfc": float(stats.spearmanr(de.loc[fc.index, "mean_expr"], fc).statistic),
        "spearman_expr_vs_lfc_p": float(stats.spearmanr(de.loc[fc.index, "mean_expr"], fc).pvalue),
    }


def share_table(fpkm: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    present = [g for g in genes if g in fpkm.index]
    total = fpkm[SAMPLE_COLS].sum(axis=0)
    part = fpkm.loc[present, SAMPLE_COLS].sum(axis=0)
    return (part / total).rename("share")


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    ensure_raw()

    fpkm = load_fpkm()
    logx = np.log2(fpkm[SAMPLE_COLS] + PSEUDO)
    mean_fpkm = fpkm[SAMPLE_COLS].mean(axis=1)
    expressed = mean_fpkm >= MIN_MEAN_FPKM

    qc_rows = []
    for sample in SAMPLE_COLS:
        qc_rows.append(
            {
                "sample": sample,
                "fpkm_sum": float(fpkm[sample].sum()),
                "n_fpkm_ge_1": int((fpkm[sample] >= 1).sum()),
                "n_fpkm_gt_0": int((fpkm[sample] > 0).sum()),
                "median_fpkm_detected": float(fpkm.loc[fpkm[sample] > 0, sample].median()),
            }
        )
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(TABLES / "sample_qc.tsv", sep="\t", index=False)

    de_map: dict[str, pd.DataFrame] = {}
    eligible_map: dict[str, pd.Index] = {}
    gene_rows = []
    global_rows = []
    program_rows = []
    direct_rows = []
    share_rows = []
    mass_rows = []
    program_mass_rows = []

    rng = np.random.default_rng(SEED)

    for contrast, treat_g, ctrl_g in CONTRASTS:
        treat = GROUPS[treat_g]
        ctrl = GROUPS[ctrl_g]
        de = welch_contrast(logx, treat, ctrl)
        de["mean_fpkm"] = mean_fpkm.reindex(de.index)
        de["mean_fpkm_treat"] = fpkm[treat].mean(axis=1)
        de["mean_fpkm_ctrl"] = fpkm[ctrl].mean(axis=1)
        de["mean_fpkm_contrast"] = fpkm[treat + ctrl].mean(axis=1)
        de["mean_expr"] = logx[treat + ctrl].mean(axis=1)
        de["delta_fpkm"] = fpkm[treat].mean(axis=1) - fpkm[ctrl].mean(axis=1)
        de_map[contrast] = de
        eligible = de.index[(de["mean_fpkm_treat"] >= MIN_GROUP_FPKM) | (de["mean_fpkm_ctrl"] >= MIN_GROUP_FPKM)]
        eligible_map[contrast] = eligible
        global_rows.append(global_row(contrast, de, eligible))

        for name, genes in SETS.items():
            program_rows.append({"contrast": contrast, **program_row(name, genes, de, eligible, rng)})
            for g in genes:
                if g not in de.index:
                    continue
                gene_rows.append(
                    {
                        "contrast": contrast,
                        "set": name,
                        "gene": g,
                        "global_mean_ge_1": bool(expressed.get(g, False)),
                        "eligible": bool(g in set(eligible)),
                        "mean_fpkm": float(mean_fpkm.get(g, np.nan)),
                        "mean_fpkm_contrast": float(de.at[g, "mean_fpkm_contrast"]),
                        "log2FC": float(de.at[g, "log2FC"]) if np.isfinite(de.at[g, "log2FC"]) else None,
                        "p": float(de.at[g, "p"]) if np.isfinite(de.at[g, "p"]) else None,
                        "q": float(de.at[g, "q"]) if np.isfinite(de.at[g, "q"]) else None,
                        "delta_fpkm": float(de.at[g, "delta_fpkm"]),
                    }
                )

        for query in ("IFN", "APM", "APM_noHLA", "C4_IFN_MHCI", "C4_IFN_MHCI_noHLA"):
            for ref in ("EPITHELIAL_noTarget", "KERATIN", "JUNCTION", "HOUSEKEEPING", "CELL_CYCLE", "RIBOSOME"):
                direct_rows.append(
                    {
                        "contrast": contrast,
                        **direct_mw(query, SETS[query], ref, SETS[ref], de, eligible),
                    }
                )

        # Compositional share of each program.
        for name, genes in SETS.items():
            if name.endswith("noHLA") or name.endswith("noTarget"):
                continue
            shares = share_table(fpkm, genes)
            a = shares[treat].to_numpy()
            b = shares[ctrl].to_numpy()
            t, p = stats.ttest_ind(a, b, equal_var=False)
            share_rows.append(
                {
                    "contrast": contrast,
                    "set": name,
                    "mean_share_treat": float(np.mean(a)),
                    "mean_share_ctrl": float(np.mean(b)),
                    "delta_share": float(np.mean(a) - np.mean(b)),
                    "welch_p": float(p),
                }
            )

        neg_total = float((-de.loc[eligible, "delta_fpkm"].clip(upper=0)).sum())
        for name, genes in SETS.items():
            present = [g for g in genes if g in eligible]
            delta = de.loc[present, "delta_fpkm"]
            neg = float((-delta.clip(upper=0)).sum())
            pos = float(delta.clip(lower=0).sum())
            program_mass_rows.append(
                {
                    "contrast": contrast,
                    "set": name,
                    "neg_delta_fpkm": neg,
                    "pos_delta_fpkm": pos,
                    "fraction_of_eligible_negative_mass": (neg / neg_total) if neg_total else None,
                }
            )

        # Where relative mass moves, among eligible genes.
        sub = de.loc[eligible].dropna(subset=["delta_fpkm"]).sort_values("delta_fpkm", ascending=False)
        pos = sub.loc[sub["delta_fpkm"] > 0, "delta_fpkm"]
        neg = sub.loc[sub["delta_fpkm"] < 0, "delta_fpkm"]
        pos_mass = float(pos.sum())
        neg_mass = float((-neg).sum())
        top = pd.concat([sub.head(25), sub.tail(25)])
        for gene, r in top.iterrows():
            mass_rows.append(
                {
                    "contrast": contrast,
                    "gene": gene,
                    "delta_fpkm": float(r["delta_fpkm"]),
                    "log2FC": float(r["log2FC"]) if np.isfinite(r["log2FC"]) else None,
                    "q": float(r["q"]) if np.isfinite(r["q"]) else None,
                    "mean_fpkm": float(r["mean_fpkm"]),
                    "side": "gainer" if r["delta_fpkm"] > 0 else "loser",
                    "pos_mass_expressed": pos_mass,
                    "neg_mass_expressed": neg_mass,
                    "top25_pos_fraction": float(sub.head(25)["delta_fpkm"].clip(lower=0).sum() / pos_mass) if pos_mass else None,
                    "top25_neg_fraction": float((-sub.tail(25)["delta_fpkm"].clip(upper=0)).sum() / neg_mass) if neg_mass else None,
                }
            )

    programs = pd.DataFrame(program_rows)
    programs.to_csv(TABLES / "program_specificity.tsv", sep="\t", index=False)
    pd.DataFrame(global_rows).to_csv(TABLES / "global_shift.tsv", sep="\t", index=False)
    pd.DataFrame(direct_rows).to_csv(TABLES / "program_vs_reference.tsv", sep="\t", index=False)
    pd.DataFrame(gene_rows).to_csv(TABLES / "program_genes.tsv", sep="\t", index=False)
    pd.DataFrame(share_rows).to_csv(TABLES / "compositional_shares.tsv", sep="\t", index=False)
    pd.DataFrame(mass_rows).to_csv(TABLES / "fpkm_mass_movers.tsv", sep="\t", index=False)
    pd.DataFrame(program_mass_rows).to_csv(TABLES / "program_fpkm_mass.tsv", sep="\t", index=False)

    # Prior-slice checkpoints (KD_noDDP medians published in PR #97).
    checkpoints = {
        "C4_IFN_MHCI": -0.5714133202056344,
        "APM": -0.4016755331913764,
        "IFN": -0.4445895968539956,
        "JUNCTION": -0.22897141224898787,
    }
    check_rows = []
    subp = programs[programs["contrast"] == "KD_noDDP"].set_index("set")
    for name, expected in checkpoints.items():
        got = float(subp.at[name, "median_log2FC_all"])
        check_rows.append({"set": name, "prior_median": expected, "this_median": got, "abs_diff": abs(got - expected)})
    pd.DataFrame(check_rows).to_csv(TABLES / "prior_median_check.tsv", sep="\t", index=False)

    summary = {
        "dataset": "GSE245459",
        "contrast_primary": "shTACSTD2 vs shNC, no cisplatin",
        "contrast_secondary": "shTACSTD2+DDP vs shNC+DDP",
        "n_per_group": 3,
        "expression_floor": "specificity tests use genes with mean FPKM >= 0.3 in the treatment group or the control group",
        "log2FC": "mean log2(FPKM+1) treat minus control",
        "n_perm": N_PERM,
        "seed": SEED,
        "fpkm_note": "FPKM sums are relative. A global crash here means a broad shift of expressed genes, not a drop in column sums.",
        "global": global_rows,
        "sample_qc": qc_rows,
        "prior_median_check": check_rows,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2))

    plot(programs, de_map, eligible_map)
    cols = ["set", "n_all", "median_log2FC_all", "n_eligible", "median_log2FC", "n_up", "n_down", "delta_median", "mw_p", "perm_p", "matched_perm_p"]
    print(pd.DataFrame(global_rows).to_string(index=False))
    print(programs[programs["contrast"] == "KD_noDDP"][cols].to_string(index=False))
    print("--- KD_DDP ---")
    print(programs[programs["contrast"] == "KD_DDP"][cols].to_string(index=False))
    print("--- direct ---")
    print(pd.DataFrame(direct_rows).to_string(index=False))
    print("--- shares ---")
    print(pd.DataFrame(share_rows).to_string(index=False))
    print("--- prior check ---")
    print(pd.DataFrame(check_rows).to_string(index=False))
    print("--- qc ---")
    print(qc.to_string(index=False))


def plot(programs: pd.DataFrame, de_map: dict[str, pd.DataFrame], eligible_map: dict[str, pd.Index]) -> None:
    order = [
        "HOUSEKEEPING",
        "RIBOSOME",
        "KERATIN",
        "EPITHELIAL_noTarget",
        "JUNCTION",
        "CELL_CYCLE",
        "MESENCHYMAL",
        "APOPTOSIS",
        "IFN",
        "APM",
        "C4_IFN_MHCI",
    ]
    labels = {
        "HOUSEKEEPING": "Housekeeping",
        "RIBOSOME": "Ribosome",
        "KERATIN": "Keratin",
        "EPITHELIAL_noTarget": "Epithelial\n(no TACSTD2/CLDN4)",
        "JUNCTION": "Junction",
        "CELL_CYCLE": "Cell cycle",
        "MESENCHYMAL": "Mesenchymal",
        "APOPTOSIS": "Apoptosis",
        "IFN": "IFN/ISG",
        "APM": "APM",
        "C4_IFN_MHCI": "C4 panel",
    }
    colors = {
        "HOUSEKEEPING": "#4C78A8",
        "RIBOSOME": "#4C78A8",
        "KERATIN": "#F58518",
        "EPITHELIAL_noTarget": "#F58518",
        "JUNCTION": "#E45756",
        "CELL_CYCLE": "#54A24B",
        "MESENCHYMAL": "#B279A2",
        "APOPTOSIS": "#9D755D",
        "IFN": "#72B7B2",
        "APM": "#72B7B2",
        "C4_IFN_MHCI": "#EECA3B",
    }
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), constrained_layout=True)
    contrasts = [("KD_noDDP", "shTACSTD2 vs shNC (no cisplatin)"), ("KD_DDP", "shTACSTD2 vs shNC on cisplatin")]
    for row, (contrast, title) in enumerate(contrasts):
        de = de_map[contrast]
        pool = de.loc[eligible_map[contrast], "log2FC"].dropna()
        ax = axes[row, 0]
        data = []
        used = []
        for name in order:
            genes = [g for g in SETS[name] if g in pool.index]
            if len(genes) < 3:
                continue
            data.append(pool.loc[genes].to_numpy())
            used.append(name)
        ax.axhline(pool.median(), color="#333333", lw=1.0, label=f"genome median {pool.median():+.2f}")
        ax.axhspan(pool.quantile(0.25), pool.quantile(0.75), color="#333333", alpha=0.08, lw=0)
        bp = ax.boxplot(
            data,
            tick_labels=[labels[n] for n in used],
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "lw": 1.2},
            whiskerprops={"color": "#444444"},
            capprops={"color": "#444444"},
        )
        for patch, name in zip(bp["boxes"], used):
            patch.set_facecolor(colors[name])
            patch.set_alpha(0.85)
        ax.axhline(0, color="#888888", lw=0.6, ls="--")
        ax.set_ylabel("log2FC  log2(FPKM+1)")
        ax.set_title(f"{title}\ngenome median {pool.median():+.2f} (n={len(pool):,})")
        ax.tick_params(axis="x", labelrotation=0, labelsize=7.5)
        ax.legend(frameon=False, fontsize=8, loc="lower left")

        axb = axes[row, 1]
        sub = programs[programs["contrast"] == contrast].set_index("set")
        ys = []
        yerr_lo = []
        yerr_hi = []
        xs = []
        cols = []
        for i, name in enumerate(used):
            med = sub.at[name, "median_log2FC"]
            ys.append(med)
            xs.append(i)
            cols.append(colors[name])
            yerr_lo.append(0)
            yerr_hi.append(0)
        axb.bar(xs, ys, color=cols, edgecolor="black", linewidth=0.4)
        axb.axhline(pool.median(), color="#333333", lw=1.0)
        axb.axhline(0, color="#888888", lw=0.6, ls="--")
        axb.set_xticks(xs)
        axb.set_xticklabels([labels[n] for n in used], fontsize=7.5)
        axb.set_ylabel("median log2FC")
        axb.set_title(f"Median log2FC (* expression-matched p<0.05)\ngenome median {pool.median():+.2f}")
        for i, name in enumerate(used):
            p = sub.at[name, "matched_perm_p"]
            if pd.notna(p) and p < 0.05:
                axb.text(i, ys[i] + (0.04 if ys[i] >= 0 else -0.08), "*", ha="center", va="bottom" if ys[i] >= 0 else "top", fontsize=12)

    fig.suptitle("GSE245459 SKOV3 shTACSTD2 — IFN/APM vs global and epithelial shift", fontsize=12)
    fig.savefig(FIGS / "ifn_vs_global_epithelial.png", dpi=160)
    fig.savefig(FIGS / "ifn_vs_global_epithelial.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
