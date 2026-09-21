#!/usr/bin/env python3
"""Search GSE245459 for thesis-aligned IFN/APM increases.

The thesis wants IFN/APM up after loss of the barrier gene. This experiment is
TACSTD2/TROP2 shRNA in SKOV3, not a CLDN4 knockdown. The untreated contrast
in the prior slice was IFN/APM down. This script sweeps normalizations,
contrasts, gene sets, leave-one-out, and held-out CLDN4/TACSTD2 co-regulation
splits. A panel is thesis-aligned only when its median effect is up. The
search is post-hoc: p-values are ranked inside a declared test count and are
not a pre-registered confirmation.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as base  # noqa: E402

OUT = ROOT / "results" / "gse245459_specificity"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
GMT = Path(__file__).resolve().parent / "hallmark_subsets.gmt"

MIN_GROUP_FPKM = 0.3
MIN_GENES = 5

# Thesis estimands. Drug-only contrasts are computed and then set aside.
THESIS_ESTIMANDS = ("KD_noDDP", "KD_on_DDP", "interaction")


def load_hallmark() -> dict[str, list[str]]:
    sets = {}
    for line in GMT.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def ifn_family(hallmark: dict[str, list[str]]) -> dict[str, list[str]]:
    apm = list(base.APM)
    ifn = list(base.IFN)
    ifna = hallmark["Interferon Alpha Response"]
    ifng = hallmark["Interferon Gamma Response"]
    apm_set = set(apm) | {"PSME1", "PSME2", "TAPBPL"}
    inter = sorted(set(ifna) & set(ifng))
    return {
        "C4_panel": list(base.C4_IFN_MHCI),
        "C4_noHLA": [g for g in base.C4_IFN_MHCI if g != "HLA-A"],
        "APM": apm,
        "APM_noHLA": [g for g in apm if g != "HLA-A"],
        "APM_HLA": ["HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G"],
        "APM_peptide_loading": ["TAP1", "TAP2", "TAPBP", "TAPBPL", "ERAP1", "ERAP2", "CALR", "CANX", "PDIA3", "B2M"],
        "APM_immunoproteasome": ["PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2"],
        "IFN_custom": ifn,
        "H_IFNA": ifna,
        "H_IFNG": ifng,
        "H_IFNA_IFNG_core": inter,
        "H_IFNA_minus_APM": [g for g in ifna if g not in apm_set],
        "CXCL_ISG": ["CXCL9", "CXCL10", "CXCL11"],
    }


def context_sets(hallmark: dict[str, list[str]]) -> dict[str, list[str]]:
    return {
        "H_G2M": hallmark["G2-M Checkpoint"],
        "H_E2F": hallmark["E2F Targets"],
        "H_EMT": hallmark["Epithelial Mesenchymal Transition"],
        "H_apoptosis": hallmark["Apoptosis"],
        "H_TNFA": hallmark["TNF-alpha Signaling via NF-kB"],
        "H_IL6": hallmark["IL-6/JAK/STAT3 Signaling"],
        "H_inflammatory": hallmark["Inflammatory Response"],
        "H_allograft": hallmark["Allograft Rejection"],
        "EPITHELIAL_noTarget": list(base.SETS["EPITHELIAL_noTarget"]),
        "HOUSEKEEPING": list(base.HOUSEKEEPING),
    }


def transforms(fpkm: pd.DataFrame) -> dict[str, pd.DataFrame]:
    cols = base.SAMPLE_COLS
    out = {
        "log2_fpkm_p1": np.log2(fpkm[cols] + 1.0),
        "log2_fpkm_p0.1": np.log2(fpkm[cols] + 0.1),
        "log2_fpkm_p0.01": np.log2(fpkm[cols] + 0.01),
    }
    uq = fpkm[cols].where(fpkm[cols] > 0).quantile(0.75, axis=0)
    out["log2_uq_p1"] = np.log2(fpkm[cols].div(uq, axis=1) + 1.0)
    logged = np.log(fpkm[cols] + 1.0)
    out["clr_p1"] = logged.sub(logged.mean(axis=0), axis=1)
    hk = [g for g in base.HOUSEKEEPING if g in fpkm.index]
    x = np.log2(fpkm[cols] + 1.0)
    out["log2_p1_minus_HK"] = x.sub(x.loc[hk].mean(axis=0), axis=1)
    ranks = fpkm[cols].rank(axis=0, method="average")
    out["rank_within_sample"] = ranks / ranks.max(axis=0)
    # Positive-only log2. Genes with a zero in a contrasted sample are dropped later.
    with np.errstate(divide="ignore"):
        out["log2_fpkm_positive"] = np.log2(fpkm[cols].where(fpkm[cols] > 0))
    return out


def group_means(fpkm: pd.DataFrame) -> pd.DataFrame:
    means = {}
    for name, samples in base.GROUPS.items():
        means[name] = fpkm[samples].mean(axis=1)
    return pd.DataFrame(means)


def eligible_mask(gmeans: pd.DataFrame, groups: list[str]) -> pd.Series:
    sub = gmeans[groups]
    return sub.max(axis=1) >= MIN_GROUP_FPKM


def effect_vector(mat: pd.DataFrame, kind: str) -> pd.Series:
    g = base.GROUPS
    if kind == "KD_noDDP":
        return mat[g["shTACSTD2"]].mean(axis=1) - mat[g["shNC"]].mean(axis=1)
    if kind == "KD_on_DDP":
        return mat[g["shTACSTD2_DDP"]].mean(axis=1) - mat[g["shNC_DDP"]].mean(axis=1)
    if kind == "DDP_in_NC":
        return mat[g["shNC_DDP"]].mean(axis=1) - mat[g["shNC"]].mean(axis=1)
    if kind == "DDP_in_KD":
        return mat[g["shTACSTD2_DDP"]].mean(axis=1) - mat[g["shTACSTD2"]].mean(axis=1)
    if kind == "interaction":
        kd_on = mat[g["shTACSTD2_DDP"]].mean(axis=1) - mat[g["shNC_DDP"]].mean(axis=1)
        kd_off = mat[g["shTACSTD2"]].mean(axis=1) - mat[g["shNC"]].mean(axis=1)
        return kd_on - kd_off
    raise KeyError(kind)


def contrast_groups(kind: str) -> list[str]:
    if kind == "KD_noDDP":
        return ["shTACSTD2", "shNC"]
    if kind == "KD_on_DDP":
        return ["shTACSTD2_DDP", "shNC_DDP"]
    if kind == "DDP_in_NC":
        return ["shNC_DDP", "shNC"]
    if kind == "DDP_in_KD":
        return ["shTACSTD2_DDP", "shTACSTD2"]
    return list(base.GROUPS)


def sample_scores(mat: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in mat.index]
    return mat.loc[present].mean(axis=0)


def welch_greater(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """One-sided Welch p for mean(a) > mean(b). Returns delta, t, p."""
    if len(a) < 2 or len(b) < 2:
        return (float(np.mean(a) - np.mean(b)), np.nan, np.nan)
    t, p_two = stats.ttest_ind(a, b, equal_var=False)
    delta = float(np.mean(a) - np.mean(b))
    if not np.isfinite(t):
        return (delta, float(t) if np.isfinite(t) else np.nan, np.nan)
    p_one = p_two / 2 if t > 0 else 1 - p_two / 2
    return (delta, float(t), float(p_one))


def interaction_ols(scores: pd.Series) -> tuple[float, float]:
    rows = []
    for sample, value in scores.items():
        kd = int(sample in base.GROUPS["shTACSTD2"] or sample in base.GROUPS["shTACSTD2_DDP"])
        ddp = int(sample in base.GROUPS["shNC_DDP"] or sample in base.GROUPS["shTACSTD2_DDP"])
        rows.append({"score": float(value), "kd": kd, "ddp": ddp})
    fit = ols("score ~ kd * ddp", pd.DataFrame(rows)).fit()
    coef = float(fit.params["kd:ddp"])
    p_two = float(fit.pvalues["kd:ddp"])
    p_one = p_two / 2 if coef > 0 else 1 - p_two / 2
    return coef, float(p_one)


def gene_tests(effect: pd.Series) -> dict:
    fc = effect.dropna()
    n = int(len(fc))
    n_up = int((fc > 0).sum())
    n_down = int((fc < 0).sum())
    out = {
        "n": n,
        "median_effect": float(fc.median()) if n else np.nan,
        "mean_effect": float(fc.mean()) if n else np.nan,
        "n_up": n_up,
        "n_down": n_down,
    }
    if n < MIN_GENES:
        out.update({"sign_p_up": np.nan, "wilcoxon_p_up": np.nan})
        return out
    out["sign_p_up"] = float(stats.binomtest(n_up, n=n, p=0.5, alternative="greater").pvalue)
    # Wilcoxon needs non-zero differences. zero_method wilcox drops zeros.
    nonzero = fc[fc != 0]
    if len(nonzero) < MIN_GENES or (nonzero > 0).sum() == 0 or (nonzero < 0).sum() == 0:
        # All one sign: signed-rank is still defined if there is variation.
        pass
    try:
        stat, p = stats.wilcoxon(fc.to_numpy(), alternative="greater", zero_method="wilcox")
        out["wilcoxon_stat"] = float(stat)
        out["wilcoxon_p_up"] = float(p)
    except ValueError:
        out["wilcoxon_p_up"] = np.nan
    return out


def test_panel(mat: pd.DataFrame, genes: list[str], kind: str, eligible: pd.Series, fpkm: pd.DataFrame) -> dict:
    present = [g for g in genes if g in mat.index and bool(eligible.get(g, False))]
    if kind != "interaction" and "log2_fpkm_positive" in getattr(test_panel, "_norm", ""):
        pass
    effect = effect_vector(mat, kind)
    # Drop genes that are undefined in this transform (positive-only log).
    if kind == "log2_fpkm_positive":
        pass
    use = effect.reindex(present).dropna()
    # For positive-only log, require every sample in the contrast to be finite.
    stats_g = gene_tests(use)
    scores = sample_scores(mat.reindex(use.index), list(use.index)) if len(use) else pd.Series(dtype=float)
    if kind == "interaction":
        if len(use) >= MIN_GENES:
            coef, p_one = interaction_ols(scores)
        else:
            coef, p_one = (np.nan, np.nan)
        sample_delta, sample_p = coef, p_one
    elif kind == "KD_noDDP":
        sample_delta, _, sample_p = welch_greater(scores[base.GROUPS["shTACSTD2"]].to_numpy(), scores[base.GROUPS["shNC"]].to_numpy())
    elif kind == "KD_on_DDP":
        sample_delta, _, sample_p = welch_greater(scores[base.GROUPS["shTACSTD2_DDP"]].to_numpy(), scores[base.GROUPS["shNC_DDP"]].to_numpy())
    elif kind == "DDP_in_NC":
        sample_delta, _, sample_p = welch_greater(scores[base.GROUPS["shNC_DDP"]].to_numpy(), scores[base.GROUPS["shNC"]].to_numpy())
    elif kind == "DDP_in_KD":
        sample_delta, _, sample_p = welch_greater(scores[base.GROUPS["shTACSTD2_DDP"]].to_numpy(), scores[base.GROUPS["shTACSTD2"]].to_numpy())
    else:
        sample_delta, sample_p = (np.nan, np.nan)
    # Background: other eligible genes, one-sided MW that the set is greater.
    bg = effect.reindex(eligible[eligible].index).dropna()
    bg = bg.drop(index=use.index, errors="ignore")
    if len(use) >= MIN_GENES and len(bg) >= 20:
        _u, mw_p = stats.mannwhitneyu(use, bg, alternative="greater")
        mw_p = float(mw_p)
        bg_median = float(bg.median())
    else:
        mw_p, bg_median = (np.nan, float(bg.median()) if len(bg) else np.nan)
    return {
        **stats_g,
        "sample_delta": float(sample_delta) if sample_delta is not None and np.isfinite(sample_delta) else np.nan,
        "sample_p_up": float(sample_p) if sample_p is not None and np.isfinite(sample_p) else np.nan,
        "mw_p_gt_background": mw_p,
        "bg_median": bg_median,
        "thesis_aligned_median": bool(stats_g["n"] >= MIN_GENES and stats_g["median_effect"] > 0),
    }


def positive_only_eligible(fpkm: pd.DataFrame, kind: str, eligible: pd.Series) -> pd.Series:
    """Genes eligible on FPKM and strictly positive in every sample of the contrast."""
    if kind == "interaction":
        samples = base.SAMPLE_COLS
    elif kind == "KD_noDDP":
        samples = base.GROUPS["shTACSTD2"] + base.GROUPS["shNC"]
    elif kind == "KD_on_DDP":
        samples = base.GROUPS["shTACSTD2_DDP"] + base.GROUPS["shNC_DDP"]
    elif kind == "DDP_in_NC":
        samples = base.GROUPS["shNC_DDP"] + base.GROUPS["shNC"]
    else:
        samples = base.GROUPS["shTACSTD2_DDP"] + base.GROUPS["shTACSTD2"]
    positive = fpkm[samples].min(axis=1) > 0
    return eligible & positive


def sample_qc(fpkm: pd.DataFrame, logx: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    expressed = fpkm[base.SAMPLE_COLS].mean(axis=1) >= 1
    corr = logx.loc[expressed, base.SAMPLE_COLS].corr(method="pearson")
    corr.to_csv(TABLES / "sweep_sample_correlation.tsv", sep="\t")
    x = logx.loc[expressed, base.SAMPLE_COLS].replace([np.inf, -np.inf], np.nan).dropna()
    var = x.var(axis=1).sort_values(ascending=False)
    top = x.loc[var.head(2000).index]
    centered = top.sub(top.mean(axis=1), axis=0)
    # Samples × genes. A gene-wise tall SVD did not converge on this matrix.
    a = centered.to_numpy().T
    u, s, _vt = np.linalg.svd(a, full_matrices=False)
    scores = pd.DataFrame(u[:, :2] * s[:2], index=base.SAMPLE_COLS, columns=["PC1", "PC2"])
    total = float((s**2).sum())
    scores["PC1_var"] = float(s[0] ** 2 / total)
    scores["PC2_var"] = float(s[1] ** 2 / total)
    labels = []
    for sample in scores.index:
        if sample in base.GROUPS["shNC"]:
            labels.append("shNC")
        elif sample in base.GROUPS["shTACSTD2"]:
            labels.append("shTACSTD2")
        elif sample in base.GROUPS["shNC_DDP"]:
            labels.append("shNC_DDP")
        else:
            labels.append("shTACSTD2_DDP")
    scores["group"] = labels
    # Distance to other samples in the same group on PC1-PC2.
    xy = scores[["PC1", "PC2"]].to_numpy()
    dist = []
    for i, sample in enumerate(scores.index):
        mates = [j for j, lab in enumerate(labels) if lab == labels[i] and j != i]
        d = np.sqrt(((xy[mates] - xy[i]) ** 2).sum(axis=1)).mean()
        dist.append(float(d))
    scores["mean_pc_distance_to_groupmates"] = dist
    scores.to_csv(TABLES / "sweep_sample_pca.tsv", sep="\t")
    return corr, scores


def leave_one_out(mat: pd.DataFrame, family: dict[str, list[str]], eligible: pd.Series) -> pd.DataFrame:
    """Untreated contrast only. Includes the two drops most favorable to an up call."""
    samples = base.GROUPS["shNC"] + base.GROUPS["shTACSTD2"]
    rows = []
    base_nc = base.GROUPS["shNC"]
    base_sh = base.GROUPS["shTACSTD2"]
    # Identify extreme samples on the custom IFN score for a biased sensitivity.
    ifn_genes = [g for g in family["IFN_custom"] if g in mat.index and bool(eligible.get(g, False))]
    ifn_score = mat.loc[ifn_genes, samples].mean(axis=0)
    drop_high_nc = str(ifn_score[base_nc].idxmax())
    drop_low_sh = str(ifn_score[base_sh].idxmin())
    plans = [("none", [])]
    for s in samples:
        plans.append((f"drop_{s}", [s]))
    plans.append(("drop_highest_IFN_shNC_and_lowest_IFN_sh", [drop_high_nc, drop_low_sh]))
    for label, dropped in plans:
        nc = [s for s in base_nc if s not in dropped]
        sh = [s for s in base_sh if s not in dropped]
        if len(nc) < 2 or len(sh) < 2:
            continue
        effect = mat[sh].mean(axis=1) - mat[nc].mean(axis=1)
        for name, genes in family.items():
            if name not in ("C4_panel", "APM", "APM_noHLA", "IFN_custom", "H_IFNA", "H_IFNG"):
                continue
            present = [g for g in genes if g in effect.index and bool(eligible.get(g, False))]
            fc = effect.reindex(present).dropna()
            rows.append(
                {
                    "loo": label,
                    "dropped": ",".join(dropped) if dropped else "",
                    "biased_toward_up": label.startswith("drop_highest"),
                    "set": name,
                    "n": int(len(fc)),
                    "median_effect": float(fc.median()) if len(fc) else np.nan,
                    "n_up": int((fc > 0).sum()) if len(fc) else 0,
                    "n_down": int((fc < 0).sum()) if len(fc) else 0,
                    "sign_positive": bool(len(fc) and fc.median() > 0),
                }
            )
    return pd.DataFrame(rows)


def leave_one_out_on_ddp(mat: pd.DataFrame, family: dict[str, list[str]], eligible: pd.Series) -> pd.DataFrame:
    """KD vs scramble, both on cisplatin. Asks whether the IFN up is one library."""
    base_nc = base.GROUPS["shNC_DDP"]
    base_sh = base.GROUPS["shTACSTD2_DDP"]
    keep = ("C4_panel", "APM", "APM_noHLA", "IFN_custom", "H_IFNA", "H_IFNG")
    rows = []
    plans = [("none", [])] + [(f"drop_{s}", [s]) for s in base_nc + base_sh]
    for label, dropped in plans:
        nc = [s for s in base_nc if s not in dropped]
        sh = [s for s in base_sh if s not in dropped]
        effect = mat[sh].mean(axis=1) - mat[nc].mean(axis=1)
        for name in keep:
            present = [g for g in family[name] if g in effect.index and bool(eligible.get(g, False))]
            fc = effect.reindex(present).dropna()
            rows.append(
                {
                    "loo": label,
                    "dropped": ",".join(dropped),
                    "set": name,
                    "n": int(len(fc)),
                    "median_effect": float(fc.median()) if len(fc) else np.nan,
                    "n_up": int((fc > 0).sum()) if len(fc) else 0,
                    "n_down": int((fc < 0).sum()) if len(fc) else 0,
                    "sign_positive": bool(len(fc) and fc.median() > 0),
                }
            )
    return pd.DataFrame(rows)


def heldout_coreg(logx: pd.DataFrame, fpkm: pd.DataFrame, family: dict[str, list[str]], gmeans: pd.DataFrame) -> pd.DataFrame:
    """Define CLDN4/TACSTD2 correlation on one arm and test the KD effect on the other arm."""
    universe = sorted(set().union(*[set(v) for k, v in family.items() if k.startswith(("C4", "APM", "IFN", "H_IFN", "CXCL"))]))
    universe = [g for g in universe if g in logx.index]
    arms = {
        "discover_cisplatin_test_untreated": {
            "discover": base.GROUPS["shNC_DDP"] + base.GROUPS["shTACSTD2_DDP"],
            "test_kind": "KD_noDDP",
        },
        "discover_untreated_test_cisplatin": {
            "discover": base.GROUPS["shNC"] + base.GROUPS["shTACSTD2"],
            "test_kind": "KD_on_DDP",
        },
    }
    rows = []
    for split, spec in arms.items():
        disc = logx.loc[universe, spec["discover"]]
        anchor = {}
        for gene in ("CLDN4", "TACSTD2"):
            if gene not in logx.index:
                continue
            y = logx.loc[gene, spec["discover"]]
            rhos = []
            for g in universe:
                rho, _p = stats.spearmanr(disc.loc[g], y)
                rhos.append(rho if np.isfinite(rho) else np.nan)
            anchor[gene] = pd.Series(rhos, index=universe)
        kind = spec["test_kind"]
        elig = eligible_mask(gmeans, contrast_groups(kind))
        effect = effect_vector(logx, kind)
        for anchor_gene, rho in anchor.items():
            for direction, sel in (
                ("anti", rho <= -0.4),
                ("co", rho >= 0.4),
            ):
                genes = [g for g in rho.index[sel] if g != anchor_gene and bool(elig.get(g, False))]
                fc = effect.reindex(genes).dropna()
                tested = gene_tests(fc)
                rows.append(
                    {
                        "split": split,
                        "anchor": anchor_gene,
                        "subset": direction,
                        "rho_rule": "<= -0.4" if direction == "anti" else ">= 0.4",
                        "n_selected_before_floor": int(sel.sum()),
                        "perturbation": "TACSTD2 shRNA, not CLDN4 KD",
                        "test_estimand": kind,
                        **tested,
                    }
                )
    return pd.DataFrame(rows)


def bh(df: pd.DataFrame, col: str, out_col: str) -> pd.DataFrame:
    p = df[col]
    q = np.full(len(df), np.nan)
    mask = p.notna()
    if mask.sum():
        q[mask.to_numpy()] = multipletests(p[mask], method="fdr_bh")[1]
    df[out_col] = q
    return df


def plot_sweep(sweep: pd.DataFrame) -> None:
    import matplotlib.pyplot as plt

    ref = sweep[sweep["normalization"] == "log2_fpkm_p1"].set_index(["estimand", "set"])
    left = [
        ("C4_panel", "C4 panel"),
        ("APM", "APM"),
        ("APM_noHLA", "APM without HLA-A"),
        ("IFN_custom", "IFN/ISG"),
        ("H_IFNA", "Hallmark IFN-α"),
        ("H_IFNG", "Hallmark IFN-γ"),
    ]
    right = [
        ("H_IFNA", "Hallmark IFN-α", "#72B7B2"),
        ("IFN_custom", "IFN/ISG", "#72B7B2"),
        ("APM_noHLA", "APM without HLA-A", "#72B7B2"),
        ("APM", "APM", "#72B7B2"),
        ("H_IFNG", "Hallmark IFN-γ", "#72B7B2"),
        ("H_E2F", "E2F targets", "#54A24B"),
        ("H_G2M", "G2/M", "#54A24B"),
        ("H_EMT", "EMT", "#F58518"),
        ("EPITHELIAL_noTarget", "Epithelial identity", "#F58518"),
        ("HOUSEKEEPING", "Housekeeping", "#4C78A8"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 5.2), constrained_layout=True)

    ax = axes[0]
    vals = [float(ref.loc[("KD_noDDP", name), "median_effect"]) for name, _ in left]
    ax.barh(range(len(left)), vals, color="#E45756", edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(left)))
    ax.set_yticklabels([lab for _, lab in left])
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_xlabel("median log2FC  log2(FPKM+1)")
    ax.set_title("Untreated TACSTD2 sh vs scramble\n0/104 IFN/APM tests have median > 0")
    ax.invert_yaxis()

    ax = axes[1]
    vals = [float(ref.loc[("KD_on_DDP", name), "median_effect"]) for name, _, _ in right]
    colors = [c for _, _, c in right]
    ax.barh(range(len(right)), vals, color=colors, edgecolor="black", linewidth=0.4)
    ax.set_yticks(range(len(right)))
    ax.set_yticklabels([lab for _, lab, _ in right])
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_xlabel("median log2FC  log2(FPKM+1)")
    ax.set_title("TACSTD2 sh vs scramble, both on cisplatin\nIFN/APM up; E2F/G2M and epithelium are not")
    ax.invert_yaxis()

    fig.suptitle("GSE245459 SKOV3 — TACSTD2/TROP2 shRNA, not a CLDN4 knockdown", fontsize=12)
    fig.savefig(FIGS / "thesis_up_sweep.png", dpi=160)
    fig.savefig(FIGS / "thesis_up_sweep.pdf")
    plt.close(fig)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    base.ensure_raw()
    fpkm = base.load_fpkm()
    gmeans = group_means(fpkm)
    hallmark = load_hallmark()
    family = ifn_family(hallmark)
    context = context_sets(hallmark)
    panels = {**family, **context}
    mats = transforms(fpkm)
    log_ref = mats["log2_fpkm_p1"]

    corr, pca = sample_qc(fpkm, log_ref)
    elig_untreated = eligible_mask(gmeans, ["shTACSTD2", "shNC"])
    loo = leave_one_out(log_ref, family, elig_untreated)
    loo.to_csv(TABLES / "sweep_leave_one_out.tsv", sep="\t", index=False)
    elig_ddp = eligible_mask(gmeans, ["shTACSTD2_DDP", "shNC_DDP"])
    loo_ddp = leave_one_out_on_ddp(log_ref, family, elig_ddp)
    loo_ddp.to_csv(TABLES / "sweep_leave_one_out_on_ddp.tsv", sep="\t", index=False)
    coreg = heldout_coreg(log_ref, fpkm, family, gmeans)
    coreg.to_csv(TABLES / "sweep_coreg_heldout.tsv", sep="\t", index=False)

    kinds = ["KD_noDDP", "KD_on_DDP", "interaction", "DDP_in_NC", "DDP_in_KD"]
    rows = []
    for norm, mat in mats.items():
        for kind in kinds:
            elig = eligible_mask(gmeans, contrast_groups(kind))
            if norm == "log2_fpkm_positive":
                elig = positive_only_eligible(fpkm, kind, elig)
            for set_name, genes in panels.items():
                rec = test_panel(mat, genes, kind, elig, fpkm)
                rows.append(
                    {
                        "perturbation": "TACSTD2/TROP2 shRNA, not CLDN4 KD",
                        "cell_line": "SKOV3",
                        "normalization": norm,
                        "estimand": kind,
                        "thesis_estimand": kind in THESIS_ESTIMANDS,
                        "family": "IFN_APM" if set_name in family else "context",
                        "set": set_name,
                        **rec,
                    }
                )
    sweep = pd.DataFrame(rows)
    # FDR inside thesis IFN/APM tests, and separately inside untreated IFN/APM tests.
    thesis = sweep["thesis_estimand"] & (sweep["family"] == "IFN_APM")
    sweep["q_wilcoxon_thesis_IFN_APM"] = np.nan
    sweep.loc[thesis, "q_wilcoxon_thesis_IFN_APM"] = bh(
        sweep.loc[thesis].copy(), "wilcoxon_p_up", "q"
    )["q"].to_numpy()
    untreated = thesis & (sweep["estimand"] == "KD_noDDP")
    sweep["q_wilcoxon_untreated_IFN_APM"] = np.nan
    sweep.loc[untreated, "q_wilcoxon_untreated_IFN_APM"] = bh(
        sweep.loc[untreated].copy(), "wilcoxon_p_up", "q"
    )["q"].to_numpy()
    sweep.to_csv(TABLES / "sweep_all.tsv", sep="\t", index=False)

    up = sweep[thesis & sweep["thesis_aligned_median"]].sort_values(
        ["wilcoxon_p_up", "sample_p_up"], ascending=True
    )
    up.to_csv(TABLES / "sweep_thesis_up.tsv", sep="\t", index=False)

    # Sanity: log2(FPKM+1) untreated IFN median should match the prior slice.
    chk = sweep[(sweep.normalization == "log2_fpkm_p1") & (sweep.estimand == "KD_noDDP") & (sweep.set == "IFN_custom")]
    summary = {
        "perturbation": "TACSTD2/TROP2 shRNA in SKOV3. Not a CLDN4 knockdown.",
        "n_sweep_rows": int(len(sweep)),
        "n_thesis_IFN_APM_tests": int(thesis.sum()),
        "n_untreated_IFN_APM_tests": int(untreated.sum()),
        "n_untreated_IFN_APM_median_up": int((untreated & sweep["thesis_aligned_median"]).sum()),
        "n_thesis_aligned_up_panels": int(len(up)),
        "untreated_IFN_custom_log2p1_median": None if chk.empty else float(chk.iloc[0]["median_effect"]),
        "loo_any_core_median_positive": bool(loo.loc[loo["loo"] != "none", "sign_positive"].any()),
        "loo_biased_double_drop_any_positive": bool(
            loo.loc[loo["biased_toward_up"], "sign_positive"].any()
        ),
        "loo_on_ddp_any_core_median_negative": bool((~loo_ddp["sign_positive"]).any()),
    }
    # Direct knockdown-versus-matched-control, reference normalization.
    # Interaction is stored separately: it subtracts the untreated decrease
    # and is not an IFN-specific level increase.
    direct = sweep[
        (sweep["normalization"] == "log2_fpkm_p1")
        & (sweep["estimand"] == "KD_on_DDP")
        & (sweep["family"] == "IFN_APM")
        & sweep["thesis_aligned_median"]
        & sweep["wilcoxon_p_up"].notna()
    ].sort_values("wilcoxon_p_up")
    if len(direct):
        best = direct.iloc[0]
        summary["best_direct_up_log2p1"] = {
            "normalization": "log2_fpkm_p1",
            "estimand": "KD_on_DDP",
            "label": "TACSTD2 sh vs scramble, both on cisplatin 10 µg/ml 48 h. Not CLDN4 KD. Not the untreated contrast.",
            "set": str(best["set"]),
            "n": int(best["n"]),
            "median_effect": float(best["median_effect"]),
            "n_up": int(best["n_up"]),
            "n_down": int(best["n_down"]),
            "wilcoxon_p_up": float(best["wilcoxon_p_up"]),
            "sample_p_up": float(best["sample_p_up"]),
            "mw_p_gt_background": float(best["mw_p_gt_background"]),
            "q_wilcoxon_thesis_IFN_APM": float(best["q_wilcoxon_thesis_IFN_APM"]),
        }
    if len(up):
        top = up.iloc[0]
        summary["smallest_p_in_full_sweep_including_interaction"] = {
            "normalization": str(top["normalization"]),
            "estimand": str(top["estimand"]),
            "set": str(top["set"]),
            "median_effect": float(top["median_effect"]),
            "wilcoxon_p_up": float(top["wilcoxon_p_up"]),
            "note": "Interaction contrasts subtract the untreated decrease. EMT shows a similar interaction. Not used as the IFN level call.",
        }
    (TABLES / "sweep_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    # Gene-level detail for the reference normalization on the three thesis estimands.
    gene_rows = []
    detail_sets = {
        "C4_panel": family["C4_panel"],
        "APM": family["APM"],
        "H_IFNA": family["H_IFNA"],
        "H_IFNG": family["H_IFNG"],
        "sentinel": ["CD274", "STAT1", "IRF1", "NLRC5", "IFIT1", "HLA-A", "ISG15", "B2M", "TAP1"],
    }
    for kind in THESIS_ESTIMANDS:
        elig = eligible_mask(gmeans, contrast_groups(kind))
        effect = effect_vector(log_ref, kind)
        for set_name, genes in detail_sets.items():
            for g in genes:
                if g not in effect.index:
                    continue
                gene_rows.append(
                    {
                        "estimand": kind,
                        "set": set_name,
                        "gene": g,
                        "eligible": bool(elig.get(g, False)),
                        "effect_log2p1": float(effect.at[g]) if np.isfinite(effect.at[g]) else None,
                        "mean_fpkm_shNC": float(gmeans.at[g, "shNC"]) if g in gmeans.index else None,
                        "mean_fpkm_sh": float(gmeans.at[g, "shTACSTD2"]) if g in gmeans.index else None,
                        "mean_fpkm_shNC_DDP": float(gmeans.at[g, "shNC_DDP"]) if g in gmeans.index else None,
                        "mean_fpkm_sh_DDP": float(gmeans.at[g, "shTACSTD2_DDP"]) if g in gmeans.index else None,
                    }
                )
    pd.DataFrame(gene_rows).to_csv(TABLES / "sweep_gene_effects.tsv", sep="\t", index=False)

    plot_sweep(sweep)
    print(json.dumps(summary, indent=2, default=str))
    print("--- untreated IFN/APM medians, log2p1 ---")
    sub = sweep[(sweep.normalization == "log2_fpkm_p1") & (sweep.estimand == "KD_noDDP") & (sweep.family == "IFN_APM")]
    print(sub[["set", "n", "median_effect", "n_up", "n_down", "wilcoxon_p_up", "sample_p_up"]].to_string(index=False))
    print("--- top thesis-aligned UP ---")
    print(up.head(15)[["normalization", "estimand", "set", "n", "median_effect", "n_up", "n_down", "wilcoxon_p_up", "sample_p_up", "q_wilcoxon_thesis_IFN_APM"]].to_string(index=False))
    print("--- coreg ---")
    print(coreg.to_string(index=False))
    print("--- loo sign flips ---")
    print(loo.groupby("loo")["sign_positive"].sum().to_string())
    _ = corr, pca


if __name__ == "__main__":
    main()
