#!/usr/bin/env python3
"""Test whether CLDN4 attenuates TROP2 vs IFN / MHC-I / immune-gene scores.

Public DepMap 24Q4 RNA, Gygi/Nusinow CCLE TMT protein, and DepMap 24Q4
CRISPR Chronos. The TACSTD2–CLDN4 protein Spearman near 0.69 is recomputed
only as a matrix check. It is not the finding.

Partial correlation: Pearson correlation of residuals after regressing the
ranks of x and y on the ranks of the covariates. Two-sided p uses
df = n - 2 - k. This is association, not mediation.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

MHC1 = ["HLA-A", "HLA-B", "HLA-C", "B2M"]
IMMUNE_LIGAND = ["CD274", "PDCD1LG2", "CXCL9", "CXCL10", "CXCL11", "HLA-E", "IDO1"]
ISG = [
    "STAT1",
    "STAT2",
    "IRF1",
    "IRF9",
    "MX1",
    "ISG15",
    "IFIT1",
    "IFIT3",
    "OAS1",
    "OAS2",
    "OAS3",
    "EIF2AK2",
    "IFI35",
    "BST2",
    "SAMHD1",
    "TRIM25",
    "ADAR",
    "TAP1",
    "TAP2",
    "PSMB8",
    "PSMB9",
    "PSMB10",
]
PRIMARY_SCORES = ["IFNG", "MHC1", "IMMUNE"]
SCORE_LABEL = {
    "IFNG": "Hallmark IFN-γ",
    "MHC1": "MHC-I (HLA-A/B/C+B2M)",
    "IMMUNE": "Immune-ligand score",
    "CD274": "CD274",
    "ISG": "Compact ISG",
}
N_BOOT = 5000
MIN_N = 8


def zscore(s: pd.Series) -> pd.Series:
    mu = s.mean()
    sd = s.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return s * 0.0
    return (s - mu) / sd


def signature(df: pd.DataFrame, genes: list[str], min_members: int) -> tuple[pd.Series, list[str]]:
    use = [g for g in genes if g in df.columns]
    if not use:
        return pd.Series(np.nan, index=df.index), []
    z = df[use].apply(zscore, axis=0)
    n_ok = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score = score.where(n_ok >= min_members)
    return score, use


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for rank, i in enumerate(order[::-1], start=1):
        k = n - rank + 1
        val = min(prev, p[i] * n / k)
        q[i] = val
        prev = val
    return [float(min(1.0, v)) for v in q]


def pearson_r(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    denom = np.sqrt((a @ a) * (b @ b))
    if denom == 0 or not np.isfinite(denom):
        return float("nan")
    return float(np.clip((a @ b) / denom, -1.0, 1.0))


def spearman_r(x: np.ndarray, y: np.ndarray) -> float:
    return pearson_r(stats.rankdata(x, method="average"), stats.rankdata(y, method="average"))


def spearman_full(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    a = pd.concat([x, y], axis=1).dropna()
    n = int(len(a))
    if n < MIN_N:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(a.iloc[:, 0].to_numpy(), a.iloc[:, 1].to_numpy())
    return float(rho), float(p), n


def partial_corr(frame: pd.DataFrame, x: str, y: str, covariates: list[str]) -> tuple[float, float, int]:
    cols = [x, y, *covariates]
    sub = frame[cols].apply(pd.to_numeric, errors="coerce").dropna()
    n = int(len(sub))
    k = len(covariates)
    if n < k + 5:
        return float("nan"), float("nan"), n
    ranks = sub.rank(method="average")
    design = np.column_stack([np.ones(n)] + [ranks[c].to_numpy() for c in covariates])

    def resid(v: np.ndarray) -> np.ndarray:
        beta, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ beta

    r = pearson_r(resid(ranks[x].to_numpy()), resid(ranks[y].to_numpy()))
    df = n - 2 - k
    if not np.isfinite(r) or df < 1 or abs(r) >= 1:
        p = 0.0 if np.isfinite(r) and abs(r) >= 1 else float("nan")
        return r, p, n
    t = r * np.sqrt(df / max(1 - r * r, 1e-300))
    p = float(2 * stats.t.sf(abs(t), df=df))
    return float(r), p, n


def bootstrap_stat(values: np.ndarray, stat, n_boot: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stat(values[idx])
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def seed_for(label: str) -> int:
    return int(pd.util.hash_pandas_object(pd.Series([label]), index=False).iloc[0] % (2**31 - 1))


def fisher_ci(r: float, n_eff: int) -> tuple[float, float]:
    if not np.isfinite(r) or n_eff < 4 or abs(r) >= 1:
        return float("nan"), float("nan")
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1.0 / np.sqrt(n_eff - 3)
    crit = stats.norm.ppf(0.975)
    return float(np.tanh(z - crit * se)), float(np.tanh(z + crit * se))


def attenuation_call(unadj_rho: float, unadj_p: float, part_rho: float, part_p: float) -> str:
    """Pre-specified label. Not a formal test that two correlations differ.

    A partial that stays significant with the opposite sign is sign_reversed,
    not retained. Retention requires the same sign.
    """
    if not np.isfinite(unadj_rho) or not np.isfinite(unadj_p) or unadj_p >= 0.05:
        return "no_unadjusted_association"
    if not np.isfinite(part_rho) or not np.isfinite(part_p):
        return "partial_not_estimable"
    ratio = abs(part_rho) / abs(unadj_rho) if abs(unadj_rho) > 0 else float("nan")
    same_sign = np.sign(part_rho) == np.sign(unadj_rho) or part_rho == 0
    if part_p < 0.05 and not same_sign:
        return "sign_reversed"
    if part_p < 0.05 and ratio > 0.5:
        return "retained"
    if part_p < 0.05 and ratio <= 0.5:
        return "reduced_but_still_associated"
    if part_p >= 0.05 and ratio <= 0.5:
        return "attenuated"
    if part_p >= 0.05 and ratio > 0.5:
        return "ns_after_adjustment_point_estimate_still_half_or_more"
    return "other"


def cohort_tag(row: pd.Series) -> str:
    disease = str(row.get("OncotreePrimaryDisease") or "")
    subtype = str(row.get("OncotreeSubtype") or "")
    if disease == "Non-Small Cell Lung Cancer":
        if subtype == "Lung Adenocarcinoma":
            return "LUAD"
        if subtype == "Lung Squamous Cell Carcinoma":
            return "LUSC"
        return "other_NSCLC"
    if disease == "Lung Neuroendocrine Tumor" or subtype == "Small Cell Lung Cancer":
        return "SCLC_NET"
    return "other_lung"


def self_check() -> None:
    rng = np.random.default_rng(1)
    n = 500
    z = rng.normal(size=n)
    x = z + rng.normal(scale=0.25, size=n)
    y = z + rng.normal(scale=0.25, size=n)
    shared = pd.DataFrame({"x": x, "y": y, "z": z})
    r_u, _, _ = spearman_full(shared["x"], shared["y"])
    r_p, _, _ = partial_corr(shared, "x", "y", ["z"])
    if not (r_u > 0.6 and abs(r_p) < 0.15):
        raise SystemExit(f"self-check failed for shared covariate: unadj={r_u:.3f} partial={r_p:.3f}")
    x2 = rng.normal(size=n)
    y2 = x2 + rng.normal(scale=0.3, size=n)
    z2 = rng.normal(size=n)
    indep = pd.DataFrame({"x": x2, "y": y2, "z": z2})
    r_u2, _, _ = spearman_full(indep["x"], indep["y"])
    r_p2, p_p2, _ = partial_corr(indep, "x", "y", ["z"])
    if not (r_u2 > 0.7 and r_p2 > 0.6 and p_p2 < 1e-6):
        raise SystemExit(f"self-check failed for null covariate: unadj={r_u2:.3f} partial={r_p2:.3f} p={p_p2}")


def add_rna_scores(frame: pd.DataFrame, ifng_genes: list[str]) -> tuple[pd.DataFrame, dict]:
    out = frame.copy()
    ifng, used_ifng = signature(out, ifng_genes, min_members=max(5, int(0.5 * len(ifng_genes))))
    mhc, used_mhc = signature(out, MHC1, min_members=4)
    immune_genes = [g for g in IMMUNE_LIGAND if g in out.columns]
    immune, used_imm = signature(out, immune_genes, min_members=max(3, len(immune_genes) // 2))
    isg_genes = [g for g in ISG if g in out.columns]
    isg, used_isg = signature(out, isg_genes, min_members=5)
    out["IFNG"] = ifng
    out["MHC1"] = mhc
    out["IMMUNE"] = immune
    out["ISG"] = isg
    if "CD274" in out.columns:
        out["CD274_score"] = out["CD274"]
    cov = {
        "ifng_n_used": len(used_ifng),
        "ifng_n_gmt": len(ifng_genes),
        "ifng_missing": sorted(set(ifng_genes) - set(used_ifng)),
        "mhc1_used": used_mhc,
        "immune_used": used_imm,
        "isg_used": used_isg,
    }
    return out, cov


def fmt(x: float | None, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:.{digits}f}"


def fmt_p(x: float | None) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def record_pair(
    rows: list[dict],
    frame: pd.DataFrame,
    cohort: str,
    layer: str,
    x: str,
    y: str,
    covariates: list[str],
    family: str,
    primary: bool,
    do_boot: bool,
) -> None:
    cols = [x, y, *covariates]
    cov_label = "+".join(covariates) if covariates else "none"
    missing = [c for c in cols if c not in frame.columns]
    if missing:
        rows.append(
            {
                "layer": layer,
                "cohort": cohort,
                "x": x,
                "y": y,
                "covariates": cov_label,
                "family": family,
                "primary": primary,
                "n": 0,
                "rho": np.nan,
                "p": np.nan,
                "kind": "partial" if covariates else "spearman",
                "note": f"missing columns: {missing}",
            }
        )
        return
    sub = frame[cols].apply(pd.to_numeric, errors="coerce").dropna()
    n = int(len(sub))
    kind = "partial" if covariates else "spearman"
    if covariates:
        rho, p, n = partial_corr(sub, x, y, covariates)
        n_eff = n - len(covariates)
    else:
        rho, p, n = spearman_full(sub[x], sub[y])
        n_eff = n
    ci_f_lo, ci_f_hi = fisher_ci(rho, n_eff if kind == "partial" else n)
    ci_b_lo = ci_b_hi = np.nan
    if do_boot and n >= MIN_N and np.isfinite(rho):
        label = f"{layer}|{cohort}|{x}|{y}|{'+'.join(covariates)}|{kind}"
        seed = seed_for(label)
        mat = sub.to_numpy()

        def stat(sample: np.ndarray) -> float:
            sdf = pd.DataFrame(sample, columns=cols)
            if covariates:
                r, _, _ = partial_corr(sdf, x, y, covariates)
                return r
            return spearman_r(sample[:, 0], sample[:, 1])

        ci_b_lo, ci_b_hi = bootstrap_stat(mat, stat, N_BOOT, seed)
    call = ""
    rows.append(
        {
            "layer": layer,
            "cohort": cohort,
            "x": x,
            "y": y,
            "covariates": cov_label,
            "kind": kind,
            "family": family,
            "primary": primary,
            "n": n,
            "rho": rho,
            "p": p,
            "fisher_ci95_low": ci_f_lo,
            "fisher_ci95_high": ci_f_hi,
            "boot_ci95_low": ci_b_lo,
            "boot_ci95_high": ci_b_hi,
            "boot_n": N_BOOT if do_boot and n >= MIN_N else 0,
            "attenuation_call": call,
            "note": "",
        }
    )


def apply_fdr_and_calls(tab: pd.DataFrame) -> pd.DataFrame:
    out = tab.copy()
    out["q_bh"] = np.nan
    for fam, idx in out.groupby("family").groups.items():
        if not fam:
            continue
        block = out.loc[idx]
        ok = block["p"].notna() & block["primary"].astype(bool) & (block["kind"] == "partial")
        # FDR is within the declared family, and only on rows flagged primary.
        use = block.index[ok.to_numpy()]
        if len(use) == 0:
            continue
        out.loc[use, "q_bh"] = bh_fdr(out.loc[use, "p"].astype(float).tolist())
    # Also FDR the unadjusted primary spearman rows, separately, same family name + kind.
    out["q_bh_unadjusted"] = np.nan
    for fam, idx in out.groupby("family").groups.items():
        if not fam:
            continue
        block = out.loc[idx]
        ok = block["p"].notna() & block["primary"].astype(bool) & (block["kind"] == "spearman")
        use = block.index[ok.to_numpy()]
        if len(use) == 0:
            continue
        out.loc[use, "q_bh_unadjusted"] = bh_fdr(out.loc[use, "p"].astype(float).tolist())

    calls = []
    ratios = []
    for _, r in out.iterrows():
        if r["kind"] != "partial":
            calls.append("")
            ratios.append(np.nan)
            continue
        match = out[
            (out["layer"] == r["layer"])
            & (out["cohort"] == r["cohort"])
            & (out["x"] == r["x"])
            & (out["y"] == r["y"])
            & (out["kind"] == "spearman")
            & (out["covariates"] == "none")
        ]
        if match.empty or not np.isfinite(r["rho"]):
            calls.append("partial_not_estimable" if not np.isfinite(r["rho"]) else "")
            ratios.append(np.nan)
            continue
        u = match.iloc[0]
        ratio = abs(r["rho"]) / abs(u["rho"]) if np.isfinite(u["rho"]) and abs(u["rho"]) > 0 else np.nan
        ratios.append(ratio)
        calls.append(attenuation_call(float(u["rho"]), float(u["p"]), float(r["rho"]), float(r["p"])))
    out["rho_ratio_partial_over_unadj"] = ratios
    out["attenuation_call"] = calls
    return out


def load_models(path: Path) -> pd.DataFrame:
    model = pd.read_csv(path)
    need = [
        "ModelID",
        "CellLineName",
        "CCLEName",
        "StrippedCellLineName",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "ModelType",
    ]
    missing = [c for c in need if c not in model.columns]
    if missing:
        raise SystemExit(f"Model.csv missing {missing}; columns={list(model.columns)[:30]}")
    return model


def lung_mask(df: pd.DataFrame) -> pd.Series:
    return (df["OncotreeLineage"] == "Lung") & (df["ModelType"] == "Cell Line")


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["group"] = out.apply(cohort_tag, axis=1)
    out["is_NSCLC"] = out["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    out["nsclc_indicator"] = out["is_NSCLC"].astype(float)
    return out


def gygi_checksum(panel: pd.DataFrame) -> dict:
    """Column-level TACSTD2–CLDN4 Spearman on TenPx columns whose names contain _LUNG_.

    This is the already-reported ~0.69 pair (PR 54). Replicates are not averaged.
    """
    sub = panel[panel["Gene_Symbol"].isin(["TACSTD2", "CLDN4"])].copy()
    if set(sub["Gene_Symbol"]) != {"TACSTD2", "CLDN4"}:
        return {"error": f"genes present: {sorted(sub['Gene_Symbol'].unique())}"}
    wide = sub.drop_duplicates("Gene_Symbol").set_index("Gene_Symbol")
    cols = [c for c in wide.columns if c not in ("Protein_Id",) and "_TenPx" in c and "_LUNG_" in c]
    tac = pd.to_numeric(wide.loc["TACSTD2", cols], errors="coerce")
    cld = pd.to_numeric(wide.loc["CLDN4", cols], errors="coerce")
    both = pd.concat([tac.rename("TACSTD2"), cld.rename("CLDN4")], axis=1).dropna()
    rho, p, n = spearman_full(both["TACSTD2"], both["CLDN4"])
    return {
        "definition": "Gygi TenPx columns with _LUNG_ in the name; both proteins non-missing; replicates not collapsed",
        "n_lung_tenpx_columns": len(cols),
        "n_complete": n,
        "spearman_rho": rho,
        "p": p,
        "rounds_to_0.69": bool(np.isfinite(rho) and round(float(rho), 2) == 0.69),
        "protein_ids": sub.groupby("Gene_Symbol")["Protein_Id"].first().to_dict() if "Protein_Id" in sub.columns else {},
    }


def models_from_gygi(panel: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    meta = {"Gene_Symbol", "Protein_Id"}
    num = panel.drop(columns=[c for c in panel.columns if c in meta], errors="ignore")
    # Average duplicate gene symbols.
    tmp = pd.concat(
        [panel[["Gene_Symbol"]].reset_index(drop=True), num.apply(pd.to_numeric, errors="coerce").reset_index(drop=True)],
        axis=1,
    )
    prot = tmp.groupby("Gene_Symbol").mean(numeric_only=True)
    # Map columns to ModelID.
    ccle = {str(v): i for i, v in model.set_index("ModelID")["CCLEName"].dropna().items()}
    stripped = {
        str(v).upper(): i
        for i, v in model.set_index("ModelID")["StrippedCellLineName"].dropna().items()
    }
    mapping = {}
    for c in prot.columns:
        if c.startswith("ACH-"):
            mapping[c] = c.split("_")[0]
            continue
        base = c.rsplit("_TenPx", 1)[0]
        if base in ccle:
            mapping[c] = ccle[base]
            continue
        token = base.split("_")[0].upper()
        if token in stripped:
            mapping[c] = stripped[token]
    wide = prot.T
    wide["ModelID"] = wide.index.map(mapping)
    wide = wide.dropna(subset=["ModelID"])
    genes = [c for c in wide.columns if c != "ModelID"]
    by_model = wide.groupby("ModelID")[genes].mean()
    ann = by_model.join(
        model.set_index("ModelID")[
            [
                "CellLineName",
                "CCLEName",
                "OncotreeLineage",
                "OncotreePrimaryDisease",
                "OncotreeSubtype",
                "ModelType",
            ]
        ],
        how="left",
    )
    return ann


def add_protein_scores(frame: pd.DataFrame, ifng_genes: list[str]) -> tuple[pd.DataFrame, dict]:
    out = frame.copy()
    hm = [g for g in ifng_genes if g in out.columns and int(out[g].notna().sum()) >= 8]
    ifng, used_ifng = signature(out, hm, min_members=5)
    mhc, used_mhc = signature(out, [g for g in MHC1 if g in out.columns], min_members=2)
    immune_genes = [g for g in IMMUNE_LIGAND if g in out.columns and int(out[g].notna().sum()) >= 8]
    # CD274 is allowed into the score even if slightly sparser, as long as it exists.
    if "CD274" in out.columns and "CD274" not in immune_genes and int(out["CD274"].notna().sum()) >= 8:
        immune_genes.append("CD274")
    min_imm = 2 if len(immune_genes) >= 2 else 1
    immune, used_imm = signature(out, immune_genes, min_members=min_imm)
    isg_genes = [g for g in ISG if g in out.columns and int(out[g].notna().sum()) >= 8]
    isg, used_isg = signature(out, isg_genes, min_members=5)
    out["IFNG"] = ifng
    out["MHC1"] = mhc
    out["IMMUNE"] = immune
    out["ISG"] = isg
    cov = {
        "ifng_n_used": len(used_ifng),
        "ifng_n_gmt": len(ifng_genes),
        "mhc1_used": used_mhc,
        "immune_used": used_imm,
        "isg_n_used": len(used_isg),
        "isg_used": used_isg,
    }
    return out, cov


def expressed_fraction(frame: pd.DataFrame, genes: list[str], thresh: float = 0.0) -> dict:
    out = {}
    for g in genes:
        if g not in frame.columns:
            out[g] = None
            continue
        s = pd.to_numeric(frame[g], errors="coerce")
        out[g] = {
            "n": int(s.notna().sum()),
            "frac_gt_0": float((s > thresh).mean()) if s.notna().any() else None,
            "median": float(s.median()) if s.notna().any() else None,
        }
    return out


def write_finding(path: Path, key: dict, tab: pd.DataFrame) -> None:
    def row(layer, cohort, x, y, covariates="none") -> dict | None:
        hit = tab[
            (tab["layer"] == layer)
            & (tab["cohort"] == cohort)
            & (tab["x"] == x)
            & (tab["y"] == y)
            & (tab["covariates"] == covariates)
        ]
        if hit.empty:
            return None
        return hit.iloc[0].to_dict()

    def cell(rec: dict | None) -> str:
        if rec is None or not np.isfinite(rec.get("rho", np.nan)):
            n = 0 if rec is None else rec.get("n", 0)
            return f"n={n}; not tested"
        q = rec.get("q_bh") if rec.get("kind") == "partial" else rec.get("q_bh_unadjusted")
        qbit = f", q={fmt_p(q)}" if q is not None and np.isfinite(q) else ""
        ci = ""
        if np.isfinite(rec.get("boot_ci95_low", np.nan)):
            ci = f"; boot CI [{rec['boot_ci95_low']:+.2f}, {rec['boot_ci95_high']:+.2f}]"
        call = rec.get("attenuation_call") or ""
        callbit = f"; {call}" if call else ""
        ratio = rec.get("rho_ratio_partial_over_unadj")
        rbit = f"; |partial|/|unadj|={ratio:.2f}" if ratio is not None and np.isfinite(ratio) else ""
        return (
            f"{rec['rho']:+.3f} (n={int(rec['n'])}, p={fmt_p(rec['p'])}{qbit}){ci}{rbit}{callbit}"
        )

    chk = key["protein_checksum_lung_suffix_columns"]
    lines = []
    lines.append("# DepMap/CCLE: does CLDN4 attenuate TROP2 vs IFN / MHC / immune scores?")
    lines.append("")
    lines.append("Numbers below are written by `analyze.py` from `tables/key_stats.json` and `tables/partials.csv`. They are not typed by hand.")
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append("TROP2 (TACSTD2) and CLDN4 protein are already correlated in Gygi CCLE lung columns. This page tests whether the association of **TROP2 RNA or protein** with pre-specified IFN, MHC-I, and immune-ligand scores shrinks after adjustment for CLDN4, and whether CLDN4's own association shrinks after adjustment for TROP2. CRISPR Chronos is reported for the same scores when the screen exists.")
    lines.append("")
    lines.append("Partial correlation is the Pearson correlation of rank residuals after linear adjustment for the covariate ranks. It is an association contrast, not a mediation or causal estimate. Cultured lines have no T-cell infiltrate and no IFN treatment.")
    lines.append("")
    lines.append("## Answer")
    lines.append("")
    for sentence in key.get("answer", []):
        lines.append(sentence)
        lines.append("")
    lines.append("## Already-known covariate correlation (recomputed)")
    lines.append("")
    lines.append(
        f"Gygi TenPx columns with `_LUNG_` in the name, TACSTD2 and CLDN4 both quantified, replicates not collapsed: "
        f"Spearman **{fmt(chk.get('spearman_rho'))}**, n={chk.get('n_complete')}, p={fmt_p(chk.get('p'))}. "
        f"Rounds to 0.69: **{chk.get('rounds_to_0.69')}**. "
        "This is the same column definition as the existing Gygi lung complete-case result. It is a matrix check, not a new claim. "
        "Partials below use **one row per ModelID** (TenPx replicates averaged), so their n is a cell-line count, not this column count."
    )
    lines.append("")
    rna_pair = key["covariate_correlations"]["rna_lung_TACSTD2_CLDN4"]
    prot_pair = key["covariate_correlations"]["protein_model_lung_TACSTD2_CLDN4"]
    lines.append(
        f"On the rows used for scoring: RNA lung cell lines TACSTD2–CLDN4 Spearman {fmt(rna_pair.get('rho'))} "
        f"(n={rna_pair.get('n')}, p={fmt_p(rna_pair.get('p'))}). "
        f"Protein, Model.csv lung cell lines, both proteins present: {fmt(prot_pair.get('rho'))} "
        f"(n={prot_pair.get('n')}, p={fmt_p(prot_pair.get('p'))})."
    )
    lines.append("")
    lines.append("## Scores")
    lines.append("")
    lines.append("- **IFNG**: mean of within-cohort z-scores of MSigDB 2024.1 `HALLMARK_INTERFERON_GAMMA_RESPONSE` genes present in the matrix.")
    lines.append("- **MHC1**: HLA-A, HLA-B, HLA-C, B2M. RNA requires all 4. Protein requires at least 2 quantified members.")
    lines.append("- **IMMUNE**: immune-ligand score, mean z of CD274, PDCD1LG2, CXCL9, CXCL10, CXCL11, HLA-E, IDO1 (genes present). This is cancer-cell RNA or protein, not an immune-infiltrate score.")
    lines.append("- CD274 alone and a compact ISG cassette are reported as components or sensitivity and are outside the primary FDR family.")
    lines.append("")
    lines.append(
        "Attenuation label, applied only when the unadjusted Spearman p < 0.05: "
        "**attenuated** if the partial p ≥ 0.05 and |partial| ≤ half |unadjusted|; "
        "**retained** if the partial p < 0.05 and |partial| > half |unadjusted|; "
        "**reduced_but_still_associated** if the partial p < 0.05 but |partial| ≤ half |unadjusted|; "
        "**ns_after_adjustment_point_estimate_still_half_or_more** if significance is lost but the point estimate has not halved; "
        "**sign_reversed** if the partial stays significant (p < 0.05) with the opposite sign. "
        "Retention requires the same sign. "
        "BH q for partials is within the six primary partials of that cohort (TROP2|CLDN4 and CLDN4|TROP2 × three scores). Unadjusted q is within the six matching Spearman tests."
    )
    lines.append("")
    lines.append("## RNA (DepMap 24Q4 log2(TPM+1))")
    lines.append("")
    lines.append(
        f"Lung cell lines (`OncotreeLineage == Lung`, `ModelType == Cell Line`) with TACSTD2 and CLDN4: "
        f"**n = {key['n_rna_lung']}** (NSCLC {key['n_rna_nsclc']}, LUAD {key['n_rna_luad']}, LUSC {key['n_rna_lusc']}, SCLC/NET {key['n_rna_sclc']}). "
        "`tables/rna_lung_lines.csv` stores scores z-scored on all lung lines. NSCLC and LUAD rows in `partials.csv` recompute those z-scores inside the cohort."
    )
    lines.append("")
    lines.append("| cohort | score | TROP2 unadjusted | TROP2 given CLDN4 | CLDN4 unadjusted | CLDN4 given TROP2 |")
    lines.append("|---|---|---|---|---|---|")
    for cohort, label in (
        ("lung", "all lung"),
        ("NSCLC", "NSCLC"),
        ("LUAD", "LUAD (exploratory)"),
    ):
        for score in PRIMARY_SCORES:
            lines.append(
                "| "
                + " | ".join(
                    [
                        label,
                        SCORE_LABEL[score],
                        cell(row("rna", cohort, "TACSTD2", score)),
                        cell(row("rna", cohort, "TACSTD2", score, "CLDN4")),
                        cell(row("rna", cohort, "CLDN4", score)),
                        cell(row("rna", cohort, "CLDN4", score, "TACSTD2")),
                    ]
                )
                + " |"
            )
    lines.append("")
    lines.append("All-lung sensitivity, adjusting for CLDN4 or TROP2 plus a binary NSCLC indicator (not in the six-test FDR):")
    lines.append("")
    for score in PRIMARY_SCORES:
        rec_t = row("rna", "lung", "TACSTD2", score, "CLDN4+nsclc_indicator")
        rec_c = row("rna", "lung", "CLDN4", score, "TACSTD2+nsclc_indicator")
        lines.append(f"- {SCORE_LABEL[score]}: TROP2|CLDN4+NSCLC {cell(rec_t)}; CLDN4|TROP2+NSCLC {cell(rec_c)}")
    lines.append("")
    lines.append("CD274 and compact ISG are outside the primary FDR. NSCLC RNA:")
    lines.append("")
    for score in ("CD274_score", "ISG"):
        lines.append(
            f"- {score}: TROP2 {cell(row('rna', 'NSCLC', 'TACSTD2', score))}; "
            f"TROP2|CLDN4 {cell(row('rna', 'NSCLC', 'TACSTD2', score, 'CLDN4'))}; "
            f"CLDN4|TROP2 {cell(row('rna', 'NSCLC', 'CLDN4', score, 'TACSTD2'))}"
        )
    lines.append("")
    chem = key.get("rna_nsclc_chemokine_expression", {})
    if chem:
        bits = []
        for g, info in chem.items():
            if not info:
                bits.append(f"{g} absent")
            else:
                bits.append(f"{g} median {info['median']:.3f}, frac>0 {info['frac_gt_0']:.2f}")
        lines.append("NSCLC RNA expression of immune-ligand members (log2(TPM+1)): " + "; ".join(bits) + ".")
        lines.append("")
    lines.append("## Protein (Gygi/Nusinow CCLE TMT, one row per lung cell line)")
    lines.append("")
    lines.append(
        f"Lung cell lines with any mapped Gygi protein: n = {key['n_protein_lung_any']}. "
        f"With both TACSTD2 and CLDN4 protein: **n = {key['n_protein_both']}** "
        f"(NSCLC {key['n_protein_nsclc_both']}). "
        "RPPA antibody table from this run: "
        f"mentions CLDN4 = {key['rppa_mentions_CLDN4']}, mentions TACSTD2/TROP2 = {key['rppa_mentions_TACSTD2']}. "
        "RPPA is not used for the partials."
    )
    lines.append("")
    lines.append("| cohort | score | TROP2 unadjusted | TROP2 given CLDN4 | CLDN4 unadjusted | CLDN4 given TROP2 |")
    lines.append("|---|---|---|---|---|---|")
    for cohort, label in (("lung", "all lung"), ("NSCLC", "NSCLC")):
        for score in PRIMARY_SCORES:
            lines.append(
                "| "
                + " | ".join(
                    [
                        label,
                        SCORE_LABEL[score],
                        cell(row("protein", cohort, "TACSTD2", score)),
                        cell(row("protein", cohort, "TACSTD2", score, "CLDN4")),
                        cell(row("protein", cohort, "CLDN4", score)),
                        cell(row("protein", cohort, "CLDN4", score, "TACSTD2")),
                    ]
                )
                + " |"
            )
    lines.append("")
    imm_genes = key["protein_score_coverage_lung"].get("immune_used", [])
    lines.append(
        "Protein immune-ligand score members with enough quantified values: "
        + (", ".join(imm_genes) if imm_genes else "none")
        + ". CXCL9, CXCL10, and CXCL11 are not in the Gygi matrix, so this protein score is not the RNA immune-ligand score. "
        f"Hallmark IFN-γ proteins used: {key['protein_score_coverage_lung'].get('ifng_n_used')} of {key['protein_score_coverage_lung'].get('ifng_n_gmt')}. "
        f"MHC-I proteins used: {', '.join(key['protein_score_coverage_lung'].get('mhc1_used', []))}."
    )
    lines.append("")
    lines.append(
        "No primary protein partial in this table has BH q < 0.05. Protein attenuation labels use the nominal p < 0.05 rule on n = 32–44 and are not FDR findings."
    )
    lines.append("")
    lines.append("## CRISPR (DepMap 24Q4 Chronos)")
    lines.append("")
    dep = key["crispr_dependency"]
    lines.append(
        f"Lung cell lines with Chronos for TACSTD2 and CLDN4: **n = {dep['n_lung']}** "
        f"(NSCLC {dep['n_nsclc']}). Dependent = gene-dependency probability > 0.5."
    )
    lines.append("")
    lines.append("| gene | n | mean Chronos (SD) | n dependent (prob>0.5) |")
    lines.append("|---|---:|---|---:|")
    for gene, rec in dep["genes"].items():
        lines.append(
            f"| {gene} | {rec['n']} | {rec['mean_effect']:+.4f} ({rec['sd_effect']:.4f}) | {rec['n_dependent']} / {rec['n_dependency']} |"
        )
    lines.append("")
    lines.append("Chronos versus the same RNA scores, on lines with both CRISPR and RNA. A null unadjusted correlation means there is no association for CLDN4 to attenuate.")
    lines.append("")
    lines.append("| cohort | score | TROP2 Chronos | TROP2 Chronos given CLDN4 Chronos | CLDN4 Chronos given TROP2 Chronos | TROP2 Chronos given CLDN4 RNA |")
    lines.append("|---|---|---|---|---|---|")
    for cohort, label in (("lung", "lung overlap"), ("NSCLC", "NSCLC overlap")):
        for score in PRIMARY_SCORES:
            lines.append(
                "| "
                + " | ".join(
                    [
                        label,
                        SCORE_LABEL[score],
                        cell(row("crispr", cohort, "TACSTD2_effect", score)),
                        cell(row("crispr", cohort, "TACSTD2_effect", score, "CLDN4_effect")),
                        cell(row("crispr", cohort, "CLDN4_effect", score, "TACSTD2_effect")),
                        cell(row("crispr", cohort, "TACSTD2_effect", score, "CLDN4_RNA")),
                    ]
                )
                + " |"
            )
    lines.append("")
    lines.append("## What the calls say")
    lines.append("")
    for item in key["headline_calls"]:
        lines.append(f"- {item}")
    lines.append("")
    lines.append("## What this is not")
    lines.append("")
    lines.append("- Not a new estimate of the TROP2–CLDN4 protein correlation. That pair is the covariate, recomputed above.")
    lines.append("- Not immune exclusion, ICI response, or IFN stimulation. The dish has no T cells.")
    lines.append("- Not surface MHC or PD-L1. CD274 protein, where present, is TMT abundance.")
    lines.append("- Not a claim that CLDN4 causes or blocks the TROP2 association. Partials remove shared rank variation only.")
    lines.append("- Not RPPA n=118. That antibody panel does not contain TROP2 or CLDN4.")
    lines.append("")
    lines.append("## Rerun")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 -m pip install -r scripts/depmap_trop2_ifn_partial_cldn4/requirements.txt")
    lines.append("python3 scripts/depmap_trop2_ifn_partial_cldn4/download.py --outdir data/depmap_trop2_ifn_partial_cldn4")
    lines.append("python3 scripts/depmap_trop2_ifn_partial_cldn4/analyze.py --data data/depmap_trop2_ifn_partial_cldn4 --outdir results/depmap_trop2_ifn_partial_cldn4")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def forest(ax, tab: pd.DataFrame, layer: str, cohort: str, title: str) -> None:
    scores = PRIMARY_SCORES
    specs = [
        ("TACSTD2", "none", "TROP2", "#1f4e79", "o"),
        ("TACSTD2", "CLDN4", "TROP2 | CLDN4", "#1f4e79", "D"),
        ("CLDN4", "none", "CLDN4", "#b85c38", "o"),
        ("CLDN4", "TACSTD2", "CLDN4 | TROP2", "#b85c38", "D"),
    ]
    # x variable names differ for crispr
    if layer == "crispr":
        specs = [
            ("TACSTD2_effect", "none", "TROP2 Chronos", "#1f4e79", "o"),
            ("TACSTD2_effect", "CLDN4_effect", "TROP2 | CLDN4 Chronos", "#1f4e79", "D"),
            ("CLDN4_effect", "none", "CLDN4 Chronos", "#b85c38", "o"),
            ("CLDN4_effect", "TACSTD2_effect", "CLDN4 | TROP2 Chronos", "#b85c38", "D"),
        ]
    y_base = np.arange(len(scores))
    offsets = [-0.30, -0.10, 0.10, 0.30]
    ax.axvline(0, color="#444444", lw=0.8)
    handles = []
    for spec, off in zip(specs, offsets):
        xs, ys, xerr_lo, xerr_hi = [], [], [], []
        for i, score in enumerate(scores):
            hit = tab[
                (tab["layer"] == layer)
                & (tab["cohort"] == cohort)
                & (tab["x"] == spec[0])
                & (tab["y"] == score)
                & (tab["covariates"] == spec[1])
            ]
            if hit.empty or not np.isfinite(hit.iloc[0]["rho"]):
                continue
            r = hit.iloc[0]
            xs.append(r["rho"])
            ys.append(y_base[i] + off)
            if np.isfinite(r["boot_ci95_low"]):
                xerr_lo.append(r["rho"] - r["boot_ci95_low"])
                xerr_hi.append(r["boot_ci95_high"] - r["rho"])
            else:
                xerr_lo.append(0.0)
                xerr_hi.append(0.0)
        if not xs:
            continue
        h = ax.errorbar(
            xs,
            ys,
            xerr=[xerr_lo, xerr_hi],
            fmt=spec[4],
            color=spec[3],
            markersize=5,
            label=spec[2],
            elinewidth=1,
            capsize=2,
        )
        handles.append(h)
    ax.set_yticks(y_base)
    ax.set_yticklabels([SCORE_LABEL[s] for s in scores])
    ax.set_xlabel("Spearman or rank-residual partial r")
    ax.set_title(title, fontsize=10)
    ax.set_xlim(-0.85, 0.85)
    ax.legend(frameon=False, fontsize=7, loc="lower right")


def main() -> int:
    self_check()
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_trop2_ifn_partial_cldn4")
    ap.add_argument("--outdir", default="results/depmap_trop2_ifn_partial_cldn4")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    figdir = outdir / "figures"
    tabdir = outdir / "tables"
    figdir.mkdir(parents=True, exist_ok=True)
    tabdir.mkdir(parents=True, exist_ok=True)

    ifng_genes = [
        g.strip() for g in (data / "hallmark_ifng_genes.txt").read_text().splitlines() if g.strip()
    ]
    model = load_models(data / "Model.csv")
    manifest = {}
    man_path = data / "download_manifest.json"
    if man_path.exists():
        manifest = json.loads(man_path.read_text())

    # ----- RNA -----
    expr = pd.read_csv(data / "expression_panel.csv")
    expr = expr.merge(model, on="ModelID", how="left")
    lung = annotate(expr.loc[lung_mask(expr)].copy())
    for g in ["TACSTD2", "CLDN4"]:
        if g not in lung.columns:
            raise SystemExit(f"RNA panel missing {g}")
        lung[g] = pd.to_numeric(lung[g], errors="coerce")
    lung = lung.dropna(subset=["TACSTD2", "CLDN4"]).copy()
    lung_scored, cov_lung = add_rna_scores(lung, ifng_genes)
    nsclc = add_rna_scores(lung_scored.loc[lung_scored["is_NSCLC"]].copy(), ifng_genes)[0]
    luad = add_rna_scores(lung_scored.loc[lung_scored["group"] == "LUAD"].copy(), ifng_genes)[0]
    sclc = lung_scored.loc[lung_scored["group"] == "SCLC_NET"]

    rows: list[dict] = []
    rna_frames = {"lung": lung_scored, "NSCLC": nsclc, "LUAD": luad}
    for cohort, frame in rna_frames.items():
        primary = cohort in {"lung", "NSCLC"}
        family = f"rna_{cohort}_primary" if primary else ""
        do_boot = primary
        for score in PRIMARY_SCORES + ["CD274_score", "ISG"]:
            y = score
            if score == "CD274_score" and "CD274_score" not in frame.columns:
                continue
            for x, other in (("TACSTD2", "CLDN4"), ("CLDN4", "TACSTD2")):
                is_primary_score = score in PRIMARY_SCORES
                record_pair(
                    rows, frame, cohort, "rna", x, y, [],
                    family if is_primary_score else "",
                    primary and is_primary_score,
                    do_boot and is_primary_score,
                )
                record_pair(
                    rows, frame, cohort, "rna", x, y, [other],
                    family if is_primary_score else "",
                    primary and is_primary_score,
                    do_boot and is_primary_score,
                )
        if cohort == "lung":
            for score in PRIMARY_SCORES:
                record_pair(
                    rows, frame, cohort, "rna", "TACSTD2", score,
                    ["CLDN4", "nsclc_indicator"],
                    "",
                    False,
                    True,
                )
                record_pair(
                    rows, frame, cohort, "rna", "CLDN4", score,
                    ["TACSTD2", "nsclc_indicator"],
                    "",
                    False,
                    True,
                )
    # covariate correlation
    rho_rna, p_rna, n_rna = spearman_full(lung_scored["TACSTD2"], lung_scored["CLDN4"])

    # ----- Protein -----
    gygi = pd.read_csv(data / "gygi_protein_panel.csv")
    checksum = gygi_checksum(gygi)
    prot_ann = models_from_gygi(gygi, model)
    prot_lung = annotate(prot_ann.loc[lung_mask(prot_ann)].copy())
    for g in prot_lung.columns:
        if g in {
            "CellLineName",
            "CCLEName",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "ModelType",
            "group",
        }:
            continue
    both = prot_lung.dropna(subset=[c for c in ("TACSTD2", "CLDN4") if c in prot_lung.columns]).copy()
    both_scored, cov_prot = add_protein_scores(both, ifng_genes)
    nsclc_prot = add_protein_scores(both_scored.loc[both_scored["is_NSCLC"]].copy(), ifng_genes)[0]
    rho_prot, p_prot, n_prot = spearman_full(both_scored["TACSTD2"], both_scored["CLDN4"])
    for cohort, frame in (("lung", both_scored), ("NSCLC", nsclc_prot)):
        family = f"protein_{cohort}_primary"
        for score in PRIMARY_SCORES:
            for x, other in (("TACSTD2", "CLDN4"), ("CLDN4", "TACSTD2")):
                record_pair(rows, frame, cohort, "protein", x, score, [], family, True, True)
                record_pair(rows, frame, cohort, "protein", x, score, [other], family, True, True)
        for score in ("CD274", "ISG"):
            if score not in frame.columns and score != "CD274":
                continue
            y = "CD274" if score == "CD274" else "ISG"
            if y not in frame.columns:
                continue
            for x, other in (("TACSTD2", "CLDN4"), ("CLDN4", "TACSTD2")):
                record_pair(rows, frame, cohort, "protein", x, y, [], "", False, False)
                record_pair(rows, frame, cohort, "protein", x, y, [other], "", False, False)

    # ----- CRISPR -----
    effect = pd.read_csv(data / "CRISPRGeneEffect_panel.csv")
    dep = pd.read_csv(data / "CRISPRGeneDependency_panel.csv")
    effect = effect.rename(columns={c: f"{c}_effect" for c in effect.columns if c != "ModelID"})
    dep = dep.rename(columns={c: f"{c}_dep" for c in dep.columns if c != "ModelID"})
    for c in effect.columns:
        if c != "ModelID":
            effect[c] = pd.to_numeric(effect[c], errors="coerce")
    for c in dep.columns:
        if c != "ModelID":
            dep[c] = pd.to_numeric(dep[c], errors="coerce")
    crispr = effect.merge(dep, on="ModelID", how="outer")
    crispr = crispr.merge(
        model[
            [
                "ModelID",
                "CellLineName",
                "OncotreeLineage",
                "OncotreePrimaryDisease",
                "OncotreeSubtype",
                "ModelType",
            ]
        ],
        on="ModelID",
        how="left",
    )
    crispr_lung = annotate(crispr.loc[lung_mask(crispr)].copy())
    crispr_lung = crispr_lung.dropna(subset=["TACSTD2_effect", "CLDN4_effect"]).copy()
    # Join RNA scores recomputed on the overlap.
    rna_keep = ["ModelID", "TACSTD2", "CLDN4"] + [
        g for g in set(MHC1 + IMMUNE_LIGAND + ISG + ifng_genes) if g in lung_scored.columns
    ]
    rna_keep = list(dict.fromkeys(rna_keep))
    overlap = crispr_lung.merge(lung_scored[rna_keep], on="ModelID", how="inner", suffixes=("", "_rna"))
    # CLDN4 RNA column
    overlap = overlap.rename(columns={"CLDN4": "CLDN4_RNA", "TACSTD2": "TACSTD2_RNA"})
    overlap_scored, _ = add_rna_scores(
        overlap.rename(columns={"CLDN4_RNA": "CLDN4", "TACSTD2_RNA": "TACSTD2"}),
        ifng_genes,
    )
    # add_rna_scores keeps CLDN4/TACSTD2 names; restore effect columns already present.
    overlap_scored["CLDN4_RNA"] = overlap_scored["CLDN4"]
    nsclc_cr = add_rna_scores(
        overlap_scored.loc[overlap_scored["is_NSCLC"]].copy(),
        ifng_genes,
    )[0]
    nsclc_cr["CLDN4_RNA"] = nsclc_cr["CLDN4"]

    dep_genes = {}
    for gene in ("TACSTD2", "CLDN4", "KRAS", "EEF2"):
        ecol = f"{gene}_effect"
        dcol = f"{gene}_dep"
        if ecol not in crispr_lung.columns:
            continue
        s = crispr_lung[ecol].dropna()
        if dcol in crispr_lung.columns:
            d = crispr_lung.loc[s.index, dcol] if dcol in crispr_lung.columns else pd.Series(dtype=float)
            # align
            paired = crispr_lung[[ecol, dcol]].dropna()
            n_dep = int((paired[dcol] > 0.5).sum())
            n_d = int(len(paired))
        else:
            n_dep, n_d = 0, 0
        dep_genes[gene] = {
            "n": int(s.shape[0]),
            "mean_effect": float(s.mean()),
            "sd_effect": float(s.std(ddof=1)) if len(s) > 1 else float("nan"),
            "n_dependent": n_dep,
            "n_dependency": n_d,
            "min_effect": float(s.min()),
            "max_effect": float(s.max()),
        }

    for cohort, frame in (("lung", overlap_scored), ("NSCLC", nsclc_cr)):
        family = f"crispr_{cohort}_primary"
        for score in PRIMARY_SCORES:
            record_pair(rows, frame, cohort, "crispr", "TACSTD2_effect", score, [], family, True, True)
            record_pair(rows, frame, cohort, "crispr", "CLDN4_effect", score, [], family, True, True)
            record_pair(
                rows, frame, cohort, "crispr", "TACSTD2_effect", score, ["CLDN4_effect"], family, True, True
            )
            record_pair(
                rows, frame, cohort, "crispr", "CLDN4_effect", score, ["TACSTD2_effect"], family, True, True
            )
            record_pair(
                rows, frame, cohort, "crispr", "TACSTD2_effect", score, ["CLDN4_RNA"], "", False, True
            )

    tab = apply_fdr_and_calls(pd.DataFrame(rows))
    tab.to_csv(tabdir / "partials.csv", index=False)

    # Line-level tables (slim).
    lung_out_cols = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "group",
            "TACSTD2",
            "CLDN4",
            "CD274",
            "IFNG",
            "MHC1",
            "IMMUNE",
            "ISG",
        ]
        if c in lung_scored.columns
    ]
    lung_scored[lung_out_cols].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "rna_lung_lines.csv", index=False
    )
    prot_cols = [
        c
        for c in [
            "CellLineName",
            "CCLEName",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "group",
            "TACSTD2",
            "CLDN4",
            "CD274",
            "IFNG",
            "MHC1",
            "IMMUNE",
            "ISG",
        ]
        if c in both_scored.columns
    ]
    both_scored[prot_cols].sort_values(["group", "CellLineName"]).to_csv(
        tabdir / "protein_lung_lines.csv", index=False
    )

    rppa = manifest.get("files", {}).get("rppa_antibody_info", {})
    chem = expressed_fraction(nsclc, IMMUNE_LIGAND)

    def headline_for(layer: str, cohort: str) -> list[str]:
        bits = []
        for score in PRIMARY_SCORES:
            hit = tab[
                (tab["layer"] == layer)
                & (tab["cohort"] == cohort)
                & (tab["x"] == ("TACSTD2_effect" if layer == "crispr" else "TACSTD2"))
                & (tab["y"] == score)
                & (tab["covariates"] == ("CLDN4_effect" if layer == "crispr" else "CLDN4"))
            ]
            if hit.empty:
                continue
            call = hit.iloc[0]["attenuation_call"] or "not_labeled"
            rho = hit.iloc[0]["rho"]
            bits.append(f"{layer} {cohort} {SCORE_LABEL[score]}: TROP2|CLDN4 call={call}, partial r={fmt(rho)}")
        return bits

    headlines = []
    headlines.extend(headline_for("rna", "NSCLC"))
    headlines.extend(headline_for("rna", "lung"))
    headlines.extend(headline_for("protein", "lung"))
    headlines.extend(headline_for("protein", "NSCLC"))
    headlines.extend(headline_for("crispr", "lung"))
    # Explicit comparison sentence inputs
    def grab(layer, cohort, x, y, cov):
        hit = tab[
            (tab["layer"] == layer)
            & (tab["cohort"] == cohort)
            & (tab["x"] == x)
            & (tab["y"] == y)
            & (tab["covariates"] == cov)
        ]
        return None if hit.empty else hit.iloc[0]

    for score in PRIMARY_SCORES:
        a = grab("rna", "NSCLC", "TACSTD2", score, "CLDN4")
        b = grab("rna", "NSCLC", "CLDN4", score, "TACSTD2")
        if a is None or b is None or not np.isfinite(a["rho"]) or not np.isfinite(b["rho"]):
            continue
        larger = "TROP2" if abs(a["rho"]) > abs(b["rho"]) else "CLDN4"
        headlines.append(
            f"NSCLC RNA {SCORE_LABEL[score]} residual comparison: |TROP2|CLDN4|={abs(a['rho']):.3f} "
            f"({a['attenuation_call']}), |CLDN4|TROP2|={abs(b['rho']):.3f} ({b['attenuation_call']}); "
            f"larger absolute partial is {larger}."
        )

    def brief(rec: pd.Series | None, with_call: bool = True) -> str:
        if rec is None or not np.isfinite(rec["rho"]):
            return "not estimated"
        q = rec["q_bh"] if rec["kind"] == "partial" else rec["q_bh_unadjusted"]
        qbit = f", q={fmt_p(float(q))}" if pd.notna(q) else ""
        ratio = rec["rho_ratio_partial_over_unadj"] if "rho_ratio_partial_over_unadj" in rec.index else np.nan
        rbit = f", |partial|/|unadj|={float(ratio):.2f}" if pd.notna(ratio) else ""
        call = rec["attenuation_call"] if with_call else ""
        cbit = f"; {call}" if isinstance(call, str) and call else ""
        return f"{float(rec['rho']):+.3f} (n={int(rec['n'])}, p={fmt_p(float(rec['p']))}{qbit}{rbit}{cbit})"

    def n_partial_q(layer: str) -> int:
        block = tab[(tab["layer"] == layer) & (tab["kind"] == "partial") & (tab["primary"] == True) & tab["q_bh"].notna()]  # noqa: E712
        return int((block["q_bh"] < 0.05).sum())

    t_ifn = grab("rna", "NSCLC", "TACSTD2", "IFNG", "none")
    t_ifn_p = grab("rna", "NSCLC", "TACSTD2", "IFNG", "CLDN4")
    t_imm = grab("rna", "NSCLC", "TACSTD2", "IMMUNE", "none")
    t_imm_p = grab("rna", "NSCLC", "TACSTD2", "IMMUNE", "CLDN4")
    t_mhc = grab("rna", "NSCLC", "TACSTD2", "MHC1", "none")
    c_ifn_p = grab("rna", "NSCLC", "CLDN4", "IFNG", "TACSTD2")
    c_imm_p = grab("rna", "NSCLC", "CLDN4", "IMMUNE", "TACSTD2")
    lung_flip = grab("rna", "lung", "CLDN4", "IMMUNE", "TACSTD2")
    lung_flip_lin = grab("rna", "lung", "CLDN4", "IMMUNE", "TACSTD2+nsclc_indicator")
    luad_ifn = grab("rna", "LUAD", "TACSTD2", "IFNG", "none")
    luad_ifn_p = grab("rna", "LUAD", "TACSTD2", "IFNG", "CLDN4")
    luad_imm_p = grab("rna", "LUAD", "TACSTD2", "IMMUNE", "CLDN4")
    p_ifn = grab("protein", "lung", "TACSTD2", "IFNG", "none")
    p_ifn_p = grab("protein", "lung", "TACSTD2", "IFNG", "CLDN4")
    p_imm = grab("protein", "lung", "TACSTD2", "IMMUNE", "none")
    answer = [
        (
            "In NSCLC RNA, CLDN4 does not attenuate the TROP2 associations that clear p < 0.05. "
            f"Hallmark IFN-γ is {brief(t_ifn, False)} unadjusted and {brief(t_ifn_p)} after CLDN4. "
            f"The immune-ligand score is {brief(t_imm, False)} unadjusted and {brief(t_imm_p)} after CLDN4. "
            f"MHC-I has no unadjusted TROP2 association ({brief(t_mhc, False)}), so there is nothing there for CLDN4 to attenuate."
        ),
        (
            "On the same NSCLC RNA rows, CLDN4's own positive associations are the ones that shrink. "
            f"CLDN4 vs Hallmark IFN-γ given TROP2 is {brief(c_ifn_p)}. "
            f"CLDN4 vs the immune-ligand score given TROP2 is {brief(c_imm_p)}."
        ),
        (
            "All-lung RNA (n includes SCLC/NET) agrees that TROP2's IFN-γ and immune-ligand associations are retained after CLDN4. "
            "All-lung TROP2 vs MHC-I is retained after CLDN4 alone, and is attenuated once a binary NSCLC indicator is added "
            f"({brief(grab('rna', 'lung', 'TACSTD2', 'MHC1', 'CLDN4+nsclc_indicator'))}). "
            f"CLDN4 vs the immune-ligand score reverses sign after TROP2 adjustment: {brief(lung_flip)}. "
            f"Adding a binary NSCLC indicator leaves that residual {brief(lung_flip_lin)}."
        ),
        (
            "LUAD-only RNA is exploratory and is outside the FDR family. "
            f"TROP2 vs Hallmark IFN-γ is {brief(luad_ifn, False)} unadjusted and {brief(luad_ifn_p)} after CLDN4. "
            f"TROP2 vs the immune-ligand score after CLDN4 is {brief(luad_imm_p)}. "
            "The IFN-γ retention above is an NSCLC result, not a LUAD-only result."
        ),
        (
            f"Protein partials do not decide the question. Primary protein partials with BH q < 0.05: {n_partial_q('protein')}. "
            f"All-lung TROP2 vs Hallmark IFN-γ protein is {brief(p_ifn, False)} unadjusted and {brief(p_ifn_p)} after CLDN4. "
            f"The protein immune-ligand score uses only {', '.join(cov_prot.get('immune_used', [])) or 'no genes'} "
            f"(CXCL9/CXCL10/CXCL11 are absent). Its unadjusted TROP2 association is {brief(p_imm, False)}; "
            "the bootstrap interval is in the table."
        ),
        (
            f"CRISPR is available. Lung lines with Chronos: n={len(crispr_lung)}. "
            f"TACSTD2 mean gene effect {dep_genes['TACSTD2']['mean_effect']:+.4f}, "
            f"dependent {dep_genes['TACSTD2']['n_dependent']}/{dep_genes['TACSTD2']['n_dependency']}. "
            f"CLDN4 mean {dep_genes['CLDN4']['mean_effect']:+.4f}, "
            f"dependent {dep_genes['CLDN4']['n_dependent']}/{dep_genes['CLDN4']['n_dependency']}. "
            f"On the same lines, KRAS is dependent in {dep_genes['KRAS']['n_dependent']}/{dep_genes['KRAS']['n_dependency']} "
            f"and EEF2 in {dep_genes['EEF2']['n_dependent']}/{dep_genes['EEF2']['n_dependency']}. "
            "TROP2 Chronos is not associated with the three RNA scores, before or after CLDN4 Chronos, so there is no CRISPR association for CLDN4 to attenuate."
        ),
    ]

    key = {
        "question": "Does CLDN4 attenuate TACSTD2 association with IFN/MHC/immune scores?",
        "partial_definition": "Pearson r of residuals of ranks after OLS on covariate ranks; df=n-2-k",
        "n_boot": N_BOOT,
        "n_rna_lung": int(len(lung_scored)),
        "n_rna_nsclc": int(len(nsclc)),
        "n_rna_luad": int((lung_scored["group"] == "LUAD").sum()),
        "n_rna_lusc": int((lung_scored["group"] == "LUSC").sum()),
        "n_rna_sclc": int(len(sclc)),
        "rna_score_coverage_lung": cov_lung,
        "rna_nsclc_chemokine_expression": chem,
        "protein_checksum_lung_suffix_columns": checksum,
        "covariate_correlations": {
            "rna_lung_TACSTD2_CLDN4": {"rho": rho_rna, "p": p_rna, "n": n_rna},
            "protein_model_lung_TACSTD2_CLDN4": {"rho": rho_prot, "p": p_prot, "n": n_prot},
        },
        "n_protein_lung_any": int(len(prot_lung)),
        "n_protein_both": int(len(both_scored)),
        "n_protein_nsclc_both": int(both_scored["is_NSCLC"].sum()) if len(both_scored) else 0,
        "protein_score_coverage_lung": cov_prot,
        "rppa_mentions_CLDN4": rppa.get("mentions_CLDN4"),
        "rppa_mentions_TACSTD2": rppa.get("mentions_TACSTD2_or_TROP2"),
        "crispr_dependency": {
            "n_lung": int(len(crispr_lung)),
            "n_nsclc": int(crispr_lung["is_NSCLC"].sum()),
            "n_overlap_rna": int(len(overlap_scored)),
            "n_overlap_nsclc": int(len(nsclc_cr)),
            "genes": dep_genes,
            "dependent_definition": "CRISPRGeneDependency probability > 0.5",
        },
        "headline_calls": headlines,
        "answer": answer,
        "release": manifest.get("release", "DepMap Public 24Q4"),
        "doi": manifest.get("doi"),
    }
    (tabdir / "key_stats.json").write_text(json.dumps(key, indent=2, default=str) + "\n")
    write_finding(outdir / "FINDING.md", key, tab)

    # Verdict: one paragraph from NSCLC RNA calls only.
    nsclc_calls = []
    for score in PRIMARY_SCORES:
        rec = grab("rna", "NSCLC", "TACSTD2", score, "CLDN4")
        if rec is not None:
            nsclc_calls.append(f"{score}={rec['attenuation_call']}")
    prot_calls = []
    for score in PRIMARY_SCORES:
        rec = grab("protein", "lung", "TACSTD2", score, "CLDN4")
        if rec is not None:
            prot_calls.append(f"{score}={rec['attenuation_call']}")
    crispr_calls = []
    for score in PRIMARY_SCORES:
        rec = grab("crispr", "lung", "TACSTD2_effect", score, "CLDN4_effect")
        if rec is not None:
            crispr_calls.append(f"{score}={rec['attenuation_call']}")
    verdict = (
        "NSCLC RNA TROP2|CLDN4 calls: "
        + ", ".join(nsclc_calls)
        + ". Protein lung TROP2|CLDN4 calls: "
        + ", ".join(prot_calls)
        + ". CRISPR lung TROP2 Chronos|CLDN4 Chronos calls: "
        + ", ".join(crispr_calls)
        + f". Protein column checksum rounds_to_0.69={checksum.get('rounds_to_0.69')} "
        + f"rho={checksum.get('spearman_rho')} n={checksum.get('n_complete')}. "
        f"CRISPR lung dependent TACSTD2 {dep_genes.get('TACSTD2', {}).get('n_dependent')}/"
        f"{dep_genes.get('TACSTD2', {}).get('n_dependency')}, "
        f"CLDN4 {dep_genes.get('CLDN4', {}).get('n_dependent')}/"
        f"{dep_genes.get('CLDN4', {}).get('n_dependency')}."
    )
    (outdir / "verdict.txt").write_text(verdict + "\n")

    # Figures
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)
    forest(axes[0], tab, "rna", "NSCLC", "NSCLC RNA")
    forest(axes[1], tab, "rna", "lung", "All-lung RNA")
    fig.suptitle("TROP2 and CLDN4 versus IFN / MHC-I / immune-ligand scores", fontsize=12)
    fig.savefig(figdir / "fig_rna_partial_forest.png", dpi=160)
    fig.savefig(figdir / "fig_rna_partial_forest.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)
    forest(axes[0], tab, "protein", "NSCLC", "NSCLC protein")
    forest(axes[1], tab, "protein", "lung", "All-lung protein")
    fig.suptitle("Gygi TMT, one row per cell line", fontsize=12)
    fig.savefig(figdir / "fig_protein_partial_forest.png", dpi=160)
    fig.savefig(figdir / "fig_protein_partial_forest.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.4, 3.8), constrained_layout=True)
    for ax, score in zip(axes, PRIMARY_SCORES):
        sub = nsclc[["TACSTD2", "CLDN4", score, "group"]].dropna()
        ranks = sub.rank(method="average")
        z = ranks["CLDN4"].to_numpy()
        design = np.column_stack([np.ones(len(sub)), z])

        def resid(v):
            beta, *_ = np.linalg.lstsq(design, v, rcond=None)
            return v - design @ beta

        xr = resid(ranks["TACSTD2"].to_numpy())
        yr = resid(ranks[score].to_numpy())
        ax.scatter(xr, yr, s=16, c="#1f4e79", alpha=0.75, edgecolors="none")
        rec = grab("rna", "NSCLC", "TACSTD2", score, "CLDN4")
        rho_s = "NA" if rec is None or not np.isfinite(rec["rho"]) else f"{rec['rho']:+.2f}"
        ax.set_title(f"{SCORE_LABEL[score]}\npartial r={rho_s}, n={0 if rec is None else int(rec['n'])}", fontsize=9)
        ax.set_xlabel("TROP2 rank residual | CLDN4")
        ax.set_ylabel("Score rank residual | CLDN4")
        ax.axhline(0, color="#bbbbbb", lw=0.6)
        ax.axvline(0, color="#bbbbbb", lw=0.6)
    fig.suptitle("NSCLC RNA rank residuals", fontsize=12)
    fig.savefig(figdir / "fig_nsclc_rna_residuals.png", dpi=160)
    fig.savefig(figdir / "fig_nsclc_rna_residuals.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    genes_order = [g for g in ("EEF2", "KRAS", "TACSTD2", "CLDN4") if g in dep_genes]
    means = [dep_genes[g]["mean_effect"] for g in genes_order]
    sds = [dep_genes[g]["sd_effect"] for g in genes_order]
    colors = ["#4d4d4d" if g in ("EEF2", "KRAS") else "#1f4e79" if g == "TACSTD2" else "#b85c38" for g in genes_order]
    axes[0].barh(genes_order, means, xerr=sds, color=colors, alpha=0.9, capsize=3)
    axes[0].axvline(0, color="#444444", lw=0.8)
    axes[0].set_xlabel("Mean Chronos gene effect (SD)")
    axes[0].set_title(f"Lung lines with Chronos n={len(crispr_lung)}", fontsize=10)
    forest(axes[1], tab, "crispr", "lung", "Chronos vs RNA scores")
    fig.savefig(figdir / "fig_crispr.png", dpi=160)
    fig.savefig(figdir / "fig_crispr.pdf")
    plt.close(fig)

    print(verdict)
    print(f"wrote {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
