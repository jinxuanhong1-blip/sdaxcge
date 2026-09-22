#!/usr/bin/env python3
"""PPT funnel: DepMap/CCLE TROP2–CLDN4/TJ coexpression + immune scores in lung lines.

Public DepMap 24Q4 RNA and Gygi/Nusinow CCLE TMT protein only.
Numbers in FINDING.md are written from tables/key_stats.json — not typed by hand.

Bridge question for the resistance/TJ PPT slide:
  1) Does TROP2 coexpress with CLDN4 protein and locked TJ gene scores in lung lines?
  2) Does TROP2 coexpress with immune-related cancer-cell scores in the same lines?
  3) Do TJ scores attenuate TROP2–immune associations (or the reverse)?

Cultured lines: no T-cell infiltrate, no IFN treatment, no ICI outcome.
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

TJ_EPITHELIAL = [
    "CLDN1",
    "CLDN3",
    "CLDN4",
    "CLDN7",
    "OCLN",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "TJP3",
    "F11R",
    "JAM2",
    "JAM3",
    "CGN",
    "CGNL1",
    "CRB3",
    "ILDR1",
    "LSR",
]
TJ_TISMO = ["CLDN3", "CLDN4", "CLDN6", "CLDN7", "CDH1", "F11R", "OCLN"]
CLDN4_TJ_EDGE = [
    "CLDN4",
    "CLDN1",
    "CLDN7",
    "CGNL1",
    "MARVELD2",
    "MARVELD3",
    "TJP1",
    "TJP2",
    "ILDR1",
    "CLDN3",
    "OCLN",
]
# CLDN4 held out of epithelial and edge scores for partials that ask about
# residual TJ beyond CLDN4 itself.
TJ_EPITHELIAL_NO_CLDN4 = [g for g in TJ_EPITHELIAL if g != "CLDN4"]
CLDN4_TJ_EDGE_NO_CLDN4 = [g for g in CLDN4_TJ_EDGE if g != "CLDN4"]
TJ_TISMO_NO_CLDN4 = [g for g in TJ_TISMO if g != "CLDN4"]

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

PRIMARY_TJ = ["CLDN4", "TJ_EPITHELIAL", "TJ_TISMO", "CLDN4_TJ_EDGE"]
PRIMARY_IMMUNE = ["IFNG", "MHC1", "IMMUNE"]
SCORE_LABEL = {
    "CLDN4": "CLDN4",
    "TJ_EPITHELIAL": "TJ epithelial (18)",
    "TJ_TISMO": "TJ TISMO (7)",
    "CLDN4_TJ_EDGE": "CLDN4 TJ edge (11)",
    "TJ_EPITHELIAL_noCLDN4": "TJ epithelial without CLDN4",
    "TJ_TISMO_noCLDN4": "TJ TISMO without CLDN4",
    "CLDN4_TJ_EDGE_noCLDN4": "CLDN4 TJ edge without CLDN4",
    "IFNG": "Hallmark IFN-γ",
    "MHC1": "MHC-I (HLA-A/B/C+B2M)",
    "IMMUNE": "Immune-ligand score",
    "ISG": "Compact ISG",
    "CD274": "CD274",
    "EPCAM": "EPCAM",
}
N_BOOT = 2000
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
        raise SystemExit(f"Model.csv missing {missing}")
    return model


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["group"] = out.apply(cohort_tag, axis=1)
    out["is_NSCLC"] = out["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"
    out["nsclc_indicator"] = out["is_NSCLC"].astype(float)
    return out


def gygi_checksum(panel: pd.DataFrame) -> dict:
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
    }


def models_from_gygi(panel: pd.DataFrame, model: pd.DataFrame) -> pd.DataFrame:
    meta = {"Gene_Symbol", "Protein_Id"}
    num = panel.drop(columns=[c for c in panel.columns if c in meta], errors="ignore")
    tmp = pd.concat(
        [
            panel[["Gene_Symbol"]].reset_index(drop=True),
            num.apply(pd.to_numeric, errors="coerce").reset_index(drop=True),
        ],
        axis=1,
    )
    prot = tmp.groupby("Gene_Symbol").mean(numeric_only=True)
    ccle = {str(v): i for i, v in model.set_index("ModelID")["CCLEName"].dropna().items()}
    stripped = {
        str(v).upper(): i for i, v in model.set_index("ModelID")["StrippedCellLineName"].dropna().items()
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


def add_scores(frame: pd.DataFrame, ifng_genes: list[str], layer: str) -> tuple[pd.DataFrame, dict]:
    out = frame.copy()
    cov: dict = {"layer": layer}

    def add_sig(name: str, genes: list[str], min_members: int, min_present: int | None = None) -> None:
        present = genes
        if min_present is not None:
            present = [
                g
                for g in genes
                if g in out.columns and int(pd.to_numeric(out[g], errors="coerce").notna().sum()) >= min_present
            ]
        score, used = signature(out, present, min_members=min_members)
        out[name] = score
        cov[name] = {"used": used, "n_requested": len(genes), "n_used": len(used)}

    if layer == "rna":
        add_sig("TJ_EPITHELIAL", TJ_EPITHELIAL, min_members=max(5, len(TJ_EPITHELIAL) // 2))
        add_sig("TJ_TISMO", TJ_TISMO, min_members=max(3, len(TJ_TISMO) // 2))
        add_sig("CLDN4_TJ_EDGE", CLDN4_TJ_EDGE, min_members=max(4, len(CLDN4_TJ_EDGE) // 2))
        add_sig("TJ_EPITHELIAL_noCLDN4", TJ_EPITHELIAL_NO_CLDN4, min_members=max(5, len(TJ_EPITHELIAL_NO_CLDN4) // 2))
        add_sig("TJ_TISMO_noCLDN4", TJ_TISMO_NO_CLDN4, min_members=max(3, len(TJ_TISMO_NO_CLDN4) // 2))
        add_sig("CLDN4_TJ_EDGE_noCLDN4", CLDN4_TJ_EDGE_NO_CLDN4, min_members=max(4, len(CLDN4_TJ_EDGE_NO_CLDN4) // 2))
        add_sig("IFNG", ifng_genes, min_members=max(5, int(0.5 * len(ifng_genes))))
        add_sig("MHC1", MHC1, min_members=4)
        add_sig("IMMUNE", [g for g in IMMUNE_LIGAND if g in out.columns], min_members=3)
        add_sig("ISG", [g for g in ISG if g in out.columns], min_members=5)
    else:
        # Protein: require ≥8 non-missing lines before allowing a gene into a score.
        add_sig("TJ_EPITHELIAL", TJ_EPITHELIAL, min_members=3, min_present=8)
        add_sig("TJ_TISMO", TJ_TISMO, min_members=2, min_present=8)
        add_sig("CLDN4_TJ_EDGE", CLDN4_TJ_EDGE, min_members=3, min_present=8)
        add_sig("TJ_EPITHELIAL_noCLDN4", TJ_EPITHELIAL_NO_CLDN4, min_members=3, min_present=8)
        add_sig("TJ_TISMO_noCLDN4", TJ_TISMO_NO_CLDN4, min_members=2, min_present=8)
        add_sig("CLDN4_TJ_EDGE_noCLDN4", CLDN4_TJ_EDGE_NO_CLDN4, min_members=3, min_present=8)
        hm = [g for g in ifng_genes if g in out.columns and int(out[g].notna().sum()) >= 8]
        add_sig("IFNG", hm, min_members=5)
        add_sig("MHC1", [g for g in MHC1 if g in out.columns], min_members=2)
        imm = [g for g in IMMUNE_LIGAND if g in out.columns and int(out[g].notna().sum()) >= 8]
        if "CD274" in out.columns and "CD274" not in imm and int(out["CD274"].notna().sum()) >= 8:
            imm.append("CD274")
        add_sig("IMMUNE", imm, min_members=2 if len(imm) >= 2 else 1)
        isg = [g for g in ISG if g in out.columns and int(out[g].notna().sum()) >= 8]
        add_sig("ISG", isg, min_members=5)
    return out, cov


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
        label = f"{layer}|{cohort}|{x}|{y}|{cov_label}|{kind}"
        seed = seed_for(label)
        mat = sub.to_numpy()

        def stat(sample: np.ndarray) -> float:
            sdf = pd.DataFrame(sample, columns=cols)
            if covariates:
                r, _, _ = partial_corr(sdf, x, y, covariates)
                return r
            return spearman_r(sample[:, 0], sample[:, 1])

        ci_b_lo, ci_b_hi = bootstrap_stat(mat, stat, N_BOOT, seed)
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
            "note": "",
        }
    )


def apply_fdr_and_calls(tab: pd.DataFrame) -> pd.DataFrame:
    out = tab.copy()
    out["q_bh"] = np.nan
    out["q_bh_unadjusted"] = np.nan
    for fam, idx in out.groupby("family").groups.items():
        if not fam:
            continue
        block = out.loc[idx]
        for kind, qcol in (("partial", "q_bh"), ("spearman", "q_bh_unadjusted")):
            ok = block["p"].notna() & block["primary"].astype(bool) & (block["kind"] == kind)
            use = block.index[ok.to_numpy()]
            if len(use) == 0:
                continue
            out.loc[use, qcol] = bh_fdr(out.loc[use, "p"].astype(float).tolist())

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
    return f"{rec['rho']:+.3f} (n={int(rec['n'])}, p={fmt_p(rec['p'])}{qbit}){ci}{rbit}{callbit}"


def pick(tab: pd.DataFrame, layer: str, cohort: str, x: str, y: str, covariates: str = "none") -> dict | None:
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


def make_figures(outdir: Path, rna: pd.DataFrame, prot: pd.DataFrame, tab: pd.DataFrame) -> None:
    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)

    # Figure 1: RNA NSCLC TROP2 vs TJ and immune scores
    nsclc = rna[rna["is_NSCLC"]].copy()
    scores = ["CLDN4", "TJ_EPITHELIAL", "TJ_TISMO", "IFNG", "MHC1", "IMMUNE"]
    fig, axes = plt.subplots(2, 3, figsize=(11, 7), constrained_layout=True)
    for ax, y in zip(axes.ravel(), scores):
        sub = nsclc[["TACSTD2", y]].dropna()
        ax.scatter(sub["TACSTD2"], sub[y], s=14, alpha=0.7, c="#1f4e79", edgecolors="none")
        rho, p, n = spearman_full(sub["TACSTD2"], sub[y])
        ax.set_title(f"{SCORE_LABEL.get(y, y)}\nρ={rho:+.2f}, n={n}, p={fmt_p(p)}", fontsize=9)
        ax.set_xlabel("TACSTD2 log2(TPM+1)")
        ax.set_ylabel(SCORE_LABEL.get(y, y))
    fig.suptitle("DepMap 24Q4 NSCLC RNA: TROP2 vs TJ and immune scores", fontsize=12)
    fig.savefig(figdir / "fig_rna_nsclc_scatter.png", dpi=160)
    fig.savefig(figdir / "fig_rna_nsclc_scatter.pdf")
    plt.close(fig)

    # Figure 2: forest of primary RNA associations
    rows = []
    for y in PRIMARY_TJ + PRIMARY_IMMUNE:
        rec = pick(tab, "rna", "NSCLC", "TACSTD2", y)
        if rec and np.isfinite(rec.get("rho", np.nan)):
            rows.append(rec)
    if rows:
        fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
        ys = np.arange(len(rows))
        rhos = [r["rho"] for r in rows]
        lo = [r.get("boot_ci95_low", r.get("fisher_ci95_low")) for r in rows]
        hi = [r.get("boot_ci95_high", r.get("fisher_ci95_high")) for r in rows]
        labels = [SCORE_LABEL.get(r["y"], r["y"]) for r in rows]
        ax.errorbar(rhos, ys, xerr=[np.array(rhos) - np.array(lo), np.array(hi) - np.array(rhos)], fmt="o", color="#1f4e79")
        ax.axvline(0, color="0.5", lw=1)
        ax.set_yticks(ys)
        ax.set_yticklabels(labels)
        ax.set_xlabel("Spearman ρ (bootstrap 95% CI)")
        ax.set_title("NSCLC RNA: TROP2 vs TJ / immune scores")
        fig.savefig(figdir / "fig_rna_nsclc_forest.png", dpi=160)
        fig.savefig(figdir / "fig_rna_nsclc_forest.pdf")
        plt.close(fig)

    # Figure 3: protein TROP2–CLDN4 + available TJ partners
    if "TACSTD2" in prot.columns and "CLDN4" in prot.columns:
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), constrained_layout=True)
        partners = ["CLDN4", "CLDN7", "CLDN1"]
        lung = prot[(prot["OncotreeLineage"] == "Lung") & (prot["ModelType"] == "Cell Line")].copy()
        for ax, g in zip(axes, partners):
            if g not in lung.columns:
                ax.set_visible(False)
                continue
            sub = lung[["TACSTD2", g]].dropna()
            ax.scatter(sub["TACSTD2"], sub[g], s=22, alpha=0.8, c="#8b4513", edgecolors="none")
            rho, p, n = spearman_full(sub["TACSTD2"], sub[g])
            ax.set_title(f"{g}\nρ={rho:+.2f}, n={n}, p={fmt_p(p)}", fontsize=9)
            ax.set_xlabel("TACSTD2 protein")
            ax.set_ylabel(f"{g} protein")
        fig.suptitle("Gygi CCLE lung protein: TROP2 vs claudins", fontsize=12)
        fig.savefig(figdir / "fig_protein_claudin_scatter.png", dpi=160)
        fig.savefig(figdir / "fig_protein_claudin_scatter.pdf")
        plt.close(fig)

    # Figure 4: bridge partials (RNA NSCLC)
    bridge_ys = ["IFNG", "IMMUNE", "MHC1"]
    fig, ax = plt.subplots(figsize=(8, 4.2), constrained_layout=True)
    labels = []
    unadj = []
    part_tj = []
    for y in bridge_ys:
        u = pick(tab, "rna", "NSCLC", "TACSTD2", y)
        p = pick(tab, "rna", "NSCLC", "TACSTD2", y, "TJ_EPITHELIAL")
        if u and p and np.isfinite(u.get("rho", np.nan)):
            labels.append(SCORE_LABEL[y])
            unadj.append(u["rho"])
            part_tj.append(p["rho"] if np.isfinite(p.get("rho", np.nan)) else np.nan)
    xpos = np.arange(len(labels))
    w = 0.35
    ax.bar(xpos - w / 2, unadj, w, label="unadjusted", color="#1f4e79")
    ax.bar(xpos + w / 2, part_tj, w, label="| TJ epithelial", color="#c4a35a")
    ax.axhline(0, color="0.5", lw=1)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Spearman / partial ρ")
    ax.set_title("NSCLC RNA: does TJ epithelial attenuate TROP2–immune?")
    ax.legend(frameon=False)
    fig.savefig(figdir / "fig_bridge_partials.png", dpi=160)
    fig.savefig(figdir / "fig_bridge_partials.pdf")
    plt.close(fig)


def write_finding(path: Path, key: dict, tab: pd.DataFrame) -> None:
    lines: list[str] = []
    lines.append("# DepMap/CCLE PPT funnel: TROP2–CLDN4/TJ protein + immune scores in lung lines")
    lines.append("")
    lines.append("Numbers below are written by `analyze.py` from `tables/key_stats.json` and `tables/correlations.csv`. They are not typed by hand.")
    lines.append("")
    lines.append("## Question")
    lines.append("")
    lines.append("For a PPT resistance / tight-junction bridge slide, do **public DepMap/CCLE lung cell lines** support:")
    lines.append("")
    lines.append("1. TROP2 coexpression with CLDN4 protein and locked TJ gene scores, and")
    lines.append("2. TROP2 coexpression with cancer-cell IFN / MHC-I / immune-ligand scores,")
    lines.append("3. with TJ scores explaining (attenuating) the TROP2–immune associations?")
    lines.append("")
    lines.append("Cultured lines have no T-cell infiltrate and no IFN treatment. Partials are association contrasts, not mediation.")
    lines.append("")
    lines.append("## Locked protein check (do not change)")
    lines.append("")
    chk = key["protein_checksum_lung_suffix_columns"]
    lines.append(
        f"Gygi TenPx `_LUNG_` columns, TACSTD2 and CLDN4 both quantified, replicates not collapsed: "
        f"**Spearman ρ = {fmt(chk.get('spearman_rho'))}** (n={chk.get('n_complete')}, p={fmt_p(chk.get('p'))}). "
        f"Rounds to 0.69 = {chk.get('rounds_to_0.69')}."
    )
    lines.append("")
    model_pair = key["protein_model_lung_TACSTD2_CLDN4"]
    lines.append(
        f"Model-mapped lung cell lines with both proteins: ρ = {fmt(model_pair.get('rho'))} "
        f"(n={model_pair.get('n')}, p={fmt_p(model_pair.get('p'))})."
    )
    lines.append("")
    rppa = key.get("rppa", {})
    lines.append(
        f"RPPA antibody info: mentions CLDN4 = {rppa.get('mentions_CLDN4')}, "
        f"TACSTD2/TROP2 = {rppa.get('mentions_TACSTD2_or_TROP2')}, Claudin-7 = {rppa.get('mentions_Claudin7')}. "
        "RPPA is not used for TROP2–CLDN4."
    )
    lines.append("")
    lines.append("## Protein: TROP2 vs TJ genes / scores")
    lines.append("")
    lines.append(
        f"Gygi lung cell lines mapped to DepMap models: n = {key['n_protein_lung']}. "
        f"With TACSTD2 protein: n = {key['n_protein_trop2']}. "
        f"Protein TJ genes used in TJ_EPITHELIAL score: {', '.join(key['protein_score_cov']['TJ_EPITHELIAL']['used']) or 'none'}."
    )
    lines.append("")
    lines.append("| cohort | partner | Spearman |")
    lines.append("|---|---|---|")
    for cohort in ("lung", "NSCLC"):
        for y in PRIMARY_TJ + ["CLDN7", "CLDN1", "OCLN", "TJP1", "F11R", "CDH1", "EPCAM"]:
            rec = pick(tab, "protein", cohort, "TACSTD2", y)
            if rec is None:
                continue
            lines.append(f"| {cohort} | {SCORE_LABEL.get(y, y)} | {cell(rec)} |")
    lines.append("")
    lines.append("## RNA: TROP2 vs TJ scores (DepMap 24Q4 lung cell lines)")
    lines.append("")
    lines.append(
        f"Lung cell lines with TACSTD2 RNA: n = {key['n_rna_lung']} "
        f"(NSCLC {key['n_rna_nsclc']}, LUAD {key['n_rna_luad']}). "
        f"TJ epithelial genes used: {key['rna_score_cov']['TJ_EPITHELIAL']['n_used']}/"
        f"{key['rna_score_cov']['TJ_EPITHELIAL']['n_requested']}."
    )
    lines.append("")
    lines.append("| cohort | TJ score | Spearman | after EPCAM | after keratins |")
    lines.append("|---|---|---|---|---|")
    for cohort in ("lung", "NSCLC", "LUAD"):
        for y in PRIMARY_TJ:
            u = pick(tab, "rna", cohort, "TACSTD2", y)
            ep = pick(tab, "rna", cohort, "TACSTD2", y, "EPCAM")
            kr = pick(tab, "rna", cohort, "TACSTD2", y, "KRT8+KRT18+KRT19")
            lines.append(
                f"| {cohort} | {SCORE_LABEL.get(y, y)} | {cell(u)} | {cell(ep)} | {cell(kr)} |"
            )
    lines.append("")
    lines.append("## RNA: TROP2 vs immune-related scores")
    lines.append("")
    lines.append(
        f"Hallmark IFN-γ genes used: {key['rna_score_cov']['IFNG']['n_used']}/"
        f"{key['rna_score_cov']['IFNG']['n_requested']}. "
        f"Immune-ligand members used: {', '.join(key['rna_score_cov']['IMMUNE']['used'])}."
    )
    lines.append("")
    lines.append("| cohort | score | TROP2 Spearman |")
    lines.append("|---|---|---|")
    for cohort in ("lung", "NSCLC", "LUAD"):
        for y in PRIMARY_IMMUNE + ["ISG", "CD274"]:
            lines.append(f"| {cohort} | {SCORE_LABEL.get(y, y)} | {cell(pick(tab, 'rna', cohort, 'TACSTD2', y))} |")
    lines.append("")
    lines.append("## Bridge: does TJ attenuate TROP2–immune (RNA NSCLC)?")
    lines.append("")
    lines.append(
        "Primary FDR family = the six NSCLC RNA partials "
        "(TROP2 vs IFNG/MHC1/IMMUNE given TJ_EPITHELIAL, and TROP2 vs TJ_EPITHELIAL/TJ_TISMO/CLDN4 given IFNG)."
    )
    lines.append("")
    lines.append("| immune score | TROP2 unadjusted | TROP2 \\| TJ epithelial | TROP2 \\| CLDN4 | TJ epithelial \\| IFNG context |")
    lines.append("|---|---|---|---|---|")
    for y in PRIMARY_IMMUNE:
        u = pick(tab, "rna", "NSCLC", "TACSTD2", y)
        tj = pick(tab, "rna", "NSCLC", "TACSTD2", y, "TJ_EPITHELIAL")
        c4 = pick(tab, "rna", "NSCLC", "TACSTD2", y, "CLDN4")
        # For MHC/immune the last column is less meaningful; still report TJ vs that score.
        tj_vs = pick(tab, "rna", "NSCLC", "TJ_EPITHELIAL", y)
        lines.append(
            f"| {SCORE_LABEL[y]} | {cell(u)} | {cell(tj)} | {cell(c4)} | {cell(tj_vs)} |"
        )
    lines.append("")
    lines.append("TROP2 vs TJ scores given IFN-γ (does immune score explain the TJ coexpression?):")
    lines.append("")
    lines.append("| TJ score | TROP2 unadjusted | TROP2 \\| IFNG |")
    lines.append("|---|---|---|")
    for y in PRIMARY_TJ:
        lines.append(
            f"| {SCORE_LABEL[y]} | {cell(pick(tab, 'rna', 'NSCLC', 'TACSTD2', y))} | "
            f"{cell(pick(tab, 'rna', 'NSCLC', 'TACSTD2', y, 'IFNG'))} |"
        )
    lines.append("")
    lines.append("## Protein immune layer (honest n)")
    lines.append("")
    lines.append(
        f"Protein immune-ligand members used: {', '.join(key['protein_score_cov']['IMMUNE']['used']) or 'none'}. "
        "CXCL9/10/11 are often absent from Gygi."
    )
    lines.append("")
    lines.append("| cohort | score | TROP2 Spearman | TROP2 \\| CLDN4 |")
    lines.append("|---|---|---|---|")
    for cohort in ("lung", "NSCLC"):
        for y in PRIMARY_IMMUNE:
            lines.append(
                f"| {cohort} | {SCORE_LABEL[y]} | {cell(pick(tab, 'protein', cohort, 'TACSTD2', y))} | "
                f"{cell(pick(tab, 'protein', cohort, 'TACSTD2', y, 'CLDN4'))} |"
            )
    lines.append("")
    lines.append("## PPT claim language (honest)")
    lines.append("")
    for line in key["ppt_claims"]:
        lines.append(f"- {line}")
    lines.append("")
    lines.append("## Not claimed")
    lines.append("")
    lines.append("- ICI resistance in patients (no response labels in DepMap).")
    lines.append("- Immune exclusion or spatial TJ barrier (cell lines only).")
    lines.append("- RPPA n≈118 TROP2–CLDN4 (no TROP2/CLDN4 antibodies).")
    lines.append("- Private KL scRNA, KD co-culture, or PDX protein.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 scripts/depmap_trop2_tj_ppt_bridge/download.py")
    lines.append("python3 scripts/depmap_trop2_tj_ppt_bridge/analyze.py")
    lines.append("```")
    lines.append("")
    path.write_text("\n".join(lines) + "\n")


def build_ppt_claims(key: dict, tab: pd.DataFrame) -> list[str]:
    claims = []
    chk = key["protein_checksum_lung_suffix_columns"]
    if chk.get("rounds_to_0.69"):
        claims.append(
            f"KEEP on PPT: Gygi lung TROP2–CLDN4 protein Spearman ρ={fmt(chk['spearman_rho'])} "
            f"(n={chk['n_complete']}, p={fmt_p(chk['p'])})."
        )
    else:
        claims.append(
            f"CHECK protein checksum: ρ={fmt(chk.get('spearman_rho'))} n={chk.get('n_complete')} "
            f"(expected ~0.69)."
        )

    def ok(layer: str, cohort: str, y: str) -> dict | None:
        rec = pick(tab, layer, cohort, "TACSTD2", y)
        if rec and np.isfinite(rec.get("rho", np.nan)) and rec.get("p", 1) < 0.05:
            return rec
        return None

    tj_hits = [y for y in PRIMARY_TJ if ok("rna", "NSCLC", y)]
    if tj_hits:
        bits = []
        for y in tj_hits:
            r = ok("rna", "NSCLC", y)
            bits.append(f"{SCORE_LABEL[y]} ρ={r['rho']:+.3f} (n={int(r['n'])}, p={fmt_p(r['p'])})")
        claims.append("KEEP on PPT (RNA NSCLC): TROP2 coexpresses with " + "; ".join(bits) + ".")
    else:
        claims.append("DO NOT claim RNA TJ coexpression in NSCLC from this run (no PRIMARY_TJ hit p<0.05).")

    imm_hits = [y for y in PRIMARY_IMMUNE if ok("rna", "NSCLC", y)]
    if imm_hits:
        bits = []
        for y in imm_hits:
            r = ok("rna", "NSCLC", y)
            bits.append(f"{SCORE_LABEL[y]} ρ={r['rho']:+.3f} (n={int(r['n'])}, p={fmt_p(r['p'])})")
        claims.append("KEEP on PPT (RNA NSCLC): TROP2 coexpresses with " + "; ".join(bits) + ".")
    else:
        claims.append("DO NOT claim RNA immune-score coexpression in NSCLC from this run.")

    # Bridge attenuation summary
    bridge_notes = []
    for y in PRIMARY_IMMUNE:
        u = pick(tab, "rna", "NSCLC", "TACSTD2", y)
        p = pick(tab, "rna", "NSCLC", "TACSTD2", y, "TJ_EPITHELIAL")
        if u and p and np.isfinite(u.get("rho", np.nan)) and u.get("p", 1) < 0.05:
            call = p.get("attenuation_call") or "NA"
            bridge_notes.append(
                f"{SCORE_LABEL[y]}: unadj {u['rho']:+.3f} → |TJ epithelial {p['rho']:+.3f} ({call})"
            )
    if bridge_notes:
        claims.append("Bridge (TROP2–immune | TJ epithelial, RNA NSCLC): " + "; ".join(bridge_notes) + ".")
    else:
        claims.append("Bridge: no unadjusted TROP2–immune association to test for TJ attenuation in NSCLC RNA.")

    # Explicit resistance caveat
    claims.append(
        "DO NOT write 'PPT resistance' as a DepMap result: these are untreated cultured lines without ICI labels."
    )
    return claims


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/depmap_trop2_tj_ppt_bridge")
    ap.add_argument("--outdir", default="results/depmap_trop2_tj_ppt_bridge")
    args = ap.parse_args()
    data = Path(args.data)
    outdir = Path(args.outdir)
    tables = outdir / "tables"
    tables.mkdir(parents=True, exist_ok=True)

    model = load_models(data / "Model.csv")
    ifng = [g for g in (data / "hallmark_ifng_genes.txt").read_text().splitlines() if g]
    expr = pd.read_csv(data / "expression_panel.csv")
    gygi = pd.read_csv(data / "gygi_protein_panel.csv")

    rna = expr.merge(model, on="ModelID", how="left")
    rna = annotate(rna)
    rna = rna[(rna["OncotreeLineage"] == "Lung") & (rna["ModelType"] == "Cell Line")].copy()
    rna, rna_cov = add_scores(rna, ifng, "rna")

    prot_all = models_from_gygi(gygi, model)
    prot_all = annotate(prot_all)
    prot = prot_all[(prot_all["OncotreeLineage"] == "Lung") & (prot_all["ModelType"] == "Cell Line")].copy()
    prot, prot_cov = add_scores(prot, ifng, "protein")

    checksum = gygi_checksum(gygi)
    both = prot[["TACSTD2", "CLDN4"]].dropna() if "TACSTD2" in prot.columns and "CLDN4" in prot.columns else pd.DataFrame()
    if len(both):
        rho_m, p_m, n_m = spearman_full(both["TACSTD2"], both["CLDN4"])
    else:
        rho_m, p_m, n_m = float("nan"), float("nan"), 0

    rows: list[dict] = []
    gene_partners_rna = PRIMARY_TJ + [
        "TJ_EPITHELIAL_noCLDN4",
        "TJ_TISMO_noCLDN4",
        "CLDN4_TJ_EDGE_noCLDN4",
        "CLDN1",
        "CLDN7",
        "OCLN",
        "TJP1",
        "F11R",
        "CDH1",
        "EPCAM",
    ]
    gene_partners_prot = [
        "CLDN4",
        "CLDN7",
        "CLDN1",
        "CLDN3",
        "OCLN",
        "TJP1",
        "TJP2",
        "F11R",
        "CDH1",
        "EPCAM",
        "TJ_EPITHELIAL",
        "TJ_TISMO",
        "CLDN4_TJ_EDGE",
        "TJ_EPITHELIAL_noCLDN4",
        "TJ_TISMO_noCLDN4",
        "CLDN4_TJ_EDGE_noCLDN4",
    ]

    # RNA unadjusted + controls
    for cohort_name, mask in (
        ("lung", rna.index == rna.index),
        ("NSCLC", rna["is_NSCLC"]),
        ("LUAD", rna["group"] == "LUAD"),
    ):
        frame = rna.loc[mask].copy()
        # re-z for cohort-specific scores where needed: keep global scores for comparability,
        # but also score inside cohort for primary immune/TJ (documented).
        frame_c, _ = add_scores(frame, ifng, "rna")
        for y in gene_partners_rna + PRIMARY_IMMUNE + ["ISG", "CD274"]:
            primary = cohort_name == "NSCLC" and y in (PRIMARY_TJ + PRIMARY_IMMUNE)
            record_pair(
                rows, frame_c, cohort_name, "rna", "TACSTD2", y, [],
                family="rna_nsclc_primary" if primary else "rna_exploratory",
                primary=primary,
                do_boot=primary or (cohort_name == "NSCLC" and y in PRIMARY_TJ + PRIMARY_IMMUNE),
            )
            if y in PRIMARY_TJ:
                for covs in (["EPCAM"], ["KRT8", "KRT18", "KRT19"], ["IFNG"], ["IMMUNE"]):
                    if all(c in frame_c.columns for c in covs):
                        record_pair(
                            rows, frame_c, cohort_name, "rna", "TACSTD2", y, covs,
                            family="rna_nsclc_primary" if (cohort_name == "NSCLC" and covs == ["IFNG"] and y in PRIMARY_TJ) else "rna_partials",
                            primary=bool(cohort_name == "NSCLC" and covs == ["IFNG"] and y in PRIMARY_TJ),
                            do_boot=cohort_name == "NSCLC",
                        )
            if y in PRIMARY_IMMUNE:
                for covs in (["CLDN4"], ["TJ_EPITHELIAL"], ["TJ_TISMO"], ["CLDN4_TJ_EDGE"], ["EPCAM"]):
                    if all(c in frame_c.columns for c in covs):
                        primary_b = cohort_name == "NSCLC" and covs == ["TJ_EPITHELIAL"]
                        record_pair(
                            rows, frame_c, cohort_name, "rna", "TACSTD2", y, covs,
                            family="rna_nsclc_primary" if primary_b else "rna_partials",
                            primary=primary_b,
                            do_boot=cohort_name == "NSCLC",
                        )
        # TJ score vs immune score
        for tj in ["TJ_EPITHELIAL", "TJ_TISMO", "CLDN4", "CLDN4_TJ_EDGE"]:
            for imm in PRIMARY_IMMUNE:
                record_pair(
                    rows, frame_c, cohort_name, "rna", tj, imm, [],
                    family="rna_tj_vs_immune",
                    primary=False,
                    do_boot=cohort_name == "NSCLC",
                )

    # Protein
    for cohort_name, mask in (
        ("lung", prot.index == prot.index),
        ("NSCLC", prot["is_NSCLC"]),
        ("LUAD", prot["group"] == "LUAD"),
    ):
        frame = prot.loc[mask].copy()
        frame_c, _ = add_scores(frame, ifng, "protein")
        for y in gene_partners_prot + PRIMARY_IMMUNE + ["ISG"]:
            primary = False
            record_pair(
                rows, frame_c, cohort_name, "protein", "TACSTD2", y, [],
                family="protein_exploratory",
                primary=primary,
                do_boot=cohort_name in ("lung", "NSCLC") and y in ("CLDN4", "TJ_EPITHELIAL", "IFNG", "IMMUNE"),
            )
            if y in PRIMARY_IMMUNE and "CLDN4" in frame_c.columns:
                record_pair(
                    rows, frame_c, cohort_name, "protein", "TACSTD2", y, ["CLDN4"],
                    family="protein_partials",
                    primary=False,
                    do_boot=cohort_name == "lung",
                )
            if y in ("TJ_EPITHELIAL", "TJ_TISMO", "CLDN4_TJ_EDGE") and "EPCAM" in frame_c.columns:
                record_pair(
                    rows, frame_c, cohort_name, "protein", "TACSTD2", y, ["EPCAM"],
                    family="protein_partials",
                    primary=False,
                    do_boot=False,
                )

    tab = apply_fdr_and_calls(pd.DataFrame(rows))
    tab.to_csv(tables / "correlations.csv", index=False)

    rna_out_cols = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "CCLEName",
            "group",
            "is_NSCLC",
            "TACSTD2",
            "CLDN4",
            "EPCAM",
            "TJ_EPITHELIAL",
            "TJ_TISMO",
            "CLDN4_TJ_EDGE",
            "IFNG",
            "MHC1",
            "IMMUNE",
            "ISG",
            "CD274",
        ]
        if c in rna.columns
    ]
    rna[rna_out_cols].to_csv(tables / "rna_lung_lines.csv", index=False)
    prot_out_cols = [
        c
        for c in [
            "CellLineName",
            "CCLEName",
            "group",
            "is_NSCLC",
            "TACSTD2",
            "CLDN4",
            "CLDN7",
            "CLDN1",
            "EPCAM",
            "TJ_EPITHELIAL",
            "TJ_TISMO",
            "CLDN4_TJ_EDGE",
            "IFNG",
            "MHC1",
            "IMMUNE",
        ]
        if c in prot.columns
    ]
    prot[prot_out_cols].to_csv(tables / "protein_lung_lines.csv", index=False)

    rppa_path = data / "CCLE_RPPA_Ab_info_20180123.csv"
    if rppa_path.exists():
        text = rppa_path.read_text(errors="replace")
        rppa = {
            "mentions_CLDN4": "CLDN4" in text or "Claudin-4" in text or "Claudin 4" in text,
            "mentions_TACSTD2_or_TROP2": "TACSTD2" in text or "TROP2" in text,
            "mentions_Claudin7": "Claudin-7" in text or "Claudin 7" in text,
        }
    else:
        rppa = {"error": "file missing"}

    key = {
        "release": "DepMap Public 24Q4",
        "protein_source": "Gygi/Nusinow CCLE TMT protein_quant_current_normalized",
        "protein_checksum_lung_suffix_columns": checksum,
        "protein_model_lung_TACSTD2_CLDN4": {"rho": rho_m, "p": p_m, "n": n_m},
        "n_rna_lung": int(len(rna)),
        "n_rna_nsclc": int(rna["is_NSCLC"].sum()),
        "n_rna_luad": int((rna["group"] == "LUAD").sum()),
        "n_protein_lung": int(len(prot)),
        "n_protein_trop2": int(prot["TACSTD2"].notna().sum()) if "TACSTD2" in prot.columns else 0,
        "n_protein_both": int(len(both)),
        "rna_score_cov": rna_cov,
        "protein_score_cov": prot_cov,
        "rppa": rppa,
        "n_boot": N_BOOT,
    }
    key["ppt_claims"] = build_ppt_claims(key, tab)
    (tables / "key_stats.json").write_text(json.dumps(key, indent=2, default=float) + "\n")

    make_figures(outdir, rna, prot, tab)
    write_finding(outdir / "FINDING.md", key, tab)

    # Verdict one-liner from computed claims
    keep = [c for c in key["ppt_claims"] if c.startswith("KEEP")]
    verdict = " | ".join(keep) if keep else "No KEEP claims from this DepMap PPT funnel run."
    (outdir / "verdict.txt").write_text(verdict + "\n")
    print(verdict, flush=True)
    print(f"wrote {outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
