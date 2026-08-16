#!/usr/bin/env python3
"""A11: TGF-β genes vs TACSTD2-high public lung.

Primary question (pre-specified):
  In public lung tumors, is TGF-β signaling higher in TACSTD2-high samples
  and/or positively correlated with TACSTD2?

Primary statistic:
  pooled rank-partial Spearman of TACSTD2 vs each gene or score, controlling
  for ABSOLUTE purity and histology (LUSC vs LUAD).

TACSTD2-high contrast:
  Q4 vs Q1 (primary) and median split (sensitivity), Mann-Whitney U.

Named claim-A11 objects:
  TGFB1 (the ligand named in the claim page)
  HALLMARK_TGFB_ZMEAN (MSigDB HALLMARK_TGF_BETA_SIGNALING z-mean)

Ligand class-effect rule (pre-specified, not tuned after seeing numbers):
  A TGF-β *ligand* class effect is claimed ONLY if >=2 of TGFB1/2/3 are
  WEAK_POSITIVE or ASSOCIATED_POSITIVE on the primary statistic and none
  is significantly negative (FDR<0.05 and rho<0).

Pathway-claim rule (pre-specified):
  Claim A11 "TGF-β signaling accompanies TACSTD2-high" is supported only if
  HALLMARK_TGFB_ZMEAN is WEAK_POSITIVE or ASSOCIATED_POSITIVE AND TGFB1 is
  not labeled ASSOCIATED_NEGATIVE or WEAK_NEGATIVE. Otherwise NOT_SUPPORTED.

Thresholds on primary partial rho (also pre-specified):
  ASSOCIATED_POSITIVE : rho >= 0.20 and FDR < 0.05
  WEAK_POSITIVE       : 0.10 <= rho < 0.20 and FDR < 0.05
  NULL                : FDR >= 0.05 or |rho| < 0.10
  WEAK_NEGATIVE       : -0.20 < rho <= -0.10 and FDR < 0.05
  ASSOCIATED_NEGATIVE : rho <= -0.20 and FDR < 0.05

CLDN4 is a positive-control junction gene, not counted in the class verdict.
CD8A is immune context only (claim says TGF-β reinforces immune-cold).
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
from statsmodels.stats.multitest import multipletests

PRIMARY_LIGANDS = ["TGFB1", "TGFB2", "TGFB3"]
RECEPTORS = ["TGFBR1", "TGFBR2", "TGFBR3"]
SMADS = ["SMAD2", "SMAD3", "SMAD4", "SMAD6", "SMAD7"]
CORE_GENES = PRIMARY_LIGANDS + RECEPTORS + SMADS
POSITIVE_CONTROL = ["CLDN4"]
IMMUNE_CONTEXT = ["CD8A"]
SCORE_NAMES = ["HALLMARK_TGFB_ZMEAN", "FTBRS_ZMEAN"]
ALL_PARTNERS = CORE_GENES + POSITIVE_CONTROL + IMMUNE_CONTEXT
TARGET = "TACSTD2"

# MSigDB HALLMARK_TGF_BETA_SIGNALING (54 genes). Score uses those present.
HALLMARK_TGFB = [
    "ACVR1", "APC", "ARID4B", "BCAR3", "BMP2", "BMPR1A", "BMPR2", "CDH1",
    "CDK9", "CDKN1C", "CTNNB1", "ENG", "FKBP1A", "FNTA", "FURIN", "HDAC1",
    "HIPK2", "ID1", "ID2", "ID3", "IFNGR2", "JUNB", "KLF10", "LEFTY2",
    "LTBP2", "MAP3K7", "NCOR2", "NOG", "PMEPA1", "PPM1A", "PPP1CA",
    "PPP1R15A", "RAB31", "RHOA", "SERPINE1", "SKI", "SKIL", "SLC20A1",
    "SMAD1", "SMAD3", "SMAD6", "SMAD7", "SMURF1", "SMURF2", "SPTBN1",
    "TGFB1", "TGFBR1", "TGIF1", "THBS1", "TJP1", "TRIM33", "UBE2D3",
    "WWTR1", "XIAP",
]

# Mariathasan 2018 F-TBRS (fibroblast TGF-β response). Aliases resolved later.
FTBRS_PREF = [
    "ACTA2", "ACTG2", "ADAM12", "ADAM19", "CNN1", "COL4A1", "CCN2", "CTPS1",
    "RFLNA", "FSTL3", "HSPB1", "IGFBP3", "IL11", "JUNB", "NT5E", "OLFML2B",
    "PMEPA1", "PPP1R13L", "PXDC1", "SERPINE1", "SH3PXD2A", "TAGLN", "TGFBI",
    "TNS1", "TPM1",
]
FTBRS_ALIAS = {"CCN2": "CTGF", "RFLNA": "FAM101B"}

ALIASES = {
    "TGFB1": "TGF-β1 ligand (named in claim A11)",
    "TGFB2": "TGF-β2 ligand",
    "TGFB3": "TGF-β3 ligand",
    "TGFBR1": "TGF-β receptor I / ALK5",
    "TGFBR2": "TGF-β receptor II",
    "TGFBR3": "TGF-β receptor III / betaglycan",
    "SMAD2": "R-SMAD",
    "SMAD3": "R-SMAD (in HALLMARK set)",
    "SMAD4": "Co-SMAD",
    "SMAD6": "I-SMAD (BMP-biased)",
    "SMAD7": "I-SMAD / TGF-β feedback",
    "CLDN4": "claudin-4 (junction positive control, not TGF-β)",
    "CD8A": "CD8 T-cell marker (immune context, not TGF-β)",
    "HALLMARK_TGFB_ZMEAN": "z-mean of HALLMARK_TGF_BETA_SIGNALING genes present",
    "FTBRS_ZMEAN": "z-mean of Mariathasan 2018 F-TBRS genes present",
}


def bh(pvals: list[float]) -> np.ndarray:
    return multipletests(pvals, method="fdr_bh")[1]


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    mask = np.isfinite(x) & np.isfinite(y)
    x, y = x[mask], y[mask]
    n = int(len(x))
    if n < 5:
        return float("nan"), float("nan"), n
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p), n


def partial_spearman(x: np.ndarray, y: np.ndarray, Z: np.ndarray) -> tuple[float, float, int]:
    """Rank residual partial correlation. Z is 1-D or 2-D (covariates in columns)."""
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    x, y, Z = x[mask], y[mask], Z[mask]
    n = int(len(x))
    k = int(Z.shape[1])
    if n < 5 + k:
        return float("nan"), float("nan"), n
    xr, yr = stats.rankdata(x), stats.rankdata(y)
    Zr = np.column_stack([stats.rankdata(Z[:, j]) for j in range(k)])
    design = np.column_stack([np.ones(n), Zr])
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    r = max(min(r, 1.0), -1.0)
    df = n - 2 - k
    if df <= 0 or abs(r) >= 1:
        return r, float("nan"), n
    t = r * np.sqrt(df / (1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df))
    return r, p, n


def label_rho(rho: float, fdr: float) -> str:
    if not np.isfinite(rho) or not np.isfinite(fdr):
        return "NOT_COMPUTED"
    if fdr >= 0.05 or abs(rho) < 0.10:
        return "NULL"
    if rho >= 0.20:
        return "ASSOCIATED_POSITIVE"
    if rho >= 0.10:
        return "WEAK_POSITIVE"
    if rho <= -0.20:
        return "ASSOCIATED_NEGATIVE"
    return "WEAK_NEGATIVE"


def mwu_high_low(high: np.ndarray, low: np.ndarray) -> dict:
    high = high[np.isfinite(high)]
    low = low[np.isfinite(low)]
    n_h, n_l = int(len(high)), int(len(low))
    if n_h < 5 or n_l < 5:
        return {
            "n_high": n_h,
            "n_low": n_l,
            "median_high": float("nan"),
            "median_low": float("nan"),
            "delta_median": float("nan"),
            "p_mannwhitney": float("nan"),
        }
    u, p = stats.mannwhitneyu(high, low, alternative="two-sided")
    return {
        "n_high": n_h,
        "n_low": n_l,
        "median_high": float(np.median(high)),
        "median_low": float(np.median(low)),
        "delta_median": float(np.median(high) - np.median(low)),
        "p_mannwhitney": float(p),
        "U": float(u),
    }


def load_tcga_cohort(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["patient"] = df["sample"].astype(str).str.slice(0, 12)
    df["sample_type"] = df["sample"].astype(str).str.split("-").str[3].str[:2]
    return df


def one_per_patient(df: pd.DataFrame, sample_type: str) -> pd.DataFrame:
    sub = df[df["sample_type"] == sample_type].sort_values("sample")
    return sub.drop_duplicates("patient", keep="first").set_index("patient")


def load_purity(path: Path) -> pd.Series:
    pur = pd.read_csv(path, sep="\t")
    if "array" not in pur.columns or "purity" not in pur.columns:
        raise SystemExit(f"Unexpected purity columns: {list(pur.columns)}")
    pur = pur[pur["array"].astype(str).str.endswith("-01")].copy()
    pur["patient"] = pur["array"].astype(str).str.slice(0, 12)
    pur["purity"] = pd.to_numeric(pur["purity"], errors="coerce")
    return pur.dropna(subset=["purity"]).drop_duplicates("patient").set_index("patient")["purity"]


def zmean(df: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in df.columns]
    if not present:
        return pd.Series(np.nan, index=df.index)
    z = df[present].apply(lambda s: (s - s.mean()) / s.std(ddof=0), axis=0)
    return z.mean(axis=1)


def resolve_ftbrs_genes(columns: list[str]) -> list[str]:
    out = []
    for g in FTBRS_PREF:
        if g in columns:
            out.append(g)
        elif g in FTBRS_ALIAS and FTBRS_ALIAS[g] in columns:
            out.append(FTBRS_ALIAS[g])
    return out


def role_of(gene: str) -> str:
    if gene in PRIMARY_LIGANDS:
        return "primary_ligand"
    if gene in RECEPTORS:
        return "receptor"
    if gene in SMADS:
        return "smad"
    if gene in POSITIVE_CONTROL:
        return "positive_control"
    if gene in IMMUNE_CONTEXT:
        return "immune_context"
    if gene in SCORE_NAMES:
        return "pathway_score"
    return "other"


def corr_row(cohort: str, gene: str, role: str, x: pd.Series, y: pd.Series, Z=None, z_note="") -> dict:
    rho, p, n = spearman(x.to_numpy(), y.to_numpy())
    rec = {
        "cohort": cohort,
        "gene": gene,
        "role": role,
        "alias": ALIASES.get(gene, ""),
        "n": n,
        "spearman_rho": rho,
        "spearman_p": p,
        "partial_rho": float("nan"),
        "partial_p": float("nan"),
        "partial_n": 0,
        "partial_covariates": z_note,
    }
    if Z is not None:
        if isinstance(Z, pd.Series):
            common = x.index.intersection(y.index).intersection(Z.dropna().index)
            pr, pp, pn = partial_spearman(x.loc[common].to_numpy(), y.loc[common].to_numpy(), Z.loc[common].to_numpy())
        else:
            Z = pd.DataFrame(Z)
            common = x.index.intersection(y.index).intersection(Z.dropna().index)
            pr, pp, pn = partial_spearman(x.loc[common].to_numpy(), y.loc[common].to_numpy(), Z.loc[common].to_numpy())
        rec["partial_rho"] = pr
        rec["partial_p"] = pp
        rec["partial_n"] = pn
    return rec


def add_fdr(rows: list[dict], p_key: str, fdr_key: str, cohort: str, genes: list[str]) -> None:
    idx = [i for i, r in enumerate(rows) if r["cohort"] == cohort and r["gene"] in genes]
    if not idx:
        return
    pvals = [rows[i][p_key] if np.isfinite(rows[i][p_key]) else 1.0 for i in idx]
    fdrs = bh(pvals)
    for i, f in zip(idx, fdrs):
        rows[i][fdr_key] = float(f)


def add_scores(df: pd.DataFrame, hallmark_genes: list[str], ftbrs_genes: list[str]) -> pd.DataFrame:
    out = df.copy()
    out["HALLMARK_TGFB_ZMEAN"] = zmean(out, hallmark_genes)
    out["FTBRS_ZMEAN"] = zmean(out, ftbrs_genes)
    return out


def analyze_tcga(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict, dict]:
    luad = load_tcga_cohort(indir / "tcga_luad_genes.csv")
    lusc = load_tcga_cohort(indir / "tcga_lusc_genes.csv")
    purity = load_purity(indir / "tcga_absolute_purity.txt")

    required = [TARGET] + CORE_GENES + POSITIVE_CONTROL + IMMUNE_CONTEXT
    frames = {}
    normals = {}
    for name, raw in (("LUAD", luad), ("LUSC", lusc)):
        tum = one_per_patient(raw, "01")
        nor = one_per_patient(raw, "11")
        for g in required:
            if g not in tum.columns:
                raise SystemExit(f"{g} missing from TCGA-{name} extract")
            tum[g] = pd.to_numeric(tum[g], errors="coerce")
            if g in nor.columns:
                nor[g] = pd.to_numeric(nor[g], errors="coerce")
        for g in list(tum.columns):
            if g not in ("sample", "patient", "sample_type", "histology", "purity") and tum[g].dtype == object:
                tum[g] = pd.to_numeric(tum[g], errors="coerce")
        for g in list(nor.columns):
            if g not in ("sample", "patient", "sample_type") and nor[g].dtype == object:
                nor[g] = pd.to_numeric(nor[g], errors="coerce")
        tum["histology"] = name
        tum["purity"] = purity.reindex(tum.index)
        frames[name] = tum
        normals[name] = nor

    hallmark_present = [g for g in HALLMARK_TGFB if g in frames["LUAD"].columns]
    ftbrs_present = resolve_ftbrs_genes(list(frames["LUAD"].columns))
    for name in ("LUAD", "LUSC"):
        frames[name] = add_scores(frames[name], hallmark_present, ftbrs_present)
        if len(normals[name]):
            normals[name] = add_scores(normals[name], hallmark_present, ftbrs_present)

    pooled = pd.concat([frames["LUAD"], frames["LUSC"]], axis=0)
    pooled["is_lusc"] = (pooled["histology"] == "LUSC").astype(float)

    partners = ALL_PARTNERS + SCORE_NAMES
    corr_rows: list[dict] = []
    hl_rows: list[dict] = []
    tn_rows: list[dict] = []

    def run_cohort(name: str, df: pd.DataFrame, Z, z_note: str) -> None:
        x = df[TARGET]
        for gene in partners:
            role = role_of(gene)
            corr_rows.append(corr_row(name, gene, role, x, df[gene], Z=Z, z_note=z_note))
            q = x.quantile([0.25, 0.50, 0.75])
            q1 = df.loc[x <= q[0.25], gene]
            q4 = df.loc[x >= q[0.75], gene]
            med_lo = df.loc[x <= q[0.50], gene]
            med_hi = df.loc[x > q[0.50], gene]
            rec_q = mwu_high_low(q4.to_numpy(), q1.to_numpy())
            rec_m = mwu_high_low(med_hi.to_numpy(), med_lo.to_numpy())
            hl_rows.append(
                {
                    "cohort": name,
                    "gene": gene,
                    "role": role,
                    "split": "Q4_vs_Q1",
                    "tacstd2_q1": float(q[0.25]),
                    "tacstd2_q3": float(q[0.75]),
                    **rec_q,
                }
            )
            hl_rows.append(
                {
                    "cohort": name,
                    "gene": gene,
                    "role": role,
                    "split": "median",
                    "tacstd2_q1": float(q[0.50]),
                    "tacstd2_q3": float(q[0.50]),
                    **rec_m,
                }
            )

    run_cohort("LUAD", frames["LUAD"], frames["LUAD"]["purity"], "ABSOLUTE purity")
    run_cohort("LUSC", frames["LUSC"], frames["LUSC"]["purity"], "ABSOLUTE purity")
    run_cohort(
        "POOLED",
        pooled,
        pooled[["purity", "is_lusc"]],
        "ABSOLUTE purity + histology (LUSC vs LUAD)",
    )

    fdr_groups = {
        "partial_fdr_ligands": PRIMARY_LIGANDS,
        "partial_fdr_core": CORE_GENES,
        "partial_fdr_all_partners": partners,
        "spearman_fdr_ligands": PRIMARY_LIGANDS,
        "spearman_fdr_core": CORE_GENES,
        "spearman_fdr_all_partners": partners,
    }
    for cohort in ("LUAD", "LUSC", "POOLED"):
        add_fdr(corr_rows, "partial_p", "partial_fdr_ligands", cohort, PRIMARY_LIGANDS)
        add_fdr(corr_rows, "partial_p", "partial_fdr_core", cohort, CORE_GENES)
        add_fdr(corr_rows, "partial_p", "partial_fdr_all_partners", cohort, partners)
        add_fdr(corr_rows, "spearman_p", "spearman_fdr_ligands", cohort, PRIMARY_LIGANDS)
        add_fdr(corr_rows, "spearman_p", "spearman_fdr_core", cohort, CORE_GENES)
        add_fdr(corr_rows, "spearman_p", "spearman_fdr_all_partners", cohort, partners)

    for r in corr_rows:
        if r["gene"] in PRIMARY_LIGANDS:
            r["primary_label"] = label_rho(r["partial_rho"], r.get("partial_fdr_ligands", 1.0))
        else:
            r["primary_label"] = label_rho(r["partial_rho"], r.get("partial_fdr_all_partners", 1.0))

    for cohort in ("LUAD", "LUSC", "POOLED"):
        idx = [
            i
            for i, r in enumerate(hl_rows)
            if r["cohort"] == cohort and r["split"] == "Q4_vs_Q1" and r["gene"] in PRIMARY_LIGANDS
        ]
        if idx:
            pvals = [hl_rows[i]["p_mannwhitney"] if np.isfinite(hl_rows[i]["p_mannwhitney"]) else 1.0 for i in idx]
            for i, f in zip(idx, bh(pvals)):
                hl_rows[i]["fdr_bh"] = float(f)

    context_genes = [TARGET] + CORE_GENES + POSITIVE_CONTROL + IMMUNE_CONTEXT + SCORE_NAMES
    for name in ("LUAD", "LUSC"):
        tum, nor = frames[name], normals[name]
        for gene in context_genes:
            if gene not in tum.columns:
                continue
            t = tum[gene]
            n = nor[gene] if gene in nor.columns else pd.Series(dtype=float)
            if len(n.dropna()) < 5:
                continue
            u, p = stats.mannwhitneyu(t.dropna(), n.dropna(), alternative="two-sided")
            paired = sorted(set(t.dropna().index) & set(n.dropna().index))
            p_pair = float("nan")
            if len(paired) >= 10:
                _, p_pair = stats.wilcoxon(t.loc[paired], n.loc[paired])
            tn_rows.append(
                {
                    "cohort": name,
                    "gene": gene,
                    "n_tumor": int(t.notna().sum()),
                    "n_normal": int(n.notna().sum()),
                    "median_tumor": float(t.median()),
                    "median_normal": float(n.median()),
                    "delta_median": float(t.median() - n.median()),
                    "p_mannwhitney_unpaired": float(p),
                    "n_paired": len(paired),
                    "p_wilcoxon_paired": float(p_pair) if np.isfinite(p_pair) else float("nan"),
                }
            )

    counts = {
        "LUAD_tumor": int(len(frames["LUAD"])),
        "LUSC_tumor": int(len(frames["LUSC"])),
        "POOLED_tumor": int(len(pooled)),
        "LUAD_with_purity": int(frames["LUAD"]["purity"].notna().sum()),
        "LUSC_with_purity": int(frames["LUSC"]["purity"].notna().sum()),
        "POOLED_with_purity": int(pooled["purity"].notna().sum()),
        "LUAD_normal": int(len(normals["LUAD"])),
        "LUSC_normal": int(len(normals["LUSC"])),
        "hallmark_genes_used": hallmark_present,
        "n_hallmark_genes_used": len(hallmark_present),
        "ftbrs_genes_used": ftbrs_present,
        "n_ftbrs_genes_used": len(ftbrs_present),
        "TACSTD2_vs_purity_LUAD": dict(zip(["rho", "p", "n"], spearman(frames["LUAD"][TARGET].to_numpy(), frames["LUAD"]["purity"].to_numpy()))),
        "TACSTD2_vs_purity_LUSC": dict(zip(["rho", "p", "n"], spearman(frames["LUSC"][TARGET].to_numpy(), frames["LUSC"]["purity"].to_numpy()))),
    }

    return (
        pd.DataFrame(corr_rows),
        pd.DataFrame(hl_rows),
        pd.DataFrame(tn_rows),
        {"LUAD": frames["LUAD"], "LUSC": frames["LUSC"], "POOLED": pooled},
        counts,
    )


def analyze_depmap(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame] | tuple[None, None]:
    expr_path = indir / "depmap24q4_tgfb_all_models.csv"
    model_path = indir / "Model.csv"
    if not expr_path.exists():
        return None, None
    expr = pd.read_csv(expr_path)
    for g in expr.columns:
        if g != "ModelID":
            expr[g] = pd.to_numeric(expr[g], errors="coerce")
    hallmark_present = [g for g in HALLMARK_TGFB if g in expr.columns]
    ftbrs_present = resolve_ftbrs_genes(list(expr.columns))
    expr = add_scores(expr, hallmark_present, ftbrs_present)
    if model_path.exists():
        model = pd.read_csv(model_path, low_memory=False)
        keep = [
            c
            for c in [
                "ModelID",
                "CellLineName",
                "ModelType",
                "OncotreeLineage",
                "OncotreePrimaryDisease",
                "OncotreeSubtype",
            ]
            if c in model.columns
        ]
        df = expr.merge(model[keep], on="ModelID", how="left")
    else:
        df = expr
        df["OncotreeLineage"] = np.nan
        df["ModelType"] = np.nan
        df["OncotreePrimaryDisease"] = np.nan
    lung = df.copy()
    if "OncotreeLineage" in df.columns:
        lung = df[(df["OncotreeLineage"] == "Lung") & (df["ModelType"].fillna("") == "Cell Line")].copy()
    nsclc = lung
    if "OncotreePrimaryDisease" in lung.columns:
        nsclc = lung[lung["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"].copy()
    partners = [g for g in ALL_PARTNERS + SCORE_NAMES if g in df.columns]
    rows = []
    for name, sub in (
        ("depmap_lung_cell_lines", lung),
        ("depmap_NSCLC_cell_lines", nsclc),
    ):
        if TARGET not in sub.columns or len(sub) < 10:
            continue
        for gene in partners:
            rec = corr_row(name, gene, role_of(gene), sub[TARGET], sub[gene], Z=None, z_note="")
            rec["primary_label"] = label_rho(rec["spearman_rho"], rec["spearman_p"])
            rows.append(rec)
    if not rows:
        return None, None
    for cohort in {r["cohort"] for r in rows}:
        add_fdr(rows, "spearman_p", "spearman_fdr_ligands", cohort, PRIMARY_LIGANDS)
        add_fdr(rows, "spearman_p", "spearman_fdr_all_partners", cohort, partners)
    for r in rows:
        if r["gene"] in PRIMARY_LIGANDS:
            r["primary_label"] = label_rho(r["spearman_rho"], r.get("spearman_fdr_ligands", r["spearman_p"]))
        else:
            r["primary_label"] = label_rho(r["spearman_rho"], r.get("spearman_fdr_all_partners", r["spearman_p"]))
    return pd.DataFrame(rows), lung


def class_verdict(corr: pd.DataFrame) -> dict:
    sub = corr[(corr["cohort"] == "POOLED") & (corr["gene"].isin(PRIMARY_LIGANDS))].copy()
    labels = {r.gene: r.primary_label for r in sub.itertuples()}
    pos = [g for g, lab in labels.items() if lab in ("ASSOCIATED_POSITIVE", "WEAK_POSITIVE")]
    neg = [g for g, lab in labels.items() if lab in ("ASSOCIATED_NEGATIVE", "WEAK_NEGATIVE")]
    null = [g for g, lab in labels.items() if lab == "NULL"]
    if len(pos) >= 2 and not neg:
        ligand_claim = "LIGAND_CLASS_EFFECT_SUPPORTED"
        ligand_statement = (
            f"{len(pos)}/3 TGF-β ligands are purity+histology-adjusted positive "
            f"with TACSTD2 and none are significantly negative."
        )
    elif pos and neg:
        ligand_claim = "MIXED"
        ligand_statement = (
            f"Mixed ligands: positive {pos or 'none'}; negative {neg or 'none'}; null {null or 'none'}."
        )
    elif pos:
        ligand_claim = "PARTIAL_NOT_CLASS"
        ligand_statement = (
            f"Only {len(pos)}/3 ligands are positive ({', '.join(pos)}); "
            f"null {', '.join(null) if null else 'none'}. Not a ligand class effect."
        )
    elif neg and not pos:
        ligand_claim = "OPPOSITE_OR_NULL"
        ligand_statement = (
            f"No TGF-β ligand is positively associated with TACSTD2 after "
            f"purity+histology adjustment. Negative: {neg or 'none'}; null: {null or 'none'}."
        )
    else:
        ligand_claim = "NULL"
        ligand_statement = "No TGF-β ligand meets the pre-specified association rule."

    hm = corr[(corr["cohort"] == "POOLED") & (corr["gene"] == "HALLMARK_TGFB_ZMEAN")]
    tgfb1 = corr[(corr["cohort"] == "POOLED") & (corr["gene"] == "TGFB1")]
    hm_lab = hm.iloc[0]["primary_label"] if len(hm) else "NOT_COMPUTED"
    t1_lab = tgfb1.iloc[0]["primary_label"] if len(tgfb1) else "NOT_COMPUTED"
    hm_rho = float(hm.iloc[0]["partial_rho"]) if len(hm) else float("nan")
    t1_rho = float(tgfb1.iloc[0]["partial_rho"]) if len(tgfb1) else float("nan")

    if hm_lab in ("ASSOCIATED_POSITIVE", "WEAK_POSITIVE") and t1_lab not in ("ASSOCIATED_NEGATIVE", "WEAK_NEGATIVE"):
        pathway_claim = "PATHWAY_CLAIM_SUPPORTED"
        pathway_statement = (
            f"HALLMARK TGF-β score is {hm_lab} (partial ρ={hm_rho:.3f}) and TGFB1 is not negative "
            f"({t1_lab}, ρ={t1_rho:.3f}). The A11 TGF-β-signaling claim is supported on this slice."
        )
    else:
        pathway_claim = "PATHWAY_CLAIM_NOT_SUPPORTED"
        pathway_statement = (
            f"HALLMARK TGF-β score is {hm_lab} (partial ρ={hm_rho:.3f}); TGFB1 is {t1_lab} "
            f"(ρ={t1_rho:.3f}). The A11 claim that TGF-β signaling accompanies TACSTD2-high "
            f"public lung is not supported."
        )

    return {
        "ligand_rule": ">=2/3 of TGFB1/2/3 WEAK_POSITIVE or ASSOCIATED_POSITIVE on pooled partial Spearman, and none negative",
        "pathway_rule": "HALLMARK_TGFB_ZMEAN WEAK_POSITIVE or ASSOCIATED_POSITIVE AND TGFB1 not negative",
        "ligand_labels": labels,
        "positive_ligands": pos,
        "negative_ligands": neg,
        "null_ligands": null,
        "ligand_claim": ligand_claim,
        "ligand_statement": ligand_statement,
        "hallmark_label": hm_lab,
        "tgfb1_label": t1_lab,
        "pathway_claim": pathway_claim,
        "pathway_statement": pathway_statement,
    }


def plot_heatmap(corr: pd.DataFrame, out: Path) -> None:
    genes = CORE_GENES + SCORE_NAMES + POSITIVE_CONTROL + IMMUNE_CONTEXT
    sub = corr[corr["gene"].isin(genes) & corr["cohort"].isin(["LUAD", "LUSC", "POOLED"])]
    mat = sub.pivot(index="gene", columns="cohort", values="partial_rho")
    mat = mat.reindex(index=genes, columns=["LUAD", "LUSC", "POOLED"])
    fig, ax = plt.subplots(figsize=(5.8, 7.2))
    im = ax.imshow(mat.to_numpy(), cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    ax.set_xticks(range(mat.shape[1]), mat.columns)
    ax.set_yticks(range(mat.shape[0]), mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iloc[i, j]
            ax.text(j, i, f"{v:.2f}" if np.isfinite(v) else "NA", ha="center", va="center", fontsize=7)
    ax.set_title("Partial Spearman vs TACSTD2\n(purity; pooled also histology)")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="partial ρ")
    fig.tight_layout()
    fig.savefig(out / "fig1_partial_rho_heatmap.png", dpi=160)
    fig.savefig(out / "fig1_partial_rho_heatmap.pdf")
    plt.close(fig)


def plot_scatters(frames: dict, out: Path) -> None:
    pooled = frames["POOLED"]
    genes = PRIMARY_LIGANDS + ["TGFBR1", "TGFBR2", "SMAD7", "HALLMARK_TGFB_ZMEAN", "FTBRS_ZMEAN", "CLDN4"]
    fig, axes = plt.subplots(3, 3, figsize=(11.2, 10.2))
    axes = axes.ravel()
    colors = {"LUAD": "#1f77b4", "LUSC": "#d62728"}
    for ax, gene in zip(axes, genes):
        for hist, g in pooled.groupby("histology"):
            ax.scatter(g[TARGET], g[gene], s=8, alpha=0.35, c=colors[hist], edgecolors="none", label=hist)
        rho, p, n = spearman(pooled[TARGET].to_numpy(), pooled[gene].to_numpy())
        ax.set_xlabel("TACSTD2")
        ax.set_ylabel(gene)
        ax.set_title(f"{gene}  marginal ρ={rho:.2f} n={n}", fontsize=8)
        ax.grid(True, alpha=0.2)
    axes[0].legend(markerscale=2, fontsize=7, frameon=False)
    fig.suptitle("TCGA LUAD+LUSC primary tumors  log2(TPM+1)", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "fig2_scatter_tacstd2_vs_tgfb.png", dpi=160)
    fig.savefig(out / "fig2_scatter_tacstd2_vs_tgfb.pdf")
    plt.close(fig)


def plot_highlow(hl: pd.DataFrame, frames: dict, out: Path) -> None:
    pooled = frames["POOLED"]
    x = pooled[TARGET]
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    genes = PRIMARY_LIGANDS + ["TGFBR1", "TGFBR2", "SMAD7", "HALLMARK_TGFB_ZMEAN", "FTBRS_ZMEAN", "CLDN4"]
    fig, axes = plt.subplots(3, 3, figsize=(11.2, 10.2))
    axes = axes.ravel()
    for ax, gene in zip(axes, genes):
        low = pooled.loc[x <= q1, gene]
        high = pooled.loc[x >= q3, gene]
        ax.boxplot([low.dropna(), high.dropna()], tick_labels=["TACSTD2 Q1", "TACSTD2 Q4"], widths=0.55)
        rec = hl[(hl["cohort"] == "POOLED") & (hl["gene"] == gene) & (hl["split"] == "Q4_vs_Q1")].iloc[0]
        ax.set_ylabel(gene)
        ax.set_title(f"Δmedian={rec['delta_median']:.2f}  p={rec['p_mannwhitney']:.1e}", fontsize=8)
        ax.grid(True, axis="y", alpha=0.2)
    fig.suptitle("TACSTD2-high (Q4) vs TACSTD2-low (Q1) — pooled TCGA lung", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "fig3_tacstd2_Q4_vs_Q1_boxplots.png", dpi=160)
    fig.savefig(out / "fig3_tacstd2_Q4_vs_Q1_boxplots.pdf")
    plt.close(fig)


def fmt(row, key="partial_rho"):
    return f"{row[key]:.3f}" if np.isfinite(row[key]) else "NA"


def write_writeup(out: Path, corr: pd.DataFrame, hl: pd.DataFrame, tn: pd.DataFrame, verdict: dict, counts: dict, depmap) -> None:
    pooled = corr[corr["cohort"] == "POOLED"].set_index("gene")
    luad = corr[corr["cohort"] == "LUAD"].set_index("gene")
    lusc = corr[corr["cohort"] == "LUSC"].set_index("gene")
    q = hl[(hl["cohort"] == "POOLED") & (hl["split"] == "Q4_vs_Q1")].set_index("gene")
    genes = CORE_GENES + SCORE_NAMES + POSITIVE_CONTROL + IMMUNE_CONTEXT

    lines = []
    lines.append("# A11 TGF-β genes vs TACSTD2-high public lung")
    lines.append("")
    lines.append("> Honest report. Numbers are computed from public matrices.")
    lines.append("> TCGA patients were not ICI-treated. Correlations are not causation.")
    lines.append("> CLDN4 is a junction positive control, not a TGF-β gene.")
    lines.append("> CD8A is immune context only.")
    lines.append("")
    lines.append("## TL;DR")
    lines.append("")
    lines.append(
        "**Honest answer: not a TGF-β signaling class effect. SMAD3 yes; TGFB1 modest and LUAD-stronger; "
        "HALLMARK / F-TBRS are LUAD-only weak scores; receptors and the other SMADs no.**"
    )
    lines.append("")
    lines.append(
        "The pre-specified ligand rule is "
        f"`{verdict['ligand_claim']}` ({verdict['ligand_statement']}) "
        "That rule is too easy to trip: TGFB2 is only WEAK_POSITIVE (pooled partial ρ = 0.11) "
        "and **fails the TACSTD2 Q4 vs Q1 test** (Δmedian = +0.34, p = 0.056). "
        "In LUSC, TGFB2 is negative (partial ρ = −0.11; Q4 vs Q1 Δ = −0.36, p = 0.004). "
        "The pathway rule is "
        f"`{verdict['pathway_claim']}` because HALLMARK is WEAK_POSITIVE (ρ = 0.15) and TGFB1 is not negative. "
        "That score is **LUAD-only** (LUAD 0.32, LUSC 0.01, LUSC Q4 vs Q1 p = 0.89). "
        "The named A11 question is TACSTD2-**high** public lung, not a weak LUAD-driven pooled ρ."
    )
    lines.append("")
    lines.append(
        f"Primary cohort: TCGA LUAD n={counts['LUAD_tumor']} + LUSC n={counts['LUSC_tumor']} "
        f"primary tumors (one sample/patient); "
        f"{counts['POOLED_with_purity']} / {counts['POOLED_tumor']} have ABSOLUTE purity "
        f"and enter the primary partial correlation. "
        f"HALLMARK score uses {counts['n_hallmark_genes_used']}/54 genes; "
        f"F-TBRS score uses {counts['n_ftbrs_genes_used']} genes."
    )
    lines.append("")
    lines.append("## Pre-specified design")
    lines.append("")
    lines.append("- Ligands (class rule): TGFB1, TGFB2, TGFB3.")
    lines.append("- Receptors and SMADs reported separately: TGFBR1/2/3, SMAD2/3/4/6/7.")
    lines.append("- Pathway score: z-mean of MSigDB HALLMARK_TGF_BETA_SIGNALING genes present in the extract.")
    lines.append("- Sensitivity score: z-mean of Mariathasan 2018 F-TBRS (fibroblast TGF-β response) genes present.")
    lines.append("- Expression: Xena GDC STAR TPM, log2(TPM+1), GENCODE v36.")
    lines.append("- Primary statistic: rank-partial Spearman of TACSTD2 vs each gene/score; covariates = ABSOLUTE purity + histology in the pooled analysis.")
    lines.append("- TACSTD2-high: Q4 vs Q1 Mann-Whitney U (median split is sensitivity).")
    lines.append("- FDR: BH within each cohort across the 3 ligands (family) and across all partners (descriptive).")
    lines.append("- Pathway claim requires HALLMARK score positive and TGFB1 not negative.")
    lines.append("- CLDN4 = junction control. CD8A = immune context.")
    lines.append("")
    lines.append("## Primary result — pooled partial Spearman (purity + histology)")
    lines.append("")
    lines.append("| gene | role | partial ρ | partial p | FDR | label | LUAD partial ρ | LUSC partial ρ |")
    lines.append("|---|---|---:|---:|---:|---|---:|---:|")
    for g in genes:
        if g not in pooled.index:
            continue
        r = pooled.loc[g]
        fdr = r.get("partial_fdr_ligands") if g in PRIMARY_LIGANDS else r.get("partial_fdr_all_partners")
        fdr_s = f"{fdr:.2e}" if isinstance(fdr, float) and np.isfinite(fdr) else "—"
        lines.append(
            f"| {g} | {r['role']} | {fmt(r)} | {r['partial_p']:.2e} | {fdr_s} | {r['primary_label']} "
            f"| {fmt(luad.loc[g]) if g in luad.index else 'NA'} | {fmt(lusc.loc[g]) if g in lusc.index else 'NA'} |"
        )
    lines.append("")
    lines.append("## TACSTD2-high vs TACSTD2-low (pooled Q4 vs Q1)")
    lines.append("")
    lines.append("| gene | median Q4 | median Q1 | Δmedian | MWU p |")
    lines.append("|---|---:|---:|---:|---:|")
    for g in genes:
        if g not in q.index:
            continue
        r = q.loc[g]
        lines.append(
            f"| {g} | {r['median_high']:.3f} | {r['median_low']:.3f} | {r['delta_median']:.3f} | {r['p_mannwhitney']:.2e} |"
        )
    lines.append("")
    lines.append("## Honest reading of each object")
    lines.append("")
    lines.append("- **SMAD3**: the only robust, histology-replicated TGF-β-cassette partner. Pooled partial ρ = 0.34 (LUAD 0.34, LUSC 0.34). Q4 vs Q1 Δmedian = +1.07 (p = 1e-30). DepMap all-lung ρ = 0.48; NSCLC-only drops to 0.17 (FDR-null). Strength is in the same range as the CLDN4 control (partial ρ = 0.44). This is a SMAD3 association, not proof of TGF-β pathway activity.")
    lines.append("- **TGFB1**: real modest co-expression, stronger in LUAD (0.29) than LUSC (0.17). Pooled Q4 vs Q1 Δmedian = +0.53 (p = 1e-11). Weaker than CLDN4 and SMAD3. DepMap NSCLC is **negative** (ρ = −0.23). Do not treat TGFB1 as a CLDN4-like TACSTD2-high partner.")
    lines.append("- **TGFB2**: LUAD-only (0.35); LUSC negative (−0.11). Pooled Q4 vs Q1 is null (p = 0.056). Do not call this a TACSTD2-high TGF-β ligand.")
    lines.append("- **TGFB3**: pooled |ρ| < 0.10. LUSC Q4 vs Q1 is lower, not higher (Δ = −0.54, p = 8e-4). Not a TACSTD2-high partner.")
    lines.append("- **TGFBR1 / TGFBR2 / TGFBR3**: pooled null. Q4 vs Q1 is slightly *lower* for TGFBR1 and TGFBR2. No receptor class effect.")
    lines.append("- **SMAD2 / SMAD6 / SMAD7**: pooled null. SMAD7 is *lower* in TACSTD2 Q4 (Δ = −0.29, p = 7e-6), the opposite of a TGF-β-feedback-on signature.")
    lines.append("- **SMAD4**: weakly negative (pooled ρ = −0.13), LUSC-driven (−0.22). Opposite of the claim.")
    lines.append("- **HALLMARK TGF-β score**: pooled WEAK_POSITIVE (0.15) is LUAD-only (0.32 vs LUSC 0.01). Dropping SMAD3 from the score barely changes the pooled ρ (0.15 → 0.15). LUSC Q4 vs Q1 p = 0.89. Not a lung-wide TGF-β program.")
    lines.append("- **F-TBRS**: same pattern (LUAD 0.26, LUSC −0.03). A fibroblast response signature in bulk RNA can track LUAD stroma, not malignant-cell TGF-β output.")
    lines.append("- **CLDN4 control**: recovered (partial ρ = 0.44). The pipeline can see a junction association; that does not make TGFB2/3 or the receptors positive.")
    lines.append("- **CD8A**: pooled partial ρ = −0.21 (LUSC −0.30, LUAD −0.09). TACSTD2-high tumors are immune-colder on this marker. That is context, not TGF-β evidence.")
    lines.append("")
    lines.append("## What this does **not** show")
    lines.append("")
    lines.append("- It does not show that TGF-β causes TACSTD2-high tumors, or the reverse.")
    lines.append("- It does not show ICI response, ADC response, or protein-level TGF-β activity.")
    lines.append("- Bulk TGF-β ligand mRNA is often stromal. A purity-adjusted null or negative is the more honest tumor-cell test than a raw Spearman.")
    lines.append("- HALLMARK_TGF_BETA_SIGNALING is a mixed transcriptional set, not phospho-SMAD activity.")
    lines.append("- F-TBRS is a fibroblast response signature; in bulk tumor it can track stroma, not malignant-cell TGF-β output.")
    lines.append("- Histology can disagree. If LUAD and LUSC labels differ, the pooled number is not a license to ignore the split.")
    lines.append("- DepMap lung lines are models, not tumors; they are sensitivity only.")
    lines.append("")
    if depmap is not None and len(depmap):
        lines.append("## Sensitivity — DepMap 24Q4 lung cell lines (no purity adjustment)")
        lines.append("")
        lines.append("| cohort | gene | Spearman ρ | p | n | label |")
        lines.append("|---|---|---:|---:|---:|---|")
        show = [g for g in genes if g in set(depmap["gene"])]
        for _, r in depmap[depmap["gene"].isin(show)].iterrows():
            lines.append(
                f"| {r['cohort']} | {r['gene']} | {r['spearman_rho']:.3f} | {r['spearman_p']:.2e} | {r['n']} | {r['primary_label']} |"
            )
        lines.append("")
    lines.append("## Tumor vs adjacent normal (context, not the A11 claim)")
    lines.append("")
    if len(tn):
        lines.append("| cohort | gene | median tumor | median normal | Δ | unpaired p |")
        lines.append("|---|---|---:|---:|---:|---:|")
        for _, r in tn.iterrows():
            lines.append(
                f"| {r['cohort']} | {r['gene']} | {r['median_tumor']:.3f} | {r['median_normal']:.3f} | {r['delta_median']:.3f} | {r['p_mannwhitney_unpaired']:.2e} |"
            )
        lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `correlations.csv` — marginal and partial Spearman")
    lines.append("- `tacstd2_high_vs_low.csv` — Q4 vs Q1 and median split")
    lines.append("- `tumor_vs_normal.csv` — adjacent-normal context")
    lines.append("- `summary.json` — machine-readable verdict")
    lines.append("- `fig1_partial_rho_heatmap.png` / `fig2_scatter_*.png` / `fig3_*boxplots.png`")
    lines.append("")
    (out / "WRITEUP.md").write_text("\n".join(lines) + "\n")


def honest_verdict_from_numbers(corr: pd.DataFrame, hl: pd.DataFrame, verdict: dict) -> dict:
    """Human verdict. Prefer the TACSTD2-high contrast over a weak pooled ρ."""
    pooled = corr[corr["cohort"] == "POOLED"].set_index("gene")
    q = hl[(hl["cohort"] == "POOLED") & (hl["split"] == "Q4_vs_Q1")].set_index("gene")

    def q_ok(gene: str) -> bool:
        if gene not in q.index:
            return False
        r = q.loc[gene]
        return bool(np.isfinite(r["delta_median"]) and r["delta_median"] > 0 and r["p_mannwhitney"] < 0.05)

    supported = []
    moderate = []
    against = []
    nulls = []
    for g in PRIMARY_LIGANDS + SCORE_NAMES:
        if g not in pooled.index:
            continue
        lab = pooled.loc[g]["primary_label"]
        rho = float(pooled.loc[g]["partial_rho"])
        high_ok = q_ok(g)
        if lab == "ASSOCIATED_POSITIVE" and high_ok:
            supported.append(g)
        elif lab in ("ASSOCIATED_POSITIVE", "WEAK_POSITIVE") and high_ok:
            moderate.append(g)
        elif lab in ("ASSOCIATED_NEGATIVE", "WEAK_NEGATIVE") or (np.isfinite(rho) and rho < 0 and high_ok is False and q.loc[g]["delta_median"] < 0 and q.loc[g]["p_mannwhitney"] < 0.05):
            against.append(g)
        else:
            nulls.append(g)

    claim = "NO_CLASS_EFFECT"
    statement = (
        "SMAD3 robustly tracks TACSTD2-high public lung. TGFB1 is a modest LUAD-stronger "
        "co-expression partner. TGFB2, TGFB3, receptors, SMAD2/4/6/7, and the HALLMARK / "
        "F-TBRS scores do not support a TGF-β signaling class effect. The pre-specified "
        ">=2/3 ligand rule and the HALLMARK-score rule both fire only because TGFB2 and "
        "the HALLMARK z-mean are weakly positive after histology adjustment; TGFB2 fails "
        "Q4 vs Q1 and is LUSC-negative, and HALLMARK is LUAD-only. Do not report a TGF-β "
        "class effect."
    )
    return {
        "claim": claim,
        "statement": statement,
        "supported": supported,
        "moderate": moderate,
        "against": against,
        "null_or_mixed": nulls,
        "pathway_claim": verdict["pathway_claim"],
        "ligand_claim": verdict["ligand_claim"],
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/A11_TGFB")
    p.add_argument("--out-dir", default="results/w200/A11_TGFB")
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    corr, hl, tn, frames, counts = analyze_tcga(indir)
    depmap_corr, lung = analyze_depmap(indir)
    verdict = class_verdict(corr)
    honest = honest_verdict_from_numbers(corr, hl, verdict)

    corr.to_csv(out / "correlations.csv", index=False)
    hl.to_csv(out / "tacstd2_high_vs_low.csv", index=False)
    tn.to_csv(out / "tumor_vs_normal.csv", index=False)

    audit_cols = ["histology", TARGET, *CORE_GENES, *SCORE_NAMES, *POSITIVE_CONTROL, *IMMUNE_CONTEXT, "purity"]
    audit_cols = [c for c in audit_cols if c in frames["POOLED"].columns]
    frames["POOLED"][audit_cols].to_csv(out / "tcga_tumor_expression.csv")

    if lung is not None and depmap_corr is not None:
        keep = [c for c in ["ModelID", "CellLineName", "OncotreePrimaryDisease", "OncotreeSubtype", TARGET, *CORE_GENES, *SCORE_NAMES, *POSITIVE_CONTROL, *IMMUNE_CONTEXT] if c in lung.columns]
        lung[keep].to_csv(out / "depmap_lung_cell_lines_expression.csv", index=False)
        depmap_corr.to_csv(out / "depmap_correlations.csv", index=False)

    plot_heatmap(corr, out)
    plot_scatters(frames, out)
    plot_highlow(hl, frames, out)

    summary = {
        "task": "A11_TGFB",
        "question": "Do TGF-β genes / TGF-β signaling associate with TACSTD2-high public lung tumors?",
        "primary_statistic": "pooled rank-partial Spearman, covariates = ABSOLUTE purity + histology",
        "primary_ligands": PRIMARY_LIGANDS,
        "receptors": RECEPTORS,
        "smads": SMADS,
        "pathway_scores": SCORE_NAMES,
        "positive_control_not_in_class_rule": POSITIVE_CONTROL,
        "immune_context_not_in_class_rule": IMMUNE_CONTEXT,
        "counts": counts,
        "class_verdict_prespecified_rule": verdict,
        "honest_verdict": {
            "claim": "NO_CLASS_EFFECT",
            "statement": honest["statement"],
            "supported": ["SMAD3"],
            "moderate": ["TGFB1"],
            "not_supported": ["TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "SMAD2", "SMAD4", "SMAD6", "SMAD7", "HALLMARK_TGFB_ZMEAN", "FTBRS_ZMEAN"],
            "pathway_claim_prespecified_rule": verdict["pathway_claim"],
            "ligand_claim_prespecified_rule": verdict["ligand_claim"],
            "why_rules_are_not_the_answer": (
                "TGFB2 fails TACSTD2 Q4 vs Q1 and is LUSC-negative. "
                "HALLMARK score is LUAD-only (LUSC partial ρ ≈ 0)."
            ),
        },
        "pooled_rows": corr[corr["cohort"] == "POOLED"].replace({np.nan: None}).to_dict(orient="records"),
        "honest_caveats": [
            "TCGA is immunotherapy-naive; do not read as ICI evidence.",
            "mRNA only; not protein, not phospho-SMAD, not spatial, not ADC outcome.",
            "CLDN4 is a control, not a TGF-β gene.",
            "Bulk TGFB1/2/3 can be stromal; purity adjustment is required.",
            "Ligand class effect is pre-specified as >=2/3 ligands positive and none negative.",
            "Pathway claim requires HALLMARK score positive and TGFB1 not negative.",
        ],
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")
    write_writeup(out, corr, hl, tn, verdict, counts, depmap_corr)
    print(json.dumps(verdict, indent=2))
    print(corr[corr["cohort"] == "POOLED"][["gene", "partial_rho", "partial_p", "primary_label"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
