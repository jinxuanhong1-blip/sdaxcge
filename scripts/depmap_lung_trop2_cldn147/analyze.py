#!/usr/bin/env python3
"""Lung-line TACSTD2 (TROP2) vs CLDN1, CLDN4, and CLDN7: protein and RNA.

Protein: Gygi/Nusinow CCLE MS (portal release name "Proteomics").
RNA: DepMap Public 24Q4 log2(TPM+1), the newest portal release with a public URL.

Slide comparison target, carried from the project CCLE protein-coexpression
slide (audited in PRs 54 and 111, not re-fit here):
  TACSTD2–CLDN4 CCLE protein, Spearman 0.69, n=118, labeled NSCLC.
No filter is chosen to hit 0.69.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

SLIDE_RHO = 0.69
SLIDE_N = 118
FOCAL = ("CLDN1", "CLDN4", "CLDN7")
CONTROLS = ("CLDN3", "CLDN5", "CLDN18")
PARTNERS = FOCAL
# Epithelial keratins used as the RNA adjustment set. Protein KRT18/KRT19
# are not in the Gygi matrix; KRT8 is a fragment with sparse lung coverage.
KERATIN = ("KRT8", "KRT18", "KRT19")
RNA_GENES = ("TACSTD2",) + FOCAL + CONTROLS + ("EPCAM",) + KERATIN
PROTEIN_REQUESTED = ("TACSTD2",) + FOCAL + CONTROLS + ("EPCAM", "KRT8", "KRT18", "KRT19")
PROTEIN_REQUIRED = ("TACSTD2", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "CLDN18", "EPCAM", "KRT8")
GENES = RNA_GENES
RNG_SEED = 0
N_BOOT = 5000


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    out = np.full(n, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out.tolist()
    idx = np.where(ok)[0]
    pv = p[idx]
    order = np.argsort(pv)
    m = len(pv)
    ranked = pv[order]
    q = ranked * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    placed = np.empty(m)
    placed[order] = q
    out[idx] = placed
    return out.tolist()


def corr_block(x, y, ci: bool = False) -> dict:
    x = pd.to_numeric(pd.Series(x), errors="coerce")
    y = pd.to_numeric(pd.Series(y), errors="coerce")
    m = x.notna() & y.notna()
    xv = x[m].to_numpy(float)
    yv = y[m].to_numpy(float)
    n = int(len(xv))
    out = {
        "n": n,
        "spearman_rho": None,
        "spearman_p": None,
        "pearson_r": None,
        "pearson_p": None,
        "spearman_ci95_low": None,
        "spearman_ci95_high": None,
        "rounds_to_slide_0.69": False,
        "abs_diff_vs_slide_0.69": None,
        "n_equals_slide_118": n == SLIDE_N,
    }
    if n < 5 or np.unique(xv).size < 2 or np.unique(yv).size < 2:
        return out
    rho, p_s = stats.spearmanr(xv, yv)
    r, p_p = stats.pearsonr(xv, yv)
    out.update(
        {
            "spearman_rho": float(rho),
            "spearman_p": float(p_s),
            "pearson_r": float(r),
            "pearson_p": float(p_p),
            "rounds_to_slide_0.69": float(f"{float(rho):.2f}") == SLIDE_RHO,
            "abs_diff_vs_slide_0.69": abs(float(rho) - SLIDE_RHO),
        }
    )
    if ci and n >= 8:
        rng = np.random.default_rng(RNG_SEED)
        boots = np.empty(N_BOOT)
        for i in range(N_BOOT):
            take = rng.integers(0, n, n)
            boots[i] = stats.spearmanr(xv[take], yv[take]).statistic
        lo, hi = np.percentile(boots, [2.5, 97.5])
        out["spearman_ci95_low"] = float(lo)
        out["spearman_ci95_high"] = float(hi)
        out["slide_0.69_inside_ci95"] = bool(lo <= SLIDE_RHO <= hi)
    return out


def load_protein(gz_path: Path) -> pd.DataFrame:
    with gzip.open(gz_path, "rt", newline="") as f:
        reader = csv.DictReader(f)
        quant = [c for c in reader.fieldnames if "_TenPx" in c and not c.endswith("_Peptides")]
        wanted = {}
        for row in reader:
            sym = row["Gene_Symbol"]
            if sym in PROTEIN_REQUESTED:
                if sym in wanted:
                    raise SystemExit(f"Duplicate Gygi gene symbol {sym}")
                wanted[sym] = row
    missing_required = [g for g in PROTEIN_REQUIRED if g not in wanted]
    if missing_required:
        raise SystemExit(f"Gygi matrix missing required genes {missing_required}")
    absent = [g for g in PROTEIN_REQUESTED if g not in wanted]
    frame = pd.DataFrame({"col": quant})
    frame["ccle"] = frame["col"].str.rsplit("_TenPx", n=1).str[0]
    for g in PROTEIN_REQUESTED:
        if g not in wanted:
            frame[g] = np.nan
            continue
        vals = []
        for c in quant:
            raw = wanted[g][c]
            vals.append(np.nan if raw in ("", "NA", "NaN", "nan") else float(raw))
        frame[g] = vals
    frame.attrs["ids"] = {g: wanted[g]["Protein_Id"] for g in wanted}
    frame.attrs["absent"] = absent
    return frame


def partial_spearman(y, x, covariates: list[pd.Series]) -> dict:
    """Spearman partial correlation: rank, then Pearson residual correlation.

    p uses df = n - 2 - k. This is not the slide Spearman and is not compared
    to 0.69 as if it were the same statistic.
    """
    pieces = [pd.to_numeric(pd.Series(y), errors="coerce").rename("y"), pd.to_numeric(pd.Series(x), errors="coerce").rename("x")]
    for i, cov in enumerate(covariates):
        pieces.append(pd.to_numeric(pd.Series(cov), errors="coerce").rename(f"c{i}"))
    df = pd.concat(pieces, axis=1).dropna()
    n = int(len(df))
    k = len(covariates)
    out = {
        "n": n,
        "partial_spearman": None,
        "partial_p": None,
        "unadjusted_spearman_same_n": None,
        "unadjusted_spearman_p_same_n": None,
    }
    if n < k + 8:
        return out
    unadj = corr_block(df["y"], df["x"], ci=False)
    out["unadjusted_spearman_same_n"] = unadj["spearman_rho"]
    out["unadjusted_spearman_p_same_n"] = unadj["spearman_p"]
    ranked = df.rank(method="average")
    z = np.column_stack([np.ones(n), ranked.iloc[:, 2:].to_numpy(float)])
    yv = ranked["y"].to_numpy(float)
    xv = ranked["x"].to_numpy(float)
    beta_y, _, _, _ = np.linalg.lstsq(z, yv, rcond=None)
    beta_x, _, _, _ = np.linalg.lstsq(z, xv, rcond=None)
    ry = yv - z @ beta_y
    rx = xv - z @ beta_x
    if np.unique(ry).size < 2 or np.unique(rx).size < 2:
        return out
    r = float(np.corrcoef(ry, rx)[0, 1])
    df_res = n - 2 - k
    if df_res <= 0 or not np.isfinite(r) or abs(r) >= 1:
        p = None
    else:
        tstat = r * np.sqrt(df_res / (1.0 - r * r))
        p = float(2 * stats.t.sf(abs(tstat), df_res))
    out["partial_spearman"] = r
    out["partial_p"] = p
    return out


def disease_group(primary: str) -> str:
    if primary == "Non-Small Cell Lung Cancer":
        return "NSCLC"
    if primary == "Lung Neuroendocrine Tumor":
        return "NET"
    if isinstance(primary, str) and primary and primary != "nan":
        return "other"
    return "unmapped"


def load_model(path: Path) -> pd.DataFrame:
    m = pd.read_csv(path)
    keep = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "ModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
        "CCLEName",
    ]
    m = m[keep].copy()
    m["CCLEName"] = m["CCLEName"].astype(str).replace({"nan": np.nan})
    return m


def add_row(rows: list[dict], **kwargs) -> dict:
    rec = dict(kwargs)
    rows.append(rec)
    return rec


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/depmap_trop2_cldn147")
    p.add_argument("--out-dir", default="results/depmap_lung_trop2_cldn147")
    p.add_argument("--protein-gz", default="")
    p.add_argument("--model-csv", default="")
    p.add_argument("--rna-csv", default="")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    protein_gz = Path(args.protein_gz) if args.protein_gz else cache / "protein_quant_current_normalized.csv.gz"
    model_csv = Path(args.model_csv) if args.model_csv else cache / "Model.csv"
    rna_csv = (
        Path(args.rna_csv)
        if args.rna_csv
        else out / "depmap24q4_selected_genes_all_models.csv"
    )
    for path in (protein_gz, model_csv, rna_csv):
        if not path.exists():
            raise SystemExit(f"Missing {path}. Run download.py first.")

    model = load_model(model_csv)
    dup_ccle = int(model["CCLEName"].dropna().duplicated().sum())
    model_ccle = model.dropna(subset=["CCLEName"]).drop_duplicates("CCLEName", keep="first")

    prot = load_protein(protein_gz)
    protein_ids = dict(prot.attrs["ids"])
    protein_absent = list(prot.attrs["absent"])
    prot = prot.merge(
        model_ccle.drop(columns=["ModelID"]).rename(columns={"CCLEName": "ccle"}),
        on="ccle",
        how="left",
    )
    prot["disease_group"] = prot["OncotreePrimaryDisease"].map(
        lambda v: disease_group("" if pd.isna(v) else str(v))
    )
    prot["is_lung_suffix"] = prot["ccle"].str.endswith("_LUNG")
    prot.to_csv(out / "gygi_all_lines_selected_genes.csv", index=False)

    rna = pd.read_csv(rna_csv)
    for g in GENES:
        rna[g] = pd.to_numeric(rna[g], errors="coerce")
    rna = rna.merge(model, on="ModelID", how="left")
    rna["disease_group"] = rna["OncotreePrimaryDisease"].map(
        lambda v: disease_group("" if pd.isna(v) else str(v))
    )
    lung_rna_mask = (rna["OncotreeLineage"] == "Lung") & (rna["ModelType"] == "Cell Line")
    lung_rna = rna.loc[lung_rna_mask].copy()
    lung_rna.to_csv(out / "depmap24q4_lung_cell_lines_selected_genes.csv", index=False)

    rows: list[dict] = []

    def protein_cohort(name: str, mask: pd.Series, note: str, primary: bool, partners, role: str) -> None:
        sub = prot.loc[mask]
        for partner in partners:
            block = corr_block(sub["TACSTD2"], sub[partner], ci=primary)
            extra = ""
            if partner in protein_absent:
                extra = " Not quantified in the Gygi CCLE MS matrix."
            elif partner == "CLDN18":
                extra = " Protein row is isoform A2 only (sp|P56856-2|CLD18_HUMAN)."
            add_row(
                rows,
                layer="protein",
                dataset="Gygi_CCLE_MS",
                cohort=name,
                partner=partner,
                primary=primary,
                role=role,
                note=note + extra,
                **block,
            )

    def rna_cohort(name: str, mask: pd.Series, note: str, primary: bool, partners, role: str) -> None:
        sub = rna.loc[mask]
        for partner in partners:
            block = corr_block(sub["TACSTD2"], sub[partner], ci=primary)
            add_row(
                rows,
                layer="RNA",
                dataset="DepMap_Public_24Q4_log2TPM",
                cohort=name,
                partner=partner,
                primary=primary,
                role=role,
                note=note,
                **block,
            )

    protein_cohort(
        "lung_suffix",
        prot["is_lung_suffix"],
        "PRIMARY protein: CCLE name ends with _LUNG. Pairwise complete cases. No imputation.",
        True,
        FOCAL,
        "focal",
    )
    protein_cohort(
        "lung_suffix",
        prot["is_lung_suffix"],
        "Specificity controls on the same _LUNG columns. Not the slide pair.",
        False,
        CONTROLS,
        "control",
    )
    protein_cohort(
        "lung_suffix_NSCLC",
        prot["is_lung_suffix"] & (prot["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"),
        "CCLE _LUNG and DepMap 24Q4 OncotreePrimaryDisease == Non-Small Cell Lung Cancer.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    protein_cohort(
        "lung_suffix_NET",
        prot["is_lung_suffix"] & (prot["OncotreePrimaryDisease"] == "Lung Neuroendocrine Tumor"),
        "CCLE _LUNG and OncotreePrimaryDisease == Lung Neuroendocrine Tumor (includes SCLC).",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    protein_cohort(
        "lung_suffix_LUAD",
        prot["is_lung_suffix"] & (prot["OncotreeSubtype"] == "Lung Adenocarcinoma"),
        "CCLE _LUNG and OncotreeSubtype == Lung Adenocarcinoma.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    protein_cohort(
        "lung_suffix_LUSC",
        prot["is_lung_suffix"] & (prot["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma"),
        "CCLE _LUNG and OncotreeSubtype == Lung Squamous Cell Carcinoma.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    protein_cohort(
        "all_lines",
        prot["TACSTD2"].notna(),
        "All Gygi lines. Context only, not the lung slide.",
        False,
        FOCAL + CONTROLS,
        "context",
    )

    rna_cohort(
        "lung_cell_lines",
        lung_rna_mask,
        "PRIMARY RNA: OncotreeLineage==Lung and ModelType==Cell Line.",
        True,
        FOCAL,
        "focal",
    )
    rna_cohort(
        "lung_cell_lines",
        lung_rna_mask,
        "Specificity controls in the same lung cell lines. Not the slide pair.",
        False,
        CONTROLS,
        "control",
    )
    rna_cohort(
        "lung_NSCLC_cell_lines",
        lung_rna_mask & (rna["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"),
        "Lung cell lines, OncotreePrimaryDisease == Non-Small Cell Lung Cancer.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    rna_cohort(
        "lung_NET_cell_lines",
        lung_rna_mask & (rna["OncotreePrimaryDisease"] == "Lung Neuroendocrine Tumor"),
        "Lung cell lines, OncotreePrimaryDisease == Lung Neuroendocrine Tumor.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    rna_cohort(
        "lung_LUAD_cell_lines",
        lung_rna_mask & (rna["OncotreeSubtype"] == "Lung Adenocarcinoma"),
        "Lung cell lines, OncotreeSubtype == Lung Adenocarcinoma.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )
    rna_cohort(
        "lung_LUSC_cell_lines",
        lung_rna_mask & (rna["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma"),
        "Lung cell lines, OncotreeSubtype == Lung Squamous Cell Carcinoma.",
        False,
        FOCAL + CONTROLS,
        "histology",
    )

    # Same CCLE names: RNA on the Gygi _LUNG lines that are in the 24Q4 matrix.
    gygi_lung_names = set(prot.loc[prot["is_lung_suffix"], "ccle"])
    rna_on_gygi = lung_rna_mask & rna["CCLEName"].isin(gygi_lung_names)
    rna_cohort(
        "gygi_lung_lines_with_rna",
        rna_on_gygi,
        "DepMap 24Q4 RNA restricted to CCLE names present in the Gygi _LUNG protein columns.",
        False,
        FOCAL,
        "sensitivity",
    )

    # Same lines, both layers quantified for that pair.
    prot_lung = prot.loc[prot["is_lung_suffix"], ["ccle", *GENES]].rename(
        columns={g: f"{g}_protein" for g in GENES}
    )
    rna_by_ccle = (
        rna.loc[lung_rna_mask, ["CCLEName", *GENES]]
        .dropna(subset=["CCLEName"])
        .drop_duplicates("CCLEName", keep="first")
        .rename(columns={"CCLEName": "ccle", **{g: f"{g}_rna" for g in GENES}})
    )
    matched = prot_lung.merge(rna_by_ccle, on="ccle", how="inner")
    matched.to_csv(out / "gygi_lung_with_depmap24q4_rna.csv", index=False)
    for partner in PARTNERS:
        both = matched["TACSTD2_protein"].notna() & matched[f"{partner}_protein"].notna()
        both = both & matched["TACSTD2_rna"].notna() & matched[f"{partner}_rna"].notna()
        sub = matched.loc[both]
        block_p = corr_block(sub["TACSTD2_protein"], sub[f"{partner}_protein"], ci=False)
        block_r = corr_block(sub["TACSTD2_rna"], sub[f"{partner}_rna"], ci=False)
        note = (
            "Same CCLE names: Gygi _LUNG lines with both proteins quantified and DepMap 24Q4 RNA. "
            f"n_gygi_lung_with_any_rna={len(matched)}."
        )
        add_row(
            rows,
            layer="protein",
            dataset="Gygi_CCLE_MS",
            cohort="matched_rna_protein_lines",
            partner=partner,
            primary=False,
            role="sensitivity",
            note=note,
            **block_p,
        )
        add_row(
            rows,
            layer="RNA",
            dataset="DepMap_Public_24Q4_log2TPM",
            cohort="matched_rna_protein_lines",
            partner=partner,
            primary=False,
            role="sensitivity",
            note=note,
            **block_r,
        )

    primary_idx = [i for i, r in enumerate(rows) if r["primary"]]
    qvals = bh([rows[i]["spearman_p"] if rows[i]["spearman_p"] is not None else np.nan for i in primary_idx])
    for i, q in zip(primary_idx, qvals):
        rows[i]["spearman_q_bh_primary6"] = q
    for r in rows:
        r.setdefault("spearman_q_bh_primary6", None)
        r.setdefault("slide_0.69_inside_ci95", None)

    corr_df = pd.DataFrame(rows)
    col_order = [
        "layer",
        "dataset",
        "cohort",
        "partner",
        "role",
        "primary",
        "n",
        "spearman_rho",
        "spearman_p",
        "spearman_q_bh_primary6",
        "spearman_ci95_low",
        "spearman_ci95_high",
        "pearson_r",
        "pearson_p",
        "rounds_to_slide_0.69",
        "abs_diff_vs_slide_0.69",
        "n_equals_slide_118",
        "slide_0.69_inside_ci95",
        "note",
    ]
    corr_df = corr_df[col_order]
    corr_df.to_csv(out / "correlations.csv", index=False)

    # Slide comparison is the lung primary protein and RNA rows.
    slide = corr_df[corr_df["primary"]].copy()
    slide.insert(0, "slide_claim", "CCLE protein coexpression TACSTD2–CLDN4 Spearman 0.69, n=118, labeled NSCLC")
    slide.to_csv(out / "slide_comparison.csv", index=False)

    lung_prot = prot.loc[prot["is_lung_suffix"]].copy()
    lung_prot.to_csv(out / "gygi_lung_selected_genes.csv", index=False)

    # Partial Spearman. Slide ρ=0.69 remains the unadjusted protein CLDN4 Spearman.
    partial_rows: list[dict] = []
    rna_adjust = [
        ("EPCAM", ["EPCAM"]),
        ("KRT8+KRT18+KRT19", list(KERATIN)),
        ("EPCAM+KRT8+KRT18+KRT19", ["EPCAM", *KERATIN]),
    ]
    protein_adjust = [("EPCAM", ["EPCAM"])]
    partial_partners = FOCAL + CONTROLS

    def add_partials(frame, layer, dataset, cohort, mask, adjustments, family: str) -> None:
        sub = frame.loc[mask]
        for partner in partial_partners:
            if partner not in sub.columns:
                continue
            for adj_name, covars in adjustments:
                if any(c not in sub.columns for c in covars):
                    continue
                if sub[covars].notna().any().sum() < len(covars):
                    partial_rows.append(
                        {
                            "layer": layer,
                            "dataset": dataset,
                            "cohort": cohort,
                            "partner": partner,
                            "adjustment": adj_name,
                            "family": family,
                            "n": 0,
                            "partial_spearman": None,
                            "partial_p": None,
                            "unadjusted_spearman_same_n": None,
                            "unadjusted_spearman_p_same_n": None,
                            "note": "Covariate not quantified in this matrix.",
                        }
                    )
                    continue
                block = partial_spearman(sub["TACSTD2"], sub[partner], [sub[c] for c in covars])
                note = (
                    "Spearman partial correlation (rank, then residual Pearson). "
                    "Not the slide statistic. Slide remains unadjusted protein TACSTD2–CLDN4 Spearman."
                )
                if layer == "protein" and partner == "CLDN18":
                    note += " CLDN18 protein is isoform A2 only."
                partial_rows.append(
                    {
                        "layer": layer,
                        "dataset": dataset,
                        "cohort": cohort,
                        "partner": partner,
                        "adjustment": adj_name,
                        "family": family,
                        "note": note,
                        **block,
                    }
                )

    add_partials(prot, "protein", "Gygi_CCLE_MS", "lung_suffix", prot["is_lung_suffix"], protein_adjust, "primary_partial")
    add_partials(
        prot,
        "protein",
        "Gygi_CCLE_MS",
        "lung_suffix_NSCLC",
        prot["is_lung_suffix"] & (prot["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"),
        protein_adjust,
        "histology_partial",
    )
    add_partials(
        prot,
        "protein",
        "Gygi_CCLE_MS",
        "lung_suffix_LUAD",
        prot["is_lung_suffix"] & (prot["OncotreeSubtype"] == "Lung Adenocarcinoma"),
        protein_adjust,
        "histology_partial",
    )
    add_partials(
        prot,
        "protein",
        "Gygi_CCLE_MS",
        "lung_suffix_LUSC",
        prot["is_lung_suffix"] & (prot["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma"),
        protein_adjust,
        "histology_partial",
    )
    add_partials(rna, "RNA", "DepMap_Public_24Q4_log2TPM", "lung_cell_lines", lung_rna_mask, rna_adjust, "primary_partial")
    add_partials(
        rna,
        "RNA",
        "DepMap_Public_24Q4_log2TPM",
        "lung_NSCLC_cell_lines",
        lung_rna_mask & (rna["OncotreePrimaryDisease"] == "Non-Small Cell Lung Cancer"),
        rna_adjust,
        "histology_partial",
    )
    add_partials(
        rna,
        "RNA",
        "DepMap_Public_24Q4_log2TPM",
        "lung_LUAD_cell_lines",
        lung_rna_mask & (rna["OncotreeSubtype"] == "Lung Adenocarcinoma"),
        rna_adjust,
        "histology_partial",
    )
    add_partials(
        rna,
        "RNA",
        "DepMap_Public_24Q4_log2TPM",
        "lung_LUSC_cell_lines",
        lung_rna_mask & (rna["OncotreeSubtype"] == "Lung Squamous Cell Carcinoma"),
        rna_adjust,
        "histology_partial",
    )
    # Protein keratin adjustment is not estimable: KRT18 and KRT19 have no Gygi row,
    # and KRT8 is a fragment. Record the coverage instead of a partial rho.
    krt8_lung = int(prot.loc[prot["is_lung_suffix"], "KRT8"].notna().sum())
    cldn4_and_krt8 = int(
        (prot["is_lung_suffix"] & prot["CLDN4"].notna() & prot["KRT8"].notna()).sum()
    )
    partial_rows.append(
        {
            "layer": "protein",
            "dataset": "Gygi_CCLE_MS",
            "cohort": "lung_suffix",
            "partner": "CLDN4",
            "adjustment": "KRT8+KRT18+KRT19",
            "family": "not_estimable",
            "n": cldn4_and_krt8,
            "partial_spearman": None,
            "partial_p": None,
            "unadjusted_spearman_same_n": None,
            "unadjusted_spearman_p_same_n": None,
            "note": (
                f"Not estimated. KRT18 and KRT19 are absent from Gygi. "
                f"KRT8 fragment is quantified in {krt8_lung}/77 _LUNG columns and in "
                f"{cldn4_and_krt8} lines that also have CLDN4. No protein keratin partial is reported."
            ),
        }
    )
    partial_df = pd.DataFrame(partial_rows)
    bh_idx = [
        i
        for i, r in enumerate(partial_rows)
        if r["family"] == "primary_partial" and r.get("partial_p") is not None and r["partner"] in ("CLDN4", "CLDN7")
    ]
    q_partial = bh([partial_rows[i]["partial_p"] for i in bh_idx])
    partial_df["partial_q_bh_cldn4_cldn7"] = np.nan
    for i, q in zip(bh_idx, q_partial):
        partial_df.loc[i, "partial_q_bh_cldn4_cldn7"] = q
    partial_df.to_csv(out / "partial_correlations.csv", index=False)

    # Figures
    disease_colors = {
        "NSCLC": "#2166ac",
        "NET": "#b2182b",
        "other": "#7f7f7f",
        "unmapped": "#c7c7c7",
    }
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 7.2))
    for col, partner in enumerate(PARTNERS):
        # protein
        ax = axes[0, col]
        sub = lung_prot.dropna(subset=["TACSTD2", partner])
        for grp, g in sub.groupby(sub["disease_group"].fillna("unmapped")):
            ax.scatter(
                g["TACSTD2"],
                g[partner],
                s=26,
                alpha=0.85,
                c=disease_colors.get(grp, "#333333"),
                label=f"{grp} n={len(g)}",
                linewidths=0,
            )
        block = corr_df[
            (corr_df.layer == "protein") & (corr_df.cohort == "lung_suffix") & (corr_df.partner == partner)
        ].iloc[0]
        ax.set_title(
            f"Protein  TACSTD2–{partner}\nSpearman ρ={block.spearman_rho:.3f}  n={int(block.n)}  p={block.spearman_p:.2e}",
            fontsize=10,
        )
        ax.set_xlabel("TACSTD2 protein (Gygi TMT)")
        ax.set_ylabel(f"{partner} protein")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, frameon=False, loc="best")

        ax = axes[1, col]
        sub = lung_rna.dropna(subset=["TACSTD2", partner])
        for grp, g in sub.groupby(sub["disease_group"].fillna("unmapped")):
            ax.scatter(
                g["TACSTD2"],
                g[partner],
                s=18,
                alpha=0.75,
                c=disease_colors.get(grp, "#333333"),
                label=f"{grp} n={len(g)}",
                linewidths=0,
            )
        block = corr_df[
            (corr_df.layer == "RNA") & (corr_df.cohort == "lung_cell_lines") & (corr_df.partner == partner)
        ].iloc[0]
        ax.set_title(
            f"RNA  TACSTD2–{partner}\nSpearman ρ={block.spearman_rho:.3f}  n={int(block.n)}  p={block.spearman_p:.2e}",
            fontsize=10,
        )
        ax.set_xlabel("TACSTD2 RNA, log2(TPM+1)")
        ax.set_ylabel(f"{partner} RNA, log2(TPM+1)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, frameon=False, loc="best")
    fig.suptitle(
        "Lung lines: TROP2 (TACSTD2) vs CLDN1 / CLDN4 / CLDN7\n"
        "Top: Gygi CCLE MS, CCLE name ends with _LUNG. Bottom: DepMap 24Q4 lung cell lines.",
        fontsize=12,
    )
    fig.tight_layout()
    fig.savefig(out / "fig_scatter_protein_rna.png", dpi=160)
    fig.savefig(out / "fig_scatter_protein_rna.pdf")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    prim = corr_df[corr_df["primary"]].copy()
    x = np.arange(len(PARTNERS))
    width = 0.36
    prot_y = [
        float(prim[(prim.layer == "protein") & (prim.partner == g)]["spearman_rho"].iloc[0]) for g in PARTNERS
    ]
    rna_y = [float(prim[(prim.layer == "RNA") & (prim.partner == g)]["spearman_rho"].iloc[0]) for g in PARTNERS]
    prot_n = [int(prim[(prim.layer == "protein") & (prim.partner == g)]["n"].iloc[0]) for g in PARTNERS]
    rna_n = [int(prim[(prim.layer == "RNA") & (prim.partner == g)]["n"].iloc[0]) for g in PARTNERS]
    b1 = ax.bar(x - width / 2, prot_y, width, color="#2166ac", label="Protein, Gygi _LUNG")
    b2 = ax.bar(x + width / 2, rna_y, width, color="#d68910", label="RNA, DepMap 24Q4 lung cell lines")
    ax.axhline(SLIDE_RHO, color="#333333", ls="--", lw=1, label="Slide protein claim ρ=0.69")
    ax.set_xticks(x)
    ax.set_xticklabels(PARTNERS)
    ax.set_ylabel("Spearman ρ with TACSTD2")
    ax.set_ylim(0, 1.05)
    ax.set_title("Lung-line coexpression vs the CCLE protein slide (ρ=0.69, n=118)")
    for bar, n in zip(b1, prot_n):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f"n={n}", ha="center", va="bottom", fontsize=8)
    for bar, n in zip(b2, rna_n):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f"n={n}", ha="center", va="bottom", fontsize=8)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out / "fig_rho_vs_slide.png", dpi=160)
    fig.savefig(out / "fig_rho_vs_slide.pdf")
    plt.close(fig)

    show_partners = ("CLDN4", "CLDN7", "CLDN1", "CLDN3", "CLDN5", "CLDN18")
    hist_specs = [
        ("protein", [("lung_suffix", "All lung"), ("lung_suffix_NSCLC", "NSCLC"), ("lung_suffix_LUAD", "LUAD"), ("lung_suffix_LUSC", "LUSC")]),
        ("RNA", [("lung_cell_lines", "All lung"), ("lung_NSCLC_cell_lines", "NSCLC"), ("lung_LUAD_cell_lines", "LUAD"), ("lung_LUSC_cell_lines", "LUSC")]),
    ]
    hist_colors = ["#2166ac", "#d68910", "#1a9850", "#762a83"]
    fig, axes = plt.subplots(2, 1, figsize=(11.4, 7.4), sharey=False)
    for ax, (layer, cohorts) in zip(axes, hist_specs):
        x = np.arange(len(show_partners))
        width = 0.18
        for j, ((cohort, label), color) in enumerate(zip(cohorts, hist_colors)):
            ys = []
            ns = []
            for partner in show_partners:
                hit = corr_df[(corr_df.layer == layer) & (corr_df.cohort == cohort) & (corr_df.partner == partner)]
                if hit.empty or pd.isna(hit.iloc[0].spearman_rho):
                    ys.append(np.nan)
                    ns.append(0)
                else:
                    ys.append(float(hit.iloc[0].spearman_rho))
                    ns.append(int(hit.iloc[0].n))
            xpos = x + (j - 1.5) * width
            bars = ax.bar(xpos, ys, width, color=color, label=label)
            for bar, n, y in zip(bars, ns, ys):
                if n and np.isfinite(y):
                    ax.text(bar.get_x() + bar.get_width() / 2, y + (0.03 if y >= 0 else -0.08), f"{n}", ha="center", va="bottom" if y >= 0 else "top", fontsize=6, color="#333333")
        ax.axhline(0, color="#666666", lw=0.6)
        if layer == "protein":
            ax.axhline(SLIDE_RHO, color="#333333", ls="--", lw=0.8, label="Slide ρ=0.69 (unadjusted protein CLDN4)")
        ax.set_xticks(x)
        ax.set_xticklabels(show_partners)
        ax.set_ylim(-0.55, 1.18)
        ax.set_ylabel("Spearman ρ with TACSTD2")
        ax.set_title(layer)
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper right")
    fig.suptitle("All-lung vs NSCLC vs LUAD vs LUSC. Numbers on bars are n. CLDN5 protein is absent from Gygi.", fontsize=11)
    fig.tight_layout()
    fig.savefig(out / "fig_histology_controls.png", dpi=160)
    fig.savefig(out / "fig_histology_controls.pdf")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharey=True)
    partial_show = ("CLDN4", "CLDN7", "CLDN1", "CLDN3")
    # Left: protein unadjusted vs EPCAM. Right: RNA unadjusted vs three adjustments.
    # Unadjusted values are the cohort Spearmans, not partial correlations.
    prot_adj_order = [("unadjusted", None), ("EPCAM", "EPCAM")]
    rna_adj_order = [
        ("unadjusted", None),
        ("EPCAM", "EPCAM"),
        ("keratin", "KRT8+KRT18+KRT19"),
        ("both", "EPCAM+KRT8+KRT18+KRT19"),
    ]
    panel_colors = ["#2166ac", "#d68910", "#1a9850", "#762a83"]
    for ax, layer, cohort, order in (
        (axes[0], "protein", "lung_suffix", prot_adj_order),
        (axes[1], "RNA", "lung_cell_lines", rna_adj_order),
    ):
        x = np.arange(len(partial_show))
        width = 0.18 if layer == "RNA" else 0.32
        for j, ((label, adj), color) in enumerate(zip(order, panel_colors)):
            ys = []
            for partner in partial_show:
                if adj is None:
                    hit = corr_df[(corr_df.layer == layer) & (corr_df.cohort == cohort) & (corr_df.partner == partner)]
                    ys.append(float(hit.iloc[0].spearman_rho) if not hit.empty and pd.notna(hit.iloc[0].spearman_rho) else np.nan)
                else:
                    hit = partial_df[
                        (partial_df.layer == layer)
                        & (partial_df.cohort == cohort)
                        & (partial_df.partner == partner)
                        & (partial_df.adjustment == adj)
                    ]
                    ys.append(float(hit.iloc[0].partial_spearman) if not hit.empty and pd.notna(hit.iloc[0].partial_spearman) else np.nan)
            ax.bar(x + (j - (len(order) - 1) / 2) * width, ys, width, color=color, label=label)
        ax.axhline(0, color="#666666", lw=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(partial_show)
        ax.set_ylim(-0.3, 1.05)
        ax.set_title(f"{layer}, {cohort}")
        ax.grid(True, axis="y", alpha=0.25)
        ax.legend(frameon=False, fontsize=7)
        ax.set_ylabel("Correlation with TACSTD2")
    fig.suptitle("Unadjusted Spearman vs partial Spearman. The slide ρ=0.69 is the unadjusted protein CLDN4 bar, not a partial.", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "fig_partial_keratin_epcam.png", dpi=160)
    fig.savefig(out / "fig_partial_keratin_epcam.pdf")
    plt.close(fig)

    # Reproduction check against the prior Gygi lung CLDN4 complete-case result.
    cldn4 = corr_df[
        (corr_df.layer == "protein") & (corr_df.cohort == "lung_suffix") & (corr_df.partner == "CLDN4")
    ].iloc[0]
    repro_ok = int(cldn4.n) == 45 and abs(float(cldn4.spearman_rho) - 0.6931488801054019) < 1e-9

    n_lung_prot = int(prot["is_lung_suffix"].sum())
    n_lung_prot_mapped = int((prot["is_lung_suffix"] & prot["CellLineName"].notna()).sum())
    summary = {
        "slide_claim": {
            "text": "CCLE protein coexpression of TROP2 (TACSTD2) with CLDN4: Spearman 0.69, n=118, labeled NSCLC",
            "rho": SLIDE_RHO,
            "n": SLIDE_N,
            "source_note": "Project slide number previously traced in PRs 54 and 111. This script does not tune filters to that number.",
        },
        "gygi_protein_ids": protein_ids,
        "n_gygi_columns": int(len(prot)),
        "n_gygi_lung_suffix": n_lung_prot,
        "n_gygi_lung_mapped_to_depmap24q4": n_lung_prot_mapped,
        "n_duplicate_ccle_names_dropped": dup_ccle,
        "n_rna_models": int(len(rna)),
        "n_rna_lung_cell_lines": int(lung_rna_mask.sum()),
        "cldn4_lung_protein_reproduces_prior_0.693_n45": bool(repro_ok),
        "protein_genes_absent": protein_absent,
        "protein_krt8_lung_n": krt8_lung,
        "protein_cldn4_and_krt8_n": cldn4_and_krt8,
        "partial_note": "Partial correlations are not the slide statistic. Slide reconciliation uses unadjusted protein TACSTD2–CLDN4 Spearman only.",
        "primary": slide[
            [
                "layer",
                "cohort",
                "partner",
                "n",
                "spearman_rho",
                "spearman_p",
                "spearman_q_bh_primary6",
                "spearman_ci95_low",
                "spearman_ci95_high",
                "pearson_r",
                "rounds_to_slide_0.69",
                "n_equals_slide_118",
            ]
        ].to_dict(orient="records"),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        "SLIDE: CCLE protein coexpression TACSTD2-CLDN4 Spearman 0.69, n=118, labeled NSCLC.",
        f"CLDN4 protein _LUNG reproduces prior Gygi complete-case result: {bool(repro_ok)} "
        f"(n={int(cldn4.n)}, rho={float(cldn4.spearman_rho):.6f}).",
        "PRIMARY protein cohort is CCLE name suffix _LUNG, pairwise complete cases.",
        "PRIMARY RNA cohort is DepMap 24Q4 OncotreeLineage==Lung and ModelType==Cell Line.",
        "No imputation. BH q is across the 6 primary unadjusted tests only.",
        "Partial correlations are a different statistic and are not used to match ρ=0.69.",
        f"Protein keratin partial not estimated: KRT18/KRT19 absent; KRT8 fragment n={krt8_lung}/77, overlap with CLDN4 n={cldn4_and_krt8}.",
        f"Protein genes absent from Gygi: {protein_absent}.",
    ]
    focus = partial_df[
        partial_df.partner.isin(["CLDN4", "CLDN7", "CLDN1", "CLDN3"])
        & partial_df.cohort.isin(["lung_suffix", "lung_cell_lines"])
    ]
    for _, rec in focus.iterrows():
        rho = rec.partial_spearman
        rho_s = "NA" if pd.isna(rho) else f"{float(rho):.3f}"
        lines.append(
            f"PARTIAL {rec.layer} {rec.cohort} {rec.partner} adj={rec.adjustment}: n={rec.n} partial={rho_s}"
        )
    for rec in summary["primary"]:
        lines.append(
            f"{rec['layer']} {rec['partner']}: n={rec['n']} Spearman={rec['spearman_rho']:.4f} "
            f"p={rec['spearman_p']:.3e} q={rec['spearman_q_bh_primary6']:.3e} "
            f"rounds_to_0.69={rec['rounds_to_slide_0.69']} n_is_118={rec['n_equals_slide_118']}"
        )
    (out / "verdict.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    if not repro_ok:
        raise SystemExit("CLDN4 lung protein correlation did not reproduce the prior Gygi n=45 / ρ=0.693 result")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
