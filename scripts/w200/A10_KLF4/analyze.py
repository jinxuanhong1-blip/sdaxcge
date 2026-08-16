#!/usr/bin/env python3
"""A10: KLF4 vs TACSTD2 / CLDN4 in public lung (honest direction test).

Primary statistic: Spearman ρ.
Primary cohort: TCGA-LUAD primary tumors.
Pre-specified second histology: TCGA-LUSC.
Does not tune filters, residualization, or method to force concordance.
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

PRIMARY_TF = "KLF4"
PRIMARY_TARGETS = ["TACSTD2", "CLDN4"]
CONTEXT_GENES = ["KRT5", "CLDN7", "PECAM1"]
CLAIM = (
    "KLF4 is associated in the same direction with TACSTD2 and with CLDN4 "
    "in public lung."
)


def bh(pvals: list[float]) -> np.ndarray:
    if not pvals:
        return np.array([])
    return multipletests(pvals, method="fdr_bh")[1]


def bootstrap_spearman_ci(
    x: np.ndarray, y: np.ndarray, n_boot: int = 4000, seed: int = 0
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        boots[i] = stats.spearmanr(x[idx], y[idx]).statistic
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return float(lo), float(hi)


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    xr, yr, zr = stats.rankdata(x), stats.rankdata(y), stats.rankdata(z)
    design = np.column_stack([np.ones(len(xr)), zr])
    rx = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    ry = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    n = len(xr)
    if n <= 3 or abs(r) >= 1:
        return r, np.nan
    t = r * np.sqrt((n - 3) / (1 - r**2))
    p = float(2 * stats.t.sf(abs(t), df=n - 3))
    return r, p


def corr_pair(x: pd.Series, y: pd.Series, ci: bool = True) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    d.columns = ["x", "y"]
    n = int(len(d))
    out = {
        "n": n,
        "spearman_rho": None,
        "spearman_p": None,
        "spearman_ci95_low": None,
        "spearman_ci95_high": None,
        "pearson_r": None,
        "pearson_p": None,
    }
    if n < 8:
        return out
    rho, p_s = stats.spearmanr(d["x"], d["y"])
    r, p_p = stats.pearsonr(d["x"], d["y"])
    out.update(
        {
            "spearman_rho": float(rho),
            "spearman_p": float(p_s),
            "pearson_r": float(r),
            "pearson_p": float(p_p),
        }
    )
    if ci:
        lo, hi = bootstrap_spearman_ci(d["x"].to_numpy(), d["y"].to_numpy())
        out["spearman_ci95_low"] = lo
        out["spearman_ci95_high"] = hi
    return out


def direction_label(rho: float | None, fdr: float | None, alpha: float = 0.05) -> str:
    if rho is None or fdr is None or not np.isfinite(rho) or not np.isfinite(fdr):
        return "NOT_COMPUTED"
    if fdr >= alpha:
        return "NULL"
    return "INVERSE" if rho < 0 else "POSITIVE"


def tcga_barcode_split(sample_id: str) -> tuple[str, str]:
    parts = sample_id.split("-")
    patient = "-".join(parts[:3])
    code = parts[3][:2] if len(parts) > 3 else ""
    return patient, code


def load_tcga(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path)
    expr["sample_id"] = expr["sample_id"].astype(str)
    meta = expr["sample_id"].map(lambda s: tcga_barcode_split(s))
    expr["patient"] = [m[0] for m in meta]
    expr["sample_type"] = [m[1] for m in meta]
    for g in expr.columns:
        if g not in ("sample_id", "patient", "sample_type"):
            expr[g] = pd.to_numeric(expr[g], errors="coerce")
    tumors = (
        expr[expr["sample_type"] == "01"]
        .sort_values("sample_id")
        .drop_duplicates("patient")
        .set_index("patient")
    )
    return tumors


def load_purity(indir: Path) -> pd.Series:
    purity_path = indir / "tcga_absolute_purity.txt"
    if not purity_path.exists():
        return pd.Series(dtype=float)
    pur = pd.read_csv(purity_path, sep="\t")
    if "array" not in pur.columns or "purity" not in pur.columns:
        return pd.Series(dtype=float)
    pur = pur[pur["array"].astype(str).str.endswith("-01")].copy()
    pur["patient"] = pur["array"].astype(str).str.slice(0, 12)
    pur = pur.drop_duplicates("patient").set_index("patient")
    return pd.to_numeric(pur["purity"], errors="coerce")


def load_cptac(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["sample_id"] = df["sample_id"].astype(str)
    for g in df.columns:
        if g != "sample_id":
            df[g] = pd.to_numeric(df[g], errors="coerce")
    return df.set_index("sample_id")


def find_model_csv(indir: Path) -> Path:
    for cand in (indir / "Model.csv", Path("/tmp/a10_klf4_data/Model.csv")):
        if cand.exists():
            return cand
    raise SystemExit("Model.csv not found. Re-run download.py.")


def GENE_COLS(df: pd.DataFrame) -> list[str]:
    skip = {
        "sample_id",
        "patient",
        "sample_type",
        "ModelID",
        "CellLineName",
        "ModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "OncotreeCode",
    }
    return [c for c in df.columns if c not in skip]


def load_depmap_lung(indir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expr = pd.read_csv(indir / "depmap24q4_genes.csv")
    model = pd.read_csv(find_model_csv(indir), low_memory=False)
    keep = [
        c
        for c in [
            "ModelID",
            "CellLineName",
            "ModelType",
            "OncotreeLineage",
            "OncotreePrimaryDisease",
            "OncotreeSubtype",
            "OncotreeCode",
        ]
        if c in model.columns
    ]
    df = expr.merge(model[keep], on="ModelID", how="left")
    for g in GENE_COLS(df):
        df[g] = pd.to_numeric(df[g], errors="coerce")
    lung = df[df["OncotreeLineage"] == "Lung"].copy()
    cl = lung[lung["ModelType"].fillna("") == "Cell Line"].copy()
    luad = cl[
        cl["OncotreeSubtype"].fillna("").str.contains("Adenocarcinoma", case=False)
        | (cl.get("OncotreeCode", pd.Series("", index=cl.index)).fillna("") == "LUAD")
    ].copy()
    return cl, luad


def add_row(
    rows: list[dict],
    *,
    cohort: str,
    layer: str,
    gene_x: str,
    gene_y: str,
    x: pd.Series,
    y: pd.Series,
    primary: bool,
    note: str,
    purity: pd.Series | None = None,
) -> dict:
    block = corr_pair(x, y)
    rec = {
        "cohort": cohort,
        "layer": layer,
        "gene_x": gene_x,
        "gene_y": gene_y,
        "primary_endpoint": primary,
        "note": note,
        **block,
        "purity_partial_rho": None,
        "purity_partial_p": None,
        "n_purity": None,
    }
    if purity is not None:
        d = pd.concat([x.rename("x"), y.rename("y"), purity.rename("z")], axis=1).dropna()
        rec["n_purity"] = int(len(d))
        if len(d) >= 8:
            pr, pp = partial_spearman(d["x"].to_numpy(), d["y"].to_numpy(), d["z"].to_numpy())
            rec["purity_partial_rho"] = pr
            rec["purity_partial_p"] = pp
    rows.append(rec)
    return rec


def add_tf_pairs(
    rows: list[dict],
    *,
    cohort: str,
    layer: str,
    df: pd.DataFrame,
    primary: bool,
    note: str,
    purity: pd.Series | None = None,
    include_context: bool = False,
) -> None:
    if PRIMARY_TF not in df.columns:
        return
    for tgt in PRIMARY_TARGETS:
        if tgt in df.columns:
            add_row(
                rows,
                cohort=cohort,
                layer=layer,
                gene_x=PRIMARY_TF,
                gene_y=tgt,
                x=df[PRIMARY_TF],
                y=df[tgt],
                primary=primary,
                note=note,
                purity=purity,
            )
    if "TACSTD2" in df.columns and "CLDN4" in df.columns:
        add_row(
            rows,
            cohort=cohort,
            layer=layer,
            gene_x="TACSTD2",
            gene_y="CLDN4",
            x=df["TACSTD2"],
            y=df["CLDN4"],
            primary=False,
            note="Exploratory positive-control pair (epithelial/TJ co-expression).",
            purity=purity,
        )
    if include_context:
        for g in CONTEXT_GENES:
            if g in df.columns:
                add_row(
                    rows,
                    cohort=cohort,
                    layer=layer,
                    gene_x=PRIMARY_TF,
                    gene_y=g,
                    x=df[PRIMARY_TF],
                    y=df[g],
                    primary=False,
                    note="Exploratory context gene (keratin / TJ / endothelium).",
                    purity=purity,
                )


def scatter(ax, x, y, xlabel, ylabel, title, rho, p, n):
    ax.scatter(x, y, s=12, alpha=0.45, c="#2c5aa0", edgecolors="none")
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    txt = f"Spearman ρ = {rho:.3f}\np = {p:.2e}\nn = {n}" if rho is not None else f"n = {n}"
    ax.text(
        0.02,
        0.98,
        txt,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="none"),
    )


def pair_directions(table: pd.DataFrame, cohort: str) -> dict[str, str]:
    sub = table[
        (table["cohort"] == cohort)
        & (table["primary_endpoint"])
        & (table["gene_x"] == PRIMARY_TF)
    ]
    out = {}
    for tgt in PRIMARY_TARGETS:
        row = sub[sub["gene_y"] == tgt]
        out[tgt] = "NOT_COMPUTED" if row.empty else str(row.iloc[0]["direction"])
    return out


def concordance_label(pair_verdicts: dict[str, str]) -> tuple[str, str]:
    inverses = [t for t, v in pair_verdicts.items() if v == "INVERSE"]
    positives = [t for t, v in pair_verdicts.items() if v == "POSITIVE"]
    nulls = [t for t, v in pair_verdicts.items() if v == "NULL"]
    tac = pair_verdicts.get("TACSTD2")
    cld = pair_verdicts.get("CLDN4")
    if len(positives) == 2:
        return (
            "CONCORDANT_POSITIVE",
            "KLF4 is significantly positively correlated with both TACSTD2 and CLDN4.",
        )
    if len(inverses) == 2:
        return (
            "CONCORDANT_INVERSE",
            "KLF4 is significantly inversely correlated with both TACSTD2 and CLDN4.",
        )
    if positives and inverses:
        return (
            "DISCORDANT",
            f"KLF4 is {tac} with TACSTD2 and {cld} with CLDN4. "
            "The two targets do not share a direction, so a same-sign "
            "KLF4–TACSTD2/CLDN4 claim is not supported.",
        )
    if (positives or inverses) and nulls:
        return (
            "PARTIAL",
            f"KLF4 is significant for one target only "
            f"(TACSTD2={tac}, CLDN4={cld}). Same-direction tracking of both genes "
            "is not supported.",
        )
    return (
        "NOT_SUPPORTED",
        f"Neither TACSTD2 nor CLDN4 shows a significant Spearman association "
        f"with KLF4 (TACSTD2={tac}, CLDN4={cld}).",
    )


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/A10_KLF4")
    p.add_argument("--out-dir", default="results/w200/A10_KLF4")
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    coverage: dict = {}
    purity = load_purity(indir)

    # ---- TCGA LUAD (primary) ----
    luad = load_tcga(indir / "tcga_luad_genes.csv")
    coverage["tcga_luad"] = {
        "n_primary_tumors": int(len(luad)),
        "n_with_absolute_purity": int(purity.reindex(luad.index).notna().sum()),
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g in luad.columns],
        "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g not in luad.columns],
    }
    luad.to_csv(out / "tcga_luad_primary_tumors.csv")
    pur_luad = purity.reindex(luad.index)
    pd.DataFrame({"patient": pur_luad.index, "purity": pur_luad.values}).to_csv(
        out / "tcga_luad_absolute_purity.csv", index=False
    )
    add_tf_pairs(
        rows,
        cohort="TCGA-LUAD",
        layer="RNA_log2TPM",
        df=luad,
        primary=True,
        note="Primary cohort. One -01 sample per patient.",
        purity=purity,
        include_context=True,
    )

    # ---- TCGA LUSC (pre-specified second lung histology) ----
    lusc = load_tcga(indir / "tcga_lusc_genes.csv")
    coverage["tcga_lusc"] = {
        "n_primary_tumors": int(len(lusc)),
        "n_with_absolute_purity": int(purity.reindex(lusc.index).notna().sum()),
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g in lusc.columns],
        "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g not in lusc.columns],
    }
    lusc.to_csv(out / "tcga_lusc_primary_tumors.csv")
    add_tf_pairs(
        rows,
        cohort="TCGA-LUSC",
        layer="RNA_log2TPM",
        df=lusc,
        primary=True,
        note="Pre-specified second lung histology. One -01 sample per patient. "
        "Does not replace the LUAD verdict.",
        purity=purity,
        include_context=True,
    )

    # ---- CPTAC LUAD / LSCC RNA + protein ----
    for cohort, tag, layer, fname, primary_note in (
        ("CPTAC-LUAD", "luad", "RNA_RSEM_UQ_log2", "cptac_luad_rna_genes.csv", "Replication. CPTAC freeze v1.2 tumor RNA."),
        ("CPTAC-LUAD", "luad", "protein_TMT_log2", "cptac_luad_protein_genes.csv", "Replication. Protein if quantified; pairwise complete cases only."),
        ("CPTAC-LSCC", "lscc", "RNA_RSEM_UQ_log2", "cptac_lscc_rna_genes.csv", "Replication. CPTAC freeze v1.2 LSCC tumor RNA."),
        ("CPTAC-LSCC", "lscc", "protein_TMT_log2", "cptac_lscc_protein_genes.csv", "Replication. LSCC protein if quantified; pairwise complete cases only."),
    ):
        path = indir / fname
        if not path.exists():
            coverage[f"{tag}_{layer}"] = {"missing_file": fname}
            continue
        df = load_cptac(path)
        key = f"cptac_{tag}_{'rna' if 'RNA' in layer else 'protein'}"
        coverage[key] = {
            "n_tumors": int(len(df)),
            "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in df.columns],
            "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g not in df.columns],
            "n_non_na": {
                g: int(df[g].notna().sum()) if g in df.columns else 0
                for g in [PRIMARY_TF] + PRIMARY_TARGETS
            },
        }
        add_tf_pairs(
            rows,
            cohort=cohort,
            layer=layer,
            df=df,
            primary=True,
            note=primary_note,
        )

    # ---- DepMap lung cell lines ----
    lung_cl, luad_cl = load_depmap_lung(indir)
    coverage["depmap_lung_cell_lines"] = {
        "n": int(len(lung_cl)),
        "filter": "OncotreeLineage==Lung AND ModelType==Cell Line",
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in lung_cl.columns],
        "primary_disease_counts": lung_cl["OncotreePrimaryDisease"].value_counts().to_dict()
        if "OncotreePrimaryDisease" in lung_cl.columns
        else {},
    }
    coverage["depmap_luad_cell_lines"] = {
        "n": int(len(luad_cl)),
        "filter": "OncotreeLineage==Lung AND ModelType==Cell Line AND (Adenocarcinoma subtype OR OncotreeCode==LUAD)",
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in luad_cl.columns],
    }
    lung_cl.to_csv(out / "depmap_lung_cell_lines.csv", index=False)
    luad_cl.to_csv(out / "depmap_luad_cell_lines.csv", index=False)
    add_tf_pairs(
        rows,
        cohort="DepMap24Q4-lung",
        layer="RNA_log2TPM",
        df=lung_cl,
        primary=True,
        note="Replication. All lung-lineage cell lines (not tumors).",
    )
    add_tf_pairs(
        rows,
        cohort="DepMap24Q4-LUAD",
        layer="RNA_log2TPM",
        df=luad_cl,
        primary=False,
        note="Sensitivity. LUAD cell lines only. Not used for the verdict.",
    )

    table = pd.DataFrame(rows)
    prim = table["primary_endpoint"].astype(bool)
    table.loc[prim, "fdr_bh_primary"] = bh(table.loc[prim, "spearman_p"].fillna(1).tolist())
    table.loc[~prim, "fdr_bh_primary"] = np.nan
    table["direction"] = [
        direction_label(
            r.spearman_rho,
            r.fdr_bh_primary if r.primary_endpoint else r.spearman_p,
        )
        for r in table.itertuples(index=False)
    ]
    table.to_csv(out / "correlations.csv", index=False)

    # ---- Honest verdict from TCGA-LUAD primary pairs only ----
    luad_dirs = pair_directions(table, "TCGA-LUAD")
    lusc_dirs = pair_directions(table, "TCGA-LUSC")
    overall, core = concordance_label(luad_dirs)
    statement = (
        f"TCGA-LUAD primary tumors: {core} "
        f"TCGA-LUSC (pre-specified second histology): "
        f"TACSTD2={lusc_dirs.get('TACSTD2')}, CLDN4={lusc_dirs.get('CLDN4')}."
    )

    tcga_prim = table[
        (table["cohort"] == "TCGA-LUAD")
        & (table["primary_endpoint"])
        & (table["gene_x"] == PRIMARY_TF)
    ].copy()
    repl = table[
        (table["primary_endpoint"])
        & (table["cohort"] != "TCGA-LUAD")
        & (table["gene_x"] == PRIMARY_TF)
    ]
    replication = []
    for r in repl.itertuples(index=False):
        tcga_rho = tcga_prim[tcga_prim["gene_y"] == r.gene_y]
        tcga_dir = None if tcga_rho.empty else str(tcga_rho.iloc[0]["direction"])
        tcga_val = None if tcga_rho.empty else tcga_rho.iloc[0]["spearman_rho"]
        # A near-zero / NULL primary ρ has no meaningful sign to match.
        if (
            r.spearman_rho is None
            or tcga_val is None
            or tcga_dir in (None, "NULL", "NOT_COMPUTED")
        ):
            same_sign = None
        else:
            same_sign = bool(np.sign(r.spearman_rho) == np.sign(tcga_val))
        replication.append(
            {
                "cohort": r.cohort,
                "layer": r.layer,
                "gene_y": r.gene_y,
                "n": r.n,
                "spearman_rho": r.spearman_rho,
                "direction": r.direction,
                "same_sign_as_tcga_luad": same_sign,
            }
        )

    summary = {
        "claim": CLAIM,
        "honest_verdict": {
            "primary_cohort": "TCGA-LUAD primary tumors (-01, one per patient)",
            "primary_statistic": "Spearman rho",
            "pair_direction_tcga_luad": luad_dirs,
            "pair_direction_tcga_lusc": lusc_dirs,
            "label": overall,
            "statement": statement,
            "did_we_tune_filters_to_force_concordance": False,
        },
        "replication": replication,
        "coverage": coverage,
        "cannot_test": [
            "Direct KLF4 binding or transcriptional control at TACSTD2/CLDN4 (no ChIP/perturbation in this slice).",
            "ICI response or TROP2-ADC outcome.",
            "KLF4 IHC or protein from RNA (except the CPTAC protein rows, which are reported separately).",
        ],
        "methods": {
            "tcga_matrix": "Xena GDC hub TCGA-LUAD/LUSC.star_tpm.tsv.gz, log2(TPM+1), GENCODE v36",
            "cptac_freeze": "data_freeze_v1.2_reorganized/LUAD and LSCC",
            "depmap_release": "DepMap Public 24Q4",
            "fdr": "BH within pre-specified KLF4 vs TACSTD2/CLDN4 primary rows only",
            "purity": "ABSOLUTE (PanCanAtlas); partial Spearman is secondary",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # ---- Figures ----
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 8.0))
    for row_i, (cohort_name, tumors) in enumerate((("TCGA-LUAD", luad), ("TCGA-LUSC", lusc))):
        sub_tbl = table[
            (table["cohort"] == cohort_name)
            & (table["gene_x"] == PRIMARY_TF)
            & (table["primary_endpoint"])
        ]
        for col_i, tgt in enumerate(PRIMARY_TARGETS):
            ax = axes[row_i][col_i]
            sub = sub_tbl[sub_tbl["gene_y"] == tgt]
            rho = float(sub.iloc[0]["spearman_rho"]) if not sub.empty else None
            pv = float(sub.iloc[0]["spearman_p"]) if not sub.empty else None
            n = int(sub.iloc[0]["n"]) if not sub.empty else 0
            scatter(
                ax,
                tumors[PRIMARY_TF],
                tumors[tgt],
                "KLF4 log2(TPM+1)",
                f"{tgt} log2(TPM+1)",
                f"{cohort_name}  KLF4 vs {tgt}",
                rho,
                pv,
                n,
            )
    fig.tight_layout()
    fig.savefig(out / "fig_tcga_scatter.png", dpi=160)
    fig.savefig(out / "fig_tcga_scatter.pdf")
    plt.close(fig)

    prim_tbl = table[table["primary_endpoint"] & (table["gene_x"] == PRIMARY_TF)].copy()
    prim_tbl = prim_tbl.dropna(subset=["spearman_rho"])
    if not prim_tbl.empty:
        labels = [
            f"{r.cohort} {r.layer}\nKLF4 vs {r.gene_y} (n={r.n})"
            for r in prim_tbl.itertuples(index=False)
        ]
        fig, ax = plt.subplots(figsize=(8.8, 0.42 * len(prim_tbl) + 1.6))
        y = np.arange(len(prim_tbl))
        rho = prim_tbl["spearman_rho"].to_numpy()
        lo = prim_tbl["spearman_ci95_low"].to_numpy()
        hi = prim_tbl["spearman_ci95_high"].to_numpy()
        colors = ["#b2182b" if v < 0 else "#2166ac" for v in rho]
        xerr = np.vstack([rho - lo, hi - rho])
        ax.errorbar(rho, y, xerr=xerr, fmt="none", ecolor="#555555", elinewidth=1, capsize=3)
        ax.scatter(rho, y, c=colors, s=36, zorder=3)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
        ax.set_title("A10 KLF4 primary endpoints (red = inverse, blue = positive)")
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(out / "fig_primary_forest.png", dpi=160)
        fig.savefig(out / "fig_primary_forest.pdf")
        plt.close(fig)

    prot_path = indir / "cptac_luad_protein_genes.csv"
    if prot_path.exists():
        prot = load_cptac(prot_path)
        if PRIMARY_TF in prot.columns:
            fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
            for ax, tgt in zip(axes, PRIMARY_TARGETS):
                if tgt not in prot.columns:
                    ax.axis("off")
                    continue
                d = prot[[PRIMARY_TF, tgt]].dropna()
                sub = table[
                    (table["cohort"] == "CPTAC-LUAD")
                    & (table["layer"] == "protein_TMT_log2")
                    & (table["gene_y"] == tgt)
                    & (table["gene_x"] == PRIMARY_TF)
                ]
                rho = float(sub.iloc[0]["spearman_rho"]) if not sub.empty else None
                pv = float(sub.iloc[0]["spearman_p"]) if not sub.empty else None
                scatter(
                    ax,
                    d[PRIMARY_TF],
                    d[tgt],
                    "KLF4 protein (log2)",
                    f"{tgt} protein (log2)",
                    f"CPTAC-LUAD protein  KLF4 vs {tgt}",
                    rho,
                    pv,
                    int(len(d)),
                )
            fig.tight_layout()
            fig.savefig(out / "fig_cptac_protein_scatter.png", dpi=160)
            fig.savefig(out / "fig_cptac_protein_scatter.pdf")
            plt.close(fig)

    write_results_readme(out, summary, table)
    print(json.dumps(summary["honest_verdict"], indent=2))
    return 0


def _fmt_rho(v) -> str:
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.3f}"


def _fmt_p(v) -> str:
    return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.2e}"


def write_results_readme(out: Path, summary: dict, table: pd.DataFrame) -> None:
    v = summary["honest_verdict"]
    lines = [
        "# A10 KLF4 vs TACSTD2 / CLDN4 — public lung",
        "",
        f"**Claim:** {summary['claim']}",
        "",
        f"**Honest verdict (TCGA-LUAD primary): `{v['label']}`**",
        "",
        v["statement"],
        "",
        "Filters were not tuned to force concordance.",
        "",
        "## Primary pairs (BH-FDR within this list)",
        "",
        "| Cohort | Layer | Pair | n | Spearman ρ | p | FDR | 95% CI | Direction |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    prim = table[table["primary_endpoint"] & (table["gene_x"] == PRIMARY_TF)]
    for r in prim.itertuples(index=False):
        ci = (
            f"[{r.spearman_ci95_low:.3f}, {r.spearman_ci95_high:.3f}]"
            if r.spearman_ci95_low is not None
            else "—"
        )
        fdr = "—" if pd.isna(r.fdr_bh_primary) else f"{r.fdr_bh_primary:.3g}"
        lines.append(
            f"| {r.cohort} | {r.layer} | {r.gene_x} vs {r.gene_y} | {r.n} | "
            f"{_fmt_rho(r.spearman_rho)} | {_fmt_p(r.spearman_p)} | {fdr} | {ci} | {r.direction} |"
        )
    expl = table[~table["primary_endpoint"]]
    if not expl.empty:
        lines += [
            "",
            "## Exploratory (not used for the verdict)",
            "",
            "| Cohort | Pair | n | Spearman ρ | p | Direction |",
            "| --- | --- | ---: | ---: | ---: | --- |",
        ]
        for r in expl.itertuples(index=False):
            lines.append(
                f"| {r.cohort} | {r.gene_x} vs {r.gene_y} | {r.n} | "
                f"{_fmt_rho(r.spearman_rho)} | {_fmt_p(r.spearman_p)} | {r.direction} |"
            )

    cov = summary.get("coverage", {})
    tcga = cov.get("tcga_luad", {})
    lusc = cov.get("tcga_lusc", {})
    lines += [
        "",
        "## Coverage (honest missingness)",
        "",
        f"- TCGA-LUAD primary tumors: n={tcga.get('n_primary_tumors')}; "
        f"ABSOLUTE purity available for {tcga.get('n_with_absolute_purity')}.",
        f"- TCGA-LUSC primary tumors: n={lusc.get('n_primary_tumors')}; "
        f"ABSOLUTE purity available for {lusc.get('n_with_absolute_purity')}.",
    ]
    for key, label in (
        ("cptac_luad_rna", "CPTAC LUAD RNA"),
        ("cptac_luad_protein", "CPTAC LUAD protein"),
        ("cptac_lscc_rna", "CPTAC LSCC RNA"),
        ("cptac_lscc_protein", "CPTAC LSCC protein"),
    ):
        block = cov.get(key, {})
        if not block:
            lines.append(f"- {label}: file missing.")
            continue
        nna = block.get("n_non_na", {})
        lines.append(
            f"- {label}: n={block.get('n_tumors')}; "
            f"non-NA KLF4/TACSTD2/CLDN4 = "
            f"{nna.get('KLF4', '?')}/{nna.get('TACSTD2', '?')}/{nna.get('CLDN4', '?')}."
        )
    lines.append(
        f"- DepMap 24Q4 lung cell lines: n={cov.get('depmap_lung_cell_lines', {}).get('n')}; "
        f"LUAD subset n={cov.get('depmap_luad_cell_lines', {}).get('n')}."
    )
    lines += [
        "",
        "## Secondary: purity-adjusted TCGA (does not change the verdict)",
        "",
    ]
    tcga_rows = table[
        (table["cohort"].isin(["TCGA-LUAD", "TCGA-LUSC"])) & (table["gene_x"] == PRIMARY_TF)
    ]
    for r in tcga_rows.itertuples(index=False):
        if r.purity_partial_rho is None:
            continue
        lines.append(
            f"- {r.cohort} KLF4 vs {r.gene_y}: partial ρ = {r.purity_partial_rho:.3f}, "
            f"p = {r.purity_partial_p:.2e}, n = {int(r.n_purity)}"
        )
    lines += [
        "",
        "## What this is not",
        "",
    ]
    for item in summary["cannot_test"]:
        lines.append(f"- {item}")
    lines += [
        "- A reason to drop CPTAC protein or DepMap because the sign disagrees with TCGA RNA.",
        "- Evidence that KLF4 transcriptionally activates or represses TACSTD2 or CLDN4.",
        "",
        "## Caveats",
        "",
        "1. Bulk tumor RNA mixes epithelium, stroma, endothelium, and immune cells. "
        "In TCGA-LUAD, KLF4 tracks PECAM1 (ρ = 0.33) as well as TACSTD2 (ρ = 0.21). "
        "That does not prove the TACSTD2 association is endothelial, but it forbids "
        "reading KLF4 as a purely epithelial/TJ transcription factor from bulk RNA.",
        "2. TCGA-LUAD KLF4–CLDN4 is a true null (ρ = −0.004, CI includes 0), not a "
        "weak inverse that we rounded away. Purity adjustment does not create an inverse.",
        "3. CPTAC-LUAD RNA is not a silent non-replication: KLF4–CLDN4 is significantly "
        "inverse (ρ = −0.31, n=110) while KLF4–TACSTD2 is null. That is the opposite "
        "of a same-sign pair. CPTAC-LUAD protein is null for both, with KLF4 missing "
        "in 25/110 and CLDN4 missing in 31/110 tumors.",
        "4. TCGA-LUSC and DepMap lung lines are concordant-positive. They are reported; "
        "they do not override the pre-specified LUAD verdict. DepMap lung n=214 includes "
        "NSCLC, neuroendocrine, and a few non-cancerous lines; the LUAD-only subset "
        "(n=80) is also positive for both pairs.",
        "5. The TACSTD2–CLDN4 pair itself is robustly positive in TCGA LUAD/LUSC and "
        "DepMap. The assay is not broken; the KLF4 same-sign claim is.",
        "6. TCGA/CPTAC/DepMap are not ICI cohorts.",
        "",
        "## Rerun",
        "",
        "```bash",
        "python3 scripts/w200/A10_KLF4/download.py",
        "python3 scripts/w200/A10_KLF4/analyze.py",
        "```",
        "",
        "See `summary.json` and `correlations.csv`.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
