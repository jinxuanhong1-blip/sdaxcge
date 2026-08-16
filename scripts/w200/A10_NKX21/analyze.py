#!/usr/bin/env python3
"""A10: NKX2-1 vs TACSTD2 / CLDN4 in public LUAD (honest direction test).

Primary statistic: Spearman ρ.
Primary cohort: TCGA-LUAD primary tumors.
Does not tune filters, residualization, or method to force an inverse.
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

PRIMARY_TF = "NKX2-1"
PRIMARY_TARGETS = ["TACSTD2", "CLDN4"]
CONTEXT_GENES = ["SFTPB", "NAPSA", "KRT5", "CLDN7"]
CLAIM = (
    "NKX2-1 is inversely correlated with TACSTD2 and with CLDN4 in public LUAD."
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


def load_tcga(indir: Path) -> tuple[pd.DataFrame, pd.Series]:
    expr = pd.read_csv(indir / "tcga_luad_genes.csv")
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
    purity_path = indir / "tcga_absolute_purity.txt"
    purity = pd.Series(dtype=float)
    if purity_path.exists():
        pur = pd.read_csv(purity_path, sep="\t")
        if "array" in pur.columns and "purity" in pur.columns:
            pur = pur[pur["array"].astype(str).str.endswith("-01")].copy()
            pur["patient"] = pur["array"].astype(str).str.slice(0, 12)
            pur = pur.drop_duplicates("patient").set_index("patient")
            purity = pd.to_numeric(pur["purity"], errors="coerce")
    return tumors, purity


def load_cptac(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["sample_id"] = df["sample_id"].astype(str)
    for g in df.columns:
        if g != "sample_id":
            df[g] = pd.to_numeric(df[g], errors="coerce")
    return df.set_index("sample_id")


def find_model_csv(indir: Path) -> Path:
    for cand in (indir / "Model.csv", Path("/tmp/a10_nkx21_data/Model.csv")):
        if cand.exists():
            return cand
    raise SystemExit("Model.csv not found. Re-run download.py.")


def load_depmap_luad(indir: Path) -> pd.DataFrame:
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
    # LUAD: adenocarcinoma subtype or OncotreeCode LUAD. Do not invent extra filters.
    luad = cl[
        cl["OncotreeSubtype"].fillna("").str.contains("Adenocarcinoma", case=False)
        | (cl.get("OncotreeCode", pd.Series("", index=cl.index)).fillna("") == "LUAD")
    ].copy()
    return luad


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


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/A10_NKX21")
    p.add_argument("--out-dir", default="results/w200/A10_NKX21")
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    coverage: dict = {}

    # ---- TCGA LUAD (primary) ----
    tumors, purity = load_tcga(indir)
    coverage["tcga_luad"] = {
        "n_primary_tumors": int(len(tumors)),
        "n_with_absolute_purity": int(purity.reindex(tumors.index).notna().sum()),
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g in tumors.columns],
        "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g not in tumors.columns],
    }
    tumors.to_csv(out / "tcga_luad_primary_tumors.csv")

    for tgt in PRIMARY_TARGETS:
        add_row(
            rows,
            cohort="TCGA-LUAD",
            layer="RNA_log2TPM",
            gene_x=PRIMARY_TF,
            gene_y=tgt,
            x=tumors[PRIMARY_TF],
            y=tumors[tgt],
            primary=True,
            note="Primary cohort. One -01 sample per patient.",
            purity=purity,
        )
    # TACSTD2–CLDN4 positive-control pair (not part of the inverse claim)
    add_row(
        rows,
        cohort="TCGA-LUAD",
        layer="RNA_log2TPM",
        gene_x="TACSTD2",
        gene_y="CLDN4",
        x=tumors["TACSTD2"],
        y=tumors["CLDN4"],
        primary=False,
        note="Exploratory positive-control pair (epithelial/TJ co-expression).",
        purity=purity,
    )
    for g in CONTEXT_GENES:
        if g in tumors.columns:
            add_row(
                rows,
                cohort="TCGA-LUAD",
                layer="RNA_log2TPM",
                gene_x=PRIMARY_TF,
                gene_y=g,
                x=tumors[PRIMARY_TF],
                y=tumors[g],
                primary=False,
                note="Exploratory context gene (alveolar/basal/TJ).",
                purity=purity,
            )

    # ---- CPTAC RNA ----
    rna = load_cptac(indir / "cptac_luad_rna_genes.csv")
    coverage["cptac_rna"] = {
        "n_tumors": int(len(rna)),
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in rna.columns],
        "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g not in rna.columns],
    }
    for tgt in PRIMARY_TARGETS:
        if PRIMARY_TF in rna.columns and tgt in rna.columns:
            add_row(
                rows,
                cohort="CPTAC-LUAD",
                layer="RNA_RSEM_UQ_log2",
                gene_x=PRIMARY_TF,
                gene_y=tgt,
                x=rna[PRIMARY_TF],
                y=rna[tgt],
                primary=True,
                note="Replication. CPTAC freeze v1.2 tumor RNA.",
            )
    if "TACSTD2" in rna.columns and "CLDN4" in rna.columns:
        add_row(
            rows,
            cohort="CPTAC-LUAD",
            layer="RNA_RSEM_UQ_log2",
            gene_x="TACSTD2",
            gene_y="CLDN4",
            x=rna["TACSTD2"],
            y=rna["CLDN4"],
            primary=False,
            note="Exploratory positive-control pair.",
        )

    # ---- CPTAC protein ----
    prot = load_cptac(indir / "cptac_luad_protein_genes.csv")
    coverage["cptac_protein"] = {
        "n_tumors": int(len(prot)),
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in prot.columns],
        "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g not in prot.columns],
        "n_non_na": {
            g: int(prot[g].notna().sum()) if g in prot.columns else 0
            for g in [PRIMARY_TF] + PRIMARY_TARGETS
        },
    }
    for tgt in PRIMARY_TARGETS:
        if PRIMARY_TF in prot.columns and tgt in prot.columns:
            add_row(
                rows,
                cohort="CPTAC-LUAD",
                layer="protein_TMT_log2",
                gene_x=PRIMARY_TF,
                gene_y=tgt,
                x=prot[PRIMARY_TF],
                y=prot[tgt],
                primary=True,
                note="Replication. Protein if quantified; pairwise complete cases only.",
            )
    if "TACSTD2" in prot.columns and "CLDN4" in prot.columns:
        add_row(
            rows,
            cohort="CPTAC-LUAD",
            layer="protein_TMT_log2",
            gene_x="TACSTD2",
            gene_y="CLDN4",
            x=prot["TACSTD2"],
            y=prot["CLDN4"],
            primary=False,
            note="Exploratory positive-control pair.",
        )

    # ---- DepMap LUAD cell lines ----
    luad_cl = load_depmap_luad(indir)
    coverage["depmap_luad_cell_lines"] = {
        "n": int(len(luad_cl)),
        "filter": "OncotreeLineage==Lung AND ModelType==Cell Line AND (Adenocarcinoma subtype OR OncotreeCode==LUAD)",
        "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in luad_cl.columns],
        "primary_disease_counts": luad_cl["OncotreePrimaryDisease"].value_counts().to_dict()
        if "OncotreePrimaryDisease" in luad_cl.columns
        else {},
    }
    luad_cl.to_csv(out / "depmap_luad_cell_lines.csv", index=False)
    for tgt in PRIMARY_TARGETS:
        if PRIMARY_TF in luad_cl.columns and tgt in luad_cl.columns:
            add_row(
                rows,
                cohort="DepMap24Q4-LUAD",
                layer="RNA_log2TPM",
                gene_x=PRIMARY_TF,
                gene_y=tgt,
                x=luad_cl[PRIMARY_TF],
                y=luad_cl[tgt],
                primary=True,
                note="Replication. In vitro cell lines, not tumors.",
            )
    if "TACSTD2" in luad_cl.columns and "CLDN4" in luad_cl.columns:
        add_row(
            rows,
            cohort="DepMap24Q4-LUAD",
            layer="RNA_log2TPM",
            gene_x="TACSTD2",
            gene_y="CLDN4",
            x=luad_cl["TACSTD2"],
            y=luad_cl["CLDN4"],
            primary=False,
            note="Exploratory positive-control pair.",
        )

    table = pd.DataFrame(rows)
    prim = table["primary_endpoint"].astype(bool)
    table.loc[prim, "fdr_bh_primary"] = bh(table.loc[prim, "spearman_p"].fillna(1).tolist())
    table.loc[~prim, "fdr_bh_primary"] = np.nan
    table["direction"] = [
        direction_label(r.spearman_rho, r.fdr_bh_primary if r.primary_endpoint else r.spearman_p)
        for r in table.itertuples(index=False)
    ]
    # For exploratory rows, direction uses raw p (labeled as such).
    table.loc[~prim, "direction"] = [
        direction_label(rho, p)
        for rho, p in zip(table.loc[~prim, "spearman_rho"], table.loc[~prim, "spearman_p"])
    ]
    table.to_csv(out / "correlations.csv", index=False)

    # ---- Honest verdict from TCGA primary pairs only ----
    tcga_prim = table[
        (table["cohort"] == "TCGA-LUAD")
        & (table["primary_endpoint"])
        & (table["gene_x"] == PRIMARY_TF)
    ].copy()
    pair_verdicts = {}
    for tgt in PRIMARY_TARGETS:
        sub = tcga_prim[tcga_prim["gene_y"] == tgt]
        if sub.empty:
            pair_verdicts[tgt] = "NOT_COMPUTED"
        else:
            pair_verdicts[tgt] = str(sub.iloc[0]["direction"])

    inverses = [t for t, v in pair_verdicts.items() if v == "INVERSE"]
    positives = [t for t, v in pair_verdicts.items() if v == "POSITIVE"]
    nulls = [t for t, v in pair_verdicts.items() if v == "NULL"]
    if len(inverses) == 2:
        overall = "SUPPORTED"
        statement = (
            "TCGA-LUAD primary tumors: NKX2-1 is significantly inversely correlated "
            "with both TACSTD2 and CLDN4 (BH-FDR < 0.05 within the pre-specified list)."
        )
    elif len(inverses) == 1:
        overall = "PARTIAL"
        statement = (
            f"TCGA-LUAD primary tumors: NKX2-1 is significantly inverse with {inverses[0]} "
            f"but not with {PRIMARY_TARGETS[1] if inverses[0] == PRIMARY_TARGETS[0] else PRIMARY_TARGETS[0]} "
            f"(that pair is {pair_verdicts[PRIMARY_TARGETS[1] if inverses[0] == PRIMARY_TARGETS[0] else PRIMARY_TARGETS[0]]}). "
            "The two-gene inverse claim is only partly supported."
        )
    elif positives and not inverses:
        overall = "OPPOSITE"
        statement = (
            "TCGA-LUAD primary tumors: NKX2-1 is significantly positively correlated "
            f"with {', '.join(positives)}. The inverse claim is not supported."
        )
    else:
        overall = "NOT_SUPPORTED"
        statement = (
            "TCGA-LUAD primary tumors: neither TACSTD2 nor CLDN4 shows a significant "
            f"inverse Spearman correlation with NKX2-1 "
            f"(TACSTD2={pair_verdicts.get('TACSTD2')}, CLDN4={pair_verdicts.get('CLDN4')})."
        )

    # Replication concordance (descriptive; does not override primary)
    repl = table[
        (table["primary_endpoint"])
        & (table["cohort"] != "TCGA-LUAD")
        & (table["gene_x"] == PRIMARY_TF)
    ]
    replication = []
    for r in repl.itertuples(index=False):
        replication.append(
            {
                "cohort": r.cohort,
                "layer": r.layer,
                "gene_y": r.gene_y,
                "n": r.n,
                "spearman_rho": r.spearman_rho,
                "direction": r.direction,
                "same_sign_as_tcga": (
                    None
                    if r.spearman_rho is None
                    or tcga_prim[tcga_prim["gene_y"] == r.gene_y].empty
                    or tcga_prim[tcga_prim["gene_y"] == r.gene_y].iloc[0]["spearman_rho"] is None
                    else bool(
                        np.sign(r.spearman_rho)
                        == np.sign(tcga_prim[tcga_prim["gene_y"] == r.gene_y].iloc[0]["spearman_rho"])
                    )
                ),
            }
        )

    summary = {
        "claim": CLAIM,
        "honest_verdict": {
            "primary_cohort": "TCGA-LUAD primary tumors (-01, one per patient)",
            "primary_statistic": "Spearman rho",
            "pair_direction_tcga": pair_verdicts,
            "label": overall,
            "statement": statement,
            "did_we_tune_filters_to_force_inverse": False,
        },
        "replication": replication,
        "coverage": coverage,
        "cannot_test": [
            "Direct NKX2-1 binding or repression at TACSTD2/CLDN4 (no ChIP/perturbation in this slice).",
            "ICI response or TROP2-ADC outcome.",
            "TTF-1 IHC (this is continuous RNA/protein, not a clinical IHC call).",
        ],
        "methods": {
            "tcga_matrix": "Xena GDC hub TCGA-LUAD.star_tpm.tsv.gz, log2(TPM+1), GENCODE v36",
            "cptac_freeze": "data_freeze_v1.2_reorganized/LUAD",
            "depmap_release": "DepMap Public 24Q4",
            "fdr": "BH within pre-specified primary NKX2-1 vs TACSTD2/CLDN4 rows only",
            "purity": "ABSOLUTE (PanCanAtlas); partial Spearman is secondary",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # ---- Figures ----
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
    for ax, tgt in zip(axes, PRIMARY_TARGETS):
        sub = tcga_prim[tcga_prim["gene_y"] == tgt]
        rho = float(sub.iloc[0]["spearman_rho"]) if not sub.empty else None
        pv = float(sub.iloc[0]["spearman_p"]) if not sub.empty else None
        n = int(sub.iloc[0]["n"]) if not sub.empty else 0
        scatter(
            ax,
            tumors[PRIMARY_TF],
            tumors[tgt],
            "NKX2-1 log2(TPM+1)",
            f"{tgt} log2(TPM+1)",
            f"TCGA-LUAD  NKX2-1 vs {tgt}",
            rho,
            pv,
            n,
        )
    fig.tight_layout()
    fig.savefig(out / "fig_tcga_scatter.png", dpi=160)
    fig.savefig(out / "fig_tcga_scatter.pdf")
    plt.close(fig)

    # Forest of primary rows
    prim_tbl = table[table["primary_endpoint"]].copy()
    prim_tbl = prim_tbl.dropna(subset=["spearman_rho"])
    if not prim_tbl.empty:
        labels = [
            f"{r.cohort} {r.layer}\nNKX2-1 vs {r.gene_y} (n={r.n})"
            for r in prim_tbl.itertuples(index=False)
        ]
        fig, ax = plt.subplots(figsize=(8.5, 0.45 * len(prim_tbl) + 1.6))
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
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
        ax.set_title("A10 primary endpoints (red = inverse, blue = positive)")
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(out / "fig_primary_forest.png", dpi=160)
        fig.savefig(out / "fig_primary_forest.pdf")
        plt.close(fig)

    if PRIMARY_TF in prot.columns and "TACSTD2" in prot.columns:
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
                "NKX2-1 protein (log2)",
                f"{tgt} protein (log2)",
                f"CPTAC-LUAD protein  NKX2-1 vs {tgt}",
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


def write_results_readme(out: Path, summary: dict, table: pd.DataFrame) -> None:
    v = summary["honest_verdict"]
    lines = [
        "# A10 NKX2-1 inverse with TACSTD2 / CLDN4 — public LUAD",
        "",
        f"**Claim:** {summary['claim']}",
        "",
        f"**Honest verdict (TCGA-LUAD primary): `{v['label']}`**",
        "",
        v["statement"],
        "",
        "Filters were not tuned to force an inverse.",
        "",
        "## Primary pairs (BH-FDR within this list)",
        "",
        "| Cohort | Layer | Pair | n | Spearman ρ | p | FDR | 95% CI | Direction |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    prim = table[table["primary_endpoint"]]
    for r in prim.itertuples(index=False):
        ci = (
            f"[{r.spearman_ci95_low:.3f}, {r.spearman_ci95_high:.3f}]"
            if r.spearman_ci95_low is not None
            else "—"
        )
        rho = "—" if r.spearman_rho is None else f"{r.spearman_rho:.3f}"
        pv = "—" if r.spearman_p is None else f"{r.spearman_p:.2e}"
        fdr = "—" if pd.isna(r.fdr_bh_primary) else f"{r.fdr_bh_primary:.3g}"
        lines.append(
            f"| {r.cohort} | {r.layer} | {r.gene_x} vs {r.gene_y} | {r.n} | {rho} | {pv} | {fdr} | {ci} | {r.direction} |"
        )
    lines += [
        "",
        "## What this is not",
        "",
    ]
    for item in summary["cannot_test"]:
        lines.append(f"- {item}")
    lines += [
        "",
        "## Rerun",
        "",
        "```bash",
        "python3 scripts/w200/A10_NKX21/download.py",
        "python3 scripts/w200/A10_NKX21/analyze.py",
        "```",
        "",
        "See `summary.json` and `correlations.csv`.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
