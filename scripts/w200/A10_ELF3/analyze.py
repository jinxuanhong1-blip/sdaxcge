#!/usr/bin/env python3
"""A10: ELF3 vs TACSTD2 / CLDN4 — public ChIP + lung co-expression.

ChIP cannot be tested in lung (no public lung ELF3 ChIP).
Co-expression primary statistic: Spearman ρ.
Primary cohorts: TCGA-LUAD and TCGA-LUSC primary tumors.
Does not tune filters to force a hoped-for sign.
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

PRIMARY_TF = "ELF3"
PRIMARY_TARGETS = ["TACSTD2", "CLDN4"]
CONTEXT_GENES = ["EPCAM", "KRT8", "KRT5", "NKX2-1"]
CLAIM = (
    "ELF3 binds TACSTD2 and CLDN4 in lung (public ChIP) and is co-expressed "
    "with them in public lung tumors."
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


def load_tcga(path: Path, purity: pd.Series | None = None) -> pd.DataFrame:
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


def load_purity(path: Path) -> pd.Series:
    if not path.exists():
        return pd.Series(dtype=float)
    pur = pd.read_csv(path, sep="\t")
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
    for cand in (indir / "Model.csv", Path("/tmp/a10_elf3_data/Model.csv")):
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
    return luad, cl


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


def fmt_ci(lo, hi) -> str:
    if lo is None or hi is None or (isinstance(lo, float) and not np.isfinite(lo)):
        return "—"
    return f"[{lo:.3f}, {hi:.3f}]"


def pair_directions(table: pd.DataFrame, cohort: str) -> dict[str, str]:
    sub = table[
        (table["cohort"] == cohort)
        & (table["primary_endpoint"])
        & (table["gene_x"] == PRIMARY_TF)
    ]
    out = {}
    for tgt in PRIMARY_TARGETS:
        row = sub[sub["gene_y"] == tgt]
        out[tgt] = str(row.iloc[0]["direction"]) if not row.empty else "NOT_COMPUTED"
    return out


def coexp_label(dirs: dict[str, str]) -> str:
    vals = [dirs.get(t, "NOT_COMPUTED") for t in PRIMARY_TARGETS]
    if vals == ["POSITIVE", "POSITIVE"]:
        return "POSITIVE"
    if vals == ["INVERSE", "INVERSE"]:
        return "INVERSE"
    if "NOT_COMPUTED" in vals:
        return "NOT_COMPUTED"
    if all(v == "NULL" for v in vals):
        return "NULL"
    return "MIXED"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--in-dir", default="results/w200/A10_ELF3")
    p.add_argument("--out-dir", default="results/w200/A10_ELF3")
    args = p.parse_args()
    indir = Path(args.in_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    chip_verdict = {}
    chip_path = indir / "chip_verdict.json"
    if chip_path.exists():
        chip_verdict = json.loads(chip_path.read_text())

    rows: list[dict] = []
    coverage: dict = {}
    purity = load_purity(indir / "tcga_absolute_purity.txt")

    # ---- TCGA LUAD / LUSC ----
    tumors = {}
    for cohort, fname, key in (
        ("TCGA-LUAD", "tcga_luad_genes.csv", "tcga_luad"),
        ("TCGA-LUSC", "tcga_lusc_genes.csv", "tcga_lusc"),
    ):
        path = indir / fname
        if not path.exists():
            coverage[key] = {"missing_file": str(path)}
            continue
        t = load_tcga(path)
        tumors[cohort] = t
        t.to_csv(out / f"{key}_primary_tumors.csv")
        pur = purity.reindex(t.index)
        pur.to_csv(out / f"{key}_absolute_purity.csv", header=["absolute_purity"])
        coverage[key] = {
            "n_primary_tumors": int(len(t)),
            "n_with_absolute_purity": int(pur.notna().sum()),
            "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g in t.columns],
            "genes_absent": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS + CONTEXT_GENES if g not in t.columns],
        }
        for tgt in PRIMARY_TARGETS:
            add_row(
                rows,
                cohort=cohort,
                layer="RNA_log2TPM",
                gene_x=PRIMARY_TF,
                gene_y=tgt,
                x=t[PRIMARY_TF],
                y=t[tgt],
                primary=True,
                note="Primary lung cohort. One -01 sample per patient.",
                purity=purity,
            )
        add_row(
            rows,
            cohort=cohort,
            layer="RNA_log2TPM",
            gene_x="TACSTD2",
            gene_y="CLDN4",
            x=t["TACSTD2"],
            y=t["CLDN4"],
            primary=False,
            note="Exploratory positive-control pair (epithelial/TJ co-expression).",
            purity=purity,
        )
        for g in CONTEXT_GENES:
            if g in t.columns:
                add_row(
                    rows,
                    cohort=cohort,
                    layer="RNA_log2TPM",
                    gene_x=PRIMARY_TF,
                    gene_y=g,
                    x=t[PRIMARY_TF],
                    y=t[g],
                    primary=False,
                    note="Exploratory context gene (epithelial/basal/lineage).",
                    purity=purity,
                )

    # ---- CPTAC RNA / protein ----
    rna_path = indir / "cptac_luad_rna_genes.csv"
    if rna_path.exists():
        rna = load_cptac(rna_path)
        coverage["cptac_rna"] = {
            "n_tumors": int(len(rna)),
            "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in rna.columns],
            "n_non_na": {
                g: int(rna[g].notna().sum()) if g in rna.columns else 0
                for g in [PRIMARY_TF] + PRIMARY_TARGETS
            },
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
    else:
        rna = pd.DataFrame()

    prot_path = indir / "cptac_luad_protein_genes.csv"
    if prot_path.exists():
        prot = load_cptac(prot_path)
        coverage["cptac_protein"] = {
            "n_tumors": int(len(prot)),
            "genes_present": [g for g in [PRIMARY_TF] + PRIMARY_TARGETS if g in prot.columns],
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
    else:
        prot = pd.DataFrame()

    # ---- DepMap ----
    depmap_expr = indir / "depmap24q4_genes.csv"
    luad_cl = pd.DataFrame()
    lung_cl = pd.DataFrame()
    if depmap_expr.exists():
        luad_cl, lung_cl = load_depmap_lung(indir)
        coverage["depmap_luad_cell_lines"] = {
            "n": int(len(luad_cl)),
            "filter": (
                "OncotreeLineage==Lung AND ModelType==Cell Line AND "
                "(Adenocarcinoma subtype OR OncotreeCode==LUAD)"
            ),
        }
        coverage["depmap_all_lung_cell_lines"] = {
            "n": int(len(lung_cl)),
            "filter": "OncotreeLineage==Lung AND ModelType==Cell Line",
        }
        luad_cl.to_csv(out / "depmap_luad_cell_lines.csv", index=False)
        lung_cl.to_csv(out / "depmap_all_lung_cell_lines.csv", index=False)
        for cohort, df, primary, note in (
            (
                "DepMap24Q4-LUAD",
                luad_cl,
                True,
                "Replication. In vitro LUAD cell lines, not tumors.",
            ),
            (
                "DepMap24Q4-Lung",
                lung_cl,
                False,
                "Exploratory. All DepMap lung cell lines (includes LUSC/SCLC/other).",
            ),
        ):
            for tgt in PRIMARY_TARGETS:
                if PRIMARY_TF in df.columns and tgt in df.columns:
                    add_row(
                        rows,
                        cohort=cohort,
                        layer="RNA_log2TPM",
                        gene_x=PRIMARY_TF,
                        gene_y=tgt,
                        x=df[PRIMARY_TF],
                        y=df[tgt],
                        primary=primary,
                        note=note,
                    )
            if "TACSTD2" in df.columns and "CLDN4" in df.columns:
                add_row(
                    rows,
                    cohort=cohort,
                    layer="RNA_log2TPM",
                    gene_x="TACSTD2",
                    gene_y="CLDN4",
                    x=df["TACSTD2"],
                    y=df["CLDN4"],
                    primary=False,
                    note="Exploratory positive-control pair.",
                )

    table = pd.DataFrame(rows)
    if table.empty:
        raise SystemExit("No correlation rows. Run download.py first.")
    prim = table["primary_endpoint"].astype(bool)
    table.loc[prim, "fdr_bh_primary"] = bh(table.loc[prim, "spearman_p"].fillna(1).tolist())
    table.loc[~prim, "fdr_bh_primary"] = np.nan
    table["direction"] = [
        direction_label(r.spearman_rho, r.fdr_bh_primary if r.primary_endpoint else r.spearman_p)
        for r in table.itertuples(index=False)
    ]
    table.loc[~prim, "direction"] = [
        direction_label(rho, p)
        for rho, p in zip(table.loc[~prim, "spearman_rho"], table.loc[~prim, "spearman_p"])
    ]
    table.to_csv(out / "correlations.csv", index=False)

    luad_dirs = pair_directions(table, "TCGA-LUAD")
    lusc_dirs = pair_directions(table, "TCGA-LUSC")
    luad_lab = coexp_label(luad_dirs)
    lusc_lab = coexp_label(lusc_dirs)

    def rho_txt(cohort: str, tgt: str) -> str:
        sub = table[
            (table["cohort"] == cohort)
            & (table["gene_x"] == PRIMARY_TF)
            & (table["gene_y"] == tgt)
            & (table["primary_endpoint"])
        ]
        if sub.empty or sub.iloc[0]["spearman_rho"] is None:
            return "NA"
        r = sub.iloc[0]
        return f"ρ={r.spearman_rho:.3f} (n={int(r.n)}, FDR={r.fdr_bh_primary:.2g})"

    coexp_statement = (
        f"TCGA-LUAD primary tumors: ELF3 vs TACSTD2 {luad_dirs.get('TACSTD2')} "
        f"({rho_txt('TCGA-LUAD', 'TACSTD2')}); ELF3 vs CLDN4 {luad_dirs.get('CLDN4')} "
        f"({rho_txt('TCGA-LUAD', 'CLDN4')}). "
        f"TCGA-LUSC: TACSTD2 {lusc_dirs.get('TACSTD2')} "
        f"({rho_txt('TCGA-LUSC', 'TACSTD2')}); CLDN4 {lusc_dirs.get('CLDN4')} "
        f"({rho_txt('TCGA-LUSC', 'CLDN4')}). "
        "Positive bulk co-expression is expected for epithelial genes and is "
        "not evidence that ELF3 binds or regulates TACSTD2/CLDN4 in lung."
    )

    chip_lung = bool(chip_verdict.get("lung_elf3_chip_public", False))
    if chip_lung:
        overall = "CHIP_PRESENT_SEE_OVERLAP"
        overall_statement = chip_verdict.get("statement", "")
    else:
        overall = "CHIP_NOT_TESTABLE_IN_LUNG"
        overall_statement = (
            "Public ELF3 ChIP does not exist in lung, so a lung binding claim "
            "at TACSTD2/CLDN4 is not testable. Non-lung ChIP (pancreas/liver/"
            "esophagus/biliary) is mixed and is not a lung result. "
            + coexp_statement
        )

    summary = {
        "claim": CLAIM,
        "honest_verdict": {
            "label": overall,
            "chip_in_lung": "NOT_TESTABLE" if not chip_lung else "PRESENT",
            "coexpression_tcga_luad": luad_lab,
            "coexpression_tcga_lusc": lusc_lab,
            "pair_direction_tcga_luad": luad_dirs,
            "pair_direction_tcga_lusc": lusc_dirs,
            "statement": overall_statement,
            "did_we_tune_filters_to_force_a_sign": False,
        },
        "chip": chip_verdict,
        "coverage": coverage,
        "cannot_test": [
            "ELF3 ChIP-seq in human lung / NSCLC / A549 (none public in ENCODE, ChIP-Atlas, or GEO survey).",
            "Whether ELF3 directly transactivates TACSTD2 or CLDN4 in lung (no lung ChIP, no ELF3 KD/KO lung RNA in this slice).",
            "ICI response or TROP2-ADC outcome.",
        ],
        "methods": {
            "tcga_matrix": "Xena GDC hub TCGA-LUAD/LUSC.star_tpm.tsv.gz, log2(TPM+1), GENCODE v36",
            "cptac_freeze": "data_freeze_v1.2_reorganized/LUAD",
            "depmap_release": "DepMap Public 24Q4",
            "fdr": "BH within pre-specified primary ELF3 vs TACSTD2/CLDN4 rows only",
            "purity": "ABSOLUTE (PanCanAtlas); partial Spearman is secondary",
            "chip": "ENCODE ENCSR770AOR ENCFF080FAU; ChIP-Atlas ELF3.5 + bed05",
        },
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # ---- Figures ----
    if "TCGA-LUAD" in tumors:
        t = tumors["TCGA-LUAD"]
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
        for ax, tgt in zip(axes, PRIMARY_TARGETS):
            sub = table[
                (table["cohort"] == "TCGA-LUAD")
                & (table["gene_x"] == PRIMARY_TF)
                & (table["gene_y"] == tgt)
                & (table["primary_endpoint"])
            ]
            rho = float(sub.iloc[0]["spearman_rho"]) if not sub.empty else None
            pv = float(sub.iloc[0]["spearman_p"]) if not sub.empty else None
            n = int(sub.iloc[0]["n"]) if not sub.empty else 0
            scatter(
                ax,
                t[PRIMARY_TF],
                t[tgt],
                "ELF3 log2(TPM+1)",
                f"{tgt} log2(TPM+1)",
                f"TCGA-LUAD  ELF3 vs {tgt}",
                rho,
                pv,
                n,
            )
        fig.tight_layout()
        fig.savefig(out / "fig_tcga_luad_scatter.png", dpi=160)
        fig.savefig(out / "fig_tcga_luad_scatter.pdf")
        plt.close(fig)

    if "TCGA-LUSC" in tumors:
        t = tumors["TCGA-LUSC"]
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
        for ax, tgt in zip(axes, PRIMARY_TARGETS):
            sub = table[
                (table["cohort"] == "TCGA-LUSC")
                & (table["gene_x"] == PRIMARY_TF)
                & (table["gene_y"] == tgt)
                & (table["primary_endpoint"])
            ]
            rho = float(sub.iloc[0]["spearman_rho"]) if not sub.empty else None
            pv = float(sub.iloc[0]["spearman_p"]) if not sub.empty else None
            n = int(sub.iloc[0]["n"]) if not sub.empty else 0
            scatter(
                ax,
                t[PRIMARY_TF],
                t[tgt],
                "ELF3 log2(TPM+1)",
                f"{tgt} log2(TPM+1)",
                f"TCGA-LUSC  ELF3 vs {tgt}",
                rho,
                pv,
                n,
            )
        fig.tight_layout()
        fig.savefig(out / "fig_tcga_lusc_scatter.png", dpi=160)
        fig.savefig(out / "fig_tcga_lusc_scatter.pdf")
        plt.close(fig)

    prim_tbl = table[table["primary_endpoint"]].copy()
    prim_tbl = prim_tbl.dropna(subset=["spearman_rho"])
    if not prim_tbl.empty:
        labels = [
            f"{r.cohort} {r.layer}\nELF3 vs {r.gene_y} (n={r.n})"
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
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Spearman ρ (95% bootstrap CI)")
        ax.set_title("A10 ELF3 primary co-expression endpoints (red = inverse, blue = positive)")
        ax.invert_yaxis()
        fig.tight_layout()
        fig.savefig(out / "fig_primary_forest.png", dpi=160)
        fig.savefig(out / "fig_primary_forest.pdf")
        plt.close(fig)

    if not prot.empty and PRIMARY_TF in prot.columns:
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
                "ELF3 protein (log2)",
                f"{tgt} protein (log2)",
                f"CPTAC-LUAD protein  ELF3 vs {tgt}",
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
    chip = summary.get("chip") or {}
    lines = [
        "# A10 ELF3 public ChIP / co-expression vs TACSTD2 / CLDN4 in lung",
        "",
        f"**Claim:** {summary['claim']}",
        "",
        f"**Honest verdict: `{v['label']}`**",
        "",
        v["statement"],
        "",
        "Filters were not tuned to force a sign. Motif scanning was not run.",
        "",
        "## 1. Public ELF3 ChIP (lung is absent)",
        "",
        chip.get(
            "statement",
            "No public human lung ELF3 ChIP-seq in ENCODE / ChIP-Atlas / GEO survey.",
        ),
        "",
        "Catalog (all `is_lung=false` except the mouse false-positive GEO text hit): "
        "`chip_experiments.csv`. Locus overlaps: `chip_locus_overlap.csv`.",
        "",
        "Non-lung overlap calls (not a lung result):",
        "",
    ]
    calls = chip.get("non_lung_calls") or {}
    if calls:
        lines.append("| Call | Result |")
        lines.append("| --- | --- |")
        for k, val in calls.items():
            lines.append(f"| {k} | {val} |")
        lines.append("")
    lines += [
        "What the non-lung peaks actually are:",
        "",
        "- **CFPAC-1 (PDAC):** MACS2 bed05 peaks recur near the TACSTD2 coding TSS "
        "(±2 kb) in several replicates, and along the CLDN4 gene span / alt TSS. "
        "A peak at the CLDN4 *coding* TSS (±2 kb) is rare (one CFPAC-1 library).",
        "- **HepG2 (ENCODE optimal IDR ENCFF080FAU):** no peak at TACSTD2 gene ±10 kb. "
        "One peak at chr7:73799582–73799998, which is the CLDN4 *alt* TSS "
        "(~31 kb upstream of the protein-coding TSS), not the coding promoter.",
        "- **ESO-26 and HBDEC2:** no bed05 peak at either locus. HBDEC2 is biliary "
        "epithelium (GSE156165), not lung.",
        "- **HEK293 (GSE280165):** additional 2024 ELF3 ChIP, not lung; peak bed not pulled.",
        "",
        "## 2. Public lung co-expression (primary pairs, BH-FDR within this list)",
        "",
        "| Cohort | Layer | Pair | n | Spearman ρ | p | FDR | 95% CI | Direction |",
        "| --- | --- | --- | ---: | ---: | ---: | ---: | --- | --- |",
    ]
    prim = table[table["primary_endpoint"]]
    for r in prim.itertuples(index=False):
        ci = fmt_ci(r.spearman_ci95_low, r.spearman_ci95_high)
        rho = "—" if r.spearman_rho is None else f"{r.spearman_rho:.3f}"
        pv = "—" if r.spearman_p is None else f"{r.spearman_p:.2e}"
        fdr = "—" if pd.isna(r.fdr_bh_primary) else f"{r.fdr_bh_primary:.3g}"
        lines.append(
            f"| {r.cohort} | {r.layer} | {r.gene_x} vs {r.gene_y} | {r.n} | {rho} | {pv} | {fdr} | {ci} | {r.direction} |"
        )
    expl = table[~table["primary_endpoint"]]
    if not expl.empty:
        lines += [
            "",
            "## Exploratory (not used for the ChIP verdict)",
            "",
            "| Cohort | Pair | n | Spearman ρ | p | Direction |",
            "| --- | --- | ---: | ---: | ---: | --- |",
        ]
        for r in expl.itertuples(index=False):
            rho = "—" if r.spearman_rho is None else f"{r.spearman_rho:.3f}"
            pv = "—" if r.spearman_p is None else f"{r.spearman_p:.2e}"
            lines.append(
                f"| {r.cohort} | {r.gene_x} vs {r.gene_y} | {r.n} | {rho} | {pv} | {r.direction} |"
            )

    cov = summary.get("coverage", {})
    lines += [
        "",
        "## Coverage",
        "",
        f"- TCGA-LUAD primary tumors: n={cov.get('tcga_luad', {}).get('n_primary_tumors')}; "
        f"ABSOLUTE purity for {cov.get('tcga_luad', {}).get('n_with_absolute_purity')}.",
        f"- TCGA-LUSC primary tumors: n={cov.get('tcga_lusc', {}).get('n_primary_tumors')}; "
        f"ABSOLUTE purity for {cov.get('tcga_lusc', {}).get('n_with_absolute_purity')}.",
        f"- CPTAC LUAD RNA: n={cov.get('cptac_rna', {}).get('n_tumors')}.",
        f"- CPTAC LUAD protein non-NA: {cov.get('cptac_protein', {}).get('n_non_na')}.",
        f"- DepMap 24Q4 LUAD cell lines: n={cov.get('depmap_luad_cell_lines', {}).get('n')}.",
        f"- DepMap 24Q4 all lung cell lines (exploratory): n={cov.get('depmap_all_lung_cell_lines', {}).get('n')}.",
        "",
        "## Secondary: purity-adjusted TCGA (does not create lung ChIP)",
        "",
    ]
    for cohort in ("TCGA-LUAD", "TCGA-LUSC"):
        sub = table[(table["cohort"] == cohort) & (table["gene_x"] == PRIMARY_TF)]
        for r in sub.itertuples(index=False):
            if r.purity_partial_rho is None:
                continue
            lines.append(
                f"- {cohort} ELF3 vs {r.gene_y}: partial ρ = {r.purity_partial_rho:.3f}, "
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
        "- A motif presence/absence argument (not run; user asked ChIP/co-expression only).",
        "- Evidence that CFPAC-1 or HepG2 peaks transfer to lung epithelium.",
        "",
        "## Caveats",
        "",
        "1. Bulk tumor RNA mixes epithelium, stroma, and immune cells. ELF3, TACSTD2, and CLDN4 are all epithelial-leaning; a positive ρ can be purity / epithelial fraction.",
        "2. Partial Spearman on ABSOLUTE purity is a check, not a cell-type deconvolution.",
        "3. DepMap lines are not tumors.",
        "4. ChIP-Atlas gene scores are MACS2 peak scores assigned to a 5 kb gene window; a non-zero score is not the same as a coding-TSS peak.",
        "5. CLDN4's Ensembl gene span starts ~31 kb upstream of the protein-coding TSS. Peaks at the alt TSS are not automatically coding-promoter binding.",
        "6. TCGA/CPTAC are not ICI cohorts.",
        "",
        "## Rerun",
        "",
        "```bash",
        "python3 scripts/w200/A10_ELF3/chip.py",
        "python3 scripts/w200/A10_ELF3/download.py",
        "python3 scripts/w200/A10_ELF3/analyze.py",
        "```",
        "",
        "See `summary.json`, `chip_verdict.json`, and `correlations.csv`.",
        "",
    ]
    (out / "README.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
