#!/usr/bin/env python3
"""TACSTD2, interferon expression, and IFN-gene dependency in DepMap lung lines."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


IFNA = """
ADAR B2M BATF2 BST2 C1S CASP1 CASP8 CCRL2 CD47 CD74 CMPK2 CMTR1 CNP CSF1
CXCL10 CXCL11 DDX60 DHX58 EIF2AK2 ELF1 EPSTI1 GBP2 GBP4 GMPR HELZ2 HERC6
HLA-C IFI27 IFI30 IFI35 IFI44 IFI44L IFIH1 IFIT2 IFIT3 IFITM1 IFITM2
IFITM3 IL15 IL4R IL7 IRF1 IRF2 IRF7 IRF9 ISG15 ISG20 LAMP3 LAP3 LGALS3BP
LPAR6 LY6E MOV10 MVB12A MX1 NCOA7 NMI NUB1 OAS1 OASL OGFR PARP12 PARP14
PARP9 PLSCR1 PNPT1 PROCR PSMA3 PSMB8 PSMB9 PSME1 PSME2 RIPK2 RNF31 RSAD2
RTP4 SAMD9 SAMD9L SELL SLC25A28 SP110 STAT2 TAP1 TDRD7 TENT5A TMEM140
TRAFD1 TRIM14 TRIM21 TRIM25 TRIM26 TRIM5 TXNIP UBA7 UBE2L6 USP18 WARS1
""".split()

IFNG = """
ADAR APOL6 ARID5B ARL4A AUTS2 B2M BANK1 BATF2 BPGM BST2 BTG1 C1R C1S CASP1
CASP3 CASP4 CASP7 CASP8 CCL2 CCL5 CCL7 CD274 CD38 CD40 CD69 CD74 CD86 CDKN1A
CFB CFH CIITA CMKLR1 CMPK2 CMTR1 CSF2RB CXCL10 CXCL11 CXCL9 DDX60 DHX58
EIF2AK2 EIF4E3 EPSTI1 FAS FCGR1A FGL2 FPR1 GBP4 GBP6 GCH1 GPR18 GZMA HELZ2
HERC6 HIF1A HLA-A HLA-B HLA-DMA HLA-DQA1 HLA-DRB1 HLA-G ICAM1 IDO1 IFI27
IFI30 IFI35 IFI44 IFI44L IFIH1 IFIT1 IFIT2 IFIT3 IFITM2 IFITM3 IFNAR2
IL10RA IL15 IL15RA IL18BP IL2RB IL4R IL6 IL7 IRF1 IRF2 IRF4 IRF5 IRF7 IRF8
IRF9 ISG15 ISG20 ISOC1 ITGB7 JAK2 KLRK1 LAP3 LATS2 LCP2 LGALS3BP LY6E
LYSMD2 MARCHF1 MT2A MTHFD2 MVP MX1 MX2 MYD88 NAMPT NCOA3 NFKB1 NFKBIA
NLRC5 NMI NOD1 NUP93 OAS2 OAS3 OASL OGFR P2RY14 PARP12 PARP14 PDE4B PELI1
PFKP PIM1 PLA2G4A PLSCR1 PML PNP PNPT1 PSMA2 PSMA3 PSMB10 PSMB2 PSMB8 PSMB9
PSME1 PSME2 PTGS2 PTPN1 PTPN2 PTPN6 RAPGEF6 RBCK1 RIGI RIPK1 RIPK2 RNF213
RNF31 RSAD2 RTP4 SAMD9L SAMHD1 SECTM1 SELP SERPING1 SLAMF7 SLC25A28 SOCS1
SOCS3 SOD2 SP110 SPPL2A SRI SSPN ST3GAL5 ST8SIA4 STAT1 STAT2 STAT3 STAT4 TAP1
TAPBP TDRD7 TMT1B TNFAIP2 TNFAIP3 TNFAIP6 TNFSF10 TOR1B TRAFD1 TRIM14 TRIM21
TRIM25 TRIM26 TXNIP UBE2L6 UPP1 USP18 VAMP5 VAMP8 VCAM1 WARS1 XAF1 XCL1 ZBP1
ZNFX1
""".split()

EXPECTED_MD5 = {
    "Model.csv": "675210d17675f3517b0ce39a3c274f16",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv":
        "71794802b750ce77c422dad0720a40af",
    "CRISPRGeneEffect.csv": "6edf7ade09b9b34199210b559d4745d3",
}

SOURCE_URLS = {
    "Model.csv": "https://ndownloader.figshare.com/files/51065297",
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv":
        "https://ndownloader.figshare.com/files/51065489",
    "CRISPRGeneEffect.csv": "https://ndownloader.figshare.com/files/51064667",
}


def md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(2**20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gene_symbol(column: str) -> str:
    return re.sub(r" \(\d+\)$", "", column)


def selected_matrix(path: Path, genes: set[str]) -> tuple[pd.DataFrame, set[str]]:
    header = pd.read_csv(path, nrows=0).columns
    id_col = header[0]
    columns = [c for c in header[1:] if gene_symbol(c) in genes]
    frame = pd.read_csv(path, usecols=[id_col, *columns])
    frame = frame.rename(columns={id_col: "ModelID", **{c: gene_symbol(c) for c in columns}})
    if frame["ModelID"].duplicated().any():
        raise ValueError(f"Duplicate model IDs in {path.name}")
    return frame.set_index("ModelID"), {gene_symbol(c) for c in columns}


def score_signature(expression: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, pd.Series]:
    available = [g for g in genes if g in expression and expression[g].std() > 0]
    z = expression[available].apply(stats.zscore, nan_policy="omit")
    z_score = z.mean(axis=1)
    rank_score = expression[available].rank(pct=True).mean(axis=1)
    return z_score, rank_score


def spearman(x: pd.Series, y: pd.Series) -> dict[str, float | int]:
    pair = pd.concat([x, y], axis=1).dropna()
    if len(pair) < 4 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
        return {"n": len(pair), "rho": np.nan, "p": np.nan}
    result = stats.spearmanr(pair.iloc[:, 0], pair.iloc[:, 1])
    return {"n": len(pair), "rho": float(result.statistic), "p": float(result.pvalue)}


def bootstrap_spearman(
    x: pd.Series, y: pd.Series, rng: np.random.Generator, iterations: int = 5000
) -> tuple[float, float]:
    pair = pd.concat([x, y], axis=1).dropna().to_numpy()
    estimates = np.empty(iterations)
    for i in range(iterations):
        sample = pair[rng.integers(0, len(pair), len(pair))]
        estimates[i] = stats.spearmanr(sample[:, 0], sample[:, 1]).statistic
    return tuple(np.nanpercentile(estimates, [2.5, 97.5]))


def partial_spearman(
    x: pd.Series, y: pd.Series, categories: pd.Series
) -> dict[str, float | int]:
    frame = pd.concat([x.rename("x"), y.rename("y"), categories.rename("category")], axis=1)
    frame = frame.dropna()
    ranked = frame[["x", "y"]].rank()
    design = pd.get_dummies(frame["category"], drop_first=True, dtype=float)
    design.insert(0, "intercept", 1.0)
    matrix = design.to_numpy()
    x_resid = ranked["x"].to_numpy() - matrix @ np.linalg.lstsq(
        matrix, ranked["x"].to_numpy(), rcond=None
    )[0]
    y_resid = ranked["y"].to_numpy() - matrix @ np.linalg.lstsq(
        matrix, ranked["y"].to_numpy(), rcond=None
    )[0]
    result = stats.pearsonr(x_resid, y_resid)
    return {"n": len(frame), "rho": float(result.statistic), "p": float(result.pvalue)}


def bh(frame: pd.DataFrame, p_col: str = "p", q_col: str = "q") -> pd.DataFrame:
    frame = frame.copy()
    frame[q_col] = np.nan
    valid = frame[p_col].notna()
    if valid.any():
        frame.loc[valid, q_col] = multipletests(
            frame.loc[valid, p_col], method="fdr_bh"
        )[1]
    return frame


def fmt_p(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.2e}" if value < 0.001 else f"{value:.3f}"


def make_scatter(
    models: pd.DataFrame, table: pd.DataFrame, output: Path, x_col: str, title: str
) -> None:
    sns.set_theme(style="whitegrid", context="talk")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
    for axis, (score, label) in zip(
        axes, [("ifna_z", "IFN-alpha score"), ("ifng_z", "IFN-gamma score")]
    ):
        data = models[[x_col, score, "DepmapModelType"]].dropna()
        sns.scatterplot(
            data=data, x=x_col, y=score, hue="DepmapModelType", legend=False,
            alpha=0.72, s=42, ax=axis,
        )
        sns.regplot(data=data, x=x_col, y=score, scatter=False, color="black", ax=axis)
        row = table.loc[table["signature"].eq(score)].iloc[0]
        axis.set_title(f"{label}\nSpearman rho={row.rho:.2f}, q={fmt_p(row.q)}")
        axis.set_xlabel("TACSTD2 log2(TPM+1)" if x_col == "tacstd2_expression"
                        else "TACSTD2 Chronos gene effect")
        axis.set_ylabel("Mean within-lung gene z-score")
    fig.suptitle(title, fontsize=18)
    fig.savefig(output, dpi=200)
    plt.close(fig)


def make_candidate_plot(candidates: pd.DataFrame, output: Path) -> None:
    data = candidates[candidates["predictor"].eq("tacstd2_expression")].copy()
    data["minus_log10_q"] = -np.log10(data["q"].clip(lower=1e-300))
    sns.set_theme(style="whitegrid", context="talk")
    fig, axis = plt.subplots(figsize=(9, 6), constrained_layout=True)
    significant = data["q"] < 0.05
    axis.scatter(
        data.loc[~significant, "rho"], data.loc[~significant, "minus_log10_q"],
        color="#808080", alpha=0.55, s=30, label="q >= 0.05",
    )
    axis.scatter(
        data.loc[significant, "rho"], data.loc[significant, "minus_log10_q"],
        color="#c23b23", alpha=0.8, s=40, label="q < 0.05",
    )
    for _, row in data.sort_values("q").head(12).iterrows():
        axis.annotate(row["gene"], (row["rho"], row["minus_log10_q"]), fontsize=8)
    axis.axvline(0, color="black", linewidth=1)
    axis.axhline(-np.log10(0.05), color="#c23b23", linestyle="--", linewidth=1)
    axis.set_xlabel("Spearman rho: TACSTD2 expression vs gene effect")
    axis.set_ylabel("-log10(BH q)")
    axis.set_title("IFN Hallmark gene dependency associations")
    axis.legend(frameon=False)
    fig.savefig(output, dpi=200)
    plt.close(fig)


def write_report(
    output: Path,
    qc: dict,
    primary: pd.DataFrame,
    dependency: pd.DataFrame,
    adjusted: pd.DataFrame,
    subgroup: pd.DataFrame,
    candidates: pd.DataFrame,
    candidate_subgroups: pd.DataFrame,
) -> None:
    def rows_as_markdown(frame: pd.DataFrame, columns: list[str]) -> str:
        shown = frame[columns].copy()
        for col in ["rho", "rho_adjusted", "ci_low", "ci_high"]:
            if col in shown:
                shown[col] = shown[col].map(lambda x: "NA" if pd.isna(x) else f"{x:.3f}")
        for col in ["p", "q", "p_adjusted", "q_adjusted"]:
            if col in shown:
                shown[col] = shown[col].map(fmt_p)
        return shown.to_markdown(index=False)

    top_expr = candidates[candidates["predictor"].eq("tacstd2_expression")].nsmallest(15, "q")
    overlap = len(set(IFNA) & set(IFNG))
    expression_candidates = candidates[
        candidates["predictor"].eq("tacstd2_expression")
    ]
    n_sig = int((expression_candidates["q"] < 0.05).sum())
    n_adjusted_sig = int((expression_candidates["q_adjusted"] < 0.05).sum())
    dependency_candidates = candidates[
        candidates["predictor"].eq("tacstd2_gene_effect")
    ]
    n_dependency_sig = int((dependency_candidates["q"] < 0.05).sum())
    subtype_caution = (
        "- The pooled adjusted candidate hit was not assumed to be histology-general; "
        "exploratory subtype estimates are reported below."
        if n_adjusted_sig
        else "- No pooled adjusted candidate hit required subtype follow-up."
    )
    lines = [
        "# DepMap lung TACSTD2 vs interferon analysis",
        "",
        "## Bottom line",
        "",
    ]
    for _, row in primary.iterrows():
        direction = "higher" if row.rho > 0 else "lower"
        label = "IFN-alpha" if row.signature == "ifna_z" else "IFN-gamma"
        lines.append(
            f"- TACSTD2 RNA was associated with **{direction} {label} score** "
            f"(rho={row.rho:.3f}, 95% bootstrap CI {row.ci_low:.3f} to "
            f"{row.ci_high:.3f}, BH q={fmt_p(row.q)}, n={int(row.n)})."
        )
    lines += [
        f"- In the candidate CRISPR screen, {n_sig} of "
        f"{len(expression_candidates)} IFN-gene associations with TACSTD2 RNA had "
        f"q<0.05; {n_adjusted_sig} remained at q<0.05 after model-type adjustment. "
        "This is an association screen, not evidence that TACSTD2 controls those genes.",
        f"- No IFN-gene effect correlated with TACSTD2 gene effect at q<0.05 "
        f"({n_dependency_sig}/{len(dependency_candidates)}).",
        subtype_caution,
        "- The analysis is observational and cross-sectional. Histology, lineage state, "
        "culture conditions, and screen quality can create correlations. No causal claim "
        "is supported.",
        "",
        "## Data and cohort",
        "",
        "- DepMap Public 24Q4 (stable Figshare archive; release DOI "
        "[10.25452/figshare.plus.27993248.v1]"
        "(https://doi.org/10.25452/figshare.plus.27993248.v1)).",
        "- Primary cohort: `OncotreeLineage == Lung`, `ModelType == Cell Line`, "
        "excluding `OncotreePrimaryDisease == Non-Cancerous`.",
        f"- Metadata lung cancer cell lines: {qc['metadata_lung_cancer_lines']}; "
        f"RNA models analyzed: {qc['expression_models']}; CRISPR-overlap models: "
        f"{qc['crispr_overlap_models']}.",
        "- RNA values are DepMap log2(TPM+1). CRISPR values are Chronos gene effects; "
        "more negative values indicate stronger dependency.",
        "",
        "## Prespecified expression tests",
        "",
        rows_as_markdown(
            primary, ["signature", "n", "rho", "ci_low", "ci_high", "p", "q"]
        ),
        "",
        "Each Hallmark score is the mean of gene-wise z-scores calculated within the "
        "lung cohort. This measures relative IFN-like transcriptional state, not IFN "
        "protein or pathway activation. The alpha and gamma sets overlap by "
        f"{overlap} genes, so the tests are not independent.",
        "",
        "## TACSTD2 dependency vs IFN scores",
        "",
        rows_as_markdown(
            dependency, ["signature", "n", "rho", "ci_low", "ci_high", "p", "q"]
        ),
        "",
        "A positive rho here means higher IFN score accompanies a less-negative "
        "(weaker) TACSTD2 dependency; a negative rho means stronger TACSTD2 dependency.",
        "",
        "## Histology sensitivity analyses",
        "",
        "Partial Spearman correlations residualize ranked variables on DepMap model type "
        "(rare types with <10 RNA-profiled models pooled as `Other`).",
        "",
        rows_as_markdown(adjusted, ["signature", "n", "rho", "p", "q"]),
        "",
        "Unadjusted major-subtype estimates:",
        "",
        rows_as_markdown(subgroup, ["cohort", "signature", "n", "rho", "p", "q"]),
        "",
        "These are sensitivity analyses, not extra discovery tests. Differences between "
        "subtypes may reflect small samples; no formal interaction test was prespecified.",
        "",
        "## IFN-gene CRISPR association screen",
        "",
        "For every available unique gene in the union of the two Hallmark sets, its "
        "Chronos gene effect was correlated with TACSTD2 RNA. BH correction is across "
        "all candidate genes for that predictor. Positive rho means high-TACSTD2 lines "
        "are less dependent on the gene; negative rho means they are more dependent. "
        "Adjusted columns are a sensitivity analysis residualizing ranks on model type.",
        "",
        rows_as_markdown(
            top_expr,
            ["gene", "sets", "n", "rho", "q", "rho_adjusted", "q_adjusted"],
        ),
        "",
        "Because pooled adjustment does not establish consistency across histologies, "
        "the following post-screen estimates show each adjusted q<0.05 hit within major "
        "subtypes. These estimates are exploratory and their p-values are intentionally "
        "not presented as confirmatory tests:",
        "",
        rows_as_markdown(candidate_subgroups, ["gene", "cohort", "n", "rho"]),
        "",
        "The full table also contains a secondary screen using TACSTD2 gene effect as "
        "the predictor, corrected as a separate family.",
        "",
        "## Robustness and limitations",
        "",
        "- Rank-average versions of both signatures are provided per model for score "
        "sensitivity, but the z-score definition is primary.",
        "- Candidate genes were fixed from MSigDB Hallmark IFN-alpha (97 genes) and "
        "IFN-gamma (200 genes); they were not selected from these data.",
        "- Hallmark sets include broad antigen-presentation, proteasome, apoptosis, and "
        "signaling genes, not only IFN-specific effectors.",
        "- Gene-effect co-variation can arise from screen/library effects and general "
        "fitness biology. Correlation does not establish synthetic lethality.",
        "- Cell lines omit immune/stromal compartments and do not model drug response. "
        "The results should not be extrapolated to patients or ADC efficacy.",
        "- No mutation, copy-number, media, growth-rate, or batch covariates were modeled. "
        "Histology adjustment is limited and cannot remove all confounding.",
        "- Current DepMap releases are newer than 24Q4, but 24Q4 was used because it has "
        "stable public URLs and file checksums. Release-specific replication is needed.",
        "- DepMap has no IFN-treated profiles. GEO IFN experiments are mostly A549, which "
        "is TACSTD2-low in DepMap, and cannot explain the basal cross-line correlation.",
        "",
        "## Reproduction",
        "",
        "See `README.md`; all source checksums are verified before analysis. Random "
        "bootstrap seed: 20260816.",
        "",
        "## Files",
        "",
        "- `tables/model_scores.csv`: model-level metadata, TACSTD2, and signatures",
        "- `tables/primary_correlations.csv`: primary expression tests",
        "- `tables/tacstd2_dependency_correlations.csv`: TACSTD2 dependency tests",
        "- `tables/histology_adjusted_correlations.csv`: partial Spearman sensitivity",
        "- `tables/subgroup_correlations.csv`: major subtype sensitivity",
        "- `tables/ifn_gene_crispr_correlations.csv`: full raw and model-type-adjusted "
        "candidate CRISPR screen",
        "- `tables/top_candidate_subgroup_correlations.csv`: exploratory subtype "
        "estimates for adjusted candidate hits",
        "- `tables/ifn_treated_contrasts.csv`: GEO IFN-versus-control TACSTD2 and ISG tests",
        "- `tables/ifn_treated_inventory.csv`: included and excluded treatment datasets",
        "- `figures/`: scatter, candidate-screen, and IFN-treatment plots",
        "- `qc_summary.json`: counts, coverage, checksums, and software-independent inputs",
    ]
    output.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).parent / "data")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    parser.add_argument("--skip-checksums", action="store_true")
    args = parser.parse_args()
    tables = args.out_dir / "tables"
    figures = args.out_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    checksums = {}
    for filename, expected in EXPECTED_MD5.items():
        path = args.data_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}; follow README.md")
        observed = md5(path)
        checksums[filename] = observed
        if not args.skip_checksums and observed != expected:
            raise ValueError(f"Checksum mismatch for {filename}: {observed} != {expected}")

    metadata = pd.read_csv(args.data_dir / "Model.csv").set_index("ModelID")
    cohort = metadata[
        metadata["OncotreeLineage"].eq("Lung")
        & metadata["ModelType"].eq("Cell Line")
        & ~metadata["OncotreePrimaryDisease"].eq("Non-Cancerous")
    ].copy()

    requested = set(IFNA) | set(IFNG) | {"TACSTD2"}
    expression, expression_genes = selected_matrix(
        args.data_dir / "OmicsExpressionProteinCodingGenesTPMLogp1.csv", requested
    )
    expression = expression.loc[expression.index.intersection(cohort.index)]
    if "TACSTD2" not in expression:
        raise ValueError("TACSTD2 absent from expression matrix")

    models = cohort.loc[expression.index, [
        "CellLineName", "DepmapModelType", "OncotreePrimaryDisease", "OncotreeSubtype"
    ]].copy()
    models["tacstd2_expression"] = expression["TACSTD2"]
    models["ifna_z"], models["ifna_rank"] = score_signature(expression, IFNA)
    models["ifng_z"], models["ifng_rank"] = score_signature(expression, IFNG)

    effects, effect_genes = selected_matrix(
        args.data_dir / "CRISPRGeneEffect.csv", requested
    )
    models["tacstd2_gene_effect"] = effects.get("TACSTD2")

    rng = np.random.default_rng(20260816)
    primary_rows = []
    for signature in ["ifna_z", "ifng_z"]:
        result = spearman(models["tacstd2_expression"], models[signature])
        result["ci_low"], result["ci_high"] = bootstrap_spearman(
            models["tacstd2_expression"], models[signature], rng
        )
        primary_rows.append({"signature": signature, **result})
    primary = bh(pd.DataFrame(primary_rows))

    dependency_rows = []
    for signature in ["ifna_z", "ifng_z"]:
        result = spearman(models["tacstd2_gene_effect"], models[signature])
        result["ci_low"], result["ci_high"] = bootstrap_spearman(
            models["tacstd2_gene_effect"], models[signature], rng
        )
        dependency_rows.append({"signature": signature, **result})
    dependency = bh(pd.DataFrame(dependency_rows))

    type_counts = models["DepmapModelType"].value_counts()
    collapsed_type = models["DepmapModelType"].where(
        models["DepmapModelType"].map(type_counts) >= 10, "Other"
    )
    adjusted_rows = []
    for signature in ["ifna_z", "ifng_z"]:
        adjusted_rows.append({
            "signature": signature,
            **partial_spearman(
                models["tacstd2_expression"], models[signature], collapsed_type
            ),
        })
    adjusted = bh(pd.DataFrame(adjusted_rows))

    subgroup_rows = []
    for model_type, count in type_counts.items():
        if count < 20:
            continue
        selected = models["DepmapModelType"].eq(model_type)
        for signature in ["ifna_z", "ifng_z"]:
            subgroup_rows.append({
                "cohort": model_type,
                "signature": signature,
                **spearman(
                    models.loc[selected, "tacstd2_expression"],
                    models.loc[selected, signature],
                ),
            })
    subgroup = bh(pd.DataFrame(subgroup_rows))

    candidate_rows = []
    for predictor in ["tacstd2_expression", "tacstd2_gene_effect"]:
        for gene in sorted((set(IFNA) | set(IFNG)) & effect_genes):
            result = spearman(models[predictor], effects[gene])
            adjusted_result = partial_spearman(
                models[predictor], effects[gene], collapsed_type
            )
            gene_sets = "+".join(
                name for name, members in [("IFNA", IFNA), ("IFNG", IFNG)] if gene in members
            )
            candidate_rows.append({
                "predictor": predictor,
                "gene": gene,
                "sets": gene_sets,
                **result,
                "n_adjusted": adjusted_result["n"],
                "rho_adjusted": adjusted_result["rho"],
                "p_adjusted": adjusted_result["p"],
            })
    candidate_families = []
    for predictor in ["tacstd2_expression", "tacstd2_gene_effect"]:
        family = pd.DataFrame(candidate_rows).loc[lambda x: x.predictor.eq(predictor)]
        family = bh(family)
        family = bh(family, p_col="p_adjusted", q_col="q_adjusted")
        candidate_families.append(family)
    candidates = pd.concat(candidate_families, ignore_index=True)
    candidates = candidates.sort_values(["predictor", "q", "p", "gene"])

    candidate_subgroup_rows = []
    adjusted_hits = candidates.loc[
        candidates["predictor"].eq("tacstd2_expression")
        & candidates["q_adjusted"].lt(0.05),
        "gene",
    ]
    for gene in adjusted_hits:
        for model_type in ["LUAD", "SCLC", "LUSC"]:
            selected = models["DepmapModelType"].eq(model_type)
            result = spearman(
                models.loc[selected, "tacstd2_expression"], effects[gene]
            )
            candidate_subgroup_rows.append({
                "gene": gene, "cohort": model_type, **result
            })
    candidate_subgroups = pd.DataFrame(
        candidate_subgroup_rows, columns=["gene", "cohort", "n", "rho", "p"]
    )

    models.index.name = "ModelID"
    models.to_csv(tables / "model_scores.csv")
    primary.to_csv(tables / "primary_correlations.csv", index=False)
    dependency.to_csv(tables / "tacstd2_dependency_correlations.csv", index=False)
    adjusted.to_csv(tables / "histology_adjusted_correlations.csv", index=False)
    subgroup.to_csv(tables / "subgroup_correlations.csv", index=False)
    candidates.to_csv(tables / "ifn_gene_crispr_correlations.csv", index=False)
    candidate_subgroups.to_csv(
        tables / "top_candidate_subgroup_correlations.csv", index=False
    )

    make_scatter(
        models, primary, figures / "tacstd2_expression_vs_ifn.png",
        "tacstd2_expression", "TACSTD2 expression and IFN state in DepMap lung lines",
    )
    make_scatter(
        models, dependency, figures / "tacstd2_dependency_vs_ifn.png",
        "tacstd2_gene_effect", "TACSTD2 dependency and IFN state in DepMap lung lines",
    )
    make_candidate_plot(candidates, figures / "ifn_gene_crispr_screen.png")

    qc = {
        "release": "DepMap Public 24Q4",
        "release_doi": "10.25452/figshare.plus.27993248.v1",
        "cohort_definition": (
            "OncotreeLineage == Lung; ModelType == Cell Line; "
            "OncotreePrimaryDisease != Non-Cancerous"
        ),
        "metadata_lung_cancer_lines": len(cohort),
        "expression_models": len(models),
        "crispr_overlap_models": int(models["tacstd2_gene_effect"].notna().sum()),
        "ifna_declared_genes": len(IFNA),
        "ifng_declared_genes": len(IFNG),
        "ifna_expression_genes": len(set(IFNA) & expression_genes),
        "ifng_expression_genes": len(set(IFNG) & expression_genes),
        "ifna_crispr_genes": len(set(IFNA) & effect_genes),
        "ifng_crispr_genes": len(set(IFNG) & effect_genes),
        "ifn_unique_crispr_genes": len((set(IFNA) | set(IFNG)) & effect_genes),
        "ifna_ifng_overlap": len(set(IFNA) & set(IFNG)),
        "checksums_md5": checksums,
        "source_urls": SOURCE_URLS,
        "bootstrap_seed": 20260816,
        "bootstrap_iterations": 5000,
    }
    (args.out_dir / "qc_summary.json").write_text(json.dumps(qc, indent=2) + "\n")
    write_report(
        args.out_dir / "REPORT.md", qc, primary, dependency, adjusted, subgroup,
        candidates, candidate_subgroups,
    )
    print(json.dumps(qc, indent=2))


if __name__ == "__main__":
    main()
