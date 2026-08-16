#!/usr/bin/env python3
"""Reproducible TACSTD2/Tacstd2 versus immune-signature correlations."""

import csv
import gzip
import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
ANALYSIS = ROOT / "analysis"
FIGURES = ROOT / "figures"
ANALYSIS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)

HUMAN = {
    "CD8_score": ["CD8A", "CD8B", "CCL5", "GZMK", "LCK", "TRAC"],
    "NK_score": ["NKG7", "KLRD1", "GNLY", "PRF1", "CTSW", "NCR1"],
    "pan_immune_score": ["PTPRC", "CD3D", "CD3E", "LST1", "CD74", "HLA-DRA"],
    "epithelial_proxy": ["EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "MUC1"],
}
MOUSE = {
    "CD8_score": ["Cd8a", "Cd8b1", "Ccl5", "Gzmk", "Lck", "Trac"],
    "NK_score": ["Nkg7", "Klrd1", "Prf1", "Ctsw", "Ncr1"],
    "pan_immune_score": ["Ptprc", "Cd3d", "Cd3e", "LST1", "Cd74", "H2-Aa"],
    "epithelial_proxy": ["Epcam", "Krt7", "Krt8", "Krt18", "Krt19", "Muc1"],
}


def read_gmt(path):
    sets = {}
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        sets[fields[0]] = fields[2:]
    return sets


def ssgsea(expression, genes, alpha=0.25):
    """ESTIMATE-compatible single-sample GSEA enrichment statistic."""
    genes = set(genes).intersection(expression.index)
    values = {}
    for sample in expression.columns:
        ranked = expression[sample].dropna().sort_values(ascending=False)
        hits = ranked.index.isin(genes)
        weights = np.abs(ranked.to_numpy()) ** alpha
        hit_walk = np.cumsum(weights * hits) / np.sum(weights * hits)
        miss_walk = np.cumsum(~hits) / np.sum(~hits)
        values[sample] = np.sum(hit_walk - miss_walk)
    return pd.Series(values)


def zscore_rows(frame):
    sd = frame.std(axis=1, ddof=0).replace(0, np.nan)
    return frame.sub(frame.mean(axis=1), axis=0).div(sd, axis=0)


def scores(expression, signatures):
    result, used = {}, {}
    for name, genes in signatures.items():
        present = [g for g in genes if g in expression.index]
        used[name] = present
        result[name] = zscore_rows(expression.loc[present]).mean(axis=0) if present else np.nan
    return pd.DataFrame(result), used


def residualize(y, controls):
    x = np.asarray(controls, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    x = np.column_stack([np.ones(len(y)), x])
    return y - x @ np.linalg.lstsq(x, y, rcond=None)[0]


def correlation(x, y, controls=None):
    frame = pd.concat([x.rename("x"), y.rename("y"), controls], axis=1).dropna()
    rho, p = stats.spearmanr(frame["x"], frame["y"])
    out = {"n": len(frame), "rho": rho, "p_value": p}
    if controls is not None:
        ranked = frame.rank(method="average")
        rx = residualize(ranked["x"].to_numpy(), ranked.iloc[:, 2:].to_numpy())
        ry = residualize(ranked["y"].to_numpy(), ranked.iloc[:, 2:].to_numpy())
        partial, pp = stats.pearsonr(rx, ry)
        out.update({"partial_rho": partial, "partial_p_value": pp})
    return out


def read_geo_matrix(path):
    metadata, rows, columns = {}, [], None
    in_table = False
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                fields = next(csv.reader([line], delimiter="\t"))
                if columns is None:
                    columns = [x.strip('"') for x in fields]
                else:
                    rows.append([x.strip('"') for x in fields])
            elif line.startswith("!Sample_"):
                fields = next(csv.reader([line], delimiter="\t"))
                metadata.setdefault(fields[0], []).append([x.strip('"') for x in fields[1:]])
    frame = pd.DataFrame(rows, columns=columns).set_index("ID_REF").apply(pd.to_numeric)
    return frame, metadata


def read_mouse_annotation(path):
    with gzip.open(path, "rt", errors="replace") as handle:
        while True:
            line = handle.readline()
            if not line:
                raise RuntimeError("GPL annotation table not found")
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                break
        table = pd.read_csv(handle, sep="\t", names=header, dtype=str, low_memory=False)
    return table.set_index("ID")["Gene symbol"]


def collapse_probes(expression, annotation):
    symbols = annotation.reindex(expression.index)
    keep = symbols.notna() & (symbols != "---") & ~symbols.str.contains("///", na=False)
    work = expression.loc[keep].copy()
    work["symbol"] = symbols.loc[keep]
    return work.groupby("symbol").median(numeric_only=True)


def read_gene_matrix(path, id_columns):
    frame = pd.read_csv(path, sep="\t", compression="gzip")
    index = frame.iloc[:, 0].astype(str)
    values = frame.drop(columns=id_columns, errors="ignore")
    values.index = index
    return values.apply(pd.to_numeric, errors="coerce").groupby(level=0).median()


def add_results(results, cohort, target, score_frame, purity_controls, adjustment):
    outcomes = ["CD8_score", "NK_score", "pan_immune_score"]
    if "ESTIMATE_ImmuneScore" in score_frame:
        outcomes.append("ESTIMATE_ImmuneScore")
    for outcome in outcomes:
        raw = correlation(target, score_frame[outcome])
        adjusted = correlation(target, score_frame[outcome], purity_controls)
        results.append({
            "cohort": cohort,
            "outcome": outcome,
            "adjustment": adjustment,
            **raw,
            "partial_rho": adjusted["partial_rho"],
            "partial_p_value": adjusted["partial_p_value"],
            "adjusted_n": adjusted["n"],
        })


def scatter_plot(cohort, target, score_frame, control, control_label):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, outcome in zip(axes, ["CD8_score", "NK_score", "pan_immune_score"]):
        frame = pd.concat([target.rename("target"), score_frame[outcome], control], axis=1).dropna()
        raw = stats.spearmanr(frame["target"], frame[outcome]).statistic
        ranked = frame.rank()
        rx = residualize(ranked["target"].to_numpy(), ranked.iloc[:, 2:].to_numpy())
        ry = residualize(ranked[outcome].to_numpy(), ranked.iloc[:, 2:].to_numpy())
        partial = stats.pearsonr(rx, ry).statistic
        ax.scatter(frame["target"], frame[outcome], s=18, alpha=0.75)
        ax.set_title(f"{outcome}\\nρ={raw:.2f}; partial={partial:.2f}")
        ax.set_xlabel("TACSTD2/Tacstd2")
        ax.set_ylabel("signature score")
    fig.suptitle(f"{cohort}: adjusted for {control_label}", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / f"{cohort}_correlations.png", dpi=180)
    plt.close(fig)


def main():
    results, genes_used = [], {}

    # OncoSG: real pathologist/computational sample-purity values are public.
    onco = pd.read_csv(DATA / "OncoSG_GIS031_selected_expression_zscores.tsv", sep="\t").set_index("sample_id").T
    purity = pd.read_csv(DATA / "OncoSG_GIS031_sample_purity.tsv", sep="\t").set_index("sample_id")["purity"]
    onco_scores, genes_used["OncoSG"] = scores(onco, HUMAN)
    onco_control = purity.rename("purity").to_frame()
    add_results(results, "OncoSG_GIS031", onco.loc["TACSTD2"], onco_scores, onco_control, "reported sample PURITY")
    scatter_plot("OncoSG_GIS031", onco.loc["TACSTD2"], onco_scores, onco_control, "reported sample purity")

    # Open durvalumab neoadjuvant NSCLC cohorts: no measured purity field, so use
    # a prespecified epithelial-content score and label it as a proxy.
    durva_specs = [
        ("GSE248378_post", "GSE248378_Durva_Post_FPKMs.txt.gz", ["gene"]),
        ("GSE253564_pre", "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz", ["gene", "Entrez.ID"]),
    ]
    estimate_sets = read_gmt(ROOT / "provenance" / "ESTIMATE_SI_geneset.gmt")
    for cohort, filename, ids in durva_specs:
        expr = np.log2(read_gene_matrix(DATA / filename, ids) + 1)
        score_frame, genes_used[cohort] = scores(expr, HUMAN)
        score_frame["ESTIMATE_ImmuneScore"] = ssgsea(expr, estimate_sets["ImmuneSignature"])
        genes_used[cohort]["ESTIMATE_ImmuneScore"] = sorted(
            set(expr.index).intersection(estimate_sets["ImmuneSignature"])
        )
        control = score_frame[["epithelial_proxy"]]
        add_results(results, cohort, expr.loc["TACSTD2"], score_frame, control, "epithelial-expression proxy")
        scatter_plot(cohort, expr.loc["TACSTD2"], score_frame, control, "epithelial-expression proxy")

    # GSE76628: this is non-neoplastic flank tissue from athymic nude mice.
    # "Purity" is undefined; adjust for epithelial content and experimental design.
    mouse_probe, metadata = read_geo_matrix(DATA / "GSE76628_series_matrix.txt.gz")
    mouse_expr = np.log2(collapse_probes(mouse_probe, read_mouse_annotation(DATA / "GPL1261-55999.txt.gz")) + 1)
    mouse_scores, genes_used["GSE76628"] = scores(mouse_expr, MOUSE)
    titles = pd.Series(metadata["!Sample_title"][0], index=mouse_expr.columns)
    day = titles.str.extract(r"_(Normal|5_Day|20_Day|60_Day)_", expand=False).fillna("unknown")
    treatment = titles.str.extract(r"_(DC101|G6|NT|Normal)_", expand=False).fillna("unknown")
    design = pd.get_dummies(pd.DataFrame({"day": day, "treatment": treatment}), drop_first=True, dtype=float)
    control = pd.concat([mouse_scores[["epithelial_proxy"]], design], axis=1)
    add_results(
        results, "GSE76628_mouse_flank", mouse_expr.loc["Tacstd2"], mouse_scores,
        control, "epithelial-expression proxy + day + treatment",
    )
    scatter_plot("GSE76628_mouse_flank", mouse_expr.loc["Tacstd2"], mouse_scores, control, "epithelial proxy + design")
    sample_meta = pd.DataFrame({"sample_id": mouse_expr.columns, "title": titles.values, "day": day.values, "treatment": treatment.values})
    sample_meta.to_csv(ANALYSIS / "GSE76628_sample_metadata.tsv", sep="\t", index=False)

    result_frame = pd.DataFrame(results)
    result_frame.to_csv(ANALYSIS / "correlations.tsv", sep="\t", index=False, float_format="%.8g")
    (ANALYSIS / "signature_genes_used.json").write_text(json.dumps(genes_used, indent=2) + "\n")
    print(result_frame.to_string(index=False))


if __name__ == "__main__":
    main()
