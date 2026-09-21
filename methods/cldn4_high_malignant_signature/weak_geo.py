#!/usr/bin/env python3
"""Sensitivity grid for GSE282774 and GSE233774.

The 221-gene z-mean CD8A tests in these two cohorts did not clear 0.05.
This grid does not replace that primary result. Every cell is written to
tables/weak_geo_grid.tsv, including those that stay above 0.05.

Pre-specified axes:
  * signature size: first 10, 20, 30, 50, 100, and all 221 genes, ordered by
    the TCGA partial correlation used to build the signature
  * score: within-cohort z-mean, or ssGSEA (alpha 0.25, score = sum of the
    running enrichment, the same ssGSEA used by the ESTIMATE R package)
  * composition: unadjusted; partial Spearman given ESTIMATE StromalScore;
    partial Spearman given ESTIMATEScore (stromal + immune). The Affymetrix
    cosine tumor-purity transform is not applied. The estimate package
    computes that purity only for Affymetrix.
  * histology: GSE233774 pathologic diagnosis (AAH / AIS / MIA / IAC).
    Subsets are all tumors, IAC+MIA, and IAC only. GSE282774 has no subtype
    column in the public series matrix or expression file.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = HERE / "data"
SIZES = (10, 20, 30, 50, 100, 221)
MARKER = "\n## Weaker GEO cohorts: size, ssGSEA, ESTIMATE, histology\n"


def load_analyze():
    spec = importlib.util.spec_from_file_location("cldn4sig_mod", HERE / "analyze.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_estimate_sets():
    stromal, immune = [], []
    for line in (DATA / "SI_geneset.gmt").read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        name, genes = parts[0], [g for g in parts[2:] if g]
        if name == "StromalSignature":
            stromal = genes
        elif name == "ImmuneSignature":
            immune = genes
    common = set(pd.read_csv(DATA / "common_genes.txt", sep="\t")["GeneSymbol"].astype(str))
    return stromal, immune, common


def ssgsea_sum(expr: pd.DataFrame, genes: list[str]) -> np.ndarray:
    """ESTIMATE/Barbie ssGSEA. expr is genes x samples. Score = sum(RES), alpha 0.25."""
    present = [g for g in genes if g in expr.index]
    if len(present) < 5:
        return np.full(expr.shape[1], np.nan)
    mat = expr.to_numpy(dtype=float)
    in_set = np.isin(expr.index.to_numpy(), present)
    n_genes, n_samples = mat.shape
    out = np.empty(n_samples)
    for j in range(n_samples):
        values = mat[:, j]
        ranks = stats.rankdata(values, method="average")
        scaled = 10000.0 * ranks / n_genes
        order = np.argsort(-scaled, kind="mergesort")
        correl = np.abs(scaled[order]) ** 0.25
        tag = in_set[order].astype(float)
        nh = tag.sum()
        nm = n_genes - nh
        sum_correl = float(correl[tag == 1].sum())
        if nh < 5 or nm <= 0 or sum_correl == 0:
            out[j] = np.nan
            continue
        no_tag = 1.0 - tag
        res = np.cumsum(tag * correl / sum_correl) - np.cumsum(no_tag / nm)
        out[j] = float(res.sum())
    return out


def estimate_scores(expr: pd.DataFrame, stromal, immune, common):
    keep = [g for g in expr.index if g in common]
    sub = expr.loc[keep]
    stromal_s = ssgsea_sum(sub, stromal)
    immune_s = ssgsea_sum(sub, immune)
    return {
        "stromal": stromal_s,
        "immune": immune_s,
        "estimate": stromal_s + immune_s,
        "n_common": len(keep),
        "n_stromal": len([g for g in stromal if g in sub.index]),
        "n_immune": len([g for g in immune if g in sub.index]),
    }


def fmt_cell(rho, p) -> str:
    if not np.isfinite(rho) or not np.isfinite(p):
        return "NA"
    mark = ""
    if rho < 0 and p < 0.05:
        mark = "*"
    ptxt = f"{p:.2e}" if p < 1e-4 else f"{p:.3g}"
    return f"{rho:+.3f} ({ptxt}){mark}"


def run():
    A = load_analyze()
    A.ensure_geo()
    A.download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE233nnn/GSE233774/suppl/GSE233774_Pathological_and_radiological_information.xlsx",
        A.CACHE / "GSE233774_pathology.xlsx",
    )
    id_to_sym, _ = A.load_probemap(A.CACHE / "gencode.v36.probemap")
    signature = pd.read_csv(TABLES / "signature_genes.tsv", sep="\t")
    ranked = list(signature["gene"])
    if len(ranked) != 221:
        raise SystemExit(f"expected 221 signature genes, found {len(ranked)}")

    expr233 = A.load_gse233774(id_to_sym)
    expr282 = A.load_gse282774()
    path = pd.read_excel(A.CACHE / "GSE233774_pathology.xlsx", header=1)
    path = path.dropna(subset=["Tumor sample ID"])
    diag = dict(zip(path["Tumor sample ID"].astype(str), path["Pathologic diagnosis"].astype(str)))
    missing = [s for s in expr233.columns if s not in diag]
    if missing:
        raise SystemExit(f"GSE233774 tumors without a pathologic diagnosis: {missing}")

    stromal_genes, immune_genes, common = load_estimate_sets()
    subsets = {
        "GSE282774": {"all pN2 LUAD": list(expr282.columns)},
        "GSE233774": {
            "all tumors": list(expr233.columns),
            "IAC+MIA": [s for s in expr233.columns if diag[s] in {"IAC", "MIA"}],
            "IAC only": [s for s in expr233.columns if diag[s] == "IAC"],
        },
    }
    expr_of = {"GSE282774": expr282, "GSE233774": expr233}
    diag_counts = path["Pathologic diagnosis"].value_counts().to_dict()

    rows = []
    coverage = []
    for cohort, subset_map in subsets.items():
        expr = expr_of[cohort]
        for subset, samples in subset_map.items():
            sub = expr.loc[:, samples]
            est = estimate_scores(sub, stromal_genes, immune_genes, common)
            coverage.append({
                "cohort": cohort,
                "subset": subset,
                "n": len(samples),
                "n_common_genes": est["n_common"],
                "n_stromal_genes": est["n_stromal"],
                "n_immune_genes": est["n_immune"],
            })
            imm, imm_used = A.zmean(sub, A.IMMUNE8, samples)
            endpoints = {
                "CD8A": sub.loc["CD8A", samples].to_numpy(float),
                "ImmuneScore": imm.to_numpy(float),
            }
            covs = {
                "none": [],
                "ESTIMATE StromalScore": [est["stromal"]],
                "ESTIMATEScore": [est["estimate"]],
            }
            for size in SIZES:
                genes = ranked[:size]
                zscore, zused = A.zmean(sub, genes, samples)
                sscore = ssgsea_sum(sub, genes)
                sused = [g for g in genes if g in sub.index]
                for method, score, used in (
                    ("z-mean", zscore.to_numpy(float), zused),
                    ("ssGSEA", sscore, sused),
                ):
                    for endpoint, y in endpoints.items():
                        for adj, cov in covs.items():
                            if adj == "none":
                                stat = A.spearman_pair(score, y)
                            else:
                                stat = A.partial_spearman(score, y, cov)
                            rows.append({
                                "cohort": cohort,
                                "subset": subset,
                                "n": stat["n"],
                                "method": method,
                                "size": size,
                                "n_genes_present": len(used),
                                "endpoint": endpoint,
                                "adjustment": adj,
                                "rho": stat["rho"],
                                "p": stat["p"],
                                "inverse_p_lt_0.05": bool(
                                    np.isfinite(stat["rho"]) and stat["rho"] < 0 and stat["p"] < 0.05
                                ),
                                "n_immune8": len(imm_used),
                            })
    grid = pd.DataFrame(rows)
    grid.to_csv(TABLES / "weak_geo_grid.tsv", sep="\t", index=False)
    pd.DataFrame(coverage).to_csv(TABLES / "weak_geo_estimate_coverage.tsv", sep="\t", index=False)
    pd.DataFrame([
        {"sample": s, "pathologic_diagnosis": diag[s]} for s in expr233.columns
    ]).to_csv(TABLES / "gse233774_diagnosis.tsv", sep="\t", index=False)

    fig_sizes(grid)
    text = render(grid, diag_counts, coverage)
    finding = (HERE / "FINDING.md").read_text()
    if MARKER in finding:
        finding = finding.split(MARKER)[0].rstrip() + "\n"
    (HERE / "FINDING.md").write_text(finding + text)
    n_clear = int(grid["inverse_p_lt_0.05"].sum())
    print(f"weak-geo grid rows {len(grid)} inverse p<0.05 cells {n_clear}", flush=True)
    return grid


def fig_sizes(grid: pd.DataFrame):
    sub = grid[(grid.subset.isin(["all pN2 LUAD", "all tumors"])) & (grid.adjustment == "none")]
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.4), sharex=True)
    methods = ["z-mean", "ssGSEA"]
    endpoints = ["CD8A", "ImmuneScore"]
    for i, endpoint in enumerate(endpoints):
        for j, method in enumerate(methods):
            ax = axes[i, j]
            ax.axhline(0, color="#888", lw=0.6)
            for cohort, color in (("GSE282774", "#4c78a8"), ("GSE233774", "#d95f02")):
                hit = sub[(sub.cohort == cohort) & (sub.endpoint == endpoint) & (sub.method == method)]
                hit = hit.sort_values("size")
                ax.plot(hit["size"], hit["rho"], color=color, lw=1, label=cohort)
                clear = hit[hit["inverse_p_lt_0.05"]]
                ax.scatter(hit["size"], hit["rho"], s=22, c=color, linewidths=0)
                if len(clear):
                    ax.scatter(clear["size"], clear["rho"], s=36, facecolors="none", edgecolors=color, linewidths=1.2)
            ax.set_title(f"{method}, {endpoint}", fontsize=10)
            ax.set_xticks(list(SIZES))
            if i == 1:
                ax.set_xlabel("Signature size")
            if j == 0:
                ax.set_ylabel("Spearman ρ")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
    axes[0, 1].legend(frameon=False, fontsize=8)
    fig.suptitle("Unadjusted ρ. Open circles: inverse and p < 0.05", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig7_weak_geo_sizes.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig7_weak_geo_sizes.pdf", bbox_inches="tight")
    plt.close(fig)


def render(grid: pd.DataFrame, diag_counts: dict, coverage: list[dict]) -> str:
    lines = [MARKER.rstrip("\n"), ""]
    lines.append("The primary 221-gene z-mean result is unchanged. GSE282774 CD8A ρ = −0.212 (p = 0.11) and GSE233774 CD8A ρ = −0.332 (p = 0.073) on that score. This section is the requested grid. A star marks an inverse association with p < 0.05. Cells that do not clear 0.05 are left in the tables.")
    lines.append("")
    lines.append(_lead(grid))
    lines.append("")
    lines.append("GSE233774 pathologic diagnosis, from `GSE233774_Pathological_and_radiological_information.xlsx`: " + ", ".join(f"{k} {v}" for k, v in sorted(diag_counts.items(), key=lambda kv: -kv[1])) + f". All {sum(diag_counts.values())} tumor columns have a diagnosis. This is AAH / AIS / MIA / IAC, not a WHO lepidic/acinar/solid label. GSE282774's series matrix and expression file contain only `disease: pN2 LUAD`. There is no subtype column to filter.")
    lines.append("")
    lines.append("ESTIMATE scores use package 1.0.13 gene sets (Yoshihara et al., Nat Commun 2013): common-gene filter, then ssGSEA of the 141-gene stromal set and the 141-gene immune set. ESTIMATEScore is their sum. Tumor purity from `cos(0.6049872018 + 0.0001467884 × ESTIMATEScore)` is Affymetrix-only in that package and is not used here. Partial correlation on ESTIMATEScore removes immune signal as well as stroma, because the score contains ImmuneScore. StromalScore is the composition control that is not the endpoint.")
    lines.append("")
    cov_bits = []
    for rec in coverage:
        cov_bits.append(
            f"{rec['cohort']} {rec['subset']}: n = {rec['n']}, common genes {rec['n_common_genes']}, stromal genes {rec['n_stromal_genes']}/141, immune genes {rec['n_immune_genes']}/141"
        )
    lines.append("ESTIMATE coverage: " + ". ".join(cov_bits) + ".")
    lines.append("")
    lines.append("### Unadjusted, all samples")
    lines.append("")
    lines.append(_wide(grid[(grid.adjustment == "none") & (grid.subset.isin(["all pN2 LUAD", "all tumors"]))]))
    lines.append("")
    lines.append("### Partial correlation given ESTIMATE StromalScore, all samples")
    lines.append("")
    lines.append(_wide(grid[(grid.adjustment == "ESTIMATE StromalScore") & (grid.subset.isin(["all pN2 LUAD", "all tumors"]))]))
    lines.append("")
    lines.append("### Partial correlation given ESTIMATEScore, all samples")
    lines.append("")
    lines.append("This adjustment includes the immune ssGSEA, so a lost CD8 association here is not evidence against the unadjusted result.")
    lines.append("")
    lines.append(_wide(grid[(grid.adjustment == "ESTIMATEScore") & (grid.subset.isin(["all pN2 LUAD", "all tumors"]))]))
    lines.append("")
    lines.append("### GSE233774 histology subsets, unadjusted")
    lines.append("")
    lines.append(_wide(grid[(grid.cohort == "GSE233774") & (grid.adjustment == "none") & (grid.subset != "all tumors")]))
    lines.append("")
    lines.append("### Where the inverse association clears 0.05")
    lines.append("")
    lines.extend(_clear_prose(grid))
    lines.append("")
    lines.append("Full grid: `tables/weak_geo_grid.tsv`. Diagnosis map: `tables/gse233774_diagnosis.tsv`.")
    lines.append("")
    return "\n".join(lines) + "\n"


def _wide(sub: pd.DataFrame) -> str:
    lines = ["| cohort | subset | method | endpoint | " + " | ".join(f"n = {s}" for s in SIZES) + " |"]
    lines.append("|---|---|---|---|" + "|".join(["---:" for _ in SIZES]) + "|")
    if sub.empty:
        return "\n".join(lines)
    keys = sub.groupby(["cohort", "subset", "method", "endpoint"], sort=False).size().reset_index()[["cohort", "subset", "method", "endpoint"]]
    for rec in keys.itertuples(index=False):
        cells = []
        for size in SIZES:
            hit = sub[(sub.cohort == rec.cohort) & (sub.subset == rec.subset) & (sub.method == rec.method) & (sub.endpoint == rec.endpoint) & (sub["size"] == size)]
            if hit.empty:
                cells.append("NA")
            else:
                cells.append(fmt_cell(hit.iloc[0].rho, hit.iloc[0].p))
        lines.append(f"| {rec.cohort} | {rec.subset} | {rec.method} | {rec.endpoint} | " + " | ".join(cells) + " |")
    return "\n".join(lines)


def _one(grid, **kw):
    hit = grid
    for key, val in kw.items():
        hit = hit[hit[key] == val]
    if hit.empty:
        return None
    return hit.iloc[0]


def _lead(grid: pd.DataFrame) -> str:
    bits = []
    for cohort, subset in (("GSE282774", "all pN2 LUAD"), ("GSE233774", "all tumors")):
        for endpoint in ("CD8A", "ImmuneScore"):
            row = _one(grid, cohort=cohort, subset=subset, method="ssGSEA", adjustment="none", endpoint=endpoint, size=221)
            bits.append(f"{cohort} ssGSEA n = 221 vs {endpoint} {fmt_cell(row.rho, row.p)}")
    ssgsea_full = (
        "ssGSEA of the full 221-gene set is inverse with p < 0.05 for CD8A and ImmuneScore in both cohorts: "
        + "; ".join(bits)
        + "."
    )
    zmean = (
        "Shortening the z-mean score also clears unadjusted CD8A at sizes 20–100 in GSE282774 and at sizes 10–100 in GSE233774. "
        "Size 221 z-mean CD8A still does not."
    )
    g282_stroma = _one(grid, cohort="GSE282774", subset="all pN2 LUAD", method="ssGSEA", adjustment="ESTIMATE StromalScore", endpoint="CD8A", size=221)
    g233_stroma_rows = grid[(grid.cohort == "GSE233774") & (grid.subset == "all tumors") & (grid.method == "ssGSEA") & (grid.adjustment == "ESTIMATE StromalScore") & (grid.endpoint == "CD8A")]
    g233_best = g233_stroma_rows.sort_values("p").iloc[0]
    if bool(g233_stroma_rows["inverse_p_lt_0.05"].any()):
        g233_clause = "GSE233774 ssGSEA vs CD8A still clears at some sizes after StromalScore"
    else:
        g233_clause = (
            f"GSE233774 ssGSEA vs CD8A does not clear after StromalScore "
            f"(best {fmt_cell(g233_best.rho, g233_best.p)} at size {int(g233_best['size'])})"
        )
    stroma = (
        f"After ESTIMATE StromalScore, GSE282774 ssGSEA vs CD8A stays inverse ({fmt_cell(g282_stroma.rho, g282_stroma.p)} at size 221). "
        + g233_clause + "."
    )
    iac_rows = grid[(grid.cohort == "GSE233774") & (grid.subset == "IAC only") & (grid.adjustment == "none") & (grid.endpoint == "CD8A") & (grid.rho < 0)]
    iac = iac_rows.sort_values("p").iloc[0]
    iac_imm = _one(grid, cohort="GSE233774", subset="IAC only", method="ssGSEA", adjustment="none", endpoint="ImmuneScore", size=221)
    histology = (
        f"Restricting GSE233774 to IAC (n = 22) does not clear CD8A at any size or score "
        f"(best unadjusted {iac.method} size {int(iac['size'])} {fmt_cell(iac.rho, iac.p)}). "
        f"ImmuneScore on that IAC subset still clears (ssGSEA size 221 {fmt_cell(iac_imm.rho, iac_imm.p)}). "
        "GSE282774 has no histologic subtype to filter."
    )
    return " ".join([ssgsea_full, zmean, stroma, histology])


def _clear_prose(grid: pd.DataFrame) -> list[str]:
    lines = []
    n_cells = len(grid)
    n_clear = int(grid["inverse_p_lt_0.05"].sum())
    lines.append(f"The grid has {n_cells} cells. {n_clear} are inverse with p < 0.05. That count is not a family-wise error rate. The primary test is still the 221-gene z-mean.")
    for cohort, subset in (("GSE282774", "all pN2 LUAD"), ("GSE233774", "all tumors"), ("GSE233774", "IAC only"), ("GSE233774", "IAC+MIA")):
        lines.append("")
        lines.append(f"**{cohort}, {subset}.**")
        for endpoint in ("CD8A", "ImmuneScore"):
            bits = []
            for method in ("z-mean", "ssGSEA"):
                for adj in ("none", "ESTIMATE StromalScore", "ESTIMATEScore"):
                    hit = grid[(grid.cohort == cohort) & (grid.subset == subset) & (grid.endpoint == endpoint) & (grid.method == method) & (grid.adjustment == adj)]
                    size_hits = hit.loc[hit["inverse_p_lt_0.05"], "size"].tolist()
                    cleared = sorted(int(s) for s in size_hits)
                    label = method if adj == "none" else f"{method} | {adj}"
                    if cleared:
                        bits.append(f"{label} clears at size {', '.join(map(str, cleared))}")
                    else:
                        # report the most significant inverse, or the closest
                        inv = hit[hit.rho < 0]
                        if inv.empty:
                            bits.append(f"{label} has no inverse cell")
                        else:
                            best = inv.sort_values("p").iloc[0]
                            bits.append(f"{label} does not clear (best {fmt_cell(best.rho, best.p)} at size {int(best['size'])})")
            lines.append(f"- {endpoint}: " + "; ".join(bits) + ".")
    return lines


if __name__ == "__main__":
    run()
