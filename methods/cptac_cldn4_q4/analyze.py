#!/usr/bin/env python3
"""CLDN4 protein Q4 vs Q1 vs ImmuneScore, GEP18, CD8A.

CPTAC LUAD and LSCC public TMT (freeze v1.2). No TACSTD2 gate.
Predictor is CLDN4 protein only. Quartiles among tumors with quantified
CLDN4 protein. Endpoints are RNA-layer (ESTIMATE ImmuneScore, Ayers GEP18,
CD8A). Honest pairwise-complete n.
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

CLDN4 = "ENSG00000189143"
CD8A = "ENSG00000153563"

# Ayers et al. GEP18 (T-cell inflamed). Same Ensembl IDs as prior CPTAC slices.
GEP18 = [
    "ENSG00000271503",  # CCL5
    "ENSG00000139193",  # CD27
    "ENSG00000120217",  # CD274
    "ENSG00000103855",  # CD276
    "ENSG00000153563",  # CD8A
    "ENSG00000174600",  # CMKLR1
    "ENSG00000138755",  # CXCL9
    "ENSG00000172215",  # CXCR6
    "ENSG00000196735",  # HLA-DQA1
    "ENSG00000196126",  # HLA-DRB1
    "ENSG00000204592",  # HLA-E
    "ENSG00000131203",  # IDO1
    "ENSG00000089692",  # LAG3
    "ENSG00000105374",  # NKG7
    "ENSG00000197646",  # PDCD1LG2
    "ENSG00000205220",  # PSMB10
    "ENSG00000115415",  # STAT1
    "ENSG00000181847",  # TIGIT
]

ENDPOINTS = ("ImmuneScore", "GEP18", "CD8A")


def find_row(index: pd.Index, prefix: str) -> str | None:
    hits = [i for i in index if str(i) == prefix or str(i).startswith(prefix + ".")]
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"multiple rows for {prefix}: {hits}")
    return hits[0]


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    return df.apply(pd.to_numeric, errors="coerce")


def extract_gene(mat: pd.DataFrame, prefix: str) -> tuple[pd.Series, str | None]:
    row = find_row(mat.index, prefix)
    if row is None:
        return pd.Series(np.nan, index=mat.columns, name=prefix), None
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = prefix
    return s, row


def mean_z(mat: pd.DataFrame, prefixes: list[str], min_genes: int = 10) -> tuple[pd.Series, list[str]]:
    rows = []
    used = []
    for p in prefixes:
        rid = find_row(mat.index, p)
        if rid is None:
            continue
        s = pd.to_numeric(mat.loc[rid], errors="coerce")
        if s.notna().sum() < 10:
            continue
        sd = s.std(ddof=0)
        if sd == 0 or not np.isfinite(sd):
            continue
        rows.append((s - s.mean()) / sd)
        used.append(rid)
    if not rows:
        return pd.Series(np.nan, index=mat.columns), used
    z = pd.concat(rows, axis=1)
    n = z.notna().sum(axis=1)
    score = z.mean(axis=1, skipna=True)
    score[n < min_genes] = np.nan
    return score, used


def load_pheno(path: Path) -> pd.DataFrame:
    ph = pd.read_csv(path, sep="\t", index_col=0)
    ph.index = ph.index.astype(str)
    ph = ph[ph.index != "data_type"]
    return ph


def immune_score(pheno: pd.DataFrame) -> pd.Series:
    if "ESTIMATE_ImmuneScore" in pheno.columns:
        s = pd.to_numeric(pheno["ESTIMATE_ImmuneScore"], errors="coerce")
        s.name = "ImmuneScore"
        return s
    hits = [c for c in pheno.columns if "immunescore" in c.lower().replace("_", "")]
    if not hits:
        raise KeyError("ESTIMATE_ImmuneScore not in phenotype")
    s = pd.to_numeric(pheno[hits[0]], errors="coerce")
    s.name = "ImmuneScore"
    return s


def quartiles(s: pd.Series) -> pd.Series:
    """Equal-count Q1–Q4 on non-NA values. rank-first so ties still fill 4 bins."""
    out = pd.Series(np.nan, index=s.index, dtype=object)
    ok = s.dropna()
    if len(ok) < 8:
        return out
    labels = pd.qcut(ok.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    out.loc[ok.index] = labels.astype(str)
    return out


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    if n < 4 or d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        return {"n": int(n), "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    return {"n": int(n), "rho": float(rho), "p": float(p)}


def mannwhitney(a: pd.Series, b: pd.Series) -> dict:
    a = pd.to_numeric(a, errors="coerce").dropna()
    b = pd.to_numeric(b, errors="coerce").dropna()
    rec = {
        "n_q4": int(len(a)),
        "n_q1": int(len(b)),
        "median_q4": float(a.median()) if len(a) else np.nan,
        "median_q1": float(b.median()) if len(b) else np.nan,
        "delta_median_q4_minus_q1": np.nan,
        "U": np.nan,
        "p": np.nan,
        "rank_biserial_q4_gt_q1": np.nan,
    }
    if len(a) and len(b):
        rec["delta_median_q4_minus_q1"] = float(a.median() - b.median())
    if len(a) < 3 or len(b) < 3:
        return rec
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # Positive r means Q4 stochastically greater than Q1.
    r = (2.0 * U) / (len(a) * len(b)) - 1.0
    rec.update({"U": float(U), "p": float(p), "rank_biserial_q4_gt_q1": float(r)})
    return rec


def cohort_files(data: Path, cohort: str) -> dict[str, Path]:
    prefix = cohort
    return {
        "protein": data
        / cohort
        / f"{prefix}_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "rna": data / cohort / f"{prefix}_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "pheno": data / cohort / f"{prefix}_phenotype.txt",
    }


def analyze_cohort(data: Path, cohort: str) -> tuple[pd.DataFrame, dict, list[dict], list[dict]]:
    files = cohort_files(data, cohort)
    prot = load_matrix(files["protein"])
    rna = load_matrix(files["rna"])
    pheno = load_pheno(files["pheno"])

    samples = sorted(set(prot.columns) & set(rna.columns) & set(pheno.index))
    cldn4, cldn4_row = extract_gene(prot, CLDN4)
    cd8a, cd8a_row = extract_gene(rna, CD8A)
    gep, gep_used = mean_z(rna, GEP18, min_genes=10)
    imm = immune_score(pheno)

    core = pd.DataFrame(index=samples)
    core["cohort"] = cohort
    core["CLDN4_protein"] = cldn4.reindex(samples)
    core["ImmuneScore"] = imm.reindex(samples)
    core["GEP18"] = gep.reindex(samples)
    core["CD8A"] = cd8a.reindex(samples)
    core["quartile"] = quartiles(core["CLDN4_protein"])
    # Explicit: TACSTD2 is never used as a gate or covariate.
    core["tacstd2_gate"] = "none"

    n_tumors = len(samples)
    n_cldn4 = int(core["CLDN4_protein"].notna().sum())
    n_na = int(core["CLDN4_protein"].isna().sum())
    q_counts = core["quartile"].value_counts(dropna=True).to_dict()
    cuts = core.loc[core["CLDN4_protein"].notna(), "CLDN4_protein"]
    q_cuts = cuts.quantile([0.25, 0.5, 0.75]).to_dict() if n_cldn4 else {}

    info = {
        "cohort": cohort,
        "n_protein_columns": int(prot.shape[1]),
        "n_rna_columns": int(rna.shape[1]),
        "n_pheno_rows": int(pheno.shape[0]),
        "n_intersect_tumors": n_tumors,
        "cldn4_protein_row": cldn4_row,
        "cd8a_rna_row": cd8a_row,
        "gep18_genes_used": len(gep_used),
        "gep18_genes_requested": len(GEP18),
        "gep18_rows": gep_used,
        "n_cldn4_protein_quantified": n_cldn4,
        "n_cldn4_protein_na": n_na,
        "quartile_counts": {k: int(v) for k, v in q_counts.items()},
        "cldn4_protein_q25": float(q_cuts.get(0.25, np.nan)),
        "cldn4_protein_q50": float(q_cuts.get(0.5, np.nan)),
        "cldn4_protein_q75": float(q_cuts.get(0.75, np.nan)),
        "tacstd2_gate": "none",
        "endpoint_non_na_among_cldn4": {
            ep: int(core.loc[core["CLDN4_protein"].notna(), ep].notna().sum()) for ep in ENDPOINTS
        },
    }

    mw_rows = []
    sp_rows = []
    for ep in ENDPOINTS:
        q4 = core.loc[core["quartile"] == "Q4", ep]
        q1 = core.loc[core["quartile"] == "Q1", ep]
        r = mannwhitney(q4, q1)
        r.update(
            {
                "cohort": cohort,
                "predictor": "CLDN4_protein",
                "endpoint": ep,
                "test": "MWU_Q4_vs_Q1",
                "tacstd2_gate": "none",
            }
        )
        mw_rows.append(r)
        s = spearman(core["CLDN4_protein"], core[ep])
        s.update(
            {
                "cohort": cohort,
                "predictor": "CLDN4_protein",
                "endpoint": ep,
                "test": "spearman_continuous",
                "tacstd2_gate": "none",
            }
        )
        sp_rows.append(s)

    return core, info, mw_rows, sp_rows


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)


def plot_boxes(cores: dict[str, pd.DataFrame], mw: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4))
    titles = {
        "ImmuneScore": "ESTIMATE ImmuneScore",
        "GEP18": "GEP18 (RNA z-mean)",
        "CD8A": "CD8A RNA",
    }
    for i, cohort in enumerate(("LUAD", "LSCC")):
        d = cores[cohort]
        sub = d[d["quartile"].isin(["Q1", "Q4"])].copy()
        for j, ep in enumerate(ENDPOINTS):
            ax = axes[i, j]
            a = sub.loc[sub["quartile"] == "Q1", ep].dropna()
            b = sub.loc[sub["quartile"] == "Q4", ep].dropna()
            ax.boxplot(
                [a, b],
                tick_labels=[f"Q1\nn={len(a)}", f"Q4\nn={len(b)}"],
                widths=0.55,
                patch_artist=True,
                boxprops=dict(facecolor="#d9e8f5", edgecolor="#2c3e50"),
                medianprops=dict(color="#c0392b", linewidth=1.6),
                whiskerprops=dict(color="#2c3e50"),
                capprops=dict(color="#2c3e50"),
                flierprops=dict(marker="o", markersize=3, markerfacecolor="#7f8c8d"),
            )
            row = mw[(mw.cohort == cohort) & (mw.endpoint == ep)].iloc[0]
            p = row["p"]
            ax.set_title(f"{cohort}  {titles[ep]}\nMWU p={p:.3g}", fontsize=9)
            ax.set_ylabel(ep, fontsize=8)
            style(ax)
    fig.suptitle(
        "CLDN4 protein Q4 vs Q1  ·  no TACSTD2 gate  ·  CPTAC public TMT",
        fontsize=11,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cptac_cldn4_q4")
    ap.add_argument("--outdir", default="methods/cptac_cldn4_q4")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    res = out / "results"
    res.mkdir(parents=True, exist_ok=True)

    cores = {}
    infos = []
    mw_all = []
    sp_all = []
    for cohort in ("LUAD", "LSCC"):
        core, info, mw_rows, sp_rows = analyze_cohort(data, cohort)
        cores[cohort] = core
        infos.append(info)
        mw_all.extend(mw_rows)
        sp_all.extend(sp_rows)

    sample = pd.concat(cores.values(), axis=0)
    sample.index.name = "case_id"
    sample.to_csv(res / "sample_scores.tsv", sep="\t")

    mw = pd.DataFrame(mw_all)
    sp = pd.DataFrame(sp_all)
    mw.to_csv(res / "q4_vs_q1.tsv", sep="\t", index=False)
    sp.to_csv(res / "spearman.tsv", sep="\t", index=False)

    n_rows = []
    for info in infos:
        n_rows.append(
            {
                "cohort": info["cohort"],
                "n_intersect_tumors": info["n_intersect_tumors"],
                "n_cldn4_protein": info["n_cldn4_protein_quantified"],
                "n_cldn4_protein_na": info["n_cldn4_protein_na"],
                "n_Q1": info["quartile_counts"].get("Q1", 0),
                "n_Q2": info["quartile_counts"].get("Q2", 0),
                "n_Q3": info["quartile_counts"].get("Q3", 0),
                "n_Q4": info["quartile_counts"].get("Q4", 0),
                "n_ImmuneScore_among_cldn4": info["endpoint_non_na_among_cldn4"]["ImmuneScore"],
                "n_GEP18_among_cldn4": info["endpoint_non_na_among_cldn4"]["GEP18"],
                "n_CD8A_among_cldn4": info["endpoint_non_na_among_cldn4"]["CD8A"],
                "gep18_genes_used": info["gep18_genes_used"],
                "tacstd2_gate": "none",
            }
        )
    ntab = pd.DataFrame(n_rows)
    ntab.to_csv(res / "n_table.tsv", sep="\t", index=False)

    plot_boxes(cores, mw, res / "fig_q4_vs_q1.png")

    summary = {
        "question": (
            "CLDN4 protein Q4 vs Q1 versus ImmuneScore, GEP18, CD8A "
            "in CPTAC LUAD and LSCC public TMT. No TACSTD2 gate."
        ),
        "predictor": "CLDN4 protein (ENSG00000189143)",
        "tacstd2_gate": "none",
        "quartile_rule": "pd.qcut(rank(method='first'), 4) among CLDN4-protein-quantified tumors",
        "endpoints": {
            "ImmuneScore": "ESTIMATE_ImmuneScore from freeze phenotype (RNA-derived)",
            "GEP18": "Ayers 18-gene T-cell inflamed; mean of per-gene z on tumor RNA",
            "CD8A": "tumor RNA ENSG00000153563 (log2 RSEM coding UQ 1500)",
        },
        "cohorts": infos,
        "q4_vs_q1": mw.to_dict(orient="records"),
        "spearman": sp.to_dict(orient="records"),
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"n": n_rows, "q4_vs_q1": mw.to_dict(orient="records")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
