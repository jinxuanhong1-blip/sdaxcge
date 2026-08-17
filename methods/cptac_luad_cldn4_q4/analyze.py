#!/usr/bin/env python3
"""LUAD extra: CLDN4 protein Q4 vs ImmuneScore and CD274.

CPTAC LUAD public TMT (freeze v1.2). LSCC Q4 is already reported in
methods/cptac_cldn4_q4 and is not re-run here. No TACSTD2 gate.
Predictor is CLDN4 protein only. Quartiles among tumors with quantified
CLDN4 protein. Primary endpoints: ESTIMATE ImmuneScore and CD274 RNA.
Supporting: CD274 protein (same TMT matrix; honest n, TMT dropout allowed).
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
CD274 = "ENSG00000120217"

PRIMARY = ("ImmuneScore", "CD274_RNA")
SUPPORTING = ("CD274_protein",)
ENDPOINTS = PRIMARY + SUPPORTING


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
    r = (2.0 * U) / (len(a) * len(b)) - 1.0
    rec.update({"U": float(U), "p": float(p), "rank_biserial_q4_gt_q1": float(r)})
    return rec


def cohort_files(data: Path) -> dict[str, Path]:
    return {
        "protein": data
        / "LUAD"
        / "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
        "rna": data / "LUAD" / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
        "pheno": data / "LUAD" / "LUAD_phenotype.txt",
    }


def analyze(data: Path) -> tuple[pd.DataFrame, dict, list[dict], list[dict]]:
    files = cohort_files(data)
    prot = load_matrix(files["protein"])
    rna = load_matrix(files["rna"])
    pheno = load_pheno(files["pheno"])

    samples = sorted(set(prot.columns) & set(rna.columns) & set(pheno.index))
    cldn4, cldn4_row = extract_gene(prot, CLDN4)
    cd274_rna, cd274_rna_row = extract_gene(rna, CD274)
    cd274_prot, cd274_prot_row = extract_gene(prot, CD274)
    imm = immune_score(pheno)

    core = pd.DataFrame(index=samples)
    core["cohort"] = "LUAD"
    core["CLDN4_protein"] = cldn4.reindex(samples)
    core["ImmuneScore"] = imm.reindex(samples)
    core["CD274_RNA"] = cd274_rna.reindex(samples)
    core["CD274_protein"] = cd274_prot.reindex(samples)
    core["quartile"] = quartiles(core["CLDN4_protein"])
    core["tacstd2_gate"] = "none"

    n_tumors = len(samples)
    n_cldn4 = int(core["CLDN4_protein"].notna().sum())
    n_na = int(core["CLDN4_protein"].isna().sum())
    q_counts = core["quartile"].value_counts(dropna=True).to_dict()
    cuts = core.loc[core["CLDN4_protein"].notna(), "CLDN4_protein"]
    q_cuts = cuts.quantile([0.25, 0.5, 0.75]).to_dict() if n_cldn4 else {}

    # Pairwise n for Q4 vs Q1 on each endpoint (honest; CD274 protein may drop).
    q4q1 = core[core["quartile"].isin(["Q1", "Q4"])]
    pairwise_q4q1 = {}
    for ep in ENDPOINTS:
        sub = q4q1[["quartile", ep]].dropna()
        pairwise_q4q1[ep] = {
            "n_q4": int((sub["quartile"] == "Q4").sum()),
            "n_q1": int((sub["quartile"] == "Q1").sum()),
        }

    info = {
        "cohort": "LUAD",
        "n_protein_columns": int(prot.shape[1]),
        "n_rna_columns": int(rna.shape[1]),
        "n_pheno_rows": int(pheno.shape[0]),
        "n_intersect_tumors": n_tumors,
        "cldn4_protein_row": cldn4_row,
        "cd274_rna_row": cd274_rna_row,
        "cd274_protein_row": cd274_prot_row,
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
        "q4_vs_q1_pairwise_n": pairwise_q4q1,
        "n_cd274_protein_quantified_all_tumors": int(core["CD274_protein"].notna().sum()),
        "n_cd274_protein_and_cldn4": int(
            core[["CLDN4_protein", "CD274_protein"]].dropna().shape[0]
        ),
    }

    mw_rows = []
    sp_rows = []
    for ep in ENDPOINTS:
        q4 = core.loc[core["quartile"] == "Q4", ep]
        q1 = core.loc[core["quartile"] == "Q1", ep]
        r = mannwhitney(q4, q1)
        r.update(
            {
                "cohort": "LUAD",
                "predictor": "CLDN4_protein",
                "endpoint": ep,
                "role": "primary" if ep in PRIMARY else "supporting",
                "test": "MWU_Q4_vs_Q1",
                "tacstd2_gate": "none",
            }
        )
        mw_rows.append(r)
        s = spearman(core["CLDN4_protein"], core[ep])
        s.update(
            {
                "cohort": "LUAD",
                "predictor": "CLDN4_protein",
                "endpoint": ep,
                "role": "primary" if ep in PRIMARY else "supporting",
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


def plot_boxes(core: pd.DataFrame, mw: pd.DataFrame, path: Path) -> None:
    titles = {
        "ImmuneScore": "ESTIMATE ImmuneScore",
        "CD274_RNA": "CD274 RNA (PD-L1)",
        "CD274_protein": "CD274 protein (PD-L1)",
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.6))
    sub = core[core["quartile"].isin(["Q1", "Q4"])].copy()
    for ax, ep in zip(axes, ENDPOINTS):
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
        row = mw[mw.endpoint == ep].iloc[0]
        p = row["p"]
        ax.set_title(f"LUAD  {titles[ep]}\nMWU p={p:.3g}", fontsize=9)
        ax.set_ylabel(ep, fontsize=8)
        style(ax)
    fig.suptitle(
        "LUAD extra  ·  CLDN4 protein Q4 vs Q1  ·  ImmuneScore / CD274  ·  no TACSTD2 gate",
        fontsize=10,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cptac_luad_cldn4_q4")
    ap.add_argument("--outdir", default="methods/cptac_luad_cldn4_q4")
    args = ap.parse_args()
    data = Path(args.data)
    out = Path(args.outdir)
    res = out / "results"
    res.mkdir(parents=True, exist_ok=True)

    core, info, mw_rows, sp_rows = analyze(data)
    core.index.name = "case_id"
    core.to_csv(res / "sample_scores.tsv", sep="\t")

    mw = pd.DataFrame(mw_rows)
    sp = pd.DataFrame(sp_rows)
    mw.to_csv(res / "q4_vs_q1.tsv", sep="\t", index=False)
    sp.to_csv(res / "spearman.tsv", sep="\t", index=False)

    ntab = pd.DataFrame(
        [
            {
                "cohort": "LUAD",
                "n_intersect_tumors": info["n_intersect_tumors"],
                "n_cldn4_protein": info["n_cldn4_protein_quantified"],
                "n_cldn4_protein_na": info["n_cldn4_protein_na"],
                "n_Q1": info["quartile_counts"].get("Q1", 0),
                "n_Q2": info["quartile_counts"].get("Q2", 0),
                "n_Q3": info["quartile_counts"].get("Q3", 0),
                "n_Q4": info["quartile_counts"].get("Q4", 0),
                "n_ImmuneScore_among_cldn4": info["endpoint_non_na_among_cldn4"]["ImmuneScore"],
                "n_CD274_RNA_among_cldn4": info["endpoint_non_na_among_cldn4"]["CD274_RNA"],
                "n_CD274_protein_among_cldn4": info["endpoint_non_na_among_cldn4"]["CD274_protein"],
                "n_Q4_Q1_ImmuneScore": "{n_q4} vs {n_q1}".format(
                    **info["q4_vs_q1_pairwise_n"]["ImmuneScore"]
                ),
                "n_Q4_Q1_CD274_RNA": "{n_q4} vs {n_q1}".format(
                    **info["q4_vs_q1_pairwise_n"]["CD274_RNA"]
                ),
                "n_Q4_Q1_CD274_protein": "{n_q4} vs {n_q1}".format(
                    **info["q4_vs_q1_pairwise_n"]["CD274_protein"]
                ),
                "tacstd2_gate": "none",
                "lscc_in_this_slice": "no (already reported in methods/cptac_cldn4_q4)",
            }
        ]
    )
    ntab.to_csv(res / "n_table.tsv", sep="\t", index=False)

    plot_boxes(core, mw, res / "fig_q4_vs_q1.png")

    summary = {
        "question": (
            "LUAD extra: CLDN4 protein Q4 vs Q1 versus ImmuneScore and CD274 "
            "in CPTAC LUAD public TMT. LSCC Q4 already reported. No TACSTD2 gate."
        ),
        "predictor": "CLDN4 protein (ENSG00000189143)",
        "tacstd2_gate": "none",
        "lscc": "not re-run; see methods/cptac_cldn4_q4",
        "quartile_rule": "pd.qcut(rank(method='first'), 4) among CLDN4-protein-quantified tumors",
        "endpoints": {
            "ImmuneScore": "ESTIMATE_ImmuneScore from freeze phenotype (RNA-derived); primary",
            "CD274_RNA": "tumor RNA ENSG00000120217 (log2 RSEM coding UQ 1500); primary",
            "CD274_protein": "same TMT matrix ENSG00000120217; supporting; honest n",
        },
        "cohort": info,
        "q4_vs_q1": mw.to_dict(orient="records"),
        "spearman": sp.to_dict(orient="records"),
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"n": ntab.to_dict(orient="records"), "q4_vs_q1": mw.to_dict(orient="records"), "spearman": sp.to_dict(orient="records")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
