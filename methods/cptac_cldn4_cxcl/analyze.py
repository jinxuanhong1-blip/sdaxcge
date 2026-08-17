#!/usr/bin/env python3
"""CLDN4 vs CXCL9 / CXCL10 protein on CPTAC LUAD and LSCC TMT, if present.

Honest pairwise-complete n. Absence after a documented row search is a result.
No TACSTD2 gate. Cohorts kept separate. Protein vs protein only.
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

GENES = {
    "CLDN4": {
        "ensembl": ["ENSG00000189143"],
        "symbols": ["CLDN4", "CLD4"],
    },
    "CXCL9": {
        "ensembl": ["ENSG00000138755"],
        "symbols": ["CXCL9", "MIG", "SCYB9"],
    },
    "CXCL10": {
        "ensembl": ["ENSG00000169245"],
        "symbols": ["CXCL10", "INP10", "IP-10", "IP10", "SCYB10"],
    },
}

PAIRS = [("CLDN4", "CXCL9"), ("CLDN4", "CXCL10")]

FREEZE_FILES = {
    "LUAD": "LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
    "LSCC": "LSCC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
}

LINKED_FILES = {
    "LUAD": "HS_CPTAC_LUAD_proteome_ratio_NArm_TUMOR.cct",
    "LSCC": "HS_CPTAC_LSCC_2020_proteome_ratio_NArm_TUMOR.cct",
}


def load_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    # LinkedOmics CCT first column is attrib_name; drop non-numeric junk rows if any.
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def find_rows(index: pd.Index, prefixes: list[str], symbols: list[str]) -> list[str]:
    hits: list[str] = []
    for i in index:
        s = str(i)
        su = s.upper()
        for p in prefixes:
            if s == p or s.startswith(p + "."):
                hits.append(s)
                break
        else:
            if su in {x.upper() for x in symbols}:
                hits.append(s)
    # unique, preserve order
    seen = set()
    out = []
    for h in hits:
        if h not in seen:
            seen.add(h)
            out.append(h)
    return out


def extract_gene(mat: pd.DataFrame, gene: str) -> dict:
    spec = GENES[gene]
    rows = find_rows(mat.index, spec["ensembl"], spec["symbols"])
    rec = {
        "gene": gene,
        "row_present": bool(rows),
        "n_rows": len(rows),
        "row_ids": ",".join(rows),
        "used_row": rows[0] if rows else "",
        "n_samples": int(mat.shape[1]),
        "n_observed": 0,
        "n_na": int(mat.shape[1]),
        "na_frac": 1.0,
        "note": "gene absent from table (not a row of NAs — never quantified / dropped)",
    }
    if not rows:
        return rec
    if len(rows) > 1:
        rec["note"] = f"multiple rows {rows}; using first"
    row = rows[0]
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    n_obs = int(s.notna().sum())
    n_na = int(s.isna().sum())
    rec.update(
        {
            "n_observed": n_obs,
            "n_na": n_na,
            "na_frac": float(n_na / len(s)) if len(s) else 1.0,
            "note": "row present" if n_obs else "row present but all NA",
        }
    )
    return rec


def gene_series(mat: pd.DataFrame, gene: str) -> tuple[pd.Series, dict]:
    rec = extract_gene(mat, gene)
    if not rec["row_present"]:
        return pd.Series(np.nan, index=mat.columns, name=gene), rec
    row = rec["used_row"]
    s = pd.to_numeric(mat.loc[row], errors="coerce")
    s.name = gene
    return s, rec


def spearman(x: pd.Series, y: pd.Series) -> dict:
    d = pd.concat([x, y], axis=1).dropna()
    n = len(d)
    rec = {
        "n_pairwise": int(n),
        "rho": np.nan,
        "p": np.nan,
        "tested": False,
        "note": "not tested",
    }
    if n == 0:
        rec["note"] = "no pairwise-complete tumors (protein missing)"
        return rec
    if n < 4:
        rec["note"] = f"n={n} < 4; Spearman not computed"
        return rec
    if d.iloc[:, 0].nunique() < 2 or d.iloc[:, 1].nunique() < 2:
        rec["note"] = f"n={n} but no variance"
        return rec
    rho, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    rec.update(
        {
            "rho": float(rho),
            "p": float(p),
            "tested": True,
            "note": "Spearman on pairwise-complete protein values",
        }
    )
    return rec


def analyze_table(mat: pd.DataFrame, cohort: str, source: str) -> tuple[list[dict], list[dict], pd.DataFrame]:
    presence = []
    series = {}
    for gene in GENES:
        s, rec = gene_series(mat, gene)
        rec.update({"cohort": cohort, "source": source})
        presence.append(rec)
        series[gene] = s

    pairs = []
    for a, b in PAIRS:
        rec = spearman(series[a], series[b])
        rec.update(
            {
                "cohort": cohort,
                "source": source,
                "predictor": f"{a}_protein",
                "endpoint": f"{b}_protein",
                "predictor_row_present": presence[[p["gene"] for p in presence].index(a)]["row_present"],
                "endpoint_row_present": presence[[p["gene"] for p in presence].index(b)]["row_present"],
                "n_predictor_observed": presence[[p["gene"] for p in presence].index(a)]["n_observed"],
                "n_endpoint_observed": presence[[p["gene"] for p in presence].index(b)]["n_observed"],
                "n_tumors": int(mat.shape[1]),
            }
        )
        pairs.append(rec)

    core = pd.DataFrame({g: series[g] for g in GENES})
    core["cohort"] = cohort
    core["source"] = source
    return presence, pairs, core


def plot_scatter(cores: dict[str, pd.DataFrame], pairs: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2))
    endpoints = ["CXCL9", "CXCL10"]
    for i, cohort in enumerate(("LUAD", "LSCC")):
        d = cores.get(cohort)
        for j, ep in enumerate(endpoints):
            ax = axes[i, j]
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            row = pairs[(pairs.cohort == cohort) & (pairs.endpoint == f"{ep}_protein")]
            if d is None or d.empty:
                ax.text(0.5, 0.5, "no table", ha="center", va="center", transform=ax.transAxes)
                ax.set_title(f"{cohort}  CLDN4 vs {ep}", fontsize=9)
                continue
            sub = d[["CLDN4", ep]].dropna()
            if sub.empty:
                ax.text(
                    0.5,
                    0.5,
                    f"{ep} protein not quantified\nn_pairwise=0",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    fontsize=9,
                )
            else:
                ax.scatter(sub["CLDN4"], sub[ep], s=18, c="#2c3e50", alpha=0.75, edgecolors="none")
                if len(row) and bool(row.iloc[0]["tested"]):
                    rho = row.iloc[0]["rho"]
                    p = row.iloc[0]["p"]
                    ax.set_title(
                        f"{cohort}  CLDN4 vs {ep} protein\nρ={rho:.3f}  p={p:.3g}  n={len(sub)}",
                        fontsize=9,
                    )
                else:
                    ax.set_title(f"{cohort}  CLDN4 vs {ep} protein\nn={len(sub)} (not tested)", fontsize=9)
                ax.set_xlabel("CLDN4 protein", fontsize=8)
                ax.set_ylabel(f"{ep} protein", fontsize=8)
                continue
            ax.set_title(f"{cohort}  CLDN4 vs {ep} protein", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle(
        "CPTAC public TMT  ·  CLDN4 vs CXCL9/CXCL10 protein if present",
        fontsize=11,
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/cptac_cldn4_cxcl")
    ap.add_argument("--outdir", default="methods/cptac_cldn4_cxcl")
    args = ap.parse_args()
    data = Path(args.data)
    res = Path(args.outdir) / "results"
    res.mkdir(parents=True, exist_ok=True)

    presence_all = []
    pairs_all = []
    cores_freeze: dict[str, pd.DataFrame] = {}
    sample_frames = []

    for cohort, fname in FREEZE_FILES.items():
        path = data / "freeze" / cohort / fname
        if not path.exists():
            presence_all.append(
                {
                    "cohort": cohort,
                    "source": "freeze_v1.2",
                    "gene": "TABLE",
                    "row_present": False,
                    "n_rows": 0,
                    "n_samples": 0,
                    "n_observed": 0,
                    "n_na": 0,
                    "na_frac": 1.0,
                    "note": f"missing file {path}",
                    "row_ids": "",
                }
            )
            continue
        mat = load_matrix(path)
        presence, pairs, core = analyze_table(mat, cohort, "freeze_v1.2")
        presence_all.extend(presence)
        pairs_all.extend(pairs)
        cores_freeze[cohort] = core
        sample_frames.append(core)

    for cohort, fname in LINKED_FILES.items():
        path = data / "linkedomics" / cohort / fname
        if not path.exists():
            presence_all.append(
                {
                    "cohort": cohort,
                    "source": "linkedomics_NArm",
                    "gene": "TABLE",
                    "row_present": False,
                    "n_rows": 0,
                    "n_samples": 0,
                    "n_observed": 0,
                    "n_na": 0,
                    "na_frac": 1.0,
                    "note": f"missing file {path}",
                    "row_ids": "",
                }
            )
            continue
        mat = load_matrix(path)
        presence, pairs, core = analyze_table(mat, cohort, "linkedomics_NArm")
        presence_all.extend(presence)
        pairs_all.extend(pairs)
        sample_frames.append(core)

    presence_df = pd.DataFrame(presence_all)
    pairs_df = pd.DataFrame(pairs_all)
    presence_df.to_csv(res / "presence.tsv", sep="\t", index=False)
    pairs_df.to_csv(res / "spearman.tsv", sep="\t", index=False)

    if sample_frames:
        sample = pd.concat(sample_frames, axis=0)
        sample.index.name = "case_id"
        sample.to_csv(res / "sample_scores.tsv", sep="\t")

    n_rows = []
    for cohort in ("LUAD", "LSCC"):
        for source in ("freeze_v1.2", "linkedomics_NArm"):
            subp = presence_df[(presence_df.cohort == cohort) & (presence_df.source == source)]
            subr = pairs_df[(pairs_df.cohort == cohort) & (pairs_df.source == source)]
            if subp.empty:
                continue
            rec = {"cohort": cohort, "source": source}
            for gene in GENES:
                g = subp[subp.gene == gene]
                if g.empty:
                    rec[f"{gene}_present"] = False
                    rec[f"n_{gene}"] = 0
                else:
                    rec[f"{gene}_present"] = bool(g.iloc[0]["row_present"]) and int(g.iloc[0]["n_observed"]) > 0
                    rec[f"n_{gene}"] = int(g.iloc[0]["n_observed"])
                    rec["n_tumors"] = int(g.iloc[0]["n_samples"])
            for _, pr in subr.iterrows():
                rec[f"n_{pr['predictor'].replace('_protein','')}_{pr['endpoint'].replace('_protein','')}"] = int(
                    pr["n_pairwise"]
                )
            n_rows.append(rec)
    ntab = pd.DataFrame(n_rows)
    ntab.to_csv(res / "n_table.tsv", sep="\t", index=False)

    freeze_pairs = pairs_df[pairs_df.source == "freeze_v1.2"].copy()
    plot_scatter(cores_freeze, freeze_pairs, res / "fig_cldn4_vs_cxcl_protein.png")

    summary = {
        "question": (
            "CLDN4 vs CXCL9/CXCL10 protein if present in CPTAC LUAD and LSCC "
            "public TMT. Honest pairwise-complete n. No TACSTD2 gate."
        ),
        "predictor": "CLDN4 protein (ENSG00000189143)",
        "endpoints": {
            "CXCL9": "CXCL9 protein (ENSG00000138755) if a row exists",
            "CXCL10": "CXCL10 protein (ENSG00000169245) if a row exists",
        },
        "primary_source": "freeze_v1.2 Ensembl gene-abundance (tumor)",
        "secondary_source": "LinkedOmics NArm CCT (gene symbols; high-missingness genes dropped)",
        "tacstd2_gate": "none",
        "presence": presence_df.to_dict(orient="records"),
        "spearman": pairs_df.to_dict(orient="records"),
        "n_table": n_rows,
    }
    (res / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"n": n_rows, "presence": presence_df.to_dict(orient="records"), "spearman": pairs_df.to_dict(orient="records")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
