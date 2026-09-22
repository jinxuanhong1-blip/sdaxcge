#!/usr/bin/env python3
"""TCGA: TJ panel genes vs CD8 / TNK scores after keratin adjustment.

Pre-specified genes: CLDN4, CLDN3, CLDN7, OCLN, F11R, CDH1.
Cohorts: locked keratin funnel (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD)
plus LUSC. Primary endpoint = keratin-partial Spearman vs CD8 score
(mean of CD8A, CD8B). Secondary = T/NK score (CD3D/E/G, CD8A, NKG7, GNLY, KLRD1).

Numbers are written by this script from Xena GDC STAR TPM. Never fabricated.
"""

from __future__ import annotations

import csv
import gzip
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from stats import bh_fdr, dl_meta_spearman, partial_spearman, spearman

ROOT = Path(__file__).resolve().parent
CACHE = Path(os.environ.get("TCGA_CACHE", "/tmp/tcga"))
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

GDC = "https://gdc-hub.s3.us-east-1.amazonaws.com/download"
PROBEMAP = CACHE / "gencode.v36.annotation.gtf.gene.probemap"

KERATIN_FUNNEL = ["LUAD", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
IMMUNE_COHORTS = ["LUAD", "LUSC"] + [c for c in KERATIN_FUNNEL if c != "LUAD"]
TJ_PANEL = ["CLDN4", "CLDN3", "CLDN7", "OCLN", "F11R", "CDH1"]
KRT = ["KRT8", "KRT18", "KRT19"]
CD8_GENES = ["CD8A", "CD8B"]
TNK_GENES = ["CD3D", "CD3E", "CD3G", "CD8A", "NKG7", "GNLY", "KLRD1"]
NEEDED = sorted(set(TJ_PANEL + KRT + CD8_GENES + TNK_GENES))


def say(msg: str) -> None:
    print(msg, flush=True)


def load_probemap() -> dict[str, str]:
    id_to_gene = {}
    with PROBEMAP.open() as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            # id, gene
            gid, gene = parts[0], parts[1]
            id_to_gene[gid] = gene.upper()
    return id_to_gene


def is_primary01(barcode: str) -> bool:
    parts = barcode.replace(".", "-").split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def patient_id(barcode: str) -> str:
    parts = barcode.replace(".", "-").split("-")
    return "-".join(parts[:3])


def load_cohort(cohort: str, id_to_gene: dict[str, str]) -> dict[str, np.ndarray]:
    path = CACHE / f"TCGA-{cohort}.star_tpm.tsv.gz"
    if not path.exists():
        raise SystemExit(f"missing {path}")
    say(f"load {cohort}")
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        keep_idx = [i for i, b in enumerate(barcodes) if is_primary01(b)]
        if len(keep_idx) < 20:
            raise RuntimeError(f"{cohort}: only {len(keep_idx)} primary tumors")
        patients = [patient_id(barcodes[i]) for i in keep_idx]
        # average replicate aliquots later
        gene_rows = {}
        for line in handle:
            tab = line.find("\t")
            gid = line[:tab]
            gene = id_to_gene.get(gid)
            if gene is None:
                # sometimes first column is already the symbol
                gene = gid.upper()
            if gene not in NEEDED:
                continue
            vals = line.rstrip("\n").split("\t")
            arr = np.array([float(vals[i + 1]) for i in keep_idx], dtype=float)
            # Xena STAR TPM is already log2(TPM+1) in this hub for these files
            gene_rows[gene] = arr
        missing = [g for g in NEEDED if g not in gene_rows]
        if missing:
            raise RuntimeError(f"{cohort} missing genes: {missing}")
    # average replicates to patient
    uniq = sorted(set(patients))
    idx_by_pat = {p: [] for p in uniq}
    for i, p in enumerate(patients):
        idx_by_pat[p].append(i)
    out = {"_patients": np.array(uniq)}
    for gene, arr in gene_rows.items():
        means = np.empty(len(uniq), dtype=float)
        for j, p in enumerate(uniq):
            means[j] = float(np.mean(arr[idx_by_pat[p]]))
        out[gene] = means
    out["_n"] = len(uniq)
    say(f"  n_patients={out['_n']}")
    return out


def score_mean(frame: dict, genes: list[str]) -> np.ndarray:
    return np.mean(np.vstack([frame[g] for g in genes]), axis=0)


def write_tsv(path: Path, rows: list[dict]):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    id_to_gene = load_probemap()
    frames = {c: load_cohort(c, id_to_gene) for c in IMMUNE_COHORTS}

    rows = []
    for cohort in IMMUNE_COHORTS:
        frame = frames[cohort]
        n = frame["_n"]
        cd8 = score_mean(frame, CD8_GENES)
        tnk = score_mean(frame, TNK_GENES)
        cov = [frame[g] for g in KRT]
        for gene in TJ_PANEL:
            x = frame[gene]
            for yname, y in (("CD8_score", cd8), ("TNK_score", tnk)):
                rho_m, p_m = spearman(x, y)
                rho_p, p_p, n_p = partial_spearman(x, y, cov)
                rows.append(
                    {
                        "cohort": cohort,
                        "in_keratin_funnel": cohort in KERATIN_FUNNEL,
                        "gene": gene,
                        "y": yname,
                        "n": n_p,
                        "rho_marginal": rho_m,
                        "p_marginal": p_m,
                        "rho_partial_krt": rho_p,
                        "p_partial_krt": p_p,
                    }
                )
    # BH within CD8 family across genes × cohorts
    cd8_idx = [i for i, r in enumerate(rows) if r["y"] == "CD8_score"]
    q = bh_fdr([rows[i]["p_partial_krt"] for i in cd8_idx])
    for i, qi in zip(cd8_idx, q):
        rows[i]["q_partial_cd8_family"] = qi
        rows[i]["q_partial_tnk_family"] = ""
    tnk_idx = [i for i, r in enumerate(rows) if r["y"] == "TNK_score"]
    q2 = bh_fdr([rows[i]["p_partial_krt"] for i in tnk_idx])
    for i, qi in zip(tnk_idx, q2):
        rows[i]["q_partial_tnk_family"] = qi
        if "q_partial_cd8_family" not in rows[i]:
            rows[i]["q_partial_cd8_family"] = ""
    write_tsv(TAB / "tcga_cohort_partial.tsv", rows)

    # Meta per gene for CD8 and TNK, funnel+LUSC and funnel-only
    meta_rows = []
    for yname in ("CD8_score", "TNK_score"):
        for cohort_set_name, cohort_set in (
            ("lung_plus_funnel", IMMUNE_COHORTS),
            ("keratin_funnel", KERATIN_FUNNEL),
        ):
            for gene in TJ_PANEL:
                sub = [r for r in rows if r["gene"] == gene and r["y"] == yname and r["cohort"] in cohort_set]
                rhos = [r["rho_partial_krt"] for r in sub]
                ns = [r["n"] for r in sub]
                signs = sum(1 for r in rhos if r < 0)
                meta = dl_meta_spearman(rhos, ns, k_cov=3)
                meta_rows.append(
                    {
                        "y": yname,
                        "cohort_set": cohort_set_name,
                        "gene": gene,
                        "n_cohorts": meta["n_cohorts"],
                        "N": int(sum(ns)),
                        "n_negative": signs,
                        "rho": meta["rho"],
                        "p": meta["p"],
                        "I2": meta["I2"],
                        "ci_lo": meta["ci_lo"],
                        "ci_hi": meta["ci_hi"],
                    }
                )
    write_tsv(TAB / "tcga_meta.tsv", meta_rows)

    # Rank for PPT: lung+funnel, CD8, keratin partial — most negative first
    rank = [r for r in meta_rows if r["y"] == "CD8_score" and r["cohort_set"] == "lung_plus_funnel"]
    rank = sorted(rank, key=lambda r: r["rho"] if math.isfinite(r["rho"]) else 999)
    rank_rows = []
    for i, r in enumerate(rank, 1):
        rank_rows.append({**r, "rank_most_negative": i, "is_lead_point_estimate": i == 1})
    write_tsv(TAB / "tcga_rank_cd8.tsv", rank_rows)

    # Forest
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    genes_ord = [r["gene"] for r in rank][::-1]
    for i, gene in enumerate(genes_ord):
        r = next(x for x in rank if x["gene"] == gene)
        color = "#b45309" if gene == "CLDN4" else "#1f2937"
        ax.errorbar(
            r["rho"],
            i,
            xerr=[[r["rho"] - r["ci_lo"]], [r["ci_hi"] - r["rho"]]],
            fmt="o",
            color=color,
            capsize=3,
        )
        ax.text(
            0.98,
            i,
            f"ρ={r['rho']:.3f}  {r['n_negative']}/{r['n_cohorts']} neg  I²={100*r['I2']:.0f}%",
            transform=ax.get_yaxis_transform(),
            va="center",
            ha="right",
            fontsize=8,
            color=color,
        )
    ax.axvline(0, color="#9ca3af", lw=1)
    ax.set_yticks(range(len(genes_ord)))
    ax.set_yticklabels(genes_ord)
    ax.set_xlabel("Keratin-partial Spearman ρ vs CD8 score")
    ax.set_title("TCGA TJ panel vs immune-cold (lung + keratin funnel)")
    fig.tight_layout()
    fig.savefig(FIG / "tcga_tj_forest.png", dpi=150)
    fig.savefig(FIG / "tcga_tj_forest.pdf")
    plt.close(fig)
    say("wrote TCGA tables + forest")
    for r in rank_rows:
        say(f"  rank{r['rank_most_negative']} {r['gene']} ρ={r['rho']:.3f}")


if __name__ == "__main__":
    main()
