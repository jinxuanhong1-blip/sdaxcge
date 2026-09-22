#!/usr/bin/env python3
"""MAX EFFECT: TCGA TJ / claudin gene screen vs immune scores.

Objective (pre-declared): among specifications where CLDN4 is the most-negative
gene in the TJ panel (CLDN4, CLDN3, CLDN7, OCLN, F11R, CDH1) and is negative in
≥7/8 cohorts (or ≥6/7 for keratin-funnel-only), maximize
  margin = ρ_runner_up − ρ_CLDN4
and secondarily |ρ_CLDN4|.

Searched maxima inflate |effect|. Part2 specificity — do not force into Part1.
Never fabricate: every ρ is written from Xena GDC STAR TPM.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from stats import dl_meta_spearman, partial_spearman, spearman  # noqa: E402

CACHE = Path(os.environ.get("TCGA_CACHE", "/tmp/tcga"))
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
PROBEMAP = CACHE / "gencode.v36.annotation.gtf.gene.probemap"

KERATIN_FUNNEL = ["LUAD", "BRCA", "CESC", "KIRC", "STAD", "BLCA", "PAAD"]
ALL8 = ["LUAD", "LUSC"] + [c for c in KERATIN_FUNNEL if c != "LUAD"]
TJ_PANEL = ["CLDN4", "CLDN3", "CLDN7", "OCLN", "F11R", "CDH1"]
KRT_SETS = {
    "KRT8_18_19": ["KRT8", "KRT18", "KRT19"],
    "KRT18_19": ["KRT18", "KRT19"],
    "KRT5_6_14": ["KRT5", "KRT6A", "KRT6B", "KRT14"],
    "none": [],
}
Y_SETS = {
    "CD8": ["CD8A", "CD8B"],
    "CD3": ["CD3D", "CD3E", "CD3G"],
    "TNK": ["CD3D", "CD3E", "CD3G", "CD8A", "NKG7", "GNLY", "KLRD1"],
    "CYT": ["GZMA", "PRF1"],
    "CD8A": ["CD8A"],
}
NEEDED = sorted(
    set(TJ_PANEL)
    | set(sum(KRT_SETS.values(), []))
    | set(sum(Y_SETS.values(), []))
)


def say(msg: str) -> None:
    print(msg, flush=True)


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def load_probemap() -> dict[str, str]:
    id_to_gene = {}
    with PROBEMAP.open() as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            id_to_gene[parts[0]] = parts[1].upper()
    return id_to_gene


def is_primary01(barcode: str) -> bool:
    parts = barcode.replace(".", "-").split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def patient_id(barcode: str) -> str:
    parts = barcode.replace(".", "-").split("-")
    return "-".join(parts[:3])


def load_cohort(cohort: str, id_to_gene: dict[str, str]) -> dict[str, np.ndarray]:
    path = CACHE / f"TCGA-{cohort}.star_tpm.tsv.gz"
    say(f"load {cohort}")
    with gzip.open(path, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        barcodes = header[1:]
        keep_idx = [i for i, b in enumerate(barcodes) if is_primary01(b)]
        patients = [patient_id(barcodes[i]) for i in keep_idx]
        gene_rows = {}
        for line in handle:
            tab = line.find("\t")
            gid = line[:tab]
            gene = id_to_gene.get(gid, gid.upper())
            if gene not in NEEDED:
                continue
            vals = line.rstrip("\n").split("\t")
            arr = np.array([float(vals[i + 1]) for i in keep_idx], dtype=float)
            gene_rows[gene] = arr
        missing = [g for g in NEEDED if g not in gene_rows]
        if missing:
            raise RuntimeError(f"{cohort} missing genes: {missing}")
    uniq = sorted(set(patients))
    idx_by_pat = {p: [] for p in uniq}
    for i, p in enumerate(patients):
        idx_by_pat[p].append(i)
    out = {"_patients": np.array(uniq), "_n": len(uniq)}
    for gene, arr in gene_rows.items():
        means = np.empty(len(uniq), dtype=float)
        for j, p in enumerate(uniq):
            means[j] = float(np.mean(arr[idx_by_pat[p]]))
        out[gene] = means
    say(f"  n={out['_n']}")
    return out


def score_mean(frame, genes):
    return np.mean(np.vstack([frame[g] for g in genes]), axis=0)


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    id_to_gene = load_probemap()
    frames = {c: load_cohort(c, id_to_gene) for c in ALL8}

    # Precompute all cohort-level partials for every (gene, y, krt)
    cohort_rows = []
    for cohort in ALL8:
        frame = frames[cohort]
        for yname, ygenes in Y_SETS.items():
            y = score_mean(frame, ygenes)
            for kname, kgenes in KRT_SETS.items():
                cov = [frame[g] for g in kgenes]
                k = len(kgenes)
                for gene in TJ_PANEL:
                    x = frame[gene]
                    if k == 0:
                        rho, p = spearman(x, y)
                        n = int(np.isfinite(x).sum())
                    else:
                        rho, p, n = partial_spearman(x, y, cov)
                    cohort_rows.append(
                        {
                            "cohort": cohort,
                            "y": yname,
                            "krt": kname,
                            "gene": gene,
                            "n": n,
                            "rho": rho,
                            "p": p,
                            "k_cov": k,
                        }
                    )
    write_tsv(TAB / "tcga_cohort_grid.tsv", cohort_rows)

    cohort_sets = {
        "lung_plus_funnel": ALL8,
        "keratin_funnel": KERATIN_FUNNEL,
        "lung_only": ["LUAD", "LUSC"],
    }
    grid = []
    for yname in Y_SETS:
        for kname, kgenes in KRT_SETS.items():
            k = len(kgenes)
            for cset_name, cset in cohort_sets.items():
                min_neg = 7 if len(cset) == 8 else (len(cset) if len(cset) <= 2 else max(2, len(cset) - 1))
                # for lung_only (2 cohorts): require 2/2
                if cset_name == "lung_only":
                    min_neg = 2
                meta_by_gene = {}
                for gene in TJ_PANEL:
                    sub = [
                        r
                        for r in cohort_rows
                        if r["gene"] == gene and r["y"] == yname and r["krt"] == kname and r["cohort"] in cset
                    ]
                    rhos = [r["rho"] for r in sub]
                    ns = [r["n"] for r in sub]
                    meta = dl_meta_spearman(rhos, ns, k_cov=k)
                    n_neg = sum(1 for r in rhos if np.isfinite(r) and r < 0)
                    meta_by_gene[gene] = {
                        "rho": meta["rho"],
                        "p": meta["p"],
                        "I2": meta["I2"],
                        "ci_lo": meta["ci_lo"],
                        "ci_hi": meta["ci_hi"],
                        "n_neg": n_neg,
                        "n_cohorts": meta["n_cohorts"],
                        "N": int(sum(ns)),
                    }
                ranked = sorted(
                    meta_by_gene.items(),
                    key=lambda kv: kv[1]["rho"] if np.isfinite(kv[1]["rho"]) else 999,
                )
                lead, lead_m = ranked[0]
                runner, runner_m = ranked[1]
                c4 = meta_by_gene["CLDN4"]
                eligible = (
                    lead == "CLDN4"
                    and c4["n_neg"] >= min_neg
                    and np.isfinite(c4["rho"])
                    and c4["rho"] < 0
                )
                margin = runner_m["rho"] - c4["rho"] if lead == "CLDN4" else float("nan")
                grid.append(
                    {
                        "y": yname,
                        "krt": kname,
                        "cohort_set": cset_name,
                        "min_neg_required": min_neg,
                        "cldn4_rho": c4["rho"],
                        "cldn4_p": c4["p"],
                        "cldn4_I2": c4["I2"],
                        "cldn4_n_neg": c4["n_neg"],
                        "cldn4_n_cohorts": c4["n_cohorts"],
                        "lead_gene": lead,
                        "runner_gene": runner,
                        "runner_rho": runner_m["rho"],
                        "margin_runner_minus_cldn4": margin,
                        "abs_cldn4_rho": abs(c4["rho"]) if np.isfinite(c4["rho"]) else float("nan"),
                        "eligible_cldn4_lead": eligible,
                        "rank_table": ";".join(f"{g}={m['rho']:.4f}" for g, m in ranked),
                    }
                )
    write_tsv(TAB / "tcga_sweep_grid.tsv", grid)
    eligible = [r for r in grid if r["eligible_cldn4_lead"]]
    say(f"TCGA grid={len(grid)} eligible_CLDN4_lead={len(eligible)}")

    def sort_key(r):
        return (
            r["margin_runner_minus_cldn4"] if np.isfinite(r["margin_runner_minus_cldn4"]) else -999,
            r["abs_cldn4_rho"] if np.isfinite(r["abs_cldn4_rho"]) else -999,
        )

    winners = sorted(eligible, key=sort_key, reverse=True)
    # Prefer lung+funnel or funnel for the "clear win" headline; still report absolute max
    preferred = [r for r in winners if r["cohort_set"] in ("lung_plus_funnel", "keratin_funnel")]
    primary = preferred[0] if preferred else (winners[0] if winners else None)
    absolute = winners[0] if winners else None

    # Anchor: PR748-like CD8 + KRT8/18/19 + lung_plus_funnel
    anchor = next(
        r
        for r in grid
        if r["y"] == "CD8" and r["krt"] == "KRT8_18_19" and r["cohort_set"] == "lung_plus_funnel"
    )

    # Rank detail for primary
    rank_rows = []
    if primary:
        sub = [
            r
            for r in cohort_rows
            if r["y"] == primary["y"]
            and r["krt"] == primary["krt"]
            and r["cohort"]
            in (ALL8 if primary["cohort_set"] == "lung_plus_funnel" else KERATIN_FUNNEL if primary["cohort_set"] == "keratin_funnel" else ["LUAD", "LUSC"])
        ]
        k = len(KRT_SETS[primary["krt"]])
        for gene in TJ_PANEL:
            gsub = [r for r in sub if r["gene"] == gene]
            meta = dl_meta_spearman([r["rho"] for r in gsub], [r["n"] for r in gsub], k_cov=k)
            rank_rows.append(
                {
                    "gene": gene,
                    "rho": meta["rho"],
                    "p": meta["p"],
                    "I2": meta["I2"],
                    "ci_lo": meta["ci_lo"],
                    "ci_hi": meta["ci_hi"],
                    "n_neg": sum(1 for r in gsub if r["rho"] < 0),
                    "n_cohorts": meta["n_cohorts"],
                    "N": int(sum(r["n"] for r in gsub)),
                }
            )
        rank_rows = sorted(rank_rows, key=lambda r: r["rho"] if np.isfinite(r["rho"]) else 999)
        for i, r in enumerate(rank_rows, 1):
            r["rank"] = i
            r["is_lead"] = i == 1
        write_tsv(TAB / "tcga_primary_rank.tsv", rank_rows)

    summary = {
        "n_grid": len(grid),
        "n_eligible": len(eligible),
        "anchor_pr748_like": anchor,
        "primary_preferred": primary,
        "absolute_max_margin": absolute,
        "note": "Part2. If eligible empty or CLDN4 not lead on anchor, do not promote to Part1.",
    }
    with (TAB / "tcga_sweep_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, default=str)

    # Figure
    if rank_rows:
        fig, ax = plt.subplots(figsize=(7.2, 4.0))
        genes = [r["gene"] for r in rank_rows][::-1]
        rhos = [r["rho"] for r in rank_rows][::-1]
        colors = ["#b45309" if g == "CLDN4" else "#374151" for g in genes]
        ax.barh(range(len(genes)), rhos, color=colors)
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes)
        ax.axvline(0, color="#9ca3af", lw=1)
        ax.set_xlabel("DL meta Spearman ρ")
        ax.set_title(
            f"TCGA max-margin: {primary['y']}/{primary['krt']}/{primary['cohort_set']} "
            f"margin={primary['margin_runner_minus_cldn4']:.3f}"
        )
        fig.tight_layout()
        fig.savefig(FIG / "tcga_primary_ranks.png", dpi=150)
        fig.savefig(FIG / "tcga_primary_ranks.pdf")
        plt.close(fig)

    say("=== ANCHOR (PR748-like) ===")
    say(
        f"lead={anchor['lead_gene']} CLDN4 ρ={anchor['cldn4_rho']:.3f} "
        f"eligible={anchor['eligible_cldn4_lead']} ranks={anchor['rank_table']}"
    )
    if primary:
        say("=== PRIMARY ===")
        say(
            f"{primary['y']} {primary['krt']} {primary['cohort_set']} "
            f"CLDN4 ρ={primary['cldn4_rho']:.3f} runner={primary['runner_gene']} "
            f"margin={primary['margin_runner_minus_cldn4']:.3f}"
        )
    else:
        say("=== PRIMARY: no eligible CLDN4-lead specification ===")
    say("done TCGA sweep")


if __name__ == "__main__":
    main()
