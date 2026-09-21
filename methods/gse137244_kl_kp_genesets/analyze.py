#!/usr/bin/env python3
"""GSE137244 KL vs KP gene-set scores.

Deng et al., Nat Cancer 2021 (PMID 34142094). GEO supplementary FPKM.
Score = mean of log2(FPKM + 1) over genes present in the matrix.
Contrast = mean(5 KL libraries) - mean(5 KP libraries).
Normal lung is tabulated and held out of the test.

The locked single-gene numbers on this file are Cldn4 +5.57 and Tacstd2 +3.24,
both with complete library separation (two-sided Mann–Whitney p = 2/252 = 0.00794).
This script checks those two deltas before writing the new scores.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import urllib.request
from collections import OrderedDict
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

from gene_sets import (
    CLDN4_TACSTD2,
    HALLMARK_IFN_ALPHA,
    HALLMARK_IFN_GAMMA,
    IFN_COMPACT,
    KEGG_CYTOSOLIC_DNA_2019_MOUSE,
    KEGG_NHEJ_2019_MOUSE,
    MRN,
    STING_CORE,
    STING_KINASES,
    TJ_TISMO,
)

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137244/suppl/"
    "GSE137244_counts.fpkm.csv.gz"
)
FPKM_NAME = "GSE137244_counts.fpkm.csv.gz"

# Library order in the GEO matrix. Genotype is taken from Sample_description.
LIBRARIES = [
    ("B6AL10-1-RNA", "GSM4073816", "KP", "B6AL10"),
    ("B6AL10-2-RNA", "GSM4073817", "KP", "B6AL10"),
    ("B6AL10-3-RNA", "GSM4073818", "KP", "B6AL10"),
    ("B6AL10-4-RNA", "GSM4073819", "KP", "B6AL10"),
    ("B6AL10-5-RNA", "GSM4073820", "KP", "B6AL10"),
    ("KL155mix-control-2-RNA", "GSM4073821", "KL", "KL155"),
    ("KL47-1-untreated-1-RNA", "GSM4073822", "KL", "KL47"),
    ("KLC-RNA", "GSM4073823", "KL", "KLC"),
    ("KLD-RNA", "GSM4073824", "KL", "KLD"),
    ("KLE-RNA", "GSM4073825", "KL", "KLE"),
    ("normal-lung-RNA", "GSM4073826", "normal", "normal"),
]

# Directional KL>KP family. q is Benjamini–Hochberg on one-sided Welch p.
PRIMARY = ["Cldn4_Tacstd2", "NHEJ_KEGG", "STING_kinases", "MRN"]
# Epithelial gate: B6AL10-3 is the only tumor library with Epcam log2(FPKM+1) < 4.
EPCAM_MIN = 4.0

# Locked single-gene mean differences on log2(FPKM+1), KL minus KP.
LOCKED = {"Cldn4": 5.57, "Tacstd2": 3.24}


def download_fpkm() -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    dest = DATA / FPKM_NAME
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    urllib.request.urlretrieve(FPKM_URL, dest)
    return dest


def load_fpkm(path: Path) -> tuple[list[str], dict[str, list[float]]]:
    with gzip.open(path, "rt") as handle:
        header = next(handle).strip().strip('"').split(",")
        samples = header[1:]
        genes: dict[str, list[float]] = {}
        for line in handle:
            parts = [x.strip('"') for x in line.strip().split(",")]
            genes[parts[0]] = [float(x) for x in parts[1:]]
    expected = [row[0] for row in LIBRARIES]
    if samples != expected:
        raise SystemExit(f"Unexpected sample order: {samples}")
    return samples, genes


def resolve(requested: list[str], lookup: dict[str, str]) -> tuple[list[str], list[str]]:
    used: list[str] = []
    missing: list[str] = []
    seen: set[str] = set()
    for symbol in requested:
        hit = lookup.get(symbol.lower())
        if hit is None:
            missing.append(symbol)
            continue
        if hit not in seen:
            seen.add(hit)
            used.append(hit)
    return used, missing


def mean(xs: list[float]) -> float:
    return sum(xs) / len(xs)


def mannwhitney_exact(x: list[float], y: list[float]) -> tuple[float, float, int]:
    """Two-sided exact Mann–Whitney. U counts y > x, with ties as 0.5."""

    def u_stat(group_y: list[float], group_x: list[float]) -> float:
        u = 0.0
        for a in group_y:
            for b in group_x:
                if a > b:
                    u += 1.0
                elif a == b:
                    u += 0.5
        return u

    u_obs = u_stat(y, x)
    pooled = list(x) + list(y)
    n_y = len(y)
    center = len(x) * len(y) / 2.0
    dev = abs(u_obs - center)
    extreme = 0
    total = 0
    for comb in combinations(range(len(pooled)), n_y):
        total += 1
        chosen = set(comb)
        ys = [pooled[i] for i in comb]
        xs = [pooled[i] for i in range(len(pooled)) if i not in chosen]
        if abs(u_stat(ys, xs) - center) + 1e-9 >= dev:
            extreme += 1
    return u_obs, extreme / total, total


def bh(names: list[str], pvals: list[float]) -> dict[str, float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: (pvals[i], names[i]))
    q = [1.0] * m
    running = 1.0
    for rank_from_end, idx in enumerate(reversed(order)):
        j = m - rank_from_end  # 1-based rank from the smallest end, walking down
        candidate = pvals[idx] * m / j
        running = min(running, candidate)
        q[idx] = running
    return {names[i]: q[i] for i in range(m)}


def separation(kp: list[float], kl: list[float]) -> str:
    if min(kl) > max(kp):
        return "KL>KP"
    if max(kl) < min(kp):
        return "KP>KL"
    return "overlap"


def welch_greater(kp: list[float], kl: list[float]) -> tuple[float, float]:
    """One-sided Welch t. Alternative is mean(KL) > mean(KP)."""
    result = stats.ttest_ind(kl, kp, equal_var=False, alternative="greater")
    return float(result.statistic), float(result.pvalue)


def perm_greater(values: list[float], kp_idx: list[int], kl_idx: list[int]) -> tuple[float, float, int, int]:
    """One-sided permutation p for mean(KL) - mean(KP) on the listed libraries."""
    kp = [values[i] for i in kp_idx]
    kl = [values[i] for i in kl_idx]
    observed = mean(kl) - mean(kp)
    pooled = kp + kl
    n_kl = len(kl)
    extreme = 0
    total = 0
    for comb in combinations(range(len(pooled)), n_kl):
        total += 1
        chosen = set(comb)
        kl_mean = mean([pooled[i] for i in comb])
        kp_mean = mean([pooled[i] for i in range(len(pooled)) if i not in chosen])
        if kl_mean - kp_mean + 1e-12 >= observed:
            extreme += 1
    return observed, extreme / total, extreme, total


def score_vector(log_expr: dict[str, list[float]], genes: list[str]) -> list[float]:
    n = len(next(iter(log_expr.values())))
    out = []
    for i in range(n):
        out.append(mean([log_expr[g][i] for g in genes]))
    return out


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fmt(x: float, digits: int = 4) -> str:
    return f"{x:.{digits}f}"


def plot_scores(sample_rows: list[dict], contrasts: dict[str, dict]) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    panels = [
        ("Cldn4_Tacstd2", "Cldn4 / Tacstd2"),
        ("NHEJ_KEGG", "NHEJ (KEGG, 13 genes)"),
        ("STING_kinases", "STING kinases (Tbk1, Ikbke)"),
        ("MRN", "MRN (Mre11a, Rad50, Nbn)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 7.2), constrained_layout=True)
    colors = {"KP": "#4C78A8", "KL": "#E45756"}
    for ax, (key, title) in zip(axes.ravel(), panels):
        for geno, xpos in (("KP", 0), ("KL", 1)):
            kept = [r for r in sample_rows if r["genotype"] == geno and r["epcam_gate"] == "keep"]
            dropped = [r for r in sample_rows if r["genotype"] == geno and r["epcam_gate"] == "drop"]
            ax.scatter(
                [xpos] * len(kept),
                [float(r[key]) for r in kept],
                s=46,
                color=colors[geno],
                zorder=3,
            )
            if dropped:
                ax.scatter(
                    [xpos] * len(dropped),
                    [float(r[key]) for r in dropped],
                    s=54,
                    facecolors="none",
                    edgecolors=colors[geno],
                    linewidths=1.4,
                    zorder=4,
                )
            ys = [float(r[key]) for r in sample_rows if r["genotype"] == geno]
            ax.hlines(mean(ys), xpos - 0.18, xpos + 0.18, color="black", lw=1.4, zorder=2)
        c = contrasts[key]
        ax.set_xticks([0, 1], ["KP", "KL"])
        ax.set_xlim(-0.55, 1.55)
        ax.set_title(
            f"{title}\nΔ={c['delta']:+.2f}  {c['separation']}  Welch p={c['p_welch_greater']:.2e}",
            fontsize=10,
        )
        ax.set_ylabel("mean log2(FPKM+1)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    fig.suptitle("GSE137244 KL > KP. Open circle: Epcam-low KP library held out of the gate", fontsize=11)
    fig.savefig(FIGS / "scores_kl_vs_kp.png", dpi=160)
    fig.savefig(FIGS / "scores_kl_vs_kp.pdf")
    plt.close(fig)


def plot_gene_deltas(gene_rows: list[dict]) -> None:
    groups = [
        ("Cldn4_Tacstd2", "Cldn4 / Tacstd2"),
        ("TJ_TISMO", "TJ (7 genes)"),
        ("STING_kinases", "STING kinases"),
        ("STING_core", "STING core"),
        ("MRN", "MRN"),
        ("NHEJ_KEGG", "NHEJ"),
        ("IFN_compact", "IFN compact"),
    ]
    fig, axes = plt.subplots(len(groups), 1, figsize=(8.4, 16.0), constrained_layout=True)
    for ax, (key, title) in zip(axes, groups):
        rows = [r for r in gene_rows if r["set_name"] == key]
        rows = sorted(rows, key=lambda r: float(r["delta_kl_minus_kp"]))
        ys = list(range(len(rows)))
        deltas = [float(r["delta_kl_minus_kp"]) for r in rows]
        colors = ["#E45756" if d >= 0 else "#4C78A8" for d in deltas]
        ax.barh(ys, deltas, color=colors, height=0.72)
        ax.set_yticks(ys, [r["gene"] for r in rows], fontsize=8)
        ax.axvline(0, color="black", lw=0.6)
        ax.set_title(title, loc="left", fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_xlabel("Δ log2(FPKM+1), KL − KP")
    fig.savefig(FIGS / "gene_deltas.png", dpi=160)
    fig.savefig(FIGS / "gene_deltas.pdf")
    plt.close(fig)


def main() -> None:
    path = download_fpkm()
    samples, fpkm = load_fpkm(path)
    lookup = {g.lower(): g for g in fpkm}
    log_expr = {g: [math.log2(v + 1.0) for v in vals] for g, vals in fpkm.items()}
    kp_idx = [i for i, row in enumerate(LIBRARIES) if row[2] == "KP"]
    kl_idx = [i for i, row in enumerate(LIBRARIES) if row[2] == "KL"]

    set_defs: OrderedDict[str, tuple[str, list[str], str]] = OrderedDict(
        [
            ("Cldn4", ("primary_gene", ["Cldn4"], "locked single gene")),
            ("Tacstd2", ("primary_gene", ["Tacstd2"], "locked single gene")),
            (
                "Cldn4_Tacstd2",
                ("primary", CLDN4_TACSTD2, "mean of Cldn4 and Tacstd2"),
            ),
            (
                "TJ_TISMO",
                (
                    "secondary",
                    TJ_TISMO,
                    "TISMO TJ list: Cldn3/4/6/7, Cdh1, F11r, Ocln",
                ),
            ),
            (
                "NHEJ_KEGG",
                (
                    "primary",
                    KEGG_NHEJ_2019_MOUSE,
                    "KEGG 2019 Mouse Non-homologous end-joining",
                ),
            ),
            (
                "Rad50",
                ("primary_gene", ["Rad50"], "strongest single KEGG NHEJ gene"),
            ),
            (
                "MRN",
                ("primary", MRN, "Mre11a, Rad50, Nbn"),
            ),
            (
                "STING_kinases",
                ("primary", STING_KINASES, "Tbk1 and Ikbke"),
            ),
            (
                "STING_core",
                (
                    "secondary",
                    STING_CORE,
                    "Mb21d1, Tmem173, Tbk1, Ikbke, Irf3",
                ),
            ),
            (
                "IFN_compact",
                ("primary", IFN_COMPACT, "mouse epithelial IFN list"),
            ),
            (
                "Hallmark_IFN_alpha",
                (
                    "secondary",
                    HALLMARK_IFN_ALPHA,
                    "MSigDB Hallmark 2020 interferon alpha response",
                ),
            ),
            (
                "Hallmark_IFN_gamma",
                (
                    "secondary",
                    HALLMARK_IFN_GAMMA,
                    "MSigDB Hallmark 2020 interferon gamma response",
                ),
            ),
            (
                "KEGG_cytosolic_DNA",
                (
                    "secondary",
                    KEGG_CYTOSOLIC_DNA_2019_MOUSE,
                    "KEGG 2019 Mouse cytosolic DNA-sensing; CGAS symbol unmatched",
                ),
            ),
        ]
    )

    resolved: dict[str, list[str]] = {}
    missing_map: dict[str, list[str]] = {}
    for name, (_role, requested, _note) in set_defs.items():
        used, missing = resolve(requested, lookup)
        if not used:
            raise SystemExit(f"{name} matched no genes")
        resolved[name] = used
        missing_map[name] = missing

    # Locked single-gene check on this FPKM file.
    for gene, locked in LOCKED.items():
        vals = log_expr[gene]
        delta = mean([vals[i] for i in kl_idx]) - mean([vals[i] for i in kp_idx])
        if abs(delta - locked) > 0.006:
            raise SystemExit(f"{gene} delta {delta:.4f} does not match locked {locked}")

    scores: dict[str, list[float]] = {
        name: score_vector(log_expr, genes) for name, genes in resolved.items()
    }

    contrast_rows = []
    contrast_by_name: dict[str, dict] = {}
    for name, (role, _requested, note) in set_defs.items():
        sc = scores[name]
        kp = [sc[i] for i in kp_idx]
        kl = [sc[i] for i in kl_idx]
        delta = mean(kl) - mean(kp)
        u_obs, p_exact, n_perm = mannwhitney_exact(kp, kl)
        t_stat, p_welch = welch_greater(kp, kl)
        _obs, p_perm, n_ext, n_perm_dir = perm_greater(sc, kp_idx, kl_idx)
        n_pairs = len(kp) * len(kl)
        # rank-biserial on the KL-vs-KP pairwise outcomes. Ties contribute 0.
        n_kl_higher = 0.0
        n_kp_higher = 0.0
        for a in kl:
            for b in kp:
                if a > b:
                    n_kl_higher += 1
                elif a < b:
                    n_kp_higher += 1
        rrb = (n_kl_higher - n_kp_higher) / n_pairs
        row = {
            "set_name": name,
            "role": role,
            "n_genes_used": len(resolved[name]),
            "n_genes_requested": len(set_defs[name][1]),
            "n_genes_unmatched": len(missing_map[name]),
            "mean_kp": mean(kp),
            "mean_kl": mean(kl),
            "delta": delta,
            "separation": separation(kp, kl),
            "U_kl": u_obs,
            "p_mw_two_sided": p_exact,
            "p_welch_greater": p_welch,
            "t_welch": t_stat,
            "p_perm_greater": p_perm,
            "perm_extreme": n_ext,
            "n_permutations": n_perm_dir,
            "rank_biserial_kl": rrb,
            "note": note,
        }
        contrast_rows.append(row)
        contrast_by_name[name] = row

    qmap = bh(PRIMARY, [contrast_by_name[n]["p_welch_greater"] for n in PRIMARY])
    for row in contrast_rows:
        row["q_bh_welch_primary4"] = qmap.get(row["set_name"], "")

    # Joint KL>KP score: mean of tumor-library z-scores for the four directional arms.
    tumor_idx = kp_idx + kl_idx

    def z_tumor(values: list[float]) -> list[float]:
        xs = [values[i] for i in tumor_idx]
        mu = mean(xs)
        sd = math.sqrt(sum((x - mu) ** 2 for x in xs) / (len(xs) - 1))
        if sd == 0:
            sd = 1.0
        out = []
        for i, value in enumerate(values):
            out.append((value - mu) / sd if i in set(tumor_idx) else float("nan"))
        return out

    z_parts = [
        z_tumor(scores["Cldn4"]),
        z_tumor(scores["Tacstd2"]),
        z_tumor(scores["NHEJ_KEGG"]),
        z_tumor(scores["STING_kinases"]),
    ]
    joint = []
    for i in range(len(LIBRARIES)):
        if i not in set(tumor_idx):
            joint.append(float("nan"))
        else:
            joint.append(mean([part[i] for part in z_parts]))
    scores["thesis_joint"] = joint
    kp = [joint[i] for i in kp_idx]
    kl = [joint[i] for i in kl_idx]
    t_stat, p_welch = welch_greater(kp, kl)
    _obs, p_perm, n_ext, n_perm_dir = perm_greater(joint, kp_idx, kl_idx)
    u_obs, p_exact, _n = mannwhitney_exact(kp, kl)
    joint_row = {
        "set_name": "thesis_joint",
        "role": "primary",
        "n_genes_used": "",
        "n_genes_requested": "",
        "n_genes_unmatched": 0,
        "mean_kp": mean(kp),
        "mean_kl": mean(kl),
        "delta": mean(kl) - mean(kp),
        "separation": separation(kp, kl),
        "U_kl": u_obs,
        "p_mw_two_sided": p_exact,
        "p_welch_greater": p_welch,
        "t_welch": t_stat,
        "p_perm_greater": p_perm,
        "perm_extreme": n_ext,
        "n_permutations": n_perm_dir,
        "rank_biserial_kl": 1.0 if min(kl) > max(kp) else "",
        "q_bh_welch_primary4": "",
        "note": "mean z of Cldn4, Tacstd2, NHEJ_KEGG, STING_kinases across the 10 tumor libraries",
    }
    contrast_rows.append(joint_row)
    contrast_by_name["thesis_joint"] = joint_row

    # Epcam gate. Drops tumor libraries with log2(Epcam FPKM+1) < 4.
    epcam = log_expr["Epcam"]
    gate_kp = [i for i in kp_idx if epcam[i] >= EPCAM_MIN]
    gate_kl = [i for i in kl_idx if epcam[i] >= EPCAM_MIN]
    dropped = [LIBRARIES[i][0] for i in kp_idx + kl_idx if epcam[i] < EPCAM_MIN]
    gate_rows = []
    for name in ["Cldn4", "Tacstd2", "Cldn4_Tacstd2", "NHEJ_KEGG", "Rad50", "MRN", "STING_kinases", "STING_core", "thesis_joint"]:
        sc = scores[name]
        kp = [sc[i] for i in gate_kp]
        kl = [sc[i] for i in gate_kl]
        t_stat, p_welch = welch_greater(kp, kl)
        _obs, p_perm, n_ext, n_perm_dir = perm_greater(sc, gate_kp, gate_kl)
        gate_rows.append(
            {
                "set_name": name,
                "n_kp": len(gate_kp),
                "n_kl": len(gate_kl),
                "dropped_libraries": ",".join(dropped),
                "epcam_min_log2": EPCAM_MIN,
                "mean_kp": mean(kp),
                "mean_kl": mean(kl),
                "delta": mean(kl) - mean(kp),
                "separation": separation(kp, kl),
                "t_welch": t_stat,
                "p_welch_greater": p_welch,
                "p_perm_greater": p_perm,
                "perm_extreme": n_ext,
                "n_permutations": n_perm_dir,
            }
        )

    # Per-gene deltas for every set that is small enough to read, plus all primary.
    gene_rows = []
    gene_sets_for_table = [
        "Cldn4",
        "Tacstd2",
        "Cldn4_Tacstd2",
        "TJ_TISMO",
        "NHEJ_KEGG",
        "Rad50",
        "MRN",
        "STING_kinases",
        "STING_core",
        "IFN_compact",
    ]
    for name in gene_sets_for_table:
        genes = resolved[name]
        deltas = []
        for gene in genes:
            vals = log_expr[gene]
            kp = [vals[i] for i in kp_idx]
            kl = [vals[i] for i in kl_idx]
            deltas.append(mean(kl) - mean(kp))
        set_delta = contrast_by_name[name]["delta"]
        for gene, gdelta in zip(genes, deltas):
            vals = log_expr[gene]
            kp = [vals[i] for i in kp_idx]
            kl = [vals[i] for i in kl_idx]
            share = ""
            if abs(set_delta) > 1e-12 and name not in ("Cldn4", "Tacstd2"):
                share = gdelta / (len(genes) * set_delta)
            gene_rows.append(
                {
                    "set_name": name,
                    "gene": gene,
                    "mean_kp": mean(kp),
                    "mean_kl": mean(kl),
                    "delta_kl_minus_kp": gdelta,
                    "separation": separation(kp, kl),
                    "share_of_set_delta": share,
                    "fpkm_all_zero": all(v == 0.0 for v in fpkm[gene]),
                }
            )

    # Leave-one-out for the four primary multi-gene scores.
    loo_rows = []
    for name in ["NHEJ_KEGG", "STING_kinases", "STING_core", "MRN", "IFN_compact", "Cldn4_Tacstd2"]:
        genes = resolved[name]
        full = contrast_by_name[name]["delta"]
        for drop in genes:
            kept = [g for g in genes if g != drop]
            sc = score_vector(log_expr, kept)
            delta = mean([sc[i] for i in kl_idx]) - mean([sc[i] for i in kp_idx])
            loo_rows.append(
                {
                    "set_name": name,
                    "dropped_gene": drop,
                    "n_genes_kept": len(kept),
                    "delta_full": full,
                    "delta_without_gene": delta,
                    "delta_change": delta - full,
                }
            )

    sample_rows = []
    for i, (library, gsm, geno, line_name) in enumerate(LIBRARIES):
        epcam_value = log_expr["Epcam"][i]
        row = {
            "library": library,
            "gsm": gsm,
            "genotype": geno,
            "line_name": line_name,
            "in_contrast": "yes" if geno in ("KP", "KL") else "no",
            "epcam_log2": f"{epcam_value:.4f}",
            "epcam_gate": "keep" if geno in ("KP", "KL") and epcam_value >= EPCAM_MIN else ("drop" if geno in ("KP", "KL") else "held_out"),
        }
        for name in scores:
            row[name] = f"{scores[name][i]:.4f}"
        sample_rows.append(row)

    genes_used_rows = []
    for name, genes in resolved.items():
        for gene in genes:
            genes_used_rows.append(
                {
                    "set_name": name,
                    "gene": gene,
                    "role": set_defs[name][0],
                    "source": set_defs[name][2],
                }
            )
    unmatched_rows = []
    for name, missing in missing_map.items():
        for gene in missing:
            unmatched_rows.append({"set_name": name, "symbol": gene})

    TABLES.mkdir(parents=True, exist_ok=True)
    score_fields = ["library", "gsm", "genotype", "line_name", "in_contrast", "epcam_log2", "epcam_gate"] + list(scores)
    write_tsv(TABLES / "sample_scores.tsv", sample_rows, score_fields)

    contrast_fields = [
        "set_name",
        "role",
        "n_genes_used",
        "n_genes_requested",
        "n_genes_unmatched",
        "mean_kp",
        "mean_kl",
        "delta",
        "separation",
        "U_kl",
        "p_mw_two_sided",
        "p_welch_greater",
        "t_welch",
        "p_perm_greater",
        "perm_extreme",
        "n_permutations",
        "rank_biserial_kl",
        "q_bh_welch_primary4",
        "note",
    ]
    # stringify with stable precision
    contrast_out = []
    float_keys = (
        "mean_kp",
        "mean_kl",
        "delta",
        "U_kl",
        "p_mw_two_sided",
        "p_welch_greater",
        "t_welch",
        "p_perm_greater",
        "rank_biserial_kl",
    )
    for row in contrast_rows:
        out = dict(row)
        for key in float_keys:
            out[key] = f"{float(row[key]):.6e}" if key.startswith("p_") or key == "t_welch" else f"{float(row[key]):.6f}"
        if row["q_bh_welch_primary4"] != "":
            out["q_bh_welch_primary4"] = f"{row['q_bh_welch_primary4']:.6e}"
        contrast_out.append(out)
    write_tsv(TABLES / "set_contrasts.tsv", contrast_out, contrast_fields)
    gate_out = []
    for row in gate_rows:
        out = dict(row)
        for key in ("mean_kp", "mean_kl", "delta", "t_welch", "p_welch_greater", "p_perm_greater", "epcam_min_log2"):
            if key.startswith("p_") or key == "t_welch":
                out[key] = f"{row[key]:.6e}"
            else:
                out[key] = f"{row[key]:.6f}"
        gate_out.append(out)
    write_tsv(
        TABLES / "epcam_gate_contrasts.tsv",
        gate_out,
        [
            "set_name",
            "n_kp",
            "n_kl",
            "dropped_libraries",
            "epcam_min_log2",
            "mean_kp",
            "mean_kl",
            "delta",
            "separation",
            "t_welch",
            "p_welch_greater",
            "p_perm_greater",
            "perm_extreme",
            "n_permutations",
        ],
    )

    gene_out = []
    for row in gene_rows:
        out = dict(row)
        for key in ("mean_kp", "mean_kl", "delta_kl_minus_kp"):
            out[key] = f"{row[key]:.6f}"
        if row["share_of_set_delta"] != "":
            out["share_of_set_delta"] = f"{row['share_of_set_delta']:.6f}"
        gene_out.append(out)
    write_tsv(
        TABLES / "gene_deltas.tsv",
        gene_out,
        [
            "set_name",
            "gene",
            "mean_kp",
            "mean_kl",
            "delta_kl_minus_kp",
            "separation",
            "share_of_set_delta",
            "fpkm_all_zero",
        ],
    )
    loo_out = []
    for row in loo_rows:
        out = dict(row)
        for key in ("delta_full", "delta_without_gene", "delta_change"):
            out[key] = f"{row[key]:.6f}"
        loo_out.append(out)
    write_tsv(
        TABLES / "leave_one_out.tsv",
        loo_out,
        ["set_name", "dropped_gene", "n_genes_kept", "delta_full", "delta_without_gene", "delta_change"],
    )
    write_tsv(
        TABLES / "genes_used.tsv",
        genes_used_rows,
        ["set_name", "gene", "role", "source"],
    )
    write_tsv(TABLES / "genes_unmatched.tsv", unmatched_rows, ["set_name", "symbol"])
    write_tsv(
        TABLES / "honest_n.tsv",
        [
            {
                "unit": "library",
                "kp": 5,
                "kl": 5,
                "normal_held_out": 1,
                "note": "Mann-Whitney uses these 10 tumor libraries. This is the unit that yields p=2/252 when separation is complete.",
            },
            {
                "unit": "line_name_on_GEO_title",
                "kp": 1,
                "kl": 5,
                "normal_held_out": 1,
                "note": "KP libraries are B6AL10-1..5 (one line prefix). KL libraries are KL155, KL47, KLC, KLD, KLE. Not mouse n.",
            },
        ],
        ["unit", "kp", "kl", "normal_held_out", "note"],
    )

    summary = {
        "accession": "GSE137244",
        "pmid": "34142094",
        "matrix": "GSE137244_counts.fpkm.csv.gz",
        "transform": "log2(FPKM+1)",
        "score": "unweighted mean across matched genes",
        "contrast": "mean KL libraries minus mean KP libraries",
        "normal_lung": "GSM4073826 held out",
        "locked_check": {
            "Cldn4_delta": contrast_by_name["Cldn4"]["delta"],
            "Tacstd2_delta": contrast_by_name["Tacstd2"]["delta"],
            "both_complete_separation": True,
            "p_mw_two_sided": contrast_by_name["Cldn4"]["p_mw_two_sided"],
            "p_welch_greater_Cldn4": contrast_by_name["Cldn4"]["p_welch_greater"],
            "p_welch_greater_STING_kinases": contrast_by_name["STING_kinases"]["p_welch_greater"],
            "p_welch_greater_joint": contrast_by_name["thesis_joint"]["p_welch_greater"],
        },
        "contrasts": contrast_rows,
        "unmatched": missing_map,
        "primary_bh": qmap,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    plot_scores([r for r in sample_rows if r["genotype"] in ("KP", "KL")], contrast_by_name)
    plot_gene_deltas(gene_rows)

    print("accession GSE137244  log2(FPKM+1)  one-sided Welch KL>KP  libraries 5 vs 5")
    print(f"{'set':22} {'delta':>8} {'sep':8} {'welch':>10} {'perm':>8}")
    for row in contrast_rows:
        if row["role"] not in ("primary", "primary_gene"):
            continue
        n_used = row["n_genes_used"]
        n_txt = f"{n_used:4d}" if n_used != "" else "    "
        print(
            f"{row['set_name']:22} {n_txt} {row['delta']:+8.3f} "
            f"{row['separation']:8} {row['p_welch_greater']:.3e} {row['p_perm_greater']:.5f}"
        )
    print("epcam gate", ",".join(dropped), f"n={len(gate_kp)} vs {len(gate_kl)}")
    for row in gate_rows:
        print(
            f"  gate {row['set_name']:18} {row['delta']:+8.3f} {row['separation']:8} "
            f"welch {row['p_welch_greater']:.3e}"
        )
    print("wrote", TABLES, "and", FIGS)


if __name__ == "__main__":
    main()
