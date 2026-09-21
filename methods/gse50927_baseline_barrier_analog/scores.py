"""Alternate scores for the baseline KO vs WT barrier-loss analog.

Sets are fixed before the comparison. Score definitions are a grid.
The reported "max" is the largest permutation z inside a predeclared
question. It is not a search over genes.

n remains 1 vs 1 GEO sample. Permutations reshuffle genes on the
deposited ranking. They are not a mouse-level test.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

N_PERM = 1000
SEED = 42
FDR_CUT = 0.05
P_CAP = 10.0

# Canonical IFN-gamma chemokines. Inflammatory set is the ELR+ CXC /
# immediate CC burst used as the VILI confound family. Neither list is
# rebuilt from this dataset's FDR hits.
CHEMOKINE_IFN = ["Cxcl9", "Cxcl10", "Cxcl11", "Ccl5"]
CHEMOKINE_INFLAM = ["Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Ccl2", "Ccl3", "Ccl4", "Ccl7"]
CHEMOKINE_ALL = [
    "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
    "Ccl11", "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl22", "Ccl24", "Ccl25",
    "Ccl27a", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10", "Cxcl11",
    "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl16", "Cxcl17", "Cx3cl1", "Xcl1",
]
# No H2 symbol. Heavy chains stay in a separate row because the line is
# mixed 129 / B6 / BALB. Chaperones are general ER genes, so a second
# row drops Calr / Canx / Pdia3.
MHC_SAFE = [
    "B2m", "Tap1", "Tap2", "Tapbp", "Tapbpl",
    "Psmb8", "Psmb9", "Psmb10", "Psme1", "Psme2",
    "Nlrc5", "Erap1", "Calr", "Canx", "Pdia3",
]
MHC_SAFE_NO_CHAPERONE = [
    "B2m", "Tap1", "Tap2", "Tapbp", "Tapbpl",
    "Psmb8", "Psmb9", "Psmb10", "Psme1", "Psme2",
    "Nlrc5", "Erap1",
]
MHC_HEAVY = ["H2-K1", "H2-D1"]
MHC_HAPLOTYPE = [
    "H2-K2", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q8",
    "H2-Q10", "H2-T23", "H2-M2", "H2-M3", "H2-Bl",
]
INJURY = ["Tnf", "Il1b", "Il6", "Egr1"]

# Pulled out of the Hallmark intersection so the IFN cassette does not
# contain the chemokine module or the haplotype-safe MHC module.
CORE_DROP = {
    "Cxcl9", "Cxcl10", "Cxcl11", "Ccl5", "Ccl7", "Xcl1", "Cxcl1", "Cxcl2",
    "B2m", "Tap1", "Tap2", "Tapbp", "Tapbpl", "Psmb8", "Psmb9", "Psmb10",
    "Psme1", "Psme2", "Nlrc5", "Cd74", "Ciita", "Erap1", "Calr", "Canx", "Pdia3",
    "Il6",
}


def load_hallmark(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        for row in rows:
            sets.setdefault(row["set"], []).append(row["gene"])
    return sets


def ifn_core(hallmark: dict[str, list[str]]) -> list[str]:
    alpha = set(hallmark["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
    gamma = set(hallmark["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
    core = sorted((alpha & gamma) - CORE_DROP)
    core = [g for g in core if not g.startswith("H2-")]
    return core


def vili_up(rec: dict) -> bool:
    return rec["logFC"] > 0 and rec["FDR"] < FDR_CUT and rec["logCPM"] > 0


def gene_value(rec: dict, kind: str, rank_pct: float) -> float:
    if kind == "logFC":
        return rec["logFC"]
    if kind == "up":
        return 1.0 if rec["logFC"] > 0 else 0.0
    if kind == "fdr_up":
        return 1.0 if rec["logFC"] > 0 and rec["FDR"] < FDR_CUT else 0.0
    if kind == "signed_p":
        p = min(max(rec["P"], 1e-300), 1.0)
        mag = min(-np.log10(p), P_CAP)
        if rec["logFC"] > 0:
            return mag
        if rec["logFC"] < 0:
            return -mag
        return 0.0
    if kind == "rank":
        return rank_pct
    raise KeyError(kind)


def eligible(table: dict, mode: str, drop_vili: dict | None) -> list[str]:
    out = []
    for symbol, rec in table.items():
        if mode == "all":
            keep = True
        elif mode == "detected":
            keep = rec["logCPM"] > 0
        elif mode == "cpm_gt_1":
            keep = rec["logCPM"] > 1
        else:
            raise KeyError(mode)
        if keep and drop_vili is not None and symbol in drop_vili and vili_up(drop_vili[symbol]):
            keep = False
        if keep:
            out.append(symbol)
    return out


DEFINITIONS = [
    ("mean_logFC_detected", "detected", "logFC", "mean", False),
    ("median_logFC_detected", "detected", "logFC", "median", False),
    ("mean_logFC_all", "all", "logFC", "mean", False),
    ("mean_logFC_cpm_gt_1", "cpm_gt_1", "logFC", "mean", False),
    ("frac_up_detected", "detected", "up", "mean", False),
    ("frac_fdr_up_detected", "detected", "fdr_up", "mean", False),
    ("mean_signed_neglog10P_detected", "detected", "signed_p", "mean", False),
    ("mean_rank_pct_detected", "detected", "rank", "mean", False),
    ("mean_logFC_detected_not_VILI_up", "detected", "logFC", "mean", True),
]


def aggregate(vals: np.ndarray, how: str) -> float:
    if how == "mean":
        return float(vals.mean())
    if how == "median":
        return float(np.median(vals))
    raise KeyError(how)


def rank_pct_map(table: dict) -> dict[str, float]:
    symbols = list(table)
    logfc = np.array([table[s]["logFC"] for s in symbols], dtype=float)
    order = np.argsort(logfc, kind="mergesort")
    sorted_fc = logfc[order]
    ranks = np.empty(len(symbols), dtype=float)
    n = len(symbols)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and sorted_fc[j + 1] == sorted_fc[i]:
            j += 1
        ranks[order[i : j + 1]] = 0.5 * ((i + 1) + (j + 1))
        i = j + 1
    return {symbols[k]: 100.0 * ranks[k] / n for k in range(n)}


def score_members(table: dict, members: list[str], kind: str, how: str, ranks: dict[str, float]) -> float | None:
    if not members:
        return None
    vals = np.array([gene_value(table[g], kind, ranks[g]) for g in members], dtype=float)
    return aggregate(vals, how)


def perm_gap(
    values: dict[str, float],
    set_a: list[str],
    set_b: list[str],
    how: str,
    rng: np.random.Generator,
) -> dict | None:
    a = [g for g in set_a if g in values]
    b = [g for g in set_b if g in values]
    if len(a) < 3 or len(b) < 3:
        return {
            "n_a": len(a),
            "n_b": len(b),
            "score_a": score_of(values, a, how),
            "score_b": score_of(values, b, how),
            "gap": "",
            "null_mean": "",
            "null_sd": "",
            "z": "",
            "perm_p_greater": "",
            "status": "n_below_3",
        }
    obs_a = score_of(values, a, how)
    obs_b = score_of(values, b, how)
    obs = obs_a - obs_b
    universe = np.array(list(values))
    pool = np.array([values[g] for g in universe])
    null = np.empty(N_PERM, dtype=float)
    na, nb = len(a), len(b)
    for i in range(N_PERM):
        pick = rng.choice(len(pool), size=na + nb, replace=False)
        null[i] = aggregate(pool[pick[:na]], how) - aggregate(pool[pick[na:]], how)
    sd = float(null.std(ddof=1))
    mu = float(null.mean())
    z = (obs - mu) / sd if sd > 0 else float("nan")
    p = (1.0 + np.sum(null >= obs)) / (N_PERM + 1.0)
    return {
        "n_a": na,
        "n_b": nb,
        "score_a": obs_a,
        "score_b": obs_b,
        "gap": obs,
        "null_mean": mu,
        "null_sd": sd,
        "z": z,
        "perm_p_greater": float(p),
        "status": "ok",
    }


def score_of(values: dict[str, float], members: list[str], how: str):
    if not members:
        return ""
    return aggregate(np.array([values[g] for g in members], dtype=float), how)


def value_map(table: dict, universe: list[str], kind: str, ranks: dict[str, float]) -> dict[str, float]:
    return {g: gene_value(table[g], kind, ranks[g]) for g in universe}


def fnum(x) -> str:
    if x == "" or x is None:
        return ""
    return f"{float(x):.9g}"


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def plot_scores(group_rows: list[dict], gap_rows: list[dict], figures: Path) -> None:
    figures.mkdir(parents=True, exist_ok=True)
    want = [
        "ifn_core",
        "chemokine_ifn",
        "chemokine_inflammatory",
        "mhc_haplotype_safe",
        "injury_mediator",
    ]
    labels = {
        "ifn_core": "IFN core",
        "chemokine_ifn": "IFN chemokines",
        "chemokine_inflammatory": "Inflammatory chemokines",
        "mhc_haplotype_safe": "MHC-I haplotype-safe",
        "injury_mediator": "Tnf/Il1b/Il6/Egr1",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.8))

    ax = axes[0]
    y = np.arange(len(want))
    base_vals, vili_vals, ylabels = [], [], []
    for name in want:
        b = next(r for r in group_rows if r["set"] == name and r["contrast"] == "naive_KO_vs_WT" and r["definition"] == "mean_logFC_detected")
        v = next(r for r in group_rows if r["set"] == name and r["contrast"] == "VILI_WT_vs_naive_WT" and r["definition"] == "mean_logFC_detected")
        base_vals.append(float(b["score"]))
        vili_vals.append(float(v["score"]))
        ylabels.append(f"{labels[name]}\n(n={b['n_scored']} / {v['n_scored']})")
    ax.hlines(y, base_vals, vili_vals, color="#d5d8dc", lw=2, zorder=1)
    ax.scatter(base_vals, y, s=42, color="#1f4e79", zorder=2, label="Baseline KO − WT")
    ax.scatter(vili_vals, y, s=42, color="#a93226", zorder=2, label="WT VILI − naive")
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_xlabel("Mean logFC, genes with logCPM > 0")
    ax.set_title("Same sets, two contrasts")
    ax.legend(frameon=False, fontsize=8, loc="best")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    q1 = [r for r in gap_rows if r["question"] == "Q1_ifn_core_minus_inflam" and r["status"] == "ok"]
    q1 = sorted(q1, key=lambda r: float(r["z"]))
    zs = [float(r["z"]) for r in q1]
    ax.barh(
        np.arange(len(q1)),
        zs,
        color=["#1f4e79" if z > 0 else "#a93226" for z in zs],
    )
    span = max(abs(z) for z in zs) if zs else 1
    for i, z in enumerate(zs):
        ax.text(
            z + (0.04 * span if z >= 0 else -0.04 * span),
            i,
            f"{z:.2f}",
            va="center",
            ha="left" if z >= 0 else "right",
            fontsize=7,
            color="#1c2833",
        )
    ax.set_yticks(np.arange(len(q1)))
    ax.set_yticklabels([r["definition"] for r in q1], fontsize=7)
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.set_xlabel("Permutation z on the baseline gap (genes, not mice)")
    ax.set_title("Core − inflammatory chemokines\n(negative: inflammatory logFC is larger)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle(
        "GSE50927 barrier-loss analog, n = 1 vs 1 (not cancer)",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(figures / "fig_score_definitions.png", dpi=160, bbox_inches="tight")
    fig.savefig(figures / "fig_score_definitions.pdf", bbox_inches="tight")
    plt.close(fig)


def run_scores(tables: dict, root: Path) -> dict:
    hallmark = load_hallmark(root / "genesets" / "hallmark_ifn_mm.tsv")
    core = ifn_core(hallmark)
    other_chem = [g for g in CHEMOKINE_ALL if g not in set(CHEMOKINE_IFN) | set(CHEMOKINE_INFLAM)]
    sets = {
        "ifn_core": core,
        "chemokine_ifn": CHEMOKINE_IFN,
        "chemokine_inflammatory": CHEMOKINE_INFLAM,
        "chemokine_other": other_chem,
        "mhc_haplotype_safe": MHC_SAFE,
        "mhc_haplotype_safe_no_chaperone": MHC_SAFE_NO_CHAPERONE,
        "mhc_heavy_chain_H2K1_H2D1": MHC_HEAVY,
        "mhc_haplotype_risk": MHC_HAPLOTYPE,
        "injury_mediator": INJURY,
    }
    baseline = tables["naive_KO_vs_WT"]
    vili = tables["VILI_WT_vs_naive_WT"]
    ranks_b = rank_pct_map(baseline)
    ranks_v = rank_pct_map(vili)

    membership = []
    for name, genes in sets.items():
        for gene in genes:
            rec_b = baseline.get(gene)
            rec_v = vili.get(gene)
            membership.append({
                "set": name,
                "gene": gene,
                "in_baseline_table": "yes" if rec_b else "no",
                "logFC_baseline": fnum(rec_b["logFC"]) if rec_b else "",
                "logCPM_baseline": fnum(rec_b["logCPM"]) if rec_b else "",
                "FDR_baseline": fnum(rec_b["FDR"]) if rec_b else "",
                "logFC_VILI_WT": fnum(rec_v["logFC"]) if rec_v else "",
                "logCPM_VILI_WT": fnum(rec_v["logCPM"]) if rec_v else "",
                "FDR_VILI_WT": fnum(rec_v["FDR"]) if rec_v else "",
                "VILI_FDR_up": "yes" if rec_v and vili_up(rec_v) else "no",
            })

    group_rows = []
    for contrast_name, table, ranks in (
        ("naive_KO_vs_WT", baseline, ranks_b),
        ("VILI_WT_vs_naive_WT", vili, ranks_v),
    ):
        for def_name, mode, kind, how, drop in DEFINITIONS:
            # Injury exclusion uses the WT VILI table even when the score
            # itself is the baseline logFC.
            universe = eligible(table, mode, vili if drop else None)
            universe_set = set(universe)
            values = value_map(table, universe, kind, ranks)
            for set_name, genes in sets.items():
                members = [g for g in genes if g in universe_set]
                sc = score_members(table, members, kind, how, ranks) if members else None
                group_rows.append({
                    "contrast": contrast_name,
                    "definition": def_name,
                    "set": set_name,
                    "n_set_nominal": len(genes),
                    "n_scored": len(members),
                    "score": fnum(sc) if sc is not None else "",
                })

    questions = [
        ("Q1_ifn_core_minus_inflam", "ifn_core", "chemokine_inflammatory"),
        ("Q2_ifn_chemokine_minus_inflam", "chemokine_ifn", "chemokine_inflammatory"),
        ("Q3_mhc_safe_minus_heavy_chain", "mhc_haplotype_safe", "mhc_heavy_chain_H2K1_H2D1"),
        ("Q4_mhc_safe_no_chaperone_minus_heavy", "mhc_haplotype_safe_no_chaperone", "mhc_heavy_chain_H2K1_H2D1"),
    ]
    # Q3/Q4 use heavy chain as the comparator so a "win" cannot come from
    # beating H2-K2. Heavy chain is only two genes, so n_b < 3 and those
    # questions will be status n_below_3 on the gap test. Absolute
    # enrichment of the haplotype-safe set is computed separately below.
    rng = np.random.default_rng(SEED)
    gap_rows = []
    for def_name, mode, kind, how, drop in DEFINITIONS:
        universe = eligible(baseline, mode, vili if drop else None)
        values = value_map(baseline, universe, kind, ranks_b)
        for q_name, a_name, b_name in questions:
            result = perm_gap(values, sets[a_name], sets[b_name], how, rng)
            gap_rows.append({
                "question": q_name,
                "definition": def_name,
                "set_a": a_name,
                "set_b": b_name,
                "contrast": "naive_KO_vs_WT",
                "n_perm": N_PERM,
                "seed": SEED,
                **{k: fnum(v) if k not in {"n_a", "n_b", "status"} else v for k, v in result.items()},
            })

    # Absolute enrichment: set mean vs random genes of the same size.
    enrich_rows = []
    for set_name in ("ifn_core", "mhc_haplotype_safe", "mhc_haplotype_safe_no_chaperone", "chemokine_ifn", "chemokine_inflammatory"):
        for def_name, mode, kind, how, drop in DEFINITIONS:
            if drop:
                continue
            universe = eligible(baseline, mode, None)
            values = value_map(baseline, universe, kind, ranks_b)
            members = [g for g in sets[set_name] if g in values]
            if len(members) < 3:
                enrich_rows.append({
                    "set": set_name,
                    "definition": def_name,
                    "n": len(members),
                    "score": fnum(score_of(values, members, how)) if members else "",
                    "z": "",
                    "perm_p_greater": "",
                    "status": "n_below_3",
                })
                continue
            obs = score_of(values, members, how)
            pool = np.array(list(values.values()))
            null = np.empty(N_PERM)
            n = len(members)
            for i in range(N_PERM):
                pick = rng.choice(len(pool), size=n, replace=False)
                null[i] = aggregate(pool[pick], how)
            sd = float(null.std(ddof=1))
            z = (obs - float(null.mean())) / sd if sd > 0 else float("nan")
            p = (1.0 + np.sum(null >= obs)) / (N_PERM + 1.0)
            enrich_rows.append({
                "set": set_name,
                "definition": def_name,
                "n": n,
                "score": fnum(obs),
                "null_mean": fnum(null.mean()),
                "z": fnum(z),
                "perm_p_greater": fnum(p),
                "status": "ok",
            })

    # Confound-adjusted gap: (baseline − VILI) on mean logFC detected.
    # Predeclared companion. Not substituted for a baseline-only win.
    adj_rows = []
    for def_name, mode, kind, how, drop in DEFINITIONS:
        if drop:
            continue
        # Values are baseline statistic minus the same statistic on the VILI table.
        # Built per gene so a permutation still moves genes, not mice.
        uni_b = set(eligible(baseline, mode, None))
        uni_v = set(eligible(vili, mode, None))
        shared = sorted(uni_b & uni_v)
        values = {}
        for g in shared:
            values[g] = gene_value(baseline[g], kind, ranks_b[g]) - gene_value(vili[g], kind, ranks_v[g])
        for q_name, a_name, b_name in (
            ("Q1_adj_baseline_minus_VILI", "ifn_core", "chemokine_inflammatory"),
            ("Q2_adj_baseline_minus_VILI", "chemokine_ifn", "chemokine_inflammatory"),
        ):
            result = perm_gap(values, sets[a_name], sets[b_name], how, rng)
            adj_rows.append({
                "question": q_name,
                "definition": def_name,
                "contrast_math": "baseline_statistic minus WT_VILI_statistic",
                "n_perm": N_PERM,
                "seed": SEED,
                **{k: fnum(v) if k not in {"n_a", "n_b", "status"} else v for k, v in result.items()},
            })

    tables_dir = root / "tables"
    write_tsv(tables_dir / "score_membership.tsv", membership)
    write_tsv(tables_dir / "score_by_set.tsv", group_rows)
    write_tsv(tables_dir / "score_gaps.tsv", gap_rows)
    write_tsv(tables_dir / "score_enrichment.tsv", enrich_rows)
    write_tsv(tables_dir / "score_confound_adjusted.tsv", adj_rows)
    plot_scores(group_rows, gap_rows, root / "figures")

    def best(rows: list[dict], question: str) -> dict | None:
        ok = [r for r in rows if r.get("question") == question and r["status"] == "ok" and r["z"] != ""]
        if not ok:
            return None
        return max(ok, key=lambda r: float(r["z"]))

    winners = {
        "Q1_baseline_ifn_core_minus_inflam": best(gap_rows, "Q1_ifn_core_minus_inflam"),
        "Q2_baseline_ifn_chemokine_minus_inflam": best(gap_rows, "Q2_ifn_chemokine_minus_inflam"),
        "Q1_confound_adjusted": best(adj_rows, "Q1_adj_baseline_minus_VILI"),
        "Q2_confound_adjusted": best(adj_rows, "Q2_adj_baseline_minus_VILI"),
    }
    # Enrichment winners are per set, max z, reported separately.
    enrich_winners = {}
    for set_name in ("ifn_core", "mhc_haplotype_safe", "mhc_haplotype_safe_no_chaperone", "chemokine_ifn"):
        ok = [r for r in enrich_rows if r["set"] == set_name and r["status"] == "ok"]
        enrich_winners[set_name] = max(ok, key=lambda r: float(r["z"])) if ok else None

    print("IFN core", len(core), "missing", [g for g in core if g not in baseline])
    for key, row in winners.items():
        print("WIN", key, row["definition"] if row else None, row["z"] if row else None, "gap", row["gap"] if row else None)
    for key, row in enrich_winners.items():
        print("ENR", key, row["definition"] if row else None, "z", row["z"] if row else None, "score", row["score"] if row else None)
    return {"winners": winners, "enrich_winners": enrich_winners, "ifn_core_n": len(core)}
