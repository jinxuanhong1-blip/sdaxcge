#!/usr/bin/env python3
"""NicheNet ligand-activity delta, CLDN4-high vs CLDN4-low, concordant-4.

Primary activity is the official classification score: regulatory potential
as the ranking, T/NK genes associated with malignant CLDN4 as the label.
activity_delta = AUROC(top-N CLDN4-high genes) - AUROC(top-N CLDN4-low genes).
AUPR-corrected is stored beside it.

Companion: background-centered mean z of each ligand's top-K prior targets.
The T/NK axis is shifted (mean z < 0), so the centered shift is the contrast
against the rest of the measured transcriptome, not against zero.

Nothing here is fit by searching labels. The grid of N and K is written in
full. The primary cell is N=200, K=50, same-sign genes, potential-ligand null.
"""

from __future__ import annotations

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PRIOR = os.environ.get("PRIOR_SCORES", "/tmp/nichenet_work/prior_scores")
INPUTS = os.environ.get("NICHENET_INPUTS", os.path.join(HERE, "inputs"))
TAB = os.path.join(HERE, "results", "tables")
FIG = os.path.join(HERE, "figures")
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

CORE_BARRIER = ["F11R", "NECTIN2", "CDH1", "LGALS9"]
CORE_IFN = ["CXCL9", "CXCL10", "CCL5", "IFNG"]
EXT_BARRIER = [
    "PVR", "NECTIN1", "NECTIN3", "CD274", "PDCD1LG2", "CEACAM1",
    "TGFB1", "CD47", "HLA-E", "MIF", "CD24", "LGALS3",
]
EXT_IFN = ["CXCL11", "CXCL16", "CCL2", "CCL4", "CXCL12", "ICAM1", "TNFSF9", "IL15", "IL18"]
PRIMARY_N = 200
PRIMARY_K = 50
N_GRID = (50, 100, 200, 400)
K_GRID = (25, 50, 100, 200)
N_PERM = 4000
RNG = np.random.default_rng(5533)


def family_of(ligand: str) -> str:
    if ligand in CORE_BARRIER:
        return "core_barrier"
    if ligand in CORE_IFN:
        return "core_ifn_recruit"
    if ligand in EXT_BARRIER:
        return "ext_barrier"
    if ligand in EXT_IFN:
        return "ext_ifn_recruit"
    return "other"


def load_matrix():
    genes = open(os.path.join(PRIOR, "genes.txt")).read().splitlines()
    ligands = open(os.path.join(PRIOR, "ligands.txt")).read().splitlines()
    rp = np.fromfile(os.path.join(PRIOR, "rp_universe.bin"), dtype=np.float64)
    rp = rp.reshape((len(genes), len(ligands)), order="F")
    z = pd.read_csv(os.path.join(PRIOR, "universe_z.tsv"), sep="\t")
    z = z.set_index("gene").loc[genes]
    return genes, ligands, rp, z


def auroc_ap(scores: np.ndarray, y: np.ndarray) -> tuple[float, float, float]:
    """scores: n_gene. y: bool. Returns auroc, aupr, aupr_corrected."""
    if y.sum() < 5 or (~y).sum() < 5:
        return np.nan, np.nan, np.nan
    auroc = float(roc_auc_score(y, scores))
    aupr = float(average_precision_score(y, scores))
    return auroc, aupr, aupr - float(y.mean())


def batch_auroc_ap(rp: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_lig = rp.shape[1]
    auroc = np.empty(n_lig)
    aupr = np.empty(n_lig)
    corr = np.empty(n_lig)
    for j in range(n_lig):
        auroc[j], aupr[j], corr[j] = auroc_ap(rp[:, j], y)
    return auroc, aupr, corr


def target_set_auroc(rp: np.ndarray, z: np.ndarray, k: int) -> np.ndarray:
    """AUROC of meta-z, positive class = top-k regulatory-potential targets.

    >0.5 means the ligand's top targets sit toward CLDN4-high T/NK genes.
    """
    n_gene, n_lig = rp.shape
    # rank z ascending so rank 1 is the most CLDN4-low gene
    order = np.argsort(z, kind="mergesort")
    rank = np.empty(n_gene, dtype=np.float64)
    rank[order] = np.arange(1, n_gene + 1, dtype=np.float64)
    out = np.empty(n_lig)
    denom = k * (n_gene - k)
    for j in range(n_lig):
        top = np.argpartition(rp[:, j], -k)[-k:]
        u = rank[top].sum() - k * (k + 1) / 2.0
        out[j] = u / denom
    return out


def perm_mean(values: np.ndarray, idx: np.ndarray, side: str) -> tuple[float, float]:
    obs = float(np.mean(values[idx]))
    n = len(idx)
    pool = len(values)
    draws = np.empty(N_PERM)
    for i in range(N_PERM):
        draws[i] = values[RNG.choice(pool, n, replace=False)].mean()
    if side == "greater":
        p = (np.sum(draws >= obs) + 1) / (N_PERM + 1)
    else:
        p = (np.sum(draws <= obs) + 1) / (N_PERM + 1)
    return obs, float(p)


def main() -> None:
    genes, ligands, rp, zdf = load_matrix()
    z = zdf["z"].to_numpy(dtype=float)
    same_pos = (zdf["n_pos"].to_numpy() == zdf["k"].to_numpy()) & (zdf["k"].to_numpy() >= 2) & (z > 0)
    same_neg = (zdf["n_neg"].to_numpy() == zdf["k"].to_numpy()) & (zdf["k"].to_numpy() >= 2) & (z < 0)
    print(f"universe {len(genes)} ligands {len(ligands)} same+ {same_pos.sum()} same- {same_neg.sum()}", flush=True)

    send = pd.read_csv(os.path.join(INPUTS, "ligand_sender_delta.tsv"), sep="\t")
    send = send[send["split"] == "delta_q4q1"].drop_duplicates("ligand")
    pot = pd.read_csv(os.path.join(INPUTS, "potential_ligands.tsv"), sep="\t")
    pot_set = set(pot.loc[pot["potential"] == True, "ligand"])  # noqa: E712
    lig_index = {g: i for i, g in enumerate(ligands)}

    # --- geneset activity grid (official direction: score = RP) ---
    grid_rows = []
    activity_primary = None
    for n in N_GRID:
        pos_idx = np.where(same_pos)[0]
        neg_idx = np.where(same_neg)[0]
        pos_idx = pos_idx[np.argsort(-z[pos_idx])[:n]]
        neg_idx = neg_idx[np.argsort(z[neg_idx])[:n]]
        y_high = np.zeros(len(genes), dtype=bool)
        y_low = np.zeros(len(genes), dtype=bool)
        y_high[pos_idx] = True
        y_low[neg_idx] = True
        print(f"N={n} high {y_high.sum()} low {y_low.sum()}", flush=True)
        ah, ap_h, ac_h = batch_auroc_ap(rp, y_high)
        al, ap_l, ac_l = batch_auroc_ap(rp, y_low)
        block = pd.DataFrame({
            "ligand": ligands,
            "n": n,
            "auroc_high": ah,
            "auroc_low": al,
            "aupr_high": ap_h,
            "aupr_low": ap_l,
            "aupr_corrected_high": ac_h,
            "aupr_corrected_low": ac_l,
            "delta_auroc": ah - al,
            "delta_aupr_corrected": ac_h - ac_l,
        })
        grid_rows.append(block)
        if n == PRIMARY_N:
            activity_primary = block
    grid = pd.concat(grid_rows, ignore_index=True)
    grid.to_csv(os.path.join(TAB, "activity_geneset_grid.tsv"), sep="\t", index=False)

    # --- target-set shift and target AUROC ---
    shift = pd.read_csv(os.path.join(PRIOR, "ligand_target_shift.tsv"), sep="\t")
    bg_mean = float(shift["bg_mean_z"].iloc[0])
    bg_sd = float(shift["bg_sd_z"].iloc[0])
    for k in K_GRID:
        col = f"mean_z_top{k}"
        shift[f"centered_z_top{k}"] = shift[col] - bg_mean
        shift[f"cohen_d_top{k}"] = shift[f"centered_z_top{k}"] / bg_sd
    tau = {}
    for k in K_GRID:
        print(f"target AUROC K={k}", flush=True)
        tau[k] = target_set_auroc(rp, z, k)
        shift[f"target_auroc_k{k}"] = tau[k]
    shift.to_csv(os.path.join(TAB, "ligand_target_shift.tsv"), sep="\t", index=False)

    # --- ligand table at the primary cell ---
    prim = activity_primary.merge(shift, on="ligand", how="left")
    prim = prim.merge(
        send[["ligand", "family", "n_patients", "mean_delta", "p_wilcoxon",
              "mean_GSE123902", "mean_GSE131907", "mean_GSE205335", "mean_GSE189357"]],
        on="ligand", how="left",
    )
    prim = prim.merge(
        pot[["ligand", "potential", "best_receptor", "best_receptor_frac", "frac_patients_pct10"]],
        on="ligand", how="left",
    )
    prim["family"] = prim["ligand"].map(family_of)
    prim["potential"] = prim["ligand"].isin(pot_set)
    # Realized outgoing activity: sender log-expression delta, gated by the
    # NicheNet receptor detection rate, times the signed target shift.
    # Units are log1p(CP10k) x detection x (AUROC - 0.5).
    prim["target_auroc"] = prim[f"target_auroc_k{PRIMARY_K}"]
    prim["centered_z"] = prim[f"centered_z_top{PRIMARY_K}"]
    prim["outgoing_expr"] = prim["mean_delta"]
    prim["outgoing_lr"] = prim["mean_delta"] * prim["best_receptor_frac"]
    prim["activity_delta_auroc"] = prim["delta_auroc"]
    prim["activity_delta_aupr"] = prim["delta_aupr_corrected"]
    # Joint: expression delta scaled by how far target AUROC sits from 0.5.
    # Positive when the ligand is higher in CLDN4-high AND its targets sit
    # toward CLDN4-high. Negative when the ligand is higher in CLDN4-low
    # AND its targets sit toward CLDN4-low (product of two negatives).
    prim["realized_activity_delta"] = prim["mean_delta"] * (prim["target_auroc"] - 0.5) * 2.0
    prim.to_csv(os.path.join(TAB, "ligand_activity_delta.tsv"), sep="\t", index=False)

    # Null distribution: potential ligands only.
    pot_mask = prim["potential"].to_numpy()
    pot_df = prim.loc[pot_mask].reset_index(drop=True)

    def idx_of(names, frame):
        m = {g: i for i, g in enumerate(frame["ligand"])}
        return np.array([m[g] for g in names if g in m])

    metrics = {
        "delta_auroc": "greater",
        "delta_aupr_corrected": "greater",
        "target_auroc": "greater",
        "centered_z": "greater",
        "realized_activity_delta": "greater",
        "outgoing_expr": "greater",
        "outgoing_lr": "greater",
    }
    # IFN tests are the opposite tail (higher activity in CLDN4-low).
    ifn_side = {
        "delta_auroc": "less",
        "delta_aupr_corrected": "less",
        "target_auroc": "less",
        "centered_z": "less",
        "realized_activity_delta": "less",
        "outgoing_expr": "less",
        "outgoing_lr": "less",
    }

    fam_rows = []
    for fam_name, members, sides in (
        ("core_barrier", CORE_BARRIER, metrics),
        ("core_ifn_recruit", CORE_IFN, ifn_side),
    ):
        present = [g for g in members if g in set(pot_df["ligand"])]
        missing = [g for g in members if g not in set(pot_df["ligand"])]
        ii = idx_of(present, pot_df)
        for metric, side in sides.items():
            vals = pot_df[metric].to_numpy(dtype=float)
            # ligands without a sender delta are dropped from that metric
            ok = np.isfinite(vals)
            vals_ok = vals[ok]
            # map idx into the filtered array
            keep_pos = {old: new for new, old in enumerate(np.where(ok)[0])}
            ii_ok = np.array([keep_pos[i] for i in ii if i in keep_pos])
            if len(ii_ok) == 0:
                obs, p = np.nan, np.nan
            else:
                obs, p = perm_mean(vals_ok, ii_ok, side)
            fam_rows.append({
                "family": fam_name,
                "metric": metric,
                "expect": "high>low" if fam_name == "core_barrier" else "low>high",
                "n_ligands_in_null": int(len(ii_ok)),
                "ligands_scored": ",".join(pot_df.iloc[ii]["ligand"].tolist()) if len(ii) else "",
                "ligands_not_potential": ",".join(missing),
                "mean": obs,
                "perm_p": p,
                "side": side,
                "n_null": int(ok.sum()),
            })
    fam = pd.DataFrame(fam_rows)
    fam.to_csv(os.path.join(TAB, "family_activity_delta.tsv"), sep="\t", index=False)

    # Grid summary for the two families, all ligands in the matrix (activity
    # does not require the expression filter; sender-gated metrics do).
    grid_sum = []
    for n, g in grid.groupby("n"):
        for fam_name, members, col, expect_pos in (
            ("core_barrier", CORE_BARRIER, "delta_auroc", True),
            ("core_barrier", CORE_BARRIER, "delta_aupr_corrected", True),
            ("core_ifn_recruit", CORE_IFN, "delta_auroc", False),
            ("core_ifn_recruit", CORE_IFN, "delta_aupr_corrected", False),
        ):
            sub = g[g.ligand.isin(members)]
            grid_sum.append({
                "family": fam_name,
                "n": int(n),
                "metric": col,
                "mean": float(sub[col].mean()),
                "min": float(sub[col].min()),
                "max": float(sub[col].max()),
                "n_ligands": int(len(sub)),
            })
    for k in K_GRID:
        for fam_name, members in (("core_barrier", CORE_BARRIER), ("core_ifn_recruit", CORE_IFN)):
            sub = shift[shift.ligand.isin(members)]
            grid_sum.append({
                "family": fam_name,
                "n": int(k),
                "metric": f"centered_z_top{k}",
                "mean": float(sub[f"centered_z_top{k}"].mean()),
                "min": float(sub[f"centered_z_top{k}"].min()),
                "max": float(sub[f"centered_z_top{k}"].max()),
                "n_ligands": int(len(sub)),
            })
            grid_sum.append({
                "family": fam_name,
                "n": int(k),
                "metric": f"target_auroc_k{k}",
                "mean": float(sub[f"target_auroc_k{k}"].mean()),
                "min": float(sub[f"target_auroc_k{k}"].min()),
                "max": float(sub[f"target_auroc_k{k}"].max()),
                "n_ligands": int(len(sub)),
            })
    pd.DataFrame(grid_sum).to_csv(os.path.join(TAB, "grid_family_summary.tsv"), sep="\t", index=False)

    # Called-ligand card
    called = CORE_BARRIER + CORE_IFN + ["CXCL11", "CCL4", "PVR", "CD274", "HLA-A", "TGFB1"]
    card = prim[prim.ligand.isin(called)].copy()
    # ranks among potential ligands
    for col in ("delta_auroc", "delta_aupr_corrected", "target_auroc", "centered_z",
                "realized_activity_delta", "outgoing_expr", "outgoing_lr"):
        card[f"rank_{col}"] = card["ligand"].map(
            pot_df[col].rank(ascending=False, method="min")
            if False else
            prim.loc[prim.potential].set_index("ligand")[col].rank(ascending=False, method="min")
        )
    card.to_csv(os.path.join(TAB, "called_ligand_activity.tsv"), sep="\t", index=False)

    # Which N / K maximizes the pre-specified family contrast.
    # Contrast = barrier mean - IFN mean, on delta_auroc (positive is on-thesis)
    # and on centered z (barrier positive, IFN negative => barrier - IFN positive).
    best = {"primary_N": PRIMARY_N, "primary_K": PRIMARY_K, "bg_mean_z": bg_mean, "bg_sd_z": bg_sd,
            "n_universe": len(genes), "n_potential": int(pot_mask.sum())}
    gsum = pd.DataFrame(grid_sum)
    dlt = gsum[gsum.metric == "delta_auroc"]
    wide = dlt.pivot(index="n", columns="family", values="mean")
    wide["contrast"] = wide["core_barrier"] - wide["core_ifn_recruit"]
    best_n = int(wide["contrast"].idxmax())
    best["max_delta_auroc_contrast_N"] = best_n
    best["max_delta_auroc_contrast"] = float(wide.loc[best_n, "contrast"])
    best["max_delta_auroc_barrier"] = float(wide.loc[best_n, "core_barrier"])
    best["max_delta_auroc_ifn"] = float(wide.loc[best_n, "core_ifn_recruit"])
    best["grid_delta_auroc"] = wide.reset_index().to_dict(orient="records")

    cz = gsum[gsum.metric.str.startswith("centered_z_top")]
    # metric name encodes k; n column is k
    rows = []
    for _, r in cz.iterrows():
        rows.append(r)
    cz = pd.DataFrame(rows)
    widek = cz.pivot(index="n", columns="family", values="mean")
    widek["contrast"] = widek["core_barrier"] - widek["core_ifn_recruit"]
    best_k = int(widek["contrast"].idxmax())
    best["max_centered_z_contrast_K"] = best_k
    best["max_centered_z_contrast"] = float(widek.loc[best_k, "contrast"])
    best["max_centered_z_barrier"] = float(widek.loc[best_k, "core_barrier"])
    best["max_centered_z_ifn"] = float(widek.loc[best_k, "core_ifn_recruit"])
    best["grid_centered_z"] = widek.reset_index().to_dict(orient="records")
    # Specificity of the geneset activity delta: most ligands are positive.
    spec_rows = []
    spec_rng = np.random.default_rng(5533)
    for n, g in grid.groupby("n"):
        sub = g[g.ligand.isin(pot_set)].reset_index(drop=True)
        vals = sub["delta_auroc"].to_numpy(dtype=float)
        b_idx = np.array([i for i, lig in enumerate(sub["ligand"]) if lig in CORE_BARRIER])
        bmean = float(vals[b_idx].mean())
        draws = np.array([
            vals[spec_rng.choice(len(vals), len(b_idx), replace=False)].mean()
            for _ in range(N_PERM)
        ])
        spec_rows.append({
            "n": int(n),
            "barrier_mean_delta_auroc": bmean,
            "potential_median_delta_auroc": float(np.median(vals)),
            "potential_mean_delta_auroc": float(np.mean(vals)),
            "excess_over_median": bmean - float(np.median(vals)),
            "perm_p_greater": float((np.sum(draws >= bmean) + 1) / (N_PERM + 1)),
        })
    spec = pd.DataFrame(spec_rows)
    spec.to_csv(os.path.join(TAB, "specificity_by_n.tsv"), sep="\t", index=False)
    best["specificity"] = spec.to_dict(orient="records")

    tail = shift.nsmallest(15, "target_auroc_k50")[
        ["ligand", "target_auroc_k50", "centered_z_top50", "mean_z_top50", "mean_rp_top50"]
    ]
    tail.to_csv(os.path.join(TAB, "low_target_auroc_tail.tsv"), sep="\t", index=False)

    with open(os.path.join(TAB, "summary.json"), "w") as f:
        json.dump(best, f, indent=2)

    print(json.dumps(best, indent=2))
    print(fam.to_string(index=False))
    cols = ["ligand", "family", "potential", "mean_delta", "best_receptor", "best_receptor_frac",
            "delta_auroc", "delta_aupr_corrected", "target_auroc", "centered_z",
            "realized_activity_delta", "outgoing_lr"]
    print(card[cols].sort_values(["family", "ligand"]).to_string(index=False))
    plot(card, pot_df, wide, widek)
    print("done", flush=True)


def plot(card: pd.DataFrame, pot_df: pd.DataFrame, wide_n: pd.DataFrame, wide_k: pd.DataFrame) -> None:
    order = CORE_BARRIER + CORE_IFN
    show = card.set_index("ligand").loc[order]
    colors = ["#1f4e79"] * 4 + ["#b85c38"] * 4
    y = np.arange(len(order))

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.6))
    ax = axes[0]
    ax.barh(y, show["outgoing_lr"].to_numpy(), color=colors, height=0.72)
    ax.axvline(0, color="#444", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("LR-gated outgoing Δ\nlog1p(CP10k) × receptor detection")
    ax.set_title("Malignant sender, Q4 − Q1")
    ax = axes[1]
    ax.barh(y, (show["target_auroc"] - 0.5).to_numpy(), color=colors, height=0.72)
    ax.axvline(0, color="#444", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("Target-set AUROC − 0.5\n>0 targets sit toward CLDN4-high")
    ax.set_title("Top 50 prior targets in T/NK")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_outgoing_and_targets.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "fig1_outgoing_and_targets.pdf"))
    plt.close()

    grid = pd.read_csv(os.path.join(TAB, "activity_geneset_grid.tsv"), sep="\t")
    pot_lig = set(pot_df["ligand"])
    med_map = {}
    for n, g in grid.groupby("n"):
        med_map[int(n)] = float(g.loc[g.ligand.isin(pot_lig), "delta_auroc"].median())

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2))
    ax = axes[0]
    ax.plot(wide_n.index, wide_n["core_barrier"], "o-", color="#1f4e79", label="barrier mean")
    ax.plot(wide_n.index, wide_n["core_ifn_recruit"], "o-", color="#b85c38", label="IFN / recruit mean")
    xs = list(wide_n.index.astype(int))
    ax.plot(xs, [med_map[int(n)] for n in xs], "s--", color="#6b7280", label="median potential ligand")
    ax.axhline(0, color="#888", lw=0.6)
    ax.set_xlabel("Geneset size N")
    ax.set_ylabel("Mean activity Δ AUROC")
    ax.set_title("Geneset grid")
    ax.legend(frameon=False, fontsize=8)
    ax = axes[1]
    ax.plot(wide_k.index, wide_k["core_barrier"], "o-", color="#1f4e79", label="barrier")
    ax.plot(wide_k.index, wide_k["core_ifn_recruit"], "o-", color="#b85c38", label="IFN / recruit")
    ax.axhline(0, color="#888", lw=0.6)
    ax.set_xlabel("Top-K prior targets")
    ax.set_ylabel("Mean centered target z")
    ax.set_title("Target-shift grid")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_grid.png"), dpi=160)
    fig.savefig(os.path.join(FIG, "fig2_grid.pdf"))
    plt.close()


if __name__ == "__main__":
    main()
