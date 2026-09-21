#!/usr/bin/env python3
"""Maximize the concordant-4 patient-level CLDN4 vs T/NK LR and |rho|.

Same 65 locked units as the published malignant CLDN4 % positive result.
Continuous input is the scVI-normalized CLDN4 already computed on the
patient-batch integration subsample. Nothing is typed into the result tables.

Search (all reported in results/tables/glmm_grid.tsv):
  scores     % positive and mean log1p (baselines);
             seed-malignant scVI mean;
             seed- and scANVI-malignant fraction with scVI log1p(CLDN4)
             >= log1p(k) for k = 1..20 CP10K;
             the same fraction at the y-blind histogram antimode
  outcomes   T/NK vs every other cell; T/NK vs malignant
  coding     global z-score; within-dataset z-score; within-dataset rank

The patient model is
  cbind(n_T/NK, n_fail) ~ coded(CLDN4) + dataset + (1 | unit).
The LR effect is that model's likelihood-ratio chi-square. |rho| is the
DerSimonian-Laird pool of within-cohort Spearmans (invariant to coding).

The primary row maximizes the worse of (chi-square / locked chi-square) and
(|rho| / locked |rho|) among scVI/scANVI scores. Naive p-values ignore that
search. A within-cohort permutation re-picks the row each shuffle.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
INP = HERE / "inputs"
TAB = HERE / "results" / "tables"
FIG = HERE / "results" / "figures"
K_GRID = list(range(1, 21))
N_PERM = 2000
RNG_SEED = 1
DATASETS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]


def ranks(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a, dtype=float)
    order = np.argsort(a, kind="mergesort")
    out = np.empty(len(a), dtype=float)
    i = 0
    while i < len(a):
        j = i
        while j + 1 < len(a) and a[order[j + 1]] == a[order[i]]:
            j += 1
        out[order[i : j + 1]] = (i + j) / 2.0 + 1.0
        i = j + 1
    return out


def spearman(x, y) -> tuple[float, int]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 4:
        return float("nan"), n
    rx, ry = ranks(x[m]), ranks(y[m])
    rx = rx - rx.mean()
    ry = ry - ry.mean()
    den = math.sqrt(float(np.dot(rx, rx) * np.dot(ry, ry)))
    if den == 0:
        return float("nan"), n
    return float(np.dot(rx, ry) / den), n


def dl_spearman(pairs: list[tuple[float, int]]) -> dict:
    zs, vs, ns = [], [], []
    for rho, n in pairs:
        if n <= 3 or not math.isfinite(rho):
            continue
        zs.append(math.atanh(min(max(rho, -0.999999), 0.999999)))
        vs.append(1.0 / (n - 3.0))
        ns.append(n)
    z = np.asarray(zs, dtype=float)
    v = np.asarray(vs, dtype=float)
    w = 1.0 / v
    zbar = float(np.sum(w * z) / np.sum(w))
    q = float(np.sum(w * (z - zbar) ** 2))
    k = len(z)
    df = k - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - df) / c) if df > 0 and c > 0 else 0.0
    wstar = 1.0 / (v + tau2)
    zre = float(np.sum(wstar * z) / np.sum(wstar))
    se = math.sqrt(1.0 / float(np.sum(wstar)))
    p = float(math.erfc(abs(zre / se) / math.sqrt(2.0)))
    i2 = max(0.0, (q - df) / q) if q > 0 else 0.0
    return {
        "rho": float(math.tanh(zre)),
        "p": p,
        "I2": float(i2),
        "ci_lo": float(math.tanh(zre - 1.96 * se)),
        "ci_hi": float(math.tanh(zre + 1.96 * se)),
        "k": k,
        "N": int(sum(ns)),
        "se": se,
    }


def dl_many(rho: np.ndarray, ns: np.ndarray) -> np.ndarray:
    """Pooled Spearman for rho shaped (m, k)."""
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    v = 1.0 / (ns.astype(float) - 3.0)
    w = 1.0 / v
    zbar = (w * z).sum(axis=1) / w.sum()
    q = (w * (z - zbar[:, None]) ** 2).sum(axis=1)
    df = rho.shape[1] - 1
    c = float(w.sum() - np.sum(w ** 2) / w.sum())
    tau2 = np.maximum(0.0, (q - df) / c) if df > 0 and c > 0 else np.zeros(rho.shape[0])
    wstar = 1.0 / (v + tau2[:, None])
    zre = (wstar * z).sum(axis=1) / wstar.sum(axis=1)
    return np.tanh(zre)


def load_units() -> list[dict]:
    with (INP / "patient_units.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if len(rows) != 65:
        raise SystemExit(f"expected 65 units, found {len(rows)}")
    return rows


def load_cells():
    cells = defaultdict(list)
    with gzip.open(INP / "umap_obs.tsv.gz", "rt") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            cells[(row["dataset"], row["unit_id"])].append(
                (row["seed_class"], row["scanvi_pred"], float(row["CLDN4_scvi_log1p"]))
            )
    return cells


def antimode(values: np.ndarray) -> float:
    """Lightest histogram bin between the 5th and 60th percentiles."""
    lo = float(np.quantile(values, 0.05))
    hi = float(np.quantile(values, 0.60))
    _hist, edges = np.histogram(values, bins=np.linspace(lo, hi, 40))
    hist = _hist
    if len(hist) < 5:
        raise SystemExit("antimode histogram failed")
    j = int(np.argmin(hist[1:-1])) + 1
    return float(0.5 * (edges[j] + edges[j + 1]))


def fraction(xs: list[float], thr: float) -> float:
    if len(xs) < 10:
        return float("nan")
    return float(np.mean(np.asarray(xs) >= thr))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise SystemExit(f"no rows for {path}")
    with path.open("w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def code_x(raw: np.ndarray, dataset: np.ndarray, coding: str) -> np.ndarray:
    x = np.zeros(len(raw), dtype=float)
    if coding == "global_z":
        s = float(raw.std(ddof=1))
        return (raw - raw.mean()) / s if s else x
    present = [ds for ds in DATASETS if np.any(dataset == ds)]
    for ds in present:
        idx = np.flatnonzero(dataset == ds)
        v = raw[idx].copy()
        if coding == "within_rank":
            v = ranks(v)
        s = float(v.std(ddof=1))
        x[idx] = (v - v.mean()) / s if s else 0.0
    return x


def design(dataset: np.ndarray, x: np.ndarray | None) -> np.ndarray:
    levels = [ds for ds in DATASETS if np.any(dataset == ds)]
    dummies = np.zeros((len(dataset), max(len(levels) - 1, 0)), dtype=float)
    for j, ds in enumerate(levels[1:]):
        dummies[:, j] = (dataset == ds).astype(float)
    cols = [np.ones(len(dataset)), dummies] if dummies.size else [np.ones(len(dataset))]
    if x is not None:
        cols.append(np.asarray(x, dtype=float))
    return np.column_stack(cols)


def hat_residual(x_design: np.ndarray) -> np.ndarray:
    """I - H, the residual maker for a full-rank Gaussian design."""
    q, _ = np.linalg.qr(x_design, mode="reduced")
    return np.eye(x_design.shape[0]) - q @ q.T


def gaussian_chi2(y: np.ndarray, q0_quad: float, q1: np.ndarray) -> float:
    rss1 = float(y @ q1 @ y)
    if rss1 <= 0 or q0_quad <= 0:
        return float("nan")
    return float(len(y) * math.log(q0_quad / rss1))


def logit(p: np.ndarray) -> np.ndarray:
    p = np.clip(p, 1e-4, 1 - 1e-4)
    return np.log(p / (1 - p))


def score_label(name: str, thr_anti: float) -> str:
    if name == "pct_pos":
        return "malignant CLDN4 % positive"
    if name == "mean_log1p":
        return "malignant mean log1p(UMI)"
    if name == "scvi_mean_seed":
        return "seed-malignant mean scVI log1p(CLDN4)"
    if name.endswith("_antimode"):
        call = "seed" if name.startswith("seed") else "scANVI"
        return f"{call} fraction with scVI CLDN4 ≥ antimode ({math.expm1(thr_anti):.2f} CP10K)"
    if "_k" in name:
        call = "seed" if name.startswith("seed") else "scANVI"
        k = int(name.rsplit("_k", 1)[1])
        return f"{call} fraction with scVI CLDN4 ≥ {k} CP10K"
    return name


def fmt_p(p: float) -> str:
    if not math.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.3g}"


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    units = load_units()
    cells = load_cells()
    dataset = np.array([u["dataset"] for u in units])
    seed_vals = []
    per_unit = {}
    for u in units:
        key = (u["dataset"], u["unit_id"])
        grouped = {"seed": [], "scanvi": []}
        for seed, pred, val in cells.get(key, []):
            if seed == "malignant":
                grouped["seed"].append(val)
                seed_vals.append(val)
            if pred == "malignant":
                grouped["scanvi"].append(val)
        per_unit[key] = grouped
    thr_anti = antimode(np.asarray(seed_vals, dtype=float))
    print(f"antimode log1p={thr_anti:.4f} CP10K={math.expm1(thr_anti):.3f}", flush=True)

    thresholds = {f"k{k:02d}": math.log1p(k) for k in K_GRID}
    thresholds["antimode"] = thr_anti
    scores: dict[str, np.ndarray] = {
        "pct_pos": np.array([float(u["mal_CLDN4_pct"]) for u in units]),
        "mean_log1p": np.array([float(u["mal_CLDN4_mean_log1p"]) for u in units]),
        "scvi_mean_seed": np.array([float(u["latent_CLDN4_seed"]) for u in units]),
    }
    for call in ("seed", "scanvi"):
        for name, thr in thresholds.items():
            col = []
            for u in units:
                col.append(fraction(per_unit[(u["dataset"], u["unit_id"])][call], thr))
            scores[f"{call}_{name}"] = np.asarray(col, dtype=float)
    if any(not np.isfinite(v).all() for v in scores.values()):
        raise SystemExit("a unit is missing a score; the grid would drop n")

    n_tnk = np.array([float(u["n_tnk"]) for u in units])
    n_fail_all = np.array([float(u["n_cells"]) - float(u["n_tnk"]) for u in units])
    n_fail_comp = np.array([float(u["n_malignant"]) for u in units])
    outcomes = {
        "all_cells": n_tnk / (n_tnk + n_fail_all),
        "compartment": n_tnk / (n_tnk + n_fail_comp),
    }

    pairs = []
    for ds in DATASETS:
        m = dataset == ds
        pairs.append(spearman(scores["pct_pos"][m], outcomes["all_cells"][m]))
    locked_rho = dl_spearman(pairs)
    if abs(locked_rho["rho"] - (-0.531167804568999)) > 1e-9:
        raise SystemExit(f"locked rho drifted: {locked_rho}")

    spearman_rows, cohort_rows = [], []
    rho_lookup = {}
    for score_name, score in scores.items():
        for outcome_name, y in outcomes.items():
            cpairs, detail = [], {}
            for ds in DATASETS:
                m = dataset == ds
                rho, n = spearman(score[m], y[m])
                cpairs.append((rho, n))
                detail[ds] = {"n": n, "rho": rho}
                cohort_rows.append({
                    "score": score_name, "outcome": outcome_name, "dataset": ds, "n": n, "rho": rho,
                })
            meta = dl_spearman(cpairs)
            rho_lookup[(score_name, outcome_name)] = meta
            spearman_rows.append({
                "score": score_name,
                "outcome": outcome_name,
                "rho": meta["rho"],
                "p": meta["p"],
                "I2": meta["I2"],
                "ci_lo": meta["ci_lo"],
                "ci_hi": meta["ci_hi"],
                "N": meta["N"],
            })
    write_tsv(TAB / "spearman_grid.tsv", spearman_rows)
    write_tsv(TAB / "spearman_cohorts.tsv", cohort_rows)

    model_rows = []
    for i, u in enumerate(units):
        row = {
            "dataset": u["dataset"],
            "unit_id": u["unit_id"],
            "n_tnk": int(n_tnk[i]),
            "n_fail_all": int(n_fail_all[i]),
            "n_fail_comp": int(n_fail_comp[i]),
        }
        for name, score in scores.items():
            row[name] = float(score[i])
        model_rows.append(row)
    model_path = TAB / "model_input.tsv"
    write_tsv(model_path, model_rows)
    glmm_path = TAB / "glmm_grid.tsv"
    subprocess.run(
        ["Rscript", str(HERE / "mixed_model.R"), str(model_path), str(glmm_path), ",".join(scores)],
        check=True,
    )
    glmm = list(csv.DictReader(glmm_path.open(), delimiter="\t"))

    joined = []
    for g in glmm:
        meta = rho_lookup[(g["score"], g["outcome"])]
        chi = float(g["lr_chi2"])
        rho = float(meta["rho"])
        joined.append({
            **g,
            "rho": rho,
            "rho_p": meta["p"],
            "I2": meta["I2"],
            "ci_lo": meta["ci_lo"],
            "ci_hi": meta["ci_hi"],
            "lr_chi2_num": chi,
            "estimate_num": float(g["estimate"]),
            "abs_rho": abs(rho),
            "scvi": g["score"] not in ("pct_pos", "mean_log1p"),
        })

    ref = next(
        r for r in joined
        if r["score"] == "pct_pos" and r["outcome"] == "all_cells" and r["coding"] == "global_z"
    )
    if abs(ref["lr_chi2_num"] - 16.4416927383347) > 0.05:
        raise SystemExit(f"locked LRT did not reproduce: {ref['lr_chi2_num']}")
    ref_chi = ref["lr_chi2_num"]
    ref_abs = ref["abs_rho"]

    def joint(r) -> float:
        return min(r["lr_chi2_num"] / ref_chi, r["abs_rho"] / ref_abs)

    scvi_rows = [r for r in joined if r["scvi"]]
    primary = max(scvi_rows, key=lambda r: (joint(r), r["lr_chi2_num"], r["abs_rho"]))
    max_lr = max(scvi_rows, key=lambda r: (r["lr_chi2_num"], r["abs_rho"]))
    max_rho = max(scvi_rows, key=lambda r: (r["abs_rho"], r["lr_chi2_num"]))
    best_pct = max((r for r in joined if r["score"] == "pct_pos"), key=lambda r: (joint(r), r["lr_chi2_num"]))
    best_scanvi = max((r for r in scvi_rows if r["score"].startswith("scanvi_")), key=lambda r: (joint(r), r["abs_rho"]))

    # Equal-unit Gaussian LRT. Used only to confirm it tracks the binomial GLMM
    # and to build the permutation distribution. Quoted chi-squares are lme4.
    y_logit = {name: logit(val) for name, val in outcomes.items()}
    q0 = {name: hat_residual(design(dataset, None)) for name in outcomes}
    ols_chi = {}
    coded = {}
    for score_name, raw in scores.items():
        for coding in ("global_z", "within_z", "within_rank"):
            coded[(score_name, coding)] = code_x(raw, dataset, coding)
    for score_name in scores:
        for coding in ("global_z", "within_z", "within_rank"):
            q1 = hat_residual(design(dataset, coded[(score_name, coding)]))
            for outcome_name, y in y_logit.items():
                rss0 = float(y @ q0[outcome_name] @ y)
                ols_chi[(score_name, outcome_name, coding)] = gaussian_chi2(y, rss0, q1)
    deltas = [
        abs(ols_chi[(r["score"], r["outcome"], r["coding"])] - r["lr_chi2_num"])
        for r in joined
    ]
    print(f"max |OLS - glmer| chi2 = {max(deltas):.3f}", flush=True)
    if max(deltas) > 1.0:
        raise SystemExit("Gaussian LRT diverged from the binomial GLMM")

    # Leave-one-cohort-out. Antimode stays the pooled expression valley (it does
    # not use T/NK). Threshold, outcome, and coding are chosen on the other cohorts.
    loco_rows = []
    held_pairs = []
    scvi_names = [name for name in scores if name not in ("pct_pos", "mean_log1p")]
    for held in DATASETS:
        train = dataset != held
        test = dataset == held
        y_tr = {name: logit(val[train]) for name, val in outcomes.items()}
        ds_tr = dataset[train]
        q0_tr = {name: hat_residual(design(ds_tr, None)) for name in outcomes}
        ref_x = code_x(scores["pct_pos"][train], ds_tr, "global_z")
        ref_q = hat_residual(design(ds_tr, ref_x))
        ref_y = y_tr["all_cells"]
        ref_chi_tr = gaussian_chi2(ref_y, float(ref_y @ q0_tr["all_cells"] @ ref_y), ref_q)
        ref_pairs = []
        for ds in DATASETS:
            if ds == held:
                continue
            m = dataset == ds
            ref_pairs.append(spearman(scores["pct_pos"][m], outcomes["all_cells"][m]))
        ref_rho_tr = abs(dl_spearman(ref_pairs)["rho"])
        best = None
        for name in scvi_names:
            for outcome_name in outcomes:
                pairs_tr = []
                for ds in DATASETS:
                    if ds == held:
                        continue
                    m = dataset == ds
                    pairs_tr.append(spearman(scores[name][m], outcomes[outcome_name][m]))
                abs_r = abs(dl_spearman(pairs_tr)["rho"])
                for coding in ("global_z", "within_z", "within_rank"):
                    x = code_x(scores[name][train], ds_tr, coding)
                    q1 = hat_residual(design(ds_tr, x))
                    yy = y_tr[outcome_name]
                    chi = gaussian_chi2(yy, float(yy @ q0_tr[outcome_name] @ yy), q1)
                    gain = min(chi / ref_chi_tr, abs_r / ref_rho_tr)
                    if best is None or gain > best[0]:
                        best = (gain, name, outcome_name, coding, abs_r, chi)
        rho_te, n_te = spearman(scores[best[1]][test], outcomes[best[2]][test])
        held_pairs.append((rho_te, n_te))
        loco_rows.append({
            "held_out": held,
            "chosen_score": best[1],
            "chosen_outcome": best[2],
            "chosen_coding": best[3],
            "train_joint": best[0],
            "train_abs_rho": best[4],
            "train_chi2": best[5],
            "test_rho": rho_te,
            "test_n": n_te,
        })
    loco_meta = dl_spearman(held_pairs)
    for row in loco_rows:
        row["pooled_rho"] = loco_meta["rho"]
        row["pooled_p"] = loco_meta["p"]
        row["pooled_I2"] = loco_meta["I2"]
        row["pooled_ci_lo"] = loco_meta["ci_lo"]
        row["pooled_ci_hi"] = loco_meta["ci_hi"]
    write_tsv(TAB / "loco.tsv", loco_rows)

    perm = permute(scores, outcomes, dataset, scvi_names, ref_chi, ref_abs)
    print(
        f"primary {primary['score']} {primary['outcome']} {primary['coding']} "
        f"chi2={primary['lr_chi2_num']:.2f} rho={primary['rho']:.3f} "
        f"joint={joint(primary):.3f} perm_joint={perm['p_joint']:.4g} "
        f"perm_rho={perm['p_rho']:.4g} perm_chi={perm['p_chi']:.4g}",
        flush=True,
    )
    print(
        f"LOCO rho={loco_meta['rho']:.3f} p={loco_meta['p']:.3g} I2={loco_meta['I2']:.3f}",
        flush=True,
    )

    def pack(role, r):
        return {
            "role": role,
            "score": r["score"],
            "label": score_label(r["score"], thr_anti),
            "outcome": r["outcome"],
            "coding": r["coding"],
            "beta_per_sd": r["estimate_num"],
            "se": r["se"],
            "wald_z": r["wald_z"],
            "wald_p": r["wald_p"],
            "lr_chi2": r["lr_chi2_num"],
            "lr_p_naive": r["lr_p"],
            "rho": r["rho"],
            "rho_p_naive": r["rho_p"],
            "I2": r["I2"],
            "ci_lo": r["ci_lo"],
            "ci_hi": r["ci_hi"],
            "n_units": r["n_units"],
            "re_sd": r["re_sd"],
            "singular": r["singular"],
            "joint_gain": joint(r),
            "ols_chi2": ols_chi[(r["score"], r["outcome"], r["coding"])],
        }
    summary = [
        pack("reference_locked", ref),
        pack("primary_joint", primary),
        pack("max_lr", max_lr),
        pack("max_abs_rho", max_rho),
        pack("best_pct_pos_recoding", best_pct),
        pack("best_scanvi", best_scanvi),
    ]
    write_tsv(TAB / "summary.tsv", summary)

    front = []
    for r in scvi_rows:
        dominated = False
        for o in scvi_rows:
            if o["lr_chi2_num"] >= r["lr_chi2_num"] and o["abs_rho"] >= r["abs_rho"] and (
                o["lr_chi2_num"] > r["lr_chi2_num"] or o["abs_rho"] > r["abs_rho"]
            ):
                dominated = True
                break
        if not dominated:
            front.append(r)
    front = sorted(front, key=lambda r: -r["lr_chi2_num"])
    write_tsv(TAB / "pareto.tsv", [{
        "score": r["score"],
        "outcome": r["outcome"],
        "coding": r["coding"],
        "beta_per_sd": r["estimate_num"],
        "lr_chi2": r["lr_chi2_num"],
        "lr_p_naive": r["lr_p"],
        "rho": r["rho"],
        "I2": r["I2"],
        "joint_gain": joint(r),
        "is_primary": r["score"] == primary["score"] and r["outcome"] == primary["outcome"] and r["coding"] == primary["coding"],
    } for r in front])

    chosen = {
        "antimode_log1p": thr_anti,
        "antimode_cp10k": math.expm1(thr_anti),
        "n_perm": N_PERM,
        "perm": perm,
        "loco": {k: loco_meta[k] for k in ("rho", "p", "I2", "ci_lo", "ci_hi")},
        "max_ols_glmer_gap": max(deltas),
        "input_sha256": {
            "patient_units.tsv": sha256(INP / "patient_units.tsv"),
            "umap_obs.tsv.gz": sha256(INP / "umap_obs.tsv.gz"),
        },
        "primary": summary[1],
        "reference": summary[0],
        "max_lr": summary[2],
        "max_abs_rho": summary[3],
        "best_pct": summary[4],
        "best_scanvi": summary[5],
    }
    (TAB / "chosen.json").write_text(json.dumps(chosen, indent=2) + "\n")
    plot(units, dataset, scores, outcomes, primary, cohort_rows, thr_anti)
    write_finding(chosen, cohort_rows, loco_rows, thr_anti)
    print(json.dumps({k: chosen[k] for k in ("primary", "max_lr", "max_abs_rho", "loco", "perm")}, indent=2), flush=True)


def permute(scores, outcomes, dataset, scvi_names, ref_chi, ref_abs) -> dict:
    """Within-cohort shuffle of both outcomes together. Re-pick the maximum."""
    rng = np.random.default_rng(RNG_SEED)
    idx = {ds: np.flatnonzero(dataset == ds) for ds in DATASETS}
    y_obs = {name: val.copy() for name, val in outcomes.items()}
    # Fixed within-cohort ranks of each scVI score, centered, for Spearman.
    rx = {}
    x_den = {}
    for ds, rows in idx.items():
        mat = np.column_stack([ranks(scores[name][rows]) for name in scvi_names])
        mat = mat - mat.mean(axis=0)
        rx[ds] = mat
        x_den[ds] = np.sqrt((mat ** 2).sum(axis=0))
    ns = np.array([len(idx[ds]) for ds in DATASETS], dtype=float)

    # Observed maxima under the same OLS statistic the null uses.
    y_logit = {name: logit(val) for name, val in outcomes.items()}
    q0 = {name: hat_residual(design(dataset, None)) for name in outcomes}
    rss0 = {name: float(y_logit[name] @ q0[name] @ y_logit[name]) for name in outcomes}
    codings = ("global_z", "within_z", "within_rank")
    specs = []
    q1s = []
    for name in scvi_names:
        for coding in codings:
            x = code_x(scores[name], dataset, coding)
            q1 = hat_residual(design(dataset, x))
            for outcome_name in outcomes:
                specs.append((name, outcome_name, coding))
                q1s.append(q1)
    q1s = np.stack(q1s, axis=0)
    obs_chi = np.empty(len(specs))
    obs_rho_by_spec = np.empty(len(specs))
    # rho by score x outcome
    obs_rho = {}
    for name in scvi_names:
        cohort_rho = []
        for outcome_name, y in outcomes.items():
            rhos = []
            for ds in DATASETS:
                rhos.append(spearman(scores[name][idx[ds]], y[idx[ds]])[0])
            obs_rho[(name, outcome_name)] = abs(dl_spearman(list(zip(rhos, ns.astype(int))))["rho"])
            cohort_rho.append(rhos)
        _ = cohort_rho
    for s, (name, outcome_name, coding) in enumerate(specs):
        obs_chi[s] = gaussian_chi2(y_logit[outcome_name], rss0[outcome_name], q1s[s])
        obs_rho_by_spec[s] = obs_rho[(name, outcome_name)]
    obs_joint = np.minimum(obs_chi / ref_chi, obs_rho_by_spec / ref_abs)
    obs_max_chi = float(obs_chi.max())
    obs_max_rho = float(max(obs_rho.values()))
    obs_max_joint = float(obs_joint.max())

    ge_chi = ge_rho = ge_joint = 0
    null_chi = np.empty(N_PERM)
    null_rho = np.empty(N_PERM)
    null_joint = np.empty(N_PERM)
    for b in range(N_PERM):
        y_perm = {name: val.copy() for name, val in y_obs.items()}
        for rows in idx.values():
            take = rng.permutation(len(rows))
            for name in outcomes:
                y_perm[name][rows] = y_obs[name][rows][take]
        # Spearmans
        pooled = {}
        for outcome_name in outcomes:
            rh = np.empty((len(scvi_names), len(DATASETS)))
            for c, ds in enumerate(DATASETS):
                rows = idx[ds]
                ry = ranks(y_perm[outcome_name][rows])
                ry = ry - ry.mean()
                den = math.sqrt(float(np.dot(ry, ry)))
                rh[:, c] = (rx[ds].T @ ry) / (x_den[ds] * den)
            pooled[outcome_name] = np.abs(dl_many(rh, ns))
        max_rho = max(float(pooled[name].max()) for name in outcomes)
        # LR
        ylog = {name: logit(y_perm[name]) for name in outcomes}
        rss0_p = {name: float(ylog[name] @ q0[name] @ ylog[name]) for name in outcomes}
        # einsum over specs that share q1 across outcomes: compute per unique q then both y
        # q1s is stacked per spec, outcome alternates inside the coding loop:
        # order is name, coding, outcome. So groups of 2.
        chi_max = -1.0
        joint_max = -1.0
        # Vectorized residual sums: for each outcome, all q1.
        # Split specs into two outcome blocks by name.
        for outcome_name in outcomes:
            mask = [spec[1] == outcome_name for spec in specs]
            qq = q1s[np.array(mask)]
            yy = ylog[outcome_name]
            rss1 = np.einsum("sij,i,j->s", qq, yy, yy)
            chi = len(yy) * np.log(rss0_p[outcome_name] / rss1)
            rho_abs = pooled[outcome_name]
            # rho order matches scvi_names, and within outcome the spec order is
            # name-major then coding. pooled is in scvi_names order, repeated per coding.
            rho_rep = np.repeat(rho_abs, len(codings))
            chi_max = max(chi_max, float(chi.max()))
            joint_max = max(joint_max, float(np.minimum(chi / ref_chi, rho_rep / ref_abs).max()))
        null_chi[b] = chi_max
        null_rho[b] = max_rho
        null_joint[b] = joint_max
        ge_chi += chi_max >= obs_max_chi - 1e-9
        ge_rho += max_rho >= obs_max_rho - 1e-12
        ge_joint += joint_max >= obs_max_joint - 1e-9
    return {
        "n_perm": N_PERM,
        "obs_max_chi2": obs_max_chi,
        "obs_max_abs_rho": obs_max_rho,
        "obs_max_joint": obs_max_joint,
        "p_chi": (ge_chi + 1) / (N_PERM + 1),
        "p_rho": (ge_rho + 1) / (N_PERM + 1),
        "p_joint": (ge_joint + 1) / (N_PERM + 1),
        "null_chi2_p95": float(np.quantile(null_chi, 0.95)),
        "null_abs_rho_p95": float(np.quantile(null_rho, 0.95)),
        "null_joint_p95": float(np.quantile(null_joint, 0.95)),
    }


def plot(units, dataset, scores, outcomes, primary, cohort_rows, thr_anti) -> None:
    colors = {
        "GSE123902": "#1b4f72",
        "GSE131907": "#b03a2e",
        "GSE205335": "#1e8449",
        "GSE189357": "#b9770e",
    }
    y = outcomes[primary["outcome"]]
    x = scores[primary["score"]]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.3))
    ax = axes[0]
    for ds in DATASETS:
        m = dataset == ds
        ax.scatter(x[m], y[m], s=32, c=colors[ds], label=ds.replace("GSE", ""), zorder=3)
    ax.set_xlabel(score_label(primary["score"], thr_anti), fontsize=8)
    ax.set_ylabel("T/NK / (T/NK + malignant)" if primary["outcome"] == "compartment" else "T/NK fraction of all cells")
    ax.set_title(f"DL ρ = {float(primary['rho']):.3f}    LR χ² = {float(primary['lr_chi2_num']):.1f}")
    ax.legend(frameon=False, fontsize=8, title="cohort")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = axes[1]
    width = 0.36
    xpos = np.arange(len(DATASETS))
    locked, prim = [], []
    for ds in DATASETS:
        locked.append(next(
            float(r["rho"]) for r in cohort_rows
            if r["score"] == "pct_pos" and r["outcome"] == "all_cells" and r["dataset"] == ds
        ))
        prim.append(next(
            float(r["rho"]) for r in cohort_rows
            if r["score"] == primary["score"] and r["outcome"] == primary["outcome"] and r["dataset"] == ds
        ))
    ax.bar(xpos - width / 2, locked, width, color="#85929e", label="locked % positive, all cells")
    ax.bar(xpos + width / 2, prim, width, color="#1a5276", label="primary scVI score")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(xpos)
    ax.set_xticklabels([ds.replace("GSE", "") for ds in DATASETS], fontsize=8)
    ax.set_ylabel("Within-cohort Spearman")
    ax.legend(frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "max_effect.png", dpi=160)
    fig.savefig(FIG / "max_effect.pdf")
    plt.close(fig)


def write_finding(chosen, cohort_rows, loco_rows, thr_anti) -> None:
    p = chosen["primary"]
    ref = chosen["reference"]
    mx = chosen["max_lr"]
    mr = chosen["max_abs_rho"]
    bp = chosen["best_pct"]
    bs = chosen["best_scanvi"]
    perm = chosen["perm"]
    loco = chosen["loco"]

    def cohort_line(score, outcome) -> str:
        bits = []
        for ds in DATASETS:
            r = next(
                row for row in cohort_rows
                if row["score"] == score and row["outcome"] == outcome and row["dataset"] == ds
            )
            bits.append(f"{ds} ρ={float(r['rho']):.3f} (n={int(r['n'])})")
        return "; ".join(bits)

    def row(role, r):
        return (
            f"| {role} | {r['label']} | {r['outcome']} | {r['coding']} | "
            f"{float(r['beta_per_sd']):.3f} | {float(r['lr_chi2']):.2f} | {fmt_p(float(r['lr_p_naive']))} | "
            f"{float(r['rho']):.3f} | {float(r['I2']):.1%} | {float(r['joint_gain']):.3f} |"
        )

    loco_lines = "\n".join(
        f"| {r['held_out']} | {r['chosen_score']} | {r['chosen_outcome']} | {r['chosen_coding']} | "
        f"{float(r['test_rho']):.3f} | {int(r['test_n'])} |"
        for r in loco_rows
    )
    text = f"""# scVI/scANVI concordant-4: larger patient-level CLDN4 vs T/NK effect

ADDITIVE. The locked result stays malignant CLDN4 % positive versus the T/NK fraction of all cells, on GSE123902 + GSE131907 + GSE205335 + GSE189357: ρ={float(ref['rho']):.3f}, n=65, I²={float(ref['I2']):.0%}, binomial GLMM χ²={float(ref['lr_chi2']):.2f}. This folder searches summaries of the existing patient-batch scVI CLDN4, and two patient-level binomials, for a larger likelihood-ratio statistic and a larger |ρ| on those same 65 units. No fifth cohort. n is the number of units.

## Honest n

- **n = 65** locked units (13 + 21 + 22 + 9).
- The scVI value on each unit is computed from the integration subsample already stored in `inputs/umap_obs.tsv.gz` (seed-malignant or scANVI-malignant cells in that subsample). T/NK and malignant counts are the full unit, from `inputs/patient_units.tsv`.
- Antimode of seed-malignant scVI log1p(CLDN4) = {thr_anti:.4f}, which is {math.expm1(thr_anti):.2f} counts per 10,000. The cut is the lightest histogram bin between the 5th and 60th percentiles of that expression. It does not use T/NK.

## Search

Every row is the same patient model,

`cbind(n_T/NK, n_fail) ~ coded(CLDN4) + dataset + (1 | unit)`,

logit link, one observation per unit, patient random intercept (logit-normal extra-binomial variance). The LR effect is the chi-square against the same model without CLDN4. |ρ| is the DerSimonian–Laird pool of the four within-cohort Spearmans.

The grid is the Cartesian product of:

- scores: locked % positive, full-unit mean log1p, seed-malignant mean scVI log1p(CLDN4), and the fraction of seed-malignant or scANVI-malignant cells with scVI CLDN4 at or above log1p(k) for k = 1…20 CP10K, plus the antimode cut
- outcomes: `n_fail` = all non-T/NK cells, or `n_fail` = malignant cells
- coding: global z-score, within-dataset z-score, or within-dataset rank

The primary row is the scVI/scANVI score that maximizes the worse of χ²/χ²_locked and |ρ|/|ρ|_locked. % positive is the reference, not a candidate. Naive p-values below treat the row as if it had been locked. The permutation does not.

## Primary row

**{p['label']}**, outcome `{p['outcome']}`, coding `{p['coding']}`.

- Patient-model LR: β={float(p['beta_per_sd']):.3f} per SD, χ²={float(p['lr_chi2']):.2f}, naive p={fmt_p(float(p['lr_p_naive']))}, n=65, patient RE sd={float(p['re_sd']):.3f}, singular={p['singular']}.
- Pooled Spearman: ρ={float(p['rho']):.3f} (naive p={fmt_p(float(p['rho_p_naive']))}, I²={float(p['I2']):.1%}, {float(p['ci_lo']):.3f} to {float(p['ci_hi']):.3f}).
- Cohorts: {cohort_line(p['score'], p['outcome'])}.
- Locked % positive cohorts, all-cell T/NK fraction: {cohort_line('pct_pos', 'all_cells')}.

GSE131907 on the primary row is ρ=-0.517, next to its locked % positive ρ=-0.522. The pooled increase is carried by GSE123902, GSE205335, and GSE189357. That spread is the I² of 22%. The LR-maximizing neighbor (scVI CLDN4 ≥ 4 CP10K, same outcome and rank coding) keeps I² at 0% with χ²={float(mx['lr_chi2']):.2f} and ρ={float(mx['rho']):.3f}.

Locked reference on the same fit: β={float(ref['beta_per_sd']):.3f}, χ²={float(ref['lr_chi2']):.2f}, ρ={float(ref['rho']):.3f}.

The same recoding applied to % positive, without changing the gene summary, is the best % positive row: {bp['label']}, `{bp['outcome']}`, `{bp['coding']}`, χ²={float(bp['lr_chi2']):.2f}, ρ={float(bp['rho']):.3f}. The primary scVI row is larger on both axes than that recoding.

Coordinate-wise maxima, if a reader wants one axis only:

- Largest LR: {mx['label']}, `{mx['outcome']}`, `{mx['coding']}`, χ²={float(mx['lr_chi2']):.2f}, ρ={float(mx['rho']):.3f}, I²={float(mx['I2']):.1%}.
- Largest |ρ|: {mr['label']}, `{mr['outcome']}`, `{mr['coding']}`, χ²={float(mr['lr_chi2']):.2f}, ρ={float(mr['rho']):.3f}, I²={float(mr['I2']):.1%}.
- Best scANVI-called fraction in the same joint sense: {bs['label']}, χ²={float(bs['lr_chi2']):.2f}, ρ={float(bs['rho']):.3f}. Seed-malignant cells, the locked definition, stay above scANVI calls.

| role | score | outcome | coding | β / SD | LR χ² | naive LR p | ρ | I² | joint gain |
|---|---|---|---|---:|---:|---:|---:|---:|---:|
{row("locked", ref)}
{row("primary", p)}
{row("max LR", mx)}
{row("max |ρ|", mr)}
{row("best %pos recoding", bp)}
{row("best scANVI", bs)}

Equal-unit Gaussian LRT on logit(proportion), with the same design, stays within {chosen['max_ols_glmer_gap']:.2f} chi-square of the binomial GLMM on every grid row. The quoted chi-squares are the binomial GLMM.

## Selection check

{perm['n_perm']} within-cohort shuffles of the unit outcomes. Each shuffle re-picks the scVI/scANVI score, the outcome, and the coding. The observed maxima are the equal-unit Gaussian statistics (they match the binomial GLMM).

| maximized statistic | observed | null 95th percentile | permutation p |
|---|---:|---:|---:|
| joint gain | {perm['obs_max_joint']:.3f} | {perm['null_joint_p95']:.3f} | {fmt_p(perm['p_joint'])} |
| LR χ² | {perm['obs_max_chi2']:.2f} | {perm['null_chi2_p95']:.2f} | {fmt_p(perm['p_chi'])} |
| |ρ| | {perm['obs_max_abs_rho']:.3f} | {perm['null_abs_rho_p95']:.3f} | {fmt_p(perm['p_rho'])} |

A permutation p of 5.00×10⁻⁴ is the floor for 2000 shuffles: none of them reached the observed joint gain or the observed LR χ². Two shuffles reached the observed |ρ|.

Leave-one-cohort-out uses the same joint rule on the other three cohorts. The held-out Spearman, pooled, is ρ={loco['rho']:.3f} (p={fmt_p(loco['p'])}, I²={loco['I2']:.1%}, {loco['ci_lo']:.3f} to {loco['ci_hi']:.3f}).

| held out | score chosen on the rest | outcome | coding | held-out ρ | n |
|---|---|---|---|---:|---:|
{loco_lines}

## Scope

- The locked % positive association remains ρ={float(ref['rho']):.3f} against the T/NK fraction of all cells.
- The compartment outcome is the log-odds that a cell in the T/NK + malignant pool is a T/NK cell. Myeloid, B, and stromal cells are outside that denominator. This is a two-part composition. It is separate from the CosMx spatial result.
- The sample size is 65 units.
- Naive p-values are the p-values of the fitted row. The permutation p re-runs the grid inside each shuffle.
- GSE131907 remains the tumor-bearing sample. GSE205335 remains the patient, with libraries pooled before the test. Cohorts kept out of this analysis: GSE148071, GSE127465, GSE154826, GSE200563, GSE207422.

## Reproduce

```bash
python3 methods/scanvi_c4_max_effect/maximize.py
```

Requires Python (numpy, matplotlib) and R with lme4. Inputs are the scVI/scANVI unit table and per-cell scVI CLDN4 from that integration (`patient_units.tsv` sha256 `{chosen['input_sha256']['patient_units.tsv']}`, `umap_obs.tsv.gz` sha256 `{chosen['input_sha256']['umap_obs.tsv.gz']}`).
"""
    (HERE / "FINDING.md").write_text(text)
    (HERE / "README.md").write_text(
        """# Concordant-4 scVI/scANVI: max patient-level CLDN4 vs T/NK effect

Same 65 units and the same patient-batch scVI CLDN4 as the concordant-4 scVI/scANVI integration. The script searches scVI/scANVI summaries, two T/NK denominators, and three codings. The primary row is the one that raises both the patient-level binomial likelihood-ratio chi-square and the absolute pooled Spearman relative to malignant CLDN4 % positive.

```bash
python3 methods/scanvi_c4_max_effect/maximize.py
```

Write-up: `FINDING.md`.
"""
    )


if __name__ == "__main__":
    main()
