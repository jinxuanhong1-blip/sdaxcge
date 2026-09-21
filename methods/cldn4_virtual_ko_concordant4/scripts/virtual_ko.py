#!/usr/bin/env python3
"""CLDN4 virtual knockout on concordant-4 malignant epithelial cells.

CellOracle's published simulation is a one-step shift from a ridge GRN:
after a regulator is set to the knockout value, each target moves by
beta_regulator * (z_knockout - z_observed). The stock package builds that
GRN from motif links among transcription factors. CLDN4 is not a TF, so a
motif base GRN has no outgoing CLDN4 edges and the package knockout is
identically zero. It cannot answer the question.

What runs here is that same one-step ridge shift, with CLDN4 added to the
regulator list next to expressed TFs. The GRN is fit inside each patient
(or sample) that has enough malignant cells. The patient is the biological
replicate. A housekeeping gene matched on detection is knocked out in the
same model as a negative control.

This is a co-expression GRN simulation. It is not a knockdown and it does
not use the private KL matrices.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

ROOT = Path(__file__).resolve().parents[1]
CACHE = Path("/tmp/cldn4_vko")
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

MIN_CELLS = 80
MAX_CELLS = 400
MIN_CLDN4_PCT = 5.0
MIN_REG_DETECT = 0.10
MIN_PROGRAM_GENES = 5
N_BOOT = 10
N_PERM = 40
ALPHAS = (1.0, 10.0, 100.0, 1000.0)
SEED = 1

PRIMARY = ("IFN", "MHC_I_APM", "chemokine", "barrier_ligands")
SECONDARY = ("chemokine_ligands", "TJ_broad")
EXPECT_UP = {"IFN", "MHC_I_APM", "chemokine", "chemokine_ligands"}
EXPECT_DOWN = {"barrier_ligands", "TJ_broad"}
COHORT_ORDER = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
COHORT_COLOR = {
    "GSE123902": "#1b4f72",
    "GSE131907": "#117a65",
    "GSE205335": "#b9770e",
    "GSE189357": "#6c3483",
}


def load_programs() -> dict[str, list[str]]:
    raw = json.loads((ROOT / "data" / "programs.json").read_text())
    return {k: [g.upper() for g in v] for k, v in raw.items() if isinstance(v, list)}


def tf_list() -> list[str]:
    return [
        g.strip().upper()
        for g in (ROOT / "data" / "human_tfs.txt").read_text().splitlines()
        if g.strip()
    ]


def seed_for(text: str) -> int:
    acc = SEED
    for ch in text:
        acc = (acc * 131 + ord(ch)) % 2_147_483_647
    return int(acc)


def fit_beta(Xz: np.ndarray, Y: np.ndarray, alpha: float) -> np.ndarray:
    """Ridge map from standardized regulators to log1p targets. Shape (p, t)."""
    y_mean = Y.mean(axis=0)
    Yc = Y - y_mean
    p = Xz.shape[1]
    xtx = Xz.T @ Xz
    xtx.flat[:: p + 1] += alpha
    xty = Xz.T @ Yc
    try:
        B = np.linalg.solve(xtx, xty)
    except np.linalg.LinAlgError:
        B = np.linalg.lstsq(xtx, xty, rcond=None)[0]
    return B


def choose_alpha(Xz: np.ndarray, Y: np.ndarray, rng: np.random.Generator) -> tuple[float, float]:
    n = Xz.shape[0]
    if n < 30:
        return 100.0, float("nan")
    folds = np.array_split(rng.permutation(n), 3)
    best_a, best_s = ALPHAS[0], -1e18
    for alpha in ALPHAS:
        fold_scores = []
        for k in range(3):
            te = folds[k]
            tr = np.concatenate([folds[i] for i in range(3) if i != k])
            if len(tr) < 10 or len(te) < 5:
                continue
            B = fit_beta(Xz[tr], Y[tr], alpha)
            pred = Xz[te] @ B + Y[tr].mean(axis=0)
            ss_res = ((Y[te] - pred) ** 2).sum(axis=0)
            ss_tot = ((Y[te] - Y[te].mean(axis=0)) ** 2).sum(axis=0)
            ok = ss_tot > 1e-8
            if not np.any(ok):
                continue
            fold_scores.append(float(np.mean(1.0 - ss_res[ok] / ss_tot[ok])))
        if not fold_scores:
            continue
        score = float(np.mean(fold_scores))
        if score > best_s:
            best_s = score
            best_a = alpha
    return float(best_a), float(best_s)


def zscore_cols(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mu = X.mean(axis=0)
    sd = X.std(axis=0, ddof=1)
    sd_safe = np.where(sd > 1e-8, sd, 1.0)
    return (X - mu) / sd_safe, mu, sd


def program_delta(beta_row: np.ndarray, gene_idx: np.ndarray, scale: float) -> tuple[float, int]:
    if gene_idx.size < MIN_PROGRAM_GENES or not np.isfinite(scale):
        return float("nan"), int(gene_idx.size)
    d = float(np.mean(beta_row[gene_idx]) * scale)
    return d, int(gene_idx.size)


def knock_scale(log_expr: np.ndarray) -> float:
    """(z of log1p=0) - (z of the observed mean). Mean z is 0, so this is z(0)."""
    sd = float(log_expr.std(ddof=1))
    mean = float(log_expr.mean())
    if sd < 1e-8:
        return float("nan")
    return (0.0 - mean) / sd


def analyze_matrix(
    logx: np.ndarray,
    genes: list[str],
    programs: dict[str, list[str]],
    tfs: set[str],
    rng: np.random.Generator,
    n_perm: int,
) -> dict | None:
    gene_index = {g: i for i, g in enumerate(genes)}
    if "CLDN4" not in gene_index:
        return None
    cldn = logx[:, gene_index["CLDN4"]]
    detect = (logx > 0).mean(axis=0)
    var = logx.var(axis=0)
    cldn_pct = float(100.0 * np.mean(cldn > 0))
    if cldn_pct < MIN_CLDN4_PCT or float(cldn.std(ddof=1)) < 1e-8:
        return None

    program_gene_set = set()
    for name in PRIMARY + SECONDARY:
        program_gene_set.update(programs[name])

    house = [
        g
        for g in programs["housekeeping_candidates"]
        if g in gene_index and g not in program_gene_set and g != "CLDN4"
    ]
    control = None
    if house:
        target_d = detect[gene_index["CLDN4"]]
        control = min(house, key=lambda g: abs(detect[gene_index[g]] - target_d))
        if detect[gene_index[control]] < MIN_REG_DETECT or var[gene_index[control]] <= 0:
            usable = [g for g in house if detect[gene_index[g]] >= MIN_REG_DETECT and var[gene_index[g]] > 0]
            control = min(usable, key=lambda g: abs(detect[gene_index[g]] - target_d)) if usable else None

    n_reg_cap = int(min(40, max(15, logx.shape[0] // 8)))
    tf_cands = []
    for g in tfs:
        j = gene_index.get(g)
        if j is None or g == "CLDN4" or g == control:
            continue
        if detect[j] >= MIN_REG_DETECT and var[j] > 0:
            tf_cands.append((var[j], g))
    tf_cands.sort(reverse=True)
    regulators = [g for _, g in tf_cands[:n_reg_cap]]
    regulators.append("CLDN4")
    if control:
        regulators.append(control)
    # unique, preserve order
    seen = set()
    regulators = [g for g in regulators if not (g in seen or seen.add(g))]

    target_names = []
    for g in genes:
        if g == "CLDN4":
            continue
        if g in program_gene_set and var[gene_index[g]] > 0:
            target_names.append(g)
    if len(regulators) < 5 or len(target_names) < MIN_PROGRAM_GENES:
        return None

    R = np.column_stack([logx[:, gene_index[g]] for g in regulators])
    Y = np.column_stack([logx[:, gene_index[g]] for g in target_names])
    Xz, _, x_sd = zscore_cols(R)
    keep_reg = x_sd > 1e-8
    if keep_reg.sum() < 5 or not keep_reg[regulators.index("CLDN4")]:
        return None
    regulators = [g for g, k in zip(regulators, keep_reg) if k]
    Xz = Xz[:, keep_reg]

    alpha, cv_r2 = choose_alpha(Xz, Y, rng)
    B_sum = np.zeros((Xz.shape[1], Y.shape[1]), dtype=np.float64)
    for _ in range(N_BOOT):
        take = rng.choice(Xz.shape[0], Xz.shape[0], replace=True)
        B_sum += fit_beta(Xz[take], Y[take], alpha)
    B = B_sum / N_BOOT
    B1 = fit_beta(Xz, Y, alpha)

    cldn_row = regulators.index("CLDN4")
    scale_cldn = knock_scale(logx[:, gene_index["CLDN4"]])
    control_row = regulators.index(control) if control in regulators else None
    scale_ctrl = knock_scale(logx[:, gene_index[control]]) if control_row is not None else float("nan")

    t_index = {g: i for i, g in enumerate(target_names)}

    def idx_for(program: str) -> np.ndarray:
        hits = [t_index[g] for g in programs[program] if g in t_index]
        return np.asarray(hits, dtype=int)

    # permutation null of the single-fit CLDN4 row, same scale
    nulls = {name: [] for name in PRIMARY + SECONDARY}
    obs1 = {}
    obs_boot = {}
    n_genes = {}
    for name in PRIMARY + SECONDARY:
        ix = idx_for(name)
        d1, n_g = program_delta(B1[cldn_row], ix, scale_cldn)
        db, _ = program_delta(B[cldn_row], ix, scale_cldn)
        obs1[name] = d1
        obs_boot[name] = db
        n_genes[name] = n_g

    if n_perm > 0:
        base_col = Xz[:, cldn_row].copy()
        for _ in range(n_perm):
            shuffled = base_col.copy()
            rng.shuffle(shuffled)
            Xn = Xz.copy()
            Xn[:, cldn_row] = shuffled
            Bp = fit_beta(Xn, Y, alpha)
            for name in PRIMARY + SECONDARY:
                d, _ = program_delta(Bp[cldn_row], idx_for(name), scale_cldn)
                nulls[name].append(d)

    perm_p = {}
    for name in PRIMARY + SECONDARY:
        arr = np.asarray(nulls[name], dtype=float)
        arr = arr[np.isfinite(arr)]
        obs = obs1[name]
        if arr.size == 0 or not np.isfinite(obs):
            perm_p[name] = float("nan")
        else:
            perm_p[name] = float((1 + np.sum(np.abs(arr) >= abs(obs))) / (1 + arr.size))

    ctrl_delta = {}
    for name in PRIMARY + SECONDARY:
        if control_row is None:
            ctrl_delta[name] = float("nan")
        else:
            d, _ = program_delta(B[control_row], idx_for(name), scale_ctrl)
            ctrl_delta[name] = d

    before = {}
    for name in PRIMARY + SECONDARY:
        genes_p = [g for g in programs[name] if g in t_index]
        if len(genes_p) < MIN_PROGRAM_GENES:
            before[name] = float("nan")
        else:
            cols = [gene_index[g] for g in genes_p]
            before[name] = float(logx[:, cols].mean())

    beta_genes = {target_names[i]: float(B[cldn_row, i]) for i in range(len(target_names))}
    return {
        "alpha": alpha,
        "cv_r2": cv_r2,
        "n_reg": len(regulators),
        "n_targets": len(target_names),
        "control": control or "",
        "cldn4_pct": cldn_pct,
        "cldn4_mean": float(cldn.mean()),
        "scale_cldn": scale_cldn,
        "before": before,
        "delta_boot": obs_boot,
        "delta_single": obs1,
        "delta_control": ctrl_delta,
        "n_genes": n_genes,
        "perm_p": perm_p,
        "beta_genes": beta_genes,
        "regulators": regulators,
    }


def self_test() -> None:
    rng = np.random.default_rng(0)
    n = 240
    g = 30
    genes = [f"TF{i}" for i in range(12)] + ["CLDN4", "GAPDH", "RPLP0"]
    ifn = [f"IFN{i}" for i in range(8)]
    mhc = [f"MHC{i}" for i in range(6)]
    genes = genes + ifn + mhc
    # rewrite a tiny programs object via the real function's expectations by
    # calling analyze_matrix with a local program dict and a local tf set.
    logx = rng.normal(0.4, 0.35, size=(n, len(genes))).astype(np.float64)
    logx = np.clip(logx, 0, None)
    gi = {g: i for i, g in enumerate(genes)}
    # plant a negative partial association: high CLDN4 -> low IFN, low barrier stand-in (MHC used as down in this test? )
    z = (logx[:, gi["CLDN4"]] - logx[:, gi["CLDN4"]].mean()) / logx[:, gi["CLDN4"]].std()
    for gname in ifn:
        logx[:, gi[gname]] = np.clip(1.2 - 0.55 * z + rng.normal(0, 0.05, n), 0, None)
    for gname in mhc:
        logx[:, gi[gname]] = np.clip(0.3 + 0.40 * z + rng.normal(0, 0.05, n), 0, None)
    programs = {
        "IFN": ifn,
        "MHC_I_APM": mhc,
        "chemokine": ifn[:5],
        "barrier_ligands": mhc,
        "chemokine_ligands": ifn[:5],
        "TJ_broad": mhc,
        "housekeeping_candidates": ["GAPDH", "RPLP0"],
    }
    tfs = {f"TF{i}" for i in range(12)}
    # analyze_matrix reads programs keys used in PRIMARY. chemokine will follow IFN.
    out = analyze_matrix(logx, genes, programs, tfs, rng, n_perm=15)
    if out is None:
        raise RuntimeError("self-test produced no model")
    ifn_d = out["delta_boot"]["IFN"]
    mhc_d = out["delta_boot"]["MHC_I_APM"]
    print(f"self-test IFN delta {ifn_d:.4f} MHC delta {mhc_d:.4f} alpha {out['alpha']}", flush=True)
    # KO removes CLDN4, so the planted negative CLDN4->IFN link must move IFN up,
    # and the planted positive CLDN4->MHC link must move MHC down.
    if not (ifn_d > 0 and mhc_d < 0):
        raise RuntimeError(f"self-test sign failure IFN={ifn_d} MHC={mhc_d}")
    print("self-test passed", flush=True)


def load_cohort(dataset: str):
    X = sparse.load_npz(CACHE / f"{dataset}.npz").tocsr()
    meta = np.load(CACHE / f"{dataset}.meta.npz", allow_pickle=True)
    genes = [str(g) for g in meta["genes"].tolist()]
    patient = np.array([str(p) for p in meta["patient"].tolist()])
    n_count = np.asarray(meta["n_count"], dtype=np.float64)
    if X.shape[0] != patient.shape[0] or X.shape[1] != len(genes):
        raise RuntimeError(f"{dataset} shape mismatch {X.shape} {len(genes)} {patient.shape}")
    return X, genes, patient, n_count


def patient_log(X, genes, n_count, idx: np.ndarray) -> np.ndarray:
    counts = np.asarray(X[idx].toarray(), dtype=np.float64)
    lib = n_count[idx]
    lib = np.where(lib > 0, lib, 1.0)
    return np.log1p(counts / lib[:, None] * 1e4)


def main() -> None:
    self_test()
    programs = load_programs()
    tfs = set(tf_list())
    if "CLDN4" in tfs:
        raise RuntimeError("CLDN4 is in the TF list; the method note would be wrong")
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)

    inv = pd.read_csv(TAB / "unit_inventory.tsv", sep="\t")
    hist = {}
    if "histology" in inv.columns:
        for rec in inv.itertuples(index=False):
            hist[(rec.dataset, str(rec.unit_id))] = "" if pd.isna(rec.histology) else str(rec.histology)

    fit_rows = []
    delta_rows = []
    beta_acc: dict[str, list[float]] = {}
    beta_n: dict[str, int] = {}

    for dataset in COHORT_ORDER:
        path = CACHE / f"{dataset}.npz"
        if not path.exists():
            print(f"missing cache {path}", flush=True)
            continue
        X, genes, patient, n_count = load_cohort(dataset)
        print(f"{dataset} cells={X.shape[0]} genes={X.shape[1]}", flush=True)
        for unit in sorted(set(patient.tolist())):
            idx_all = np.flatnonzero(patient == unit)
            n_all = int(idx_all.size)
            reason = "ok"
            if n_all < MIN_CELLS:
                reason = f"below_min_cells_{MIN_CELLS}"
            rng = np.random.default_rng(seed_for(f"{dataset}:{unit}"))
            if n_all > MAX_CELLS:
                idx = rng.choice(idx_all, MAX_CELLS, replace=False)
            else:
                idx = idx_all
            logx = patient_log(X, genes, n_count, idx) if reason == "ok" else None
            result = None
            if reason == "ok":
                cldn_i = genes.index("CLDN4")
                pct = float(100.0 * np.mean(logx[:, cldn_i] > 0))
                if pct < MIN_CLDN4_PCT:
                    reason = f"cldn4_pct_below_{MIN_CLDN4_PCT}"
                else:
                    result = analyze_matrix(logx, genes, programs, tfs, rng, N_PERM)
                    if result is None:
                        reason = "grn_not_fit"
            if result is None:
                fit_rows.append(
                    {
                        "dataset": dataset,
                        "unit_id": unit,
                        "histology": hist.get((dataset, unit), ""),
                        "n_cells_qc": n_all,
                        "n_cells_fit": 0 if reason != "ok" else int(idx.size),
                        "status": reason,
                    }
                )
                print(f"  {dataset} {unit} n={n_all} {reason}", flush=True)
                continue
            stable = {}
            for name in PRIMARY + SECONDARY:
                a = result["delta_boot"][name]
                b = result["delta_single"][name]
                stable[name] = bool(np.isfinite(a) and np.isfinite(b) and (a == 0 or b == 0 or np.sign(a) == np.sign(b)))
            fit_rows.append(
                {
                    "dataset": dataset,
                    "unit_id": unit,
                    "histology": hist.get((dataset, unit), ""),
                    "n_cells_qc": n_all,
                    "n_cells_fit": int(idx.size),
                    "status": "fit",
                    "alpha": result["alpha"],
                    "cv_r2": result["cv_r2"],
                    "n_reg": result["n_reg"],
                    "n_targets": result["n_targets"],
                    "control_gene": result["control"],
                    "cldn4_pct": result["cldn4_pct"],
                    "cldn4_mean_log1p": result["cldn4_mean"],
                    "scale_cldn": result["scale_cldn"],
                }
            )
            for name in PRIMARY + SECONDARY:
                before = result["before"][name]
                delta = result["delta_boot"][name]
                delta_rows.append(
                    {
                        "dataset": dataset,
                        "unit_id": unit,
                        "histology": hist.get((dataset, unit), ""),
                        "program": name,
                        "tier": "primary" if name in PRIMARY else "secondary",
                        "expect": "up" if name in EXPECT_UP else "down",
                        "n_cells_fit": int(idx.size),
                        "n_genes": result["n_genes"][name],
                        "before": before,
                        "after": before + delta if np.isfinite(before) and np.isfinite(delta) else float("nan"),
                        "delta_cldn4": delta,
                        "delta_cldn4_single": result["delta_single"][name],
                        "delta_control": result["delta_control"][name],
                        "perm_p": result["perm_p"][name],
                        "stable_sign": stable[name],
                        "control_gene": result["control"],
                        "kd_like": bool(
                            np.isfinite(delta)
                            and ((name in EXPECT_UP and delta > 0) or (name in EXPECT_DOWN and delta < 0))
                        ),
                    }
                )
            for g, b in result["beta_genes"].items():
                beta_acc.setdefault(g, []).append(b)
                beta_n[g] = beta_n.get(g, 0) + 1
            print(
                f"  {dataset} {unit} n={n_all} fit={idx.size} regs={result['n_reg']} "
                f"IFN={result['delta_boot']['IFN']:.4f} MHC={result['delta_boot']['MHC_I_APM']:.4f} "
                f"chem={result['delta_boot']['chemokine']:.4f} barrier={result['delta_boot']['barrier_ligands']:.4f}",
                flush=True,
            )

    seen = {(r["dataset"], str(r["unit_id"])) for r in fit_rows}
    for rec in inv.itertuples(index=False):
        key = (rec.dataset, str(rec.unit_id))
        if key in seen:
            continue
        fit_rows.append(
            {
                "dataset": rec.dataset,
                "unit_id": str(rec.unit_id),
                "histology": "" if pd.isna(rec.histology) else str(rec.histology),
                "n_cells_qc": int(rec.n_malignant_qc) if pd.notna(rec.n_malignant_qc) else 0,
                "n_cells_fit": 0,
                "status": str(rec.status),
            }
        )

    fit_df = pd.DataFrame(fit_rows)
    delta_df = pd.DataFrame(delta_rows)
    fit_df.to_csv(TAB / "grn_fit_qc.tsv", sep="\t", index=False)
    delta_df.to_csv(TAB / "patient_program_deltas.tsv", sep="\t", index=False)

    gene_rows = []
    program_of = {}
    for name in PRIMARY + SECONDARY:
        for g in programs[name]:
            program_of.setdefault(g, []).append(name)
    for g, vals in sorted(beta_acc.items()):
        arr = np.asarray(vals, dtype=float)
        gene_rows.append(
            {
                "gene": g,
                "programs": ",".join(program_of.get(g, [])),
                "n_patients": int(arr.size),
                "median_beta_cldn4": float(np.median(arr)),
                "mean_beta_cldn4": float(np.mean(arr)),
            }
        )
    pd.DataFrame(gene_rows).to_csv(TAB / "gene_median_beta.tsv", sep="\t", index=False)

    summary_rows = []
    for name in PRIMARY + SECONDARY:
        sub = delta_df[(delta_df["program"] == name) & (delta_df["stable_sign"])].copy()
        cohorts = [("ALL", sub)] + [(ds, sub[sub["dataset"] == ds]) for ds in COHORT_ORDER]
        # ADC-only sensitivity: drop non-ADC GSE205335 rows when histology is known
        adc = sub[~((sub["dataset"] == "GSE205335") & (~sub["histology"].isin(["ADC", "LUAD", ""])))]
        # empty histology on the other cohorts stays in; GSE205335 non-ADC drops
        cohorts.append(("ALL_drop_GSE205335_nonADC", adc))
        for label, part in cohorts:
            summary_rows.append(summarize_block(name, label, part))
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(TAB / "program_summary.tsv", sep="\t", index=False)

    plot_before_after(delta_df)
    plot_deltas(delta_df)
    write_summary_json(summary, fit_df, delta_df)
    print(summary[summary["block"].isin(["ALL"] + COHORT_ORDER)].to_string(index=False), flush=True)


def summarize_block(program: str, block: str, part: pd.DataFrame) -> dict:
    d = part["delta_cldn4"].to_numpy(dtype=float)
    d = d[np.isfinite(d)]
    before = part["before"].to_numpy(dtype=float)
    after = part["after"].to_numpy(dtype=float)
    ctrl = part["delta_control"].to_numpy(dtype=float)
    expect = "up" if program in EXPECT_UP else "down"
    if d.size == 0:
        return {
            "program": program,
            "block": block,
            "expect": expect,
            "n": 0,
            "median_before": np.nan,
            "median_after": np.nan,
            "median_delta": np.nan,
            "mean_delta": np.nan,
            "n_kd_like": 0,
            "frac_kd_like": np.nan,
            "wilcoxon_p": np.nan,
            "median_perm_p": np.nan,
            "n_perm_p_lt_0.05": 0,
            "median_delta_control": np.nan,
            "wilcoxon_p_cldn4_minus_control": np.nan,
        }
    kd = (d > 0) if expect == "up" else (d < 0)
    p_w = float("nan")
    if d.size >= 8 and np.any(d != 0):
        try:
            p_w = float(stats.wilcoxon(d, alternative="two-sided", zero_method="wilcox").pvalue)
        except ValueError:
            p_w = float("nan")
    diff = part["delta_cldn4"].to_numpy(dtype=float) - part["delta_control"].to_numpy(dtype=float)
    diff = diff[np.isfinite(diff)]
    p_diff = float("nan")
    if diff.size >= 8 and np.any(diff != 0):
        try:
            p_diff = float(stats.wilcoxon(diff, alternative="two-sided", zero_method="wilcox").pvalue)
        except ValueError:
            p_diff = float("nan")
    perm = part["perm_p"].to_numpy(dtype=float)
    perm = perm[np.isfinite(perm)]
    return {
        "program": program,
        "block": block,
        "expect": expect,
        "n": int(d.size),
        "median_before": float(np.nanmedian(before)),
        "median_after": float(np.nanmedian(after)),
        "median_delta": float(np.median(d)),
        "mean_delta": float(np.mean(d)),
        "n_kd_like": int(kd.sum()),
        "frac_kd_like": float(kd.mean()),
        "wilcoxon_p": p_w,
        "median_perm_p": float(np.median(perm)) if perm.size else float("nan"),
        "n_perm_p_lt_0.05": int(np.sum(perm <= 0.05)) if perm.size else 0,
        "median_delta_control": float(np.nanmedian(ctrl)),
        "wilcoxon_p_cldn4_minus_control": p_diff,
    }


def plot_before_after(delta_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.6), sharey=False)
    for ax, name in zip(axes, PRIMARY):
        sub = delta_df[(delta_df["program"] == name) & (delta_df["stable_sign"])]
        for ds, part in sub.groupby("dataset"):
            color = COHORT_COLOR.get(ds, "gray")
            for rec in part.itertuples(index=False):
                ax.plot([0, 1], [rec.before, rec.after], color=color, alpha=0.35, lw=0.8)
            ax.scatter(np.zeros(len(part)), part["before"], s=12, color=color, label=ds, zorder=3)
            ax.scatter(np.ones(len(part)), part["after"], s=12, color=color, zorder=3)
        ax.set_xticks([0, 1], ["before", "CLDN4 virtual KO"])
        ax.set_title(name, fontsize=10)
        ax.set_xlim(-0.25, 1.25)
    axes[0].set_ylabel("mean log1p CP10k")
    handles, labels = axes[0].get_legend_handles_labels()
    # one handle per cohort
    uniq = dict(zip(labels, handles))
    fig.legend(uniq.values(), uniq.keys(), loc="upper center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle("Malignant-cell program score, patient / sample unit", y=1.05, fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "before_after_programs.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "before_after_programs.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_deltas(delta_df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(12.5, 3.8), sharey=False)
    rng = np.random.default_rng(1)
    for ax, name in zip(axes, PRIMARY):
        sub = delta_df[(delta_df["program"] == name) & (delta_df["stable_sign"])]
        for i, ds in enumerate(COHORT_ORDER):
            part = sub[sub["dataset"] == ds]["delta_cldn4"].to_numpy(dtype=float)
            part = part[np.isfinite(part)]
            if part.size == 0:
                continue
            jitter = rng.normal(0, 0.06, size=part.size)
            ax.scatter(np.full(part.size, i) + jitter, part, s=16, color=COHORT_COLOR[ds], alpha=0.85)
            ax.hlines(np.median(part), i - 0.25, i + 0.25, color="black", lw=1.4)
        ax.axhline(0, color="gray", lw=0.6)
        ax.set_xticks(range(len(COHORT_ORDER)), ["123902", "131907", "205335", "189357"], rotation=0)
        expect = "expect UP" if name in EXPECT_UP else "expect DOWN"
        ax.set_title(f"{name}\n{expect}", fontsize=9)
    axes[0].set_ylabel("virtual-KO delta (log1p CP10k)")
    fig.tight_layout()
    fig.savefig(FIG / "delta_by_cohort.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIG / "delta_by_cohort.pdf", bbox_inches="tight")
    plt.close(fig)


def write_summary_json(summary: pd.DataFrame, fit_df: pd.DataFrame, delta_df: pd.DataFrame) -> None:
    payload = {
        "method": "celloracle_equivalent_one_step_ridge",
        "celloracle_package_executed": False,
        "reason_not_stock_celloracle": (
            "CLDN4 is absent from human_tfs_from_priors.txt. "
            "A motif base GRN assigns CLDN4 no targets, so the stock KO delta is zero."
        ),
        "min_cells": MIN_CELLS,
        "max_cells": MAX_CELLS,
        "min_cldn4_pct": MIN_CLDN4_PCT,
        "n_boot": N_BOOT,
        "n_perm": N_PERM,
        "n_units_inventory": int(fit_df.shape[0]) if len(fit_df) else 0,
        "n_units_fit": int((fit_df["status"] == "fit").sum()) if len(fit_df) else 0,
        "fit_status_counts": fit_df["status"].value_counts().to_dict() if len(fit_df) else {},
        "primary_all": summary[(summary["block"] == "ALL") & (summary["program"].isin(PRIMARY))]
        .to_dict(orient="records"),
    }
    (TAB / "summary.json").write_text(json.dumps(payload, indent=2, default=float))
    print(json.dumps(payload["fit_status_counts"], indent=2), flush=True)


if __name__ == "__main__":
    main()
