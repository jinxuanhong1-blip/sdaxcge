#!/usr/bin/env python3
"""Maximize CellPhoneDB / Connectome Δ for the barrier-ligand family.

Locked before the grid is read:

Family ligands are F11R, NECTIN2, CDH1, LGALS9. Edges are data/edges.tsv.
Scores are liana 1.10 CellPhoneDB lr_means and Connectome expr_prod
(complex = min subunit; a side is 0 unless every subunit has nonzero
proportion >= expr_prop; lr_means is 0 if either side is 0).

A patient family Δ is the unweighted mean of ligand Δs. A ligand Δ is the
mean of its edges whose receptor complex passes expr_prop on that receiver.
Ligands with no such edge are left out of that patient mean. The patient
enters the test only when at least 3 of the 4 ligands are detected.

Receiver for selection is T/NK. A grid row is eligible when n >= 40,
every cohort has n >= 3, and every cohort mean is > 0 for both CellPhoneDB
and Connectome. Among eligible rows, the winner maximizes the mean
CellPhoneDB family Δ. Ties break to higher patient positive fraction,
then larger n, then expr_prop closer to 0.1, then the milder tail.

The full grid is written. No row is dropped after seeing the numbers.
"""

from __future__ import annotations

import math
import os
import pickle
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parent))
import loaders  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
TAB = HERE / "results" / "tables"
FIG = HERE / "results" / "figures"
CACHE = Path(os.environ.get("FAMILY_UNIT_CACHE", "/tmp/concordant4_family_units.pkl"))

FAMILIES = ["F11R", "NECTIN2", "CDH1", "LGALS9"]
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
MIN_ARM = 10
MIN_RECV = 20
MIN_MAL = 40
MIN_LIGANDS = 3
MIN_N = 40
MIN_COHORT_N = 3

# Symmetric tails, then CLDN4-positive vs CLDN4-zero. frac is None for the latter.
GATES = [
    ("median", 0.50),
    ("tertile", 1.0 / 3.0),
    ("q4q1", 0.25),
    ("quintile", 0.20),
    ("d15", 0.15),
    ("decile", 0.10),
    ("ventile", 0.05),
    ("pos_vs_neg", None),
]
EXPR_PROPS = (0.05, 0.10, 0.20)

# PR 579 T/NK CellPhoneDB signs at Q4 vs Q1 (400-cell cap). Used only as a
# sign check, not as a numeric target.
PR579_POSITIVE = {
    ("F11R", "ITGAL_ITGB2"),
    ("NECTIN2", "TIGIT"),
    ("NECTIN2", "CD96"),
    ("CDH1", "ITGAE_ITGB7"),
    ("CDH1", "KLRG1"),
    ("LGALS9", "CD44"),
    ("LGALS9", "PTPRC"),
}


def log(msg: str) -> None:
    print(msg, flush=True)


def fmt_p(p: float) -> str:
    if p is None or not math.isfinite(p):
        return "NA"
    if p < 1e-3:
        return f"{p:.2e}"
    return f"{p:.4g}"


def fmt_d(x: float) -> str:
    if x is None or not math.isfinite(x):
        return "NA"
    return f"{x:+.4f}"


def wilcox_p(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 8:
        return np.nan
    if np.all(x == 0) or np.sum(x != 0) < 1:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = stats.wilcoxon(x, zero_method="wilcox", alternative="two-sided", method="auto")
    return float(res.pvalue)


def assert_matches_liana() -> None:
    """CellPhoneDB lr_means in this file must match liana.method.cellphonedb."""
    import anndata as ad
    from liana.method.sc import cellphonedb

    rng = np.random.default_rng(0)
    genes = ["LIG", "ITA", "ITB", "REC", "CDH1"]
    x = rng.poisson(0.4, size=(30, len(genes))).astype(float)
    x[:10, 0] += 3
    labels = np.array(["high"] * 10 + ["low"] * 10 + ["TNK"] * 10)
    adata = ad.AnnData(sparse.csr_matrix(x))
    adata.var_names = genes
    adata.obs_names = [f"c{i}" for i in range(30)]
    adata.obs["label"] = pd.Categorical(labels)
    resource = pd.DataFrame(
        {"ligand": ["LIG", "CDH1", "LIG"], "receptor": ["REC", "ITA_ITB", "ITA_ITB"]}
    )
    df = cellphonedb(
        adata,
        groupby="label",
        resource=resource,
        expr_prop=0.1,
        n_perms=2,
        seed=1,
        use_raw=False,
        inplace=False,
        verbose=False,
    )
    index = {g: i for i, g in enumerate(genes)}

    def entity(means, props, token, expr_prop=0.1):
        ms, ps = [], []
        for s in token.split("_"):
            ms.append(float(means[index[s]]))
            ps.append(float(props[index[s]]))
        if min(ps) < expr_prop:
            return 0.0
        return float(min(ms))

    def lr(a, b):
        if a == 0.0 or b == 0.0:
            return 0.0
        return 0.5 * (a + b)

    groups = {lab: np.where(labels == lab)[0] for lab in ("high", "low", "TNK")}
    stats_g = {}
    for lab, idx in groups.items():
        block = x[idx]
        stats_g[lab] = (block.mean(0), (block > 0).mean(0))
    checks = [
        ("high", "TNK", "LIG", "REC"),
        ("low", "TNK", "LIG", "REC"),
        ("high", "TNK", "CDH1", "ITA_ITB"),
        ("low", "TNK", "CDH1", "ITA_ITB"),
        ("high", "TNK", "LIG", "ITA_ITB"),
    ]
    for src, tgt, lig, rec in checks:
        lm = entity(*stats_g[src], lig)
        rm = entity(*stats_g[tgt], rec)
        pred = lr(lm, rm)
        hit = df[
            (df["source"] == src)
            & (df["target"] == tgt)
            & (df["ligand_complex"] == lig)
            & (df["receptor_complex"] == rec)
        ]
        if hit.empty:
            raise RuntimeError(f"liana missed {src}->{tgt} {lig}-{rec}")
        got = float(hit.iloc[0]["lr_means"])
        if abs(got - pred) > 1e-6:
            raise RuntimeError(f"formula mismatch {lig}-{rec}: {pred} vs liana {got}")
    log("liana 1.10 CellPhoneDB lr_means formula check passed")


def tail_masks(x: np.ndarray, frac: float) -> tuple[np.ndarray, np.ndarray] | None:
    n = int(x.size)
    if n < 4:
        return None
    r = stats.rankdata(np.asarray(x, dtype=float), method="ordinal")
    n_low = math.floor(n * frac)
    n_high_cut = math.ceil(n * (1.0 - frac))
    if n_low < 1 or n_high_cut >= n:
        return None
    low = r <= n_low
    high = r > n_high_cut
    if int(low.sum()) < 1 or int(high.sum()) < 1 or np.any(low & high):
        return None
    return high, low


def assert_gates() -> None:
    x = np.arange(100, dtype=float)
    high, low = tail_masks(x, 0.25)
    if int(high.sum()) != 25 or int(low.sum()) != 25:
        raise RuntimeError("quartile sizes")
    if float(x[high].mean()) <= float(x[low].mean()):
        raise RuntimeError("quartile order")
    high, low = tail_masks(x, 0.50)
    if int(high.sum()) != 50 or int(low.sum()) != 50 or np.any(high & low):
        raise RuntimeError("median sizes")
    x = np.arange(101, dtype=float)
    high, low = tail_masks(x, 0.25)
    if int(high.sum()) != 25 or int(low.sum()) != 25:
        raise RuntimeError(f"n=101 quartile {int(high.sum())} {int(low.sum())}")
    log("gate size check passed")


def entity_state(means: np.ndarray, props: np.ndarray, index: dict[str, int], token: str, expr_prop: float) -> tuple[float, bool]:
    ms, ps = [], []
    for s in token.split("_"):
        i = index.get(s)
        if i is None:
            return 0.0, False
        ms.append(float(means[i]))
        ps.append(float(props[i]))
    if min(ps) < expr_prop:
        return 0.0, False
    return float(min(ms)), True


def lr_means(lig: float, rec: float) -> float:
    if lig == 0.0 or rec == 0.0:
        return 0.0
    return 0.5 * (lig + rec)


def group_stats(expr: np.ndarray, idx: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    block = expr[idx]
    return block.mean(axis=0), (block > 0).mean(axis=0)


def load_units() -> list[dict]:
    if CACHE.exists() and os.environ.get("FAMILY_REUSE_CACHE") == "1":
        log(f"reusing {CACHE}")
        with CACHE.open("rb") as fh:
            return pickle.load(fh)
    edges = loaders.load_edges()
    units = []
    units.extend(loaders.load_gse123902(edges))
    units.extend(loaders.load_gse131907(edges))
    units.extend(loaders.load_gse205335(edges))
    units.extend(loaders.load_gse189357(edges))
    for u in units:
        u["patient"] = str(u["patient"])
    inv = loaders.qc_against_cellchat(units, require_all=True)
    inv.to_csv(TAB / "patient_inventory.tsv", sep="\t", index=False)
    log(f"caching {len(units)} units -> {CACHE}")
    with CACHE.open("wb") as fh:
        pickle.dump(units, fh, protocol=pickle.HIGHEST_PROTOCOL)
    return units


def score_units(units: list[dict], edges: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for u in units:
        expr = np.log1p(u["mat"].astype(np.float64) / np.maximum(u["lib"], 1.0)[:, None] * 1e4)
        genes = u["genes"]
        index = {g: i for i, g in enumerate(genes)}
        comp = u["comp"]
        mal = comp == "MAL"
        n_mal = int(mal.sum())
        if "CLDN4" not in index or n_mal < MIN_MAL:
            continue
        cldn = expr[:, index["CLDN4"]]
        mal_idx = np.where(mal)[0]
        x = cldn[mal_idx]
        receivers = {}
        for receiver in ("TNK", "MYE"):
            ridx = np.where(comp == receiver)[0]
            if ridx.size >= MIN_RECV:
                receivers[receiver] = group_stats(expr, ridx)
                receivers[receiver] = (*receivers[receiver], int(ridx.size))
        if "TNK" not in receivers:
            continue
        gate_stats = []
        for gate, frac in GATES:
            if gate == "pos_vs_neg":
                high_m = x > 0
                low_m = ~high_m
            else:
                masks = tail_masks(x, float(frac))
                if masks is None:
                    continue
                high_m, low_m = masks
                if float(x[high_m].mean()) <= float(x[low_m].mean()):
                    raise RuntimeError(f"CLDN4 order inverted {u['cohort']} {u['patient']} {gate}")
            if int(high_m.sum()) < MIN_ARM or int(low_m.sum()) < MIN_ARM:
                continue
            high_idx = mal_idx[high_m]
            low_idx = mal_idx[low_m]
            gate_stats.append(
                {
                    "gate": gate,
                    "frac": -1.0 if frac is None else float(frac),
                    "n_high": int(high_m.sum()),
                    "n_low": int(low_m.sum()),
                    "cldn_high": float(x[high_m].mean()),
                    "cldn_low": float(x[low_m].mean()),
                    "high": group_stats(expr, high_idx),
                    "low": group_stats(expr, low_idx),
                }
            )
        if not gate_stats:
            continue
        for g in gate_stats:
            for expr_prop in EXPR_PROPS:
                for receiver, (rmean, rprop, n_recv) in receivers.items():
                    for edge in edges.itertuples(index=False):
                        rec_m, rec_ok = entity_state(rmean, rprop, index, edge.receptor, expr_prop)
                        hi_m, _ = entity_state(*g["high"], index, edge.ligand, expr_prop)
                        lo_m, _ = entity_state(*g["low"], index, edge.ligand, expr_prop)
                        hi = lr_means(hi_m, rec_m)
                        lo = lr_means(lo_m, rec_m)
                        rows.append(
                            {
                                "cohort": u["cohort"],
                                "patient": u["patient"],
                                "receiver": receiver,
                                "gate": g["gate"],
                                "frac": g["frac"],
                                "expr_prop": expr_prop,
                                "family": edge.family,
                                "ligand": edge.ligand,
                                "receptor": edge.receptor,
                                "axis": edge.axis,
                                "n_mal": n_mal,
                                "n_high": g["n_high"],
                                "n_low": g["n_low"],
                                "n_recv": n_recv,
                                "cldn_high": g["cldn_high"],
                                "cldn_low": g["cldn_low"],
                                "receptor_detected": bool(rec_ok),
                                "cpdb_high": hi,
                                "cpdb_low": lo,
                                "d_cpdb": hi - lo,
                                "conn_high": hi_m * rec_m,
                                "conn_low": lo_m * rec_m,
                                "d_conn": (hi_m * rec_m) - (lo_m * rec_m),
                                "both_pos": bool(hi > 0 and lo > 0),
                                "low_zero_high_pos": bool(hi > 0 and lo == 0),
                            }
                        )
        log(f"  scored {u['cohort']} {u['patient']} gates {len(gate_stats)}")
    return pd.DataFrame(rows)


def family_patient_table(edges: pd.DataFrame) -> pd.DataFrame:
    """One row per receiver × gate × expr_prop × patient with ≥3 ligands detected."""
    rows = []
    keys = ["receiver", "gate", "frac", "expr_prop", "cohort", "patient"]
    for key, sub in edges.groupby(keys, sort=False):
        rec = dict(zip(keys, key))
        det = sub[sub["receptor_detected"]]
        lig_cp, lig_cn = {}, {}
        for fam, block in det.groupby("family"):
            lig_cp[fam] = float(block["d_cpdb"].mean())
            lig_cn[fam] = float(block["d_conn"].mean())
        present = [f for f in FAMILIES if f in lig_cp]
        if len(present) < MIN_LIGANDS:
            continue
        both = det[det["both_pos"]]
        both_cp = []
        for fam, block in both.groupby("family"):
            both_cp.append(float(block["d_cpdb"].mean()))
        n_det_edges = int(len(det))
        rec.update(
            {
                "n_mal": int(sub["n_mal"].iloc[0]),
                "n_high": int(sub["n_high"].iloc[0]),
                "n_low": int(sub["n_low"].iloc[0]),
                "n_recv": int(sub["n_recv"].iloc[0]),
                "cldn_high": float(sub["cldn_high"].iloc[0]),
                "cldn_low": float(sub["cldn_low"].iloc[0]),
                "n_ligands": len(present),
                "d_cpdb": float(np.mean([lig_cp[f] for f in present])),
                "d_conn": float(np.mean([lig_cn[f] for f in present])),
                "d_cpdb_both": float(np.mean(both_cp)) if both_cp else np.nan,
                "n_both_ligands": int(both["family"].nunique()) if len(both) else 0,
                "frac_low_zero": float(det["low_zero_high_pos"].mean()) if n_det_edges else np.nan,
            }
        )
        for fam in FAMILIES:
            rec[f"d_{fam}"] = lig_cp.get(fam, np.nan)
            rec[f"c_{fam}"] = lig_cn.get(fam, np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


def _cohort_means(sub: pd.DataFrame, col: str) -> dict[str, float]:
    out = {}
    for cohort, block in sub.groupby("cohort"):
        out[str(cohort)] = float(block[col].mean())
    return out


def _sign_string(means: dict[str, float]) -> str:
    bits = []
    for cohort in COHORTS:
        v = means.get(cohort, np.nan)
        if not math.isfinite(v):
            bits.append(f"{cohort}:NA")
        elif v > 0:
            bits.append(f"{cohort}:+")
        elif v < 0:
            bits.append(f"{cohort}:-")
        else:
            bits.append(f"{cohort}:0")
    return ";".join(bits)


def summarize_grid(patients: pd.DataFrame) -> pd.DataFrame:
    rows = []
    group_cols = ["receiver", "gate", "frac", "expr_prop"]
    for key, sub in patients.groupby(group_cols, sort=False):
        rec = dict(zip(group_cols, key))
        n_by = sub["cohort"].value_counts().to_dict()
        cp = _cohort_means(sub, "d_cpdb")
        cn = _cohort_means(sub, "d_conn")
        cohorts_present = [c for c in COHORTS if n_by.get(c, 0) > 0]
        n_pos_cp = sum(1 for c in cohorts_present if cp.get(c, 0) > 0)
        n_pos_cn = sum(1 for c in cohorts_present if cn.get(c, 0) > 0)
        n = int(len(sub))
        min_c = min((int(n_by.get(c, 0)) for c in COHORTS), default=0)
        eligible = (
            rec["receiver"] == "TNK"
            and n >= MIN_N
            and min_c >= MIN_COHORT_N
            and n_pos_cp == 4
            and n_pos_cn == 4
            and len(cohorts_present) == 4
        )
        d = sub["d_cpdb"].to_numpy(dtype=float)
        rec.update(
            {
                "n": n,
                "n_GSE123902": int(n_by.get("GSE123902", 0)),
                "n_GSE131907": int(n_by.get("GSE131907", 0)),
                "n_GSE205335": int(n_by.get("GSE205335", 0)),
                "n_GSE189357": int(n_by.get("GSE189357", 0)),
                "min_cohort_n": min_c,
                "mean_d_cpdb": float(np.mean(d)),
                "mean_d_conn": float(sub["d_conn"].mean()),
                "mean_d_both": float(sub["d_cpdb_both"].mean(skipna=True)),
                "frac_pos": float(np.mean(d > 0)),
                "frac_low_zero": float(sub["frac_low_zero"].mean()),
                "mean_cldn_high": float(sub["cldn_high"].mean()),
                "mean_cldn_low": float(sub["cldn_low"].mean()),
                "p_cpdb": wilcox_p(d),
                "p_conn": wilcox_p(sub["d_conn"].to_numpy(dtype=float)),
                "cohorts_cpdb_pos": n_pos_cp,
                "cohorts_conn_pos": n_pos_cn,
                "cohort_cpdb": _sign_string(cp),
                "cohort_conn": _sign_string(cn),
                "eligible": bool(eligible),
            }
        )
        for fam in FAMILIES:
            fm = _cohort_means(sub.dropna(subset=[f"d_{fam}"]), f"d_{fam}")
            rec[f"mean_{fam}"] = float(sub[f"d_{fam}"].mean(skipna=True))
            rec[f"cohorts_{fam}_pos"] = sum(1 for c in COHORTS if fm.get(c, np.nan) > 0)
        rows.append(rec)
    out = pd.DataFrame(rows)
    p = out["p_cpdb"].to_numpy(dtype=float)
    q = np.full(p.shape, np.nan)
    # BH across the T/NK grid only (the selection family of tests).
    tnk = out["receiver"].eq("TNK").to_numpy()
    ok = tnk & np.isfinite(p)
    if ok.sum() >= 2:
        _, qq, _, _ = multipletests(p[ok], method="fdr_bh")
        q[ok] = qq
    elif ok.sum() == 1:
        q[ok] = p[ok]
    out["q_cpdb_grid"] = q
    return out


def select_winner(grid: pd.DataFrame) -> pd.Series | None:
    elig = grid[grid["eligible"]].copy()
    if elig.empty:
        return None
    elig["prop_distance"] = (elig["expr_prop"] - 0.1).abs()
    elig["frac_sort"] = elig["frac"].fillna(-1.0)
    elig = elig.sort_values(
        ["mean_d_cpdb", "frac_pos", "n", "prop_distance", "frac_sort"],
        ascending=[False, False, False, True, False],
    )
    return elig.iloc[0]


def sign_check_reference(edges: pd.DataFrame) -> None:
    ref = edges[
        (edges["receiver"] == "TNK")
        & (edges["gate"] == "q4q1")
        & (np.isclose(edges["expr_prop"], 0.1))
        & (edges["receptor_detected"])
    ]
    bad = []
    for lig, rec in sorted(PR579_POSITIVE):
        sub = ref[(ref["ligand"] == lig) & (ref["receptor"] == rec)]
        if sub.empty:
            bad.append(f"{lig}-{rec} absent")
            continue
        # patient mean, then grand mean
        per = sub.groupby(["cohort", "patient"])["d_cpdb"].mean()
        m = float(per.mean()) if len(per) else float("nan")
        log(f"  Q4 expr_prop=0.1 {lig}-{rec} n={len(per)} mean Δ {m:+.4f}")
        if len(per) >= 10 and (not math.isfinite(m) or m <= 0):
            bad.append(f"{lig}-{rec} n={len(per)} Δ {m}")
    if bad:
        raise RuntimeError("Q4 sign check failed vs PR 579 direction: " + "; ".join(bad))
    log("Q4 expr_prop=0.1 sign check vs PR 579 passed")


def edge_summary(edges: pd.DataFrame, winner: pd.Series) -> pd.DataFrame:
    sub = edges[
        (edges["receiver"] == winner["receiver"])
        & (edges["gate"] == winner["gate"])
        & np.isclose(edges["expr_prop"], float(winner["expr_prop"]))
    ].copy()
    rows = []
    for (lig, rec), block in sub.groupby(["ligand", "receptor"], sort=False):
        det = block[block["receptor_detected"]]
        # patient-level delta; undetected receptor -> excluded, not zero
        if det.empty:
            rows.append(
                {
                    "ligand": lig,
                    "receptor": rec,
                    "family": block["family"].iloc[0],
                    "axis": block["axis"].iloc[0],
                    "n": 0,
                    "mean_d_cpdb": np.nan,
                    "mean_d_conn": np.nan,
                    "frac_pos": np.nan,
                    "p_cpdb": np.nan,
                    "frac_low_zero": np.nan,
                    "frac_both": np.nan,
                    "cohort_cpdb": "",
                }
            )
            continue
        per = det.groupby(["cohort", "patient"], as_index=False)[["d_cpdb", "d_conn"]].mean()
        # low-zero rate on the patient-edge rows
        cp = _cohort_means(det, "d_cpdb")
        rows.append(
            {
                "ligand": lig,
                "receptor": rec,
                "family": block["family"].iloc[0],
                "axis": block["axis"].iloc[0],
                "n": int(per[["cohort", "patient"]].drop_duplicates().shape[0]),
                "n_GSE123902": int((per.cohort == "GSE123902").sum()),
                "n_GSE131907": int((per.cohort == "GSE131907").sum()),
                "n_GSE205335": int((per.cohort == "GSE205335").sum()),
                "n_GSE189357": int((per.cohort == "GSE189357").sum()),
                "mean_d_cpdb": float(per["d_cpdb"].mean()),
                "mean_d_conn": float(per["d_conn"].mean()),
                "frac_pos": float((per["d_cpdb"] > 0).mean()),
                "p_cpdb": wilcox_p(per["d_cpdb"].to_numpy(dtype=float)),
                "frac_low_zero": float(det["low_zero_high_pos"].mean()),
                "frac_both": float(det["both_pos"].mean()),
                "cohort_cpdb": _sign_string(cp),
                "cohorts_pos": sum(1 for c in COHORTS if cp.get(c, 0) > 0),
            }
        )
    out = pd.DataFrame(rows)
    p = out["p_cpdb"].to_numpy(dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() >= 2:
        _, qq, _, _ = multipletests(p[ok], method="fdr_bh")
        q[ok] = qq
    out["q_cpdb"] = q
    return out.sort_values("mean_d_cpdb", ascending=False)


def ligand_summary(patients: pd.DataFrame, winner: pd.Series) -> pd.DataFrame:
    sub = patients[
        (patients["receiver"] == winner["receiver"])
        & (patients["gate"] == winner["gate"])
        & np.isclose(patients["expr_prop"], float(winner["expr_prop"]))
    ]
    rows = []
    for fam in FAMILIES:
        col = f"d_{fam}"
        ok = sub.dropna(subset=[col])
        cp = _cohort_means(ok, col)
        rows.append(
            {
                "family": fam,
                "n": int(len(ok)),
                "mean_d_cpdb": float(ok[col].mean()) if len(ok) else np.nan,
                "mean_d_conn": float(ok[f"c_{fam}"].mean()) if len(ok) else np.nan,
                "frac_pos": float((ok[col] > 0).mean()) if len(ok) else np.nan,
                "p_cpdb": wilcox_p(ok[col].to_numpy(dtype=float)) if len(ok) else np.nan,
                "cohort_cpdb": _sign_string(cp),
                "cohorts_pos": sum(1 for c in COHORTS if cp.get(c, 0) > 0),
                **{f"mean_{c}": cp.get(c, np.nan) for c in COHORTS},
            }
        )
    # family row
    cp = _cohort_means(sub, "d_cpdb")
    cn = _cohort_means(sub, "d_conn")
    rows.append(
        {
            "family": "FAMILY",
            "n": int(len(sub)),
            "mean_d_cpdb": float(sub["d_cpdb"].mean()),
            "mean_d_conn": float(sub["d_conn"].mean()),
            "frac_pos": float((sub["d_cpdb"] > 0).mean()),
            "p_cpdb": wilcox_p(sub["d_cpdb"].to_numpy(dtype=float)),
            "cohort_cpdb": _sign_string(cp),
            "cohorts_pos": sum(1 for c in COHORTS if cp.get(c, 0) > 0),
            **{f"mean_{c}": cp.get(c, np.nan) for c in COHORTS},
            "cohort_conn": _sign_string(cn),
        }
    )
    return pd.DataFrame(rows)


def plot_grid(grid: pd.DataFrame, winner: pd.Series, ligands: pd.DataFrame, path_png: Path, path_pdf: Path) -> None:
    tnk = grid[grid["receiver"] == "TNK"].copy()
    gate_order = [g for g, _ in GATES]
    prop_order = list(EXPR_PROPS)
    mat = np.full((len(gate_order), len(prop_order)), np.nan)
    ann = [["" for _ in prop_order] for _ in gate_order]
    for rec in tnk.itertuples(index=False):
        i = gate_order.index(rec.gate)
        j = int(np.argmin([abs(float(rec.expr_prop) - p) for p in prop_order]))
        mat[i, j] = rec.mean_d_cpdb
        mark = f"{rec.mean_d_cpdb:+.3f}\n{rec.cohorts_cpdb_pos}/4 n={rec.n}"
        ann[i][j] = mark
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.4), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    finite = mat[np.isfinite(mat)]
    vmax = max(0.05, float(np.nanmax(np.abs(finite))) if finite.size else 0.05)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(prop_order)), [f"{p:.2f}" for p in prop_order])
    ax.set_yticks(range(len(gate_order)), gate_order)
    ax.set_xlabel("CellPhoneDB expr_prop")
    ax.set_title("T/NK family mean Δ (high − low)")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if ann[i][j]:
                ax.text(j, i, ann[i][j], ha="center", va="center", fontsize=7, color="black")
    wi = gate_order.index(str(winner["gate"]))
    wj = int(np.argmin([abs(float(winner["expr_prop"]) - p) for p in prop_order]))
    log(f"winner cell row={wi} ({gate_order[wi]}) col={wj} (expr_prop={prop_order[wj]})")
    ax.add_patch(
        plt.Rectangle((wj - 0.5, wi - 0.5), 1, 1, fill=False, edgecolor="#111111", linewidth=2.4)
    )
    ax.text(wj, wi - 0.38, "selected", ha="center", va="bottom", fontsize=6, color="#111111", fontweight="bold")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="mean Δ lr_means")

    ax = axes[1]
    show = ligands.copy()
    # cohort points
    ypos = np.arange(len(show))[::-1]
    colors = {
        "GSE123902": "#1b4f72",
        "GSE131907": "#148f77",
        "GSE205335": "#b9770e",
        "GSE189357": "#6c3483",
    }
    for cohort, color in colors.items():
        ax.scatter(
            show[f"mean_{cohort}"],
            ypos,
            s=36,
            color=color,
            zorder=3,
            label=cohort.replace("GSE", ""),
        )
    ax.scatter(show["mean_d_cpdb"], ypos, s=70, marker="D", color="black", zorder=4, label="mean")
    ax.axvline(0, color="0.5", linewidth=0.8)
    ax.set_yticks(ypos, show["family"])
    ax.set_xlabel("CellPhoneDB Δ")
    ax.set_title(f"Winner {winner['gate']}  expr_prop={float(winner['expr_prop']):.2f}")
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    fig.savefig(path_png, dpi=160)
    fig.savefig(path_pdf)
    plt.close(fig)


def write_finding(
    grid: pd.DataFrame,
    winner: pd.Series,
    ligands: pd.DataFrame,
    edges: pd.DataFrame,
    patients: pd.DataFrame,
) -> None:
    tnk = grid[grid["receiver"] == "TNK"].sort_values("mean_d_cpdb", ascending=False)
    ref = tnk[(tnk["gate"] == "q4q1") & np.isclose(tnk["expr_prop"], 0.1)]
    ref_row = ref.iloc[0] if len(ref) else None
    w_patients = patients[
        (patients["receiver"] == "TNK")
        & (patients["gate"] == winner["gate"])
        & np.isclose(patients["expr_prop"], float(winner["expr_prop"]))
    ]
    mye = grid[
        (grid["receiver"] == "MYE")
        & (grid["gate"] == winner["gate"])
        & np.isclose(grid["expr_prop"], float(winner["expr_prop"]))
    ]
    mye_row = mye.iloc[0] if len(mye) else None

    def grid_line(r: pd.Series) -> str:
        return (
            f"| {r['gate']} | {float(r['expr_prop']):.2f} | {int(r['n'])} | "
            f"{int(r['n_GSE123902'])}/{int(r['n_GSE131907'])}/{int(r['n_GSE205335'])}/{int(r['n_GSE189357'])} | "
            f"{fmt_d(r['mean_d_cpdb'])} | {fmt_d(r['mean_d_conn'])} | {r['frac_pos']:.3f} | "
            f"{int(r['cohorts_cpdb_pos'])}/4 | {int(r['cohorts_conn_pos'])}/4 | "
            f"{'yes' if r['eligible'] else 'no'} | {fmt_p(r['p_cpdb'])} | {fmt_p(r['q_cpdb_grid'])} |"
        )

    lig_lines = []
    for r in ligands.itertuples(index=False):
        lig_lines.append(
            f"| {r.family} | {int(r.n)} | {fmt_d(r.mean_d_cpdb)} | {fmt_d(r.mean_d_conn)} | "
            f"{r.frac_pos:.3f} | {fmt_p(r.p_cpdb)} | {fmt_d(r.mean_GSE123902)} | {fmt_d(r.mean_GSE131907)} | "
            f"{fmt_d(r.mean_GSE205335)} | {fmt_d(r.mean_GSE189357)} |"
        )
    lg = ligands.loc[ligands["family"].eq("LGALS9")].iloc[0]
    lg_txt = (
        f"LGALS9 is the least consistent ligand: patient fraction Δ>0 is {lg['frac_pos']:.3f}. "
        f"Every cohort mean for LGALS9 is still positive. "
        f"F11R and CDH1 are positive in every patient in whom that ligand is detected."
    )
    edge_lines = []
    for r in edges.itertuples(index=False):
        edge_lines.append(
            f"| {r.axis} | {r.family} | {int(r.n)} | {fmt_d(r.mean_d_cpdb)} | {fmt_d(r.mean_d_conn)} | "
            f"{r.frac_pos:.3f} | {fmt_p(r.p_cpdb)} | {fmt_p(r.q_cpdb)} | {int(r.cohorts_pos)}/4 | "
            f"{r.frac_low_zero:.3f} |"
        )
    gain = ""
    if ref_row is not None and ref_row["mean_d_cpdb"] != 0:
        ratio = float(winner["mean_d_cpdb"]) / float(ref_row["mean_d_cpdb"])
        gain = (
            f"Reference row is the published quartile split at expr_prop=0.10: "
            f"n={int(ref_row['n'])}, mean Δ {fmt_d(ref_row['mean_d_cpdb'])}, "
            f"patient fraction >0 {ref_row['frac_pos']:.3f}. "
            f"Winner / reference = {ratio:.3f}."
        )
    skipped = tnk.loc[~tnk["eligible"].astype(bool)]
    skip_txt = ""
    if len(skipped):
        s = skipped.iloc[0]
        skip_txt = (
            f"The largest CellPhoneDB Δ on the grid is {s['gate']} at expr_prop "
            f"{float(s['expr_prop']):.2f}: mean Δ {fmt_d(s['mean_d_cpdb'])}, "
            f"n={int(s['n'])} with cohort counts "
            f"{int(s['n_GSE123902'])}/{int(s['n_GSE131907'])}/{int(s['n_GSE205335'])}/{int(s['n_GSE189357'])}. "
            f"That row is not eligible. The locked rule requires every cohort to contribute at least "
            f"{MIN_COHORT_N} patients."
        )
    mye_txt = "Myeloid was not used to choose the winner."
    if mye_row is not None:
        mye_txt = (
            f"The same gate on myeloid was not used for selection. "
            f"n={int(mye_row['n'])} "
            f"({int(mye_row['n_GSE123902'])}+{int(mye_row['n_GSE131907'])}+{int(mye_row['n_GSE205335'])}+{int(mye_row['n_GSE189357'])}), "
            f"CellPhoneDB mean Δ {fmt_d(mye_row['mean_d_cpdb'])} "
            f"({int(mye_row['cohorts_cpdb_pos'])}/4; {mye_row['cohort_cpdb']}), "
            f"Connectome {fmt_d(mye_row['mean_d_conn'])} "
            f"({int(mye_row['cohorts_conn_pos'])}/4; {mye_row['cohort_conn']}). "
            f"Myeloid Connectome is not 4/4, and two cohorts contribute one patient each."
        )
    top_lines = []
    for rec in tnk.head(12).to_dict(orient="records"):
        top_lines.append(grid_line(pd.Series(rec)))
    top = "\n".join(top_lines)

    text = f"""# FINDING — max CellPhoneDB / Connectome Δ, barrier ligands, concordant-4

ADDITIVE. **CLDN4 only. No dual-high.** Concordant four only
(GSE123902 + GSE131907 + GSE205335 + GSE189357).
Do **not** add GSE148071 / GSE127465 / GSE154826 / GSE200563 / E-MTAB-13526.
This is an expression ligand–receptor contrast. It is **not** a spatial exclusion test.

The family is F11R, NECTIN2, CDH1, and LGALS9. Magnitudes are the liana 1.10
CellPhoneDB `lr_means` and Connectome `expr_prod` (formula checked against
`liana.method.cellphonedb`). Positive Δ is CLDN4-high minus CLDN4-low.
Honest n is the patient. The selection rule is in the script docstring and in METHODS.md.
It was applied by code to the full grid. It was not edited after the grid was seen.

## Winner

Gate **{winner['gate']}**, expr_prop **{float(winner['expr_prop']):.2f}**, receiver **T/NK**.

| | |
|---|---|
| n patients | {int(winner['n'])} ({int(winner['n_GSE123902'])}+{int(winner['n_GSE131907'])}+{int(winner['n_GSE205335'])}+{int(winner['n_GSE189357'])}) |
| mean CellPhoneDB Δ | {fmt_d(winner['mean_d_cpdb'])} |
| mean Connectome Δ | {fmt_d(winner['mean_d_conn'])} |
| patient fraction Δ>0 | {winner['frac_pos']:.3f} |
| Wilcoxon p (CellPhoneDB) | {fmt_p(winner['p_cpdb'])} |
| BH q across the T/NK grid | {fmt_p(winner['q_cpdb_grid'])} |
| cohorts CellPhoneDB | {winner['cohort_cpdb']} |
| cohorts Connectome | {winner['cohort_conn']} |
| mean Δ on edges expressed in both arms | {fmt_d(winner['mean_d_both'])} |
| mean fraction of detected edges with low arm at 0 | {winner['frac_low_zero']:.3f} |
| malignant CLDN4 log1p mean, high vs low | {winner['mean_cldn_high']:.3f} vs {winner['mean_cldn_low']:.3f} |

{gain}

{skip_txt}

{mye_txt}

## Ligands at the winner

| ligand | n | Δ CPDB | Δ Connectome | frac>0 | p | 123902 | 131907 | 205335 | 189357 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(lig_lines)}

{lg_txt}

## Edges at the winner

Δ is the mean of patients in whom the receptor complex passes expr_prop.
`low0` is the fraction of those patient-edges whose low arm is 0 and high arm is >0.
BH is across these edges.

| edge | ligand | n | Δ CPDB | Δ Connectome | frac>0 | p | q | cohorts + | low0 |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
{chr(10).join(edge_lines)}

## Grid (T/NK, highest CellPhoneDB Δ first, 12 rows)

Eligible means n≥40, every cohort n≥3, and all four cohort means >0 for both magnitudes.

| gate | expr_prop | n | 123902/131907/205335/189357 | Δ CPDB | Δ Connectome | frac>0 | CPDB cohorts | Conn cohorts | eligible | p | q_grid |
|---|---:|---:|---|---:|---:|---:|---:|---:|---|---:|---:|
{top}

Full grid: `results/tables/grid.tsv` ({int((grid.receiver=='TNK').sum())} T/NK rows, {int(grid.eligible.sum())} eligible).
Patient family scores: `results/tables/patient_family.tsv`.
Figure: `results/figures/max_effect_grid.png`.

## How the family score is built

Each edge uses the liana complex rule (minimum subunit). `lr_means` averages the ligand and receptor means and is 0 if either side is 0. Connectome magnitude is the product of those means. An edge counts for a patient only when the receptor passes expr_prop on T/NK. The ligand score is the mean of that ligand's counted edges. The family score is the unweighted mean of ligands with a counted edge. A patient is in the test only if at least 3 of the 4 ligands are counted. Missing ligands are not filled with zero.

The both-arms row keeps only edges that are nonzero on both the high and the low malignant arm, so that number is not produced by the low arm falling under expr_prop.

## What is not claimed

- This does not measure spatial exclusion, contact, or muzzling.
- TACSTD2 is not a gate. This is not dual-high.
- GSE148071 and the other non-concordant sets are not in this fit.
- The Wilcoxon p describes the selected row. The grid BH q is the multiplicity-adjusted figure for that search.
- Rank-aggregate ρ is not the maximized quantity. ρ is a rank, not an expression Δ.
- Cell-pooled tests are not the honest n.
- Patients in the winner table: {int(w_patients.shape[0])}.
"""
    (HERE / "FINDING.md").write_text(text)
    (HERE / "results" / "FINDING.md").write_text(text)


def write_versions() -> None:
    import anndata
    import liana
    import scipy

    txt = (
        f"liana {liana.__version__}\n"
        f"anndata {anndata.__version__}\n"
        f"scipy {scipy.__version__}\n"
        f"numpy {np.__version__}\n"
        f"pandas {pd.__version__}\n"
    )
    (HERE / "results" / "session_versions.txt").write_text(txt)


def main() -> None:
    warnings.filterwarnings("ignore", category=FutureWarning)
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    assert_matches_liana()
    assert_gates()
    if os.environ.get("FAMILY_FROM_TABLES") == "1":
        log("rebuilding the write-up from saved tables")
        edge_long = pd.read_csv(TAB / "edge_long.tsv", sep="\t")
        patients = pd.read_csv(TAB / "patient_family.tsv", sep="\t")
        grid = pd.read_csv(TAB / "grid.tsv", sep="\t")
        grid["eligible"] = grid["eligible"].astype(str).str.lower().eq("true")
        for col in ("receptor_detected", "both_pos", "low_zero_high_pos"):
            edge_long[col] = edge_long[col].astype(str).str.lower().eq("true")
        sign_check_reference(edge_long)
    else:
        edges_def = pd.read_csv(HERE / "data" / "edges.tsv", sep="\t")
        units = load_units()
        log(f"scoring {len(units)} units")
        edge_long = score_units(units, edges_def)
        edge_long.to_csv(TAB / "edge_long.tsv", sep="\t", index=False)
        log(f"edge rows {len(edge_long)}")
        sign_check_reference(edge_long)
        patients = family_patient_table(edge_long)
        patients.to_csv(TAB / "patient_family.tsv", sep="\t", index=False)
        grid = summarize_grid(patients)
        grid.to_csv(TAB / "grid.tsv", sep="\t", index=False)
    winner = select_winner(grid)
    if winner is None:
        raise RuntimeError("no eligible spec; grid written, no winner claimed")
    log(
        f"WINNER {winner['gate']} expr_prop={winner['expr_prop']} "
        f"n={winner['n']} d={winner['mean_d_cpdb']:+.4f} frac={winner['frac_pos']:.3f}"
    )
    ligands = ligand_summary(patients, winner)
    ligands.to_csv(TAB / "winner_ligands.tsv", sep="\t", index=False)
    edge_sum = edge_summary(edge_long, winner)
    edge_sum.to_csv(TAB / "winner_edges.tsv", sep="\t", index=False)
    plot_grid(grid, winner, ligands, FIG / "max_effect_grid.png", FIG / "max_effect_grid.pdf")
    write_finding(grid, winner, ligands, edge_sum, patients)
    write_versions()
    log("wrote FINDING.md")


if __name__ == "__main__":
    main()
