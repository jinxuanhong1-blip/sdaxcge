#!/usr/bin/env python3
"""Last sweep of GSE22493 after the HGNC remap stayed discordant.

Pre-specified rule, fixed before looking at these extra results:

A processing is thesis-aligned (IFN/APM up after CLDN4 siRNA) only when
BOTH the IFN_IMMUNE set and the APM set pass on the same processing, and
  * CLDN4 is negative on every array where it is measured (siRNA < OE)
  * the processing is an experiment-level summary (not one array alone)
  * the set effect is positive
  * the one-sided up-test p < 0.05

Gene/probe tests use a one-sided Mann-Whitney against the other genes or
probes. ssGSEA uses a one-sided permutation p on the gene-median ranking
(seed 42, 999 permutations) plus a positive median per-array ES.

Flipping GSM558701 back to the raw Ch2/Ch1 sign is recorded and does NOT
count: the deposited matrix is already that array's normalized log ratio
with the sign reversed (dye-swap correction; Spearman of deposited vs
-Ch2 N Log Ratio is +1). Undoing it makes CLDN4 positive on that array.

Nothing here is invented. If no processing passes, the script writes
FINAL_DISCORDANT.
"""

from __future__ import annotations

import gzip
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze as A  # noqa: E402

STING = [
    "CGAS",
    "STING1",
    "TBK1",
    "IKBKE",
    "IRF3",
    "MAVS",
    "DDX58",
    "IFIH1",
    "TREX1",
    "IFI16",
    "ZBP1",
    "IRF7",
    "CXCL10",
    "IFNB1",
]
NHEJ = [
    "XRCC6",
    "XRCC5",
    "PRKDC",
    "DCLRE1C",
    "LIG4",
    "XRCC4",
    "NHEJ1",
    "PAXX",
    "POLL",
    "POLM",
    "PNKP",
    "APTX",
]
PANELS = {
    "PRIORITY": A.PRIORITY,
    "APM": A.APM,
    "IFN_IMMUNE": A.IFN_IMMUNE,
    "TJ_CORE": A.TJ_CORE,
    "STING": STING,
    "NHEJ": NHEJ,
}
# Thesis success requires these two, both up.
REQUIRED = ("IFN_IMMUNE", "APM")
UP_PANELS = {"PRIORITY", "APM", "IFN_IMMUNE", "STING", "NHEJ"}
N_PERM = 999
SEED = 42


def canon_list(genes: list[str], hgnc: dict) -> list[str]:
    out = []
    seen = set()
    for g in genes:
        got, _how = A.resolve_token(g, hgnc)
        sym = got or g
        if sym not in seen:
            seen.add(sym)
            out.append(sym)
    return out


def load_raw_probe(path: Path) -> dict[str, dict]:
    with gzip.open(path, "rt", errors="replace") as fh:
        lines = fh.read().splitlines()
    header = None
    start = None
    for i, ln in enumerate(lines):
        if ln.startswith("BEGIN DATA"):
            header = lines[i + 1].split("\t")
            start = i + 2
            break
    if header is None or start is None:
        raise SystemExit(f"No data block in {path}")
    idx = {c: j for j, c in enumerate(header)}
    out: dict[str, dict] = {}
    for ln in lines[start:]:
        if ln.startswith("END DATA") or not ln.strip():
            break
        fields = ln.split("\t")
        if not fields or not fields[0].isdigit():
            continue
        probe = fields[idx["Index"]]

        def num(col: str) -> float | None:
            try:
                v = float(fields[idx[col]])
            except (ValueError, KeyError):
                return None
            if not math.isfinite(v):
                return None
            return v

        ch1 = num("Ch1 Median - B")
        ch2 = num("Ch2 Median - B")
        log_ratio = None
        if ch1 is not None and ch2 is not None and ch1 > 0 and ch2 > 0:
            log_ratio = math.log2(ch2 / ch1)
        out[probe] = {
            "log2": log_ratio,
            "flag": fields[idx["Flags"]].strip(),
            "nlog": num("Ch2 N Log Ratio"),
        }
    return out


def matrix_from_columns(columns: list[dict[str, float | None]]) -> dict[str, list[float | None]]:
    probes = set()
    for col in columns:
        probes.update(col)
    out = {}
    for probe in probes:
        out[probe] = [col.get(probe) for col in columns]
    return out


def effects_from_df(df: pd.DataFrame, cols: list[str], min_n: int) -> dict[str, float]:
    effects = {}
    for _, row in df.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if v is not None and isinstance(v, (int, float)) and math.isfinite(float(v)):
                vals.append(float(v))
        if len(vals) >= min_n:
            effects[str(row["symbol"])] = float(np.median(vals))
    return effects


def cldn4_values(df: pd.DataFrame, cols: list[str]) -> list[float]:
    hit = df[df["symbol"] == "CLDN4"]
    if hit.empty:
        return []
    row = hit.iloc[0]
    vals = []
    for c in cols:
        v = row[c]
        if v is not None and isinstance(v, (int, float)) and math.isfinite(float(v)):
            vals.append(float(v))
    return vals


def orientation_ok(cldn_vals: list[float]) -> bool:
    if not cldn_vals:
        return False
    return all(v < 0 for v in cldn_vals)


def mw_up(panel_vals: list[float], background: list[float]) -> float:
    return A.safe_mw(panel_vals, background, "greater")


def ssgsea_es(scores: np.ndarray, in_set: np.ndarray, alpha: float = 0.25) -> float:
    """Barbie / GSVA ssGSEA. Positive ES: set genes rank toward high scores."""
    scores = np.asarray(scores, dtype=float)
    in_set = np.asarray(in_set, dtype=bool)
    n = scores.size
    n_hit = int(in_set.sum())
    if n_hit < 5 or n_hit >= n:
        return float("nan")
    order = np.argsort(-scores, kind="mergesort")
    inset = in_set[order]
    ranks = np.empty(n, dtype=float)
    ranks[order] = np.arange(n, 0, -1, dtype=float)
    weights = np.power(ranks[order], alpha)
    hit_w = np.where(inset, weights, 0.0)
    hit_sum = float(hit_w.sum())
    if hit_sum <= 0:
        return float("nan")
    phit = np.cumsum(hit_w) / hit_sum
    pmiss = np.cumsum(~inset) / float(n - n_hit)
    return float(np.sum(phit - pmiss))


def perm_p_up(scores: np.ndarray, in_set: np.ndarray, observed: float, rng: np.random.Generator) -> float:
    if not math.isfinite(observed):
        return float("nan")
    null = np.empty(N_PERM, dtype=float)
    for i in range(N_PERM):
        null[i] = ssgsea_es(rng.permutation(scores), in_set)
    # One-sided: enrichment at the high-score end.
    return float((np.sum(null >= observed) + 1) / (N_PERM + 1))


def gsva_n3(expr: np.ndarray, in_set: np.ndarray) -> float:
    """Across-sample rank GSVA-style KS, then the median over the 3 samples.

    With three arrays each gene has only three rank levels, so this is a
    coarse check, not a substitute for the ratio test. Positive means the
    set is enriched among genes that are high in that sample relative to
    the other samples.
    """
    # expr: genes x samples, finite
    n_gene, n_s = expr.shape
    if n_s < 2 or int(in_set.sum()) < 5:
        return float("nan")
    ranks = pd.DataFrame(expr).rank(axis=1, method="average").to_numpy()
    centered = (ranks - (n_s + 1) / 2.0) / (n_s / 4.0)
    sample_es = []
    n_hit = int(in_set.sum())
    for s in range(n_s):
        order = np.argsort(-centered[:, s], kind="mergesort")
        inset = in_set[order]
        hit = np.where(inset, 1.0 / n_hit, 0.0)
        miss = np.where(~inset, 1.0 / (n_gene - n_hit), 0.0)
        sample_es.append(float(np.sum(np.cumsum(hit - miss))))
    return float(np.median(sample_es))


def quantile_normalize(mat: np.ndarray) -> np.ndarray:
    """Column-wise quantile normalization. NaNs stay missing."""
    out = mat.copy()
    grid = np.linspace(0.0, 1.0, 1000)
    refs = []
    masks = []
    for j in range(mat.shape[1]):
        col = mat[:, j]
        mask = np.isfinite(col)
        masks.append(mask)
        vals = np.sort(col[mask])
        if vals.size < 10:
            refs.append(None)
            continue
        refs.append(np.interp(grid, np.linspace(0.0, 1.0, vals.size), vals))
    good = [r for r in refs if r is not None]
    if not good:
        return out
    ref = np.mean(np.vstack(good), axis=0)
    for j in range(mat.shape[1]):
        mask = masks[j]
        if int(mask.sum()) < 2 or refs[j] is None:
            continue
        ranks = stats.rankdata(mat[mask, j], method="average")
        pct = (ranks - 1.0) / max(int(mask.sum()) - 1, 1)
        out[mask, j] = np.interp(pct, grid, ref)
    return out


def df_to_matrix(df: pd.DataFrame) -> tuple[list[str], np.ndarray]:
    use = df.dropna(subset=A.GSMS, how="all").copy()
    symbols = use["symbol"].tolist()
    mat = use[A.GSMS].to_numpy(dtype=float)
    return symbols, mat


def collapse_probe_matrix(probes: list[dict], matrix: dict[str, list[float | None]]) -> pd.DataFrame:
    return A.collapse_genes(probes, matrix, "symbol")


def probe_effects(probes: list[dict], matrix: dict[str, list[float | None]], cols_idx: list[int], min_n: int) -> dict[str, list[float]]:
    """symbol -> list of probe medians (one per probe)."""
    by: dict[str, list[float]] = {}
    for rec in probes:
        sym = rec["symbol"]
        if not sym:
            continue
        row = matrix.get(rec["probe"])
        if row is None:
            continue
        vals = []
        for j in cols_idx:
            v = row[j]
            if v is not None and math.isfinite(v):
                vals.append(v)
        if len(vals) >= min_n:
            by.setdefault(sym, []).append(float(np.median(vals)))
    return by


def main() -> None:
    hgnc = A.load_hgnc(A.CACHE / "hgnc_complete_set.txt")
    platform = A.load_platform(A.CACHE / "GPL10555_family.soft.gz")
    _samples, deposited = A.load_matrix(A.CACHE / "GSE22493_series_matrix.txt.gz")
    probes = A.annotate_probes(platform, hgnc)
    panel_ids = {name: canon_list(genes, hgnc) for name, genes in PANELS.items()}

    raw_cols = []
    flag_cols = []
    for gsm in A.GSMS:
        raw = load_raw_probe(A.CACHE / "raw" / f"{gsm}.txt.gz")
        raw_cols.append({pid: rec["log2"] for pid, rec in raw.items()})
        flag_cols.append({pid: rec["log2"] if rec["flag"] == "3" else None for pid, rec in raw.items()})
    raw_matrix = matrix_from_columns(raw_cols)
    flag_matrix = matrix_from_columns(flag_cols)

    # Author dye-swap: deposited GSM558701 ranks match -Ch2/Ch1.
    flipped = {}
    for pid, row in raw_matrix.items():
        flipped[pid] = [
            row[0],
            None if row[1] is None else -row[1],
            row[2],
        ]
    flipped_flag = {}
    for pid, row in flag_matrix.items():
        flipped_flag[pid] = [
            row[0],
            None if row[1] is None else -row[1],
            row[2],
        ]

    def center_matrix(matrix: dict[str, list[float | None]]) -> dict[str, list[float | None]]:
        cols = []
        for j in range(3):
            vals = [row[j] for row in matrix.values() if row[j] is not None and math.isfinite(row[j])]
            cols.append(float(np.median(vals)) if vals else 0.0)
        out = {}
        for pid, row in matrix.items():
            out[pid] = [
                None if row[j] is None or not math.isfinite(row[j]) else row[j] - cols[j]
                for j in range(3)
            ]
        return out

    gene_tables = {
        "deposited": collapse_probe_matrix(probes, deposited),
        "deposited_median_center": collapse_probe_matrix(probes, center_matrix(deposited)),
        "raw_ch2_over_ch1": collapse_probe_matrix(probes, raw_matrix),
        "raw_flag3": collapse_probe_matrix(probes, flag_matrix),
        "raw_701_flipped_to_deposited_sign": collapse_probe_matrix(probes, flipped),
        "raw_flag3_701_flipped": collapse_probe_matrix(probes, flipped_flag),
        "raw_flag3_701_flipped_median_center": collapse_probe_matrix(probes, center_matrix(flipped_flag)),
    }

    # Quantile-normalize the deposited gene x array matrix, then re-median.
    symbols, mat = df_to_matrix(gene_tables["deposited"])
    qn = quantile_normalize(mat)
    qn_df = pd.DataFrame(qn, columns=A.GSMS)
    qn_df.insert(0, "symbol", symbols)
    qn_df["n_probes"] = 1
    qn_df["n_arrays"] = np.isfinite(qn).sum(axis=1)
    qn_df["tokens"] = ""
    qn_df["median_log2"] = np.nanmedian(qn, axis=1)
    qn_df["mean_log2"] = np.nanmean(qn, axis=1)
    gene_tables["deposited_quantile_norm"] = qn_df

    rows: list[dict] = []

    def add_effect_row(
        method: str,
        scope: str,
        unit: str,
        effects: dict[str, float],
        cldn_vals: list[float],
        panel: str,
        panel_syms: list[str],
        note: str,
    ) -> None:
        vals = [effects[g] for g in panel_syms if g in effects and math.isfinite(effects[g])]
        bg = [v for s, v in effects.items() if s not in set(panel_syms) and math.isfinite(v)]
        med = float(np.median(vals)) if vals else float("nan")
        n_up = sum(v > 0 for v in vals)
        n_down = sum(v < 0 for v in vals)
        p = mw_up(vals, bg)
        ok = orientation_ok(cldn_vals)
        # TJ thesis is down; this sweep's pass flag is the IFN/APM up rule only.
        want_up = panel in UP_PANELS
        passed = bool(
            scope == "experiment"
            and ok
            and want_up
            and math.isfinite(med)
            and med > 0
            and math.isfinite(p)
            and p < 0.05
        )
        rows.append(
            {
                "method": method,
                "scope": scope,
                "unit": unit,
                "panel": panel,
                "n": len(vals),
                "n_up": n_up,
                "n_down": n_down,
                "effect": med,
                "effect_kind": "median_log2",
                "p_up": p,
                "cldn4_values": ";".join(f"{v:.4f}" for v in cldn_vals),
                "cldn4_median": float(np.median(cldn_vals)) if cldn_vals else float("nan"),
                "orientation_ok": "yes" if ok else "no",
                "pass_panel_up": "yes" if passed else "no",
                "note": note,
            }
        )

    experiment_methods = [
        ("deposited", "deposited series-matrix VALUE; primary analysis repeated"),
        ("deposited_median_center", "deposited ratios, each array median-centered"),
        ("deposited_quantile_norm", "deposited gene matrix quantile-normalized across arrays"),
        ("raw_ch2_over_ch1", "ScanArray log2(Ch2/Ch1), both channels > 0; GSM558701 NOT flipped"),
        ("raw_flag3", "same raw ratio, ScanArray flag 3 only; GSM558701 NOT flipped"),
        (
            "raw_701_flipped_to_deposited_sign",
            "raw log2(Ch2/Ch1) with GSM558701 negated so it matches the deposited dye-swap sign",
        ),
        ("raw_flag3_701_flipped", "flag-3 raw ratios with GSM558701 negated"),
        (
            "raw_flag3_701_flipped_median_center",
            "flag-3 raw ratios, GSM558701 negated, then each array median-centered",
        ),
    ]
    # Unflipped raw does not preserve per-array CLDN4 sign. The orientation
    # function catches that. Flipped raw is a reconstruction of the deposit,
    # not a new biological contrast.
    for method, note in experiment_methods:
        df = gene_tables[method]
        eff = effects_from_df(df, A.GSMS, min_n=2)
        cldn = cldn4_values(df, A.GSMS)
        for panel, syms in panel_ids.items():
            add_effect_row(method, "experiment", "gene", eff, cldn, panel, syms, note)

    # Leave-one-out. Dropping one array was partly reported (700+702). All three are here.
    for drop_i, drop in enumerate(A.GSMS):
        keep = [g for g in A.GSMS if g != drop]
        df = gene_tables["deposited"]
        eff = effects_from_df(df, keep, min_n=2)
        cldn = cldn4_values(df, keep)
        note = f"leave-one-out drop {drop}; both remaining arrays required"
        for panel, syms in panel_ids.items():
            add_effect_row(f"loo_drop_{drop}", "experiment", "gene", eff, cldn, panel, syms, note)

    # Single arrays. Scope replicate: cannot carry the thesis by itself.
    df = gene_tables["deposited"]
    for gsm in A.GSMS:
        eff = effects_from_df(df, [gsm], min_n=1)
        cldn = cldn4_values(df, [gsm])
        for panel, syms in panel_ids.items():
            add_effect_row(
                f"single_{gsm}",
                "replicate",
                "gene",
                eff,
                cldn,
                panel,
                syms,
                "one array only; not an experiment-level call",
            )

    # Probe-level on the deposited matrix. Each probe is a unit.
    probe_by = probe_effects(probes, deposited, [0, 1, 2], min_n=2)
    # CLDN4 probe values
    cldn_df = gene_tables["deposited"]
    cldn = cldn4_values(cldn_df, A.GSMS)
    for panel, syms in panel_ids.items():
        vals = []
        for g in syms:
            vals.extend(probe_by.get(g, []))
        bg = []
        in_set = set(syms)
        for g, plist in probe_by.items():
            if g not in in_set:
                bg.extend(plist)
        med = float(np.median(vals)) if vals else float("nan")
        p = mw_up(vals, bg)
        ok = orientation_ok(cldn)
        passed = bool(
            ok and panel in UP_PANELS and math.isfinite(med) and med > 0 and math.isfinite(p) and p < 0.05
        )
        rows.append(
            {
                "method": "deposited_probe_level",
                "scope": "experiment",
                "unit": "probe",
                "panel": panel,
                "n": len(vals),
                "n_up": sum(v > 0 for v in vals),
                "n_down": sum(v < 0 for v in vals),
                "effect": med,
                "effect_kind": "median_of_probe_medians",
                "p_up": p,
                "cldn4_values": ";".join(f"{v:.4f}" for v in cldn),
                "cldn4_median": float(np.median(cldn)) if cldn else float("nan"),
                "orientation_ok": "yes" if ok else "no",
                "pass_panel_up": "yes" if passed else "no",
                "note": "no gene collapse; each mapped probe is one observation",
            }
        )

    # Quantile thresholds on |gene median log2|.
    dep_eff = effects_from_df(gene_tables["deposited"], A.GSMS, min_n=2)
    abs_all = np.array([abs(v) for v in dep_eff.values() if math.isfinite(v)])
    cldn_abs = abs(float(np.median(cldn))) if cldn else float("nan")
    thresholds = {
        "q50_abs": float(np.quantile(abs_all, 0.50)),
        "q75_abs": float(np.quantile(abs_all, 0.75)),
        "q90_abs": float(np.quantile(abs_all, 0.90)),
        "abs_ge_CLDN4": cldn_abs,
    }
    for thr_name, thr in thresholds.items():
        kept = {s: v for s, v in dep_eff.items() if math.isfinite(v) and abs(v) >= thr}
        for panel, syms in panel_ids.items():
            vals = [kept[g] for g in syms if g in kept]
            # Background is other large-effect genes, so the test asks whether
            # the large IFN/APM effects are shifted up relative to other large effects.
            bg = [v for s, v in kept.items() if s not in set(syms)]
            med = float(np.median(vals)) if vals else float("nan")
            p = mw_up(vals, bg)
            n_up = sum(v > 0 for v in vals)
            n_nz = sum(v != 0 for v in vals)
            sign_p = (
                float(stats.binomtest(n_up, n_nz, 0.5, alternative="greater").pvalue) if n_nz else float("nan")
            )
            ok = orientation_ok(cldn)
            # Pass if the large-effect subset is up vs other large-effect genes,
            # or, when MW is unstable because n is small, do not pass.
            passed = bool(
                ok and panel in UP_PANELS and math.isfinite(med) and med > 0 and math.isfinite(p) and p < 0.05
            )
            rows.append(
                {
                    "method": f"deposited_abs_threshold_{thr_name}",
                    "scope": "experiment",
                    "unit": "gene",
                    "panel": panel,
                    "n": len(vals),
                    "n_up": n_up,
                    "n_down": sum(v < 0 for v in vals),
                    "effect": med,
                    "effect_kind": "median_log2_among_|effect|>=threshold",
                    "p_up": p,
                    "cldn4_values": ";".join(f"{v:.4f}" for v in cldn),
                    "cldn4_median": float(np.median(cldn)) if cldn else float("nan"),
                    "orientation_ok": "yes" if ok else "no",
                    "pass_panel_up": "yes" if passed else "no",
                    "note": f"threshold={thr:.4f}; sign-test p_up={sign_p:.4g}; genes below the bar are ignored",
                }
            )

    # ssGSEA on deposited per-array rankings, significance on the gene-median ranking.
    rng = np.random.default_rng(SEED)
    dep_df = gene_tables["deposited"]
    # Per-array scores: genes finite on that array.
    for panel, syms in panel_ids.items():
        per_es = []
        for gsm in A.GSMS:
            sub = dep_df[np.isfinite(dep_df[gsm].to_numpy(dtype=float))]
            scores = sub[gsm].to_numpy(dtype=float)
            inset = sub["symbol"].isin(syms).to_numpy()
            per_es.append(ssgsea_es(scores, inset))
        # Experiment ranking: gene median across arrays with >=2 values.
        scores = np.array([dep_eff[s] for s in dep_eff], dtype=float)
        names = list(dep_eff)
        inset = np.array([s in set(syms) for s in names])
        obs = ssgsea_es(scores, inset)
        p = perm_p_up(scores, inset, obs, rng)
        med_es = float(np.nanmedian(per_es))
        ok = orientation_ok(cldn)
        n_pos = sum(math.isfinite(e) and e > 0 for e in per_es)
        passed = bool(
            ok and panel in UP_PANELS and math.isfinite(med_es) and med_es > 0 and math.isfinite(p) and p < 0.05
        )
        rows.append(
            {
                "method": "ssgsea_deposited",
                "scope": "experiment",
                "unit": "enrichment",
                "panel": panel,
                "n": int(inset.sum()),
                "n_up": n_pos,
                "n_down": sum(math.isfinite(e) and e < 0 for e in per_es),
                "effect": med_es,
                "effect_kind": "median_per_array_ssgsea_ES",
                "p_up": p,
                "cldn4_values": ";".join(f"{v:.4f}" for v in cldn),
                "cldn4_median": float(np.median(cldn)) if cldn else float("nan"),
                "orientation_ok": "yes" if ok else "no",
                "pass_panel_up": "yes" if passed else "no",
                "note": (
                    "alpha=0.25; per-array ES "
                    + ",".join("NA" if not math.isfinite(e) else f"{e:.3f}" for e in per_es)
                    + f"; perm p is on the gene-median ranking (n_perm={N_PERM}, seed={SEED}); "
                    + f"ranking ES={obs:.3f}"
                ),
            }
        )

    # n=3 GSVA-style score on the deposited gene matrix (complete cases).
    complete = dep_df.dropna(subset=A.GSMS)
    expr = complete[A.GSMS].to_numpy(dtype=float)
    complete_symbols = complete["symbol"].tolist()
    for panel, syms in panel_ids.items():
        inset = np.array([s in set(syms) for s in complete_symbols])
        es = gsva_n3(expr, inset)
        # Permute gene labels on the same coarse score.
        null = np.empty(N_PERM, dtype=float)
        n_hit = int(inset.sum())
        for i in range(N_PERM):
            pick = np.zeros(len(complete_symbols), dtype=bool)
            choose = rng.choice(len(complete_symbols), size=n_hit, replace=False)
            pick[choose] = True
            null[i] = gsva_n3(expr, pick)
        p = float((np.sum(null >= es) + 1) / (N_PERM + 1)) if math.isfinite(es) else float("nan")
        ok = orientation_ok(cldn)
        passed = bool(ok and panel in UP_PANELS and math.isfinite(es) and es > 0 and math.isfinite(p) and p < 0.05)
        rows.append(
            {
                "method": "gsva_n3_deposited",
                "scope": "experiment",
                "unit": "enrichment",
                "panel": panel,
                "n": n_hit,
                "n_up": "" if not math.isfinite(es) else (1 if es > 0 else 0),
                "n_down": "" if not math.isfinite(es) else (1 if es < 0 else 0),
                "effect": es,
                "effect_kind": "median_across_samples_of_KS",
                "p_up": p,
                "cldn4_values": ";".join(f"{v:.4f}" for v in cldn),
                "cldn4_median": float(np.median(cldn)) if cldn else float("nan"),
                "orientation_ok": "yes" if ok else "no",
                "pass_panel_up": "yes" if passed else "no",
                "note": "n=3 across-sample ranks only; three rank levels; not a rescue by itself unless p<0.05 and ES>0",
            }
        )

    # Within-array mean rank (1 = lowest ratio). High mean rank = up.
    rank_effects: dict[str, float] = {}
    rank_mat = {}
    for gsm in A.GSMS:
        col = dep_df[["symbol", gsm]].dropna()
        # rank 0..1
        r = stats.rankdata(col[gsm].to_numpy(dtype=float), method="average")
        r = (r - 1) / max(len(r) - 1, 1)
        for sym, rv in zip(col["symbol"], r):
            rank_mat.setdefault(sym, []).append(float(rv))
    for sym, vals in rank_mat.items():
        if len(vals) >= 2:
            rank_effects[sym] = float(np.median(vals))
    for panel, syms in panel_ids.items():
        vals = [rank_effects[g] for g in syms if g in rank_effects]
        bg = [v for s, v in rank_effects.items() if s not in set(syms)]
        med = float(np.median(vals)) if vals else float("nan")
        p = mw_up(vals, bg)
        # Rank median > 0.5 means the set sits above the middle of the array.
        ok = orientation_ok(cldn)
        passed = bool(ok and panel in UP_PANELS and math.isfinite(med) and med > 0.5 and math.isfinite(p) and p < 0.05)
        rows.append(
            {
                "method": "deposited_within_array_rank",
                "scope": "experiment",
                "unit": "gene",
                "panel": panel,
                "n": len(vals),
                "n_up": sum(v > 0.5 for v in vals),
                "n_down": sum(v < 0.5 for v in vals),
                "effect": med,
                "effect_kind": "median_of_within_array_percentile",
                "p_up": p,
                "cldn4_values": ";".join(f"{v:.4f}" for v in cldn),
                "cldn4_median": float(np.median(cldn)) if cldn else float("nan"),
                "orientation_ok": "yes" if ok else "no",
                "pass_panel_up": "yes" if passed else "no",
                "note": "percentile within each array, then median; up means percentile > 0.5",
            }
        )

    sweep = pd.DataFrame(rows)
    # A method is thesis-aligned only if IFN_IMMUNE and APM both pass.
    aligned_methods = []
    for method, sub in sweep.groupby("method"):
        got = {}
        for panel in REQUIRED:
            hit = sub[sub["panel"] == panel]
            if hit.empty:
                got[panel] = False
            else:
                got[panel] = bool((hit["pass_panel_up"] == "yes").all())
        if all(got.values()):
            aligned_methods.append(method)
    verdict = "FINAL_DISCORDANT" if not aligned_methods else "THESIS_ALIGNED"
    sweep["final_verdict"] = verdict
    sweep["method_aligns_ifn_and_apm"] = sweep["method"].isin(aligned_methods).map({True: "yes", False: "no"})

    out = A.TAB
    out.mkdir(parents=True, exist_ok=True)
    sweep.to_csv(out / "sweep_final.tsv", sep="\t", index=False, float_format="%.6g")

    # Compact experiment-level IFN/APM/PRIORITY/TJ/STING/NHEJ table for the writeup.
    keep_methods_prefix_ok = sweep["scope"].eq("experiment")
    compact = sweep[keep_methods_prefix_ok & sweep["panel"].isin(["IFN_IMMUNE", "APM", "PRIORITY", "TJ_CORE", "STING", "NHEJ"])]
    compact.to_csv(out / "sweep_compact.tsv", sep="\t", index=False, float_format="%.6g")

    # Figure: IFN and APM effect for experiment-level gene/enrichment rows.
    fig_rows = sweep[
        sweep["panel"].isin(["IFN_IMMUNE", "APM"])
        & sweep["scope"].eq("experiment")
        & (
            sweep["effect_kind"].astype(str).str.contains("log2")
            | sweep["effect_kind"].eq("median_of_probe_medians")
        )
    ].copy()
    method_order = list(dict.fromkeys(fig_rows["method"]))
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    y_base = np.arange(len(method_order))
    for panel, color, shift in (("IFN_IMMUNE", "#d95f0e", -0.15), ("APM", "#2c7fb8", 0.15)):
        xs = []
        ys = []
        for i, method in enumerate(method_order):
            hit = fig_rows[(fig_rows["method"] == method) & (fig_rows["panel"] == panel)]
            if hit.empty:
                continue
            xs.append(float(hit.iloc[0]["effect"]))
            ys.append(i + shift)
        ax.scatter(xs, ys, c=color, s=28, label=panel, zorder=3)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y_base)
    ax.set_yticklabels(method_order, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("set median log2 (siRNA / CLDN4-OE)")
    ax.set_title("Last sweep: no processing puts both IFN and APM up")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    A.FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(A.FIG / "fig_sweep_ifn_apm.png", dpi=160)
    fig.savefig(A.FIG / "fig_sweep_ifn_apm.pdf")
    plt.close(fig)

    # Resolved set membership for the audit.
    membership = []
    for panel, syms in panel_ids.items():
        for s in syms:
            membership.append({"panel": panel, "hgnc_symbol": s, "in_deposited_ge2": "yes" if s in dep_eff else "no"})
    pd.DataFrame(membership).to_csv(out / "sweep_set_membership.tsv", sep="\t", index=False)

    print("VERDICT", verdict)
    print("aligned_methods", aligned_methods)
    show = sweep[(sweep["panel"].isin(REQUIRED)) & (sweep["scope"] == "experiment")][
        ["method", "panel", "n", "n_up", "n_down", "effect", "p_up", "orientation_ok", "pass_panel_up"]
    ]
    print(show.to_string(index=False))
    # Near misses: positive effect, orientation ok, experiment, but p>=0.05
    near = sweep[
        (sweep["panel"].isin(REQUIRED))
        & (sweep["scope"] == "experiment")
        & (sweep["orientation_ok"] == "yes")
        & (sweep["effect"] > 0)
    ]
    print("--- positive IFN/APM effects with orientation ok ---")
    if near.empty:
        print("(none)")
    else:
        print(near[["method", "panel", "n", "n_up", "n_down", "effect", "p_up", "pass_panel_up"]].to_string(index=False))


if __name__ == "__main__":
    main()
