#!/usr/bin/env python3
"""Score each downloaded cohort against the pre-specified TF network.

Hypothesis (fixed before looking at numbers)
-------------------------------------------
In public human and mouse lung RNA, ELF3, GRHL1, KLF4 and TFAP2A co-vary with
TACSTD2 and CLDN4, and NKX2-1 varies in the opposite direction.

This is a co-expression hunt, not a test of direct binding or regulation.

Pair rule
  * positive pair: Spearman ρ >= 0.30 and two-sided p < 0.05
  * NKX2-1 anti pair: Spearman ρ <= -0.20 and two-sided p < 0.05
  * n_used >= 12 after library-size and missingness filters

Cohort call
  * SUPPORTED          : >= 6 of 8 positive pairs AND both NKX2-1 pairs anti
  * PARTIAL            : >= 4 of 8 positive pairs (NKX2-1 not required)
  * NOT_SUPPORTED      : otherwise
  * UNINFORMATIVE      : targets or TFs undetectable, or n_used < 12

The same rules are applied to EPCAM-residual (rank-partial) correlations.
A SUPPORTED call that dies after EPCAM residualization is reported as
composition-confounded, not as a validated network.

Library sizes below 1e6 counts are dropped. A gene is treated as missing in a
sample if its CPM is 0; a gene used in a pair must be detected in >= 30% of
kept samples. log2(CPM+1) is the expression scale.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

TFS_POS = ["ELF3", "GRHL1", "KLF4", "TFAP2A"]
TARGETS = ["TACSTD2", "CLDN4"]
TF_NEG = "NKX2-1"
POS_PAIRS = [(tf, tgt) for tf in TFS_POS for tgt in TARGETS]
NEG_PAIRS = [(TF_NEG, tgt) for tgt in TARGETS]
COMP_MARKERS = ["EPCAM", "SFTPC", "SCGB1A1", "KRT5", "PTPRC", "COL1A1"]


def read_counts(path: Path) -> tuple[pd.Series, pd.DataFrame]:
    df = pd.read_csv(path, sep="\t", index_col=0)
    lib = df.loc["__libsize__"].astype(float)
    expr = df.drop(index="__libsize__").astype(float)
    return lib, expr


def log_cpm(expr: pd.DataFrame, lib: pd.Series) -> pd.DataFrame:
    return np.log2(expr.div(lib.clip(lower=1.0), axis=1) * 1e6 + 1.0)


def rank_resid(y: np.ndarray, z: np.ndarray) -> np.ndarray:
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    rz = np.column_stack([np.ones(len(rz)), rz])
    beta, *_ = np.linalg.lstsq(rz, ry, rcond=None)
    return ry - rz @ beta


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    if len(x) < 6 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan, np.nan
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    if len(x) < 8 or np.std(z) == 0:
        return np.nan, np.nan
    return spearman(rank_resid(x, z), rank_resid(y, z))


def detected_frac(cpm_like: pd.Series) -> float:
    # log2(CPM+1) == 0 iff CPM == 0
    return float((cpm_like > 0).mean())


def load_tcga_map(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    out: dict[str, str] = {}
    with path.open() as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for rec in rdr:
            st = rec.get("sample_type") or ""
            if st == "Solid Tissue Normal":
                label = "normal"
            elif "Tumor" in st or st == "Primary Tumor" or st == "Recurrent Tumor":
                label = "tumor"
            else:
                label = "other"
            out[rec["external_id"]] = label
            out[rec.get("tcga_barcode") or ""] = label
            # gene_sums columns sometimes carry a trailing .1
            out[rec["external_id"] + ".1"] = label
    out.pop("", None)
    return out


def tcga_sample_type(columns: list[str], mapping: dict[str, str]) -> dict[str, str]:
    out = {}
    for c in columns:
        if c in mapping:
            out[c] = mapping[c]
            continue
        key = c.split(".")[0]
        if key in mapping:
            out[c] = mapping[key]
            continue
        # last-resort barcode parse: TCGA-XX-XXXX-11A...
        toks = c.replace(".", "-").split("-")
        code = None
        for t in toks:
            if len(t) >= 2 and t[:2].isdigit():
                code = t[:2]
                break
        if code == "11":
            out[c] = "normal"
        elif code in {"01", "02", "03", "05", "06"}:
            out[c] = "tumor"
        else:
            out[c] = "unknown"
    return out


def split_cohorts(cid: str, source: str, lib: pd.Series, tcga_map: dict[str, str]) -> list[tuple[str, str, pd.Index]]:
    cols = lib.index
    if source == "tcga":
        types = tcga_sample_type(list(cols), tcga_map)
        tumor = [c for c in cols if types[c] == "tumor"]
        normal = [c for c in cols if types[c] == "normal"]
        out = []
        if len(tumor) >= 12:
            out.append((f"{cid}_tumor", "tumor", pd.Index(tumor)))
        if len(normal) >= 8:
            out.append((f"{cid}_normal", "adjacent_normal", pd.Index(normal)))
        if not out:
            out.append((cid, "tumor", cols))
        return out
    return [(cid, "tissue" if cid == "GTEx_LUNG" else "unknown", cols)]


def score_block(logcpm: pd.DataFrame, min_detect: float = 0.30) -> dict:
    n = logcpm.shape[1]
    detect = {g: detected_frac(logcpm.loc[g]) if g in logcpm.index else 0.0 for g in TFS_POS + TARGETS + [TF_NEG] + COMP_MARKERS}

    def pair_stats(a: str, b: str, covariate: str | None = None) -> dict:
        if a not in logcpm.index or b not in logcpm.index:
            return {"rho": np.nan, "p": np.nan, "n": 0, "detect_a": detect.get(a, 0.0), "detect_b": detect.get(b, 0.0)}
        if detect[a] < min_detect or detect[b] < min_detect:
            return {"rho": np.nan, "p": np.nan, "n": n, "detect_a": detect[a], "detect_b": detect[b]}
        x = logcpm.loc[a].to_numpy()
        y = logcpm.loc[b].to_numpy()
        if covariate and covariate in logcpm.index and detect[covariate] >= min_detect:
            z = logcpm.loc[covariate].to_numpy()
            rho, p = partial_spearman(x, y, z)
        else:
            rho, p = spearman(x, y)
        return {"rho": rho, "p": p, "n": n, "detect_a": detect[a], "detect_b": detect[b]}

    raw_pos = {f"{a}__{b}": pair_stats(a, b) for a, b in POS_PAIRS}
    raw_neg = {f"{a}__{b}": pair_stats(a, b) for a, b in NEG_PAIRS}
    epcam_pos = {f"{a}__{b}": pair_stats(a, b, "EPCAM") for a, b in POS_PAIRS}
    epcam_neg = {f"{a}__{b}": pair_stats(a, b, "EPCAM") for a, b in NEG_PAIRS}
    sftpc_pos = {f"{a}__{b}": pair_stats(a, b, "SFTPC") for a, b in POS_PAIRS}
    sftpc_neg = {f"{a}__{b}": pair_stats(a, b, "SFTPC") for a, b in NEG_PAIRS}

    def call(pos, neg) -> str:
        n_ok = sum(
            1
            for v in pos.values()
            if v["rho"] == v["rho"] and v["rho"] >= 0.30 and v["p"] < 0.05
        )
        n_anti = sum(
            1
            for v in neg.values()
            if v["rho"] == v["rho"] and v["rho"] <= -0.20 and v["p"] < 0.05
        )
        n_finite = sum(1 for v in pos.values() if v["rho"] == v["rho"])
        if n_finite < 4 or n < 12:
            return "UNINFORMATIVE"
        if n_ok >= 6 and n_anti == 2:
            return "SUPPORTED"
        if n_ok >= 4:
            return "PARTIAL"
        return "NOT_SUPPORTED"

    def mean_rho(d) -> float:
        xs = [v["rho"] for v in d.values() if v["rho"] == v["rho"]]
        return float(np.mean(xs)) if xs else float("nan")

    score = mean_rho(raw_pos) - mean_rho(raw_neg)

    # permutation of sample labels for the five TFs (4 pos + NKX2-1)
    rng = np.random.default_rng(1)
    tfs = [g for g in TFS_POS + [TF_NEG] if g in logcpm.index]
    if tfs and n >= 12:
        mat = logcpm.loc[tfs].to_numpy()
        tgt = {t: logcpm.loc[t].to_numpy() for t in TARGETS if t in logcpm.index}
        null = []
        for _ in range(300):
            perm = rng.permutation(n)
            ppos, pneg = [], []
            for tf in TFS_POS:
                if tf not in logcpm.index:
                    continue
                x = mat[tfs.index(tf), perm]
                for t, y in tgt.items():
                    rho, _ = spearman(x, y)
                    if rho == rho:
                        ppos.append(rho)
            if TF_NEG in logcpm.index:
                x = mat[tfs.index(TF_NEG), perm]
                for t, y in tgt.items():
                    rho, _ = spearman(x, y)
                    if rho == rho:
                        pneg.append(rho)
            if ppos and pneg:
                null.append(float(np.mean(ppos) - np.mean(pneg)))
        if null and score == score:
            p_perm = (1.0 + sum(s >= score for s in null)) / (1.0 + len(null))
            null_mean = float(np.mean(null))
            null_sd = float(np.std(null, ddof=1)) if len(null) > 1 else float("nan")
        else:
            p_perm = null_mean = null_sd = float("nan")
    else:
        p_perm = null_mean = null_sd = float("nan")

    comp = {}
    for marker in COMP_MARKERS:
        for gene in TFS_POS + TARGETS + [TF_NEG]:
            st = pair_stats(gene, marker)
            comp[f"{gene}__{marker}"] = st["rho"]

    return {
        "n_used": n,
        "call_raw": call(raw_pos, raw_neg),
        "call_epcam": call(epcam_pos, epcam_neg),
        "call_sftpc": call(sftpc_pos, sftpc_neg),
        "mean_rho_pos": mean_rho(raw_pos),
        "mean_rho_nkx": mean_rho(raw_neg),
        "score": score,
        "p_perm": p_perm,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "n_pos_pass": sum(1 for v in raw_pos.values() if v["rho"] == v["rho"] and v["rho"] >= 0.30 and v["p"] < 0.05),
        "n_neg_pass": sum(1 for v in raw_neg.values() if v["rho"] == v["rho"] and v["rho"] <= -0.20 and v["p"] < 0.05),
        "n_pos_pass_epcam": sum(1 for v in epcam_pos.values() if v["rho"] == v["rho"] and v["rho"] >= 0.30 and v["p"] < 0.05),
        "n_neg_pass_epcam": sum(1 for v in epcam_neg.values() if v["rho"] == v["rho"] and v["rho"] <= -0.20 and v["p"] < 0.05),
        "pairs_raw": {k: v for k, v in {**raw_pos, **raw_neg}.items()},
        "pairs_epcam": {k: v for k, v in {**epcam_pos, **epcam_neg}.items()},
        "pairs_sftpc": {k: v for k, v in {**sftpc_pos, **sftpc_neg}.items()},
        "detect": detect,
        "composition_rho": comp,
    }


def flatten_pairs(prefix: str, pairs: dict) -> dict:
    out = {}
    for k, v in pairs.items():
        out[f"{prefix}_{k}_rho"] = v["rho"]
        out[f"{prefix}_{k}_p"] = v["p"]
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohorts", default="results/hunt_tf/tables/cohorts_to_download.tsv")
    ap.add_argument("--matrix-dir", default="results/hunt_tf/matrices")
    ap.add_argument("--out-dir", default="results/hunt_tf/tables")
    args = ap.parse_args()

    cohorts = list(csv.DictReader(open(args.cohorts), delimiter="\t"))
    matrix_dir = Path(args.matrix_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    tcga_map = load_tcga_map(out_dir / "tcga_sample_map.tsv")

    summary_rows = []
    pair_rows = []
    long_json = []

    for rec in cohorts:
        cid = rec["cohort_id"]
        path = matrix_dir / f"{cid}.counts.tsv"
        if not path.exists():
            continue
        lib, expr = read_counts(path)
        keep = lib[lib >= 1_000_000].index
        if rec["source"] == "sra" and rec["selected_runs"]:
            nominated = set(rec["selected_runs"].split(","))
            keep = pd.Index([c for c in keep if c.split(".")[0] in nominated or c in nominated])
        if len(keep) < 8:
            keep = lib[lib >= 200_000].index
        expr = expr.loc[:, keep]
        lib = lib.loc[keep]
        logcpm = log_cpm(expr, lib)

        material_default = rec["material"]
        splits = split_cohorts(cid, rec["source"], lib, tcga_map)
        # split_cohorts for SRA should just return the cohort as-is
        if rec["source"] == "sra":
            splits = [(cid, material_default, lib.index)]
        elif rec["source"] == "gtex":
            splits = [(cid, "tissue", lib.index)]

        for split_id, material, cols in splits:
            block = logcpm.loc[:, cols]
            if block.shape[1] < 8:
                continue
            res = score_block(block)
            row = {
                "cohort_id": split_id,
                "parent_id": cid,
                "source": rec["source"],
                "organism": rec["organism"],
                "material": material,
                "study_title": rec["study_title"],
                "n_nominated": rec["n_nominated"],
                **{k: res[k] for k in (
                    "n_used", "call_raw", "call_epcam", "call_sftpc",
                    "mean_rho_pos", "mean_rho_nkx", "score", "p_perm",
                    "null_mean", "null_sd", "n_pos_pass", "n_neg_pass",
                    "n_pos_pass_epcam", "n_neg_pass_epcam",
                )},
            }
            for g, d in res["detect"].items():
                row[f"detect_{g}"] = d
            for k, v in res["composition_rho"].items():
                row[f"comp_{k}"] = v
            summary_rows.append(row)

            for scale, pairs in (("raw", res["pairs_raw"]), ("epcam", res["pairs_epcam"]), ("sftpc", res["pairs_sftpc"])):
                for pair, st in pairs.items():
                    a, b = pair.split("__")
                    pair_rows.append({
                        "cohort_id": split_id,
                        "organism": rec["organism"],
                        "material": material,
                        "scale": scale,
                        "tf": a,
                        "target": b,
                        "rho": st["rho"],
                        "p": st["p"],
                        "n": st["n"],
                        "detect_tf": st["detect_a"],
                        "detect_target": st["detect_b"],
                    })
            long_json.append({"cohort": row, "pairs": res})

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary = summary.sort_values(["organism", "score"], ascending=[True, False])
    summary.to_csv(out_dir / "cohort_scores.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(pair_rows).to_csv(out_dir / "pair_correlations.tsv", sep="\t", index=False, float_format="%.4f")
    (out_dir / "cohort_scores.json").write_text(json.dumps(long_json, indent=2, default=str) + "\n")

    if not summary.empty:
        print(summary.groupby(["organism", "call_raw"]).size().to_string())
        print("--- EPCAM-residual ---")
        print(summary.groupby(["organism", "call_epcam"]).size().to_string())


if __name__ == "__main__":
    main()
