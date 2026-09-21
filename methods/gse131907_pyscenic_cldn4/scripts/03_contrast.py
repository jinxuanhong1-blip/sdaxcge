#!/usr/bin/env python3
"""Patient-level CLDN4-high vs low contrast of pySCENIC AUCell.

Writes results/gse131907_pyscenic_cldn4/RESULTS.md from the computed tables.
Observational signs are reported next to the KD-like expectation and next to
the two public CLDN4-loss references. Those references are not fit here.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

PREP = Path("/tmp/gse131907_pyscenic/prepared")
ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "gse131907_pyscenic_cldn4"
TAB = RES / "tables"
FIG = RES / "figures"

IFN_TFS = {"STAT1", "STAT2", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", "IRF7", "IRF8", "IRF9"}
MHC_TFS = {"NLRC5", "CIITA", "RFX5", "RFXAP", "RFXANK"}
TJ_TFS = {"ELF3", "GRHL1", "GRHL2", "GRHL3", "KLF4", "KLF5", "OVOL1", "OVOL2", "TFAP2A", "TFAP2C"}
GIVEN_TFS = {"ELF3"}

# KD-like expectation used by prior AUCell work: IFN/MHC lower in CLDN4-high.
# GSE207704 CLDN4 CRISPR: IFN trends down after loss, which predicts the
# opposite observational sign (higher IFN in CLDN4-high).
# GSE50927 whole-lung Cldn4 KO, n=1: IFN/MHC up after loss, which predicts
# lower IFN in CLDN4-high. Not a cancer-cell replicate.


def fmt(x, sig: int = 3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not np.isfinite(v):
        return "NA"
    if v == 0:
        return "0"
    if abs(v) >= 0.001:
        return f"{v:.4g}"
    return f"{v:.2e}"


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    n = int(ok.sum())
    if n == 0:
        return out
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out[order] = np.clip(q, 0, 1)
    return out


def wilcoxon_paired(delta: np.ndarray) -> dict:
    d = np.asarray(delta, dtype=float)
    d = d[np.isfinite(d)]
    n = int(d.size)
    if n < 5:
        return {"n": n, "stat": np.nan, "p": np.nan, "delta_median": float(np.median(d)) if n else np.nan,
                "n_pos": int((d > 0).sum()) if n else 0, "n_neg": int((d < 0).sum()) if n else 0}
    med = float(np.median(d))
    n_pos = int((d > 0).sum())
    n_neg = int((d < 0).sum())
    if np.allclose(d, 0):
        return {"n": n, "stat": np.nan, "p": 1.0, "delta_median": 0.0, "n_pos": n_pos, "n_neg": n_neg}
    try:
        w = stats.wilcoxon(d, zero_method="wilcox", alternative="two-sided", method="auto")
        p = float(w.pvalue)
        stat = float(w.statistic)
    except ValueError:
        p, stat = np.nan, np.nan
    return {"n": n, "stat": stat, "p": p, "delta_median": med, "n_pos": n_pos, "n_neg": n_neg}


def assign_split(obs: pd.DataFrame) -> pd.DataFrame:
    """Within-sample CLDN4 split. Index matches obs.index."""
    high = pd.Series(np.nan, index=obs.index, dtype=float)
    paired = pd.Series(False, index=obs.index)
    rule = pd.Series("", index=obs.index, dtype=object)
    sample_rows = []
    for sample, idx in obs.groupby("sample").groups.items():
        idx = pd.Index(idx)
        if len(idx) < 20:
            sample_rows.append({
                "sample": sample, "n": int(len(idx)), "n_high": 0, "n_low": 0,
                "paired": False, "rule": "lt20", "median_cldn4": np.nan,
                "patient_id": obs.loc[idx, "patient_id"].iloc[0] if len(idx) else "",
                "origin": obs.loc[idx, "Sample_Origin"].iloc[0] if len(idx) else "",
            })
            continue
        x = obs.loc[idx, "cldn4_log1p_cp10k"].to_numpy(dtype=float)
        med = float(np.nanmedian(x))
        if np.isfinite(med) and med > 0:
            h = x >= med
            how = "within_sample_median"
        else:
            h = x > 0
            how = "gt0_because_median_0"
        finite = np.isfinite(x)
        h = h & finite
        low = (~h) & finite
        n_high = int(h.sum())
        n_low = int(low.sum())
        is_paired = n_high >= 20 and n_low >= 20
        high.loc[idx] = np.where(finite, h.astype(float), np.nan)
        rule.loc[idx] = how
        if is_paired:
            paired.loc[idx] = True
        sample_rows.append({
            "sample": sample,
            "n": int(finite.sum()),
            "n_high": n_high,
            "n_low": n_low,
            "paired": is_paired,
            "rule": how,
            "median_cldn4": med,
            "patient_id": obs.loc[idx, "patient_id"].iloc[0],
            "origin": obs.loc[idx, "Sample_Origin"].mode().iloc[0] if len(obs.loc[idx, "Sample_Origin"].mode()) else "",
            "mean_cldn4_high": float(np.nanmean(x[h])) if n_high else np.nan,
            "mean_cldn4_low": float(np.nanmean(x[low])) if n_low else np.nan,
        })
    out = obs.copy()
    out["high"] = high
    out["sample_paired"] = paired
    out["split_rule"] = rule
    return out, pd.DataFrame(sample_rows)


def patient_table(obs: pd.DataFrame, values: np.ndarray) -> pd.DataFrame:
    sub = obs.loc[obs["sample_paired"]].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["value"] = values[sub.index.to_numpy()]
    rows = []
    for pid, g in sub.groupby("patient_id"):
        hi = g.loc[g["high"] == 1, "value"].to_numpy(dtype=float)
        lo = g.loc[g["high"] == 0, "value"].to_numpy(dtype=float)
        if hi.size < 20 or lo.size < 20:
            continue
        # Cell-weighted means. Sample weight is recorded so a huge sample is visible.
        vc = g["sample"].value_counts(normalize=True)
        rows.append({
            "patient_id": pid,
            "n_high": int(hi.size),
            "n_low": int(lo.size),
            "n_samples": int(g["sample"].nunique()),
            "max_sample_fraction": float(vc.iloc[0]),
            "origins": ",".join(sorted(g["Sample_Origin"].dropna().unique())),
            "mean_high": float(np.nanmean(hi)),
            "mean_low": float(np.nanmean(lo)),
            "delta": float(np.nanmean(hi) - np.nanmean(lo)),
        })
    return pd.DataFrame(rows)


def sample_deltas(obs: pd.DataFrame, values: np.ndarray, sample_info: pd.DataFrame) -> np.ndarray:
    paired_samples = set(sample_info.loc[sample_info["paired"], "sample"])
    deltas = []
    for sample, g in obs.groupby("sample"):
        if sample not in paired_samples:
            continue
        v = values[g.index.to_numpy()]
        h = g["high"].to_numpy()
        if np.nanmean(h == 1) == 0 or np.nansum(h == 0) == 0:
            continue
        deltas.append(float(np.nanmean(v[h == 1]) - np.nanmean(v[h == 0])))
    return np.asarray(deltas, dtype=float)


def relation(delta: float, expected: str) -> str:
    if expected in {"null", "unspecified"} or not np.isfinite(delta):
        return "NA"
    if delta == 0:
        return "ZERO"
    higher = delta > 0
    if expected == "lower_in_high":
        return "OPPOSITE" if higher else "SAME"
    if expected == "higher_in_high":
        return "OPPOSITE" if not higher else "SAME"
    return "NA"


def axis_of_program(name: str) -> str:
    if "INTERFERON" in name or name.endswith("IFN_UNION") or name == "IFN_UNION":
        return "IFN"
    if name == "MHC1_APM":
        return "MHC"
    if name == "TJ_STRUCT":
        return "TJ"
    if name == "KERATIN":
        return "KERATIN"
    if name == "CTRL_RIBO":
        return "CTRL"
    return "OTHER"


def md_table(df: pd.DataFrame, cols: list[str]) -> str:
    use = df.loc[:, [c for c in cols if c in df.columns]].copy()
    header = "| " + " | ".join(use.columns) + " |"
    sep = "| " + " | ".join("---" for _ in use.columns) + " |"
    lines = [header, sep]
    for rec in use.itertuples(index=False):
        cells = []
        for v in rec:
            if isinstance(v, float):
                cells.append(fmt(v))
            else:
                cells.append(str(v).replace("|", "/"))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def plain_direction(delta: float) -> str:
    if not np.isfinite(delta):
        return "undefined"
    if delta > 0:
        return "higher in CLDN4-high"
    if delta < 0:
        return "higher in CLDN4-low"
    return "no difference"


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    obs = pd.read_csv(PREP / "obs_malignant.tsv.gz", sep="\t")
    auc = np.load(PREP / "auc.npy")
    cols = (PREP / "auc_columns.txt").read_text().splitlines()
    regs = pd.read_csv(PREP / "regulons.tsv", sep="\t")
    prep = json.loads((PREP / "prepare_summary.json").read_text())
    run = json.loads((PREP / "pyscenic_run.json").read_text()) if (PREP / "pyscenic_run.json").exists() else {}
    if len(obs) != auc.shape[0] or auc.shape[1] != len(cols):
        raise SystemExit(f"shape mismatch obs={len(obs)} auc={auc.shape} cols={len(cols)}")
    obs = obs.reset_index(drop=True)

    programs_used = {}
    pjson = PREP / "program_genes_used.json"
    if pjson.exists():
        programs_used = json.loads(pjson.read_text())
    ifn_genes = set(programs_used.get("IFN_UNION", []))
    if not ifn_genes:
        ifn_genes = set(programs_used.get("HALLMARK_INTERFERON_GAMMA_RESPONSE", [])) | set(
            programs_used.get("HALLMARK_INTERFERON_ALPHA_RESPONSE", [])
        )
    universe = set((PREP / "genes.txt").read_text().splitlines())
    universe_n = len(universe)

    obs_s, sample_info = assign_split(obs)
    # tLung-only split uses only tLung cells, so medians are not borrowed from metastases.
    tlung_obs = obs.loc[obs["Sample_Origin"] == "tLung"].copy()
    tlung_labeled, tlung_samples = assign_split(tlung_obs)

    # Overlap of each cisTarget regulon with the Hallmark IFN union.
    overlap_rows = []
    reg_by_name = {}
    for rec in regs.itertuples(index=False):
        targets = str(rec.targets).split(";") if isinstance(rec.targets, str) else []
        targets = [g for g in targets if g and g != "CLDN4" and g in universe]
        k = len(set(targets) & ifn_genes)
        n_set = len(set(targets))
        K = len(ifn_genes)
        if n_set and K and universe_n > K:
            # P(X >= k)
            p_ov = float(stats.hypergeom.sf(k - 1, universe_n, K, n_set)) if k else 1.0
        else:
            p_ov = np.nan
        overlap_rows.append({"name": rec.name, "tf": rec.tf, "n_targets_nocldn4": n_set,
                             "ifn_overlap": k, "ifn_overlap_p": p_ov, "nes": rec.nes,
                             "contains_CLDN4": rec.contains_CLDN4})
        reg_by_name[rec.name] = rec
    ov = pd.DataFrame(overlap_rows)
    ov["ifn_overlap_fdr"] = bh(ov["ifn_overlap_p"].to_numpy())
    ov.to_csv(TAB / "regulon_ifn_overlap.tsv", sep="\t", index=False)

    def family(col: str) -> str:
        if col.startswith("PROGRAM::"):
            return "program_AUCell"
        if col.endswith("|CLDN4out"):
            return "cistarget_CLDN4_held_out"
        return "cistarget_as_returned"

    def meta_for(col: str) -> dict:
        if col.startswith("PROGRAM::"):
            program = col.split("::", 1)[1]
            axis = axis_of_program(program)
            return {"tf": "", "axis": axis, "kind": "program_AUCell_not_cistarget",
                    "elf3_given": False, "nes": np.nan, "n_targets": len(programs_used.get(program, [])),
                    "base_name": program}
        base = col.replace("|CLDN4out", "")
        rec = reg_by_name.get(base)
        tf = rec.tf if rec is not None else base.split("(")[0]
        ov_row = ov.loc[ov["name"] == base]
        ifn_hit = False
        if len(ov_row):
            ifn_hit = bool(ov_row["ifn_overlap"].iloc[0] >= 5 and ov_row["ifn_overlap_fdr"].iloc[0] < 0.05)
        if tf in IFN_TFS:
            axis = "IFN"
        elif tf in MHC_TFS:
            axis = "MHC"
        elif tf in TJ_TFS:
            axis = "TJ"
        elif ifn_hit:
            axis = "IFN_by_targets"
        else:
            axis = "other"
        return {
            "tf": tf,
            "axis": axis,
            "kind": "cistarget_regulon",
            "elf3_given": tf in GIVEN_TFS,
            "nes": float(rec.nes) if rec is not None else np.nan,
            "n_targets": int(rec.n_targets) if rec is not None else np.nan,
            "base_name": base,
        }

    # Score every AUCell column at patient level (primary) and sample level.
    contrast_rows = []
    patient_long = []
    for j, col in enumerate(cols):
        info = meta_for(col)
        # Primary contrast uses CLDN4-held-out regulons and program sets.
        # The as-returned regulon is stored so a sign flip caused by CLDN4 itself is visible.
        pt = patient_table(obs_s, auc[:, j])
        st = wilcoxon_paired(pt["delta"].to_numpy() if len(pt) else np.array([]))
        sd = sample_deltas(obs_s, auc[:, j], sample_info)
        st_s = wilcoxon_paired(sd)
        # tLung
        if len(tlung_labeled):
            # Map tLung labels back is already on original index.
            pt_t = patient_table(tlung_labeled, auc[:, j])
            st_t = wilcoxon_paired(pt_t["delta"].to_numpy() if len(pt_t) else np.array([]))
        else:
            st_t = {"n": 0, "p": np.nan, "delta_median": np.nan, "n_pos": 0, "n_neg": 0}
        expected = {
            "IFN": "lower_in_high",
            "IFN_by_targets": "lower_in_high",
            "MHC": "lower_in_high",
            "TJ": "higher_in_high",
            "KERATIN": "higher_in_high",
            "CTRL": "null",
            "other": "unspecified",
        }[info["axis"]]
        delta = st["delta_median"]
        row = {
            "score": col,
            "tf": info["tf"],
            "axis": info["axis"],
            "kind": info["kind"],
            "family": family(col),
            "elf3_given": info["elf3_given"],
            "nes": info["nes"],
            "n_targets": info["n_targets"],
            "n_patients": st["n"],
            "delta_median_high_minus_low": delta,
            "direction": plain_direction(delta),
            "n_patients_delta_pos": st["n_pos"],
            "n_patients_delta_neg": st["n_neg"],
            "p_patient": st["p"],
            "n_samples": st_s["n"],
            "delta_median_samples": st_s["delta_median"],
            "p_sample": st_s["p"],
            "n_patients_tLung": st_t["n"],
            "delta_median_tLung": st_t["delta_median"],
            "direction_tLung": plain_direction(st_t["delta_median"]),
            "p_tLung": st_t["p"],
            "kd_like_expectation": expected,
            "vs_kd_like": relation(delta, expected),
            "vs_GSE207704": relation(delta, "higher_in_high") if info["axis"] in {"IFN", "IFN_by_targets", "MHC"} else "NA",
            "vs_GSE50927": relation(delta, "lower_in_high") if info["axis"] in {"IFN", "IFN_by_targets", "MHC"} else "NA",
        }
        contrast_rows.append(row)
        if len(pt):
            tmp = pt.copy()
            tmp.insert(0, "score", col)
            patient_long.append(tmp)
    contrast = pd.DataFrame(contrast_rows)
    # FDR within the primary family: cisTarget held-out, and within programs.
    contrast["fdr_patient"] = np.nan
    for fam in ["cistarget_CLDN4_held_out", "program_AUCell"]:
        m = contrast["family"] == fam
        contrast.loc[m, "fdr_patient"] = bh(contrast.loc[m, "p_patient"].to_numpy())
    contrast.to_csv(TAB / "contrasts.tsv", sep="\t", index=False)
    if patient_long:
        pd.concat(patient_long, ignore_index=True).to_csv(TAB / "patient_deltas.tsv.gz", sep="\t", index=False)
    sample_info.to_csv(TAB / "sample_split.tsv", sep="\t", index=False)
    regs.to_csv(TAB / "regulons.tsv", sep="\t", index=False)

    # Split check on CLDN4 itself.
    cldn4_pt = patient_table(obs_s, obs_s["cldn4_log1p_cp10k"].to_numpy(dtype=float))
    cldn4_w = wilcoxon_paired(cldn4_pt["delta"].to_numpy() if len(cldn4_pt) else np.array([]))

    # Between-patient companion: CLDN4+ fraction vs mean AUCell, all QC malignant cells.
    # Different question from the within-sample split.
    between_rows = []
    pat_all = []
    for pid, g in obs.groupby("patient_id"):
        if len(g) < 20:
            continue
        pat_all.append({
            "patient_id": pid,
            "n": int(len(g)),
            "frac_cldn4_pos": float((g["cldn4_count"] > 0).mean()),
            "mean_cldn4": float(g["cldn4_log1p_cp10k"].mean()),
            "idx": g.index.to_numpy(),
        })
    focus_cols = contrast.loc[
        contrast["family"].isin(["cistarget_CLDN4_held_out", "program_AUCell"])
        & contrast["axis"].isin(["IFN", "IFN_by_targets", "MHC", "TJ", "KERATIN", "CTRL"])
    , "score"].tolist()
    # Also the full program list even if axis filter missed one.
    for col in cols:
        if col.startswith("PROGRAM::") and col not in focus_cols:
            focus_cols.append(col)
    for col in focus_cols:
        j = cols.index(col)
        xs, ys = [], []
        for rec in pat_all:
            xs.append(rec["frac_cldn4_pos"])
            ys.append(float(np.nanmean(auc[rec["idx"], j])))
        if len(xs) >= 5:
            rho, p = stats.spearmanr(xs, ys)
        else:
            rho, p = np.nan, np.nan
        between_rows.append({
            "score": col, "n_patients": len(xs), "spearman_rho_frac_cldn4_pos": float(rho),
            "p": float(p), "question": "between_patient_not_the_within_sample_split",
        })
    between = pd.DataFrame(between_rows)
    between.to_csv(TAB / "between_patient_spearman.tsv", sep="\t", index=False)

    primary = contrast.loc[contrast["family"].isin(["cistarget_CLDN4_held_out", "program_AUCell"])].copy()
    # Primary IFN board = STAT/IRF cisTarget regulons + Hallmark/union program AUCell.
    # IFN_by_targets is Hallmark-gene overlap (FOS, JUNB, ATF3, …). It is not an IFN TF.
    ifn = primary.loc[primary["axis"] == "IFN"].sort_values(
        ["kind", "p_patient"], kind="mergesort"
    )
    overlap = primary.loc[primary["axis"] == "IFN_by_targets"].sort_values(
        ["p_patient"], kind="mergesort"
    )
    other_axes = primary.loc[primary["axis"].isin(["MHC", "TJ", "KERATIN", "CTRL"])].sort_values(
        ["axis", "kind", "p_patient"]
    )
    cist = primary.loc[primary["family"] == "cistarget_CLDN4_held_out"].sort_values("p_patient")

    # Figures
    plot_sign_board(
        ifn, FIG / "fig_ifn_sign_board",
        title=("Primary IFN board — observational AUCell, not a knockdown\n"
               "STAT/IRF cisTarget regulons and Hallmark IFN program AUCell only.\n"
               "Red = OPPOSITE the KD-like expectation (IFN lower in CLDN4-high). Navy = SAME."),
    )
    plot_sign_board(
        overlap, FIG / "fig_hallmark_overlap_not_ifn_tf",
        title=("Hallmark-IFN target overlap — these TFs are not interferon regulators\n"
               "Red / navy compare the observational Δ to the IFN KD-like direction only,\n"
               "so an opposite sign is not dropped. The color does not rename the TF."),
    )
    plot_forest(cist, FIG / "fig_cistarget_forest")
    plot_paired(obs_s, auc, cols, FIG / "fig_patient_paired_ifn_tj")
    plot_n(obs, sample_info, cldn4_pt, FIG / "fig_honest_n")

    n_patients_paired = int(cldn4_pt["patient_id"].nunique()) if len(cldn4_pt) else 0
    n_samples_paired = int(sample_info["paired"].sum())
    n_gt0 = int((sample_info["rule"] == "gt0_because_median_0").sum())
    n_samples_ge20 = int((sample_info["n"] >= 20).sum())
    origins = obs.groupby("Sample_Origin").size().to_dict()
    # patients in the paired set
    paired_patients = sorted(cldn4_pt["patient_id"].unique()) if len(cldn4_pt) else []
    n_single = int((cldn4_pt["n_samples"] == 1).sum()) if len(cldn4_pt) and "n_samples" in cldn4_pt else 0
    n_dom = int((cldn4_pt["max_sample_fraction"] > 0.7).sum()) if len(cldn4_pt) else 0

    def sign_line(rec, expectation: str) -> str:
        return (
            f"- **{rec.score}** ({rec.kind}): observational median Δ = {fmt(rec.delta_median_high_minus_low)} "
            f"({rec.direction}); {int(rec.n_patients_delta_pos)}/{int(rec.n_patients)} patients Δ>0; "
            f"Wilcoxon p={fmt(rec.p_patient)}; tLung median Δ = {fmt(rec.delta_median_tLung)} "
            f"({rec.direction_tLung}, p={fmt(rec.p_tLung)}). "
            f"Versus KD-like expectation ({expectation}): **{rec.vs_kd_like}**. "
            f"Versus GSE207704 (IFN/MHC down after CLDN4 loss, predicts Δ>0): **{rec.vs_GSE207704}**. "
            f"Versus GSE50927 (IFN/MHC up after Cldn4 KO, n=1 lung, predicts Δ<0): **{rec.vs_GSE50927}**."
        )

    # Headline counts use the primary IFN board only.
    n_opp = int((ifn["vs_kd_like"] == "OPPOSITE").sum()) if len(ifn) else 0
    n_same = int((ifn["vs_kd_like"] == "SAME").sum()) if len(ifn) else 0
    n_ov_opp = int((overlap["vs_kd_like"] == "OPPOSITE").sum()) if len(overlap) else 0
    n_ov_same = int((overlap["vs_kd_like"] == "SAME").sum()) if len(overlap) else 0
    ifn_lines = [sign_line(rec, "IFN lower in CLDN4-high") for rec in ifn.itertuples(index=False)]
    if not ifn_lines:
        ifn_lines.append(
            "- No Hallmark-IFN program and no STAT/IRF cisTarget regulon was scored. "
            "That absence is reported here rather than filled in."
        )
    overlap_lines = [
        sign_line(rec, "IFN lower in CLDN4-high, bookkeeping only — this TF is not an IFN regulator")
        for rec in overlap.itertuples(index=False)
    ]
    if not overlap_lines:
        overlap_lines.append("- No cisTarget regulon met the Hallmark-IFN overlap rule (k≥5 and FDR<0.05).")

    cist_ifn = ifn.loc[ifn["kind"] == "cistarget_regulon"] if len(ifn) else ifn
    prog_ifn = ifn.loc[ifn["kind"] == "program_AUCell_not_cistarget"] if len(ifn) else ifn

    def row_by_score(frame: pd.DataFrame, score: str):
        hit = frame.loc[frame["score"] == score]
        return None if hit.empty else hit.iloc[0]

    def brief(rec) -> str:
        if rec is None:
            return "was not scored"
        return (
            f"is {rec.direction} (median Δ={fmt(rec.delta_median_high_minus_low)}, "
            f"{int(rec.n_patients_delta_pos)}/{int(rec.n_patients)} patients Δ>0, "
            f"p={fmt(rec.p_patient)}; tLung Δ={fmt(rec.delta_median_tLung)}, "
            f"{rec.direction_tLung}, p={fmt(rec.p_tLung)}); "
            f"**{rec.vs_kd_like}** versus the KD-like expectation, "
            f"**{rec.vs_GSE207704}** versus GSE207704, **{rec.vs_GSE50927}** versus GSE50927"
        )

    def between_of(score: str) -> str:
        hit = between.loc[between["score"] == score]
        if hit.empty:
            return "not computed"
        rec = hit.iloc[0]
        return (
            f"ρ={fmt(rec.spearman_rho_frac_cldn4_pos)}, p={fmt(rec.p)}, "
            f"n={int(rec.n_patients)} patients"
        )

    def tf_bits(frame: pd.DataFrame) -> str:
        if frame is None or len(frame) == 0:
            return "none recovered at NES ≥ 3.0"
        parts = []
        for rec in frame.sort_values("p_patient").itertuples(index=False):
            given = " (A10-given, not a discovery)" if bool(rec.elf3_given) else ""
            flip = ""
            if rec.direction != rec.direction_tLung and rec.direction_tLung not in {"undefined", "NA"}:
                flip = (
                    f"; all-site direction and tLung direction disagree "
                    f"(tLung {rec.direction_tLung}, p={fmt(rec.p_tLung)})"
                )
            parts.append(
                f"{rec.tf}{given} {rec.direction} (Δ={fmt(rec.delta_median_high_minus_low)}, "
                f"p={fmt(rec.p_patient)}, **{rec.vs_kd_like}** KD-like{flip})"
            )
        return "; ".join(parts)

    mhc_rows = other_axes.loc[other_axes["axis"] == "MHC"] if len(other_axes) else other_axes
    tj_rows = other_axes.loc[other_axes["axis"] == "TJ"] if len(other_axes) else other_axes
    ker_rows = other_axes.loc[other_axes["axis"] == "KERATIN"] if len(other_axes) else other_axes
    ctrl_rows = other_axes.loc[other_axes["axis"] == "CTRL"] if len(other_axes) else other_axes
    mhc_prog = row_by_score(other_axes, "PROGRAM::MHC1_APM")
    tj_prog = row_by_score(other_axes, "PROGRAM::TJ_STRUCT")
    ker_prog = row_by_score(other_axes, "PROGRAM::KERATIN")
    ribo = row_by_score(other_axes, "PROGRAM::CTRL_RIBO")
    ifn_union = row_by_score(ifn, "PROGRAM::IFN_UNION")
    ifng = row_by_score(ifn, "PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE")
    ifna = row_by_score(ifn, "PROGRAM::HALLMARK_INTERFERON_ALPHA_RESPONSE")

    # Held-out vs as-returned sign, primary axes only.
    held = contrast.loc[contrast["family"] == "cistarget_CLDN4_held_out",
                        ["score", "tf", "axis", "direction", "delta_median_high_minus_low"]].copy()
    held["base"] = held["score"].str.replace("|CLDN4out", "", regex=False)
    raw = contrast.loc[contrast["family"] == "cistarget_as_returned",
                       ["score", "direction", "delta_median_high_minus_low"]].copy()
    raw = raw.rename(columns={"score": "base", "direction": "direction_raw",
                              "delta_median_high_minus_low": "delta_raw"})
    merged_sign = held.merge(raw, on="base", how="left")
    merged_sign = merged_sign.loc[merged_sign["axis"].isin(["IFN", "MHC", "TJ"])]
    sign_flips = merged_sign.loc[merged_sign["direction"] != merged_sign["direction_raw"]]

    def opening() -> str:
        paras = []
        paras.append(
            "The primary IFN board is STAT/IRF cisTarget regulons plus Hallmark IFN program AUCell "
            "(IFN-α, IFN-γ, and their union). Regulons that only overlap Hallmark IFN genes are "
            "**not** on this board and are **not** counted in the SAME / OPPOSITE totals below. "
            "They are printed in full in the next subsection so those signs are not dropped."
        )
        paras.append(
            f"**PROGRAM::IFN_UNION** {brief(ifn_union)}. "
            f"**PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE** {brief(ifng)}. "
            f"**PROGRAM::HALLMARK_INTERFERON_ALPHA_RESPONSE** {brief(ifna)}. "
            "KD-like expectation for these rows: IFN **lower** in CLDN4-high (Δ < 0), the frame in which "
            "CLDN4 loss would raise IFN. That expectation was **not** fit on this matrix. "
            "GSE207704 (public CLDN4 CRISPR in T47D/MCF7) has IFN trending down after loss and predicts Δ > 0. "
            "GSE50927 (public whole-lung Cldn4 KO, n=1, not a cancer-cell test) has IFN/MHC up after loss "
            "and predicts Δ < 0. A row can be SAME versus one reference and OPPOSITE versus the other."
        )
        if len(cist_ifn):
            same = ", ".join(cist_ifn.loc[cist_ifn["vs_kd_like"] == "SAME", "tf"].astype(str)) or "none"
            opp = ", ".join(cist_ifn.loc[cist_ifn["vs_kd_like"] == "OPPOSITE", "tf"].astype(str)) or "none"
            paras.append(
                f"cisTarget STAT/IRF regulons at NES ≥ 3.0, CLDN4 held out ({len(cist_ifn)} regulons). "
                f"Observational sign **OPPOSITE** the KD-like expectation (higher in CLDN4-high): {opp}. "
                f"Observational sign **SAME** as the KD-like expectation (higher in CLDN4-low): {same}. "
                "Both lists are on the board. A non-significant p does not remove the sign."
            )
        else:
            paras.append("cisTarget at NES ≥ 3.0 recovered no STAT/IRF regulon.")
        paras.append(
            f"**PROGRAM::MHC1_APM** {brief(mhc_prog)}. "
            f"KD-like expectation: MHC lower in CLDN4-high. "
            f"cisTarget MHC TFs: {tf_bits(mhc_rows.loc[mhc_rows['kind'] == 'cistarget_regulon'])}. "
            "If the MHC program and a cisTarget MHC regulon disagree, both stay in this paragraph and in the MHC table."
        )
        paras.append(
            f"**PROGRAM::TJ_STRUCT** {brief(tj_prog)}. "
            f"**PROGRAM::KERATIN** {brief(ker_prog)}. "
            "KD-like expectation: TJ and keratin higher in CLDN4-high. "
            f"cisTarget TJ-associated TFs: {tf_bits(tj_rows.loc[tj_rows['kind'] == 'cistarget_regulon'])}. "
            f"**PROGRAM::CTRL_RIBO** {brief(ribo)} (null control; SAME/OPPOSITE is not applied)."
        )
        paras.append(
            f"Between-patient Spearman of malignant CLDN4+ fraction versus mean AUCell "
            f"(a different question from the within-sample split): "
            f"IFN_UNION {between_of('PROGRAM::IFN_UNION')}; "
            f"STAT1 {between_of('STAT1(+)|CLDN4out')}; "
            f"MHC1_APM {between_of('PROGRAM::MHC1_APM')}; "
            f"RFXANK {between_of('RFXANK(+)|CLDN4out')}; "
            f"TJ_STRUCT {between_of('PROGRAM::TJ_STRUCT')}. "
            "A null between-patient IFN correlation does not replace the within-sample patient Δ, "
            "and the within-sample Δ does not replace the between-patient correlation."
        )
        if len(sign_flips):
            flip_txt = "; ".join(
                f"{r.tf} held-out {r.direction} vs as-returned {r.direction_raw}"
                for r in sign_flips.itertuples(index=False)
            )
            paras.append(
                "Removing CLDN4 from the regulon changed the observational direction for: " + flip_txt + "."
            )
        else:
            paras.append(
                "Removing CLDN4 from IFN, MHC, and TJ cisTarget regulons did not change any observational direction "
                "relative to the as-returned regulon."
            )
        return "\n\n".join(paras)

    # Regulon inventory for the methods table.
    reg_show = regs.copy()
    reg_show["targets_short"] = reg_show["targets"].fillna("").map(
        lambda s: ";".join(str(s).split(";")[:12]) + ("…" if str(s).count(";") >= 12 else "")
    )

    n_table = pd.DataFrame([
        {"item": "cells_in_matrix", "n": prep["n_cells_matrix"], "note": "GSE131907 UMI header"},
        {"item": "author_malignant", "n": prep["n_author_malignant"],
         "note": "Cell_subtype in {Malignant cells, tS1, tS2, tS3}"},
        {"item": "malignant_UMI_ge_200", "n": prep["n_malignant_umi_ge_200"], "note": "QC cells scored by AUCell"},
        {"item": "genes_aucell_universe", "n": prep["n_genes_universe"], "note": "detected in ≥1% of QC malignant cells"},
        {"item": "grn_cells_subsample", "n": prep["n_grn_cells"], "note": f"stratified by sample, seed {prep['seed']}"},
        {"item": "grn_genes", "n": prep["n_grn_genes"], "note": "HVG + TFs detected in ≥5% + program genes"},
        {"item": "tf_regulators_grnboost2", "n": prep["n_tf_regulators"], "note": "Lambert/Aerts TF list, detection ≥0.05"},
        {"item": "cistarget_regulons_NES3", "n": int(len(regs)), "note": "activating modules kept by pySCENIC defaults"},
        {"item": "samples_ge20_malignant", "n": n_samples_ge20, "note": "before the high/low tail rule"},
        {"item": "samples_paired_ge20_each_tail", "n": n_samples_paired, "note": "unit inside each patient"},
        {"item": "samples_split_by_gt0", "n": n_gt0, "note": "median CLDN4 log1p was 0"},
        {"item": "patients_paired", "n": n_patients_paired, "note": "primary Wilcoxon unit"},
        {"item": "paired_patients_single_sample", "n": n_single, "note": "one sample supplied both tails"},
        {"item": "paired_patients_one_sample_gt_70pct_cells", "n": n_dom, "note": "cell-weighted mean can follow that sample"},
        {"item": "patients_tLung_paired", "n": int(contrast["n_patients_tLung"].max() if len(contrast) else 0),
         "note": "primary-tumor cells only; split recomputed"},
    ])
    n_table.to_csv(TAB / "honest_n.tsv", sep="\t", index=False)

    show_cols = ["score", "tf", "axis", "kind", "n_patients", "delta_median_high_minus_low",
                 "direction", "n_patients_delta_pos", "n_patients_delta_neg", "p_patient", "fdr_patient",
                 "vs_kd_like", "vs_GSE207704", "vs_GSE50927",
                 "delta_median_tLung", "direction_tLung", "p_tLung"]
    reg_cols = ["name", "tf", "nes", "n_targets", "contains_CLDN4"]

    text = f"""# RESULTS — pySCENIC on GSE131907 malignant cells, CLDN4-high vs low

Cohort: **GSE131907** (Kim et al., *Nat Commun* 2020, PMID 32385277), author-malignant cells.
Method that was actually run: **pySCENIC 0.12.1** (GRNBoost2 → cisTarget motif pruning → AUCell).
SCENIC+ was not run: this GEO deposit is scRNA-seq only (no matched scATAC).
Concordant-4 was not co-embedded. This machine has 15 GiB RAM; cisTarget rankings plus a four-cohort malignant matrix do not fit. GSE131907 is the largest single member of that concordant set.

This file separates three things that are easy to mix:

1. **Observational** regulon / program AUCell on tumor cells (this run).
2. **KD-like expectation** used by earlier AUCell work in this project: IFN and MHC activity **lower** in CLDN4-high cells (Δ high−low < 0), TJ **higher** in CLDN4-high. That is the direction in which CLDN4 loss would raise IFN. It is **not** estimated from these cells.
3. **Public CLDN4-loss transcriptomes**, which do not agree with each other. **GSE207704** (CLDN4 CRISPR in T47D/MCF7) has IFN trending **down** after loss, which predicts observational Δ > 0. **GSE50927** (whole-lung Cldn4 KO, n=1) has IFN/MHC **up** after loss (NES +1.51 in the earlier public summary), which predicts observational Δ < 0. GSE50927 is not a cancer-cell replicate.

A10 ELF3–CLDN4 is taken as given. ELF3 rows are labeled `elf3_given` and are not a new discovery.

## Observational IFN sign (read this before the rest)

{opening()}

Primary IFN board, observational sign **OPPOSITE** the KD-like expectation: **{n_opp}**. Primary IFN board, **SAME** sign as that expectation: **{n_same}**. These two counts are only the STAT/IRF and Hallmark/union rows below. They do not include Hallmark-overlap regulons.

{chr(10).join(ifn_lines)}

The KD-like column is the pre-specified “IFN lower in CLDN4-high” direction. The GSE207704 column is often the opposite reference, because that public CRISPR series did not open IFN. A row can be SAME versus one reference and OPPOSITE versus the other. That disagreement is the result, not a reason to drop a row.

Program AUCell rows are **not** cisTarget regulons. They use the same recovery curve (top 5% of the detected-gene ranking) on Hallmark IFN sets with CLDN4 removed. cisTarget rows on this board are STAT1, STAT2, or IRF1–IRF9 regulons that passed motif enrichment at NES ≥ 3.0, rescored after removing CLDN4 so the split gene cannot sit inside the signature.

### Primary IFN table

{md_table(ifn, show_cols) if len(ifn) else "_No IFN rows._"}

### Hallmark-IFN target overlap (not IFN transcription factors)

These cisTarget regulons are **not** STAT or IRF. Each has ≥5 targets in the Hallmark IFN union and a hypergeometric FDR < 0.05. The SAME / OPPOSITE labels use the IFN KD-like direction so an opposite observational sign is still printed. That column does not rename FOS, JUNB, ATF3, or any other TF here as an interferon regulator.

Overlap rows if the IFN KD-like direction is used as the bookkeeping reference: **{n_ov_opp}** OPPOSITE, **{n_ov_same}** SAME. These counts are **not** added to the primary IFN totals above.

{chr(10).join(overlap_lines)}

{md_table(overlap, show_cols) if len(overlap) else "_No overlap rows._"}

### MHC, TJ, keratin, ribosome control

KD-like expectation: MHC lower in CLDN4-high; TJ and keratin higher in CLDN4-high; ribosome null. GSE207704 / GSE50927 columns are filled only for MHC (same public IFN/MHC notes as above). A cisTarget row that disagrees with its program is left in this table.

{md_table(other_axes, show_cols) if len(other_axes) else "_No comparator rows._"}

## Split check

Within-sample CLDN4 log1p(CP10k), same patient aggregate as the regulon test: median Δ = {fmt(cldn4_w["delta_median"])}, {cldn4_w["n_pos"]}/{cldn4_w["n"]} patients positive, Wilcoxon p={fmt(cldn4_w["p"])}. The high group is the upper tail of CLDN4 by construction. This row only checks that the labels were applied.

## Honest n

Primary unit = **patient**. Cells from a sample enter that patient only when the sample itself has ≥20 CLDN4-high and ≥20 CLDN4-low malignant cells. High/low is the within-sample median of CLDN4 log1p(CP10k). If that median is 0, high = CLDN4 > 0. Patient means are cell-weighted across that patient's paired samples. Cell counts are inventory, not the test.

{md_table(n_table, ["item", "n", "note"])}

Author-malignant cells by `Sample_Origin`: {", ".join(f"{k} {v}" for k, v in sorted(origins.items(), key=lambda kv: -kv[1]))}.

mBrain can dominate a pooled cell test. The patient aggregate stops one metastasis from being hundreds of pseudo-replicates, but a patient whose cells are mostly one sample still follows that sample ({n_dom} paired patients have one sample with >70% of the paired cells).

Paired patients: {", ".join(paired_patients) if paired_patients else "none"}.

## cisTarget regulons (CLDN4 held out)

GRNBoost2 was fit on a stratified subsample of {prep["n_grn_cells"]} malignant cells and {prep["n_grn_genes"]} genes (log1p CP10k; library size from all genes; seed {prep["seed"]}). The subsample was not restricted to CLDN4-high cells. cisTarget used hg38 v10 cluster rankings for 500 bp upstream / 100 bp downstream and for ±10 kb, motif table `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`, rank threshold 1500, NES ≥ 3.0, pySCENIC default activating modules. AUCell was then run on **all** QC malignant cells, ranking genes detected in ≥1% of those cells. Within a cell, raw-count ranks match log1p(CP10k) ranks.

{int((regs["contains_CLDN4"] == True).sum()) if "contains_CLDN4" in regs else "NA"} regulons contained CLDN4 before it was removed for the contrast. Full target lists: `tables/regulons.tsv`.

{md_table(cist, ["score", "tf", "axis", "nes", "n_targets", "elf3_given", "n_patients", "delta_median_high_minus_low", "direction", "n_patients_delta_pos", "p_patient", "fdr_patient", "vs_kd_like", "delta_median_samples", "p_sample", "delta_median_tLung", "direction_tLung", "p_tLung"]) if len(cist) else "_No cisTarget regulon._"}

FDR is Benjamini–Hochberg within the CLDN4-held-out cisTarget family. It is a descriptive multiplicity adjustment, not a binding claim.

### As-returned regulons (CLDN4 still inside, if it was a target)

Kept so a sign that appears only when CLDN4 is inside the signature is visible. Not the primary contrast.

{md_table(contrast.loc[contrast["family"] == "cistarget_as_returned", ["score", "tf", "n_patients", "delta_median_high_minus_low", "direction", "p_patient", "vs_kd_like"]], ["score", "tf", "n_patients", "delta_median_high_minus_low", "direction", "p_patient", "vs_kd_like"])}

## tLung sensitivity

Same rules, restricted to `Sample_Origin == tLung` before the median split. Columns `delta_median_tLung`, `direction_tLung`, and `p_tLung` are in the tables above. A sign that flips between the all-site patient test and tLung is a real disagreement and is not averaged away.

## Between-patient companion (not the high vs low split)

Spearman of the patient's malignant CLDN4+ fraction versus the mean AUCell of all that patient's malignant cells. This is the “CLDN4-high patients” question. It is not the within-sample high versus low contrast, and the two can disagree.

{md_table(between, ["score", "n_patients", "spearman_rho_frac_cldn4_pos", "p", "question"]) if len(between) else "_No between-patient rows._"}

## What this is / is not

- **Is** a pySCENIC run (GRNBoost2 + cisTarget + AUCell) on GSE131907 author-malignant cells, with the CLDN4-high versus low contrast at patient level.
- **Is** an explicit report of observational IFN sign against the KD-like expectation and against GSE207704 and GSE50927, including opposite signs.
- **Is not** SCENIC+. There is no scATAC in this accession.
- **Is not** a concordant-4 integrated GRN. Memory stopped that merge. The other three cohorts were not scored here.
- **Is not** a knockdown, a coculture, or an estimate of the private CLDN4 KD. Those data are not in this public repository and were not re-fit.
- **Is not** ChIP, and cisTarget motif recovery is not proof that the TF binds the target in these tumors.
- **Is not** an ICI or MPR result. GSE131907 is treatment-naive.
- **Is not** a new ELF3–CLDN4 discovery.
- **Is not** a cell-level p-value. EBUS and brain-metastasis samples can dominate cell counts; the patient is the unit.
- numpy shim: `np.object = object` before importing pySCENIC 0.12.1, because NumPy 2 removed that alias. No GRN equation was changed.

## Figures

- `figures/fig_ifn_sign_board.png` — primary IFN board only (STAT/IRF cisTarget + Hallmark program AUCell), colored by SAME versus OPPOSITE the KD-like expectation
- `figures/fig_hallmark_overlap_not_ifn_tf.png` — regulons with Hallmark-IFN target overlap; not IFN transcription factors; signs kept
- `figures/fig_cistarget_forest.png` — cisTarget regulons, CLDN4 held out; red = STAT/IRF, gold = MHC TF, green = TJ-associated TF
- `figures/fig_patient_paired_ifn_tj.png` — per-patient high versus low for the IFN union and the TJ program
- `figures/fig_honest_n.png` — cell, sample, and patient counts

## Reproduce

```bash
python3 methods/gse131907_pyscenic_cldn4/scripts/01_prepare.py
python3 methods/gse131907_pyscenic_cldn4/scripts/02_pyscenic.py
python3 methods/gse131907_pyscenic_cldn4/scripts/03_contrast.py
```

Downloads (not committed): GEO UMI + annotation + series matrix; Aerts lab `allTFs_hg38.txt`, `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`, and the two hg38 v10 genes-vs-motifs ranking feathers; MSigDB Hallmark GMT `h.all.v2023.2.Hs.symbols.gmt`.

Run flags: cistarget_run={run.get("cistarget_run", True)}; scenicplus_run={run.get("scenicplus_run", False)}; pyscenic={run.get("pyscenic", "0.12.1")}.
"""
    (RES / "RESULTS.md").write_text(text)
    # Compact machine-readable headline for the PR. Numbers only, already in the tables.
    ifn_cols = ["score", "kind", "delta_median_high_minus_low", "direction", "p_patient",
                "vs_kd_like", "vs_GSE207704", "vs_GSE50927"]
    headline = {
        "n_patients_paired": n_patients_paired,
        "n_samples_paired": n_samples_paired,
        "n_regulons": int(len(regs)),
        "primary_ifn_board": "STAT_IRF_cistarget_plus_hallmark_program_AUCell",
        "ifn_rows": ifn[ifn_cols].to_dict(orient="records") if len(ifn) else [],
        "n_ifn_opposite_kd_like": n_opp,
        "n_ifn_same_kd_like": n_same,
        "hallmark_overlap_not_ifn_tf_n": int(len(overlap)),
        "hallmark_overlap_opposite_if_compared_to_ifn_kd_like": n_ov_opp,
        "hallmark_overlap_same_if_compared_to_ifn_kd_like": n_ov_same,
        "note": "Overlap counts are not IFN transcription factors and are not added to n_ifn_*.",
    }
    (TAB / "headline.json").write_text(json.dumps(headline, indent=2) + "\n")
    print(json.dumps(headline, indent=2), flush=True)
    print(f"[wrote] {RES / 'RESULTS.md'}", flush=True)


def plot_sign_board(ifn: pd.DataFrame, stem: Path, title: str) -> None:
    fig, ax = plt.subplots(figsize=(11, max(3.2, 0.42 * max(len(ifn), 1) + 1.8)))
    if ifn.empty:
        ax.text(0.5, 0.5, "No rows", ha="center")
    else:
        plot = ifn.iloc[::-1]
        y = np.arange(len(plot))
        colors = []
        for flag in plot["vs_kd_like"]:
            if flag == "OPPOSITE":
                colors.append("#9B2226")
            elif flag == "SAME":
                colors.append("#1D3557")
            else:
                colors.append("#6C757D")
        ax.barh(y, plot["delta_median_high_minus_low"], color=colors, edgecolor="none")
        ax.axvline(0, color="black", lw=0.8)
        labels = []
        for rec in plot.itertuples(index=False):
            kind = "program" if str(rec.kind).startswith("program") else "cisTarget"
            labels.append(f"{rec.score}  [{kind}]  {rec.vs_kd_like} vs KD-like")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel("Patient median Δ AUCell (CLDN4-high − CLDN4-low)")
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=160)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_forest(cist: pd.DataFrame, stem: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, max(3.5, 0.32 * max(len(cist), 1) + 1.2)))
    if cist.empty:
        ax.text(0.5, 0.5, "No cisTarget regulons", ha="center")
    else:
        # Show all if <= 40, else the 30 smallest p plus every IFN/MHC/TJ row.
        if len(cist) > 40:
            keep_axis = cist["axis"].isin(["IFN", "MHC", "TJ"])
            top = cist.nsmallest(30, "p_patient")
            plot = pd.concat([cist.loc[keep_axis], top]).drop_duplicates("score")
            plot = plot.sort_values("delta_median_high_minus_low")
        else:
            plot = cist.sort_values("delta_median_high_minus_low")
        y = np.arange(len(plot))
        colors = []
        for a in plot["axis"]:
            if a == "IFN":
                colors.append("#9B2226")
            elif a == "MHC":
                colors.append("#BB3E03")
            elif a == "TJ":
                colors.append("#2A9D8F")
            else:
                colors.append("#495057")
        ax.barh(y, plot["delta_median_high_minus_low"], color=colors)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.tf} {r.axis}" for r in plot.itertuples(index=False)], fontsize=7)
        ax.set_xlabel("Patient median Δ (CLDN4 held out of the regulon)")
    ax.set_title("cisTarget regulons, observational CLDN4-high − low\nRed = STAT/IRF, rust = MHC TF, green = TJ-associated TF, grey = other")
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=160)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_paired(obs: pd.DataFrame, auc: np.ndarray, cols: list[str], stem: Path) -> None:
    want = []
    for name in ["PROGRAM::IFN_UNION", "PROGRAM::HALLMARK_INTERFERON_GAMMA_RESPONSE", "PROGRAM::TJ_STRUCT", "PROGRAM::MHC1_APM"]:
        if name in cols:
            want.append(name)
    if not want:
        return
    fig, axes = plt.subplots(1, len(want), figsize=(3.4 * len(want), 4.2), sharey=False)
    if len(want) == 1:
        axes = [axes]
    rng = np.random.default_rng(0)
    for ax, name in zip(axes, want):
        j = cols.index(name)
        pt = patient_table(obs, auc[:, j])
        if pt.empty:
            ax.set_title(name)
            continue
        for rec in pt.itertuples(index=False):
            jitter = float(rng.normal(0, 0.03))
            ax.plot([0 + jitter, 1 + jitter], [rec.mean_low, rec.mean_high], color="#adb5bd", lw=0.7)
        ax.scatter(np.zeros(len(pt)), pt["mean_low"], s=16, c="#1D3557", zorder=3)
        ax.scatter(np.ones(len(pt)), pt["mean_high"], s=16, c="#9B2226", zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["CLDN4-low", "CLDN4-high"])
        short = name.replace("PROGRAM::", "")
        ax.set_title(short, fontsize=9)
        ax.set_ylabel("Patient-mean AUCell")
    fig.suptitle("Observational patient means. Not a knockdown.", fontsize=11)
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=160)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


def plot_n(obs: pd.DataFrame, sample_info: pd.DataFrame, cldn4_pt: pd.DataFrame, stem: Path) -> None:
    labels = ["QC malignant\ncells", "samples ≥20", "samples paired\n≥20 / tail", "patients paired"]
    vals = [
        len(obs),
        int((sample_info["n"] >= 20).sum()),
        int(sample_info["paired"].sum()),
        int(cldn4_pt["patient_id"].nunique()) if len(cldn4_pt) else 0,
    ]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.bar(labels, vals, color=["#ced4da", "#adb5bd", "#495057", "#1D3557"])
    for i, v in enumerate(vals):
        ax.text(i, v, str(v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Count")
    ax.set_title("Honest n — the test is the patient bar")
    fig.tight_layout()
    fig.savefig(stem.with_suffix(".png"), dpi=160)
    fig.savefig(stem.with_suffix(".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
