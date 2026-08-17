#!/usr/bin/env python3
"""Build the patient-level CLDN4 + T/NK/B composition table.

Inputs are already-extracted public tables (PR #283 / #284 / #290 / #318).
No new GEO matrices are downloaded. CLDN4 only — TACSTD2 is never a gate.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "data" / "extracted"
OUT = ROOT / "results"
OUT.mkdir(parents=True, exist_ok=True)


def _median_split(score: pd.Series) -> tuple[pd.Series, float]:
    s = pd.to_numeric(score, errors="coerce")
    cut = float(s.median())
    high = pd.Series(np.nan, index=s.index, dtype=float)
    ok = s.notna()
    high.loc[ok] = (s.loc[ok] >= cut).astype(float)
    return high, cut


def gse207422() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Hu marker lineages + A3-malignant CLDN4 (PR #318). Post-tx primary."""
    lin = pd.read_csv(EXT / "GSE207422_lineage_counts.tsv", sep="\t")
    pp = pd.read_csv(EXT / "GSE207422_per_patient.tsv", sep="\t")
    wide = lin.pivot_table(index="Sample", columns="lineage", values="n_cells", aggfunc="sum", fill_value=0)
    for col in ["T", "NK", "B", "plasma", "epithelial", "myeloid", "neutrophil", "fibroblast", "endothelial", "mast", "other"]:
        if col not in wide.columns:
            wide[col] = 0
    wide["TNK"] = wide["T"] + wide["NK"]
    wide["TLS"] = wide["B"] + wide["plasma"]
    wide["Stromal"] = wide["fibroblast"] + wide["endothelial"]
    wide["Myeloid"] = wide["myeloid"] + wide["neutrophil"] + wide["mast"]
    wide["Other"] = wide["other"]
    wide = wide.reset_index().rename(columns={"Sample": "sample_id"})
    meta = pp.rename(columns={"Sample": "sample_id"})
    df = wide.merge(meta, on="sample_id", how="left")
    df["cohort"] = "GSE207422"
    df["subcohort"] = "Hu2023"
    df["annotation"] = "hu_markers"
    df["cldn4_score"] = df["mal_cldn4_mean_log1p_cp10k"]
    df["cldn4_n"] = df["n_malig_a3"]
    df["cldn4_definition"] = "A3-malignant mean log1p(CP10k); empty if 0 A3-malignant cells"
    df["n_malignant"] = df["n_malig_a3"]
    df["n_epithelial"] = df["n_epithelial"]
    df["n_T"] = df["T"]
    df["n_NK"] = df["NK"]
    df["n_TNK"] = df["TNK"]
    df["n_B"] = df["B"]
    df["n_plasma"] = df["plasma"]
    df["n_B_plasma"] = df["TLS"]
    df["n_cells_comp"] = df[["epithelial", "TNK", "B", "Myeloid", "Stromal", "Other"]].sum(axis=1)
    df["response"] = df["path_response"]
    df["patient_id"] = df["patient"]
    # primary: post-treatment, finite A3-malignant CLDN4
    df["eligible"] = (df["timing"] == "post") & df["cldn4_score"].notna() & (df["n_malignant"] >= 1)
    # P06 has 1 A3-malignant cell with CLDN4=0 — kept, noisy (same as PR #318)
    high, cut = _median_split(df.loc[df["eligible"], "cldn4_score"])
    df["cldn4_high"] = np.nan
    df.loc[high.index, "cldn4_high"] = high
    df["cldn4_median_cut"] = cut
    df["slice"] = "GSE207422_post_A3mal"
    df["primary"] = df["eligible"]
    df["provenance"] = "PR318 lineage_counts + per_patient (A3-malignant CLDN4)"
    long_parts = []
    for ct, src in [
        ("Epithelial", "epithelial"),
        ("TNK", "TNK"),
        ("B", "B"),
        ("Myeloid", "Myeloid"),
        ("Stromal", "Stromal"),
        ("Other", "Other"),
    ]:
        tmp = df[["sample_id", "patient_id", "cohort", "slice", src]].copy()
        tmp["celltype"] = ct
        tmp["n"] = tmp[src]
        long_parts.append(tmp[["cohort", "slice", "sample_id", "patient_id", "celltype", "n"]])
    return df, pd.concat(long_parts, ignore_index=True)


def gse207422_drmref() -> tuple[pd.DataFrame, pd.DataFrame]:
    """DRMref public labels (PR #283 counts + PR #290 CLDN4). Sensitivity, not primary."""
    long = pd.read_csv(EXT / "composition_long_pr283.tsv", sep="\t")
    sub = long[(long["cohort"] == "GSE207422") & (long["annotation"] == "drmref_collapsed")].copy()
    wide = sub.pivot_table(index="sample_id", columns="celltype", values="n", aggfunc="sum", fill_value=0)
    # expected: TLS, TNK, stromal, epithelial, myeloid, mast
    for col in ["TLS", "TNK", "stromal", "epithelial", "myeloid", "mast"]:
        if col not in wide.columns:
            wide[col] = 0
    wide = wide.reset_index()
    drm = pd.read_csv(EXT / "GSE207422_drmref_patients.tsv", sep="\t")
    drm = drm.rename(columns={"sample": "sample_id"})
    df = wide.merge(drm, on="sample_id", how="left")
    pp = pd.read_csv(EXT / "GSE207422_per_patient.tsv", sep="\t").rename(columns={"Sample": "sample_id"})
    df = df.merge(pp[["sample_id", "timing", "path_response", "patient"]], on="sample_id", how="left")
    df["cohort"] = "GSE207422"
    df["subcohort"] = "DRMref"
    df["annotation"] = "drmref_collapsed"
    df["cldn4_score"] = df["malig_CLDN4_mean"]
    df["cldn4_n"] = df["n_malignant"]
    df["cldn4_definition"] = "DRMref malignant mean (PR290); not A3"
    df["n_T"] = np.nan
    df["n_NK"] = np.nan
    df["n_TNK"] = df["TNK"]
    df["n_B"] = np.nan
    df["n_plasma"] = np.nan
    df["n_B_plasma"] = df["TLS"]
    df["n_epithelial"] = df["epithelial"]
    df["n_cells_comp"] = df[["epithelial", "TNK", "TLS", "myeloid", "stromal", "mast"]].sum(axis=1)
    df["response"] = df["response"].fillna(df["path_response"])
    df["patient_id"] = df["patient"]
    df["eligible"] = df["cldn4_score"].notna() & (df["n_malignant"] >= 1)
    high, cut = _median_split(df.loc[df["eligible"], "cldn4_score"])
    df["cldn4_high"] = np.nan
    df.loc[high.index, "cldn4_high"] = high
    df["cldn4_median_cut"] = cut
    df["slice"] = "GSE207422_drmref"
    df["primary"] = False
    df["provenance"] = "PR283 drmref_collapsed + PR290 drmref CLDN4"
    df["B"] = df["TLS"]  # TLS is the B-lineage part in this collapse
    df["Myeloid"] = df["myeloid"]
    df["Stromal"] = df["stromal"]
    df["Other"] = df["mast"]
    df["Epithelial"] = df["epithelial"]
    long_parts = []
    for ct, src in [("Epithelial", "epithelial"), ("TNK", "TNK"), ("B", "TLS"), ("Myeloid", "myeloid"), ("Stromal", "stromal"), ("Other", "mast")]:
        tmp = df[["sample_id", "patient_id", "cohort", "slice", src]].copy()
        tmp["celltype"] = ct
        tmp["n"] = tmp[src]
        long_parts.append(tmp[["cohort", "slice", "sample_id", "patient_id", "celltype", "n"]])
    return df, pd.concat(long_parts, ignore_index=True)


def _gse241934(which: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Author major.cell.type + epithelial CLDN4 (PR #283 + #290)."""
    long = pd.read_csv(EXT / "composition_long_pr283.tsv", sep="\t")
    if which == "IIT":
        cohort_key = "GSE241934_IIT"
        scores = pd.read_csv(EXT / "GSE241934_IIT_patients.tsv", sep="\t")
        slice_name = "GSE241934_IIT"
        annotation = "author_major"
    else:
        cohort_key = "GSE241934_RWC"
        scores = pd.read_csv(EXT / "GSE241934_Real_patients.tsv", sep="\t")
        slice_name = "GSE241934_RWC"
        annotation = "author_major"
    sub = long[(long["cohort"] == cohort_key) & (long["annotation"] == "author_major")].copy()
    wide = sub.pivot_table(index="sample_id", columns="celltype", values="n", aggfunc="sum", fill_value=0)
    for col in ["T", "NK", "B", "epithelial", "myeloid", "fibroblast", "endothelial", "mast"]:
        if col not in wide.columns:
            wide[col] = 0
    wide["TNK"] = wide["T"] + wide["NK"]
    wide["Stromal"] = wide["fibroblast"] + wide["endothelial"]
    wide["Myeloid"] = wide["myeloid"] + wide["mast"]
    wide["Other"] = 0
    wide = wide.reset_index()
    scores = scores.rename(columns={"patient": "sample_id"})
    df = wide.merge(scores, on="sample_id", how="left")
    df["cohort"] = "GSE241934"
    df["subcohort"] = which
    df["annotation"] = annotation
    df["cldn4_score"] = df["CLDN4_log1p_cp10k"]
    df["cldn4_n"] = df["n_epi"]
    df["cldn4_definition"] = "author epithelial mean log1p(CP10k); CopyKAT IDs not public"
    df["n_malignant"] = df["n_epi"]
    df["n_epithelial"] = df["n_epi"]
    df["n_T"] = df["T"]
    df["n_NK"] = df["NK"]
    df["n_TNK"] = df["TNK"]
    df["n_B"] = df["B"]
    df["n_plasma"] = np.nan
    df["n_B_plasma"] = df["B"]
    df["n_cells_comp"] = df[["epithelial", "TNK", "B", "Myeloid", "Stromal"]].sum(axis=1)
    df["response"] = df["pathology"]
    df["patient_id"] = df["sample_id"]
    df["timing"] = "post"
    # honest eligibility: finite CLDN4 and n_epi >= 10 (same floor as PR290 eligible intent)
    df["eligible"] = df["cldn4_score"].notna() & (df["n_epi"] >= 10)
    high, cut = _median_split(df.loc[df["eligible"], "cldn4_score"])
    df["cldn4_high"] = np.nan
    df.loc[high.index, "cldn4_high"] = high
    df["cldn4_median_cut"] = cut
    df["slice"] = slice_name
    df["primary"] = df["eligible"]
    df["provenance"] = f"PR283 {cohort_key} author_major + PR290 {which} epithelial CLDN4"
    df["Epithelial"] = df["epithelial"]
    long_parts = []
    for ct, src in [("Epithelial", "epithelial"), ("TNK", "TNK"), ("B", "B"), ("Myeloid", "Myeloid"), ("Stromal", "Stromal")]:
        tmp = df[["sample_id", "patient_id", "cohort", "slice", src]].copy()
        tmp["celltype"] = ct
        tmp["n"] = tmp[src]
        long_parts.append(tmp[["cohort", "slice", "sample_id", "patient_id", "celltype", "n"]])
    return df, pd.concat(long_parts, ignore_index=True)


def gse291670() -> tuple[pd.DataFrame, pd.DataFrame]:
    comp = pd.read_csv(EXT / "composition_all_pr284.tsv", sep="\t")
    sub = comp[comp["dataset"] == "GSE291670"].copy()
    scores = pd.read_csv(EXT / "GSE291670_patients.tsv", sep="\t")
    scores["join"] = scores["short"]
    sub["join"] = sub["sample_id"]
    df = sub.merge(scores, on="join", how="left", suffixes=("", "_sc"))
    df["cohort"] = "GSE291670"
    df["subcohort"] = "Camrelizumab"
    df["annotation"] = "marker_hierarchy"
    df["cldn4_score"] = df["mal_CLDN4_mean_log1p_cp10k"]
    df["cldn4_n"] = df["mal_CLDN4_n"]
    df["cldn4_definition"] = "marker-malignant mean log1p(CP10k)"
    df["n_malignant"] = df["n_malignant"]
    df["n_epithelial"] = df["Epithelial"]
    df["n_T"] = np.nan
    df["n_NK"] = np.nan
    df["n_TNK"] = df["T_NK"]
    df["n_B"] = df["B"]
    df["n_plasma"] = np.nan
    df["n_B_plasma"] = df["B"]
    df["n_cells_comp"] = df[["Epithelial", "T_NK", "B", "Myeloid", "Stromal", "Other"]].sum(axis=1)
    df["response"] = df["pathologic_response"]
    df["patient_id"] = df["patient_id"]
    df["timing"] = "post"
    df["eligible"] = df["cldn4_score"].notna() & (df["cldn4_n"] >= 10)
    high, cut = _median_split(df.loc[df["eligible"], "cldn4_score"])
    df["cldn4_high"] = np.nan
    df.loc[high.index, "cldn4_high"] = high
    df["cldn4_median_cut"] = cut
    df["slice"] = "GSE291670"
    df["primary"] = df["eligible"]
    df["provenance"] = "PR284 6-part counts + PR290 malignant CLDN4"
    df["TNK"] = df["T_NK"]
    long_parts = []
    for ct, src in [("Epithelial", "Epithelial"), ("TNK", "T_NK"), ("B", "B"), ("Myeloid", "Myeloid"), ("Stromal", "Stromal"), ("Other", "Other")]:
        tmp = df[["sample_id", "patient_id", "cohort", "slice", src]].copy()
        tmp["celltype"] = ct
        tmp["n"] = tmp[src]
        long_parts.append(tmp[["cohort", "slice", "sample_id", "patient_id", "celltype", "n"]])
    return df, pd.concat(long_parts, ignore_index=True)


def gse205335() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = pd.read_csv(EXT / "GSE205335_patients.tsv", sep="\t")
    df["cohort"] = "GSE205335"
    df["subcohort"] = "ICI_RECIST"
    df["annotation"] = "author_malig_plus_counts"
    df["sample_id"] = df["patient"]
    df["patient_id"] = df["patient"]
    df["cldn4_score"] = df["mal_CLDN4_mean"]
    df["cldn4_n"] = df["n_malignant"]
    df["cldn4_definition"] = "author-malignant mean (PR290 extract)"
    df["n_T"] = np.nan
    df["n_NK"] = np.nan
    df["n_TNK"] = df["n_tnk"]
    df["n_B"] = np.nan
    df["n_plasma"] = np.nan
    df["n_B_plasma"] = df["n_b_plasma"]
    df["n_epithelial"] = df["n_malignant"]
    df["n_cells_comp"] = df["n_cells"]
    df["B"] = df["n_b_plasma"]
    df["TNK"] = df["n_tnk"]
    df["Epithelial"] = df["n_malignant"]
    df["Other"] = (df["n_cells"] - df["n_tnk"] - df["n_b_plasma"] - df["n_malignant"]).clip(lower=0)
    df["Myeloid"] = np.nan
    df["Stromal"] = np.nan
    df["timing"] = "on/post_ICI"
    df["eligible"] = df["cldn4_score"].notna() & (df["n_malignant"] >= 10)
    high, cut = _median_split(df.loc[df["eligible"], "cldn4_score"])
    df["cldn4_high"] = np.nan
    df.loc[high.index, "cldn4_high"] = high
    df["cldn4_median_cut"] = cut
    df["slice"] = "GSE205335"
    df["primary"] = df["eligible"]
    df["provenance"] = "PR290 GSE205335_patients (malignant CLDN4; B = B+plasma only)"
    long_parts = []
    for ct, src in [("Epithelial", "Epithelial"), ("TNK", "TNK"), ("B", "B"), ("Other", "Other")]:
        tmp = df[["sample_id", "patient_id", "cohort", "slice", src]].copy()
        tmp["celltype"] = ct
        tmp["n"] = tmp[src]
        long_parts.append(tmp[["cohort", "slice", "sample_id", "patient_id", "celltype", "n"]])
    return df, pd.concat(long_parts, ignore_index=True)


def gse253013() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Tumor samples only. TNK from PR290; B from PR290 TLS table (marker_coarse TLS was 0)."""
    pat = pd.read_csv(EXT / "GSE253013_patients.tsv", sep="\t")
    tls = pd.read_csv(EXT / "tls_patient_scores.tsv", sep="\t")
    tls = tls[tls["cohort"] == "GSE253013"][["patient", "n_B", "n_plasma", "cldn4_mal_mean_log1p_cp10k", "n_malignant"]].copy()
    tls = tls.rename(columns={"n_malignant": "n_malignant_tls", "cldn4_mal_mean_log1p_cp10k": "cldn4_tls"})
    df = pat[pat["tissue"] == "Tumor"].copy()
    df = df.merge(tls, on="patient", how="left")
    df["cohort"] = "GSE253013"
    df["subcohort"] = "tumor"
    df["annotation"] = "marker_malig_plus_TLS_B"
    df["sample_id"] = df["patient"]
    df["patient_id"] = df["patient"]
    # Prefer the same-table malignant-like CLDN4 (paired with n_malignant_like)
    df["cldn4_score"] = df["CLDN4_mean_log1p_cp10k"]
    df["cldn4_n"] = df["n_malignant_like"]
    df["cldn4_definition"] = "marker-malignant-like mean log1p(CP10k); B counts from TLS extract"
    df["n_malignant"] = df["n_malignant_like"]
    df["n_epithelial"] = df["n_epithelial"]
    df["n_T"] = np.nan
    df["n_NK"] = np.nan
    df["n_TNK"] = df["n_tnk"]
    df["n_B"] = df["n_B"]
    df["n_plasma"] = df["n_plasma"]
    df["n_B_plasma"] = df["n_B"].fillna(0) + df["n_plasma"].fillna(0)
    df["n_cells_comp"] = df["n_cells"]
    df["TNK"] = df["n_tnk"]
    df["B"] = df["n_B"].fillna(0)
    df["Epithelial"] = df["n_epithelial"]
    df["Other"] = (df["n_cells"] - df["n_tnk"] - df["B"] - df["n_epithelial"]).clip(lower=0)
    df["response"] = np.nan
    df["timing"] = "tumor"
    # noisy malignant-like n<10 kept out of primary
    df["eligible"] = df["cldn4_score"].notna() & (df["n_malignant_like"] >= 10)
    high, cut = _median_split(df.loc[df["eligible"], "cldn4_score"])
    df["cldn4_high"] = np.nan
    df.loc[high.index, "cldn4_high"] = high
    df["cldn4_median_cut"] = cut
    df["slice"] = "GSE253013_tumor"
    df["primary"] = df["eligible"]
    df["provenance"] = "PR290 GSE253013_patients (mal CLDN4/TNK) + TLS n_B; 9.3GB RDS not downloaded"
    long_parts = []
    for ct, src in [("Epithelial", "Epithelial"), ("TNK", "TNK"), ("B", "B"), ("Other", "Other")]:
        tmp = df[["sample_id", "patient_id", "cohort", "slice", src]].copy()
        tmp["celltype"] = ct
        tmp["n"] = tmp[src]
        long_parts.append(tmp[["cohort", "slice", "sample_id", "patient_id", "celltype", "n"]])
    return df, pd.concat(long_parts, ignore_index=True)


KEEP = [
    "slice",
    "cohort",
    "subcohort",
    "patient_id",
    "sample_id",
    "annotation",
    "timing",
    "response",
    "n_cells_comp",
    "n_epithelial",
    "n_malignant",
    "n_T",
    "n_NK",
    "n_TNK",
    "n_B",
    "n_plasma",
    "n_B_plasma",
    "frac_TNK",
    "frac_B",
    "frac_B_plasma",
    "cldn4_score",
    "cldn4_n",
    "cldn4_definition",
    "cldn4_median_cut",
    "cldn4_high",
    "eligible",
    "primary",
    "provenance",
]


def finalize(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    n = pd.to_numeric(out["n_cells_comp"], errors="coerce")
    out["frac_TNK"] = pd.to_numeric(out["n_TNK"], errors="coerce") / n
    b = pd.to_numeric(out["n_B"], errors="coerce")
    bp = pd.to_numeric(out["n_B_plasma"], errors="coerce")
    out["frac_B"] = b / n
    out["frac_B_plasma"] = bp / n
    for c in KEEP:
        if c not in out.columns:
            out[c] = np.nan
    return out[KEEP]


def main() -> None:
    frames = []
    longs = []
    builders = [
        gse207422,
        gse207422_drmref,
        lambda: _gse241934("IIT"),
        lambda: _gse241934("RWC"),
        gse291670,
        gse205335,
        gse253013,
    ]
    for fn in builders:
        wide, long = fn()
        frames.append(finalize(wide))
        longs.append(long)
    table = pd.concat(frames, ignore_index=True)
    long = pd.concat(longs, ignore_index=True)
    table.to_csv(OUT / "composition_table.tsv", sep="\t", index=False)
    long.to_csv(OUT / "composition_long.tsv", sep="\t", index=False)

    honest = []
    for sl, sub in table.groupby("slice", sort=False):
        elig = sub[sub["eligible"].astype(bool)]
        honest.append(
            {
                "slice": sl,
                "n_extracted": int(sub["patient_id"].nunique()),
                "n_cldn4_eligible": int(elig["patient_id"].nunique()),
                "n_high": int((elig["cldn4_high"] == 1).sum()),
                "n_low": int((elig["cldn4_high"] == 0).sum()),
                "B_available": bool(elig["n_B"].notna().any() or elig["n_B_plasma"].notna().any()),
                "B_is_B_plasma_only": bool(elig["n_B"].isna().all() and elig["n_B_plasma"].notna().any()),
                "primary": bool(sub["primary"].any()),
                "cldn4_definition": elig["cldn4_definition"].iloc[0] if len(elig) else "",
                "annotation": sub["annotation"].iloc[0],
                "median_cut": float(elig["cldn4_median_cut"].iloc[0]) if len(elig) else float("nan"),
            }
        )
    hon = pd.DataFrame(honest)
    hon.to_csv(OUT / "honest_n.tsv", sep="\t", index=False)
    (OUT / "build_summary.json").write_text(
        json.dumps(
            {
                "n_rows": int(len(table)),
                "slices": hon.to_dict(orient="records"),
                "note": "n is patients. CLDN4-high is a within-slice median split. TACSTD2 is not used.",
            },
            indent=2,
        )
    )
    print(hon.to_string(index=False))
    print(f"wrote {OUT / 'composition_table.tsv'} rows={len(table)}")


if __name__ == "__main__":
    main()
