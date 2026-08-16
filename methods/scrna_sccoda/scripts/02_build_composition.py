#!/usr/bin/env python3
"""Build sample-level composition and malignant TACSTD2 tables (patient unit)."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from lib.annotations import (  # noqa: E402
    AUTHOR_MAJOR_CANON,
    HU_LINEAGES,
    NORMAL_LUNG,
    collapse_label,
    is_tls_fine,
    is_tnk_fine,
)

PAPER_GROUP_207422 = {
    "BD_immune01": "TN",
    "BD_immune02": "NMPR",
    "BD_immune03": "MPR",
    "BD_immune04": "NMPR",
    "BD_immune05": "TN",
    "BD_immune06": "MPR",
    "BD_immune07": "NMPR",
    "BD_immune08": "TN",
    "BD_immune09": "NMPR",
    "BD_immune10": "NMPR",
    "BD_immune11": "MPR",
    "BD_immune12": "NMPR",
    "BD_immune13": "NMPR",
    "BD_immune14": "MPR",
    "BD_immune15": "NMPR",
}


def _score(expr: dict[str, np.ndarray], genes: list[str], n: int) -> np.ndarray:
    mats = [np.log1p(expr[g].astype(np.float32)) for g in genes if g in expr]
    if not mats:
        return np.zeros(n, dtype=np.float32)
    return np.mean(np.vstack(mats), axis=0)


def assign_hu_lineage(expr: dict[str, np.ndarray], n: int) -> np.ndarray:
    names = list(HU_LINEAGES)
    scores = np.vstack([_score(expr, HU_LINEAGES[k], n) for k in names])
    best = scores.argmax(axis=0)
    best_val = scores.max(axis=0)
    second = np.partition(scores, -2, axis=0)[-2]
    assigned = np.array(names, dtype=object)[best]
    assigned[(best_val < 0.15) | ((best_val - second) < 0.05)] = "other"
    return assigned


def malignant_like(expr: dict[str, np.ndarray], lineage: np.ndarray) -> np.ndarray:
    n = len(lineage)
    normal = _score(expr, NORMAL_LUNG, n)
    return (lineage == "epithelial") & (normal < 0.05)


def load_npz(path: Path) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray | None]:
    z = np.load(path, allow_pickle=True)
    cells = z["cells"]
    genes = list(z["genes"])
    mat = z["mat"]
    expr = {g: mat[i] for i, g in enumerate(genes)}
    total = z["total"] if "total" in z.files else None
    return cells, expr, total


def _add_counts(rows: list[dict], cohort, sample, annotation, counter) -> None:
    for lab, n in counter.items():
        rows.append(
            {
                "cohort": cohort,
                "sample_id": sample,
                "annotation": annotation,
                "celltype": lab,
                "n": int(n),
            }
        )


def build_gse207422(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cells, expr, total = load_npz(root / "data" / "GSE207422" / "extracted_markers.npz")
    n = len(cells)
    sample_ids = np.array([str(c).rsplit("_", 1)[0] for c in cells])
    lineage = assign_hu_lineage(expr, n)
    malig = malignant_like(expr, lineage)
    tac = expr["TACSTD2"].astype(np.float64)
    lib = total.astype(np.float64) if total is not None else np.ones(n)

    meta_xlsx = root / "data" / "GSE207422" / "GSE207422_NSCLC_scRNAseq_metadata.xlsx"
    raw = pd.read_excel(meta_xlsx)
    raw = raw.dropna(subset=["Sample"])
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)].copy()
    raw["sample_id"] = raw["Sample"].astype(str)
    raw["paper_group"] = raw["sample_id"].map(PAPER_GROUP_207422)
    raw["mpr"] = raw["Pathologic Response"].replace({"pCR": "MPR"})
    raw.loc[raw["paper_group"] == "TN", "mpr"] = "NA_pre"
    raw["timing"] = np.where(raw["paper_group"] == "TN", "pre", "post")

    cov_rows = []
    comp_rows = []
    for sid, sub_idx in _groupby_indices(sample_ids):
        lin = lineage[sub_idx]
        msk_mal = malig[sub_idx]
        tnk = np.isin(lin, ["T", "NK"])
        hu_counts = pd.Series(lin).value_counts().to_dict()
        _add_counts(comp_rows, "GSE207422", sid, "hu_markers", hu_counts)
        collapsed = {collapse_label("hu_collapsed", k): 0 for k in set(hu_counts) | {"T", "NK", "B"}}
        for k, v in hu_counts.items():
            collapsed[collapse_label("hu_collapsed", k)] = collapsed.get(collapse_label("hu_collapsed", k), 0) + v
        _add_counts(comp_rows, "GSE207422", sid, "hu_collapsed", collapsed)
        coarse = {
            "TNK": int(tnk.sum()),
            "TLS": int(np.isin(lin, ["B", "plasma"]).sum()),
            "epithelial": int((lin == "epithelial").sum()),
            "other": int((~np.isin(lin, ["T", "NK", "B", "plasma", "epithelial"])).sum()),
        }
        _add_counts(comp_rows, "GSE207422", sid, "marker_coarse", coarse)

        n_mal = int(msk_mal.sum())
        if n_mal >= 1:
            umi = tac[sub_idx][msk_mal]
            lib_m = lib[sub_idx][msk_mal]
            cp = np.where(lib_m > 0, umi / lib_m * 1e4, np.nan)
            mal_log1p = float(np.mean(np.log1p(umi)))
            mal_cp = float(np.nanmean(np.log1p(cp)))
            mal_pos = float(np.mean(umi > 0))
        else:
            mal_log1p = mal_cp = mal_pos = np.nan
        meta = raw.loc[raw["sample_id"] == sid]
        cov_rows.append(
            {
                "cohort": "GSE207422",
                "subcohort": "Hu2023",
                "sample_id": sid,
                "n_cells": int(len(sub_idx)),
                "n_epithelial": int((lin == "epithelial").sum()),
                "n_malignant_like": n_mal,
                "n_tnk": int(tnk.sum()),
                "tacstd2_malig_mean_log1p": mal_log1p,
                "tacstd2_malig_mean_log1p_cp10k": mal_cp,
                "tacstd2_malig_pct_pos": mal_pos,
                "mpr": (meta["mpr"].iloc[0] if len(meta) else "unknown"),
                "timing": (meta["timing"].iloc[0] if len(meta) else "unknown"),
                "histology": (str(meta["Pathology"].iloc[0]) if len(meta) else ""),
                "eligible_tacstd2": int(n_mal >= 10),
                "eligible_mpr": int(len(meta) > 0 and meta["timing"].iloc[0] == "post"),
            }
        )

    # DRMref public-derived (sibling open extract; not author CopyKAT)
    drm_path = root / "data" / "public_derived" / "GSE207422_drmref_celltype_counts.tsv"
    drm_tac = root / "data" / "public_derived" / "GSE207422_drmref_malig_tacstd2.tsv"
    if drm_path.exists():
        drm = pd.read_csv(drm_path, sep="\t")
        for (patient, resp), g in drm.groupby(["patient", "response"]):
            sid = {
                "P02": "BD_immune02",
                "P03": "BD_immune03",
                "P04": "BD_immune04",
                "P06": "BD_immune06",
                "P07": "BD_immune07",
                "P09": "BD_immune09",
                "P10": "BD_immune10",
                "P11": "BD_immune11",
                "P12": "BD_immune12",
                "P13": "BD_immune13",
                "P14": "BD_immune14",
                "P15": "BD_immune15",
            }.get(str(patient))
            if sid is None:
                continue
            raw_counts = dict(zip(g["celltype"], g["n"]))
            _add_counts(comp_rows, "GSE207422", sid, "drmref", raw_counts)
            collapsed = {}
            for k, v in raw_counts.items():
                ck = collapse_label("drmref_collapsed", k)
                collapsed[ck] = collapsed.get(ck, 0) + int(v)
            _add_counts(comp_rows, "GSE207422", sid, "drmref_collapsed", collapsed)
        if drm_tac.exists():
            dt = pd.read_csv(drm_tac, sep="\t")
            cov = pd.DataFrame(cov_rows).set_index("sample_id")
            for _, r in dt.iterrows():
                sid = r["sample"]
                if sid in cov.index:
                    cov.loc[sid, "tacstd2_drmref_malig_log1p_cp10k"] = r["mean_TACSTD2_log1p_cp10k_mal"]
                    cov.loc[sid, "n_drmref_malignant"] = r["n_malignant"]
            cov_rows = cov.reset_index().to_dict("records")

    return pd.DataFrame(comp_rows), pd.DataFrame(cov_rows)


def _groupby_indices(ids: np.ndarray) -> list[tuple[str, np.ndarray]]:
    out = []
    for sid in pd.unique(ids):
        out.append((str(sid), np.flatnonzero(ids == sid)))
    return out


def _norm_mpr(val: str) -> str:
    v = str(val).strip()
    if v in {"MPR", "pCR"}:
        return "MPR"
    if v in {"non-MPR", "NMPR"}:
        return "NMPR"
    return v


def build_gse241934(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = root / "data" / "GSE241934"
    parts = [
        ("IIT", d / "GSE241934_IIT_Meta.txt.gz", d / "extracted_IIT_markers.npz"),
        ("RWC", d / "GSE241934_Real_Meta.txt.gz", d / "extracted_RWC_markers.npz"),
    ]
    comp_rows, cov_rows = [], []
    for sub, meta_path, npz_path in parts:
        if not meta_path.exists():
            continue
        with gzip.open(meta_path, "rt") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            idx = {h: i for i, h in enumerate(header)}
            cells = []
            for line in fh:
                p = line.rstrip("\n").split("\t")
                cells.append(
                    {
                        "cellID": p[idx["cellID"]],
                        "sample_id": p[idx["sampleID"]],
                        "major": AUTHOR_MAJOR_CANON.get(p[idx["major.cell.type"]], p[idx["major.cell.type"]]),
                        "fine": p[idx["cell.type"]],
                        "mpr": _norm_mpr(p[idx["Pathological Response"]]),
                        "rate": p[idx["Pathological Response Rate"]],
                        "histology": p[idx["Histology"]],
                        "egfr": p[idx["EGFR"]],
                        "nCount": float(p[idx["nCount_RNA"]]),
                    }
                )
        meta = pd.DataFrame(cells)
        expr = None
        cells_npz = None
        if npz_path.exists():
            cells_npz, expr, _ = load_npz(npz_path)
            order = {c: i for i, c in enumerate(cells_npz)}
            meta["_i"] = meta["cellID"].map(order)
            if meta["_i"].isna().any():
                nmiss = int(meta["_i"].isna().sum())
                print(f"  warn {sub}: {nmiss} meta barcodes missing from MTX extract", flush=True)
        for sid, g in meta.groupby("sample_id"):
            major_counts = g["major"].value_counts().to_dict()
            _add_counts(comp_rows, f"GSE241934_{sub}", sid, "author_major", major_counts)
            collapsed = {}
            for k, v in major_counts.items():
                ck = collapse_label("author_collapsed", k)
                collapsed[ck] = collapsed.get(ck, 0) + int(v)
            _add_counts(comp_rows, f"GSE241934_{sub}", sid, "author_collapsed", collapsed)
            fine_counts = defaultdict(int)
            for lab in g["fine"]:
                if is_tnk_fine(lab) or lab in {"T", "NK"}:
                    fine_counts["TNK"] += 1
                elif is_tls_fine(lab) or lab == "B":
                    fine_counts["TLS"] += 1
                elif lab in {"NA", "nan", ""}:
                    # unassigned fine labels: fall back to major via this cell later
                    pass
                else:
                    fine_counts["other_immune"] += 1
            # cells with NA fine label keep major collapse
            na = g["fine"].isin(["NA", "nan", ""])
            for k, v in g.loc[na, "major"].value_counts().items():
                ck = collapse_label("author_collapsed", k)
                fine_counts[ck] = fine_counts.get(ck, 0) + int(v)
            _add_counts(comp_rows, f"GSE241934_{sub}", sid, "author_fine_immune", dict(fine_counts))

            if expr is not None and g["_i"].notna().any():
                ii_all = g.loc[g["_i"].notna(), "_i"].astype(int).to_numpy()
                expr_s = {k: v[ii_all] for k, v in expr.items()}
                hu = assign_hu_lineage(expr_s, len(ii_all))
                hu_counts = pd.Series(hu).value_counts().to_dict()
                _add_counts(comp_rows, f"GSE241934_{sub}", sid, "hu_markers", hu_counts)
                collapsed_h = {}
                for k, v in hu_counts.items():
                    ck = collapse_label("hu_collapsed", k)
                    collapsed_h[ck] = collapsed_h.get(ck, 0) + int(v)
                _add_counts(comp_rows, f"GSE241934_{sub}", sid, "hu_collapsed", collapsed_h)

            n_epi = int((g["major"] == "epithelial").sum())
            mal_log1p = mal_cp = mal_pos = np.nan
            n_mal = 0
            if expr is not None and "TACSTD2" in expr and g["_i"].notna().any():
                ii = g.loc[g["_i"].notna() & (g["major"] == "epithelial"), "_i"].astype(int).to_numpy()
                n_mal = int(len(ii))
                if n_mal:
                    umi = expr["TACSTD2"][ii].astype(np.float64)
                    lib = g.loc[g["_i"].notna() & (g["major"] == "epithelial"), "nCount"].to_numpy()
                    mal_log1p = float(np.mean(np.log1p(umi)))
                    cp = np.where(lib > 0, umi / lib * 1e4, np.nan)
                    mal_cp = float(np.nanmean(np.log1p(cp)))
                    mal_pos = float(np.mean(umi > 0))
            cov_rows.append(
                {
                    "cohort": f"GSE241934_{sub}",
                    "subcohort": sub,
                    "sample_id": sid,
                    "n_cells": int(len(g)),
                    "n_epithelial": n_epi,
                    "n_malignant_like": n_mal,
                    "n_tnk": int(g["major"].isin(["T", "NK"]).sum()),
                    "tacstd2_malig_mean_log1p": mal_log1p,
                    "tacstd2_malig_mean_log1p_cp10k": mal_cp,
                    "tacstd2_malig_pct_pos": mal_pos,
                    "mpr": g["mpr"].iloc[0],
                    "timing": "post",
                    "histology": g["histology"].iloc[0],
                    "egfr": g["egfr"].iloc[0],
                    "eligible_tacstd2": int(n_mal >= 10 or n_epi >= 10),
                    "eligible_mpr": 1,
                }
            )
    comp = pd.DataFrame(comp_rows)
    cov = pd.DataFrame(cov_rows)
    if not comp.empty:
        pooled_c = comp.copy()
        pooled_c["cohort"] = "GSE241934"
        pooled_v = cov.copy()
        pooled_v["cohort"] = "GSE241934"
        comp = pd.concat([comp, pooled_c], ignore_index=True)
        cov = pd.concat([cov, pooled_v], ignore_index=True)
    return comp, cov


def build_gse253013(root: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Public-derived patient tables only. 9.3 GB RDS is not re-downloaded."""
    p = root / "data" / "public_derived" / "GSE253013_per_patient_metrics.tsv"
    a = root / "data" / "public_derived" / "GSE253013_author_label_per_patient.tsv"
    if not p.exists():
        return pd.DataFrame(), pd.DataFrame()
    df = pd.read_csv(p, sep="\t")
    tumor = df[df["tissue"] == "Tumor"].copy()
    comp_rows, cov_rows = [], []
    author = pd.read_csv(a, sep="\t") if a.exists() else None
    for _, r in tumor.iterrows():
        sid = r["patient"]
        n_cells = int(r["n_cells"])
        n_epi = int(r["n_epithelial"])
        n_mal = int(r["n_malignant_like"])
        n_tnk = int(r["n_tnk"])
        n_other = max(n_cells - n_epi - n_tnk, 0)
        n_other_epi = max(n_epi - n_mal, 0)
        _add_counts(
            comp_rows,
            "GSE253013",
            sid,
            "marker_coarse",
            {"TNK": n_tnk, "TLS": 0, "epithelial": n_epi, "other": n_other},
        )
        _add_counts(
            comp_rows,
            "GSE253013",
            sid,
            "marker_malig_split",
            {
                "TNK": n_tnk,
                "malignant_like": n_mal,
                "other_epithelial": n_other_epi,
                "other": n_other,
            },
        )
        if author is not None and sid in set(author["patient"]):
            ar = author.loc[author["patient"] == sid].iloc[0]
            n_ae = int(ar["n_author_epithelial"])
            n_at = int(ar["n_author_tcells"])
            # author table does not give a full composition; residual is unknown
            _add_counts(
                comp_rows,
                "GSE253013",
                sid,
                "author_garnett_limited",
                {"TNK": n_at, "epithelial": n_ae, "other": max(n_cells - n_ae - n_at, 0)},
            )
        cov_rows.append(
            {
                "cohort": "GSE253013",
                "subcohort": "Sze2024_tumor",
                "sample_id": sid,
                "n_cells": n_cells,
                "n_epithelial": n_epi,
                "n_malignant_like": n_mal,
                "n_tnk": n_tnk,
                "tacstd2_malig_mean_log1p": float(r["TACSTD2_mean_log1p"]),
                "tacstd2_malig_mean_log1p_cp10k": float(r["TACSTD2_mean_log1p_cp10k"]),
                "tacstd2_malig_pct_pos": float(r["TACSTD2_pct_pos"]) / 100.0,
                "mpr": "NA_no_public_label",
                "timing": "treatment_naive",
                "histology": "LUAD",
                "eligible_tacstd2": int(bool(r["eligible_malig"])),
                "eligible_mpr": 0,
                "provenance": "public_derived_sibling_extract_of_GEO_RDS",
            }
        )
    return pd.DataFrame(comp_rows), pd.DataFrame(cov_rows)


def median_split(cov: pd.DataFrame) -> pd.DataFrame:
    out = cov.copy()
    out["tacstd2_high"] = np.nan
    for cohort, g in out.groupby("cohort"):
        elig = g["eligible_tacstd2"].astype(bool) & g["tacstd2_malig_mean_log1p"].notna()
        if elig.sum() < 4:
            continue
        med = g.loc[elig, "tacstd2_malig_mean_log1p"].median()
        out.loc[g.index[elig], "tacstd2_high"] = (
            g.loc[elig, "tacstd2_malig_mean_log1p"] > med
        ).astype(float)
        out.loc[g.index[elig], "tacstd2_median_cut"] = med
    # z-score within cohort among eligible
    out["tacstd2_z"] = np.nan
    for cohort, g in out.groupby("cohort"):
        elig = g["eligible_tacstd2"].astype(bool) & g["tacstd2_malig_mean_log1p"].notna()
        x = g.loc[elig, "tacstd2_malig_mean_log1p"]
        if elig.sum() < 3 or x.std(ddof=1) == 0:
            continue
        out.loc[g.index[elig], "tacstd2_z"] = (x - x.mean()) / x.std(ddof=1)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args()
    out = args.root / "results"
    out.mkdir(parents=True, exist_ok=True)

    comps, covs = [], []
    for fn, name in (
        (build_gse207422, "GSE207422"),
        (build_gse241934, "GSE241934"),
        (build_gse253013, "GSE253013"),
    ):
        try:
            c, v = fn(args.root)
        except FileNotFoundError as exc:
            print(f"SKIP {name}: {exc}", flush=True)
            continue
        print(f"{name}: composition rows={len(c)} covariate rows={len(v)}", flush=True)
        comps.append(c)
        covs.append(v)
    composition = pd.concat(comps, ignore_index=True) if comps else pd.DataFrame()
    covariates = median_split(pd.concat(covs, ignore_index=True)) if covs else pd.DataFrame()
    composition.to_csv(out / "composition_long.tsv", sep="\t", index=False)
    covariates.to_csv(out / "sample_covariates.tsv", sep="\t", index=False)
    audit = {
        "n_composition_rows": int(len(composition)),
        "n_samples": int(covariates[["cohort", "sample_id"]].drop_duplicates().shape[0]) if len(covariates) else 0,
        "cohorts": sorted(covariates["cohort"].unique().tolist()) if len(covariates) else [],
        "annotations": sorted(composition["annotation"].unique().tolist()) if len(composition) else [],
        "gse253013_rds_downloaded": False,
        "gse253013_note": "9.3 GB RDS over 2 GB budget; used public-derived patient tables",
    }
    (out / "build_audit.json").write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
