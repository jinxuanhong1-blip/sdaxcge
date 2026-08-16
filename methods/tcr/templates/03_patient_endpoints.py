#!/usr/bin/env python3
"""Patient-level E1 endpoints from an author TCR table or scirpy cells.tsv.

GSE243013-class default: join
  GSE243013_T_with_TCR_annotation.csv.gz
  GSE243013_NSCLC_immune_scRNA_metadata.csv.gz
and write patient_endpoints.tsv + e2_status.txt.

Does not invent TACSTD2 scores. If --tumor-scores is omitted, E2 is
'not estimable'.

Example:
  python3 03_patient_endpoints.py \\
      --tcr GSE243013_T_with_TCR_annotation.csv.gz \\
      --meta GSE243013_NSCLC_immune_scRNA_metadata.csv.gz \\
      --out workdir/endpoints
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

CD8_TEX = ("CD8T_Tex_CXCL13",)
CD8_TEX_SENS = ("CD8T_Tex_CXCL13", "CD8T_terminal_Tex_LAYN")
CD8_LABEL_SUBSTR = ("CD8T", "CD8")


def _read_table(path: Path) -> pd.DataFrame:
    if str(path).endswith(".gz"):
        return pd.read_csv(path, low_memory=False)
    if path.suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t", low_memory=False)
    return pd.read_csv(path, low_memory=False)


def _bin_response(x: str, rule: str) -> str:
    s = str(x).strip()
    sl = s.lower().replace(" ", "")
    if sl in {"nan", "none", "", "unknowm", "unknown"}:
        return "unknown"
    is_pcr = sl in {"pcr", "pcr"}
    is_mpr = sl == "mpr"
    is_non = sl in {"non-mpr", "nonmpr", "nmpr"}
    if rule == "pcr_or_mpr_vs_nonmpr":
        if is_pcr or is_mpr:
            return "MPR"
        if is_non:
            return "non-MPR"
    elif rule == "pcr_vs_nonpcr":
        if is_pcr:
            return "pCR"
        if is_mpr or is_non:
            return "non-pCR"
    elif rule == "mpr_only_vs_nonmpr":
        if is_pcr:
            return "drop"
        if is_mpr:
            return "MPR"
        if is_non:
            return "non-MPR"
    return "unknown"


def _is_cd8(label: str) -> bool:
    s = str(label)
    return any(k in s for k in CD8_LABEL_SUBSTR)


def endpoints_from_author_tcr(
    tcr: pd.DataFrame,
    sample_meta: pd.DataFrame,
    cxcl13_clusters: tuple[str, ...],
    min_tcr: int,
    min_cd8: int,
) -> pd.DataFrame:
    df = tcr.copy()
    df["clone_id"] = df["TRA_cdr3"].fillna("") + "|" + df["TRB_cdr3"].fillna("")
    # drop cells with no CDR3 on either arm
    df = df[df["clone_id"] != "|"].copy()
    df["is_cxcl13"] = df["sub_cell_type"].isin(cxcl13_clusters)
    df["is_cd8"] = df["sub_cell_type"].map(_is_cd8)

    size = df.groupby(["sampleID", "clone_id"], observed=True)["cellID"].transform("size")
    df["clone_size"] = size
    df["expanded"] = df["clone_size"] >= 2

    rows = []
    for sid, g in df.groupby("sampleID", observed=True):
        n_tcr = len(g)
        n_cd8 = int(g["is_cd8"].sum())
        n_clones = g["clone_id"].nunique()
        n_exp_clones = g.loc[g["expanded"], "clone_id"].nunique()
        n_cxcl13 = int(g["is_cxcl13"].sum())
        n_exp_cxcl13_cells = int((g["is_cxcl13"] & g["expanded"]).sum())
        exp_cxcl13_clones = g.loc[g["is_cxcl13"] & g["expanded"], "clone_id"]
        n_exp_cxcl13_clones = int(exp_cxcl13_clones.nunique())
        rows.append(
            {
                "sample_id": sid,
                "n_tcr": n_tcr,
                "n_cd8_tcr": n_cd8,
                "n_clones": n_clones,
                "n_exp_clones": n_exp_clones,
                "frac_exp_cells": (g["expanded"].sum() / n_tcr) if n_tcr else pd.NA,
                "n_cxcl13": n_cxcl13,
                "n_exp_cxcl13_cells": n_exp_cxcl13_cells,
                "n_exp_cxcl13_clones": n_exp_cxcl13_clones,
                "exp_cxcl13_per_k": (n_exp_cxcl13_clones / n_cd8 * 1000.0) if n_cd8 else pd.NA,
                "frac_cxcl13_that_are_exp": (n_exp_cxcl13_cells / n_cxcl13) if n_cxcl13 else pd.NA,
                "pass_depth": n_tcr >= min_tcr and n_cd8 >= min_cd8,
            }
        )
    out = pd.DataFrame(rows)
    out = out.merge(sample_meta, on="sample_id", how="left")
    return out


def sample_meta_from_gse243013(meta: pd.DataFrame, rule: str) -> pd.DataFrame:
    cols = [
        "sampleID",
        "cancer_type",
        "pathological_response",
        "pathological_response_rate",
        "radiological_response",
        "anti-PD1_therapy",
        "chemotherapy",
    ]
    have = [c for c in cols if c in meta.columns]
    s = meta[have].drop_duplicates("sampleID").rename(columns={"sampleID": "sample_id"})
    s["path_response_bin"] = s["pathological_response"].map(lambda x: _bin_response(x, rule))
    return s


def e2_status(tumor_scores: Path | None) -> str:
    if tumor_scores is None or not Path(tumor_scores).exists():
        return (
            "not estimable: no paired tumor-epithelial TACSTD2/CLDN4 table "
            "for these patient_ids (tacstd2_source=none)."
        )
    tbl = _read_table(Path(tumor_scores))
    need = {"patient_id", "tacstd2"}
    if not need.issubset(set(c.lower() for c in tbl.columns)):
        return (
            "not estimable: tumor-score table present but missing "
            "patient_id and/or tacstd2 columns."
        )
    return f"estimable: {tumor_scores}"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tcr", type=Path, help="Author TCR annotation CSV (GSE243013-class)")
    p.add_argument("--meta", type=Path, help="Author cell metadata CSV")
    p.add_argument("--scirpy-cells", type=Path, help="Alternative: cells.tsv from 01_scirpy_clonotypes.py")
    p.add_argument("--tumor-scores", type=Path, default=None)
    p.add_argument("--path-rule", default="pcr_or_mpr_vs_nonmpr")
    p.add_argument("--min-tcr", type=int, default=50)
    p.add_argument("--min-cd8", type=int, default=20)
    p.add_argument("--include-layn", action="store_true")
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    clusters = CD8_TEX_SENS if args.include_layn else CD8_TEX

    if args.tcr is None and args.scirpy_cells is None:
        raise SystemExit("provide --tcr or --scirpy-cells")

    if args.tcr is not None:
        tcr = _read_table(args.tcr)
        if args.meta is None:
            raise SystemExit("--meta is required with --tcr")
        meta = _read_table(args.meta)
        sample_meta = sample_meta_from_gse243013(meta, args.path_rule)
        end = endpoints_from_author_tcr(tcr, sample_meta, clusters, args.min_tcr, args.min_cd8)
    else:
        cells = _read_table(args.scirpy_cells)
        # Expected columns from template 01: sample_id, clone_id, clone_size, expanded
        if "expanded" not in cells.columns:
            cells["clone_size"] = cells.groupby(["sample_id", "clone_id"], observed=True)[
                "clone_id"
            ].transform("size")
            cells["expanded"] = cells["clone_size"] >= 2
        if "is_cxcl13" not in cells.columns:
            cells["is_cxcl13"] = False
        if "is_cd8" not in cells.columns:
            cells["is_cd8"] = True
        rows = []
        for sid, g in cells.groupby("sample_id", observed=True):
            n_tcr = len(g)
            n_cd8 = int(g["is_cd8"].sum())
            n_cxcl13 = int(g["is_cxcl13"].sum())
            n_exp_cxcl13_clones = g.loc[g["is_cxcl13"] & g["expanded"], "clone_id"].nunique()
            rows.append(
                {
                    "sample_id": sid,
                    "n_tcr": n_tcr,
                    "n_cd8_tcr": n_cd8,
                    "n_clones": g["clone_id"].nunique(),
                    "n_exp_clones": g.loc[g["expanded"], "clone_id"].nunique(),
                    "frac_exp_cells": g["expanded"].mean() if n_tcr else pd.NA,
                    "n_cxcl13": n_cxcl13,
                    "n_exp_cxcl13_cells": int((g["is_cxcl13"] & g["expanded"]).sum()),
                    "n_exp_cxcl13_clones": int(n_exp_cxcl13_clones),
                    "exp_cxcl13_per_k": (n_exp_cxcl13_clones / n_cd8 * 1000.0) if n_cd8 else pd.NA,
                    "frac_cxcl13_that_are_exp": (
                        (g["is_cxcl13"] & g["expanded"]).sum() / n_cxcl13 if n_cxcl13 else pd.NA
                    ),
                    "pass_depth": n_tcr >= args.min_tcr and n_cd8 >= args.min_cd8,
                }
            )
        end = pd.DataFrame(rows)

    end.to_csv(args.out / "patient_endpoints.tsv", sep="\t", index=False)
    spec = {
        "primary_endpoint": "exp_cxcl13_per_k",
        "path_response_bin_rule": args.path_rule,
        "cxcl13_clusters": list(clusters),
        "min_tcr_cells": args.min_tcr,
        "min_cd8_tcr_cells": args.min_cd8,
        "n_samples": int(end.shape[0]),
        "n_pass_depth": int(end["pass_depth"].sum()) if "pass_depth" in end else None,
        "unit": "patient/sample",
        "note": "No inference is run in this template. Tests belong in a locked spec.",
    }
    (args.out / "e1_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    (args.out / "e2_status.txt").write_text(e2_status(args.tumor_scores) + "\n")
    print(f"wrote {args.out}  samples={end.shape[0]}")


if __name__ == "__main__":
    main()
