#!/usr/bin/env python3
"""Flatten computed stats into one honest CSV + markdown table."""
import json
from pathlib import Path

import pandas as pd

OUT = Path(__file__).resolve().parents[2] / "results" / "noskip" / "zenodo"


def fmt_p(p):
    if p is None or (isinstance(p, float) and (p != p)):
        return ""
    if isinstance(p, str):
        return p
    p = float(p)
    if p == 0.0:
        return "<1e-300"
    if p < 1e-6:
        return f"{p:.3e}"
    return f"{p:.6g}"


def row(**kw):
    keys = [
        "dataset", "doi", "modality", "subset", "n",
        "TACSTD2_present", "CLDN4_present",
        "TACSTD2_pct", "CLDN4_pct", "rho", "p", "note",
    ]
    out = {k: kw.get(k) for k in keys}
    out["p"] = fmt_p(out.get("p"))
    if out.get("rho") is not None:
        out["rho"] = float(out["rho"])
    return out


def main():
    rows = []
    # 10731914 human
    s = json.loads((OUT / "10731914_stats.json").read_text())
    for key, subset in [("spearman_all_cells", "all_cells_cp10k"), ("spearman_epithelial", "epithelial_cp10k")]:
        sp = s[key]
        extra = {}
        if subset.startswith("epithelial"):
            extra = {"TACSTD2_pct": 73.48, "CLDN4_pct": 79.97}  # from luad_expression_by_celltype
        else:
            extra = {"TACSTD2_pct": 26.88, "CLDN4_pct": 27.53}
        rows.append(row(
            dataset="8-LUAD human scRNA-seq",
            doi="10.5281/zenodo.10731914",
            modality="scRNA-seq counts, CP10K for ρ",
            subset=subset,
            n=sp["n"],
            TACSTD2_present=True, CLDN4_present=True,
            TACSTD2_pct=extra["TACSTD2_pct"], CLDN4_pct=extra["CLDN4_pct"],
            rho=sp["rho"], p=("<1e-300" if sp["p"] == 0.0 else sp["p"]),
            note="p=0.0 from scipy underflow reported as <1e-300",
        ))
    # visium
    vis = pd.read_csv(OUT / "8417887_visium_per_sample.csv")
    for _, r in vis.iterrows():
        rows.append(row(
            dataset="early LUAD Visium",
            doi="10.5281/zenodo.8417887",
            modality="10x Visium filtered_feature_bc_matrix.h5, raw counts",
            subset=r["sample"], n=int(r["n"]),
            TACSTD2_present=True, CLDN4_present=True,
            TACSTD2_pct=r["TACSTD2_pct"], CLDN4_pct=r["CLDN4_pct"],
            rho=r["rho"], p=r["p"],
        ))
    vp = json.loads((OUT / "8417887_visium_pooled_spearman.json").read_text())
    rows.append(row(
        dataset="early LUAD Visium", doi="10.5281/zenodo.8417887",
        modality="10x Visium raw counts", subset="all_spots_pooled",
        n=vp["n"], TACSTD2_present=True, CLDN4_present=True,
        TACSTD2_pct=vp["TACSTD2_pct"], CLDN4_pct=vp["CLDN4_pct"],
        rho=vp["rho"], p=vp["p"],
    ))
    # 8417887 scrna tables
    scr = json.loads((OUT / "8417887_scrna_tables.json").read_text())
    for lab in ("malignant", "immune"):
        r = scr[lab]
        sp = r["spearman"]
        rows.append(row(
            dataset=f"early LUAD scRNA {lab}",
            doi="10.5281/zenodo.8417887",
            modality="dense gene x cell txt, raw",
            subset=lab, n=r["n"],
            TACSTD2_present=True, CLDN4_present=True,
            TACSTD2_pct=r["TACSTD2_pct"], CLDN4_pct=r["CLDN4_pct"],
            rho=sp["rho"], p=sp["p"],
        ))
    # 11205626
    paired = pd.read_csv(OUT / "11205626_cell_protocol_per_sample.csv")
    for _, r in paired.iterrows():
        rows.append(row(
            dataset="paired normal-LUAD 10x (Cell protocol only)",
            doi="10.5281/zenodo.11205626",
            modality="10x h5 raw counts",
            subset=f"{r['tissue']} patient {r['patient']}",
            n=int(r["n"]), TACSTD2_present=True, CLDN4_present=True,
            TACSTD2_pct=r["TACSTD2_pct"], CLDN4_pct=r["CLDN4_pct"],
            rho=r["rho"], p=r["p"],
        ))
    pool = json.loads((OUT / "11205626_pooled_spearman.json").read_text())
    for k in ("normal", "tumor", "all"):
        r = pool[k]
        rows.append(row(
            dataset="paired normal-LUAD 10x (Cell protocol only)",
            doi="10.5281/zenodo.11205626",
            modality="10x h5 raw counts",
            subset=r["subset"], n=r["n"],
            TACSTD2_present=True, CLDN4_present=True,
            TACSTD2_pct=r["TACSTD2_pct"], CLDN4_pct=r["CLDN4_pct"],
            rho=r["rho"], p=r["p"],
            note="SN and immune-depleted libraries excluded to avoid double-counting",
        ))
    # mouse
    m = json.loads((OUT / "10731914_mouse_stats.json").read_text())
    rows.append(row(
        dataset="mouse LUAD scRNA-seq",
        doi="10.5281/zenodo.10731914",
        modality="scRNA-seq raw UMI (symbols Tacstd2/Cldn4)",
        subset="all_cells", n=m["n_cells"],
        TACSTD2_present=True, CLDN4_present=True,
        TACSTD2_pct=m["Tacstd2_pct"], CLDN4_pct=m["Cldn4_pct"],
        rho=m["spearman_all"]["rho"], p=m["spearman_all"]["p"],
        note="Cldn4 detected in ~1/64000 cells; ρ is not interpretable",
    ))
    # tcells
    tc = json.loads((OUT / "13947395_tcell_expression.json").read_text())
    for r in tc:
        rows.append(row(
            dataset="Fig6 processed T-cell objects",
            doi="10.5281/zenodo.13947395",
            modality="h5ad X (cells x genes)",
            subset=r["file"], n=r.get("n"),
            TACSTD2_present=r["TACSTD2_present"], CLDN4_present=r["CLDN4_present"],
            TACSTD2_pct=r.get("TACSTD2_pct"), CLDN4_pct=r.get("CLDN4_pct"),
            rho=r.get("rho"), p=r.get("p"),
            note="T-cell-only objects; both genes nearly absent (1–2% cells)",
        ))
    # negatives
    rows.append(row(
        dataset="Nanostring nCounter NSCLC anti-PD1",
        doi="10.5281/zenodo.2635194",
        modality="nCounter 784-gene immune panel",
        subset="learn+validate n=55+36", n=91,
        TACSTD2_present=False, CLDN4_present=False,
        note="Neither gene is on the panel. ρ not computed.",
    ))
    rows.append(row(
        dataset="I3LUNG clinical/radiomics",
        doi="10.64898/2026.01.16.25342913",
        modality="clinical + genomics drivers + radiomics",
        subset="I3LUNG_DATA.zip", n=None,
        TACSTD2_present=False, CLDN4_present=False,
        note="No expression matrix.",
    ))
    rows.append(row(
        dataset="anti-PD1 NSCLC IMC",
        doi="10.5281/zenodo.8041882",
        modality="41-marker IMC protein panel",
        subset="all cells", n=99659,
        TACSTD2_present=False, CLDN4_present=False,
        note="Protein panel; TACSTD2/CLDN4 not among antibodies. ρ not computed.",
    ))
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "master_stats.csv", index=False)
    print(df.to_string(index=False))
    print("\nWrote", OUT / "master_stats.csv")


if __name__ == "__main__":
    main()
