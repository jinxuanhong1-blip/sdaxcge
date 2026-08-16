#!/usr/bin/env python3
"""Leftover GeoMx: GSE265899, GSE289483, GSE334014, GSE326968."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import B_GENES, BROAD_EPI, EPI, T_GENES, dump_json, present, score_mean, spearman  # noqa: E402

OUT = ROOT / "results" / "leftover_spatial"


def load_gene_x_aoi(path, index_col=0):
    df = pd.read_csv(path, index_col=index_col)
    df.index = df.index.astype(str).str.strip().str.upper()
    # drop empty columns
    df = df.loc[:, df.columns.notna()]
    df = df.apply(pd.to_numeric, errors="coerce")
    return df


def scores_from_gx(gx):
    genes = list(gx.index)
    aoi = gx.T
    aoi.columns = [str(c).upper() for c in aoi.columns]
    return pd.DataFrame(
        {
            "CLDN4": aoi["CLDN4"] if "CLDN4" in aoi.columns else np.nan,
            "TACSTD2": aoi["TACSTD2"] if "TACSTD2" in aoi.columns else np.nan,
            "epi": score_mean(aoi, present(EPI, genes)),
            "broad": score_mean(aoi, present(BROAD_EPI, genes)),
            "T": score_mean(aoi, present(T_GENES, genes)),
            "B": score_mean(aoi, present(B_GENES, genes)),
            "TB": score_mean(aoi, present(T_GENES + B_GENES, genes)),
            "PTPRC": aoi["PTPRC"] if "PTPRC" in aoi.columns else np.nan,
        },
        index=aoi.index,
    )


def corr_block(sc, label):
    out = {"subset": label, "n": int(len(sc))}
    for a, b in [
        ("CLDN4", "TB"),
        ("TACSTD2", "TB"),
        ("epi", "TB"),
        ("CLDN4", "T"),
        ("TACSTD2", "T"),
        ("CLDN4", "B"),
        ("TACSTD2", "B"),
        ("CLDN4", "PTPRC"),
        ("TACSTD2", "PTPRC"),
        ("broad", "TB"),
    ]:
        out[f"{a}_vs_{b}"] = spearman(sc[a], sc[b])
    return out


def run_gse265899():
    gx = load_gene_x_aoi(ROOT / "data" / "GSE265899" / "GSE265899_Q3Norm.csv.gz")
    # index may be 'GENES' header already used
    sc = scores_from_gx(gx)
    sc["kind"] = np.where(sc.index.str.lower().str.startswith("tumor"), "tumor", "immune")
    sc["pair"] = sc.index.str.extract(r"(\d+)$", expand=False)
    recs = [corr_block(sc, "all_AOI"), corr_block(sc[sc.kind == "tumor"], "tumor"), corr_block(sc[sc.kind == "immune"], "immune")]
    # paired tumor vs immune for target genes
    paired = sc.dropna(subset=["pair"]).pivot_table(index="pair", columns="kind", values=["CLDN4", "TACSTD2", "TB", "T", "B"])
    paired_stats = {}
    for gene in ["CLDN4", "TACSTD2", "TB"]:
        if (gene, "tumor") in paired.columns and (gene, "immune") in paired.columns:
            t = paired[(gene, "tumor")]
            i = paired[(gene, "immune")]
            m = t.notna() & i.notna()
            from scipy import stats as st

            if m.sum() >= 5:
                w = st.wilcoxon(t[m], i[m], alternative="two-sided")
                paired_stats[gene] = {
                    "n_pairs": int(m.sum()),
                    "median_tumor": float(t[m].median()),
                    "median_immune": float(i[m].median()),
                    "wilcoxon_p": float(w.pvalue),
                    "n_tumor_gt_immune": int((t[m] > i[m]).sum()),
                }
    return {"series": "GSE265899", "platform": "GeoMx WTA", "corrs": recs, "paired": paired_stats, "n_tumor": int((sc.kind == "tumor").sum()), "n_immune": int((sc.kind == "immune").sum())}


def run_gse289483():
    gx = load_gene_x_aoi(ROOT / "data" / "GSE289483" / "GSE289483_processed_q3norm_gene_expr.csv.gz")
    sc = scores_from_gx(gx)
    ann_path = ROOT / "results" / "leftover_spatial" / "tables" / "gse289483_sample_annotations.tsv"
    if not ann_path.exists():
        ann_path = ROOT / "data" / "GSE289483" / "sample_annotations.tsv"
    ann = pd.read_csv(ann_path, sep="\t")
    ann["dsp_id"] = ann["dsp_id"].astype(str)
    sc = sc.copy()
    sc["dsp_id"] = sc.index.astype(str)
    # matrix IDs already DSP_...
    merged = sc.merge(ann, on="dsp_id", how="left")
    recs = [corr_block(merged, "all_AOI")]
    for col, val in [
        ("segmentation", "CD45-"),
        ("segmentation", "CD45+"),
        ("component", "Sarcomatoid"),
        ("component", "Adenocarcinomas"),
        ("component", "Normal"),
    ]:
        sub = merged[merged[col] == val]
        if len(sub) >= 8:
            recs.append(corr_block(sub, f"{col}={val}"))
    # CD45- vs CD45+ target levels
    from scipy import stats as st

    levels = {}
    a = merged[merged.segmentation == "CD45-"]
    b = merged[merged.segmentation == "CD45+"]
    for gene in ["CLDN4", "TACSTD2", "TB"]:
        if len(a) >= 5 and len(b) >= 5:
            u = st.mannwhitneyu(a[gene].dropna(), b[gene].dropna(), alternative="two-sided")
            levels[gene] = {
                "n_CD45minus": int(a[gene].notna().sum()),
                "n_CD45plus": int(b[gene].notna().sum()),
                "median_CD45minus": float(a[gene].median()),
                "median_CD45plus": float(b[gene].median()),
                "mw_p": float(u.pvalue),
            }
    return {
        "series": "GSE289483",
        "platform": "GeoMx WTA",
        "n_aoi_matrix": int(len(sc)),
        "n_aoi_annotated": int(merged["component"].notna().sum()),
        "corrs": recs,
        "CD45_levels": levels,
    }


def run_gse334014():
    df = pd.read_excel(ROOT / "data" / "GSE334014" / "GSE334014_CountMatrix_ALL.xlsx")
    df = df.rename(columns={df.columns[0]: "gene"})
    df["gene"] = df["gene"].astype(str).str.upper()
    gx = df.set_index("gene")
    gx = gx.apply(pd.to_numeric, errors="coerce")
    sc = scores_from_gx(gx)
    tail = sc.index.to_series().astype(str).str.split("|").str[-1].str.strip()
    sc["compartment"] = np.where(
        tail.str.contains(r"PanCK", case=False),
        "PanCK",
        np.where(tail.str.contains(r"stroma", case=False), "stroma", "other"),
    )
    recs = [corr_block(sc, "all_AOI"), corr_block(sc[sc.compartment == "PanCK"], "PanCK"), corr_block(sc[sc.compartment == "stroma"], "stroma")]
    from scipy import stats as st

    a = sc[sc.compartment == "PanCK"]
    b = sc[sc.compartment == "stroma"]
    levels = {}
    for gene in ["CLDN4", "TACSTD2", "TB"]:
        u = st.mannwhitneyu(a[gene].dropna(), b[gene].dropna(), alternative="two-sided")
        levels[gene] = {
            "n_PanCK": int(a[gene].notna().sum()),
            "n_stroma": int(b[gene].notna().sum()),
            "median_PanCK": float(a[gene].median()),
            "median_stroma": float(b[gene].median()),
            "mw_p": float(u.pvalue),
        }
    return {"series": "GSE334014", "platform": "GeoMx WTA", "n_aoi": int(len(sc)), "corrs": recs, "compartment_levels": levels}


def run_gse326968():
    gx = load_gene_x_aoi(ROOT / "data" / "GSE326968" / "GSE326968_geomx_Q3_log2_normalized_matrix.csv.gz")
    sc = scores_from_gx(gx)
    recs = [corr_block(sc, "all_AOI")]
    return {"series": "GSE326968", "platform": "GeoMx WTA", "n_aoi": int(len(sc)), "corrs": recs, "note": "CLAD / allograft; no public CD45 segment key in the Q3 matrix"}


def flatten_corrs(payload):
    rows = []
    for c in payload["corrs"]:
        rec = {"series": payload["series"], "subset": c["subset"], "n": c["n"]}
        for k, v in c.items():
            if isinstance(v, dict) and "rho" in v:
                rec[f"{k}_n"] = v["n"]
                rec[f"{k}_rho"] = v["rho"]
                rec[f"{k}_p"] = v["p"]
        rows.append(rec)
    return rows


def main():
    (OUT / "tables").mkdir(parents=True, exist_ok=True)
    payloads = [
        run_gse265899(),
        run_gse289483(),
        run_gse334014(),
        run_gse326968(),
    ]
    rows = []
    for p in payloads:
        rows.extend(flatten_corrs(p))
    pd.DataFrame(rows).to_csv(OUT / "tables" / "geomx_correlations.csv", index=False)
    dump_json(OUT / "tables" / "geomx_summary.json", payloads)
    print("wrote geomx tables")


if __name__ == "__main__":
    main()
