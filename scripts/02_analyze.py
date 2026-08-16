"""GSE154826 / Leader 2021: TACSTD2 vs LCAM and CD8.

Uses the authors' published cell->cluster map (input_tables/cell_metadata.csv)
and the LCAM definitions from their own scripts:

  figure_5abcd_s5a.R  cell-frequency LCAM
      hi = T_activated + IgG + MoMac-II
      lo = B + AM + cDC2 + AZU1_mac + Tcm/naive_II + cDC1
      after within-lineage (T / B&plasma / MNP / lin_neg) normalisation
  get_LCAM_scores.R   bulk-RNA gene LCAM (applied here to immune-cell pbulk)

TACSTD2 protein is not on the CITE-seq panel. EPCAM / CD8 / CD45 ADTs are.
Most libraries are CD45+ bead or FACS enriched; epithelial cells sit in the
authors' "epi_endo_fibro_doublet" gate. Only two patients (695 LUAD, 706 LUSC)
were processed without CD45 enrichment (digest / dead-cell kit).
"""

from __future__ import annotations

import json
import os
import sys
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_gse154826 as L

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)

MIN_IMMUNE = 80
MIN_GATE = 15
EPS = 1e-2
CD8_CLUSTERS = ["CD8 Trm", "T_GZMK", "T_Nklike"]
CD4_CLUSTERS = ["CD4 Trm", "Treg", "Tcm/naive_I"]
EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
ENDO_GENES = ["PECAM1", "VWF", "CLDN5", "CDH5"]
FIBRO_GENES = ["COL1A1", "COL1A2", "DCN", "LUM"]


def _cpm(num, den):
    den = float(den)
    if den <= 0:
        return np.nan
    return 1e4 * float(num) / den


def load_annotated_cells():
    cols_always = [
        "barcode", "gex_umi", "n_genes", "mito_umi", "amp_batch_ID",
        "sample_ID", "cluster_ID", "annotated",
        "g_TACSTD2", "g_EPCAM", "g_PTPRC", "g_CD8A", "g_CD8B", "g_CD3D",
        "g_CD3E", "g_CD4", "g_KRT8", "g_KRT18", "g_KRT19", "g_KRT7",
        "g_CDH1", "g_CLDN4", "g_PECAM1", "g_VWF", "g_CLDN5", "g_CDH5",
        "g_COL1A1", "g_COL1A2", "g_DCN", "g_LUM", "g_NKG7", "g_MS4A1",
    ]
    frames = []
    for f in sorted(os.listdir(L.BATCHDIR)):
        if not (f.startswith("cellstats_") and f.endswith(".parquet")):
            continue
        path = os.path.join(L.BATCHDIR, f)
        raw = pd.read_parquet(path)
        raw = raw[raw.annotated].copy()
        keep = [c for c in cols_always if c in raw.columns]
        adt = [c for c in raw.columns if c.startswith("adt_")]
        frames.append(raw[keep + adt])
    cells = pd.concat(frames, ignore_index=True)
    ann = L.load_annots_list()
    geo = L.load_sample_annots()
    cells["sample_ID"] = cells["sample_ID"].astype(int)
    cells = cells.merge(ann, left_on="cluster_ID", right_on="cluster", how="left")
    cells = cells.merge(geo, on="sample_ID", how="left")
    cells["sub_lineage"] = cells["sub_lineage"].fillna("")
    cells["is_gate"] = cells["lineage"] == L.GATE_LINEAGE
    cells["is_immune"] = ~cells["is_gate"]
    cells["pt"] = cells["patient_ID"].astype(str) + " " + cells["tissue"]
    # epithelial-like score inside the gate (UMI sums)
    cells["epi_umi"] = cells[[f"g_{g}" for g in EPI_GENES]].sum(axis=1)
    cells["endo_umi"] = cells[[f"g_{g}" for g in ENDO_GENES]].sum(axis=1)
    cells["fibro_umi"] = cells[[f"g_{g}" for g in FIBRO_GENES]].sum(axis=1)
    cells["is_epi_like"] = cells.is_gate & (cells.epi_umi > cells.endo_umi) & (
        cells.epi_umi > cells.fibro_umi
    ) & (cells.epi_umi >= 1)
    return cells


def lineage_norm_table(cells, group_col):
    """Author LCAM table: drop gate, row-normalise, then within-lineage re-normalise."""
    imm = cells[cells.sub_lineage != ""]
    tab = pd.crosstab(imm[group_col], imm["sub_lineage"]).astype(float)
    tab = tab.div(tab.sum(axis=1), axis=0)
    ng = (
        L.load_annots_list()
        .drop_duplicates("sub_lineage")
        .set_index("sub_lineage")["norm_group_fig5"]
    )
    for gname in L.NORM_GROUPS:
        cols = [c for c in tab.columns if ng.get(c) == gname]
        if not cols:
            continue
        s = tab[cols].sum(axis=1).replace(0, np.nan)
        tab[cols] = tab[cols].div(s, axis=0)
    return tab


def lcam_from_table(tab):
    hi = [c for c in L.LCAM_HI_CLUSTERS if c in tab.columns]
    lo = [c for c in L.LCAM_LO_CLUSTERS if c in tab.columns]
    out = pd.DataFrame(index=tab.index)
    out["lcam_hi"] = np.log(tab[hi] + EPS).sum(axis=1)
    out["lcam_lo"] = np.log(tab[lo] + EPS).sum(axis=1)
    out["lcam_diff"] = out["lcam_hi"] - out["lcam_lo"]
    return out


def _first(s):
    return s.iloc[0]


def aggregate_unit(df):
    """One row of sample- or patient-level metrics from a cell frame."""
    imm = df[df.is_immune]
    gate = df[df.is_gate]
    epi = df[df.is_epi_like]
    t = df[df.lineage == "T"]
    row = {
        "n_cells": len(df),
        "n_immune": len(imm),
        "n_gate": len(gate),
        "n_epi_like": len(epi),
        "n_T": len(t),
        "frac_gate": len(gate) / max(len(df), 1),
        "frac_epi_like": len(epi) / max(len(df), 1),
        "tac_cpm_all": _cpm(df.g_TACSTD2.sum(), df.gex_umi.sum()),
        "tac_cpm_immune": _cpm(imm.g_TACSTD2.sum(), imm.gex_umi.sum()),
        "tac_cpm_gate": _cpm(gate.g_TACSTD2.sum(), gate.gex_umi.sum()),
        "tac_cpm_epi": _cpm(epi.g_TACSTD2.sum(), epi.gex_umi.sum()),
        "epcam_cpm_all": _cpm(df.g_EPCAM.sum(), df.gex_umi.sum()),
        "epcam_cpm_gate": _cpm(gate.g_EPCAM.sum(), gate.gex_umi.sum()),
        "cd8a_cpm_all": _cpm(df.g_CD8A.sum(), df.gex_umi.sum()),
        "cd8a_cpm_immune": _cpm(imm.g_CD8A.sum(), imm.gex_umi.sum()),
        "cd8a_cpm_T": _cpm(t.g_CD8A.sum(), t.gex_umi.sum()),
        "frac_tac_pos": (df.g_TACSTD2 > 0).mean(),
        "frac_tac_pos_gate": (gate.g_TACSTD2 > 0).mean() if len(gate) else np.nan,
        "cd8_trm_of_all": (df.sub_lineage == "CD8 Trm").mean(),
        "cd8_trm_of_T": (t.sub_lineage == "CD8 Trm").mean() if len(t) else np.nan,
        "cd8_broad_of_T": t.sub_lineage.isin(CD8_CLUSTERS).mean() if len(t) else np.nan,
        "t_activated_of_T": (t.sub_lineage == "T_activated").mean() if len(t) else np.nan,
        "igg_of_B": (
            (df.sub_lineage == "IgG").sum() / max((df.lineage == "B&plasma").sum(), 1)
        ),
        "patient_ID": _first(df.patient_ID),
        "tissue": _first(df.tissue),
        "disease": _first(df.disease),
        "prep_s1": _first(df.prep_s1),
        "prep": _first(df.prep),
        "library_chemistry": _first(df.library_chemistry),
        "Species": _first(df.Species),
    }
    adt_cols = [c for c in df.columns if c.startswith("adt_")]
    has = bool(adt_cols) and bool(df[adt_cols].notna().any().any())
    row["has_adt"] = has
    for name in ("EPCAM", "CD8", "CD45", "CD3", "CD4", "CD19"):
        col = f"adt_{name}"
        if has and col in df.columns and df[col].notna().any():
            row[f"adt_{name}_mean"] = float(df[col].mean())
            row[f"adt_{name}_imm"] = float(imm[col].mean()) if len(imm) else np.nan
            row[f"adt_{name}_gate"] = float(gate[col].mean()) if len(gate) else np.nan
        else:
            row[f"adt_{name}_mean"] = np.nan
            row[f"adt_{name}_imm"] = np.nan
            row[f"adt_{name}_gate"] = np.nan
    return pd.Series(row)


def bulk_lcam_from_pbulk():
    """Immune-cell pbulk LCAM gene score (get_LCAM_scores.R, per sample)."""
    want = list(dict.fromkeys(L.LCAM_BULK_GENES + ["TACSTD2"]))
    rows = []
    seen_genes = set()
    ann = L.load_annots_list().set_index("cluster")
    for f in sorted(os.listdir(L.BATCHDIR)):
        if not f.startswith("pbulk_"):
            continue
        z = np.load(os.path.join(L.BATCHDIR, f), allow_pickle=True)
        genes = z["genes"].astype(str)
        seen_genes.update(genes.tolist())
        gidx = {g: i for i, g in enumerate(genes) if g in want}
        is_imm = np.array(
            [ann.loc[int(c), "lineage"] != L.GATE_LINEAGE for c in z["cluster_ID"]]
        )
        counts = z["counts"]
        for sid in np.unique(z["sample_ID"]):
            m = (z["sample_ID"] == sid) & is_imm
            if not m.any():
                continue
            sub = counts[m]
            rows.append((
                int(sid),
                int(z["n_cells"][m].sum()),
                float(sub.sum()),
                {g: float(sub[:, i].sum()) for g, i in gidx.items()},
            ))
    if not rows:
        return pd.DataFrame()
    present = [g for g in L.LCAM_BULK_GENES if g in seen_genes]
    missing = [g for g in L.LCAM_BULK_GENES if g not in seen_genes]
    sids = [r[0] for r in rows]
    ncells = np.array([r[1] for r in rows])
    lib = np.array([r[2] for r in rows], float)
    Xg = np.array([[r[3].get(g, 0.0) for g in present] for r in rows], dtype=float)
    cpm = 1e6 * Xg / np.where(lib[:, None] == 0, np.nan, lib[:, None])
    zmat = np.log10(1e-6 + cpm)
    zmat = (zmat - np.nanmean(zmat, axis=0)) / (np.nanstd(zmat, axis=0) + 1e-12)
    subtypes = [L.LCAM_BULK_SUBTYPE[g] for g in present]
    uniq = list(dict.fromkeys(subtypes))
    sub_scores = np.zeros((len(sids), len(uniq)))
    for j, st in enumerate(uniq):
        idx = [k for k, s in enumerate(subtypes) if s == st]
        sub_scores[:, j] = zmat[:, idx].mean(axis=1)
    sub_scores = (sub_scores - sub_scores.mean(axis=0)) / (sub_scores.std(axis=0, ddof=0) + 1e-12)
    hi = [i for i, s in enumerate(uniq) if s in L.LCAM_HI_SUBTYPES]
    lo = [i for i, s in enumerate(uniq) if s in L.LCAM_LO_SUBTYPES]
    tac = np.array([r[3].get("TACSTD2", 0.0) for r in rows])
    out = pd.DataFrame({
        "sample_ID": sids,
        "n_immune_pbulk": ncells,
        "bulk_lcam_hi": sub_scores[:, hi].mean(axis=1),
        "bulk_lcam_lo": sub_scores[:, lo].mean(axis=1),
        "bulk_lcam_diff": sub_scores[:, hi].mean(axis=1) - sub_scores[:, lo].mean(axis=1),
        "tac_cpm_immune_pbulk": 1e4 * tac / np.where(lib == 0, np.nan, lib),
    })
    out.attrs["missing_lcam_genes"] = missing
    out.attrs["present_lcam_genes"] = present
    return out


def spearman_rows(df, pairs, cohort, extra=None):
    rows = []
    for x, y in pairs:
        sub = df[[x, y]].dropna()
        n = len(sub)
        if n < 5:
            r, p = np.nan, np.nan
        else:
            r, p = spearmanr(sub[x], sub[y])
        rec = {"cohort": cohort, "x": x, "y": y, "n": n, "rho": r, "p": p}
        if extra:
            rec.update(extra)
        rows.append(rec)
    return rows


def partial_spearman(df, x, y, z):
    """Spearman of residuals after rank-regressing x and y on z."""
    sub = df[[x, y, z]].dropna()
    n = len(sub)
    if n < 6:
        return n, np.nan, np.nan
    rx = sub[x].rank()
    ry = sub[y].rank()
    rz = sub[z].rank()
    # linear residual on ranks
    def resid(a, b):
        b = np.asarray(b, float)
        a = np.asarray(a, float)
        b = np.column_stack([np.ones(len(b)), b])
        coef, *_ = np.linalg.lstsq(b, a, rcond=None)
        return a - b @ coef
    r, p = spearmanr(resid(rx, rz), resid(ry, rz))
    return n, r, p


def save_fig(fig, name):
    path = os.path.join(L.FIGURES, name)
    fig.savefig(path, dpi=160, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def main():
    print("loading annotated cells", flush=True)
    cells = load_annotated_cells()
    print(f"  {len(cells)} cells, {cells.sample_ID.nunique()} samples, "
          f"{cells.patient_ID.nunique()} patients", flush=True)

    # lineage detectability
    det = (
        cells.groupby("lineage")
        .agg(
            n_cells=("barcode", "size"),
            tac_mean_umi=("g_TACSTD2", "mean"),
            tac_detect=("g_TACSTD2", lambda s: (s > 0).mean()),
            epcam_mean_umi=("g_EPCAM", "mean"),
            cd8a_mean_umi=("g_CD8A", "mean"),
        )
        .reset_index()
    )
    det.to_csv(os.path.join(L.TABLES, "tacstd2_by_lineage.tsv"), sep="\t", index=False)

    # sample-level
    print("aggregating samples", flush=True)
    sample = cells.groupby("sample_ID").apply(aggregate_unit, include_groups=False).reset_index()
    tab_s = lineage_norm_table(cells, "sample_ID")
    sample = sample.merge(lcam_from_table(tab_s), left_on="sample_ID", right_index=True, how="left")
    bulk = bulk_lcam_from_pbulk()
    missing_genes = list(bulk.attrs.get("missing_lcam_genes", []))
    present_genes = list(bulk.attrs.get("present_lcam_genes", []))
    sample = sample.merge(bulk, on="sample_ID", how="left")
    sample["v2_beads_tumor"] = (
        (sample.library_chemistry == "V2")
        & (sample.prep_s1 == "beads")
        & (sample.tissue == "Tumor")
        & (sample.n_immune >= MIN_IMMUNE)
    )
    sample["usable_lcam"] = sample.n_immune >= MIN_IMMUNE
    sample.to_csv(os.path.join(L.TABLES, "sample_metrics.tsv"), sep="\t", index=False)

    # patient x tissue (paper unit)
    print("aggregating patient-tissue", flush=True)
    pt = cells.groupby("pt").apply(aggregate_unit, include_groups=False).reset_index()
    tab_p = lineage_norm_table(cells, "pt")
    pt = pt.merge(lcam_from_table(tab_p), left_on="pt", right_index=True, how="left")
    # patient-level bulk LCAM: pool samples of the same patient-tissue
    geo = L.load_sample_annots()
    sid_to_pt = (
        geo.assign(pt=geo.patient_ID.astype(str) + " " + geo.tissue)
        .set_index("sample_ID")["pt"]
    )
    if len(bulk):
        b2 = bulk.copy()
        b2["pt"] = b2.sample_ID.map(sid_to_pt)
        # recompute from sample bulk scores as n-weighted mean (approximation);
        # better: already have sample-level. For patient, take n-weighted mean.
        w = b2.groupby("pt").apply(
            lambda d: pd.Series({
                "bulk_lcam_diff": np.average(d.bulk_lcam_diff, weights=d.n_immune_pbulk),
                "bulk_lcam_hi": np.average(d.bulk_lcam_hi, weights=d.n_immune_pbulk),
                "bulk_lcam_lo": np.average(d.bulk_lcam_lo, weights=d.n_immune_pbulk),
            }),
            include_groups=False,
        )
        pt = pt.merge(w, left_on="pt", right_index=True, how="left")
    pt["v2_beads_tumor"] = (
        (pt.library_chemistry == "V2")
        & (pt.prep_s1 == "beads")
        & (pt.tissue == "Tumor")
        & (pt.n_immune >= MIN_IMMUNE)
    )
    pt["usable_lcam"] = pt.n_immune >= MIN_IMMUNE
    pt.to_csv(os.path.join(L.TABLES, "patient_tissue_metrics.tsv"), sep="\t", index=False)

    # ---- correlations -------------------------------------------------------
    pairs = [
        ("tac_cpm_all", "lcam_diff"),
        ("tac_cpm_gate", "lcam_diff"),
        ("tac_cpm_epi", "lcam_diff"),
        ("tac_cpm_immune", "lcam_diff"),
        ("frac_gate", "lcam_diff"),
        ("tac_cpm_all", "cd8_trm_of_T"),
        ("tac_cpm_gate", "cd8_trm_of_T"),
        ("tac_cpm_epi", "cd8_trm_of_T"),
        ("tac_cpm_all", "cd8_broad_of_T"),
        ("tac_cpm_gate", "cd8_broad_of_T"),
        ("tac_cpm_epi", "cd8_broad_of_T"),
        ("tac_cpm_all", "cd8a_cpm_immune"),
        ("tac_cpm_gate", "cd8a_cpm_immune"),
        ("tac_cpm_epi", "cd8a_cpm_immune"),
        ("tac_cpm_all", "bulk_lcam_diff"),
        ("tac_cpm_gate", "bulk_lcam_diff"),
        ("frac_gate", "cd8_trm_of_T"),
        ("epcam_cpm_all", "lcam_diff"),
        ("epcam_cpm_all", "cd8_trm_of_T"),
    ]
    corr_rows = []

    def add_cohort(df, name, extra=None):
        d = df[df.usable_lcam].copy()
        corr_rows.extend(spearman_rows(d, pairs, name, extra))

    tum = pt[pt.tissue == "Tumor"]
    add_cohort(tum, "patient Tumor (all preps)")
    add_cohort(tum[tum.v2_beads_tumor], "patient Tumor V2 beads (paper LCAM set)")
    add_cohort(tum[tum.disease == "LUAD"], "patient Tumor LUAD")
    add_cohort(tum[tum.disease == "LUSC"], "patient Tumor LUSC")
    add_cohort(tum[tum.v2_beads_tumor & (tum.disease == "LUAD")], "patient Tumor V2 beads LUAD")
    add_cohort(tum[tum.v2_beads_tumor & (tum.disease == "LUSC")], "patient Tumor V2 beads LUSC")
    add_cohort(tum[tum.prep_s1 == "digest/dead cell"], "patient Tumor digest (n=2, underpowered)")
    add_cohort(pt[pt.tissue == "Normal"], "patient Normal (all preps)")
    add_cohort(
        pt[(pt.tissue == "Normal") & (pt.prep_s1 == "beads") & (pt.library_chemistry == "V2")],
        "patient Normal V2 beads",
    )
    # sample-level sensitivity
    st = sample[(sample.tissue == "Tumor") & sample.usable_lcam]
    add_cohort(st, "sample Tumor (all preps)")
    add_cohort(st[st.v2_beads_tumor], "sample Tumor V2 beads")
    add_cohort(st[st.has_adt == True], "sample Tumor with ADT (CITE-seq)")

    # CITE-seq ADT pairs
    adt_pairs = [
        ("tac_cpm_all", "adt_CD8_imm"),
        ("tac_cpm_gate", "adt_CD8_imm"),
        ("tac_cpm_epi", "adt_CD8_imm"),
        ("adt_EPCAM_mean", "lcam_diff"),
        ("adt_EPCAM_mean", "adt_CD8_imm"),
        ("adt_EPCAM_gate", "lcam_diff"),
        ("frac_gate", "adt_CD8_imm"),
    ]
    cite_pt = tum[(tum.has_adt == True) & tum.usable_lcam]
    corr_rows.extend(spearman_rows(cite_pt, adt_pairs, "patient Tumor CITE-seq ADT"))
    cite_s = st[st.has_adt == True]
    corr_rows.extend(spearman_rows(cite_s, adt_pairs, "sample Tumor CITE-seq ADT"))

    corr = pd.DataFrame(corr_rows)
    # BH within the primary cohort only; also global for inspection
    def bh(p):
        p = np.asarray(p, float)
        out = np.full_like(p, np.nan, dtype=float)
        ok = np.isfinite(p)
        if ok.sum() == 0:
            return out
        pv = p[ok]
        n = len(pv)
        order = np.argsort(pv)
        ranks = np.empty(n, int)
        ranks[order] = np.arange(1, n + 1)
        q = pv * n / ranks
        q_sorted = np.minimum.accumulate(q[order][::-1])[::-1]
        q_adj = np.empty(n)
        q_adj[order] = np.clip(q_sorted, 0, 1)
        out[ok] = q_adj
        return out

    corr["q_global"] = bh(corr["p"].values)
    corr.to_csv(os.path.join(L.TABLES, "spearman.tsv"), sep="\t", index=False)

    # partial: all-cell TACSTD2 vs LCAM / CD8 after residualising gate fraction
    partial_rows = []
    for cohort_name, d in [
        ("patient Tumor V2 beads (paper LCAM set)", tum[tum.v2_beads_tumor]),
        ("patient Tumor (all preps)", tum[tum.usable_lcam]),
        ("sample Tumor V2 beads", st[st.v2_beads_tumor]),
    ]:
        for x, y, z in [
            ("tac_cpm_all", "lcam_diff", "frac_gate"),
            ("tac_cpm_all", "cd8_trm_of_T", "frac_gate"),
            ("tac_cpm_gate", "lcam_diff", "frac_gate"),
            ("tac_cpm_epi", "lcam_diff", "frac_gate"),
        ]:
            n, r, p = partial_spearman(d, x, y, z)
            partial_rows.append({
                "cohort": cohort_name, "x": x, "y": y, "z": z,
                "n": n, "rho_partial": r, "p": p,
            })
    partial = pd.DataFrame(partial_rows)
    partial.to_csv(os.path.join(L.TABLES, "partial_spearman.tsv"), sep="\t", index=False)

    # cell-level CITE-seq sanity: TACSTD2 RNA vs CD8 ADT (expected negative: different lineages)
    cite_cells = cells[cells.columns.intersection(["g_TACSTD2", "gex_umi", "adt_CD8", "adt_EPCAM", "is_gate", "lineage", "tissue", "sample_ID"])]
    cell_level = []
    if "adt_CD8" in cells.columns:
        cc = cells[cells.adt_CD8.notna() & (cells.tissue == "Tumor")].copy()
        if len(cc):
            r, p = spearmanr(cc.g_TACSTD2, cc.adt_CD8)
            cell_level.append({
                "test": "cell Tumor TACSTD2 UMI vs CD8 ADT (lineage sanity)",
                "n": int(len(cc)), "rho": float(r), "p": float(p),
                "note": "Expected negative: TACSTD2 is epithelial, CD8 is T cells. Not a sample-level test.",
            })
            r, p = spearmanr(cc.g_TACSTD2, cc.adt_EPCAM) if "adt_EPCAM" in cc.columns else (np.nan, np.nan)
            cell_level.append({
                "test": "cell Tumor TACSTD2 UMI vs EPCAM ADT",
                "n": int(len(cc)), "rho": float(r) if np.isfinite(r) else None,
                "p": float(p) if np.isfinite(p) else None,
                "note": "Same-lineage positive control if Trop-2 RNA tracks EpCAM protein.",
            })
    pd.DataFrame(cell_level).to_csv(os.path.join(L.TABLES, "cite_cell_level_sanity.tsv"), sep="\t", index=False)

    # digest two-patient descriptive
    digest = tum[tum.prep_s1 == "digest/dead cell"][
        ["pt", "patient_ID", "disease", "n_cells", "n_gate", "n_epi_like",
         "tac_cpm_all", "tac_cpm_gate", "tac_cpm_epi", "lcam_diff",
         "cd8_trm_of_T", "cd8_broad_of_T", "cd8a_cpm_immune", "frac_gate"]
    ]
    digest.to_csv(os.path.join(L.TABLES, "digest_two_patients.tsv"), sep="\t", index=False)

    # ---- key stats ----------------------------------------------------------
    def grab(cohort, x, y):
        hit = corr[(corr.cohort == cohort) & (corr.x == x) & (corr.y == y)]
        if hit.empty:
            return {}
        r = hit.iloc[0]
        return {"n": int(r.n), "rho": None if pd.isna(r.rho) else float(r.rho),
                "p": None if pd.isna(r.p) else float(r.p)}

    primary = "patient Tumor V2 beads (paper LCAM set)"
    key = {
        "dataset": "GSE154826 / Leader et al. Cancer Cell 2021",
        "n_annotated_cells": int(len(cells)),
        "n_samples_annotated": int(cells.sample_ID.nunique()),
        "n_patients_annotated": int(cells.patient_ID.nunique()),
        "n_patient_tumor": int((pt.tissue == "Tumor").sum()),
        "n_patient_tumor_v2_beads": int(tum.v2_beads_tumor.sum()),
        "n_digest_tumor_patients": int((tum.prep_s1 == "digest/dead cell").sum()),
        "trop2_adt_on_panel": False,
        "epcam_cd8_cd45_adt": True,
        "cd45_enrichment": "61/119 GEO libraries CD45+ beads; 10 CD45+ FACS; 22 digest/dead-cell (2 patients); 4 CD2+",
        "lcam_hi_clusters": L.LCAM_HI_CLUSTERS,
        "lcam_lo_clusters": L.LCAM_LO_CLUSTERS,
        "missing_bulk_lcam_genes": missing_genes,
        "present_bulk_lcam_genes": present_genes,
        "primary_cohort": primary,
        "primary_tac_gate_vs_lcam": grab(primary, "tac_cpm_gate", "lcam_diff"),
        "primary_tac_epi_vs_lcam": grab(primary, "tac_cpm_epi", "lcam_diff"),
        "primary_tac_all_vs_lcam": grab(primary, "tac_cpm_all", "lcam_diff"),
        "primary_tac_gate_vs_cd8trm": grab(primary, "tac_cpm_gate", "cd8_trm_of_T"),
        "primary_tac_all_vs_cd8trm": grab(primary, "tac_cpm_all", "cd8_trm_of_T"),
        "primary_frac_gate_vs_lcam": grab(primary, "frac_gate", "lcam_diff"),
        "all_tumor_tac_gate_vs_lcam": grab("patient Tumor (all preps)", "tac_cpm_gate", "lcam_diff"),
        "all_tumor_tac_all_vs_lcam": grab("patient Tumor (all preps)", "tac_cpm_all", "lcam_diff"),
        "cite_patient_tac_gate_vs_cd8adt": grab("patient Tumor CITE-seq ADT", "tac_cpm_gate", "adt_CD8_imm"),
        "cite_patient_tac_all_vs_cd8adt": grab("patient Tumor CITE-seq ADT", "tac_cpm_all", "adt_CD8_imm"),
        "cite_patient_frac_gate_vs_cd8adt": grab("patient Tumor CITE-seq ADT", "frac_gate", "adt_CD8_imm"),
        "cite_sample_tac_all_vs_cd8adt": grab("sample Tumor CITE-seq ADT", "tac_cpm_all", "adt_CD8_imm"),
        "cite_sample_tac_gate_vs_cd8adt": grab("sample Tumor CITE-seq ADT", "tac_cpm_gate", "adt_CD8_imm"),
        "cite_sample_frac_gate_vs_cd8adt": grab("sample Tumor CITE-seq ADT", "frac_gate", "adt_CD8_imm"),
        "cite_patient_epcamadt_vs_lcam": grab("patient Tumor CITE-seq ADT", "adt_EPCAM_mean", "lcam_diff"),
        "user_expected": "negative TACSTD2 vs LCAM / CD8",
        "verdict_note": (
            "CD45+ CITE-seq. TACSTD2 protein not measured. "
            "Primary test is TACSTD2 RNA inside the authors' epi/endo/fibro gate "
            "versus the authors' own lineage-normalised LCAM difference, "
            "patient-level V2 beads Tumor."
        ),
    }
    # attach partials
    key["partial"] = partial.to_dict(orient="records")
    key["digest_patients"] = digest.to_dict(orient="records")
    key["cell_level_cite"] = cell_level
    with open(os.path.join(L.TABLES, "key_stats.json"), "w") as fh:
        json.dump(key, fh, indent=2, default=str)

    # ---- figures ------------------------------------------------------------
    print("figures", flush=True)
    prim = tum[tum.v2_beads_tumor].copy()
    prim["hist"] = prim.disease.fillna("NA")

    def scatter(ax, d, x, y, title, xlab, ylab):
        colors = {"LUAD": "#1f77b4", "LUSC": "#d62728"}
        plotted = d[[x, y, "hist"]].dropna(subset=[x, y])
        for h, g in plotted.groupby("hist"):
            ax.scatter(g[x], g[y], s=36, alpha=0.85, c=colors.get(h, "#888"),
                       edgecolors="k", linewidths=0.3, label=f"{h} n={len(g)}")
        sub = plotted[[x, y]]
        if len(sub) >= 5:
            r, p = spearmanr(sub[x], sub[y])
            ax.set_title(f"{title}\nSpearman ρ={r:.2f}  p={p:.3g}  n={len(sub)}", fontsize=10)
        else:
            ax.set_title(f"{title}\nn={len(sub)} (no Spearman)", fontsize=10)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        if len(plotted):
            ax.legend(frameon=False, fontsize=8)
        ax.axhline(0, color="#aaa", lw=0.6)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0))
    scatter(axes[0], prim, "tac_cpm_gate", "lcam_diff",
            "Primary: gate TACSTD2 vs LCAM",
            "TACSTD2 CPM in epi/endo/fibro gate", "LCAM difference (hi − lo)")
    scatter(axes[1], prim, "tac_cpm_all", "lcam_diff",
            "All-cell TACSTD2 vs LCAM (composition-mixed)",
            "TACSTD2 CPM, all annotated cells", "LCAM difference (hi − lo)")
    scatter(axes[2], prim, "frac_gate", "lcam_diff",
            "Gate fraction vs LCAM (confound check)",
            "Fraction of cells in epi/endo/fibro gate", "LCAM difference (hi − lo)")
    fig.tight_layout()
    save_fig(fig, "fig1_primary_tacstd2_vs_lcam.png")

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.0))
    scatter(axes[0], prim, "tac_cpm_gate", "cd8_trm_of_T",
            "Gate TACSTD2 vs CD8 Trm / T",
            "TACSTD2 CPM in gate", "CD8 Trm fraction of T cells")
    scatter(axes[1], prim, "tac_cpm_epi", "cd8_broad_of_T",
            "Epi-like TACSTD2 vs broad CD8 / T",
            "TACSTD2 CPM in epi-like gate cells", "CD8 Trm+GZMK+NKlike / T")
    scatter(axes[2], prim, "tac_cpm_gate", "cd8a_cpm_immune",
            "Gate TACSTD2 vs immune CD8A",
            "TACSTD2 CPM in gate", "CD8A CPM in immune cells")
    fig.tight_layout()
    save_fig(fig, "fig2_tacstd2_vs_cd8.png")

    # forest of primary-cohort rhos
    forest_keys = [
        ("tac_cpm_gate", "lcam_diff", "gate TACSTD2 vs LCAM"),
        ("tac_cpm_epi", "lcam_diff", "epi-like TACSTD2 vs LCAM"),
        ("tac_cpm_all", "lcam_diff", "all-cell TACSTD2 vs LCAM"),
        ("tac_cpm_immune", "lcam_diff", "immune TACSTD2 vs LCAM"),
        ("frac_gate", "lcam_diff", "gate fraction vs LCAM"),
        ("tac_cpm_gate", "cd8_trm_of_T", "gate TACSTD2 vs CD8 Trm/T"),
        ("tac_cpm_all", "cd8_trm_of_T", "all-cell TACSTD2 vs CD8 Trm/T"),
        ("tac_cpm_gate", "cd8_broad_of_T", "gate TACSTD2 vs broad CD8/T"),
        ("tac_cpm_gate", "cd8a_cpm_immune", "gate TACSTD2 vs immune CD8A"),
        ("epcam_cpm_all", "lcam_diff", "all-cell EPCAM vs LCAM"),
    ]
    fdf = []
    for x, y, lab in forest_keys:
        hit = corr[(corr.cohort == primary) & (corr.x == x) & (corr.y == y)]
        if hit.empty:
            continue
        r = hit.iloc[0]
        fdf.append({"lab": lab, "rho": r.rho, "p": r.p, "n": r.n})
    fdf = pd.DataFrame(fdf).iloc[::-1]
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    yidx = np.arange(len(fdf))
    cols = ["#b2182b" if (pd.notna(r) and r < 0) else "#2166ac" for r in fdf.rho]
    ax.axvline(0, color="#444", lw=0.8)
    ax.scatter(fdf.rho, yidx, c=cols, s=42, zorder=3, edgecolors="k", linewidths=0.3)
    for i, row in enumerate(fdf.itertuples()):
        if pd.isna(row.rho):
            continue
        ax.text(row.rho + (0.03 if row.rho >= 0 else -0.03), i,
                f"ρ={row.rho:.2f} p={row.p:.2g}",
                va="center", ha="left" if row.rho >= 0 else "right", fontsize=8)
    ax.set_yticks(yidx)
    ax.set_yticklabels(fdf.lab)
    ax.set_xlabel(f"Spearman ρ  ·  {primary}  ·  n≈{int(prim.shape[0])}")
    ax.set_xlim(-1.05, 1.05)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_title("User expected negative. Red = negative ρ.")
    fig.tight_layout()
    save_fig(fig, "fig3_forest_primary.png")

    # TACSTD2 by lineage
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    order = ["epi_endo_fibro_doublet", "MNP", "B&plasma", "T", "NK", "mast", "pDC"]
    dat, labs = [], []
    for lin in order:
        v = cells.loc[cells.lineage == lin, "g_TACSTD2"]
        dat.append(np.log1p(v.values))
        labs.append(f"{lin}\nn={len(v):,}")
    ax.boxplot(dat, tick_labels=labs, showfliers=False, medianprops={"color": "#b2182b"})
    ax.set_ylabel("log1p TACSTD2 UMI / cell")
    ax.set_title("TACSTD2 is almost restricted to the authors' epi/endo/fibro gate")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save_fig(fig, "fig4_tacstd2_by_lineage.png")

    # all tumor patients: beads vs digest
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    all_t = tum[tum.usable_lcam].copy()
    all_t["hist"] = all_t.disease
    markers = {"beads": "o", "digest/dead cell": "D", "sort": "s", "CD2": "^"}
    colors = {"LUAD": "#1f77b4", "LUSC": "#d62728"}
    for (prep, h), g in all_t.groupby(["prep_s1", "hist"]):
        ax.scatter(g.tac_cpm_gate, g.lcam_diff, s=50 if prep != "beads" else 32,
                   marker=markers.get(prep, "o"), c=colors.get(h, "#888"),
                   edgecolors="k", linewidths=0.3, alpha=0.85,
                   label=f"{h} {prep} n={len(g)}")
    sub = all_t[["tac_cpm_gate", "lcam_diff"]].dropna()
    r, p = spearmanr(sub.tac_cpm_gate, sub.lcam_diff)
    ax.set_title(f"All tumor patients, gate TACSTD2 vs LCAM\nρ={r:.2f} p={p:.3g} n={len(sub)}")
    ax.set_xlabel("TACSTD2 CPM in epi/endo/fibro gate")
    ax.set_ylabel("LCAM difference (hi − lo)")
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    save_fig(fig, "fig5_all_tumor_preps.png")

    # CITE-seq ADT
    if len(cite_pt):
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
        d = cite_pt.copy()
        d["hist"] = d.disease
        scatter(axes[0], d, "tac_cpm_gate", "adt_CD8_imm",
                "CITE-seq tumor: gate TACSTD2 vs CD8 ADT",
                "TACSTD2 CPM in gate", "Mean CD8 ADT (immune cells)")
        scatter(axes[1], d, "tac_cpm_all", "adt_CD8_imm",
                "CITE-seq: all-cell TACSTD2 vs CD8 ADT\n(composition-mixed)",
                "TACSTD2 CPM, all annotated cells", "Mean CD8 ADT (immune cells)")
        fig.tight_layout()
        save_fig(fig, "fig6_cite_adt.png")

        # composition artifact: gate fraction vs CD8 protein
        fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
        cite_s2 = cite_s.copy()
        cite_s2["hist"] = cite_s2.disease
        scatter(axes[0], cite_s2, "frac_gate", "adt_CD8_imm",
                "Sample CITE-seq: gate fraction vs CD8 ADT",
                "Fraction of cells in epi/endo/fibro gate", "Mean CD8 ADT (immune cells)")
        scatter(axes[1], cite_s2, "tac_cpm_gate", "adt_CD8_imm",
                "Sample CITE-seq: gate TACSTD2 vs CD8 ADT",
                "TACSTD2 CPM in gate", "Mean CD8 ADT (immune cells)")
        fig.tight_layout()
        save_fig(fig, "fig7_cite_composition.png")

    # coverage / design schematic table
    cov = pd.DataFrame({
        "item": [
            "GEO accession",
            "Paper",
            "Annotated cells used",
            "Patients with annotated cells",
            "Tumor patient-tissue units",
            "V2 beads Tumor (paper LCAM set)",
            "Digest/dead-cell Tumor patients",
            "Trop-2 / TACSTD2 ADT",
            "EPCAM / CD8 / CD45 ADT",
            "CD45 enrichment",
            "LCAM cell-frequency definition",
            "User expected sign",
        ],
        "value": [
            "GSE154826",
            "Leader, Grout et al. Cancer Cell 2021",
            f"{len(cells):,}",
            str(cells.patient_ID.nunique()),
            str(int((pt.tissue == 'Tumor').sum())),
            str(int(tum.v2_beads_tumor.sum())),
            "2 (695 LUAD, 706 LUSC) — correlation not computed",
            "NOT on the panel",
            "Yes, in hashed CITE-seq libraries",
            "Yes for beads/sort; no for digest (2 patients)",
            "T_activated + IgG + MoMac-II  minus  B+AM+cDC2+AZU1_mac+Tcm/naive_II+cDC1",
            "negative",
        ],
    })
    cov.to_csv(os.path.join(L.TABLES, "design_coverage.tsv"), sep="\t", index=False)

    focus = corr[corr.cohort.isin([
        primary,
        "patient Tumor (all preps)",
        "patient Tumor CITE-seq ADT",
        "sample Tumor CITE-seq ADT",
    ]) & corr.x.isin([
        "tac_cpm_gate", "tac_cpm_epi", "tac_cpm_all", "frac_gate",
    ]) & corr.y.isin([
        "lcam_diff", "cd8_trm_of_T", "cd8_broad_of_T", "cd8a_cpm_immune",
        "adt_CD8_imm", "bulk_lcam_diff",
    ])].copy()
    focus.to_csv(os.path.join(L.TABLES, "primary_focus.tsv"), sep="\t", index=False)

    print("wrote", L.TABLES)
    print("primary gate TACSTD2 vs LCAM", key["primary_tac_gate_vs_lcam"])
    print("primary all-cell TACSTD2 vs LCAM", key["primary_tac_all_vs_lcam"])
    print("primary gate TACSTD2 vs CD8 Trm", key["primary_tac_gate_vs_cd8trm"])


if __name__ == "__main__":
    main()
