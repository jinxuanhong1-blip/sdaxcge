"""A11 single-cell analyses on GSE131907 (Kim et al. 2020, LUAD).

Patient is the unit of analysis. Site is not mixed in the primary test:
tLung tS1/tS2/tS3 are primary-tumour epithelium; 'Malignant cells' in
Kim's annotation are metastatic sites only. Pooling them with brain mets
inflates CD47–TACSTD2 because both genes drop at that site.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as C
from stats_utils import spearman, bh_fdr

DERIVED = C.DATA / "derived"
TAB = C.OUT / "tables"
TAB.mkdir(parents=True, exist_ok=True)

MALIGNANT = {"Malignant cells", "tS1", "tS2", "tS3"}
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE"}


def load():
    mat = pd.read_parquet(DERIVED / "sc_GSE131907_panel.parquet")
    ann = pd.read_parquet(DERIVED / "sc_GSE131907_annot.parquet")
    return mat, ann


def s1_celltype(mat, ann):
    tumor = ann[ann.Sample_Origin.isin(TUMOR_ORIGINS)]
    genes = [g for g in [C.TARGET] + C.primary_panel() if g in mat.columns]
    rows = []
    for ct, idx in tumor.groupby("Cell_type").groups.items():
        sub = mat.loc[idx, genes]
        n = len(idx)
        n_pat = tumor.loc[idx, "Sample"].nunique()
        for g in genes:
            v = sub[g].values
            rows.append(dict(
                cell_type=ct, gene=g, n_cells=n, n_samples=n_pat,
                mean_log2tpm=float(v.mean()),
                pct_expr=float((v > 0).mean() * 100),
                median=float(np.median(v)),
            ))
    return pd.DataFrame(rows)


def _pseudobulk(mat, ann, mask, min_cells):
    sub = ann.loc[mask]
    genes = [g for g in [C.TARGET] + C.primary_panel() if g in mat.columns]
    pb = mat.loc[sub.index, genes].groupby(sub["Sample"]).mean()
    n_cells = sub.groupby("Sample").size().rename("n_cells")
    origin = sub.groupby("Sample")["Sample_Origin"].agg(
        lambda s: s.value_counts().index[0])
    pb = pb.join(n_cells).join(origin)
    return pb[pb.n_cells >= min_cells]


def s2_slices(mat, ann):
    """Within-malignant-cell TACSTD2 correlations, one row per gene x slice."""
    slices = {
        "primary_tLung_tS": (
            ann.Cell_subtype.isin(MALIGNANT) & ann.Sample_Origin.eq("tLung"), 20),
        "mets_mLN_tLB": (
            ann.Cell_subtype.eq("Malignant cells")
            & ann.Sample_Origin.isin({"mLN", "tL/B"}), 20),
        "brain_mets": (
            ann.Cell_subtype.eq("Malignant cells") & ann.Sample_Origin.eq("mBrain"), 20),
        "tumor_sites_mixed_no_brain": (
            ann.Cell_subtype.isin(MALIGNANT)
            & ann.Sample_Origin.isin({"tLung", "tL/B", "mLN", "PE"}), 30),
    }
    rows, pbs = [], []
    for name, (mask, min_c) in slices.items():
        pb = _pseudobulk(mat, ann, mask, min_c)
        pb = pb.assign(slice=name)
        pbs.append(pb.reset_index())
        tac = pb[C.TARGET]
        for g in [c for c in pb.columns if c not in
                  (C.TARGET, "n_cells", "Sample_Origin", "slice")]:
            r = spearman(tac, pb[g])
            rows.append(dict(
                slice=name, gene=g, axis=C.gene_to_axis().get(g, ""),
                rho=r["rho"], lo=r["lo"], hi=r["hi"], p=r["p"], n=r["n"],
                n_cells_total=int(pb.n_cells.sum()),
            ))
    df = pd.DataFrame(rows)
    # FDR within each slice, primary panel only.
    df["q"] = np.nan
    for sl in df.slice.unique():
        m = df.slice.eq(sl)
        df.loc[m, "q"] = bh_fdr(df.loc[m, "p"].values)
    return df, pd.concat(pbs, ignore_index=True)


def s3_composition(mat, ann):
    tumor = ann[ann.Sample_Origin.eq("tLung")]
    frac = (tumor.groupby(["Sample", "Cell_type"]).size()
            .unstack(fill_value=0))
    frac = frac.div(frac.sum(axis=1), axis=0)

    mal = ann[ann.Cell_subtype.isin(MALIGNANT) & ann.Sample_Origin.eq("tLung")]
    tac = (mat.loc[mal.index, C.TARGET].groupby(mal["Sample"]).mean()
           .rename("TACSTD2_mal"))
    n_mal = mal.groupby("Sample").size()
    keep = n_mal[n_mal >= 20].index
    tac = tac.reindex(keep)
    frac = frac.reindex(keep)

    rows = []
    for ct in frac.columns:
        r = spearman(tac, frac[ct])
        rows.append(dict(cell_type=ct, rho=r["rho"], lo=r["lo"], hi=r["hi"],
                         p=r["p"], n=r["n"]))
    df = pd.DataFrame(rows)
    df["q"] = bh_fdr(df["p"].values)
    return df


def main():
    mat, ann = load()
    print("sc matrix", mat.shape)
    s1 = s1_celltype(mat, ann)
    s1.to_csv(TAB / "20_sc_celltype_expression.csv", index=False)

    s2, pb = s2_slices(mat, ann)
    s2.to_csv(TAB / "21_sc_malignant_pseudobulk.csv", index=False)
    pb.to_csv(TAB / "21_sc_malignant_pseudobulk_samples.csv", index=False)

    s3 = s3_composition(mat, ann)
    s3.to_csv(TAB / "22_sc_tme_composition.csv", index=False)

    focus = s2[s2.gene.isin(["CD47", "SIRPA", "THBS1", "NECTIN4", "LGALS3"])]
    print(focus.to_string(index=False))
    print("\nS1 TACSTD2/CD47")
    print(s1[s1.gene.isin(["TACSTD2", "CD47"])]
          .sort_values(["gene", "mean_log2tpm"], ascending=[True, False])
          .to_string(index=False))
    print("wrote sc tables")


if __name__ == "__main__":
    main()
