#!/usr/bin/env python3
"""UCSC Xena TCGA-LUAD: TACSTD2 and CLDN4 vs immune features after
ESTIMATE and ABSOLUTE tumor-purity adjustment.

Primary matrix: Xena legacy HiSeqV2 (log2(RSEM norm_count+1)), which is the
same RNAseqV2 freeze MD Anderson used for published ESTIMATE scores.
Sensitivity matrix: Xena GDC STAR log2(TPM+1).

Outputs -> results/w200/Xena_LUAD/
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import warnings

import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.join(ROOT, "data", "w200", "raw")
OUT = os.path.join(ROOT, "results", "w200", "Xena_LUAD")
TABLES = os.path.join(OUT, "tables")
FIG = os.path.join(OUT, "figures")
os.makedirs(TABLES, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

# Yoshihara et al. 2013 Nature Communications: ESTIMATE score -> purity.
# Affymetrix-calibrated; applied here to RNAseqV2 scores as MD Anderson does.
ESTIMATE_A = 0.6049872018
ESTIMATE_B = 0.0001467884

PREDICTORS = ["TACSTD2", "CLDN4"]

MARKER_GENES = [
    "CD8A", "CD8B", "CD3E", "CD2",
    "GZMA", "GZMB", "PRF1", "NKG7", "IFNG",
    "CXCL9", "CXCL10",
    "CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "IDO1",
    "FOXP3", "MS4A1", "CD68",
]

# Ayers et al. 2017 JCI: 18-gene T-cell-inflamed GEP (pembrolizumab).
GEP18 = [
    "CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6",
    "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2",
    "PSMB10", "STAT1", "TIGIT",
]

WOLF_SIGNATURES = {
    "LIexpression_score": "Wolf_LymphocyteInfiltration",
    "IFNG_score_21050467": "Wolf_IFNgamma",
    "TGFB_score_21050467": "Wolf_TGFbeta",
    "CSF1_response": "Wolf_MacrophageCSF1",
    "CHANG_CORE_SERUM_RESPONSE_UP": "Wolf_WoundHealing_CSR",
}

CIBERSORT_CELLS = {
    "T.cells.CD8": "CIBERSORT_CD8_T",
    "T.cells.regulatory..Tregs.": "CIBERSORT_Tregs",
    "NK.cells.activated": "CIBERSORT_NK_activated",
    "Macrophages.M1": "CIBERSORT_Macrophage_M1",
    "Macrophages.M2": "CIBERSORT_Macrophage_M2",
}

# Immune_score is a term in ESTIMATE_score; adjusting Immune_score for
# ESTIMATE purity is circular by construction. Still reported, flagged.
CIRCULAR_WITH_ESTIMATE = {"ESTIMATE_Immune_score", "ESTIMATE_Stromal_score"}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sample01(barcode: str) -> str | None:
    """Any TCGA barcode -> patient-level primary-tumor id (TCGA-XX-XXXX-01), or None."""
    b = barcode.replace(".", "-")
    parts = b.split("-")
    if len(parts) < 4 or not parts[3].startswith("01"):
        return None
    return "-".join(parts[:3]) + "-01"


def estimate_purity_from_score(score: np.ndarray | pd.Series) -> np.ndarray:
    pur = np.cos(ESTIMATE_A + ESTIMATE_B * np.asarray(score, dtype=float))
    return np.clip(pur, 0.0, 1.0)


def partial_spearman(x, y, *covariates):
    """Spearman of x,y controlling for one or more covariates (Pearson on rank residuals)."""
    xr = stats.rankdata(x).astype(float)
    yr = stats.rankdata(y).astype(float)
    zcols = [stats.rankdata(z).astype(float) for z in covariates]
    zc = np.column_stack([np.ones(len(xr))] + zcols)
    bx, *_ = np.linalg.lstsq(zc, xr, rcond=None)
    by, *_ = np.linalg.lstsq(zc, yr, rcond=None)
    rx = xr - zc @ bx
    ry = yr - zc @ by
    if np.std(rx) < 1e-12 or np.std(ry) < 1e-12:
        return np.nan, np.nan
    r, _ = stats.pearsonr(rx, ry)
    n = len(x)
    df = n - 2 - len(covariates)
    if df <= 0:
        return float(r), np.nan
    t = r * np.sqrt(df / max(1e-12, 1 - r * r))
    p = 2 * stats.t.sf(abs(t), df)
    return float(r), float(p)


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    out = np.empty(n)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, pv[idx] * n / rank)
        out[idx] = val
        prev = val
    q[ok] = out
    return q


def spearman_row(x, y):
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10:
        return np.nan, np.nan, int(mask.sum())
    r, p = stats.spearmanr(x[mask], y[mask])
    return float(r), float(p), int(mask.sum())


def load_hiseqv2(genes: list[str]) -> pd.DataFrame:
    path = os.path.join(RAW, "TCGA.LUAD.HiSeqV2.gz")
    expr = pd.read_csv(path, sep="\t", index_col=0)
    keep = [c for c in expr.columns if c.endswith("-01")]
    expr = expr.loc[:, keep]
    missing = [g for g in genes if g not in expr.index]
    if missing:
        raise SystemExit(f"HiSeqV2 missing genes: {missing}")
    # one column per -01 already; average if any exact dups
    out = expr.loc[genes].T.groupby(level=0).mean()
    out.index.name = "sample"
    return out


def load_star_tpm(genes: list[str]) -> pd.DataFrame:
    pm = pd.read_csv(os.path.join(RAW, "gencode.v36.annotation.gtf.gene.probemap"), sep="\t")
    sub = pm[pm["gene"].isin(genes)][["id", "gene"]]
    found = set(sub["gene"])
    missing = [g for g in genes if g not in found]
    if missing:
        raise SystemExit(f"probemap missing genes: {missing}")
    # if a symbol has multiple Ensembl rows, keep the first (canonical gencode id)
    id_to_gene = {}
    for _, row in sub.iterrows():
        id_to_gene.setdefault(row["id"], row["gene"])
    want = set(id_to_gene)
    collected = {}
    path = os.path.join(RAW, "TCGA-LUAD.star_tpm.tsv.gz")
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        samples = header[1:]
        for line in f:
            gid = line.split("\t", 1)[0]
            if gid not in want:
                continue
            vals = np.array(line.rstrip("\n").split("\t")[1:], dtype=float)
            collected[id_to_gene[gid]] = vals
    missing2 = [g for g in genes if g not in collected]
    if missing2:
        raise SystemExit(f"STAR TPM missing genes: {missing2}")
    mat = pd.DataFrame(collected, index=samples)
    # primary tumors only; collapse 01A/01B to patient -01 (prefer 01A by averaging)
    sid = [sample01(s) for s in mat.index]
    mat = mat.loc[[s is not None for s in sid]].copy()
    mat.index = [sample01(s) for s in mat.index]
    mat = mat.groupby(level=0).mean()
    mat.index.name = "sample"
    return mat


def add_derived(eg: pd.DataFrame) -> pd.DataFrame:
    out = eg.copy()
    out["CYT_GZMA_PRF1"] = out[["GZMA", "PRF1"]].mean(axis=1)
    out["CD8_score"] = out[["CD8A", "CD8B"]].mean(axis=1)
    gep = [g for g in GEP18 if g in out.columns]
    out["GEP18"] = out[gep].mean(axis=1)
    return out


def load_purity() -> pd.DataFrame:
    est = pd.read_csv(os.path.join(RAW, "MDACC_estimate_LUAD_RNAseqV2.txt"), sep="\t")
    est["sample"] = est["ID"].map(sample01)
    est = est.dropna(subset=["sample"]).groupby("sample").mean(numeric_only=True)
    est["estimate_purity"] = estimate_purity_from_score(est["ESTIMATE_score"])
    est = est.rename(columns={
        "Immune_score": "ESTIMATE_Immune_score",
        "Stromal_score": "ESTIMATE_Stromal_score",
        "ESTIMATE_score": "ESTIMATE_score",
    })

    ab = pd.read_csv(os.path.join(RAW, "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"), sep="\t")
    ab = ab[ab["array"].str.endswith("-01", na=False)][["array", "purity"]].dropna()
    ab["sample"] = ab["array"].map(sample01)
    abs_p = ab.dropna(subset=["sample"]).groupby("sample")["purity"].mean().rename("absolute_purity")

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unknown extension")
        aran = pd.read_excel(
            os.path.join(RAW, "Aran2015_purity_S1.xlsx"),
            sheet_name=0, header=3,
        )
    aran = aran.rename(columns={"Sample ID": "sid", "Cancer type": "ctype"})
    aran = aran[aran["ctype"] == "LUAD"].copy()
    aran["sample"] = aran["sid"].map(sample01)
    aran = aran.dropna(subset=["sample"]).groupby("sample")[["ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"]].mean()
    aran = aran.rename(columns={
        "ESTIMATE": "aran_estimate_purity",
        "ABSOLUTE": "aran_absolute_purity",
        "LUMP": "aran_lump_purity",
        "IHC": "aran_ihc_purity",
        "CPE": "aran_cpe_purity",
    })

    pur = est.join(abs_p, how="outer").join(aran, how="left")
    pur.index.name = "sample"
    return pur


def load_leuk() -> pd.Series:
    lf = pd.read_csv(
        os.path.join(RAW, "TCGA_all_leuk_estimate.masked.20170107.tsv"),
        sep="\t", header=None, names=["cancer", "aliquot", "leuk"],
    )
    lf = lf[lf["cancer"] == "LUAD"].copy()
    lf["sample"] = lf["aliquot"].map(sample01)
    return lf.dropna(subset=["sample"]).groupby("sample")["leuk"].mean().rename(
        "LeukocyteFraction_methylation"
    )


def load_cibersort() -> pd.DataFrame:
    cs = pd.read_csv(os.path.join(RAW, "TCGA.Kallisto.fullIDs.cibersort.relative.tsv"), sep="\t")
    cs = cs[cs["CancerType"] == "LUAD"].copy()
    cs["sample"] = cs["SampleID"].map(sample01)
    out = cs.dropna(subset=["sample"]).groupby("sample")[list(CIBERSORT_CELLS)].mean()
    out = out.rename(columns=CIBERSORT_CELLS)
    return out


def load_wolf(valid_samples: pd.Index) -> pd.DataFrame:
    sig = pd.read_csv(os.path.join(RAW, "Scores_160_Signatures.tsv.gz"), sep="\t")
    sig = sig[sig["SetName"].isin(WOLF_SIGNATURES)].set_index("SetName")
    sig = sig.drop(columns=["Source"]).T
    sig.index = [sample01(b) or "" for b in sig.index]
    sig = sig[sig.index != ""]
    sig = sig.loc[sig.index.isin(valid_samples)]
    sig = sig.astype(float).groupby(level=0).mean()
    sig.columns = [WOLF_SIGNATURES[c] for c in sig.columns]
    return sig


def feature_class(feat: str) -> str:
    if feat in CIRCULAR_WITH_ESTIMATE:
        return "estimate_score_circular"
    if feat == "LeukocyteFraction_methylation":
        return "leukocyte_fraction_DNA"
    if feat.startswith("CIBERSORT"):
        return "CIBERSORT"
    if feat.startswith("Wolf"):
        return "Wolf_signature"
    if feat in ("CYT_GZMA_PRF1", "CD8_score", "GEP18"):
        return "gene_score"
    return "gene_marker"


def correlate_block(df: pd.DataFrame, predictor: str, purity_col: str, features: list[str],
                    matrix: str) -> pd.DataFrame:
    rows = []
    x_all = df[predictor].to_numpy(dtype=float)
    z_all = df[purity_col].to_numpy(dtype=float)
    for feat in features:
        y_all = df[feat].to_numpy(dtype=float)
        mask = np.isfinite(x_all) & np.isfinite(y_all) & np.isfinite(z_all)
        n = int(mask.sum())
        if n < 10:
            rho = p = prho = pp = np.nan
        else:
            rho, p = stats.spearmanr(x_all[mask], y_all[mask])
            prho, pp = partial_spearman(x_all[mask], y_all[mask], z_all[mask])
        rows.append(dict(
            matrix=matrix,
            predictor=predictor,
            purity_method=purity_col,
            feature=feat,
            feature_class=feature_class(feat),
            n=n,
            spearman_rho=rho,
            spearman_p=p,
            partial_rho=prho,
            partial_p=pp,
            circular_with_this_purity=bool(
                feat in CIRCULAR_WITH_ESTIMATE and purity_col == "estimate_purity"
            ),
        ))
    res = pd.DataFrame(rows)
    res["spearman_fdr"] = bh_fdr(res["spearman_p"])
    res["partial_fdr"] = bh_fdr(res["partial_p"])
    res["sig_partial_fdr05"] = res["partial_fdr"] < 0.05
    return res


def purity_context(df: pd.DataFrame, matrix: str) -> pd.DataFrame:
    rows = []
    pairs = [
        ("TACSTD2", "CLDN4"),
        ("TACSTD2", "estimate_purity"),
        ("TACSTD2", "absolute_purity"),
        ("CLDN4", "estimate_purity"),
        ("CLDN4", "absolute_purity"),
        ("estimate_purity", "absolute_purity"),
        ("estimate_purity", "aran_estimate_purity"),
        ("absolute_purity", "aran_absolute_purity"),
        ("estimate_purity", "LeukocyteFraction_methylation"),
        ("absolute_purity", "LeukocyteFraction_methylation"),
        ("ESTIMATE_Immune_score", "estimate_purity"),
        ("ESTIMATE_Immune_score", "absolute_purity"),
        ("CYT_GZMA_PRF1", "estimate_purity"),
        ("CYT_GZMA_PRF1", "absolute_purity"),
        ("CD8_score", "absolute_purity"),
        ("GEP18", "absolute_purity"),
    ]
    for a, b in pairs:
        if a not in df.columns or b not in df.columns:
            continue
        r, p, n = spearman_row(df[a].to_numpy(float), df[b].to_numpy(float))
        rows.append(dict(matrix=matrix, a=a, b=b, n=n, spearman_rho=r, spearman_p=p))
    return pd.DataFrame(rows)


def make_figures(df: pd.DataFrame, res: pd.DataFrame) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    prim = res[(res["matrix"] == "HiSeqV2") & (res["purity_method"].isin(
        ["estimate_purity", "absolute_purity"]
    ))].copy()

    def bar_panel(predictor: str, path: str, title: str) -> None:
        sub = prim[prim["predictor"] == predictor]
        feats = (
            sub[sub["purity_method"] == "absolute_purity"]
            .sort_values("partial_rho")["feature"].tolist()
        )
        fig, ax = plt.subplots(figsize=(9, 0.28 * len(feats) + 1.8))
        y = np.arange(len(feats))
        h = 0.38
        for shift, method, color, label in (
            (-h / 2, "estimate_purity", "#1f77b4", "partial | ESTIMATE purity"),
            (+h / 2, "absolute_purity", "#d62728", "partial | ABSOLUTE purity"),
        ):
            m = sub[sub["purity_method"] == method].set_index("feature").loc[feats]
            ax.barh(y + shift, m["partial_rho"], height=h, color=color, alpha=0.85, label=label)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(feats, fontsize=7)
        ax.set_xlabel("Partial Spearman rho")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc="lower right")
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)

    bar_panel(
        "TACSTD2",
        os.path.join(FIG, "bar_tacstd2_partial_rho.png"),
        "TCGA-LUAD HiSeqV2: TACSTD2 vs immune features\n"
        "partial Spearman after ESTIMATE vs ABSOLUTE purity",
    )
    bar_panel(
        "CLDN4",
        os.path.join(FIG, "bar_cldn4_partial_rho.png"),
        "TCGA-LUAD HiSeqV2: CLDN4 vs immune features\n"
        "partial Spearman after ESTIMATE vs ABSOLUTE purity",
    )

    plot_feats = [
        "LeukocyteFraction_methylation", "ESTIMATE_Immune_score",
        "CYT_GZMA_PRF1", "CD8_score", "GEP18", "CD274",
    ]
    for pred, fname in (("TACSTD2", "scatter_tacstd2_vs_immune.png"),
                        ("CLDN4", "scatter_cldn4_vs_immune.png")):
        fig, axes = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)
        sc = None
        for ax, feat in zip(axes.ravel(), plot_feats):
            mask = df[feat].notna() & df["absolute_purity"].notna()
            sc = ax.scatter(
                df.loc[mask, pred], df.loc[mask, feat],
                c=df.loc[mask, "absolute_purity"], cmap="viridis", s=10, alpha=0.7,
            )
            row = prim[(prim["predictor"] == pred) & (prim["feature"] == feat)
                       & (prim["purity_method"] == "absolute_purity")]
            if len(row):
                r = row.iloc[0]
                ax.set_title(
                    f"{feat}\nrho={r.spearman_rho:.2f}; "
                    f"|ABS rho={r.partial_rho:.2f} (FDR={r.partial_fdr:.1e}, n={int(r.n)})",
                    fontsize=8,
                )
            ax.set_xlabel(f"{pred} log2(RSEM+1)")
            ax.set_ylabel(feat, fontsize=8)
        if sc is not None:
            fig.colorbar(sc, ax=axes, label="ABSOLUTE purity", shrink=0.55)
        fig.suptitle(f"TCGA-LUAD primary tumors: {pred} vs immune features", fontsize=12)
        fig.savefig(os.path.join(FIG, fname), dpi=150)
        plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), constrained_layout=True)
    m = df["absolute_purity"].notna() & df["estimate_purity"].notna()
    axes[0].scatter(df.loc[m, "estimate_purity"], df.loc[m, "absolute_purity"], s=10, alpha=0.6, c="#333")
    r, p, n = spearman_row(df.loc[m, "estimate_purity"].to_numpy(float),
                           df.loc[m, "absolute_purity"].to_numpy(float))
    axes[0].set_xlabel("ESTIMATE purity (Yoshihara formula)")
    axes[0].set_ylabel("ABSOLUTE purity")
    axes[0].set_title(f"ESTIMATE vs ABSOLUTE\nrho={r:.2f} p={p:.1e} n={n}")
    for ax, pred in zip(axes[1:], PREDICTORS):
        sc = ax.scatter(df.loc[m, pred], df.loc[m, "absolute_purity"],
                        c=df.loc[m, "estimate_purity"], cmap="coolwarm", s=10, alpha=0.7)
        r1, p1, _ = spearman_row(df.loc[m, pred].to_numpy(float),
                                 df.loc[m, "absolute_purity"].to_numpy(float))
        r2, p2, _ = spearman_row(df.loc[m, pred].to_numpy(float),
                                 df.loc[m, "estimate_purity"].to_numpy(float))
        ax.set_xlabel(f"{pred} log2(RSEM+1)")
        ax.set_ylabel("ABSOLUTE purity")
        ax.set_title(f"{pred} vs purity\n|ABS rho={r1:.2f}; |EST rho={r2:.2f}")
    fig.colorbar(sc, ax=axes[1:], label="ESTIMATE purity", shrink=0.8)
    fig.savefig(os.path.join(FIG, "purity_context.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    m = df["absolute_purity"].notna()
    sc = ax.scatter(df.loc[m, "TACSTD2"], df.loc[m, "CLDN4"],
                    c=df.loc[m, "absolute_purity"], cmap="viridis", s=12, alpha=0.75)
    r, p, n = spearman_row(df.loc[m, "TACSTD2"].to_numpy(float),
                           df.loc[m, "CLDN4"].to_numpy(float))
    ax.set_xlabel("TACSTD2 log2(RSEM+1)")
    ax.set_ylabel("CLDN4 log2(RSEM+1)")
    ax.set_title(f"TACSTD2 vs CLDN4 in TCGA-LUAD\nSpearman rho={r:.2f} p={p:.1e} n={n}")
    fig.colorbar(sc, ax=ax, label="ABSOLUTE purity")
    fig.savefig(os.path.join(FIG, "scatter_tacstd2_vs_cldn4.png"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    genes = sorted(set(PREDICTORS + MARKER_GENES + GEP18))
    print("Loading HiSeqV2 ...")
    hiseq = add_derived(load_hiseqv2(genes))
    print(f"  primary tumors: {hiseq.shape[0]}")

    print("Loading STAR TPM (sensitivity) ...")
    star = add_derived(load_star_tpm(genes))
    print(f"  primary tumors: {star.shape[0]}")

    print("Loading purity / immune tables ...")
    pur = load_purity()
    leuk = load_leuk()
    cs = load_cibersort()

    def assemble(eg: pd.DataFrame) -> pd.DataFrame:
        df = eg.join(pur, how="left")
        df = df.join(leuk, how="left")
        df = df.join(cs, how="left")
        wolf = load_wolf(df.index)
        df = df.join(wolf, how="left")
        return df

    df_h = assemble(hiseq)
    df_s = assemble(star)

    # Primary analysis set: HiSeqV2 + both purity estimates
    primary = df_h.dropna(subset=["TACSTD2", "CLDN4", "estimate_purity", "absolute_purity"]).copy()
    print(f"Primary n (HiSeqV2 + ESTIMATE + ABSOLUTE): {len(primary)}")

    features = (
        ["LeukocyteFraction_methylation", "ESTIMATE_Immune_score", "ESTIMATE_Stromal_score"]
        + [c for c in WOLF_SIGNATURES.values()]
        + list(CIBERSORT_CELLS.values())
        + ["CYT_GZMA_PRF1", "CD8_score", "GEP18"]
        + MARKER_GENES
    )

    blocks = []
    for matrix, df in (("HiSeqV2", df_h), ("STAR_TPM", df_s)):
        for pred in PREDICTORS:
            for pur_col in ("estimate_purity", "absolute_purity"):
                work = df.dropna(subset=[pred, pur_col])
                blocks.append(correlate_block(work, pred, pur_col, features, matrix))
    res = pd.concat(blocks, ignore_index=True)

    ctx = pd.concat(
        [purity_context(df_h, "HiSeqV2"), purity_context(df_s, "STAR_TPM")],
        ignore_index=True,
    )

    # Two-covariate partial on the primary set (both purities at once)
    both_rows = []
    for pred in PREDICTORS:
        x = primary[pred].to_numpy(float)
        z1 = primary["estimate_purity"].to_numpy(float)
        z2 = primary["absolute_purity"].to_numpy(float)
        for feat in features:
            y = primary[feat].to_numpy(float)
            mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z1) & np.isfinite(z2)
            n = int(mask.sum())
            if n < 10:
                prho = pp = np.nan
            else:
                prho, pp = partial_spearman(x[mask], y[mask], z1[mask], z2[mask])
            both_rows.append(dict(
                matrix="HiSeqV2", predictor=pred, purity_method="estimate_and_absolute",
                feature=feat, feature_class=feature_class(feat), n=n,
                spearman_rho=np.nan, spearman_p=np.nan,
                partial_rho=prho, partial_p=pp,
                circular_with_this_purity=feat in CIRCULAR_WITH_ESTIMATE,
            ))
    both = pd.DataFrame(both_rows)
    both["spearman_fdr"] = np.nan
    both["partial_fdr"] = bh_fdr(both["partial_p"])
    both["sig_partial_fdr05"] = both["partial_fdr"] < 0.05
    res = pd.concat([res, both], ignore_index=True)

    res.to_csv(os.path.join(TABLES, "correlations.tsv"), sep="\t", index=False)
    ctx.to_csv(os.path.join(TABLES, "purity_context.tsv"), sep="\t", index=False)
    primary.reset_index().to_csv(os.path.join(TABLES, "analysis_table_hiseqv2.tsv"), sep="\t", index=False)
    star_prim = df_s.dropna(subset=["TACSTD2", "CLDN4", "estimate_purity", "absolute_purity"])
    star_prim.reset_index().to_csv(os.path.join(TABLES, "analysis_table_star_tpm.tsv"), sep="\t", index=False)

    make_figures(primary, res)

    # Provenance
    inputs = {
        "expression_hiseqv2": {
            "file": "TCGA.LUAD.HiSeqV2.gz",
            "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap/HiSeqV2.gz",
            "desc": "UCSC Xena TCGA-LUAD HiSeqV2, log2(RSEM normalized_count + 1)",
        },
        "expression_star_tpm": {
            "file": "TCGA-LUAD.star_tpm.tsv.gz",
            "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_tpm.tsv.gz",
            "desc": "UCSC Xena GDC hub TCGA-LUAD STAR log2(TPM+1)",
        },
        "probemap": {
            "file": "gencode.v36.annotation.gtf.gene.probemap",
            "url": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/gencode.v36.annotation.gtf.gene.probemap",
            "desc": "GENCODE v36 gene symbol map for STAR matrix",
        },
        "estimate_mdacc": {
            "file": "MDACC_estimate_LUAD_RNAseqV2.txt",
            "url": "https://bioinformatics.mdanderson.org/estimate/tables/lung_adenocarcinoma_RNAseqV2.txt",
            "desc": "MD Anderson precomputed ESTIMATE scores (RNAseqV2 LUAD)",
        },
        "absolute": {
            "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
            "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
            "desc": "PanCanAtlas ABSOLUTE purity/ploidy (DNA/SNP-array)",
        },
        "aran2015": {
            "file": "Aran2015_purity_S1.xlsx",
            "url": "https://static-content.springer.com/esm/art%3A10.1038%2Fncomms9971/MediaObjects/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx",
            "desc": "Aran 2015 Nat Commun Supp Data 1 (ESTIMATE/ABSOLUTE/LUMP/IHC/CPE)",
        },
        "leukocyte_fraction": {
            "file": "TCGA_all_leuk_estimate.masked.20170107.tsv",
            "url": "https://api.gdc.cancer.gov/data/6f75c9d7-5134-4ed1-b8f3-72856c98a4e8",
            "desc": "PanImmune DNA-methylation leukocyte fraction (RNA-independent)",
        },
        "cibersort": {
            "file": "TCGA.Kallisto.fullIDs.cibersort.relative.tsv",
            "url": "https://api.gdc.cancer.gov/data/b3df502e-3594-46ef-9f94-d041a20a0b9a",
            "desc": "PanImmune CIBERSORT relative fractions",
        },
        "signatures": {
            "file": "Scores_160_Signatures.tsv.gz",
            "url": "https://api.gdc.cancer.gov/data/80a82092-161d-4615-9d96-e858f113618d",
            "desc": "PanImmune 160 signature scores (Wolf subset used)",
        },
    }
    for v in inputs.values():
        p = os.path.join(RAW, v["file"])
        v["bytes"] = os.path.getsize(p)
        v["sha256"] = sha256(p)

    def ctx_val(matrix, a, b, field):
        hit = ctx[(ctx["matrix"] == matrix) & (ctx["a"] == a) & (ctx["b"] == b)]
        if hit.empty:
            return None
        return float(hit.iloc[0][field]) if pd.notna(hit.iloc[0][field]) else None

    prov = {
        "question": (
            "In UCSC Xena TCGA-LUAD primary tumors, are TACSTD2 and CLDN4 "
            "associated with immune features after ESTIMATE and/or ABSOLUTE purity adjustment?"
        ),
        "primary_matrix": "HiSeqV2 (matched to MDACC ESTIMATE RNAseqV2)",
        "sensitivity_matrix": "GDC STAR log2(TPM+1)",
        "n_hiseqv2_primary": int(hiseq.shape[0]),
        "n_star_primary": int(star.shape[0]),
        "n_primary_both_purities": int(len(primary)),
        "n_star_both_purities": int(len(star_prim)),
        "methods": {
            "unadjusted": "Spearman",
            "adjusted": (
                "partial Spearman (Pearson on rank residuals after regressing out "
                "ranked purity); t-test with n-2-k df"
            ),
            "estimate_purity": (
                f"cos({ESTIMATE_A} + {ESTIMATE_B} * ESTIMATE_score), Yoshihara 2013; "
                "scores from MD Anderson RNAseqV2 LUAD table"
            ),
            "absolute_purity": "PanCanAtlas ABSOLUTE consensus calls (DNA/SNP-array)",
            "multiple_testing": (
                "BH-FDR within each (matrix, predictor, purity_method) block"
            ),
            "duplicates": "aliquots averaged to patient-level -01 sample",
            "circularity": (
                "ESTIMATE Immune/Stromal scores are terms in ESTIMATE_score; "
                "partialling ESTIMATE purity out of those scores is circular. "
                "ESTIMATE purity itself is RNA-derived, so adjusting RNA immune "
                "genes for ESTIMATE purity is only partly independent. "
                "ABSOLUTE (DNA) and methylation leukocyte fraction are the "
                "non-circular readouts."
            ),
        },
        "key_context_spearman_hiseqv2": {
            "TACSTD2_vs_CLDN4": ctx_val("HiSeqV2", "TACSTD2", "CLDN4", "spearman_rho"),
            "TACSTD2_vs_ESTIMATE_purity": ctx_val("HiSeqV2", "TACSTD2", "estimate_purity", "spearman_rho"),
            "TACSTD2_vs_ABSOLUTE_purity": ctx_val("HiSeqV2", "TACSTD2", "absolute_purity", "spearman_rho"),
            "CLDN4_vs_ESTIMATE_purity": ctx_val("HiSeqV2", "CLDN4", "estimate_purity", "spearman_rho"),
            "CLDN4_vs_ABSOLUTE_purity": ctx_val("HiSeqV2", "CLDN4", "absolute_purity", "spearman_rho"),
            "ESTIMATE_vs_ABSOLUTE_purity": ctx_val("HiSeqV2", "estimate_purity", "absolute_purity", "spearman_rho"),
            "ESTIMATE_purity_vs_leukocyte_fraction": ctx_val(
                "HiSeqV2", "estimate_purity", "LeukocyteFraction_methylation", "spearman_rho"
            ),
            "ABSOLUTE_purity_vs_leukocyte_fraction": ctx_val(
                "HiSeqV2", "absolute_purity", "LeukocyteFraction_methylation", "spearman_rho"
            ),
        },
        "inputs": inputs,
    }
    with open(os.path.join(OUT, "provenance.json"), "w") as f:
        json.dump(prov, f, indent=2)

    # Console summary for the report
    print("\n===== PURITY CONTEXT (HiSeqV2) =====")
    print(ctx[ctx["matrix"] == "HiSeqV2"].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\n===== PRIMARY: TACSTD2 | ABSOLUTE (HiSeqV2) =====")
    show = res[(res["matrix"] == "HiSeqV2") & (res["predictor"] == "TACSTD2")
               & (res["purity_method"] == "absolute_purity")].sort_values("partial_p")
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\n===== PRIMARY: CLDN4 | ABSOLUTE (HiSeqV2) =====")
    show = res[(res["matrix"] == "HiSeqV2") & (res["predictor"] == "CLDN4")
               & (res["purity_method"] == "absolute_purity")].sort_values("partial_p")
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\n===== PRIMARY: TACSTD2 | ESTIMATE (HiSeqV2) =====")
    show = res[(res["matrix"] == "HiSeqV2") & (res["predictor"] == "TACSTD2")
               & (res["purity_method"] == "estimate_purity")].sort_values("partial_p")
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("\n===== PRIMARY: CLDN4 | ESTIMATE (HiSeqV2) =====")
    show = res[(res["matrix"] == "HiSeqV2") & (res["predictor"] == "CLDN4")
               & (res["purity_method"] == "estimate_purity")].sort_values("partial_p")
    print(show.to_string(index=False, float_format=lambda v: f"{v:.3g}"))
    print("done")


if __name__ == "__main__":
    main()
