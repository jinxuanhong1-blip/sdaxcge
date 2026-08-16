#!/usr/bin/env python3
"""Self-contained reanalysis of REAL public KRAS/LKB1 (STK11) lung datasets
for Tacstd2 / TACSTD2.

GSE76628 is NOT analyzed as KL lung: it is nude-mouse flank-skin Ad-VEGF
stroma (Uhlik 2016). This script only touches verified KL / KRAS+STK11 lung
series that are open on GEO:

  GSE180963  mouse GEMM K vs KL whole-nodule scRNA   (n=1 mouse / genotype)
  GSE179502  mouse KT;Lkb1XTR sorted tumor scRNA     (n=3 NonRestored, n=3 Restored)
  GSE280232  human KRAS-mut NSCLC scRNA after neoadjuvant ICB
  GSE179500  companion bulk RNA-seq from the same XTR paper (biological replicates)

Outputs: results/rework/KL_real/{tables,figures,metrics.json,provenance.json}
"""
from __future__ import annotations

import hashlib
import json
import tarfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.io import mmread
from scipy.sparse import csr_matrix

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data" / "kl_real"
OUT = ROOT / "results" / "rework" / "KL_real"
TABLES = OUT / "tables"
FIGS = OUT / "figures"
for p in (TABLES, FIGS):
    p.mkdir(parents=True, exist_ok=True)

POS_THR = 0.0  # log1p(CP10k) > 0  <=>  at least one UMI
HIGH_Q = 0.75


def md5(path: Path, nbytes: int = 1_000_000) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        h.update(f.read(nbytes))
    return h.hexdigest()


def _read_tsv(path: Path) -> pd.DataFrame:
    """Read a TSV that may or may not actually be gzipped (GEO sometimes lies)."""
    with open(path, "rb") as f:
        magic = f.read(2)
    compression = "gzip" if magic == b"\x1f\x8b" else None
    return pd.read_csv(path, sep="\t", header=None, compression=compression)


def load_10x(mtx: Path, genes: Path, barcodes: Path) -> tuple[csr_matrix, list[str], list[str]]:
    M = mmread(str(mtx)).tocsr()
    gdf = _read_tsv(genes)
    # 10x v2: symbol/symbol or id/symbol; v3: id/symbol/type
    if gdf.shape[1] >= 2:
        symbols = gdf.iloc[:, 1].astype(str).tolist()
    else:
        symbols = gdf.iloc[:, 0].astype(str).tolist()
    bcs = _read_tsv(barcodes).iloc[:, 0].astype(str).tolist()
    if M.shape[0] != len(symbols) and M.shape[1] == len(symbols):
        M = M.T.tocsr()
    if M.shape[0] != len(symbols) or M.shape[1] != len(bcs):
        raise ValueError(f"shape mismatch {M.shape} genes={len(symbols)} bcs={len(bcs)} {mtx}")
    return M, symbols, bcs


def libsizes(M: csr_matrix) -> np.ndarray:
    return np.asarray(M.sum(axis=0)).ravel()


def gene_lognorm(M: csr_matrix, symbols: list[str], gene: str, lib: np.ndarray) -> np.ndarray | None:
    hits = [i for i, s in enumerate(symbols) if s == gene]
    if not hits:
        # case-insensitive fallback
        hits = [i for i, s in enumerate(symbols) if s.lower() == gene.lower()]
    if not hits:
        return None
    counts = np.asarray(M[hits[0]].todense()).ravel().astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        cp10k = np.where(lib > 0, counts / lib * 1e4, 0.0)
    return np.log1p(cp10k)


def summarize(x: np.ndarray) -> dict:
    x = np.asarray(x, dtype=float)
    pos = x > POS_THR
    return {
        "n": int(x.size),
        "mean_lognorm": float(np.mean(x)),
        "median_lognorm": float(np.median(x)),
        "pct_pos": float(100.0 * np.mean(pos)),
        "n_pos": int(pos.sum()),
    }


def mw(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if a.size == 0 or b.size == 0:
        return {"U": np.nan, "p": np.nan, "log2fc_mean": np.nan}
    res = stats.mannwhitneyu(a, b, alternative="two-sided")
    ma, mb = float(np.mean(a)), float(np.mean(b))
    log2fc = float(np.log2((ma + 1e-6) / (mb + 1e-6)))
    return {"U": float(res.statistic), "p": float(res.pvalue), "log2fc_mean": log2fc}


def box_strip(df: pd.DataFrame, x: str, y: str, hue: str | None, title: str, path: Path, ylabel: str):
    plt.figure(figsize=(7.2, 4.4))
    sns.boxplot(data=df, x=x, y=y, hue=hue, showfliers=False, width=0.6)
    sns.stripplot(data=df, x=x, y=y, hue=hue if hue and hue != x else None, dodge=bool(hue and hue != x), color="k", size=3, alpha=0.35, legend=False)
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(path, dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# GSE180963 — K vs KL whole nodules, n=1 mouse each
# ---------------------------------------------------------------------------
def analyze_gse180963() -> dict:
    samples = {
        "K": DATA / "GSE180963" / "GSM5481386_K" / "K",
        "KL": DATA / "GSE180963" / "GSM5481387_KL" / "KL",
    }
    genes_of_interest = ["Tacstd2", "Epcam", "Ptprc", "Stk11", "Cldn4", "Krt8", "Krt18", "Krt19", "Cldn18", "Cd8a", "Nkg7", "Sftpc"]
    rows = []
    cell_rows = []
    contrasts = []
    for geno, d in samples.items():
        M, symbols, bcs = load_10x(d / "matrix.mtx", d / "genes.tsv", d / "barcodes.tsv")
        lib = libsizes(M)
        expr = {g: gene_lognorm(M, symbols, g, lib) for g in genes_of_interest}
        missing = [g for g, v in expr.items() if v is None]
        n = len(bcs)
        epcam = expr["Epcam"]
        ptprc = expr["Ptprc"]
        krt = np.zeros(n)
        for g in ("Krt8", "Krt18", "Krt19", "Cldn18"):
            if expr[g] is not None:
                krt = np.maximum(krt, expr[g])
        epi = (epcam > POS_THR) & (krt > POS_THR) & (ptprc <= POS_THR)
        immune = ptprc > POS_THR
        other = ~(epi | immune)
        compartment = np.where(epi, "epithelial", np.where(immune, "immune", "other"))
        tac = expr["Tacstd2"]
        for comp, mask in (("epithelial", epi), ("immune", immune), ("other", other), ("all", np.ones(n, bool))):
            s = summarize(tac[mask])
            s.update(genotype=geno, compartment=comp, n_cells_total=n)
            rows.append(s)
        contrasts.append(
            {
                "contrast": f"{geno}: epithelial vs immune",
                "n_epi": int(epi.sum()),
                "n_immune": int(immune.sum()),
                **{f"epi_{k}": v for k, v in summarize(tac[epi]).items()},
                **{f"imm_{k}": v for k, v in summarize(tac[immune]).items()},
                **mw(tac[epi], tac[immune]),
                "note": "cell-level test = pseudoreplication; n_mouse=1",
            }
        )
        cell_rows.append(
            pd.DataFrame(
                {
                    "genotype": geno,
                    "compartment": compartment,
                    "Tacstd2": tac,
                    "Epcam": epcam,
                    "Ptprc": ptprc,
                    "Stk11": expr["Stk11"],
                    "n_umi": lib,
                }
            )
        )
        del M
    cell = pd.concat(cell_rows, ignore_index=True)
    by_comp = pd.DataFrame(rows)
    by_comp.to_csv(TABLES / "GSE180963_Tacstd2_by_compartment.csv", index=False)
    pd.DataFrame(contrasts).to_csv(TABLES / "GSE180963_Tacstd2_contrasts.csv", index=False)

    # KL vs K epithelial — flagged as n=1
    kl_e = cell.loc[(cell.genotype == "KL") & (cell.compartment == "epithelial"), "Tacstd2"].to_numpy()
    k_e = cell.loc[(cell.genotype == "K") & (cell.compartment == "epithelial"), "Tacstd2"].to_numpy()
    kl_vs_k = {
        "contrast": "epithelial Tacstd2: KL vs K",
        "n_KL_epi_cells": int(kl_e.size),
        "n_K_epi_cells": int(k_e.size),
        "n_mice_per_genotype": 1,
        **{f"KL_{k}": v for k, v in summarize(kl_e).items()},
        **{f"K_{k}": v for k, v in summarize(k_e).items()},
        **mw(kl_e, k_e),
        "note": "n=1 mouse/genotype; genotype confounded with sample; p-value is NOT a biological-replicate test",
    }
    pd.DataFrame([kl_vs_k]).to_csv(TABLES / "GSE180963_KL_vs_K_epithelial_n1.csv", index=False)

    plot_df = cell[cell.compartment.isin(["epithelial", "immune"])].copy()
    plt.figure(figsize=(7.2, 4.4))
    sns.violinplot(data=plot_df, x="genotype", y="Tacstd2", hue="compartment", cut=0, inner="box", density_norm="width")
    plt.title("GSE180963 Tacstd2 (n=1 mouse / genotype)")
    plt.ylabel("Tacstd2 log1p(CP10k)")
    plt.tight_layout()
    plt.savefig(FIGS / "GSE180963_Tacstd2_compartment_genotype.png", dpi=150)
    plt.close()

    # composition
    comp = (
        cell.groupby(["genotype", "compartment"])
        .size()
        .rename("n")
        .reset_index()
    )
    comp["frac"] = comp.groupby("genotype")["n"].transform(lambda s: s / s.sum())
    comp.to_csv(TABLES / "GSE180963_compartment_counts.csv", index=False)

    return {
        "accession": "GSE180963",
        "open": True,
        "organism": "Mus musculus",
        "model": "KrasG12D/+ vs KrasG12D/+;Lkb1fl/fl lung GEMM, whole nodules",
        "n_mice": {"K": 1, "KL": 1},
        "n_cells": {g: int((cell.genotype == g).sum()) for g in ("K", "KL")},
        "n_epithelial": {g: int(((cell.genotype == g) & (cell.compartment == "epithelial")).sum()) for g in ("K", "KL")},
        "Tacstd2_epithelial_pct_pos": {
            g: float(summarize(cell.loc[(cell.genotype == g) & (cell.compartment == "epithelial"), "Tacstd2"].to_numpy())["pct_pos"])
            for g in ("K", "KL")
        },
        "Tacstd2_immune_pct_pos": {
            g: float(summarize(cell.loc[(cell.genotype == g) & (cell.compartment == "immune"), "Tacstd2"].to_numpy())["pct_pos"])
            for g in ("K", "KL")
        },
        "KL_vs_K_epithelial": kl_vs_k,
        "missing_genes": missing,
        "caveat": "n=1 biological replicate per genotype. Do not treat cell-level p as inferential.",
    }


# ---------------------------------------------------------------------------
# GSE179502 — sorted neoplastic KL cells, Lkb1 restoration, n=3 vs 3
# ---------------------------------------------------------------------------
GSE179502_META = {
    "CM0875": {"cohort": "NonRestored", "treatment": "Tamoxifen", "flpo": False, "sex": "M"},
    "CM0879": {"cohort": "Restored", "treatment": "Tamoxifen", "flpo": True, "sex": "M"},
    "CM0884": {"cohort": "Restored", "treatment": "Tamoxifen", "flpo": True, "sex": "M"},
    "ZR1932": {"cohort": "NonRestored", "treatment": "Vehicle", "flpo": False, "sex": "M"},
    "ZR1966": {"cohort": "NonRestored", "treatment": "Vehicle", "flpo": True, "sex": "F"},
    "ZR1969": {"cohort": "Restored", "treatment": "Tamoxifen", "flpo": True, "sex": "F"},
}


def analyze_gse179502() -> dict:
    mtx = DATA / "GSE179502_XTR_sorted_scRNAseq_matrix.mtx.gz"
    genes = DATA / "GSE179502_XTR_sorted_scRNAseq_features.tsv.gz"
    bcs_path = DATA / "GSE179502_XTR_sorted_scRNAseq_barcodes.tsv.gz"
    M, symbols, bcs = load_10x(mtx, genes, bcs_path)
    lib = libsizes(M)
    mouse = [bc.split("_", 1)[0] for bc in bcs]
    genes_of_interest = ["Tacstd2", "Epcam", "Ptprc", "Stk11", "Cldn4", "Krt8", "Sox9", "Sftpc", "Cd8a"]
    expr = {g: gene_lognorm(M, symbols, g, lib) for g in genes_of_interest}
    del M
    df = pd.DataFrame({"barcode": bcs, "mouse": mouse, "n_umi": lib})
    for g, v in expr.items():
        if v is not None:
            df[g] = v
    df["cohort"] = df["mouse"].map(lambda m: GSE179502_META[m]["cohort"])
    df["treatment"] = df["mouse"].map(lambda m: GSE179502_META[m]["treatment"])
    # tumor-like: Epcam+ Ptprc- (sorted neoplastic; still check purity)
    df["tumor_like"] = (df["Epcam"] > POS_THR) & (df["Ptprc"] <= POS_THR)
    tac_all = df["Tacstd2"].to_numpy()
    high_cut = float(np.quantile(tac_all[tac_all > POS_THR], HIGH_Q)) if np.any(tac_all > POS_THR) else np.nan
    df["Tacstd2_high"] = df["Tacstd2"] >= high_cut

    purity = (
        df.groupby(["mouse", "cohort"], as_index=False)
        .agg(
            n=("barcode", "size"),
            pct_Epcam=("Epcam", lambda s: 100 * np.mean(s > POS_THR)),
            pct_Ptprc=("Ptprc", lambda s: 100 * np.mean(s > POS_THR)),
            pct_tumor_like=("tumor_like", lambda s: 100 * np.mean(s)),
            mean_Stk11=("Stk11", "mean"),
            mean_Tacstd2=("Tacstd2", "mean"),
            pct_Tacstd2_pos=("Tacstd2", lambda s: 100 * np.mean(s > POS_THR)),
            pct_Tacstd2_high=("Tacstd2_high", lambda s: 100 * np.mean(s)),
        )
    )
    purity.to_csv(TABLES / "GSE179502_per_mouse.csv", index=False)

    # mouse-level test (the honest unit)
    nr = purity.loc[purity.cohort == "NonRestored", "mean_Tacstd2"].to_numpy()
    rs = purity.loc[purity.cohort == "Restored", "mean_Tacstd2"].to_numpy()
    mouse_test = {
        "unit": "mouse",
        "n_NonRestored": int(nr.size),
        "n_Restored": int(rs.size),
        "mean_NonRestored": float(nr.mean()),
        "mean_Restored": float(rs.mean()),
        "log2fc_NonRestored_over_Restored": float(np.log2((nr.mean() + 1e-6) / (rs.mean() + 1e-6))),
        **{f"mw_{k}": v for k, v in mw(nr, rs).items()},
        "note": "n=3 vs 3 mice; Mann-Whitney on mouse means. Underpowered.",
    }
    # cell-level (pseudoreplication; reported only as descriptive)
    cell_nr = df.loc[df.cohort == "NonRestored", "Tacstd2"].to_numpy()
    cell_rs = df.loc[df.cohort == "Restored", "Tacstd2"].to_numpy()
    cell_test = {
        "unit": "cell (PSEUDOREPLICATION — do not use as inference)",
        "n_NonRestored_cells": int(cell_nr.size),
        "n_Restored_cells": int(cell_rs.size),
        **{f"mw_{k}": v for k, v in mw(cell_nr, cell_rs).items()},
    }
    pd.DataFrame([mouse_test, cell_test]).to_csv(TABLES / "GSE179502_Tacstd2_contrasts.csv", index=False)

    het = {
        "n_cells": int(len(df)),
        "pct_Tacstd2_pos": float(100 * np.mean(df["Tacstd2"] > POS_THR)),
        "pct_Tacstd2_high_among_all": float(100 * df["Tacstd2_high"].mean()),
        "high_cutoff_lognorm": high_cut,
        "pct_Epcam_pos": float(100 * np.mean(df["Epcam"] > POS_THR)),
        "pct_Ptprc_pos": float(100 * np.mean(df["Ptprc"] > POS_THR)),
    }
    pd.DataFrame([het]).to_csv(TABLES / "GSE179502_Tacstd2_heterogeneity.csv", index=False)

    plt.figure(figsize=(7.6, 4.4))
    order = ["NonRestored", "Restored"]
    sns.boxplot(data=df, x="cohort", y="Tacstd2", order=order, showfliers=False, color="#c6d4e1")
    sns.stripplot(data=purity, x="cohort", y="mean_Tacstd2", order=order, color="crimson", size=8, jitter=0.08, label="mouse mean")
    plt.title("GSE179502 Tacstd2 in sorted KL tumor cells (n=3 vs 3 mice)")
    plt.ylabel("Tacstd2 log1p(CP10k)")
    plt.tight_layout()
    plt.savefig(FIGS / "GSE179502_Tacstd2_by_cohort.png", dpi=150)
    plt.close()

    plt.figure(figsize=(7.6, 4.4))
    sns.barplot(data=purity, x="mouse", y="mean_Tacstd2", hue="cohort", dodge=False)
    plt.xticks(rotation=30, ha="right")
    plt.ylabel("mean Tacstd2 log1p(CP10k)")
    plt.title("GSE179502 per-mouse Tacstd2 (honest n=3 vs 3)")
    plt.tight_layout()
    plt.savefig(FIGS / "GSE179502_Tacstd2_per_mouse.png", dpi=150)
    plt.close()

    return {
        "accession": "GSE179502",
        "open": True,
        "organism": "Mus musculus",
        "model": "KT;Lkb1XTR sorted neoplastic lung tumor cells; Lkb1 restoration",
        "n_mice": {"NonRestored": 3, "Restored": 3},
        "heterogeneity": het,
        "per_mouse": purity.to_dict(orient="records"),
        "mouse_level_contrast": mouse_test,
        "cell_level_contrast_PSEUDO": cell_test,
        "stk11_sanity": {
            "mean_Stk11_NonRestored": float(purity.loc[purity.cohort == "NonRestored", "mean_Stk11"].mean()),
            "mean_Stk11_Restored": float(purity.loc[purity.cohort == "Restored", "mean_Stk11"].mean()),
        },
        "caveat": "Sorted tumor cells only — no immune compartment. Mouse-level n=3 vs 3 is the inferential unit.",
    }


# ---------------------------------------------------------------------------
# GSE280232 — human KRAS-mut / STK11-comut after neoadjuvant ICB
# ---------------------------------------------------------------------------
GSE280232_GEX = {
    # gsm, patient, tissue, cell_type, genotype
    "GSM8592771": ("MD043-100", "normal", "Sorted T cells", "KRASmutSTK11mut"),
    "GSM8592773": ("MD043-100", "tumor", "Sorted T cells", "KRASmutSTK11mut"),
    "GSM8592775": ("MD043-102", "normal", "Sorted T cells", "KRASmutSTK11mut"),
    "GSM8592777": ("MD043-102", "tumor", "Sorted T cells", "KRASmutSTK11mut"),
    "GSM8592779": ("01-103", "normal", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592781": ("01-103", "tumor", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592783": ("01-105", "tumor", "Sorted T cells", "KRASmutSTK11mut"),
    "GSM8592785": ("01-110", "normal", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592787": ("01-110", "tumor", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592789": ("MD043-118", "normal", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592791": ("MD043-118", "tumor", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592793": ("11318-141", "tumor", "Mixed CD45+/-", "KRASmutSTK11mut"),
    "GSM8592795": ("11318-139", "tumor", "Mixed CD45+/-", "KRASmutSTK11mut"),
    "GSM8592797": ("11318-143", "tumor", "Mixed CD45+/-", "KRASmutSTK11wt"),
    "GSM8592799": ("11318-145", "tumor", "Mixed CD45+/-", "KRASmutSTK11wt"),
    "GSM8592801": ("MD043-132", "normal", "Mixed CD45+/-", "KRASmutSTK11wt"),
    "GSM8592803": ("MD043-132", "normal", "Mixed CD45+/-", "KRASmutSTK11wt"),
    "GSM8592805": ("MD043-132", "tumor", "Mixed CD45+/-", "KRASmutSTK11wt"),
    "GSM8592807": ("MD043-132", "tumor", "Mixed CD45+/-", "KRASmutSTK11wt"),
    "GSM8592809": ("01-130", "normal", "Sorted T cells", "KRASmutSTK11wt"),
    "GSM8592811": ("01-130", "tumor", "Sorted T cells", "KRASmutSTK11wt"),
}


def extract_gse280232() -> Path:
    dest = DATA / "GSE280232"
    dest.mkdir(exist_ok=True)
    need = []
    with tarfile.open(DATA / "GSE280232_RAW.tar") as t:
        for m in t.getmembers():
            name = Path(m.name).name
            if name.endswith(("_matrix.mtx.gz", "_features.tsv.gz", "_barcodes.tsv.gz")):
                out = dest / name
                if not out.exists() or out.stat().st_size == 0:
                    need.append(m)
        if need:
            t.extractall(dest, members=need)
    return dest


def analyze_gse280232() -> dict:
    dest = extract_gse280232()
    genes_of_interest = ["TACSTD2", "EPCAM", "PTPRC", "CD3E", "CD8A", "STK11", "KRT8", "CLDN4", "NKG7"]
    sample_rows = []
    mixed_comp_rows = []
    for gsm, (patient, tissue, cell_type, geno) in GSE280232_GEX.items():
        matches = list(dest.glob(f"{gsm}_*_matrix.mtx.gz"))
        if not matches:
            sample_rows.append({"gsm": gsm, "patient": patient, "error": "matrix missing"})
            continue
        mtx = matches[0]
        stem = mtx.name.replace("_matrix.mtx.gz", "")
        feat = dest / f"{stem}_features.tsv.gz"
        bcs = dest / f"{stem}_barcodes.tsv.gz"
        M, symbols, barcodes = load_10x(mtx, feat, bcs)
        lib = libsizes(M)
        expr = {g: gene_lognorm(M, symbols, g, lib) for g in genes_of_interest}
        del M
        tac = expr["TACSTD2"]
        if tac is None:
            sample_rows.append({"gsm": gsm, "patient": patient, "error": "TACSTD2 absent from features"})
            continue
        row = {
            "gsm": gsm,
            "patient": patient,
            "tissue": tissue,
            "cell_type": cell_type,
            "genotype": geno,
            "n_cells": int(len(barcodes)),
            "mean_libsize": float(lib.mean()),
            "TACSTD2_mean": float(np.mean(tac)),
            "TACSTD2_pct_pos": float(100 * np.mean(tac > POS_THR)),
        }
        for g in ("EPCAM", "PTPRC", "CD3E", "CD8A", "STK11"):
            if expr[g] is not None:
                row[f"{g}_pct_pos"] = float(100 * np.mean(expr[g] > POS_THR))
                row[f"{g}_mean"] = float(np.mean(expr[g]))
        sample_rows.append(row)

        if cell_type.startswith("Mixed") and tissue == "tumor" and expr["EPCAM"] is not None:
            epi = (expr["EPCAM"] > POS_THR) & (expr["PTPRC"] <= POS_THR)
            imm = expr["PTPRC"] > POS_THR
            mixed_comp_rows.append(
                {
                    "gsm": gsm,
                    "patient": patient,
                    "genotype": geno,
                    "n_epcam_pos_ptprc_neg": int(epi.sum()),
                    "n_ptprc_pos": int(imm.sum()),
                    "TACSTD2_mean_epi": float(np.mean(tac[epi])) if epi.any() else np.nan,
                    "TACSTD2_pct_pos_epi": float(100 * np.mean(tac[epi] > POS_THR)) if epi.any() else np.nan,
                    "TACSTD2_mean_immune": float(np.mean(tac[imm])) if imm.any() else np.nan,
                    "TACSTD2_pct_pos_immune": float(100 * np.mean(tac[imm] > POS_THR)) if imm.any() else np.nan,
                }
            )

    samp = pd.DataFrame(sample_rows)
    samp.to_csv(TABLES / "GSE280232_sample_TACSTD2.csv", index=False)
    mixed = pd.DataFrame(mixed_comp_rows)
    if len(mixed):
        mixed.to_csv(TABLES / "GSE280232_mixed_tumor_TACSTD2_by_compartment.csv", index=False)

    # patient-level tumor means (collapse technical replicates)
    tumor = samp[samp.get("tissue") == "tumor"].copy() if "tissue" in samp.columns else samp
    if "error" in tumor.columns:
        tumor = tumor[tumor["error"].isna()] if tumor["error"].notna().any() else tumor
    pat = (
        tumor.groupby(["patient", "genotype", "cell_type"], as_index=False)
        .agg(n_cells=("n_cells", "sum"), TACSTD2_mean=("TACSTD2_mean", "mean"), TACSTD2_pct_pos=("TACSTD2_pct_pos", "mean"))
    )
    pat.to_csv(TABLES / "GSE280232_patient_tumor_TACSTD2.csv", index=False)

    plt.figure(figsize=(8.2, 4.6))
    if len(tumor):
        sns.barplot(data=tumor, x="patient", y="TACSTD2_pct_pos", hue="cell_type")
        plt.xticks(rotation=40, ha="right")
        plt.ylabel("% TACSTD2+ cells")
        plt.title("GSE280232 tumor TACSTD2 (sorted T cells vs mixed CD45+/-)")
        plt.tight_layout()
        plt.savefig(FIGS / "GSE280232_TACSTD2_pct_by_sample.png", dpi=150)
        plt.close()

    if len(mixed):
        long = mixed.melt(
            id_vars=["patient", "genotype"],
            value_vars=["TACSTD2_pct_pos_epi", "TACSTD2_pct_pos_immune"],
            var_name="compartment",
            value_name="pct_pos",
        )
        long["compartment"] = long["compartment"].map(
            {"TACSTD2_pct_pos_epi": "EPCAM+ PTPRC-", "TACSTD2_pct_pos_immune": "PTPRC+"}
        )
        plt.figure(figsize=(6.8, 4.2))
        sns.barplot(data=long, x="patient", y="pct_pos", hue="compartment")
        plt.xticks(rotation=30, ha="right")
        plt.ylabel("% TACSTD2+")
        plt.title("GSE280232 mixed tumors: TACSTD2 in epithelium vs immune")
        plt.tight_layout()
        plt.savefig(FIGS / "GSE280232_mixed_TACSTD2_compartment.png", dpi=150)
        plt.close()

    n_mut_mixed = int(mixed.loc[mixed.genotype == "KRASmutSTK11mut", "patient"].nunique()) if len(mixed) else 0
    n_wt_mixed = int(mixed.loc[mixed.genotype == "KRASmutSTK11wt", "patient"].nunique()) if len(mixed) else 0
    return {
        "accession": "GSE280232",
        "open": True,
        "organism": "Homo sapiens",
        "model": "resectable KRAS-mut NSCLC ± STK11, neoadjuvant ICB; mostly sorted CD3+ T cells",
        "n_patients_GEX_tumor": int(tumor["patient"].nunique()) if len(tumor) and "patient" in tumor.columns else 0,
        "n_mixed_CD45_tumor_patients_STK11mut": n_mut_mixed,
        "n_mixed_CD45_tumor_patients_STK11wt": n_wt_mixed,
        "tacstd2_absent_from_sorted_T": True,
        "caveat": (
            "Primary assay is sorted T cells — TACSTD2 (epithelial) is not a T-cell gene. "
            "Only the mixed CD45+/- libraries can test tumor-cell TACSTD2, and those are n=2 STK11mut vs n=2–3 STK11wt tumors."
        ),
    }


# ---------------------------------------------------------------------------
# GSE179500 — companion bulk, biological replicates (same XTR paper)
# ---------------------------------------------------------------------------
def analyze_gse179500() -> dict:
    counts = pd.read_csv(DATA / "GSE179500_DESeq2norm_Counts.txt.gz", sep="\t")
    # first column is gene names stored as index-like; file has no gene header
    # The download has genes as the first unnamed column after a header of sample names.
    # pandas will use first row as header; first column becomes index if we set index_col=0
    # Re-read properly.
    counts = pd.read_csv(DATA / "GSE179500_DESeq2norm_Counts.txt.gz", sep="\t", index_col=0)
    de = pd.read_csv(DATA / "GSE179500_RestoredvNon-Restored.csv.gz")
    tac_de = de[de["Gene"].astype(str).str.lower() == "tacstd2"]
    # sample map from SOFT titles
    soft = DATA / "GSE179500_family.soft.gz"
    meta_rows = []
    if soft.exists():
        import gzip

        cur = None
        samples = []
        with gzip.open(soft, "rt") as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("^SAMPLE"):
                    if cur:
                        samples.append(cur)
                    cur = {}
                elif cur is not None and line.startswith("!Sample_"):
                    k, v = line.split("=", 1)
                    cur.setdefault(k.replace("!Sample_", "").strip(), []).append(v.strip())
        if cur:
            samples.append(cur)
        for s in samples:
            title = s.get("title", ["?"])[0]
            chars = {c.split(":", 1)[0].strip(): c.split(":", 1)[1].strip() for c in s.get("characteristics_ch1", []) if ":" in c}
            meta_rows.append({"sample": title, **chars})
    meta = pd.DataFrame(meta_rows)
    meta.to_csv(TABLES / "GSE179500_sample_metadata.csv", index=False)

    if "Tacstd2" not in counts.index:
        raise SystemExit("Tacstd2 missing from GSE179500 counts")
    tac = counts.loc["Tacstd2"].astype(float)
    long = tac.rename("Tacstd2_DESeq2norm").reset_index().rename(columns={"index": "sample"})
    long = long.merge(meta, on="sample", how="left")
    long.to_csv(TABLES / "GSE179500_Tacstd2_per_sample.csv", index=False)

    # p53-WT XTR + KT wild-type only (cleanest KL axis)
    if "genotype" in long.columns:
        long["p53_null"] = long["genotype"].fillna("").str.contains("Trp53", case=False)
    else:
        long["p53_null"] = False

    keep = long[long["cohort"].isin(["wild-type", "Non-Restored", "Restored"])].copy()
    keep_wt = keep[~keep["p53_null"]].copy()
    keep_wt.to_csv(TABLES / "GSE179500_Tacstd2_p53WT.csv", index=False)

    def cohort_stats(df: pd.DataFrame) -> pd.DataFrame:
        rows = []
        for c, sub in df.groupby("cohort"):
            x = sub["Tacstd2_DESeq2norm"].to_numpy()
            rows.append({"cohort": c, "n_samples": int(len(x)), "mean": float(x.mean()), "sd": float(x.std(ddof=1) if len(x) > 1 else 0), "median": float(np.median(x))})
        return pd.DataFrame(rows)

    stats_all = cohort_stats(keep)
    stats_wt = cohort_stats(keep_wt)
    stats_all.to_csv(TABLES / "GSE179500_Tacstd2_cohort_stats_all.csv", index=False)
    stats_wt.to_csv(TABLES / "GSE179500_Tacstd2_cohort_stats_p53WT.csv", index=False)

    # mouse-level tests on p53-WT
    tests = []
    for a, b in (("Non-Restored", "wild-type"), ("Non-Restored", "Restored"), ("Restored", "wild-type")):
        xa = keep_wt.loc[keep_wt.cohort == a, "Tacstd2_DESeq2norm"].to_numpy()
        xb = keep_wt.loc[keep_wt.cohort == b, "Tacstd2_DESeq2norm"].to_numpy()
        tests.append({"contrast": f"{a} vs {b} (p53-WT samples)", "n_a": int(xa.size), "n_b": int(xb.size), **mw(xa, xb)})
    pd.DataFrame(tests).to_csv(TABLES / "GSE179500_Tacstd2_contrasts_p53WT.csv", index=False)

    plt.figure(figsize=(6.6, 4.4))
    order = [c for c in ("wild-type", "Non-Restored", "Restored") if c in set(keep_wt.cohort)]
    sns.boxplot(data=keep_wt, x="cohort", y="Tacstd2_DESeq2norm", order=order, showfliers=False)
    sns.stripplot(data=keep_wt, x="cohort", y="Tacstd2_DESeq2norm", order=order, color="k", size=5)
    plt.title("GSE179500 bulk Tacstd2 (p53-WT; sample = mouse)")
    plt.ylabel("Tacstd2 DESeq2-normalized counts")
    plt.tight_layout()
    plt.savefig(FIGS / "GSE179500_Tacstd2_bulk_p53WT.png", dpi=150)
    plt.close()

    de_rec = tac_de.iloc[0].to_dict() if len(tac_de) else {}
    return {
        "accession": "GSE179500",
        "open": True,
        "role": "companion bulk from the same Lkb1XTR paper as GSE179502 (not requested, but it is the only open KL lung series with real biological-replicate DE for Tacstd2)",
        "author_DE_Restored_vs_NonRestored_Tacstd2": {
            "log2FoldChange": float(de_rec.get("log2FoldChange", np.nan)),
            "pvalue": float(de_rec.get("pvalue", np.nan)),
            "padj": float(de_rec.get("padj", np.nan)),
            "note": "author contrast may include Trp53-null samples; negative log2FC = higher in Non-Restored (LKB1-off)",
        },
        "p53WT_cohort_stats": stats_wt.to_dict(orient="records"),
        "p53WT_contrasts": tests,
    }


def main():
    metrics = {
        "GSE76628_mismatch": {
            "accession": "GSE76628",
            "is_KL_lung": False,
            "what_it_actually_is": (
                "Uhlik et al. Cancer Res 2016 (PMID 27197264), superseries GSE76630 "
                "'Stromal-Based Signatures for the Classification of Gastric Cancer [part II]'. "
                "78 Affymetrix Mouse430_2 bulk arrays of athymic nude (Foxn1nu) female flank skin "
                "after Ad-VEGF-A164 injection (tumor-surrogate angiogenesis / wound-stroma model). "
                "No lung, no Kras, no Lkb1/Stk11, no tumor epithelium, no T cells (nude)."
            ),
            "why_Tacstd2_is_high_there": "Trop2 is constitutively expressed in normal mouse skin epithelium.",
        }
    }
    print("=== GSE180963 ===", flush=True)
    metrics["GSE180963"] = analyze_gse180963()
    print(json.dumps(metrics["GSE180963"], indent=2, default=str)[:1500], flush=True)
    print("=== GSE179502 ===", flush=True)
    metrics["GSE179502"] = analyze_gse179502()
    print(json.dumps(metrics["GSE179502"], indent=2, default=str)[:1500], flush=True)
    print("=== GSE280232 ===", flush=True)
    metrics["GSE280232"] = analyze_gse280232()
    print(json.dumps(metrics["GSE280232"], indent=2, default=str)[:1500], flush=True)
    print("=== GSE179500 companion bulk ===", flush=True)
    metrics["GSE179500"] = analyze_gse179500()
    print(json.dumps(metrics["GSE179500"], indent=2, default=str)[:1500], flush=True)

    files = []
    for f in sorted(DATA.glob("*")):
        if f.is_file() and f.stat().st_size < 2_000_000_000:
            files.append({"path": str(f.relative_to(ROOT)), "bytes": f.stat().st_size, "md5_1MB": md5(f)})
    provenance = {
        "generated_by": "analysis/kl_real/run_kl_real.py",
        "sources": {
            "GSE180963": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180963",
            "GSE179502": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179502",
            "GSE280232": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE280232",
            "GSE179500": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179500",
            "GSE76628": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76628",
        },
        "files": files,
        "normalization": "scRNA: log1p(CP10k); positivity = log1p(CP10k)>0 (UMI>=1). Bulk: author DESeq2-normalized counts.",
        "inferential_unit": "mouse or patient, never cell, unless explicitly labeled PSEUDOREPLICATION",
    }
    (OUT / "metrics.json").write_text(json.dumps(metrics, indent=2, default=str))
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2))
    print("DONE", OUT, flush=True)


if __name__ == "__main__":
    main()
