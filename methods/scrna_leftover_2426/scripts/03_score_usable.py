#!/usr/bin/env python3
"""Score leftover public lung scRNA: epithelial/malignant TACSTD2/CLDN4 vs T/NK.

Public GEO processed matrices only. Patient (or sample) is the inferential unit.
Cell-level p-values are not reported as a patient-level claim.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats

DATA = Path("/tmp/scrna_leftover_2426")
OUT = Path(__file__).resolve().parents[1] / "tables"
OUT.mkdir(parents=True, exist_ok=True)

# Marker panels (symbol). Lineage = argmax of mean log1p(CP10k) scores.
MARKERS = {
    "epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CDH1"],
    "t": ["CD3D", "CD3E", "CD3G", "CD2", "TRAC", "CD8A", "CD4"],
    "nk": ["NKG7", "GNLY", "KLRD1", "KLRF1", "NCAM1", "KLRB1"],
    "b": ["MS4A1", "CD79A", "CD79B", "CD19"],
    "myeloid": ["LYZ", "CD68", "CD14", "CSF1R", "FCGR3A"],
    "fibroblast": ["COL1A1", "DCN", "LUM", "COL1A2"],
    "endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
}
NORMAL_LUNG = ["SFTPA2", "SFTPC", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]
TARGETS = ["TACSTD2", "CLDN4", "CD274"]

# GSE267108 PD-L1 groups inferred by matching published per-group cell counts
# (Zhang 2025 J Cancer; 38244 / 14231 / 6335) to GEO barcode counts (exact).
GSE267108_PDL1 = {
    "CHW": "neg",
    "HRZ": "neg",
    "LYS": "neg",
    "WXQ": "neg",
    "WWM": "pos",
    "ZXM": "pos",
    "LJM": "high",
    "CFR": "high",
}

GSE274595_META = {
    "GSM8453374": {"patient": "PT002", "histology": "LUAD", "assay": "tissue"},
    "GSM8453375": {"patient": "PT004", "histology": "LUAD", "assay": "tissue"},
    "GSM8453376": {"patient": "PT005", "histology": "LUSC", "assay": "tissue"},
    "GSM8453377": {"patient": "PT007", "histology": "NSCLC", "assay": "nuclei"},
    "GSM8453378": {"patient": "PT008", "histology": "NSCLC", "assay": "nuclei"},
    "GSM8453379": {"patient": "PT009", "histology": "NSCLC", "assay": "nuclei"},
    "GSM8453380": {"patient": "PT010", "histology": "NSCLC", "assay": "nuclei"},
    "GSM8453381": {"patient": "PT012", "histology": "NSCLC", "assay": "nuclei"},
}


def _read_features(path: Path) -> list[str]:
    genes = []
    import gzip

    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and not parts[1].startswith("ENSG"):
                genes.append(parts[1])
            else:
                genes.append(parts[0])
    return genes


def _read_barcodes(path: Path) -> list[str]:
    import gzip

    with gzip.open(path, "rt") as f:
        return [line.rstrip("\n") for line in f]


def load_mtx(prefix: Path):
    """Load a 10x-style MTX folder/prefix. Returns genes, barcodes, csr cells x genes."""
    import gzip

    feat = prefix.parent / (prefix.name + "_features.tsv.gz")
    if not feat.exists():
        feat = prefix.parent / (prefix.name + "_genes.tsv.gz")
    bar = prefix.parent / (prefix.name + "_barcodes.tsv.gz")
    mtx = prefix.parent / (prefix.name + "_matrix.mtx.gz")
    genes = _read_features(feat)
    barcodes = _read_barcodes(bar)
    # MTX header
    with gzip.open(mtx, "rt") as f:
        line = f.readline()
        while line.startswith("%"):
            line = f.readline()
        n_row, n_col, n_nz = [int(x) for x in line.split()]
        rows, cols, data = [], [], []
        for line in f:
            a, b, c = line.split()
            rows.append(int(a) - 1)
            cols.append(int(b) - 1)
            data.append(float(c))
    # Cell Ranger: genes x cells
    if n_row == len(genes) and n_col == len(barcodes):
        X = sparse.csr_matrix((data, (cols, rows)), shape=(len(barcodes), len(genes)))
    elif n_row == len(barcodes) and n_col == len(genes):
        X = sparse.csr_matrix((data, (rows, cols)), shape=(len(barcodes), len(genes)))
    else:
        raise ValueError(
            f"shape mismatch {prefix.name}: mtx {n_row}x{n_col} genes={len(genes)} bars={len(barcodes)}"
        )
    return genes, barcodes, X


def gene_index(genes: list[str]) -> dict[str, int]:
    idx = {}
    for i, g in enumerate(genes):
        g = g.split(".")[0]
        if g not in idx:
            idx[g] = i
    return idx


def col(X, idx, name):
    j = idx.get(name)
    if j is None:
        return None
    return np.asarray(X[:, j].todense()).ravel()


def mean_score(X_log, idx, names):
    cols = [idx[n] for n in names if n in idx]
    if not cols:
        return np.zeros(X_log.shape[0])
    return np.asarray(X_log[:, cols].mean(axis=1)).ravel()


def annotate(X, genes, min_umi=200, min_genes=100):
    idx = gene_index(genes)
    umi = np.asarray(X.sum(axis=1)).ravel()
    ngenes = np.asarray((X > 0).sum(axis=1)).ravel()
    keep = (umi >= min_umi) & (ngenes >= min_genes)
    X = X[keep]
    umi = umi[keep]
    # CP10k log1p
    scale = 1e4 / np.maximum(umi, 1.0)
    X_log = X.multiply(scale[:, None])
    X_log.data = np.log1p(X_log.data)
    X_log = X_log.tocsr()
    scores = {k: mean_score(X_log, idx, v) for k, v in MARKERS.items()}
    score_mat = np.vstack([scores[k] for k in MARKERS])
    lineage = np.array(list(MARKERS.keys()))[score_mat.argmax(axis=0)]
    # T/NK combined
    tnk = (lineage == "t") | (lineage == "nk")
    epi = lineage == "epithelial"
    normal = mean_score(X_log, idx, NORMAL_LUNG)
    mal = epi & (normal < np.median(normal[epi]) if epi.any() else False)
    # if few epi, skip median split and treat all epi as epi (mal proxy = epi)
    if epi.sum() < 20:
        mal = epi
    out = {
        "n_qc": int(keep.sum()),
        "n_input": int(keep.size if hasattr(keep, "size") else len(keep)),
        "lineage": lineage,
        "tnk": tnk,
        "epi": epi,
        "mal": mal,
        "X_log": X_log,
        "idx": idx,
        "umi": umi,
    }
    for g in TARGETS:
        v = col(X_log, idx, g)
        out[g] = v if v is not None else np.full(X_log.shape[0], np.nan)
    return out


def sample_metrics(ann, sample_id, extra=None):
    n = ann["n_qc"]
    epi = ann["epi"]
    mal = ann["mal"]
    tnk = ann["tnk"]
    row = {
        "sample": sample_id,
        "n_qc": n,
        "n_epi": int(epi.sum()),
        "n_mal": int(mal.sum()),
        "n_tnk": int(tnk.sum()),
        "frac_epi": epi.mean() if n else np.nan,
        "frac_mal": mal.mean() if n else np.nan,
        "frac_tnk": tnk.mean() if n else np.nan,
        "epi_TACSTD2": float(np.nanmean(ann["TACSTD2"][epi])) if epi.any() else np.nan,
        "epi_CLDN4": float(np.nanmean(ann["CLDN4"][epi])) if epi.any() else np.nan,
        "mal_TACSTD2": float(np.nanmean(ann["TACSTD2"][mal])) if mal.any() else np.nan,
        "mal_CLDN4": float(np.nanmean(ann["CLDN4"][mal])) if mal.any() else np.nan,
        "tnk_TACSTD2": float(np.nanmean(ann["TACSTD2"][tnk])) if tnk.any() else np.nan,
        "tnk_CLDN4": float(np.nanmean(ann["CLDN4"][tnk])) if tnk.any() else np.nan,
        "epi_CD274": float(np.nanmean(ann["CD274"][epi])) if epi.any() else np.nan,
    }
    if extra:
        row.update(extra)
    return row


def spearman(df, x, y, min_n=4, extra_filter=None):
    d = df.copy()
    if extra_filter is not None:
        d = d.loc[extra_filter(d)]
    d = d.dropna(subset=[x, y])
    n = len(d)
    if n < min_n:
        return {"n": n, "rho": np.nan, "p": np.nan, "note": f"n<{min_n}"}
    rho, p = stats.spearmanr(d[x], d[y])
    return {"n": n, "rho": float(rho), "p": float(p), "note": ""}


def mwu(df, value, group, a, b):
    d = df.dropna(subset=[value, group])
    xa = d.loc[d[group] == a, value]
    xb = d.loc[d[group] == b, value]
    if len(xa) < 2 or len(xb) < 2:
        return {
            "n_a": int(len(xa)),
            "n_b": int(len(xb)),
            "median_a": float(xa.median()) if len(xa) else np.nan,
            "median_b": float(xb.median()) if len(xb) else np.nan,
            "U": np.nan,
            "p": np.nan,
            "note": "n<2 in a group",
        }
    U, p = stats.mannwhitneyu(xa, xb, alternative="two-sided")
    return {
        "n_a": int(len(xa)),
        "n_b": int(len(xb)),
        "median_a": float(xa.median()),
        "median_b": float(xb.median()),
        "U": float(U),
        "p": float(p),
        "note": "",
    }


def score_gse274595():
    rows = []
    d = DATA / "GSE274595"
    for gsm, meta in GSE274595_META.items():
        feats = list(d.glob(f"{gsm}_*_features.tsv.gz"))
        if not feats:
            continue
        prefix = Path(str(feats[0]).replace("_features.tsv.gz", ""))
        genes, bars, X = load_mtx(prefix)
        ann = annotate(X, genes)
        rows.append(sample_metrics(ann, gsm, {**meta, "dataset": "GSE274595"}))
        print("GSE274595", gsm, meta["patient"], rows[-1]["n_qc"], rows[-1]["n_epi"], rows[-1]["n_tnk"])
    return pd.DataFrame(rows)


def score_gse267108():
    rows = []
    d = DATA / "GSE267108"
    for feat in sorted(d.glob("*_features.tsv.gz")):
        sid = feat.name.replace("_features.tsv.gz", "")
        pid = re.sub(r"1$", "", sid)
        genes, bars, X = load_mtx(feat.parent / sid)
        ann = annotate(X, genes)
        extra = {
            "dataset": "GSE267108",
            "patient": pid,
            "pdl1_inferred": GSE267108_PDL1.get(pid, "unknown"),
            "histology": "LUAD",
            "assay": "BD_Rhapsody_tumor",
        }
        rows.append(sample_metrics(ann, sid, extra))
        print("GSE267108", sid, rows[-1]["n_qc"], rows[-1]["n_epi"], rows[-1]["n_tnk"], extra["pdl1_inferred"])
    return pd.DataFrame(rows)


def score_gse337519():
    d = DATA / "GSE337519"
    prefix = d / "GSM9856929"
    genes, bars, X = load_mtx(prefix)
    ann = annotate(X, genes)
    row = sample_metrics(
        ann,
        "GSM9856929",
        {
            "dataset": "GSE337519",
            "patient": "luad1",
            "histology": "LUAD",
            "assay": "paired_pre_post_merged_or_single",
            "n_barcodes_raw": len(bars),
        },
    )
    print("GSE337519", row["n_qc"], row["n_epi"], row["n_tnk"])
    return pd.DataFrame([row])


def score_gse186446_bulk():
    """Bulk featureCounts companion of the 3-pt ICB T-cell scRNA series. Not scRNA."""
    p = Path("/workspace/methods/scrna_leftover_2426/metadata/GSE186446_fcount_aggr.txt.gz")
    df = pd.read_csv(p, sep="\t")
    df.columns = [c.strip() for c in df.columns]
    gene = df["gene"].astype(str)
    mat = df.drop(columns=["id", "gene"]).apply(pd.to_numeric, errors="coerce")
    mat.index = gene.values
    # collapse duplicate symbols
    mat = mat.groupby(level=0).sum()
    # CPM-like log1p
    lib = mat.sum(axis=0).replace(0, np.nan)
    cpm = mat.divide(lib, axis=1) * 1e6
    logcpm = np.log1p(cpm)
    samples = []
    meta_path = Path("/workspace/methods/scrna_leftover_2426/metadata/GSE186446_samples.json")
    meta_by_title = {r["title"].replace(" ", ""): r for r in json.loads(meta_path.read_text())}
    for col in logcpm.columns:
        key = col.replace(" ", "")
        pid = key.split("_")[0]
        info = meta_by_title.get(key)
        if info is None and key.startswith("1344_C"):
            # fcount uses C1–C8; GEO titles are 1344_R1–R8, all viable tumor
            info = {
                "patient number": "MSK 1344",
                "involved by tumor?": "Yes",
                "region characterization": "Viable Tumor (C-core; GEO titles are R1–R8)",
            }
        rec = {
            "sample": col,
            "dataset": "GSE186446_bulk_fcount",
            "patient": (info or {}).get("patient number", f"MSK {pid}"),
            "viable_tumor": (info or {}).get("involved by tumor?", ""),
            "region": (info or {}).get("region characterization", ""),
        }
        for g in ["TACSTD2", "CLDN4", "EPCAM", "KRT8", "CD3D", "CD3E", "CD8A", "NKG7", "GNLY", "PTPRC"]:
            rec[g] = float(logcpm.loc[g, col]) if g in logcpm.index else np.nan
        rec["tnk_score"] = np.nanmean([rec["CD3D"], rec["CD8A"], rec["NKG7"]])
        samples.append(rec)
    return pd.DataFrame(samples)


def stats_block(name, df, x_cols, y="frac_tnk", group=None, group_a=None, group_b=None, min_epi=10, min_tnk=20):
    elig = df[(df["n_epi"] >= min_epi) & (df["n_tnk"] >= min_tnk)].copy()
    rows = []
    for x in x_cols:
        s = spearman(elig, x, y)
        rows.append(
            {
                "dataset": name,
                "analysis": "spearman_patient",
                "comparison": f"{x} vs {y}",
                "metric": x,
                "n": s["n"],
                "stat": "spearman_rho",
                "value": s["rho"],
                "p_value": s["p"],
                "note": s["note"]
                + f"; eligible n_epi>={min_epi} and n_tnk>={min_tnk}; unit=row",
            }
        )
    if group and group_a and group_b:
        for x in x_cols:
            m = mwu(elig, x, group, group_a, group_b)
            rows.append(
                {
                    "dataset": name,
                    "analysis": f"mwu_{group}",
                    "comparison": f"{x} {group_a} vs {group_b}",
                    "metric": x,
                    "n": f"{m['n_a']}+{m['n_b']}",
                    "stat": "mannwhitney_u",
                    "value": m["U"],
                    "p_value": m["p"],
                    "note": f"{group_a} med={m['median_a']:.4f} (n={m['n_a']}); {group_b} med={m['median_b']:.4f} (n={m['n_b']}); {m['note']}",
                }
            )
    return rows, elig


def main():
    stats_rows = []

    g595 = score_gse274595()
    g595.to_csv(OUT / "GSE274595_per_sample.tsv", sep="\t", index=False)
    r, elig = stats_block("GSE274595", g595, ["epi_TACSTD2", "epi_CLDN4", "mal_TACSTD2", "mal_CLDN4"])
    stats_rows.extend(r)
    # tissue-only sensitivity
    r2, _ = stats_block(
        "GSE274595_tissue",
        g595[g595["assay"] == "tissue"],
        ["epi_TACSTD2", "epi_CLDN4"],
        min_epi=5,
        min_tnk=10,
    )
    for x in r2:
        x["note"] += "; tissue lanes only (nuclei excluded)"
    stats_rows.extend(r2)

    g267 = score_gse267108()
    g267.to_csv(OUT / "GSE267108_per_sample.tsv", sep="\t", index=False)
    g267 = g267.copy()
    g267["pdl1_bin"] = g267["pdl1_inferred"].map(lambda x: "neg" if x == "neg" else "pos_or_high")
    r, _ = stats_block(
        "GSE267108",
        g267,
        ["epi_TACSTD2", "epi_CLDN4", "mal_TACSTD2", "mal_CLDN4"],
        group="pdl1_bin",
        group_a="neg",
        group_b="pos_or_high",
    )
    stats_rows.extend(r)

    g337 = score_gse337519()
    g337.to_csv(OUT / "GSE337519_per_sample.tsv", sep="\t", index=False)
    stats_rows.append(
        {
            "dataset": "GSE337519",
            "analysis": "descriptive_n1",
            "comparison": "epithelial TACSTD2 / CLDN4 vs T/NK (n=1 patient)",
            "metric": "epi_TACSTD2",
            "n": 1,
            "stat": "mean_log1p_cp10k",
            "value": float(g337.iloc[0]["epi_TACSTD2"]),
            "p_value": np.nan,
            "note": f"n=1; cannot test response or patient-level ρ. n_epi={int(g337.iloc[0]['n_epi'])} n_tnk={int(g337.iloc[0]['n_tnk'])} frac_tnk={g337.iloc[0]['frac_tnk']:.4f} epi_CLDN4={g337.iloc[0]['epi_CLDN4']:.4f}",
        }
    )

    bulk = score_gse186446_bulk()
    bulk.to_csv(OUT / "GSE186446_bulk_per_sample.tsv", sep="\t", index=False)
    # viable tumor regions only; NOT scRNA
    vt = bulk[bulk["viable_tumor"].astype(str).str.lower().eq("yes")].copy()
    s = spearman(vt, "TACSTD2", "tnk_score", min_n=4)
    stats_rows.append(
        {
            "dataset": "GSE186446_bulk_fcount",
            "analysis": "spearman_region_bulk",
            "comparison": "bulk TACSTD2 vs T/NK score (viable tumor regions)",
            "metric": "TACSTD2",
            "n": s["n"],
            "stat": "spearman_rho",
            "value": s["rho"],
            "p_value": s["p"],
            "note": "NOT scRNA. Public processed file is bulk featureCounts of 32 regions / 3 ICB patients. scRNA in this SuperSeries is T-cell only.",
        }
    )
    s = spearman(vt, "CLDN4", "tnk_score", min_n=4)
    stats_rows.append(
        {
            "dataset": "GSE186446_bulk_fcount",
            "analysis": "spearman_region_bulk",
            "comparison": "bulk CLDN4 vs T/NK score (viable tumor regions)",
            "metric": "CLDN4",
            "n": s["n"],
            "stat": "spearman_rho",
            "value": s["rho"],
            "p_value": s["p"],
            "note": "NOT scRNA. Region unit; 3 patients so patient-unit n=3 is underpowered.",
        }
    )
    # patient means
    pat = (
        vt.groupby("patient")[["TACSTD2", "CLDN4", "tnk_score"]]
        .mean()
        .reset_index()
    )
    pat.to_csv(OUT / "GSE186446_bulk_per_patient.tsv", sep="\t", index=False)
    s = spearman(pat, "TACSTD2", "tnk_score", min_n=3)
    stats_rows.append(
        {
            "dataset": "GSE186446_bulk_fcount",
            "analysis": "spearman_patient_bulk",
            "comparison": "bulk TACSTD2 vs T/NK (patient mean of viable tumor)",
            "metric": "TACSTD2",
            "n": s["n"],
            "stat": "spearman_rho",
            "value": s["rho"],
            "p_value": s["p"],
            "note": f"n=3 patients. {s['note']}",
        }
    )
    s = spearman(pat, "CLDN4", "tnk_score", min_n=3)
    stats_rows.append(
        {
            "dataset": "GSE186446_bulk_fcount",
            "analysis": "spearman_patient_bulk",
            "comparison": "bulk CLDN4 vs T/NK (patient mean of viable tumor)",
            "metric": "CLDN4",
            "n": s["n"],
            "stat": "spearman_rho",
            "value": s["rho"],
            "p_value": s["p"],
            "note": f"n=3 patients. {s['note']}",
        }
    )

    stats_df = pd.DataFrame(stats_rows)
    stats_df.to_csv(OUT / "stats.tsv", sep="\t", index=False)
    print(stats_df.to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
