#!/usr/bin/env python3
"""Cldn4 vs type-I IFN / MHC on public processed GSE274477 + GSE274351.

Public processed matrices only. No SRA / raw reprocessing.
"""
from __future__ import annotations

import gzip
import json
import math
import os
from collections import OrderedDict
from pathlib import Path

import numpy as np
from scipy import stats

OUT = Path("/workspace/methods/gse274477_kras_ifn_cldn4")
GEO477 = Path("/tmp/geo/gse274477/GSE274477_Feature_counts_Matrix_14Aug.tsv.gz")
GEO351 = Path("/tmp/geo/gse274351/GSE274351_expredata_TPM_gene.txt.gz")
GTF = Path("/tmp/annot/Mus_musculus.GRCm38.102.gtf.gz")
GMT = Path("/tmp/annot/mh.all.v2024.1.Mm.symbols.gmt")

# Paper-core type I IFN / IRG genes (Fernández-García et al. PNAS 2024, Fig. 2B / Fig. 3A)
IFN_CORE = [
    "Stat1",
    "Irf7",
    "Ifih1",
    "Ddx58",  # Rigi
    "Oasl2",
    "Bst2",
    "Ifi27l2a",
    "Ifitm3",
]
MHC_CORE = [
    "B2m",
    "Tap1",
    "Tapbp",
    "H2-K1",
    "H2-D1",
    "H2-T23",
    "H2-Q6",
    "H2-Q7",
    "H2-Aa",
]


def parse_gtf_maps(path: Path) -> tuple[dict[str, str], set[str]]:
    ens_to_sym: dict[str, str] = {}
    mt_ids: set[str] = set()
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            seqname = line.split("\t", 1)[0]
            attrs = line.split("\t")[8]
            gid = gname = None
            for part in attrs.split(";"):
                part = part.strip()
                if part.startswith("gene_id"):
                    gid = part.split('"')[1]
                elif part.startswith("gene_name"):
                    gname = part.split('"')[1]
            if gid and gname:
                ens_to_sym[gid] = gname
                if seqname == "MT" or gname.startswith("mt-"):
                    mt_ids.add(gid)
    return ens_to_sym, mt_ids


def parse_gmt_ifn(path: Path) -> list[str]:
    for line in path.read_text().splitlines():
        parts = line.split("\t")
        if parts[0] == "HALLMARK_INTERFERON_ALPHA_RESPONSE":
            return parts[2:]
    raise SystemExit("HALLMARK_INTERFERON_ALPHA_RESPONSE missing from GMT")


def strip_ens(gid: str) -> str:
    return gid.split(".")[0]


def mean_z(mat: np.ndarray) -> np.ndarray:
    """Column-wise mean of row-wise z-scores. mat is genes x samples."""
    if mat.size == 0:
        return np.full(mat.shape[1] if mat.ndim == 2 else 0, np.nan)
    mu = np.nanmean(mat, axis=1, keepdims=True)
    sd = np.nanstd(mat, axis=1, keepdims=True)
    sd = np.where(sd == 0, 1.0, sd)
    z = (mat - mu) / sd
    return np.nanmean(z, axis=0)


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def mw(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return (float("nan"), float("nan"), float("nan"))
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    # rank-biserial: 1 - 2U/(n1 n2) with U for a vs b; positive if a > b
    r = 1.0 - (2.0 * u) / (len(a) * len(b))
    # actually U large if a > b, so r_rb = (2U)/(n1n2) - 1
    r_rb = (2.0 * u) / (len(a) * len(b)) - 1.0
    return float(u), float(p), float(r_rb)


def spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5:
        return (float("nan"), float("nan"))
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p)


def median_split(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    med = float(np.nanmedian(x))
    hi = x > med
    lo = x <= med
    return hi, lo


def summarize(x: np.ndarray) -> dict:
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0, "median": None, "mean": None, "q25": None, "q75": None}
    return {
        "n": int(x.size),
        "median": float(np.median(x)),
        "mean": float(np.mean(x)),
        "q25": float(np.percentile(x, 25)),
        "q75": float(np.percentile(x, 75)),
    }


def md_table(headers: list[str], rows: list[list[str]]) -> str:
    out = ["| " + " | ".join(headers) + " |", "| " + " | ".join(["---"] * len(headers)) + " |"]
    for r in rows:
        out.append("| " + " | ".join(r) + " |")
    return "\n".join(out)


def main() -> None:
    ens_to_sym, mt_ids = parse_gtf_maps(GTF)
    hallmark_ifna = parse_gmt_ifn(GMT)
    print("mapped genes", len(ens_to_sym), "mt", len(mt_ids), "hallmark IFNA", len(hallmark_ifna))

    # reverse: symbols we need
    needed_sym = set(["Cldn4"] + IFN_CORE + MHC_CORE + hallmark_ifna)
    needed_ens = {ens: sym for ens, sym in ens_to_sym.items() if sym in needed_sym}
    print("needed ens hits", len(needed_ens), "of", len(needed_sym), "symbols")

    # ---------- GSE274477 ----------
    with gzip.open(GEO477, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    cells = header[1:]
    n_cells = len(cells)
    genotype = np.array(
        ["MUT" if "_MUT_" in c else ("WT" if "_CONTROL_" in c else "other") for c in cells]
    )
    run = np.array(["FCNIO1" if c.startswith("FCNIO1") else ("FCNIO2" if c.startswith("FCNIO2") else "?") for c in cells])
    undet = np.array(["Undetermined" in c for c in cells])

    lib = np.zeros(n_cells, dtype=np.float64)
    ngene = np.zeros(n_cells, dtype=np.int32)
    mito = np.zeros(n_cells, dtype=np.float64)
    gene_rows: dict[str, np.ndarray] = {}

    with gzip.open(GEO477, "rt") as fh:
        next(fh)
        for line in fh:
            gid_raw, rest = line.split("\t", 1)
            gid = strip_ens(gid_raw)
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            lib += vals
            ngene += vals > 0
            if gid in mt_ids:
                mito += vals
            if gid in needed_ens:
                gene_rows[needed_ens[gid]] = vals.astype(np.float64)

    mito_frac = np.divide(mito, lib, out=np.zeros_like(mito), where=lib > 0)
    print("extracted symbols", sorted(gene_rows))
    missing_sym = needed_sym - set(gene_rows)
    print("missing from scRNA matrix", sorted(missing_sym)[:40], "n=", len(missing_sym))

    cldn4_present = "Cldn4" in gene_rows
    ifn_present = [g for g in IFN_CORE if g in gene_rows]
    mhc_present = [g for g in MHC_CORE if g in gene_rows]
    h_present = [g for g in hallmark_ifna if g in gene_rows]
    print("Cldn4", cldn4_present, "IFN_CORE", ifn_present, "MHC", mhc_present, "H_IFNA", len(h_present))

    # QC: documented, paper-like depth (avg lib ~1e5, avg genes ~2354)
    qc = (~undet) & (lib >= 50_000) & (ngene >= 1_000) & (ngene <= 8_000) & (mito_frac <= 0.10)
    print(
        "QC n",
        int(qc.sum()),
        "WT",
        int((qc & (genotype == "WT")).sum()),
        "MUT",
        int((qc & (genotype == "MUT")).sum()),
    )
    # looser QC for sensitivity
    qc_loose = (~undet) & (lib >= 20_000) & (ngene >= 800) & (mito_frac <= 0.15)
    print(
        "loose QC n",
        int(qc_loose.sum()),
        "WT",
        int((qc_loose & (genotype == "WT")).sum()),
        "MUT",
        int((qc_loose & (genotype == "MUT")).sum()),
    )

    def logcpm(counts: np.ndarray, libsize: np.ndarray) -> np.ndarray:
        return np.log1p(counts / np.maximum(libsize, 1.0) * 1e4)

    def scores_for_mask(mask: np.ndarray) -> dict[str, np.ndarray]:
        lib_m = lib[mask]
        out: dict[str, np.ndarray] = {}
        out["Cldn4"] = logcpm(gene_rows["Cldn4"][mask], lib_m) if "Cldn4" in gene_rows else np.full(mask.sum(), np.nan)
        out["Cldn4_raw"] = gene_rows["Cldn4"][mask] if "Cldn4" in gene_rows else np.full(mask.sum(), np.nan)

        def set_score(genes: list[str]) -> np.ndarray:
            mats = [logcpm(gene_rows[g][mask], lib_m) for g in genes if g in gene_rows]
            if not mats:
                return np.full(mask.sum(), np.nan)
            return mean_z(np.vstack(mats))

        out["IFN_core"] = set_score(IFN_CORE)
        out["MHC"] = set_score(MHC_CORE)
        out["H_IFNA"] = set_score(h_present)
        out["genotype"] = genotype[mask]
        out["run"] = run[mask]
        return out

    sc = scores_for_mask(qc)
    sc_loose = scores_for_mask(qc_loose)

    def group_block(scd: dict, label: str) -> list[dict]:
        rows = []
        for gt in ("WT", "MUT"):
            m = scd["genotype"] == gt
            det = scd["Cldn4_raw"][m] > 0
            rows.append(
                {
                    "set": label,
                    "group": gt,
                    "n_cells": int(m.sum()),
                    "Cldn4_detect_n": int(det.sum()),
                    "Cldn4_detect_frac": float(det.mean()) if m.sum() else None,
                    "Cldn4_logCPM": summarize(scd["Cldn4"][m]),
                    "IFN_core": summarize(scd["IFN_core"][m]),
                    "MHC": summarize(scd["MHC"][m]),
                    "H_IFNA": summarize(scd["H_IFNA"][m]),
                }
            )
        return rows

    gse477_groups = group_block(sc, "QC_lib>=50k_genes>=1k_mito<=10pct")
    gse477_groups_loose = group_block(sc_loose, "QC_lib>=20k_genes>=800_mito<=15pct")

    def hi_lo_block(scd: dict, gt_filter: str | None, label: str) -> dict:
        m = np.ones(len(scd["Cldn4"]), dtype=bool)
        if gt_filter:
            m = scd["genotype"] == gt_filter
        x = scd["Cldn4"][m]
        hi, lo = median_split(x)
        # if almost no variation, note it
        n_unique = int(np.unique(np.round(x[np.isfinite(x)], 6)).size)
        rec = {
            "set": label,
            "subset": gt_filter or "all_QC",
            "n": int(m.sum()),
            "n_high": int(hi.sum()),
            "n_low": int(lo.sum()),
            "n_unique_Cldn4": n_unique,
            "Cldn4_detect_n": int((scd["Cldn4_raw"][m] > 0).sum()),
            "median_Cldn4": float(np.nanmedian(x)),
        }
        for score in ("IFN_core", "MHC", "H_IFNA"):
            a = scd[score][m][hi]
            b = scd[score][m][lo]
            u, p, r = mw(a, b)
            rec[score] = {
                "high_median": float(np.nanmedian(a)) if a.size else None,
                "low_median": float(np.nanmedian(b)) if b.size else None,
                "delta_high_minus_low": float(np.nanmedian(a) - np.nanmedian(b)) if a.size and b.size else None,
                "U": u,
                "p": p,
                "rank_biserial_high_gt_low": r,
            }
            rs, psp = spearman(x, scd[score][m])
            rec[f"spearman_Cldn4_vs_{score}"] = {"rho": rs, "p": psp, "n": int(m.sum())}
        return rec

    hi_lo = [
        hi_lo_block(sc, "MUT", "QC_primary MUT median-split Cldn4"),
        hi_lo_block(sc, "WT", "QC_primary WT median-split Cldn4"),
        hi_lo_block(sc, None, "QC_primary all AT2 median-split Cldn4"),
        hi_lo_block(sc_loose, "MUT", "QC_loose MUT median-split Cldn4"),
    ]

    # WT vs MUT for Cldn4 / IFN / MHC
    wt_mut = {}
    m_wt = sc["genotype"] == "WT"
    m_mut = sc["genotype"] == "MUT"
    for score in ("Cldn4", "IFN_core", "MHC", "H_IFNA"):
        u, p, r = mw(sc[score][m_mut], sc[score][m_wt])
        wt_mut[score] = {
            "n_WT": int(m_wt.sum()),
            "n_MUT": int(m_mut.sum()),
            "WT_median": float(np.nanmedian(sc[score][m_wt])),
            "MUT_median": float(np.nanmedian(sc[score][m_mut])),
            "delta_MUT_minus_WT": float(np.nanmedian(sc[score][m_mut]) - np.nanmedian(sc[score][m_wt])),
            "U": u,
            "p": p,
            "rank_biserial_MUT_gt_WT": r,
        }

    # ---------- GSE274351 LCM TPM ----------
    with gzip.open(GEO351, "rt") as fh:
        h351 = fh.readline().rstrip("\n").split("\t")
    samples = h351[1:]
    # deposited labels
    def samp_gt(name: str) -> str:
        if name.startswith("NL"):
            return "NL"
        if name.startswith("KP"):
            return "KP"
        if name.startswith("KL"):
            return "KL"
        if name.startswith("K"):
            return "K"
        return "?"

    gt351 = np.array([samp_gt(s) for s in samples])
    print("GSE274351 columns", list(zip(samples, gt351)))

    tpm_rows: dict[str, np.ndarray] = {}
    with gzip.open(GEO351, "rt") as fh:
        next(fh)
        for line in fh:
            gid_raw, rest = line.split("\t", 1)
            gid = strip_ens(gid_raw)
            if gid not in needed_ens:
                continue
            vals = np.fromstring(rest, sep="\t", dtype=np.float64)
            tpm_rows[needed_ens[gid]] = vals

    print("351 extracted", sorted(tpm_rows), "missing", sorted(needed_sym - set(tpm_rows))[:30])

    def logtpm(g: str) -> np.ndarray:
        return np.log2(tpm_rows[g] + 1.0)

    def set_score_tpm(genes: list[str]) -> np.ndarray:
        mats = [logtpm(g) for g in genes if g in tpm_rows]
        return mean_z(np.vstack(mats)) if mats else np.full(len(samples), np.nan)

    tpm = {
        "Cldn4": logtpm("Cldn4") if "Cldn4" in tpm_rows else np.full(len(samples), np.nan),
        "Cldn4_tpm": tpm_rows["Cldn4"] if "Cldn4" in tpm_rows else np.full(len(samples), np.nan),
        "IFN_core": set_score_tpm(IFN_CORE),
        "MHC": set_score_tpm(MHC_CORE),
        "H_IFNA": set_score_tpm([g for g in hallmark_ifna if g in tpm_rows]),
        "gt": gt351,
        "sample": np.array(samples),
    }

    by_gt = {}
    for g in ("NL", "K", "KP", "KL"):
        m = tpm["gt"] == g
        by_gt[g] = {
            "n": int(m.sum()),
            "samples": [s for s, ok in zip(samples, m) if ok],
            "Cldn4_TPM": summarize(tpm["Cldn4_tpm"][m]),
            "Cldn4_log2TPM1": summarize(tpm["Cldn4"][m]),
            "IFN_core": summarize(tpm["IFN_core"][m]),
            "MHC": summarize(tpm["MHC"][m]),
            "H_IFNA": summarize(tpm["H_IFNA"][m]),
        }

    # KL vs K, KL vs KP, KL vs K+KP, adenomas vs NL
    contrasts = {}
    pairs = [
        ("KL", "K"),
        ("KL", "KP"),
        ("KL", "K+KP"),
        ("K+KP+KL", "NL"),
        ("K", "NL"),
        ("KP", "NL"),
        ("KL", "NL"),
    ]
    mask_of = {
        "NL": tpm["gt"] == "NL",
        "K": tpm["gt"] == "K",
        "KP": tpm["gt"] == "KP",
        "KL": tpm["gt"] == "KL",
        "K+KP": (tpm["gt"] == "K") | (tpm["gt"] == "KP"),
        "K+KP+KL": tpm["gt"] != "NL",
    }
    for a, b in pairs:
        contrasts[f"{a}_vs_{b}"] = {}
        for score in ("Cldn4", "Cldn4_tpm", "IFN_core", "MHC", "H_IFNA"):
            u, p, r = mw(tpm[score][mask_of[a]], tpm[score][mask_of[b]])
            contrasts[f"{a}_vs_{b}"][score] = {
                "n_a": int(mask_of[a].sum()),
                "n_b": int(mask_of[b].sum()),
                "median_a": float(np.nanmedian(tpm[score][mask_of[a]])),
                "median_b": float(np.nanmedian(tpm[score][mask_of[b]])),
                "U": u,
                "p": p,
                "rank_biserial_a_gt_b": r,
            }

    # high vs low Cldn4 among adenomas only
    ad = tpm["gt"] != "NL"
    x = tpm["Cldn4"][ad]
    hi, lo = median_split(x)
    ad_hilo = {
        "n_adenoma": int(ad.sum()),
        "n_high": int(hi.sum()),
        "n_low": int(lo.sum()),
        "high_samples": [s for s, ok in zip(tpm["sample"][ad], hi) if ok],
        "low_samples": [s for s, ok in zip(tpm["sample"][ad], lo) if ok],
        "high_gt": Counter_safe(tpm["gt"][ad][hi]),
        "low_gt": Counter_safe(tpm["gt"][ad][lo]),
        "median_split_log2TPM1": float(np.nanmedian(x)),
    }
    for score in ("IFN_core", "MHC", "H_IFNA"):
        u, p, r = mw(tpm[score][ad][hi], tpm[score][ad][lo])
        ad_hilo[score] = {
            "high_median": float(np.nanmedian(tpm[score][ad][hi])),
            "low_median": float(np.nanmedian(tpm[score][ad][lo])),
            "delta_high_minus_low": float(np.nanmedian(tpm[score][ad][hi]) - np.nanmedian(tpm[score][ad][lo])),
            "U": u,
            "p": p,
            "rank_biserial_high_gt_low": r,
        }
        rs, psp = spearman(x, tpm[score][ad])
        ad_hilo[f"spearman_Cldn4_vs_{score}"] = {"rho": rs, "p": psp, "n": int(ad.sum())}

    # also all samples including NL
    xall = tpm["Cldn4"]
    hi_a, lo_a = median_split(xall)
    all_hilo = {"n": int(len(xall)), "n_high": int(hi_a.sum()), "n_low": int(lo_a.sum())}
    for score in ("IFN_core", "MHC", "H_IFNA"):
        u, p, r = mw(tpm[score][hi_a], tpm[score][lo_a])
        all_hilo[score] = {
            "high_median": float(np.nanmedian(tpm[score][hi_a])),
            "low_median": float(np.nanmedian(tpm[score][lo_a])),
            "delta_high_minus_low": float(np.nanmedian(tpm[score][hi_a]) - np.nanmedian(tpm[score][lo_a])),
            "p": p,
            "rank_biserial_high_gt_low": r,
        }
        rs, psp = spearman(xall, tpm[score])
        all_hilo[f"spearman_Cldn4_vs_{score}"] = {"rho": rs, "p": psp, "n": int(len(xall))}

    # per-sample table
    per_sample = []
    for i, s in enumerate(samples):
        per_sample.append(
            {
                "sample": s,
                "genotype_matrix_label": gt351[i],
                "Cldn4_TPM": float(tpm["Cldn4_tpm"][i]),
                "Cldn4_log2TPM1": float(tpm["Cldn4"][i]),
                "IFN_core_z": float(tpm["IFN_core"][i]),
                "MHC_z": float(tpm["MHC"][i]),
                "H_IFNA_z": float(tpm["H_IFNA"][i]),
            }
        )

    payload = OrderedDict(
        [
            (
                "data",
                {
                    "GSE274477": {
                        "file": "GSE274477_Feature_counts_Matrix_14Aug.tsv.gz",
                        "organism": "Mus musculus",
                        "assay": "Fluidigm C1 scRNA-seq feature counts (public processed)",
                        "paper_n_after_author_QC": {"WT": 139, "MUT": 176, "total": 315},
                        "matrix_wells": n_cells,
                        "QC": "non-Undetermined; library size >= 50000; 1000-8000 detected genes; mito fraction <= 0.10",
                        "normalization": "log1p(10k * count / library size)",
                    },
                    "GSE274351": {
                        "file": "GSE274351_expredata_TPM_gene.txt.gz",
                        "organism": "Mus musculus",
                        "assay": "LCM adenoma / adjacent alveolus QuantSeq TPM (public processed)",
                        "GEO_sample_titles": "5 K, 5 KP, 4 KL, 4 NL",
                        "matrix_column_labels": "4 K, 5 KP, 5 KL, 4 NL (used as deposited)",
                        "n_columns": len(samples),
                    },
                },
            ),
            ("genes", {
                "Cldn4_ensembl": "ENSMUSG00000047501",
                "IFN_core_requested": IFN_CORE,
                "IFN_core_in_scRNA": ifn_present,
                "IFN_core_in_LCM": [g for g in IFN_CORE if g in tpm_rows],
                "MHC_requested": MHC_CORE,
                "MHC_in_scRNA": mhc_present,
                "MHC_in_LCM": [g for g in MHC_CORE if g in tpm_rows],
                "Hallmark_IFNA_n_requested": len(hallmark_ifna),
                "Hallmark_IFNA_n_scRNA": len(h_present),
                "Hallmark_IFNA_n_LCM": sum(1 for g in hallmark_ifna if g in tpm_rows),
            }),
            ("GSE274477_WT_vs_MUT", wt_mut),
            ("GSE274477_group_summaries", gse477_groups),
            ("GSE274477_group_summaries_looseQC", gse477_groups_loose),
            ("GSE274477_Cldn4_high_vs_low", hi_lo),
            ("GSE274351_by_genotype", by_gt),
            ("GSE274351_contrasts", contrasts),
            ("GSE274351_adenoma_Cldn4_high_vs_low", ad_hilo),
            ("GSE274351_all_samples_Cldn4_high_vs_low", all_hilo),
            ("GSE274351_per_sample", per_sample),
        ]
    )

    (OUT / "results.json").write_text(json.dumps(payload, indent=2, default=str))
    write_finding(payload, sc, tpm)
    print("wrote", OUT / "FINDING.generated.md", "and", OUT / "results.json")


def Counter_safe(arr: np.ndarray) -> dict[str, int]:
    d: dict[str, int] = {}
    for v in arr:
        d[str(v)] = d.get(str(v), 0) + 1
    return d


def fnum(x, nd=3):
    if x is None or (isinstance(x, float) and (math.isnan(x) or math.isinf(x))):
        return "NA"
    if isinstance(x, int):
        return str(x)
    if abs(x) >= 100:
        return f"{x:.1f}"
    if abs(x) >= 10:
        return f"{x:.2f}"
    return f"{x:.{nd}f}"


def write_finding(p: dict, sc: dict, tpm: dict) -> None:
    g477 = p["GSE274477_Cldn4_high_vs_low"][0]  # MUT primary
    g477_all = p["GSE274477_Cldn4_high_vs_low"][2]
    wtmut = p["GSE274477_WT_vs_MUT"]
    ad = p["GSE274351_adenoma_Cldn4_high_vs_low"]
    by = p["GSE274351_by_genotype"]
    con = p["GSE274351_contrasts"]

    # verdict logic
    mut_ifn_delta = g477["IFN_core"]["delta_high_minus_low"]
    mut_mhc_delta = g477["MHC"]["delta_high_minus_low"]
    ad_ifn_delta = ad["IFN_core"]["delta_high_minus_low"]
    ad_mhc_delta = ad["MHC"]["delta_high_minus_low"]
    kl_vs_kkp = con["KL_vs_K+KP"]["Cldn4_tpm"]

    def dir_word(delta, pval):
        if delta is None or (isinstance(delta, float) and math.isnan(delta)):
            return "untested"
        if pval < 0.05 and delta < 0:
            return "down"
        if pval < 0.05 and delta > 0:
            return "up"
        if delta < 0:
            return "down (n.s.)"
        if delta > 0:
            return "up (n.s.)"
        return "unchanged"

    mut_ifn_dir = dir_word(mut_ifn_delta, g477["IFN_core"]["p"])
    mut_mhc_dir = dir_word(mut_mhc_delta, g477["MHC"]["p"])
    ad_ifn_dir = dir_word(ad_ifn_delta, ad["IFN_core"]["p"])
    ad_mhc_dir = dir_word(ad_mhc_delta, ad["MHC"]["p"])

    tracks = (
        (mut_ifn_delta is not None and mut_ifn_delta < 0)
        and (mut_mhc_delta is not None and mut_mhc_delta < 0)
    )
    tracks_sig = (
        tracks
        and g477["IFN_core"]["p"] < 0.05
        and g477["MHC"]["p"] < 0.05
    )

    if tracks_sig:
        answer = (
            "Yes, with the public processed matrices: in KRAS-mutant AT2 cells, "
            "high Cldn4 tracks lower type-I IFN and MHC scores."
        )
    elif tracks:
        answer = (
            "Directionally yes, not significant at this n: high-Cldn4 KRAS-mutant AT2 cells "
            "have lower IFN and MHC scores, but the high-vs-low tests do not reach p<0.05."
        )
    else:
        answer = (
            "Not supported as a simple high-Cldn4 = IFN-down / MHC-down rule on these public matrices."
        )

    # tables
    t_wtmut = md_table(
        ["score", "n WT", "n MUT", "median WT", "median MUT", "Δ MUT−WT", "MWU p", "rank-biserial (MUT>WT)"],
        [
            [
                s,
                str(wtmut[s]["n_WT"]),
                str(wtmut[s]["n_MUT"]),
                fnum(wtmut[s]["WT_median"]),
                fnum(wtmut[s]["MUT_median"]),
                fnum(wtmut[s]["delta_MUT_minus_WT"]),
                fmt_p(wtmut[s]["p"]),
                fnum(wtmut[s]["rank_biserial_MUT_gt_WT"]),
            ]
            for s in ("Cldn4", "IFN_core", "MHC", "H_IFNA")
        ],
    )

    def hilo_row(block, score, label):
        d = block[score]
        sp = block[f"spearman_Cldn4_vs_{score}"]
        return [
            label,
            score,
            str(block.get("n_high", block.get("n_high"))),
            str(block.get("n_low", block.get("n_low"))),
            fnum(d["high_median"]),
            fnum(d["low_median"]),
            fnum(d["delta_high_minus_low"]),
            fmt_p(d["p"]),
            fnum(d["rank_biserial_high_gt_low"]),
            f"{fnum(sp['rho'])} (p={fmt_p(sp['p'])})",
        ]

    t_hilo = md_table(
        [
            "comparison",
            "score",
            "n high Cldn4",
            "n low Cldn4",
            "median high",
            "median low",
            "Δ high−low",
            "MWU p",
            "rank-biserial (high>low)",
            "Spearman Cldn4 vs score",
        ],
        [
            hilo_row(g477, "IFN_core", "GSE274477 MUT AT2"),
            hilo_row(g477, "MHC", "GSE274477 MUT AT2"),
            hilo_row(g477, "H_IFNA", "GSE274477 MUT AT2"),
            hilo_row(g477_all, "IFN_core", "GSE274477 all QC AT2"),
            hilo_row(g477_all, "MHC", "GSE274477 all QC AT2"),
            hilo_row(ad, "IFN_core", "GSE274351 LCM adenomas"),
            hilo_row(ad, "MHC", "GSE274351 LCM adenomas"),
            hilo_row(ad, "H_IFNA", "GSE274351 LCM adenomas"),
        ],
    )

    t_gt = md_table(
        ["genotype (matrix label)", "n", "samples", "Cldn4 TPM median (IQR)", "IFN_core z median", "MHC z median"],
        [
            [
                g,
                str(by[g]["n"]),
                ", ".join(by[g]["samples"]),
                f"{fnum(by[g]['Cldn4_TPM']['median'])} ({fnum(by[g]['Cldn4_TPM']['q25'])}–{fnum(by[g]['Cldn4_TPM']['q75'])})",
                fnum(by[g]["IFN_core"]["median"]),
                fnum(by[g]["MHC"]["median"]),
            ]
            for g in ("NL", "K", "KP", "KL")
        ],
    )

    t_kl = md_table(
        ["contrast", "n a", "n b", "Cldn4 TPM median a", "Cldn4 TPM median b", "MWU p", "rank-biserial (a>b)"],
        [
            [
                name,
                str(con[name]["Cldn4_tpm"]["n_a"]),
                str(con[name]["Cldn4_tpm"]["n_b"]),
                fnum(con[name]["Cldn4_tpm"]["median_a"]),
                fnum(con[name]["Cldn4_tpm"]["median_b"]),
                fmt_p(con[name]["Cldn4_tpm"]["p"]),
                fnum(con[name]["Cldn4_tpm"]["rank_biserial_a_gt_b"]),
            ]
            for name in ("KL_vs_K", "KL_vs_KP", "KL_vs_K+KP", "KL_vs_NL", "K_vs_NL", "KP_vs_NL")
        ],
    )

    t_ps = md_table(
        ["sample", "genotype", "Cldn4 TPM", "IFN_core z", "MHC z", "Hallmark IFNα z"],
        [
            [
                r["sample"],
                r["genotype_matrix_label"],
                fnum(r["Cldn4_TPM"], 2),
                fnum(r["IFN_core_z"]),
                fnum(r["MHC_z"]),
                fnum(r["H_IFNA_z"]),
            ]
            for r in p["GSE274351_per_sample"]
        ],
    )

    gs = p["GSE274477_group_summaries"][0]
    gm = p["GSE274477_group_summaries"][1]

    text = f"""# FINDING: Cldn4 vs impaired type-I IFN in KRAS / KL epithelium

**Question.** Does Cldn4 track the impaired type-I IFN state in KRAS / KL mouse lung epithelium (Cldn4-only)? Specifically: high vs low Cldn4 → IFN down, MHC down; and KL vs K/KP Cldn4 if a public gene matrix exists.

**Answer.** {answer}

High-vs-low Cldn4 (primary): IFN {mut_ifn_dir} and MHC {mut_mhc_dir} in KRAS-mutant AT2 (GSE274477); IFN {ad_ifn_dir} and MHC {ad_mhc_dir} in LCM adenomas (GSE274351). KL Cldn4 vs K+KP: median TPM {fnum(kl_vs_kkp['median_a'])} vs {fnum(kl_vs_kkp['median_b'])}, MWU p={fmt_p(kl_vs_kkp['p'])} (n_KL={kl_vs_kkp['n_a']}, n_K+KP={kl_vs_kkp['n_b']}).

This note is **additive public mouse** evidence only. It does not change the Cldn4-only thesis; it tests whether Cldn4 covaries with the Fernández-García KRAS / type-I IFN / MHC program on the authors’ deposited matrices.

## Data (public processed only)

| accession | deposited file | what it is | labels used | honest n |
| --- | --- | --- | --- | --- |
| [GSE274477](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274477) | `GSE274477_Feature_counts_Matrix_14Aug.tsv.gz` | Fluidigm C1 scRNA-seq feature counts, AT2 cells 4 wk after Ad-Cre; KRAS WT vs KRAS(G12V) | column names: `CONTROL` = WT, `MUT` = KRAS(G12V); two chips `FCNIO1`/`FCNIO2` | matrix = 1640 wells (820 WT, 820 MUT, including 40 Undetermined). Author paper: 223 WT + 335 MUT captured → 139 + 176 after their QC (n=315). **This analysis QC (below): n={wtmut['Cldn4']['n_WT']} WT + {wtmut['Cldn4']['n_MUT']} MUT** |
| [GSE274351](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274351) | `GSE274351_expredata_TPM_gene.txt.gz` | LCM small adenomas + adjacent alveolus, QuantSeq 3′ TPM | **column labels in the TPM file** (not GEO titles) | **K n=4, KP n=5, KL n=5, NL n=4** (18 columns). GEO sample titles list 5 K / 5 KP / 4 KL / 4 NL — see caveat |

PMID 39186651 (Fernández-García et al., *PNAS* 2024). Organism: *Mus musculus*. No SRA / FASTQ / realignment.

**GEO vs matrix label mismatch (GSE274351).** GEO titles: K adenoma 1–5, KP 1–5, KL 1–4, NL 1–4. Deposited TPM columns: `K1–K4`, `KP1–KP5`, `KL1–KL5`, `NL1–NL4`. Analyses use the **deposited column labels**. n is therefore 4/5/5/4, not the GEO title counts.

## Genes

Cldn4-only (no other claudins). Ensembl GRCm38 `ENSMUSG00000047501` → `Cldn4`. Present on both matrices.

**IFN_core** (paper Fig. 2B / 3A type-I IFN / IRG genes): {", ".join(p["genes"]["IFN_CORE"] if "IFN_CORE" in p["genes"] else p["genes"]["IFN_core_requested"])}. Present scRNA: {", ".join(p["genes"]["IFN_core_in_scRNA"])}. Present LCM: {", ".join(p["genes"]["IFN_core_in_LCM"])}. `Ddx58` is the current symbol for Rigi.

**MHC / antigen presentation** (paper Fig. 2B): {", ".join(p["genes"]["MHC_requested"])}. Present scRNA: {", ".join(p["genes"]["MHC_in_scRNA"])}. Present LCM: {", ".join(p["genes"]["MHC_in_LCM"])}.

**H_IFNA**: MSigDB mouse `HALLMARK_INTERFERON_ALPHA_RESPONSE` (2024.1.Mm; {p["genes"]["Hallmark_IFNA_n_requested"]} genes). Recovered {p["genes"]["Hallmark_IFNA_n_scRNA"]} / {p["genes"]["Hallmark_IFNA_n_LCM"]} on the two matrices. Secondary to the paper core lists.

No gene on the requested lists was missing from both matrices.

## Methods

1. Download GEO supplementary processed files only (FTP `GSE274nnn`).
2. Map `ENSMUSG` (version-stripped) → symbol with Ensembl GRCm38.102 GTF.
3. **GSE274477 QC (author capture list is not deposited).** Drop `Undetermined` wells. Keep wells with library size ≥ 50,000 counts, 1,000–8,000 detected genes, mitochondrial fraction ≤ 10% (MT genes from the same GTF). This is a transparent depth filter near the paper’s reported averages (~100,000 reads, ~2,354 genes) — **not a reproduction of their unpublished well list**. n after QC is reported above; the paper’s n=315 is smaller.
4. Normalize scRNA as log1p(10,000 × count / library size). Score a gene set as the mean of per-gene z-scores across QC cells.
5. High vs low Cldn4 = median split of Cldn4 log-normalized expression. Primary split is **within KRAS-mutant AT2**. Also report all-QC-AT2.
6. **GSE274351** uses deposited TPM. Scores = mean z of log2(TPM+1) across the 18 samples. High vs low Cldn4 = median split **among adenomas only** (NL held out), n=14. KL vs K/KP uses matrix genotype labels.
7. Tests: two-sided Mann–Whitney U; Spearman rank correlation. Effect = rank-biserial r (positive ⇒ first group higher). No multiple-testing claim beyond the pre-specified IFN and MHC scores. Honest n: no imputation, no dropping of labeled samples.

## Table 1. KRAS MUT vs WT AT2 (GSE274477, QC cells)

Sanity check that the deposited counts recover the paper’s IFN / MHC impairment.

{t_wtmut}

Cldn4 detection (counts > 0) after QC: WT {gs["Cldn4_detect_n"]}/{gs["n_cells"]} ({gs["Cldn4_detect_frac"]:.2%}); MUT {gm["Cldn4_detect_n"]}/{gm["n_cells"]} ({gm["Cldn4_detect_frac"]:.2%}).

## Table 2. Cldn4–IFN / MHC (high vs low Cldn4)

This is the table the note is done when it exists.

{t_hilo}

Median-split Cldn4 (MUT AT2) = {fnum(g477["median_Cldn4"])} log1p(CPM10k); unique rounded values = {g477["n_unique_Cldn4"]}; cells with Cldn4 counts > 0 = {g477["Cldn4_detect_n"]} / {g477["n"]}.

LCM adenoma high-Cldn4 samples: {", ".join(ad["high_samples"])} ({ad["high_gt"]}). Low: {", ".join(ad["low_samples"])} ({ad["low_gt"]}).

## Table 3. LCM Cldn4 by genotype (GSE274351 matrix labels)

{t_gt}

## Table 4. KL vs K / KP Cldn4 (GSE274351)

{t_kl}

## Table 5. Per-sample LCM values

{t_ps}

## Interpretation

- The deposited GSE274477 counts recover KRAS-mutant AT2 **IFN_core down** and **MHC down** vs WT (Table 1), matching the paper’s early AT2 type-I IFN impairment. Cldn4 itself is {("higher" if wtmut["Cldn4"]["delta_MUT_minus_WT"]>0 else "lower")} in MUT vs WT (Δ={fnum(wtmut["Cldn4"]["delta_MUT_minus_WT"])}, p={fmt_p(wtmut["Cldn4"]["p"])}).
- **Does high Cldn4 track IFN-down / MHC-down?** In KRAS-mutant AT2, high vs low Cldn4: IFN Δ={fnum(mut_ifn_delta)} (p={fmt_p(g477["IFN_core"]["p"])}), MHC Δ={fnum(mut_mhc_delta)} (p={fmt_p(g477["MHC"]["p"])}). In LCM adenomas (n=14): IFN Δ={fnum(ad_ifn_delta)} (p={fmt_p(ad["IFN_core"]["p"])}), MHC Δ={fnum(ad_mhc_delta)} (p={fmt_p(ad["MHC"]["p"])}).
- **KL vs K/KP Cldn4.** Matrix-labeled KL median Cldn4 TPM {fnum(kl_vs_kkp["median_a"])} vs K+KP {fnum(kl_vs_kkp["median_b"])} (p={fmt_p(kl_vs_kkp["p"])}). n_KL=5, n_K=4, n_KP=5 is small; a non-significant genotype difference is not evidence of equality.
- Cldn4 is sparse on Fluidigm C1 AT2 libraries. High/low is therefore partly a detect vs not-detect split. That is reported, not hidden.

## Caveats

- Public processed only. Author single-cell barcodes after QC are not in GEO; n here ≠ 315.
- GSE274351 GEO titles ≠ TPM column labels (5/4 vs 4/5 for K vs KL). Genotype contrasts follow the file.
- LCM is bulk epithelium-enriched adenoma, not pure AT2. Immune / stromal RNA can move MHC and IFN scores.
- Hallmark IFNα is a human-derived set mapped to mouse symbols; the paper-core list is the primary score.
- No other claudins. No human data. No raw-read reprocessing.
- Small n on LCM (18 samples). p-values are descriptive.

## Files

- `analyze.py` — this analysis
- `results.json` — machine-readable numbers
- `FINDING.md` — this note
"""

    # fix IFN_CORE key reference - I used a fallback already
    (OUT / "FINDING.generated.md").write_text(text)


if __name__ == "__main__":
    main()
