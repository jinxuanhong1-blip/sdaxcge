#!/usr/bin/env python3
"""Score open CLDN4 knockdown/knockout matrices for NHEJ down, STING up, IFN up, APM up.

ArrayExpress/BioStudies search date: 2026-09-21.
Processed matrices for the two ArrayExpress RNA experiments are the GEO deposits
of the same accessions (E-GEOD-50927 = GSE50927, E-GEOD-22493 = GSE22493).
BioStudies file HTTP returned 404 and ftp.ebi.ac.uk TLS is blocked in this
environment. GSE207704 is the human CLDN4 CRISPR matrix linked from BioStudies
paper S-EPMC10105442; it has no E-GEOD record.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import urllib.request
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "ae_cldn4_kd_wave"
CACHE = Path(os.environ.get("AE_CACHE", "/tmp/ae_dl/geo"))
OUT.mkdir(parents=True, exist_ok=True)
CACHE.mkdir(parents=True, exist_ok=True)

GEO = {
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "GSE50927_VILIwtGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtGenes.csv.gz",
    "GSE50927_VILIwtkoloGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
    "GSE50927_VILIwtkohiGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
    "GSE22493_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/matrix/GSE22493_series_matrix.txt.gz",
    "GPL10555_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10555/soft/GPL10555_family.soft.gz",
    "GSE207704_CLDN4_RNAseq.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz",
}

# Alias groups. First symbol present in the matrix is used; aliases are not averaged.
NHEJ = [
    ("XRCC6",),
    ("XRCC5",),
    ("PRKDC",),
    ("LIG4",),
    ("XRCC4",),
    ("NHEJ1",),
    ("PAXX",),
    ("DCLRE1C",),
    ("POLL",),
    ("POLM",),
    ("PNKP",),
    ("APTX",),
]
NHEJ_ANCHOR = [("TP53BP1", "TRP53BP1")]
STING = [
    ("CGAS", "MB21D1"),
    ("STING1", "TMEM173"),
    ("TBK1",),
    ("IKBKE",),
    ("IRF3",),
    ("IRF7",),
    ("IFNB1",),
    ("CXCL10",),
    ("CCL5",),
]
# Sensitivity split: sensors/adapters only. Chemokine outputs stay in the primary STING set.
STING_SENSOR = [
    ("CGAS", "MB21D1"),
    ("STING1", "TMEM173"),
    ("TBK1",),
    ("IKBKE",),
    ("IRF3",),
    ("IRF7",),
]
IFN = [
    ("ISG15",),
    ("MX1",),
    ("MX2",),
    ("OAS1",),
    ("OAS2",),
    ("OAS3",),
    ("OASL",),
    ("IFIT1",),
    ("IFIT2",),
    ("IFIT3",),
    ("IFI6",),
    ("IFI27",),
    ("IFI44",),
    ("IFI44L",),
    ("RSAD2",),
    ("BST2",),
    ("STAT1",),
    ("STAT2",),
    ("IRF1",),
    ("IRF9",),
    ("IFITM1",),
    ("IFITM3",),
    ("USP18",),
    ("HERC5",),
    ("IFIH1",),
    ("DDX58",),
    ("CXCL9",),
    ("CXCL11",),
    ("ISG20",),
]
APM_HUMAN = [
    ("HLA-A",),
    ("HLA-B",),
    ("HLA-C",),
    ("B2M",),
    ("TAP1",),
    ("TAP2",),
    ("TAPBP",),
    ("PSMB8",),
    ("PSMB9",),
    ("PSMB10",),
    ("NLRC5",),
    ("ERAP1",),
    ("ERAP2",),
    ("CALR",),
    ("PDIA3",),
    ("HLA-E",),
]
APM_MOUSE = [
    ("H2-K1",),
    ("H2-D1",),
    ("B2M",),
    ("TAP1",),
    ("TAP2",),
    ("TAPBP",),
    ("PSMB8",),
    ("PSMB9",),
    ("PSMB10",),
    ("NLRC5",),
    ("ERAP1",),
    ("CALR",),
    ("PDIA3",),
    ("H2-T23",),
    ("H2-M3",),
]

AXES = {
    "NHEJ": {"genes": NHEJ, "expect": "down"},
    "STING": {"genes": STING, "expect": "up"},
    "IFN": {"genes": IFN, "expect": "up"},
    "APM": {"genes": None, "expect": "up"},  # species-specific
}


def fetch(name: str) -> Path:
    dest = CACHE / name
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    url = GEO[name]
    req = urllib.request.Request(url, headers={"User-Agent": "research"})
    with urllib.request.urlopen(req, timeout=120) as r, dest.open("wb") as out:
        out.write(r.read())
    return dest


def mouse_symbol(human: str) -> str:
    special = {
        "TP53BP1": "Trp53bp1",
        "TRP53BP1": "Trp53bp1",
        "HLA-A": "H2-K1",
        "HLA-B": "H2-D1",
        "HLA-E": "H2-T23",
        "H2-K1": "H2-K1",
        "H2-D1": "H2-D1",
        "H2-T23": "H2-T23",
        "H2-M3": "H2-M3",
    }
    if human in special:
        return special[human]
    return human[0] + human[1:].lower()


def species_aliases(groups, species: str):
    out = []
    for group in groups:
        if species == "mouse":
            mapped = []
            for g in group:
                mapped.append(mouse_symbol(g))
            # unique, keep order
            seen = []
            for g in mapped:
                if g not in seen:
                    seen.append(g)
            out.append(tuple(seen))
        else:
            out.append(group)
    return out


def bh(pvals):
    p = np.asarray(pvals, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return q
    pv = p[ok]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]
    adj = ranked * n / (np.arange(1, n + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    q_ok = np.empty(n)
    q_ok[order] = adj
    q[ok] = q_ok
    return q


def index_rows(rows, symbol_key):
    """Keep the highest-logCPM row when a symbol is duplicated."""
    best = {}
    for row in rows:
        sym = (row.get(symbol_key) or "").strip()
        if not sym or sym == "NA":
            continue
        cpm = float(row["logCPM"]) if row.get("logCPM") not in (None, "", "NA") else -1e9
        prev = best.get(sym.lower())
        if prev is None or cpm > prev[0]:
            best[sym.lower()] = (
                cpm,
                {
                    "symbol": sym,
                    "log2fc": float(row["logFC"]),
                    "pval": float(row["PValue"]),
                    "fdr": float(row["FDR"]),
                },
            )
    return {k: v[1] for k, v in best.items()}


def load_edger(path: Path):
    with gzip.open(path, "rt", newline="") as f:
        rows = list(csv.DictReader(f))
    key = "Marker.Symbol" if "Marker.Symbol" in rows[0] else "GeneSymbol"
    table = index_rows(rows, key)
    return table, len(rows)


def load_gpl(path: Path):
    id_to_gene = {}
    with gzip.open(path, "rt", errors="replace") as f:
        in_table = False
        header = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(f).rstrip("\n").split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.rstrip("\n").split("\t")
            rec = dict(zip(header, parts))
            gene = (rec.get("ORF") or "").strip()
            if gene:
                id_to_gene[rec["ID"]] = gene
    return id_to_gene


def load_gse22493(matrix_path: Path, id_to_gene: dict):
    """VALUE is log2(Cy5/Cy3) = log2(CLDN4-siRNA / CLDN4-overexpression)."""
    probes = []
    with gzip.open(matrix_path, "rt", errors="replace") as f:
        in_table = False
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                header = next(f).rstrip("\n").split("\t")
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(parts) < 4 or parts[0] == "ID_REF":
                continue
            vals = []
            for x in parts[1:4]:
                if x in ("", "null", "NA"):
                    vals.append(np.nan)
                else:
                    vals.append(float(x))
            probes.append((parts[0], np.array(vals, dtype=float)))

    by_gene = defaultdict(list)
    for pid, vals in probes:
        gene = id_to_gene.get(pid)
        if gene:
            by_gene[gene.upper()].append(vals)

    table = {}
    gene_ps = []
    genes = []
    for gene, arrs in by_gene.items():
        stacked = np.vstack(arrs)  # probes x 3
        stacked = stacked[np.isfinite(stacked).any(axis=1)]
        if stacked.size == 0:
            continue
        with np.errstate(all="ignore"):
            per_array = np.nanmedian(stacked, axis=0)
        if not np.isfinite(per_array).any():
            continue
        mean = float(np.nanmean(per_array))
        finite = per_array[np.isfinite(per_array)]
        if len(finite) >= 2 and np.nanstd(finite, ddof=1) > 0:
            p = float(stats.ttest_1samp(finite, 0.0).pvalue)
        elif len(finite) >= 1 and np.allclose(finite, 0):
            p = 1.0
        else:
            p = np.nan
        table[gene.lower()] = {
            "symbol": gene,
            "log2fc": mean,
            "pval": p,
            "fdr": np.nan,
            "n_probes": len(arrs),
            "replicates": per_array.tolist(),
        }
        genes.append(gene.lower())
        gene_ps.append(p)
    q = bh(gene_ps)
    for g, qq in zip(genes, q):
        table[g]["fdr"] = float(qq) if np.isfinite(qq) else np.nan
    return table, len(by_gene)


def load_gse207704(path: Path):
    """Collapsed FPKM. log2FC uses pseudocount 0.1. No replicate-level p."""
    lines = {"MCF7": {}, "T47D": {}}
    with gzip.open(path, "rt", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        # columns: gene_short_name, MCF7_CLDN4KO_FPKM, MCF7_WT_FPKM, T47D_CLDN4KO_FPKM, T47D_WT_FPKM
        idx = {name.split(" ")[0]: i for i, name in enumerate(header)}
        name_i = header.index("gene_short_name")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= name_i:
                continue
            gene = parts[name_i].strip().strip('"')
            if not gene or gene == "gene_short_name":
                continue
            try:
                mko = float(parts[idx["MCF7_CLDN4KO_FPKM"]])
                mwt = float(parts[idx["MCF7_WT_FPKM"]])
                tko = float(parts[idx["T47D_CLDN4KO_FPKM"]])
                twt = float(parts[idx["T47D_WT_FPKM"]])
            except (KeyError, ValueError):
                continue
            key = gene.upper()
            expr = mko + mwt + tko + twt
            prev = lines["MCF7"].get(key)
            if prev is not None and prev["expr"] >= expr:
                continue
            lines["MCF7"][key] = {
                "symbol": gene,
                "log2fc": math.log2((mko + 0.1) / (mwt + 0.1)),
                "pval": np.nan,
                "fdr": np.nan,
                "expr": expr,
                "ko": mko,
                "wt": mwt,
            }
            lines["T47D"][key] = {
                "symbol": gene,
                "log2fc": math.log2((tko + 0.1) / (twt + 0.1)),
                "pval": np.nan,
                "fdr": np.nan,
                "expr": expr,
                "ko": tko,
                "wt": twt,
            }
    mean = {}
    for key, m in lines["MCF7"].items():
        t = lines["T47D"][key]
        mean[key.lower()] = {
            "symbol": m["symbol"],
            "log2fc": (m["log2fc"] + t["log2fc"]) / 2.0,
            "pval": np.nan,
            "fdr": np.nan,
            "mcf7": m["log2fc"],
            "t47d": t["log2fc"],
        }
    mcf7 = {k.lower(): {kk: vv for kk, vv in v.items() if kk != "expr"} for k, v in lines["MCF7"].items()}
    t47d = {k.lower(): {kk: vv for kk, vv in v.items() if kk != "expr"} for k, v in lines["T47D"].items()}
    return mcf7, t47d, mean


def lookup(table, aliases):
    for a in aliases:
        hit = table.get(a.lower())
        if hit is not None and np.isfinite(hit["log2fc"]):
            return a, hit
    return None, None


def score_axis(table, groups, expect: str):
    used = []
    missing = []
    for group in groups:
        alias, hit = lookup(table, group)
        if hit is None:
            missing.append("|".join(group))
            continue
        used.append(
            {
                "gene": hit["symbol"],
                "alias_queried": alias,
                "log2fc": hit["log2fc"],
                "pval": hit.get("pval", np.nan),
                "fdr": hit.get("fdr", np.nan),
                "n_probes": hit.get("n_probes", 1),
            }
        )
    vals = np.array([u["log2fc"] for u in used], dtype=float)
    n = len(vals)
    out = {
        "n_measured": n,
        "n_missing": len(missing),
        "missing": ";".join(missing),
        "median_log2fc": float(np.median(vals)) if n else np.nan,
        "mean_log2fc": float(np.mean(vals)) if n else np.nan,
        "n_up": int(np.sum(vals > 0)) if n else 0,
        "n_down": int(np.sum(vals < 0)) if n else 0,
        "n_zero": int(np.sum(vals == 0)) if n else 0,
        "genes": used,
    }
    if n >= 6:
        diff = vals[vals != 0]
        alt = "less" if expect == "down" else "greater"
        alt_opp = "greater" if expect == "down" else "less"
        try:
            if len(diff) >= 6:
                w = stats.wilcoxon(diff, alternative=alt, zero_method="wilcox")
                out["wilcoxon_p_expected"] = float(w.pvalue)
                wopp = stats.wilcoxon(diff, alternative=alt_opp, zero_method="wilcox")
                out["wilcoxon_p_opposite"] = float(wopp.pvalue)
            else:
                out["wilcoxon_p_expected"] = np.nan
                out["wilcoxon_p_opposite"] = np.nan
        except ValueError:
            out["wilcoxon_p_expected"] = np.nan
            out["wilcoxon_p_opposite"] = np.nan
        set_syms = {u["gene"].lower() for u in used}
        bg = np.array(
            [v["log2fc"] for k, v in table.items() if k not in set_syms and np.isfinite(v["log2fc"])],
            dtype=float,
        )
        if len(bg) > 20 and n >= 6:
            mw = stats.mannwhitneyu(vals, bg, alternative="two-sided")
            out["mw_p_vs_background"] = float(mw.pvalue)
            out["background_median_log2fc"] = float(np.median(bg))
        else:
            out["mw_p_vs_background"] = np.nan
            out["background_median_log2fc"] = np.nan
    else:
        out["wilcoxon_p_expected"] = np.nan
        out["wilcoxon_p_opposite"] = np.nan
        out["mw_p_vs_background"] = np.nan
        out["background_median_log2fc"] = np.nan
    fdrs = np.array([u["fdr"] for u in used], dtype=float)
    lfc = vals
    if expect == "down":
        hit_mask = (lfc < 0) & np.isfinite(fdrs) & (fdrs < 0.05)
        frac = (out["n_down"] / n) if n else np.nan
    else:
        hit_mask = (lfc > 0) & np.isfinite(fdrs) & (fdrs < 0.05)
        frac = (out["n_up"] / n) if n else np.nan
    out["frac_expected"] = float(frac) if n else np.nan
    out["frac_opposite"] = float((out["n_up"] if expect == "down" else out["n_down"]) / n) if n else np.nan
    out["n_fdr05_expected"] = int(np.sum(hit_mask)) if n else 0
    out["call"] = call_axis(
        out["median_log2fc"],
        n,
        out["frac_expected"],
        out["frac_opposite"],
        out["wilcoxon_p_expected"],
        out.get("wilcoxon_p_opposite", np.nan),
        expect,
    )
    return out


def call_axis(median, n, frac_exp, frac_opp, p_exp, p_opp, expect) -> str:
    """Gene-set call. Not a sample-level replicate test.

    supports / opposite require |median log2FC| >= 0.10 and a one-sided
    Wilcoxon p < 0.05. directional / opposite_directional are the same
    magnitude and a 2/3 fraction, without that p-value.
    """
    if n < 6 or not np.isfinite(median):
        return "insufficient"
    mag = abs(median) >= 0.10
    med_expected = median < 0 if expect == "down" else median > 0
    if med_expected:
        if mag and np.isfinite(p_exp) and p_exp < 0.05:
            return "supports"
        if mag and np.isfinite(frac_exp) and frac_exp >= 0.67:
            return "directional"
        return "null"
    if mag and np.isfinite(p_opp) and p_opp < 0.05:
        return "opposite"
    if mag and np.isfinite(frac_opp) and frac_opp >= 0.67:
        return "opposite_directional"
    return "null"


def fmt(x, nd=4):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return ""
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def main():
    for name in GEO:
        fetch(name)

    uninj, n_un = load_edger(fetch("GSE50927_Cldn4lungWTvsKOgenes.csv.gz"))
    vlow, _ = load_edger(fetch("GSE50927_VILIwtkoloGenes.csv.gz"))
    vhigh, _ = load_edger(fetch("GSE50927_VILIwtkohiGenes.csv.gz"))
    vwt, _ = load_edger(fetch("GSE50927_VILIwtGenes.csv.gz"))
    cldn4_base = uninj["cldn4"]["log2fc"]
    if not (cldn4_base < -5):
        raise SystemExit(f"QC fail: uninjured Cldn4 logFC {cldn4_base} is not a KO-down contrast")

    id_to_gene = load_gpl(fetch("GPL10555_family.soft.gz"))
    g224, n224 = load_gse22493(fetch("GSE22493_series_matrix.txt.gz"), id_to_gene)
    mcf7, t47d, mean207 = load_gse207704(fetch("GSE207704_CLDN4_RNAseq.txt.gz"))

    contrasts = [
        {
            "contrast_id": "E-GEOD-50927_uninjured_KOvsWT",
            "accession": "E-GEOD-50927",
            "geo": "GSE50927",
            "organism": "mouse",
            "species_key": "mouse",
            "perturbation": "germline Cldn4 KO vs WT, whole lung, no VILI",
            "assay": "RNA-seq author edgeR",
            "n_note": "1 vs 1; deposited edgeR logFC/FDR, not a new replicate test",
            "table": uninj,
            "in_arrayexpress": "yes",
            "primary": "yes",
        },
        {
            "contrast_id": "E-GEOD-50927_VILIlow_KOvsWT",
            "accession": "E-GEOD-50927",
            "geo": "GSE50927",
            "organism": "mouse",
            "species_key": "mouse",
            "perturbation": "Cldn4 KO VILI-low vs WT VILI",
            "assay": "RNA-seq author edgeR",
            "n_note": "1 vs 1; Cldn4 logFC strongly negative so treated as KO vs WT",
            "table": vlow,
            "in_arrayexpress": "yes",
            "primary": "no",
        },
        {
            "contrast_id": "E-GEOD-50927_VILIhigh_KOvsWT",
            "accession": "E-GEOD-50927",
            "geo": "GSE50927",
            "organism": "mouse",
            "species_key": "mouse",
            "perturbation": "Cldn4 KO VILI-high vs WT VILI",
            "assay": "RNA-seq author edgeR",
            "n_note": "1 vs 1; Cldn4 logFC strongly negative so treated as KO vs WT",
            "table": vhigh,
            "in_arrayexpress": "yes",
            "primary": "no",
        },
        {
            "contrast_id": "E-GEOD-22493_siRNA_vs_OE",
            "accession": "E-GEOD-22493",
            "geo": "GSE22493",
            "organism": "human",
            "species_key": "human",
            "perturbation": "SKOV-3-IP-Luc CLDN4 siRNA (Cy5) vs CLDN4 overexpression control (Cy3)",
            "assay": "two-color array, 3 arrays, log2(KD/OE)",
            "n_note": "n=3 arrays; control is overexpression, not scramble/WT",
            "table": g224,
            "in_arrayexpress": "yes",
            "primary": "yes",
        },
        {
            "contrast_id": "GSE207704_MCF7_KOvsWT",
            "accession": "GSE207704",
            "geo": "GSE207704",
            "organism": "human",
            "species_key": "human",
            "perturbation": "MCF7 CLDN4 CRISPR KO vs WT, collapsed FPKM",
            "assay": "RNA-seq FPKM, replicates already pooled",
            "n_note": "no replicate p; BioStudies paper S-EPMC10105442, no E-GEOD",
            "table": mcf7,
            "in_arrayexpress": "no",
            "primary": "no",
        },
        {
            "contrast_id": "GSE207704_T47D_KOvsWT",
            "accession": "GSE207704",
            "geo": "GSE207704",
            "organism": "human",
            "species_key": "human",
            "perturbation": "T47D CLDN4 CRISPR KO vs WT, collapsed FPKM",
            "assay": "RNA-seq FPKM, replicates already pooled",
            "n_note": "no replicate p; BioStudies paper S-EPMC10105442, no E-GEOD",
            "table": t47d,
            "in_arrayexpress": "no",
            "primary": "yes",
        },
        {
            "contrast_id": "GSE207704_mean_MCF7_T47D",
            "accession": "GSE207704",
            "geo": "GSE207704",
            "organism": "human",
            "species_key": "human",
            "perturbation": "mean log2FC of MCF7 and T47D CLDN4 CRISPR KO vs WT",
            "assay": "RNA-seq FPKM, two lines",
            "n_note": "lines are not replicates of one model; mean is a summary only",
            "table": mean207,
            "in_arrayexpress": "no",
            "primary": "yes",
        },
    ]

    score_rows = []
    gene_rows = []
    summaries = []
    for c in contrasts:
        apm = APM_MOUSE if c["species_key"] == "mouse" else APM_HUMAN
        axis_defs = {
            "NHEJ": (species_aliases(NHEJ, c["species_key"]), "down"),
            "STING": (species_aliases(STING, c["species_key"]), "up"),
            "STING_sensor": (species_aliases(STING_SENSOR, c["species_key"]), "up"),
            "IFN": (species_aliases(IFN, c["species_key"]), "up"),
            "APM": (species_aliases(apm, c["species_key"]), "up"),
        }
        calls = {}
        for axis, (groups, expect) in axis_defs.items():
            sc = score_axis(c["table"], groups, expect)
            calls[axis] = sc["call"]
            score_rows.append(
                {
                    "contrast_id": c["contrast_id"],
                    "accession": c["accession"],
                    "geo": c["geo"],
                    "in_arrayexpress": c["in_arrayexpress"],
                    "organism": c["organism"],
                    "perturbation": c["perturbation"],
                    "primary": c["primary"],
                    "axis": axis,
                    "expected": expect,
                    "call": sc["call"],
                    "n_measured": sc["n_measured"],
                    "n_missing": sc["n_missing"],
                    "n_up": sc["n_up"],
                    "n_down": sc["n_down"],
                    "frac_expected": sc["frac_expected"],
                    "median_log2fc": sc["median_log2fc"],
                    "mean_log2fc": sc["mean_log2fc"],
                    "wilcoxon_p_expected_direction": sc["wilcoxon_p_expected"],
                    "wilcoxon_p_opposite_direction": sc.get("wilcoxon_p_opposite", np.nan),
                    "mw_p_vs_background": sc["mw_p_vs_background"],
                    "background_median_log2fc": sc["background_median_log2fc"],
                    "n_fdr05_expected_direction": sc["n_fdr05_expected"],
                    "missing_genes": sc["missing"],
                }
            )
            for g in sc["genes"]:
                matches = (g["log2fc"] < 0) if expect == "down" else (g["log2fc"] > 0)
                gene_rows.append(
                    {
                        "contrast_id": c["contrast_id"],
                        "accession": c["accession"],
                        "axis": axis,
                        "expected": expect,
                        "gene": g["gene"],
                        "log2fc_KOorKD_vs_control": g["log2fc"],
                        "p": g["pval"],
                        "fdr": g["fdr"],
                        "n_probes": g["n_probes"],
                        "matches_expected": "yes" if matches else "no",
                    }
                )
        # anchor gene, not part of the four-axis call
        anchor_groups = species_aliases(NHEJ_ANCHOR, c["species_key"])
        alias, hit = lookup(c["table"], anchor_groups[0])
        if hit:
            gene_rows.append(
                {
                    "contrast_id": c["contrast_id"],
                    "accession": c["accession"],
                    "axis": "NHEJ_anchor_TP53BP1",
                    "expected": "down",
                    "gene": hit["symbol"],
                    "log2fc_KOorKD_vs_control": hit["log2fc"],
                    "p": hit.get("pval", np.nan),
                    "fdr": hit.get("fdr", np.nan),
                    "n_probes": hit.get("n_probes", 1),
                    "matches_expected": "yes" if hit["log2fc"] < 0 else "no",
                }
            )
        pattern = " ".join(f"{ax}:{calls[ax]}" for ax in ("NHEJ", "STING", "IFN", "APM"))
        n_support = sum(calls[ax] == "supports" for ax in ("NHEJ", "STING", "IFN", "APM"))
        summaries.append({**c, "pattern": pattern, "n_axes_support": n_support, "calls": calls})

    # perturbation QC
    qc = []

    def add_qc(contrast, gene, table):
        hit = table.get(gene.lower())
        qc.append(
            {
                "contrast_id": contrast,
                "gene": gene,
                "log2fc": None if hit is None else hit["log2fc"],
                "p": None if hit is None else hit.get("pval"),
                "fdr": None if hit is None else hit.get("fdr"),
                "n_probes": None if hit is None else hit.get("n_probes", 1),
            }
        )

    for gene in ("Cldn4", "Cldn3", "Cldn7"):
        add_qc("E-GEOD-50927_uninjured_KOvsWT", gene, uninj)
        add_qc("E-GEOD-50927_VILIlow_KOvsWT", gene, vlow)
        add_qc("E-GEOD-50927_VILIhigh_KOvsWT", gene, vhigh)
        add_qc("E-GEOD-50927_VILI_in_WT_NOT_a_KO_contrast", gene, vwt)
    add_qc("E-GEOD-22493_siRNA_vs_OE", "CLDN4", g224)
    add_qc("E-GEOD-22493_siRNA_vs_OE", "TACSTD2", g224)
    add_qc("GSE207704_MCF7_KOvsWT", "CLDN4", mcf7)
    add_qc("GSE207704_T47D_KOvsWT", "CLDN4", t47d)
    add_qc("GSE207704_mean_MCF7_T47D", "CLDN4", mean207)

    inventory = [
        {
            "accession": "E-GEOD-50927",
            "repository": "ArrayExpress (GEO mirror GSE50927)",
            "organism": "Mus musculus",
            "perturbation": "germline Cldn4 knockout",
            "assay": "RNA-seq, whole lung",
            "qualifies": "yes",
            "scored": "yes",
            "why": "Direct Cldn4 KO transcriptome. Uninjured KO vs WT is the primary contrast. VILI-low and VILI-high are secondary KO vs WT (Cldn4 logFC -12.7 and -10.2). VILIwtGenes is injury in WT (Cldn4 logFC +3.96) and is not scored.",
        },
        {
            "accession": "E-GEOD-22493",
            "repository": "ArrayExpress (GEO mirror GSE22493)",
            "organism": "Homo sapiens",
            "perturbation": "CLDN4 siRNA vs CLDN4 overexpression",
            "assay": "two-color microarray, n=3",
            "qualifies": "yes_with_caveat",
            "scored": "yes",
            "why": "Direct CLDN4 siRNA, but the control channel is CLDN4 overexpression, not scramble or parental WT. Ovarian SKOV-3-IP-Luc, not lung.",
        },
        {
            "accession": "GSE207704",
            "repository": "GEO only; paper record S-EPMC10105442 in BioStudies has a docx supplement, no expression matrix and no E-GEOD accession",
            "organism": "Homo sapiens",
            "perturbation": "CLDN4 CRISPR knockout, MCF7 and T47D",
            "assay": "RNA-seq collapsed FPKM",
            "qualifies": "linked_not_arrayexpress",
            "scored": "yes",
            "why": "Only public human CLDN4 CRISPR transcriptome. Scored from the GEO FPKM table so the ArrayExpress-only result is not mistaken for the human cancer-cell test. Replicates are already pooled.",
        },
        {
            "accession": "E-GEOD-60885",
            "repository": "ArrayExpress",
            "organism": "Homo sapiens",
            "perturbation": "none",
            "assay": "DNA methylation",
            "qualifies": "no",
            "scored": "no",
            "why": "Methylation study that names Claudin-4 as a trophoblast invasion gene. Not a CLDN4 knockdown or knockout transcriptome.",
        },
        {
            "accession": "E-GEOD-84742",
            "repository": "ArrayExpress",
            "organism": "Mus musculus",
            "perturbation": "Hopx/Klf4/Tcf7l2/Hnf4a axis, not Cldn4",
            "assay": "transcription profiling",
            "qualifies": "no",
            "scored": "no",
            "why": "Claudin expression gradients during colonic differentiation. Not a Cldn4 genetic perturbation.",
        },
        {
            "accession": "E-GEOD-48443",
            "repository": "ArrayExpress",
            "organism": "Mus musculus",
            "perturbation": "Cldn18 knockout",
            "assay": "transcription profiling",
            "qualifies": "no",
            "scored": "no",
            "why": "Claudin-18 knockout lung. Not CLDN4.",
        },
        {
            "accession": "E-GEOD-26055",
            "repository": "ArrayExpress",
            "organism": "Homo sapiens",
            "perturbation": "CLDN7 overexpression",
            "assay": "transcription profiling",
            "qualifies": "no",
            "scored": "no",
            "why": "Claudin-7, not CLDN4 knockdown.",
        },
        {
            "accession": "S-BSST1967",
            "repository": "BioStudies",
            "organism": "Homo sapiens",
            "perturbation": "CPP-S4 peptide against CLDN4 palmitoylation; CLDN4 knockdown is described but not the deposited RNA-seq contrast",
            "assay": "RNA-seq DESeq2 file read_count.cnt_vs_ctp plus proteomics ExpAll.xlsx",
            "qualifies": "no",
            "scored": "no",
            "why": "Xu et al. Cell Reports Medicine 2025 (PMC12281411). Deposited bulk transcriptome is MHCC97H CPP-S4 vs mock, and proteomics are lenvatinib-vs-mock and palmitoyl-proteome. Raw RNA-seq is CNGB CNP0006650. Not a genetic KD/KO/CRISPR/siRNA matrix. ftp.ebi.ac.uk TLS failed here, so the xls/xlsx were not opened; the contrast is taken from the paper.",
        },
        {
            "accession": "S-EPMC8988515",
            "repository": "BioStudies (PMC author manuscript)",
            "organism": "Homo sapiens",
            "perturbation": "CLDN4 shRNA functional assays, not a deposited KD transcriptome",
            "assay": "NHEJ reporter and supplements",
            "qualifies": "no",
            "scored": "no",
            "why": "Yamamoto et al. Mol Cancer Ther 2022. Table S1 is TCGA CLDN4-high vs low (1,582 transcripts), not shCLDN4 RNA-seq. Supplement xlsx download was blocked by PMC recaptcha after the HTML identified the table. No GEO accession.",
        },
        {
            "accession": "S-EPMC12603150",
            "repository": "BioStudies (Scientific Reports 2025)",
            "organism": "Homo sapiens",
            "perturbation": "CLDN4 CRISPRi in OVCA429 and OVCAR3",
            "assay": "functional IFN/STING assays; no public matrix",
            "qualifies": "no",
            "scored": "no",
            "why": "Data availability statement: raw data available on request. No GEO, ArrayExpress, or PRIDE accession. Not scored.",
        },
        {
            "accession": "E-MTAB / E-MEXP / E-PROT / PXD",
            "repository": "ArrayExpress and PRIDE via BioStudies",
            "organism": "",
            "perturbation": "CLDN4 or Cldn4 or claudin-4",
            "assay": "any",
            "qualifies": "no",
            "scored": "no",
            "why": "2026-09-21 queries returned 0 E-MTAB, 0 E-MEXP, 0 E-PROT, and 0 PXD studies. PRIDE keyword search for CLDN4, claudin-4, and Cldn4 returned 0 projects.",
        },
    ]

    def write_tsv(path, rows):
        if not rows:
            return
        fields = list(rows[0].keys())
        with path.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow({k: fmt(row[k]) if not isinstance(row[k], str) else row[k] for k in fields})

    write_tsv(OUT / "signature_scores.tsv", score_rows)
    write_tsv(OUT / "de_genes.tsv", gene_rows)
    write_tsv(OUT / "qc_perturbation.tsv", qc)
    write_tsv(OUT / "inventory.tsv", inventory)

    # figure: median log2FC by axis for primary contrasts
    primary_ids = [
        "E-GEOD-50927_uninjured_KOvsWT",
        "E-GEOD-22493_siRNA_vs_OE",
        "GSE207704_T47D_KOvsWT",
        "GSE207704_MCF7_KOvsWT",
    ]
    short = {
        "E-GEOD-50927_uninjured_KOvsWT": "Cldn4 KO lung",
        "E-GEOD-22493_siRNA_vs_OE": "siRNA vs OE",
        "GSE207704_T47D_KOvsWT": "T47D CRISPR",
        "GSE207704_MCF7_KOvsWT": "MCF7 CRISPR",
    }
    axes = ["NHEJ", "STING", "IFN", "APM"]
    fig, ax = plt.subplots(figsize=(8.2, 4.4))
    x = np.arange(len(axes))
    width = 0.18
    for i, cid in enumerate(primary_ids):
        meds = []
        for axis in axes:
            rec = next(r for r in score_rows if r["contrast_id"] == cid and r["axis"] == axis)
            meds.append(rec["median_log2fc"])
        xpos = x + (i - 1.5) * width
        ax.bar(xpos, meds, width, label=short[cid])
        for xp, med in zip(xpos, meds):
            ax.text(xp, med + (0.03 if med >= 0 else -0.08), f"{med:+.2f}", ha="center", va="bottom" if med >= 0 else "top", fontsize=6, rotation=90)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_ylim(-0.85, 1.05)
    ax.set_xticks(x)
    ax.set_xticklabels(["NHEJ\n(expect down)", "STING\n(expect up)", "IFN\n(expect up)", "APM\n(expect up)"])
    ax.set_ylabel("Median log2FC (KD/KO vs control)")
    ax.set_title("CLDN4 loss: NHEJ, STING, IFN, APM")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_signature_medians.png", dpi=160)
    plt.close()

    payload = {
        "n_uninjured_genes": n_un,
        "n_gse22493_genes": n224,
        "cldn4_uninjured_log2fc": cldn4_base,
        "cldn4_vili_wt_log2fc": vwt["cldn4"]["log2fc"],
        "cldn4_vili_low_log2fc": vlow["cldn4"]["log2fc"],
        "cldn4_vili_high_log2fc": vhigh["cldn4"]["log2fc"],
        "summaries": [
            {
                "contrast_id": s["contrast_id"],
                "pattern": s["pattern"],
                "n_axes_support": s["n_axes_support"],
                "primary": s["primary"],
            }
            for s in summaries
        ],
    }
    (OUT / "key_stats.json").write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))
    print("--- scores ---")
    for r in score_rows:
        print(
            f"{r['contrast_id']:36} {r['axis']:5} {r['call']:12} med={r['median_log2fc']:+.3f} "
            f"up/down={r['n_up']}/{r['n_down']} n={r['n_measured']} "
            f"W={r['wilcoxon_p_expected_direction']} MW={r['mw_p_vs_background']} fdr={r['n_fdr05_expected_direction']}"
        )


if __name__ == "__main__":
    main()
