#!/usr/bin/env python3
"""GSE309652 NSCLC NanoString Metabolic Pathways: CLDN4 vs CD8A / ImmuneScore / CD274 / IFN.

Additive public bulk slice. Unit is the array. CLDN4-only (no dual-high).
Downloads stay under $GSE309652_CLDN4_DATA (default /tmp/gse309652_cldn4)
and are not committed.

GPL31904 is the 768-gene nCounter Human Metabolic Pathways Panel.
CLDN4 is not on that panel. The script still inventories honest n,
histology, ICI labels, and the immune genes that *are* deposited.
"""

from __future__ import annotations

import gzip
import io
import json
import os
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
DATA = Path(os.environ.get("GSE309652_CLDN4_DATA", "/tmp/gse309652_cldn4"))
SEED = 20260817
N_BOOT = 2000

GEO_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE309nnn/GSE309652/"
    "matrix/GSE309652_series_matrix.txt.gz"
)
GPL_SOFT = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL31nnn/GPL31904/"
    "soft/GPL31904_family.soft.gz"
)

# Locked lists from the same-repo CLDN4 slices (OncoSG / Ayers / GSE4573).
# Do not silently substitute missing members.
PRIMARY_GENES = ["CLDN4", "CD8A", "CD274", "IFNG"]
A1_IMMUNE = ["CD8A", "GZMA", "GZMB", "IFNG", "EOMES", "CXCL9", "CXCL10", "TBX21"]
AYERS6_IFN = ["IFNG", "STAT1", "CXCL9", "CXCL10", "IDO1", "HLA-DRA"]
EPITHELIAL = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "KRT7"]
COMPANION = ["TACSTD2"]
CLAUDINS = ["CLDN1", "CLDN3", "CLDN4", "CLDN7"]

HOLDS_N = 40


def dl(dest: Path, url: str) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    print(f"[get] {url}", flush=True)
    import urllib.request

    req = urllib.request.Request(url, headers={"User-Agent": "gse309652-cldn4/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r, open(dest, "wb") as out:
        out.write(r.read())
    return dest


def parse_series_meta(path: Path) -> tuple[pd.DataFrame, dict]:
    series: dict[str, list[str]] = {}
    sample_fields: dict[str, list[list[str]]] = defaultdict(list)
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!Series_"):
                k = line.split("\t", 1)[0][8:]
                v = line.split("\t", 1)[1].strip().strip('"') if "\t" in line else ""
                series.setdefault(k, []).append(v)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][8:]
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                sample_fields[key].append(vals)
    gsms = sample_fields["geo_accession"][0]
    titles = sample_fields["title"][0]
    sources = sample_fields["source_name_ch1"][0]
    rows = []
    for i, gsm in enumerate(gsms):
        rec = {
            "gsm": gsm,
            "title": titles[i],
            "source_name": sources[i],
        }
        for char_row in sample_fields.get("characteristics_ch1", []):
            val = char_row[i]
            if val and ":" in val:
                k, v = val.split(":", 1)
                rec[k.strip().lower()] = v.strip()
        rows.append(rec)
    meta = pd.DataFrame(rows).set_index("gsm")
    meta.index.name = "gsm"
    return meta, {k: " | ".join(v) for k, v in series.items()}


def parse_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as fh:
        text = fh.read()
    lines = text.splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_begin")) + 1
    end = next(i for i, l in enumerate(lines) if l.startswith("!series_matrix_table_end"))
    df = pd.read_csv(io.StringIO("\n".join(lines[start:end])), sep="\t")
    id_col = df.columns[0]
    df = df.set_index(id_col)
    df.index = df.index.astype(str).str.strip().str.strip('"')
    df.columns = [str(c).strip().strip('"') for c in df.columns]
    return df.apply(pd.to_numeric, errors="coerce")


def gpl_has_symbol(path: Path, symbol: str) -> bool:
    """True only if the GPL platform table contains the symbol."""
    needle = symbol.upper()
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table and needle in line.upper():
                return True
    return False


def mean_z(expr: pd.DataFrame, genes: list[str]) -> tuple[pd.Series, list[str]]:
    present = [g for g in genes if g in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns, name="mean_z"), present
    z = expr.loc[present].astype(float).apply(lambda r: (r - r.mean()) / (r.std(ddof=0) or np.nan), axis=1)
    return z.mean(axis=0).rename("mean_z"), present


def spearman_ci(x, y, n_boot: int = N_BOOT, seed: int = SEED):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = int(x.size)
    if n < 8:
        return dict(n=n, rho=np.nan, p=np.nan, ci_low=np.nan, ci_high=np.nan)
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        if np.unique(x[idx]).size < 3 or np.unique(y[idx]).size < 3:
            boots[i] = np.nan
            continue
        boots[i] = stats.spearmanr(x[idx], y[idx])[0]
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return dict(n=n, rho=float(rho), p=float(p), ci_low=float(lo), ci_high=float(hi))


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_rho(r) -> str:
    if r is None or not np.isfinite(r):
        return "—"
    return f"{r:+.3f}"


def histo_map(tissue: str) -> str:
    t = str(tissue).strip()
    if t == "Adenocarcinoma":
        return "LUAD"
    if t == "Squamous cell carcinoma":
        return "LUSC"
    return "Other"


def present_missing(genes: list[str], index) -> tuple[list[str], list[str]]:
    present = [g for g in genes if g in index]
    missing = [g for g in genes if g not in index]
    return present, missing


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)

    matrix_path = dl(DATA / "GSE309652_series_matrix.txt.gz", GEO_MATRIX)
    gpl_path = dl(DATA / "GPL31904_family.soft.gz", GPL_SOFT)

    meta, series = parse_series_meta(matrix_path)
    expr = parse_matrix(matrix_path)
    expr = expr.reindex(columns=meta.index)

    n_array = int(expr.shape[1])
    n_genes = int(expr.shape[0])
    n_missing = int(expr.isna().sum().sum())
    assert n_array == 72, f"expected 72 arrays, got {n_array}"
    assert n_genes == 768, f"expected 768 genes, got {n_genes}"

    meta["histology"] = meta["tissue"].map(histo_map)
    meta["pdl1_num"] = pd.to_numeric(meta.get("pdl1"), errors="coerce")
    meta["suvmax_num"] = pd.to_numeric(meta.get("suvmax"), errors="coerce")
    meta["response"] = meta.get("response", pd.Series(index=meta.index, dtype=str)).astype(str)
    n_luad = int((meta["histology"] == "LUAD").sum())
    n_lusc = int((meta["histology"] == "LUSC").sum())
    n_other = int((meta["histology"] == "Other").sum())
    n_r = int((meta["response"] == "R").sum())
    n_nr = int((meta["response"] == "NR").sum())
    n_pdl1 = int(meta["pdl1_num"].notna().sum())

    cldn4_on_matrix = "CLDN4" in expr.index
    cldn4_on_gpl = gpl_has_symbol(gpl_path, "CLDN4")
    # also scan the deposited symbol list (matrix IDs *are* gene symbols)
    any_cldn = [g for g in expr.index if str(g).upper().startswith("CLDN")]

    a1_present, a1_missing = present_missing(A1_IMMUNE, expr.index)
    ifn_present, ifn_missing = present_missing(AYERS6_IFN, expr.index)
    epi_present, epi_missing = present_missing(EPITHELIAL, expr.index)

    immune_score, a1_used = mean_z(expr, A1_IMMUNE)
    ifn_score, ifn_used = mean_z(expr, AYERS6_IFN)
    epi_score, epi_used = mean_z(expr, EPITHELIAL)

    # Inventory scores live on the 72 columns even though CLDN4 does not.
    scores = pd.DataFrame(
        {
            "ImmuneScore_A1_7of8": immune_score,
            "IFN_Ayers_4of6": ifn_score,
            "epithelial_mean_z": epi_score,
        },
        index=meta.index,
    )
    for g in ["CD8A", "CD274", "IFNG", "PTPRC", "STAT1", "CXCL9", "IDO1"]:
        if g in expr.index:
            scores[g] = expr.loc[g].astype(float)

    annot = meta.join(scores)
    annot.to_csv(TABLES / "sample_annotation.tsv", sep="\t")

    # Gene coverage (requested + locked lists).
    coverage_rows = []
    for label, genes in [
        ("primary", PRIMARY_GENES),
        ("ImmuneScore_A1", A1_IMMUNE),
        ("IFN_Ayers6", AYERS6_IFN),
        ("epithelial", EPITHELIAL),
        ("companion", COMPANION),
        ("claudin_family", CLAUDINS),
    ]:
        present, missing = present_missing(genes, expr.index)
        coverage_rows.append(
            {
                "set": label,
                "n_listed": len(genes),
                "n_present": len(present),
                "present": ",".join(present),
                "missing": ",".join(missing),
            }
        )
    cov = pd.DataFrame(coverage_rows)
    cov.to_csv(TABLES / "gene_coverage.tsv", sep="\t", index=False)

    probe_rows = []
    for g in PRIMARY_GENES + A1_IMMUNE + AYERS6_IFN + EPITHELIAL + COMPANION + CLAUDINS:
        on = g in expr.index
        n_finite = int(np.isfinite(expr.loc[g].astype(float)).sum()) if on else 0
        probe_rows.append(
            {
                "gene": g,
                "on_GPL31904_matrix": int(on),
                "n_finite": n_finite,
                "mean_count": float(expr.loc[g].astype(float).mean()) if on else np.nan,
                "median_count": float(expr.loc[g].astype(float).median()) if on else np.nan,
            }
        )
    probe = pd.DataFrame(probe_rows).drop_duplicates("gene")
    probe.to_csv(TABLES / "probe_confirm.tsv", sep="\t", index=False)

    # Label inventory — public vs missing.
    inv = pd.DataFrame(
        [
            {"item": "arrays in series matrix", "public": "yes", "n": n_array, "note": "GSM9271195–GSM9271266; 768 genes × 72 GSM"},
            {"item": "unique GSM / unique titles", "public": "yes", "n": n_array, "note": "all unique pathology IDs"},
            {"item": "source / tissue", "public": "yes", "n": n_array, "note": "archival tumor; source_name = histology; no adjacent-normal rows"},
            {"item": "histology LUAD (Adenocarcinoma)", "public": "yes", "n": n_luad, "note": "GEO tissue characteristic, not inferred"},
            {"item": "histology LUSC (Squamous cell carcinoma)", "public": "yes", "n": n_lusc, "note": "GEO tissue characteristic, not inferred"},
            {"item": "histology Other", "public": "yes", "n": n_other, "note": "GEO code Other; not recoded to LUAD/LUSC"},
            {"item": "anti-PD-(L)1 ICI cohort", "public": "yes", "n": n_array, "note": "stage IV NSCLC; pretreatment archival tumor"},
            {"item": "ICI response R / NR", "public": "yes", "n": n_r + n_nr, "note": f"{n_r} R / {n_nr} NR"},
            {"item": "PD-L1 IHC numeric", "public": "yes", "n": n_pdl1, "note": "GEO pdl1 percent"},
            {"item": "SUVmax numeric", "public": "yes", "n": int(meta["suvmax_num"].notna().sum()), "note": "PET; not this claim"},
            {"item": "EGFR / ALK / p53 / PI3K", "public": "yes", "n": n_array, "note": "binary/coded on GEO; not this claim"},
            {"item": "age <65 / >=65", "public": "yes", "n": n_array, "note": "binned only"},
            {"item": "stage (beyond IV enrollment)", "public": "no", "n": 0, "note": "cohort is stage IV; no TNM on GEO"},
            {"item": "tumor % / ESTIMATE purity", "public": "no", "n": 0, "note": "not deposited; epithelial 6-gene set is 0/6"},
            {"item": "OS / PFS time", "public": "no", "n": 0, "note": "design text mentions survival; times are not GEO characteristics"},
            {"item": "CLDN4 on matrix / GPL31904", "public": "no", "n": 0, "note": "0/768 symbols; 0 CLDN* genes; GPL SOFT has no CLDN4"},
            {"item": "CD8A finite", "public": "yes", "n": n_array, "note": "symbol CD8A"},
            {"item": "CD274 finite", "public": "yes", "n": n_array, "note": "symbol CD274"},
            {"item": "IFNG finite", "public": "yes", "n": n_array, "note": "symbol IFNG"},
            {"item": "ImmuneScore A1 8-gene", "public": "partial", "n": n_array, "note": f"{len(a1_used)}/8 present; missing {','.join(a1_missing) or 'none'}"},
            {"item": "IFN Ayers-6", "public": "partial", "n": n_array, "note": f"{len(ifn_used)}/6 present; missing {','.join(ifn_missing) or 'none'}"},
            {"item": "epithelial mean-z (6-gene)", "public": "no", "n": 0, "note": f"0/6; missing {','.join(epi_missing)}; KRT1 is on panel and is not used as a substitute"},
            {"item": "TACSTD2 companion", "public": "no", "n": 0, "note": "not on GPL31904; not substituted"},
            {"item": "primary pairwise n (CLDN4 + any partner)", "public": "no", "n": 0, "note": "CLDN4 row does not exist"},
        ]
    )
    inv.to_csv(TABLES / "label_inventory.tsv", sep="\t", index=False)

    # CLDN4 pairwise table: honest ABSENT rows. No fabricated ρ.
    partners = [
        ("CD8A", "CD8A" in expr.index),
        ("ImmuneScore_A1_7of8", len(a1_used) > 0),
        ("CD274", "CD274" in expr.index),
        ("IFN_Ayers_4of6", len(ifn_used) > 0),
        ("IFNG", "IFNG" in expr.index),
    ]
    subsets = [
        ("all_tumors", meta.index),
        ("LUAD", meta.index[meta["histology"] == "LUAD"]),
        ("LUSC", meta.index[meta["histology"] == "LUSC"]),
        ("Other", meta.index[meta["histology"] == "Other"]),
    ]
    spearman_rows = []
    for subset_name, idx in subsets:
        for partner, partner_ok in partners:
            spearman_rows.append(
                {
                    "gene_x": "CLDN4",
                    "gene_y": partner,
                    "subset": subset_name,
                    "n_arrays": int(len(idx)),
                    "n": 0,
                    "rho": np.nan,
                    "p": np.nan,
                    "ci_low": np.nan,
                    "ci_high": np.nan,
                    "rho_adj_epithelial": np.nan,
                    "p_adj_epithelial": np.nan,
                    "partner_on_panel": int(partner_ok),
                    "cldn4_on_panel": 0,
                    "epithelial_genes": f"{len(epi_used)}/6",
                    "verdict": "ABSENT",
                    "role": "primary",
                    "note": "CLDN4 not on GPL31904 Metabolic Pathways Panel",
                }
            )
    stats_df = pd.DataFrame(spearman_rows)
    stats_df.to_csv(TABLES / "spearman_cldn4_vs_cd8a_immunescore_cd274_ifn.tsv", sep="\t", index=False)

    # Panel-inventory Spearman among genes that *are* present (not the claim).
    inventory_pairs = []
    if "CD8A" in expr.index and "CD274" in expr.index:
        inventory_pairs.append(("CD8A", "CD274"))
    if "CD8A" in expr.index and "IFNG" in expr.index:
        inventory_pairs.append(("CD8A", "IFNG"))
    if "CD274" in expr.index and "IFNG" in expr.index:
        inventory_pairs.append(("CD274", "IFNG"))
    if "CD8A" in expr.index:
        inventory_pairs.append(("CD8A", "ImmuneScore_A1_7of8"))
    if "CD274" in expr.index:
        inventory_pairs.append(("CD274", "ImmuneScore_A1_7of8"))
    if n_pdl1 >= 8 and "CD274" in expr.index:
        inventory_pairs.append(("CD274", "pdl1_ihc"))

    inv_rows = []
    for a, b in inventory_pairs:
        if a == "pdl1_ihc" or b == "pdl1_ihc":
            x = annot["CD274"] if a != "pdl1_ihc" else annot["pdl1_num"]
            y = annot["pdl1_num"] if b == "pdl1_ihc" else annot[a if a != "CD274" else "CD274"]
            if b == "pdl1_ihc":
                x = annot[a]
                y = annot["pdl1_num"]
        else:
            x = annot[a] if a in annot.columns else expr.loc[a]
            y = annot[b] if b in annot.columns else expr.loc[b]
        s = spearman_ci(x.values, y.values)
        inv_rows.append({"gene_x": a, "gene_y": b, "subset": "all_tumors", "role": "panel_inventory_not_claim", **s})
    inv_stats = pd.DataFrame(inv_rows)
    inv_stats.to_csv(TABLES / "spearman_panel_inventory.tsv", sep="\t", index=False)

    # Histology contrast on deposited immune genes (inventory).
    histo_rows = []
    for g in ["CD8A", "CD274", "IFNG", "ImmuneScore_A1_7of8", "IFN_Ayers_4of6"]:
        if g not in annot.columns:
            continue
        a = annot.loc[annot["histology"] == "LUAD", g].astype(float)
        b = annot.loc[annot["histology"] == "LUSC", g].astype(float)
        if a.notna().sum() >= 5 and b.notna().sum() >= 5:
            u, p = stats.mannwhitneyu(a.dropna(), b.dropna(), alternative="two-sided")
        else:
            u, p = np.nan, np.nan
        histo_rows.append(
            {
                "gene": g,
                "n_LUAD": int(a.notna().sum()),
                "n_LUSC": int(b.notna().sum()),
                "median_LUAD": float(a.median()),
                "median_LUSC": float(b.median()),
                "mwu_U": float(u) if np.isfinite(u) else np.nan,
                "mwu_p": float(p) if np.isfinite(p) else np.nan,
                "role": "panel_inventory_not_claim",
            }
        )
    histo_df = pd.DataFrame(histo_rows)
    histo_df.to_csv(TABLES / "histology_contrast.tsv", sep="\t", index=False)

    # Empty high/low — no CLDN4 to split.
    hl = pd.DataFrame(
        [
            {
                "axis": axis,
                "n_Q4": 0,
                "n_Q1": 0,
                "mwu_p": np.nan,
                "verdict": "ABSENT",
                "note": "no CLDN4 values to quartile",
            }
            for axis in ["CD8A", "ImmuneScore_A1_7of8", "CD274", "IFN_Ayers_4of6"]
        ]
    )
    hl.to_csv(TABLES / "highlow_cldn4.tsv", sep="\t", index=False)

    one_row = pd.DataFrame(
        [
            {
                "dataset": "GSE309652",
                "platform": "GPL31904 NanoString nCounter Human Metabolic Pathways Panel (768 genes)",
                "histology": f"NSCLC mixed (LUAD n={n_luad}; LUSC n={n_lusc}; Other n={n_other})",
                "n_arrays": n_array,
                "n_CLDN4": 0,
                "n_CLDN4_CD8A": 0,
                "n_CLDN4_ImmuneScore": 0,
                "n_CLDN4_CD274": 0,
                "n_CLDN4_IFN": 0,
                "CLDN4": "ABSENT",
                "CD8A": "present (768-gene symbol CD8A; n=72 finite)",
                "ImmuneScore": f"A1 8-gene {len(a1_used)}/8 (missing {','.join(a1_missing)}); ESTIMATE not computable",
                "CD274": "present (symbol CD274; n=72 finite)",
                "IFN": f"Ayers-6 {len(ifn_used)}/6 (missing {','.join(ifn_missing)}); IFNG present n=72",
                "epithelial_residual": f"NOT_COMPUTABLE ({len(epi_used)}/6; {','.join(epi_missing)})",
                "CLDN4_CD8A_rho": np.nan,
                "CLDN4_CD8A_p": np.nan,
                "CLDN4_ImmuneScore_rho": np.nan,
                "CLDN4_ImmuneScore_p": np.nan,
                "CLDN4_CD274_rho": np.nan,
                "CLDN4_CD274_p": np.nan,
                "CLDN4_IFN_rho": np.nan,
                "CLDN4_IFN_p": np.nan,
                "CLDN4_CD8A_verdict": "ABSENT",
                "CLDN4_ImmuneScore_verdict": "ABSENT",
                "CLDN4_CD274_verdict": "ABSENT",
                "CLDN4_IFN_verdict": "ABSENT",
                "headline_verdict": "ABSENT",
                "ICI": f"stage IV anti-PD-(L)1; response {n_r} R / {n_nr} NR (CLDN4 vs response not testable)",
                "dual_high": "not used",
                "note": "CLDN4 is not on GPL31904; GPL SOFT CLDN4 hit=0; no CLDN* symbols in the 768-gene matrix",
            }
        ]
    )
    one_row.to_csv(TABLES / "one_row.tsv", sep="\t", index=False)

    # ---- figures ----
    # fig1: coverage of requested / locked genes
    fig_genes = [
        ("CLDN4", 1 if "CLDN4" in expr.index else 0, "1/1" if "CLDN4" in expr.index else "0/1"),
        ("CD8A", 1 if "CD8A" in expr.index else 0, "1/1" if "CD8A" in expr.index else "0/1"),
        ("CD274", 1 if "CD274" in expr.index else 0, "1/1" if "CD274" in expr.index else "0/1"),
        ("IFNG", 1 if "IFNG" in expr.index else 0, "1/1" if "IFNG" in expr.index else "0/1"),
        ("ImmuneScore\n(A1 8-gene)", len(a1_used) / 8, f"{len(a1_used)}/8"),
        ("IFN\n(Ayers-6)", len(ifn_used) / 6, f"{len(ifn_used)}/6"),
        ("epithelial\n(6-gene)", len(epi_used) / 6, f"{len(epi_used)}/6"),
        ("TACSTD2", 1 if "TACSTD2" in expr.index else 0, "1/1" if "TACSTD2" in expr.index else "0/1"),
    ]
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    names = [a for a, _, _ in fig_genes]
    heights = np.array([h for _, h, _ in fig_genes], float)
    colors = []
    for h in heights:
        if h >= 1:
            colors.append("#2c5f8a")
        elif h > 0:
            colors.append("#c48a2a")
        else:
            colors.append("#b33a3a")
    ax.bar(range(len(names)), heights, color=colors, edgecolor="none")
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names)
    ax.set_ylim(-0.05, 1.25)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["absent", "partial", "present"])
    ax.set_title(
        f"GSE309652  GPL31904 Metabolic Pathways  n={n_array} arrays\n"
        f"CLDN4 ABSENT  ·  A1 ImmuneScore {len(a1_used)}/8  ·  Ayers-6 {len(ifn_used)}/6  ·  epi {len(epi_used)}/6"
    )
    for i, (_, h, lab) in enumerate(fig_genes):
        ax.text(i, min(1.05, h + 0.06), lab, ha="center", va="bottom", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig1_gene_coverage.png", dpi=160)
    fig.savefig(FIGURES / "fig1_gene_coverage.pdf")
    plt.close(fig)

    # fig2: deposited immune genes by histology (inventory; not a CLDN4 test)
    fig, axes = plt.subplots(1, 3, figsize=(10.2, 4.0))
    order = ["LUAD", "LUSC", "Other"]
    palette = {"LUAD": "#2c5f8a", "LUSC": "#6b3a2a", "Other": "#6a6a6a"}
    for ax, gene, ylab in zip(
        axes,
        ["CD8A", "CD274", "IFNG"],
        ["CD8A  nSolver count", "CD274  nSolver count", "IFNG  nSolver count"],
    ):
        data = [annot.loc[annot["histology"] == h, gene].astype(float).values for h in order]
        parts = ax.boxplot(data, tick_labels=order, patch_artist=True, widths=0.55, showfliers=True)
        for patch, h in zip(parts["boxes"], order):
            patch.set_facecolor(palette[h])
            patch.set_alpha(0.55)
        ax.set_ylabel(ylab)
        ax.set_title(f"{gene}  (panel inventory)")
    fig.suptitle(
        "GSE309652  immune genes that *are* on GPL31904\n"
        "Not a CLDN4 test — CLDN4 row does not exist",
        fontsize=10,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_immune_by_histology.png", dpi=160)
    fig.savefig(FIGURES / "fig2_immune_by_histology.pdf")
    plt.close(fig)

    # fig3: one-row visual
    fig, ax = plt.subplots(figsize=(11.2, 3.6))
    ax.axis("off")
    col_labels = [
        "dataset",
        "n arrays",
        "LUAD/LUSC/Other",
        "CLDN4",
        "CD8A",
        "ImmuneScore",
        "CD274",
        "IFN",
        "epi residual",
        "verdict",
    ]
    cell = [
        [
            "GSE309652",
            str(n_array),
            f"{n_luad}/{n_lusc}/{n_other}",
            "ABSENT",
            "present",
            f"{len(a1_used)}/8",
            "present",
            f"{len(ifn_used)}/6",
            f"{len(epi_used)}/6",
            "ABSENT",
        ]
    ]
    table = ax.table(cellText=cell, colLabels=col_labels, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1.05, 2.0)
    for (r, c), cell_obj in table.get_celld().items():
        cell_obj.set_edgecolor("#cccccc")
        if r == 0:
            cell_obj.set_facecolor("#2c5f8a")
            cell_obj.set_text_props(color="white", weight="bold")
        elif c in (3, 9):
            cell_obj.set_facecolor("#f2d6d6")
            cell_obj.set_text_props(weight="bold")
        else:
            cell_obj.set_facecolor("#f7f7f7")
    ax.set_title(
        "GSE309652 one-row  ·  CLDN4-only  ·  no dual-high\n"
        "NanoString Metabolic Pathways Panel does not measure CLDN4",
        pad=8,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_one_row.png", dpi=160)
    fig.savefig(FIGURES / "fig3_one_row.pdf")
    plt.close(fig)

    summary = {
        "dataset": "GSE309652",
        "pmid_on_geo": None,
        "title": series.get("title", ""),
        "platform": "GPL31904 NanoString nCounter Human Metabolic Pathways Panel",
        "processing": (
            "nSolver positive-control + housekeeping normalization as deposited. "
            "Values are linear nSolver counts. Spearman is rank-based."
        ),
        "n_arrays": n_array,
        "n_genes": n_genes,
        "n_missing_values": n_missing,
        "n_LUAD": n_luad,
        "n_LUSC": n_lusc,
        "n_Other": n_other,
        "n_response_R": n_r,
        "n_response_NR": n_nr,
        "n_pdl1_numeric": n_pdl1,
        "cldn4_on_matrix": cldn4_on_matrix,
        "cldn4_on_gpl_soft": cldn4_on_gpl,
        "cldn_family_on_matrix": any_cldn,
        "ImmuneScore_A1_present": a1_used,
        "ImmuneScore_A1_missing": a1_missing,
        "IFN_Ayers6_present": ifn_used,
        "IFN_Ayers6_missing": ifn_missing,
        "epithelial_present": epi_used,
        "epithelial_missing": epi_missing,
        "primary_verdict": "ABSENT",
        "dual_high": False,
        "holds_rule": "n>=40, partial Spearman ρ<0, p<0.05 (PR 229 / GSE4573). Not applicable: CLDN4 n=0.",
        "panel_inventory_spearman": inv_stats.to_dict(orient="records"),
        "histology_contrast": histo_df.to_dict(orient="records"),
        "geo_matrix_url": GEO_MATRIX,
        "gpl_soft_url": GPL_SOFT,
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    print("n_arrays", n_array, "LUAD", n_luad, "LUSC", n_lusc, "Other", n_other)
    print("CLDN4 on matrix", cldn4_on_matrix, "on GPL SOFT", cldn4_on_gpl, "CLDN* ", any_cldn)
    print("A1", a1_used, "missing", a1_missing)
    print("IFN", ifn_used, "missing", ifn_missing)
    print("epi", epi_used, "missing", epi_missing)
    print("response R/NR", n_r, n_nr, "pdl1", n_pdl1, "missing_expr", n_missing)
    print(inv_stats.to_string(index=False))
    print("[done]", TABLES)


if __name__ == "__main__":
    main()
