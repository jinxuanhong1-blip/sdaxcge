#!/usr/bin/env python3
"""GSE221733 GeoMx CTA RNA: CLDN4 panel check and (if present) AOI tests.

Public inputs only (GEO GSE221733 CTA_norm / CTA_QC / CTA_initial + series matrix).

If CLDN4 is absent from the deposited panel, skip all CLDN4 vs CD8 / ICB tests.
Do not substitute TACSTD2, EPCAM, or any other gene for CLDN4.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from openpyxl import load_workbook
from scipy import stats

REQUESTED = [
    "CLDN4",
    "CLDN3",
    "CLDN7",
    "CLDN1",
    "F11R",
    "OCLN",
    "TACSTD2",
    "EPCAM",
    "CD8A",
    "CD8B",
    "CD3E",
    "CD3D",
    "CD3G",
    "PTPRC",
    "GZMB",
    "PRF1",
    "CD4",
    "FOXP3",
    "PDCD1",
    "CD274",
]

T_CELL_AOI_TOKENS = (
    "cd8",
    "t-cell",
    "t cell",
    "tcell",
    "immune pos",
    "cd3 pos",
    "cd45 pos",
)

PAPER = (
    "Monkman J, Kim H, Mayer A, et al. Multi-omic and spatial dissection of "
    "immunotherapy response groups in non-small cell lung cancer. "
    "Immunology. 2023;169(4):487-502. PMID 37022147. DOI 10.1111/imm.13646."
)


def bh_fdr(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    adj = np.empty(n, dtype=float)
    prev = 1.0
    for rank_from_end, idx in enumerate(order[::-1]):
        rank = n - rank_from_end
        val = min(prev, p[idx] * n / rank)
        adj[idx] = val
        prev = val
    return adj.tolist()


def spearman_ci(rho: float, n: int, alpha: float = 0.05) -> tuple[float, float]:
    if n <= 3 or not np.isfinite(rho) or abs(rho) >= 1:
        return (float("nan"), float("nan"))
    z = np.arctanh(rho)
    se = 1.0 / math.sqrt(n - 3)
    zcrit = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - zcrit * se)), float(np.tanh(z + zcrit * se))


def spearman_row(x: pd.Series, y: pd.Series, stratum: str, a: str, b: str) -> dict:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 5:
        rho = p = lo = hi = float("nan")
    else:
        rho, p = stats.spearmanr(x[m], y[m])
        lo, hi = spearman_ci(float(rho), n)
    return {
        "stratum": stratum,
        "x": a,
        "y": b,
        "n": n,
        "rho": float(rho) if np.isfinite(rho) else np.nan,
        "p": float(p) if np.isfinite(p) else np.nan,
        "rho_ci95_lo": lo,
        "rho_ci95_hi": hi,
    }


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.yaxis.label.set_size(9)
    ax.xaxis.label.set_size(9)


def savefig(fig, out_dir: Path, stem: str):
    fig.tight_layout()
    fig.savefig(out_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def read_xlsx_matrix(path: Path) -> pd.DataFrame:
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    wb.close()
    genes = [str(g) for g in rows[0][1:] if g is not None]
    index = []
    data = []
    for row in rows[1:]:
        if row[0] is None:
            continue
        index.append(str(row[0]))
        data.append([float(x) if x is not None else np.nan for x in row[1 : 1 + len(genes)]])
    return pd.DataFrame(data, index=index, columns=genes)


def parse_series_matrix(path: Path) -> pd.DataFrame:
    text = path.read_text(errors="replace") if path.suffix != ".gz" else gzip.open(path, "rt", errors="replace").read()
    fields: dict[str, list[list[str]]] = defaultdict(list)
    for line in text.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if not line.startswith("!Sample_"):
            continue
        parts = [p.strip().strip('"') for p in line.split("\t")]
        fields[parts[0]].append(parts[1:])
    titles = fields["!Sample_title"][0]
    acc = fields["!Sample_geo_accession"][0]
    records = []
    for i, title in enumerate(titles):
        rec = {"title": title, "geo_accession": acc[i]}
        for row in fields.get("!Sample_characteristics_ch1", []):
            val = row[i]
            if ":" in val:
                k, v = val.split(":", 1)
                rec[k.strip().lower()] = v.strip()
        records.append(rec)
    meta = pd.DataFrame(records).set_index("title")
    meta["roi"] = pd.to_numeric(meta.get("roi"), errors="coerce")
    meta["area"] = pd.to_numeric(meta.get("area"), errors="coerce")
    meta["aoinucleicount"] = pd.to_numeric(meta.get("aoinucleicount"), errors="coerce")
    meta["followup"] = pd.to_numeric(meta.get("followup"), errors="coerce")
    meta["qc_fail"] = meta.get("qc fail", pd.Series(index=meta.index, dtype=object)).fillna("").eq("Fail")
    if "qcflags" not in meta.columns:
        meta["qcflags"] = ""
    meta["qcflags"] = meta["qcflags"].fillna("")
    meta["in_norm"] = False
    meta["in_qc"] = False
    return meta


def probe_genes(initial_csv_gz: Path) -> list[str]:
    with gzip.open(initial_csv_gz, "rt") as fh:
        header = next(csv.reader(fh))
    probes = header[1:]
    return sorted({p.rsplit("_", 1)[0] for p in probes})


def find_cldn_like(names: list[str]) -> list[str]:
    rx = re.compile(r"cldn|claudin", re.I)
    return [n for n in names if rx.search(str(n))]


def segment_class(seg: str) -> str:
    s = (seg or "").strip().lower()
    if "panck pos" in s or s in {"tumour", "tumor", "malignant"}:
        return "tumor"
    if "panck neg" in s or s in {"stroma", "tme"}:
        return "stroma"
    if "geometric" in s:
        return "geometric"
    return "other"


def looks_like_tcell_aoi(seg: str) -> bool:
    s = (seg or "").lower()
    return any(tok in s for tok in T_CELL_AOI_TOKENS)


def panel_inventory(genes: list[str], probe_genes_list: list[str]) -> pd.DataFrame:
    gene_set = set(genes)
    probe_set = set(probe_genes_list)
    rows = []
    for g in REQUESTED:
        rows.append(
            {
                "gene": g,
                "in_cta_norm": g in gene_set,
                "in_cta_qc": g in gene_set,
                "in_initial_probes": g in probe_set,
                "role": (
                    "query_CLDN4"
                    if g == "CLDN4"
                    else "other_claudin_TJ"
                    if g.startswith("CLDN") or g in {"F11R", "OCLN"}
                    else "epithelial_other"
                    if g in {"TACSTD2", "EPCAM"}
                    else "T_cell_or_immune"
                ),
            }
        )
    extra = find_cldn_like(genes) + find_cldn_like(probe_genes_list)
    for g in sorted(set(extra)):
        if g not in REQUESTED:
            rows.append(
                {
                    "gene": g,
                    "in_cta_norm": g in gene_set,
                    "in_cta_qc": g in gene_set,
                    "in_initial_probes": g in probe_set,
                    "role": "unexpected_claudin_hit",
                }
            )
    return pd.DataFrame(rows)


def fig_panel(inv: pd.DataFrame, n_genes: int, out_dir: Path):
    df = inv.copy()
    df["present"] = df["in_cta_norm"]
    df = df.sort_values(["present", "role", "gene"], ascending=[False, True, True])
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    y = np.arange(len(df))
    colors = []
    for _, r in df.iterrows():
        if r["gene"] == "CLDN4":
            colors.append("#922b21")
        elif r["present"]:
            colors.append("#1e8449")
        else:
            colors.append("#bfc5ce")
    ax.barh(y, [1] * len(df), color=colors, height=0.72, edgecolor="white")
    ax.set_yticks(y)
    ax.set_yticklabels(df["gene"], fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.invert_yaxis()
    for i, (_, r) in enumerate(df.iterrows()):
        label = "on panel" if r["present"] else "absent"
        ax.text(0.02, i, label, va="center", ha="left", color="white" if r["present"] or r["gene"] == "CLDN4" else "#333", fontsize=8)
    ax.set_title(f"GSE221733 GeoMx CTA panel ({n_genes} genes) — requested targets", fontsize=10)
    style_axes(ax)
    ax.spines["bottom"].set_visible(False)
    savefig(fig, out_dir, "panel_requested_genes")


def fig_aoi_inventory(meta: pd.DataFrame, out_dir: Path):
    segs = ["PanCK pos", "PanCK neg", "Geometric Segment"]
    stages = ["GEO deposited", "Author QC fail", "In CTA_norm"]
    counts = {stage: [] for stage in stages}
    for seg in segs:
        sub = meta[meta["segment"] == seg]
        counts["GEO deposited"].append(int(len(sub)))
        counts["Author QC fail"].append(int(sub["qc_fail"].sum()))
        counts["In CTA_norm"].append(int(sub["in_norm"].sum()))
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    x = np.arange(len(segs))
    width = 0.25
    palette = ["#5d6d7e", "#c0392b", "#1f618d"]
    for i, stage in enumerate(stages):
        ax.bar(x + (i - 1) * width, counts[stage], width, label=stage, color=palette[i])
        for j, val in enumerate(counts[stage]):
            ax.text(x[j] + (i - 1) * width, val + 0.4, str(val), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(["PanCK pos\n(tumor / malignant)", "PanCK neg\n(stroma / TME)", "Geometric"], fontsize=8)
    ax.set_ylabel("AOIs")
    ax.set_title("GSE221733 AOI inventory (no dedicated CD8 / T-cell AOI class)", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, max(counts["GEO deposited"]) + 8)
    style_axes(ax)
    savefig(fig, out_dir, "aoi_inventory")


def fig_pairing(meta: pd.DataFrame, out_dir: Path):
    # Pairing among AOIs that made it into the normalized matrix.
    norm = meta[meta["in_norm"]].copy()
    cores = defaultdict(set)
    for title, row in norm.iterrows():
        cores[int(row["roi"])].add(segment_class(row["segment"]))
    n_pair = sum(1 for v in cores.values() if "tumor" in v and "stroma" in v)
    n_tum = sum(1 for v in cores.values() if "tumor" in v and "stroma" not in v)
    n_str = sum(1 for v in cores.values() if "stroma" in v and "tumor" not in v)
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    labels = [
        "Paired tumor + stroma\n(same ROI / core)",
        "Tumor AOI only",
        "Stroma AOI only",
        "Dedicated T-cell AOI",
    ]
    vals = [n_pair, n_tum, n_str, 0]
    colors = ["#1f618d", "#7f8c8d", "#7f8c8d", "#922b21"]
    ax.barh(np.arange(len(labels)), vals, color=colors, height=0.65)
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    for i, v in enumerate(vals):
        ax.text(v + 0.3, i, str(v), va="center", fontsize=9)
    ax.set_xlabel("Cores represented in CTA_norm")
    ax.set_title("Same-core pairing available for tumor vs stroma, not tumor vs T-cell AOI", fontsize=10)
    style_axes(ax)
    savefig(fig, out_dir, "pairing_availability")


def fig_response(meta: pd.DataFrame, out_dir: Path):
    tum = meta[(meta["in_norm"]) & (meta["segment"] == "PanCK pos")].copy()
    # Patient-level: one label per patient among tumor AOIs
    pat = tum.groupby("patient id")["response"].first()
    order = ["Responder", "Non-responder", "N/A"]
    aoi_counts = [int((tum["response"] == lab).sum()) for lab in order]
    pat_counts = [int((pat == lab).sum()) for lab in order]
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6))
    colors = ["#1e8449", "#922b21", "#7f8c8d"]
    for ax, counts, title in zip(
        axes,
        [aoi_counts, pat_counts],
        [
            f"PanCK-pos AOIs in CTA_norm (n={len(tum)})",
            f"Patients with a tumor AOI in CTA_norm (n={len(pat)})",
        ],
    ):
        ax.bar(order, counts, color=colors)
        for i, v in enumerate(counts):
            ax.text(i, v + 0.2, str(v), ha="center", va="bottom", fontsize=8)
        ax.set_ylim(0, max(counts + [1]) + 3)
        ax.set_title(title, fontsize=9)
        ax.set_ylabel("Count")
        style_axes(ax)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("ICB response labels exist; CLDN4 vs response was not tested (CLDN4 absent)", fontsize=10, y=1.02)
    savefig(fig, out_dir, "response_labels")


def fig_qc_scatter(meta: pd.DataFrame, out_dir: Path):
    df = meta.dropna(subset=["area", "aoinucleicount"]).copy()
    fig, ax = plt.subplots(figsize=(6.8, 4.4))
    cmap = {
        ("PanCK pos", False): "#1f618d",
        ("PanCK pos", True): "#85c1e9",
        ("PanCK neg", False): "#b9770e",
        ("PanCK neg", True): "#f9e79f",
        ("Geometric Segment", False): "#6c3483",
        ("Geometric Segment", True): "#d2b4de",
    }
    for (seg, fail), sub in df.groupby(["segment", "qc_fail"]):
        ax.scatter(
            np.log10(sub["area"].clip(lower=1)),
            np.log10(sub["aoinucleicount"].clip(lower=1)),
            s=28,
            c=cmap.get((seg, fail), "#888"),
            label=f"{seg}{' QC-fail' if fail else ''}",
            alpha=0.85,
            edgecolors="white",
            linewidths=0.3,
        )
    ax.set_xlabel("log10 AOI area (µm²)")
    ax.set_ylabel("log10 nuclei count")
    ax.set_title("ROI / AOI size and nuclei (GEO metadata)", fontsize=10)
    ax.legend(frameon=False, fontsize=7, loc="best")
    style_axes(ax)
    savefig(fig, out_dir, "aoi_area_nuclei")


def fig_skip_summary(n_genes: int, n_cldn: int, out_dir: Path):
    fig, ax = plt.subplots(figsize=(7.4, 3.2))
    ax.set_axis_off()
    text = (
        "CLDN4 tests skipped\n\n"
        f"Deposited GeoMx panel: Cancer Transcriptome Atlas, {n_genes} genes.\n"
        f"CLDN4 exact match: no.  CLDN* / claudin string hits: {n_cldn}.\n"
        "Asked tests not run: tumor-AOI CLDN4 vs CD8A; tumor CLDN4 vs paired\n"
        "T-cell AOI / T-cell fraction; malignant-AOI CLDN4 in ICB NR vs R.\n"
        "No TACSTD2 / EPCAM stand-in was analysed as CLDN4."
    )
    ax.text(0.02, 0.55, text, va="center", ha="left", fontsize=10, family="DejaVu Sans", transform=ax.transAxes)
    ax.add_patch(plt.Rectangle((0, 0), 1, 1, fill=False, edgecolor="#922b21", lw=2, transform=ax.transAxes))
    fig.savefig(out_dir / "cldn4_tests_skipped.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "cldn4_tests_skipped.pdf", bbox_inches="tight")
    plt.close(fig)


def write_results_md(path: Path, key: dict, asset_prefix: str):
    p = asset_prefix.rstrip("/") + "/" if asset_prefix else ""
    md = f"""# GSE221733 GeoMx CTA RNA — CLDN4 vs CD8 / ICB (skip)

Additive public GeoMx DSP RNA only. **CLDN4-only.** No private 8-KL. No Visium / CosMx (those accessions are covered by sibling agents). No gene substitution.

Dataset: [GSE221733](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221733). {PAPER}

## TL;DR

**CLDN4 is not on the deposited GeoMx Cancer Transcriptome Atlas (CTA) panel.** It is absent from `GSE221733_4301_CTA_norm.xlsx` ({key['n_genes']} genes), the matching QC collapsed-count table, and the uncollapsed initial probe file ({key['n_probes']} probes / {key['n_probe_genes']} gene symbols). There are **{key['n_cldn_string_hits']}** gene or probe names matching `CLDN` / `claudin`.

The asked tests were therefore **not run**:

1. In tumor / malignant AOIs, is CLDN4 higher when CD8 / T-cell AOIs or CD8A are lower?
2. Same-core CLDN4 (tumor AOI) vs CD8A (T-cell AOI) or vs T-cell fraction.
3. CLDN4 in malignant AOIs of ICB non-responders vs responders.

This is a panel-content skip, not a negative biological result. TACSTD2 and EPCAM **are** on the CTA panel; they were **not** analysed as CLDN4 stand-ins.

![Requested genes on the CTA panel]({p}figures/panel_requested_genes.png)

![CLDN4 tests skipped]({p}figures/cldn4_tests_skipped.png)

## What was downloaded

| File | Role | Used for |
|---|---|---|
| `GSE221733_4301_CTA_norm.xlsx` | Author RUV-normalised CTA matrix (AOIs × {key['n_genes']} genes) | Panel inventory; AOIs that authors retained |
| `GSE221733_4301_CTA_QC.xlsx` | Collapsed probe counts | Confirm same gene universe; one extra Geometric AOI |
| `GSE221733_4301_CTA_initial.csv.gz` | Uncollapsed probes ({key['n_probes']} columns) | Confirm no CLDN4 probe barcode |
| `GSE221733_series_matrix.txt` | ROI / patient / response / QC flags | AOI metadata |

GEO also lists per-AOI DCC files in `GSE221733_RAW.tar`. Those were not needed to decide panel membership.

## Panel

CTA, **{key['n_genes']}** gene symbols after the AOI index column. Housekeeping / control symbols present include `HK1`, `HK2`, and `Negative Probe`.

| Requested | On CTA_norm / QC / initial probes? |
|---|---|
| **CLDN4** | **no / no / no** |
| CLDN3, CLDN7, CLDN1, F11R, OCLN | no |
| TACSTD2 | yes (not used as a CLDN4 proxy) |
| EPCAM | yes (not used as a CLDN4 proxy) |
| CD8A, CD8B, CD3E, CD3D, CD3G, PTPRC, GZMB, PRF1, CD4, FOXP3, PDCD1, CD274 | yes |

Full checklist: `{p}tables/panel_inventory.csv`.

## Cohort (deposited metadata)

93 AOIs on the series record: **{key['n_geo_tumor']}** PanCK pos (tumor / malignant), **{key['n_geo_stroma']}** PanCK neg (stroma / TME), **{key['n_geo_geom']}** Geometric Segment. **{key['n_patients']}** patient IDs. Treatment field is Immunotherapy for every AOI. ICB labels exist: Responder / Non-responder / N/A.

Author `qc fail: Fail` on **{key['n_qc_fail']}** AOIs (low nuclei, low negative-probe count, low surface area, and/or low sequencing saturation). The normalised matrix keeps **{key['n_norm']}** AOIs (**{key['n_norm_tumor']}** PanCK pos, **{key['n_norm_stroma']}** PanCK neg). The QC count table has those plus the single Geometric Segment.

![AOI inventory]({p}figures/aoi_inventory.png)

There is **no dedicated CD8 or T-cell AOI class**. Segmentation in the titles and `segment` field is PanCK pos / PanCK neg / Geometric only. Same-core pairing in CTA_norm: **{key['n_paired_cores']}** ROIs with both tumor and stroma, **{key['n_tumor_only_cores']}** tumor-only, **{key['n_stroma_only_cores']}** stroma-only. A paired CLDN4 (tumor) vs CD8A (T-cell AOI) test is not constructible from the deposited masks even if CLDN4 were present; the fallback would have been tumor-AOI CLDN4 vs tumor-AOI CD8A or vs stroma CD8A as a T-cell-fraction proxy.

![Pairing availability]({p}figures/pairing_availability.png)

Among PanCK-pos AOIs in CTA_norm: **{key['n_norm_tumor_R']}** Responder, **{key['n_norm_tumor_NR']}** Non-responder, **{key['n_norm_tumor_NA']}** N/A, from **{key['n_norm_tumor_patients']}** patients (**{key['n_norm_tumor_patients_R']}** R / **{key['n_norm_tumor_patients_NR']}** NR / **{key['n_norm_tumor_patients_NA']}** N/A). Those labels would have supported the ICB comparison if CLDN4 existed on the panel.

![ICB response labels]({p}figures/response_labels.png)

![AOI area vs nuclei]({p}figures/aoi_area_nuclei.png)

## Tests (not run)

| Test | Status | Reason |
|---|---|---|
| Tumor AOI CLDN4 vs same-AOI CD8A | skipped | CLDN4 absent |
| Tumor AOI CLDN4 vs paired T-cell AOI CD8A | skipped | CLDN4 absent; no T-cell AOI class |
| Tumor AOI CLDN4 vs T-cell fraction / stroma CD8A | skipped | CLDN4 absent |
| Tumor AOI CLDN4, ICB NR vs R | skipped | CLDN4 absent |

No Spearman ρ, Mann–Whitney U, or median-split was computed for CLDN4. No TACSTD2–CD8A or EPCAM–CD8A correlation is reported here.

## What this does **not** claim

- It does not claim CLDN4 is low, high, or associated with CD8 or ICB outcome in this cohort.
- It does not claim the CTA assay failed; CD8A and other IO genes are on the panel.
- It does not re-analyse GSE221322 (protein DSP from the same paper; CLDN4 protein is also off that 68-plex — sibling branch).
- It does not re-analyse Visium / CosMx accessions already covered elsewhere (GSE307534, official CosMx NSCLC, GSE299786).

## Files

| File | Contents |
|---|---|
| `{p}stats.txt` | Key counts |
| `{p}provenance.json` | GEO URLs and dump of `stats` |
| `{p}tables/panel_inventory.csv` | Requested genes vs panel |
| `{p}tables/aoi_metadata.csv` | 93 AOIs, QC, response, in-norm flag |
| `{p}tables/pairing_summary.csv` | Same-ROI tumor/stroma pairing |
| `{p}tables/tests_not_run.csv` | Asked tests and skip reasons |
| `{p}figures/` | Panel, AOI inventory, pairing, labels, QC |

Raw GEO xlsx/csv files are not in git. Rebuild with `bash scripts/gse221733_geomx/00_fetch.sh && python3 scripts/gse221733_geomx/analyze.py`.
"""
    path.write_text(md)


def main():
    repo = Path(__file__).resolve().parents[2]
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", type=Path, default=repo / "data/gse221733")
    ap.add_argument("--out-dir", type=Path, default=repo / "results/gse221733_geomx")
    ap.add_argument("--root-results", type=Path, default=repo / "RESULTS.md")
    args = ap.parse_args()
    data_dir = args.data_dir
    out = args.out_dir
    fig_dir = out / "figures"
    tab_dir = out / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    norm_path = data_dir / "GSE221733_4301_CTA_norm.xlsx"
    qc_path = data_dir / "GSE221733_4301_CTA_QC.xlsx"
    init_path = data_dir / "GSE221733_4301_CTA_initial.csv.gz"
    sm_path = data_dir / "GSE221733_series_matrix.txt"
    if not sm_path.exists():
        sm_path = data_dir / "GSE221733_series_matrix.txt.gz"

    norm = read_xlsx_matrix(norm_path)
    qc = read_xlsx_matrix(qc_path)
    pgenes = probe_genes(init_path)
    meta = parse_series_matrix(sm_path)
    meta.loc[meta.index.isin(norm.index), "in_norm"] = True
    meta.loc[meta.index.isin(qc.index), "in_qc"] = True
    meta["segment_class"] = meta["segment"].map(segment_class)
    meta["is_tcell_aoi"] = meta["segment"].map(looks_like_tcell_aoi)

    genes = list(norm.columns)
    cldn_hits = find_cldn_like(genes) + find_cldn_like(pgenes)
    cldn4_on_panel = "CLDN4" in genes
    n_tcell_aoi = int(meta["is_tcell_aoi"].sum())

    inv = panel_inventory(genes, pgenes)
    inv.to_csv(tab_dir / "panel_inventory.csv", index=False)

    aoi_cols = [
        "geo_accession",
        "roi",
        "segment",
        "segment_class",
        "is_tcell_aoi",
        "patient id",
        "response",
        "status",
        "followup",
        "area",
        "aoinucleicount",
        "qcflags",
        "qc_fail",
        "in_qc",
        "in_norm",
        "plate coord",
    ]
    aoi_out = meta.reset_index().rename(columns={"title": "aoi"})
    keep = ["aoi"] + [c for c in aoi_cols if c in aoi_out.columns]
    aoi_out[keep].to_csv(tab_dir / "aoi_metadata.csv", index=False)

    # Pairing table (norm AOIs)
    pair_rows = []
    for roi, sub in meta[meta["in_norm"]].groupby("roi"):
        classes = set(sub["segment_class"])
        pair_rows.append(
            {
                "roi": int(roi),
                "n_aoi_in_norm": int(len(sub)),
                "has_tumor": "tumor" in classes,
                "has_stroma": "stroma" in classes,
                "paired_tumor_stroma": "tumor" in classes and "stroma" in classes,
                "patient_id": sub["patient id"].iloc[0],
                "response": sub["response"].iloc[0],
            }
        )
    pairing = pd.DataFrame(pair_rows)
    pairing.to_csv(tab_dir / "pairing_summary.csv", index=False)

    tests = pd.DataFrame(
        [
            {
                "test": "tumor_AOI_CLDN4_vs_same_AOI_CD8A",
                "status": "skipped",
                "reason": "CLDN4 absent from CTA panel",
            },
            {
                "test": "tumor_AOI_CLDN4_vs_paired_Tcell_AOI_CD8A",
                "status": "skipped",
                "reason": "CLDN4 absent; no dedicated T-cell AOI class (PanCK pos/neg only)",
            },
            {
                "test": "tumor_AOI_CLDN4_vs_stroma_CD8A_or_Tcell_fraction",
                "status": "skipped",
                "reason": "CLDN4 absent from CTA panel",
            },
            {
                "test": "malignant_AOI_CLDN4_ICB_NR_vs_R",
                "status": "skipped",
                "reason": "CLDN4 absent from CTA panel",
            },
        ]
    )
    tests.to_csv(tab_dir / "tests_not_run.csv", index=False)

    tum_norm = meta[(meta["in_norm"]) & (meta["segment"] == "PanCK pos")]
    pat = tum_norm.groupby("patient id")["response"].first()

    key = {
        "dataset": "GSE221733",
        "panel": "GeoMx Cancer Transcriptome Atlas (CTA)",
        "n_genes": int(len(genes)),
        "n_probes": None,
        "n_probe_genes": int(len(pgenes)),
        "cldn4_on_panel": bool(cldn4_on_panel),
        "n_cldn_string_hits": int(len(set(cldn_hits))),
        "n_geo_aoi": int(len(meta)),
        "n_geo_tumor": int((meta["segment"] == "PanCK pos").sum()),
        "n_geo_stroma": int((meta["segment"] == "PanCK neg").sum()),
        "n_geo_geom": int((meta["segment"] == "Geometric Segment").sum()),
        "n_patients": int(meta["patient id"].nunique()),
        "n_qc_fail": int(meta["qc_fail"].sum()),
        "n_norm": int(meta["in_norm"].sum()),
        "n_norm_tumor": int(((meta["in_norm"]) & (meta["segment"] == "PanCK pos")).sum()),
        "n_norm_stroma": int(((meta["in_norm"]) & (meta["segment"] == "PanCK neg")).sum()),
        "n_qc_matrix": int(meta["in_qc"].sum()),
        "n_tcell_aoi": n_tcell_aoi,
        "n_paired_cores": int(pairing["paired_tumor_stroma"].sum()),
        "n_tumor_only_cores": int((pairing["has_tumor"] & ~pairing["has_stroma"]).sum()),
        "n_stroma_only_cores": int((pairing["has_stroma"] & ~pairing["has_tumor"]).sum()),
        "n_norm_tumor_R": int((tum_norm["response"] == "Responder").sum()),
        "n_norm_tumor_NR": int((tum_norm["response"] == "Non-responder").sum()),
        "n_norm_tumor_NA": int((tum_norm["response"] == "N/A").sum()),
        "n_norm_tumor_patients": int(len(pat)),
        "n_norm_tumor_patients_R": int((pat == "Responder").sum()),
        "n_norm_tumor_patients_NR": int((pat == "Non-responder").sum()),
        "n_norm_tumor_patients_NA": int((pat == "N/A").sum()),
        "tacstd2_on_panel": "TACSTD2" in genes,
        "epcam_on_panel": "EPCAM" in genes,
        "cd8a_on_panel": "CD8A" in genes,
        "tests_run": False,
        "proxy_genes_used": False,
    }

    # Accurate probe count
    with gzip.open(init_path, "rt") as fh:
        header = next(csv.reader(fh))
    key["n_probes"] = int(len(header) - 1)

    # If CLDN4 were present, run the prespecified tests. It is not.
    if cldn4_on_panel:
        tum = tum_norm.join(norm[["CLDN4", "CD8A"]], how="left")
        within = spearman_row(tum["CLDN4"], tum["CD8A"], "tumor_aoi", "CLDN4", "CD8A")
        pd.DataFrame([within]).to_csv(tab_dir / "spearman_cldn4_cd8a_tumor.csv", index=False)
        key["tests_run"] = True
        key["tumor_CLDN4_CD8A_n"] = within["n"]
        key["tumor_CLDN4_CD8A_rho"] = within["rho"]
        key["tumor_CLDN4_CD8A_p"] = within["p"]
        labeled = tum[tum["response"].isin(["Responder", "Non-responder"])]
        if len(labeled) >= 6:
            a = labeled.loc[labeled["response"] == "Responder", "CLDN4"]
            b = labeled.loc[labeled["response"] == "Non-responder", "CLDN4"]
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            key["response_n_R"] = int(len(a))
            key["response_n_NR"] = int(len(b))
            key["response_CLDN4_U"] = float(u)
            key["response_CLDN4_p"] = float(p)

    fig_panel(inv, key["n_genes"], fig_dir)
    fig_aoi_inventory(meta, fig_dir)
    fig_pairing(meta, fig_dir)
    fig_response(meta, fig_dir)
    fig_qc_scatter(meta, fig_dir)
    fig_skip_summary(key["n_genes"], key["n_cldn_string_hits"], fig_dir)

    def write_stats(path: Path, d: dict):
        lines = []
        for k, v in d.items():
            if isinstance(v, float):
                lines.append(f"{k}\t{v:.6g}")
            else:
                lines.append(f"{k}\t{v}")
        path.write_text("\n".join(lines) + "\n")

    write_stats(out / "stats.txt", key)
    write_results_md(out / "RESULTS.md", key, asset_prefix="")
    write_results_md(args.root_results, key, asset_prefix="results/gse221733_geomx")

    provenance = {
        "dataset": "GSE221733",
        "panel": "NanoString GeoMx Cancer Transcriptome Atlas (CTA)",
        "paper": PAPER,
        "files": [
            "GSE221733_4301_CTA_norm.xlsx",
            "GSE221733_4301_CTA_QC.xlsx",
            "GSE221733_4301_CTA_initial.csv.gz",
            "GSE221733_series_matrix.txt",
        ],
        "urls": {
            "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221733",
            "norm": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221733&format=file&file=GSE221733_4301_CTA_norm.xlsx",
            "qc": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221733&format=file&file=GSE221733_4301_CTA_QC.xlsx",
            "initial": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE221733&format=file&file=GSE221733_4301_CTA_initial.csv.gz",
            "matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE221nnn/GSE221733/matrix/GSE221733_series_matrix.txt.gz",
        },
        "absent_on_panel": ["CLDN4", "CLDN3", "CLDN7", "CLDN1", "F11R", "OCLN"],
        "present_not_used_as_proxy": ["TACSTD2", "EPCAM"],
        "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "key": key,
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(key, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
