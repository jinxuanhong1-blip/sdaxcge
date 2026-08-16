#!/usr/bin/env python3
"""Hunt: Tacstd2 (Trop2) in GSE76628 -- public data only.

Requested contrast: "Tacstd2-high tumor subset vs CD8/NK".

Reality check (documented in results/hunt_gse76628/README.md): GSE76628 is a
BULK Affymetrix Mouse 430 2.0 (GPL1261) series of 78 whole-tissue samples from
the Ad-VEGF-A164 tumor-surrogate angiogenesis model (mouse ear/flank), with
anti-VEGFR2 (DC101) / anti-VEGF (G6-31) treatment arms. It contains no tumor
cells, no sorted or single-cell populations, and no KL (Kras;Lkb1) mice. The
requested cell-level contrast is therefore infeasible in this dataset. This
script runs the closest honest, feasible analyses instead:

  1. Tacstd2 expression by experimental group (is it expressed at all, and
     does it change with the angiogenic response / treatment?).
  2. Bulk-level Spearman correlation of Tacstd2 with CD8 T and NK cell marker
     scores across all 78 samples. This is a tissue-level co-abundance proxy,
     NOT a tumor-cell-vs-lymphocyte expression comparison.
  3. Expression-above-background assessment (within-array percentile ranks).

Inputs are downloaded from NCBI GEO FTP if not already cached in data/:
  - GSE76628_series_matrix.txt.gz (processed matrix, linear-scale intensities)
  - GPL1261.annot.gz (probe -> gene symbol annotation)

Outputs: results/hunt_gse76628/ (CSV tables, PNG figures, provenance.json).
"""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
OUT = REPO / "results" / "hunt_gse76628"

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76628/matrix/"
    "GSE76628_series_matrix.txt.gz"
)
ANNOT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL1nnn/GPL1261/annot/"
    "GPL1261.annot.gz"
)

TARGET = "Tacstd2"
# Marker panels (mouse symbols). Bulk-tissue proxies only.
PANELS = {
    "CD8_T": ["Cd8a", "Cd8b1", "Cd3e", "Cd3d", "Cd3g"],
    "NK": ["Ncr1", "Klrb1c", "Klrk1", "Eomes", "Prf1", "Gzmb", "Nkg7", "Il2rb"],
    "Epithelial": ["Epcam", "Krt8", "Krt18", "Krt19"],
}

GROUP_ORDER = [
    "Normal",
    "5_Day_NT", "5_Day_DC101", "5_Day_G6",
    "20_Day_NT", "20_Day_DC101", "20_Day_G6",
    "60_Day_NT", "60_Day_DC101", "60_Day_G6",
]


def fetch(url: str, dest: Path) -> dict:
    if not dest.exists():
        tmp = dest.with_suffix(dest.suffix + ".part")
        with urllib.request.urlopen(url) as r, open(tmp, "wb") as f:
            shutil.copyfileobj(r, f)
        tmp.rename(dest)
    md5 = hashlib.md5(dest.read_bytes()).hexdigest()
    return {"url": url, "path": str(dest.relative_to(REPO)), "md5": md5,
            "bytes": dest.stat().st_size}


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (expression probes x samples, sample metadata)."""
    meta_lines: dict[str, list[list[str]]] = {}
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                in_table = False
                continue
            if in_table:
                table_lines.append(line)
            elif line.startswith("!Sample_"):
                key = line.split("\t", 1)[0][len("!Sample_"):]
                vals = [v.strip('"') for v in line.split("\t")[1:]]
                meta_lines.setdefault(key, []).append(vals)

    header = [c.strip('"') for c in table_lines[0].split("\t")]
    rows = [ln.split("\t") for ln in table_lines[1:] if ln]
    expr = pd.DataFrame(
        [r[1:] for r in rows],
        index=[r[0].strip('"') for r in rows],
        columns=header[1:],
        dtype=float,
    )

    meta = pd.DataFrame({
        "gsm": meta_lines["geo_accession"][0],
        "title": meta_lines["title"][0],
    }).set_index("gsm")
    return expr, meta


def assign_group(title: str) -> str:
    # Titles look like can_Ad_VEGF_20_Day_DC101_072 or can_Ad_VEGF_Normal_097
    core = title.replace("can_Ad_VEGF_", "")
    parts = core.split("_")
    if parts[0] == "Normal":
        return "Normal"
    return "_".join(parts[:-1])  # drop trailing animal number


def parse_annot(path: Path) -> pd.Series:
    """Return probe ID -> gene symbol."""
    ids, syms = [], []
    in_table = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                i_id = header.index("ID")
                i_sym = header.index("Gene symbol")
                continue
            if line.startswith("!platform_table_end"):
                break
            if in_table:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(i_id, i_sym) and f[i_sym]:
                    ids.append(f[i_id])
                    syms.append(f[i_sym])
    return pd.Series(syms, index=ids, name="symbol")


def main() -> None:
    DATA.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    prov_files = [
        fetch(MATRIX_URL, DATA / "GSE76628_series_matrix.txt.gz"),
        fetch(ANNOT_URL, DATA / "GPL1261.annot.gz"),
    ]

    expr_lin, meta = parse_series_matrix(DATA / "GSE76628_series_matrix.txt.gz")
    meta["group"] = meta["title"].map(assign_group)
    assert set(meta["group"]) == set(GROUP_ORDER), sorted(set(meta["group"]))

    # Series matrix stores linear-scale intensities (MAS5-like); log2 for stats.
    expr = np.log2(expr_lin.clip(lower=1.0) + 1.0)

    sym = parse_annot(DATA / "GPL1261.annot.gz")
    wanted = [TARGET] + [g for panel in PANELS.values() for g in panel]
    probes = sym[sym.isin(wanted)]

    # Per-array percentile rank of each probe (linear scale ranks are the same
    # as log ranks); used for "expressed above background?" assessment.
    pct = expr.rank(axis=0, pct=True) * 100.0

    # ---- per-sample marker table (all probes for genes of interest) ----
    long_rows = []
    for probe, gene in probes.items():
        for gsm in expr.columns:
            long_rows.append({
                "gene": gene, "probe": probe, "gsm": gsm,
                "group": meta.loc[gsm, "group"],
                "log2_intensity": round(float(expr.loc[probe, gsm]), 4),
                "within_array_percentile": round(float(pct.loc[probe, gsm]), 2),
            })
    per_sample = pd.DataFrame(long_rows)
    per_sample.to_csv(OUT / "marker_expression_per_sample.csv", index=False)

    # ---- gene-level values: for each gene use the probe with highest mean ----
    def best_probe(gene: str) -> str:
        cand = probes[probes == gene].index
        return expr.loc[cand].mean(axis=1).idxmax()

    genes_found = sorted(set(probes.values))
    gene_expr = pd.DataFrame(
        {g: expr.loc[best_probe(g)] for g in genes_found}
    )
    gene_expr["group"] = meta.loc[gene_expr.index, "group"].values

    missing = sorted(set(wanted) - set(genes_found))

    # ---- Tacstd2 group stats (all Tacstd2 probes reported) ----
    tac_probes = probes[probes == TARGET].index.tolist()
    grp_rows = []
    for probe in tac_probes:
        vals = expr.loc[probe]
        for grp in GROUP_ORDER:
            gsms = meta.index[meta["group"] == grp]
            v = vals[gsms]
            grp_rows.append({
                "probe": probe, "group": grp, "n": len(v),
                "mean_log2": round(v.mean(), 3), "sd_log2": round(v.std(), 3),
                "median_log2": round(v.median(), 3),
                "mean_within_array_percentile": round(
                    float(pct.loc[probe, gsms].mean()), 1),
            })
    grp_stats = pd.DataFrame(grp_rows)

    kw_rows = []
    for probe in tac_probes:
        groups = [expr.loc[probe, meta.index[meta["group"] == g]].values
                  for g in GROUP_ORDER]
        kw = stats.kruskal(*groups)
        kw_rows.append({"probe": probe, "kruskal_H": round(kw.statistic, 3),
                        "kruskal_p": kw.pvalue})
    grp_stats = grp_stats.merge(pd.DataFrame(kw_rows), on="probe")
    grp_stats.to_csv(OUT / "tacstd2_group_stats.csv", index=False)

    # ---- bulk correlations: Tacstd2 vs individual markers and panel scores ----
    tac = gene_expr[TARGET]
    corr_rows = []
    for panel, genes in PANELS.items():
        present = [g for g in genes if g in gene_expr.columns]
        for g in present:
            rho, p = stats.spearmanr(tac, gene_expr[g])
            corr_rows.append({"comparator": g, "panel": panel,
                              "type": "single_gene",
                              "spearman_rho": round(rho, 3), "p_value": p})
        if present:
            # panel score: mean of z-scored member genes
            z = (gene_expr[present] - gene_expr[present].mean()) / gene_expr[present].std()
            score = z.mean(axis=1)
            rho, p = stats.spearmanr(tac, score)
            corr_rows.append({"comparator": f"{panel}_score", "panel": panel,
                              "type": "panel_score",
                              "spearman_rho": round(rho, 3), "p_value": p})
    corr = pd.DataFrame(corr_rows)
    corr["note"] = ("bulk tissue co-abundance across 78 samples; "
                    "NOT a cell-level tumor-vs-lymphocyte comparison")
    corr.to_csv(OUT / "tacstd2_vs_immune_correlations.csv", index=False)

    # ---- figures ----
    fig, axes = plt.subplots(len(tac_probes), 1,
                             figsize=(10, 3.2 * len(tac_probes)), squeeze=False)
    for ax, probe in zip(axes.ravel(), tac_probes):
        data = [expr.loc[probe, meta.index[meta["group"] == g]].values
                for g in GROUP_ORDER]
        ax.boxplot(data, tick_labels=GROUP_ORDER)
        med = float(expr.median().median())
        ax.axhline(med, color="grey", ls="--", lw=0.8,
                   label=f"array-wide median (~{med:.1f})")
        ax.set_ylabel("log2 intensity")
        ax.set_title(f"Tacstd2 ({probe}) by group -- GSE76628 bulk arrays")
        ax.tick_params(axis="x", rotation=30)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "tacstd2_by_group.png", dpi=150)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = {g: c for g, c in zip(
        GROUP_ORDER, plt.cm.tab10(np.linspace(0, 1, len(GROUP_ORDER))))}
    for ax, panel in zip(axes, ["CD8_T", "NK"]):
        present = [g for g in PANELS[panel] if g in gene_expr.columns]
        z = (gene_expr[present] - gene_expr[present].mean()) / gene_expr[present].std()
        score = z.mean(axis=1)
        for grp in GROUP_ORDER:
            m = gene_expr["group"] == grp
            ax.scatter(score[m], tac[m], s=18, color=colors[grp], label=grp)
        rho, p = stats.spearmanr(tac, score)
        ax.set_xlabel(f"{panel} marker score (mean z of {', '.join(present)})")
        ax.set_ylabel(f"{TARGET} log2 intensity (best probe)")
        ax.set_title(f"{TARGET} vs {panel} score\n"
                     f"Spearman rho={rho:.2f}, p={p:.1e} (bulk proxy)")
    axes[1].legend(fontsize=6, loc="best", ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / "tacstd2_vs_cd8_nk_scores.png", dpi=150)
    plt.close(fig)

    # ---- provenance ----
    prov = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "accession": "GSE76628",
        "platform": "GPL1261 (Affymetrix Mouse Genome 430 2.0, bulk array)",
        "n_samples": int(expr.shape[1]),
        "n_probes": int(expr.shape[0]),
        "data_scope": "public GEO files only; no private/KL-mouse data used or available",
        "processing": "series-matrix linear intensities, log2(x+1); "
                      "gene level = highest-mean probe per gene",
        "genes_requested": wanted,
        "genes_not_on_array_annotation": missing,
        "files": prov_files,
    }
    (OUT / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")

    print("groups:", meta["group"].value_counts().to_dict())
    print("Tacstd2 probes:", tac_probes)
    print("missing genes:", missing)
    print(grp_stats.to_string(index=False))
    print(corr.drop(columns="note").to_string(index=False))


if __name__ == "__main__":
    main()
