#!/usr/bin/env python3
"""GSE131907 LUAD atlas (Kim et al. Nat Commun 2020, PMID 32385277).

Not ICI-response labeled. Used only as an epithelial/malignant TACSTD2–CLDN4
vs T/NK reference. Author cell annotations. Processed raw UMI txt (<2 GB);
the 2.9 GB log2TPM txt is skipped per instructions.

Extracts panel-gene rows only (does not scan all genes for library size).
Metrics: mean log1p(UMI) and percent positive. Sample is the unit for Spearman.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/tmp/scrna_data/GSE131907")
MATRIX = DATA / "raw_UMI_matrix.txt.gz"
ANNOT = DATA / "cell_annotation.txt.gz"
PANEL = Path("/workspace/scripts/scrna/gene_panel.tsv")
RES = Path("/workspace/results/scrna")
RES.mkdir(parents=True, exist_ok=True)

WANT = ["TACSTD2", "CLDN4", "EPCAM", "PTPRC", "CD3D", "CD8A", "NKG7", "MS4A1"]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "mBrain", "PE"}
MIN_EPI = 20
MIN_TNK = 20


def spear(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 4:
        return np.nan, np.nan, int(m.sum())
    r, p = stats.spearmanr(x[m], y[m])
    return float(r), float(p), int(m.sum())


def fmt_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    return f"{p:.2e}" if p < 0.001 else f"{p:.4f}"


def fmt_r(r):
    if r is None or not np.isfinite(r):
        return "NA"
    return f"{r:.3f}"


def main() -> None:
    panel = set(pd.read_csv(PANEL, sep="\t")["gene"]) | set(WANT)
    print("streaming panel rows from GSE131907 raw UMI...", flush=True)
    with gzip.open(MATRIX, "rb") as fh:
        header = fh.readline().decode("ascii").rstrip("\n")
        cells = header.split("\t")[1:]
        n_cells = len(cells)
        print(f"cells {n_cells}", flush=True)
        kept: dict[str, np.ndarray] = {}
        n_genes = 0
        for raw in fh:
            n_genes += 1
            tab = raw.find(b"\t")
            gene = raw[:tab].decode("ascii")
            if gene not in panel:
                continue
            arr = np.fromstring(raw[tab + 1 :].rstrip(b"\r\n"), sep="\t", dtype=np.float32)
            if arr.size != n_cells:
                raise ValueError(f"{gene}: {arr.size} != {n_cells}")
            kept[gene] = arr
            print(f"  kept {gene} ({len(kept)})", flush=True)
    print(f"scanned {n_genes} genes; kept {sorted(kept)}", flush=True)
    missing = sorted(panel - set(kept))
    print("missing from matrix:", missing)

    expr = pd.DataFrame({g: kept[g] for g in sorted(kept)})
    expr.index = cells
    expr.index.name = "Index"

    ann = pd.read_csv(ANNOT, sep="\t")
    df = ann.merge(expr, left_on="Index", right_index=True, how="inner")
    print(f"merged cells {len(df)} / annot {len(ann)}", flush=True)
    if len(df) != len(ann):
        print("WARNING: barcode mismatch", len(ann) - len(df))

    df["is_epi"] = df["Cell_type"] == "Epithelial cells"
    df["is_mal"] = df["Cell_subtype"] == "Malignant cells"
    df["is_tnk"] = df["Cell_type"].isin(["T lymphocytes", "NK cells"])
    df["is_cd8"] = df["Cell_subtype"].isin(
        ["Exhausted CD8+ T", "Naive CD8+ T", "Cytotoxic CD8+ T", "CD8 low T"]
    )
    df["is_b"] = df["Cell_type"] == "B lymphocytes"
    df["is_immune"] = df["Cell_type"].isin(
        ["T lymphocytes", "NK cells", "B lymphocytes", "Myeloid cells", "MAST cells"]
    )
    df["tumor_site"] = df["Sample_Origin"].isin(TUMOR_ORIGINS)

    stats_rows = []

    def add(analysis, comparison, metric, n, stat_name, value, p, note=""):
        stats_rows.append(dict(
            dataset="GSE131907", analysis=analysis, comparison=comparison,
            metric=metric, n=n, stat=stat_name, value=value, p_value=p, note=note,
        ))

    # cell-level compartment (exploratory; huge n)
    for gene in ("TACSTD2", "CLDN4"):
        if gene not in df.columns:
            continue
        logv = np.log1p(df[gene].to_numpy(float))
        df[f"{gene}_log1p"] = logv
        for label, mask_a, mask_b in (
            ("epithelial vs T/NK", df.is_epi, df.is_tnk),
            ("malignant vs T/NK", df.is_mal, df.is_tnk),
            ("epithelial vs T/NK (tLung)", df.is_epi & (df.Sample_Origin == "tLung"),
             df.is_tnk & (df.Sample_Origin == "tLung")),
        ):
            a = logv[mask_a.to_numpy()]
            b = logv[mask_b.to_numpy()]
            if len(a) < 20 or len(b) < 20:
                continue
            u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            add("compartment_cell", label, f"{gene}_log1p_umi",
                f"{len(a)}+{len(b)}", "mannwhitney_u", float(u), float(p),
                f"exploratory cell-level; {gene} epi/mal med={np.median(a):.3f} "
                f"T/NK med={np.median(b):.3f}; %pos {100*(a>0).mean():.1f} vs {100*(b>0).mean():.1f}")
            print(f"  {gene} {label}: n={len(a)}+{len(b)} med={np.median(a):.3f} vs {np.median(b):.3f} "
                  f"%pos {100*(a>0).mean():.1f}/{100*(b>0).mean():.1f} p={fmt_p(p)}")

    # sanity EPCAM / PTPRC
    if "EPCAM" in df.columns and "PTPRC" in df.columns:
        print(f"  sanity EPCAM %pos epi={100*(df.loc[df.is_epi,'EPCAM']>0).mean():.1f} "
              f"T/NK={100*(df.loc[df.is_tnk,'EPCAM']>0).mean():.1f}")
        print(f"  sanity PTPRC %pos epi={100*(df.loc[df.is_epi,'PTPRC']>0).mean():.1f} "
              f"T/NK={100*(df.loc[df.is_tnk,'PTPRC']>0).mean():.1f}")

    # sample-level
    rows = []
    for sample, g in df.groupby("Sample"):
        n = len(g)
        n_epi = int(g.is_epi.sum())
        n_mal = int(g.is_mal.sum())
        n_tnk = int(g.is_tnk.sum())
        n_cd8 = int(g.is_cd8.sum())
        n_b = int(g.is_b.sum())
        rec = dict(
            sample=sample,
            origin=g["Sample_Origin"].iloc[0],
            n_cells=n,
            n_epithelial=n_epi,
            n_malignant=n_mal,
            n_tnk=n_tnk,
            n_cd8=n_cd8,
            n_b=n_b,
            frac_tnk=n_tnk / n,
            frac_cd8=n_cd8 / n,
            frac_b=n_b / n,
            frac_immune=int(g.is_immune.sum()) / n,
        )
        for gene in ("TACSTD2", "CLDN4"):
            if f"{gene}_log1p" not in g.columns:
                continue
            epi = g.loc[g.is_epi]
            mal = g.loc[g.is_mal]
            rec[f"epi_{gene}_mean"] = float(epi[f"{gene}_log1p"].mean()) if n_epi >= MIN_EPI else np.nan
            rec[f"epi_{gene}_pct"] = float((epi[gene] > 0).mean() * 100) if n_epi >= MIN_EPI else np.nan
            rec[f"mal_{gene}_mean"] = float(mal[f"{gene}_log1p"].mean()) if n_mal >= MIN_EPI else np.nan
            rec[f"mal_{gene}_pct"] = float((mal[gene] > 0).mean() * 100) if n_mal >= MIN_EPI else np.nan
        rows.append(rec)
    samp = pd.DataFrame(rows)
    samp.to_csv(RES / "gse131907_sample_table.tsv", sep="\t", index=False)
    print(samp.head(8).to_string(index=False))

    cohorts = {
        "tLung": samp.loc[samp.origin == "tLung"],
        "tumor_sites": samp.loc[samp.origin.isin(TUMOR_ORIGINS)],
        "nLung": samp.loc[samp.origin == "nLung"],
    }
    for cname, sub in cohorts.items():
        for gene in ("TACSTD2", "CLDN4"):
            for expo, ilabel in (
                (f"epi_{gene}_mean", "T/NK fraction"),
                (f"epi_{gene}_mean", "CD8 fraction"),
                (f"epi_{gene}_mean", "B fraction"),
            ):
                immune = {"T/NK fraction": "frac_tnk", "CD8 fraction": "frac_cd8",
                          "B fraction": "frac_b"}[ilabel]
                ok = sub.dropna(subset=[expo, immune])
                r, p, n = spear(ok[expo], ok[immune])
                add("spearman_sample", f"{cname} epithelial {gene} vs {ilabel}",
                    expo, n, "spearman_rho", r, p,
                    "not ICI-labeled; author annotations; log1p UMI; >=20 epi cells")
                print(f"  {cname} epi {gene} vs {ilabel}: n={n} rho={fmt_r(r)} p={fmt_p(p)}")
        # malignant-only (mostly mets)
        for gene in ("TACSTD2", "CLDN4"):
            ok = sub.dropna(subset=[f"mal_{gene}_mean", "frac_tnk"])
            r, p, n = spear(ok[f"mal_{gene}_mean"], ok["frac_tnk"])
            add("spearman_sample", f"{cname} malignant {gene} vs T/NK fraction",
                f"mal_{gene}_mean", n, "spearman_rho", r, p,
                "author Malignant cells subtype; mostly mets/PE")
            print(f"  {cname} mal {gene} vs T/NK: n={n} rho={fmt_r(r)} p={fmt_p(p)}")

    # paired epithelial vs T/NK TACSTD2 at sample level (tumor sites)
    pair_rows = []
    for sample, g in df.loc[df.tumor_site].groupby("Sample"):
        epi = g.loc[g.is_epi, "TACSTD2_log1p"] if "TACSTD2_log1p" in g.columns else pd.Series(dtype=float)
        tnk = g.loc[g.is_tnk, "TACSTD2_log1p"] if "TACSTD2_log1p" in g.columns else pd.Series(dtype=float)
        if len(epi) >= MIN_EPI and len(tnk) >= MIN_TNK:
            pair_rows.append((float(epi.mean()), float(tnk.mean())))
    if len(pair_rows) >= 4:
        a, b = zip(*pair_rows)
        w, p = stats.wilcoxon(a, b)
        add("compartment_paired", "tumor-site epithelial vs T/NK TACSTD2",
            "epi_TACSTD2_mean", len(pair_rows), "wilcoxon_w", float(w), float(p),
            f"mal/epi med={np.median(a):.3f}; T/NK med={np.median(b):.3f}")
        print(f"  paired tumor-site epi vs T/NK TACSTD2 n={len(pair_rows)} W={w} p={fmt_p(p)}")

    # figures
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))
    for ax, gene in zip(axes, ("TACSTD2", "CLDN4")):
        col = f"{gene}_log1p"
        data = [
            df.loc[df.is_epi, col].to_numpy(),
            df.loc[df.is_tnk, col].to_numpy(),
        ]
        ax.boxplot(data, tick_labels=["epithelial", "T/NK"], showfliers=False)
        ax.set_ylabel(f"{gene} log1p UMI")
        ax.set_title(f"{gene} (cell-level, exploratory)")
    fig.suptitle("GSE131907 atlas: TACSTD2/CLDN4 are epithelial-restricted (no ICI labels)", fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "gse131907_compartment.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 4.0))
    tlung = samp.loc[samp.origin == "tLung"].dropna(subset=["epi_TACSTD2_mean"])
    tumor = samp.loc[samp.origin.isin(TUMOR_ORIGINS)].dropna(subset=["epi_TACSTD2_mean"])
    for ax, sub, title in (
        (axes[0], tlung, "tLung"),
        (axes[1], tumor, "tumor sites"),
    ):
        ax.scatter(sub["epi_TACSTD2_mean"], sub["frac_tnk"], s=36, c="#4C78A8", edgecolor="k", linewidth=0.4)
        r, p, n = spear(sub["epi_TACSTD2_mean"], sub["frac_tnk"])
        ax.set_title(f"{title} n={n}  ρ={fmt_r(r)}  p={fmt_p(p)}")
        ax.set_xlabel("epithelial TACSTD2 (mean log1p UMI)")
        ax.set_ylabel("T/NK fraction")
    fig.suptitle("GSE131907: epithelial TACSTD2 vs T/NK (not ICI-labeled)", fontsize=10)
    fig.tight_layout()
    fig.savefig(RES / "gse131907_epi_tacstd2_vs_tnk.png", dpi=160)
    plt.close(fig)

    # cell-type %pos bar
    order = ["Epithelial cells", "T lymphocytes", "NK cells", "B lymphocytes",
             "Myeloid cells", "MAST cells", "Fibroblasts", "Endothelial cells"]
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    pcts = [100 * (df.loc[df.Cell_type == ct, "TACSTD2"] > 0).mean() if (df.Cell_type == ct).any() else 0
            for ct in order]
    ax.bar(range(len(order)), pcts, color="#4C78A8")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=35, ha="right")
    ax.set_ylabel("TACSTD2 % positive (UMI>0)")
    ax.set_title("GSE131907: TACSTD2 detection by author cell type")
    fig.tight_layout()
    fig.savefig(RES / "gse131907_tacstd2_by_celltype.png", dpi=160)
    plt.close(fig)

    new_stats = pd.DataFrame(stats_rows)
    prev = RES / "stats.tsv"
    if prev.exists():
        old = pd.read_csv(prev, sep="\t")
        old = old.loc[old.dataset != "GSE131907"]
        out = pd.concat([old, new_stats], ignore_index=True)
    else:
        out = new_stats
    out.to_csv(prev, sep="\t", index=False)
    print(f"wrote {prev} rows={len(out)}")
    print(new_stats.to_string(index=False))


if __name__ == "__main__":
    main()
