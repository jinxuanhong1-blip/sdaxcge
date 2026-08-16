#!/usr/bin/env python3
"""GSE271689 GeoMx WTA and E-MTAB-13530 Visium WTA: TACSTD2/CLDN4 vs immune genes."""

from __future__ import annotations

import gzip
import json
import tarfile
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.sparse import csc_matrix
from statsmodels.stats.multitest import multipletests

from config import (
    CLDN4_SYMBOL,
    DATA,
    FIGDIR,
    RESULTS,
    TACSTD2_SYMBOL,
    VISIUM_IMMUNE_GENES,
)

sns.set_theme(style="ticks", context="talk")
plt.rcParams["pdf.fonttype"] = 42

GEOMX_GENES = [TACSTD2_SYMBOL, CLDN4_SYMBOL] + VISIUM_IMMUNE_GENES


def load_pkc_rts_map(pkc_path: Path, genes: list[str]) -> dict[str, str]:
    obj = json.loads(pkc_path.read_text())
    wanted = {g.upper() for g in genes}
    rts_to_gene = {}
    panel_has = {g: False for g in genes}
    for t in obj["Targets"]:
        name = str(t.get("DisplayName") or "")
        if name.upper() in wanted:
            panel_has[name] = True
            for p in t.get("Probes") or []:
                rts = p.get("RTS_ID")
                if rts:
                    rts_to_gene[rts] = name
    return rts_to_gene, panel_has, obj.get("Name"), len(obj["Targets"])


def parse_dcc_counts(path: Path, rts_map: dict[str, str]) -> dict[str, float]:
    counts = {g: 0.0 for g in set(rts_map.values())}
    aligned = raw = None
    opener = gzip.open if path.suffix == ".gz" else open
    in_sum = False
    with opener(path, "rt", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith("Raw,"):
                raw = float(line.split(",", 1)[1])
            elif line.startswith("Aligned,"):
                aligned = float(line.split(",", 1)[1])
            if line == "<Code_Summary>":
                in_sum = True
                continue
            if line == "</Code_Summary>":
                break
            if not in_sum or "," not in line:
                continue
            rts, _, rest = line.partition(",")
            if rts in rts_map:
                try:
                    counts[rts_map[rts]] = float(rest)
                except ValueError:
                    pass
    counts["_aligned"] = aligned if aligned is not None else float("nan")
    counts["_raw"] = raw if raw is not None else float("nan")
    return counts


def parse_series_meta(path: Path) -> pd.DataFrame:
    """Parse GSE series matrix sample characteristics (no OS fields present)."""
    with gzip.open(path, "rt", errors="ignore") as f:
        lines = f.readlines()
    titles = chars = geos = None
    extra = []
    for line in lines:
        if line.startswith("!Sample_title"):
            titles = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_geo_accession"):
            geos = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        elif line.startswith("!Sample_characteristics_ch1"):
            extra.append([x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]])
    n = len(geos)
    df = pd.DataFrame({"geo": geos, "title": titles})
    for row in extra:
        if len(row) != n:
            continue
        keys = []
        vals = []
        for item in row:
            if ": " in item:
                k, v = item.split(": ", 1)
            else:
                k, v = "characteristic", item
            keys.append(k)
            vals.append(v)
        key = pd.Series(keys).mode().iloc[0]
        df[key] = vals
    return df


def geomx_analysis() -> dict:
    pkc = DATA / "geomx_pkc/Hs_R_NGS_WTA_v1.0.pkc"
    rts_map, panel_has, panel_name, n_targets = load_pkc_rts_map(pkc, GEOMX_GENES)
    tar_path = DATA / "gse271689/GSE271689_RAW.tar"
    rows = []
    with tarfile.open(tar_path, "r") as tar:
        for m in tar.getmembers():
            if not m.name.endswith(".dcc.gz"):
                continue
            f = tar.extractfile(m)
            tmp = Path("/tmp") / Path(m.name).name
            tmp.write_bytes(f.read())
            counts = parse_dcc_counts(tmp, rts_map)
            gsm = Path(m.name).name.split("_")[0]
            rec = {"geo": gsm, "dcc": Path(m.name).name}
            rec.update(counts)
            rows.append(rec)
            tmp.unlink(missing_ok=True)
    expr = pd.DataFrame(rows)
    meta = parse_series_meta(DATA / "gse271689/GSE271689_series_matrix.txt.gz")
    df = expr.merge(meta, on="geo", how="left")
    # drop NTC
    if "spotid" in df.columns:
        df = df[~df["spotid"].astype(str).str.contains("No Template", case=False, na=False)]
    if "cell type" in df.columns:
        df = df[df["cell type"].astype(str).str.upper() != "NA"]
    aligned = pd.to_numeric(df["_aligned"], errors="coerce")
    for g in [c for c in df.columns if c in GEOMX_GENES]:
        raw_c = pd.to_numeric(df[g], errors="coerce")
        df[f"{g}_cpm"] = np.where(aligned > 0, raw_c / aligned * 1e6, np.nan)
    df.to_csv(RESULTS / "gse271689_geomx_target_counts.tsv", sep="\t", index=False)

    # compartment comparison
    comp_rows = []
    if "cell type" in df.columns:
        for gene in (TACSTD2_SYMBOL, CLDN4_SYMBOL):
            for ct, sub in df.groupby("cell type"):
                x = pd.to_numeric(sub[gene], errors="coerce")
                xc = pd.to_numeric(sub[f"{gene}_cpm"], errors="coerce")
                comp_rows.append(
                    {
                        "gene": gene,
                        "cell_type": ct,
                        "n_AOI": int(x.notna().sum()),
                        "median_count": float(x.median()),
                        "mean_count": float(x.mean()),
                        "median_cpm": float(xc.median()),
                        "pct_AOI_count_gt0": float((x > 0).mean() * 100),
                    }
                )
    pd.DataFrame(comp_rows).to_csv(RESULTS / "gse271689_compartment_summary.tsv", sep="\t", index=False)

    corr_rows = []
    for ct, sub in df.groupby(df.get("cell type", pd.Series(["all"] * len(df)))):
        for x in (TACSTD2_SYMBOL, CLDN4_SYMBOL):
            for y in VISIUM_IMMUNE_GENES:
                xc, yc = f"{x}_cpm", f"{y}_cpm"
                if xc not in sub.columns or yc not in sub.columns:
                    continue
                pair = sub[[xc, yc]].apply(pd.to_numeric, errors="coerce").dropna()
                if len(pair) < 8:
                    rho = p = np.nan
                else:
                    rho, p = stats.spearmanr(pair[xc], pair[yc])
                corr_rows.append(
                    {
                        "cell_type": ct,
                        "x": x,
                        "y": y,
                        "n": int(len(pair)),
                        "spearman_rho_cpm": rho,
                        "p": p,
                        "normalization": "count_per_million_aligned",
                    }
                )
    corr = pd.DataFrame(corr_rows)
    if not corr.empty and corr["p"].notna().any():
        mask = corr["p"].notna()
        corr.loc[mask, "fdr_bh"] = multipletests(corr.loc[mask, "p"], method="fdr_bh")[1]
    corr.to_csv(RESULTS / "gse271689_gene_immune_correlations.tsv", sep="\t", index=False)

    # figures
    if "cell type" in df.columns:
        fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.6))
        for ax, gene in zip(axes, (TACSTD2_SYMBOL, CLDN4_SYMBOL)):
            col = f"{gene}_cpm"
            plot_df = df[["cell type", col]].copy()
            plot_df[col] = np.log1p(pd.to_numeric(plot_df[col], errors="coerce"))
            sns.boxplot(data=plot_df, x="cell type", y=col, ax=ax, color="#4c6ef5")
            sns.stripplot(data=plot_df, x="cell type", y=col, ax=ax, color="k", size=3, alpha=0.35)
            ax.set_ylabel(f"log1p({gene} CPM)")
            ax.set_title(f"GSE271689 GeoMx WTA {gene}")
        fig.tight_layout()
        fig.savefig(FIGDIR / "gse271689_compartment_box.png")
        fig.savefig(FIGDIR / "gse271689_compartment_box.pdf")
        plt.close(fig)

    os_fields = [
        c
        for c in df.columns
        if any(k in c.lower() for k in ("surviv", "os", "pfs", "dead", "vital"))
    ]
    return {
        "panel_name": panel_name,
        "n_WTA_targets": n_targets,
        "panel_has_TACSTD2": bool(panel_has.get(TACSTD2_SYMBOL)),
        "panel_has_CLDN4": bool(panel_has.get(CLDN4_SYMBOL)),
        "rts_TACSTD2": [k for k, v in rts_map.items() if v == TACSTD2_SYMBOL],
        "rts_CLDN4": [k for k, v in rts_map.items() if v == CLDN4_SYMBOL],
        "n_AOI_after_NTC_filter": int(len(df)),
        "cell_types": sorted(df["cell type"].dropna().unique().tolist()) if "cell type" in df.columns else [],
        "treatment_values": sorted(df["treatment"].dropna().astype(str).unique().tolist())
        if "treatment" in df.columns
        else [],
        "OS_fields_in_GEO": os_fields,
        "OS_available": False,
        "note": "GEO series matrix/SOFT have spotid, tissue, cell type, treatment=immunotherapy only. No OS/PFS columns.",
    }


def read_visium_h5(path: Path, genes: list[str]) -> pd.DataFrame:
    with h5py.File(path, "r") as f:
        names = f["matrix/features/name"][:].astype(str)
        data = f["matrix/data"][:]
        indices = f["matrix/indices"][:]
        indptr = f["matrix/indptr"][:]
        shape = tuple(int(x) for x in f["matrix/shape"][:])
        # shape is (n_genes, n_barcodes)
        mat = csc_matrix((data, indices, indptr), shape=shape)
    name_to_i = {n: i for i, n in enumerate(names)}
    out = {}
    for g in genes:
        if g in name_to_i:
            out[g] = np.asarray(mat[name_to_i[g], :].todense()).ravel()
        else:
            out[g] = np.full(shape[1], np.nan)
    return pd.DataFrame(out)


def section_class(name: str) -> str:
    stem = name.split("-")[0]
    if stem.startswith("D"):
        return "healthy_donor"
    if "_T" in stem:
        return "tumor"
    if "_B" in stem:
        return "adjacent_nontumor"
    return "other"


def visium_analysis() -> dict:
    genes = [TACSTD2_SYMBOL, CLDN4_SYMBOL] + VISIUM_IMMUNE_GENES
    h5s = sorted((DATA / "emtab13530").glob("*-filtered_feature_bc_matrix.h5"))
    sdrf = pd.read_csv(DATA / "emtab13530/E-MTAB-13530.sdrf.txt", sep="\t")
    sdrf_sec = sdrf.drop_duplicates("Source Name")[
        [
            "Source Name",
            "Characteristics[individual]",
            "Characteristics[disease]",
            "Characteristics[sampling site]",
            "Characteristics[sex]",
            "Characteristics[age]",
            "Factor Value[disease]",
            "Factor Value[sampling site]",
        ]
    ].rename(columns={"Source Name": "section"})

    sec_rows = []
    spot_corr = []
    for p in h5s:
        section = p.name.replace("-filtered_feature_bc_matrix.h5", "")
        df = read_visium_h5(p, genes)
        klass = section_class(p.name)
        rec = {
            "section": section,
            "class": klass,
            "n_spots": int(len(df)),
        }
        for g in genes:
            x = pd.to_numeric(df[g], errors="coerce")
            rec[f"{g}_mean"] = float(x.mean())
            rec[f"{g}_detect_pct"] = float((x > 0).mean() * 100)
        sec_rows.append(rec)
        if klass == "tumor":
            for xg in (TACSTD2_SYMBOL, CLDN4_SYMBOL):
                for yg in VISIUM_IMMUNE_GENES:
                    pair = df[[xg, yg]].apply(pd.to_numeric, errors="coerce").dropna()
                    # require some expression
                    if len(pair) < 30:
                        rho = p = np.nan
                    else:
                        rho, p = stats.spearmanr(pair[xg], pair[yg])
                    spot_corr.append(
                        {
                            "section": section,
                            "x": xg,
                            "y": yg,
                            "n_spots": int(len(pair)),
                            "spearman_rho": rho,
                            "p": p,
                        }
                    )
    sec = pd.DataFrame(sec_rows).merge(sdrf_sec, on="section", how="left")
    sec.to_csv(RESULTS / "emtab13530_section_summary.tsv", sep="\t", index=False)
    corr = pd.DataFrame(spot_corr)
    if not corr.empty and corr["p"].notna().any():
        mask = corr["p"].notna()
        corr.loc[mask, "fdr_bh"] = multipletests(corr.loc[mask, "p"], method="fdr_bh")[1]
    corr.to_csv(RESULTS / "emtab13530_tumor_spot_correlations.tsv", sep="\t", index=False)

    # class-level comparison of detection / mean
    class_cmp = []
    for g in (TACSTD2_SYMBOL, CLDN4_SYMBOL):
        for klass, sub in sec.groupby("class"):
            class_cmp.append(
                {
                    "gene": g,
                    "class": klass,
                    "n_sections": int(len(sub)),
                    "median_detect_pct": float(sub[f"{g}_detect_pct"].median()),
                    "median_mean_count": float(sub[f"{g}_mean"].median()),
                }
            )
    pd.DataFrame(class_cmp).to_csv(RESULTS / "emtab13530_class_summary.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8))
    for ax, g in zip(axes, (TACSTD2_SYMBOL, CLDN4_SYMBOL)):
        sns.boxplot(data=sec, x="class", y=f"{g}_detect_pct", ax=ax, color="#2a9d8f")
        sns.stripplot(data=sec, x="class", y=f"{g}_detect_pct", ax=ax, color="k", size=4, alpha=0.5)
        ax.set_ylabel(f"{g} detection (% spots >0)")
        ax.set_title(f"E-MTAB-13530 Visium {g}")
        ax.tick_params(axis="x", rotation=20)
    fig.tight_layout()
    fig.savefig(FIGDIR / "emtab13530_detection_by_class.png")
    fig.savefig(FIGDIR / "emtab13530_detection_by_class.pdf")
    plt.close(fig)

    # pooled tumor-section mean rho
    if not corr.empty:
        fig, ax = plt.subplots(figsize=(8.6, 4.8))
        sub = corr[corr["x"] == TACSTD2_SYMBOL]
        order = VISIUM_IMMUNE_GENES
        sns.boxplot(data=sub, x="y", y="spearman_rho", ax=ax, order=order, color="#457b9d")
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylabel("Per-section Spearman ρ (spots)")
        ax.set_xlabel("Immune gene")
        ax.set_title("E-MTAB-13530 tumor sections: TACSTD2 vs immune genes")
        ax.tick_params(axis="x", rotation=40)
        fig.tight_layout()
        fig.savefig(FIGDIR / "emtab13530_tacstd2_spot_corr.png")
        fig.savefig(FIGDIR / "emtab13530_tacstd2_spot_corr.pdf")
        plt.close(fig)

    return {
        "n_h5_sections": int(len(h5s)),
        "both_genes_in_WTA_matrix": True,
        "n_tumor_sections": int((sec["class"] == "tumor").sum()),
        "n_adjacent_sections": int((sec["class"] == "adjacent_nontumor").sum()),
        "n_healthy_sections": int((sec["class"] == "healthy_donor").sum()),
        "ICI_labels": False,
        "note": "Visium WTA processed Space Ranger h5. Study is NSCLC lesions + adjacent + donors; no ICI labels in SDRF.",
    }


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    geomx = geomx_analysis()
    visium = visium_analysis()
    (RESULTS / "spatial_run_summary.json").write_text(
        json.dumps({"GSE271689": geomx, "E-MTAB-13530": visium}, indent=2)
    )
    print(json.dumps({"GSE271689": geomx, "E-MTAB-13530": visium}, indent=2))


if __name__ == "__main__":
    main()
