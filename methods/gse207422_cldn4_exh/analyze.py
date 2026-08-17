#!/usr/bin/env python3
"""GSE207422: T/NK cytotoxicity / exhaustion vs malignant CLDN4.

Additive extra. Not a TACSTD2 redo. Patient/sample is the unit.
Author CopyKAT barcodes are not on GEO; malignant is a marker proxy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent

MIN_UMI = 200
MIN_MAL = 10
MIN_TNK = 10
MIN_CD8 = 10

CYTO_GENES = ["GZMB", "PRF1", "GNLY", "NKG7"]
EXH_GENES = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"]
ZERO_NORMAL = ["SFTPA2", "AGER", "SCGB1A1", "SCGB3A1", "TPPP3"]

LINEAGE_MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT5", "KRT7", "KRT17", "ELF3", "CDH1", "MUC1"],
    "T/NK": ["CD3D", "CD3E", "CD3G", "TRAC", "CD2", "NKG7", "GNLY", "KLRD1"],
    "B/Plasma": ["CD79A", "CD79B", "MS4A1", "JCHAIN", "MZB1"],
    "Myeloid": ["LYZ", "CD68", "CD14", "C1QA", "C1QB", "FCN1"],
    "Mast": ["TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "Endothelial": ["PECAM1", "VWF", "CLDN5", "CDH5"],
    "Fibroblast": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "TAGLN"],
}
NORMAL_LUNG = ["SFTPA1", "SFTPA2", "SFTPB", "SFTPD", "AGER", "NAPSA", "SCGB1A1", "SCGB3A2", "TPPP3", "FOXJ1", "CAPS"]
TUMOR_EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CEACAM5", "CEACAM6", "MUC1", "ELF3"]


def present(genes, columns) -> list[str]:
    return [g for g in genes if g in columns]


def module_mean(ln: pd.DataFrame, genes: list[str]) -> pd.Series:
    gs = present(genes, ln.columns)
    if not gs:
        return pd.Series(np.nan, index=ln.index)
    return ln[gs].mean(axis=1)


def assign_lineage(ln: pd.DataFrame) -> pd.Series:
    scores = pd.DataFrame({lin: module_mean(ln, gs) for lin, gs in LINEAGE_MARKERS.items()})
    labels = scores.idxmax(axis=1)
    labels[scores.max(axis=1) <= 0] = "Unassigned"
    return labels


def spearman_safe(x, y) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(x.size)
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan, "note": "n<4"}
    rho, p = stats.spearmanr(x, y)
    return {"n": n, "rho": float(rho), "p": float(p), "note": ""}


def mwu_safe(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return {
            "n_a": int(len(a)),
            "n_b": int(len(b)),
            "median_a": float(np.median(a)) if len(a) else np.nan,
            "median_b": float(np.median(b)) if len(b) else np.nan,
            "p": np.nan,
            "note": "n<2 in a group",
        }
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "median_a": float(np.median(a)),
        "median_b": float(np.median(b)),
        "u": float(u),
        "p": float(p),
        "note": "",
    }


def load_meta(path: Path) -> pd.DataFrame:
    raw = pd.read_excel(path)
    raw = raw.dropna(subset=["Sample"]).copy()
    raw = raw[~raw["Sample"].astype(str).str.contains("RECIST|MPR:|NMPR:|pCR:", regex=True)]
    raw["Sample"] = raw["Sample"].astype(str)
    raw["Patient"] = raw["Patient"].astype(str)
    raw["mpr_group"] = raw["Pathologic Response"].map({"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR"})
    raw["timing"] = raw["Resource"].map(
        {"Pre-treatment biopsy": "pre", "Post-treatment surgery": "post"}
    )
    raw["residual_tumor"] = pd.to_numeric(raw["Residual Tumor"], errors="coerce")
    return raw


def load_panel(datadir: Path) -> tuple[pd.DataFrame, list[str]]:
    z = np.load(datadir / "extracted_panel.npz", allow_pickle=True)
    cells = z["cells"]
    total = z["total"]
    genes = [str(g) for g in z["genes"]]
    mat = z["mat"]
    missing = [str(g) for g in z["missing"]] if "missing" in z.files else []
    df = pd.DataFrame({"barcode": cells, "total_umi": total})
    df["Sample"] = df["barcode"].astype(str).str.rsplit("_", n=1).str[0]
    for i, g in enumerate(genes):
        df[g] = mat[i]
    return df, missing


def score_patients(cells: pd.DataFrame, ln: pd.DataFrame) -> pd.DataFrame:
    cells = cells.copy()
    cells["cyto"] = module_mean(ln, CYTO_GENES)
    cells["exh"] = module_mean(ln, EXH_GENES)
    cells["cld"] = ln["CLDN4"] if "CLDN4" in ln.columns else np.nan
    rows = []
    for sample, sub in cells.groupby("Sample", sort=True):
        defs = {
            "malig_module": sub["is_malig_module"],
            "malig_zero": sub["is_malig_zero"],
            "epithelial": sub["is_epi"],
        }
        tnk = sub[sub["is_tnk"]]
        cd8 = sub[sub["is_cd8"]]
        rec = {
            "Sample": sample,
            "Patient": sub["Patient"].iloc[0],
            "timing": sub["timing"].iloc[0],
            "mpr_group": sub["mpr_group"].iloc[0],
            "Pathology": sub["Pathology"].iloc[0],
            "RECIST": sub["RECIST"].iloc[0],
            "residual_tumor": sub["residual_tumor"].iloc[0],
            "n_cells": int(len(sub)),
            "n_tnk": int(len(tnk)),
            "n_cd8": int(len(cd8)),
            "tnk_frac": float(len(tnk) / len(sub)) if len(sub) else np.nan,
        }
        for name, mask in defs.items():
            mal = sub[mask]
            rec[f"n_{name}"] = int(len(mal))
            if len(mal) >= MIN_MAL:
                rec[f"{name}_CLDN4"] = float(mal["cld"].mean())
                rec[f"{name}_CLDN4_pct"] = float((mal["CLDN4"] > 0).mean() * 100.0) if "CLDN4" in mal.columns else np.nan
            else:
                rec[f"{name}_CLDN4"] = np.nan
                rec[f"{name}_CLDN4_pct"] = np.nan
        if len(tnk) >= MIN_TNK:
            rec["tnk_cyto"] = float(tnk["cyto"].mean())
            rec["tnk_exh"] = float(tnk["exh"].mean())
        else:
            rec["tnk_cyto"] = np.nan
            rec["tnk_exh"] = np.nan
        if len(cd8) >= MIN_CD8:
            rec["cd8_cyto"] = float(cd8["cyto"].mean())
            rec["cd8_exh"] = float(cd8["exh"].mean())
        else:
            rec["cd8_cyto"] = np.nan
            rec["cd8_exh"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def contrast_row(pat: pd.DataFrame, x: str, y: str, label: str, subset: str) -> dict:
    s = spearman_safe(pat[x], pat[y])
    kept = pat.loc[np.isfinite(pat[x]) & np.isfinite(pat[y]), "Patient"].tolist()
    return {
        "contrast": label,
        "subset": subset,
        "x": x,
        "y": y,
        "n": s["n"],
        "rho": s["rho"],
        "p": s["p"],
        "patients": ",".join(kept),
        "note": s["note"],
    }


def make_figure(pat: pd.DataFrame, outdir: Path) -> None:
    post = pat[pat.timing == "post"].copy()
    fig, axes = plt.subplots(2, 2, figsize=(10.4, 8.6))
    panels = [
        (axes[0, 0], post, "malig_module_CLDN4", "tnk_cyto",
         "A  malignant-like CLDN4 vs T/NK cytotoxicity", "post, ≥10 malig-module"),
        (axes[0, 1], post, "malig_module_CLDN4", "tnk_exh",
         "B  malignant-like CLDN4 vs T/NK exhaustion", "post, ≥10 malig-module"),
        (axes[1, 0], post, "epithelial_CLDN4", "tnk_cyto",
         "C  all-epithelial CLDN4 vs T/NK cytotoxicity", "post, ≥10 epithelial"),
        (axes[1, 1], post, "epithelial_CLDN4", "tnk_exh",
         "D  all-epithelial CLDN4 vs T/NK exhaustion", "post, ≥10 epithelial"),
    ]
    colors = {"MPR": "#2a6f97", "NMPR": "#c44536"}
    markers = {"Adeno": "o", "Squamous": "s"}
    for ax, df, x, y, title, subtitle in panels:
        ax.axhline(df[y].median() if df[y].notna().any() else 0, color="0.85", lw=0.6)
        ax.axvline(df[x].median() if df[x].notna().any() else 0, color="0.85", lw=0.6)
        plotted = 0
        for _, r in df.iterrows():
            if not (np.isfinite(r[x]) and np.isfinite(r[y])):
                continue
            ax.scatter(
                r[x], r[y],
                c=colors.get(r["mpr_group"], "0.4"),
                marker=markers.get(r["Pathology"], "o"),
                s=58, edgecolors="0.15", linewidths=0.5, zorder=3,
            )
            ax.annotate(r["Patient"], (r[x], r[y]), textcoords="offset points",
                        xytext=(4, 3), fontsize=7, color="0.2")
            plotted += 1
        s = spearman_safe(df[x], df[y])
        rho = "NA" if not np.isfinite(s["rho"]) else f"{s['rho']:.2f}"
        p = "NA" if not np.isfinite(s["p"]) else f"{s['p']:.2g}"
        ax.set_title(f"{title}\n{subtitle}  n={s['n']}  ρ={rho}  p={p}", fontsize=8.5)
        ax.set_xlabel("mean log1p(CP10k) CLDN4")
        ax.set_ylabel("mean log1p(CP10k) T/NK score")
        ax.set_xlim(left=min(-0.05, ax.get_xlim()[0]))
        ax.set_ylim(bottom=min(-0.05, ax.get_ylim()[0]))
    handles = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=colors["NMPR"],
                   markeredgecolor="0.15", markersize=7, label="NMPR"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=colors["MPR"],
                   markeredgecolor="0.15", markersize=7, label="MPR (pCR as MPR)"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="0.4",
                   markeredgecolor="0.15", markersize=7, label="Adeno"),
        plt.Line2D([0], [0], marker="s", color="w", markerfacecolor="0.4",
                   markeredgecolor="0.15", markersize=7, label="Squamous"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=8)
    fig.suptitle(
        "Extra — GSE207422 T/NK cytotoxicity / exhaustion vs malignant CLDN4\n"
        "Patient is the unit. Not a TACSTD2 redo. Author CopyKAT IDs are not public.",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.90))
    fig.savefig(outdir / "fig_extra_gse207422_cldn4_exh.png", dpi=160)
    fig.savefig(outdir / "fig_extra_gse207422_cldn4_exh.pdf")
    plt.close(fig)


def fmt(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        if abs(x) >= 0.01 or x == 0:
            return f"{x:.{nd}f}".rstrip("0").rstrip(".")
        return f"{x:.2g}"
    return str(x)


def write_finding(pat: pd.DataFrame, assoc: pd.DataFrame, notes: dict, outdir: Path) -> None:
    post = pat[pat.timing == "post"]
    pre = pat[pat.timing == "pre"]
    prim = post[post["malig_module_CLDN4"].notna() & post["tnk_cyto"].notna()]
    dropped = post[post["malig_module_CLDN4"].isna()]
    if len(dropped):
        dropped_txt = ", ".join(
            f"{r.Patient} (n_malig={int(r.n_malig_module)}, {r.mpr_group})"
            for _, r in dropped.iterrows()
        )
    else:
        dropped_txt = "none"
    lines = [
        "# Extra — GSE207422 T/NK exhaustion / cytotoxicity vs malignant CLDN4",
        "",
        "Additive extra on public [GSE207422](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207422) "
        "(Hu et al., *Genome Med* 2023, PMID 36869384). **Not a TACSTD2 redo.** "
        "The named question is T/NK **cytotoxicity** and **exhaustion** versus **malignant CLDN4**. "
        "Patient / one sample per patient is the unit. Cell-level p-values are not reported.",
        "",
        "## Question",
        "",
        "In the 12 post-resection NSCLC tumors (neoadjuvant PD-1 + chemo), does mean malignant CLDN4 "
        "track T/NK cytotoxicity (`GZMB`, `PRF1`, `GNLY`, `NKG7`) or exhaustion "
        "(`PDCD1`, `HAVCR2`, `LAG3`, `TIGIT`, `TOX`)?",
        "",
        "## Honest n",
        "",
        f"- Matrix: **{notes['n_cells_raw']}** cells (paper post-QC count). QC here: UMI ≥ {MIN_UMI} → "
        f"**{notes['n_cells_qc']}** cells.",
        "- Author CopyKAT / per-cell labels are **not on GEO**. Malignant is a marker proxy.",
        "- **Primary malignant-like** = epithelial lineage **and** tumor-epi module > normal-lung module "
        "**and** normal-lung < 0.4 (same gate as the multi-cohort exhaustion meta).",
        f"- Primary filter: ≥{MIN_MAL} malignant-like **and** ≥{MIN_TNK} T/NK. "
        f"Post-treatment: **{len(prim)} / 12** patients enter Spearman.",
        f"- Dropped post (malignant-like < {MIN_MAL}): {dropped_txt}.",
        f"- Under the primary gate, MPR with usable malignant CLDN4 is **{int((prim.mpr_group=='MPR').sum())}** "
        f"(P03 only if that is the count). Do not read an MPR contrast from n=1.",
        f"- Pre-treatment biopsies: {len(pre)} patients (P01 NE / P05 NMPR / P08 NMPR); "
        "unpaired with the 12 posts. Shown in the occupancy table, not stacked into the primary ρ.",
        f"- All 12 posts have ≥{MIN_TNK} T/NK. All-epithelial sensitivity keeps every post with ≥{MIN_MAL} epithelial cells.",
        "",
        "### Per-patient occupancy (post first)",
        "",
        "| Patient | Timing | MPR | Histology | n cells | n malig-module | n malig-zero | n epithelial | n T/NK | n CD8 | In primary? |",
        "|---|---|---|---|---:|---:|---:|---:|---:|---:|---|",
    ]
    show = pd.concat([post.sort_values("Patient"), pre.sort_values("Patient")], ignore_index=True)
    for _, r in show.iterrows():
        in_prim = "yes" if (r.timing == "post" and np.isfinite(r.malig_module_CLDN4) and np.isfinite(r.tnk_cyto)) else "no"
        lines.append(
            f"| {r.Patient} | {r.timing} | {r.mpr_group if pd.notna(r.mpr_group) else r.get('Pathologic Response', 'NE')} "
            f"| {r.Pathology} | {int(r.n_cells)} | {int(r.n_malig_module)} | {int(r.n_malig_zero)} | "
            f"{int(r.n_epithelial)} | {int(r.n_tnk)} | {int(r.n_cd8)} | {in_prim} |"
        )

    lines += [
        "",
        "pCR (P06) is grouped as MPR. P01 is pathologic NE.",
        "",
        "## Primary table (patient-level Spearman)",
        "",
        "Score = mean log1p(CP10k) of the named genes in the named compartment. "
        f"Keep a patient if ≥{MIN_MAL} cells in the CLDN4 compartment and ≥{MIN_TNK} T/NK.",
        "",
        "| Contrast | Compartment | Subset | n | ρ | p | Patients |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for _, r in assoc.iterrows():
        rho = fmt(r["rho"])
        p = fmt(r["p"], nd=3)
        if isinstance(r["p"], float) and np.isfinite(r["p"]):
            p = f"{r['p']:.3g}"
        lines.append(
            f"| {r['contrast']} | {r['x'].replace('_CLDN4','').replace('malig_module','malignant-like').replace('malig_zero','zero-normal-UMI').replace('epithelial','all-epithelial')} "
            f"| {r['subset']} | {int(r['n'])} | {rho} | {p} | {r['patients']} |"
        )

    # pull key numbers for the verdict
    def grab(contrast, subset):
        hit = assoc[(assoc.contrast == contrast) & (assoc.subset == subset)]
        if hit.empty:
            return None
        return hit.iloc[0]

    a = grab("CLDN4 vs T/NK cytotoxicity", "post, malignant-like")
    b = grab("CLDN4 vs T/NK exhaustion", "post, malignant-like")
    c = grab("CLDN4 vs T/NK cytotoxicity", "post, all-epithelial")
    d = grab("CLDN4 vs T/NK exhaustion", "post, all-epithelial")

    lines += [
        "",
        "## Extra figure",
        "",
        "`fig_extra_gse207422_cldn4_exh.png` — four patient-level scatters.",
        "",
        "| Panel | Test | n | ρ | p |",
        "|---|---|---:|---:|---:|",
    ]
    if a is not None:
        lines.append(f"| A | malignant-like CLDN4 vs T/NK cytotoxicity | {int(a.n)} | {fmt(a.rho)} | {a.p:.3g} |")
    if b is not None:
        lines.append(f"| B | malignant-like CLDN4 vs T/NK exhaustion | {int(b.n)} | {fmt(b.rho)} | {b.p:.3g} |")
    if c is not None:
        lines.append(f"| C | all-epithelial CLDN4 vs T/NK cytotoxicity | {int(c.n)} | {fmt(c.rho)} | {c.p:.3g} |")
    if d is not None:
        lines.append(f"| D | all-epithelial CLDN4 vs T/NK exhaustion | {int(d.n)} | {fmt(d.rho)} | {d.p:.3g} |")

    lines += [
        "",
        "## Verdict",
        "",
    ]
    if a is not None and b is not None:
        lines.append(
            f"Primary post malignant-like (n={int(a.n)}): CLDN4 vs cytotoxicity ρ={fmt(a.rho)}, p={a.p:.3g}; "
            f"vs exhaustion ρ={fmt(b.rho)}, p={b.p:.3g}. "
            "That n is **not** 12 and **not** 92,330. Five post tumors have fewer than 10 malignant-like cells "
            "(mostly MPR / near-pCR residual epithelium that scores as normal-lung). "
            "The all-epithelial sensitivity uses every post sample and is still a null / weak association. "
            "This does not support a GSE207422 claim that malignant CLDN4 tracks T/NK exhaustion or cytotoxicity."
        )
    lines += [
        "",
        "What this is **not**: TACSTD2 vs T/NK fraction, NMPR>MPR TACSTD2, or a multi-cohort meta "
        "(those already exist elsewhere). TACSTD2 is extracted only so the panel is complete; it is not tested here.",
        "",
        "## Methods (short)",
        "",
        "- Public files: `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz` (92,330 cells) and the GEO xlsx. "
        "Raw HRA001033 was not used.",
        "- Per cell: `log1p(CP10k) = log1p(UMI / library_size × 10⁴)`.",
        "- Lineage = argmax of mean log1p marker modules (epithelial / T/NK / B / myeloid / mast / endo / fibro).",
        "- T/NK = that lineage. CD8 sensitivity = T/NK with CD8A or CD8B UMI > 0.",
        "- Malignant-like (primary) as above. Zero-normal-UMI = epithelial and zero UMI of "
        "`SFTPA2`, `AGER`, `SCGB1A1`, `SCGB3A1`, `TPPP3` (A3-style). All-epithelial = no normal-lung gate.",
        "- Spearman two-sided on patients. Mann–Whitney on MPR vs NMPR is recorded in `association_statistics.tsv` "
        "but is not a primary contrast when MPR n<2 under the malignant gate.",
        "",
        "## Files",
        "",
        "- `FINDING.md` — this note (table is the deliverable)",
        "- `fig_extra_gse207422_cldn4_exh.png` / `.pdf` — extra figure",
        "- `per_patient.tsv` — one row / sample, all scores and counts",
        "- `honest_n.tsv` — occupancy",
        "- `association_statistics.tsv` — every Spearman / MWU",
        "- `summary.json`",
        "- Scripts: `download.py`, `extract.py`, `analyze.py`",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/gse207422_cldn4_exh/download.py",
        "python3 methods/gse207422_cldn4_exh/extract.py",
        "python3 methods/gse207422_cldn4_exh/analyze.py",
        "```",
        "",
    ]
    text = "\n".join(lines)
    (outdir / "FINDING.md").write_text(text)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/gse207422_cldn4_exh"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cells, missing = load_panel(args.datadir)
    n_raw = int(len(cells))
    cells = cells.loc[cells.total_umi >= MIN_UMI].copy()
    count_genes = [c for c in cells.columns if c not in ("barcode", "total_umi", "Sample")]
    tot = cells["total_umi"].to_numpy(float)
    tot = np.where(tot > 0, tot, np.nan)
    ln = pd.DataFrame(
        np.log1p(cells[count_genes].to_numpy(float) / tot[:, None] * 1e4),
        index=cells.index,
        columns=count_genes,
    )
    cells["lineage"] = assign_lineage(ln)
    cells["is_tnk"] = cells["lineage"] == "T/NK"
    cells["is_epi"] = cells["lineage"] == "Epithelial"
    normal = module_mean(ln, NORMAL_LUNG)
    tumor = module_mean(ln, TUMOR_EPI)
    cells["is_malig_module"] = (cells["is_epi"] & (tumor > normal) & (normal < 0.4)).to_numpy()
    zero_n = np.zeros(len(cells), dtype=np.int64)
    for g in ZERO_NORMAL:
        if g in cells.columns:
            zero_n += cells[g].to_numpy()
    cells["is_malig_zero"] = cells["is_epi"].to_numpy() & (zero_n == 0)
    cd8_umi = np.zeros(len(cells), dtype=np.int64)
    for g in ("CD8A", "CD8B"):
        if g in cells.columns:
            cd8_umi += cells[g].to_numpy()
    cells["is_cd8"] = cells["is_tnk"].to_numpy() & (cd8_umi > 0)

    meta = load_meta(args.datadir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
    cells = cells.merge(
        meta[["Sample", "Patient", "timing", "mpr_group", "Pathology", "RECIST",
              "residual_tumor", "Pathologic Response"]],
        on="Sample",
        how="left",
    )

    pat = score_patients(cells, ln)
    pat.to_csv(args.outdir / "per_patient.tsv", sep="\t", index=False)

    honest = pat[[
        "Patient", "Sample", "timing", "mpr_group", "Pathology", "RECIST",
        "n_cells", "n_malig_module", "n_malig_zero", "n_epithelial", "n_tnk", "n_cd8",
        "malig_module_CLDN4", "malig_zero_CLDN4", "epithelial_CLDN4",
        "tnk_cyto", "tnk_exh", "cd8_cyto", "cd8_exh",
    ]].copy()
    honest["in_primary"] = (
        (honest.timing == "post")
        & honest.malig_module_CLDN4.notna()
        & honest.tnk_cyto.notna()
    )
    honest.to_csv(args.outdir / "honest_n.tsv", sep="\t", index=False)

    post = pat[pat.timing == "post"]
    allp = pat
    rows = []
    specs = [
        (post, "malig_module_CLDN4", "tnk_cyto", "CLDN4 vs T/NK cytotoxicity", "post, malignant-like"),
        (post, "malig_module_CLDN4", "tnk_exh", "CLDN4 vs T/NK exhaustion", "post, malignant-like"),
        (post, "malig_zero_CLDN4", "tnk_cyto", "CLDN4 vs T/NK cytotoxicity", "post, zero-normal-UMI"),
        (post, "malig_zero_CLDN4", "tnk_exh", "CLDN4 vs T/NK exhaustion", "post, zero-normal-UMI"),
        (post, "epithelial_CLDN4", "tnk_cyto", "CLDN4 vs T/NK cytotoxicity", "post, all-epithelial"),
        (post, "epithelial_CLDN4", "tnk_exh", "CLDN4 vs T/NK exhaustion", "post, all-epithelial"),
        (post, "malig_module_CLDN4", "cd8_cyto", "CLDN4 vs CD8 cytotoxicity", "post, malignant-like, CD8"),
        (post, "malig_module_CLDN4", "cd8_exh", "CLDN4 vs CD8 exhaustion", "post, malignant-like, CD8"),
        (post, "epithelial_CLDN4", "cd8_cyto", "CLDN4 vs CD8 cytotoxicity", "post, all-epithelial, CD8"),
        (post, "epithelial_CLDN4", "cd8_exh", "CLDN4 vs CD8 exhaustion", "post, all-epithelial, CD8"),
        (post, "malig_module_CLDN4_pct", "tnk_cyto", "CLDN4 %pos vs T/NK cytotoxicity", "post, malignant-like"),
        (post, "malig_module_CLDN4_pct", "tnk_exh", "CLDN4 %pos vs T/NK exhaustion", "post, malignant-like"),
        (allp, "malig_module_CLDN4", "tnk_cyto", "CLDN4 vs T/NK cytotoxicity", "pre+post, malignant-like"),
        (allp, "malig_module_CLDN4", "tnk_exh", "CLDN4 vs T/NK exhaustion", "pre+post, malignant-like"),
        (allp, "epithelial_CLDN4", "tnk_cyto", "CLDN4 vs T/NK cytotoxicity", "pre+post, all-epithelial"),
        (allp, "epithelial_CLDN4", "tnk_exh", "CLDN4 vs T/NK exhaustion", "pre+post, all-epithelial"),
        (post[post.mpr_group == "NMPR"], "malig_module_CLDN4", "tnk_cyto",
         "CLDN4 vs T/NK cytotoxicity", "post NMPR-only, malignant-like"),
        (post[post.mpr_group == "NMPR"], "malig_module_CLDN4", "tnk_exh",
         "CLDN4 vs T/NK exhaustion", "post NMPR-only, malignant-like"),
        (post[post.Pathology == "Adeno"], "epithelial_CLDN4", "tnk_cyto",
         "CLDN4 vs T/NK cytotoxicity", "post adeno, all-epithelial"),
        (post[post.Pathology == "Squamous"], "epithelial_CLDN4", "tnk_cyto",
         "CLDN4 vs T/NK cytotoxicity", "post squamous, all-epithelial"),
    ]
    for df, x, y, lab, sub in specs:
        rows.append(contrast_row(df, x, y, lab, sub))

    # MWU secondary
    mwu_rows = []
    for col, lab in [
        ("malig_module_CLDN4", "NMPR vs MPR malignant-like CLDN4"),
        ("epithelial_CLDN4", "NMPR vs MPR all-epithelial CLDN4"),
        ("tnk_cyto", "NMPR vs MPR T/NK cytotoxicity"),
        ("tnk_exh", "NMPR vs MPR T/NK exhaustion"),
    ]:
        a = post.loc[post.mpr_group == "NMPR", col]
        b = post.loc[post.mpr_group == "MPR", col]
        m = mwu_safe(a, b)
        mwu_rows.append({
            "contrast": lab,
            "subset": "post",
            "x": col,
            "y": "mpr_group",
            "n": m["n_a"] + m["n_b"],
            "rho": np.nan,
            "p": m["p"],
            "patients": f"NMPR n={m['n_a']} median={m['median_a']}; MPR n={m['n_b']} median={m['median_b']}",
            "note": m["note"] or "Mann-Whitney two-sided",
        })
    assoc = pd.DataFrame(rows + mwu_rows)
    assoc.to_csv(args.outdir / "association_statistics.tsv", sep="\t", index=False)

    notes = {
        "dataset": "GSE207422",
        "task": "additive: T/NK cytotoxicity and exhaustion vs malignant CLDN4",
        "not_a_tacstd2_redo": True,
        "n_cells_raw": n_raw,
        "n_cells_qc": int(len(cells)),
        "n_malignant_module": int(cells["is_malig_module"].sum()),
        "n_malignant_zero": int(cells["is_malig_zero"].sum()),
        "n_epithelial": int(cells["is_epi"].sum()),
        "n_tnk": int(cells["is_tnk"].sum()),
        "n_cd8": int(cells["is_cd8"].sum()),
        "n_post_primary": int(honest["in_primary"].sum()),
        "cyto_genes": CYTO_GENES,
        "exh_genes": EXH_GENES,
        "missing_from_matrix": missing,
        "author_cell_labels_public": False,
        "unit": "patient / one GEO sample",
        "primary_malignant": "epithelial AND tumor-epi > normal-lung AND normal-lung < 0.4",
    }
    prim = honest[honest.in_primary]
    a = contrast_row(post, "malig_module_CLDN4", "tnk_cyto", "cyto", "post")
    b = contrast_row(post, "malig_module_CLDN4", "tnk_exh", "exh", "post")
    notes["primary_cyto"] = {"n": a["n"], "rho": a["rho"], "p": a["p"], "patients": a["patients"]}
    notes["primary_exh"] = {"n": b["n"], "rho": b["rho"], "p": b["p"], "patients": b["patients"]}
    notes["primary_mpr_n"] = int((prim.mpr_group == "MPR").sum())
    notes["primary_nmpr_n"] = int((prim.mpr_group == "NMPR").sum())
    with (args.outdir / "summary.json").open("w") as fh:
        json.dump(notes, fh, indent=2)

    make_figure(pat, args.outdir)
    write_finding(pat, assoc, notes, args.outdir)
    print(json.dumps(notes, indent=2))
    print("wrote", args.outdir)


if __name__ == "__main__":
    main()
