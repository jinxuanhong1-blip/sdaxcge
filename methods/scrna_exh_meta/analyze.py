#!/usr/bin/env python3
"""Patient-level T/NK cytotoxicity and exhaustion vs malignant TACSTD2/CLDN4.

Public lung ICI / neoadjuvant scRNA only. Patient is the unit.
Meta-analysis is Fisher-z inverse-variance (FE + DerSimonian–Laird RE).
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent

CYTO_GENES = ["GZMB", "PRF1", "GNLY", "NKG7"]
EXH_GENES = ["PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX"]
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
MIN_UMI = 200
MIN_MAL = 10
MIN_TNK = 10
MIN_N_META = 5


def present(genes, columns) -> list[str]:
    return [g for g in genes if g in columns]


def log1p_cp10k(counts: pd.DataFrame, total: pd.Series) -> pd.DataFrame:
    tot = np.asarray(total, float)
    tot = np.where(tot > 0, tot, np.nan)
    return pd.DataFrame(
        np.log1p(np.asarray(counts, float) / tot[:, None] * 1e4),
        index=counts.index,
        columns=counts.columns,
    )


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


def malignant_like(ln: pd.DataFrame, lineage: pd.Series) -> pd.Series:
    is_epi = lineage == "Epithelial"
    normal = module_mean(ln, NORMAL_LUNG)
    tumor = module_mean(ln, TUMOR_EPI)
    return is_epi & (tumor > normal) & (normal < 0.4)


def spearman_safe(x, y) -> dict:
    x, y = np.asarray(x, float), np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    n = int(x.size)
    if n < 4:
        return {"n": n, "rho": np.nan, "p": np.nan}
    rho, p = stats.spearmanr(x, y)
    return {"n": n, "rho": float(rho), "p": float(p)}


def fisher_z_meta(rows: list[dict]) -> dict:
    usable = [r for r in rows if r["n"] >= MIN_N_META and np.isfinite(r.get("rho", np.nan))]
    if len(usable) < 2:
        return {
            "k": len(usable),
            "n_total": int(sum(r["n"] for r in usable)),
            "rho_fe": np.nan, "p_fe": np.nan, "ci_fe_lo": np.nan, "ci_fe_hi": np.nan,
            "rho_re": np.nan, "p_re": np.nan, "ci_re_lo": np.nan, "ci_re_hi": np.nan,
            "Q": np.nan, "I2": np.nan, "tau2": np.nan,
        }
    z = np.array([math.atanh(max(min(r["rho"], 0.999999), -0.999999)) for r in usable])
    n = np.array([r["n"] for r in usable], float)
    se = 1.0 / np.sqrt(n - 3.0)
    w = 1.0 / se ** 2
    z_fe = float(np.sum(w * z) / np.sum(w))
    se_fe = float(1.0 / math.sqrt(np.sum(w)))
    Q = float(np.sum(w * (z - z_fe) ** 2))
    df = len(usable) - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (Q - df) / c) if c > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q) * 100.0 if Q > 0 else 0.0
    w_re = 1.0 / (se ** 2 + tau2)
    z_re = float(np.sum(w_re * z) / np.sum(w_re))
    se_re = float(1.0 / math.sqrt(np.sum(w_re)))

    def pack(zval, seval, prefix):
        lo, hi = zval - 1.96 * seval, zval + 1.96 * seval
        p = float(2 * stats.norm.sf(abs(zval / seval))) if seval > 0 else np.nan
        return {
            f"rho_{prefix}": float(math.tanh(zval)),
            f"p_{prefix}": p,
            f"ci_{prefix}_lo": float(math.tanh(lo)),
            f"ci_{prefix}_hi": float(math.tanh(hi)),
        }

    out = {"k": len(usable), "n_total": int(n.sum()), "Q": Q, "I2": I2, "tau2": tau2}
    out.update(pack(z_fe, se_fe, "fe"))
    out.update(pack(z_re, se_re, "re"))
    return out


def score_patients(cells: pd.DataFrame, ln: pd.DataFrame, patient_col: str, extra_cols: list[str]) -> pd.DataFrame:
    cells = cells.copy()
    cells["cyto"] = module_mean(ln, CYTO_GENES)
    cells["exh"] = module_mean(ln, EXH_GENES)
    cells["tac"] = ln["TACSTD2"] if "TACSTD2" in ln.columns else np.nan
    cells["cld"] = ln["CLDN4"] if "CLDN4" in ln.columns else np.nan
    rows = []
    for pid, sub in cells.groupby(patient_col, sort=True):
        mal = sub[sub["is_malignant"]]
        tnk = sub[sub["is_tnk"]]
        rec = {
            "patient": pid,
            "n_cells": int(len(sub)),
            "n_malignant": int(len(mal)),
            "n_tnk": int(len(tnk)),
            "tnk_frac": float(len(tnk) / len(sub)) if len(sub) else np.nan,
            "mal_TACSTD2": float(mal["tac"].mean()) if len(mal) else np.nan,
            "mal_CLDN4": float(mal["cld"].mean()) if len(mal) else np.nan,
            "tnk_cyto": float(tnk["cyto"].mean()) if len(tnk) else np.nan,
            "tnk_exh": float(tnk["exh"].mean()) if len(tnk) else np.nan,
        }
        for c in extra_cols:
            if c in sub.columns:
                rec[c] = sub[c].iloc[0]
        if rec["n_malignant"] < MIN_MAL:
            rec["mal_TACSTD2"] = np.nan
            rec["mal_CLDN4"] = np.nan
        if rec["n_tnk"] < MIN_TNK:
            rec["tnk_cyto"] = np.nan
            rec["tnk_exh"] = np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def cohort_spearmans(cohort: str, pat: pd.DataFrame) -> list[dict]:
    rows = []
    for x, y, xlab, ylab in [
        ("mal_TACSTD2", "tnk_cyto", "TACSTD2", "cytotoxicity"),
        ("mal_TACSTD2", "tnk_exh", "TACSTD2", "exhaustion"),
        ("mal_CLDN4", "tnk_cyto", "CLDN4", "cytotoxicity"),
        ("mal_CLDN4", "tnk_exh", "CLDN4", "exhaustion"),
    ]:
        s = spearman_safe(pat[x], pat[y])
        rows.append({
            "cohort": cohort,
            "malignant_gene": xlab,
            "tnk_score": ylab,
            "n": s["n"],
            "rho": s["rho"],
            "p": s["p"],
            "n_malignant_median": float(pat.loc[pat[x].notna() & pat[y].notna(), "n_malignant"].median())
            if s["n"] else np.nan,
            "n_tnk_median": float(pat.loc[pat[x].notna() & pat[y].notna(), "n_tnk"].median())
            if s["n"] else np.nan,
        })
    return rows


def gene_note(ln: pd.DataFrame) -> dict:
    return {
        "cyto_present": present(CYTO_GENES, ln.columns),
        "exh_present": present(EXH_GENES, ln.columns),
        "targets_present": present(["TACSTD2", "CLDN4"], ln.columns),
    }


def analyze_gse207422(datadir: Path) -> tuple[pd.DataFrame, list[dict], dict]:
    cells = pd.read_csv(datadir / "gse207422_panel_cells.tsv.gz", sep="\t")
    cells["sample"] = cells["barcode"].str.rsplit("_", n=1).str[0]
    cells = cells.loc[cells.total_umi >= MIN_UMI].copy()
    count_genes = [c for c in cells.columns if c not in ("barcode", "total_umi", "sample")]
    ln = log1p_cp10k(cells[count_genes], cells["total_umi"])
    cells["lineage"] = assign_lineage(ln)
    cells["is_tnk"] = cells["lineage"] == "T/NK"
    cells["is_malignant"] = malignant_like(ln, cells["lineage"]).to_numpy()
    meta = pd.read_excel(datadir / "GSE207422_NSCLC_scRNAseq_metadata.xlsx").iloc[:15]
    meta = meta.rename(columns={"Sample": "sample"})
    meta["mpr_group"] = meta["Pathologic Response"].map({"MPR": "MPR", "pCR": "MPR", "NMPR": "NMPR"})
    meta["timing"] = meta["Resource"].map({
        "Pre-treatment biopsy": "pre", "Post-treatment surgery": "post",
    })
    cells = cells.merge(meta[["sample", "Patient", "timing", "mpr_group", "Pathology"]], on="sample", how="left")
    # one post-tx sample per patient is the ICI-exposed unit
    post = cells[cells.timing == "post"].copy()
    pat = score_patients(post, ln.loc[post.index], "Patient", ["mpr_group", "Pathology", "timing"])
    pat["cohort"] = "GSE207422"
    notes = gene_note(ln)
    notes.update({
        "label_source": "marker lineage; malignant-like = epithelial AND tumor-epi > normal-lung AND normal-lung < 0.4",
        "unit": "post-treatment surgery sample (1/patient)",
        "n_cells_qc": int(len(cells)),
        "n_malignant_post": int(post["is_malignant"].sum()),
        "n_tnk_post": int(post["is_tnk"].sum()),
        "cyto_genes": CYTO_GENES,
        "exh_genes": EXH_GENES,
    })
    return pat, cohort_spearmans("GSE207422", pat), notes


def analyze_gse205335(datadir: Path) -> tuple[pd.DataFrame, list[dict], dict]:
    cells = pd.read_csv(datadir / "gse205335_panel_cells.tsv.gz", sep="\t")
    ident = pd.read_csv(datadir / "GSE205335_Lung_IO_CellIdentity.txt.gz", sep="\t")
    samp = pd.read_csv(HERE / "gse205335_sample_table.tsv", sep="\t")
    cells = cells.merge(ident, on="barcode", how="inner")
    cells = cells.merge(samp[["orig.ident", "patient", "tissue", "recist", "cancer_subtype"]], on="orig.ident")
    cells = cells.loc[cells.total_umi >= MIN_UMI].copy()
    count_genes = [c for c in cells.columns if c not in (
        "barcode", "total_umi", "orig.ident", "core.patient", "cluster.total",
        "lineage.total", "lineage.sub", "cluster.sub", "celltype",
        "patient", "tissue", "recist", "cancer_subtype",
    )]
    ln = log1p_cp10k(cells[count_genes], cells["total_umi"])
    cells["is_malignant"] = cells["lineage.sub"] == "Malignant cells"
    cells["is_tnk"] = cells["lineage.total"] == "T/NK cells"
    tumor = cells[~cells["tissue"].astype(str).str.startswith("Normal")].copy()
    tumor["response"] = tumor["recist"].map({"PR": "R", "SD": "NR", "PD": "NR"})
    pat = score_patients(tumor, ln.loc[tumor.index], "patient", ["response", "cancer_subtype"])
    pat["cohort"] = "GSE205335"
    notes = gene_note(ln)
    notes.update({
        "label_source": "author lineage.sub / lineage.total",
        "unit": "mean of non-normal-tissue samples per patient",
        "n_cells_qc": int(len(cells)),
        "n_malignant_tumor": int(tumor["is_malignant"].sum()),
        "n_tnk_tumor": int(tumor["is_tnk"].sum()),
    })
    return pat, cohort_spearmans("GSE205335", pat), notes


def analyze_gse241934(datadir: Path, split: str) -> tuple[pd.DataFrame, list[dict], dict]:
    tag = "iit" if split == "IIT" else "rwc"
    cells = pd.read_csv(datadir / f"gse241934_{tag}_panel_cells.tsv.gz", sep="\t")
    meta = pd.read_csv(datadir / f"gse241934_{tag}_meta_join.tsv.gz", sep="\t")
    cells = cells.merge(meta, on="barcode", how="left")
    cells = cells.loc[cells.total_umi >= MIN_UMI].copy()
    skip = {
        "barcode", "total_umi", "orig.ident", "sampleID", "major.cell.type",
        "major_cell_type", "Pathological Response", "EGFR", "PD1", "Histology",
        "nCount_RNA", "nFeature_RNA", "percent.mt",
    }
    count_genes = [c for c in cells.columns if c not in skip]
    ln = log1p_cp10k(cells[count_genes], cells["total_umi"])
    maj = cells["major.cell.type"].fillna(cells.get("major_cell_type"))
    cells["is_malignant"] = maj == "Epi"
    cells["is_tnk"] = maj.isin(["T", "NK"])
    cells["mpr_group"] = cells["Pathological Response"].replace({"pCR": "MPR", "non-MPR": "NMPR"})
    cohort = f"GSE241934_{split}"
    pat = score_patients(cells, ln, "sampleID", ["mpr_group", "Histology", "PD1"])
    pat["cohort"] = cohort
    notes = gene_note(ln)
    notes.update({
        "label_source": "author major.cell.type; Epi used as the malignant compartment (no CopyKAT on GEO)",
        "unit": "one resected tumor / patient",
        "n_cells_qc": int(len(cells)),
        "n_epi": int(cells["is_malignant"].sum()),
        "n_tnk": int(cells["is_tnk"].sum()),
        "split": split,
    })
    return pat, cohort_spearmans(cohort, pat), notes


def analyze_marker_10x(path: Path, cohort: str, extra_cols: list[str]) -> tuple[pd.DataFrame, list[dict], dict]:
    cells = pd.read_csv(path, sep="\t")
    cells = cells.loc[cells.total_umi >= MIN_UMI].copy()
    skip = {"barcode", "total_umi", "sample", "patient", "response"}
    count_genes = [c for c in cells.columns if c not in skip]
    ln = log1p_cp10k(cells[count_genes], cells["total_umi"])
    cells["lineage"] = assign_lineage(ln)
    cells["is_tnk"] = cells["lineage"] == "T/NK"
    cells["is_malignant"] = malignant_like(ln, cells["lineage"]).to_numpy()
    pat = score_patients(cells, ln, "patient", extra_cols)
    pat["cohort"] = cohort
    notes = gene_note(ln)
    notes.update({
        "label_source": "marker lineage; malignant-like = epithelial AND tumor-epi > normal-lung AND normal-lung < 0.4",
        "unit": "one sample / patient",
        "n_cells_qc": int(len(cells)),
        "n_malignant": int(cells["is_malignant"].sum()),
        "n_tnk": int(cells["is_tnk"].sum()),
    })
    return pat, cohort_spearmans(cohort, pat), notes


def forest_panel(ax, rows: list[dict], meta: dict, title: str) -> None:
    labels = []
    rhos = []
    los = []
    his = []
    for r in rows:
        if not np.isfinite(r.get("rho", np.nan)) or r["n"] < 4:
            continue
        z = math.atanh(max(min(r["rho"], 0.999999), -0.999999))
        se = 1.0 / math.sqrt(r["n"] - 3) if r["n"] > 3 else np.nan
        lo, hi = math.tanh(z - 1.96 * se), math.tanh(z + 1.96 * se)
        labels.append(f"{r['cohort']}  n={r['n']}")
        rhos.append(r["rho"])
        los.append(lo)
        his.append(hi)
    if meta.get("k", 0) >= 2 and np.isfinite(meta.get("rho_re", np.nan)):
        labels.append(f"RE meta  k={meta['k']}  N={meta['n_total']}")
        rhos.append(meta["rho_re"])
        los.append(meta["ci_re_lo"])
        his.append(meta["ci_re_hi"])
    y = np.arange(len(labels))[::-1]
    ax.axvline(0, color="0.6", lw=0.8)
    for yi, rho, lo, hi, lab in zip(y, rhos, los, his, labels):
        is_meta = lab.startswith("RE")
        ax.plot([lo, hi], [yi, yi], color="#1f4e79" if is_meta else "#4c78a8", lw=2 if is_meta else 1.4)
        ax.plot(rho, yi, "D" if is_meta else "o", color="#1f4e79" if is_meta else "#4c78a8", ms=6 if is_meta else 5)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlim(-1.05, 1.05)
    ax.set_xlabel("Spearman ρ")
    p = meta.get("p_re", np.nan)
    rho = meta.get("rho_re", np.nan)
    i2 = meta.get("I2", np.nan)
    ax.set_title(
        f"{title}\nRE ρ={rho:.2f} (95% {meta.get('ci_re_lo', np.nan):.2f},{meta.get('ci_re_hi', np.nan):.2f})  "
        f"p={p:.2g}  I²={i2:.0f}%",
        fontsize=8,
    )


def write_methods_snippet(outdir: Path, notes: dict, spear: pd.DataFrame, meta_rows: list[dict]) -> None:
    lines = [
        "# Extra figure — T/NK cytotoxicity and exhaustion vs malignant TACSTD2 / CLDN4",
        "",
        "Public lung ICI / neoadjuvant scRNA with both a malignant (or author epithelial) compartment and T/NK.",
        "Patient is the unit. Scores are the mean of log1p(CP10k) over the named genes in the named compartment.",
        "",
        "- Cytotoxicity in T/NK: GZMB, PRF1, GNLY, NKG7",
        "- Exhaustion in T/NK: PDCD1, HAVCR2, LAG3, TIGIT, TOX",
        "- Malignant TACSTD2 and CLDN4: mean log1p(CP10k) in malignant / author-Epi cells",
        "- Keep a patient if ≥10 malignant and ≥10 T/NK cells",
        "- Spearman on the patient table; Fisher-z IVW fixed-effect and DerSimonian–Laird random-effect meta (n≥5 / cohort)",
        "",
        "## Cohorts",
        "",
        "| Cohort | Setting | Labels | Patients in Spearman |",
        "|---|---|---|---|",
    ]
    n_by = spear.groupby("cohort")["n"].max()
    setting = {
        "GSE207422": "neoadjuvant PD-1 + chemo (Hu 2023); post-tx; marker malignant-like",
        "GSE205335": "palliative ICI atlas (Park/Ahn/Lee); author malignant + T/NK",
        "GSE241934_IIT": "NEOTIDE EGFR-mut neoadjuvant sintilimab + chemo; author Epi + T/NK",
        "GSE241934_RWC": "real-world neoadjuvant PD-1 + chemo; author Epi + T/NK",
        "GSE291670": "neoadjuvant anlotinib + camrelizumab; marker malignant-like",
        "GSE233203": "leftover: pre-ABCP pleural effusion; marker compartments",
    }
    for c, n in n_by.items():
        lines.append(f"| {c} | {setting.get(c, '')} | see METHODS.md | {int(n)} |")
    lines += ["", "## Patient-level Spearman (honest n / ρ / p)", "",
              "| Contrast | Cohort | n | ρ | p |",
              "|---|---|---:|---:|---:|"]
    for _, r in spear.sort_values(["malignant_gene", "tnk_score", "cohort"]).iterrows():
        rho = "NA" if not np.isfinite(r["rho"]) else f"{r['rho']:.3f}"
        p = "NA" if not np.isfinite(r["p"]) else f"{r['p']:.3g}"
        lines.append(f"| {r['malignant_gene']} vs T/NK {r['tnk_score']} | {r['cohort']} | {int(r['n'])} | {rho} | {p} |")
    lines += ["", "## Random-effect meta (Fisher z, n≥5)", "",
              "| Contrast | k | N | ρ_RE | 95% CI | p_RE | I² |",
              "|---|---:|---:|---:|---|---:|---:|"]
    for m in meta_rows:
        if not np.isfinite(m.get("rho_re", np.nan)):
            lines.append(f"| {m['contrast']} | {m['k']} | {m['n_total']} | NA | — | NA | NA |")
            continue
        lines.append(
            f"| {m['contrast']} | {m['k']} | {m['n_total']} | {m['rho_re']:.3f} | "
            f"{m['ci_re_lo']:.3f} to {m['ci_re_hi']:.3f} | {m['p_re']:.3g} | {m['I2']:.0f}% |"
        )
    lines += [
        "",
        "## Notes on n",
        "",
        "- GSE207422: 12 post-tx samples; **7** have ≥10 malignant-like cells (P02/P06/P11/P13/P14 drop).",
        "- GSE241934 RWC: 34 tumors; **29** have ≥10 author-Epi cells.",
        "- GSE233203 leftover: 7 PE samples; **6** have ≥10 T/NK (NCCLu_397 has 3 T/NK).",
        "- GSE291670 TACSTD2 vs cytotoxicity is the only cohort-level p<0.05 (**n=6, ρ=0.943, p=0.0048**). That point is why I²=68% on that contrast; dropping it gives RE ρ=−0.03, p=0.88 (N=75).",
        "- Named-only (drop leftover): TACSTD2–cyto RE ρ=0.16, p=0.58, N=75. Named neoadjuvant-only TACSTD2–exhaustion RE ρ=−0.30, p=0.045, N=53, I²=0% (post-hoc slice; not the primary).",
        "",
        "GSE146100 is a leftover with both compartments but n=1 patient and is not in the meta.",
        "GSE243013 / GSE266035 / T-sorted series lack a malignant compartment and are not scored.",
        "",
    ]
    (outdir / "WRITEUP.md").write_text("\n".join(lines))
    (Path("/workspace/paper") / "extra_scrna_exh_meta.md").write_text("\n".join(lines))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=Path("/tmp/scrna_exh_meta"))
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    figdir = Path("/workspace/results/scrna_exh_meta")
    figdir.mkdir(parents=True, exist_ok=True)

    all_pat = []
    all_spear = []
    notes = {}

    pat, sp, nt = analyze_gse207422(args.datadir)
    all_pat.append(pat); all_spear += sp; notes["GSE207422"] = nt
    print("GSE207422", pat[["patient", "n_malignant", "n_tnk", "mal_TACSTD2", "tnk_cyto", "tnk_exh"]].to_string(), flush=True)

    pat, sp, nt = analyze_gse205335(args.datadir)
    all_pat.append(pat); all_spear += sp; notes["GSE205335"] = nt
    print("GSE205335 n", len(pat), flush=True)

    for split in ("IIT", "RWC"):
        pat, sp, nt = analyze_gse241934(args.datadir, split)
        all_pat.append(pat); all_spear += sp; notes[f"GSE241934_{split}"] = nt
        print(f"GSE241934_{split} n", len(pat), flush=True)

    pat, sp, nt = analyze_marker_10x(args.datadir / "gse291670_panel_cells.tsv.gz", "GSE291670", ["response"])
    all_pat.append(pat); all_spear += sp; notes["GSE291670"] = nt
    print("GSE291670", pat.to_string(), flush=True)

    leftover_path = args.datadir / "gse233203_panel_cells.tsv.gz"
    leftover_included = False
    if leftover_path.exists():
        pat, sp, nt = analyze_marker_10x(leftover_path, "GSE233203", ["response"])
        both = int(((pat.n_malignant >= MIN_MAL) & (pat.n_tnk >= MIN_TNK)).sum())
        nt["patients_with_both"] = both
        notes["GSE233203"] = nt
        print("GSE233203 both-compartment patients", both, flush=True)
        if both >= MIN_N_META:
            all_pat.append(pat)
            all_spear += sp
            leftover_included = True
        else:
            pat.to_csv(outdir / "gse233203_patients_not_in_meta.tsv", sep="\t", index=False)
            notes["GSE233203"]["meta"] = "excluded: fewer than 5 patients with both compartments"

    patients = pd.concat(all_pat, ignore_index=True)
    spear = pd.DataFrame(all_spear)
    patients.to_csv(outdir / "patient_scores.tsv", sep="\t", index=False)
    spear.to_csv(outdir / "cohort_spearman.tsv", sep="\t", index=False)

    meta_rows = []
    contrasts = [
        ("TACSTD2", "cytotoxicity"),
        ("TACSTD2", "exhaustion"),
        ("CLDN4", "cytotoxicity"),
        ("CLDN4", "exhaustion"),
    ]
    for gene, score in contrasts:
        sub = spear[(spear.malignant_gene == gene) & (spear.tnk_score == score)].to_dict("records")
        m = fisher_z_meta(sub)
        m["contrast"] = f"{gene} vs T/NK {score}"
        m["malignant_gene"] = gene
        m["tnk_score"] = score
        meta_rows.append(m)
    meta_df = pd.DataFrame(meta_rows)
    meta_df.to_csv(outdir / "meta_spearman.tsv", sep="\t", index=False)

    # extra figure
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.4))
    for ax, (gene, score), mrow in zip(axes.ravel(), contrasts, meta_rows):
        sub = spear[(spear.malignant_gene == gene) & (spear.tnk_score == score)].to_dict("records")
        forest_panel(ax, sub, mrow, f"malignant {gene} vs T/NK {score}")
    fig.suptitle(
        "EXTRA  public lung ICI/neoadjuvant scRNA\n"
        "patient-level Spearman, T/NK cytotoxicity (GZMB/PRF1/GNLY/NKG7) "
        "and exhaustion (PDCD1/HAVCR2/LAG3/TIGIT/TOX)",
        fontsize=11,
    )
    fig.tight_layout()
    for dest in (outdir, figdir):
        fig.savefig(dest / "fig_extra_scrna_exh_meta.png", dpi=160)
        fig.savefig(dest / "fig_extra_scrna_exh_meta.pdf")
    plt.close(fig)

    # scatter companion (same extra figure family)
    fig, axes = plt.subplots(2, 2, figsize=(10.6, 8.2))
    colors = {
        "GSE207422": "#1f77b4",
        "GSE205335": "#ff7f0e",
        "GSE241934_IIT": "#2ca02c",
        "GSE241934_RWC": "#98df8a",
        "GSE291670": "#d62728",
        "GSE233203": "#9467bd",
    }
    pairs = [
        ("mal_TACSTD2", "tnk_cyto", "malignant TACSTD2", "T/NK cytotoxicity"),
        ("mal_TACSTD2", "tnk_exh", "malignant TACSTD2", "T/NK exhaustion"),
        ("mal_CLDN4", "tnk_cyto", "malignant CLDN4", "T/NK cytotoxicity"),
        ("mal_CLDN4", "tnk_exh", "malignant CLDN4", "T/NK exhaustion"),
    ]
    for ax, (x, y, xl, yl) in zip(axes.ravel(), pairs):
        for cohort, sub in patients.groupby("cohort"):
            m = sub[[x, y]].dropna()
            ax.scatter(m[x], m[y], s=28, alpha=0.85, c=colors.get(cohort, "0.4"), label=f"{cohort} n={len(m)}")
        ax.set_xlabel(xl + "  mean log1p(CP10k)")
        ax.set_ylabel(yl + "  mean log1p(CP10k)")
        ax.legend(fontsize=6, loc="best")
    fig.suptitle("EXTRA  patient-level points (same scores as the forest)", fontsize=11)
    fig.tight_layout()
    for dest in (outdir, figdir):
        fig.savefig(dest / "fig_extra_scrna_exh_meta_scatter.png", dpi=160)
    plt.close(fig)

    # documented sensitivities (not the primary)
    sens_rows = []
    subsets = {
        "primary_all_six": None,
        "named_only_drop_GSE233203": spear.cohort != "GSE233203",
        "drop_GSE291670": spear.cohort != "GSE291670",
        "named_neoadjuvant": spear.cohort.isin(
            ["GSE207422", "GSE241934_IIT", "GSE241934_RWC", "GSE291670"]
        ),
    }
    for name, mask in subsets.items():
        src = spear if mask is None else spear[mask]
        for gene, score in contrasts:
            sub = src[(src.malignant_gene == gene) & (src.tnk_score == score)].to_dict("records")
            m = fisher_z_meta(sub)
            m["subset"] = name
            m["contrast"] = f"{gene} vs T/NK {score}"
            sens_rows.append(m)
    pd.DataFrame(sens_rows).to_csv(outdir / "sensitivity_meta.tsv", sep="\t", index=False)

    write_methods_snippet(outdir, notes, spear, meta_rows)
    summary = {
        "contrasts": meta_rows,
        "cohort_notes": notes,
        "leftover_GSE233203_in_meta": leftover_included,
        "min_malignant_cells": MIN_MAL,
        "min_tnk_cells": MIN_TNK,
        "min_n_for_meta": MIN_N_META,
        "cyto_genes": CYTO_GENES,
        "exh_genes": EXH_GENES,
    }
    with open(outdir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(meta_df.to_string(index=False), flush=True)
    print("wrote", outdir)


if __name__ == "__main__":
    main()
