#!/usr/bin/env python3
"""CLDN4-high malignant → T/NK ligand–receptor on TISCH2 GSE148071.

Primary: documented CellPhoneDB-style mean score on TISCH2 log2(TPM/10+1).
Secondary: LIANA mt.cellphonedb if importable.
CellChat is not run (no R). Cells are not treated as replicates.
"""

from __future__ import annotations

import argparse
import json
import traceback
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]


def log(msg: str) -> None:
    print(msg, flush=True)


def _as_str_array(arr) -> np.ndarray:
    out = []
    for x in arr:
        if isinstance(x, bytes):
            out.append(x.decode())
        else:
            out.append(str(x))
    return np.array(out, dtype=object)


def load_tisch_genes(path: Path, wanted: list[str]) -> tuple[pd.DataFrame, list[str], dict]:
    audit = {"path": str(path), "keys": [], "shape": None, "n_genes": None, "n_cells": None}
    with h5py.File(path, "r") as f:
        audit["keys"] = list(f.keys())
        grp = None
        if all(k in f for k in ("data", "indices", "indptr")):
            grp = f
        else:
            for key in f.keys():
                g = f[key]
                if isinstance(g, h5py.Group) and all(k in g for k in ("data", "indices", "indptr")):
                    grp = g
                    break
        if grp is None:
            raise ValueError(f"no sparse matrix in {path}: {audit['keys']}")

        def _names(candidates):
            for c in candidates:
                node = grp.get(c) if hasattr(grp, "get") else None
                if node is None:
                    node = f.get(c)
                if node is None:
                    continue
                if isinstance(node, h5py.Dataset):
                    return _as_str_array(node[:])
                if isinstance(node, h5py.Group):
                    for sub in ("name", "id", "gene_names"):
                        if sub in node:
                            return _as_str_array(node[sub][:])
            return None

        genes = _names(["features", "gene_names", "genes", "rownames"])
        barcodes = _names(["barcodes", "cell_names", "colnames"])
        shape = tuple(int(x) for x in grp["shape"][:]) if "shape" in grp else None
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]

    if shape is None:
        raise ValueError(f"{path}: missing shape")
    n0, n1 = int(shape[0]), int(shape[1])
    if genes is not None and barcodes is not None and len(genes) == n1 and len(barcodes) == n0:
        n_genes, n_cells = n1, n0
        mat = sp.csc_matrix((data, indices, indptr), shape=(n_cells, n_genes))

        def col(i: int) -> np.ndarray:
            return np.asarray(mat[:, i].todense()).ravel()
    else:
        n_genes, n_cells = n0, n1
        mat = sp.csc_matrix((data, indices, indptr), shape=(n_genes, n_cells))

        def col(i: int) -> np.ndarray:
            return np.asarray(mat.getrow(i).todense()).ravel()

    if genes is None or barcodes is None:
        raise ValueError(f"{path}: missing gene or barcode names")
    if len(genes) != n_genes or len(barcodes) != n_cells:
        raise ValueError(
            f"{path}: dim mismatch genes={len(genes)} barcodes={len(barcodes)} "
            f"shape=({n_genes},{n_cells})"
        )

    audit["shape"] = [n_genes, n_cells]
    audit["n_genes"] = n_genes
    audit["n_cells"] = n_cells
    gene_set = set(genes.tolist())
    present = [g for g in wanted if g in gene_set]
    audit["n_wanted"] = len(wanted)
    audit["n_present"] = len(present)
    audit["n_absent"] = len(wanted) - len(present)
    gene_index = {g: int(np.where(genes == g)[0][0]) for g in present}
    table = {"barcode": barcodes}
    for g, i in gene_index.items():
        table[g] = col(i).astype(np.float32)
    return pd.DataFrame(table).set_index("barcode"), present, audit


def pathway_of(ligand: str, receptor: str, cfg: dict) -> str:
    lig_units = set(str(ligand).split("+"))
    rec_units = set(str(receptor).split("+"))
    for name, block in cfg["pathways"].items():
        if lig_units & set(block["ligands"]) and rec_units & set(block["receptors"]):
            return name
    for name, block in cfg["pathways"].items():
        if lig_units & set(block["ligands"]) or rec_units & set(block["receptors"]):
            return name + "_partial"
    return "other"


def group_gene_stats(expr: pd.DataFrame, mask: np.ndarray, genes: list[str]) -> tuple[dict, dict, int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    sub = expr.loc[mask, genes]
    means = {g: float(sub[g].mean()) for g in genes}
    fracs = {g: float((sub[g] > 0).mean()) for g in genes}
    return means, fracs, n


def partner_from_stats(units: list[str], means: dict, fracs: dict) -> tuple[float, float]:
    m, f = [], []
    for g in units:
        if g not in means:
            return np.nan, np.nan
        m.append(means[g])
        f.append(fracs[g])
    return float(np.min(m)), float(np.min(f))


def score_pairs(pairs: pd.DataFrame, expr: pd.DataFrame, sender: np.ndarray, receiver: np.ndarray, genes: list[str], expr_prop: float) -> pd.DataFrame:
    s_mean, s_frac, n_s = group_gene_stats(expr, sender, genes)
    r_mean, r_frac, n_r = group_gene_stats(expr, receiver, genes)
    rows = []
    for rec in pairs.itertuples(index=False):
        lig_u = str(rec.ligand).split("+")
        rec_u = str(rec.receptor).split("+")
        l_mean, l_frac = partner_from_stats(lig_u, s_mean, s_frac)
        rec_m, rec_f = partner_from_stats(rec_u, r_mean, r_frac)
        if not np.isfinite(l_mean) or not np.isfinite(rec_m):
            continue
        rows.append(
            {
                "ligand": rec.ligand,
                "receptor": rec.receptor,
                "pathway": rec.pathway,
                "n_sender": n_s,
                "n_receiver": n_r,
                "ligand_mean": l_mean,
                "receptor_mean": rec_m,
                "ligand_frac": l_frac,
                "receptor_frac": rec_f,
                "cpdb_mean_score": 0.5 * (l_mean + rec_m),
                "product_score": l_mean * rec_m,
                "pass_expr_prop": bool((l_frac >= expr_prop) and (rec_f >= expr_prop)),
            }
        )
    return pd.DataFrame(rows)


def wilcoxon_safe(a, b) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 3 or np.allclose(a, b):
        return np.nan
    try:
        return float(stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided").pvalue)
    except ValueError:
        return np.nan


def md_cell(x) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        if abs(x) >= 0.01 or x == 0:
            return f"{x:.3f}"
        return f"{x:.3g}"
    return str(x)


def run_liana(expr: pd.DataFrame, labels: pd.Series, out_csv: Path, n_perms: int, max_per_group: int, seed: int) -> str:
    try:
        import anndata as ad
        import liana as li
    except Exception as exc:
        return f"LIANA_IMPORT_FAILED: {exc}"
    try:
        rng = np.random.default_rng(seed)
        keep_idx = []
        for _lab, idx in labels.groupby(labels, sort=False).groups.items():
            idx = np.asarray(list(idx))
            if len(idx) > max_per_group:
                idx = rng.choice(idx, size=max_per_group, replace=False)
            keep_idx.append(idx)
        keep_idx = np.concatenate(keep_idx)
        sub = expr.loc[keep_idx]
        adata = ad.AnnData(sub.to_numpy(dtype=np.float32))
        adata.obs_names = sub.index.astype(str)
        adata.var_names = sub.columns.astype(str)
        adata.obs["group"] = labels.loc[keep_idx].astype(str).values
        counts = adata.obs["group"].value_counts().to_dict()
        li.mt.cellphonedb(
            adata,
            groupby="group",
            resource_name="cellphonedb",
            expr_prop=0.10,
            n_perms=n_perms,
            use_raw=False,
            verbose=True,
            key_added="liana_res",
        )
        res = adata.uns["liana_res"].copy()
        res.to_csv(out_csv, index=False)
        return f"LIANA_OK n_edges={len(res)} downsampled={counts} file={out_csv.name}"
    except Exception as exc:
        return f"LIANA_RUN_FAILED: {exc}\n{traceback.format_exc()}"


def write_finding(
    outdir: Path,
    summary: dict,
    ntab: pd.DataFrame,
    pooled: pd.DataFrame,
    ranks: pd.DataFrame,
    liana_note: str,
) -> None:
    top_n = int(summary["top_table_n"])
    passed = pooled[pooled["pass_expr_prop"]].sort_values("cpdb_mean_score", ascending=False)
    top = passed.head(top_n)
    focus_paths = ["T_recruit", "checkpoint", "MHC_I"]
    focus = passed[passed["pathway"].isin(focus_paths)].sort_values(
        ["pathway", "cpdb_mean_score"], ascending=[True, False]
    )

    def table(df: pd.DataFrame, extra: list[str] | None = None) -> list[str]:
        cols = ["pathway", "ligand", "receptor", "ligand_frac", "receptor_frac", "cpdb_mean_score"]
        if extra:
            cols += extra
        lines = [
            "| " + " | ".join(cols) + " |",
            "|" + "|".join(["---"] * 3 + ["---:"] * (len(cols) - 3)) + "|",
        ]
        for r in df.itertuples(index=False):
            vals = []
            for c in cols:
                v = getattr(r, c)
                if c in {"ligand_frac", "receptor_frac", "cpdb_mean_score", "median_delta"}:
                    vals.append(md_cell(float(v)))
                elif c in {"n_patients", "n_sender", "n_receiver"}:
                    vals.append(str(int(v)))
                elif c in {"pval", "padj"}:
                    vals.append(f"{float(v):.3g}" if np.isfinite(v) else "NA")
                else:
                    vals.append(str(v))
            lines.append("| " + " | ".join(vals) + " |")
        return lines

    n_pass = int(passed.shape[0])
    n_scored = int(pooled.shape[0])
    lines = [
        "# FINDING — GSE148071 LIANA/LR from CLDN4-high malignant to T/NK",
        "",
        f"**Object:** TISCH2 `NSCLC_GSE148071` (Wu et al. 2021 *Nat Commun*, PMID 33953163; GEO GSE148071).",
        f"**Sender:** TISCH2 `Malignant` cells with CLDN4 ≥ global malignant median ({summary['cldn4_threshold']:.3f} on TISCH2 log2(TPM/10+1)).",
        f"**Receiver:** TISCH2 T/NK major-lineage. In this object that is **CD8T + Tprolif only** — **0 NK / CD4T / Treg / NKT cells** are labeled.",
        f"**Unit of inference:** patient (n={summary['n_patients_meta']}; one Sample each). Cells are not replicates.",
        "",
        "## Honest n",
        "",
        "| Item | n | Note |",
        "| --- | ---: | --- |",
        f"| Patients in TISCH2 object | **{summary['n_patients_meta']}** | Wu et al. advanced NSCLC biopsies |",
        f"| Cells aligned to h5 | {summary['n_cells']:,} | all Tissue=Tumor |",
        f"| Malignant | {summary['n_malignant']:,} | TISCH2 major-lineage |",
        f"| T/NK total | **{summary['n_tnk']:,}** | CD8T={summary['n_cd8t']:,}; Tprolif={summary['n_tprolif']:,}; NK=0 |",
        f"| CLDN4-high / low malignant | {summary['n_cldn4_high']:,} / {summary['n_cldn4_low']:,} | median split among malignant |",
        f"| Patients with ≥{summary['min_malig_high']} CLDN4-high malignant **and** ≥{summary['min_tnk']} T/NK | **{summary['n_patients_usable']}** | used for per-patient scores |",
        f"| Patients also having ≥{summary['min_malig_low']} CLDN4-low malignant | **{summary['n_patients_paired']}** | used for high−low Wilcoxon |",
        f"| Patients dropped (too few T/NK or senders) | {summary['n_patients_dropped']} | see `n_cells_patients.tsv` |",
        f"| CPDB pairs scored / passing 10% expr_prop | {n_scored:,} / **{n_pass:,}** | pooled CLDN4-high → T/NK |",
        f"| LIANA | {summary['liana_status_short']} | {liana_note.split(chr(10))[0]} |",
        f"| CellChat | not run | R unavailable |",
        "",
        "Per-patient counts: `results/n_cells_patients.tsv`. Do not treat cell counts as the sample size.",
        "",
        f"Usable outgoing patients (n={summary['n_patients_usable']}): {', '.join(summary['usable_patients'])}.",
        f"Paired high vs low (n={summary['n_patients_paired']}): {', '.join(summary['paired_patients'])}.",
        "",
        "## What was run",
        "",
        "Primary score is the documented CellPhoneDB mean (Efremova 2020; Garcia-Alonso 2022) on the public TISCH2/MAESTRO matrix **log2(TPM/10+1)**: partner expression = min(subunit means); pair score = mean of the two partner means. A pair passes `expr_prop` when both partners are detected (>0) in ≥10% of cells in their group. This is **not** a CellChat probability and **not** a raw-UMI reprocess.",
        "",
        f"LIANA `mt.cellphonedb` status: `{summary['liana_status_short']}`. LIANA p-values (if present) are within-object specificity after downsampling, not patient-level tests.",
        "",
        "## LR table — pooled CLDN4-high malignant → T/NK",
        "",
        f"Top {min(top_n, len(top))} pairs by CPDB mean score among those passing 10% expression in **both** partners. Sender n={summary['n_cldn4_high']:,} cells; receiver n={summary['n_tnk']:,} cells. Pooled ranks are descriptive.",
        "",
    ]
    if top.empty:
        lines += ["No pair passed the 10% expression filter.", ""]
    else:
        lines += table(top) + [""]

    lines += [
        "### Focus axes that cleared 10% (T-recruit / checkpoint / MHC-I)",
        "",
    ]
    if focus.empty:
        lines += [
            "No T-recruit, checkpoint, or MHC-I pair passed 10% expression from CLDN4-high malignant to T/NK in this object.",
            "",
        ]
    else:
        lines += table(focus) + [""]

    if ranks is not None and len(ranks):
        lines += [
            "## Patient-level high vs low (honest paired n)",
            "",
            f"Paired Wilcoxon across **{summary['n_patients_paired']}** patients with ≥{summary['min_malig_high']} CLDN4-high, ≥{summary['min_malig_low']} CLDN4-low malignant, and ≥{summary['min_tnk']} T/NK. Δ = high − low CPDB score. FDR is BH within this contrast. Negative Δ = weaker from the CLDN4-high state.",
            "",
        ]
        show = ranks[ranks["pathway"].isin(focus_paths + ["T_recruit_partial", "checkpoint_partial", "MHC_I_partial"])].copy()
        if show.empty:
            show = ranks.head(15)
        else:
            show = show.sort_values(["pathway", "median_delta"])
        if len(show):
            lines += [
                "| pathway | ligand | receptor | n_patients | median_delta | pval | padj |",
                "| --- | --- | --- | ---: | ---: | ---: | ---: |",
            ]
            for r in show.itertuples(index=False):
                p = f"{r.pval:.3g}" if np.isfinite(r.pval) else "NA"
                q = f"{r.padj:.3g}" if np.isfinite(r.padj) else "NA"
                lines.append(
                    f"| {r.pathway} | {r.ligand} | {r.receptor} | {int(r.n_patients)} | {r.median_delta:+.3f} | {p} | {q} |"
                )
            n_neg = int((show["median_delta"] < 0).sum())
            n_sig = int((show["padj"] < 0.05).sum()) if show["padj"].notna().any() else 0
            lines += [
                "",
                f"Among the {len(show)} focus/partial pairs with ≥3 paired patients: {n_neg} have median Δ < 0; **{n_sig} reach FDR < 0.05**.",
                "",
            ]

    liana_csv = outdir / "liana_cellphonedb.csv"
    if liana_csv.exists() and summary.get("liana_status_short") == "LIANA_OK":
        li = pd.read_csv(liana_csv)
        src_col = next((c for c in li.columns if c.lower() in {"source", "ligand_complex"}), None)
        # LIANA 1.x columns: source, target, ligand_complex, receptor_complex, lr_means, cellphone_pvals
        if {"source", "target"}.issubset(li.columns):
            sub = li[(li["source"] == "Malig_CLDN4high") & (li["target"] == "TNK")].copy()
            score_col = "lr_means" if "lr_means" in sub.columns else None
            if score_col:
                sub = sub.sort_values(score_col, ascending=False).head(15)
            lines += [
                "## LIANA CellPhoneDB (secondary, pooled / downsampled)",
                "",
                "Edges with `source=Malig_CLDN4high` and `target=TNK` after ≤2,000 cells/group and 50 permutations.",
                "These p-values are within-object specificity, not patient-level tests.",
                "",
            ]
            if sub.empty:
                lines += ["No LIANA edges from CLDN4-high malignant to TNK passed `expr_prop=0.10`.", ""]
            else:
                lig_c = "ligand_complex" if "ligand_complex" in sub.columns else "ligand"
                rec_c = "receptor_complex" if "receptor_complex" in sub.columns else "receptor"
                p_c = "cellphone_pvals" if "cellphone_pvals" in sub.columns else None
                hdr = ["ligand", "receptor", "lr_means"] + (["cellphone_pvals"] if p_c else [])
                lines += [
                    "| " + " | ".join(hdr) + " |",
                    "|" + "|".join(["---"] * 2 + ["---:"] * (len(hdr) - 2)) + "|",
                ]
                for r in sub.itertuples(index=False):
                    row = [str(getattr(r, lig_c)), str(getattr(r, rec_c)), md_cell(float(getattr(r, score_col)))]
                    if p_c:
                        pv = getattr(r, p_c)
                        row.append(f"{float(pv):.3g}" if np.isfinite(pv) else "NA")
                    lines.append("| " + " | ".join(row) + " |")
                lines += [""]

    lines += [
        "## Readout",
        "",
        summary["trend_sentence"],
        "",
        "## Limits",
        "",
        "- TISCH2 labels and log-normalized matrix, not a GEO raw-UMI / CellRanger reprocess.",
        "- **No NK cells** are present under TISCH2 major-lineage in GSE148071; T/NK = CD8T + Tprolif.",
        "- Tprolif is a cycling T-cell bin, not a separate lineage proof.",
        "- Malignant is the TISCH2 call, not a re-run of inferCNV/CopyKAT.",
        "- Pooled LR ranks mix patients; cite the paired-n table for inference.",
        "- Chemokine dropout is high; read `ligand_frac` / `n_patients` with every rank.",
        "- CellChat was not run.",
        "",
        "## Files",
        "",
        "| File | Role |",
        "| --- | --- |",
        "| `results/n_cells_patients.tsv` | Honest per-patient cell counts |",
        "| `results/lr_cldn4high_to_tnk.tsv` | Full pooled LR table |",
        "| `results/lr_cldn4high_to_tnk_pass.tsv` | Pairs passing 10% expr_prop |",
        "| `results/patient_cldn4_outgoing.tsv` | Per-patient high vs low scores |",
        "| `results/ranks_cldn4_outgoing.tsv` | Patient-level Wilcoxon ranks |",
        "| `results/liana_cellphonedb.csv` | LIANA output if the import ran |",
        "| `results/summary.json` | Machine-readable n and method flags |",
        "| `results/figures/` | n-cell bars and top-pair plot |",
        "",
        f"Generated {summary['generated_at']}.",
        "",
    ]
    (outdir / "FINDING.md").write_text("\n".join(lines) + "\n")
    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    ap.add_argument("--outdir", type=Path, default=HERE / "results")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    P = cfg["params"]
    pairs = pd.read_csv(HERE / "resources" / "cellphonedb_v5_lr_pairs.tsv", sep="\t")
    # Re-label with this package's pathway overlay (adds checkpoint).
    pairs["pathway"] = [pathway_of(a, b, cfg) for a, b in zip(pairs["ligand"], pairs["receptor"])]

    wanted = sorted(
        set(Path(HERE / "resources" / "lr_genes.txt").read_text().split())
        | set(cfg["state_genes"])
        | {"CLDN4", "TACSTD2", "EPCAM", "PTPRC", "CD3D", "CD3E", "NKG7", "GNLY"}
    )
    # complexes in overlay
    for col in ("ligand", "receptor"):
        for s in pairs[col].astype(str):
            wanted.extend(s.split("+"))
    wanted = sorted(set(wanted))

    meta = pd.read_csv(args.datadir / "NSCLC_GSE148071_CellMetainfo_table.tsv", sep="\t", low_memory=False)
    cell_col = meta.columns[0]
    meta[cell_col] = meta[cell_col].astype(str)
    meta = meta.set_index(cell_col)
    lin_col = next(c for c in meta.columns if "major-lineage" in c.lower())
    meta["_lineage"] = meta[lin_col].astype(str).str.strip()
    meta["_patient"] = meta["Patient"].astype(str)
    tnk_set = set(cfg["tnk_lineages"])
    malig_set = set(cfg["malignant_lineages"])
    meta["_is_malig"] = meta["_lineage"].isin(malig_set)
    meta["_is_tnk"] = meta["_lineage"].isin(tnk_set)

    h5_path = args.datadir / "NSCLC_GSE148071_expression.h5"
    expr, present, h5_audit = load_tisch_genes(h5_path, wanted)
    expr.index = expr.index.astype(str)
    common = meta.index.intersection(expr.index)
    log(f"aligned cells={len(common)} / meta={len(meta)} / h5={len(expr)}")
    meta = meta.loc[common]
    expr = expr.loc[common]
    genes = [g for g in present if g in expr.columns]

    if "CLDN4" not in expr.columns:
        raise SystemExit("CLDN4 absent from TISCH2 h5")

    is_malig = meta["_is_malig"].to_numpy()
    is_tnk = meta["_is_tnk"].to_numpy()
    cldn4 = expr["CLDN4"].to_numpy(dtype=np.float32)
    thr = float(np.median(cldn4[is_malig]))
    is_high = is_malig & (cldn4 >= thr)
    is_low = is_malig & (cldn4 < thr)

    lineage_counts = meta["_lineage"].value_counts().to_dict()
    n_cd8t = int(lineage_counts.get("CD8T", 0))
    n_tprolif = int(lineage_counts.get("Tprolif", 0))
    n_nk = int(sum(v for k, v in lineage_counts.items() if k in {"NK", "NKT"}))

    log(
        f"malignant={is_malig.sum()} T/NK={is_tnk.sum()} "
        f"CLDN4-high={is_high.sum()} low={is_low.sum()} thr={thr:.3f}"
    )

    pooled = score_pairs(pairs, expr, is_high, is_tnk, genes, P["expr_prop"])
    pooled = pooled.sort_values("cpdb_mean_score", ascending=False)
    pooled.to_csv(args.outdir / "lr_cldn4high_to_tnk.tsv", sep="\t", index=False)
    passed = pooled[pooled["pass_expr_prop"]].copy()
    passed.to_csv(args.outdir / "lr_cldn4high_to_tnk_pass.tsv", sep="\t", index=False)
    log(f"pooled pairs scored={len(pooled)} pass={len(passed)}")

    # per-patient
    ntab_rows = []
    patient_high = []
    patient_low = []
    for pat, sub in meta.groupby("_patient", sort=True):
        idx = sub.index
        m_mal = meta.loc[idx, "_is_malig"].to_numpy()
        m_tnk = meta.loc[idx, "_is_tnk"].to_numpy()
        loc = meta.index.get_indexer(idx)
        hi = is_high[loc]
        lo = is_low[loc]
        ntab_rows.append(
            {
                "Patient": pat,
                "Sample": str(sub["Sample"].iloc[0]) if "Sample" in sub.columns else pat,
                "n_cells": int(len(sub)),
                "n_malignant": int(m_mal.sum()),
                "n_tnk": int(m_tnk.sum()),
                "n_CD8T": int((sub["_lineage"] == "CD8T").sum()),
                "n_Tprolif": int((sub["_lineage"] == "Tprolif").sum()),
                "n_NK": int(sub["_lineage"].isin(["NK", "NKT"]).sum()),
                "n_cldn4_high": int(hi.sum()),
                "n_cldn4_low": int(lo.sum()),
                "cldn4_mean_malignant": float(expr.loc[idx[m_mal], "CLDN4"].mean()) if m_mal.any() else np.nan,
                "usable_outgoing": bool(hi.sum() >= P["min_malig_high"] and m_tnk.sum() >= P["min_tnk"]),
                "usable_paired": bool(
                    hi.sum() >= P["min_malig_high"]
                    and lo.sum() >= P["min_malig_low"]
                    and m_tnk.sum() >= P["min_tnk"]
                ),
            }
        )
        full_mask_hi = np.zeros(len(meta), dtype=bool)
        full_mask_lo = np.zeros(len(meta), dtype=bool)
        full_mask_tnk = np.zeros(len(meta), dtype=bool)
        full_mask_hi[loc] = hi
        full_mask_lo[loc] = lo
        full_mask_tnk[loc] = m_tnk
        if hi.sum() >= P["min_malig_high"] and m_tnk.sum() >= P["min_tnk"]:
            hdf = score_pairs(pairs, expr, full_mask_hi, full_mask_tnk, genes, P["expr_prop"])
            hdf["Patient"] = pat
            patient_high.append(hdf)
        if (
            hi.sum() >= P["min_malig_high"]
            and lo.sum() >= P["min_malig_low"]
            and m_tnk.sum() >= P["min_tnk"]
        ):
            ldf = score_pairs(pairs, expr, full_mask_lo, full_mask_tnk, genes, P["expr_prop"])
            ldf["Patient"] = pat
            patient_low.append(ldf)

    ntab = pd.DataFrame(ntab_rows).sort_values("Patient")
    ntab.to_csv(args.outdir / "n_cells_patients.tsv", sep="\t", index=False)
    n_usable = int(ntab["usable_outgoing"].sum())
    n_paired = int(ntab["usable_paired"].sum())
    log(f"patients usable outgoing={n_usable} paired={n_paired}")

    ranks = pd.DataFrame()
    if patient_high and patient_low:
        high_all = pd.concat(patient_high, ignore_index=True)
        low_all = pd.concat(patient_low, ignore_index=True)
        high_all.to_csv(args.outdir / "patient_cldn4_high_outgoing.tsv", sep="\t", index=False)
        low_all.to_csv(args.outdir / "patient_cldn4_low_outgoing.tsv", sep="\t", index=False)
        key = ["Patient", "ligand", "receptor", "pathway"]
        merged = high_all.merge(
            low_all,
            on=key,
            suffixes=("_high", "_low"),
        )
        merged["delta_high_minus_low"] = merged["cpdb_mean_score_high"] - merged["cpdb_mean_score_low"]
        merged["pass_either"] = merged["pass_expr_prop_high"] | merged["pass_expr_prop_low"]
        merged.to_csv(args.outdir / "patient_cldn4_outgoing.tsv", sep="\t", index=False)

        rows = []
        for (lig, recp, path), g in merged.groupby(["ligand", "receptor", "pathway"], sort=False):
            use = g[g["pass_either"]]
            if len(use) < 3:
                continue
            pval = wilcoxon_safe(use["cpdb_mean_score_high"], use["cpdb_mean_score_low"])
            rows.append(
                {
                    "ligand": lig,
                    "receptor": recp,
                    "pathway": path,
                    "n_patients": int(len(use)),
                    "median_delta": float(use["delta_high_minus_low"].median()),
                    "mean_delta": float(use["delta_high_minus_low"].mean()),
                    "pval": pval,
                }
            )
        ranks = pd.DataFrame(rows)
        if len(ranks):
            mask = ranks["pval"].notna()
            ranks["padj"] = np.nan
            if mask.sum():
                ranks.loc[mask, "padj"] = multipletests(ranks.loc[mask, "pval"], method="fdr_bh")[1]
            ranks = ranks.sort_values(["padj", "pval", "median_delta"], na_position="last")
            ranks.to_csv(args.outdir / "ranks_cldn4_outgoing.tsv", sep="\t", index=False)

    # LIANA secondary
    liana_note = "not attempted"
    labels = pd.Series(index=meta.index, dtype=object)
    labels[is_high] = "Malig_CLDN4high"
    labels[is_low] = "Malig_CLDN4low"
    labels[is_tnk] = "TNK"
    keep = labels.notna()
    liana_note = run_liana(
        expr.loc[keep, genes],
        labels.loc[keep],
        args.outdir / "liana_cellphonedb.csv",
        n_perms=int(P["liana_n_perms"]),
        max_per_group=int(P["liana_max_cells_per_group"]),
        seed=int(P["random_seed"]),
    )
    log(liana_note)
    liana_short = liana_note.split()[0]

    # figures
    figdir = args.outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(9.2, 3.8))
    x = np.arange(len(ntab))
    ax.bar(x - 0.2, ntab["n_cldn4_high"], width=0.4, label="CLDN4-high malignant", color="#4C72B0")
    ax.bar(x + 0.2, ntab["n_tnk"], width=0.4, label="T/NK (CD8T+Tprolif)", color="#DD8452")
    ax.set_xticks(x)
    ax.set_xticklabels(ntab["Patient"], rotation=90, fontsize=6)
    ax.set_ylabel("cells")
    ax.set_title("GSE148071 TISCH2: CLDN4-high malignant and T/NK per patient")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(figdir / "n_cells_by_patient.png", dpi=150)
    fig.savefig(figdir / "n_cells_by_patient.pdf")
    plt.close(fig)

    if len(passed):
        top = passed.head(20).iloc[::-1]
        fig, ax = plt.subplots(figsize=(8.0, 5.6))
        ax.barh(np.arange(len(top)), top["cpdb_mean_score"], color="#4C72B0")
        ax.set_yticks(np.arange(len(top)))
        ax.set_yticklabels([f"{a}–{b}" for a, b in zip(top["ligand"], top["receptor"])], fontsize=7)
        ax.set_xlabel("CPDB mean score (TISCH2 log2(TPM/10+1))")
        ax.set_title("Top pooled LR: CLDN4-high malignant → T/NK")
        fig.tight_layout()
        fig.savefig(figdir / "top_lr_cldn4high_to_tnk.png", dpi=150)
        fig.savefig(figdir / "top_lr_cldn4high_to_tnk.pdf")
        plt.close(fig)

    # trend sentence from data
    n_pass = int(len(passed))
    top_names = ", ".join(f"{a}–{b}" for a, b in zip(passed.head(5)["ligand"], passed.head(5)["receptor"])) if n_pass else "none"
    focus_pass = passed[passed["pathway"].isin(["T_recruit", "checkpoint", "MHC_I"])]
    if n_paired >= 3 and len(ranks):
        foc_r = ranks[ranks["pathway"].isin(["T_recruit", "checkpoint", "MHC_I"])]
        n_sig = int((foc_r["padj"] < 0.05).sum()) if len(foc_r) and foc_r["padj"].notna().any() else 0
        trend = (
            f"Pooled, {n_pass} pairs pass 10% expression from CLDN4-high malignant to T/NK "
            f"(top: {top_names}). Patient-level high vs low uses n={n_paired} patients; "
            f"{n_sig} focus-axis pairs reach FDR < 0.05. T/NK in this object is CD8T+Tprolif "
            f"(NK=0). This is the rank in public TISCH2 GSE148071, not a CellChat run."
        )
    else:
        trend = (
            f"Pooled, {n_pass} pairs pass 10% expression from CLDN4-high malignant to T/NK "
            f"(top: {top_names}). Paired high vs low n={n_paired} is too small for a stable "
            f"Wilcoxon on most pairs. T/NK = CD8T+Tprolif (NK=0)."
        )

    summary = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "accession": "GSE148071",
        "tisch_id": "NSCLC_GSE148071",
        "pmid": "33953163",
        "normalization": "TISCH2/MAESTRO log2(TPM/10+1)",
        "malignant_definition": "TISCH2 major-lineage == Malignant",
        "tnk_definition": "TISCH2 major-lineage in CD8T/Tprolif (NK absent)",
        "n_cells": int(len(meta)),
        "n_patients_meta": int(meta["_patient"].nunique()),
        "n_malignant": int(is_malig.sum()),
        "n_tnk": int(is_tnk.sum()),
        "n_cd8t": n_cd8t,
        "n_tprolif": n_tprolif,
        "n_nk": n_nk,
        "n_cldn4_high": int(is_high.sum()),
        "n_cldn4_low": int(is_low.sum()),
        "cldn4_threshold": thr,
        "min_malig_high": int(P["min_malig_high"]),
        "min_malig_low": int(P["min_malig_low"]),
        "min_tnk": int(P["min_tnk"]),
        "n_patients_usable": n_usable,
        "n_patients_paired": n_paired,
        "n_patients_dropped": int((~ntab["usable_outgoing"]).sum()),
        "n_pairs_scored": int(len(pooled)),
        "n_pairs_pass": n_pass,
        "liana_status": liana_note,
        "liana_status_short": liana_short,
        "cellchat_status": "not run (R unavailable)",
        "h5": h5_audit,
        "lineage_counts": lineage_counts,
        "top_table_n": int(P["top_table_n"]),
        "trend_sentence": trend,
        "usable_patients": ntab.loc[ntab["usable_outgoing"], "Patient"].tolist(),
        "paired_patients": ntab.loc[ntab["usable_paired"], "Patient"].tolist(),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2))
    write_finding(args.outdir, summary, ntab, pooled, ranks, liana_note)
    log(f"wrote FINDING.md and {args.outdir}")


if __name__ == "__main__":
    main()
