"""Common analysis engine shared by the per-dataset scripts."""
import json
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import gene_panels as GP
import stats_utils as SU

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results" / "fable_tacstd2_kdko"
DATA = RES / "data"
RES.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- mygene mapping
def map_ensembl_to_symbol(ensembl_ids, species, cache_name):
    """Map Ensembl gene IDs -> symbol via mygene.info (batched POST). Cached."""
    cache = DATA / cache_name
    if cache.exists():
        return json.loads(cache.read_text())
    ids = [e.split(".")[0] for e in ensembl_ids]
    out = {}
    url = "https://mygene.info/v3/query"
    for i in range(0, len(ids), 900):
        chunk = ids[i:i + 900]
        body = ("q=" + ",".join(chunk) +
                "&scopes=ensembl.gene&fields=symbol&species=" + species).encode()
        req = urllib.request.Request(url, data=body,
                                     headers={"Content-Type": "application/x-www-form-urlencoded"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=60) as r:
                    res = json.loads(r.read().decode())
                break
            except Exception:
                if attempt == 3:
                    raise
                time.sleep(2 ** attempt)
        for rec in res:
            q = rec.get("query")
            sym = rec.get("symbol")
            if q and sym and q not in out:
                out[q] = sym
        time.sleep(0.3)
    cache.write_text(json.dumps(out))
    return out


# ------------------------------------------------------- gene-level collapse
def collapse_to_gene(symbols, logmat):
    """Collapse multiple rows per symbol -> one row (highest mean expression)."""
    symbols = np.asarray([str(s) for s in symbols])
    mean_expr = logmat.mean(axis=1)
    order = np.argsort(-mean_expr)
    seen = {}
    keep_idx = []
    for idx in order:
        s = symbols[idx].upper()
        if s in ("", "NAN", "NONE"):
            continue
        if s not in seen:
            seen[s] = idx
            keep_idx.append(idx)
    keep_idx = np.array(sorted(keep_idx))
    return symbols[keep_idx], logmat[keep_idx], keep_idx


# ------------------------------------------------------- panel helpers
def _panel_mask(symbols_upper, panel):
    pset = set(GP.upper(x) for x in panel)
    return np.array([s in pset for s in symbols_upper])


def run_expression_dataset(dataset_id, description, symbols, logmat,
                           case_idx, ctrl_idx, value_kind, notes=""):
    """Full analysis for a two-group expression matrix (log scale).

    Writes:
      <dataset>_pergene_panel.tsv   panel + perturbation gene per-gene stats
      <dataset>_setlevel.tsv        gene-set competitive & self-contained tests
    Returns a summary dict.
    """
    symbols_up = np.array([GP.upper(s) for s in symbols])
    stats = SU.per_gene_two_group(logmat, case_idx, ctrl_idx)

    # ranking stat for competitive test: signed -log10(p) * sign(log2FC)
    signed = np.sign(stats["log2FC"]) * -np.log10(np.clip(stats["t_pvalue"], 1e-300, 1))

    # ---- per-gene panel table ----
    panel_of = {}
    for pname, genes in GP.PANELS.items():
        for g in genes:
            panel_of.setdefault(GP.upper(g), []).append(pname)
    for g in GP.PERTURBATION:
        panel_of.setdefault(GP.upper(g), []).append("perturbation")
    panel_of.setdefault("CLDN4", []).append("headline")

    rows = []
    wanted = set(panel_of)
    for i, s in enumerate(symbols_up):
        if s in wanted:
            rows.append({
                "dataset": dataset_id,
                "gene": s,
                "panels": ";".join(sorted(set(panel_of[s]))),
                "log2FC": stats["log2FC"][i],
                "mean_case": stats["mean_case"][i],
                "mean_ctrl": stats["mean_ctrl"][i],
                "t_pvalue": stats["t_pvalue"][i],
                "t_fdr": stats["t_fdr"][i],
                "mwu_pvalue": stats["mwu_pvalue"][i],
                "cohens_d": stats["cohens_d"][i],
            })
    pergene = pd.DataFrame(rows).sort_values(["panels", "gene"]).reset_index(drop=True)
    pergene.to_csv(RES / f"{dataset_id}_pergene_panel.tsv", sep="\t", index=False,
                   float_format="%.5g")

    # ---- set-level tests (competitive on genome-wide signed stat; self-contained on log2FC) ----
    set_rows = []
    for pname, genes in list(GP.PANELS.items()) + list(GP.AXES.items()):
        mask = _panel_mask(symbols_up, genes)
        comp = SU.competitive_set_test(signed, mask)
        sc = SU.selfcontained_set_test(stats["log2FC"][mask])
        set_rows.append({
            "dataset": dataset_id,
            "gene_set": pname,
            "n_detected": int(mask.sum()),
            "comp_direction": comp["direction"],
            "comp_rank_biserial": comp["rank_biserial"],
            "comp_p": comp["p_two_sided"],
            "sc_mean_log2FC": sc.get("mean_log2FC"),
            "sc_median_log2FC": sc.get("median_log2FC"),
            "sc_n_up": sc.get("n_up"),
            "sc_n_down": sc.get("n_down"),
            "sc_wilcoxon_p": sc.get("wilcoxon_p"),
        })
    setdf = pd.DataFrame(set_rows)
    setdf.to_csv(RES / f"{dataset_id}_setlevel.tsv", sep="\t", index=False,
                 float_format="%.4g")

    # ---- QC: perturbation knockdown ----
    qc = {}
    for g in GP.PERTURBATION:
        gi = np.where(symbols_up == g)[0]
        if len(gi):
            i = gi[0]
            qc[g] = {"log2FC": float(stats["log2FC"][i]),
                     "t_pvalue": float(stats["t_pvalue"][i]),
                     "t_fdr": float(stats["t_fdr"][i])}

    summary = {
        "dataset": dataset_id,
        "description": description,
        "value_kind": value_kind,
        "n_case": len(case_idx),
        "n_ctrl": len(ctrl_idx),
        "n_genes_tested": int(logmat.shape[0]),
        "notes": notes,
        "perturbation_qc": qc,
        "cldn4": next((r for r in rows if r["gene"] == "CLDN4"), None),
    }
    (RES / f"{dataset_id}_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return summary, pergene, setdf


def run_deseq2_table(dataset_id, description, df, symbol_col, l2fc_col,
                     p_col, padj_col, stat_col, basemean_col, notes=""):
    """Analysis for an author-provided DESeq2 result table (contrast = KO vs WT)."""
    df = df.copy()
    df["_sym"] = df[symbol_col].astype(str).str.upper().str.replace("_", "-")
    # collapse duplicate symbols -> most expressed (max baseMean)
    df = df.sort_values(basemean_col, ascending=False).drop_duplicates("_sym")
    symbols_up = df["_sym"].values
    l2fc = df[l2fc_col].astype(float).values
    pval = df[p_col].astype(float).fillna(1.0).values
    padj = df[padj_col].astype(float).fillna(1.0).values
    stat = df[stat_col].astype(float).fillna(0.0).values

    panel_of = {}
    for pname, genes in GP.PANELS.items():
        for g in genes:
            panel_of.setdefault(GP.upper(g), []).append(pname)
    for g in GP.PERTURBATION:
        panel_of.setdefault(GP.upper(g), []).append("perturbation")
    panel_of.setdefault("CLDN4", []).append("headline")

    rows = []
    for i, s in enumerate(symbols_up):
        if s in panel_of:
            rows.append({
                "dataset": dataset_id, "gene": s,
                "panels": ";".join(sorted(set(panel_of[s]))),
                "log2FC": l2fc[i], "pvalue": pval[i], "padj": padj[i],
                "wald_stat": stat[i],
            })
    pergene = pd.DataFrame(rows).sort_values(["panels", "gene"]).reset_index(drop=True)
    pergene.to_csv(RES / f"{dataset_id}_pergene_panel.tsv", sep="\t", index=False,
                   float_format="%.5g")

    set_rows = []
    for pname, genes in list(GP.PANELS.items()) + list(GP.AXES.items()):
        mask = _panel_mask(symbols_up, genes)
        comp = SU.competitive_set_test(stat, mask)
        sc = SU.selfcontained_set_test(l2fc[mask])
        set_rows.append({
            "dataset": dataset_id, "gene_set": pname, "n_detected": int(mask.sum()),
            "comp_direction": comp["direction"], "comp_rank_biserial": comp["rank_biserial"],
            "comp_p": comp["p_two_sided"], "sc_mean_log2FC": sc.get("mean_log2FC"),
            "sc_median_log2FC": sc.get("median_log2FC"), "sc_n_up": sc.get("n_up"),
            "sc_n_down": sc.get("n_down"), "sc_wilcoxon_p": sc.get("wilcoxon_p"),
        })
    setdf = pd.DataFrame(set_rows)
    setdf.to_csv(RES / f"{dataset_id}_setlevel.tsv", sep="\t", index=False, float_format="%.4g")

    qc = {}
    for g in GP.PERTURBATION:
        gi = np.where(symbols_up == g)[0]
        if len(gi):
            i = gi[0]
            qc[g] = {"log2FC": float(l2fc[i]), "pvalue": float(pval[i]), "padj": float(padj[i])}

    summary = {
        "dataset": dataset_id, "description": description, "value_kind": "author DESeq2",
        "n_genes_tested": int(len(df)), "notes": notes, "perturbation_qc": qc,
        "cldn4": next((r for r in rows if r["gene"] == "CLDN4"), None),
    }
    (RES / f"{dataset_id}_summary.json").write_text(json.dumps(summary, indent=2, default=float))
    return summary, pergene, setdf
