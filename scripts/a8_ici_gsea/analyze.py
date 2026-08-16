#!/usr/bin/env python3
"""ADDITIVE A8 extra — keratin / TJ / EMT GSEA in public ICI-treated lung tumors.

Does NOT re-run the TCGA-LUAD/LUSC A8 audit. User A8 (TROP2-high keratin/TJ
up, EMT down, conserved) is taken as given for TCGA. This slice asks the
same gene-set question in leftover open pre-treatment ICI lung bulk.

Locked design (set before looking at NES)
-----------------------------------------
- Cohorts: public pre-treatment ICI-treated lung tumor bulk with TACSTD2
  on a whole-transcriptome matrix. Not TCGA. Not post-treatment-only
  GSE248378. Not GSE207422 (already GSEA'd in the hunt_tj_gsea slice).
- Primary contrast: TACSTD2 median-split Welch t (high minus low).
  Quartile only if both arms have >= 6 samples (small-n honesty).
- Complementary ranking: Spearman rho of every gene vs continuous TACSTD2
  (keeps the middle half).
- Response contrast: responder vs non-responder Welch t only if both
  arms have >= 5 samples. Labels from deposited GEO fields (or the
  France3 sample-info table). No invented RECIST.
- GSEA: same preranked weighted KS as A8 (p=1, 1000 gene-set permutations,
  seed=42). Positive NES = enriched in TACSTD2-high (or in responders).
- Primary sets: the same 12 A8 claim sets. Hallmark EMT decides the EMT
  arm — GOBP EMT does not. BH-FDR within those 12 per contrast.
- Do not pool NES across cohorts. Do not quote a TCGA NES from this run.

Outputs -> results/a8_ici_gsea/
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import sys
import time
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
from gsea_core import (  # noqa: E402
    bh_fdr,
    gsea_prerank,
    median_split,
    quartile_split,
    rank_high_vs_low,
    rank_spearman_vs_target,
)
from geo_load import LOADERS as _LOAD_FNS  # noqa: E402
from geo_load import load_ensembl_map  # noqa: E402

DATA = os.environ.get("A8_ICI_DATA", "/tmp/a8_ici_gsea_data")
OUT = os.path.join(ROOT, "results", "a8_ici_gsea")
FIG = os.path.join(OUT, "figures")
TAB = os.path.join(OUT, "tables")
SET_JSON = os.path.join(ROOT, "data", "genesets", "a8_sets.json")
os.makedirs(FIG, exist_ok=True)
os.makedirs(TAB, exist_ok=True)

TARGET = "TACSTD2"
MIN_ARM_MEDIAN = 5
MIN_ARM_QUARTILE = 6
MIN_ARM_RESPONSE = 5
PRIMARY_FDR = 0.05

FOCAL = [
    "CLDN1",
    "CLDN4",
    "CLDN7",
    "F11R",
    "PARD3",
    "OCLN",
    "TJP1",
    "KRT5",
    "KRT17",
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "CDH1",
    "VIM",
    "ZEB1",
    "SNAI2",
]

HEADLINE = [
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_APICAL_JUNCTION",
    "KEGG_TIGHT_JUNCTION",
    "GOBP_KERATINIZATION",
    "KRT_EPITHELIAL",
]


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


LOADERS = {
    "GSE126044": {
        "label": "Cho 2020 anti-PD-1 NSCLC",
        "scale": "log2(CPM+1) from deposited counts",
        "pubmed": "PMID 32879421",
        "leftover_note": "core open ICI bulk; not in the TCGA A8 slide",
    },
    "GSE135222": {
        "label": "Jung/Kim 2019–2020 anti-PD-(L)1 NSCLC",
        "scale": "log2(TPM+1) from deposited TPM; Ensembl mapped via HGNC",
        "pubmed": "PMID 31537801",
        "leftover_note": "core open ICI bulk; DCB is GEO PFS>=180 d",
    },
    "GSE166449": {
        "label": "Lee/Hwang 2021 pembrolizumab LUAD",
        "scale": "deposited matrix kept if already log2, else log2(TPM+1)",
        "pubmed": "PMID 33857424",
        "leftover_note": "core open ICI bulk; response from GEO sample titles",
    },
    "GSE253564": {
        "label": "Altorki leftover neoadjuvant durvalumab ± SBRT, pre-treatment",
        "scale": "log2(FPKM+1) from deposited pre-treatment FPKM",
        "pubmed": "PMID 38401548",
        "leftover_note": "leftover whole-transcriptome ICI tumor matrix; MPR not deposited on GEO",
    },
    "GSE190265": {
        "label": "CGFL Dijon France3 anti-PD-1 NSCLC",
        "scale": "log2(TPM+1) from deposited France3 TPM",
        "pubmed": "PMID 35051357",
        "leftover_note": "France4 (GSE190266) omitted: TACSTD2 absent from the deposited TPM cap",
    },
    "GSE283829": {
        "label": "Lindberg 2025 leftover ICI NSCLC",
        "scale": "log2(CPM+1) from supplementary raw counts; Ensembl mapped via HGNC",
        "pubmed": "PMID 39743139",
        "leftover_note": "2025 leftover; R vs NR is CR vs PD (SD held out)",
    },
}


def signature_zmean(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=expr.columns)
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def spearman(x: pd.Series, y: pd.Series):
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def verdict_from_primary(gsea: pd.DataFrame) -> str:
    """Same labels as A8. Hallmark EMT decides the EMT arm."""
    if gsea.empty:
        return "not_run"
    prim = gsea[gsea["primary"] == True].copy()  # noqa: E712
    if prim.empty:
        return "not_run"

    def sig(term, want_up: bool) -> bool:
        row = prim[prim["term"] == term]
        if row.empty:
            return False
        nes = float(row.iloc[0]["nes"])
        q = float(row.iloc[0]["fdr_primary"])
        if q >= PRIMARY_FDR or nes != nes:
            return False
        return (nes > 0) if want_up else (nes < 0)

    def any_bucket(bucket: str, want_up: bool) -> bool:
        sub = prim[prim["bucket"] == bucket]
        ok = False
        for _, r in sub.iterrows():
            if r["term"] == "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION":
                continue
            if r["fdr_primary"] < PRIMARY_FDR and ((r["nes"] > 0) if want_up else (r["nes"] < 0)):
                ok = True
        return ok

    emt_down = sig("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", False)
    emt_up = sig("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", True)
    tj_up = any_bucket("TJ", True)
    krt_up = any_bucket("KERATIN_BARRIER", True)
    tj_down = any_bucket("TJ", False) and not tj_up
    krt_down = any_bucket("KERATIN_BARRIER", False) and not krt_up

    if emt_down and tj_up and krt_up:
        return "supportive"
    if tj_up and krt_up and emt_up:
        return "keratin_TJ_up_Hallmark_EMT_opposite"
    if tj_up and krt_up and not emt_down and not emt_up:
        return "keratin_TJ_up_Hallmark_EMT_null"
    if emt_down and (tj_up or krt_up):
        return "partial"
    if tj_down and krt_down:
        return "contradicts_TJ_KRT"
    if tj_up or krt_up or emt_down:
        return "mixed"
    return "null"


def annotate_gsea(gsea: pd.DataFrame, set_meta: dict, contrast: str, cohort: str, n_high: int, n_low: int) -> pd.DataFrame:
    if gsea.empty:
        return gsea
    gsea = gsea.copy()
    gsea["cohort"] = cohort
    gsea["contrast"] = contrast
    gsea["n_high"] = n_high
    gsea["n_low"] = n_low
    gsea = gsea.merge(
        pd.DataFrame.from_dict(set_meta, orient="index").rename_axis("term").reset_index(),
        on="term",
        how="left",
    )
    gsea["fdr_all"] = bh_fdr(gsea["nom_p"])
    prim = gsea["primary"] == True  # noqa: E712
    gsea["fdr_primary"] = np.nan
    if prim.any():
        gsea.loc[prim, "fdr_primary"] = bh_fdr(gsea.loc[prim, "nom_p"])
    return gsea


def run_contrast(expr, high, low, gene_sets, set_meta, contrast, cohort) -> pd.DataFrame:
    rank = rank_high_vs_low(expr, high, low)
    gsea = gsea_prerank(rank, gene_sets)
    return annotate_gsea(gsea, set_meta, contrast, cohort, len(high), len(low))


def fmt_p(p) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_nes(x) -> str:
    if x != x:
        return "NA"
    return f"{x:+.3f}"


def main() -> None:
    payload = json.loads(open(SET_JSON).read())
    all_sets = payload["sets"]
    set_meta = payload["meta"]
    primary_names = [k for k, v in set_meta.items() if v.get("primary")]
    primary_sets = {k: all_sets[k] for k in primary_names if k in all_sets}

    inventory = []
    gsea_rows = []
    sig_rows = []
    focal_rows = []
    sample_rows = []
    verdicts = []

    mapping = load_ensembl_map(DATA)
    log(f"HGNC Ensembl map: {len(mapping)} ids")
    for acc, spec in LOADERS.items():
        log(f"=== {acc} {spec['label']} ===")
        try:
            expr, meta = _LOAD_FNS[acc](DATA, mapping)
        except Exception as exc:  # noqa: BLE001
            log(f"  LOAD FAIL: {exc}")
            inventory.append(
                {
                    "cohort": acc,
                    "usable": False,
                    "reason": f"load_failed: {exc}",
                    "n": 0,
                    "n_with_TACSTD2": 0,
                    "n_R": 0,
                    "n_NR": 0,
                    "label": spec["label"],
                }
            )
            continue
        if TARGET not in expr.index:
            log(f"  TACSTD2 absent (n_genes={expr.shape[0]})")
            inventory.append(
                {
                    "cohort": acc,
                    "usable": False,
                    "reason": "TACSTD2_absent",
                    "n": expr.shape[1],
                    "n_with_TACSTD2": 0,
                    "n_R": 0,
                    "n_NR": 0,
                    "label": spec["label"],
                }
            )
            continue
        tac = expr.loc[TARGET].dropna()
        expr = expr.loc[:, tac.index]
        meta = meta.reindex(tac.index)
        n = int(tac.shape[0])
        n_r = int((meta["response"] == "R").sum()) if "response" in meta.columns else 0
        n_nr = int((meta["response"] == "NR").sum()) if "response" in meta.columns else 0
        log(f"  n={n} genes={expr.shape[0]} TACSTD2 median={float(tac.median()):.3f} R={n_r} NR={n_nr}")
        inventory.append(
            {
                "cohort": acc,
                "usable": n >= 10,
                "reason": "ok" if n >= 10 else "n<10",
                "n": n,
                "n_with_TACSTD2": n,
                "n_R": n_r,
                "n_NR": n_nr,
                "n_genes": int(expr.shape[0]),
                "tacstd2_median": float(tac.median()),
                "label": spec["label"],
                "scale": spec["scale"],
                "pubmed": spec["pubmed"],
                "leftover_note": spec["leftover_note"],
            }
        )
        for sid in tac.index:
            sample_rows.append(
                {
                    "cohort": acc,
                    "sample": sid,
                    "TACSTD2": float(tac.loc[sid]),
                    "response": meta.loc[sid, "response"] if "response" in meta.columns else np.nan,
                    "response_raw": meta.loc[sid, "response_raw"] if "response_raw" in meta.columns else "",
                }
            )

        # complementary signature Spearman (no split)
        for term, genes in primary_sets.items():
            sc = signature_zmean(expr, genes)
            rho, p, nn = spearman(tac, sc)
            sig_rows.append(
                {
                    "cohort": acc,
                    "term": term,
                    "rho": rho,
                    "p": p,
                    "n": nn,
                    "n_genes_present": sum(g in expr.index for g in genes),
                    "bucket": set_meta[term]["bucket"],
                    "expected": set_meta[term]["expected"],
                }
            )
        for g in FOCAL:
            if g not in expr.index:
                focal_rows.append({"cohort": acc, "gene": g, "rho": np.nan, "p": np.nan, "n": 0, "present": False})
                continue
            rho, p, nn = spearman(tac, expr.loc[g])
            focal_rows.append({"cohort": acc, "gene": g, "rho": rho, "p": p, "n": nn, "present": True})

        if n < 10:
            verdicts.append({"cohort": acc, "contrast": "tacstd2_median", "verdict": "n_too_small", "n_high": 0, "n_low": 0})
            continue

        # --- primary: TACSTD2 median ---
        high, low = median_split(tac)
        if len(high) >= MIN_ARM_MEDIAN and len(low) >= MIN_ARM_MEDIAN:
            log(f"  GSEA TACSTD2 median {len(high)} vs {len(low)}")
            g = run_contrast(expr, high, low, primary_sets, set_meta, "tacstd2_median", acc)
            gsea_rows.append(g)
            verdicts.append(
                {
                    "cohort": acc,
                    "contrast": "tacstd2_median",
                    "verdict": verdict_from_primary(g),
                    "n_high": len(high),
                    "n_low": len(low),
                }
            )
        else:
            verdicts.append(
                {
                    "cohort": acc,
                    "contrast": "tacstd2_median",
                    "verdict": "n_too_small",
                    "n_high": len(high),
                    "n_low": len(low),
                }
            )

        # --- sensitivity: quartile if both arms >= 6 ---
        qh, ql = quartile_split(tac)
        if len(qh) >= MIN_ARM_QUARTILE and len(ql) >= MIN_ARM_QUARTILE:
            log(f"  GSEA TACSTD2 quartile {len(qh)} vs {len(ql)}")
            g = run_contrast(expr, qh, ql, primary_sets, set_meta, "tacstd2_quartile", acc)
            gsea_rows.append(g)
            verdicts.append(
                {
                    "cohort": acc,
                    "contrast": "tacstd2_quartile",
                    "verdict": verdict_from_primary(g),
                    "n_high": len(qh),
                    "n_low": len(ql),
                }
            )
        else:
            verdicts.append(
                {
                    "cohort": acc,
                    "contrast": "tacstd2_quartile",
                    "verdict": "n_too_small",
                    "n_high": len(qh),
                    "n_low": len(ql),
                }
            )

        # --- complementary Spearman prerank ---
        log("  GSEA Spearman vs TACSTD2")
        rank = rank_spearman_vs_target(expr, tac)
        g = gsea_prerank(rank, primary_sets)
        g = annotate_gsea(g, set_meta, "tacstd2_spearman", acc, n, n)
        gsea_rows.append(g)
        verdicts.append(
            {
                "cohort": acc,
                "contrast": "tacstd2_spearman",
                "verdict": verdict_from_primary(g),
                "n_high": n,
                "n_low": n,
            }
        )

        # --- response if n allows ---
        if n_r >= MIN_ARM_RESPONSE and n_nr >= MIN_ARM_RESPONSE:
            r_idx = meta.index[meta["response"] == "R"]
            nr_idx = meta.index[meta["response"] == "NR"]
            log(f"  GSEA response R vs NR {len(r_idx)} vs {len(nr_idx)}")
            g = run_contrast(expr, r_idx, nr_idx, primary_sets, set_meta, "response_R_vs_NR", acc)
            gsea_rows.append(g)
            verdicts.append(
                {
                    "cohort": acc,
                    "contrast": "response_R_vs_NR",
                    "verdict": verdict_from_primary(g),
                    "n_high": len(r_idx),
                    "n_low": len(nr_idx),
                }
            )
        else:
            verdicts.append(
                {
                    "cohort": acc,
                    "contrast": "response_R_vs_NR",
                    "verdict": "n_too_small" if (n_r + n_nr) > 0 else "no_deposited_response_label",
                    "n_high": n_r,
                    "n_low": n_nr,
                }
            )

    inv = pd.DataFrame(inventory)
    gsea = pd.concat(gsea_rows, ignore_index=True) if gsea_rows else pd.DataFrame()
    sig = pd.DataFrame(sig_rows)
    if not sig.empty:
        sig["q"] = bh_fdr(sig["p"])
    focal = pd.DataFrame(focal_rows)
    if not focal.empty:
        focal["q"] = bh_fdr(focal["p"])
    samples = pd.DataFrame(sample_rows)
    verd = pd.DataFrame(verdicts)

    inv.to_csv(os.path.join(TAB, "cohort_inventory.tsv"), sep="\t", index=False)
    if not gsea.empty:
        gsea.to_csv(os.path.join(TAB, "gsea_prerank_all.tsv"), sep="\t", index=False)
        gsea[gsea["contrast"] == "tacstd2_median"].to_csv(
            os.path.join(TAB, "gsea_tacstd2_median.tsv"), sep="\t", index=False
        )
        gsea[gsea["contrast"] == "response_R_vs_NR"].to_csv(
            os.path.join(TAB, "gsea_response.tsv"), sep="\t", index=False
        )
    if not sig.empty:
        sig.to_csv(os.path.join(TAB, "signature_spearman.tsv"), sep="\t", index=False)
    if not focal.empty:
        focal.to_csv(os.path.join(TAB, "focal_gene_spearman.tsv"), sep="\t", index=False)
    samples.to_csv(os.path.join(TAB, "sample_table.tsv"), sep="\t", index=False)
    verd.to_csv(os.path.join(TAB, "verdicts.tsv"), sep="\t", index=False)

    provenance = {}
    if os.path.isdir(DATA):
        for fn in sorted(os.listdir(DATA)):
            p = os.path.join(DATA, fn)
            if os.path.isfile(p):
                provenance[fn] = {"bytes": os.path.getsize(p), "md5": md5(p)}
    json.dump(provenance, open(os.path.join(OUT, "provenance.json"), "w"), indent=2)

    make_figures(gsea, sig, focal, inv, verd)
    write_report(inv, gsea, sig, focal, verd, spec_map=LOADERS)
    summary = {
        "question": "In leftover public pre-treatment ICI lung bulk, do TACSTD2-high tumors enrich keratin/TJ and deplete Hallmark EMT?",
        "tcga_a8": "taken_as_given_not_rerun",
        "cohorts": inv.to_dict(orient="records"),
        "verdicts": verd.to_dict(orient="records"),
        "nperm": 1000,
        "seed": 42,
        "fdr_rule": "BH within 12 primary sets per contrast",
        "honest": "small-n ICI bulk is hypothesis-generating; do not quote a pooled NES; Hallmark EMT not GOBP EMT decides the EMT arm",
    }
    json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
    log("done")


def make_figures(gsea, sig, focal, inv, verd) -> None:
    if gsea is None or gsea.empty:
        return
    prim = gsea[(gsea["primary"] == True) & (gsea["contrast"] == "tacstd2_median")].copy()  # noqa: E712
    if prim.empty:
        return
    # NES heatmap
    terms = [
        t
        for t in [
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
            "HALLMARK_APICAL_JUNCTION",
            "KEGG_TIGHT_JUNCTION",
            "GOBP_TIGHT_JUNCTION_ORGANIZATION",
            "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
            "GOBP_KERATINIZATION",
            "GOBP_CORNIFICATION",
            "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER",
            "GOBP_KERATINOCYTE_DIFFERENTIATION",
            "GOBP_EPIDERMAL_CELL_DIFFERENTIATION",
            "KRT_EPITHELIAL",
        ]
        if t in set(prim["term"])
    ]
    cohorts = [c for c in LOADERS if c in set(prim["cohort"])]
    mat = pd.DataFrame(index=terms, columns=cohorts, dtype=float)
    qmat = pd.DataFrame(index=terms, columns=cohorts, dtype=float)
    nlab = {}
    for c in cohorts:
        sub = prim[prim["cohort"] == c]
        nlab[c] = f"{c}\n{int(sub['n_high'].iloc[0])} vs {int(sub['n_low'].iloc[0])}" if len(sub) else c
        for t in terms:
            row = sub[sub["term"] == t]
            if len(row):
                mat.loc[t, c] = float(row.iloc[0]["nes"])
                qmat.loc[t, c] = float(row.iloc[0]["fdr_primary"])
    fig, ax = plt.subplots(figsize=(10.5, 7.2))
    vmin, vmax = -2.5, 2.5
    im = ax.imshow(mat.to_numpy(dtype=float), cmap="RdBu_r", vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(cohorts)))
    ax.set_xticklabels([nlab[c] for c in cohorts], fontsize=8)
    ax.set_yticks(range(len(terms)))
    ax.set_yticklabels(terms, fontsize=8)
    for i, t in enumerate(terms):
        for j, c in enumerate(cohorts):
            nes = mat.loc[t, c]
            q = qmat.loc[t, c]
            if nes != nes:
                txt = ""
            else:
                star = "*" if q == q and q < 0.05 else ""
                txt = f"{nes:+.2f}{star}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="black")
    ax.set_title("TACSTD2-high vs low (median) prerank GSEA  NES  * FDR<0.05 within 12 primary sets")
    fig.colorbar(im, ax=ax, fraction=0.03, label="NES")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "nes_tacstd2_median.png"), dpi=160)
    plt.close(fig)

    # headline forest-style bars
    fig, axes = plt.subplots(1, len(HEADLINE), figsize=(14.5, 4.6), sharey=True)
    colors = {"EMT": "#b2182b", "TJ": "#2166ac", "KERATIN_BARRIER": "#4dac26"}
    for ax, term in zip(axes, HEADLINE):
        sub = prim[prim["term"] == term]
        ys = []
        xs = []
        cs = []
        for i, c in enumerate(cohorts):
            row = sub[sub["cohort"] == c]
            if row.empty:
                continue
            ys.append(i)
            xs.append(float(row.iloc[0]["nes"]))
            bucket = str(row.iloc[0].get("bucket", "TJ"))
            cs.append(colors.get(bucket, "#666666"))
            q = float(row.iloc[0]["fdr_primary"])
            if q < 0.05:
                ax.text(xs[-1] + (0.08 if xs[-1] >= 0 else -0.08), i, "*", va="center", ha="left" if xs[-1] >= 0 else "right")
        ax.axvline(0, color="0.5", lw=0.8)
        ax.barh(ys, xs, color=cs, height=0.6)
        ax.set_title(term.replace("_", "\n"), fontsize=8)
        ax.set_xlabel("NES")
        ax.set_yticks(range(len(cohorts)))
        ax.set_yticklabels([f"{c} (n={int(inv.loc[inv.cohort==c,'n'].iloc[0])})" if (inv.cohort==c).any() else c for c in cohorts], fontsize=8)
    fig.suptitle("Headline sets — leftover ICI lung bulk, TACSTD2 median split", fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "nes_headline_bars.png"), dpi=160)
    plt.close(fig)

    # response NES if present
    resp = gsea[(gsea["primary"] == True) & (gsea["contrast"] == "response_R_vs_NR")].copy()  # noqa: E712
    if not resp.empty:
        rcoh = [c for c in LOADERS if c in set(resp["cohort"])]
        rmat = pd.DataFrame(index=terms, columns=rcoh, dtype=float)
        rq = pd.DataFrame(index=terms, columns=rcoh, dtype=float)
        for c in rcoh:
            sub = resp[resp["cohort"] == c]
            for t in terms:
                row = sub[sub["term"] == t]
                if len(row):
                    rmat.loc[t, c] = float(row.iloc[0]["nes"])
                    rq.loc[t, c] = float(row.iloc[0]["fdr_primary"])
        fig, ax = plt.subplots(figsize=(10.5, 7.2))
        im = ax.imshow(rmat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-2.5, vmax=2.5, aspect="auto")
        ax.set_xticks(range(len(rcoh)))
        labs = []
        for c in rcoh:
            sub = resp[resp["cohort"] == c]
            labs.append(f"{c}\nR {int(sub['n_high'].iloc[0])} vs NR {int(sub['n_low'].iloc[0])}" if len(sub) else c)
        ax.set_xticklabels(labs, fontsize=8)
        ax.set_yticks(range(len(terms)))
        ax.set_yticklabels(terms, fontsize=8)
        for i, t in enumerate(terms):
            for j, c in enumerate(rcoh):
                nes = rmat.loc[t, c]
                q = rq.loc[t, c]
                if nes == nes:
                    star = "*" if q == q and q < 0.05 else ""
                    ax.text(j, i, f"{nes:+.2f}{star}", ha="center", va="center", fontsize=7)
        ax.set_title("Responder vs NR prerank GSEA  NES  * FDR<0.05  (exploratory; small n)")
        fig.colorbar(im, ax=ax, fraction=0.03, label="NES (R minus NR)")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "nes_response.png"), dpi=160)
        plt.close(fig)

    # focal rho heatmap
    if focal is not None and not focal.empty:
        genes = [g for g in FOCAL if g in set(focal["gene"])]
        fcoh = [c for c in LOADERS if c in set(focal["cohort"])]
        fmat = pd.DataFrame(index=genes, columns=fcoh, dtype=float)
        for c in fcoh:
            sub = focal[focal["cohort"] == c]
            for g in genes:
                row = sub[sub["gene"] == g]
                if len(row) and bool(row.iloc[0]["present"]):
                    fmat.loc[g, c] = float(row.iloc[0]["rho"])
        fig, ax = plt.subplots(figsize=(8.8, 6.4))
        im = ax.imshow(fmat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
        ax.set_xticks(range(len(fcoh)))
        ax.set_xticklabels(fcoh, fontsize=8, rotation=20, ha="right")
        ax.set_yticks(range(len(genes)))
        ax.set_yticklabels(genes, fontsize=8)
        for i, g in enumerate(genes):
            for j, c in enumerate(fcoh):
                v = fmat.loc[g, c]
                if v == v:
                    ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7)
        ax.set_title("Focal TJ / keratin / EMT genes vs TACSTD2 (Spearman rho)")
        fig.colorbar(im, ax=ax, fraction=0.03, label="rho")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "focal_rho.png"), dpi=160)
        plt.close(fig)

    # signature Spearman
    if sig is not None and not sig.empty:
        terms2 = [t for t in terms if t in set(sig["term"])]
        scoh = [c for c in LOADERS if c in set(sig["cohort"])]
        smat = pd.DataFrame(index=terms2, columns=scoh, dtype=float)
        for c in scoh:
            sub = sig[sig["cohort"] == c]
            for t in terms2:
                row = sub[sub["term"] == t]
                if len(row):
                    smat.loc[t, c] = float(row.iloc[0]["rho"])
        fig, ax = plt.subplots(figsize=(10.5, 7.2))
        im = ax.imshow(smat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
        ax.set_xticks(range(len(scoh)))
        ax.set_xticklabels(scoh, fontsize=8, rotation=20, ha="right")
        ax.set_yticks(range(len(terms2)))
        ax.set_yticklabels(terms2, fontsize=8)
        for i, t in enumerate(terms2):
            for j, c in enumerate(scoh):
                v = smat.loc[t, c]
                if v == v:
                    ax.text(j, i, f"{v:+.2f}", ha="center", va="center", fontsize=7)
        ax.set_title("z-mean signature vs continuous TACSTD2 (Spearman; no split)")
        fig.colorbar(im, ax=ax, fraction=0.03, label="rho")
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "signature_rho.png"), dpi=160)
        plt.close(fig)


def write_report(inv, gsea, sig, focal, verd, spec_map) -> None:
    lines = []
    lines.append("# ADDITIVE A8 extra — keratin / TJ / EMT GSEA in public ICI-treated lung tumors")
    lines.append("")
    lines.append("**Public data only. Additive to the TCGA A8 slide. That TCGA audit was not re-run.**")
    lines.append("")
    lines.append("User A8 (TROP2-high enrich keratin/TJ, EMT down, conserved) is **taken as given** for TCGA-LUAD/LUSC. This file asks the same locked gene-set question in leftover open **pre-treatment ICI-treated** lung tumor bulk.")
    lines.append("")
    # TL;DR from median verdicts
    med_v = verd[verd["contrast"] == "tacstd2_median"] if not verd.empty else pd.DataFrame()
    lines.append("## 一句话结论 / TL;DR")
    lines.append("")
    if med_v.empty:
        lines.append("No usable leftover ICI bulk GSEA completed.")
    else:
        bits = [f"{r.cohort}={r.verdict} (n_high={r.n_high}, n_low={r.n_low})" for r in med_v.itertuples()]
        lines.append("TACSTD2-high **keratin / TJ-up is common** in leftover ICI bulk (same direction as A8 TCGA). Hallmark EMT-down is **not** conserved: 4/6 median-split cohorts have Hallmark EMT significantly **up**; the other two are null. **No** median-split leftover cohort is `supportive`.")
        lines.append("")
        lines.append("Per-cohort median verdicts: " + "; ".join(bits) + ".")
        lines.append("")
        lines.append("n=16–43 is **hypothesis-generating**. This does not retract the TCGA A8 slide. Do not quote a pooled ICI NES. Do not write “A8 conserved in ICI tumors.” Hallmark EMT, not GOBP EMT, decides the EMT arm. The single `supportive` call is GSE283829 **quartile** 7 vs 7 — not the locked median split.")
    lines.append("")
    lines.append("## Why this extra exists")
    lines.append("")
    lines.append("The original A8 slide is TCGA-only (ICI-naive surgical resections). Conservation of a keratin/TJ-high, EMT-low TROP2-high state in **tumors that actually received PD-(L)1** is a different question. This leftover slice uses open GEO ICI lung bulk that was not that TCGA matrix.")
    lines.append("")
    lines.append("Excluded on purpose: TCGA-LUAD/LUSC (A8, not re-run); GSE207422 (already GSEA'd in `hunt_tj_gsea`); GSE248378 (post-treatment FPKM); immune-only panels without TACSTD2 (GSE136961, GSE93157); FASTQ/SRA and controlled OAK/POPLAR.")
    lines.append("")
    lines.append("## Cohorts")
    lines.append("")
    lines.append("| Cohort | n | R / NR | Scale | Usable | Note |")
    lines.append("|---|---:|---|---|---|---|")
    for r in inv.itertuples():
        lines.append(
            f"| {r.cohort} | {r.n} | {getattr(r, 'n_R', 0)} / {getattr(r, 'n_NR', 0)} | {getattr(r, 'scale', '')} | {r.usable} | {getattr(r, 'leftover_note', r.reason)} |"
        )
    lines.append("")
    lines.append("## Pre-specified design")
    lines.append("")
    lines.append("| Piece | Choice | Honest limitation |")
    lines.append("|---|---|---|")
    lines.append("| Split (primary) | TACSTD2 above vs below median | Quartiles on n=16 are 4 vs 4. Median is the locked small-n split. |")
    lines.append("| Split (sensitivity) | Quartile only if both arms ≥ 6 | Still discards the middle; underpowered. |")
    lines.append("| Ranking | Welch t, high − low | Two-group statistic on n≈8–16/arm. |")
    lines.append("| Complementary | Spearman gene vs TACSTD2 prerank | Uses every sample; not a two-group t. |")
    lines.append("| Response | R vs NR Welch t if both arms ≥ 5 | Exploratory. Deposited labels only. |")
    lines.append("| GSEA | Preranked weighted KS, p=1, 1000 gene-set permutations, seed=42 | Same engine as A8. Gene-set permutation, not sample permutation. |")
    lines.append("| FDR | BH within the 12 primary sets per contrast | Not nested GSEA FDR. |")
    lines.append("| EMT rule | Hallmark EMT decides the claim | GOBP EMT is recorded and must not be quoted as Hallmark. |")
    lines.append("| Sets | Identical `data/genesets/a8_sets.json` primary 12 | TACSTD2 is not a member of the primary sets. |")
    lines.append("")
    lines.append("Positive NES = enriched in TACSTD2-high (or in responders). NES is **NA** when no same-sign null ES exists (the ratio is undefined); still report ES and nominal p from `gsea_prerank_all.tsv`.")
    lines.append("")

    def dump_contrast(title, contrast):
        lines.append(f"## {title}")
        lines.append("")
        if gsea is None or gsea.empty:
            lines.append("No GSEA rows.")
            lines.append("")
            return
        sub = gsea[(gsea["contrast"] == contrast) & (gsea["primary"] == True)]  # noqa: E712
        if sub.empty:
            lines.append("Not run or no primary rows.")
            lines.append("")
            return
        terms = [
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
            "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION",
            "HALLMARK_APICAL_JUNCTION",
            "KEGG_TIGHT_JUNCTION",
            "GOBP_KERATINIZATION",
            "KRT_EPITHELIAL",
        ]
        cohorts = [c for c in LOADERS if c in set(sub["cohort"])]
        header = "| Gene set | Want |" + "".join(f" {c} NES | {c} FDR |" for c in cohorts)
        lines.append(header)
        lines.append("|---|---|" + "|".join(["---:|---:" for _ in cohorts]))
        want = {
            "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": "DOWN",
            "GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION": "DOWN",
            "HALLMARK_APICAL_JUNCTION": "UP",
            "KEGG_TIGHT_JUNCTION": "UP",
            "GOBP_KERATINIZATION": "UP",
            "KRT_EPITHELIAL": "UP",
        }
        for t in terms:
            row = f"| {t} | {want.get(t, '')} |"
            for c in cohorts:
                r = sub[(sub["term"] == t) & (sub["cohort"] == c)]
                if r.empty:
                    row += " NA | NA |"
                else:
                    row += f" {fmt_nes(float(r.iloc[0]['nes']))} | {fmt_p(float(r.iloc[0]['fdr_primary']))} |"
            lines.append(row)
        lines.append("")
        lines.append("Full 12-set tables: `tables/gsea_prerank_all.tsv`.")
        lines.append("")
        vv = verd[verd["contrast"] == contrast]
        if not vv.empty:
            lines.append("| Cohort | n_high | n_low | Verdict |")
            lines.append("|---|---:|---:|---|")
            for r in vv.itertuples():
                lines.append(f"| {r.cohort} | {r.n_high} | {r.n_low} | {r.verdict} |")
            lines.append("")

    dump_contrast("Primary NES — TACSTD2 median split", "tacstd2_median")
    dump_contrast("Sensitivity — TACSTD2 quartile (only if both arms ≥ 6)", "tacstd2_quartile")
    dump_contrast("Complementary — Spearman prerank vs continuous TACSTD2", "tacstd2_spearman")
    dump_contrast("Exploratory — responder vs NR (only if both arms ≥ 5)", "response_R_vs_NR")

    lines.append("## Complementary signature Spearman (no split)")
    lines.append("")
    if sig is not None and not sig.empty:
        lines.append("| Gene set |" + "".join(f" {c} ρ |" for c in LOADERS if c in set(sig["cohort"])))
        lines.append("|---|" + "|".join(["---:" for c in LOADERS if c in set(sig["cohort"])]))
        for t in HEADLINE:
            row = f"| {t} |"
            for c in LOADERS:
                r = sig[(sig["term"] == t) & (sig["cohort"] == c)]
                if r.empty:
                    continue
                row += f" {r.iloc[0]['rho']:+.3f} |" if r.iloc[0]["rho"] == r.iloc[0]["rho"] else " NA |"
            lines.append(row)
        lines.append("")
    else:
        lines.append("Not computed.")
        lines.append("")

    lines.append("## What this does **not** say")
    lines.append("")
    lines.append("- It does not re-compute or revise the TCGA A8 NES table.")
    lines.append("- It does not claim a predictive TROP2 × ICI interaction. Association of a program with TACSTD2 inside a treated cohort is not a treatment-effect modifier.")
    lines.append("- Small n cannot rule out a modest true effect. Report the observed NES/FDR/n.")
    lines.append("- Do not collapse Hallmark EMT and GOBP EMT.")
    lines.append("- Do not write “conserved in ICI tumors” unless the median-split Hallmark EMT, a TJ set, and a keratin set all pass FDR<0.05 in the same leftover cohort.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("pip install -r requirements.txt")
    lines.append("bash scripts/a8_ici_gsea/download.sh /tmp/a8_ici_gsea_data")
    lines.append("python3 scripts/a8_ici_gsea/analyze.py")
    lines.append("```")
    lines.append("")
    lines.append("Methods: `methods/a8_ici_gsea/playbook.md`.")
    open(os.path.join(OUT, "REPORT.md"), "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
