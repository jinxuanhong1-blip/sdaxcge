#!/usr/bin/env python3
"""CLDN4-only combo of GSE253564 (bulk FPKM) + GSE148071 (scRNA, given ρ).

This is **not CellChat**. GSE253564 public processed data are bulk FPKM
(no scRNA / GeoMx / fraction table). Combo is a patient-level Fisher-z
forest of the two *given* CLDN4 vs T/NK Spearman rhos, plus a bulk
chemokine / ligand–receptor product proxy on the FPKM matrix.

Given rhos are not re-audited.
"""
from __future__ import annotations

import json
import math
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from lib_stats import (
    bh_fdr,
    fisher_combine,
    fisher_z,
    fisher_z_var,
    implied_spearman_p,
    random_effects_dl,
    spearman,
    stouffer,
)

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
OUT = HERE / "results"
FPKM_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE253nnn/GSE253564/suppl/"
    "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"
)
FPKM_NAME = "GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz"

# Given. Do not re-estimate.
GIVEN = [
    {
        "cohort": "GSE253564",
        "assay": "bulk RNA-seq FPKM (pre-treatment)",
        "pmid": "38401548",
        "title": "Altorki et al. Nat Commun 2024; neoadjuvant durvalumab ± SBRT",
        "n": 32,
        "rho": -0.48,
        "n_note": "n=32 public pre-treatment FPKM columns",
        "rho_note": "given CLDN4 vs T/NK; not re-audited",
        "cellchat_eligible": False,
    },
    {
        "cohort": "GSE148071",
        "assay": "scRNA-seq (Wu et al.; 42 biopsies, partial usable)",
        "pmid": "33953163",
        "title": "Wu et al. Nat Commun 2021; advanced NSCLC scRNA",
        "n": 25,
        "rho": -0.49,
        "n_note": "n=25 partial of 42 (given)",
        "rho_note": "given CLDN4 vs T/NK; not re-audited",
        "cellchat_eligible": False,  # cannot combo-CellChat with bulk FPKM
    },
]

# Alias: GEO FPKM uses PVRL2, not NECTIN2. CXCL12 absent from this matrix.
ALIASES = {"NECTIN2": "PVRL2"}

T_RECRUIT_LIGANDS = [
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "CCL5",
    "CCL3",
    "CCL4",
    "CXCL16",
    "CCL19",
    "CCL21",
    "CXCL13",
    "CX3CL1",
    "CXCL12",
]
T_RECRUIT_RECEPTORS = [
    "CXCR3",
    "CCR5",
    "CCR1",
    "CXCR6",
    "CCR7",
    "CXCR5",
    "CX3CR1",
    "CXCR4",
]
IFN_GENES = ["IFNG", "IFNB1", "IFNGR1", "IFNGR2", "STAT1"]
CHECKPOINT_GENES = [
    "CD274",
    "PDCD1",
    "CD276",
    "PVR",
    "NECTIN2",
    "TIGIT",
    "LGALS9",
    "HAVCR2",
    "CEACAM1",
    "HLA-E",
    "HLA-A",
    "HLA-B",
    "HLA-C",
]
# Document-only T/NK markers (not used to replace the given ρ).
TNK_MARKERS = ["CD3D", "CD3E", "CD2", "CD8A", "CD8B", "NKG7", "GNLY", "KLRD1"]

# Curated Mal↔T/NK pairs used as bulk product proxies (CellChat-style labels only).
CURATED_PAIRS = [
    ("CXCL9", "CXCR3", "T_recruit"),
    ("CXCL10", "CXCR3", "T_recruit"),
    ("CXCL11", "CXCR3", "T_recruit"),
    ("CCL5", "CCR5", "T_recruit"),
    ("CCL3", "CCR5", "T_recruit"),
    ("CCL4", "CCR5", "T_recruit"),
    ("CCL5", "CCR1", "T_recruit"),
    ("CXCL16", "CXCR6", "T_recruit"),
    ("CCL19", "CCR7", "T_recruit"),
    ("CCL21", "CCR7", "T_recruit"),
    ("CXCL13", "CXCR5", "T_recruit"),
    ("CX3CL1", "CX3CR1", "T_recruit"),
    ("CXCL12", "CXCR4", "T_recruit"),
    ("IFNG", "IFNGR1", "IFN"),
    ("IFNG", "IFNGR2", "IFN"),
    ("CD274", "PDCD1", "checkpoint"),
    ("PVR", "TIGIT", "checkpoint"),
    ("NECTIN2", "TIGIT", "checkpoint"),
    ("PVR", "CD226", "checkpoint"),
    ("NECTIN2", "CD226", "checkpoint"),
    ("LGALS9", "HAVCR2", "checkpoint"),
    ("HLA-E", "KLRD1", "MHC_I"),
    ("HLA-A", "CD8A", "MHC_I"),
    ("HLA-B", "CD8A", "MHC_I"),
    ("HLA-C", "CD8A", "MHC_I"),
]


def _ci_rho(rho: float, n: int) -> tuple[float, float]:
    z = fisher_z(rho)
    se = math.sqrt(fisher_z_var(n))
    return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))


def ensure_fpkm() -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    dest = DATA / FPKM_NAME
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    cached = Path("/tmp/geo_combo") / FPKM_NAME
    if cached.exists() and cached.stat().st_size > 0:
        dest.write_bytes(cached.read_bytes())
        return dest
    print("downloading", FPKM_URL)
    urllib.request.urlretrieve(FPKM_URL, dest)
    return dest


def load_log2_fpkm(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t")
    gene = raw["gene"].astype(str)
    mat = raw.drop(columns=["gene", "Entrez.ID"]).apply(pd.to_numeric, errors="coerce")
    mat.index = gene
    # keep first occurrence if duplicated symbols
    mat = mat[~mat.index.duplicated(keep="first")]
    return np.log2(mat + 1.0)


def resolve(mat: pd.DataFrame, symbol: str) -> str | None:
    if symbol in mat.index:
        return symbol
    alt = ALIASES.get(symbol)
    if alt and alt in mat.index:
        return alt
    return None


def zmean(mat: pd.DataFrame, symbols: list[str]) -> tuple[pd.Series, list[str]]:
    used = []
    cols = []
    for s in symbols:
        key = resolve(mat, s)
        if key is None:
            continue
        x = mat.loc[key].astype(float)
        sd = float(x.std(ddof=1))
        if not math.isfinite(sd) or sd == 0:
            continue
        cols.append((x - float(x.mean())) / sd)
        used.append(s if key == s else f"{s}->{key}")
    if not cols:
        return pd.Series(np.nan, index=mat.columns), used
    return pd.concat(cols, axis=1).mean(axis=1), used


def forest(rows: pd.DataFrame, re: dict, path: Path) -> None:
    labels, rhos, los, his, ps = [], [], [], [], []
    for _, r in rows.iterrows():
        lo, hi = _ci_rho(float(r["rho"]), int(r["n"]))
        labels.append(f"{r['cohort']}  n={int(r['n'])}")
        rhos.append(float(r["rho"]))
        los.append(lo)
        his.append(hi)
        ps.append(float(r["p_implied"]))
    labels.append(f"RE pooled  N={re['n_patients_total']}")
    rhos.append(re["pooled_rho"])
    los.append(re["ci95_rho"][0])
    his.append(re["ci95_rho"][1])
    ps.append(re["p"])

    fig, ax = plt.subplots(figsize=(8.8, 2.4))
    y = np.arange(len(labels))
    for i, (rho, lo, hi) in enumerate(zip(rhos, los, his)):
        color = "#1f4e79" if i == len(labels) - 1 else "#4a4a4a"
        marker = "D" if i == len(labels) - 1 else "o"
        ax.plot([lo, hi], [y[i], y[i]], color=color, lw=1.6, solid_capstyle="round")
        ax.plot(rho, y[i], marker=marker, color=color, ms=7 if i == len(labels) - 1 else 6)
        ax.text(
            1.02,
            y[i],
            f"ρ={rho:+.2f}  p={ps[i]:.3g}",
            va="center",
            ha="left",
            fontsize=8,
            family="monospace",
        )
    ax.axvline(0, color="#888888", lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Spearman ρ  (CLDN4 vs T/NK; given cohort rhos, not re-audited)")
    ax.set_xlim(-1.05, 1.05)
    ax.set_title(
        f"CLDN4 vs T/NK  ·  RE ρ={re['pooled_rho']:+.3f} "
        f"[{re['ci95_rho'][0]:+.3f}, {re['ci95_rho'][1]:+.3f}]  "
        f"p={re['p']:.3g}  I²={re['I2']:.0f}%  N={re['n_patients_total']}",
        fontsize=10,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    rows = []
    for g in GIVEN:
        p = implied_spearman_p(g["rho"], g["n"])
        lo, hi = _ci_rho(g["rho"], g["n"])
        rows.append(
            {
                **g,
                "gene": "CLDN4",
                "endpoint": "T/NK",
                "p_implied": p,
                "ci95_lo": lo,
                "ci95_hi": hi,
                "p_source": "Spearman t from given ρ and n (not a matrix re-audit)",
            }
        )
    cohort_df = pd.DataFrame(rows)
    rhos = cohort_df["rho"].tolist()
    ns = cohort_df["n"].astype(int).tolist()
    ps = cohort_df["p_implied"].tolist()
    re = random_effects_dl(rhos, ns)
    st = stouffer(rhos, ps, ns)
    fi = fisher_combine(ps)

    combo_row = {
        "cohort": "combo_RE_GSE253564_GSE148071",
        "assay": "patient-level Fisher-z DerSimonian–Laird (not CellChat)",
        "pmid": "38401548+33953163",
        "title": "CLDN4-only combo of given rhos",
        "n": re["n_patients_total"],
        "rho": re["pooled_rho"],
        "n_note": "N = 32 + 25 (independent patients)",
        "rho_note": "RE pooled; given rhos not re-audited",
        "cellchat_eligible": False,
        "gene": "CLDN4",
        "endpoint": "T/NK",
        "p_implied": re["p"],
        "ci95_lo": re["ci95_rho"][0],
        "ci95_hi": re["ci95_rho"][1],
        "p_source": re["method"],
        "I2": re["I2"],
        "tau2": re["tau2"],
        "Q": re["Q"],
        "stouffer_z": st["z"],
        "stouffer_p": st["p"],
        "fisher_p": fi["p"],
    }
    combo_df = pd.concat([cohort_df, pd.DataFrame([combo_row])], ignore_index=True)
    combo_df.to_csv(OUT / "combo_cldn4_tnk.tsv", sep="\t", index=False)
    forest(cohort_df, re, OUT / "forest_CLDN4_tnk.png")

    # --- bulk chemokine / LR-proxy (GSE253564 only) ---
    mat = load_log2_fpkm(ensure_fpkm())
    n_bulk = int(mat.shape[1])
    if n_bulk != 32:
        raise SystemExit(f"expected 32 FPKM columns, got {n_bulk}")
    cldn4 = mat.loc["CLDN4"].astype(float)

    gene_rows = []
    panels = {
        "T_recruit_ligand": T_RECRUIT_LIGANDS,
        "T_recruit_receptor": T_RECRUIT_RECEPTORS,
        "IFN": IFN_GENES,
        "checkpoint": CHECKPOINT_GENES,
        "TNK_marker_doc_only": TNK_MARKERS,
    }
    for panel, genes in panels.items():
        for g in genes:
            key = resolve(mat, g)
            if key is None:
                gene_rows.append(
                    {
                        "panel": panel,
                        "gene": g,
                        "matrix_id": "",
                        "n": n_bulk,
                        "rho_vs_CLDN4": np.nan,
                        "p": np.nan,
                        "status": "absent_from_FPKM",
                    }
                )
                continue
            rho, p, n = spearman(cldn4, mat.loc[key])
            gene_rows.append(
                {
                    "panel": panel,
                    "gene": g,
                    "matrix_id": key,
                    "n": n,
                    "rho_vs_CLDN4": rho,
                    "p": p,
                    "status": "ok",
                }
            )
    gene_df = pd.DataFrame(gene_rows)
    # FDR within non-doc panels only (the hypothesis-facing chemokine/LR genes)
    face = gene_df["panel"].isin(
        ["T_recruit_ligand", "T_recruit_receptor", "IFN", "checkpoint"]
    ) & (gene_df["status"] == "ok")
    q = np.full(len(gene_df), np.nan)
    q[face.to_numpy()] = bh_fdr(gene_df.loc[face, "p"].tolist())
    gene_df["q_BH"] = q
    gene_df.sort_values(["panel", "rho_vs_CLDN4"], inplace=True)
    gene_df.to_csv(OUT / "gse253564_cldn4_vs_genes.tsv", sep="\t", index=False)

    pair_rows = []
    for lig, rec, pathway in CURATED_PAIRS:
        lk, rk = resolve(mat, lig), resolve(mat, rec)
        if lk is None or rk is None:
            pair_rows.append(
                {
                    "ligand": lig,
                    "receptor": rec,
                    "pathway": pathway,
                    "n": n_bulk,
                    "rho_vs_CLDN4": np.nan,
                    "p": np.nan,
                    "status": "missing_subunit",
                    "ligand_id": lk or "",
                    "receptor_id": rk or "",
                }
            )
            continue
        # geometric-mean proxy of (FPKM+1): mean of log2(FPKM+1)
        score = 0.5 * (mat.loc[lk].astype(float) + mat.loc[rk].astype(float))
        rho, p, n = spearman(cldn4, score)
        pair_rows.append(
            {
                "ligand": lig,
                "receptor": rec,
                "pathway": pathway,
                "n": n,
                "rho_vs_CLDN4": rho,
                "p": p,
                "status": "ok",
                "ligand_id": lk,
                "receptor_id": rk,
            }
        )
    pair_df = pd.DataFrame(pair_rows)
    okp = pair_df["status"] == "ok"
    q = np.full(len(pair_df), np.nan)
    q[okp.to_numpy()] = bh_fdr(pair_df.loc[okp, "p"].tolist())
    pair_df["q_BH"] = q
    pair_df.sort_values(["pathway", "rho_vs_CLDN4"], inplace=True)
    pair_df.to_csv(OUT / "gse253564_cldn4_lr_proxy.tsv", sep="\t", index=False)

    score_rows = []
    for name, genes in {
        "T_recruit_ligand": T_RECRUIT_LIGANDS,
        "T_recruit_receptor": T_RECRUIT_RECEPTORS,
        "IFN": IFN_GENES,
        "checkpoint": CHECKPOINT_GENES,
    }.items():
        s, used = zmean(mat, genes)
        rho, p, n = spearman(cldn4, s)
        score_rows.append(
            {
                "score": name,
                "n": n,
                "n_genes_used": len(used),
                "genes_used": ",".join(used),
                "rho_vs_CLDN4": rho,
                "p": p,
            }
        )
    score_df = pd.DataFrame(score_rows)
    score_df["q_BH"] = bh_fdr(score_df["p"].tolist())
    score_df.to_csv(OUT / "gse253564_cldn4_vs_scores.tsv", sep="\t", index=False)

    summary = {
        "not_cellchat": True,
        "reason": (
            "GSE253564 public processed file is bulk FPKM "
            "(GSE253564_Pre-treatment_Samples_Pubs_FPKMs.txt.gz), "
            "not scRNA / GeoMx / bulk-with-fractions. "
            "Cannot merge with GSE148071 scRNA for CellChat."
        ),
        "given_rhos_not_re_audited": True,
        "dual_high": False,
        "gene": "CLDN4",
        "combo": {
            "k": re["k"],
            "N": re["n_patients_total"],
            "pooled_rho": re["pooled_rho"],
            "ci95": re["ci95_rho"],
            "p_RE": re["p"],
            "I2": re["I2"],
            "tau2": re["tau2"],
            "Q": re["Q"],
            "fixed_rho": re["fixed_rho"],
            "stouffer_z": st["z"],
            "stouffer_p": st["p"],
            "fisher_p": fi["p"],
            "method": re["method"],
        },
        "cohorts": combo_df.to_dict(orient="records"),
        "gse253564": {
            "n_genes": int(mat.shape[0]),
            "n_samples": n_bulk,
            "file": FPKM_NAME,
            "transform": "log2(FPKM+1)",
            "lr_proxy": "Spearman(CLDN4, mean(log2(L+1), log2(R+1)))",
        },
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print(json.dumps(summary["combo"], indent=2))
    print("scores:")
    print(score_df.to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
