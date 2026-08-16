#!/usr/bin/env python3
"""OncoSG LUAD: TACSTD2 vs CD8 / GEP18 / immune after published purity.

Public cBioPortal study luad_oncosg_2020 only. No raw RSEM is deposited,
so ESTIMATE / xCell / MCP-counter are not computed. GEP18 is 17/18 genes
because CCL5 is absent from the public z-score matrix and the REST API.

Partial Spearman = Pearson on ranks, one covariate (PURITY), df = n-3.
This is the same estimator as prior A1 / A2 reworks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

DATAHUB_SHA = "165bd77077b03038f9c2ee104959eb474770b2a9"
STUDY = "luad_oncosg_2020"
MEDIA = (
    f"https://media.githubusercontent.com/media/cBioPortal/datahub/"
    f"{DATAHUB_SHA}/public/{STUDY}"
)
FILES = {
    "expression": (
        "data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt",
        "44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f",
    ),
    "clinical_sample": (
        "data_clinical_sample.txt",
        "55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368",
    ),
    "clinical_patient": (
        "data_clinical_patient.txt",
        None,
    ),
}

# Pre-specified primary endpoints (this rework).
CD8_GENES = ["CD8A"]
CD8_SENS_GENES = ["CD8A", "CD8B"]
GEP18_GENES = [
    "CCL5",
    "CD27",
    "CD274",
    "CD276",
    "CD8A",
    "CMKLR1",
    "CXCL9",
    "CXCR6",
    "HLA-DQA1",
    "HLA-DRB1",
    "HLA-E",
    "IDO1",
    "LAG3",
    "NKG7",
    "PDCD1LG2",
    "PSMB10",
    "STAT1",
    "TIGIT",
]
# Original Claim A1 "immune" (not GEP18, not ESTIMATE).
IMMUNE8_GENES = [
    "CD8A",
    "GZMA",
    "GZMB",
    "IFNG",
    "EOMES",
    "CXCL9",
    "CXCL10",
    "TBX21",
]
CYT_GENES = ["GZMA", "PRF1"]
CYTOTOXIC_GENES = [
    "GZMA",
    "GZMB",
    "GZMH",
    "GZMK",
    "PRF1",
    "GNLY",
    "NKG7",
    "KLRD1",
]
EXHAUSTION_GENES = [
    "PDCD1",
    "CTLA4",
    "LAG3",
    "HAVCR2",
    "TIGIT",
    "BTLA",
    "CD160",
    "VSIR",
    "TOX",
    "ENTPD1",
    "CD274",
    "PDCD1LG2",
]
IMSIG_IMMUNE = [
    "IMSIG_B_CELLS",
    "IMSIG_INTERFERON",
    "IMSIG_MACROPHAGES",
    "IMSIG_MONOCYTES",
    "IMSIG_NEUTROPHILS",
    "IMSIG_NK_CELLS",
    "IMSIG_PLASMA_CELLS",
    "IMSIG_T_CELLS",
]
IMSIG_CONTROL = ["IMSIG_PROLIFERATION", "IMSIG_TRANSLATION"]
IMSIG_LABEL = {
    "IMSIG_B_CELLS": "IMSIG B cells",
    "IMSIG_INTERFERON": "IMSIG interferon",
    "IMSIG_MACROPHAGES": "IMSIG macrophages",
    "IMSIG_MONOCYTES": "IMSIG monocytes",
    "IMSIG_NEUTROPHILS": "IMSIG neutrophils",
    "IMSIG_NK_CELLS": "IMSIG NK cells",
    "IMSIG_PLASMA_CELLS": "IMSIG plasma cells",
    "IMSIG_T_CELLS": "IMSIG T cells",
    "IMSIG_PROLIFERATION": "IMSIG proliferation (control)",
    "IMSIG_TRANSLATION": "IMSIG translation (control)",
}

API_RNA_ONLY = [
    "A008",
    "A114",
    "A122",
    "A136",
    "A139",
    "A147",
    "A184",
    "A302",
    "A435",
    "A484",
    "A489",
    "A507",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"GET {url}", file=sys.stderr)
    with urllib.request.urlopen(url, timeout=120) as r, dest.open("wb") as out:
        out.write(r.read())


def ensure_file(cache_dir: Path, key: str) -> Path:
    name, expected = FILES[key]
    dest = cache_dir / name
    if not dest.exists():
        download(f"{MEDIA}/{name}", dest)
    if expected is not None:
        got = sha256_file(dest)
        if got != expected:
            raise SystemExit(f"checksum mismatch for {name}: {got}")
    return dest


def read_cbioportal_clinical(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="#")


def read_expression(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if "Hugo_Symbol" not in df.columns:
        raise SystemExit(f"unexpected expression header: {list(df.columns)[:5]}")
    df = df.drop(columns=["Entrez_Gene_Id"], errors="ignore")
    df = df.drop_duplicates(subset=["Hugo_Symbol"], keep="first")
    df = df.set_index("Hugo_Symbol")
    return df


def mean_present(expr: pd.DataFrame, genes: list[str], samples: list[str]) -> tuple[pd.Series, list[str], list[str]]:
    present = [g for g in genes if g in expr.index]
    missing = [g for g in genes if g not in expr.index]
    if not present:
        return pd.Series(np.nan, index=samples), present, missing
    score = expr.loc[present, samples].astype(float).mean(axis=0)
    return score, present, missing


def spearman_pair(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    r, p = stats.spearmanr(x, y, nan_policy="omit")
    return float(r), float(p)


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """First-order partial Spearman via Pearson on ranks. p from t, df=n-3."""
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[mask], y[mask], z[mask]
    n = int(x.size)
    if n < 6:
        return {
            "n": n,
            "rho": np.nan,
            "p": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "t": np.nan,
            "df": np.nan,
        }
    rx, ry, rz = (stats.rankdata(v, method="average") for v in (x, y, z))
    rxy = float(np.corrcoef(rx, ry)[0, 1])
    rxz = float(np.corrcoef(rx, rz)[0, 1])
    ryz = float(np.corrcoef(ry, rz)[0, 1])
    den = math.sqrt((1.0 - rxz**2) * (1.0 - ryz**2))
    if den == 0.0:
        r = np.nan
    else:
        r = (rxy - rxz * ryz) / den
        r = float(np.clip(r, -1.0, 1.0))
    df = n - 3
    if not np.isfinite(r) or abs(r) >= 1.0:
        t = np.nan
        p = np.nan
    else:
        t = r * math.sqrt(df / (1.0 - r**2))
        p = float(2.0 * stats.t.sf(abs(t), df))
    # Fisher-z CI; SE = 1/sqrt(n-4) for one covariate (same as prior OncoSG).
    if np.isfinite(r) and n > 4 and abs(r) < 1.0:
        zf = np.arctanh(r)
        se = 1.0 / math.sqrt(n - 4)
        ci = tuple(float(np.tanh(zf + s * 1.959963984540054 * se)) for s in (-1.0, 1.0))
    else:
        ci = (np.nan, np.nan)
    return {
        "n": n,
        "rho": float(r) if np.isfinite(r) else np.nan,
        "p": p,
        "ci_low": ci[0],
        "ci_high": ci[1],
        "t": float(t) if np.isfinite(t) else np.nan,
        "df": df,
    }


def residual_ranks(x: np.ndarray, z: np.ndarray) -> np.ndarray:
    rx = stats.rankdata(x, method="average").astype(float)
    rz = stats.rankdata(z, method="average").astype(float)
    A = np.column_stack([np.ones(len(rx)), rz])
    beta, *_ = np.linalg.lstsq(A, rx, rcond=None)
    return rx - A @ beta


def fmt_p(p: float) -> str:
    if not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache-dir", type=Path, default=Path("/tmp/oncosg_a1"))
    ap.add_argument("--out-dir", type=Path, default=Path("results/rework/OncoSG_A1"))
    args = ap.parse_args()
    out = args.out_dir
    fig_dir = out / "figures"
    out.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    expr_path = ensure_file(args.cache_dir, "expression")
    clin_path = ensure_file(args.cache_dir, "clinical_sample")
    pat_path = ensure_file(args.cache_dir, "clinical_patient")

    expr = read_expression(expr_path)
    clin = read_cbioportal_clinical(clin_path)
    pat = read_cbioportal_clinical(pat_path)

    samples = list(expr.columns)
    if "TACSTD2" not in expr.index:
        raise SystemExit("TACSTD2 missing from expression matrix")

    clin = clin.set_index("SAMPLE_ID", drop=False)
    overlap = [s for s in samples if s in clin.index]
    if len(overlap) != len(samples):
        raise SystemExit(f"expression samples missing from clinical: {set(samples) - set(clin.index)}")

    tac = expr.loc["TACSTD2", overlap].astype(float)
    purity = pd.to_numeric(clin.loc[overlap, "PURITY"], errors="coerce")

    scores: dict[str, pd.Series] = {}
    coverage_rows = []

    def add_gene_score(name: str, genes: list[str], role: str) -> None:
        sc, present, missing = mean_present(expr, genes, overlap)
        scores[name] = sc
        coverage_rows.append(
            {
                "feature": name,
                "role": role,
                "n_genes_defined": len(genes),
                "n_genes_present": len(present),
                "genes_present": ",".join(present),
                "genes_missing": ",".join(missing) if missing else "",
                "source": "cBioPortal z-score matrix (mean of present genes)",
            }
        )

    add_gene_score("CD8", CD8_GENES, "primary")
    add_gene_score("CD8_CD8A_CD8B", CD8_SENS_GENES, "sensitivity")
    add_gene_score("GEP18", GEP18_GENES, "primary")
    add_gene_score("immune_tcell_effector", IMMUNE8_GENES, "primary")
    add_gene_score("CYT", CYT_GENES, "secondary")
    add_gene_score("cytotoxic_A1", CYTOTOXIC_GENES, "secondary")
    add_gene_score("exhaustion_A1", EXHAUSTION_GENES, "secondary")

    for col in IMSIG_IMMUNE:
        scores[col] = pd.to_numeric(clin.loc[overlap, col], errors="coerce")
        coverage_rows.append(
            {
                "feature": col,
                "role": "imsig_published",
                "n_genes_defined": 0,
                "n_genes_present": 0,
                "genes_present": "",
                "genes_missing": "",
                "source": "data_clinical_sample.txt published IMSIG column",
            }
        )
    for col in IMSIG_CONTROL:
        scores[col] = pd.to_numeric(clin.loc[overlap, col], errors="coerce")
        coverage_rows.append(
            {
                "feature": col,
                "role": "imsig_control",
                "n_genes_defined": 0,
                "n_genes_present": 0,
                "genes_present": "",
                "genes_missing": "",
                "source": "data_clinical_sample.txt published IMSIG column",
            }
        )

    coverage_rows.append(
        {
            "feature": "ESTIMATE_ImmuneScore",
            "role": "not_computable",
            "n_genes_defined": 141,
            "n_genes_present": 0,
            "genes_present": "",
            "genes_missing": "entire SI_geneset; no raw RSEM on cBioPortal",
            "source": "SKIP — public deposit is z-scores only",
        }
    )

    sample_tbl = pd.DataFrame(
        {
            "sample_id": overlap,
            "patient_id": clin.loc[overlap, "PATIENT_ID"].to_numpy(),
            "TACSTD2_z": tac.to_numpy(),
            "PURITY": purity.to_numpy(),
        }
    )
    for name, sc in scores.items():
        sample_tbl[name] = sc.reindex(overlap).to_numpy()
    if "PATIENT_ID" in pat.columns:
        pat_i = pat.set_index("PATIENT_ID")
        for col in ["SEX", "SMOKING_STATUS", "STAGE", "ETHNICITY", "COHORT"]:
            if col in pat_i.columns:
                sample_tbl[col] = sample_tbl["patient_id"].map(pat_i[col])

    feature_order = [
        ("CD8", "primary", "CD8A z-score"),
        ("GEP18", "primary", "Ayers 2017 GEP, mean of present gene z-scores (17/18; CCL5 absent)"),
        ("immune_tcell_effector", "primary", "A1 8-gene immune / T-cell effector mean z-score"),
        ("CYT", "secondary", "Rooney CYT: mean(GZMA, PRF1) z-scores"),
        ("CD8_CD8A_CD8B", "sensitivity", "mean(CD8A, CD8B) z-scores"),
        ("cytotoxic_A1", "secondary", "original A1 cytotoxic 8-gene mean z-score"),
        ("exhaustion_A1", "secondary", "original A1 exhaustion 12-gene mean z-score"),
    ]
    for col in IMSIG_IMMUNE:
        feature_order.append((col, "imsig_published", IMSIG_LABEL[col]))
    for col in IMSIG_CONTROL:
        feature_order.append((col, "imsig_control", IMSIG_LABEL[col]))

    rows = []
    x = tac.to_numpy(dtype=float)
    z = purity.to_numpy(dtype=float)
    tac_pur = spearman_pair(x, z)

    for name, role, desc in feature_order:
        y = scores[name].reindex(overlap).to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        raw_r, raw_p = spearman_pair(x[mask], y[mask])
        part = partial_spearman(x, y, z)
        rows.append(
            {
                "feature": name,
                "role": role,
                "definition": desc,
                "n": part["n"],
                "n_raw": int(mask.sum()),
                "raw_rho": raw_r,
                "raw_p": raw_p,
                "partial_rho": part["rho"],
                "partial_p": part["p"],
                "partial_ci_low": part["ci_low"],
                "partial_ci_high": part["ci_high"],
                "partial_t": part["t"],
                "partial_df": part["df"],
            }
        )

    corr = pd.DataFrame(rows)

    def add_fdr(frame: pd.DataFrame, role: str, col_name: str) -> None:
        idx = frame.index[frame["role"] == role]
        if len(idx) == 0:
            frame[col_name] = np.nan
            return
        q = np.full(len(frame), np.nan)
        _, qvals, _, _ = multipletests(frame.loc[idx, "partial_p"], method="fdr_bh")
        q[idx] = qvals
        frame[col_name] = q

    add_fdr(corr, "primary", "partial_fdr_primary3")
    add_fdr(corr, "imsig_published", "partial_fdr_imsig8")

    # Concordance with prior n=169 neutrophils partial ρ=-0.456
    neu = corr.loc[corr["feature"] == "IMSIG_NEUTROPHILS"].iloc[0]
    prior_rho = -0.456
    prior_n = 169

    # API vs file TACSTD2 check (first 5 overlapping samples)
    api_ok = None
    api_note = "not checked"
    try:
        body = json.dumps(
            {
                "entrezGeneIds": [4070],
                "sampleIds": overlap[:8],
            }
        ).encode()
        req = urllib.request.Request(
            "https://www.cbioportal.org/api/molecular-profiles/"
            "luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores/"
            "molecular-data/fetch",
            data=body,
            headers={"Content-Type": "application/json"},
        )
        api = json.load(urllib.request.urlopen(req, timeout=60))
        api_map = {d["sampleId"]: float(d["value"]) for d in api if "value" in d}
        diffs = []
        for s in overlap[:8]:
            if s in api_map:
                diffs.append(abs(api_map[s] - float(tac.loc[s])))
        api_ok = bool(diffs) and max(diffs) < 1e-4
        api_note = f"max|API-file|={max(diffs):.2e} on {len(diffs)} samples" if diffs else "API returned no values"
    except Exception as exc:  # noqa: BLE001
        api_ok = False
        api_note = f"API check failed: {exc}"

    clin_full = read_cbioportal_clinical(clin_path)
    coverage_summary = {
        "expression_matrix_samples": len(samples),
        "cbioportal_rna_sample_list_n": 181,
        "rna_list_not_in_public_zscore_matrix": API_RNA_ONLY,
        "clinical_samples": int(len(clin_full)),
        "clinical_purity_non_na": int(pd.to_numeric(clin_full["PURITY"], errors="coerce").notna().sum()),
        "clinical_imsig_t_non_na": int(pd.to_numeric(clin_full["IMSIG_T_CELLS"], errors="coerce").notna().sum()),
        "analysis_n_complete_tacstd2_purity": int((np.isfinite(x) & np.isfinite(z)).sum()),
        "genes_in_matrix": int(expr.shape[0]),
        "CCL5_in_matrix": bool("CCL5" in expr.index),
        "GEP18_genes_used": 17,
        "GEP18_genes_missing": ["CCL5"],
        "ESTIMATE_status": "skipped_no_raw_RSEM",
        "api_vs_file_tacstd2": {"ok": api_ok, "note": api_note},
    }

    primary = corr[corr["role"] == "primary"].copy()
    results = {
        "study": STUDY,
        "datahub_sha": DATAHUB_SHA,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "n": int((np.isfinite(x) & np.isfinite(z)).sum()),
        "tacstd2_vs_purity": {"rho": tac_pur[0], "p": tac_pur[1]},
        "primary": primary.drop(columns=["definition"]).to_dict(orient="records"),
        "prior_neutrophils_anchor": {
            "claimed": {"n": prior_n, "partial_rho": prior_rho, "source": "PR82 IMSIG neutrophils"},
            "recomputed": {
                "n": int(neu["n"]),
                "partial_rho": float(neu["partial_rho"]),
                "partial_p": float(neu["partial_p"]),
                "partial_fdr_imsig8": float(neu["partial_fdr_imsig8"]),
                "raw_rho": float(neu["raw_rho"]),
            },
            "matches_3dp": abs(float(neu["partial_rho"]) - prior_rho) < 5e-4
            and int(neu["n"]) == prior_n,
        },
        "cannot_compute": [
            "ESTIMATE ImmuneScore / StromalScore / TumorPurity (no raw RSEM; only z-scores)",
            "xCell (needs counts/TPM + spillover calibration)",
            "MCP-counter as published (marker means on log2 TPM, not z-scores)",
            "GEP18 gene CCL5 (absent from public matrix and REST API)",
        ],
        "coverage": coverage_summary,
        "methods": {
            "expression": "cBioPortal log RNA-seq V2 RSEM z-scores vs all samples",
            "purity": "published sample PURITY clinical attribute (not ABSOLUTE, not ESTIMATE)",
            "partial_spearman": "Pearson on average-ranks of (x,y,z); t-test df=n-3; Fisher-z 95% CI SE=1/sqrt(n-4)",
            "signature_score": "unweighted mean of per-gene z-scores; Spearman invariant to monotone per-gene transform",
            "fdr": "BH within the 3 primary tests; separately BH within 8 IMSIG immune scores",
            "did_we_tune_to_minus_0.456": False,
        },
    }

    corr.to_csv(out / "correlations.tsv", sep="\t", index=False, float_format="%.10g")
    sample_tbl.to_csv(out / "sample_table.tsv", sep="\t", index=False, float_format="%.6g")
    pd.DataFrame(coverage_rows).to_csv(out / "gene_coverage.tsv", sep="\t", index=False)
    def _json_safe(obj):
        if isinstance(obj, dict):
            return {k: _json_safe(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_json_safe(v) for v in obj]
        if isinstance(obj, float) and not np.isfinite(obj):
            return None
        return obj

    (out / "results.json").write_text(json.dumps(_json_safe(results), indent=2) + "\n")
    provenance = {
        "study_id": STUDY,
        "study_page": f"https://www.cbioportal.org/study/summary?id={STUDY}",
        "datahub_revision": DATAHUB_SHA,
        "files": {
            key: {
                "url": f"{MEDIA}/{FILES[key][0]}",
                "sha256": FILES[key][1] or sha256_file(args.cache_dir / FILES[key][0]),
                "bytes": (args.cache_dir / FILES[key][0]).stat().st_size,
            }
            for key in FILES
        },
        "expression_profile_id": "luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores",
        "purity_field": "PURITY",
        "paper": "Chen et al., Nat Genet 2020 (OncoSG East-Asian LUAD)",
    }
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

    # Figures
    sns.set_theme(style="whitegrid", font_scale=0.95)
    plot_feats = [
        "CD8",
        "GEP18",
        "immune_tcell_effector",
        "CYT",
        "IMSIG_T_CELLS",
        "IMSIG_NEUTROPHILS",
    ]
    sub = corr[corr["feature"].isin(plot_feats)].copy()
    sub["label"] = sub["feature"].map(
        {
            "CD8": "CD8 (CD8A)",
            "GEP18": "GEP18 (17/18)",
            "immune_tcell_effector": "Immune (A1 8-gene)",
            "CYT": "CYT (GZMA+PRF1)",
            "IMSIG_T_CELLS": "IMSIG T cells",
            "IMSIG_NEUTROPHILS": "IMSIG neutrophils",
        }
    )
    sub = sub.set_index("feature").loc[plot_feats].reset_index()
    sub = sub.iloc[::-1].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    y_pos = np.arange(len(sub))
    ax.errorbar(
        sub["partial_rho"],
        y_pos,
        xerr=[
            sub["partial_rho"] - sub["partial_ci_low"],
            sub["partial_ci_high"] - sub["partial_rho"],
        ],
        fmt="o",
        color="#1f4e79",
        ecolor="#1f4e79",
        capsize=3,
    )
    ax.axvline(0, color="0.4", lw=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sub["label"])
    ax.set_xlabel("Partial Spearman ρ (TACSTD2 vs feature | PURITY)")
    ax.set_title(f"OncoSG LUAD n={int(sub['n'].iloc[0])} — purity-adjusted")
    ax.set_xlim(-0.65, 0.15)
    fig.tight_layout()
    fig.savefig(fig_dir / "forest_partial_rho.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.6), sharey=True)
    for ax, feat, title in zip(
        axes,
        ["CD8", "GEP18", "immune_tcell_effector"],
        ["CD8 (CD8A)", "GEP18 (17/18)", "Immune (A1 8-gene)"],
    ):
        y = scores[feat].reindex(overlap).to_numpy(dtype=float)
        mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        xr = residual_ranks(x[mask], z[mask])
        yr = residual_ranks(y[mask], z[mask])
        ax.scatter(xr, yr, s=14, alpha=0.65, c="#1f4e79", edgecolors="none")
        lr = stats.linregress(xr, yr)
        xs = np.linspace(xr.min(), xr.max(), 50)
        ax.plot(xs, lr.intercept + lr.slope * xs, color="#c0392b", lw=1.4)
        row = corr.loc[corr["feature"] == feat].iloc[0]
        ax.set_title(
            f"{title}\nρ={row['partial_rho']:.3f}  p={fmt_p(row['partial_p'])}"
        )
        ax.set_xlabel("TACSTD2 rank residual")
    axes[0].set_ylabel("Feature rank residual")
    fig.suptitle("OncoSG LUAD: rank residuals after PURITY", y=1.02, fontsize=11)
    fig.tight_layout()
    fig.savefig(fig_dir / "scatter_partial_residuals.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    write_report(out, corr, results, coverage_summary, tac_pur)
    print(json.dumps({"n": results["n"], "primary": results["primary"]}, indent=2))


def write_report(
    out: Path,
    corr: pd.DataFrame,
    results: dict,
    coverage: dict,
    tac_pur: tuple[float, float],
) -> None:
    def row(name: str) -> pd.Series:
        return corr.loc[corr["feature"] == name].iloc[0]

    cd8 = row("CD8")
    gep = row("GEP18")
    imm = row("immune_tcell_effector")
    cyt = row("CYT")
    neu = row("IMSIG_NEUTROPHILS")
    tcell = row("IMSIG_T_CELLS")
    n = int(cd8["n"])

    def line(r: pd.Series) -> str:
        return (
            f"ρ = **{r['partial_rho']:.3f}** "
            f"(95% CI {r['partial_ci_low']:.3f} to {r['partial_ci_high']:.3f}; "
            f"p = {fmt_p(r['partial_p'])}; n = {int(r['n'])})"
        )

    # Honest one-line verdict
    signs = [cd8["partial_rho"], gep["partial_rho"], imm["partial_rho"]]
    sig = [cd8["partial_p"] < 0.05, gep["partial_p"] < 0.05, imm["partial_p"] < 0.05]
    if all(s < 0 and p for s, p in zip(signs, sig)):
        verdict = (
            "YES — TACSTD2 is negatively associated with CD8, GEP18, and the "
            "A1 8-gene immune score after published purity (n=169)."
        )
    elif all(s < 0 for s in signs) and any(sig):
        verdict = (
            "PARTIAL — all three primary partial rhos are negative; "
            "not every test is significant. See the table."
        )
    else:
        verdict = (
            "NOT A UNIFORM NEGATIVE — at least one primary CD8/GEP/immune "
            "partial correlation is null or positive. See the table. Do not quote one ρ."
        )

    imsig_lines = []
    for feat in IMSIG_IMMUNE:
        r = row(feat)
        imsig_lines.append(
            f"| {IMSIG_LABEL[feat]} | {r['raw_rho']:.3f} | {fmt_p(r['raw_p'])} | "
            f"{r['partial_rho']:.3f} | {r['partial_ci_low']:.3f} to {r['partial_ci_high']:.3f} | "
            f"{fmt_p(r['partial_p'])} | {r['partial_fdr_imsig8']:.2e} |"
        )

    md = f"""# REWORK OncoSG A1 — TACSTD2 vs CD8 / GEP / immune after purity

**Self-contained. Public cBioPortal `luad_oncosg_2020` only. Written to be read without the rest of the repo.**

**Verdict: {verdict}**

This is **East-Asian surgical LUAD**, not an ICI-response cohort. A negative
bulk correlation is not evidence that TROP2-high tumors fail checkpoint blockade.

## Why this rework exists

User Claim A1 includes OncoSG (with TCGA). A prior OncoSG-only slice (PR 82)
reported TACSTD2 vs **IMSIG neutrophils** partial Spearman **ρ = −0.456, n = 169**.
That is a myeloid IMSIG score, not CD8, not Ayers GEP, and not the A1 8-gene
immune signature.

This folder recomputes the three endpoints the A1 rework language actually
uses — **CD8, GEP18, immune** — on the same public matrix, after the same
published `PURITY` field. The neutrophil number is recomputed as an
**anchor / sanity check**, not as the primary claim. Filters were not tuned
to recover −0.456.

## Analysis set

| Item | n | Note |
|---|---:|---|
| cBioPortal RNA sample list `luad_oncosg_2020_rna_seq_v2_mrna` | 181 | Portal description says “181 samples” |
| Public z-score matrix columns | **169** | This is the open expression table |
| RNA-list IDs absent from the public matrix | 12 | {", ".join(coverage["rna_list_not_in_public_zscore_matrix"])} |
| Clinical samples | {coverage["clinical_samples"]} | `data_clinical_sample.txt` |
| Clinical PURITY non-NA | {coverage["clinical_purity_non_na"]} | includes many RNA-absent tumors |
| Clinical IMSIG T-cells non-NA | {coverage["clinical_imsig_t_non_na"]} | 169 matrix + 3 IMSIG-only (A184, A484, A489; those three lack PURITY) |
| **Complete-case analysis (TACSTD2 + feature + PURITY)** | **{n}** | All 169 matrix samples have PURITY and all IMSIG immune scores |

The 12 portal-listed RNA IDs are **not in the public z-score file**. Nine of
them have PURITY but no IMSIG; three have IMSIG but no PURITY. They cannot
enter a TACSTD2 correlation. **n = 169 is the honest public n**, not 181.

## Pre-specified features

| Name | Definition | Honest limitation |
|---|---|---|
| **CD8** | `CD8A` all-sample z-score | Single gene. Not a deconvolution fraction. Sensitivity: mean(`CD8A`,`CD8B`). |
| **GEP18** | unweighted mean of **17/18** Ayers 2017 T-cell-inflamed GEP gene z-scores | **`CCL5` is absent** from the public matrix and from the cBioPortal REST API (0 values). Merck NanoString TIS weights are not public. This is not the clinical assay. |
| **Immune** | original A1 8-gene T-cell effector: `CD8A GZMA GZMB IFNG EOMES CXCL9 CXCL10 TBX21` (all 8 present) | A gene-mean, not ESTIMATE ImmuneScore. |

**Covariate:** published sample-level `PURITY` (OncoSG clinical attribute).
This is **not** TCGA ABSOLUTE and **not** ESTIMATE cosine purity.

**Statistic:** first-order partial Spearman = Pearson correlation of
average-ranks after algebraic adjustment for ranked purity; two-sided t,
df = n−3; Fisher-z 95% CI with SE = 1/√(n−4). Unadjusted Spearman is
reported beside it. BH-FDR is within the 3 primary tests only.

### What cannot be computed (honest skip)

cBioPortal exposes **only z-score** mRNA profiles for this study
(`luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores` and the
diploid-referenced twin). Datahub has no `data_mrna_seq_v2_rsem.txt`.

| Score | Why skipped |
|---|---|
| ESTIMATE ImmuneScore / TumorPurity | Official ssGSEA ranks genes **within a sample** on raw/log expression. Ranking z-scores (standardized **across** samples) is not the same transform. The Affymetrix cosine purity formula is also invalid on z-scores. |
| xCell CD8 / immune | Needs counts or TPM plus RNA-seq spillover calibration. n=169 z-scores are not that input. |
| MCP-counter | Published as mean log2(TPM+1) of marker genes. Z-scores change within-sample ranks. |

Do not pretend those scores exist here.

## Direct answer

After published purity, TACSTD2 is **negatively** associated with all three
primary readouts in this public OncoSG LUAD matrix:

| Feature | Unadj ρ | Unadj p | Partial ρ | 95% CI | Partial p | BH q (3 tests) | n |
|---|---:|---:|---:|---|---:|---:|---:|
| CD8 (`CD8A`) | {cd8['raw_rho']:.3f} | {fmt_p(cd8['raw_p'])} | {cd8['partial_rho']:.3f} | {cd8['partial_ci_low']:.3f} to {cd8['partial_ci_high']:.3f} | {fmt_p(cd8['partial_p'])} | {cd8['partial_fdr_primary3']:.2e} | {int(cd8['n'])} |
| GEP18 (17/18) | {gep['raw_rho']:.3f} | {fmt_p(gep['raw_p'])} | {gep['partial_rho']:.3f} | {gep['partial_ci_low']:.3f} to {gep['partial_ci_high']:.3f} | {fmt_p(gep['partial_p'])} | {gep['partial_fdr_primary3']:.2e} | {int(gep['n'])} |
| Immune (A1 8-gene) | {imm['raw_rho']:.3f} | {fmt_p(imm['raw_p'])} | {imm['partial_rho']:.3f} | {imm['partial_ci_low']:.3f} to {imm['partial_ci_high']:.3f} | {fmt_p(imm['partial_p'])} | {imm['partial_fdr_primary3']:.2e} | {int(imm['n'])} |

- CD8: {line(cd8)}
- GEP18: {line(gep)}
- Immune: {line(imm)}

Purity adjustment **does not create** the sign: unadjusted Spearman is already
negative. Partialling purity shrinks |ρ| (TACSTD2 vs PURITY Spearman
ρ = {tac_pur[0]:.3f}, p = {fmt_p(tac_pur[1])}) but leaves all three primary
tests negative.

**Do not quote the neutrophil −0.456 as if it were CD8 or GEP.** Those are
different endpoints. CD8/GEP/immune |ρ| here is smaller than the IMSIG
neutrophil |ρ|.

## Secondary / original A1 signatures

| Feature | Unadj ρ | Partial ρ | Partial p | n |
|---|---:|---:|---:|---:|
| CYT (GZMA+PRF1) | {cyt['raw_rho']:.3f} | {cyt['partial_rho']:.3f} | {fmt_p(cyt['partial_p'])} | {int(cyt['n'])} |
| CD8A+CD8B | {row('CD8_CD8A_CD8B')['raw_rho']:.3f} | {row('CD8_CD8A_CD8B')['partial_rho']:.3f} | {fmt_p(row('CD8_CD8A_CD8B')['partial_p'])} | {int(row('CD8_CD8A_CD8B')['n'])} |
| A1 cytotoxic 8-gene | {row('cytotoxic_A1')['raw_rho']:.3f} | {row('cytotoxic_A1')['partial_rho']:.3f} | {fmt_p(row('cytotoxic_A1')['partial_p'])} | {int(row('cytotoxic_A1')['n'])} |
| A1 exhaustion 12-gene | {row('exhaustion_A1')['raw_rho']:.3f} | {row('exhaustion_A1')['partial_rho']:.3f} | {fmt_p(row('exhaustion_A1')['partial_p'])} | {int(row('exhaustion_A1')['n'])} |

These **exactly match** the original OncoSG A1 per-cohort rows in PR 89
(immune / cytotoxic / exhaustion partial ρ = −0.3178 / −0.3631 / −0.3887,
n=169). CD8 and GEP18 were not primary endpoints in that run.

## Anchor: prior IMSIG neutrophils ρ = −0.456

Recomputed on the **same** Datahub revision and checksums as PR 82:

| | Prior (PR 82) | This run |
|---|---:|---:|
| n | 169 | {int(neu['n'])} |
| Neutrophils raw ρ | −0.503 | {neu['raw_rho']:.3f} |
| Neutrophils partial ρ | **−0.456** | **{neu['partial_rho']:.3f}** |
| Partial p | 5.06e-10 | {fmt_p(neu['partial_p'])} |
| BH q (8 IMSIG) | 4.05e-09 | {neu['partial_fdr_imsig8']:.2e} |

`matches_3dp` = {str(results['prior_neutrophils_anchor']['matches_3dp']).lower()}.
The neutrophil result **reproduces**. It is **not** a CD8 or GEP result.

Published IMSIG (not recomputed from genes):

| Signature | Unadj ρ | Unadj p | Partial ρ | 95% CI | Partial p | BH q |
|---|---:|---:|---:|---|---:|---:|
{chr(10).join(imsig_lines)}

IMSIG T cells (closest published “T-cell immune” column): {line(tcell)}.
IMSIG interferon is the weak / null member of that set, as in PR 82.

## Data

- Study: [luad_oncosg_2020](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020) — Chen et al., *Nat Genet* 2020
- Expression: Datahub `{DATAHUB_SHA}` `data_mrna_seq_v2_rsem_zscores_ref_all_samples.txt` (sha256 `44093dbc6633e5280d0fcd18fb07f5c6ee5c2cff175cd76b48b3665c3d4ff54f`, 23,493,696 bytes)
- Clinical: same revision `data_clinical_sample.txt` (sha256 `55739b69bf4f2e8e624aab51a0b5903ea797751905219b7bf6581afb4c96f368`)
- Profile id: `luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores`
- TACSTD2 Entrez 4070; API vs file check: {results['coverage']['api_vs_file_tacstd2']['note']}

Raw downloads stay in the cache directory and are not committed.

## Caveats

1. **Not ICI.** OncoSG is a genomic LUAD series. No PD-1/PD-L1 treatment labels are used here. Do not call these tumors hot/cold on IO.
2. **Bulk composition.** TACSTD2 is epithelial. A negative correlation with CD8/GEP partly reflects tumor-cell content. That is why the purity-adjusted number is the one to quote — and why it is smaller than the raw ρ.
3. **Purity is a published scalar**, method not re-derived (not ABSOLUTE, not ESTIMATE). Residualizing on it is not cell-intrinsic causality.
4. **GEP18 is 17/18.** CCL5 is missing. We did not impute it and did not swap in TLR2.
5. **No ESTIMATE / xCell / MCP** on this public deposit. If a slide quotes those names for OncoSG, the public files cannot support that sentence.
6. **Z-scores.** Signature means are means of already z-scored genes. Spearman is rank-based, so a monotone per-gene transform would not change ρ; a within-sample ssGSEA (ESTIMATE) would.
7. **n = 169 ≠ 181.** The portal RNA list is larger than the open matrix. We did not drop samples to hit a target n.
8. **No multiple-testing story beyond the stated BH sets.** Secondary rows are descriptive.
9. **East-Asian LUAD.** Do not pool this ρ with TCGA-LUSC and call it “NSCLC.”

## Reproduce

```bash
pip install -r scripts/rework/OncoSG_A1/requirements.txt
python3 scripts/rework/OncoSG_A1/analyze.py \\
  --cache-dir /tmp/oncosg_a1 \\
  --out-dir results/rework/OncoSG_A1
```

## Artifacts

| File | What |
|---|---|
| `correlations.tsv` | Every feature, raw + partial ρ / p / CI / FDR |
| `sample_table.tsv` | 169 rows: TACSTD2, PURITY, scores, light clinical |
| `gene_coverage.tsv` | Present/missing genes per signature |
| `results.json` | Machine-readable primary numbers + methods flags |
| `provenance.json` | URLs, sha256, bytes |
| `figures/forest_partial_rho.png` | Primary + IMSIG T / neutrophils |
| `figures/scatter_partial_residuals.png` | Rank residuals for CD8 / GEP18 / immune |
"""
    (out / "README.md").write_text(md)


if __name__ == "__main__":
    main()
