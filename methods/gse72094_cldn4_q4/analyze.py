#!/usr/bin/env python3
"""ADDITIVE GSE72094: CLDN4 Q4 vs Q1 vs CD8A / CD274.

Residual CLDN4 vs CD8A after ESTIMATE is already NS in
methods/gse72094_cldn4 (PR 302; n=442, partial ρ=−0.051, p=0.282).
This page does not re-claim that residual. Extra cut is Q4 vs Q1
(honest n = quartile arms, not 442) plus continuous CLDN4 vs CD274.

Public GEO series matrix + GPL15048 GeneSymbol only.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
CACHE = Path("/tmp/gse72094_cldn4_q4")
FIG = HERE / "figures"
TAB = HERE / "tables"
for d in (CACHE, FIG, TAB):
    d.mkdir(parents=True, exist_ok=True)

MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72094/"
    "matrix/GSE72094_series_matrix.txt.gz"
)
GPL_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
    "?acc=GPL15048&targ=self&form=text&view=data"
)
GENES = ["CLDN4", "CD8A", "CD274", "TACSTD2"]

# Known residual from methods/gse72094_cldn4 (not re-claimed).
KNOWN_RESIDUAL = {
    "test": "CLDN4 vs CD8A after ESTIMATEScore",
    "n": 442,
    "partial_rho": -0.051325397585842876,
    "partial_p": 0.281610353415495,
    "unadj_rho": -0.21902110384004267,
    "unadj_p": 3.3489187131883584e-06,
    "source": "methods/gse72094_cldn4 (PR 302)",
}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> Path:
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"fetch {url}")
    urllib.request.urlretrieve(url, tmp)
    tmp.replace(dest)
    print(f"  wrote {dest} ({dest.stat().st_size} bytes)")
    return dest


def fmt_p(p) -> str:
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def fmt_r(x, digits: int = 3) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    return f"{x:+.{digits}f}" if digits == 3 else f"{x:.{digits}f}"


def residualize(y, Z):
    x = np.column_stack([np.ones(len(y)), Z])
    b, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ b


def spearman_pair(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m = np.isfinite(x) & np.isfinite(y)
    n = int(m.sum())
    if n < 6 or np.nanstd(x[m]) == 0 or np.nanstd(y[m]) == 0:
        return {"n": n, "rho": np.nan, "p": np.nan}
    r, p = stats.spearmanr(x[m], y[m])
    return {"n": n, "rho": float(r), "p": float(p)}


def partial_spearman(x, y, z):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    z = np.asarray(z, float)
    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    n = int(m.sum())
    if n < 8:
        return {"n": n, "rho": np.nan, "p": np.nan}
    xr = stats.rankdata(x[m])
    yr = stats.rankdata(y[m])
    zr = stats.rankdata(z[m])
    rx = residualize(xr, zr.reshape(-1, 1))
    ry = residualize(yr, zr.reshape(-1, 1))
    r, _ = stats.pearsonr(rx, ry)
    r = float(np.clip(r, -0.999999, 0.999999))
    df = n - 3
    t = r * math.sqrt(df / max(1e-12, 1 - r * r))
    p = float(2 * stats.t.sf(abs(t), df))
    return {"n": n, "rho": r, "p": p}


def load_gpl(path: Path) -> pd.Series:
    skip = 0
    with open(path, "r", errors="replace") as f:
        for i, line in enumerate(f):
            if line.startswith("ID\t"):
                skip = i
                break
    df = pd.read_csv(path, sep="\t", skiprows=skip, dtype=str, low_memory=False)
    if "platform_table_end" in str(df.iloc[-1, 0]):
        df = df.iloc[:-1]
    s = df.set_index("ID")["GeneSymbol"].astype(str)
    s = s.replace({"nan": np.nan, "None": np.nan, "": np.nan})
    s = s.dropna()
    s = s.map(lambda x: str(x).split("///")[0].strip())
    return s[s.str.len() > 0]


def parse_meta_and_extract(matrix: Path, keep_probes: set[str]):
    meta_rows: dict[str, list[str]] = {}
    samples = None
    expr_start = None
    with gzip.open(matrix, "rt", errors="replace") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("!Sample_geo_accession"):
            samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_title"):
            meta_rows["title"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_source_name"):
            meta_rows["source"] = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
        if line.startswith("!Sample_characteristics"):
            vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            key = None
            for v in vals:
                if v and ":" in v:
                    key = v.split(":", 1)[0].strip().lower()
                    break
            if key is None:
                key = f"char_{len(meta_rows)}"
            cleaned = []
            for v in vals:
                cleaned.append(v.split(":", 1)[1].strip() if ":" in v else v)
            k = key
            n = 2
            while k in meta_rows:
                k = f"{key}_{n}"
                n += 1
            meta_rows[k] = cleaned
        if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
            expr_start = i
            break
    if samples is None or expr_start is None:
        raise SystemExit(f"could not parse GEO matrix {matrix}")
    meta = pd.DataFrame(meta_rows, index=samples)
    header = [x.strip().strip('"') for x in lines[expr_start].rstrip("\n").split("\t")]
    rows = []
    idx = []
    for line in lines[expr_start + 1 :]:
        if line.startswith("!"):
            break
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        pid = parts[0].strip().strip('"')
        if pid not in keep_probes:
            continue
        idx.append(pid)
        rows.append([float(x) if x not in ("", "NA", "null") else np.nan for x in parts[1:]])
    expr = pd.DataFrame(rows, index=idx, columns=header[1:])
    expr = expr.loc[:, samples]
    return meta, expr


def collapse_maxmean(expr: pd.DataFrame, probe2gene: pd.Series) -> pd.DataFrame:
    common = expr.index.intersection(probe2gene.index)
    g = probe2gene.loc[common]
    e = expr.loc[common]
    means = e.mean(axis=1)
    pick = means.groupby(g).idxmax()
    out = e.loc[pick.values]
    out.index = pick.index
    return out


def quartiles(s: pd.Series) -> pd.Series:
    out = pd.Series(index=s.index, dtype=object)
    ok = s.dropna()
    labels = pd.qcut(ok.rank(method="first"), 4, labels=["Q1", "Q2", "Q3", "Q4"])
    out.loc[ok.index] = labels.astype(str)
    return out


def mwu(a, b) -> dict:
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    n4, n1 = len(a), len(b)
    rec = {
        "n_q4": n4,
        "n_q1": n1,
        "median_q4": float(np.median(a)) if n4 else np.nan,
        "median_q1": float(np.median(b)) if n1 else np.nan,
        "delta_median": np.nan,
        "U": np.nan,
        "rank_biserial_q4_gt_q1": np.nan,
        "p": np.nan,
    }
    if n4 < 2 or n1 < 2:
        return rec
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = 2.0 * float(U) / (n4 * n1) - 1.0
    rec.update(
        {
            "median_q4": float(np.median(a)),
            "median_q1": float(np.median(b)),
            "delta_median": float(np.median(a) - np.median(b)),
            "U": float(U),
            "rank_biserial_q4_gt_q1": float(r),
            "p": float(p),
        }
    )
    return rec


def save_box(path: Path, q4, q1, title: str, ylab: str) -> None:
    fig, ax = plt.subplots(figsize=(3.6, 4.0))
    data = [np.asarray(q1, float), np.asarray(q4, float)]
    data = [d[~np.isnan(d)] for d in data]
    bp = ax.boxplot(
        data,
        tick_labels=[f"Q1\nn={len(data[0])}", f"Q4\nn={len(data[1])}"],
        widths=0.55,
        showfliers=False,
        patch_artist=True,
    )
    for patch, c in zip(bp["boxes"], ["#c9d6df", "#f6b26b"]):
        patch.set_facecolor(c)
        patch.set_alpha(0.85)
    rng = np.random.default_rng(0)
    for i, d in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.08, 0.08, size=len(d)), d, s=10, c="k", alpha=0.45, zorder=3)
    ax.set_ylabel(ylab)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def main() -> None:
    matrix = fetch(MATRIX_URL, CACHE / "GSE72094_series_matrix.txt.gz")
    plat = fetch(GPL_URL, CACHE / "GPL15048_datatable.txt")

    p2g = load_gpl(plat)
    keep_probes = set(p2g.index[p2g.isin(GENES)])
    print(f"GPL probes for {GENES}: {len(keep_probes)}")
    for g in GENES:
        print(f"  {g}: {(p2g == g).sum()} probes")

    print("parsing matrix (target probes only)")
    meta, probes = parse_meta_and_extract(matrix, keep_probes)
    print(f"  samples={meta.shape[0]} extracted_probes={probes.shape[0]}")
    print(f"  source unique: {sorted(meta['source'].astype(str).str.lower().unique())}")

    keep = meta["source"].astype(str).str.lower().str.contains("adenocarcinoma")
    n_luad = int(keep.sum())
    print(f"  LUAD filter: {n_luad} / {len(keep)}")
    meta_t = meta.loc[keep]
    probes_t = probes.loc[:, meta_t.index]

    genes = collapse_maxmean(probes_t, p2g.loc[p2g.isin(GENES)])
    print(f"  collapsed genes: {genes.index.tolist()}")
    missing = [g for g in GENES if g not in genes.index]
    if missing:
        raise SystemExit(f"missing genes after collapse: {missing}")

    st = pd.DataFrame({g: genes.loc[g] for g in GENES})
    st.index.name = "sample"
    for col in ["patient_id", "title", "kras_status", "egfr_status", "stk11_status", "tp53_status", "gender", "stage"]:
        if col in meta_t.columns:
            st[col] = meta_t[col]

    n_patient = int(st["patient_id"].nunique()) if "patient_id" in st.columns else n_luad
    n_dup_patient = int((st["patient_id"].value_counts() > 1).sum()) if "patient_id" in st.columns else 0
    print(f"  unique patient_id={n_patient} multi-sample patients={n_dup_patient}")

    harvested_candidates = [
        HERE / "harvested" / "estimate_scores.tsv",
        CACHE / "harvested_samples.tsv",
    ]
    harvested = next((p for p in harvested_candidates if p.exists()), None)
    if harvested is not None:
        hv = pd.read_csv(harvested, sep="\t").set_index("sample")
        est_cols = [c for c in ["ESTIMATE_Score", "ESTIMATE_ImmuneScore"] if c in hv.columns]
        joined = st.join(hv[est_cols], how="inner")
        if len(joined) == n_luad and "ESTIMATE_Score" in joined.columns:
            print(f"  harvest join n={len(joined)} from {harvested}")
            st["ESTIMATE_Score"] = joined["ESTIMATE_Score"]
            if "ESTIMATE_ImmuneScore" in joined.columns:
                st["ESTIMATE_ImmuneScore"] = joined["ESTIMATE_ImmuneScore"]
        else:
            print(f"  harvest join incomplete: {len(joined)} / {n_luad}")

    st["CLDN4_quartile"] = quartiles(st["CLDN4"])
    qcounts = st["CLDN4_quartile"].value_counts().to_dict()
    print(f"  CLDN4 quartiles: {qcounts}")

    # Continuous
    cont_rows = []
    for feat in ["CD8A", "CD274", "TACSTD2"]:
        rec = spearman_pair(st["CLDN4"], st[feat])
        rec.update({"predictor": "CLDN4", "endpoint": feat, "kind": "continuous Spearman"})
        if "ESTIMATE_Score" in st.columns:
            pr = partial_spearman(st["CLDN4"], st[feat], st["ESTIMATE_Score"])
            rec["n_partial"] = pr["n"]
            rec["partial_rho"] = pr["rho"]
            rec["partial_p"] = pr["p"]
            rec["covariate"] = "ESTIMATEScore (harvested ssGSEA, PR 302)"
        cont_rows.append(rec)
    cont = pd.DataFrame(cont_rows)

    # Q4 vs Q1
    q4_mask = st["CLDN4_quartile"] == "Q4"
    q1_mask = st["CLDN4_quartile"] == "Q1"
    q_rows = []
    for feat in ["CD8A", "CD274", "TACSTD2"]:
        rec = mwu(st.loc[q4_mask, feat], st.loc[q1_mask, feat])
        rec.update(
            {
                "predictor": "CLDN4 quartile",
                "endpoint": feat,
                "kind": "MWU Q4 vs Q1",
            }
        )
        q_rows.append(rec)
    qtab = pd.DataFrame(q_rows)

    # reverse CD8A quartiles on CLDN4 (supporting, not a second discovery)
    st["CD8A_quartile"] = quartiles(st["CD8A"])
    rev = mwu(
        st.loc[st["CD8A_quartile"] == "Q4", "CLDN4"],
        st.loc[st["CD8A_quartile"] == "Q1", "CLDN4"],
    )
    rev.update({"predictor": "CD8A quartile", "endpoint": "CLDN4", "kind": "MWU Q4 vs Q1 (reverse)"})
    qtab = pd.concat([qtab, pd.DataFrame([rev])], ignore_index=True)

    cont.to_csv(TAB / "continuous.tsv", sep="\t", index=False)
    qtab.to_csv(TAB / "q4_vs_q1.tsv", sep="\t", index=False)
    st.to_csv(TAB / "samples.tsv", sep="\t")

    n_q4 = int(qcounts.get("Q4", 0))
    n_q1 = int(qcounts.get("Q1", 0))
    n_q2 = int(qcounts.get("Q2", 0))
    n_q3 = int(qcounts.get("Q3", 0))
    nna = {g: int(st[g].notna().sum()) for g in GENES}

    ntab = pd.DataFrame(
        [
            {"item": "GEO series matrix columns", "n": int(meta.shape[0])},
            {"item": "source_name contains adenocarcinoma", "n": n_luad},
            {"item": "unique patient_id", "n": n_patient},
            {"item": "patients with >1 array", "n": n_dup_patient},
            {"item": "CLDN4 non-NA", "n": nna["CLDN4"]},
            {"item": "CD8A non-NA", "n": nna["CD8A"]},
            {"item": "CD274 non-NA", "n": nna["CD274"]},
            {"item": "TACSTD2 non-NA", "n": nna["TACSTD2"]},
            {"item": "CLDN4 Q1", "n": n_q1},
            {"item": "CLDN4 Q2", "n": n_q2},
            {"item": "CLDN4 Q3", "n": n_q3},
            {"item": "CLDN4 Q4", "n": n_q4},
            {"item": "Q4 vs Q1 used", "n": f"{n_q4} vs {n_q1}"},
            {"item": "Do not write this n for Q4 vs Q1", "n": n_luad},
        ]
    )
    ntab.to_csv(TAB / "n_table.tsv", sep="\t", index=False)

    def crow(endpoint):
        return cont.loc[cont.endpoint == endpoint].iloc[0]

    def qrow(endpoint, predictor="CLDN4 quartile"):
        return qtab.loc[(qtab.endpoint == endpoint) & (qtab.predictor == predictor)].iloc[0]

    c_cd8 = crow("CD8A")
    c_pdl1 = crow("CD274")
    c_tac = crow("TACSTD2")
    q_cd8 = qrow("CD8A")
    q_pdl1 = qrow("CD274")
    q_tac = qrow("TACSTD2")
    q_rev = qrow("CLDN4", "CD8A quartile")

    # Residual reprint: prefer this-run partial if ESTIMATE joined, else known.
    if "partial_rho" in c_cd8.index and np.isfinite(c_cd8.partial_rho):
        res_rho, res_p, res_n = float(c_cd8.partial_rho), float(c_cd8.partial_p), int(c_cd8.n_partial)
    else:
        res_rho, res_p, res_n = KNOWN_RESIDUAL["partial_rho"], KNOWN_RESIDUAL["partial_p"], KNOWN_RESIDUAL["n"]

    def holds_q(row, alpha=0.05):
        return "yes" if row.p < alpha else "no"

    verdict = pd.DataFrame(
        [
            {
                "cut": "CLDN4 vs CD8A after ESTIMATEScore (known residual)",
                "n": res_n,
                "metric": "partial Spearman ρ",
                "effect": res_rho,
                "p": res_p,
                "holds": "no (already NS)",
            },
            {
                "cut": "CLDN4 vs CD8A unadjusted (known, not extra)",
                "n": int(c_cd8.n),
                "metric": "Spearman ρ",
                "effect": float(c_cd8.rho),
                "p": float(c_cd8.p),
                "holds": "already known",
            },
            {
                "cut": "CLDN4 vs CD274 (continuous extra)",
                "n": int(c_pdl1.n),
                "metric": "Spearman ρ",
                "effect": float(c_pdl1.rho),
                "p": float(c_pdl1.p),
                "holds": "yes" if c_pdl1.p < 0.05 else "no",
            },
            {
                "cut": "CLDN4 Q4 vs Q1, endpoint CD8A",
                "n": f"{int(q_cd8.n_q4)} vs {int(q_cd8.n_q1)}",
                "metric": "rank-biserial (Q4>Q1)",
                "effect": float(q_cd8.rank_biserial_q4_gt_q1),
                "p": float(q_cd8.p),
                "holds": holds_q(q_cd8),
            },
            {
                "cut": "CLDN4 Q4 vs Q1, endpoint CD274",
                "n": f"{int(q_pdl1.n_q4)} vs {int(q_pdl1.n_q1)}",
                "metric": "rank-biserial (Q4>Q1)",
                "effect": float(q_pdl1.rank_biserial_q4_gt_q1),
                "p": float(q_pdl1.p),
                "holds": holds_q(q_pdl1),
            },
        ]
    )
    verdict.to_csv(TAB / "verdict.tsv", sep="\t", index=False)

    save_box(
        FIG / "cldn4_q4q1_cd8a.png",
        st.loc[q4_mask, "CD8A"],
        st.loc[q1_mask, "CD8A"],
        f"GSE72094 CLDN4 Q4 vs Q1 · CD8A\nn={n_q4} vs {n_q1}",
        "CD8A (IRON log2)",
    )
    save_box(
        FIG / "cldn4_q4q1_cd274.png",
        st.loc[q4_mask, "CD274"],
        st.loc[q1_mask, "CD274"],
        f"GSE72094 CLDN4 Q4 vs Q1 · CD274\nn={n_q4} vs {n_q1}",
        "CD274 (IRON log2)",
    )

    # CD274 residual if ESTIMATE present (supporting; not the Q4 extra)
    pdl1_partial_note = ""
    if "partial_rho" in c_pdl1.index and np.isfinite(c_pdl1.get("partial_rho", np.nan)):
        pdl1_partial_note = (
            f"Supporting residual (not the Q4 extra): after harvested ESTIMATEScore, "
            f"CLDN4 vs CD274 flips to partial ρ={fmt_r(c_pdl1.partial_rho)} "
            f"p={fmt_p(c_pdl1.partial_p)} (n={int(c_pdl1.n_partial)}). "
            "Unadjusted and Q4 vs Q1 CD274 stay null, so this page does not upgrade "
            "the residual sign-flip into a PD-L1 claim."
        )

    holds_cd8_q = q_cd8.p < 0.05
    holds_pdl1_q = q_pdl1.p < 0.05
    holds_pdl1_c = c_pdl1.p < 0.05

    what_holds = []
    if holds_cd8_q:
        what_holds.append(
            f"CLDN4-high (Q4) tumors have lower **CD8A** than CLDN4-low (Q1), "
            f"n=**{n_q4} vs {n_q1}**, Δ median {fmt_r(q_cd8.delta_median)}, "
            f"rank-biserial {fmt_r(q_cd8.rank_biserial_q4_gt_q1)}, p={fmt_p(q_cd8.p)}. "
            "That is the same unadjusted immune-low direction as the known continuous CD8A row, "
            "not a residual-independent axis."
        )
    if holds_pdl1_q or holds_pdl1_c:
        bits = []
        if holds_pdl1_c:
            bits.append(
                f"continuous CLDN4 vs **CD274** n={int(c_pdl1.n)} ρ={fmt_r(c_pdl1.rho)} p={fmt_p(c_pdl1.p)}"
            )
        if holds_pdl1_q:
            bits.append(
                f"Q4 vs Q1 CD274 n=**{n_q4} vs {n_q1}** r={fmt_r(q_pdl1.rank_biserial_q4_gt_q1)} p={fmt_p(q_pdl1.p)}"
            )
        what_holds.append("Extra CD274: " + "; ".join(bits) + ".")
    if not what_holds:
        what_holds.append("Neither Q4 vs Q1 extra cut (CD8A, CD274) is significant at the honest quartile n.")

    what_not = [
        f"The already-reported residual **CLDN4 vs CD8A after ESTIMATEScore remains not significant** "
        f"(n={res_n}, partial ρ={fmt_r(res_rho)}, p={fmt_p(res_p)}). This PR does not upgrade that residual."
    ]
    if not holds_pdl1_c:
        what_not.append(
            f"Continuous CLDN4 vs **CD274** is null at n={int(c_pdl1.n)} "
            f"(ρ={fmt_r(c_pdl1.rho)}, p={fmt_p(c_pdl1.p)})."
        )
    if not holds_pdl1_q:
        what_not.append(
            f"CLDN4 Q4 vs Q1 vs **CD274** is not significant at n={n_q4} vs {n_q1} "
            f"(r={fmt_r(q_pdl1.rank_biserial_q4_gt_q1)}, p={fmt_p(q_pdl1.p)})."
        )
    what_not.append(
        f"Do not write n={n_luad} for the Q4 vs Q1 tests. Those tests drop Q2+Q3 and use **{n_q4} vs {n_q1}**."
    )

    finding = f"""# Finding — GSE72094 CLDN4 Q4 vs Q1 vs CD8A / CD274

**Additive extra cut only.** Schabath / Moffitt **GSE72094** (PMID 26477306): resected **lung adenocarcinoma**, **GPL15048** HuRSTA, author IRON-normalized series matrix. All matrix samples are `source_name = lung adenocarcinoma`. Patient = array: **{n_patient} unique `patient_id` / {n_luad} arrays** (multi-sample patients = {n_dup_patient}).

**Already known (PR 302 / `methods/gse72094_cldn4`, not re-claimed):** continuous CLDN4 vs CD8A, n=442, Spearman ρ=**{fmt_r(c_cd8.rho)}**, p=**{fmt_p(c_cd8.p)}**. After ESTIMATEScore residual the same pair is **not significant** (partial ρ=**{fmt_r(res_rho)}**, p=**{fmt_p(res_p)}**).

Two extra questions:

1. Continuous **CLDN4 vs CD274** (CD274 was not a primary row on the residual page).
2. **Q4 vs Q1** — CLDN4 quartiles on **CD8A** and **CD274**. Honest n is the quartile arms, not 442.

Reproduce: `python3 methods/gse72094_cldn4_q4/analyze.py` (GEO files cached under `/tmp/gse72094_cldn4_q4/`).

---

## Verdict

| Cut | n | Metric | Effect | p | Holds? |
|---|---:|---|---:|---:|---|
| CLDN4 vs CD8A after ESTIMATEScore (known residual) | {res_n} | partial Spearman ρ | {fmt_r(res_rho)} | {fmt_p(res_p)} | **no (already NS)** |
| CLDN4 vs CD8A unadjusted (known, not extra) | {int(c_cd8.n)} | Spearman ρ | {fmt_r(c_cd8.rho)} | {fmt_p(c_cd8.p)} | already known |
| **CLDN4 vs CD274** (continuous extra) | {int(c_pdl1.n)} | Spearman ρ | {fmt_r(c_pdl1.rho)} | {fmt_p(c_pdl1.p)} | **{"yes" if holds_pdl1_c else "no"}** |
| **CLDN4 Q4 vs Q1, endpoint CD8A** | **{n_q4} vs {n_q1}** | rank-biserial (Q4>Q1) | {fmt_r(q_cd8.rank_biserial_q4_gt_q1)} | {fmt_p(q_cd8.p)} | **{holds_q(q_cd8)}** |
| **CLDN4 Q4 vs Q1, endpoint CD274** | **{n_q4} vs {n_q1}** | rank-biserial (Q4>Q1) | {fmt_r(q_pdl1.rank_biserial_q4_gt_q1)} | {fmt_p(q_pdl1.p)} | **{holds_q(q_pdl1)}** |

**What holds.** {" ".join(what_holds)}

**What does not hold.** {" ".join(what_not)}

---

## Honest n

| Item | n |
|---|---:|
| GEO series / matrix columns | **{int(meta.shape[0])}** |
| `source_name` = lung adenocarcinoma | **{n_luad}** |
| Unique `patient_id` | **{n_patient}** |
| Patients with >1 array | **{n_dup_patient}** |
| CLDN4 / CD8A / CD274 / TACSTD2 non-NA | **{nna["CLDN4"]} / {nna["CD8A"]} / {nna["CD274"]} / {nna["TACSTD2"]}** |
| CLDN4 Q1 / Q2 / Q3 / Q4 | **{n_q1} / {n_q2} / {n_q3} / {n_q4}** |
| Q4 vs Q1 used | **{n_q4} vs {n_q1}** |

Quartiles are `pd.qcut(rank(method="first"), 4)` among the {n_luad} tumors with non-NA CLDN4. Do not write n={n_luad} for the quartile tests.

---

## Extra continuous — CLDN4 vs CD274

Spearman on all {n_luad} tumors. Author IRON log2.

| Endpoint | n | ρ | p | Kind |
|---|---:|---:|---:|---|
| CD8A | {int(c_cd8.n)} | {fmt_r(c_cd8.rho)} | {fmt_p(c_cd8.p)} | known residual-page unadjusted |
| **CD274** | {int(c_pdl1.n)} | **{fmt_r(c_pdl1.rho)}** | **{fmt_p(c_pdl1.p)}** | extra |
| TACSTD2 (companion) | {int(c_tac.n)} | {fmt_r(c_tac.rho)} | {fmt_p(c_tac.p)} | companion only |

{pdl1_partial_note if pdl1_partial_note else "ESTIMATE residual for CD274 was not recomputed in this run; CD8 residual is the PR 302 reprint."}

---

## Extra Q4 vs Q1

Two-sided Mann–Whitney U. Rank-biserial r = 2U/(n4 n1) − 1 (positive = Q4 higher).

### CLDN4 quartiles

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| **CD8A** | **{int(q_cd8.n_q4)} / {int(q_cd8.n_q1)}** | {q_cd8.median_q4:.3f} | {q_cd8.median_q1:.3f} | {fmt_r(q_cd8.delta_median)} | {q_cd8.U:.0f} | **{fmt_r(q_cd8.rank_biserial_q4_gt_q1)}** | **{fmt_p(q_cd8.p)}** |
| **CD274** | **{int(q_pdl1.n_q4)} / {int(q_pdl1.n_q1)}** | {q_pdl1.median_q4:.3f} | {q_pdl1.median_q1:.3f} | {fmt_r(q_pdl1.delta_median)} | {q_pdl1.U:.0f} | **{fmt_r(q_pdl1.rank_biserial_q4_gt_q1)}** | **{fmt_p(q_pdl1.p)}** |
| TACSTD2 (companion) | {int(q_tac.n_q4)} / {int(q_tac.n_q1)} | {q_tac.median_q4:.3f} | {q_tac.median_q1:.3f} | {fmt_r(q_tac.delta_median)} | {q_tac.U:.0f} | {fmt_r(q_tac.rank_biserial_q4_gt_q1)} | {fmt_p(q_tac.p)} |

### CD8A quartiles (endpoint = CLDN4; supporting)

| Endpoint | n Q4 / Q1 | median Q4 | median Q1 | Δ (Q4−Q1) | U | r | p |
|---|---:|---:|---:|---:|---:|---:|---:|
| CLDN4 | **{int(q_rev.n_q4)} / {int(q_rev.n_q1)}** | {q_rev.median_q4:.3f} | {q_rev.median_q1:.3f} | {fmt_r(q_rev.delta_median)} | {q_rev.U:.0f} | {fmt_r(q_rev.rank_biserial_q4_gt_q1)} | {fmt_p(q_rev.p)} |

The reverse cut is the same sign as the known continuous CD8A row. Do not treat it as a second independent discovery.

---

## Methods (this slice)

- **Matrix:** GEO `GSE72094_series_matrix.txt.gz` (author-processed IRON + RNA-quality batch correction). Probe × sample values used as published.
- **Annotation:** GPL15048 `GeneSymbol` (NCBI platform table). First symbol if `///`. **Max-mean** probe collapse to HUGO, same rule as `scripts/gse72094_cldn4.py`.
- **CD8** = `CD8A`. **CD274** = PD-L1 RNA (not protein). TACSTD2 is a companion only.
- **Quartiles:** `pd.qcut(rank(method="first"), 4)` on the {n_luad} non-NA CLDN4 values.
- **Continuous test:** two-sided Spearman.
- **Quartile test:** two-sided Mann–Whitney U, rank-biserial as above.
- **Residual reprint:** ESTIMATEScore from `harvested/estimate_scores.tsv` (PR 302 ssGSEA). This page does not re-fit ESTIMATE.
- Surgical LUAD. No ICI-response labels on this series.

---

## What this does not claim

- It does not re-open the PR 302 residual (already NS) as a new finding.
- It does not claim CLDN4 vs PD-L1 **protein**. CD274 is microarray RNA.
- It does not upgrade the supporting residual sign-flip (CLDN4 vs CD274) into a PD-L1 claim. Unadjusted and Q4 vs Q1 CD274 are null.
- It does not treat n={n_luad} as the Q4 vs Q1 n.
- Q4 vs Q1 CD8A is a coarsened version of the known unadjusted continuous test, not a new cohort and not a purity-residual claim.

---

## Files

- `FINDING.md` — this page
- `analyze.py` — download-once, score, write tables/figures
- `tables/q4_vs_q1.tsv` — MWU extra cuts
- `tables/continuous.tsv` — Spearman extra + known CD8A reprint
- `tables/n_table.tsv` — honest n
- `tables/verdict.tsv`
- `tables/samples.tsv` — per-array genes + quartiles
- `harvested/estimate_scores.tsv` — PR 302 ESTIMATEScore (residual reprint only)
- `figures/cldn4_q4q1_cd8a.png`, `cldn4_q4q1_cd274.png`
"""
    (HERE / "FINDING.md").write_text(finding)

    summary = {
        "cohort": "GSE72094",
        "pmid": "26477306",
        "n_luad": n_luad,
        "n_patient": n_patient,
        "n_q4": n_q4,
        "n_q1": n_q1,
        "known_residual_cd8": KNOWN_RESIDUAL,
        "this_run_residual_cd8": {"rho": res_rho, "p": res_p, "n": res_n},
        "continuous": cont.to_dict(orient="records"),
        "q4_vs_q1": qtab.to_dict(orient="records"),
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "GSE72094_series_matrix.txt.gz": {"bytes": matrix.stat().st_size, "md5": md5(matrix), "url": MATRIX_URL},
            "GPL15048_datatable.txt": {"bytes": plat.stat().st_size, "md5": md5(plat), "url": GPL_URL},
        },
    }
    with open(TAB / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)

    print(finding)
    print(f"Wrote {HERE}")


if __name__ == "__main__":
    main()
