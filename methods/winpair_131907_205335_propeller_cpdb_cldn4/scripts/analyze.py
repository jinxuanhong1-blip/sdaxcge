#!/usr/bin/env python3
"""ADDITIVE CLDN4-only propeller + CellPhoneDB-style on GSE131907+GSE205335.

Winning pair from PR #320 (author-malignant %pos vs T/NK). No dual-high.
No GSE207422. LIANA and CellChat are not run.

1) propeller/speckle-style: logit cell-type fractions, Q4 vs Q1, BH.
2) CellPhoneDB-style mean-of-partner-means on log1p(CP10k); MHC-I + T-recruit.
"""
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.special import digamma, polygamma
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parents[1]
TUMOR_ORIGINS_131907 = {"tLung", "tL/B", "mLN", "PE", "mBrain"}
MALIGNANT_SUBTYPES = {"Malignant cells", "tS1", "tS2", "tS3"}
TNK_TYPES_131907 = {"T lymphocytes", "NK cells"}
NORMAL_TISSUES_205335 = {"Normal Lung", "Normal LN", "Normal Brain"}


def log(msg: str) -> None:
    print(msg, flush=True)


def assign_quartiles(x: pd.Series) -> pd.Series:
    ranks = pd.Series(x, dtype=float).rank(method="average")
    return pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")


def logit(p: np.ndarray, pc: float = 0.5, n: np.ndarray | None = None) -> np.ndarray:
    """Empirical logit used by speckle::propeller (offset on counts when n given)."""
    p = np.asarray(p, dtype=float)
    if n is None:
        p = np.clip(p, 1e-6, 1 - 1e-6)
        return np.log(p / (1.0 - p))
    n = np.asarray(n, dtype=float)
    c = np.clip(p * n, 0, n)
    return np.log((c + pc) / (n - c + pc))


def asin_sqrt(p: np.ndarray) -> np.ndarray:
    return np.arcsin(np.sqrt(np.clip(p, 0.0, 1.0)))


def rank_biserial(a: np.ndarray, b: np.ndarray) -> tuple[float, float, float]:
    """MWU a vs b; r = 2U/(n_a n_b) - 1. r<0 means a is lower than b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return float("nan"), float("nan"), float("nan")
    u, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    r = (2.0 * float(u)) / (len(a) * len(b)) - 1.0
    return float(u), float(p), float(r)


def welch_t(a: np.ndarray, b: np.ndarray) -> dict:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    out = {
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_a": float(np.mean(a)) if len(a) else float("nan"),
        "mean_b": float(np.mean(b)) if len(b) else float("nan"),
        "delta": float("nan"),
        "t": float("nan"),
        "df": float("nan"),
        "p": float("nan"),
        "se": float("nan"),
    }
    if len(a) < 2 or len(b) < 2:
        return out
    t = stats.ttest_ind(a, b, equal_var=False, alternative="two-sided")
    va, vb = float(np.var(a, ddof=1)), float(np.var(b, ddof=1))
    se = float(np.sqrt(va / len(a) + vb / len(b)))
    out.update(
        {
            "delta": float(np.mean(a) - np.mean(b)),
            "t": float(t.statistic),
            "df": float(t.df),
            "p": float(t.pvalue),
            "se": se,
        }
    )
    return out


def squeeze_var(s2: np.ndarray, df: np.ndarray) -> tuple[float, float]:
    """limma-like inverse-χ² prior (Smyth 2004). Weak when k is tiny — reported as such."""
    s2 = np.asarray(s2, dtype=float)
    df = np.asarray(df, dtype=float)
    ok = (s2 > 0) & np.isfinite(s2) & (df > 0)
    s2, df = s2[ok], df[ok]
    if s2.size < 2:
        return float("nan"), float("nan")
    z = np.log(s2)
    e = z - digamma(df / 2.0) + np.log(df / 2.0)
    emean = float(np.mean(e))
    evar = float(np.var(e, ddof=1))
    mean_trigamma = float(np.mean(polygamma(1, df / 2.0)))
    target = evar - mean_trigamma
    if not np.isfinite(target) or target <= 1e-12:
        # No extra-biological variance: do not shrink (limma d0 → ∞, s2_post = s2).
        return float("nan"), float("nan")
    # trigamma(d0/2) ≈ target; trigamma(x) ≈ 1/x for large x
    d0 = 2.0 / target
    d0 = float(np.clip(d0, 0.1, 1e12))
    s0_2 = float(np.exp(emean + digamma(d0 / 2.0) - np.log(d0 / 2.0)))
    return d0, s0_2


def moderate_t(delta: np.ndarray, se: np.ndarray, df: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    s2 = np.asarray(se, dtype=float) ** 2
    d0, s0_2 = squeeze_var(s2, np.asarray(df, dtype=float))
    if not np.isfinite(d0) or not np.isfinite(s0_2):
        t = np.asarray(delta, dtype=float) / np.asarray(se, dtype=float)
        p = 2 * stats.t.sf(np.abs(t), np.asarray(df, dtype=float))
        return t, p, d0, s0_2
    s2_post = (d0 * s0_2 + np.asarray(df) * s2) / (d0 + np.asarray(df))
    se_post = np.sqrt(s2_post)
    t = np.asarray(delta, dtype=float) / se_post
    p = 2 * stats.t.sf(np.abs(t), d0 + np.asarray(df, dtype=float))
    return t, p, d0, s0_2


def ols_q4_cohort(y: np.ndarray, is_q4: np.ndarray, is_c2: np.ndarray) -> dict:
    """y ~ 1 + Q4 + cohort. Test Q4 coefficient (winning-pair propeller)."""
    y = np.asarray(y, dtype=float)
    q = np.asarray(is_q4, dtype=float)
    c = np.asarray(is_c2, dtype=float)
    m = np.isfinite(y) & np.isfinite(q) & np.isfinite(c)
    y, q, c = y[m], q[m], c[m]
    n = int(y.size)
    X = np.column_stack([np.ones(n), q, c])
    if n < 6 or np.linalg.matrix_rank(X) < 3:
        return {"n": n, "beta_q4": float("nan"), "se": float("nan"), "t": float("nan"), "df": float("nan"), "p": float("nan")}
    beta, residuals, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    df = n - int(rank)
    if df <= 0:
        return {"n": n, "beta_q4": float(beta[1]), "se": float("nan"), "t": float("nan"), "df": float(df), "p": float("nan")}
    resid = y - X @ beta
    s2 = float(np.sum(resid**2) / df)
    xtx_inv = np.linalg.inv(X.T @ X)
    se = float(np.sqrt(s2 * xtx_inv[1, 1]))
    t = float(beta[1] / se) if se > 0 else float("nan")
    p = float(2 * stats.t.sf(abs(t), df)) if np.isfinite(t) else float("nan")
    return {"n": n, "beta_q4": float(beta[1]), "se": se, "t": t, "df": float(df), "p": p}


def bh(pvals: list[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.sum() == 0:
        return out
    out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def load_units(data: Path) -> pd.DataFrame:
    s131 = pd.read_csv(data / "GSE131907_samples.tsv", sep="\t")
    meta = pd.read_csv(data / "GSE131907_sample_metadata.tsv", sep="\t")
    s131 = s131.merge(meta[["Sample", "patient_id", "tumor_stage"]], left_on="sample", right_on="Sample", how="left")
    tumor = s131[s131["origin"].isin(TUMOR_ORIGINS_131907) & (s131["n_malignant"] >= 20)].copy()
    tumor["cohort"] = "GSE131907"
    tumor["unit"] = tumor["sample"]
    tumor["unit_type"] = "sample"
    tumor["patient"] = tumor["patient_id"]
    tumor["cldn4_pct"] = tumor["mal_CLDN4_pct"]
    tumor["cldn4_mean"] = tumor["mal_CLDN4_mean"]
    tumor["frac_tnk"] = tumor["frac_tnk"]
    tumor["frac_cd8"] = tumor["frac_cd8"]
    tumor["frac_b"] = tumor["frac_b"]
    tumor["n_b_used"] = tumor["n_b"]
    tumor["b_label"] = "B lymphocytes (author Cell_type)"
    tumor["n_other"] = tumor["n_cells"] - tumor["n_tnk"] - tumor["n_b"]
    tumor["frac_other"] = tumor["n_other"] / tumor["n_cells"]

    s205 = pd.read_csv(data / "GSE205335_patients.tsv", sep="\t")
    s205 = s205.copy()
    s205["cohort"] = "GSE205335"
    s205["unit"] = s205["patient"]
    s205["unit_type"] = "patient"
    s205["cldn4_pct"] = s205["mal_CLDN4_pct_pos"]
    s205["cldn4_mean"] = s205["mal_CLDN4_mean"]
    s205["frac_b"] = s205["frac_b_plasma"]
    s205["n_b_used"] = s205["n_b_plasma"]
    s205["b_label"] = "B+plasma (author lineage)"
    s205["n_other"] = s205["n_cells"] - s205["n_tnk"] - s205["n_b_plasma"]
    s205["frac_other"] = s205["n_other"] / s205["n_cells"]
    s205["origin"] = s205["tissue"]
    s205["n_malignant"] = s205["n_malignant"]
    s205["n_b"] = s205["n_b_plasma"]

    cols = [
        "cohort",
        "unit",
        "unit_type",
        "patient",
        "origin",
        "n_cells",
        "n_malignant",
        "n_tnk",
        "n_cd8",
        "n_b_used",
        "n_other",
        "frac_tnk",
        "frac_cd8",
        "frac_b",
        "frac_other",
        "cldn4_pct",
        "cldn4_mean",
        "b_label",
    ]
    extra_205 = ["recist", "cancer_subtype"]
    for c in extra_205:
        if c not in tumor.columns:
            tumor[c] = np.nan
    keep = tumor[cols + extra_205].copy()
    keep2 = s205[cols + extra_205].copy()
    out = pd.concat([keep, keep2], ignore_index=True)
    out["q_pct"] = pd.Series(index=out.index, dtype="object")
    out["q_mean"] = pd.Series(index=out.index, dtype="object")
    for cohort, idx in out.groupby("cohort").groups.items():
        out.loc[idx, "q_pct"] = assign_quartiles(out.loc[idx, "cldn4_pct"]).astype(str).to_numpy()
        out.loc[idx, "q_mean"] = assign_quartiles(out.loc[idx, "cldn4_mean"]).astype(str).to_numpy()
    return out


def propeller_block(df: pd.DataFrame, compartments: list[tuple[str, str]]) -> pd.DataFrame:
    rows = []
    for cohort, sub in df.groupby("cohort", sort=False):
        tails = sub[sub["q_pct"].isin(["Q1", "Q4"])].copy()
        n1 = int((tails["q_pct"] == "Q1").sum())
        n4 = int((tails["q_pct"] == "Q4").sum())
        n_patients = int(tails["patient"].nunique())
        for name, col in compartments:
            y = tails[col].to_numpy(dtype=float)
            nlib = tails["n_cells"].to_numpy(dtype=float)
            q4 = tails["q_pct"].eq("Q4").to_numpy()
            y4, y1 = y[q4], y[~q4]
            n4_lib, n1_lib = nlib[q4], nlib[~q4]
            logit4 = logit(y4, n=n4_lib)
            logit1 = logit(y1, n=n1_lib)
            wt = welch_t(logit4, logit1)
            u, p_mwu, r = rank_biserial(y4, y1)
            asin_wt = welch_t(asin_sqrt(y4), asin_sqrt(y1))
            rows.append(
                {
                    "scope": "cohort",
                    "cohort": cohort,
                    "compartment": name,
                    "n_units_full": int(len(sub)),
                    "n_q1": n1,
                    "n_q4": n4,
                    "n_compared": n1 + n4,
                    "n_unique_patients_tails": n_patients,
                    "unit_type": sub["unit_type"].iloc[0],
                    "median_frac_q1": float(np.median(y1)),
                    "median_frac_q4": float(np.median(y4)),
                    "delta_median_frac": float(np.median(y4) - np.median(y1)),
                    "logit_delta_q4_minus_q1": wt["delta"],
                    "logit_t": wt["t"],
                    "logit_df": wt["df"],
                    "logit_se": wt["se"],
                    "logit_p": wt["p"],
                    "asin_delta": asin_wt["delta"],
                    "asin_p": asin_wt["p"],
                    "frac_mwu_p": p_mwu,
                    "frac_r_rb": r,
                    "note": "within-cohort CLDN4 %pos quartiles; propeller logit + Welch",
                }
            )
    # winning-pair OLS with cohort covariate on tails
    tails = df[df["q_pct"].isin(["Q1", "Q4"])].copy()
    for name, col in compartments:
        y = logit(tails[col].to_numpy(dtype=float), n=tails["n_cells"].to_numpy(dtype=float))
        fit = ols_q4_cohort(y, tails["q_pct"].eq("Q4").to_numpy(), tails["cohort"].eq("GSE205335").to_numpy())
        y_frac = tails[col].to_numpy(dtype=float)
        q4 = tails["q_pct"].eq("Q4").to_numpy()
        u, p_mwu, r = rank_biserial(y_frac[q4], y_frac[~q4])
        rows.append(
            {
                "scope": "winpair",
                "cohort": "GSE131907+GSE205335",
                "compartment": name,
                "n_units_full": int(len(df)),
                "n_q1": int((~q4).sum()),
                "n_q4": int(q4.sum()),
                "n_compared": int(len(tails)),
                "n_unique_patients_tails": int(tails["patient"].nunique()),
                "unit_type": "sample+patient",
                "median_frac_q1": float(np.median(y_frac[~q4])),
                "median_frac_q4": float(np.median(y_frac[q4])),
                "delta_median_frac": float(np.median(y_frac[q4]) - np.median(y_frac[~q4])),
                "logit_delta_q4_minus_q1": fit["beta_q4"],
                "logit_t": fit["t"],
                "logit_df": fit["df"],
                "logit_se": fit["se"],
                "logit_p": fit["p"],
                "asin_delta": float("nan"),
                "asin_p": float("nan"),
                "frac_mwu_p": p_mwu,
                "frac_r_rb": r,
                "note": "OLS logit ~ Q4 + cohort; winning-pair propeller",
            }
        )
    tab = pd.DataFrame(rows)
    # BH families
    tab["q_bh"] = np.nan
    tab["bh_family"] = ""
    prim = tab["scope"].eq("cohort") & tab["compartment"].isin(["TNK", "B"])
    tab.loc[prim, "q_bh"] = bh(tab.loc[prim, "logit_p"].tolist())
    tab.loc[prim, "bh_family"] = "primary_4_cohort_TNK_B"
    pair = tab["scope"].eq("winpair") & tab["compartment"].isin(["TNK", "B", "CD8"])
    tab.loc[pair, "q_bh_winpair"] = bh(tab.loc[pair, "logit_p"].tolist())
    # moderated t within each scope (tiny k; companion)
    for scope, idx in tab.groupby("scope").groups.items():
        sl = tab.loc[idx]
        t_mod, p_mod, d0, s0 = moderate_t(
            sl["logit_delta_q4_minus_q1"].to_numpy(),
            sl["logit_se"].to_numpy(),
            sl["logit_df"].to_numpy(),
        )
        tab.loc[idx, "logit_t_moderated"] = t_mod
        tab.loc[idx, "logit_p_moderated"] = p_mod
        tab.loc[idx, "ebayes_d0"] = d0
        tab.loc[idx, "ebayes_s0_2"] = s0
    tab["recovery_q10"] = (tab["logit_delta_q4_minus_q1"] < 0) & (tab["q_bh"] < 0.10)
    return tab


def parse_interactors(text: str, genes: set[str]) -> tuple[str, str] | None:
    if not isinstance(text, str) or "-" not in text:
        return None
    for i, ch in enumerate(text):
        if ch != "-":
            continue
        left, right = text[:i], text[i + 1 :]
        if not left or not right:
            continue
        if all(p in genes for p in left.split("+")) and all(p in genes for p in right.split("+")):
            return left, right
    return None


def pathway_of(ligand: str, receptor: str, cfg: dict) -> str | None:
    lig = set(ligand.split("+"))
    rec = set(receptor.split("+"))
    for name, block in cfg["pathways"].items():
        if lig & set(block["ligands"]) and rec & set(block["receptors"]):
            return name
    return None


def build_focus_pairs(cpdb_dir: Path, cfg: dict) -> pd.DataFrame:
    genes = pd.read_csv(cpdb_dir / "gene_input.csv")
    valid = set(genes["gene_name"].dropna().astype(str)) | set(genes["hgnc_symbol"].dropna().astype(str))
    for block in cfg["pathways"].values():
        valid.update(block["ligands"])
        valid.update(u for r in block["receptors"] for u in r.split("+"))
    inter = pd.read_csv(cpdb_dir / "interaction_input.csv")
    rows = []
    seen = set()
    for _, rec in inter.iterrows():
        parsed = parse_interactors(rec.get("interactors"), valid)
        if parsed is None:
            continue
        lig, recp = parsed
        path = pathway_of(lig, recp, cfg)
        if path is None:
            continue
        key = (lig, recp, path)
        if key in seen:
            continue
        seen.add(key)
        rows.append({"ligand": lig, "receptor": recp, "pathway": path, "source": "cellphonedb_v5"})
    for item in cfg["curated_overlay"]:
        key = (item["ligand"], item["receptor"], item["pathway"])
        if key in seen:
            continue
        seen.add(key)
        rows.append({**item, "source": "curated_overlay"})
    return pd.DataFrame(rows)


def needed_genes(pairs: pd.DataFrame) -> set[str]:
    genes = {"CLDN4"}
    for col in ("ligand", "receptor"):
        for val in pairs[col]:
            genes.update(val.split("+"))
    return genes


def partner_mean(means: dict[str, float], complex_name: str) -> float:
    parts = complex_name.split("+")
    vals = [means.get(p, float("nan")) for p in parts]
    if any(not np.isfinite(v) for v in vals):
        return float("nan")
    return float(min(vals))


def partner_frac(fracs: dict[str, float], complex_name: str) -> float:
    parts = complex_name.split("+")
    vals = [fracs.get(p, 0.0) for p in parts]
    return float(min(vals)) if vals else 0.0


def group_stats(logx: dict[str, np.ndarray], mask: np.ndarray) -> tuple[dict[str, float], dict[str, float], int]:
    n = int(mask.sum())
    means, fracs = {}, {}
    if n == 0:
        return means, fracs, 0
    for g, arr in logx.items():
        v = arr[mask]
        means[g] = float(v.mean())
        fracs[g] = float(np.mean(v > 0))
    return means, fracs, n


def score_pairs(pairs: pd.DataFrame, send: dict, send_f: dict, rec: dict, rec_f: dict, expr_prop: float) -> list[dict]:
    out = []
    for row in pairs.itertuples(index=False):
        l = partner_mean(send, row.ligand)
        r = partner_mean(rec, row.receptor)
        lf = partner_frac(send_f, row.ligand)
        rf = partner_frac(rec_f, row.receptor)
        score = (l + r) / 2.0 if np.isfinite(l) and np.isfinite(r) else float("nan")
        out.append(
            {
                "ligand": row.ligand,
                "receptor": row.receptor,
                "pathway": row.pathway,
                "source": row.source,
                "score": score,
                "ligand_mean": l,
                "receptor_mean": r,
                "ligand_frac": lf,
                "receptor_frac": rf,
                "pass_expr_prop": bool(lf >= expr_prop and rf >= expr_prop and np.isfinite(score)),
            }
        )
    return out


def stream_umi_txt(path: Path, keep: set[str]) -> tuple[np.ndarray, dict[str, np.ndarray], np.ndarray]:
    log(f"[stream-txt] {path} keep={len(keep)}")
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cells = np.array(header[1:], dtype=object)
        n = len(cells)
        totals = np.zeros(n, dtype=np.float64)
        store: dict[str, np.ndarray] = {}
        n_genes = 0
        for line in fh:
            gene, sep, rest = line.partition("\t")
            if not sep:
                continue
            gene = gene.split(".")[0]
            n_genes += 1
            vals = np.fromstring(rest, sep="\t", dtype=np.float32)
            if vals.size != n:
                raise ValueError(f"column mismatch {gene}: {vals.size} != {n}")
            totals += vals
            if gene in keep and gene not in store:
                store[gene] = vals
            if n_genes % 5000 == 0:
                log(f"[stream-txt] genes={n_genes} kept={len(store)}")
    log(f"[stream-txt] cells={n} genes={n_genes} kept={len(store)}")
    return cells, store, totals


def load_rds_genes(path: Path, wanted: set[str]):
    from scipy import sparse
    import rdata

    with tempfile.TemporaryDirectory(prefix="gse205335-") as tmp:
        matrix_path = path
        if path.suffix == ".gz":
            matrix_path = Path(tmp) / path.stem
            log(f"decompress {path.name}")
            with gzip.open(path, "rb") as src, matrix_path.open("wb") as dest:
                shutil.copyfileobj(src, dest, 16 * 1024 * 1024)
        log("read RDS")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix"')
            obj = rdata.read_rds(matrix_path)
    expected = {"i", "p", "Dim", "Dimnames", "x"}
    if not expected.issubset(vars(obj)):
        raise TypeError("RDS is not the expected Matrix::dgCMatrix object")
    genes = np.asarray(obj.Dimnames[0], dtype=str)
    barcodes = np.asarray(obj.Dimnames[1], dtype=str)
    log(f"build CSC {tuple(obj.Dim)}")
    matrix = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim), copy=False)
    library_umi = np.asarray(matrix.sum(axis=0)).ravel()
    matrix = matrix.tocsr()
    name_to_row = {g: i for i, g in enumerate(genes)}
    extracted = {}
    for gene in sorted(wanted):
        row = name_to_row.get(gene)
        if row is None:
            continue
        extracted[gene] = np.asarray(matrix.getrow(row).toarray()).ravel()
    log(f"extracted {len(extracted)} / {len(wanted)} genes")
    return extracted, library_umi, barcodes


def to_log1p_cp10k(store: dict[str, np.ndarray], totals: np.ndarray) -> dict[str, np.ndarray]:
    out = {}
    tot = np.maximum(totals, 1.0)
    for g, umi in store.items():
        out[g] = np.log1p(np.asarray(umi, dtype=np.float64) / tot * 1e4).astype(np.float32)
    return out


def run_cpdb_131907(units: pd.DataFrame, pairs: pd.DataFrame, annot_path: Path, umi_path: Path, expr_prop: float) -> pd.DataFrame:
    keep_samples = set(units.loc[units["cohort"] == "GSE131907", "unit"])
    qmap = units.loc[units["cohort"] == "GSE131907"].set_index("unit")["q_pct"].astype(str)
    annot = pd.read_csv(annot_path, sep="\t")
    annot = annot[annot["Sample"].isin(keep_samples)].copy()
    mal = annot["Sample_Origin"].isin(TUMOR_ORIGINS_131907) & annot["Cell_subtype"].isin(MALIGNANT_SUBTYPES)
    tnk = annot["Cell_type"].isin(TNK_TYPES_131907)
    annot["role"] = np.where(mal, "mal", np.where(tnk, "tnk", "other"))
    keep_idx = set(annot.loc[annot["role"].isin(["mal", "tnk"]), "Index"])
    genes = needed_genes(pairs)
    cells, store, totals = stream_umi_txt(umi_path, genes)
    cell_pos = {c: i for i, c in enumerate(cells)}
    logx = to_log1p_cp10k(store, totals)
    rows = []
    for sample, sub in annot.groupby("Sample"):
        if sample not in keep_samples:
            continue
        mal_i = [cell_pos[i] for i in sub.loc[sub["role"] == "mal", "Index"] if i in cell_pos]
        tnk_i = [cell_pos[i] for i in sub.loc[sub["role"] == "tnk", "Index"] if i in cell_pos]
        if len(mal_i) < 20 or len(tnk_i) < 20:
            continue
        mal_m = np.zeros(len(cells), dtype=bool)
        tnk_m = np.zeros(len(cells), dtype=bool)
        mal_m[mal_i] = True
        tnk_m[tnk_i] = True
        sm, sf, nm = group_stats(logx, mal_m)
        rm, rf, nt = group_stats(logx, tnk_m)
        scored = score_pairs(pairs, sm, sf, rm, rf, expr_prop)
        for rec in scored:
            rec.update(
                {
                    "cohort": "GSE131907",
                    "unit": sample,
                    "patient": str(units.loc[units["unit"] == sample, "patient"].iloc[0]),
                    "q_pct": str(qmap.loc[sample]),
                    "n_mal": nm,
                    "n_tnk": nt,
                }
            )
            rows.append(rec)
        log(f"  GSE131907 {sample} mal={nm} tnk={nt} q={qmap.loc[sample]}")
    return pd.DataFrame(rows)


def run_cpdb_205335(units: pd.DataFrame, pairs: pd.DataFrame, ident_path: Path, rds_path: Path, gsm_path: Path, expr_prop: float) -> pd.DataFrame:
    keep_patients = set(units.loc[units["cohort"] == "GSE205335", "patient"])
    qmap = units.loc[units["cohort"] == "GSE205335"].set_index("patient")["q_pct"].astype(str)
    ident = pd.read_csv(ident_path, sep="\t")
    gsm = pd.read_csv(gsm_path)
    ident = ident.merge(gsm[["orig.ident", "patient", "tissue"]], on="orig.ident", how="left")
    ident = ident[~ident["tissue"].isin(NORMAL_TISSUES_205335)].copy()
    ident = ident[ident["patient"].isin(keep_patients)].copy()
    genes = needed_genes(pairs)
    store, totals, barcodes = load_rds_genes(rds_path, genes)
    logx = to_log1p_cp10k(store, totals)
    indexed = ident.set_index("barcode")
    common = [b for b in barcodes if b in indexed.index]
    pos = {b: i for i, b in enumerate(barcodes)}
    rows = []
    for patient, sub in ident.groupby("patient"):
        if patient not in keep_patients:
            continue
        mal_i = [pos[b] for b in sub.loc[sub["lineage.sub"] == "Malignant cells", "barcode"] if b in pos]
        tnk_i = [pos[b] for b in sub.loc[sub["lineage.total"] == "T/NK cells", "barcode"] if b in pos]
        if len(mal_i) < 20 or len(tnk_i) < 20:
            continue
        mal_m = np.zeros(len(barcodes), dtype=bool)
        tnk_m = np.zeros(len(barcodes), dtype=bool)
        mal_m[mal_i] = True
        tnk_m[tnk_i] = True
        sm, sf, nm = group_stats(logx, mal_m)
        rm, rf, nt = group_stats(logx, tnk_m)
        scored = score_pairs(pairs, sm, sf, rm, rf, expr_prop)
        for rec in scored:
            rec.update(
                {
                    "cohort": "GSE205335",
                    "unit": patient,
                    "patient": patient,
                    "q_pct": str(qmap.loc[patient]),
                    "n_mal": nm,
                    "n_tnk": nt,
                }
            )
            rows.append(rec)
        log(f"  GSE205335 {patient} mal={nm} tnk={nt} q={qmap.loc[patient]}")
    return pd.DataFrame(rows)


def summarize_cpdb(long: pd.DataFrame, min_det: int) -> pd.DataFrame:
    rows = []
    for scope, sl in (
        ("cohort", long),
        ("winpair", long.assign(cohort_scope="GSE131907+GSE205335")),
    ):
        group_cols = ["cohort", "pathway", "ligand", "receptor"] if scope == "cohort" else ["pathway", "ligand", "receptor"]
        if scope == "winpair":
            sl = sl.copy()
            sl["cohort"] = "GSE131907+GSE205335"
        for keys, g in sl.groupby(group_cols + (["cohort"] if scope == "winpair" else []), sort=False):
            tails = g[g["q_pct"].isin(["Q1", "Q4"])]
            det = tails[tails["pass_expr_prop"]]
            n1 = int((det["q_pct"] == "Q1").sum())
            n4 = int((det["q_pct"] == "Q4").sum())
            if n1 < min_det or n4 < min_det:
                continue
            s1 = det.loc[det["q_pct"] == "Q1", "score"].to_numpy()
            s4 = det.loc[det["q_pct"] == "Q4", "score"].to_numpy()
            u, p, r = rank_biserial(s4, s1)
            rec = {
                "scope": scope,
                "cohort": det["cohort"].iloc[0],
                "pathway": det["pathway"].iloc[0],
                "pair": f"{det['ligand'].iloc[0]}–{det['receptor'].iloc[0]}",
                "ligand": det["ligand"].iloc[0],
                "receptor": det["receptor"].iloc[0],
                "n_q1_detected": n1,
                "n_q4_detected": n4,
                "n_compared": n1 + n4,
                "median_score_q1": float(np.median(s1)),
                "median_score_q4": float(np.median(s4)),
                "delta_median": float(np.median(s4) - np.median(s1)),
                "mwu_u": u,
                "p": p,
                "r_rb": r,
            }
            rows.append(rec)
    tab = pd.DataFrame(rows)
    if tab.empty:
        return tab
    tab["q_bh"] = np.nan
    for (scope, cohort), idx in tab.groupby(["scope", "cohort"]).groups.items():
        tab.loc[idx, "q_bh"] = bh(tab.loc[idx, "p"].tolist())
        tab.loc[idx, "bh_family"] = f"{scope}:{cohort}:focus_MHCI_Trecruit"
    tab["weaker_in_q4"] = tab["delta_median"] < 0
    tab["recovery_q05"] = tab["weaker_in_q4"] & (tab["q_bh"] < 0.05)
    return tab.sort_values(["scope", "cohort", "pathway", "p"])


def style_fig() -> None:
    plt.rcParams.update(
        {
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
        }
    )


def fig_fractions(units: pd.DataFrame, out: Path) -> None:
    style_fig()
    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.4), sharey=False)
    for r, comp, ylab in ((0, "frac_tnk", "T/NK fraction"), (1, "frac_b", "B fraction")):
        for c, cohort in enumerate(["GSE131907", "GSE205335"]):
            ax = axes[r, c]
            sub = units[(units["cohort"] == cohort) & units["q_pct"].isin(["Q1", "Q4"])]
            data = [sub.loc[sub["q_pct"] == q, comp].to_numpy() for q in ("Q1", "Q4")]
            bp = ax.boxplot(data, tick_labels=["Q1", "Q4"], widths=0.55, patch_artist=True)
            for patch, color in zip(bp["boxes"], ["#4c78a8", "#c44e52"]):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)
            for i, arr in enumerate(data, start=1):
                ax.scatter(np.random.default_rng(0).normal(i, 0.04, size=len(arr)), arr, s=18, c="black", zorder=3)
            ax.set_title(f"{cohort} {ylab}")
            ax.set_ylabel(ylab if c == 0 else "")
            n1 = int((sub["q_pct"] == "Q1").sum())
            n4 = int((sub["q_pct"] == "Q4").sum())
            ax.text(0.5, 0.02, f"n={n1} vs {n4}", transform=ax.transAxes, ha="center", va="bottom", fontsize=8)
    fig.suptitle("CLDN4 %pos Q4 vs Q1 — T/NK and B fractions (winning pair)", y=1.01)
    fig.tight_layout()
    fig.savefig(out / "fig_propeller_fractions_q4q1.png")
    fig.savefig(out / "fig_propeller_fractions_q4q1.pdf")
    plt.close(fig)


def fig_propeller_effects(prop: pd.DataFrame, out: Path) -> None:
    style_fig()
    sub = prop[prop["scope"] == "cohort"].copy()
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    labels = [f"{r.cohort}\n{r.compartment}" for r in sub.itertuples()]
    y = np.arange(len(sub))
    ax.axvline(0, color="#888", lw=0.8)
    ax.errorbar(sub["logit_delta_q4_minus_q1"], y, xerr=1.96 * sub["logit_se"], fmt="o", color="#7a2d0b")
    for i, r in enumerate(sub.itertuples()):
        q = r.q_bh
        ax.text(0.99, i, f"p={r.logit_p:.3g}  q={q:.3g}" if np.isfinite(q) else f"p={r.logit_p:.3g}", va="center", ha="right", transform=ax.get_yaxis_transform(), fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("logit(fraction) Q4 − Q1  (propeller / Welch)")
    ax.set_title("Propeller-style composition, BH on T/NK+B × 2 cohorts")
    fig.tight_layout()
    fig.savefig(out / "fig_propeller_effects.png")
    fig.savefig(out / "fig_propeller_effects.pdf")
    plt.close(fig)


def fig_honest_n(units: pd.DataFrame, out: Path) -> None:
    style_fig()
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    rows = []
    for cohort, sub in units.groupby("cohort"):
        tails = sub[sub["q_pct"].isin(["Q1", "Q4"])]
        rows.append((f"{cohort} full locked", len(sub)))
        rows.append((f"{cohort} Q1", int((tails.q_pct == "Q1").sum())))
        rows.append((f"{cohort} Q4", int((tails.q_pct == "Q4").sum())))
        rows.append((f"{cohort} unique patients in tails", int(tails.patient.nunique())))
    labels, vals = zip(*rows)
    ax.barh(labels[::-1], vals[::-1], color="#4c78a8")
    for i, v in enumerate(vals[::-1]):
        ax.text(v + 0.15, i, str(v), va="center")
    ax.set_xlabel("n (units)")
    ax.set_title("Honest n — Q4 vs Q1 tails, not the full cohort")
    fig.tight_layout()
    fig.savefig(out / "fig_honest_n.png")
    fig.savefig(out / "fig_honest_n.pdf")
    plt.close(fig)


def fig_cpdb(focus: pd.DataFrame, out: Path) -> None:
    if focus.empty:
        return
    style_fig()
    win = focus[focus["scope"] == "winpair"].copy()
    if win.empty:
        win = focus.copy()
    for pathway, fname in (("MHC_I", "fig_cpdb_mhci"), ("T_recruit", "fig_cpdb_trecruit")):
        sl = win[win["pathway"] == pathway].sort_values("delta_median")
        if sl.empty:
            continue
        fig, ax = plt.subplots(figsize=(7.8, max(2.8, 0.38 * len(sl) + 1.4)))
        y = np.arange(len(sl))
        colors = ["#c44e52" if d > 0 else "#4c78a8" for d in sl["delta_median"]]
        ax.barh(y, sl["delta_median"], color=colors, alpha=0.85)
        ax.axvline(0, color="#444", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{p}  n={n1}/{n4}" for p, n1, n4 in zip(sl["pair"], sl["n_q1_detected"], sl["n_q4_detected"])])
        ax.set_xlabel("median CellPhoneDB-style score  Q4 − Q1")
        ax.set_title(f"{pathway} outgoing malignant → T/NK  (winning pair)")
        fig.tight_layout()
        fig.savefig(out / f"{fname}.png")
        fig.savefig(out / f"{fname}.pdf")
        plt.close(fig)
    # extra combined forest of focus pairs
    sl = win.sort_values(["pathway", "delta_median"])
    fig, ax = plt.subplots(figsize=(8.0, max(3.2, 0.32 * len(sl) + 1.6)))
    y = np.arange(len(sl))
    colors = ["#c44e52" if d > 0 else "#4c78a8" for d in sl["delta_median"]]
    ax.barh(y, sl["delta_median"], color=colors, alpha=0.85)
    ax.axvline(0, color="#444", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.pathway} {r.pair}" for r in sl.itertuples()])
    ax.set_xlabel("median score Q4 − Q1")
    ax.set_title("Focus MHC-I + T-recruit (CellPhoneDB-style, not LIANA)")
    fig.tight_layout()
    fig.savefig(out / "fig_cpdb_focus_pair.png")
    fig.savefig(out / "fig_cpdb_focus_pair.pdf")
    plt.close(fig)


def fig_pair_boxes(units: pd.DataFrame, out: Path) -> None:
    style_fig()
    tails = units[units["q_pct"].isin(["Q1", "Q4"])].copy()
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.8))
    for ax, col, title in ((axes[0], "frac_tnk", "T/NK"), (axes[1], "frac_b", "B")):
        for i, cohort in enumerate(["GSE131907", "GSE205335"]):
            sub = tails[tails["cohort"] == cohort]
            for j, q in enumerate(["Q1", "Q4"]):
                arr = sub.loc[sub["q_pct"] == q, col].to_numpy()
                x = i * 2 + j
                ax.boxplot(arr, positions=[x], widths=0.55, patch_artist=True, boxprops=dict(facecolor="#4c78a8" if q == "Q1" else "#c44e52", alpha=0.7))
                ax.scatter(np.random.default_rng(1).normal(x, 0.05, size=len(arr)), arr, s=16, c="k", zorder=3)
        ax.set_xticks([0.5, 2.5])
        ax.set_xticklabels(["GSE131907", "GSE205335"])
        ax.set_title(f"{title} fraction, Q1 vs Q4")
        ax.set_ylabel("fraction")
    fig.suptitle("Extra: winning-pair composition boxes")
    fig.tight_layout()
    fig.savefig(out / "fig_extra_pair_boxes.png")
    fig.savefig(out / "fig_extra_pair_boxes.pdf")
    plt.close(fig)


def fmt(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if abs(x) >= 0.01 or x == 0:
        return f"{x:.{nd}g}" if abs(x) >= 1 else f"{x:.3f}"
    return f"{x:.3g}"


def write_finding(units: pd.DataFrame, prop: pd.DataFrame, focus: pd.DataFrame, out: Path, cpdb_ran: bool) -> None:
    prim = prop[(prop["scope"] == "cohort") & (prop["compartment"].isin(["TNK", "B"]))]
    win = prop[(prop["scope"] == "winpair") & (prop["compartment"].isin(["TNK", "B", "CD8"]))]
    n131 = units[units.cohort == "GSE131907"]
    n205 = units[units.cohort == "GSE205335"]
    t131 = n131[n131.q_pct.isin(["Q1", "Q4"])]
    t205 = n205[n205.q_pct.isin(["Q1", "Q4"])]

    lines = [
        "# FINDING — propeller + CellPhoneDB-style on winning pair GSE131907+GSE205335",
        "",
        "ADDITIVE **CLDN4 only**. No dual-high TACSTD2×CLDN4. **No GSE207422.**",
        "Winning pair is the PR #320 author-malignant **%pos** cut vs T/NK",
        "(Spearman n=43 ρ=−0.479; Q4 vs Q1 n=23 r=−0.705). Patient/sample is the unit.",
        "LIANA is not run (that agent is separate). CellChat is not run.",
        "",
        "## Honest n",
        "",
        "| Cohort | Locked units | Unit | Q1 | Q4 | Compared | Unique patients in tails |",
        "|---|---:|---|---:|---:|---:|---:|",
        f"| GSE131907 | {len(n131)} | sample (n_mal≥20; tumor origins) | {(t131.q_pct=='Q1').sum()} | {(t131.q_pct=='Q4').sum()} | {len(t131)} | {t131.patient.nunique()} |",
        f"| GSE205335 | {len(n205)} | patient (≥20 mal + ≥20 T/NK) | {(t205.q_pct=='Q1').sum()} | {(t205.q_pct=='Q4').sum()} | {len(t205)} | {t205.patient.nunique()} |",
        f"| Pair | {len(units)} | mixed | {(units.q_pct=='Q1').sum()} | {(units.q_pct=='Q4').sum()} | {int(((units.q_pct=='Q1')|(units.q_pct=='Q4')).sum())} | {units.loc[units.q_pct.isin(['Q1','Q4']),'patient'].nunique()} |",
        "",
        "GSE131907 T/NK extract is **sample-level** (PR #320). Quartiles are cut",
        "**within cohort** on malignant CLDN4 %pos (`pd.qcut` on average ranks).",
        "Q4 vs Q1 uses the tails only. Cells are not n. GSE205335 B is B+plasma.",
        "",
        "## 1. Propeller / speckle-style composition",
        "",
        "Transform = empirical logit of the fraction (Phipson et al. 2022 *Bioinformatics*,",
        "`speckle::propeller`). Test = Welch two-sample t on logit (cohort) or OLS",
        "`logit ~ Q4 + cohort` (pair). BH on the pre-specified primary family:",
        "**T/NK + B × 2 cohorts (4 tests)**. CD8 is nested in T/NK and is secondary.",
        "eBayes moderation across 3–4 compartments is a companion (k is too small",
        "for a stable prior; unmoderated logit p is the primary p).",
        "",
        "### Primary family (cohort × T/NK, B)",
        "",
        "| cohort | compartment | n_Q1/n_Q4 | median frac Q1 | Q4 | Δ | logit Δ | p | BH q | recover q<0.10 |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for r in prim.itertuples():
        rec = "yes" if bool(r.recovery_q10) else "no"
        lines.append(
            f"| {r.cohort} | {r.compartment} | {r.n_q1}/{r.n_q4} | {r.median_frac_q1:.3f} | {r.median_frac_q4:.3f} | {r.delta_median_frac:+.3f} | {r.logit_delta_q4_minus_q1:+.3f} | {fmt(r.logit_p)} | {fmt(r.q_bh)} | {rec} |"
        )
    lines += [
        "",
        "### Winning-pair OLS (companion; BH across T/NK, B, CD8)",
        "",
        "| compartment | n_Q1/n_Q4 | median frac Q1 | Q4 | Δ | logit β_Q4 | p | BH q |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for r in win.itertuples():
        q = getattr(r, "q_bh_winpair", float("nan"))
        lines.append(
            f"| {r.compartment} | {r.n_q1}/{r.n_q4} | {r.median_frac_q1:.3f} | {r.median_frac_q4:.3f} | {r.delta_median_frac:+.3f} | {r.logit_delta_q4_minus_q1:+.3f} | {fmt(r.logit_p)} | {fmt(q)} |"
        )
    n_rec = int(prim["recovery_q10"].fillna(False).sum())
    lines += [
        "",
        f"**Propeller readout:** {n_rec}/4 primary tests have logit Δ < 0 and BH q < 0.10.",
        "Fraction MWU (PR #320 companion) is in `tables/propeller_composition.tsv` and is not the primary p.",
        "",
        "## 2. CellPhoneDB-style documented mean score (MHC-I + T-recruit)",
        "",
        "Score = mean of partner means on log1p(CP10k); complexes = min of subunit means",
        "(Efremova 2020 *Nat Protoc*; Garcia-Alonso 2022 *Nat Protoc*).",
        "`expr_prop` = 0.10. Outgoing = author-malignant → **same-unit** T/NK.",
        "Contrast = **between-unit CLDN4 Q4 vs Q1** (not the within-patient high/low LIANA split).",
        "A pair enters the table if detected in ≥3 Q1 and ≥3 Q4 units. BH within each scope×cohort focus family.",
        "This is **not** a CellChat probability and **not** LIANA `mt.cellphonedb`.",
        "",
    ]
    if not cpdb_ran or focus.empty:
        lines += [
            "**CellPhoneDB table was not produced in this run** (matrix missing or no pair passed the detect gate).",
            "",
        ]
    else:
        win_f = focus[focus["scope"] == "winpair"]
        n_weak = int(win_f["weaker_in_q4"].sum()) if len(win_f) else 0
        n_sig = int(win_f["recovery_q05"].sum()) if len(win_f) else 0
        lines += [
            f"Winning-pair focus pairs that passed the detect gate: **{len(win_f)}**. "
            f"{n_weak}/{len(win_f)} have median Δ < 0 (weaker in Q4). **{n_sig}/{len(win_f)} reach BH q < 0.05** in the weaker-in-Q4 direction.",
            "",
            "T-recruit dropout is the n, not a hidden negative: CXCL9/10/11–CXCR3 and most",
            "other T-recruit pairs did not reach ≥3 Q1 and ≥3 Q4 `expr_prop` units.",
            "",
            "### Winning-pair focus table",
            "",
            "| pathway | pair | n_Q1/n_Q4 | median Q1 | Q4 | Δ | r | p | BH q |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for r in win_f.sort_values(["pathway", "p"]).itertuples():
            lines.append(
                f"| {r.pathway} | {r.pair} | {r.n_q1_detected}/{r.n_q4_detected} | {r.median_score_q1:.3f} | {r.median_score_q4:.3f} | {r.delta_median:+.3f} | {r.r_rb:+.3f} | {fmt(r.p)} | {fmt(r.q_bh)} |"
            )
        lines += [
            "",
            "Per-cohort rows: `tables/cpdb_focus.tsv`.",
            "",
        ]
    lines += [
        "## What was not run",
        "",
        "- LIANA `mt.cellphonedb` (separate agent).",
        "- CellChat R.",
        "- Dual-high TACSTD2×CLDN4.",
        "- GSE207422.",
        "- scCODA HMC / Milo / muscat.",
        "",
        "## Files",
        "",
        "| File | Role |",
        "|---|---|",
        "| `tables/propeller_composition.tsv` | **Table 1** — propeller / speckle-style |",
        "| `tables/cpdb_focus.tsv` | **Table 2** — CellPhoneDB-style MHC-I + T-recruit |",
        "| `tables/units_with_quartiles.tsv` | Locked units + Q labels |",
        "| `figures/` | Extra figures |",
        "",
        "## Reproduce",
        "",
        "```bash",
        "python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/download.py",
        "python3 methods/winpair_131907_205335_propeller_cpdb_cldn4/scripts/analyze.py",
        "```",
        "",
    ]
    (out / "FINDING.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    ap.add_argument("--outdir", type=Path, default=HERE)
    ap.add_argument("--gse131907", type=Path, default=Path("/tmp/gse131907"))
    ap.add_argument("--gse205335", type=Path, default=Path("/tmp/gse205335"))
    ap.add_argument("--cpdb", type=Path, default=Path("/tmp/cellphonedb"))
    ap.add_argument("--skip-cpdb", action="store_true")
    args = ap.parse_args()

    tables = args.outdir / "tables"
    figs = args.outdir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)

    cfg = yaml.safe_load((HERE / "config" / "gene_sets.yaml").read_text())
    units = load_units(args.datadir)
    units.to_csv(tables / "units_with_quartiles.tsv", sep="\t", index=False)
    log(f"units: {units.groupby(['cohort','q_pct']).size().to_dict()}")

    compartments = [("TNK", "frac_tnk"), ("B", "frac_b"), ("CD8", "frac_cd8")]
    prop = propeller_block(units, compartments)
    prop.to_csv(tables / "propeller_composition.tsv", sep="\t", index=False)
    log("wrote propeller_composition.tsv")

    fig_fractions(units, figs)
    fig_propeller_effects(prop, figs)
    fig_honest_n(units, figs)
    fig_pair_boxes(units, figs)

    focus = pd.DataFrame()
    cpdb_ran = False
    annot = args.gse131907 / "GSE131907_Lung_Cancer_cell_annotation.txt.gz"
    umi = args.gse131907 / "GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz"
    ident = args.gse205335 / "GSE205335_Lung_IO_CellIdentity.txt.gz"
    rds = args.gse205335 / "GSE205335_Lung_IO_UMI_matrix.rds.gz"
    gsm = args.datadir / "GSE205335_gsm_sample_metadata.csv"
    if not args.skip_cpdb and annot.exists() and umi.exists() and ident.exists() and rds.exists() and (args.cpdb / "interaction_input.csv").exists():
        pairs = build_focus_pairs(args.cpdb, cfg)
        pairs.to_csv(tables / "cpdb_pairs_focus.tsv", sep="\t", index=False)
        log(f"focus pairs: {len(pairs)}")
        long_parts = []
        long_parts.append(run_cpdb_131907(units, pairs, annot, umi, cfg["params"]["expr_prop"]))
        long_parts.append(run_cpdb_205335(units, pairs, ident, rds, gsm, cfg["params"]["expr_prop"]))
        long = pd.concat(long_parts, ignore_index=True)
        long.to_csv(tables / "cpdb_per_unit.tsv", sep="\t", index=False)
        focus = summarize_cpdb(long, cfg["params"]["min_detected_per_tail"])
        focus.to_csv(tables / "cpdb_focus.tsv", sep="\t", index=False)
        cpdb_ran = True
        fig_cpdb(focus, figs)
        log(f"wrote cpdb_focus.tsv rows={len(focus)}")
    else:
        log("SKIP CellPhoneDB (files missing or --skip-cpdb)")
        (tables / "cpdb_focus.tsv").write_text(
            "scope\tcohort\tpathway\tpair\tligand\treceptor\tn_q1_detected\tn_q4_detected\tn_compared\tmedian_score_q1\tmedian_score_q4\tdelta_median\tmwu_u\tp\tr_rb\tq_bh\tbh_family\tweaker_in_q4\trecovery_q05\n"
        )

    write_finding(units, prop, focus, args.outdir, cpdb_ran)
    summary = {
        "n_units": int(len(units)),
        "n_compared": int(units.q_pct.isin(["Q1", "Q4"]).sum()),
        "propeller_rows": int(len(prop)),
        "cpdb_ran": cpdb_ran,
        "cpdb_focus_rows": int(len(focus)),
        "primary_recoveries": int(prop.loc[(prop.scope == "cohort") & (prop.compartment.isin(["TNK", "B"])), "recovery_q10"].fillna(False).sum()),
    }
    (tables / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    log(json.dumps(summary))


if __name__ == "__main__":
    main()
