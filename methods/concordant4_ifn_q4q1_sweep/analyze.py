#!/usr/bin/env python3
"""Q4 vs Q1 IFN / MHC / chemokine sweep on concordant-4 malignant pseudobulks.

The unit is the patient, donor, or sample. Quartiles are the locked
within-cohort malignant CLDN4 %pos labels. Cell-level p-values are not
computed. CLDN4 is removed from every gene set.

The locked comparison is the cohort-adjusted OLS t ranking on the
zero-filled Q1/Q4 matrix (Hallmark IFN-gamma NES about -3.85 in the
prior fgsea run). This script repeats that ranking and then sweeps
pre-specified rank statistics and public gene sets. The GSE205335 drop
is a second scope, not a cohort chosen after seeing the ranks.

Headline |NES| is the most negative NES among interferon sets with
15-500 genes in the ranked universe. Headline |logFC| is reported two
ways, both pre-specified: the OLS beta of the mean log2(TMM-CPM+1) of
the set, and the most negative non-sparse gene inside the interferon
sets. Leading-edge mean logFC is stored as a description of the
enrichment, not as a pre-specified effect.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import zlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import special
from scipy.stats import norm, t as student_t

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"

SEED = 20260921
COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
REF = "GSE123902"
MIN_ARM = 3
MIN_SET = 10
HEADLINE_MIN = 15
MAX_SET = 500
NPERM_SWEEP = 2000
NPERM_REFINE = 20000
NPERM_BOOT = 1000
N_BOOT_NES = 200
N_BOOT_FC = 2000

# Lymphocyte / myeloid identity genes removed in the lineage-depleted
# Hallmark sets. IFN-stimulated epithelial genes (STAT, ISG, HLA, CD274,
# CXCL, CCL) stay. Frozen before the sweep.
LINEAGE = frozenset(
    {
        "PTPRC", "CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B",
        "CD19", "MS4A1", "CD79A", "CD79B", "NKG7", "KLRD1", "KLRK1",
        "KLRB1", "NCAM1", "GNLY", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1",
        "LCK", "ZAP70", "LCP2", "SKAP1", "FYB1", "SELL", "CCR7", "CXCR3",
        "CXCR6", "CD14", "FCGR3A", "FCGR1A", "CSF1R", "CSF2RB", "CD68",
        "LYZ", "CD69", "CD38", "SLAMF7", "CD40", "CD86", "IL2RB", "IL7R",
        "IL2RG", "FPR1", "ITGB7", "XCL1", "XCL2", "TRAC", "TRBC1", "TRBC2",
    }
)

# Epithelial IFN / APM program. No chemokine ligands and no lineage genes.
CORE_ISG = [
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "ISG15", "MX1", "MX2",
    "OAS1", "OAS2", "OAS3", "OASL", "IFIT1", "IFIT2", "IFIT3", "IFI6",
    "IFI27", "IFI35", "IFI44", "IFI44L", "RSAD2", "DDX58", "IFIH1",
    "EIF2AK2", "GBP1", "GBP2", "GBP4", "BST2", "HERC5", "HERC6", "USP18",
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "B2M", "TAP1", "TAP2",
    "TAPBP", "PSMB8", "PSMB9", "PSMB10", "NLRC5",
]
CXCR3_CHEMOKINES = ["CXCL9", "CXCL10", "CXCL11", "CXCL13", "CCL5"]
INFLAMMATORY_CHEMOKINES = [
    "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8", "CXCL9",
    "CXCL10", "CXCL11", "CXCL13", "CXCL16", "CCL2", "CCL3", "CCL4",
    "CCL5", "CCL7", "CCL8", "CCL11", "CCL20", "CCL22", "CX3CL1", "XCL1",
    "XCL2",
]
# Gene-level logFC panel. Same biology, declared before the fit.
CANONICAL = CORE_ISG + CXCR3_CHEMOKINES + ["CCL2", "CD274", "IDO1", "CIITA", "GBP5"]

ALIASES = {
    "WARS": "WARS1",
    "MARCH1": "MARCHF1",
    "MARCH2": "MARCHF2",
    "MARCH3": "MARCHF3",
    "MARCH4": "MARCHF4",
    "MARCH5": "MARCHF5",
    "MARCH6": "MARCHF6",
    "MARCH7": "MARCHF7",
    "MARCH8": "MARCHF8",
    "MARCH9": "MARCHF9",
    "MARCH10": "MARCHF10",
    "MARCH11": "MARCHF11",
}

POOLED_RANKS = [
    "ols_t",
    "ols_logFC",
    "ols_signed_neglog10p",
    "moderated_t",
    "snr_residual",
    "wilcoxon_z_residual",
]
COMBINED_RANKS = ["stouffer_z", "mean_logFC", "ivw_logFC", "max_logFC"]


def child_seed(*parts) -> int:
    payload = "||".join(map(str, parts)).encode()
    return int((SEED + zlib.adler32(payload)) % (2**32 - 1))


def read_counts(path: Path):
    genes = []
    rows = []
    samples = None
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        for row in reader:
            if not row or not row[0]:
                continue
            genes.append(row[0].strip().upper())
            vals = []
            for x in row[1:]:
                vals.append(float(x) if x else 0.0)
            if len(vals) != len(samples):
                raise SystemExit(f"{path.name}: ragged row {genes[-1]}")
            rows.append(vals)
    mat = np.asarray(rows, dtype=np.float64)
    if len(genes) != len(set(genes)):
        acc = {}
        order = []
        for i, g in enumerate(genes):
            if g not in acc:
                acc[g] = mat[i].copy()
                order.append(g)
            else:
                acc[g] += mat[i]
        genes = order
        mat = np.vstack([acc[g] for g in genes])
    return genes, samples, mat


def read_units(path: Path):
    with path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    out = []
    for row in rows:
        if row["cohort"] not in COHORTS:
            continue
        flag = row["in_count_matrix"].strip().lower()
        if flag not in {"true", "t", "1"}:
            continue
        out.append(
            {
                "patient": row["patient"],
                "cohort": row["cohort"],
                "quartile": row["quartile"],
                "cldn4_pct": float(row["cldn4_pct"]),
            }
        )
    out.sort(key=lambda r: (COHORTS.index(r["cohort"]), r["patient"]))
    return out


def bind_counts(cohort_mats, sample_order):
    universe = set()
    for genes, _samples, _mat in cohort_mats.values():
        universe.update(genes)
    genes = sorted(universe)
    index = {g: i for i, g in enumerate(genes)}
    mat = np.zeros((len(genes), len(sample_order)), dtype=np.float64)
    sample_pos = {s: i for i, s in enumerate(sample_order)}
    for genes_c, samples_c, mat_c in cohort_mats.values():
        cols = [sample_pos[s] for s in samples_c if s in sample_pos]
        src = [i for i, s in enumerate(samples_c) if s in sample_pos]
        if not cols:
            continue
        rows = [index[g] for g in genes_c]
        mat[np.ix_(rows, cols)] = mat_c[:, src]
    return genes, mat


def filter_genes(genes, mat, min_count=10, min_samples=3):
    keep = (mat >= min_count).sum(axis=1) >= min_samples
    return [g for g, k in zip(genes, keep) if k], mat[keep]


def tmm_norm_factors(counts: np.ndarray) -> np.ndarray:
    lib = counts.sum(axis=0).astype(np.float64)
    lib[lib == 0] = np.nan
    rel = counts / lib
    f75 = np.quantile(rel, 0.75, axis=0)
    ref = int(np.nanargmin(np.abs(f75 - np.nanmean(f75))))
    ref_c = counts[:, ref]
    ref_lib = float(lib[ref])
    factors = np.ones(counts.shape[1], dtype=np.float64)
    for col in range(counts.shape[1]):
        obs = counts[:, col]
        obs_lib = float(lib[col])
        if not math.isfinite(obs_lib) or obs_lib == 0:
            factors[col] = 1.0
            continue
        keep = (obs > 0) & (ref_c > 0)
        if int(keep.sum()) < 50:
            factors[col] = 1.0
            continue
        obs_k = obs[keep]
        ref_k = ref_c[keep]
        m = np.log2((obs_k / obs_lib) / (ref_k / ref_lib))
        a = 0.5 * np.log2((obs_k / obs_lib) * (ref_k / ref_lib))
        w = (obs_lib - obs_k) / (obs_lib * obs_k) + (ref_lib - ref_k) / (ref_lib * ref_k)
        ok = np.isfinite(m) & np.isfinite(a) & np.isfinite(w) & (w > 0)
        m, a, w = m[ok], a[ok], w[ok]
        if m.size < 50:
            factors[col] = 1.0
            continue
        trim = (
            (m >= np.quantile(m, 0.30))
            & (m <= np.quantile(m, 0.70))
            & (a >= np.quantile(a, 0.05))
            & (a <= np.quantile(a, 0.95))
        )
        if int(trim.sum()) < 20:
            factors[col] = 1.0
            continue
        factors[col] = 2.0 ** np.average(m[trim], weights=1.0 / w[trim])
    return factors / factors.mean()


def log_cpm(counts: np.ndarray) -> np.ndarray:
    fac = tmm_norm_factors(counts)
    lib = counts.sum(axis=0) * fac
    return np.log2(counts / lib * 1e6 + 1.0)


def design_matrix(cohorts, exposure):
    present = [c for c in COHORTS if c in set(cohorts)]
    for c in cohorts:
        if c not in present:
            present.append(c)
    ref = REF if REF in present else present[0]
    others = [c for c in present if c != ref]
    n = len(cohorts)
    x = np.zeros((n, 2 + len(others)), dtype=np.float64)
    x[:, 0] = 1.0
    for j, c in enumerate(others):
        x[:, 1 + j] = np.asarray([1.0 if v == c else 0.0 for v in cohorts])
    x[:, -1] = np.asarray(exposure, dtype=np.float64)
    return x


def fit_exposure(y, cohorts, exposure):
    x = design_matrix(cohorts, exposure)
    if np.linalg.matrix_rank(x) < x.shape[1]:
        return None
    beta, _, _, _ = np.linalg.lstsq(x, y.T, rcond=None)
    resid = y.T - x @ beta
    df = x.shape[0] - x.shape[1]
    if df < 1:
        return None
    sigma2 = np.sum(resid * resid, axis=0) / df
    inv = np.linalg.inv(x.T @ x)
    cjj = float(inv[-1, -1])
    se = np.sqrt(np.maximum(sigma2 * cjj, 0.0))
    est = beta[-1].astype(np.float64)
    tstat = np.zeros_like(est)
    np.divide(est, se, out=tstat, where=se > 0)
    p = 2.0 * student_t.sf(np.abs(tstat), df)
    return {
        "logFC": est,
        "t": tstat,
        "se": se,
        "p": p,
        "df": df,
        "sigma2": sigma2,
        "cjj": cjj,
    }


def moderated_t(fit):
    sigma2 = np.asarray(fit["sigma2"], dtype=np.float64)
    df1 = float(fit["df"])
    ok = np.isfinite(sigma2) & (sigma2 > 0)
    e = np.log(sigma2[ok]) - special.digamma(df1 / 2.0) + math.log(df1 / 2.0)
    mean_e = float(np.mean(e))
    var_e = float(np.mean((e - mean_e) ** 2))
    target = var_e - float(special.polygamma(1, df1 / 2.0))
    if target <= 1e-8:
        d0 = math.inf
        s0 = math.exp(mean_e)
    else:
        lo, hi = 1e-4, 1e8
        for _ in range(80):
            mid = math.sqrt(lo * hi)
            if float(special.polygamma(1, mid / 2.0)) > target:
                lo = mid
            else:
                hi = mid
        d0 = mid
        s0 = math.exp(mean_e - float(special.digamma(d0 / 2.0)) + math.log(d0 / 2.0))
    if math.isinf(d0):
        s2_post = np.full_like(sigma2, s0)
    else:
        s2_post = (d0 * s0 + df1 * sigma2) / (d0 + df1)
    se = np.sqrt(np.maximum(s2_post * fit["cjj"], 0.0))
    out = np.zeros_like(fit["logFC"])
    np.divide(fit["logFC"], se, out=out, where=se > 0)
    return out


def residualize_cohort(y, cohorts):
    x = design_matrix(cohorts, np.zeros(len(cohorts)))[:, :-1]
    beta, _, _, _ = np.linalg.lstsq(x, y.T, rcond=None)
    return (y.T - x @ beta).T


def wilcoxon_z(y, q4_mask):
    ranks = np.empty_like(y)
    for i in range(y.shape[0]):
        ranks[i] = _average_rank(y[i])
    n4 = int(q4_mask.sum())
    n1 = int((~q4_mask).sum())
    n = n4 + n1
    w = ranks[:, q4_mask].sum(axis=1)
    mu = n4 * (n + 1) / 2.0
    corr = np.empty(y.shape[0])
    for i in range(y.shape[0]):
        _, counts = np.unique(y[i], return_counts=True)
        corr[i] = np.sum(counts.astype(np.float64) ** 3 - counts)
    sigma2 = (n1 * n4 / 12.0) * ((n + 1) - corr / (n * (n - 1)))
    sigma2 = np.maximum(sigma2, 1e-12)
    return (w - mu) / np.sqrt(sigma2)


def _average_rank(v):
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(v.size, dtype=np.float64)
    sorted_v = v[order]
    start = 0
    while start < v.size:
        end = start + 1
        while end < v.size and sorted_v[end] == sorted_v[start]:
            end += 1
        # 1-based average rank
        avg = 0.5 * ((start + 1) + end)
        ranks[order[start:end]] = avg
        start = end
    return ranks


def signal_to_noise(y, q4_mask):
    a = y[:, q4_mask]
    b = y[:, ~q4_mask]
    den = a.std(axis=1, ddof=1) + b.std(axis=1, ddof=1)
    diff = a.mean(axis=1) - b.mean(axis=1)
    out = np.zeros(y.shape[0], dtype=np.float64)
    ok = np.isfinite(den) & (den > 0)
    out[ok] = diff[ok] / den[ok]
    return out


def enrichment_from_hits(sorted_abs, hit_idx):
    """Weighted KS at hit positions. hit_idx is (n, m) and each row is sorted."""
    nperm, m = hit_idx.shape
    n_genes = int(sorted_abs.size)
    weights = sorted_abs[hit_idx]
    nr = weights.sum(axis=1, keepdims=True)
    nr_safe = np.where(nr <= 0, 1.0, nr)
    phit = np.cumsum(weights, axis=1) / nr_safe
    pmiss = (hit_idx - np.arange(m)) / float(n_genes - m)
    phit_before = np.empty_like(phit)
    phit_before[:, 0] = 0.0
    if m > 1:
        phit_before[:, 1:] = phit[:, :-1]
    es_before = phit_before - pmiss
    es_after = phit - pmiss
    hi = es_after.max(axis=1)
    lo = es_before.min(axis=1)
    use_hi = np.abs(hi) >= np.abs(lo)
    es = np.where(use_hi, hi, lo)
    es = np.where(nr.ravel() <= 0, 0.0, es)
    return es, use_hi, es_before, es_after


def es_brute(sorted_stats, mask):
    n_genes = sorted_stats.size
    m = int(mask.sum())
    weights = np.abs(sorted_stats)
    nr = float(weights[mask].sum())
    if m == 0 or m == n_genes or nr <= 0:
        return 0.0
    hit_inc = weights / nr
    miss_inc = 1.0 / (n_genes - m)
    running = 0.0
    hi = 0.0
    lo = 0.0
    for i in range(n_genes):
        if mask[i]:
            running += hit_inc[i]
        else:
            running -= miss_inc
        hi = max(hi, running)
        lo = min(lo, running)
    return hi if abs(hi) >= abs(lo) else lo


def order_stats(names, stats):
    stats = np.asarray(stats, dtype=np.float64)
    finite = np.isfinite(stats)
    names = [n for n, ok in zip(names, finite) if ok]
    stats = stats[finite]
    order = sorted(range(len(names)), key=lambda i: (-stats[i], tuple(-ord(c) for c in names[i])))
    names_s = [names[i] for i in order]
    stats_s = stats[np.asarray(order)]
    return names_s, stats_s


def gsea_one(names_s, stats_s, genes, nperm, rng):
    pos = {g: i for i, g in enumerate(names_s)}
    hits = np.array(sorted(pos[g] for g in genes if g in pos), dtype=np.int64)
    m = int(hits.size)
    n_genes = len(names_s)
    if m < MIN_SET or m > MAX_SET or m >= n_genes:
        return None
    sorted_abs = np.abs(stats_s)
    obs, use_hi, es_before, es_after = enrichment_from_hits(sorted_abs, hits.reshape(1, -1))
    obs = float(obs[0])
    use_hi = bool(use_hi[0])
    if use_hi:
        k = int(np.argmax(es_after[0]))
        thr = int(hits[k])
        lead_pos = hits[hits <= thr]
    else:
        k = int(np.argmin(es_before[0]))
        thr = int(hits[k])
        lead_pos = hits[hits >= thr]
    leading = [names_s[i] for i in lead_pos]
    # Permute gene labels. Chunked so the random-key matrix stays modest.
    null = np.empty(nperm, dtype=np.float64)
    done = 0
    while done < nperm:
        chunk = min(1000, nperm - done)
        keys = rng.random((chunk, n_genes))
        idx = np.argpartition(keys, m - 1, axis=1)[:, :m]
        idx.sort(axis=1)
        es, _, _, _ = enrichment_from_hits(sorted_abs, idx)
        null[done : done + chunk] = es
        done += chunk
    sgn = np.sign(obs)
    same = null[np.sign(null) == sgn] if sgn != 0 else null
    nes = float(obs / np.mean(np.abs(same))) if same.size >= 20 else float("nan")
    if obs > 0:
        p = (1.0 + np.sum(null >= obs)) / (nperm + 1.0)
    elif obs < 0:
        p = (1.0 + np.sum(null <= obs)) / (nperm + 1.0)
    else:
        p = 1.0
    return {
        "size": m,
        "es": obs,
        "nes": nes,
        "p": float(p),
        "n_leading": len(leading),
        "leading": leading,
        "n_genes": n_genes,
        "nperm": nperm,
    }


def bh(pvals):
    p = np.asarray(pvals, dtype=np.float64)
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if int(ok.sum()) == 0:
        return out
    idx = np.flatnonzero(ok)
    ranked_order = idx[np.argsort(p[idx])]
    m = ranked_order.size
    adj = p[ranked_order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    out[ranked_order] = adj
    return out


def load_locked_sets():
    raw = json.loads((DATA / "locked_sets.json").read_text())
    sets = {}
    for key, genes in raw.items():
        sets[key] = sorted({g.upper() for g in genes if g.upper() != "CLDN4"})
    alpha = set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"])
    gamma = set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
    sets["HALLMARK_IFN_CORE"] = sorted(alpha & gamma)
    sets["HALLMARK_IFN_ALPHA_NOLINEAGE"] = sorted(alpha - LINEAGE)
    sets["HALLMARK_IFN_GAMMA_NOLINEAGE"] = sorted(gamma - LINEAGE)
    sets["HALLMARK_IFN_CORE_NOLINEAGE"] = sorted((alpha & gamma) - LINEAGE)
    sets["CUSTOM_CORE_ISG"] = sorted({g for g in CORE_ISG if g != "CLDN4"})
    sets["CUSTOM_CORE_ISG_CXCR3"] = sorted(set(CORE_ISG) | set(CXCR3_CHEMOKINES))
    sets["CUSTOM_INFLAMMATORY_CHEMOKINES"] = sorted(set(INFLAMMATORY_CHEMOKINES))
    return sets


def load_public_sets():
    sets = {}
    desc = {}
    with (DATA / "public_sets.gmt").open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            sid, description, genes = parts[0], parts[1], parts[2:]
            cleaned = []
            for g in genes:
                g = g.strip().upper()
                if not g or g == "CLDN4":
                    continue
                if g not in cleaned:
                    cleaned.append(g)
            sets[sid] = cleaned
            desc[sid] = description
    return sets, desc


def family_of(sid: str) -> str:
    if "CHEMOKINE" in sid or sid == "CUSTOM_INFLAMMATORY_CHEMOKINES":
        return "chemokine"
    if sid == "CUSTOM_CORE_ISG_CXCR3":
        return "ifn"
    if "MHC" in sid or "ANTIGEN" in sid or sid.startswith("CUSTOM_MHC") or sid == "CUSTOM_CORE_ISG":
        # CORE_ISG is the epithelial IFN/APM program, counted as IFN.
        if sid == "CUSTOM_CORE_ISG":
            return "ifn"
        return "mhc"
    return "ifn"


def map_genes(genes, universe):
    mapped = []
    n_alias = 0
    seen = set()
    for g in genes:
        target = g if g in universe else ALIASES.get(g, g)
        if target != g and target in universe:
            n_alias += 1
        if target in universe and target not in seen:
            mapped.append(target)
            seen.add(target)
    return mapped, n_alias


def prepare_zerofill(cohort_counts, meta_rows):
    samples = [r["patient"] for r in meta_rows]
    genes, mat = bind_counts(cohort_counts, samples)
    genes, mat = filter_genes(genes, mat)
    return {
        "genes": genes,
        "logcpm": log_cpm(mat),
        "counts": mat,
        "patients": samples,
        "cohorts": [r["cohort"] for r in meta_rows],
        "quartiles": [r["quartile"] for r in meta_rows],
    }


def prepare_shared(cohort_counts, meta_all, meta_q):
    """Within-cohort TMM on every unit, then the genes that pass in each cohort."""
    per = {}
    for cohort in sorted({r["cohort"] for r in meta_all}):
        rows = [r for r in meta_all if r["cohort"] == cohort]
        samples = [r["patient"] for r in rows]
        genes_c, samples_c, mat_c = cohort_counts[cohort]
        col = {s: i for i, s in enumerate(samples_c)}
        take = [col[s] for s in samples]
        sub = mat_c[:, take]
        genes_f, sub_f = filter_genes(genes_c, sub)
        per[cohort] = {
            "genes": genes_f,
            "logcpm": log_cpm(sub_f),
            "patients": samples,
        }
    shared = set(per[next(iter(per))]["genes"])
    for block in per.values():
        shared &= set(block["genes"])
    shared = sorted(shared)
    q_patients = [r["patient"] for r in meta_q]
    q_index = {p: i for i, p in enumerate(q_patients)}
    y = np.zeros((len(shared), len(q_patients)), dtype=np.float64)
    for cohort, block in per.items():
        loc = {g: i for i, g in enumerate(block["genes"])}
        rows = [loc[g] for g in shared]
        for j, patient in enumerate(block["patients"]):
            if patient not in q_index:
                continue
            y[:, q_index[patient]] = block["logcpm"][rows, j]
    return {
        "genes": shared,
        "logcpm": y,
        "counts": None,
        "patients": q_patients,
        "cohorts": [r["cohort"] for r in meta_q],
        "quartiles": [r["quartile"] for r in meta_q],
    }


def cohort_q4_matrices(cohort_counts, meta_q):
    """Per-cohort Q1/Q4 TMM used by the combined ranks and their bootstrap."""
    prepared = {}
    for cohort in COHORTS:
        rows = [r for r in meta_q if r["cohort"] == cohort]
        n1 = sum(r["quartile"] == "Q1" for r in rows)
        n4 = sum(r["quartile"] == "Q4" for r in rows)
        if cohort not in cohort_counts or n1 < MIN_ARM or n4 < MIN_ARM:
            continue
        genes_c, samples_c, mat_c = cohort_counts[cohort]
        col = {s: i for i, s in enumerate(samples_c)}
        sub = mat_c[:, [col[r["patient"]] for r in rows]]
        genes_f, sub_f = filter_genes(genes_c, sub)
        prepared[cohort] = {
            "genes": genes_f,
            "logcpm": log_cpm(sub_f),
            "quartiles": [r["quartile"] for r in rows],
        }
    return prepared


def cohort_blocks(cohort_counts, meta_q):
    """Q1/Q4-only TMM inside each cohort with both arms at least MIN_ARM."""
    blocks = []
    for cohort in COHORTS:
        rows = [r for r in meta_q if r["cohort"] == cohort]
        n1 = sum(r["quartile"] == "Q1" for r in rows)
        n4 = sum(r["quartile"] == "Q4" for r in rows)
        if n1 < MIN_ARM or n4 < MIN_ARM:
            continue
        genes_c, samples_c, mat_c = cohort_counts[cohort]
        col = {s: i for i, s in enumerate(samples_c)}
        samples = [r["patient"] for r in rows]
        sub = mat_c[:, [col[s] for s in samples]]
        genes_f, sub_f = filter_genes(genes_c, sub)
        y = log_cpm(sub_f)
        exposure = [1.0 if r["quartile"] == "Q4" else 0.0 for r in rows]
        fit = fit_exposure(y, [cohort] * len(rows), exposure)
        if fit is None:
            continue
        blocks.append(
            {
                "cohort": cohort,
                "genes": genes_f,
                "fit": fit,
                "n": len(rows),
                "n_q1": n1,
                "n_q4": n4,
            }
        )
    return blocks


def align_mean(blocks, field, weights=None):
    """Equal-weight or IVW mean across cohort vectors. Genes in >=2 cohorts, or all cohorts if only two exist."""
    tables = []
    for block in blocks:
        tables.append({g: float(v) for g, v in zip(block["genes"], block["fit"][field])})
    if weights == "ivw":
        ses = []
        for block in blocks:
            ses.append({g: float(v) for g, v in zip(block["genes"], block["fit"]["se"])})
    need = 2 if len(blocks) >= 2 else 1
    universe = {}
    for table in tables:
        for g in table:
            universe[g] = universe.get(g, 0) + 1
    genes = sorted(g for g, k in universe.items() if k >= need)
    out = np.zeros(len(genes))
    for i, g in enumerate(genes):
        vals = []
        wts = []
        for b, table in enumerate(tables):
            if g not in table:
                continue
            vals.append(table[g])
            if weights == "ivw":
                se = ses[b].get(g, np.nan)
                wts.append(0.0 if not np.isfinite(se) or se <= 0 else 1.0 / (se * se))
        vals = np.asarray(vals, dtype=np.float64)
        if weights == "ivw":
            wts = np.asarray(wts, dtype=np.float64)
            if wts.sum() <= 0:
                out[i] = float(np.mean(vals))
            else:
                out[i] = float(np.average(vals, weights=wts))
        else:
            out[i] = float(np.mean(vals))
    return genes, out


def align_max(blocks):
    tables = []
    for block in blocks:
        tables.append({g: float(v) for g, v in zip(block["genes"], block["fit"]["logFC"])})
    need = 2 if len(blocks) >= 2 else 1
    universe = {}
    for table in tables:
        for g in table:
            universe[g] = universe.get(g, 0) + 1
    genes = sorted(g for g, k in universe.items() if k >= need)
    out = np.zeros(len(genes))
    for i, g in enumerate(genes):
        out[i] = max(table[g] for table in tables if g in table)
    return genes, out


def stouffer_from_blocks(blocks):
    tables_z = []
    for block in blocks:
        tstat = block["fit"]["t"]
        p = block["fit"]["p"]
        z = np.sign(tstat) * norm.isf(np.clip(p, 1e-300, 1) / 2.0)
        tables_z.append({g: float(v) for g, v in zip(block["genes"], z)})
    need = 2 if len(blocks) >= 2 else 1
    universe = {}
    for table in tables_z:
        for g in table:
            universe[g] = universe.get(g, 0) + 1
    genes = sorted(g for g, k in universe.items() if k >= need)
    out = np.zeros(len(genes))
    k = None
    for i, g in enumerate(genes):
        vals = [table[g] for table in tables_z if g in table]
        out[i] = float(np.sum(vals) / math.sqrt(len(vals)))
    return genes, out


def pooled_rank_vector(expr, rank):
    q4 = np.asarray([q == "Q4" for q in expr["quartiles"]])
    fit = fit_exposure(expr["logcpm"], expr["cohorts"], q4.astype(float))
    if fit is None:
        raise SystemExit(f"rank-deficient design for {rank}")
    if rank == "ols_t":
        stat = fit["t"]
    elif rank == "ols_logFC":
        stat = fit["logFC"]
    elif rank == "ols_signed_neglog10p":
        stat = np.sign(fit["logFC"]) * -np.log10(np.clip(fit["p"], 1e-300, 1))
    elif rank == "moderated_t":
        stat = moderated_t(fit)
    elif rank == "snr_residual":
        resid = residualize_cohort(expr["logcpm"], expr["cohorts"])
        stat = signal_to_noise(resid, q4)
    elif rank == "wilcoxon_z_residual":
        resid = residualize_cohort(expr["logcpm"], expr["cohorts"])
        stat = wilcoxon_z(resid, q4)
    else:
        raise KeyError(rank)
    return expr["genes"], np.asarray(stat, dtype=np.float64), fit


def write_tsv(path: Path, rows, fieldnames):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def fmt_p(p):
    if p is None or not math.isfinite(p):
        return "NA"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_n(x, d=2):
    if x is None or not math.isfinite(x):
        return "NA"
    return f"{x:.{d}f}"


def build_ranks(expr_by_universe, blocks):
    ranks = {}
    fits = {}
    for universe, expr in expr_by_universe.items():
        for rank in POOLED_RANKS:
            print(f"  rank {universe} {rank}", flush=True)
            genes, stat, fit = pooled_rank_vector(expr, rank)
            ranks[(universe, rank)] = (genes, stat)
            fits[(universe, rank)] = fit
    print("  rank within_cohort combined", flush=True)
    genes, stat = stouffer_from_blocks(blocks)
    ranks[("within_cohort", "stouffer_z")] = (genes, stat)
    genes, stat = align_mean(blocks, "logFC")
    ranks[("within_cohort", "mean_logFC")] = (genes, stat)
    genes, stat = align_mean(blocks, "logFC", weights="ivw")
    ranks[("within_cohort", "ivw_logFC")] = (genes, stat)
    genes, stat = align_max(blocks)
    ranks[("within_cohort", "max_logFC")] = (genes, stat)
    return ranks, fits


def evaluate_grid(scope, ranks, gene_sets, nperm):
    rows = []
    n_jobs = len(ranks) * len(gene_sets)
    done = 0
    for (universe, rank), (genes, stat) in ranks.items():
        names_s, stats_s = order_stats(genes, stat)
        universe_set = set(names_s)
        abs_check = None
        for sid, members in gene_sets.items():
            done += 1
            mapped, n_alias = map_genes(members, universe_set)
            if len(mapped) < MIN_SET or len(mapped) > MAX_SET:
                rows.append(
                    {
                        "scope": scope,
                        "universe": universe,
                        "rank": rank,
                        "set": sid,
                        "family": family_of(sid),
                        "size": len(mapped),
                        "n_alias": n_alias,
                        "n_genes_ranked": len(names_s),
                        "es": "",
                        "nes": "",
                        "p": "",
                        "n_leading": "",
                        "leading_top": "",
                        "status": "size_out_of_range",
                    }
                )
                continue
            rng = np.random.default_rng(child_seed(scope, universe, rank, sid, nperm))
            result = gsea_one(names_s, stats_s, mapped, nperm, rng)
            if abs_check is None and sid == "HALLMARK_INTERFERON_GAMMA_RESPONSE":
                mask = np.isin(np.asarray(names_s), mapped)
                brute = es_brute(stats_s, mask)
                fast = result["es"]
                if abs(brute - fast) > 1e-8:
                    raise SystemExit(f"ES mismatch brute {brute} fast {fast}")
                abs_check = True
            lead = result["leading"]
            # Most extreme end of the leading edge (bottom of the list for a down set).
            if result["es"] < 0:
                top = list(reversed(lead))[:20]
            else:
                top = lead[:20]
            rows.append(
                {
                    "scope": scope,
                    "universe": universe,
                    "rank": rank,
                    "set": sid,
                    "family": family_of(sid),
                    "size": result["size"],
                    "n_alias": n_alias,
                    "n_genes_ranked": result["n_genes"],
                    "es": result["es"],
                    "nes": result["nes"],
                    "p": result["p"],
                    "n_leading": result["n_leading"],
                    "leading_top": ",".join(top),
                    "leading_genes": ",".join(lead),
                    "status": "ok",
                }
            )
            if done % 25 == 0:
                print(f"    {scope} {done}/{n_jobs} {sid} NES {result['nes']:.3f}", flush=True)
    return rows


def refine_rows(rows, ranks, gene_sets, keys, nperm):
    """Recompute NES for selected (scope, universe, rank, set) keys."""
    index = {(r["scope"], r["universe"], r["rank"], r["set"]): r for r in rows}
    for key in keys:
        row = index.get(key)
        if row is None or row["status"] != "ok":
            continue
        scope, universe, rank, sid = key
        genes, stat = ranks[scope][(universe, rank)]
        names_s, stats_s = order_stats(genes, stat)
        mapped, _ = map_genes(gene_sets[sid], set(names_s))
        rng = np.random.default_rng(child_seed("refine", scope, universe, rank, sid, nperm))
        result = gsea_one(names_s, stats_s, mapped, nperm, rng)
        row["es"] = result["es"]
        row["nes"] = result["nes"]
        row["p"] = result["p"]
        row["n_leading"] = result["n_leading"]
        row["nperm"] = nperm
        lead = result["leading"]
        if result["es"] < 0:
            top = list(reversed(lead))[:20]
        else:
            top = lead[:20]
        row["leading_top"] = ",".join(top)
        row["leading_genes"] = ",".join(lead)
        print(
            f"  refine {scope} {universe} {rank} {sid} NES {result['nes']:.3f} ES {result['es']:.4f}",
            flush=True,
        )


def attach_fdr(rows):
    groups = {}
    families = {}
    for i, row in enumerate(rows):
        if row["status"] != "ok" or row["p"] == "":
            continue
        groups.setdefault((row["scope"], row["universe"], row["rank"]), []).append(i)
        families.setdefault((row["scope"], row["universe"], row["rank"], row["family"]), []).append(i)
    for idx in groups.values():
        adj = bh([float(rows[i]["p"]) for i in idx])
        for j, i in enumerate(idx):
            rows[i]["fdr_rank"] = float(adj[j])
    for idx in families.values():
        adj = bh([float(rows[i]["p"]) for i in idx])
        for j, i in enumerate(idx):
            rows[i]["fdr_family"] = float(adj[j])


def leading_logfc(row, logfc_map):
    genes = [g for g in str(row.get("leading_genes", "")).split(",") if g]
    vals = [logfc_map[g] for g in genes if g in logfc_map]
    if not vals:
        return float("nan"), float("nan"), float("nan")
    arr = np.asarray(vals)
    lineage_n = sum(g in LINEAGE for g in genes)
    return float(np.mean(arr)), float(np.median(arr)), lineage_n / len(genes)


def signature_table(expr, gene_sets, scope):
    q4 = np.asarray([q == "Q4" for q in expr["quartiles"]], dtype=float)
    fit = fit_exposure(expr["logcpm"], expr["cohorts"], q4)
    logfc = {g: float(v) for g, v in zip(expr["genes"], fit["logFC"])}
    counts = expr["counts"]
    rows = []
    gene_rows = []
    ifn_union = set()
    for sid, members in gene_sets.items():
        if family_of(sid) == "ifn":
            ifn_union.update(members)
    # Program-level logFC and bootstrap CI.
    rng = np.random.default_rng(child_seed("sig", scope))
    cohorts = expr["cohorts"]
    quartiles = expr["quartiles"]
    for sid, members in gene_sets.items():
        mapped, n_alias = map_genes(members, set(expr["genes"]))
        if len(mapped) < 5:
            continue
        ix = [expr["genes"].index(g) for g in mapped]
        score = expr["logcpm"][ix].mean(axis=0)
        beta, se, p, df = _ols_beta(score, cohorts, q4)
        boots = np.empty(N_BOOT_FC)
        for b in range(N_BOOT_FC):
            take = _stratified_indices(cohorts, quartiles, rng)
            boots[b] = _ols_beta(score[take], [cohorts[i] for i in take], q4[take])[0]
        boots = boots[np.isfinite(boots)]
        lo, hi = np.quantile(boots, [0.025, 0.975]) if boots.size >= 100 else (np.nan, np.nan)
        gene_fc = np.array([logfc[g] for g in mapped])
        rows.append(
            {
                "scope": scope,
                "set": sid,
                "family": family_of(sid),
                "size": len(mapped),
                "n_alias": n_alias,
                "signature_logFC": beta,
                "signature_se": se,
                "signature_p": p,
                "df": df,
                "signature_ci_lo": float(lo),
                "signature_ci_hi": float(hi),
                "mean_gene_logFC": float(np.mean(gene_fc)),
                "median_gene_logFC": float(np.median(gene_fc)),
                "n": len(quartiles),
                "n_q1": int((q4 == 0).sum()),
                "n_q4": int((q4 == 1).sum()),
            }
        )
    # Canonical and IFN-union genes.
    nonzero = None
    if counts is not None:
        nonzero = (counts > 0).sum(axis=1)
    nz_map = {g: int(v) for g, v in zip(expr["genes"], nonzero)} if nonzero is not None else {}
    wanted = set(CANONICAL) | (ifn_union & set(expr["genes"]))
    for g in sorted(wanted):
        if g not in logfc:
            continue
        nz = nz_map.get(g, "")
        gene_rows.append(
            {
                "scope": scope,
                "gene": g,
                "logFC": logfc[g],
                "t": float(fit["t"][expr["genes"].index(g)]),
                "p": float(fit["p"][expr["genes"].index(g)]),
                "n_nonzero": nz,
                "sparse": int(nz < 8) if nz != "" else "",
                "canonical": int(g in set(CANONICAL)),
                "in_ifn_union": int(g in ifn_union),
            }
        )
    return rows, gene_rows, logfc, fit


def _ols_beta(y, cohorts, exposure):
    fit_y = np.asarray(y, dtype=np.float64).reshape(1, -1)
    fit = fit_exposure(fit_y, cohorts, exposure)
    if fit is None:
        return float("nan"), float("nan"), float("nan"), float("nan")
    return float(fit["logFC"][0]), float(fit["se"][0]), float(fit["p"][0]), int(fit["df"])


def _stratified_indices(cohorts, quartiles, rng):
    chunks = []
    pairs = sorted(set(zip(cohorts, quartiles)))
    for cohort, quartile in pairs:
        w = [i for i, (c, q) in enumerate(zip(cohorts, quartiles)) if c == cohort and q == quartile]
        if not w:
            continue
        chunks.append(rng.choice(np.asarray(w), size=len(w), replace=True))
    return np.concatenate(chunks)


def bootstrap_nes(expr, genes_in_set, nperm, rng):
    """Patient bootstrap of NES for one set on the OLS t ranking. TMM stays fixed."""
    y = expr["logcpm"]
    cohorts = expr["cohorts"]
    quartiles = expr["quartiles"]
    out = np.empty(N_BOOT_NES)
    for b in range(N_BOOT_NES):
        take = _stratified_indices(cohorts, quartiles, rng)
        y_b = y[:, take]
        co_b = [cohorts[i] for i in take]
        q_b = np.asarray([1.0 if quartiles[i] == "Q4" else 0.0 for i in take])
        fit = fit_exposure(y_b, co_b, q_b)
        if fit is None:
            out[b] = np.nan
            continue
        names_s, stats_s = order_stats(expr["genes"], fit["t"])
        result = gsea_one(
            names_s,
            stats_s,
            genes_in_set,
            nperm,
            np.random.default_rng(rng.integers(0, 2**31 - 1)),
        )
        out[b] = result["nes"] if result else np.nan
        if (b + 1) % 25 == 0:
            print(f"    bootstrap {b+1}/{N_BOOT_NES}", flush=True)
    ok = out[np.isfinite(out)]
    if ok.size < 20:
        return {"n_boot_ok": int(ok.size), "nes_boot_median": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan")}
    lo, hi = np.quantile(ok, [0.025, 0.975])
    return {
        "n_boot_ok": int(ok.size),
        "nes_boot_median": float(np.median(ok)),
        "ci_lo": float(lo),
        "ci_hi": float(hi),
    }


def eligible(row, family):
    if row["status"] != "ok" or row["family"] != family:
        return False
    try:
        size = int(row["size"])
        nes = float(row["nes"])
    except (TypeError, ValueError):
        return False
    return HEADLINE_MIN <= size <= MAX_SET and math.isfinite(nes)


def pick_most_negative(rows, family):
    pool = [r for r in rows if eligible(r, family)]
    if not pool:
        return None
    pool.sort(key=lambda r: (float(r["nes"]), -int(r["size"]), r["set"], r["rank"], r["universe"]))
    return pool[0]


def forest(path, labels, estimates, lo, hi, colors, xlabel, title):
    fig_h = max(3.2, 0.38 * len(labels) + 1.2)
    fig, ax = plt.subplots(figsize=(8.2, fig_h))
    y = np.arange(len(labels))[::-1]
    for i, (est, a, b, color) in enumerate(zip(estimates, lo, hi, colors)):
        if est is None or not math.isfinite(est):
            continue
        if a is not None and b is not None and math.isfinite(a) and math.isfinite(b):
            ax.plot([a, b], [y[i], y[i]], color=color, lw=1.6, solid_capstyle="round")
        ax.plot(est, y[i], "o", color=color, ms=5.5)
    ax.axvline(0, color="#888888", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    fig.savefig(path.with_suffix(".pdf"))
    plt.close(fig)


def short_label(sid):
    return (
        sid.replace("HALLMARK_INTERFERON_", "Hallmark IFN ")
        .replace("HALLMARK_IFN_", "Hallmark IFN ")
        .replace("REACTOME_INTERFERON_", "Reactome IFN ")
        .replace("REACTOME_", "Reactome ")
        .replace("KEGG_", "KEGG ")
        .replace("WP_INTERFERON_", "WP IFN ")
        .replace("WP_", "WP ")
        .replace("GO_", "GO ")
        .replace("CUSTOM_", "")
        .replace("_", " ")
    )


def write_finding(path, context):
    b = context["baseline"]
    win = context["ifn_max"]
    win_drop = context["ifn_max_drop_same"]
    drop_best = context["ifn_max_on_drop"]
    mhc = context["mhc_max"]
    chem = context["chem_max"]
    sig = context["sig_max"]
    gene = context["gene_max"]
    gene_drop = context["gene_max_drop"]
    lines = []
    a = lines.append
    a("# Concordant-4 malignant pseudobulk: Q4 vs Q1 IFN / MHC / chemokine sweep")
    a("")
    a("ADDITIVE. **CLDN4-only.** Cohorts are the locked concordant four:")
    a("GSE123902 + GSE131907 + GSE205335 + GSE189357. Quartiles are the locked")
    a("within-cohort malignant CLDN4 %pos labels. The unit is the patient,")
    a("donor, or sample. P4001 is already out of the count matrix. Expression")
    a(f"n is **{context['n_all']}**. Q1 n=**{context['n_q1']}**, Q4 n=**{context['n_q4']}**.")
    a("GSE189357 Q4 has 2 units, so that cohort is skipped for within-cohort")
    a("ranks and kept in the pooled model. CLDN4 is removed from every set.")
    a("")
    a("This does not replace the locked fgsea result. It asks whether another")
    a("pre-specified ranking or interferon set is more negative than Hallmark")
    a("IFN-gamma on the cohort-adjusted t, and whether that result is still")
    a("negative after GSE205335 is removed.")
    a("")
    a("## Reproduction of the locked ranking")
    a("")
    a("Zero-filled Q1/Q4 counts, one TMM, OLS `log2(CPM+1) ~ cohort + Q4`.")
    a(f"Genes ranked: **{context['n_genes']}**. CLDN4 logFC **{fmt_n(context['cldn4_logfc'], 3)}**")
    a(f"(t {fmt_n(context['cldn4_t'], 2)}, p {fmt_p(context['cldn4_p'])}).")
    a("The prior fgsea point was ES −0.687 and NES −3.85 for Hallmark IFN-gamma.")
    a(f"This run: ES **{fmt_n(b['es'], 3)}**, NES **{fmt_n(float(b['nes']), 2)}**")
    a(f"(permutation p {fmt_p(float(b['p']))}, {b.get('nperm', NPERM_REFINE)} permutations).")
    if context["baseline_boot"]:
        boot = context["baseline_boot"]
        a(f"Patient bootstrap 95% (B={boot['n_boot_ok']}): "
          f"**{fmt_n(boot['ci_lo'], 2)} to {fmt_n(boot['ci_hi'], 2)}**.")
    a("Positive NES would mean the set is higher in CLDN4-high. The sign here is down.")
    a("")
    a("## What was swept")
    a("")
    a("Pooled ranks, each on the zero-filled Q4 matrix and on the shared-gene")
    a("within-cohort TMM matrix: OLS t, OLS logFC, signed −log10 p, limma-style")
    a("moderated t, cohort-residualized signal-to-noise, cohort-residualized")
    a("Wilcoxon z. Cohort-combined ranks, equal weight, arms with n≥3")
    a("(GSE189357 binary dropped): Stouffer z, mean logFC, inverse-variance")
    a("logFC, and the max within-cohort logFC (a gene has to be down in every")
    a("tested cohort to sit at the bottom). Gene sets are the locked Hallmark")
    a("IFN lists, their intersection, the same lists after a frozen lymphocyte/")
    a("myeloid identity filter, a frozen epithelial ISG/APM list, Enrichr")
    a("Reactome 2022 / KEGG 2021 / WikiPathways 2021 / GO Biological Process")
    a("2023 interferon, MHC-I, and chemokine sets, and a frozen inflammatory")
    a("chemokine-ligand list. Headline sets have 15–500 genes after")
    a("intersection with the ranked universe. GSEA weight is |stat|^1.")
    a("NES = ES / mean |null ES| of the same sign. The null is gene-set")
    a("permutation. FDR is Benjamini–Hochberg inside one rank, not across the")
    a("grid. The grid maximum is descriptive.")
    a("")
    a("## Most negative interferon NES")
    a("")
    a(f"Full concordant-4 grid ({context['n_ifn_full']} eligible IFN tests):")
    a(f"**{win['set']}**, rank `{win['rank']}`, universe `{win['universe']}`,")
    a(f"size {win['size']}, ES {fmt_n(float(win['es']), 3)}, NES **{fmt_n(float(win['nes']), 2)}**,")
    a(f"p {fmt_p(float(win['p']))}, FDR within rank {fmt_p(float(win.get('fdr_rank', float('nan'))))},")
    a(f"FDR within IFN sets on that rank {fmt_p(float(win.get('fdr_family', float('nan'))))}.")
    a(f"Leading edge {win['n_leading']} genes, mean OLS logFC {fmt_n(float(win['lead_mean']), 3)},")
    a(f"lineage-identity fraction {fmt_n(float(win['lineage_frac']), 2)}.")
    if context["ifn_boot"]:
        boot = context["ifn_boot"]
        a(f"Patient bootstrap of this rank, TMM held fixed, 95% "
          f"(B={boot['n_boot_ok']}): **{fmt_n(boot['ci_lo'], 2)} to {fmt_n(boot['ci_hi'], 2)}**.")
    delta = float(win["nes"]) - float(b["nes"])
    if win["set"] == b["set"] and win["rank"] == b["rank"] and win["universe"] == b["universe"]:
        a("That row is the grid maximum. No other pre-specified rank or interferon")
        a("set was more negative.")
    else:
        a(f"Difference from the reproduced Hallmark IFN-gamma OLS-t NES: **{delta:+.2f}**.")
    a("")
    a("Same configuration with GSE205335 removed and the matrix re-fit:")
    if win_drop is None:
        a("not estimable.")
    else:
        a(f"NES **{fmt_n(float(win_drop['nes']), 2)}** (size {win_drop['size']}, "
          f"p {fmt_p(float(win_drop['p']))}).")
    a("")
    a("Most negative IFN NES on the GSE205335-excluded grid itself:")
    a(f"**{drop_best['set']}**, `{drop_best['rank']}` / `{drop_best['universe']}`,")
    a(f"size {drop_best['size']}, NES **{fmt_n(float(drop_best['nes']), 2)}**,")
    a(f"p {fmt_p(float(drop_best['p']))}.")
    if context["drop_boot"]:
        boot = context["drop_boot"]
        a(f"Bootstrap 95% (B={boot['n_boot_ok']}): **{fmt_n(boot['ci_lo'], 2)} to {fmt_n(boot['ci_hi'], 2)}**.")
    a("")
    a("Hallmark IFN-gamma on the locked OLS-t ranking, GSE205335 removed:")
    a(f"NES **{fmt_n(float(context['baseline_drop']['nes']), 2)}**.")
    a("The prior leave-one-out value was −2.37. GSE131907 alone was flat in that run;")
    a("a more negative pooled NES that shrinks after dropping GSE205335 is carried")
    a("in part by that cohort.")
    a("")
    a("## MHC and chemokine")
    a("")
    a(f"Most negative MHC NES in the full grid: **{mhc['set']}**, `{mhc['rank']}` / "
      f"`{mhc['universe']}`, size {mhc['size']}, NES **{fmt_n(float(mhc['nes']), 2)}**.")
    mhc_d = context["mhc_same_drop"]
    if mhc_d is not None:
        a(f"Same configuration without GSE205335: NES **{fmt_n(float(mhc_d['nes']), 2)}**.")
    a(f"Most negative chemokine NES in the full grid: **{chem['set']}**, `{chem['rank']}` / "
      f"`{chem['universe']}`, size {chem['size']}, NES **{fmt_n(float(chem['nes']), 2)}**.")
    chem_d = context["chem_same_drop"]
    if chem_d is not None:
        a(f"Same configuration without GSE205335: NES **{fmt_n(float(chem_d['nes']), 2)}**.")
    a("")
    a("## logFC of the IFN program")
    a("")
    a("Signature logFC is the cohort-adjusted OLS beta of the mean")
    a("`log2(TMM-CPM+1)` across genes in the set, on the zero-filled Q4 matrix.")
    a("It is in log2 units. Mean and median gene logFC are the unweighted")
    a("summaries of the gene-level betas. Leading-edge means are not used to")
    a("pick the maximum.")
    a("")
    a(f"Most negative IFN signature logFC: **{sig['set']}** "
      f"**{fmt_n(float(sig['signature_logFC']), 3)}** "
      f"(95% {fmt_n(float(sig['signature_ci_lo']), 3)} to {fmt_n(float(sig['signature_ci_hi']), 3)}, "
      f"p {fmt_p(float(sig['signature_p']))}, size {sig['size']}).")
    a(f"Median gene logFC in that set: {fmt_n(float(sig['median_gene_logFC']), 3)}. "
      f"Mean gene logFC: {fmt_n(float(sig['mean_gene_logFC']), 3)}.")
    sig_d = context["sig_same_drop"]
    if sig_d is not None:
        a(f"Same set without GSE205335: signature logFC **{fmt_n(float(sig_d['signature_logFC']), 3)}** "
          f"(95% {fmt_n(float(sig_d['signature_ci_lo']), 3)} to {fmt_n(float(sig_d['signature_ci_hi']), 3)}).")
    a("")
    a(f"Most negative non-sparse canonical IFN/MHC/chemokine gene "
      f"(nonzero in ≥8 of the Q4/Q1 units): **{gene['gene']}** "
      f"logFC **{fmt_n(float(gene['logFC']), 3)}** (p {fmt_p(float(gene['p']))}, "
      f"nonzero {gene['n_nonzero']}).")
    if gene_drop is not None:
        a(f"Without GSE205335 the most negative non-sparse canonical gene is "
          f"**{gene_drop['gene']}** logFC **{fmt_n(float(gene_drop['logFC']), 3)}**.")
    a("")
    a("## Read this as a sweep")
    a("")
    a("The baseline OLS-t Hallmark IFN-gamma test is the pre-specified comparison")
    a("to the prior NES of about −3.85. The grid records how far |NES| and |logFC|")
    a("move under the other declared ranks and sets. A gain that disappears when")
    a("GSE205335 is removed is not a four-cohort result. Within-cohort Wilcoxon")
    a("and Stouffer ranks skip GSE189357 because Q4 has n=2. Normalization factors")
    a("are not re-estimated inside the patient bootstrap.")
    a("")
    a("Tables: `results/tables/nes_sweep.tsv`, `signature_logfc.tsv`, `gene_logfc.tsv`, `headline.tsv`.")
    a("")
    path.write_text("\n".join(lines) + "\n")


def plot_all(context):
    # IFN sets at their most negative rank, full scope, with drop-GSE205335 overlay when the same key exists.
    full = [r for r in context["rows"] if r["scope"] == "full" and eligible(r, "ifn")]
    best = {}
    for row in full:
        cur = best.get(row["set"])
        if cur is None or float(row["nes"]) < float(cur["nes"]):
            best[row["set"]] = row
    ordered = sorted(best.values(), key=lambda r: float(r["nes"]))
    drop_index = {
        (r["universe"], r["rank"], r["set"]): r
        for r in context["rows"]
        if r["scope"] == "drop_GSE205335" and r["status"] == "ok"
    }
    labels, est, lo, hi, colors = [], [], [], [], []
    for row in ordered:
        twin = drop_index.get((row["universe"], row["rank"], row["set"]))
        labels.append(f"{short_label(row['set'])}  [{row['rank']}]")
        est.append(float(row["nes"]))
        lo.append(float("nan"))
        hi.append(float("nan"))
        colors.append("#1f4e79")
        labels.append("    without GSE205335")
        if twin is None or twin["nes"] == "":
            est.append(float("nan"))
        else:
            est.append(float(twin["nes"]))
        lo.append(float("nan"))
        hi.append(float("nan"))
        colors.append("#c47b2b")
    forest(
        FIG / "forest_ifn_nes_best_rank.png",
        labels, est, lo, hi, colors,
        "NES, Q4 vs Q1 (negative = down in CLDN4-high)",
        "Most negative IFN NES per set, then the same rank without GSE205335",
    )

    # Hallmark IFN-gamma across ranks.
    gamma = [
        r for r in context["rows"]
        if r["set"] == "HALLMARK_INTERFERON_GAMMA_RESPONSE" and r["status"] == "ok" and r["scope"] == "full"
    ]
    gamma.sort(key=lambda r: (r["universe"], r["rank"]))
    labels, est, lo, hi, colors = [], [], [], [], []
    for row in gamma:
        twin = drop_index.get((row["universe"], row["rank"], row["set"]))
        labels.append(f"{row['universe']} · {row['rank']}")
        est.append(float(row["nes"]))
        lo.append(float("nan"))
        hi.append(float("nan"))
        colors.append("#1f4e79")
        labels.append("    without GSE205335")
        est.append(float(twin["nes"]) if twin is not None and twin["nes"] != "" else float("nan"))
        lo.append(float("nan"))
        hi.append(float("nan"))
        colors.append("#c47b2b")
    forest(
        FIG / "forest_ifng_by_rank.png",
        labels, est, lo, hi, colors,
        "NES, Hallmark IFN-gamma",
        "Hallmark IFN-gamma NES across the rank grid",
    )

    sig_rows = [r for r in context["sig_rows"] if r["scope"] == "full" and r["family"] == "ifn" and int(r["size"]) >= HEADLINE_MIN]
    sig_rows.sort(key=lambda r: float(r["signature_logFC"]))
    drop_sig = {(r["set"]): r for r in context["sig_rows"] if r["scope"] == "drop_GSE205335"}
    labels, est, lo, hi, colors = [], [], [], [], []
    for row in sig_rows:
        labels.append(short_label(row["set"]))
        est.append(float(row["signature_logFC"]))
        lo.append(float(row["signature_ci_lo"]))
        hi.append(float(row["signature_ci_hi"]))
        colors.append("#1f4e79")
        twin = drop_sig.get(row["set"])
        labels.append("    without GSE205335")
        if twin is None:
            est.append(float("nan"))
            lo.append(float("nan"))
            hi.append(float("nan"))
        else:
            est.append(float(twin["signature_logFC"]))
            lo.append(float(twin["signature_ci_lo"]))
            hi.append(float(twin["signature_ci_hi"]))
        colors.append("#c47b2b")
    forest(
        FIG / "forest_signature_logfc.png",
        labels, est, lo, hi, colors,
        "Signature log2 fold change, Q4 vs Q1",
        "Mean-logCPM IFN signature, cohort-adjusted OLS",
    )

    genes = [
        r for r in context["gene_rows"]
        if r["scope"] == "full" and int(r["canonical"]) == 1 and int(r["sparse"]) == 0
    ]
    genes.sort(key=lambda r: float(r["logFC"]))
    fig_h = max(3.2, 0.28 * len(genes) + 1.0)
    fig, ax = plt.subplots(figsize=(7.4, fig_h))
    y = np.arange(len(genes))[::-1]
    vals = [float(r["logFC"]) for r in genes]
    colors = ["#1f4e79" if v < 0 else "#b04a3a" for v in vals]
    ax.axvline(0, color="#888888", lw=0.8)
    ax.hlines(y, 0, vals, color="#b7c3d0", lw=1)
    ax.scatter(vals, y, c=colors, s=22, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([r["gene"] for r in genes], fontsize=8)
    ax.set_xlabel("OLS logFC, Q4 vs Q1")
    ax.set_title("Canonical IFN / MHC / chemokine genes")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "lollipop_canonical_logfc.png", dpi=160)
    fig.savefig(FIG / "lollipop_canonical_logfc.pdf")
    plt.close(fig)


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    print("Loading counts and units", flush=True)
    units = read_units(DATA / "tnk_units.tsv")
    cohort_counts = {}
    for cohort in COHORTS:
        genes, samples, mat = read_counts(DATA / f"{cohort}_malignant_counts.tsv.gz")
        need = [r["patient"] for r in units if r["cohort"] == cohort]
        missing = sorted(set(need) - set(samples))
        if missing:
            raise SystemExit(f"{cohort} missing {missing}")
        cohort_counts[cohort] = (genes, samples, mat)
        print(f"  {cohort} genes {len(genes)} samples {len(samples)}", flush=True)

    meta_q = [r for r in units if r["quartile"] in {"Q1", "Q4"}]
    print(f"Q4 vs Q1 n={len(meta_q)} Q1={sum(r['quartile']=='Q1' for r in meta_q)} "
          f"Q4={sum(r['quartile']=='Q4' for r in meta_q)}", flush=True)

    locked = load_locked_sets()
    public, public_desc = load_public_sets()
    gene_sets = {}
    gene_sets.update(locked)
    gene_sets.update(public)
    # Write the exact sets used, including derived Hallmark subsets.
    with (TAB / "gene_sets_used.gmt").open("w") as handle:
        for sid, genes in gene_sets.items():
            handle.write(sid + "\t" + family_of(sid) + "\t" + "\t".join(genes) + "\n")

    def scope_inputs(drop):
        rows_all = [r for r in units if not (drop and r["cohort"] == "GSE205335")]
        rows_q = [r for r in rows_all if r["quartile"] in {"Q1", "Q4"}]
        counts = {c: cohort_counts[c] for c in COHORTS if not (drop and c == "GSE205335")}
        # bind_counts expects every cohort key that appears; pass only kept cohorts.
        expr_zero = prepare_zerofill(counts, rows_q)
        expr_shared = prepare_shared(counts, rows_all, rows_q)
        blocks = cohort_blocks(counts, rows_q)
        return rows_all, rows_q, expr_zero, expr_shared, blocks

    print("Preparing full scope", flush=True)
    _all, q_full, zero_full, shared_full, blocks_full = scope_inputs(False)
    print(f"  zerofill genes {len(zero_full['genes'])} shared genes {len(shared_full['genes'])} "
          f"within-cohort blocks {len(blocks_full)}", flush=True)

    fit0 = fit_exposure(
        zero_full["logcpm"],
        zero_full["cohorts"],
        np.asarray([1.0 if q == "Q4" else 0.0 for q in zero_full["quartiles"]]),
    )
    cldn_i = zero_full["genes"].index("CLDN4")
    cldn_fc = float(fit0["logFC"][cldn_i])
    cldn_t = float(fit0["t"][cldn_i])
    cldn_p = float(fit0["p"][cldn_i])
    print(f"CLDN4 logFC {cldn_fc:.4f} t {cldn_t:.3f} p {cldn_p:.4g} n_genes {len(zero_full['genes'])} df {fit0['df']}", flush=True)
    if len(zero_full["genes"]) != 21604:
        raise SystemExit(f"Gene count {len(zero_full['genes'])} != 21604; TMM/filter does not match the locked run.")
    if not math.isfinite(cldn_fc) or cldn_fc <= 0 or abs(cldn_fc - 1.673) > 0.02:
        raise SystemExit(f"CLDN4 logFC {cldn_fc} does not match the locked 1.673.")

    print("Preparing drop-GSE205335 scope", flush=True)
    _all_d, q_drop, zero_drop, shared_drop, blocks_drop = scope_inputs(True)
    print(f"  zerofill genes {len(zero_drop['genes'])} shared {len(shared_drop['genes'])}", flush=True)

    print("Building ranks", flush=True)
    ranks = {}
    ranks["full"], fits_full = build_ranks(
        {"zerofill_q4": zero_full, "shared_within": shared_full}, blocks_full
    )
    ranks["drop_GSE205335"], _fits_drop = build_ranks(
        {"zerofill_q4": zero_drop, "shared_within": shared_drop}, blocks_drop
    )

    # Deterministic ES check against the published Hallmark IFN-gamma ES.
    base_genes, base_stat = ranks["full"][("zerofill_q4", "ols_t")]
    names_s, stats_s = order_stats(base_genes, base_stat)
    gamma_genes, _ = map_genes(gene_sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"], set(names_s))
    mask = np.isin(np.asarray(names_s), gamma_genes)
    es_gamma = es_brute(stats_s, mask)
    print(f"Hallmark IFN-gamma ES {es_gamma:.6f} (locked fgsea ES -0.687153)", flush=True)
    if abs(es_gamma - (-0.687153053169553)) > 0.01:
        raise SystemExit("IFN-gamma ES does not reproduce the locked fgsea ranking.")

    print("Sweep", flush=True)
    rows = []
    rows.extend(evaluate_grid("full", ranks["full"], gene_sets, NPERM_SWEEP))
    rows.extend(evaluate_grid("drop_GSE205335", ranks["drop_GSE205335"], gene_sets, NPERM_SWEEP))

    def top_keys(scope, family, k=8):
        pool = [r for r in rows if r["scope"] == scope and eligible(r, family)]
        pool.sort(key=lambda r: float(r["nes"]))
        return [(r["scope"], r["universe"], r["rank"], r["set"]) for r in pool[:k]]

    refine_keys = set()
    refine_keys.update(top_keys("full", "ifn"))
    refine_keys.update(top_keys("drop_GSE205335", "ifn"))
    refine_keys.update(top_keys("full", "mhc", 4))
    refine_keys.update(top_keys("full", "chemokine", 4))
    refine_keys.add(("full", "zerofill_q4", "ols_t", "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
    refine_keys.add(("drop_GSE205335", "zerofill_q4", "ols_t", "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
    refine_keys.add(("full", "zerofill_q4", "ols_t", "HALLMARK_INTERFERON_ALPHA_RESPONSE"))
    print(f"Refining {len(refine_keys)} headline candidates at {NPERM_REFINE} perms", flush=True)
    # refine_rows looks up ranks[scope]; pass the nested dict via a small adapter.
    class RankView(dict):
        pass

    refine_rows(rows, ranks, gene_sets, refine_keys, NPERM_REFINE)
    attach_fdr(rows)

    print("Signature and gene logFC", flush=True)
    sig_full, gene_full, logfc_full, _ = signature_table(zero_full, gene_sets, "full")
    sig_drop, gene_drop, logfc_drop, _ = signature_table(zero_drop, gene_sets, "drop_GSE205335")
    sig_rows = sig_full + sig_drop
    gene_rows = gene_full + gene_drop

    for row in rows:
        if row["status"] != "ok":
            row["lead_mean"] = ""
            row["lead_median"] = ""
            row["lineage_frac"] = ""
            continue
        fcmap = logfc_full if row["scope"] == "full" else logfc_drop
        mean_fc, med_fc, frac = leading_logfc(row, fcmap)
        row["lead_mean"] = mean_fc
        row["lead_median"] = med_fc
        row["lineage_frac"] = frac

    ifn_max = pick_most_negative([r for r in rows if r["scope"] == "full"], "ifn")
    mhc_max = pick_most_negative([r for r in rows if r["scope"] == "full"], "mhc")
    chem_max = pick_most_negative([r for r in rows if r["scope"] == "full"], "chemokine")
    ifn_drop_best = pick_most_negative([r for r in rows if r["scope"] == "drop_GSE205335"], "ifn")

    def find_row(scope, universe, rank, sid):
        for row in rows:
            if (row["scope"], row["universe"], row["rank"], row["set"]) == (scope, universe, rank, sid):
                return row
        return None

    baseline = find_row("full", "zerofill_q4", "ols_t", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
    baseline_drop = find_row("drop_GSE205335", "zerofill_q4", "ols_t", "HALLMARK_INTERFERON_GAMMA_RESPONSE")
    ifn_same_drop = find_row("drop_GSE205335", ifn_max["universe"], ifn_max["rank"], ifn_max["set"])
    mhc_same = find_row("drop_GSE205335", mhc_max["universe"], mhc_max["rank"], mhc_max["set"])
    chem_same = find_row("drop_GSE205335", chem_max["universe"], chem_max["rank"], chem_max["set"])

    ifn_sig = [r for r in sig_full if r["family"] == "ifn" and int(r["size"]) >= HEADLINE_MIN]
    ifn_sig.sort(key=lambda r: float(r["signature_logFC"]))
    sig_max = ifn_sig[0]
    sig_same = next(r for r in sig_drop if r["set"] == sig_max["set"])

    canon = [r for r in gene_full if int(r["canonical"]) == 1 and str(r["sparse"]) == "0"]
    canon.sort(key=lambda r: float(r["logFC"]))
    gene_max = canon[0]
    canon_d = [r for r in gene_drop if int(r["canonical"]) == 1 and str(r["sparse"]) == "0"]
    canon_d.sort(key=lambda r: float(r["logFC"]))
    gene_max_drop = canon_d[0] if canon_d else None

    print("Bootstrapping headline NES", flush=True)
    cohort_prep = {
        "full": cohort_q4_matrices(cohort_counts, q_full),
        "drop_GSE205335": cohort_q4_matrices(
            {c: v for c, v in cohort_counts.items() if c != "GSE205335"}, q_drop
        ),
    }

    def boot_rank(scope_name, expr, universe, rank, sid):
        # Pooled ranks resample columns of the matrix that defined the rank.
        # Within-cohort combined ranks resample inside each cohort's own
        # Q1/Q4 TMM. Those are different matrices; do not substitute one.
        rng = np.random.default_rng(child_seed("boot", scope_name, universe, rank, sid))
        if universe == "within_cohort":
            universe_genes = set()
            for block in expr.values():
                universe_genes.update(block["genes"])
            members, _ = map_genes(gene_sets[sid], universe_genes)
            return _boot_combined(expr, rank, members, rng)
        members, _ = map_genes(gene_sets[sid], set(expr["genes"]))
        return _boot_pooled(expr, rank, members, rng)

    def _boot_pooled(expr, rank, members, rng):
        y = expr["logcpm"]
        cohorts = expr["cohorts"]
        quartiles = expr["quartiles"]
        out = np.empty(N_BOOT_NES)
        for b in range(N_BOOT_NES):
            take = _stratified_indices(cohorts, quartiles, rng)
            expr_b = {
                "genes": expr["genes"],
                "logcpm": y[:, take],
                "cohorts": [cohorts[i] for i in take],
                "quartiles": [quartiles[i] for i in take],
            }
            genes, stat, _fit = pooled_rank_vector(expr_b, rank)
            names_s, stats_s = order_stats(genes, stat)
            result = gsea_one(
                names_s, stats_s, members, NPERM_BOOT,
                np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
            )
            out[b] = result["nes"] if result else np.nan
            if (b + 1) % 20 == 0:
                print(f"    pooled bootstrap {b+1}/{N_BOOT_NES}", flush=True)
        return _boot_summary(out)

    def _boot_combined(expr, rank, members, rng):
        # Resample inside each cohort's own Q1/Q4 TMM, then rebuild the
        # combined rank. expr here is a dict cohort -> genes, logcpm, quartiles.
        # Passing the pooled matrix instead changes the normalization and is
        # not this rank's bootstrap.
        if not isinstance(expr, dict) or "logcpm" in expr:
            raise SystemExit("within-cohort bootstrap needs per-cohort TMM blocks")
        out = np.empty(N_BOOT_NES)
        for b in range(N_BOOT_NES):
            blocks = []
            for cohort, block in expr.items():
                quartiles = block["quartiles"]
                y = block["logcpm"]
                take = []
                ok_arm = True
                for q in ("Q1", "Q4"):
                    ix = np.flatnonzero(np.asarray(quartiles) == q)
                    if ix.size < MIN_ARM:
                        ok_arm = False
                        break
                    take.append(rng.choice(ix, size=ix.size, replace=True))
                if not ok_arm:
                    continue
                take = np.concatenate(take)
                exposure = (np.asarray(quartiles)[take] == "Q4").astype(float)
                fit = fit_exposure(y[:, take], [cohort] * take.size, exposure)
                if fit is None:
                    continue
                blocks.append({"genes": block["genes"], "fit": fit, "cohort": cohort})
            if len(blocks) < 2:
                out[b] = np.nan
                continue
            if rank == "stouffer_z":
                genes, stat = stouffer_from_blocks(blocks)
            elif rank == "mean_logFC":
                genes, stat = align_mean(blocks, "logFC")
            elif rank == "ivw_logFC":
                genes, stat = align_mean(blocks, "logFC", weights="ivw")
            elif rank == "max_logFC":
                genes, stat = align_max(blocks)
            else:
                out[b] = np.nan
                continue
            names_s, stats_s = order_stats(genes, stat)
            result = gsea_one(
                names_s, stats_s, members, NPERM_BOOT,
                np.random.default_rng(int(rng.integers(0, 2**31 - 1))),
            )
            out[b] = result["nes"] if result else np.nan
            if (b + 1) % 20 == 0:
                print(f"    combined bootstrap {b+1}/{N_BOOT_NES}", flush=True)
        return _boot_summary(out)

    boot_targets = []
    boot_targets.append(("baseline", "full", zero_full, "zerofill_q4", "ols_t", "HALLMARK_INTERFERON_GAMMA_RESPONSE"))
    boot_targets.append(("ifn_max", "full", zero_full if ifn_max["universe"] != "shared_within" else shared_full,
                         ifn_max["universe"], ifn_max["rank"], ifn_max["set"]))
    boot_targets.append((
        "ifn_drop_best", "drop_GSE205335",
        zero_drop if ifn_drop_best["universe"] != "shared_within" else shared_drop,
        ifn_drop_best["universe"], ifn_drop_best["rank"], ifn_drop_best["set"],
    ))
    # Deduplicate identical jobs.
    seen_jobs = set()
    boots = {}
    for name, scope_name, expr, universe, rank, sid in boot_targets:
        job = (scope_name, universe, rank, sid)
        if job in seen_jobs:
            boots[name] = boots[[k for k, v in boots.items() if False]]
        # Always compute; cache by job.
        cache_key = job
        if cache_key not in boots:
            print(f"Bootstrap {name}: {job}", flush=True)
            expr_use = expr
            if universe == "within_cohort":
                expr_use = cohort_prep[scope_name]
            boots[cache_key] = boot_rank(scope_name, expr_use, universe, rank, sid)
        boots[name] = boots[cache_key]

    # Fix the clumsy cache: boots[name] and boots[job] mixed. Rebuild cleanly below if needed.
    # The loop above stores both. Retrieve by name. The line `boots[name] = boots[[...]]` is dead
    # because seen_jobs is never filled. Fill cache properly by using a side dict.
    # (left as a straightforward recompute-safe map)

    headline = []

    def add_head(kind, row, boot=None):
        if row is None:
            return
        item = {
            "kind": kind,
            "scope": row.get("scope", ""),
            "set": row.get("set", row.get("gene", "")),
            "universe": row.get("universe", "zerofill_q4"),
            "rank": row.get("rank", "ols_logFC"),
            "size": row.get("size", ""),
            "nes": row.get("nes", ""),
            "es": row.get("es", ""),
            "p": row.get("p", row.get("signature_p", "")),
            "fdr_family": row.get("fdr_family", ""),
            "signature_logFC": row.get("signature_logFC", row.get("logFC", "")),
            "ci_lo": "" if not boot else boot["ci_lo"],
            "ci_hi": "" if not boot else boot["ci_hi"],
            "n_boot_ok": "" if not boot else boot["n_boot_ok"],
        }
        headline.append(item)

    add_head("baseline_ifng_ols_t", baseline, boots.get("baseline"))
    add_head("baseline_ifng_ols_t_drop_GSE205335", baseline_drop, None)
    add_head("ifn_nes_max", ifn_max, boots.get("ifn_max"))
    add_head("ifn_nes_max_same_config_drop_GSE205335", ifn_same_drop, None)
    add_head("ifn_nes_max_on_drop_grid", ifn_drop_best, boots.get("ifn_drop_best"))
    add_head("mhc_nes_max", mhc_max, None)
    add_head("chemokine_nes_max", chem_max, None)
    add_head("ifn_signature_logFC_max", {**sig_max, "scope": "full", "universe": "zerofill_q4", "rank": "signature_mean_logCPM", "nes": "", "es": ""}, None)
    add_head("canonical_gene_logFC_max", {
        "scope": "full", "set": gene_max["gene"], "universe": "zerofill_q4", "rank": "ols_logFC",
        "size": 1, "nes": "", "es": "", "p": gene_max["p"], "signature_logFC": gene_max["logFC"],
    }, None)

    fields = [
        "scope", "universe", "rank", "set", "family", "size", "n_alias", "n_genes_ranked",
        "es", "nes", "p", "fdr_rank", "fdr_family", "n_leading", "lead_mean", "lead_median",
        "lineage_frac", "leading_top", "status", "nperm",
    ]
    for row in rows:
        row.setdefault("fdr_rank", "")
        row.setdefault("fdr_family", "")
        row.setdefault("nperm", NPERM_SWEEP if row["status"] == "ok" else "")
        row.pop("leading_genes", None)
    write_tsv(TAB / "nes_sweep.tsv", rows, fields)
    write_tsv(
        TAB / "signature_logfc.tsv",
        sig_rows,
        ["scope", "set", "family", "size", "n_alias", "n", "n_q1", "n_q4", "signature_logFC",
         "signature_se", "signature_p", "df", "signature_ci_lo", "signature_ci_hi",
         "mean_gene_logFC", "median_gene_logFC"],
    )
    write_tsv(
        TAB / "gene_logfc.tsv",
        gene_rows,
        ["scope", "gene", "logFC", "t", "p", "n_nonzero", "sparse", "canonical", "in_ifn_union"],
    )
    write_tsv(
        TAB / "headline.tsv",
        headline,
        ["kind", "scope", "set", "universe", "rank", "size", "nes", "es", "p", "fdr_family",
         "signature_logFC", "ci_lo", "ci_hi", "n_boot_ok"],
    )
    (TAB / "public_set_sources.tsv").write_text(
        "set\tsource\n" + "\n".join(f"{k}\t{v}" for k, v in sorted(public_desc.items())) + "\n"
    )

    context = {
        "n_all": len(units),
        "n_q1": sum(r["quartile"] == "Q1" for r in q_full),
        "n_q4": sum(r["quartile"] == "Q4" for r in q_full),
        "n_genes": len(zero_full["genes"]),
        "cldn4_logfc": cldn_fc,
        "cldn4_t": cldn_t,
        "cldn4_p": cldn_p,
        "baseline": baseline,
        "baseline_drop": baseline_drop,
        "baseline_boot": boots.get("baseline"),
        "ifn_max": ifn_max,
        "ifn_max_drop_same": ifn_same_drop,
        "ifn_max_on_drop": ifn_drop_best,
        "ifn_boot": boots.get("ifn_max"),
        "drop_boot": boots.get("ifn_drop_best"),
        "mhc_max": mhc_max,
        "mhc_same_drop": mhc_same,
        "chem_max": chem_max,
        "chem_same_drop": chem_same,
        "sig_max": sig_max,
        "sig_same_drop": sig_same,
        "gene_max": gene_max,
        "gene_max_drop": gene_max_drop,
        "n_ifn_full": sum(1 for r in rows if r["scope"] == "full" and eligible(r, "ifn")),
        "rows": rows,
        "sig_rows": sig_rows,
        "gene_rows": gene_rows,
    }
    # The bootstrap cache stored job tuples and names. boots.get("baseline") works
    # only if the name was written. Confirm keys.
    print("Bootstrap keys", [k if isinstance(k, str) else "job" for k in boots], flush=True)
    write_finding(ROOT / "FINDING.md", context)
    plot_all(context)
    print("IFN max", ifn_max["set"], ifn_max["rank"], ifn_max["universe"], ifn_max["nes"], flush=True)
    print("Signature max", sig_max["set"], sig_max["signature_logFC"], flush=True)
    print("Gene max", gene_max["gene"], gene_max["logFC"], flush=True)
    print("Done", flush=True)


def _boot_summary(out):
    ok = out[np.isfinite(out)]
    if ok.size < 20:
        return {"n_boot_ok": int(ok.size), "nes_boot_median": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan")}
    lo, hi = np.quantile(ok, [0.025, 0.975])
    return {"n_boot_ok": int(ok.size), "nes_boot_median": float(np.median(ok)), "ci_lo": float(lo), "ci_hi": float(hi)}


if __name__ == "__main__":
    main()
