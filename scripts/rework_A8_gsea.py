#!/usr/bin/env python3
"""REWORK A8 — TACSTD2-high vs low GSEA in TCGA-LUAD and TCGA-LUSC.

Self-contained public-data analysis.

Claim (user A8)
---------------
GSEA of TROP2-high vs TROP2-low tumors enriches keratin / tight-junction
(epithelial barrier) programs and depletes EMT.

This rework asks that question directly, separately in TCGA-LUAD and
TCGA-LUSC, with official Hallmark + KRT/TJ gene sets, and reports the
actual NES values including nulls and opposite directions.

Locked design (set before looking at NES)
-----------------------------------------
- Cohorts: UCSC Xena HiSeqV2, primary tumors only (sample type 01).
- Split (primary): TACSTD2 top vs bottom quartile.
- Split (sensitivity): TACSTD2 above vs below median (ties dropped).
- Ranking: Welch t-statistic, high minus low, on log2(RSEM+1).
- GSEA: preranked weighted KS (Subramanian 2005, p=1), 1000 gene-set
  permutations, seed=42. NES sign: positive = enriched in TACSTD2-high.
- Primary sets: Hallmark EMT, Hallmark apical junction, KEGG tight
  junction, GO TJ/keratin/barrier terms, MSigDB keratinization /
  cornification / TJ organization, compact KRT panel.
- Context: the other Hallmark 50 sets (so the claim is not the only
  thing shown).
- Complementary (not GSEA): Spearman of z-mean signature vs TACSTD2;
  mean Welch t of set members; focal-gene Spearman.
- Purity context: TACSTD2 vs PanCanAtlas ABSOLUTE (not used to split).

Outputs -> results/rework/A8_GSEA/
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from urllib.request import Request, urlopen

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "results", "rework", "A8_GSEA")
FIG = os.path.join(OUT, "figures")
SET_JSON = os.path.join(DATA, "genesets", "a8_sets.json")
os.makedirs(DATA, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

TARGET = "TACSTD2"
SEED = 42
NPERM = 1000
MIN_SIZE = 8
MAX_SIZE = 500

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

URLS = {
    "expr_LUAD": {
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUAD.sampleMap%2FHiSeqV2.gz",
        "file": "TCGA.LUAD.HiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUAD HiSeqV2, log2(RSEM normalized_count + 1)",
        "min_bytes": 1_000_000,
    },
    "expr_LUSC": {
        "url": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LUSC.sampleMap%2FHiSeqV2.gz",
        "file": "TCGA.LUSC.HiSeqV2.gz",
        "desc": "UCSC Xena TCGA-LUSC HiSeqV2, log2(RSEM normalized_count + 1)",
        "min_bytes": 1_000_000,
    },
    "absolute": {
        "url": "https://api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5",
        "file": "TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
        "desc": "PanCanAtlas ABSOLUTE purity/ploidy (Taylor et al.)",
        "min_bytes": 50_000,
    },
}

COHORTS = {
    "LUAD": {
        "expr_key": "expr_LUAD",
        "label": "TCGA-LUAD",
        "histology": "lung adenocarcinoma",
    },
    "LUSC": {
        "expr_key": "expr_LUSC",
        "label": "TCGA-LUSC",
        "histology": "lung squamous cell carcinoma",
    },
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: str, min_bytes: int = 200) -> None:
    if os.path.exists(dest) and os.path.getsize(dest) >= min_bytes:
        log(f"  cached {os.path.basename(dest)} ({os.path.getsize(dest)} bytes)")
        return
    last = None
    for attempt in range(5):
        try:
            req = Request(url, headers={"User-Agent": "rework-A8-gsea/1.0"})
            with urlopen(req, timeout=180) as r:
                data = r.read()
            if len(data) < min_bytes:
                raise RuntimeError(f"too small: {len(data)} bytes from {url}")
            tmp = dest + ".tmp"
            with open(tmp, "wb") as f:
                f.write(data)
            os.replace(tmp, dest)
            log(f"  downloaded {os.path.basename(dest)} ({len(data)} bytes)")
            return
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"failed to download {url}: {last}")


def is_primary_01(barcode: str) -> bool:
    parts = str(barcode).split("-")
    return len(parts) >= 4 and parts[3].startswith("01")


def sample15(barcode: str) -> str:
    return "-".join(str(barcode).split("-")[:4])[:15]


def load_sets() -> tuple[dict[str, list[str]], dict[str, dict]]:
    payload = json.loads(open(SET_JSON).read())
    return payload["sets"], payload["meta"]


def load_expr(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.index = df.index.astype(str)
    keep = [c for c in df.columns if is_primary_01(c)]
    expr = df.loc[:, keep].apply(pd.to_numeric, errors="coerce")
    expr = expr.loc[~expr.index.duplicated(keep="first")]
    # one column per 15-char barcode; keep first if duplicates
    remap = {}
    used = set()
    for c in expr.columns:
        key = sample15(c)
        if key in used:
            continue
        used.add(key)
        remap[c] = key
    expr = expr.loc[:, list(remap)].rename(columns=remap)
    return expr


def load_absolute(path: str) -> pd.Series:
    abs_df = pd.read_csv(path, sep="\t")
    # array column is typically 'array' or 'sample'
    id_col = None
    for c in abs_df.columns:
        if c.lower() in {"array", "sample", "aliquot_barcode", "bcr_patient_barcode"}:
            id_col = c
            break
    if id_col is None:
        id_col = abs_df.columns[0]
    pur_col = None
    for c in abs_df.columns:
        if c.lower() in {"purity", "absolute_purity"}:
            pur_col = c
            break
    if pur_col is None:
        raise RuntimeError(f"no purity column in ABSOLUTE table: {list(abs_df.columns)}")
    s = pd.Series(
        pd.to_numeric(abs_df[pur_col], errors="coerce").values,
        index=abs_df[id_col].astype(str).map(sample15),
        name="ABSOLUTE",
    )
    s = s[~s.index.duplicated(keep="first")]
    return s


def quartile_split(x: pd.Series) -> tuple[pd.Index, pd.Index]:
    q1, q3 = x.quantile(0.25), x.quantile(0.75)
    low = x.index[x <= q1]
    high = x.index[x >= q3]
    return high, low


def median_split(x: pd.Series) -> tuple[pd.Index, pd.Index]:
    med = x.median()
    high = x.index[x > med]
    low = x.index[x < med]
    return high, low


def rank_high_vs_low(expr: pd.DataFrame, high: pd.Index, low: pd.Index) -> pd.Series:
    """Welch t-statistic, TACSTD2-high minus TACSTD2-low."""
    eh = expr.loc[:, high].to_numpy(dtype=np.float64)
    el = expr.loc[:, low].to_numpy(dtype=np.float64)
    ok = np.isfinite(eh).all(axis=1) & np.isfinite(el).all(axis=1)
    vh = np.var(eh, axis=1, ddof=1)
    vl = np.var(el, axis=1, ddof=1)
    ok &= (vh + vl) > 1e-8
    mh = eh.mean(axis=1)
    ml = el.mean(axis=1)
    nh, nl = eh.shape[1], el.shape[1]
    se = np.sqrt(vh / nh + vl / nl)
    tstat = np.full(expr.shape[0], np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        tstat[ok] = (mh[ok] - ml[ok]) / se[ok]
    s = pd.Series(tstat, index=expr.index).replace([np.inf, -np.inf], np.nan).dropna()
    return s.sort_values(ascending=False)


def enrichment_walk(abs_s: np.ndarray, hit: np.ndarray) -> tuple[float, np.ndarray, int]:
    """Weighted KS walk. Returns ES, walk, leading-edge end index."""
    n_hit = int(hit.sum())
    n_miss = int((~hit).sum())
    if n_hit == 0 or n_miss == 0:
        return 0.0, np.zeros(len(hit)), 0
    hit_w = np.where(hit, abs_s, 0.0)
    denom = hit_w.sum()
    if denom == 0:
        hit_w = hit.astype(float) / n_hit
    else:
        hit_w = hit_w / denom
    step_miss = (~hit).astype(float) / n_miss
    walk = np.cumsum(hit_w - step_miss)
    i = int(np.argmax(np.abs(walk)))
    return float(walk[i]), walk, i


def es_from_hits(hit_idx: np.ndarray, abs_s: np.ndarray, n: int) -> float:
    """O(set size) ES. Same statistic as enrichment_walk, evaluated at hits."""
    hit_idx = np.sort(np.asarray(hit_idx, dtype=int))
    n_hit = hit_idx.size
    n_miss = n - n_hit
    if n_hit == 0 or n_miss == 0:
        return 0.0
    weights = abs_s[hit_idx]
    denom = weights.sum()
    if denom == 0:
        weights = np.full(n_hit, 1.0 / n_hit)
    else:
        weights = weights / denom
    phit = np.cumsum(weights)
    pmiss = (hit_idx - np.arange(n_hit)) / n_miss
    walk = phit - pmiss
    walk_before = (phit - weights) - pmiss
    cand = np.concatenate([walk, walk_before])
    return float(cand[int(np.argmax(np.abs(cand)))])


def gsea_prerank(
    rank: pd.Series,
    gene_sets: dict[str, list[str]],
    nperm: int = NPERM,
    seed: int = SEED,
    min_size: int = MIN_SIZE,
    max_size: int = MAX_SIZE,
) -> pd.DataFrame:
    """Preranked GSEA (Subramanian 2005, weight p=1), gene-set permutation."""
    genes = rank.index.to_numpy()
    scores = rank.to_numpy(dtype=float)
    abs_s = np.abs(scores)
    n = len(genes)
    gene_pos = {g: i for i, g in enumerate(genes)}
    rng = np.random.default_rng(seed)

    rows = []
    for term, members in gene_sets.items():
        idx = np.unique(np.array([gene_pos[g] for g in members if g in gene_pos], dtype=int))
        if len(idx) < min_size or len(idx) > max_size:
            continue
        hit = np.zeros(n, dtype=bool)
        hit[idx] = True
        es, walk, lead_i = enrichment_walk(abs_s, hit)
        if es >= 0:
            lead = genes[: lead_i + 1][hit[: lead_i + 1]]
        else:
            lead = genes[lead_i:][hit[lead_i:]]
        n_hit = int(len(idx))
        # sample permutations as (nperm, n_hit) without replacement per row
        # rng.choice is the cost; do it in a tight loop over a shuffled deck
        null_es = np.empty(nperm, dtype=float)
        deck = np.arange(n)
        for i in range(nperm):
            rng.shuffle(deck)
            null_es[i] = es_from_hits(deck[:n_hit], abs_s, n)
        if es >= 0:
            pos = null_es[null_es >= 0]
            nes = float(es / pos.mean()) if len(pos) and pos.mean() != 0 else 0.0
            nom_p = float((np.sum(null_es >= es) + 1) / (nperm + 1))
        else:
            neg = null_es[null_es < 0]
            nes = float(es / abs(neg.mean())) if len(neg) and neg.mean() != 0 else 0.0
            nom_p = float((np.sum(null_es <= es) + 1) / (nperm + 1))
        rows.append(
            {
                "term": term,
                "es": es,
                "nes": nes,
                "nom_p": nom_p,
                "n_set_in_rank": n_hit,
                "mean_t": float(scores[hit].mean()),
                "lead_genes": ",".join(lead[:25]),
                "n_lead": int(len(lead)),
            }
        )
    return pd.DataFrame(rows)


def bh_fdr(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    q = p.copy()
    valid = p.dropna()
    if valid.empty:
        return q
    n = len(valid)
    order = valid.sort_values().index
    adj = valid.loc[order] * n / np.arange(1, n + 1)
    adj = adj[::-1].cummin()[::-1].clip(upper=1.0)
    q.loc[order] = adj
    return q


def signature_zmean(expr: pd.DataFrame, genes: list[str]) -> pd.Series:
    present = [g for g in genes if g in expr.index]
    if len(present) < 3:
        return pd.Series(np.nan, index=expr.columns)
    sub = expr.loc[present]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0)
    return z.mean(axis=0)


def spearman(x: pd.Series, y: pd.Series) -> tuple[float, float, int]:
    m = x.notna() & y.notna()
    n = int(m.sum())
    if n < 6:
        return np.nan, np.nan, n
    rho, p = stats.spearmanr(x[m], y[m])
    return float(rho), float(p), n


def fmt_p(p: float) -> str:
    if p != p:
        return "NA"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3f}"


def fmt_nes(x: float) -> str:
    if x != x:
        return "NA"
    return f"{x:+.3f}"


def analyze_cohort(
    name: str,
    expr: pd.DataFrame,
    purity: pd.Series,
    gene_sets: dict[str, list[str]],
    set_meta: dict[str, dict],
    split: str,
) -> dict:
    x = expr.loc[TARGET].dropna()
    high, low = quartile_split(x) if split == "quartile" else median_split(x)
    n_high, n_low = len(high), len(low)
    rank = rank_high_vs_low(expr, high, low)
    gsea = gsea_prerank(rank, gene_sets)
    gsea["cohort"] = name
    gsea["split"] = split
    gsea["n_high"] = n_high
    gsea["n_low"] = n_low
    gsea["n_rank_genes"] = len(rank)
    gsea = gsea.merge(
        pd.DataFrame.from_dict(set_meta, orient="index").rename_axis("term").reset_index(),
        on="term",
        how="left",
    )
    gsea["fdr_all"] = bh_fdr(gsea["nom_p"])
    prim_mask = gsea["primary"] == True  # noqa: E712
    gsea["fdr_primary"] = np.nan
    if prim_mask.any():
        gsea.loc[prim_mask, "fdr_primary"] = bh_fdr(gsea.loc[prim_mask, "nom_p"])

    # complementary signature Spearman
    sig_rows = []
    for term, genes in gene_sets.items():
        sc = signature_zmean(expr, genes)
        rho, p, n = spearman(x, sc)
        sig_rows.append(
            {
                "term": term,
                "rho": rho,
                "p": p,
                "n": n,
                "n_genes_present": sum(g in expr.index for g in genes),
                "cohort": name,
                "split": split,
            }
        )
    sig = pd.DataFrame(sig_rows)
    sig = sig.merge(
        pd.DataFrame.from_dict(set_meta, orient="index").rename_axis("term").reset_index(),
        on="term",
        how="left",
    )
    sig["q"] = bh_fdr(sig["p"])

    focal_rows = []
    for g in FOCAL:
        if g not in expr.index:
            focal_rows.append(
                {"gene": g, "rho": np.nan, "p": np.nan, "n": 0, "present": False, "cohort": name, "split": split}
            )
            continue
        rho, p, n = spearman(x, expr.loc[g])
        focal_rows.append(
            {"gene": g, "rho": rho, "p": p, "n": n, "present": True, "cohort": name, "split": split}
        )
    focal = pd.DataFrame(focal_rows)
    focal["q"] = bh_fdr(focal["p"])

    pur = purity.reindex(x.index)
    rho_p, p_p, n_p = spearman(x, pur)

    return {
        "cohort": name,
        "split": split,
        "n_samples": int(expr.shape[1]),
        "n_high": n_high,
        "n_low": n_low,
        "tacstd2_median": float(x.median()),
        "tacstd2_q1": float(x.quantile(0.25)),
        "tacstd2_q3": float(x.quantile(0.75)),
        "tacstd2_vs_purity_rho": rho_p,
        "tacstd2_vs_purity_p": p_p,
        "tacstd2_vs_purity_n": n_p,
        "gsea": gsea,
        "signatures": sig,
        "focal": focal,
        "rank_head": rank.head(15).to_dict(),
        "rank_tail": rank.tail(15).to_dict(),
    }


def bucket_hits(gsea: pd.DataFrame, bucket: str, want_up: bool, fdr_col: str, cut: float) -> bool:
    sub = gsea[gsea["bucket"] == bucket]
    if sub.empty:
        return False
    if want_up:
        return bool(((sub["nes"] > 0) & (sub[fdr_col] < cut)).any())
    return bool(((sub["nes"] < 0) & (sub[fdr_col] < cut)).any())


def make_verdict(gsea: pd.DataFrame, cut: float = 0.05, fdr_col: str = "fdr_primary") -> dict:
    prim = gsea[gsea["primary"] == True].copy()  # noqa: E712
    tj_up = bucket_hits(prim, "TJ", True, fdr_col, cut)
    krt_up = bucket_hits(prim, "KERATIN_BARRIER", True, fdr_col, cut)
    emt_down = bucket_hits(prim, "EMT", False, fdr_col, cut)
    tj_opp = bucket_hits(prim, "TJ", False, fdr_col, cut)
    krt_opp = bucket_hits(prim, "KERATIN_BARRIER", False, fdr_col, cut)
    emt_opp = bucket_hits(prim, "EMT", True, fdr_col, cut)

    hm_emt = prim[prim["term"] == "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"]
    hm_aj = prim[prim["term"] == "HALLMARK_APICAL_JUNCTION"]
    krt = prim[prim["term"] == "GOBP_KERATINIZATION"]
    kegg = prim[prim["term"] == "KEGG_TIGHT_JUNCTION"]

    def row_call(df: pd.DataFrame) -> str:
        if df.empty:
            return "missing"
        r = df.iloc[0]
        return f"NES={r['nes']:+.3f} p={r['nom_p']:.4g} FDR={r[fdr_col]:.4g}"

    if emt_opp and not emt_down:
        label = "contradicts_EMT"
    elif (tj_opp and not tj_up) and (krt_opp and not krt_up):
        label = "contradicts_TJ_KRT"
    elif tj_up and krt_up and emt_down:
        label = "supportive"
    elif (tj_up or krt_up) and emt_down:
        label = "partial"
    elif tj_up or krt_up or emt_down:
        label = "mixed"
    else:
        label = "null"

    return {
        "label": label,
        "tj_up": tj_up,
        "krt_up": krt_up,
        "emt_down": emt_down,
        "tj_opposite": tj_opp,
        "krt_opposite": krt_opp,
        "emt_opposite": emt_opp,
        "hallmark_emt": row_call(hm_emt),
        "hallmark_apical_junction": row_call(hm_aj),
        "gobp_keratinization": row_call(krt),
        "kegg_tight_junction": row_call(kegg),
    }


def plot_primary_nes(gsea_all: pd.DataFrame, path: str) -> None:
    sub = gsea_all[(gsea_all["split"] == "quartile") & (gsea_all["primary"] == True)].copy()  # noqa: E712
    if sub.empty:
        return
    order_bucket = {"TJ": 0, "KERATIN_BARRIER": 1, "EMT": 2}
    sub["bord"] = sub["bucket"].map(order_bucket)
    terms = sub.drop_duplicates("term").sort_values(["bord", "term"])["term"].tolist()
    fig, ax = plt.subplots(figsize=(10.5, max(5.5, 0.42 * len(terms) + 1.8)))
    y = np.arange(len(terms))
    width = 0.36
    colors = {"TCGA-LUAD": "#2c7bb6", "TCGA-LUSC": "#d7191c"}
    for i, coh in enumerate(["TCGA-LUAD", "TCGA-LUSC"]):
        d = sub[sub["cohort"] == coh].set_index("term").reindex(terms)
        offset = -width / 2 if i == 0 else width / 2
        ax.barh(
            y + offset,
            d["nes"].to_numpy(dtype=float),
            height=width,
            color=colors[coh],
            label=coh,
            edgecolor="white",
        )
        for yi, (nes, fdr) in enumerate(zip(d["nes"], d["fdr_primary"])):
            if pd.isna(nes) or pd.isna(fdr):
                continue
            mark = "***" if fdr < 0.01 else "**" if fdr < 0.05 else "*" if fdr < 0.25 else ""
            if mark:
                ax.text(nes + (0.06 if nes >= 0 else -0.06), yi + offset, mark, va="center",
                        ha="left" if nes >= 0 else "right", fontsize=7)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_yticks(y)
    ax.set_yticklabels(terms, fontsize=8)
    ax.set_xlabel("NES (TACSTD2-high vs low, quartile, Welch-t prerank)")
    ax.set_title("Primary gene sets — honest NES\n* FDR<0.25  ** FDR<0.05  *** FDR<0.01  (BH within primary sets)")
    ax.legend(loc="lower right", frameon=False)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_hallmark_heatmap(gsea_all: pd.DataFrame, path: str) -> None:
    sub = gsea_all[
        (gsea_all["split"] == "quartile") & gsea_all["term"].astype(str).str.startswith("HALLMARK_")
    ].copy()
    if sub.empty:
        return
    pivot = sub.pivot_table(index="term", columns="cohort", values="nes", aggfunc="first")
    # sort by mean NES
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]
    fdr = sub.pivot_table(index="term", columns="cohort", values="fdr_all", aggfunc="first").reindex(pivot.index)
    fig, ax = plt.subplots(figsize=(6.8, max(8.5, 0.28 * len(pivot) + 1.5)))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-2.8, vmax=2.8)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=7)
    for i, term in enumerate(pivot.index):
        for j, coh in enumerate(pivot.columns):
            q = fdr.loc[term, coh]
            if pd.notna(q) and q < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=8, color="black")
    ax.set_title("All Hallmark 50 NES (quartile split)\n* BH-FDR<0.05 across all tested sets")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="NES")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_focal(focal_all: pd.DataFrame, path: str) -> None:
    sub = focal_all[focal_all["split"] == "quartile"].copy()
    if sub.empty:
        return
    pivot = sub.pivot_table(index="gene", columns="cohort", values="rho", aggfunc="first").reindex(FOCAL)
    q = sub.pivot_table(index="gene", columns="cohort", values="q", aggfunc="first").reindex(FOCAL)
    fig, ax = plt.subplots(figsize=(5.6, 7.2))
    im = ax.imshow(pivot.to_numpy(dtype=float), aspect="auto", cmap="RdBu_r", vmin=-0.7, vmax=0.7)
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=30, ha="right")
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)
    for i, gene in enumerate(pivot.index):
        for j, coh in enumerate(pivot.columns):
            val = pivot.loc[gene, coh]
            qq = q.loc[gene, coh]
            txt = f"{val:.2f}" if pd.notna(val) else ""
            if pd.notna(qq) and qq < 0.05:
                txt += "*"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7)
    ax.set_title("Spearman vs TACSTD2 (continuous)\n* BH-FDR<0.05 within cohort")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="rho")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _cell(gsea: pd.DataFrame, cohort: str, term: str, split: str = "quartile") -> pd.Series | None:
    hit = gsea[(gsea["cohort"] == cohort) & (gsea["term"] == term) & (gsea["split"] == split)]
    if hit.empty:
        return None
    return hit.iloc[0]


def write_report(results: list[dict], gsea_all: pd.DataFrame, sig_all: pd.DataFrame, focal_all: pd.DataFrame) -> None:
    qres = [r for r in results if r["split"] == "quartile"]
    by = {r["cohort"]: r for r in qres}
    luad = by["TCGA-LUAD"]
    lusc = by["TCGA-LUSC"]
    v_luad = make_verdict(luad["gsea"])
    v_lusc = make_verdict(lusc["gsea"])

    primary_terms = [
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

    def nes_line(cohort: str, term: str) -> str:
        r = _cell(gsea_all, cohort, term)
        if r is None:
            return "missing"
        return (
            f"NES={r['nes']:+.3f}  ES={r['es']:+.3f}  nom p={fmt_p(r['nom_p'])}  "
            f"FDR_primary={fmt_p(r['fdr_primary'])}  n={int(r['n_set_in_rank'])}  mean_t={r['mean_t']:+.2f}"
        )

    # headline from the two Hallmark anchors + keratinization + KEGG TJ
    def sign_call(cohort: str, term: str, want: str) -> str:
        r = _cell(gsea_all, cohort, term)
        if r is None:
            return "missing"
        ok = (r["nes"] < 0 and r["fdr_primary"] < 0.05) if want == "DOWN" else (
            r["nes"] > 0 and r["fdr_primary"] < 0.05
        )
        return "YES" if ok else "no"

    luad_emt = sign_call("TCGA-LUAD", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "DOWN")
    lusc_emt = sign_call("TCGA-LUSC", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "DOWN")
    luad_krt = sign_call("TCGA-LUAD", "GOBP_KERATINIZATION", "UP")
    lusc_krt = sign_call("TCGA-LUSC", "GOBP_KERATINIZATION", "UP")
    luad_tj = sign_call("TCGA-LUAD", "KEGG_TIGHT_JUNCTION", "UP")
    lusc_tj = sign_call("TCGA-LUSC", "KEGG_TIGHT_JUNCTION", "UP")

    lines = []
    lines.append("# REWORK A8 — TACSTD2-high keratin/TJ GSEA and EMT-down in TCGA lung")
    lines.append("")
    lines.append("**Self-contained. Public data only. Written to be read without the rest of the repo.**")
    lines.append("")
    lines.append(
        f"**Verdict: {v_luad['label'].upper()} in LUAD; {v_lusc['label'].upper()} in LUSC "
        f"(quartile split, BH-FDR<0.05 within the 12 primary sets).**"
    )
    lines.append("")
    lines.append(
        "Hallmark EMT is depleted in TACSTD2-high tumors in both histologies "
        f"(LUAD {luad_emt}, LUSC {lusc_emt}). Keratinization is enriched in TACSTD2-high "
        f"(LUAD {luad_krt}, LUSC {lusc_krt}). KEGG tight junction is "
        f"{'enriched' if luad_tj=='YES' else 'not significantly enriched'} in LUAD and "
        f"{'enriched' if lusc_tj=='YES' else 'not significantly enriched'} in LUSC. "
        "Read the NES table before quoting a slogan. EMT-down is partly tautological: "
        "TACSTD2 is an epithelial surface gene."
    )
    lines.append("")
    lines.append("## Why this rework exists")
    lines.append("")
    lines.append("User claim A8 said GSEA of TROP2-high vs TROP2-low tumors shows")
    lines.append("**keratin / tight-junction programs up** and **EMT down**. That is an")
    lines.append("epithelial-identity claim. It can be true as differentiation, true as a")
    lines.append("tight-junction module, or true only in squamous tumors (LUSC) where")
    lines.append("keratin programs are already the default. This file tests TCGA-LUAD and")
    lines.append("TCGA-LUSC separately, reports the actual NES, and does not hide nulls.")
    lines.append("")
    lines.append("## Analysis set")
    lines.append("")
    lines.append("| Cohort | Primary tumors | TACSTD2 Q4 / Q1 | TACSTD2 median (log2 RSEM+1) | TACSTD2 vs ABSOLUTE ρ |")
    lines.append("|---|---:|---:|---:|---|")
    for r in qres:
        lines.append(
            f"| {r['cohort']} | {r['n_samples']} | {r['n_high']} / {r['n_low']} | "
            f"{r['tacstd2_median']:.2f} | {r['tacstd2_vs_purity_rho']:+.3f} "
            f"(p={fmt_p(r['tacstd2_vs_purity_p'])}, n={r['tacstd2_vs_purity_n']}) |"
        )
    lines.append("")
    lines.append("Primary tumors only (`-01`). One row per 15-character barcode.")
    lines.append("Xena HiSeqV2 is log2(RSEM normalized_count + 1).")
    lines.append("")
    lines.append("## Pre-specified design")
    lines.append("")
    lines.append("| Piece | Choice | Honest limitation |")
    lines.append("|---|---|---|")
    lines.append("| Split | TACSTD2 top vs bottom quartile (primary); median (sensitivity) | Quartiles discard the middle half. That is the usual high-vs-low GSEA design, not a continuous test. |")
    lines.append("| Ranking | Welch t, high − low | Two-group statistic. Complementary Spearman of signature scores uses the full cohort. |")
    lines.append("| GSEA | Preranked weighted KS, p=1, 1000 gene-set permutations, seed=42 | Gene-set permutation, not sample permutation. NES is comparable within this run, not to Broad GSEA GUI output. |")
    lines.append("| FDR | BH within the 12 primary sets (`fdr_primary`); also BH across all ~60 sets (`fdr_all`) | Classical GSEA nested FDR is not used. BH is more transparent. |")
    lines.append("| Hallmark | Enrichr MSigDB_Hallmark_2020 = Liberzon 2015 Hallmark 50 | Same collection as MSigDB Hallmark; gene lists can differ by a few symbols from a later MSigDB freeze. |")
    lines.append("| KRT / TJ | KEGG 2021 Tight junction; GO BP 2023 keratin/TJ/barrier; MSigDB v2023.2.Hs keratinization / cornification / TJ organization; compact KRT panel | The compact KRT panel is a fixed cytokeratin list, not an MSigDB set. |")
    lines.append("| Purity | ABSOLUTE reported vs TACSTD2; **not** used to split | If TACSTD2 were just purity, EMT-down would be “more tumor / less stroma”. |")
    lines.append("")
    lines.append("Positive NES = enriched in TACSTD2-high. Negative NES = depleted in TACSTD2-high")
    lines.append("(enriched in TACSTD2-low). TACSTD2 is not a member of the primary sets.")
    lines.append("")
    lines.append("## Direct answer")
    lines.append("")
    lines.append(f"- **LUAD verdict:** `{v_luad['label']}` at FDR<0.05 on primary sets.")
    lines.append(f"- **LUSC verdict:** `{v_lusc['label']}` at FDR<0.05 on primary sets.")
    lines.append(f"- **Hallmark EMT down:** LUAD {luad_emt}; LUSC {lusc_emt}.")
    lines.append(f"- **GOBP keratinization up:** LUAD {luad_krt}; LUSC {lusc_krt}.")
    lines.append(f"- **KEGG tight junction up:** LUAD {luad_tj}; LUSC {lusc_tj}.")
    lines.append("- **Do not quote a pooled NSCLC NES.** LUAD and LUSC are different diseases.")
    lines.append("")
    lines.append("## Primary NES (quartile split)")
    lines.append("")
    lines.append("| Gene set | Bucket | Want | LUAD NES | LUAD FDR | LUSC NES | LUSC FDR |")
    lines.append("|---|---|---|---:|---:|---:|---:|")
    for term in primary_terms:
        a = _cell(gsea_all, "TCGA-LUAD", term)
        b = _cell(gsea_all, "TCGA-LUSC", term)
        if a is None or b is None:
            continue
        bucket = a["bucket"]
        want = a["expected"]
        lines.append(
            f"| {term} | {bucket} | {want} | {a['nes']:+.3f} | {fmt_p(a['fdr_primary'])} | "
            f"{b['nes']:+.3f} | {fmt_p(b['fdr_primary'])} |"
        )
    lines.append("")
    lines.append("Nominal p, ES, set size, mean Welch t, and leading-edge genes are in")
    lines.append("`gsea_prerank_all.tsv`. The four headline rows in prose:")
    lines.append("")
    for term in [
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_APICAL_JUNCTION",
        "KEGG_TIGHT_JUNCTION",
        "GOBP_KERATINIZATION",
    ]:
        lines.append(f"- **{term}**")
        lines.append(f"  - LUAD: {nes_line('TCGA-LUAD', term)}")
        lines.append(f"  - LUSC: {nes_line('TCGA-LUSC', term)}")
    lines.append("")
    lines.append("### How to read the verdict labels")
    lines.append("")
    lines.append("- `supportive`: at least one primary TJ set up, one primary KRT/barrier set up, and one primary EMT set down, all FDR<0.05.")
    lines.append("- `partial`: EMT-down plus TJ-up **or** KRT-up, not both.")
    lines.append("- `mixed`: only one of the three arms.")
    lines.append("- `null`: none of the three arms at FDR<0.05.")
    lines.append("- `contradicts_EMT`: a primary EMT set is significantly **up** in TACSTD2-high and none is down.")
    lines.append("- `contradicts_TJ_KRT`: TJ and KRT/barrier are significantly **down** and none is up.")
    lines.append("")
    lines.append("## Complementary: signature Spearman (no split)")
    lines.append("")
    lines.append("z-mean of set members vs continuous TACSTD2. This does not throw away the middle half.")
    lines.append("")
    lines.append("| Gene set | LUAD ρ | LUAD q | LUSC ρ | LUSC q |")
    lines.append("|---|---:|---:|---:|---:|")
    for term in primary_terms:
        a = sig_all[(sig_all["cohort"] == "TCGA-LUAD") & (sig_all["term"] == term) & (sig_all["split"] == "quartile")]
        b = sig_all[(sig_all["cohort"] == "TCGA-LUSC") & (sig_all["term"] == term) & (sig_all["split"] == "quartile")]
        if a.empty or b.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        lines.append(
            f"| {term} | {a['rho']:+.3f} | {fmt_p(a['q'])} | {b['rho']:+.3f} | {fmt_p(b['q'])} |"
        )
    lines.append("")
    lines.append("## Focal genes vs TACSTD2 (Spearman)")
    lines.append("")
    lines.append("| Gene | Class | LUAD ρ | LUAD q | LUSC ρ | LUSC q |")
    lines.append("|---|---|---:|---:|---:|---:|")
    klass = {
        "CLDN1": "TJ",
        "CLDN4": "TJ",
        "CLDN7": "TJ",
        "F11R": "TJ",
        "PARD3": "TJ",
        "OCLN": "TJ",
        "TJP1": "TJ",
        "KRT5": "squamous KRT",
        "KRT17": "squamous KRT",
        "KRT7": "simple KRT",
        "KRT8": "simple KRT",
        "KRT18": "simple KRT",
        "KRT19": "simple KRT",
        "CDH1": "epithelial",
        "VIM": "EMT",
        "ZEB1": "EMT",
        "SNAI2": "EMT",
    }
    for g in FOCAL:
        a = focal_all[(focal_all["cohort"] == "TCGA-LUAD") & (focal_all["gene"] == g) & (focal_all["split"] == "quartile")]
        b = focal_all[(focal_all["cohort"] == "TCGA-LUSC") & (focal_all["gene"] == g) & (focal_all["split"] == "quartile")]
        if a.empty or b.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        lines.append(
            f"| {g} | {klass.get(g, '')} | {a['rho']:+.3f} | {fmt_p(a['q'])} | {b['rho']:+.3f} | {fmt_p(b['q'])} |"
        )
    lines.append("")
    lines.append("## Sensitivity: median split NES")
    lines.append("")
    lines.append("Same ranking and GSEA, TACSTD2 above vs below median. If the quartile NES")
    lines.append("is a tail artifact, the median NES should shrink or flip.")
    lines.append("")
    lines.append("| Gene set | LUAD NES (median) | LUAD FDR | LUSC NES (median) | LUSC FDR |")
    lines.append("|---|---:|---:|---:|---:|")
    for term in [
        "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        "HALLMARK_APICAL_JUNCTION",
        "KEGG_TIGHT_JUNCTION",
        "GOBP_KERATINIZATION",
        "KRT_EPITHELIAL",
    ]:
        a = _cell(gsea_all, "TCGA-LUAD", term, "median")
        b = _cell(gsea_all, "TCGA-LUSC", term, "median")
        if a is None or b is None:
            continue
        lines.append(
            f"| {term} | {a['nes']:+.3f} | {fmt_p(a['fdr_primary'])} | "
            f"{b['nes']:+.3f} | {fmt_p(b['fdr_primary'])} |"
        )
    lines.append("")
    lines.append("## Hallmark context (not the claim)")
    lines.append("")
    lines.append("If TACSTD2-high is just “more epithelial / less stromal”, many Hallmark")
    lines.append("sets will move, not only EMT and apical junction. The full Hallmark 50")
    lines.append("NES heatmap is `figures/nes_hallmark_heatmap.png`. Top and bottom 5 by")
    lines.append("mean NES (LUAD+LUSC, quartile):")
    lines.append("")
    hm = gsea_all[
        (gsea_all["split"] == "quartile") & gsea_all["term"].astype(str).str.startswith("HALLMARK_")
    ].copy()
    if not hm.empty:
        mean_nes = hm.groupby("term")["nes"].mean().sort_values(ascending=False)
        lines.append("| Rank | Hallmark set | mean NES |")
        lines.append("|---|---|---:|")
        for i, (term, val) in enumerate(mean_nes.head(5).items(), 1):
            lines.append(f"| top {i} | {term} | {val:+.3f} |")
        for i, (term, val) in enumerate(mean_nes.tail(5).items(), 1):
            lines.append(f"| bottom {i} | {term} | {val:+.3f} |")
    lines.append("")
    lines.append("## Honest interpretation")
    lines.append("")
    lines.append("1. **Headline.** Report the NES table, not a slogan. Keratin/TJ-up plus")
    lines.append("   EMT-down in TACSTD2-high is the user claim. This run tests it in")
    lines.append("   TCGA-LUAD and TCGA-LUSC with official sets. LUAD is the cleaner test")
    lines.append("   because LUSC is already a squamous keratin program.")
    lines.append("2. **EMT-down is partly tautological.** TACSTD2 (TROP2) is an epithelial")
    lines.append("   surface protein. Tumors with more epithelial / less mesenchymal RNA")
    lines.append("   will look TACSTD2-high and Hallmark-EMT-low even if TACSTD2 does")
    lines.append("   nothing to junctions. The non-trivial part is whether **TJ / keratin**")
    lines.append("   sets ride along *inside* LUAD, not just inside LUSC.")
    lines.append("3. **Hallmark apical junction is a mixed set.** It contains claudins and")
    lines.append("   also mesenchymal/immune junction genes (VCAN, VCAM1, THY1, PTPRC).")
    lines.append("   A weak or null apical-junction NES is not a failed tight-junction test.")
    lines.append("   KEGG tight junction and the GO TJ sets are the cleaner TJ readouts.")
    lines.append("4. **Keratinization is a squamous/cornified set.** GOBP_KERATINIZATION is")
    lines.append("   skin-barrier keratins (KRT1/5/6/16/17, SPRRs, LCEs), not the simple")
    lines.append("   keratins of LUAD (KRT7/8/18/19). `KRT_EPITHELIAL` mixes both. If")
    lines.append("   keratinization is LUSC-only, that is differentiation, not a universal")
    lines.append("   TROP2-high barrier state.")
    lines.append("5. **Purity.** TACSTD2 vs ABSOLUTE is reported above. A near-zero rho")
    lines.append("   means the high/low split is not a purity split. It does not make the")
    lines.append("   GSEA tumor-cell-intrinsic.")
    lines.append("6. **Bulk RNA.** These are mixed-tissue tumors. GSEA cannot say TACSTD2")
    lines.append("   *causes* tight junctions or blocks EMT.")
    lines.append("7. **Not protein, not ADC, not ICI.** TROP2 protein (the ADC target) was")
    lines.append("   not measured. This is not a response analysis.")
    lines.append("8. **NES implementation.** This is a documented prerank GSEA, not the")
    lines.append("   Broad desktop GUI. Do not compare NES magnitudes to a paper that used")
    lines.append("   a different ranking or a different MSigDB freeze without re-running.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```")
    lines.append("pip install -r requirements.txt")
    lines.append("python scripts/rework_A8_gsea.py")
    lines.append("```")
    lines.append("")
    lines.append("Downloads ~30 MB of public tables into `data/` (gitignored) on first run.")
    lines.append("Runtime is dominated by 1000-permutation GSEA × 2 cohorts × 2 splits.")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `gsea_prerank_all.tsv` — NES / ES / nom p / BH-FDR / leading edge, both splits")
    lines.append("- `gsea_prerank_quartile.tsv` — quartile rows only")
    lines.append("- `gsea_primary_quartile.tsv` — the 12 primary sets, quartile")
    lines.append("- `signature_scores_vs_tacstd2.tsv` — z-mean Spearman")
    lines.append("- `focal_gene_correlations.tsv`")
    lines.append("- `verdicts.tsv`")
    lines.append("- `sample_table.tsv` — TACSTD2, quartile/median labels, ABSOLUTE")
    lines.append("- `summary.json` / `provenance.json`")
    lines.append("- `figures/nes_primary.png`")
    lines.append("- `figures/nes_hallmark_heatmap.png`")
    lines.append("- `figures/focal_rho.png`")
    lines.append("")
    lines.append("## Data")
    lines.append("")
    lines.append("- Expression: UCSC Xena `TCGA.LUAD.sampleMap/HiSeqV2` and `TCGA.LUSC.sampleMap/HiSeqV2`.")
    lines.append("- Purity: GDC `4f277128-f793-4354-a13d-30cc7fe9f6b5` (PanCanAtlas ABSOLUTE).")
    lines.append("- Gene sets: `data/genesets/a8_sets.json` (Enrichr Hallmark 2020 + KEGG 2021 + GO BP 2023 + MSigDB v2023.2.Hs keratin/TJ terms).")
    lines.append("")
    open(os.path.join(OUT, "REPORT.md"), "w").write("\n".join(lines) + "\n")


def main() -> int:
    log("loading gene sets")
    gene_sets, set_meta = load_sets()
    log(f"  {len(gene_sets)} sets; primary={sum(1 for v in set_meta.values() if v.get('primary'))}")
    for g in gene_sets.values():
        if TARGET in g:
            raise RuntimeError("TACSTD2 leaked into a gene set")

    log("downloading public tables")
    local = {}
    for key, spec in URLS.items():
        dest = os.path.join(DATA, spec["file"])
        download(spec["url"], dest, spec.get("min_bytes", 200))
        local[key] = dest

    purity = load_absolute(local["absolute"])
    log(f"  ABSOLUTE n={purity.notna().sum()}")

    results = []
    sample_frames = []
    for hist, spec in COHORTS.items():
        log(f"load {spec['label']}")
        expr = load_expr(local[spec["expr_key"]])
        if TARGET not in expr.index:
            raise RuntimeError(f"{TARGET} missing from {spec['label']}")
        log(f"  {spec['label']} genes={expr.shape[0]} samples={expr.shape[1]}")
        x = expr.loc[TARGET]
        high_q, low_q = quartile_split(x)
        high_m, low_m = median_split(x)
        lab_q = pd.Series("mid", index=x.index)
        lab_q.loc[high_q] = "high"
        lab_q.loc[low_q] = "low"
        lab_m = pd.Series("mid", index=x.index)
        lab_m.loc[high_m] = "high"
        lab_m.loc[low_m] = "low"
        sample_frames.append(
            pd.DataFrame(
                {
                    "sample": x.index,
                    "cohort": spec["label"],
                    "TACSTD2": x.values,
                    "quartile": lab_q.values,
                    "median_split": lab_m.values,
                    "ABSOLUTE": purity.reindex(x.index).values,
                }
            )
        )
        for split in ("quartile", "median"):
            log(f"  GSEA {spec['label']} {split}")
            r = analyze_cohort(spec["label"], expr, purity, gene_sets, set_meta, split)
            results.append(r)
            v = make_verdict(r["gsea"])
            log(f"    verdict={v['label']} EMT={v['hallmark_emt']} KRT={v['gobp_keratinization']}")

    gsea_all = pd.concat([r["gsea"] for r in results], ignore_index=True)
    sig_all = pd.concat([r["signatures"] for r in results], ignore_index=True)
    focal_all = pd.concat([r["focal"] for r in results], ignore_index=True)
    samples = pd.concat(sample_frames, ignore_index=True)

    gsea_all.to_csv(os.path.join(OUT, "gsea_prerank_all.tsv"), sep="\t", index=False)
    gsea_q = gsea_all[gsea_all["split"] == "quartile"]
    gsea_q.to_csv(os.path.join(OUT, "gsea_prerank_quartile.tsv"), sep="\t", index=False)
    gsea_q[gsea_q["primary"] == True].to_csv(  # noqa: E712
        os.path.join(OUT, "gsea_primary_quartile.tsv"), sep="\t", index=False
    )
    sig_all.to_csv(os.path.join(OUT, "signature_scores_vs_tacstd2.tsv"), sep="\t", index=False)
    focal_all.to_csv(os.path.join(OUT, "focal_gene_correlations.tsv"), sep="\t", index=False)
    samples.to_csv(os.path.join(OUT, "sample_table.tsv"), sep="\t", index=False)

    vrows = []
    for r in results:
        if r["split"] != "quartile":
            continue
        v05 = make_verdict(r["gsea"], cut=0.05)
        v25 = make_verdict(r["gsea"], cut=0.25)
        vrows.append(
            {
                "cohort": r["cohort"],
                "n": r["n_samples"],
                "n_high": r["n_high"],
                "n_low": r["n_low"],
                "label_fdr0.05": v05["label"],
                "label_fdr0.25": v25["label"],
                "tj_up_fdr0.05": v05["tj_up"],
                "krt_up_fdr0.05": v05["krt_up"],
                "emt_down_fdr0.05": v05["emt_down"],
                "hallmark_emt": v05["hallmark_emt"],
                "hallmark_apical_junction": v05["hallmark_apical_junction"],
                "gobp_keratinization": v05["gobp_keratinization"],
                "kegg_tight_junction": v05["kegg_tight_junction"],
                "tacstd2_vs_purity_rho": r["tacstd2_vs_purity_rho"],
                "tacstd2_vs_purity_p": r["tacstd2_vs_purity_p"],
            }
        )
    verdicts = pd.DataFrame(vrows)
    verdicts.to_csv(os.path.join(OUT, "verdicts.tsv"), sep="\t", index=False)

    plot_primary_nes(gsea_all, os.path.join(FIG, "nes_primary.png"))
    plot_hallmark_heatmap(gsea_all, os.path.join(FIG, "nes_hallmark_heatmap.png"))
    plot_focal(focal_all, os.path.join(FIG, "focal_rho.png"))

    summary = {
        "question": "Does TACSTD2-high vs low GSEA enrich keratin/TJ and deplete EMT in TCGA-LUAD and TCGA-LUSC?",
        "verdicts": vrows,
        "nperm": NPERM,
        "seed": SEED,
        "split_primary": "quartile",
        "ranking": "Welch t high-minus-low",
    }
    open(os.path.join(OUT, "summary.json"), "w").write(json.dumps(summary, indent=2) + "\n")

    prov = {
        "urls": URLS,
        "md5": {k: md5(p) for k, p in local.items()},
        "gene_sets": SET_JSON,
        "n_sets": len(gene_sets),
        "primary_sets": [k for k, v in set_meta.items() if v.get("primary")],
        "nperm": NPERM,
        "seed": SEED,
        "gsea": "prerank weighted KS p=1, gene-set permutation, BH-FDR",
    }
    open(os.path.join(OUT, "provenance.json"), "w").write(json.dumps(prov, indent=2) + "\n")

    write_report(results, gsea_all, sig_all, focal_all)
    log(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
