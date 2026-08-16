#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_gse334497.py
================
Public mouse Trop2 (Tacstd2) knockout RNA-seq (verified GEO accession GSE334497).

Title  : TROP2/claudin program mediates immune exclusion to impede checkpoint
         blockade in breast cancer
Model  : 4T1 mouse mammary tumors, Trop2 KO vs WT, n=5 vs 5, in vivo
Suppl  : GSE334497_normalized_counts.csv.gz  (~1.2 MB; public FTP)
IDs    : mouse Ensembl gene IDs, mapped to MGI symbols via mygene.info

Sample mapping (from GEO series matrix !Sample_description / !Sample_title):
    KO : KO162, KO164, KO165, KO172, RESUB-KO163R
    WT : control170, RESUB-171R, RESUB-170R, RESUB-169R, RESUB-168R

CAVEAT: the supplement is *normalized* (non-integer) counts, not raw counts, so
this is NOT a DESeq2/edgeR run. Welch t-test on log2(norm+1) is exploratory.
n=5 is better than the n=2-3 pitfall but still a processed matrix.

Opposite-gene question (mouse symbols):
    Tacstd2 KO → does Cldn4 go UP?
Gene-set GSEA uses human Hallmark/KEGG after toupper(MGI) — documented
limitation; see playbook.md §9. Fine for these well-conserved sets.
"""
from __future__ import annotations
import gzip
import json
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
DATA = HERE / "data"
RES = HERE / "results"
GMT = HERE.parent / "gene_sets" / "playbook_sets.gmt"
RES.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

SUPPL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334497/suppl/GSE334497_normalized_counts.csv.gz"
MATRIX = DATA / "GSE334497_normalized_counts.csv.gz"
MAP_CACHE = DATA / "ensembl_to_mgi.json"

KO = ["KO162", "KO164", "KO165", "KO172", "RESUB-KO163R"]
WT = ["control170", "RESUB-171R", "RESUB-170R", "RESUB-169R", "RESUB-168R"]

FOCUS_MOUSE = [
    "Tacstd2", "Cldn4", "Epcam",
    "Cldn1", "Cldn3", "Cldn7",
    "Tjp1", "Tjp2", "Ocln", "F11r",
    "Cdh1", "Vim", "Snai1", "Snai2", "Zeb1",
    "Stat1", "Stat2", "Irf1", "Isg15", "Mx1",
]


def log(msg: str) -> None:
    print(msg, flush=True)


def ensure_matrix() -> None:
    if MATRIX.exists():
        log(f"[data] cached {MATRIX.name} ({MATRIX.stat().st_size:,} B)")
        return
    log(f"[data] downloading {SUPPL}")
    urllib.request.urlretrieve(SUPPL, MATRIX)


def load_counts() -> pd.DataFrame:
    with gzip.open(MATRIX, "rt") as fh:
        df = pd.read_csv(fh, index_col=0)
    df.index.name = "ensembl"
    # strip possible version suffix
    df.index = df.index.astype(str).str.replace(r"\.\d+$", "", regex=True)
    return df


def map_symbols(ensembl_ids: list[str]) -> dict[str, str]:
    if MAP_CACHE.exists():
        return json.loads(MAP_CACHE.read_text())
    mapping = {}
    # mygene.info accepts batches of ~1000
    for i in range(0, len(ensembl_ids), 800):
        chunk = ensembl_ids[i:i + 800]
        body = json.dumps({"ids": chunk, "fields": "symbol", "species": "mouse"}).encode()
        req = urllib.request.Request(
            "https://mygene.info/v3/gene", data=body,
            headers={"Content-Type": "application/json", "User-Agent": "Mozilla/5.0"},
        )
        recs = json.loads(urllib.request.urlopen(req, timeout=60).read())
        for rec in recs:
            q = rec.get("query")
            sym = rec.get("symbol")
            if q and sym:
                mapping[q] = sym
        log(f"[map] {min(i+800, len(ensembl_ids))}/{len(ensembl_ids)}")
    MAP_CACHE.write_text(json.dumps(mapping, indent=0))
    return mapping


def bh_fallback(p: np.ndarray) -> np.ndarray:
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return out
    order = np.argsort(pv)
    q = pv[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n)
    tmp[order] = q
    out[ok] = tmp
    return out


def de(df: pd.DataFrame) -> pd.DataFrame:
    t = np.log2(df[KO].to_numpy(float) + 1.0)
    c = np.log2(df[WT].to_numpy(float) + 1.0)
    lfc = t.mean(axis=1) - c.mean(axis=1)
    tstat, p = stats.ttest_ind(t, c, axis=1, equal_var=False, nan_policy="omit")
    out = pd.DataFrame({
        "ensembl": df.index,
        "gene": df["gene"].to_numpy(),
        "log2FC": lfc,
        "t": tstat,
        "pvalue": p,
        "mean_KO": df[KO].mean(axis=1).to_numpy(),
        "mean_WT": df[WT].mean(axis=1).to_numpy(),
    })
    out["padj"] = bh_fallback(out["pvalue"].to_numpy())
    return out.sort_values("log2FC", ascending=False).reset_index(drop=True)


def run_gsea(fc: pd.DataFrame):
    import gseapy as gp
    # Hallmark/KEGG GMT is human-symbol; map mouse -> HUMAN-LIKE via toupper.
    rnk = fc.dropna(subset=["gene", "log2FC"]).copy()
    rnk = rnk[np.isfinite(rnk["log2FC"])]
    rnk["hs"] = rnk["gene"].astype(str).str.upper()
    rnk = rnk.drop_duplicates("hs")
    rnk = rnk.sort_values("log2FC", ascending=False)
    rnk_path = RES / "KO_vs_WT.rnk"
    rnk[["hs", "log2FC"]].to_csv(rnk_path, sep="\t", header=False, index=False)
    if not GMT.exists():
        return None
    pre = gp.prerank(
        rnk=rnk_path.as_posix(), gene_sets=GMT.as_posix(),
        min_size=5, max_size=1000, permutation_num=1000, seed=42,
        outdir=None, no_plot=True, verbose=False,
    )
    res = pre.res2d.copy()
    res.to_csv(RES / "KO_vs_WT_gsea.tsv", sep="\t", index=False)
    return res


def main() -> None:
    ensure_matrix()
    raw = load_counts()
    log(f"[data] {raw.shape[0]} genes x {raw.shape[1]} samples")
    mapping = map_symbols(list(raw.index))
    raw["gene"] = [mapping.get(i, "") for i in raw.index]
    n_mapped = int((raw["gene"] != "").sum())
    log(f"[map] {n_mapped}/{len(raw)} Ensembl IDs -> MGI symbols (mygene.info)")

    fc = de(raw)
    fc.to_csv(RES / "KO_vs_WT_de.tsv.gz", sep="\t", index=False)
    focus = fc[fc["gene"].isin(FOCUS_MOUSE)].copy()
    focus["gene"] = pd.Categorical(focus["gene"], categories=FOCUS_MOUSE, ordered=True)
    focus = focus.sort_values("gene")
    focus.to_csv(RES / "focus_genes.tsv", sep="\t", index=False)

    lines = [
        "=" * 72,
        "GSE334497  4T1 Trop2 KO vs WT  (mouse, n=5, public normalized counts)",
        "=" * 72,
        "NOT DESeq2: normalized (non-integer) counts + Welch t on log2(x+1).",
        "Symbols: MGI via mygene.info; GSEA uses toupper(MGI) vs human GMT.",
        "",
        f"{'gene':<9}{'log2FC':>9}{'padj':>12}{'mean_KO':>12}{'mean_WT':>12}",
    ]
    for _, r in focus.iterrows():
        lines.append(f"{r['gene']:<9}{r['log2FC']:>9.2f}{r['padj']:>12.3g}"
                     f"{r['mean_KO']:>12.1f}{r['mean_WT']:>12.1f}")
    t2 = focus[focus["gene"] == "Tacstd2"]
    c4 = focus[focus["gene"] == "Cldn4"]
    if len(t2):
        r = t2.iloc[0]
        lines.append(f"[on-target] Tacstd2 log2FC={r['log2FC']:+.2f} padj={r['padj']:.3g}")
    if len(c4):
        r = c4.iloc[0]
        got = "UP" if r["log2FC"] > 0 else "DOWN"
        lines.append(f"[opposite ] Cldn4   log2FC={r['log2FC']:+.2f} padj={r['padj']:.3g} -> {got} "
                     f"(hypothesis was UP on Tacstd2 loss)")
    try:
        gsea = run_gsea(fc)
        if gsea is not None:
            lines.append("")
            lines.append("GSEA prerank (positive NES = up in KO; human GMT via toupper):")
            for _, r in gsea.iterrows():
                lines.append(f"  {r['Term']}: NES={r['NES']:.2f}  FDR={r['FDR q-val']:.3g}")
    except Exception as e:
        lines.append(f"[gsea] skipped ({e!r})")

    text = "\n".join(lines)
    (RES / "run_summary.txt").write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
