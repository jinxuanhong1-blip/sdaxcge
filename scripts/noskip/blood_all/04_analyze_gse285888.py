#!/usr/bin/env python3
"""Stream GSE285888 PBMC scRNA matrix; patient-level TACSTD2/CLDN4 vs CR/DR/PD."""
from __future__ import annotations

import gzip
import json
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

DATA = Path("/tmp/geo_blood/gse285888")
OUT = Path("results/noskip/blood_all")
OUT.mkdir(parents=True, exist_ok=True)

PANEL = ["TACSTD2", "CLDN4", "EPCAM", "KRT8", "KRT18", "PTPRC", "CD3D", "CD8A", "GZMB", "PRF1", "IL1B", "CXCL8"]


def main():
    meta_path = DATA / "GSM8712033_metadata.csv.gz"
    mat_path = DATA / "GSM8712033_matrix.txt.gz"
    meta = pd.read_csv(meta_path)
    print("meta columns", list(meta.columns)[:30], "n", len(meta), flush=True)
    # guess cell id / patient / response
    print(meta.head(2).to_string(), flush=True)
    # common Seurat export: Unnamed: 0 is cell barcode
    cell_col = meta.columns[0]
    # find patient and response-like columns
    print("dtypes ok", flush=True)
    for c in meta.columns:
        nunq = meta[c].nunique(dropna=True)
        if nunq < 40:
            print(f"  col {c} nunique={nunq} vals={meta[c].dropna().astype(str).value_counts().head(8).to_dict()}")

    # Stream matrix: first line header = cell barcodes
    with gzip.open(mat_path, "rt") as f:
        header = [x.strip().strip('"') for x in f.readline().rstrip("\n").split("\t")]
        # first token may be gene or empty
        if header[0] in ("", "gene", "GENE", "index"):
            cells = header[1:]
            gene_first = True
        else:
            # maybe first is a cell
            cells = header
            gene_first = False
        print("n cells in matrix header", len(cells), "gene_first", gene_first, "ex", cells[:3], flush=True)
        want = set(PANEL)
        found = {}
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            gene = parts[0].strip().strip('"')
            if gene in want:
                vals = np.array(parts[1:], dtype=float)
                found[gene] = vals
                print("got", gene, "nz", int((vals > 0).sum()), "mean", float(vals.mean()), flush=True)
            if len(found) == len(want):
                break
            if i % 5000 == 0 and i:
                print("scanned", i, "found", list(found), flush=True)

    # align meta to cells
    # metadata row order often matches matrix columns
    if len(meta) == len(cells):
        meta = meta.copy()
        meta["_cell"] = cells
    else:
        # try match
        meta = meta.copy()
        meta["_cell"] = meta[cell_col].astype(str)

    # patient / response
    pat_col = None
    for c in ("orig.ident", "patient", "Patient", "donor", "sample"):
        if c in meta.columns:
            pat_col = c
            break
    if pat_col is None:
        for c in meta.columns:
            if meta[c].nunique() < 50 and meta[c].nunique() > 5:
                pat_col = c
                break
    resp_col = None
    for c in meta.columns:
        cl = c.lower()
        if any(k in cl for k in ("response", "subtype", "group", "outcome")):
            resp_col = c
            break
    irae_col = None
    for c in meta.columns:
        if "irae" in c.lower() or "irAE" in c:
            irae_col = c
            break
    print("pat", pat_col, "resp", resp_col, "irae", irae_col, flush=True)

    # nCount
    ncount = None
    for c in meta.columns:
        if c.lower() in ("ncount_rna", "ncount", "n_counts"):
            ncount = meta[c].to_numpy(float)
            break

    cell_to_i = {c: i for i, c in enumerate(cells)}
    # if meta _cell matches
    order = []
    for c in meta["_cell"].astype(str):
        order.append(cell_to_i.get(c))
    if sum(x is not None for x in order) < 0.9 * len(meta):
        # assume same order
        order = list(range(min(len(meta), len(cells))))
        meta = meta.iloc[: len(order)]

    # patient pseudobulk
    recs = []
    detect_rows = []
    for gene, vec in found.items():
        det = float((vec > 0).mean())
        detect_rows.append({"dataset": "GSE285888", "gene": gene, "n_cells": int(len(vec)), "frac_detected": det, "mean_count": float(vec.mean())})
        # per patient
        by_pat = defaultdict(lambda: {"counts": 0.0, "n": 0, "nz": 0, "lib": 0.0})
        for mi, oi in enumerate(order):
            if oi is None:
                continue
            pat = str(meta.iloc[mi][pat_col]) if pat_col else "NA"
            by_pat[pat]["counts"] += float(vec[oi])
            by_pat[pat]["n"] += 1
            by_pat[pat]["nz"] += int(vec[oi] > 0)
            if ncount is not None:
                by_pat[pat]["lib"] += float(ncount[mi])
        for pat, d in by_pat.items():
            # response from majority of cells
            sub = meta[meta[pat_col].astype(str) == pat] if pat_col else meta
            resp = str(sub[resp_col].mode().iloc[0]) if resp_col and len(sub) else ""
            irae = str(sub[irae_col].mode().iloc[0]) if irae_col and len(sub) else ""
            cp10k = (d["counts"] / d["lib"] * 1e4) if d["lib"] > 0 else d["counts"] / max(d["n"], 1)
            recs.append(
                {
                    "dataset": "GSE285888",
                    "compartment": "PBMC_scRNA",
                    "patient": pat,
                    "response_raw": resp,
                    "irae": irae,
                    "gene": gene,
                    "n_cells": d["n"],
                    "frac_detected": d["nz"] / max(d["n"], 1),
                    "cp10k": float(cp10k),
                    "sum_count": d["counts"],
                }
            )

    pb = pd.DataFrame(recs)
    if pb.empty:
        raise SystemExit(f"no genes extracted; found={list(found)}")
    pb.to_csv(OUT / "GSE285888_patient_pseudobulk.tsv", sep="\t", index=False)
    pd.DataFrame(detect_rows).to_csv(OUT / "GSE285888_detectability.tsv", sep="\t", index=False)

    # contrasts
    stats_rows = []
    def groups(df, gene, col, pos, neg):
        sub = df[df.gene == gene]
        a = sub[sub[col].astype(str).str.contains(pos, case=False, na=False)]["cp10k"]
        b = sub[sub[col].astype(str).str.contains(neg, case=False, na=False)]["cp10k"]
        if len(a) < 2 or len(b) < 2:
            return None
        U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return {"n1": int(len(a)), "n2": int(len(b)), "median1": float(a.median()), "median2": float(b.median()), "p": float(p)}

    # classify response
    def norm_resp(s):
        s = str(s).upper()
        if "CR" in s.split() or s.startswith("CR") or "COMPLETE" in s:
            return "CR"
        if "DR" in s or "DURABLE" in s:
            return "DR"
        if "PD" in s or "PROGRESS" in s:
            return "PD"
        return s

    pb["resp3"] = pb["response_raw"].map(norm_resp)
    print(pb.drop_duplicates("patient")["resp3"].value_counts().to_dict(), flush=True)
    print(pb.drop_duplicates("patient")["response_raw"].value_counts().head(15).to_dict(), flush=True)

    for gene in ("TACSTD2", "CLDN4", "IL1B", "CXCL8", "PTPRC", "EPCAM"):
        if gene not in found:
            continue
        sub = pb[pb.gene == gene]
        # CR vs PD
        a = sub[sub.resp3 == "CR"]["cp10k"]
        b = sub[sub.resp3 == "PD"]["cp10k"]
        r = sub[sub.resp3.isin(["CR", "DR"])]["cp10k"]
        if len(a) >= 2 and len(b) >= 2:
            U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            stats_rows.append({"dataset": "GSE285888", "gene": gene, "contrast": "CR vs PD", "metric": "patient_cp10k", "n1": int(len(a)), "n2": int(len(b)), "median1": float(a.median()), "median2": float(b.median()), "p": float(p)})
        if len(r) >= 2 and len(b) >= 2:
            U, p = stats.mannwhitneyu(r, b, alternative="two-sided")
            stats_rows.append({"dataset": "GSE285888", "gene": gene, "contrast": "CR+DR vs PD", "metric": "patient_cp10k", "n1": int(len(r)), "n2": int(len(b)), "median1": float(r.median()), "median2": float(b.median()), "p": float(p)})
        # detection rate CR vs PD
        a = sub[sub.resp3 == "CR"]["frac_detected"]
        b = sub[sub.resp3 == "PD"]["frac_detected"]
        if len(a) >= 2 and len(b) >= 2:
            U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            stats_rows.append({"dataset": "GSE285888", "gene": gene, "contrast": "CR vs PD", "metric": "patient_frac_detected", "n1": int(len(a)), "n2": int(len(b)), "median1": float(a.median()), "median2": float(b.median()), "p": float(p)})

    pd.DataFrame(stats_rows).to_csv(OUT / "GSE285888_response_stats.tsv", sep="\t", index=False)
    # plot detectability
    fig, ax = plt.subplots(figsize=(7, 4))
    genes = [r["gene"] for r in detect_rows]
    fracs = [r["frac_detected"] * 100 for r in detect_rows]
    colors = ["#c0392b" if g in ("TACSTD2", "CLDN4") else "#2980b9" for g in genes]
    ax.bar(genes, fracs, color=colors)
    ax.set_ylabel("% cells with count > 0")
    ax.set_title("GSE285888 baseline PBMC scRNA detectability")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(OUT / "GSE285888_detectability.png", dpi=140)
    plt.close(fig)
    print(json.dumps(detect_rows, indent=2))
    print(pd.DataFrame(stats_rows).to_string(index=False))


if __name__ == "__main__":
    main()
