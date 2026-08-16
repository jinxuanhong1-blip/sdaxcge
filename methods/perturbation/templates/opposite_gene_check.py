#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
opposite_gene_check.py  -  test a directional "compensatory / opposite gene"
hypothesis after a knockdown/knockout.

MOTIVATING QUESTION (from the playbook)
    "If we knock down CLDN4, does TACSTD2 (TROP2) go UP?"
    More generally: after perturbing gene A, does a functionally related /
    paralogous / antagonistic gene B move in the expected direction, with
    statistical support?

This does NOT prove a mechanism; it turns a hand-wavy expectation into a
falsifiable, reported result:
    * effect size  : log2 fold change of B
    * significance : padj of B (from DESeq2/edgeR)
    * context      : where B sits in the genome-wide LFC distribution (rank /
                     percentile) so you can say "top 3% up-regulated" instead
                     of just "went up".

USAGE
    python opposite_gene_check.py de_results.tsv --target CLDN4 \
        --opposite TACSTD2 --expect up --lfc-col log2FoldChange --padj-col padj
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("de_results")
    ap.add_argument("--gene-col", default="gene")
    ap.add_argument("--lfc-col", default="log2FoldChange")
    ap.add_argument("--padj-col", default="padj")
    ap.add_argument("--target", default="CLDN4", help="the perturbed gene")
    ap.add_argument("--opposite", default="TACSTD2", help="the hypothesized responder")
    ap.add_argument("--expect", choices=["up", "down"], default="up")
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()

    df = pd.read_csv(args.de_results, sep="\t")
    if args.gene_col not in df.columns:      # allow index-as-gene tables
        df = df.reset_index().rename(columns={df.index.name or "index": args.gene_col})
    df = df.dropna(subset=[args.lfc_col])
    df = df.sort_values(args.lfc_col, ascending=False).reset_index(drop=True)
    n = len(df)

    def row(sym):
        r = df[df[args.gene_col] == sym]
        return r.iloc[0] if len(r) else None

    def describe(sym, role):
        r = row(sym)
        if r is None:
            print(f"[{role}] {sym}: NOT DETECTED in results.")
            return None
        rank = df.index[df[args.gene_col] == sym][0] + 1
        pct = 100 * rank / n
        padj = r.get(args.padj_col, np.nan)
        print(f"[{role}] {sym}: log2FC={r[args.lfc_col]:+.3f}  "
              f"padj={padj:.3g}  rank={rank}/{n} (top {pct:.1f}% up)")
        return r

    print("=" * 60)
    describe(args.target, "target ")   # on-target: should move strongly (usually down)
    opp = describe(args.opposite, "responder")
    print("-" * 60)
    if opp is not None:
        lfc = opp[args.lfc_col]
        padj = opp.get(args.padj_col, np.nan)
        got = "up" if lfc > 0 else "down"
        sig = (not np.isnan(padj)) and (padj < args.alpha)
        verdict = "SUPPORTED" if (got == args.expect and sig) else (
            "wrong direction" if got != args.expect else "not significant")
        print(f"Hypothesis: {args.opposite} goes {args.expect.upper()} when "
              f"{args.target} is perturbed.")
        print(f"Observed  : {got.upper()} (log2FC={lfc:+.3f}), "
              f"significant={sig} (alpha={args.alpha}).")
        print(f"VERDICT   : {verdict}")
    print("=" * 60)


if __name__ == "__main__":
    main()
