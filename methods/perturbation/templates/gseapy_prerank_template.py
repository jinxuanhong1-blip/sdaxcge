#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
gseapy_prerank_template.py  -  pre-ranked GSEA in Python (gseapy).

Ranks genes by an upstream DE metric (shrunken LFC or signed -log10 p) and tests
Hallmark IFN-alpha / IFN-gamma / EMT + KEGG Tight junction.

    pip install gseapy pandas numpy

Gene sets are fetched once from Enrichr and cached to a local .gmt so the run is
reproducible offline afterwards. For MOUSE data, either map symbols to human
first (see ortholog_map_mouse_human.R) or fetch the mouse Enrichr libraries.

USAGE
    python gseapy_prerank_template.py ranked.rnk out_dir [human|mouse]
"""
from __future__ import annotations
import sys
import urllib.request
from pathlib import Path
import numpy as np
import pandas as pd
import gseapy as gp

ENRICHR = "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName={lib}"
# (library, term-without-parenthetical-size) pairs to keep.
HUMAN_SETS = [
    ("MSigDB_Hallmark_2020", "Interferon Alpha Response"),
    ("MSigDB_Hallmark_2020", "Interferon Gamma Response"),
    ("MSigDB_Hallmark_2020", "Epithelial Mesenchymal Transition"),
    ("KEGG_2021_Human", "Tight junction"),
]
MOUSE_SETS = [
    ("MSigDB_Hallmark_2020", "Interferon Alpha Response"),      # symbols are human-cased
    ("MSigDB_Hallmark_2020", "Interferon Gamma Response"),
    ("MSigDB_Hallmark_2020", "Epithelial Mesenchymal Transition"),
    ("KEGG_2019_Mouse", "Tight junction"),
]


def build_gmt(sets, gmt_path: Path) -> Path:
    if gmt_path.exists():
        return gmt_path
    libs = {lib: self_fetch(lib) for lib in {s[0] for s in sets}}
    wanted = {(lib, term.lower()) for lib, term in sets}
    lines = []
    for lib, text in libs.items():
        for row in text.strip().split("\n"):
            p = row.split("\t")
            base = p[0].split(" (")[0].strip().lower()
            if (lib, base) in wanted:
                genes = [g.split(",")[0].strip() for g in p[2:] if g.strip()]
                lines.append("\t".join([f"{lib}::{p[0].replace(' ', '_')}", ""] + genes))
    gmt_path.write_text("\n".join(lines) + "\n")
    return gmt_path


def self_fetch(lib: str) -> str:
    req = urllib.request.Request(ENRICHR.format(lib=lib),
                                 headers={"User-Agent": "Mozilla/5.0"})
    return urllib.request.urlopen(req, timeout=60).read().decode("utf-8", "replace")


def main(rnk_file, out_dir, species="human"):
    out = Path(out_dir); out.mkdir(parents=True, exist_ok=True)
    sets = MOUSE_SETS if species == "mouse" else HUMAN_SETS
    gmt = build_gmt(sets, out / "gene_sets.gmt")

    rnk = pd.read_csv(rnk_file, sep="\t", header=None, names=["gene", "metric"])
    rnk = rnk[np.isfinite(rnk["metric"])].drop_duplicates("gene")
    rnk = rnk.sort_values("metric", ascending=False)

    pre = gp.prerank(rnk=rnk, gene_sets=gmt.as_posix(), min_size=5, max_size=1000,
                     permutation_num=1000, seed=42, outdir=None, no_plot=True)
    res = pre.res2d.sort_values("NES")
    res.to_csv(out / "gsea_prerank_report.tsv", sep="\t", index=False)
    print(res[["Term", "NES", "NOM p-val", "FDR q-val"]].to_string(index=False))
    print("Done ->", out)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0] if a else "ranked.rnk",
         a[1] if len(a) > 1 else "gseapy_out",
         a[2] if len(a) > 2 else "human")
