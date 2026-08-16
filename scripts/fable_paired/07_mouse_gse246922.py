"""GSE246922 — mouse lung ICB acquired-resistance models (Memon et al. 2024).

KP (Kras/p53) and LLC1 lung-cancer lines: parental vs IFNγ-exposed vs ICB-resistant
derivatives. This is the closest *public processed* mouse-lung ICB analog to the
TISMO 'Tacstd2 rises after IO' claim (TISMO itself is a Shiny portal; bulk
download endpoints were not exposed without the UI).

Expression is already VST (approx. log2). Genes via Ensembl:
  Tacstd2 ENSMUSG00000051397
  Cldn4   ENSMUSG00000047501
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import config as C
from common import compare_groups, savefig, PALETTE, strip_box
import matplotlib.pyplot as plt

MOUSE = {"Tacstd2": "ENSMUSG00000051397", "Cldn4": "ENSMUSG00000047501"}


def _load(path, model: str) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    recs = []
    for gene, ens in MOUSE.items():
        if ens not in df.index:
            print(f"[warn] {ens} ({gene}) missing from {model}")
            continue
        for col, val in df.loc[ens].items():
            recs.append({"model": model, "sample": col, "gene": gene, "expr": float(val)})
    return pd.DataFrame(recs)


def _group(sample: str, model: str) -> str:
    s = sample.lower()
    if model == "KP":
        if s.startswith("kpy"):
            return "IFNg"
        if s.startswith("res"):
            return "ICB_resistant"
        if s.startswith("kp"):
            return "parental"
    if model == "LLC1":
        if s.startswith("llc1y"):
            return "IFNg"
        if s.startswith("res"):
            return "ICB_resistant"
        if s.startswith("llc1"):
            return "parental"
    return "other"


def main() -> None:
    kp = _load(C.raw_path("GSE246922_KP"), "KP")
    llc = _load(C.raw_path("GSE246922_LLC1"), "LLC1")
    long = pd.concat([kp, llc], ignore_index=True)
    long["group"] = [ _group(s, m) for s, m in zip(long["sample"], long["model"]) ]
    long.to_csv(C.TABLES_DIR / "gse246922_mouse_long.csv", index=False)
    print(long.groupby(["model", "group", "gene"]).size())

    rows = []
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    order = ["parental", "IFNg", "ICB_resistant"]
    colors = {"parental": PALETTE["pre"], "IFNg": PALETTE["on"],
              "ICB_resistant": PALETTE["post"]}
    for i, model in enumerate(["KP", "LLC1"]):
        for j, gene in enumerate(["Tacstd2", "Cldn4"]):
            sub = long[(long.model == model) & (long.gene == gene)]
            groups = {g: sub.loc[sub.group == g, "expr"].values for g in order
                      if (sub.group == g).any()}
            strip_box(axes[i, j], groups, colors,
                      f"{gene} (VST)", f"{model} {gene}")
            if "parental" in groups and "ICB_resistant" in groups:
                rows.append(compare_groups(
                    f"GSE246922_{model}", gene, "ICB_resistant_vs_parental",
                    groups["ICB_resistant"], groups["parental"],
                    "ICB_resistant", "parental").as_row())
            if "parental" in groups and "IFNg" in groups:
                rows.append(compare_groups(
                    f"GSE246922_{model}", gene, "IFNg_vs_parental",
                    groups["IFNg"], groups["parental"],
                    "IFNg", "parental").as_row())
    fig.suptitle("GSE246922 mouse lung ICB models (KP, LLC1) — Tacstd2 / Cldn4", fontsize=11)
    savefig(fig, C.FIGURES_DIR / "fig6_gse246922_mouse_lung_icb.png")
    out = pd.DataFrame(rows)
    out.to_csv(C.TABLES_DIR / "gse246922_mouse_stats.csv", index=False)
    print(out.to_string())


if __name__ == "__main__":
    main()
