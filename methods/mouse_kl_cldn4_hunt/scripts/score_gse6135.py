#!/usr/bin/env python3
"""Mouse/tumor-level Cldn4 vs T/NK and IFN/MHC in GSE6135 (Ji et al. KL microarray)."""

from __future__ import annotations

import gzip
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from gene_sets import CLDN4, IFN_MHC, T_NK  # noqa: E402
from score_utils import association_table, mean_score, present_genes, zscore_rows  # noqa: E402

RAW = Path("/tmp/geo_dl")
OUT = Path("/workspace/analysis/GSE6135")


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    titles = None
    geo = None
    chars = None
    table_lines = []
    in_table = False
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                titles = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_geo_accession"):
                geo = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1") and chars is None:
                chars = [x.strip().strip('"') for x in line.rstrip().split("\t")[1:]]
            elif line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            elif line.startswith("!series_matrix_table_end"):
                break
            elif in_table:
                table_lines.append(line)
    from io import StringIO

    expr = pd.read_csv(StringIO("".join(table_lines)), sep="\t")
    expr = expr.set_index("ID_REF")
    expr.columns = [c.strip('"') for c in expr.columns]
    meta = pd.DataFrame({"geo": geo, "title": titles, "characteristics": chars})
    return expr, meta


def load_annot(path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(path, "rt", errors="replace") as f:
        go = False
        header = None
        for line in f:
            if line.startswith("!platform_table_begin"):
                go = True
                header = next(f).rstrip().split("\t")
                continue
            if line.startswith("!platform_table_end"):
                break
            if go:
                parts = line.rstrip("\n").split("\t")
                rows.append(parts)
    ann = pd.DataFrame(rows, columns=header)
    return ann[["ID", "Gene symbol"]].rename(columns={"ID": "probe", "Gene symbol": "symbol"})


def classify_genotype(title: str, characteristics: str) -> str:
    t = f"{title} {characteristics}".lower()
    if "lkb" in t or "stk11" in t:
        return "KL"
    if "p53" in t:
        return "KP"
    if "p16" in t or "ink4" in t:
        return "K_Ink4a"
    if "kras" in t or "k-ras" in t:
        return "K"
    return "other"


def mouse_from_title(title: str) -> str:
    # e.g. "592T1 || K-ras Lkb1 L/L || primary lung tumor || Ad"
    head = title.split("||")[0].strip()
    head = head.replace(" ", "")
    m = re.match(r"([0-9]+)(?:\([^)]+\))?[-_]?T?\d*", head)
    if m:
        return m.group(1)
    return head


def collapse_symbol_matrix(expr: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    ann = ann.copy()
    # Multi-symbol probes: keep if any wanted symbol is exact later; split on ///
    probe_to_syms = {}
    for _, r in ann.iterrows():
        raw = str(r["symbol"]) if pd.notna(r["symbol"]) else ""
        if not raw or raw == "nan":
            continue
        probe_to_syms[r["probe"]] = [s.strip() for s in raw.split("///") if s.strip()]
    wanted = set(CLDN4 + T_NK + IFN_MHC)
    gene_rows = {g: [] for g in wanted}
    for probe, syms in probe_to_syms.items():
        if probe not in expr.index:
            continue
        for s in syms:
            if s in gene_rows:
                gene_rows[s].append(probe)
    out = {}
    for g, probes in gene_rows.items():
        if probes:
            out[g] = expr.loc[probes].mean(axis=0)
    return pd.DataFrame(out).T


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    expr, meta = parse_series_matrix(RAW / "GSE6135-GPL8321_series_matrix.txt.gz")
    ann = load_annot(RAW / "GPL8321.annot.gz")
    gene_expr = collapse_symbol_matrix(expr, ann)
    meta = meta.set_index("geo")
    meta["genotype"] = [classify_genotype(t, c) for t, c in zip(meta["title"], meta["characteristics"])]
    meta["mouse_id"] = [
        f"{g}_{mouse_from_title(t)}"
        for g, t in zip(meta["genotype"], meta["title"])
    ]
    meta["is_metastasis"] = meta["title"].str.contains(r"\bMet\b", case=False, regex=True)

    present = {
        "Cldn4": present_genes(gene_expr.index, CLDN4),
        "T_NK": present_genes(gene_expr.index, T_NK),
        "IFN_MHC": present_genes(gene_expr.index, IFN_MHC),
    }
    z = zscore_rows(gene_expr)
    tumor = meta.copy()
    tumor["cldn4"] = gene_expr.loc["Cldn4"] if "Cldn4" in gene_expr.index else np.nan
    tumor["tnk_score"] = mean_score(z, present["T_NK"])
    tumor["ifn_mhc_score"] = mean_score(z, present["IFN_MHC"])
    tumor["unit"] = "tumor"
    tumor = tumor.reset_index().rename(columns={"geo": "sample_id"})

    # Mouse-level: mean tumors from the same parsed mouse id.
    g = tumor.groupby("mouse_id", sort=False)
    mouse = pd.DataFrame(
        {
            "mouse_id": g.size().index,
            "genotype": g["genotype"].first(),
            "n_tumors": g.size().values,
            "cldn4_mean": g["cldn4"].mean().values,
            "tnk_score": g["tnk_score"].mean().values,
            "ifn_mhc_score": g["ifn_mhc_score"].mean().values,
            "sample_ids": g["sample_id"].apply(lambda s: ",".join(s)).values,
            "titles": g["title"].apply(lambda s: " | ".join(s)).values,
        }
    )

    assoc_tumor = association_table(tumor, "cldn4", ["tnk_score", "ifn_mhc_score"])
    assoc_tumor.insert(0, "unit", "tumor")
    assoc_tumor.insert(1, "subset", "all_mouse_tumors")
    assoc_mouse = association_table(mouse, "cldn4_mean", ["tnk_score", "ifn_mhc_score"])
    assoc_mouse.insert(0, "unit", "mouse")
    assoc_mouse.insert(1, "subset", "all_mice")
    parts = [assoc_tumor, assoc_mouse]
    for geno, sub in tumor.groupby("genotype"):
        a = association_table(sub, "cldn4", ["tnk_score", "ifn_mhc_score"])
        a.insert(0, "unit", "tumor")
        a.insert(1, "subset", f"genotype={geno}")
        parts.append(a)
    assoc = pd.concat(parts, ignore_index=True)

    genes_used = pd.DataFrame(
        {
            "set": ["Cldn4", "T_NK", "IFN_MHC"],
            "requested": [",".join(CLDN4), ",".join(T_NK), ",".join(IFN_MHC)],
            "present": [
                ",".join(present["Cldn4"]),
                ",".join(present["T_NK"]),
                ",".join(present["IFN_MHC"]),
            ],
            "n_present": [len(present["Cldn4"]), len(present["T_NK"]), len(present["IFN_MHC"])],
        }
    )

    # KL vs other Cldn4 (tumor-level). Honest n.
    kl = tumor.loc[tumor["genotype"] == "KL", "cldn4"]
    other = tumor.loc[tumor["genotype"] != "KL", "cldn4"]
    from scipy import stats as _stats

    kl_row = {
        "contrast": "KL_vs_other_Cldn4",
        "n_KL_tumors": int(len(kl)),
        "n_other_tumors": int(len(other)),
        "mean_Cldn4_KL": float(kl.mean()) if len(kl) else float("nan"),
        "mean_Cldn4_other": float(other.mean()) if len(other) else float("nan"),
        "mannwhitney_p": float("nan"),
    }
    if len(kl) >= 2 and len(other) >= 2:
        _, p = _stats.mannwhitneyu(kl, other, alternative="two-sided")
        kl_row["mannwhitney_p"] = float(p)
    pd.DataFrame([kl_row]).to_csv(OUT / "kl_vs_other_cldn4.tsv", sep="\t", index=False)

    tumor.to_csv(OUT / "tumor_level_scores.tsv", sep="\t", index=False)
    mouse.to_csv(OUT / "mouse_level_scores.tsv", sep="\t", index=False)
    assoc.to_csv(OUT / "high_vs_low_summary.tsv", sep="\t", index=False)
    genes_used.to_csv(OUT / "genes_used.tsv", sep="\t", index=False)

    print(
        f"GSE6135: n_tumors={len(tumor)} n_mice={len(mouse)} "
        f"Cldn4={bool(present['Cldn4'])} T/NK={present['T_NK']} IFN/MHC={present['IFN_MHC']}"
    )
    print(tumor.groupby("genotype").size().to_string())
    print(assoc.to_string(index=False))


if __name__ == "__main__":
    main()
