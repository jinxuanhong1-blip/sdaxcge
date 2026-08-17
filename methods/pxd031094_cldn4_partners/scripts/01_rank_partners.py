#!/usr/bin/env python3
"""Rank PXD031094 CLDN4 CoIP partners from MaxQuant proteinGroups_Cldn4.txt.

Only proteins present in the public SEARCH table are reported.
Gene symbols come from MaxQuant 'Gene names' or UniProt GN= in FASTA headers.
No proteins are invented.
"""
from __future__ import annotations

import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import t as tdist
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"

LFQ_BAIT = [
    "LFQ intensity Cldn4_1",
    "LFQ intensity Cldn4_2",
    "LFQ intensity Cldn4_3",
    "LFQ intensity Cldn4_4",
]
LFQ_CTRL = [
    "LFQ intensity GFPctrl_1",
    "LFQ intensity GFPctrl_2",
    "LFQ intensity GFPctrl_3",
    "LFQ intensity GFPctrl_4",
]

# Category membership is annotation only. Partners still have to exist in MaxQuant.
IFN_SETS = [
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "REACTOME_INTERFERON_SIGNALING",
    "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
]
# Peptide-loading / core MHC only. The large Reactome "class I processing"
# set (381 genes) is not used — it tags ubiquitin/proteasome machinery as MHC.
MHC_SETS = [
    "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
]
TJ_SETS = [
    "REACTOME_TIGHT_JUNCTION_INTERACTIONS",
    "KEGG_TIGHT_JUNCTION",
    "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
    "GOBP_TIGHT_JUNCTION_ORGANIZATION",
    "HALLMARK_APICAL_JUNCTION",
]
TRAFFIC_SETS = [
    "REACTOME_VESICLE_MEDIATED_TRANSPORT",
    "REACTOME_MEMBRANE_TRAFFICKING",
    "REACTOME_RAB_REGULATION_OF_TRAFFICKING",
    "REACTOME_CLATHRIN_MEDIATED_ENDOCYTOSIS",
    "REACTOME_ER_TO_GOLGI_ANTEROGRADE_TRANSPORT",
    "REACTOME_LYSOSOME_VESICLE_BIOGENESIS",
]
IMMUNE_SETS = [
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_ALLOGRAFT_REJECTION",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
]


def parse_gns(header: object) -> list[str]:
    if not isinstance(header, str):
        return []
    out: list[str] = []
    for g in re.findall(r"GN=([^\s;]+)", header):
        if g not in out:
            out.append(g)
    return out


def parse_protein_name(header: object) -> str:
    if not isinstance(header, str):
        return ""
    first = header.split(";")[0]
    m = re.search(r"^[>\w.|-]+ (.+?) OS=", first)
    return m.group(1) if m else ""


def first_uniprot(majority: object) -> str:
    if not isinstance(majority, str) or not majority.strip():
        return ""
    return majority.split(";")[0]


def gene_symbol(row: pd.Series) -> str:
    gn = row.get("Gene names")
    if isinstance(gn, str) and gn.strip():
        return gn.strip().split(";")[0]
    gns = row["gns"]
    if gns:
        return gns[0]
    name = str(row.get("prot_name") or "")
    if re.search(r"DLA class I", name, re.I):
        return "DLA88"
    return ""


def load_category_sets() -> dict[str, set[str]]:
    blob = json.loads((DATA / "category_sets.json").read_text())
    sets = {k: set(v) for k, v in blob["sets"].items()}
    curated = blob.get("curated", {})
    ifn = set().union(*(sets[k] for k in IFN_SETS))
    mhc = set().union(*(sets[k] for k in MHC_SETS))
    mhc |= set(curated.get("MHC_CORE", []))
    tj = set().union(*(sets[k] for k in TJ_SETS))
    tj |= set(curated.get("TJ_CORE", []))
    traffic = set().union(*(sets[k] for k in TRAFFIC_SETS))
    immune = set().union(*(sets[k] for k in IMMUNE_SETS))
    immune |= ifn | mhc
    return {"IFN": ifn, "MHC": mhc, "TJ": tj, "trafficking": traffic, "immune": immune}


def annotate_categories(gene: str, cats: dict[str, set[str]]) -> str:
    if not gene:
        return ""
    g = gene.upper()
    aliases = {gene, g, g.replace("_", "-")}
    # canine MHC class I
    if g.startswith("DLA") or g in {"DLA88", "DLA-88"}:
        aliases.add("DLA88")
    hits = []
    for label in ("IFN", "MHC", "TJ", "trafficking", "immune"):
        s = cats[label]
        if any(a in s for a in aliases):
            hits.append(label)
            continue
        # RAB family: only if the symbol itself is a RAB* gene present in MaxQuant
        if label == "trafficking" and re.fullmatch(r"RAB\d+[A-Z]?|RABL\d+|RAB\d+GAP\d+|RAB\d+FIP\d+", g):
            hits.append(label)
    return ";".join(hits)


def moderated_p(log2fc: pd.Series, logX: pd.DataFrame, bait: list[str], ctrl: list[str]) -> pd.Series:
    """Simple variance-moderated two-sample t (not limma). df_prior=10."""
    s2c = logX[bait].var(axis=1, ddof=1)
    s2g = logX[ctrl].var(axis=1, ddof=1)
    s2 = (3 * s2c + 3 * s2g) / 6
    s2_prior = float(s2.median())
    df_prior = 10.0
    s2_post = (6 * s2 + df_prior * s2_prior) / (6 + df_prior)
    se = np.sqrt(s2_post * (0.25 + 0.25))
    tmod = log2fc / se
    return pd.Series(2 * tdist.sf(np.abs(tmod.to_numpy()), 16), index=log2fc.index)


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    src = DATA / "proteinGroups_Cldn4.txt"
    df = pd.read_csv(src, sep="\t", low_memory=False)
    n_raw = len(df)

    df["gns"] = df["Fasta headers"].map(parse_gns)
    df["prot_name"] = df["Fasta headers"].map(parse_protein_name)
    df["uniprot"] = df["Majority protein IDs"].map(first_uniprot)
    df["gene"] = df.apply(gene_symbol, axis=1)

    keep = (
        df["Reverse"].isna()
        & df["Potential contaminant"].isna()
        & df["Only identified by site"].isna()
        & (df["Unique peptides"] >= 2)
    )
    d = df.loc[keep].copy()
    n_filt = len(d)

    X = d[LFQ_BAIT + LFQ_CTRL].astype(float).replace(0, np.nan)
    pos = X.to_numpy(dtype=float)
    pos = pos[np.isfinite(pos) & (pos > 0)]
    impute = float(pos.min() / 2.0)
    logX = np.log2(X.fillna(impute))

    d["n_cldn4"] = (d[LFQ_BAIT] > 0).sum(axis=1).astype(int)
    d["n_gfp"] = (d[LFQ_CTRL] > 0).sum(axis=1).astype(int)
    d["mean_lfq_cldn4"] = d[LFQ_BAIT].replace(0, np.nan).mean(axis=1)
    d["mean_lfq_gfp"] = d[LFQ_CTRL].replace(0, np.nan).mean(axis=1)
    d["log2fc"] = logX[LFQ_BAIT].mean(axis=1) - logX[LFQ_CTRL].mean(axis=1)

    welch = []
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Precision loss occurred")
        for i in d.index:
            a = logX.loc[i, LFQ_BAIT].to_numpy()
            b = logX.loc[i, LFQ_CTRL].to_numpy()
            _t, pv = stats.ttest_ind(a, b, equal_var=False)
            welch.append(pv)
    d["p_welch"] = welch
    d["p_mod"] = moderated_p(d["log2fc"], logX, LFQ_BAIT, LFQ_CTRL)
    d["fdr_welch"] = multipletests(d["p_welch"].fillna(1.0), method="fdr_bh")[1]
    d["fdr_mod"] = multipletests(d["p_mod"].fillna(1.0), method="fdr_bh")[1]

    cats = load_category_sets()
    d["categories"] = d["gene"].map(lambda g: annotate_categories(g, cats))
    d["is_bait"] = d["gene"].eq("CLDN4") | d["uniprot"].eq("O14493")
    d["exclusive"] = (d["n_cldn4"] >= 2) & (d["n_gfp"] == 0)
    d["stringent"] = (d["n_cldn4"] >= 2) & (d["fdr_mod"] < 0.05) & (d["log2fc"] > 1)

    # Rank partners: bait-detected, unique peptides already filtered.
    partners = d.loc[d["n_cldn4"] >= 2].copy()
    partners["enrich_score"] = partners["log2fc"]
    partners = partners.sort_values(
        ["is_bait", "stringent", "exclusive", "log2fc", "fdr_mod"],
        ascending=[False, False, False, False, True],
    )
    partners["rank"] = np.arange(1, len(partners) + 1)

    cols = [
        "rank",
        "gene",
        "uniprot",
        "prot_name",
        "is_bait",
        "exclusive",
        "stringent",
        "n_cldn4",
        "n_gfp",
        "mean_lfq_cldn4",
        "mean_lfq_gfp",
        "log2fc",
        "p_mod",
        "fdr_mod",
        "p_welch",
        "fdr_welch",
        "Unique peptides",
        "Razor + unique peptides",
        "categories",
        "Majority protein IDs",
    ]
    out = partners[cols].rename(
        columns={
            "prot_name": "protein_name",
            "Unique peptides": "unique_peptides",
            "Razor + unique peptides": "razor_unique_peptides",
            "Majority protein IDs": "majority_protein_ids",
        }
    )
    # Do not fill empty gene with a made-up symbol.
    out["gene"] = out["gene"].replace("", np.nan)

    ranked_path = RESULTS / "cldn4_coip_partners_ranked.tsv"
    out.to_csv(ranked_path, sep="\t", index=False, float_format="%.6g")

    primary = out.loc[out["stringent"] | out["exclusive"]].copy()
    primary["rank"] = np.arange(1, len(primary) + 1)
    primary_path = RESULTS / "cldn4_coip_partners_primary.tsv"
    primary.to_csv(primary_path, sep="\t", index=False, float_format="%.6g")

    mapped = out.loc[out["categories"].fillna("").str.len() > 0].copy()
    mapped["rank"] = np.arange(1, len(mapped) + 1)
    mapped_path = RESULTS / "cldn4_coip_category_mapped.tsv"
    mapped.to_csv(mapped_path, sep="\t", index=False, float_format="%.6g")

    cat_primary = primary.loc[primary["categories"].fillna("").str.len() > 0].copy()
    cat_primary["rank"] = np.arange(1, len(cat_primary) + 1)
    cat_primary_path = RESULTS / "cldn4_coip_category_hits_primary.tsv"
    cat_primary.to_csv(cat_primary_path, sep="\t", index=False, float_format="%.6g")

    # Category hits among primary (non-bait) partners
    prim_nb = primary.loc[~primary["is_bait"]]
    cat_counts = {}
    for lab in ("IFN", "MHC", "TJ", "trafficking", "immune"):
        cat_counts[lab] = int(prim_nb["categories"].fillna("").str.contains(rf"(?:^|;){lab}(?:;|$)").sum())

    summary = {
        "accession": "PXD031094",
        "source_file": "proteinGroups_Cldn4.txt",
        "source_url": URL_NOTE,
        "n_protein_groups_raw": int(n_raw),
        "n_after_reverse_contaminant_site_unique2": int(n_filt),
        "n_detected_in_ge2_cldn4": int(len(partners)),
        "impute_lfq": impute,
        "n_stringent_fdr_mod_0.05_log2fc_gt1": int(partners["stringent"].sum()),
        "n_exclusive_ge2_cldn4_0_gfp": int(partners["exclusive"].sum()),
        "n_primary_stringent_or_exclusive": int(len(primary)),
        "n_primary_nonbait": int((~primary["is_bait"]).sum()),
        "n_primary_with_category": int(len(cat_primary)),
        "bait_recovered": bool(partners["is_bait"].any()),
        "bait_log2fc": float(partners.loc[partners["is_bait"], "log2fc"].iloc[0]) if partners["is_bait"].any() else None,
        "primary_category_counts_nonbait": cat_counts,
        "filters": {
            "remove": ["Reverse", "Potential contaminant", "Only identified by site"],
            "min_unique_peptides": 2,
            "min_cldn4_lfq_replicates": 2,
            "stringent": "moderated two-sample t FDR<0.05 and log2FC>1",
            "exclusive": "LFQ>0 in >=2 Cldn4 and 0 GFP",
        },
        "not_done": [
            "RAW re-search",
            "abundance vs ICI",
            "author Table S1 (Cell Reports supplement blocked)",
        ],
    }
    (RESULTS / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(json.dumps(summary, indent=2))
    print(f"wrote {ranked_path}")
    print(f"wrote {primary_path}")
    print(f"wrote {mapped_path}")
    print(f"wrote {cat_primary_path}")


URL_NOTE = (
    "https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/10/PXD031094/"
    "proteinGroups_Cldn4.txt"
)


if __name__ == "__main__":
    main()
