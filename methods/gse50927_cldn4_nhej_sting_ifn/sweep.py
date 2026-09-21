#!/usr/bin/env python3
"""Last sweep on GSE50927 deposited EdgeR logFC.

Looks for a NHEJ-down / STING-up mRNA bridge after the core panels were flat:
H2afx and other DNA-damage transcripts, micronucleus-envelope and chromosome-
segregation proxies, cytosolic DNA sensors other than Cgas/Sting, HR versus
NHEJ, and a descriptive injury x genotype contrast-of-contrasts.

GSVA is not run. The series matrix has no expression rows, and the 10 FASTQ
files are two lane-splits of five libraries (one SRX per GSM). Recounting
would still be 1 library vs 1 library.

A genotype contrast is called a NHEJ-down / STING-up slice only when the
NHEJ-core mean deposited logFC is <= -0.25 and the STING-core mean is >= +0.25,
or when at least three NHEJ genes and three STING genes each move by |logFC|
>= 0.5 in that direction. One gene on each side is reported and does not
meet the slice rule. Wilcoxon and GSEA remain gene-set tests on one contrast,
not mouse-level tests.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze import CONTRASTS, compact_score, load_edger, to_mouse
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank

HERE = Path(__file__).resolve().parent
DATA = Path("/tmp/gse50927_nhej")
TAB = HERE / "tables"
FIG = HERE / "figures"
TAB.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

# Pre-specified before reading sweep output. Panel means are KO minus WT.
SLICE_NHEJ_MAX = -0.25
SLICE_STING_MIN = 0.25
GENE_ABS = 0.5
FLAT_MEAN = 0.25

FILES = {
    "naive": DATA / "GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "vili_wt": DATA / "GSE50927_VILIwtGenes.csv.gz",
    "vili_kohi": DATA / "GSE50927_VILIwtkohiGenes.csv.gz",
    "vili_kolo": DATA / "GSE50927_VILIwtkoloGenes.csv.gz",
}

# Sets scored in this sweep. Original 12-set FDR table is left unchanged.
SWEEP_ORDER = [
    "H2AFX_MRNA",
    "DDR_TXN_CORE",
    "HALLMARK_P53_PATHWAY",
    "GO_P53_MEDIATOR_0030330",
    "KEGG_P53_mmu04115",
    "HALLMARK_DNA_REPAIR",
    "DSB_CHROMATIN_CORE",
    "GO_DSB_CHROMATIN_0140861",
    "REACTOME_DSB_SENSING_R-HSA-5693548",
    "HR_CORE",
    "KEGG_HR_mmu03440",
    "GO_HR_0000724",
    "REACTOME_HDR_R-HSA-5685942",
    "NHEJ_CORE",
    "STING_CORE",
    "CYTOSOLIC_SENSORS_BEYOND",
    "KEGG_CYTOSOLIC_BEYOND_CGAS_STING",
    "AIM2_CORE",
    "ZBP1_CORE",
    "POLR3_CYTOSOLIC",
    "MICRONUCLEUS_ENVELOPE",
    "GO_HOMOLOGOUS_CHROMOSOME_SEGREGATION_0045143",
]

LABEL = {
    "H2AFX_MRNA": "H2afx mRNA",
    "DDR_TXN_CORE": "DDR transcript core",
    "HALLMARK_P53_PATHWAY": "Hallmark p53",
    "GO_P53_MEDIATOR_0030330": "GO p53 mediator",
    "KEGG_P53_mmu04115": "KEGG p53",
    "HALLMARK_DNA_REPAIR": "Hallmark DNA repair",
    "DSB_CHROMATIN_CORE": "DSB chromatin core",
    "GO_DSB_CHROMATIN_0140861": "GO DSB chromatin",
    "REACTOME_DSB_SENSING_R-HSA-5693548": "Reactome DSB sensing",
    "HR_CORE": "HR core",
    "KEGG_HR_mmu03440": "KEGG HR",
    "GO_HR_0000724": "GO HR",
    "REACTOME_HDR_R-HSA-5685942": "Reactome HDR",
    "NHEJ_CORE": "NHEJ core",
    "STING_CORE": "STING core",
    "CYTOSOLIC_SENSORS_BEYOND": "Sensors beyond cGAS/STING",
    "KEGG_CYTOSOLIC_BEYOND_CGAS_STING": "KEGG cytosolic DNA minus cGAS/STING/IFN",
    "AIM2_CORE": "AIM2 core",
    "ZBP1_CORE": "ZBP1 core",
    "POLR3_CYTOSOLIC": "RNA pol III cytosolic",
    "MICRONUCLEUS_ENVELOPE": "Micronucleus envelope proxy",
    "GO_HOMOLOGOUS_CHROMOSOME_SEGREGATION_0045143": "GO chromosome segregation",
}


def load_sweep_sets(universe: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    raw = json.loads((HERE / "genesets" / "sweep_panels.json").read_text())
    base = json.loads((HERE / "genesets" / "panels.json").read_text())
    raw["NHEJ_CORE"] = base["NHEJ_CORE"]
    raw["STING_CORE"] = base["STING_CORE"]
    sets: dict[str, list[str]] = {}
    rows = []
    for name in SWEEP_ORDER:
        used, missing = [], []
        seen: set[str] = set()
        for g in raw[name]:
            m = to_mouse(g, universe)
            if m is None:
                missing.append(g)
                continue
            if m not in seen:
                seen.add(m)
                used.append(m)
        sets[name] = used
        rows.append(
            {
                "set": name,
                "label": LABEL[name],
                "n_listed": len(raw[name]),
                "n_in_rank": len(used),
                "n_missing": len(missing),
                "missing": ",".join(missing),
                "symbols_used": ",".join(used) if len(used) <= 24 else "",
            }
        )
    return sets, pd.DataFrame(rows)


def mean_fc(df: pd.DataFrame, genes: list[str]) -> float:
    return float(compact_score(df, genes, flip=False)["mean_logfc_deposited"])


def plot_primary(compact: pd.DataFrame) -> None:
    sub = compact[compact["contrast"] == "naive_KO_vs_WT"].copy()
    sub = sub[sub["set"] != "H2AFX_MRNA"]
    sub["label"] = sub["set"].map(LABEL)
    sub = sub.iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.2, 6.4))
    colors = ["#3C6E9F" if v >= 0 else "#B85C38" for v in sub["mean_logfc_deposited"]]
    ax.barh(sub["label"], sub["mean_logfc_deposited"], color=colors)
    ax.axvline(0, color="black", lw=0.6)
    ax.axvline(SLICE_STING_MIN, color="#888888", lw=0.6, ls="--")
    ax.axvline(SLICE_NHEJ_MAX, color="#888888", lw=0.6, ls="--")
    ax.set_xlabel("Mean deposited logFC (KO − WT) on naive lung")
    ax.set_title("Sweep sets. Dashed lines are the ±0.25 slice cutoffs.")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_sweep_naive_means.png", dpi=160)
    fig.savefig(FIG / "fig4_sweep_naive_means.pdf")
    plt.close(fig)


def plot_interaction(inter: pd.DataFrame) -> None:
    show = [
        "NHEJ_CORE",
        "STING_CORE",
        "HR_CORE",
        "DDR_TXN_CORE",
        "DSB_CHROMATIN_CORE",
        "CYTOSOLIC_SENSORS_BEYOND",
        "MICRONUCLEUS_ENVELOPE",
        "HALLMARK_DNA_REPAIR",
    ]
    sub = inter[inter["set"].isin(show)].copy()
    sub["label"] = sub["set"].map(LABEL)
    # one row per set already
    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    y = np.arange(len(show))
    hi = sub.set_index("set").loc[show, "interaction_high_mean"]
    lo = sub.set_index("set").loc[show, "interaction_low_mean"]
    ax.barh(y + 0.15, hi, height=0.3, color="#D4762C", label="KOhigh effect − naive KO effect")
    ax.barh(y - 0.15, lo, height=0.3, color="#3C6E9F", label="KOlow effect − naive KO effect")
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([LABEL[s] for s in show], fontsize=8)
    ax.set_xlabel("Contrast-of-contrasts (logFC units, 1 library vs 1)")
    ax.legend(frameon=False, fontsize=8)
    ax.set_title("Descriptive injury × genotype difference")
    fig.tight_layout()
    fig.savefig(FIG / "fig5_interaction.png", dpi=160)
    fig.savefig(FIG / "fig5_interaction.pdf")
    plt.close(fig)


def main() -> None:
    frames = {k: load_edger(p) for k, p in FILES.items()}
    universe = set(frames["naive"].index)
    sets, coverage = load_sweep_sets(universe)
    coverage.to_csv(TAB / "sweep_coverage.tsv", sep="\t", index=False)

    gsea_rows = []
    compact_rows = []
    for name, meta in CONTRASTS.items():
        df = frames[meta["key"]]
        rank = df["logFC"].replace([np.inf, -np.inf], np.nan).dropna().sort_values(ascending=False)
        gsea_sets = {k: v for k, v in sets.items() if len(v) >= 6}
        res = gsea_prerank(rank, gsea_sets, nperm=NPERM, seed=SEED, min_size=6)
        res["contrast"] = name
        res["fdr"] = bh_fdr(res["nom_p"])
        gsea_rows.append(res)
        for panel, genes in sets.items():
            row = compact_score(df, genes, flip=meta["flip_to_high_minus_low"])
            row.update({"contrast": name, "role": meta["role"], "set": panel, "label": LABEL[panel]})
            compact_rows.append(row)
    gsea = pd.concat(gsea_rows, ignore_index=True)
    gsea["label"] = gsea["term"].map(LABEL)
    compact = pd.DataFrame(compact_rows)
    gsea.to_csv(TAB / "sweep_gsea.tsv", sep="\t", index=False)
    compact.to_csv(TAB / "sweep_compact.tsv", sep="\t", index=False)

    # Gene-level for the bridge genes plus the single-gene and small-set members.
    focal = []
    for key in (
        "NHEJ_CORE",
        "STING_CORE",
        "HR_CORE",
        "DDR_TXN_CORE",
        "DSB_CHROMATIN_CORE",
        "CYTOSOLIC_SENSORS_BEYOND",
        "AIM2_CORE",
        "ZBP1_CORE",
        "MICRONUCLEUS_ENVELOPE",
        "H2AFX_MRNA",
    ):
        focal.extend(sets[key])
    focal = list(dict.fromkeys(focal))
    gene_rows = []
    for symbol in focal:
        row = {"symbol": symbol}
        for name, meta in CONTRASTS.items():
            df = frames[meta["key"]]
            if symbol not in df.index:
                row[f"{name}_logFC"] = np.nan
                row[f"{name}_logCPM"] = np.nan
                row[f"{name}_FDR"] = np.nan
                continue
            row[f"{name}_logFC"] = float(df.loc[symbol, "logFC"])
            row[f"{name}_logCPM"] = float(df.loc[symbol, "logCPM"])
            row[f"{name}_FDR"] = float(df.loc[symbol, "FDR"])
        row["interaction_high"] = row["VILI_KOhigh_vs_WT_logFC"] - row["naive_KO_vs_WT_logFC"]
        row["interaction_low"] = row["VILI_KOlow_vs_WT_logFC"] - row["naive_KO_vs_WT_logFC"]
        gene_rows.append(row)
    genes = pd.DataFrame(gene_rows)
    genes.to_csv(TAB / "sweep_gene_level.tsv", sep="\t", index=False)

    # Panel interaction on the intersection of genes with logCPM >= 0 in both contrasts.
    inter_rows = []
    for panel, members in sets.items():
        if panel == "H2AFX_MRNA":
            continue
        row = {"set": panel, "label": LABEL[panel]}
        for tag, key in (("high", "vili_kohi"), ("low", "vili_kolo")):
            a = frames[key]
            b = frames["naive"]
            both = []
            for g in members:
                if g not in a.index or g not in b.index:
                    continue
                if float(a.loc[g, "logCPM"]) < 0 or float(b.loc[g, "logCPM"]) < 0:
                    continue
                both.append(float(a.loc[g, "logFC"]) - float(b.loc[g, "logFC"]))
            row[f"interaction_{tag}_n"] = len(both)
            row[f"interaction_{tag}_mean"] = float(np.mean(both)) if both else np.nan
            if len(both) >= 6 and np.any(np.array(both) != 0):
                from scipy.stats import wilcoxon

                try:
                    _, p = wilcoxon(both, zero_method="wilcox", alternative="two-sided")
                except ValueError:
                    p = np.nan
            else:
                p = np.nan
            row[f"interaction_{tag}_wilcoxon_p"] = float(p) if pd.notna(p) else np.nan
        naive_mean = mean_fc(frames["naive"], members)
        hr = mean_fc(frames["naive"], sets["HR_CORE"])
        nhej = mean_fc(frames["naive"], sets["NHEJ_CORE"])
        row["naive_mean_KO_minus_WT"] = naive_mean
        row["naive_HR_minus_NHEJ"] = hr - nhej
        inter_rows.append(row)
    inter = pd.DataFrame(inter_rows)
    inter.to_csv(TAB / "sweep_interaction.tsv", sep="\t", index=False)

    # Balance on every genotype contrast.
    bal_rows = []
    for name, meta in CONTRASTS.items():
        if name == "VILI_WT_vs_naive_WT":
            continue
        df = frames[meta["key"]]
        hr = mean_fc(df, sets["HR_CORE"])
        nhej = mean_fc(df, sets["NHEJ_CORE"])
        bal_rows.append(
            {
                "contrast": name,
                "HR_core_mean_KO_minus_WT": hr,
                "NHEJ_core_mean_KO_minus_WT": nhej,
                "HR_minus_NHEJ": hr - nhej,
            }
        )
    bal = pd.DataFrame(bal_rows)
    bal.to_csv(TAB / "sweep_hr_nhej_balance.tsv", sep="\t", index=False)

    def panel_mean(contrast: str, panel: str) -> float:
        hit = compact[(compact["contrast"] == contrast) & (compact["set"] == panel)]
        return float(hit["mean_logfc_deposited"].iloc[0])

    slice_rows = []
    for name in ("naive_KO_vs_WT", "VILI_KOhigh_vs_WT", "VILI_KOlow_vs_WT"):
        nhej = panel_mean(name, "NHEJ_CORE")
        sting = panel_mean(name, "STING_CORE")
        slice_rows.append(
            {
                "contrast": name,
                "kind": "genotype_panel",
                "NHEJ_mean_KO_minus_WT": nhej,
                "STING_mean_KO_minus_WT": sting,
                "meets_NHEJ_down_STING_up": bool(nhej <= SLICE_NHEJ_MAX and sting >= SLICE_STING_MIN),
            }
        )
    nhej_genes = sets["NHEJ_CORE"]
    sting_genes = sets["STING_CORE"]
    for tag, col in (
        ("naive_KO_vs_WT", "naive_KO_vs_WT_logFC"),
        ("VILI_KOhigh_vs_WT", "VILI_KOhigh_vs_WT_logFC"),
        ("VILI_KOlow_vs_WT", "VILI_KOlow_vs_WT_logFC"),
        ("interaction_high", "interaction_high"),
        ("interaction_low", "interaction_low"),
    ):
        ng = genes[genes["symbol"].isin(nhej_genes)]
        sg = genes[genes["symbol"].isin(sting_genes)]
        n_down = ng.loc[ng[col] <= -GENE_ABS, "symbol"].tolist()
        s_up = sg.loc[sg[col] >= GENE_ABS, "symbol"].tolist()
        # A gene slice requires three genes on each side. One NHEJ gene and one
        # STING gene at |logFC|>=0.5 are listed and do not flip the panel call.
        slice_rows.append(
            {
                "contrast": tag,
                "kind": "genes_|logFC|>=0.5",
                "NHEJ_genes_down": ",".join(n_down),
                "n_NHEJ_down": len(n_down),
                "STING_genes_up": ",".join(s_up),
                "n_STING_up": len(s_up),
                "meets_gene_slice": bool(len(n_down) >= 3 and len(s_up) >= 3),
            }
        )
    # Panel interaction slice
    nhej_i = inter[inter["set"] == "NHEJ_CORE"].iloc[0]
    sting_i = inter[inter["set"] == "STING_CORE"].iloc[0]
    for tag in ("high", "low"):
        nv = float(nhej_i[f"interaction_{tag}_mean"])
        sv = float(sting_i[f"interaction_{tag}_mean"])
        slice_rows.append(
            {
                "contrast": f"interaction_{tag}",
                "kind": "panel_interaction",
                "NHEJ_interaction_mean": nv,
                "STING_interaction_mean": sv,
                "meets_NHEJ_down_STING_up": bool(nv <= SLICE_NHEJ_MAX and sv >= SLICE_STING_MIN),
            }
        )
    slices = pd.DataFrame(slice_rows)
    slices.to_csv(TAB / "sweep_slice.tsv", sep="\t", index=False)

    # Which primary sets are outside the flat window.
    primary = compact[compact["contrast"] == "naive_KO_vs_WT"].set_index("set")
    gprim = gsea[gsea["contrast"] == "naive_KO_vs_WT"].set_index("term")
    nonflat = []
    for panel in SWEEP_ORDER:
        mean = float(primary.loc[panel, "mean_logfc_deposited"])
        p = primary.loc[panel, "wilcoxon_p"]
        fdr = float(gprim.loc[panel, "fdr"]) if panel in gprim.index else np.nan
        mean_move = abs(mean) >= FLAT_MEAN
        p_move = pd.notna(p) and float(p) < 0.05
        fdr_move = pd.notna(fdr) and fdr < 0.05
        if mean_move or p_move or fdr_move:
            nonflat.append(
                {
                    "set": panel,
                    "mean_KO_minus_WT": mean,
                    "wilcoxon_p": None if pd.isna(p) else float(p),
                    "gsea_fdr": None if pd.isna(fdr) else fdr,
                }
            )

    h2 = genes[genes["symbol"] == "H2afx"].iloc[0]
    panel_slices = [r for r in slice_rows if r.get("kind") == "genotype_panel"]
    inter_slices = [r for r in slice_rows if r.get("kind") == "panel_interaction"]
    gene_slices = [r for r in slice_rows if r.get("kind") == "genes_|logFC|>=0.5" and r.get("meets_gene_slice")]
    bridge_flat = (
        not any(r["meets_NHEJ_down_STING_up"] for r in panel_slices)
        and not any(r["meets_NHEJ_down_STING_up"] for r in inter_slices)
        and not gene_slices
    )
    gsva = pd.DataFrame(
        [
            {"item": "series_matrix_bytes", "value": "2328", "note": "FTP GSE50927_series_matrix.txt.gz"},
            {"item": "expression_rows", "value": "0", "note": "header ID_REF plus five GSM columns only"},
            {"item": "author_tables", "value": "4 EdgeR logFC tables", "note": "no per-sample counts"},
            {"item": "sra_runs", "value": "10", "note": "SRR988118–SRR988127"},
            {"item": "experiments", "value": "5", "note": "SRX352050–SRX352054, one per GSM"},
            {"item": "runs_per_experiment", "value": "2", "note": "same sample_accession; lane splits, not a second mouse"},
            {"item": "gsva_computed", "value": "no", "note": "a recount would remain 1 library vs 1 library"},
            {"item": "custom_sets_scored_on", "value": "deposited logFC", "note": "sweep_compact.tsv and sweep_gsea.tsv"},
        ]
    )
    gsva.to_csv(TAB / "gsva_status.tsv", sep="\t", index=False)

    verdict = {
        "n": "1 vs 1 GSM",
        "slice_rule": f"NHEJ mean<={SLICE_NHEJ_MAX} and STING mean>={SLICE_STING_MIN} on KO-minus-WT logFC",
        "gene_rule": f"at least 3 NHEJ genes and 3 STING genes with |logFC|>={GENE_ABS} in the matching direction",
        "bridge_slice_found": (not bridge_flat),
        "panel_slices": panel_slices,
        "interaction_slices": inter_slices,
        "gene_slices_meeting_rule": gene_slices,
        "naive_nonflat_sets": nonflat,
        "H2afx_naive_logFC": float(h2["naive_KO_vs_WT_logFC"]),
        "H2afx_naive_FDR": float(h2["naive_KO_vs_WT_FDR"]),
        "H2afx_naive_logCPM": float(h2["naive_KO_vs_WT_logCPM"]),
        "final_call": (
            "FINAL: non-cancer lung supports IFN after Cldn4 loss but not NHEJ–STING mRNA bridge"
            if bridge_flat
            else "A pre-specified NHEJ-down / STING-up slice is present; do not mark the bridge flat"
        ),
    }
    (TAB / "sweep_verdict.json").write_text(json.dumps(verdict, indent=2) + "\n")
    plot_primary(compact)
    plot_interaction(inter)
    print(verdict["final_call"])
    print("nonflat", json.dumps(nonflat, indent=2))
    print(bal.to_string(index=False))
    print(slices.to_string(index=False))
    cols = ["contrast", "set", "n_used", "mean_logfc_deposited", "wilcoxon_p"]
    print(compact.loc[compact["contrast"] == "naive_KO_vs_WT", cols].to_string(index=False))


if __name__ == "__main__":
    main()
