#!/usr/bin/env python3
"""Sweep DNA-damage, micronuclei, cytosolic-DNA, and NHEJ panels.

Same contrast as analyze.py: GSE50927 naive Cldn4 KO vs WT, no VILI.
logFC is KO minus WT. n = 1 vs 1 GSM. There is no genes-by-sample matrix
(series matrix Sample_data_row_count = 0), so GSVA and ssGSEA are not run.
The enrichment that this object supports is preranked GSEA on the deposited
logFC rank, plus mean logFC, gene-label permutation, IFN-gene removal, and
leave-one-out.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from gsea_prerank import SEED, bh_fdr, gsea_prerank

ROOT = Path(__file__).resolve().parent
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"
SETDIR = ROOT / "genesets"
DATA = Path(os.environ.get("GSE50927_NHEJ_STING_DATA", "/tmp/gse50927_nhej_sting"))
URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/"
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
)
SOURCE = "GSE50927_Cldn4lungWTvsKOgenes.csv.gz"
NPERM_MEAN = 5000
MIN_LOGCPM = 0.0

# KEGG / old symbols that are not the mm9 symbol in this table.
ALIASES = {
    "CGAS": "MB21D1",
    "STING1": "TMEM173",
    "PAXX": "1110057K04RIK",
    "H2AX": "H2AFX",
}

# Explicit panels. Membership is written to tables/sweep_set_members.tsv.
CUSTOM = {
    "custom_NHEJ_core": (
        "repair",
        ["Xrcc6", "Xrcc5", "Prkdc", "Dclre1c", "Nhej1", "Xrcc4", "Lig4", "Poll", "Polm", "1110057K04Rik"],
    ),
    "custom_NHEJ_Ku_DNAPKcs": (
        "repair",
        ["Xrcc6", "Xrcc5", "Prkdc"],
    ),
    "custom_NHEJ_ligation": (
        "repair",
        ["Xrcc4", "Lig4", "Nhej1", "1110057K04Rik"],
    ),
    "custom_altNHEJ_MMEJ": (
        "repair",
        ["Parp1", "Lig3", "Xrcc1", "Polq", "Mre11a", "Rad50", "Nbn", "Fen1", "Lig1"],
    ),
    "custom_nuclear_envelope": (
        "structure",
        [
            "Lmna", "Lmnb1", "Lmnb2", "Banf1", "Emd", "Lemd2", "Lemd3", "Tmpo",
            "Sun1", "Sun2", "Syne1", "Syne2", "Nup153", "Nup50", "Nup98", "Nup210",
            "Nup85", "Nup107", "Aaas", "Tor1a",
        ],
    ),
    "custom_chromosome_missegregation": (
        "proliferation",
        [
            "Aurka", "Aurkb", "Plk1", "Bub1", "Bub1b", "Mad2l1", "Mad2l2", "Espl1",
            "Incenp", "Cdca8", "Cenpa", "Cenpe", "Cenpf", "Kif2c", "Sgo1", "Pttg1",
            "Knl1", "Ndc80", "Nusap1", "Top2a",
        ],
    ),
    "custom_cytosolic_DNA_sensors": (
        "sensor",
        [
            "Mb21d1", "Tmem173", "Tbk1", "Irf3", "Zbp1", "Aim2", "Ddx41",
            "Ifi204", "Ifi203", "Ifi202b", "Pycard", "Casp1", "Dhx9", "Dhx36",
        ],
    ),
    "custom_cytosolic_DNA_negative_regulators": (
        "sensor",
        ["Trex1", "Samhd1", "Rnaseh2a", "Rnaseh2b", "Rnaseh2c", "Adar"],
    ),
    "custom_damage_arrest_transcripts": (
        "damage_signal",
        [
            "Cdkn1a", "Gadd45a", "Gadd45b", "Bbc3", "Bax", "Ddb2", "Xpc", "Polk",
            "Mgmt", "Trp53", "Mdm2", "Ccng1", "Sesn1", "Sesn2", "Fas",
        ],
    ),
    "custom_proliferation_markers": (
        "proliferation",
        [
            "Mki67", "Top2a", "Pcna", "Mcm2", "Mcm3", "Mcm4", "Mcm5", "Mcm6",
            "Cdk1", "Ccnb1", "Ccnb2", "Ccna2", "Birc5", "Ube2c",
        ],
    ),
    "custom_phospho_proxy_ISG": (
        "output",
        [
            "Ccl5", "Cxcl9", "Cxcl10", "Isg15", "Ifit1", "Ifit2", "Ifit3",
            "Rsad2", "Stat1", "Irf7", "Oasl2", "Gbp4", "Mx1",
        ],
    ),
}

HALLMARK_CLASS = {
    "HALLMARK_DNA_REPAIR": "repair",
    "HALLMARK_P53_PATHWAY": "damage_signal",
    "HALLMARK_UV_RESPONSE_UP": "damage_signal",
    "HALLMARK_UV_RESPONSE_DN": "damage_signal",
    "HALLMARK_APOPTOSIS": "damage_signal",
    "HALLMARK_G2M_CHECKPOINT": "proliferation",
    "HALLMARK_MITOTIC_SPINDLE": "proliferation",
    "HALLMARK_E2F_TARGETS": "proliferation",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "output",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "output",
    "HALLMARK_INFLAMMATORY_RESPONSE": "output",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": "output",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING": "output",
}

KEGG_CLASS = {
    "KEGG_Non-homologous_end-joining": "repair",
    "KEGG_Homologous_recombination": "repair",
    "KEGG_Base_excision_repair": "repair",
    "KEGG_Mismatch_repair": "repair",
    "KEGG_Nucleotide_excision_repair": "repair",
    "KEGG_DNA_replication": "proliferation",
    "KEGG_Cell_cycle": "proliferation",
    "KEGG_p53_signaling_pathway": "damage_signal",
    "KEGG_Apoptosis": "damage_signal",
    "KEGG_Cytosolic_DNA-sensing_pathway": "sensor",
}

# Classes allowed to compete as an upstream bridge. Proliferation is reported
# separately because an IFN-high whole lung is full of cycling immune cells.
BRIDGE_CLASSES = {"repair", "damage_signal", "sensor", "structure"}


def download() -> Path:
    DATA.mkdir(parents=True, exist_ok=True)
    dest = DATA / SOURCE
    if not dest.exists() or dest.stat().st_size < 1000:
        urllib.request.urlretrieve(URL, dest)
    return dest


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_edger(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.rename(columns={"Marker.Symbol": "symbol"})
    df["symbol"] = df["symbol"].astype(str)
    for col in ("logFC", "logCPM", "PValue", "FDR"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.sort_values("logCPM", ascending=False).drop_duplicates("symbol")
    df = df.set_index("symbol", drop=False)
    return df


def resolve_members(members: list[str], upper: dict[str, str]) -> tuple[list[str], list[str]]:
    found: list[str] = []
    missing: list[str] = []
    seen: set[str] = set()
    for raw in members:
        key = ALIASES.get(raw.upper(), raw.upper())
        hit = upper.get(key)
        if hit is None:
            missing.append(raw)
            continue
        if hit in seen:
            continue
        seen.add(hit)
        found.append(hit)
    return found, missing


def perm_p(values: np.ndarray, universe: np.ndarray, nperm: int, seed: int) -> float:
    if len(values) == 0 or len(universe) <= len(values):
        return math.nan
    obs = float(values.mean())
    rng = np.random.default_rng(seed)
    exceed = 0
    k = len(values)
    for _ in range(nperm):
        draw = rng.choice(universe, size=k, replace=False)
        if abs(float(draw.mean())) >= abs(obs) - 1e-15:
            exceed += 1
    return (exceed + 1) / (nperm + 1)


def leave_one_out(symbols: list[str], logfc: pd.Series) -> dict:
    if len(symbols) < 2:
        return {
            "loo_min_mean": math.nan,
            "loo_max_mean": math.nan,
            "loo_worst_gene": "",
            "loo_sign_flip": "",
            "loo_n": 0,
        }
    full = float(logfc.loc[symbols].mean())
    means = []
    for gene in symbols:
        others = [g for g in symbols if g != gene]
        means.append((gene, float(logfc.loc[others].mean())))
    worst = min(means, key=lambda item: item[1] if full >= 0 else -item[1])
    flip = any((full > 0 and m <= 0) or (full < 0 and m >= 0) or full == 0 for _, m in means)
    return {
        "loo_min_mean": min(m for _, m in means),
        "loo_max_mean": max(m for _, m in means),
        "loo_worst_gene": worst[0],
        "loo_sign_flip": "yes" if flip else "no",
        "loo_n": len(means),
    }


def json_safe(obj):
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, dict):
        return {k: json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [json_safe(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        return json_safe(float(obj))
    if isinstance(obj, (np.integer,)):
        return int(obj)
    return obj


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    df = load_edger(download())
    anchor = float(df.loc["Cldn4", "logFC"])
    if not (-6.2 < anchor < -5.9):
        raise SystemExit(f"Cldn4 logFC {anchor} is not KO minus WT")
    rank = df["logFC"].sort_values(ascending=False)
    upper = {}
    for symbol in df.index:
        upper.setdefault(symbol.upper(), symbol)

    hallmark = read_gmt(SETDIR / "mh.all.v2023.2.Mm.symbols.gmt")
    kegg = read_gmt(SETDIR / "kegg_mouse_selected.gmt")

    catalog: list[tuple[str, str, list[str]]] = []
    for name, cls in HALLMARK_CLASS.items():
        catalog.append((name, cls, hallmark[name]))
    for name, cls in KEGG_CLASS.items():
        catalog.append((name, cls, kegg[name]))
    for name, (cls, genes) in CUSTOM.items():
        catalog.append((name, cls, genes))

    ifn_sources = [
        hallmark["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        hallmark["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        CUSTOM["custom_phospho_proxy_ISG"][1],
    ]
    ifn_genes: set[str] = set()
    for src in ifn_sources:
        resolved, _ = resolve_members(src, upper)
        ifn_genes.update(resolved)

    member_rows = []
    resolved_sets: dict[str, list[str]] = {}
    classes: dict[str, str] = {}
    for name, cls, genes in catalog:
        found, missing = resolve_members(genes, upper)
        resolved_sets[name] = found
        classes[name] = cls
        for gene in found:
            member_rows.append(
                {
                    "set": name,
                    "set_class": cls,
                    "symbol": gene,
                    "in_rank": "yes",
                    "logFC": float(df.loc[gene, "logFC"]),
                    "logCPM": float(df.loc[gene, "logCPM"]),
                    "FDR": float(df.loc[gene, "FDR"]),
                    "ifn_overlap": "yes" if gene in ifn_genes else "no",
                    "expressed_logCPM_ge_0": "yes" if float(df.loc[gene, "logCPM"]) >= MIN_LOGCPM else "no",
                }
            )
        for gene in missing:
            member_rows.append(
                {
                    "set": name,
                    "set_class": cls,
                    "symbol": gene,
                    "in_rank": "no",
                    "logFC": math.nan,
                    "logCPM": math.nan,
                    "FDR": math.nan,
                    "ifn_overlap": "",
                    "expressed_logCPM_ge_0": "",
                }
            )

    # Reproduction lock: IFN-γ NES on this rank is the number from the earlier pass.
    ifn_only = gsea_prerank(
        rank,
        {"HALLMARK_INTERFERON_GAMMA_RESPONSE": resolved_sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"]},
        nperm=1000,
        seed=SEED,
    )
    ifn_nes = float(ifn_only.loc[0, "nes"])
    if abs(ifn_nes - 1.508101) > 0.02:
        raise SystemExit(f"IFN-γ NES {ifn_nes} drifted from the locked 1.508")

    gsea = gsea_prerank(rank, resolved_sets, nperm=1000, seed=SEED)
    gsea_by_term = {row.term: row for row in gsea.itertuples(index=False)} if len(gsea) else {}

    decont_sets = {
        name: [g for g in genes if g not in ifn_genes]
        for name, genes in resolved_sets.items()
        if classes[name] != "output"
    }
    gsea_de = gsea_prerank(rank, decont_sets, nperm=1000, seed=SEED)
    gsea_de_by = {row.term: row for row in gsea_de.itertuples(index=False)} if len(gsea_de) else {}

    universe = rank.to_numpy(dtype=float)
    summary_rows = []
    loo_rows = []
    for name, genes in resolved_sets.items():
        expressed = [g for g in genes if float(df.loc[g, "logCPM"]) >= MIN_LOGCPM]
        values = df.loc[expressed, "logFC"].to_numpy(dtype=float) if expressed else np.array([])
        de_genes = [g for g in expressed if g not in ifn_genes] if classes[name] != "output" else list(expressed)
        # Output sets are the IFN result; decontamination would empty them.
        if classes[name] == "output":
            de_genes = list(expressed)
        de_values = df.loc[de_genes, "logFC"].to_numpy(dtype=float) if de_genes else np.array([])
        seed_i = SEED + (sum(ord(ch) for ch in name) % 100000)
        full_p = perm_p(values, universe, NPERM_MEAN, seed_i)
        de_p = perm_p(de_values, universe, NPERM_MEAN, seed_i + 17)
        loo = leave_one_out(de_genes, df["logFC"])
        grow = gsea_by_term.get(name)
        drow = gsea_de_by.get(name)
        n_ifn = sum(1 for g in expressed if g in ifn_genes)
        summary_rows.append(
            {
                "set": name,
                "set_class": classes[name],
                "bridge_eligible": "yes" if classes[name] in BRIDGE_CLASSES else "no",
                "n_in_rank": len(genes),
                "n_expressed": len(expressed),
                "n_ifn_overlap_expressed": n_ifn,
                "n_decontaminated_expressed": len(de_genes),
                "mean_logFC": float(values.mean()) if len(values) else math.nan,
                "median_logFC": float(np.median(values)) if len(values) else math.nan,
                "n_up": int((values > 0).sum()) if len(values) else 0,
                "n_down": int((values < 0).sum()) if len(values) else 0,
                "mean_perm_p": full_p,
                "decontam_mean_logFC": float(de_values.mean()) if len(de_values) else math.nan,
                "decontam_median_logFC": float(np.median(de_values)) if len(de_values) else math.nan,
                "decontam_n_up": int((de_values > 0).sum()) if len(de_values) else 0,
                "decontam_n_down": int((de_values < 0).sum()) if len(de_values) else 0,
                "decontam_mean_perm_p": de_p,
                "nes": float(grow.nes) if grow is not None else math.nan,
                "nes_nom_p": float(grow.nom_p) if grow is not None else math.nan,
                "es": float(grow.es) if grow is not None else math.nan,
                "n_gsea": int(grow.n_set_in_rank) if grow is not None else 0,
                "lead_genes": grow.lead_genes if grow is not None else "",
                "decontam_nes": float(drow.nes) if drow is not None else math.nan,
                "decontam_nes_nom_p": float(drow.nom_p) if drow is not None else math.nan,
                "decontam_lead_genes": drow.lead_genes if drow is not None else "",
                "loo_min_mean": loo["loo_min_mean"],
                "loo_worst_gene": loo["loo_worst_gene"],
                "loo_sign_flip": loo["loo_sign_flip"],
            }
        )
        if de_genes:
            full_de_mean = float(df.loc[de_genes, "logFC"].mean())
            for gene in de_genes:
                others = [g for g in de_genes if g != gene]
                loo_mean = float(df.loc[others, "logFC"].mean()) if others else math.nan
                loo_rows.append(
                    {
                        "set": name,
                        "set_class": classes[name],
                        "left_out": gene,
                        "left_out_logFC": float(df.loc[gene, "logFC"]),
                        "left_out_FDR": float(df.loc[gene, "FDR"]),
                        "mean_without": loo_mean,
                        "delta_vs_full": loo_mean - full_de_mean,
                    }
                )

    summary = pd.DataFrame(summary_rows)
    # BH within the sweep, separately for full NES p and for bridge-eligible decontam NES p.
    summary["nes_fdr_in_sweep"] = bh_fdr(summary["nes_nom_p"])
    bridge_mask = summary["bridge_eligible"].eq("yes")
    summary["decontam_nes_fdr_bridge_family"] = np.nan
    summary.loc[bridge_mask, "decontam_nes_fdr_bridge_family"] = bh_fdr(
        summary.loc[bridge_mask, "decontam_nes_nom_p"]
    ).to_numpy()

    def rank_key(frame: pd.DataFrame) -> pd.Series:
        # Prefer a positive decontaminated NES, then its nominal p, then the mean.
        nes = frame["decontam_nes"]
        p = frame["decontam_nes_nom_p"]
        mean = frame["decontam_mean_logFC"]
        score = np.where(nes.notna(), nes, mean.fillna(-999))
        # Sort later; this column is the score (higher is a stronger KO-up enrichment).
        return pd.Series(score, index=frame.index)

    summary["bridge_score"] = rank_key(summary)
    eligible = summary[summary["bridge_eligible"].eq("yes")].copy()
    eligible = eligible.sort_values(
        ["decontam_nes", "decontam_mean_logFC"],
        ascending=[False, False],
        na_position="last",
    )
    prolif = summary[summary["set_class"].eq("proliferation")].copy()
    prolif = prolif.sort_values(["nes", "mean_logFC"], ascending=[False, False], na_position="last")

    def top_stable(frame: pd.DataFrame) -> dict | None:
        if frame.empty:
            return None
        stable = frame[frame["loo_sign_flip"].ne("yes") & frame["decontam_mean_logFC"].gt(0)]
        pick = stable.iloc[0] if len(stable) else frame.iloc[0]
        return pick.to_dict()

    best_bridge = top_stable(eligible)
    best_prolif = top_stable(prolif)
    # A bridge is called only when the decontaminated enrichment is KO-up,
    # leave-one-out does not flip the sign, and the bridge-family BH FDR is < 0.05.
    called = False
    if best_bridge is not None:
        fdr = best_bridge.get("decontam_nes_fdr_bridge_family")
        nes = best_bridge.get("decontam_nes")
        called = (
            best_bridge.get("loo_sign_flip") == "no"
            and isinstance(best_bridge.get("decontam_mean_logFC"), float)
            and best_bridge["decontam_mean_logFC"] > 0
            and isinstance(fdr, float)
            and not math.isnan(fdr)
            and fdr < 0.05
            and isinstance(nes, float)
            and nes > 0
        )

    blocked = [
        {
            "method": "GSVA",
            "ran": "no",
            "reason": "Needs a genes-by-sample matrix. GSE50927 series matrix Sample_data_row_count is 0 for all 5 GSM. Deposited object is one EdgeR logFC per gene, n=1 vs 1.",
        },
        {
            "method": "ssGSEA",
            "ran": "no",
            "reason": "ssGSEA scores each sample. This contrast has no sample expression column. A score built from the logFC rank would be preranked GSEA with a different name.",
        },
        {
            "method": "preranked GSEA",
            "ran": "yes",
            "reason": "Subramanian 2005 weighted KS on the deposited logFC rank. 1000 gene-set permutations, seed 42, independent stream per set. Positive NES = KO-up. Gene-set permutation, not a mouse-level test.",
        },
    ]

    summary = summary.sort_values(["set_class", "set"])
    summary.to_csv(TABLES / "sweep_summary.tsv", sep="\t", index=False)
    pd.DataFrame(member_rows).to_csv(TABLES / "sweep_set_members.tsv", sep="\t", index=False)
    pd.DataFrame(loo_rows).to_csv(TABLES / "sweep_leave_one_out.tsv", sep="\t", index=False)
    pd.DataFrame(blocked).to_csv(TABLES / "sweep_methods.tsv", sep="\t", index=False)
    eligible.to_csv(TABLES / "sweep_bridge_rank.tsv", sep="\t", index=False)

    call = {
        "contrast": "naive Cldn4 KO vs WT, no VILI",
        "n": "1 vs 1",
        "logFC": "KO minus WT",
        "Cldn4_logFC": anchor,
        "ifn_gamma_nes_lock": ifn_nes,
        "n_ifn_genes_used_for_decontamination": len(ifn_genes),
        "gsva": "not run",
        "ssgsea": "not run",
        "bridge_called_at_fdr_0.05": called,
        "best_non_proliferation_set": None if best_bridge is None else best_bridge["set"],
        "best_non_proliferation_decontam_nes": None if best_bridge is None else best_bridge["decontam_nes"],
        "best_non_proliferation_decontam_fdr": None
        if best_bridge is None
        else best_bridge["decontam_nes_fdr_bridge_family"],
        "best_non_proliferation_decontam_mean": None if best_bridge is None else best_bridge["decontam_mean_logFC"],
        "best_non_proliferation_loo_flip": None if best_bridge is None else best_bridge["loo_sign_flip"],
        "best_proliferation_set": None if best_prolif is None else best_prolif["set"],
        "best_proliferation_nes": None if best_prolif is None else best_prolif["nes"],
        "best_proliferation_mean": None if best_prolif is None else best_prolif["mean_logFC"],
        "rule": (
            "Upstream bridge candidates are repair, damage_signal, sensor, and "
            "nuclear-envelope sets after Hallmark IFN-α/γ and the phospho-proxy "
            "ISG list are removed. A call requires decontaminated NES > 0, "
            "leave-one-out mean staying positive, and BH FDR < 0.05 inside that "
            "family. Proliferation sets are reported beside the bridge, not as the bridge."
        ),
    }
    (TABLES / "sweep_call.json").write_text(json.dumps(json_safe(call), indent=2) + "\n")

    plot_sweep(summary, FIGURES / "fig2_sweep_nes.png", FIGURES / "fig2_sweep_nes.pdf")
    print(json.dumps(json_safe(call), indent=2))
    show = summary.sort_values("nes", ascending=False, na_position="last")
    cols = [
        "set",
        "set_class",
        "n_expressed",
        "n_ifn_overlap_expressed",
        "mean_logFC",
        "nes",
        "nes_nom_p",
        "nes_fdr_in_sweep",
        "decontam_mean_logFC",
        "decontam_nes",
        "decontam_nes_nom_p",
        "decontam_nes_fdr_bridge_family",
        "loo_sign_flip",
        "loo_worst_gene",
    ]
    print(show[cols].to_string(index=False))


def plot_sweep(summary: pd.DataFrame, path_png: Path, path_pdf: Path) -> None:
    frame = summary.dropna(subset=["nes"]).sort_values("nes")
    colors = {
        "output": "#A33B3B",
        "proliferation": "#3C6E8F",
        "repair": "#5E6B73",
        "damage_signal": "#C47B2B",
        "sensor": "#2F6F4E",
        "structure": "#6B4C7A",
    }
    fig, ax = plt.subplots(figsize=(10.2, 8.2))
    y = np.arange(len(frame))
    bar_colors = [colors.get(c, "#888888") for c in frame["set_class"]]
    ax.barh(y, frame["nes"], color=bar_colors, height=0.72)
    ax.set_yticks(y)
    labels = []
    for row in frame.itertuples(index=False):
        label = row.set.replace("HALLMARK_", "").replace("KEGG_", "KEGG ").replace("custom_", "")
        labels.append(label)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.axvline(0, color="#222222", linewidth=0.8)
    ax.set_xlabel("prerank NES  (positive = enriched in Cldn4 KO)")
    ax.set_title("GSE50927 naive lung, no VILI  ·  n = 1 vs 1  ·  not GSVA / ssGSEA", loc="left", fontsize=11)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    handles = [
        plt.Rectangle((0, 0), 1, 1, color=colors[k])
        for k in ("output", "proliferation", "repair", "damage_signal", "sensor", "structure")
    ]
    ax.legend(
        handles,
        ["IFN / phospho-proxy output", "proliferation", "repair / NHEJ", "damage signal", "cytosolic DNA", "nuclear envelope"],
        frameon=False,
        fontsize=8,
        loc="lower right",
    )
    fig.tight_layout()
    fig.savefig(path_png, dpi=160)
    fig.savefig(path_pdf)
    plt.close(fig)


if __name__ == "__main__":
    main()
