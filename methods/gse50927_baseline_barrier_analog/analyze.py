#!/usr/bin/env python3
"""GSE50927 baseline Cldn4 KO vs WT (no VILI): IFN / chemokine / MHC.

Barrier-loss analog only. Not a cancer, tumour, or ICI contrast.

Primary contrast is the author edgeR table for naive whole lung
(Cldn4 KO minus WT, no ventilator). The three VILI tables are used only
to separate injury confounds. They are not the KD-match contrast.

Honest n is 1 vs 1 GEO sample per contrast. The series text says
duplicates; no per-sample count matrix is deposited. Author FDR is an
edgeR dispersion call on that unreplicated design, not a mouse-level test.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binomtest, wilcoxon

ROOT = Path(__file__).resolve().parent
DATA = Path(os.environ.get("GSE50927_BASELINE_DATA", "/tmp/gse50927_baseline_barrier"))
TABLES = ROOT / "tables"
FIGURES = ROOT / "figures"

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl"
FILES = {
    "naive_KO_vs_WT": "GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "VILI_WT_vs_naive_WT": "GSE50927_VILIwtGenes.csv.gz",
    "VILI_KOhigh_vs_WT": "GSE50927_VILIwtkohiGenes.csv.gz",
    "VILI_KOlow_vs_WT": "GSE50927_VILIwtkoloGenes.csv.gz",
}

# IFN/ISG list is the public panel already used for this accession, with
# Cxcl9/Cxcl10/Cxcl11 moved into the chemokine module so the two calls
# do not double-count. MHC list is the same public panel, then split
# below into classical MHC-I, MHC-II, and haplotype-risk bins.
IFN_GENES = [
    "Ifnb1", "Ifna1", "Ifna2", "Ifna4", "Ifna5", "Ifng", "Ifnar1", "Ifnar2",
    "Stat1", "Stat2", "Stat3", "Irf1", "Irf3", "Irf7", "Irf9",
    "Isg15", "Isg20", "Mx1", "Mx2",
    "Oas1a", "Oas1b", "Oas1g", "Oas2", "Oas3", "Oasl1", "Oasl2",
    "Ifit1", "Ifit2", "Ifit3", "Ifitm1", "Ifitm2", "Ifitm3", "Rsad2",
    "Gbp2", "Gbp3", "Gbp4", "Gbp5", "Gbp6", "Gbp7",
    "Ifi44", "Ifi204", "Ifi27", "Ifi35", "Usp18", "Adar", "Dhx58",
    "Ddx58", "Ifih1", "Tlr3", "Tlr7", "Tlr9", "Nlrc5",
    "Ifi47", "Ubd", "Irg1", "Irgm1", "Irgm2", "Igtp", "Tgtp1", "Tgtp2",
    "Iigp1", "Gbp2b",
]
CHEMOKINE_GENES = [
    "Ccl1", "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
    "Ccl11", "Ccl12", "Ccl17", "Ccl19", "Ccl20", "Ccl22", "Ccl24", "Ccl25",
    "Ccl27a", "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10", "Cxcl11",
    "Cxcl12", "Cxcl13", "Cxcl14", "Cxcl16", "Cxcl17", "Cx3cl1", "Xcl1",
]
MHC_CLASSICAL = [
    "B2m", "H2-K1", "H2-D1", "Tap1", "Tap2", "Tapbp",
    "Psmb8", "Psmb9", "Psmb10", "Nlrc5",
]
MHC_II = [
    "Ciita", "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "H2-Ea",
    "H2-DMb1", "H2-DMa", "H2-Oa", "H2-Ob", "H2-DMb2",
]
MHC_HAPLOTYPE = [
    "H2-K2", "H2-Q1", "H2-Q2", "H2-Q4", "H2-Q6", "H2-Q7", "H2-Q8",
    "H2-Q10", "H2-T23", "H2-M2", "H2-M3", "H2-Bl",
]
INJURY_MEDIATORS = ["Tnf", "Il1b", "Il6", "Egr1"]

PANELS = {
    "IFN_ISG": IFN_GENES,
    "chemokine": CHEMOKINE_GENES,
    "MHC_I_classical": MHC_CLASSICAL,
    "MHC_II": MHC_II,
    "MHC_haplotype_risk": MHC_HAPLOTYPE,
    "injury_mediator": INJURY_MEDIATORS,
}

# Primary hit rule. logCPM > 0 drops the edgeR floor (about -2.07) and
# genes that are called only because a near-zero count jumped.
FDR_CUT = 0.05
CPM_CUT = 0.0
LOW_CPM = 1.0  # flagged, still eligible for a hit if logCPM > 0


def download() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for name in FILES.values():
        dest = DATA / name
        if dest.exists() and dest.stat().st_size > 1000:
            continue
        url = f"{FTP}/{name}"
        print(f"GET {url}")
        urllib.request.urlretrieve(url, dest)


def load_table(path: Path) -> dict[str, dict]:
    with gzip.open(path, "rt", newline="") as handle:
        rows = list(csv.DictReader(handle))
    by_symbol: dict[str, dict] = {}
    for row in rows:
        symbol = (row.get("Marker.Symbol") or row.get("GeneSymbol") or "").strip()
        if not symbol or symbol in by_symbol:
            continue
        by_symbol[symbol] = {
            "entrez": row["EntrezID"],
            "logFC": float(row["logFC"]),
            "logCPM": float(row["logCPM"]),
            "P": float(row["PValue"]),
            "FDR": float(row["FDR"]),
        }
    return by_symbol


def hit(rec: dict, direction: str) -> bool:
    if rec["logCPM"] <= CPM_CUT or rec["FDR"] >= FDR_CUT:
        return False
    if direction == "up":
        return rec["logFC"] > 0
    return rec["logFC"] < 0


def classify(base: dict, vili: dict, module: str) -> str:
    if module == "MHC_haplotype_risk":
        return "haplotype_risk_excluded"
    base_up = hit(base, "up")
    base_down = hit(base, "down")
    vili_up = hit(vili, "up")
    base_up_low = (
        base["logFC"] > 0
        and base["FDR"] < FDR_CUT
        and base["logCPM"] <= CPM_CUT
    )
    # A baseline decrease that flips up only after ventilation is an injury
    # confound with the opposite resting sign (Egr1). Do not file it as
    # injury-only, which would hide the baseline decrease.
    if base_down and vili_up:
        return "baseline_down_and_injury_up"
    if base_up and vili_up:
        return "injury_shared"
    if base_up and not vili_up:
        return "baseline_only"
    if base_up_low and vili_up:
        return "injury_shared_low_baseline_cpm"
    if base_up_low and not vili_up:
        return "baseline_low_cpm_only"
    if vili_up and not base_up:
        return "injury_only"
    if base_down:
        return "baseline_down"
    return "not_hit"


def gene_set_test(logfcs: list[float]) -> dict:
    """Gene-wise direction test. Not a mouse-level p."""
    usable = [v for v in logfcs if v != 0 and math.isfinite(v)]
    n_up = sum(1 for v in usable if v > 0)
    out = {
        "n_detected": len(logfcs),
        "n_nonzero": len(usable),
        "n_up": n_up,
        "median_logFC": median(logfcs) if logfcs else float("nan"),
        "wilcoxon_stat": "",
        "wilcoxon_p": "",
        "binom_p": "",
    }
    if len(usable) >= 5:
        stat, p = wilcoxon(usable, alternative="two-sided", zero_method="wilcox")
        out["wilcoxon_stat"] = float(stat)
        out["wilcoxon_p"] = float(p)
        out["binom_p"] = float(binomtest(n_up, len(usable), 0.5, alternative="two-sided").pvalue)
    return out


def median(vals: list[float]) -> float:
    s = sorted(vals)
    n = len(s)
    if n == 0:
        return float("nan")
    mid = n // 2
    if n % 2:
        return s[mid]
    return 0.5 * (s[mid - 1] + s[mid])


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def plot(genes: list[dict]) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    colors = {
        "baseline_only": "#1f4e79",
        "injury_shared": "#c47b09",
        "injury_shared_low_baseline_cpm": "#e0a100",
        "injury_only": "#a93226",
        "baseline_down": "#5d6d7e",
        "baseline_down_and_injury_up": "#117a65",
        "baseline_low_cpm_only": "#85929e",
        "not_hit": "#d5d8dc",
        "haplotype_risk_excluded": "#6c3483",
    }
    focus = [
        g for g in genes
        if g["module"] in {"IFN_ISG", "chemokine", "MHC_I_classical"}
        and g["class"] in {
            "baseline_only", "injury_shared", "injury_shared_low_baseline_cpm", "injury_only",
        }
        and (
            g["class"] != "baseline_only" or float(g["logCPM_baseline"]) >= LOW_CPM
            or g["module"] != "IFN_ISG"
        )
    ]
    # Keep the scatter readable: IFN baseline-only is large, so plot every
    # non-grey class but only label the chemokines, MHC, and a few IFN genes.
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 5.6), gridspec_kw={"width_ratios": [1.05, 1.15]})

    # Left: class counts
    modules = ["IFN_ISG", "chemokine", "MHC_I_classical"]
    module_labels = ["IFN / ISG", "Chemokine", "MHC-I classical"]
    class_order = ["baseline_only", "injury_shared", "injury_only"]
    class_labels = ["baseline only", "injury-shared", "injury only"]
    ax = axes[0]
    import numpy as np
    x = np.arange(len(modules))
    width = 0.24
    for i, (cls, lab) in enumerate(zip(class_order, class_labels)):
        counts = []
        for module in modules:
            counts.append(sum(1 for g in genes if g["module"] == module and g["class"] == cls))
        # injury_shared bar also shows the low-baseline-CPM shared genes as a hatch overlay? 
        # Keep the primary class only; low-CPM shared is annotated in the finding.
        ax.bar(x + (i - 1) * width, counts, width, label=lab, color=colors[cls], zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(module_labels)
    ax.set_ylabel("Genes (author FDR < 0.05, logCPM > 0)")
    ax.set_title("Baseline KO call, split from VILI")
    ax.legend(frameon=False, fontsize=8)
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, color="#eeeeee")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Right: baseline logFC vs WT-VILI logFC
    ax = axes[1]
    # Named points are the decision genes. The remaining points stay unlabeled
    # so the injury axis is readable. Full names are in tables/gene_level.tsv.
    label_genes = {
        "Ccl5", "Cx3cl1", "Cxcl1", "Ccl2", "B2m", "H2-K1", "Isg15", "Egr1", "Il6",
    }
    offsets = {
        "Il6": (8, 4),
        "Ccl2": (8, -2),
        "Cxcl1": (-34, 8),
        "Ccl5": (8, 4),
        "Cx3cl1": (-38, 6),
        "Isg15": (8, -12),
        "B2m": (-26, -12),
        "H2-K1": (-36, -2),
        "Egr1": (8, 6),
    }
    plotted = [
        g for g in genes
        if g["module"] in {"IFN_ISG", "chemokine", "MHC_I_classical", "injury_mediator"}
        and g["in_table"] == "yes"
        and (float(g["logCPM_baseline"]) > CPM_CUT or float(g["logCPM_vili"]) > CPM_CUT)
        and (g["class"] != "not_hit" or g["gene"] in label_genes)
    ]
    for g in plotted:
        ax.scatter(
            float(g["logFC_baseline"]),
            float(g["logFC_VILI_WT"]),
            s=36 if g["gene"] in label_genes else 16,
            color=colors.get(g["class"], "#bbbbbb"),
            zorder=3 if g["gene"] in label_genes else 2,
            linewidths=0,
        )
    for g in plotted:
        if g["gene"] not in label_genes:
            continue
        ax.annotate(
            g["gene"],
            (float(g["logFC_baseline"]), float(g["logFC_VILI_WT"])),
            textcoords="offset points",
            xytext=offsets.get(g["gene"], (4, 3)),
            fontsize=8,
            color="#1c2833",
            arrowprops={"arrowstyle": "-", "color": "#7f8c8d", "lw": 0.6},
        )
    ax.axvline(0, color="#bbbbbb", lw=0.6)
    ax.axhline(0, color="#bbbbbb", lw=0.6)
    ax.set_xlabel("Baseline logFC (KO − WT, no VILI)")
    ax.set_ylabel("Injury logFC (WT VILI − naive WT)")
    ax.set_title("Same gene, two contrasts")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.suptitle(
        "GSE50927 barrier-loss analog (not cancer): naive Cldn4 KO lung",
        fontsize=12,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(FIGURES / "fig_baseline_vs_injury.png", dpi=160, bbox_inches="tight")
    fig.savefig(FIGURES / "fig_baseline_vs_injury.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    download()
    tables = {key: load_table(DATA / name) for key, name in FILES.items()}
    naive = tables["naive_KO_vs_WT"]
    cldn4 = naive["Cldn4"]
    if abs(cldn4["logFC"] + 6.061296179) > 1e-6:
        raise SystemExit(f"Cldn4 logFC drifted: {cldn4['logFC']}")
    if cldn4["FDR"] > 1e-20:
        raise SystemExit(f"Cldn4 FDR unexpected: {cldn4['FDR']}")

    gene_rows = []
    for module, symbols in PANELS.items():
        for symbol in symbols:
            base = naive.get(symbol)
            vili = tables["VILI_WT_vs_naive_WT"].get(symbol)
            kohi = tables["VILI_KOhigh_vs_WT"].get(symbol)
            kolo = tables["VILI_KOlow_vs_WT"].get(symbol)
            if base is None or vili is None or kohi is None or kolo is None:
                gene_rows.append({
                    "gene": symbol,
                    "module": module,
                    "in_table": "no",
                    "entrez": "",
                    "logFC_baseline": "",
                    "logCPM_baseline": "",
                    "P_baseline": "",
                    "FDR_baseline": "",
                    "logFC_VILI_WT": "",
                    "FDR_VILI_WT": "",
                    "logCPM_vili": "",
                    "logFC_KOhigh": "",
                    "FDR_KOhigh": "",
                    "logFC_KOlow": "",
                    "FDR_KOlow": "",
                    "class": "absent",
                    "low_abundance_baseline": "",
                })
                continue
            label = classify(base, vili, module)
            gene_rows.append({
                "gene": symbol,
                "module": module,
                "in_table": "yes",
                "entrez": base["entrez"],
                "logFC_baseline": f"{base['logFC']:.9g}",
                "logCPM_baseline": f"{base['logCPM']:.9g}",
                "P_baseline": f"{base['P']:.9g}",
                "FDR_baseline": f"{base['FDR']:.9g}",
                "logFC_VILI_WT": f"{vili['logFC']:.9g}",
                "FDR_VILI_WT": f"{vili['FDR']:.9g}",
                "logCPM_vili": f"{vili['logCPM']:.9g}",
                "logFC_KOhigh": f"{kohi['logFC']:.9g}",
                "FDR_KOhigh": f"{kohi['FDR']:.9g}",
                "logFC_KOlow": f"{kolo['logFC']:.9g}",
                "FDR_KOlow": f"{kolo['FDR']:.9g}",
                "class": label,
                "low_abundance_baseline": "yes" if base["logCPM"] < LOW_CPM else "no",
            })

    # Nlrc5 is on both IFN and classical MHC. Keep both rows (module membership)
    # but the MHC call is the one used for the MHC sentence. IFN Wilcoxon
    # includes it once, inside IFN_ISG, which matches the prior IFN panel.

    summary_rows = []
    tests = {}
    for module in PANELS:
        members = [g for g in gene_rows if g["module"] == module and g["in_table"] == "yes"]
        detected = [g for g in members if float(g["logCPM_baseline"]) > CPM_CUT]
        counts = {}
        for g in members:
            counts[g["class"]] = counts.get(g["class"], 0) + 1
        all_test = gene_set_test([float(g["logFC_baseline"]) for g in detected])
        held_out = {
            "injury_shared",
            "injury_shared_low_baseline_cpm",
            "injury_only",
            "haplotype_risk_excluded",
        }
        clean = [g for g in detected if g["class"] not in held_out]
        clean_test = gene_set_test([float(g["logFC_baseline"]) for g in clean])
        tests[module] = {"all_detected": all_test, "excluding_injury_genes": clean_test, "counts": counts}
        summary_rows.append({
            "module": module,
            "n_panel": len(PANELS[module]),
            "n_in_table": len(members),
            "n_absent": sum(1 for g in gene_rows if g["module"] == module and g["in_table"] == "no"),
            "n_detected_logCPM_gt_0": len(detected),
            "n_baseline_only": counts.get("baseline_only", 0),
            "n_baseline_only_logCPM_ge_1": sum(
                1 for g in members
                if g["class"] == "baseline_only" and float(g["logCPM_baseline"]) >= LOW_CPM
            ),
            "n_injury_shared": counts.get("injury_shared", 0),
            "n_injury_shared_low_baseline_cpm": counts.get("injury_shared_low_baseline_cpm", 0),
            "n_injury_only": counts.get("injury_only", 0),
            "n_baseline_down": counts.get("baseline_down", 0),
            "n_baseline_down_and_injury_up": counts.get("baseline_down_and_injury_up", 0),
            "n_baseline_only_down_in_VILI": sum(
                1
                for g in members
                if g["class"] == "baseline_only"
                and float(g["logFC_VILI_WT"]) < 0
                and float(g["FDR_VILI_WT"]) < FDR_CUT
            ),
            "n_baseline_low_cpm_only": counts.get("baseline_low_cpm_only", 0),
            "n_haplotype_risk_excluded": counts.get("haplotype_risk_excluded", 0),
            "median_logFC_detected": all_test["median_logFC"],
            "n_up_detected": all_test["n_up"],
            "n_nonzero_detected": all_test["n_nonzero"],
            "wilcoxon_p_detected": all_test["wilcoxon_p"],
            "binom_p_detected": all_test["binom_p"],
            "median_logFC_excluding_injury": clean_test["median_logFC"],
            "wilcoxon_p_excluding_injury": clean_test["wilcoxon_p"],
            "role": (
                "haplotype bin, not a match numerator"
                if module == "MHC_haplotype_risk"
                else "injury confound control, not a match module"
                if module == "injury_mediator"
                else "baseline KO vs WT match module"
            ),
        })

    fields = list(gene_rows[0].keys())
    write_tsv(TABLES / "gene_level.tsv", gene_rows, fields)
    write_tsv(TABLES / "module_summary.tsv", summary_rows, list(summary_rows[0].keys()))

    qc_rows = []
    for key, table in tables.items():
        rec = table["Cldn4"]
        qc_rows.append({
            "contrast": key,
            "role": "primary_baseline" if key == "naive_KO_vs_WT" else "injury_confound_not_the_analog",
            "gene": "Cldn4",
            "entrez": rec["entrez"],
            "logFC": f"{rec['logFC']:.9g}",
            "logCPM": f"{rec['logCPM']:.9g}",
            "P": f"{rec['P']:.9g}",
            "FDR": f"{rec['FDR']:.9g}",
            "sign": "KO minus WT" if key != "VILI_WT_vs_naive_WT" else "VILI minus naive WT",
        })
    write_tsv(TABLES / "cldn4_qc.tsv", qc_rows, list(qc_rows[0].keys()))

    inventory = [{
        "dataset": "GSE50927",
        "label": "barrier-loss analog only",
        "cancer": "no",
        "ICI": "no",
        "tissue": "mouse whole lung",
        "primary_contrast": "naive Cldn4 KO vs WT (no VILI)",
        "honest_n": "1 vs 1 GSM",
        "design_text": "duplicates claimed; not deposited as two processed columns",
        "n_edgeR_tables": 4,
        "n_expression_matrix_rows": 0,
        "Cldn4_logFC_baseline": f"{cldn4['logFC']:.9g}",
        "Cldn4_FDR_baseline": f"{cldn4['FDR']:.9g}",
        "source": "author edgeR CSV, GEO supplementary",
    }]
    write_tsv(TABLES / "label_inventory.tsv", inventory, list(inventory[0].keys()))

    payload = {
        "label": "barrier-loss analog only",
        "not": ["cancer", "tumour", "LUAD", "ICI", "CLDN4-high tumor geography"],
        "primary_contrast": "naive_KO_vs_WT",
        "honest_n": "1 vs 1 GSM",
        "Cldn4_baseline": cldn4,
        "modules": tests,
        "rules": {
            "hit": "FDR < 0.05 and logCPM > 0 and logFC in the stated direction",
            "baseline_only": "up at naive KO vs WT, not up at FDR < 0.05 in WT VILI vs naive WT",
            "injury_shared": "up in both the baseline genotype contrast and WT VILI",
            "injury_only": "up in WT VILI, not a baseline hit",
            "haplotype_risk": "nonclassical H2 genes on the mixed 129/B6/BALB line; excluded from the MHC-I numerator",
            "gene_set_p": "Wilcoxon / binomial on gene logFC. Genes are not independent. Not a mouse-level p.",
        },
    }
    (TABLES / "summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    plot(gene_rows)

    print("Cldn4 baseline", cldn4["logFC"], cldn4["FDR"])
    for row in summary_rows:
        print(
            row["module"],
            "base", row["n_baseline_only"],
            "base>=1", row["n_baseline_only_logCPM_ge_1"],
            "shared", row["n_injury_shared"],
            "shared_low", row["n_injury_shared_low_baseline_cpm"],
            "injury_only", row["n_injury_only"],
            "down", row["n_baseline_down"],
            "low_only", row["n_baseline_low_cpm_only"],
            "median", row["median_logFC_detected"],
            "W", row["wilcoxon_p_detected"],
            "W_clean", row["wilcoxon_p_excluding_injury"],
        )
    print("wrote", TABLES, FIGURES)


if __name__ == "__main__":
    main()
