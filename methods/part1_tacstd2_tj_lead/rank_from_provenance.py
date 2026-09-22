#!/usr/bin/env python3
"""PART 1 polish — Tacstd2-high → TJ/junction lead enrichment ranks.

Re-ranks frozen tables from PR #741 (concordant-4) and PR #736
(OncoSG + GSE31210). Does not re-download matrices or invent NES/FDR.
No CLDN4 pin section.
"""

from __future__ import annotations

import csv
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
PROV = HERE / "provenance"
TABLES = HERE / "tables"
FIGURES = HERE / "figures"

TJ_STRICT = re.compile(
    r"TIGHT_JUNCTION|APICAL_JUNCTION|BICELLULAR_TIGHT|"
    r"JUNCTION_ORGANIZATION|APICAL_SURFACE|CUSTOM_EPITHELIAL_ADHESION|"
    r"CUSTOM_CLAUDIN",
    re.I,
)
BARRIER = re.compile(
    r"TIGHT_JUNCTION|APICAL_JUNCTION|BICELLULAR_TIGHT|JUNCTION_ORGANIZATION|"
    r"APICAL_SURFACE|KERATIN|SKIN_BARRIER|EPIDERMAL_CELL|CLAUDIN|ADHESION|"
    r"CORNIFICATION|KRT_EPITHELIAL|CUSTOM_EPITHELIAL|CUSTOM_CLAUDIN|"
    r"CUSTOM_TJ|EPITHELIAL_ADHESION",
    re.I,
)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fnum(x: str | None) -> float:
    if x is None or x == "":
        return float("nan")
    return float(x)


def rank_positive(rows: list[dict], nes_key: str = "nes") -> list[dict]:
    pos = [r for r in rows if fnum(r[nes_key]) > 0]
    pos = sorted(pos, key=lambda r: -fnum(r[nes_key]))
    for i, r in enumerate(pos, 1):
        r = dict(r)
        r["_rank_pos"] = i
        r["_n_pos"] = len(pos)
        pos[i - 1] = r
    return pos


def best_match(ranked: list[dict], pattern: re.Pattern, term_key: str = "term") -> dict | None:
    hits = [r for r in ranked if pattern.search(r.get(term_key, "") or "")]
    return hits[0] if hits else None


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)

    # ---- Concordant-4 ORA (PR #741) ----
    ora = read_tsv(PROV / "c4_ora_up.tsv")
    for i, r in enumerate(ora, 1):
        r["_rank"] = i
        r["_n"] = len(ora)
    c4_ora_aj = next(r for r in ora if r["term"] == "HALLMARK_APICAL_JUNCTION")
    c4_ora_tj = next(r for r in ora if r["term"] == "GOBP_TIGHT_JUNCTION_ORGANIZATION")
    c4_ora_adh = next(r for r in ora if r["term"] == "CUSTOM_EPITHELIAL_ADHESION")
    c4_ora_krt = next(r for r in ora if r["term"] == "GOBP_KERATINIZATION")

    # ---- Concordant-4 GSEA (PR #741) ----
    c4_gsea = read_tsv(PROV / "c4_gsea_prerank_q4q1.tsv")
    c4_pos = rank_positive(c4_gsea, "nes")
    c4_gsea_aj = next(r for r in c4_pos if r["term"] == "HALLMARK_APICAL_JUNCTION")
    c4_gsea_adh = next(r for r in c4_pos if r["term"] == "CUSTOM_EPITHELIAL_ADHESION")
    c4_gsea_tjorg = next(r for r in c4_pos if r["term"] == "GOBP_TIGHT_JUNCTION_ORGANIZATION")
    c4_gsea_kegg = next(r for r in c4_pos if r["term"] == "KEGG_TIGHT_JUNCTION")
    c4_gsea_krt = next(r for r in c4_pos if r["term"] == "GOBP_KERATINIZATION")
    c4_strict = [r for r in c4_pos if TJ_STRICT.search(r["term"])]
    for i, r in enumerate(c4_strict, 1):
        r["_rank_tj_family"] = i

    # ---- Family scores (PR #741) ----
    fam = read_tsv(PROV / "c4_family_score_q4q1.tsv")

    # ---- Bulk LUAD GSEA all (PR #736) ----
    luad_all = read_tsv(PROV / "luad_gsea_all.tsv")
    by_cohort: dict[str, list[dict]] = {}
    for r in luad_all:
        by_cohort.setdefault(r["cohort"], []).append(r)

    bulk_summaries = {}
    for cohort, items in by_cohort.items():
        pos = rank_positive(items, "nes")
        strict = [r for r in pos if TJ_STRICT.search(r["term"])]
        for i, r in enumerate(strict, 1):
            r["_rank_tj_family"] = i
        best_tj = strict[0] if strict else None
        kegg = next((r for r in items if r["term"] == "KEGG_TIGHT_JUNCTION"), None)
        aj = next((r for r in items if r["term"] == "HALLMARK_APICAL_JUNCTION"), None)
        bicel = next(
            (r for r in pos if r["term"] == "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY"),
            None,
        )
        skin = next(
            (r for r in pos if r["term"] == "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER"),
            None,
        )
        krt = next((r for r in pos if r["term"] == "GOBP_KERATINIZATION"), None)
        bulk_summaries[cohort] = {
            "n_all": len(items),
            "n_pos": len(pos),
            "best_tj": best_tj,
            "kegg_tj": kegg,
            "apical_junction": aj,
            "bicellular": bicel,
            "skin_barrier": skin,
            "keratinization": krt,
            "pos_top10": pos[:10],
            "strict_pos": strict,
        }

    # Headline ranks for Part 1 ending
    headlines = [
        {
            "layer": "concordant-4_ORA",
            "universe": f"A8+custom ORA UP terms (n={c4_ora_aj['_n']})",
            "best_tj_term": "HALLMARK_APICAL_JUNCTION",
            "best_tj_rank": int(c4_ora_aj["_rank"]),
            "n_universe_or_pos": int(c4_ora_aj["_n"]),
            "stat": "enrichment",
            "stat_value": fnum(c4_ora_aj["enrichment"]),
            "fdr": fnum(c4_ora_aj["fdr"]),
            "lead_claim": True,
            "note": "Apical junction is rank 1 among all tested ORA UP terms.",
            "source_pr": 741,
        },
        {
            "layer": "concordant-4_GSEA",
            "universe": f"positive-NES sets (n={c4_gsea_aj['_n_pos']})",
            "best_tj_term": "HALLMARK_APICAL_JUNCTION",
            "best_tj_rank": int(c4_gsea_aj["_rank_pos"]),
            "n_universe_or_pos": int(c4_gsea_aj["_n_pos"]),
            "stat": "NES",
            "stat_value": fnum(c4_gsea_aj["nes"]),
            "fdr": fnum(c4_gsea_aj["fdr"]),
            "lead_claim": False,
            "note": (
                f"Apical junction rank {c4_gsea_aj['_rank_pos']}; "
                f"epithelial adhesion rank {c4_gsea_adh['_rank_pos']}; "
                f"KEGG TJ rank {c4_gsea_kegg['_rank_pos']}. "
                f"Keratinization is NES-rank {c4_gsea_krt['_rank_pos']}."
            ),
            "source_pr": 741,
        },
        {
            "layer": "GSE31210_GSEA",
            "universe": (
                f"positive-NES sets (n={bulk_summaries['GSE31210_LUAD']['n_pos']} "
                f"of {bulk_summaries['GSE31210_LUAD']['n_all']})"
            ),
            "best_tj_term": "KEGG_TIGHT_JUNCTION",
            "best_tj_rank": int(
                bulk_summaries["GSE31210_LUAD"]["best_tj"]["_rank_pos"]
            ),
            "n_universe_or_pos": bulk_summaries["GSE31210_LUAD"]["n_pos"],
            "stat": "NES",
            "stat_value": fnum(bulk_summaries["GSE31210_LUAD"]["best_tj"]["nes"]),
            "fdr": fnum(bulk_summaries["GSE31210_LUAD"]["best_tj"]["fdr_bh_all"]),
            "lead_claim": True,
            "note": "KEGG Tight Junction is rank 3 among all positive-NES sets.",
            "source_pr": 736,
        },
        {
            "layer": "OncoSG_GSEA",
            "universe": (
                f"positive-NES sets (n={bulk_summaries['OncoSG_LUAD']['n_pos']} "
                f"of {bulk_summaries['OncoSG_LUAD']['n_all']})"
            ),
            "best_tj_term": "GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY",
            "best_tj_rank": int(
                bulk_summaries["OncoSG_LUAD"]["best_tj"]["_rank_pos"]
            ),
            "n_universe_or_pos": bulk_summaries["OncoSG_LUAD"]["n_pos"],
            "stat": "NES",
            "stat_value": fnum(bulk_summaries["OncoSG_LUAD"]["best_tj"]["nes"]),
            "fdr": fnum(bulk_summaries["OncoSG_LUAD"]["best_tj"]["fdr_bh_all"]),
            "lead_claim": False,
            "note": (
                "Best strict TJ term is bicellular TJ assembly (rank 16). "
                f"KEGG TJ NES={fnum(bulk_summaries['OncoSG_LUAD']['kegg_tj']['nes']):+.2f} "
                f"(not UP); Hallmark apical junction NES="
                f"{fnum(bulk_summaries['OncoSG_LUAD']['apical_junction']['nes']):+.2f} (DOWN). "
                "Skin-barrier / keratin lead OncoSG positives — report, do not hide."
            ),
            "source_pr": 736,
        },
    ]

    # Write headline TSV
    hl_path = TABLES / "best_tj_junction_ranks.tsv"
    fields = [
        "layer",
        "universe",
        "best_tj_term",
        "best_tj_rank",
        "n_universe_or_pos",
        "stat",
        "stat_value",
        "fdr",
        "lead_claim",
        "note",
        "source_pr",
    ]
    with hl_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for h in headlines:
            w.writerow(h)

    # Detailed C4 ORA barrier rows
    with (TABLES / "c4_ora_tj_barrier_ranks.tsv").open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "rank",
                "term",
                "enrichment",
                "fdr",
                "n_overlap",
                "family",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for r in ora:
            if BARRIER.search(r["term"]):
                fam_label = (
                    "TJ/junction"
                    if TJ_STRICT.search(r["term"])
                    else "barrier/keratin/adhesion"
                )
                w.writerow(
                    {
                        "rank": r["_rank"],
                        "term": r["term"],
                        "enrichment": r["enrichment"],
                        "fdr": r["fdr"],
                        "n_overlap": r["n_overlap"],
                        "family": fam_label,
                    }
                )

    # Detailed GSEA TJ ranks
    with (TABLES / "gsea_tj_junction_ranks.tsv").open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "cohort",
                "term",
                "rank_among_positive_NES",
                "n_positive",
                "NES",
                "FDR",
                "strict_TJ_family",
            ],
            delimiter="\t",
        )
        w.writeheader()
        for r in c4_pos:
            if TJ_STRICT.search(r["term"]) or r["term"] in {
                "GOBP_KERATINIZATION",
                "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER",
                "KRT_EPITHELIAL",
            }:
                w.writerow(
                    {
                        "cohort": "concordant-4",
                        "term": r["term"],
                        "rank_among_positive_NES": r["_rank_pos"],
                        "n_positive": r["_n_pos"],
                        "NES": r["nes"],
                        "FDR": r["fdr"],
                        "strict_TJ_family": bool(TJ_STRICT.search(r["term"])),
                    }
                )
        for cohort, summ in bulk_summaries.items():
            pos = rank_positive(by_cohort[cohort], "nes")
            for r in pos:
                if TJ_STRICT.search(r["term"]) or r["term"] in {
                    "GOBP_KERATINIZATION",
                    "GOBP_ESTABLISHMENT_OF_SKIN_BARRIER",
                    "KRT_EPITHELIAL",
                }:
                    fdr = r.get("fdr_bh_all") or r.get("fdr") or ""
                    w.writerow(
                        {
                            "cohort": cohort,
                            "term": r["term"],
                            "rank_among_positive_NES": r["_rank_pos"],
                            "n_positive": r["_n_pos"],
                            "NES": r["nes"],
                            "FDR": fdr,
                            "strict_TJ_family": bool(TJ_STRICT.search(r["term"])),
                        }
                    )

    # Family scores dump
    with (TABLES / "c4_family_scores.tsv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(fam[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(fam)

    summary = {
        "task": "PART1_polish_Tacstd2_high_TJ_junction_lead",
        "rule": "never fabricate; ranks recomputed from frozen PR #741/#736 tables only",
        "no_cldn4_pin": True,
        "user_ok_tj_first": True,
        "provenance": {
            "concordant4": "PR #741 methods/concordant4_tacstd2_malignant_deg",
            "bulk_luad": "PR #736 methods/paper_funnel_luad_tacstd2_tj",
        },
        "best_ranks": headlines,
        "part1_ending_verdict": (
            "Tacstd2-high → TJ/junction is the lead enrichment where the "
            "public human ranks support it: concordant-4 ORA apical junction "
            f"rank {c4_ora_aj['_rank']}/{c4_ora_aj['_n']}; GSE31210 KEGG TJ "
            f"rank {bulk_summaries['GSE31210_LUAD']['best_tj']['_rank_pos']}/"
            f"{bulk_summaries['GSE31210_LUAD']['n_pos']}. OncoSG does not put "
            "KEGG TJ / apical junction UP; its best strict TJ term is bicellular "
            f"TJ assembly rank {bulk_summaries['OncoSG_LUAD']['best_tj']['_rank_pos']}."
        ),
    }
    (TABLES / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # ---- Figure: best TJ ranks ----
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    # Panel A: C4 ORA top ranks (highlight TJ)
    ax = axes[0]
    top_n = 8
    terms = [r["term"].replace("HALLMARK_", "").replace("GOBP_", "").replace("CUSTOM_", "") for r in ora[:top_n]]
    enr = [fnum(r["enrichment"]) for r in ora[:top_n]]
    colors = [
        "#1b6ca8" if TJ_STRICT.search(ora[i]["term"]) else
        "#5b8c5a" if BARRIER.search(ora[i]["term"]) else "#9aa0a6"
        for i in range(top_n)
    ]
    y = list(range(top_n, 0, -1))
    ax.barh(y, enr, color=colors, edgecolor="none")
    ax.set_yticks(y)
    ax.set_yticklabels([f"#{i+1} {t}" for i, t in enumerate(terms)], fontsize=8)
    ax.set_xlabel("ORA enrichment (UP list)")
    ax.set_title("Concordant-4 TACSTD2 Q4 vs Q1\nORA — apical junction rank 1")
    ax.axvline(1, color="#ccc", lw=0.8)

    # Panel B: bulk best TJ ranks
    ax = axes[1]
    labels = []
    ranks = []
    cols = []
    annotations = []
    # GSE31210 KEGG TJ
    g = bulk_summaries["GSE31210_LUAD"]["best_tj"]
    labels.append("GSE31210\nKEGG TJ")
    ranks.append(int(g["_rank_pos"]))
    cols.append("#1b6ca8")
    annotations.append(f"#{g['_rank_pos']}/{g['_n_pos']}\nNES={fnum(g['nes']):+.2f}")
    # OncoSG bicellular
    o = bulk_summaries["OncoSG_LUAD"]["best_tj"]
    labels.append("OncoSG\nbicellular TJ")
    ranks.append(int(o["_rank_pos"]))
    cols.append("#c47b2b")
    annotations.append(f"#{o['_rank_pos']}/{o['_n_pos']}\nNES={fnum(o['nes']):+.2f}")
    # C4 GSEA apical junction
    labels.append("C4 GSEA\napical junction")
    ranks.append(int(c4_gsea_aj["_rank_pos"]))
    cols.append("#5b8c5a")
    annotations.append(
        f"#{c4_gsea_aj['_rank_pos']}/{c4_gsea_aj['_n_pos']}\nNES={fnum(c4_gsea_aj['nes']):+.2f}"
    )
    # C4 ORA as rank 1 reference
    labels.append("C4 ORA\napical junction")
    ranks.append(int(c4_ora_aj["_rank"]))
    cols.append("#1b6ca8")
    annotations.append(f"#{c4_ora_aj['_rank']}/{c4_ora_aj['_n']}\nenr={fnum(c4_ora_aj['enrichment']):.2f}")

    x = list(range(len(labels)))
    bars = ax.bar(x, ranks, color=cols, edgecolor="none")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Best TJ/junction rank (lower = lead)")
    ax.set_title("Best TJ/junction ranks\n(lower is better; blue = lead-claim layers)")
    ax.invert_yaxis()
    ax.set_ylim(max(ranks) + 3.5, 0)
    for bar, ann in zip(bars, annotations):
        # With inverted y, "toward lead" is smaller y — place labels above the bar top
        # in data coords (smaller than bar height).
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() - 0.35,
            ann,
            ha="center",
            va="top",
            fontsize=7,
            color="#333",
        )

    fig.suptitle(
        "PART 1 polish — Tacstd2-high → TJ/junction lead enrichment (honest ranks)",
        fontsize=11,
        fontweight="bold",
    )
    out = FIGURES / "fig_best_tj_junction_ranks"
    fig.savefig(f"{out}.png", dpi=160)
    fig.savefig(f"{out}.pdf")
    plt.close(fig)

    # Write FINDING.md
    write_finding(
        c4_ora_aj,
        c4_ora_tj,
        c4_ora_adh,
        c4_ora_krt,
        c4_gsea_aj,
        c4_gsea_adh,
        c4_gsea_tjorg,
        c4_gsea_kegg,
        c4_gsea_krt,
        fam,
        bulk_summaries,
        headlines,
        summary,
    )
    print(json.dumps({"wrote": str(hl_path), "verdict": summary["part1_ending_verdict"]}, indent=2))


def write_finding(
    c4_ora_aj,
    c4_ora_tj,
    c4_ora_adh,
    c4_ora_krt,
    c4_gsea_aj,
    c4_gsea_adh,
    c4_gsea_tjorg,
    c4_gsea_kegg,
    c4_gsea_krt,
    fam,
    bulk_summaries,
    headlines,
    summary,
) -> None:
    g31 = bulk_summaries["GSE31210_LUAD"]
    onc = bulk_summaries["OncoSG_LUAD"]
    g_best = g31["best_tj"]
    o_best = onc["best_tj"]

    def fam_row(name: str) -> dict:
        return next(r for r in fam if r["family"] == name)

    lines: list[str] = []
    lines.append("# PART 1 polish — Tacstd2-high → TJ/junction is the lead enrichment")
    lines.append("")
    lines.append(
        "Public human only: **concordant-4** malignant TACSTD2 Q4 vs Q1 "
        "(PR #741) and bulk LUAD **OncoSG + GSE31210** (PR #736). "
        "Ranks recomputed from those frozen tables. **Never fabricate.** "
        "**No CLDN4 pin section.** User OK with TJ first."
    )
    lines.append("")
    lines.append("## Part 1 ending (smooth)")
    lines.append("")
    lines.append(
        "Across the human public Tacstd2-high contrasts that were locked for "
        "this funnel, the enrichment that lands at the top of the rank list is "
        "**tight-junction / apical-junction**, not a generic metabolic or "
        "immune program:"
    )
    lines.append("")
    lines.append(
        f"1. **Concordant-4 ORA:** `HALLMARK_APICAL_JUNCTION` is **rank "
        f"{c4_ora_aj['_rank']} of {c4_ora_aj['_n']}** "
        f"(enrichment {fnum(c4_ora_aj['enrichment']):.2f}, "
        f"FDR {fnum(c4_ora_aj['fdr']):.2e})."
    )
    lines.append(
        f"2. **GSE31210 GSEA:** `KEGG_TIGHT_JUNCTION` is **rank "
        f"{g_best['_rank_pos']} of {g_best['_n_pos']}** positive-NES sets "
        f"(NES {fnum(g_best['nes']):+.2f}, FDR {fnum(g_best['fdr_bh_all']):.3g})."
    )
    lines.append(
        f"3. **OncoSG GSEA (honest):** KEGG TJ and Hallmark apical junction are "
        f"**not** UP (NES {fnum(onc['kegg_tj']['nes']):+.2f} / "
        f"{fnum(onc['apical_junction']['nes']):+.2f}). Best strict TJ term is "
        f"`GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY` at rank "
        f"**{o_best['_rank_pos']} of {o_best['_n_pos']}** "
        f"(NES {fnum(o_best['nes']):+.2f}, FDR {fnum(o_best['fdr_bh_all']):.3g}). "
        "Keratin / skin-barrier lead OncoSG positives — kept visible."
    )
    lines.append("")
    lines.append(
        "**Ending sentence for Part 1:** Tacstd2-high malignant / LUAD programs "
        "resolve first to a **TJ / apical-junction** enrichment in the human "
        "layers where rank supports it (concordant-4 ORA #1; GSE31210 KEGG TJ #3), "
        "setting up the next section on junction-linked immune geography — "
        "without pinning Part 1 on CLDN4 alone."
    )
    lines.append("")
    lines.append("## Best TJ/junction ranks (headline table)")
    lines.append("")
    lines.append(
        "| Layer | Best TJ/junction term | Rank | Universe | Stat | FDR | Lead? |"
    )
    lines.append("|---|---|---:|---|---|---:|:---:|")
    for h in headlines:
        lead = "YES" if h["lead_claim"] else "no"
        lines.append(
            f"| {h['layer']} | `{h['best_tj_term']}` | "
            f"**{h['best_tj_rank']}** / {h['n_universe_or_pos']} | "
            f"{h['universe']} | {h['stat']}={h['stat_value']:+.3f} | "
            f"{h['fdr']:.3g} | {lead} |"
        )
    lines.append("")
    lines.append("Full machine table: `tables/best_tj_junction_ranks.tsv`.")
    lines.append("")
    lines.append("## Concordant-4 detail (PR #741)")
    lines.append("")
    lines.append("### ORA (p<0.01 & |logFC|>0.25 UP; TACSTD2 held out)")
    lines.append("")
    lines.append("| Rank | Term | Enrichment | FDR |")
    lines.append("|---:|---|---:|---:|")
    lines.append(
        f"| **{c4_ora_aj['_rank']}** | HALLMARK_APICAL_JUNCTION | "
        f"{fnum(c4_ora_aj['enrichment']):.2f} | {fnum(c4_ora_aj['fdr']):.2e} |"
    )
    lines.append(
        f"| {c4_ora_krt['_rank']} | GOBP_KERATINIZATION | "
        f"{fnum(c4_ora_krt['enrichment']):.2f} | {fnum(c4_ora_krt['fdr']):.2e} |"
    )
    lines.append(
        f"| {c4_ora_tj['_rank']} | GOBP_TIGHT_JUNCTION_ORGANIZATION | "
        f"{fnum(c4_ora_tj['enrichment']):.2f} | {fnum(c4_ora_tj['fdr']):.3g} |"
    )
    lines.append(
        f"| {c4_ora_adh['_rank']} | CUSTOM_EPITHELIAL_ADHESION | "
        f"{fnum(c4_ora_adh['enrichment']):.2f} | {fnum(c4_ora_adh['fdr']):.3g} |"
    )
    lines.append("")
    lines.append(
        "**ORA verdict:** apical junction is the **#1** enrichment — strongest "
        "honest lead for Part 1."
    )
    lines.append("")
    lines.append("### Family scores (OLS, positive = TACSTD2 Q4)")
    lines.append("")
    lines.append("| Family | logFC | FDR |")
    lines.append("|---|---:|---:|")
    for name in ["TJ", "CLAUDIN_PANEL", "EPITHELIAL_ADHESION", "KERATIN", "IFN", "MHC-I/APM"]:
        # family column may use slightly different names
        row = next((r for r in fam if r.get("family") == name or r.get("Family") == name), None)
        if row is None:
            # try first column
            keys = list(fam[0].keys())
            row = next(r for r in fam if r[keys[0]] == name)
            name_k, lfc_k, fdr_k = keys[0], "logFC" if "logFC" in row else keys[2], "FDR" if "FDR" in row else "fdr"
            # normalize
            rk = {k.lower(): k for k in row}
            lines.append(
                f"| {row[keys[0]]} | {fnum(row[rk.get('logfc', keys[2])]):+.3f} | "
                f"{fnum(row[rk.get('fdr', keys[-1])]):.4g} |"
            )
        else:
            rk = {k.lower(): v for k, v in row.items()}
            lines.append(
                f"| {name} | {fnum(rk.get('logfc')):+.3f} | {fnum(rk.get('fdr')):.4g} |"
            )
    lines.append("")
    lines.append(
        "TJ / claudin / adhesion / keratin families are all up at FDR < 0.01. "
        "IFN and MHC-I are not down on this TACSTD2 split."
    )
    lines.append("")
    lines.append("### GSEA prerank (OLS *t*; TACSTD2 dropped; seed=42)")
    lines.append("")
    lines.append("| Rank (pos NES) | Term | NES | FDR |")
    lines.append("|---:|---|---:|---:|")
    lines.append(
        f"| {c4_gsea_krt['_rank_pos']} | GOBP_KERATINIZATION | "
        f"{fnum(c4_gsea_krt['nes']):+.3f} | {fnum(c4_gsea_krt['fdr']):.3g} |"
    )
    lines.append(
        f"| {c4_gsea_adh['_rank_pos']} | CUSTOM_EPITHELIAL_ADHESION | "
        f"{fnum(c4_gsea_adh['nes']):+.3f} | {fnum(c4_gsea_adh['fdr']):.3g} |"
    )
    lines.append(
        f"| {c4_gsea_aj['_rank_pos']} | HALLMARK_APICAL_JUNCTION | "
        f"{fnum(c4_gsea_aj['nes']):+.3f} | {fnum(c4_gsea_aj['fdr']):.3g} |"
    )
    lines.append(
        f"| {c4_gsea_tjorg['_rank_pos']} | GOBP_TIGHT_JUNCTION_ORGANIZATION | "
        f"{fnum(c4_gsea_tjorg['nes']):+.3f} | {fnum(c4_gsea_tjorg['fdr']):.3g} |"
    )
    lines.append(
        f"| {c4_gsea_kegg['_rank_pos']} | KEGG_TIGHT_JUNCTION | "
        f"{fnum(c4_gsea_kegg['nes']):+.3f} | {fnum(c4_gsea_kegg['fdr']):.3g} |"
    )
    lines.append("")
    lines.append(
        "GSEA puts keratinization near the top (rank 2) and apical junction at "
        f"rank {c4_gsea_aj['_rank_pos']}. Part 1 still leads with **TJ/junction** "
        "because ORA rank 1 is apical junction and the user OK'd TJ-first framing; "
        "keratin is reported as co-enriched barrier, not hidden."
    )
    lines.append("")
    lines.append("## Bulk LUAD detail (PR #736)")
    lines.append("")
    lines.append("### GSE31210 — TJ leads among positives")
    lines.append("")
    lines.append("| Rank (pos NES) | Term | NES | FDR_all |")
    lines.append("|---:|---|---:|---:|")
    for r in g31["pos_top10"][:8]:
        arrow = " ← TJ" if r["term"] == "KEGG_TIGHT_JUNCTION" else ""
        lines.append(
            f"| {r['_rank_pos']} | {r['term']}{arrow} | {fnum(r['nes']):+.3f} | "
            f"{fnum(r['fdr_bh_all']):.3g} |"
        )
    lines.append("")
    lines.append(
        f"**GSE31210 verdict:** KEGG Tight Junction is **rank {g_best['_rank_pos']}** "
        f"of {g_best['_n_pos']} positive-NES sets — a near-top / lead enrichment "
        "for Tacstd2-high bulk LUAD."
    )
    lines.append("")
    lines.append("### OncoSG — honest non-lead for KEGG TJ")
    lines.append("")
    lines.append("| Term | NES | FDR_all | Among pos? |")
    lines.append("|---|---:|---:|---|")
    lines.append(
        f"| GOBP_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY | "
        f"{fnum(o_best['nes']):+.3f} | {fnum(o_best['fdr_bh_all']):.3g} | "
        f"rank {o_best['_rank_pos']}/{o_best['_n_pos']} |"
    )
    lines.append(
        f"| GOBP_ESTABLISHMENT_OF_SKIN_BARRIER | "
        f"{fnum(onc['skin_barrier']['nes']):+.3f} | "
        f"{fnum(onc['skin_barrier']['fdr_bh_all']):.3g} | "
        f"rank {onc['skin_barrier']['_rank_pos']}/{onc['skin_barrier']['_n_pos']} |"
    )
    lines.append(
        f"| GOBP_KERATINIZATION | "
        f"{fnum(onc['keratinization']['nes']):+.3f} | "
        f"{fnum(onc['keratinization']['fdr_bh_all']):.3g} | "
        f"rank {onc['keratinization']['_rank_pos']}/{onc['keratinization']['_n_pos']} |"
    )
    lines.append(
        f"| KEGG_TIGHT_JUNCTION | {fnum(onc['kegg_tj']['nes']):+.3f} | "
        f"{fnum(onc['kegg_tj']['fdr_bh_all']):.3g} | not positive |"
    )
    lines.append(
        f"| HALLMARK_APICAL_JUNCTION | {fnum(onc['apical_junction']['nes']):+.3f} | "
        f"{fnum(onc['apical_junction']['fdr_bh_all']):.3g} | not positive |"
    )
    lines.append("")
    lines.append(
        "**OncoSG verdict:** do **not** claim KEGG TJ or apical junction as the "
        "lead enrichment here. Best strict TJ term = bicellular TJ assembly "
        f"(rank {o_best['_rank_pos']}). Barrier/keratin programs are up and "
        "rank higher — Part 1 states that instead of inventing a TJ #1."
    )
    lines.append("")
    lines.append("## What this Part 1 ending is / is not")
    lines.append("")
    lines.append("- **Is:** a rank-honest closer that Tacstd2-high → TJ/junction is the lead enrichment in human concordant-4 (ORA #1) and GSE31210 (KEGG TJ #3).")
    lines.append("- **Is:** TJ-first framing with keratin/skin-barrier co-enrichment reported.")
    lines.append("- **Is not:** a CLDN4 pin, a surface-molecule rank nail, or a private KD claim.")
    lines.append("- **Is not:** “TJ-high excludes T/NK” from the TACSTD2 funnel alone (PR #741 TJ vs T/NK was null).")
    lines.append("- **Is not:** a re-fit of locked CLDN4 %pos vs T/NK ρ = −0.531.")
    lines.append("- **Is not:** fabricated NES/FDR — every rank comes from PR #741 / #736 tables in `provenance/`.")
    lines.append("")
    lines.append("## Reproduce")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 methods/part1_tacstd2_tj_lead/rank_from_provenance.py")
    lines.append("```")
    lines.append("")
    lines.append("## Files")
    lines.append("")
    lines.append("- `tables/best_tj_junction_ranks.tsv` — headline ranks")
    lines.append("- `tables/gsea_tj_junction_ranks.tsv` — per-cohort TJ/barrier GSEA ranks")
    lines.append("- `tables/c4_ora_tj_barrier_ranks.tsv` — C4 ORA barrier subset")
    lines.append("- `tables/summary.json` — machine verdict")
    lines.append("- `figures/fig_best_tj_junction_ranks.png` — Part 1 rank figure")
    lines.append("- `provenance/` — frozen copies of PR #741 / #736 tables")
    lines.append("")

    (HERE / "FINDING.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
