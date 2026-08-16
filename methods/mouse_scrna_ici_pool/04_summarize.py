#!/usr/bin/env python3
"""Write combined direction table, figures, and bilingual WRITEUP."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results/mouse_scrna_ici_pool"
FIG = OUT / "figures"


def fmt_p(p):
    if p is None or (isinstance(p, float) and (np.isnan(p) or p != p)):
        return "NA (n<2)"
    p = float(p)
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.3g}"


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    samp = pd.read_csv(OUT / "sample_inventory.tsv", sep="\t")
    cells = pd.read_csv(OUT / "compartment_scores.tsv", sep="\t")
    contr = pd.read_csv(OUT / "contrasts.tsv", sep="\t")
    both = pd.read_csv(OUT / "epi_vs_tnk_per_sample.tsv", sep="\t")

    # direction table: ICB-relevant epithelial contrasts + epi vs T/NK
    keep = contr[
        contr.contrast.str.contains(
            "aPD1_vs|CA170_vs|SCLC_aPD1|CD45neg_vs|epithelial_vs_tnk|RTICI_vs"
        )
    ].copy()
    keep.to_csv(OUT / "direction_table.tsv", sep="\t", index=False)

    # figure: epithelial Tacstd2 ICB-ish deltas
    epi_icb = contr[(contr.compartment == "epithelial") & (contr.gene == "Tacstd2")]
    if not epi_icb.empty:
        fig, ax = plt.subplots(figsize=(8, 3.6))
        y = np.arange(len(epi_icb))
        colors = ["#b2182b" if d > 0 else "#2166ac" for d in epi_icb.delta_log]
        ax.barh(y, epi_icb.delta_log, color=colors)
        ax.set_yticks(y)
        ax.set_yticklabels([f"{r.gse} {r.contrast}" for r in epi_icb.itertuples()], fontsize=8)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_xlabel("Tacstd2 Δ mean log (treated − reference)")
        ax.set_title("Epithelial/tumor Tacstd2 (n=1 vs 1 unless noted)")
        fig.tight_layout()
        fig.savefig(FIG / "epithelial_tacstd2_deltas.png", dpi=140)
        plt.close()

    if not both.empty:
        fig, ax = plt.subplots(figsize=(5.2, 4.2))
        for gene, g in both.groupby("gene"):
            ax.scatter(g.tnk_mean, g.epi_mean, label=gene, s=40)
        lim = max(both.epi_mean.max(), both.tnk_mean.max()) * 1.1
        ax.plot([0, lim], [0, lim], "k--", lw=0.8)
        ax.set_xlabel("T/NK mean log")
        ax.set_ylabel("Epithelial/tumor mean log")
        ax.legend()
        ax.set_title("Paired samples with both compartments (≥20 cells)")
        fig.tight_layout()
        fig.savefig(FIG / "epi_vs_tnk_scatter.png", dpi=140)
        plt.close()

    # WRITEUP
    lines = []
    lines.append("# Mouse lung ICI scRNA pool — Tacstd2 / Cldn4")
    lines.append("")
    lines.append("**Additive.** User A4 TISMO (49/64 Tacstd2 up after ICB, Wilcoxon p=5.8e-5) is **taken as given** and is not re-audited. This slice is public **mouse lung ICI scRNA**, not TISMO bulk.")
    lines.append("")
    lines.append("## Verdict")
    lines.append("")
    lines.append("No public mouse-lung ICI scRNA series in this pool has **n≥2 biological libraries per arm** of epithelial/tumor Tacstd2 or Cldn4 after ICB vs a true no-ICB control. Several series **lack epithelial/tumor cells** (CD45+/CD3+ sorts). GSE283827 (SCLC ± aPD-1) cannot be scored: MTX has 32,589 genes and **no features.tsv**.")
    lines.append("")
    lines.append("Paired same-library epithelial vs T/NK: Tacstd2 is higher in epithelial in **10/13** samples (Wilcoxon p=0.068 — not a firm restriction claim). Cldn4 is higher in epithelial in **11/13** (p=9.8e-4). GSE157881 deposited CD45− Tacstd2 is **lower** than CD45+ T/NK. This is not an ICB-response claim and does **not** reproduce TISMO 49/64.")
    lines.append("")
    lines.append("## Catalog (scored vs documented)")
    lines.append("")
    lines.append("| Series | Model | ICB in *these* libraries | Epithelial/tumor | T/NK | Tacstd2 / Cldn4 |")
    lines.append("|---|---|---|---|---|---|")
    lines.append("| GSE157881 | HKP1 lung, 0 vs 4 Gy | **No** (RT only) | CD45− present | CD45+ present | scored |")
    lines.append("| GSE157882 | HKP1, club-cell DT after 4 Gy | **No** | **ABSENT** (CD45+ only) | present | scored in T/NK |")
    lines.append("| GSE176091 | VC-LUAD, CA170 vs PBS | CA170 (VISTA/PD-L1) | **ABSENT** (CD45+ TIL) | present; n=2 vs 2 | scored in T/NK |")
    lines.append("| GSE267557 | MC38-bearing lung, young vs aged | both aPD-1 | **ABSENT** (CD45+) | present | not downloaded (1.2 GB); no control arm |")
    lines.append("| GSE268525 | LLC1-sgLkb1 lung met | ICI vs RT+ICI (**no untreated**) | scored if Epcam/Cdh1+Krt8 | scored | n=1 vs 1 |")
    lines.append("| GSE283827 | Rb/p53 SCLC ± aPD-1 ± ERBB2i | yes (Ly16 aPD-1, Ly29 vehicle) | authors: epithelial | authors: immune | **not scored — no gene names** |")
    lines.append("| GSE303943 | CMT167R subcutaneous | **No** (PKCi vs solvent) | scored | scored | not ICB; subcutaneous |")
    lines.append("| GSE133604 leftover | KP ± aPD-1 | yes | scored | scored | n=1 vs 1; features = mm10-3.0 31053 |")
    lines.append("| GSE129297 leftover | SCLC GEMM ± aPD-1 | yes | scored | scored | n=1 vs 1; unfiltered 10x, UMI≥500 |")
    lines.append("| GSE232730 leftover | KPL-3M CD45+ | yes | **ABSENT** | present | n=1 vs 1 |")
    lines.append("| GSE222158 leftover | FVB lung CD45+ | yes | **ABSENT** | present | n=1 vs 1 |")
    lines.append("| GSE275877 leftover | LKR13 ICI R vs acquired-NR | yes | unknown | unknown | **skip**: 31,053 × 33.97M unfiltered barcodes |")
    lines.append("")
    lines.append("## Sample-level epithelial / T/NK Tacstd2")
    lines.append("")
    epi = cells[cells.compartment == "epithelial"][["gse", "sample", "arm", "n_cells", "Tacstd2_mean", "Tacstd2_pct", "Cldn4_mean", "Cldn4_pct"]]
    tnk = cells[cells.compartment == "tnk"][["gse", "sample", "n_cells", "Tacstd2_mean", "Tacstd2_pct"]]
    lines.append("### Epithelial / tumor")
    lines.append("")
    lines.append("| Series | Sample | Arm | n cells | Tacstd2 mean | Tacstd2 %pos | Cldn4 mean | Cldn4 %pos |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for r in epi.itertuples():
        if r.n_cells < 1:
            continue
        lines.append(
            f"| {r.gse} | {r.sample} | {r.arm} | {int(r.n_cells)} | {r.Tacstd2_mean:.3f} | {r.Tacstd2_pct:.1f} | {r.Cldn4_mean:.3f} | {r.Cldn4_pct:.1f} |"
        )
    lines.append("")
    lines.append("### T/NK")
    lines.append("")
    lines.append("| Series | Sample | n cells | Tacstd2 mean | Tacstd2 %pos |")
    lines.append("|---|---|---|---|---|")
    for r in tnk.itertuples():
        if r.n_cells < 1:
            continue
        lines.append(f"| {r.gse} | {r.sample} | {int(r.n_cells)} | {r.Tacstd2_mean:.3f} | {r.Tacstd2_pct:.1f} |")
    lines.append("")
    lines.append("## Contrasts (honest n / p)")
    lines.append("")
    lines.append("| Series | Contrast | Compartment | Gene | n_a | n_b | Δ | Welch p | MWU p | Note |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|")
    for r in contr.itertuples():
        da = f"{r.delta_log:.3f}" if pd.notna(r.delta_log) else "NA"
        note = "" if pd.isna(r.note) else str(r.note)
        lines.append(
            f"| {r.gse} | {r.contrast} | {r.compartment} | {r.gene} | {int(r.n_a)} | {int(r.n_b)} | {da} | {fmt_p(r.welch_p)} | {fmt_p(r.mwu_p)} | {note} |"
        )
    lines.append("")
    lines.append("## Combined direction (epithelial ICB-ish + epi vs T/NK)")
    lines.append("")
    lines.append("Only **GSE176091** has n=2 vs 2, and that is **T/NK in a CD45+ sort** (epithelial absent). All epithelial ICB-vs-control tests are **n=1 vs 1** (leftover GSE133604, GSE129297). Sign is reported; p is not.")
    lines.append("")
    lines.append("| Series | Contrast | Compartment | Gene | Direction (treated − ref) | n | p |")
    lines.append("|---|---|---|---|---|---|---|")
    lines.append("| GSE133604 leftover | KP aPD-1 vs IgG | epithelial | Tacstd2 | down (tiny) | 1 vs 1 | NA |")
    lines.append("| GSE133604 leftover | KP aPD-1 vs IgG | epithelial | Cldn4 | down | 1 vs 1 | NA |")
    lines.append("| GSE129297 leftover | SCLC aPD-1 vs Ctrl | epithelial | Tacstd2 | up (tiny) | 1 vs 1 | NA |")
    lines.append("| GSE129297 leftover | SCLC aPD-1 vs Ctrl | epithelial | Cldn4 | down | 1 vs 1 | NA |")
    lines.append("| GSE176091 | CA170 vs PBS | T/NK (CD45+; epi **ABSENT**) | Tacstd2 | down | 2 vs 2 | Welch 0.17; MWU 0.33 |")
    lines.append("| GSE176091 | CA170 vs PBS | T/NK | Cldn4 | tie | 2 vs 2 | Welch 0.98; MWU 1 |")
    lines.append("| GSE268525 | RT+ICI vs ICI (no untreated) | epithelial | Tacstd2 | tie/down | 1 vs 1 | NA |")
    lines.append("| GSE268525 | RT+ICI vs ICI | epithelial | Cldn4 | up | 1 vs 1 | NA |")
    lines.append("| GSE157881 | CD45− vs CD45+ (RT, not ICB) | epi − T/NK | Tacstd2 | **down** (epi < T/NK) | 1 vs 1 | NA |")
    lines.append("| POOL | same-library epi vs T/NK | epithelial − T/NK | Tacstd2 | 10/13 epi > T/NK | 13 paired | Wilcoxon 0.068 |")
    lines.append("| POOL | same-library epi vs T/NK | epithelial − T/NK | Cldn4 | 11/13 epi > T/NK | 13 paired | Wilcoxon 9.8e-4 |")
    lines.append("")
    if not both.empty:
        n = both[both.gene == "Tacstd2"].shape[0]
        n_up = int((both[both.gene == "Tacstd2"].delta_epi_minus_tnk > 0).sum()) if n else 0
        n_c = both[both.gene == "Cldn4"].shape[0]
        n_c_up = int((both[both.gene == "Cldn4"].delta_epi_minus_tnk > 0).sum()) if n_c else 0
        lines.append(f"Paired epithelial vs T/NK Tacstd2: **{n_up}/{n}** samples epithelial > T/NK (samples with ≥20 cells in both); Wilcoxon signed-rank p=0.068. Cldn4: **{n_c_up}/{n_c}**, p=9.8e-4.")
    lines.append("")
    lines.append("## Methods (short)")
    lines.append("")
    lines.append("- Public GEO processed MTX/TSV only. Tacstd2 = ENSMUSG00000051397; Cldn4 = ENSMUSG00000047501.")
    lines.append("- Unsorted 10x: epithelial/tumor = Epcam+ or (Cdh1+ and Krt8+) or Ascl1/Chga/Insm1+; T/NK = Cd3d/e, Cd8a, Nkg7, Ncr1 and not epithelial.")
    lines.append("- GSE157881/882 deposited log-like values (used as-is). 10x MTX: log1p(UMI).")
    lines.append("- GSE129297/GSE133604 gene names: Cell Ranger mm10-3.0 31,053-gene table from GSE275877 features (same n_genes).")
    lines.append("- Welch / MWU only if n≥2 samples/arm. No invented R vs NR labels.")
    lines.append("")
    lines.append("## 中文摘要")
    lines.append("")
    lines.append("TISMO 用户 A4（49/64 Tacstd2 在 ICB 后升高，p=5.8e-5）**视为已知，不重算**。本切片是公开小鼠肺 ICI **单细胞**，不是 TISMO bulk。")
    lines.append("候选集中：GSE157881/882 是放疗/club 细胞清除而非 ICB；GSE176091 / GSE267557 / GSE232730 / GSE222158 为 CD45+ 或 CD3+，**无上皮/肿瘤细胞**；GSE283827 有 aPD-1 但 **未提供 features.tsv**，无法对 Tacstd2/Cldn4 诚实计分；GSE268525 无未治疗对照；GSE303943 为皮下 PKCi。")
    lines.append("有上皮且有 ICB vs 对照的只剩 leftover GSE133604（KP）与 GSE129297（SCLC），均为 **n=1 vs 1**，只报方向，不算 p。配对上皮 vs T/NK：Tacstd2 10/13 上皮更高（Wilcoxon p=0.068）；Cldn4 11/13（p=9.8e-4）。GSE157881 CD45− Tacstd2 反而低于 CD45+ T/NK。不能外推 TISMO 49/64。")
    lines.append("")

    (OUT / "WRITEUP.md").write_text("\n".join(lines), encoding="utf-8")
    (OUT / "summary_counts.json").write_text(json.dumps({
        "n_samples_scored": int((samp.status == "ok").sum()) if "status" in samp else int(len(samp)),
        "n_contrasts": int(len(contr)),
        "n_paired_epi_tnk": int(len(both)),
    }, indent=2), encoding="utf-8")
    print("wrote", OUT / "WRITEUP.md")


if __name__ == "__main__":
    main()
