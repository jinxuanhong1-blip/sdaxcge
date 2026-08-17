#!/usr/bin/env python3
"""Write FINDING.md from analyze.py outputs (honest n + ligand table)."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def fmt_n(x) -> str:
    return f"{int(x):,}"


def n_row(n_cells: list[dict], group: str) -> dict | None:
    for r in n_cells:
        if r["group"] == group:
            return r
    return None


def pair_rows(df: pd.DataFrame, direction: str, high: bool, k: int = 10) -> list[str]:
    sub = df[df["direction"] == direction].copy()
    if sub.empty:
        return ["| *(none)* | | | | | |"]
    if high:
        sub = sub[sub["delta_prob"] > 0].sort_values("delta_prob", ascending=False)
    else:
        sub = sub[sub["delta_prob"] < 0].sort_values("delta_prob")
    if sub.empty:
        return ["| *(none)* | | | | | |"]
    lines = []
    for rec in sub.head(k).itertuples(index=False):
        pmin = min(float(rec.pval_high), float(rec.pval_low))
        lines.append(
            f"| {rec.interaction_name} | {rec.pathway_name} | {rec.ligand_class} | "
            f"{rec.delta_prob:+.3f} | {rec.prob_high:.3f} | {rec.prob_low:.3f} | {pmin:.3f} |"
        )
    return lines


def main() -> None:
    out = ROOT / "results"
    header = json.loads((out / "summary.json").read_text())
    kept = header.get("kept") or {}
    per_sample = pd.read_csv(out / "per_sample_post.tsv", sep="\t")
    ligand_path = out / "ligand_table.tsv"
    ligand = pd.read_csv(ligand_path, sep="\t") if ligand_path.exists() else pd.DataFrame()

    epi_total = int(per_sample["n_epithelial"].sum())
    tnk_total = int(per_sample["n_TNK"].sum())
    n_pts = int(per_sample["sample"].nunique())
    n_mpr = int((per_sample["response"] == "MPR").sum())
    n_nmpr = int((per_sample["response"] == "NMPR").sum())
    top = per_sample.sort_values("n_epithelial", ascending=False).iloc[0]
    top_frac = 100.0 * float(top.n_epithelial) / epi_total if epi_total else 0.0
    nmpr_epi = int(per_sample.loc[per_sample.response == "NMPR", "n_epithelial"].sum())
    top_nmpr_frac = 100.0 * float(top.n_epithelial) / nmpr_epi if (top.response == "NMPR" and nmpr_epi) else None

    run_lines = [
        "| Split | Mal high | Mal low | TNK | Detected LR | Sig LR (out / in) | Cold/barrier | Recruit-up in high |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for s in header.get("runs", []):
        cells = {r["group"]: r for r in s["n_cells"]}
        if s["mode"] == "combo_mpr":
            hi = f"{cells.get('Mal_high_NMPR', {}).get('n_cells', '?')}+{cells.get('Mal_high_MPR', {}).get('n_cells', '?')}"
            lo = f"{cells.get('Mal_low_NMPR', {}).get('n_cells', '?')}+{cells.get('Mal_low_MPR', {}).get('n_cells', '?')}"
            tnk = f"{cells.get('TNK_NMPR', {}).get('n_cells', '?')}+{cells.get('TNK_MPR', {}).get('n_cells', '?')}"
        else:
            hi = cells.get("Mal_high", {}).get("n_cells", "?")
            lo = cells.get("Mal_low", {}).get("n_cells", "?")
            tnk = cells.get("TNK", {}).get("n_cells", "?")
        cb = s["cold_barrier"]
        mark = "**" if s["mode"] == header.get("kept_split") else ""
        run_lines.append(
            f"| {mark}{s['mode']}{mark} | {mark}{hi}{mark} | {mark}{lo}{mark} | {mark}{tnk}{mark} | "
            f"{mark}{s['n_detected']}{mark} | {mark}{s['n_significant']} "
            f"({s['n_sig_outgoing_mal_to_tnk']} / {s['n_sig_incoming_tnk_to_mal']}){mark} | "
            f"{mark}{cb['score']}{mark} | {mark}{cb['n_recruit_up_in_high']}{mark} |"
        )

    kept_cells = {r["group"]: r for r in kept.get("n_cells", [])}
    mal_hi = kept_cells.get("Mal_high", {})
    mal_lo = kept_cells.get("Mal_low", {})
    tnk = kept_cells.get("TNK", {})

    n_sig = int(kept.get("n_significant", 0))
    n_det = int(kept.get("n_detected", 0))
    n_out = int(kept.get("n_sig_outgoing_mal_to_tnk", 0))
    n_in = int(kept.get("n_sig_incoming_tnk_to_mal", 0))
    n_tested = int(kept.get("n_lr_tested", 0))

    if not ligand.empty:
        n_lig_out = int((ligand.direction == "outgoing").sum())
        n_lig_in = int((ligand.direction == "incoming").sum())
        n_high_out = int(((ligand.direction == "outgoing") & (ligand.delta_prob > 0)).sum())
        n_low_out = int(((ligand.direction == "outgoing") & (ligand.delta_prob < 0)).sum())
    else:
        n_lig_out = n_lig_in = n_high_out = n_low_out = 0

    sample_lines = [
        "| Sample | Response | Epithelial | T+NK | Mean CLDN4 (epi) |",
        "| --- | --- | ---: | ---: | ---: |",
    ]
    for rec in per_sample.sort_values(["response", "sample"]).itertuples(index=False):
        mean_c = rec.mean_CLDN4_epithelial
        mean_s = "NA" if pd.isna(mean_c) else f"{float(mean_c):.2f}"
        sample_lines.append(
            f"| {rec.sample} | {rec.response} | {int(rec.n_epithelial):,} | {int(rec.n_TNK):,} | {mean_s} |"
        )

    kept_mode = header.get("kept_split")
    cb = kept.get("cold_barrier", {})
    reason = header.get("kept_reason", "")

    md = f"""# FINDING — CellChat-style CLDN4-high vs low vs T/NK (GSE207422)

Additive to prior TACSTD2 CellChat (`methods/scrna_cellchat`). That folder is taken as given and was **not** re-run. **CLDN4 only.**

Public UMI only (Hu et al., *Genome Medicine* 2023, PMID 36869384). CellChat R and LIANA were not run. Pairs below use Hill probability + 100 high/low permutations in `scripts/analyze.py`. A pair is significant if detected (`expr_prop ≥ 0.10` on both sides), `P > 0`, and permutation `p < 0.05`. Smallest possible p with 100 permutations is 1/101 = 0.0099.

No edge is drawn unless it meets that rule. The ligand table is `results/ligand_table.tsv` (also `ligand_table_kept.tsv`). Extra figure: `results/fig_extra_ligand_table.png`. Honest per-sample n: `results/per_sample_post.tsv` and `fig_n_per_sample.png`.

## English

### Data and n

| Item | n | Note |
| --- | ---: | --- |
| Genes × barcodes in the GEO UMI | {fmt_n(header['matrix_genes'])} × {fmt_n(header['matrix_cells'])} | one matrix |
| CellChatDB v2 protein pairs with every subunit in the matrix | {fmt_n(header['lr_pairs_in_matrix'])} | of 2,239 protein pairs |
| Post-treatment samples used for splits | {n_pts} | MPR+pCR = {n_mpr}; NMPR = {n_nmpr} |
| Post-treatment epithelial cells | {fmt_n(epi_total)} | malignant proxy (CopyKAT IDs not public) |
| Post-treatment T+NK | {fmt_n(tnk_total)} | T and NK merged as TNK |
| Pre-treatment biopsies | excluded | 3 samples |

**{top['sample']}** ({top['response']}) is **{fmt_n(top.n_epithelial)} / {fmt_n(epi_total)} ({top_frac:.1f}%)** of post-treatment epithelial cells{f" and **{fmt_n(top.n_epithelial)} / {fmt_n(nmpr_epi)} ({top_nmpr_frac:.1f}%)** of NMPR epithelial cells" if top_nmpr_frac is not None else ""}. Group means are **cell-pooled**, so this sample can dominate NMPR and pooled contrasts. Unit of inference for a patient claim would be n={n_mpr} vs {n_nmpr}; that test is **not** what the permutation does.

### Per-sample n (post-treatment)

{chr(10).join(sample_lines)}

### Splits compared (CLDN4 `log1p(CP10k)` on epithelium)

Tertile high/low drops the middle third.

{chr(10).join(run_lines)}

Combo score adds NMPR and MPR arms and is not comparable to a single contrast.

### Kept split

**{kept_mode}.** {reason}

Kept n: CLDN4-high epithelial **{fmt_n(mal_hi.get('n_cells', 0))}** cells / **{mal_hi.get('n_samples', 0)}** patients (mean CLDN4 {mal_hi.get('mean_CLDN4_log1p_cp10k', float('nan')):.3f}); CLDN4-low **{fmt_n(mal_lo.get('n_cells', 0))}** / **{mal_lo.get('n_samples', 0)}** (mean {mal_lo.get('mean_CLDN4_log1p_cp10k', float('nan')):.3f}); T/NK **{fmt_n(tnk.get('n_cells', 0))}** / **{tnk.get('n_samples', 0)}**.

Cold/barrier score = {cb.get('score')} (barrier-up {cb.get('n_barrier_up_in_high')}, inhib-up {cb.get('n_inhib_up_in_high')}, recruit-up {cb.get('n_recruit_up_in_high')}, recruit-down {cb.get('n_recruit_down_in_high')}, attack-in-down {cb.get('n_attack_in_down_in_high')}, barrier-down {cb.get('n_barrier_down_in_high')}).

### Ligand table (kept split, significant ΔP only)

Among {fmt_n(n_tested)} directed tests: detected {fmt_n(n_det)}; **significant {fmt_n(n_sig)}** (outgoing Mal→TNK {n_out}; incoming TNK→Mal {n_in}). Ligand table rows: **{len(ligand)}** differential pairs ({n_lig_out} outgoing, {n_lig_in} incoming; outgoing higher in CLDN4-high = {n_high_out}, higher in low = {n_low_out}).

Full table: [`results/ligand_table.tsv`](results/ligand_table.tsv). Extra figure: [`results/fig_extra_ligand_table.png`](results/fig_extra_ligand_table.png).

#### Outgoing Mal → T/NK, higher in CLDN4-high

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
{chr(10).join(pair_rows(ligand, "outgoing", True, 12))}

#### Outgoing Mal → T/NK, higher in CLDN4-low

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
{chr(10).join(pair_rows(ligand, "outgoing", False, 12))}

#### Incoming T/NK → Mal, higher into CLDN4-high

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
{chr(10).join(pair_rows(ligand, "incoming", True, 8))}

#### Incoming T/NK → Mal, higher into CLDN4-low

| Pair | Pathway | Class | ΔP (high−low) | P high | P low | p min |
| --- | --- | --- | ---: | ---: | ---: | ---: |
{chr(10).join(pair_rows(ligand, "incoming", False, 8))}

### What is not claimed

- These are permutation tests on cell-pooled truncated means, not a patient-level mixed model. Patient n = {n_mpr} MPR vs {n_nmpr} NMPR.
- {top['sample']} supplies a large share of post epithelial UMIs.
- CellChat R visualizations were not generated. Figures show only pairs that pass p < 0.05.
- Author cell-type and CopyKAT objects were not used.
- TACSTD2 was not used to define high/low. This is not a dual-high or TACSTD2 re-analysis.
- GSE253013 / GSE241934 were not used.

### Files

`results/ligand_table.tsv`, `ligand_table_kept.tsv`, `key_pairs_kept.tsv`, `per_sample_post.tsv`, `summary.json`, `fig_extra_ligand_table.png`, `fig_n_per_sample.png`, `fig_n_*.png`, `fig_nsig_*.png`, `fig_top_*.png`.

---

## 中文

**这是加性分析，只切 CLDN4。** 先前的 TACSTD2 CellChat 当作已给，不重跑。公开 GSE207422 一张 UMI（{fmt_n(header['matrix_genes'])} 基因 × {fmt_n(header['matrix_cells'])} 细胞）。术后 **MPR {n_mpr}（含 pCR）/ NMPR {n_nmpr}**。上皮（恶性代理）{fmt_n(epi_total)}，T+NK {fmt_n(tnk_total)}。**{top['sample']}** 占术后上皮 {top_frac:.1f}%。均值是细胞池化，不是患者混合模型。

保留切分 **{kept_mode}**：CLDN4-high {fmt_n(mal_hi.get('n_cells', 0))} / {mal_hi.get('n_samples', 0)} 例，low {fmt_n(mal_lo.get('n_cells', 0))} / {mal_lo.get('n_samples', 0)} 例，T/NK {fmt_n(tnk.get('n_cells', 0))}。显著差异配体–受体 **{len(ligand)}** 行，见 `results/ligand_table.tsv`。额外图 `fig_extra_ligand_table.png`。
"""
    (ROOT / "FINDING.md").write_text(md)
    print(f"wrote {ROOT / 'FINDING.md'} kept={kept_mode} ligand_rows={len(ligand)}")


if __name__ == "__main__":
    main()
