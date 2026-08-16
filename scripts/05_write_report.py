"""Assemble results/claim_A10/REPORT.md from the tables. No peeking while computing."""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402


def load(name):
    p = os.path.join(C.TABLES, name)
    return pd.read_csv(p) if os.path.exists(p) else None


def fmt_rho(r):
    if r is None or (isinstance(r, float) and not np.isfinite(r)):
        return "NA"
    return f"{r:.3f}"


def pair_row(df, ds, a, b):
    if df is None:
        return None
    r = df[(df.dataset == ds) & (((df.gene_a == a) & (df.gene_b == b)) |
                                 ((df.gene_a == b) & (df.gene_b == a)))]
    return r.iloc[0] if len(r) else None


def verdict_from_rules(direction_ok, fdr, abs_rho, bg_pct, n_datasets_ok, n_datasets):
    """Apply the pre-registered decision rules as closely as a script can.

    Human review still owns the overall wording; this only proposes a label.
    """
    if n_datasets_ok == 0:
        return "NOT SUPPORTED"
    if direction_ok and fdr < 0.05 and abs_rho >= 0.3 and (
            bg_pct >= 90 or bg_pct <= 10) and n_datasets_ok >= 2:
        return "SUPPORTED"
    if direction_ok and n_datasets_ok >= 1:
        return "PARTIALLY SUPPORTED"
    return "NOT SUPPORTED"


def main():
    det = load("human_bulk_gene_detection.csv")
    pairs = load("human_bulk_pairwise_spearman.csv")
    coh = load("human_bulk_module_coherence.csv")
    part = load("human_bulk_partial_correlation.csv")
    mod = load("human_bulk_module_score_assoc.csv")
    strata = load("human_tcga_nkx2-1_low_vs_high.csv")
    tn = load("human_tcga_tumour_vs_normal.csv")
    pert = load("perturbation_contrasts.csv")
    hmeans = load("census_human_celltype_means.csv")
    hwithin = load("census_human_within_type_spearman.csv")
    mmeans = load("census_mouse_celltype_means.csv")
    mwithin = load("census_mouse_within_type_spearman.csv")

    lines = []
    A = lines.append
    A("# Claim A10 — results (honest)")
    A("")
    A("This file is generated from `results/claim_A10/tables/` by `scripts/05_write_report.py`.")
    A("The claim, predictions, and decision rules were pre-registered in")
    A("`claims/claim_A10.md` before any of these numbers were inspected.")
    A("")
    A("**Claim (operationalised):** ELF3, GRHL1, KLF4 and TFAP2A form a coherent")
    A("transcription-factor module that is positively coupled to TACSTD2 and CLDN4,")
    A("in a lung-epithelial state where NKX2-1 is down. Expected in human and mouse.")
    A("")
    A("## Overall verdict")
    A("")
    A("_Filled after the per-prediction sections below._")
    A("")

    # ---------- detection ----------
    A("## Data actually used")
    A("")
    A("| dataset | role | n | notes |")
    A("|---|---|---|---|")
    if det is not None:
        for ds in det.dataset.unique():
            n = "see pairwise table"
            if pairs is not None and (pairs.dataset == ds).any():
                n = int(pairs.loc[pairs.dataset == ds, "n"].iloc[0])
            A(f"| {ds} (recount3) | human bulk observational | {n} | log2CPM of coverage sums |")
    if hmeans is not None:
        A(f"| CELLxGENE Census human lung epithelium | P6 / cell type | "
          f"{int(hmeans.n_cells.sum())} cells | primary data, normal + LUAD |")
    if mmeans is not None:
        A(f"| CELLxGENE Census mouse lung | P5/P6 | {int(mmeans.n_cells.sum())} cells | "
          f"all primary lung cells in Census |")
    if pert is not None and len(pert):
        for ds in pert.dataset.unique():
            A(f"| {ds} | perturbation / model | — | see `perturbation_contrasts.csv` |")
    A("")
    A("### Claim-gene detection (human bulk)")
    A("")
    if det is not None:
        A("| dataset | gene | mean log2CPM | fraction samples > 1 |")
        A("|---|---|---:|---:|")
        for _, r in det[det.gene.isin(C.CLAIM_GENES_H)].iterrows():
            A(f"| {r.dataset} | {r.gene} | {r.mean_log2cpm:.2f} | {r.frac_samples_gt1_log2cpm:.2f} |")
    A("")

    # ---------- P1 ----------
    A("## P1 — the four TFs are a module")
    A("")
    A("Pre-registered: mean pairwise Spearman among ELF3/GRHL1/KLF4/TFAP2A exceeds")
    A("expression-matched random 4-gene sets, and |ρ| ≥ 0.3.")
    A("")
    if coh is not None:
        A("| dataset | mean pairwise ρ | matched-null mean | matched-null p95 | empirical p | random-pair percentile |")
        A("|---|---:|---:|---:|---:|---:|")
        for _, r in coh.iterrows():
            A(f"| {r.dataset} | {r.mean_pairwise_rho:.3f} | {r.matched_null_mean:.3f} | "
              f"{r.matched_null_p95:.3f} | {r.empirical_p_vs_matched_null:.4g} | "
              f"{r.percentile_vs_random_pairs:.1f} |")
    A("")
    A("Pairwise TF–TF (human bulk):")
    A("")
    A("| dataset | pair | ρ | 95% CI | FDR | bg percentile |")
    A("|---|---|---:|---|---:|---:|")
    if pairs is not None:
        for ds in ["GTEx-lung", "TCGA-LUAD"]:
            for i, a in enumerate(C.MODULE_H):
                for b in C.MODULE_H[i + 1:]:
                    r = pair_row(pairs, ds, a, b)
                    if r is None:
                        continue
                    A(f"| {ds} | {a}–{b} | {r.spearman_rho:.3f} | "
                      f"[{r.ci_low:.3f}, {r.ci_high:.3f}] | {r.fdr:.2e} | {r.bg_percentile:.1f} |")
    A("")

    # ---------- P2 ----------
    A("## P2 — module / each TF correlates positively with TACSTD2 and CLDN4")
    A("")
    A("| dataset | pair | ρ | 95% CI | FDR | bg percentile |")
    A("|---|---|---:|---|---:|---:|")
    if pairs is not None:
        for ds in ["GTEx-lung", "TCGA-LUAD"]:
            for a in C.MODULE_H:
                for b in C.TARGET_H:
                    r = pair_row(pairs, ds, a, b)
                    if r is None:
                        continue
                    A(f"| {ds} | {a}–{b} | {r.spearman_rho:.3f} | "
                      f"[{r.ci_low:.3f}, {r.ci_high:.3f}] | {r.fdr:.2e} | {r.bg_percentile:.1f} |")
    A("")
    if mod is not None:
        A("Module score vs targets:")
        A("")
        A("| dataset | y | ρ | 95% CI | FDR |")
        A("|---|---|---:|---|---:|")
        for _, r in mod[mod.y.isin(C.TARGET_H)].iterrows():
            A(f"| {r.dataset} | {r.y} | {r.spearman_rho:.3f} | "
              f"[{r.ci_low:.3f}, {r.ci_high:.3f}] | {r.fdr:.2e} |")
        A("")

    # ---------- P3 ----------
    A("## P3 — NKX2-1 is anti-correlated with the module and the targets")
    A("")
    A("| dataset | pair | ρ | 95% CI | FDR | bg percentile |")
    A("|---|---|---:|---|---:|---:|")
    if pairs is not None:
        for ds in ["GTEx-lung", "TCGA-LUAD"]:
            for a in C.MODULE_H + C.TARGET_H:
                r = pair_row(pairs, ds, a, C.NKX_H)
                if r is None:
                    continue
                A(f"| {ds} | {a}–NKX2-1 | {r.spearman_rho:.3f} | "
                  f"[{r.ci_low:.3f}, {r.ci_high:.3f}] | {r.fdr:.2e} | {r.bg_percentile:.1f} |")
    A("")
    if strata is not None and len(strata):
        A("TCGA-LUAD NKX2-1-low vs NKX2-1-high tumours (quartile split):")
        A("")
        A("| gene | n_low | n_high | log2FC (low−high) | Cliff's δ | FDR |")
        A("|---|---:|---:|---:|---:|---:|")
        for _, r in strata[strata.gene.isin(C.CLAIM_GENES_H + ["EPCAM", "SFTPC", "SCGB1A1", "KRT5"])].iterrows():
            A(f"| {r.gene} | {int(r.n_low)} | {int(r.n_high)} | {r['log2fc_low_minus_high']:.3f} | "
              f"{r.cliffs_delta:.3f} | {r.fdr:.2e} |")
        A("")

    # ---------- composition ----------
    A("## Composition check (bulk stand-in for P6)")
    A("")
    A("If the bulk associations are only 'this biopsy has more airway epithelium',")
    A("they should collapse after residualising on EPCAM / SFTPC / SCGB1A1 / lineage markers.")
    A("")
    if part is not None:
        A("| dataset | covariates | pair | partial ρ | FDR |")
        A("|---|---|---|---:|---:|")
        focus_pairs = [(a, b) for a in C.MODULE_H for b in C.TARGET_H]
        focus_pairs += [(a, C.NKX_H) for a in ["ELF3", "TACSTD2", "CLDN4"]]
        for ds in ["GTEx-lung", "TCGA-LUAD"]:
            for cov in ["none", "EPCAM", "EPCAM+SFTPC+SCGB1A1", "all_lineage_markers"]:
                for a, b in focus_pairs:
                    r = part[(part.dataset == ds) & (part.covariates == cov) &
                             (((part.gene_a == a) & (part.gene_b == b)) |
                              ((part.gene_a == b) & (part.gene_b == a)))]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    A(f"| {ds} | {cov} | {a}–{b} | {r.partial_rho:.3f} | {r.fdr:.2e} |")
        A("")

    # ---------- P6 scRNA ----------
    A("## P6 — within-cell-type scRNA (Census)")
    A("")
    if hmeans is not None and len(hmeans):
        A("Cell-type means (normal human lung epithelium, claim genes):")
        A("")
        cols = ["cell_type", "n_cells"] + C.CLAIM_GENES_H
        sub = hmeans[hmeans.disease == "normal"][cols].copy()
        if len(sub):
            A("| " + " | ".join(cols) + " |")
            A("|" + "|".join(["---"] * 2 + ["---:"] * (len(cols) - 2)) + "|")
            for _, r in sub.sort_values("n_cells", ascending=False).head(12).iterrows():
                A("| " + " | ".join(
                    [str(r.cell_type), str(int(r.n_cells))] +
                    [f"{r[g]:.2f}" if np.isfinite(r[g]) else "NA" for g in C.CLAIM_GENES_H]
                ) + " |")
        A("")
    if hwithin is not None and len(hwithin):
        A("Within-type Spearman for the pairs the claim cares about:")
        A("")
        A("| disease | cell type | n | pair | ρ | FDR |")
        A("|---|---|---:|---|---:|---:|")
        focus = [("ELF3", "TACSTD2"), ("ELF3", "CLDN4"), ("GRHL1", "TACSTD2"),
                 ("KLF4", "TACSTD2"), ("TFAP2A", "TACSTD2"),
                 ("NKX2-1", "TACSTD2"), ("NKX2-1", "CLDN4"), ("NKX2-1", "ELF3"),
                 ("TACSTD2", "CLDN4")]
        types = ["basal cell", "club cell", "ciliated cell", "lung secretory cell",
                 "pulmonary alveolar type 2 cell", "pulmonary alveolar type 1 cell",
                 "malignant cell", "epithelial cell of lung"]
        for dis in ["normal", "lung adenocarcinoma"]:
            for ct in types:
                for a, b in focus:
                    r = hwithin[(hwithin.disease == dis) & (hwithin.cell_type == ct) &
                                (((hwithin.gene_a == a) & (hwithin.gene_b == b)) |
                                 ((hwithin.gene_a == b) & (hwithin.gene_b == a)))]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    A(f"| {dis} | {ct} | {int(r.n_cells)} | {a}–{b} | "
                      f"{r.spearman_rho:.3f} | {r.fdr:.2e} |")
        A("")

    # ---------- P4 ----------
    A("## P4 — experimental NKX2-1 loss raises MODULE and TARGET")
    A("")
    A("Positive log2FC = higher in the NKX-low / KD arm (claim direction).")
    A("")
    if pert is not None and len(pert):
        A("| dataset | comparison | gene | n_low | n_high | log2FC | Cliff's δ | FDR |")
        A("|---|---|---|---:|---:|---:|---:|---:|")
        for _, r in pert.iterrows():
            nlo = r.get("n_nkx_low", r.get("n_low", ""))
            nhi = r.get("n_nkx_high", r.get("n_high", ""))
            d = r.get("cliffs_delta", np.nan)
            fdr = r.get("fdr", np.nan)
            lfc = r.get("log2fc_low_minus_high", np.nan)
            A(f"| {r.dataset} | {r.comparison} | {r.gene} | {nlo} | {nhi} | "
              f"{lfc if pd.notna(lfc) else 'NA'} | "
              f"{'' if pd.isna(d) else f'{d:.3f}'} | "
              f"{'' if pd.isna(fdr) else f'{fdr:.2e}'} |")
    else:
        A("No perturbation contrast table was produced. See `logs/perturbation_notes.txt`.")
    A("")

    # ---------- P5 mouse scRNA ----------
    A("## P5 — mouse")
    A("")
    if mmeans is not None and len(mmeans):
        cols = ["cell_type", "disease", "n_cells"] + [g for g in C.CLAIM_GENES_M if g in mmeans.columns]
        A("| " + " | ".join(cols) + " |")
        A("|" + "|".join(["---"] * 3 + ["---:"] * (len(cols) - 3)) + "|")
        show = mmeans.sort_values("n_cells", ascending=False).head(15)
        for _, r in show.iterrows():
            A("| " + " | ".join(
                [str(r.cell_type), str(r.disease), str(int(r.n_cells))] +
                [f"{r[g]:.2f}" if g in r and np.isfinite(r[g]) else "NA" for g in C.CLAIM_GENES_M if g in mmeans.columns]
            ) + " |")
        A("")
        A("Census mouse lung is dominated by endothelium, fibroblasts and AT2; airway")
        A("basal/club/ciliated cells are sparse or unlabelled. That is a power limitation")
        A("for the mouse observational arm, not evidence of absence.")
        A("")
    if mwithin is not None and len(mwithin):
        A(f"Mouse within-type pairs computed: {len(mwithin)}. See `census_mouse_within_type_spearman.csv`.")
        A("")

    A("## What this is not")
    A("")
    A("- Not a ChIP / motif / reporter assay. No claim that ELF3 (etc.) *bind* TACSTD2/CLDN4.")
    A("- Not TTF-1 IHC. NKX2-1 here is RNA.")
    A("- Not an ICI or TROP2-ADC outcome analysis.")
    A("- Not a test of any claim other than A10.")
    A("")
    A("## Limitations (pre-declared, restated with what actually happened)")
    A("")
    A("1. Bulk lung is composition-confounded. P6 exists to gate P1–P3.")
    A("2. ELF3/GRHL1/KLF4/TFAP2A are generic epithelial TFs; a positive P1/P2 can be")
    A("   'epithelium present' rather than a specific program. Matched-null and")
    A("   within-type tests are the guardrails.")
    A("3. Perturbation datasets are small; per-gene FDR is weak. Consistency matters more.")
    A("4. Cell-line KD is not primary lung.")
    A("5. Mouse Census lung has little labelled airway epithelium.")
    A("6. recount3 gene sums are coverage, not strict counts; rank statistics are invariant")
    A("   to per-gene length factors.")
    A("")
    A("## Files")
    A("")
    A("- Tables: `results/claim_A10/tables/`")
    A("- Figures: `results/claim_A10/figures/`")
    A("- Logs: `results/claim_A10/logs/`")
    A("- Spec: `claims/claim_A10.md`")
    A("")

    text = "\n".join(lines) + "\n"
    path = os.path.join(C.OUT, "REPORT.md")
    with open(path, "w") as fh:
        fh.write(text)
    print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
