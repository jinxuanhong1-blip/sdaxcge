# Claim A10 — pre-registered specification

**Status:** written and committed *before* any analysis output was inspected.
**Scope:** Claim A10 only. No other claim is addressed by this work.

## 1. The claim as received

The claim was supplied in shorthand:

> `ELF3/GRHL1/KLF4/TFAP2A vs TACSTD2/CLDN4; NKX2-1 down. Public human+mouse lung.`

No longer-form statement of "Claim A10" exists in this repository, so the shorthand is
operationalised below. **If the intended claim differs from this operationalisation, the
conclusions in `results/claim_A10/REPORT.md` may not apply.** The operationalisation is stated
explicitly so that it can be corrected and the analysis re-run.

## 2. Operationalised claim

> In lung epithelium there is a coherent transcription-factor module — **ELF3, GRHL1, KLF4,
> TFAP2A** — whose activity is coordinately associated with a cell-surface / tight-junction
> program — **TACSTD2** (TROP2) and **CLDN4** — and this coupled TF→target program marks a cell
> state in which the lung lineage factor **NKX2-1** is **down**-regulated. The relationship is
> expected to hold in both human and mouse lung.

Gene sets:

| role | human symbols | mouse symbols |
|---|---|---|
| TF module (`MODULE`) | ELF3, GRHL1, KLF4, TFAP2A | Elf3, Grhl1, Klf4, Tfap2a |
| target program (`TARGET`) | TACSTD2, CLDN4 | Tacstd2, Cldn4 |
| antagonist (`NKX`) | NKX2-1 | Nkx2-1 |

## 3. Testable predictions

| id | prediction | rationale |
|---|---|---|
| **P1** | The four `MODULE` TFs are positively inter-correlated across lung samples/cells, more than a random gene set of matched expression. | Required for "module" to be a meaningful unit rather than four unrelated genes. |
| **P2** | `MODULE` (and each TF individually) correlates positively with `TACSTD2` and `CLDN4`. | The "vs" in the claim is read as "drives / is coupled to". |
| **P3** | `NKX2-1` correlates **negatively** with `MODULE` and with `TARGET`. | "NKX2-1 down" in the same cell state. |
| **P4** | Experimental **loss of NKX2-1** *increases* `MODULE` and `TARGET` expression. | Directional/causal version of P3. Correlation alone cannot establish it. |
| **P5** | P1–P4 reproduce in **mouse** lung. | The claim asserts human+mouse. |
| **P6** | The associations are **not solely** a by-product of cell-type composition (airway vs alveolar mixing) — they survive analysis *within* an epithelial cell type in single-cell data. | Bulk co-expression across a tissue with variable epithelial composition can manufacture all of P1–P3 trivially. |

## 4. Decision rules (fixed in advance)

Statistics: Spearman rank correlation (robust to the non-normality of expression data);
Benjamini–Hochberg FDR across all gene-pair tests reported per dataset.

Every observed correlation is additionally converted to a **background percentile**: the fraction
of 200,000 random gene–gene pairs (drawn from expressed genes in the same dataset, so the
percentile absorbs whatever global co-expression structure exists) with a correlation at least as
extreme. This guards against calling `rho = 0.35` "strong" in a dataset where the median
gene-pair correlation is already 0.3.

Per-prediction verdicts:

* **SUPPORTED** — effect in the predicted direction, FDR < 0.05, |rho| >= 0.3 (or |log2FC| >= 0.5
  for perturbation), and background percentile in the top/bottom 10% of random pairs;
  reproduced in **>= 2 independent datasets per species claimed**.
* **PARTIALLY SUPPORTED** — predicted direction holds for some genes/datasets but not others, or
  effect is significant but small (|rho| < 0.3), or holds in one species only.
* **NOT SUPPORTED** — no consistent effect in the predicted direction (effects near zero, or
  direction inconsistent across datasets).
* **REFUTED** — significant effect in the **opposite** direction, consistently, across
  independent datasets.

The overall Claim A10 verdict is the conjunction: the claim is only **SUPPORTED** if P1, P2, P3
and P5 are supported; P4 upgrades it from associational to causal; failure of P6 downgrades any
support to "explained by cell composition".

## 5. Data (all public, no controlled access)

Selected before results were inspected. Chosen for (a) size sufficient for correlation, (b) direct
NKX2-1 perturbation for P4, (c) both species, (d) single-cell for P6.

**Human, bulk, observational**
* GTEx v8 lung (normal), uniformly reprocessed by **recount3** — n ~ 655.
* TCGA-LUAD (lung adenocarcinoma + adjacent normal), **recount3** — n ~ 600.

**Human, perturbation**
* **GSE129340** — NKX2-1/TTF-1 siRNA knockdown in NCI-H441 (NKX2-1-positive LUAD) and NCI-H209.
* **GSE229541** — NKX2-1 dosage manipulation in human LUAD lines.

**Mouse, perturbation / observational bulk**
* **GSE129583** — conditional *Nkx2-1* deletion in mouse lung alveolar epithelium (Little et al.).
* **GSE115899**, **GSE145152**, **GSE188435** — *Nkx2-1*-negative vs -positive mouse lung
  adenocarcinoma models (Snyder/Camolotto lineage-switch series).

**Single cell, both species (P6, and the mouse observational arm)**
* **CELLxGENE Census** (stable release), `tissue_general == "lung"`, primary data only:
  human ~6M cells, mouse ~2.2e5 cells.

## 6. Pre-declared threats to validity

1. **Composition confounding.** All bulk correlations are confounded by epithelial fraction. P6
   exists specifically to test this and its result gates the interpretation of P1–P3.
2. **Circularity of "module".** ELF3, GRHL1, KLF4, TFAP2A are all broadly epithelial TFs; any
   positive P1/P2 result may reflect generic "epithelium present" rather than a specific program.
   The random-background percentile partially addresses this; an *epithelium-matched* background
   is also reported.
3. **Direction.** P1–P3, P5 are correlational. Only P4 speaks to causation, and only in the
   NKX2-1 → module direction, not module → target.
4. **Cell lines are not lung.** GSE129340/GSE229541 are cancer cell lines; concordance with
   primary tissue is not assumed.
5. **Perturbation datasets are small** (n = 2–5 per arm is typical); per-gene FDR from such
   designs is weak evidence in isolation, hence the requirement for consistency across datasets.
6. **Orthology.** Mouse/human one-to-one orthology is assumed for these seven genes, which is the
   case, but mouse *Tfap2a* and *Grhl1* expression in lung is low, so mouse tests for those two
   may be underpowered — this will be reported as a power limitation, not as evidence of absence.
