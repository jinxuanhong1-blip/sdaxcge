# A11 — CD47 / Galectin / Nectin / TGFB vs TACSTD2-high in public lung cancer

**Headline.** After ABSOLUTE purity adjustment, **CD47 tracks TACSTD2 in TCGA-LUAD (partial ρ = 0.30) and does not track it in TCGA-LUSC (partial ρ = 0.01).** The LUAD association sits at the pre-specified |ρ| ≥ 0.3 threshold, is unchanged by purity, survives a post-hoc housekeeping residualization (ρ = 0.29), and is consistent inside every purity tertile. It is **absent, with a confidence interval that excludes even a small effect, in the pre-registered LUSC replication.** Under the rules written before any result was inspected, CD47 vs TACSTD2 is **not established as a general NSCLC finding.**

Purity was the confounder this analysis was built to kill. It is not the confounder. TACSTD2 is essentially uncorrelated with ABSOLUTE purity in both histologies (LUAD ρ = 0.04, p = 0.34; LUSC ρ = 0.01, p = 0.90). Adjusting for it moves almost no estimate off the diagonal. The thing that actually decides the claim is histology, not stroma.

---

## What was asked, and what was actually tested

TACSTD2 (TROP2) is the target of sacituzumab govitecan and datopotamab deruxtecan. The question is whether TROP2-high lung tumours also run a distinct CD47 / galectin / nectin / TGF-β immune-evasion programme — the kind of co-expression that would justify a specific ADC + innate-checkpoint combination.

The analysis was pre-registered in `00_analysis_plan.md` (commit `655308f`) before any result was inspected. The decision rule, copied here because it is doing real work:

> An axis counts as genuinely associated with TACSTD2 only if it is (i) significant after purity adjustment in bulk LUAD, (ii) directionally consistent in LUSC, and (iii) directionally consistent in within-malignant-cell single-cell pseudobulk. Anything less is reported as **not established**. |ρ| ≥ 0.3 is the threshold for "worth acting on".

Bulk is A1–A4. Single-cell is S1–S3 on GSE131907 (Kim et al. 2020, LUAD; 208,506 cells). GSE127465 was planned as an sc replication and was not run.

## Data (all public, all real, none simulated)

| Role | Source | Used |
|---|---|---|
| Discovery | TCGA-LUAD STAR TPM, UCSC Xena GDC hub | 528 primary tumours (515 with ABSOLUTE) |
| Replication | TCGA-LUSC, same pipeline | 501 primary tumours (493 with ABSOLUTE) |
| Purity | TCGA ABSOLUTE, GDC PanCanAtlas | 10,787 calls; matched on the 15-character barcode |
| Adjacent normal (context only) | TCGA-LUAD / LUSC type-11 | 59 / 51 |
| scRNA, LUAD | GSE131907, Kim et al. 2020, author cell labels | 208,506 cells; 10 primary + 11 LN/bronchus-met + 10 brain-met patients with ≥20 malignant cells |

Expression is the Xena `log2(TPM+1)` matrix. Primary tumours only (`-01`). LUAD and LUSC were never pooled. Gene symbols from GENCODE v36. The CD47 axis is CD47, SIRPA, THBS1. The other three axes and a 10-gene housekeeping panel were run on the same pipeline so that CD47 is interpreted against a real noise floor, not against zero.

## The confounder that did not show up

TACSTD2 is an epithelial gene, so the prior was that "TACSTD2-high" in bulk is mostly "more tumour / more epithelium." That prior is **false for TACSTD2 in TCGA lung**:

| | LUAD ρ (95% CI) | LUSC ρ (95% CI) |
|---|---|---|
| TACSTD2 vs ABSOLUTE purity | 0.04 (−0.04, 0.13) | 0.01 (−0.08, 0.09) |
| TACSTD2 vs RNA epithelial score | 0.46 (0.39, 0.53) | 0.16 (0.08, 0.25) |
| TACSTD2 vs 10-gene housekeeping score | 0.32 (0.24, 0.40) | 0.01 (−0.08, 0.10) |

DNA purity and TACSTD2 are independent. RNA epithelial score does track TACSTD2 in LUAD — that is differentiation / epithelial transcription, not tumour fraction. The housekeeping score also tracks TACSTD2 in LUAD only, which is why A4 exists and why a post-hoc HK residualization is reported below, labelled as post-hoc.

Figure `03` is the whole purity story in one panel: every gene in the primary panel sits on the y = x line of naive ρ vs purity-adjusted ρ. Purity is not doing anything here.

## CD47 vs TACSTD2 — the result

### LUAD (discovery)

| Estimator | ρ or δ | 95% CI | n | q (BH, 30-gene panel) |
|---|---|---|---|---|
| Spearman, naive | 0.303 | 0.224, 0.379 | 528 | < 10⁻¹¹ |
| Partial Spearman \| ABSOLUTE purity | 0.304 | 0.223, 0.381 | 515 | < 10⁻¹¹ |
| Partial Spearman \| epithelial score | 0.281 | — | 528 | — |
| Partial Spearman \| housekeeping score (post-hoc) | 0.284 | 0.203, 0.360 | 528 | — |
| Partial Spearman \| purity + HK (post-hoc) | 0.288 | 0.207, 0.366 | 515 | — |
| Tertile high vs low, Cliff's δ | 0.382 | — | 176 vs 176 | < 10⁻¹¹ |

Inside ABSOLUTE-purity tertiles (A3, pre-specified):

| Stratum | Purity range | ρ (95% CI) | n |
|---|---|---|---|
| Low | 0.09–0.35 | 0.26 (0.11, 0.39) | 172 |
| Mid | 0.36–0.53 | 0.34 (0.20, 0.47) | 168 |
| High | 0.54–1.00 | 0.33 (0.19, 0.45) | 175 |

The LUAD association is not a low-purity artefact, is not an epithelial-content artefact, and is not a housekeeping-programme artefact. It is a real among-tumour correlation of two epithelial genes, of modest size, just at the pre-specified actionability line.

### LUSC (pre-registered replication)

| Estimator | ρ or δ | 95% CI | n | q |
|---|---|---|---|---|
| Spearman, naive | 0.007 | −0.080, 0.095 | 501 | 0.90 |
| Partial Spearman \| purity | 0.006 | −0.083, 0.094 | 493 | 0.93 |
| Partial Spearman \| purity + HK | 0.004 | −0.085, 0.092 | 493 | — |
| Tertile high vs low, Cliff's δ | −0.004 | — | 167 vs 167 | 0.95 |

The LUSC confidence interval excludes |ρ| = 0.10. This is not "underpowered replication." A 0.30 effect of the LUAD size would have been seen at essentially p = 0. The association is absent.

A3 in LUSC is more interesting than the pooled null: low-purity ρ = 0.14, mid 0.09, **high-purity ρ = −0.18 (95% CI −0.32, −0.03; p = 0.019, q = 0.057).** The FDR-adjusted q is 0.057, so this is a hint, not a finding. Directionally it is the opposite of LUAD, and it appears in the stratum where tumour-cell RNA should dominate. If anything, high-purity LUSC argues against a shared CD47–TACSTD2 programme.

### The rest of the CD47 axis

SIRPA (the myeloid receptor) is weakly positive in LUAD after purity (ρ = 0.17, q = 0.0003) and **collapses after housekeeping residualization (ρ = 0.02, p = 0.59).** THBS1 is null in both histologies. There is no evidence for a co-regulated CD47–SIRPA–THBS1 module. The LUAD signal is CD47 the ligand, alone.

### Tumour vs adjacent normal (context, not the question)

CD47 is **lower** in primary tumour than in adjacent lung: LUAD log2FC = −0.71, Cliff's δ = −0.56; LUSC log2FC = −1.23, δ = −0.84. TACSTD2 itself is almost unenriched in LUAD tumour vs adjacent (log2FC = +0.12) and only moderately enriched in LUSC (+0.80). "TACSTD2-high LUAD has more CD47" is a ranking among tumours. It is not "tumours express more CD47 than lung." Any combination rationale that starts from CD47 being a tumour-enriched don't-eat-me signal in TROP2-high NSCLC is starting from the wrong baseline.

## Housekeeping calibration (A4) — why LUAD still needs a grain of salt

In LUAD, TACSTD2 correlates with ACTB (purity-adjusted ρ = 0.40) and UBC (0.39) **at least as strongly as with CD47 (0.30).** Several other housekeepers are also positive (PPIA 0.21, B2M 0.17, GAPDH 0.16). In LUSC the same housekeepers sit around zero, matching CD47.

That is why the post-hoc HK residualization was run. CD47 survives it (0.29). SIRPA does not. So CD47 is not *just* the LUAD TACSTD2-high expression programme — but TACSTD2-high LUAD is a more transcriptionally intense / well-differentiated state, and any single-gene ρ of 0.30 in that cohort has to be read against that floor. We do not have a genome-wide null; ten housekeepers are a floor, not a proper empirical null. This is a real limit, not a caveat for the appendix.

## Single cell (GSE131907, LUAD) — criterion (iii)

Two facts from S1 are solid and do not depend on patient n:

- **TACSTD2 is epithelial.** In tumour-tissue cells, mean log2TPM is 2.09 in epithelium (82% expressing) and ≤ 0.26 in every other compartment.
- **CD47 is highest in epithelium too** (mean 1.55, 76% expressing), with a real myeloid second place (0.87, 55%). The bulk LUAD correlation *can* be a tumour-cell relationship. It is not required to be stroma.

S2 is the test that matters, and it is smaller than the first draft of this analysis made it look.

Kim annotates primary-tumour epithelium as tS1/tS2/tS3, and uses the label "Malignant cells" only at metastatic sites. Brain mets run low for both genes (patient-mean TACSTD2 1.12, CD47 0.80) relative to primary tLung (2.31, 1.44). Pooling sites therefore correlates **site** as well as biology. The mixed-site ρ = 0.72 (n = 21, q = 0.006) is real as a number and is not the confirmatory test.

| Slice | n patients | n cells | CD47 ρ (95% CI) | p | q (30-gene) |
|---|---|---|---|---|---|
| Primary tLung, tS1/2/3 | 10 | 6,347 | 0.55 (−0.12, 0.88) | 0.098 | 0.71 |
| LN / bronchus mets | 11 | 9,361 | 0.75 (0.28, 0.93) | 0.007 | 0.19 |
| Brain mets | 10 | 15,423 | 0.61 (−0.03, 0.90) | 0.060 | 0.72 |
| Mixed tumour sites, no brain | 21 | 15,708 | 0.72 (0.42, 0.88) | 2×10⁻⁴ | 0.006 |

Direction is the same as LUAD bulk in every slice. That satisfies the weak reading of criterion (iii). The strong reading — a within-primary, FDR-significant, |ρ| ≥ 0.3 confirmation — is **not met**. n = 10 primary patients is not a confirmation, it is a compatible underpowered point estimate. SIRPA is null in every slice. NECTIN4, the bulk replicator, is underpowered here (primary ρ = −0.10; mixed 0.28).

S3 (tLung cell-type fractions vs malignant TACSTD2, n = 10) is too small to interpret and is in the table for completeness.

GSE127465 (7 patients) was the planned sc replication. It is not in this folder. Seven patients would not have rescued the primary-n problem.

## The other three axes, briefly and honestly

These were pre-specified. They are not the focus of this update.

**Nectin / NECTIN4 is the only gene that clearly replicates.** Purity-adjusted ρ with TACSTD2: LUAD 0.49, LUSC 0.53. Both CIs sit well above 0.3. HK residualization does not touch it (0.49 / 0.53). NECTIN4 is a known epithelial / squamous-adjacent nectin; co-expression with TACSTD2 is biologically unsurprising and is the one result in this folder that would survive a sceptical review on its own. It is not CD47.

**LGALS3** is directionally consistent (LUAD 0.44, LUSC 0.28). LUAD is above the 0.3 line; LUSC is just under. After HK: 0.37 / 0.28. Closest runner-up to a replicated ligand association. Still fails the strict |ρ| ≥ 0.3-in-both rule.

**TGFB1** is positive in both (LUAD 0.32, LUSC 0.17) but the LUSC effect is small. TGFB2 flips sign (LUAD +0.33, LUSC −0.11). The F-TBRS activity score is +0.23 in LUAD and −0.12 in LUSC. There is no shared TGF-β programme.

**Immune-cell receptors** (TIGIT, CD226, CD96, LAG3, HAVCR2) are near-zero in LUAD and negative in LUSC, tracking the LUSC TACSTD2-high / immune-low composition pattern (LUSC TACSTD2 vs immune score ρ = −0.26). Those are TME-fraction effects, which is why the ligand-only modules were pre-specified.

## Verdict against the pre-registered rule

| Criterion | CD47 | NECTIN4 | LGALS3 |
|---|---|---|---|
| (i) LUAD, purity-adjusted, q < 0.05 and \|ρ\| ≥ 0.3 | yes (0.30) | yes (0.49) | yes (0.44) |
| (ii) LUSC, same direction | **no (0.01)** | yes (0.53) | yes (0.28), below 0.3 |
| (iii) malignant-cell scRNA, same direction | yes, every site slice; primary n=10 p=0.10 | no (primary −0.10; mixed 0.28) | weak (mixed 0.34, p=0.14) |
| **Call** | **not established** as NSCLC; LUAD-restricted bulk is real and sc-compatible | **bulk-replicated**; sc does not confirm (underpowered, not the question) | not established |

CD47 vs TACSTD2 is a **LUAD-restricted, modest, purity-independent bulk correlation that is compatible with malignant-cell co-expression and is absent in LUSC.** It is not a NSCLC finding. It is not a CD47-axis finding (SIRPA/THBS1 do not travel with it). It is not tumour-vs-normal enrichment.

A LUAD-only claim is the most the public data will currently carry. It was not the pre-registered claim.

## What this does and does not say about combinations

It does **not** say "give a CD47 blocker with a TROP2 ADC in NSCLC." LUSC — a large fraction of the TROP2-ADC eligible population — has no CD47–TACSTD2 relationship at all, and CD47 RNA is lower in tumour than in adjacent lung in both histologies.

It does say that if anyone is going to look for a CD47 × TROP2 combination window, the only public bulk signal is in adenocarcinoma, it is ρ ≈ 0.3, and it has not been shown to live in the malignant cell. That is a hypothesis-generating observation, not a translational rationale.

NECTIN4 × TACSTD2 is the association in this dataset that actually looks like biology rather than cohort structure. That is a different drug axis.

## Limitations (the ones that can change the answer)

1. **Bulk RNA cannot assign the CD47 transcript to the malignant cell.** Purity independence, HK residualization, and S1 (CD47 highest in epithelium) argue against a pure stroma story. S2 in primary tumours is directionally consistent and underpowered (n = 10). That is not the same as proven tumour-cell co-expression.
2. **No protein, no spatial, no signalling.** RNA ρ = 0.3 is not CD47 protein on the TROP2-high cell, and is not "don't-eat-me" activity.
3. **LUSC is a different disease, not a failed experiment.** A histology-restricted LUAD association can be true. It is not what was asked, and it is not what the decision rule accepts.
4. **The LUAD housekeeping floor is real.** Two housekeepers out-correlate CD47. Residual technical or differentiation structure remains possible. We did not run a genome-wide empirical null.
5. **ABSOLUTE is missing for 13 LUAD and 8 LUSC primaries.** Naive and adjusted n differ slightly; the naive LUAD ρ on the full 528 is 0.303 vs 0.304 on the 515, so this is not material.
6. **One pipeline, one quantification (STAR TPM, GENCODE v36).** No RNA-seq batch of TCGA was re-processed here.
7. **Adjacent-normal n is small** (59 / 51). The tumour-vs-normal CD47 depletion is large enough that n is not the issue; the comparison is still same-patient-unpaired and is context only.
8. **No LUAD transcriptional subtype** (TRU / PI / PP) was available in the Xena clinical table we used. If the LUAD CD47–TACSTD2 correlation is just "TRU tumours are high for both," that would be a different interpretation and is untested.
9. **GSE131907 is LUAD.** It cannot speak to the LUSC null. Kim's "Malignant cells" label is metastatic; primary epithelium is tS1/tS2/tS3. Analyses that ignore that split will report an inflated CD47–TACSTD2 ρ.
10. **GSE127465 was not analysed.** Planned sc replication, 7 patients, would not have fixed the n problem.

## How to reproduce

```
.venv/bin/python src/prep_bulk.py
.venv/bin/python src/analyze_bulk.py
.venv/bin/python src/analyze_bulk_sensitivity.py
.venv/bin/python src/extract_sc_panel.py
.venv/bin/python src/analyze_sc.py
.venv/bin/python src/plot_a11.py
```

Inputs are not in git (`data/` is gitignored; ~320 MB bulk + 3.4 GB sc). URLs are in `data/bulk/urls.txt` and the GEO suppl paths in the plan. Outputs live in `results/w200/A11/tables/` and `figures/`.

## Files

- `00_analysis_plan.md` — pre-registered, committed before results
- `tables/10_bulk_gene_level.csv` — every primary-panel and housekeeping gene, naive and adjusted
- `tables/11_bulk_module_level.csv` — ligand modules, F-TBRS, compartment scores, purity
- `tables/12_bulk_purity_strata.csv` — A3
- `tables/13_bulk_tumor_vs_normal.csv`
- `tables/14_bulk_cd47_sensitivity.csv` — post-hoc HK residualization
- `figures/01_cd47_vs_tacstd2_scatter.{png,pdf}` — the money plot
- `figures/02_forest_naive_vs_purity.{png,pdf}`
- `figures/03_purity_shift_and_housekeepers.{png,pdf}`
- `figures/04_cd47_purity_strata.{png,pdf}`
- `figures/05_tumor_vs_normal.{png,pdf}`
- `tables/20_sc_celltype_expression.csv` — S1
- `tables/21_sc_malignant_pseudobulk.csv` — S2, site-stratified
- `tables/21_sc_malignant_pseudobulk_samples.csv` — per-patient means
- `tables/22_sc_tme_composition.csv` — S3
- `figures/06_sc_malignant_pseudobulk.{png,pdf}`
