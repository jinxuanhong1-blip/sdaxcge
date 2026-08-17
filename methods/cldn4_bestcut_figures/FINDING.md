# FINDING — CLDN4 best-cut extra figure pack

**Additive redraw only.** This folder does **not** re-audit, re-download, or re-estimate any public extra. It takes already-significant public CLDN4 extras and redraws them as a single high-contrast slide pack.

Locked numbers: `locked_numbers.tsv`. Figures: `figures/`. Reproduce: `python3 methods/cldn4_bestcut_figures/make_figures.py`.

No p-value was invented. Display rounding on slides matches the source write-ups. Exact floats stay in the TSV.

Fisher-z 95% whiskers appear only on the forest for rows that had no published CI (GSE218989, PR290 combos, GSE265899). Those intervals are a display transform of the **already published** n and ρ. They are not a new test.

---

## Extra figures (9 PNG / PDF pairs)

| File | Extra | Quoted n / effect / p |
|---|---|---|
| `fig01_forest_spearman` | One forest of Spearman extras | all rows below |
| `fig02_gse50927_ifn_nes` | GSE50927 IFN NES | n=1 vs 1; IFN-γ NES **+1.508** FDR **0.005**; IFN-α **+1.587** FDR **0.010**; MHC-I **+1.423** FDR **0.050**; Cldn4 log2FC **−6.061** |
| `fig03_cptac_lscc_protein` | CPTAC LSCC protein −0.43 | n=**78**; ImmuneScore ρ=**−0.432** p=**7.9e-05**; GEP18 **−0.461** p=**2.1e-05**; CD8A **−0.437** p=**6.4e-05** |
| `fig04_gse218989_cd8` | GSE218989 CLDN4 vs CD8 | n=**355**; ρ=**−0.172** p=**0.0011**. Response is null (AUC 0.47, MWU p=0.40) and is shown as the published null, not a hit |
| `fig05_pr290_combos` | PR290 combos | T/NK author %pos n=**43** ρ=**−0.479** p=**0.00152**; CXCL13+ n=**60** ρ=**−0.425** p=**0.00121** |
| `fig06_gse285029_ifn` | GSE285029 CLDN4 vs IFN | n=**234**; IFN-compact ρ=**+0.204** p=**0.00175**; MHC-I ρ=**+0.263** p=**4.6e-05** |
| `fig07_gse265899_spatial` | spatial GSE265899 | all AOI n=**95** ρ=**−0.821** p=**2.4e-24**; tumor AOI n=**48** ρ=**−0.399** p=**0.00492** |
| `fig08_gse312098_cldn4_down` | GSE312098 CLDN4 down | n=3 vs 3; CLDN4 log2FC **−0.86** p=**1.7e-05** q=**0.0043**; IFN MW p=**2.6e-08** |
| `fig09_multipanel_extras` | One multi-panel extra figure set | A–G of the seven extras, each quoting n / ρ or NES / p |

`fig01` is the one forest. `fig09` is the one multi-panel extra figure set. `fig02`–`fig08` are the individual high-contrast slides.

---

## Source extras (taken as given)

| Extra | Source PR | Source table |
|---|---|---|
| GSE50927 IFN NES | [#287](https://github.com/jinxuanhong1-blip/sdaxcge/pull/287) | `methods/cldn4_ko_gsea/tables/gsea_headline.tsv` |
| CPTAC LSCC protein | [#289](https://github.com/jinxuanhong1-blip/sdaxcge/pull/289) | `methods/cptac_lusc_cldn4_protein/tables/associations.tsv` |
| GSE218989 CLDN4 vs CD8 | [#292](https://github.com/jinxuanhong1-blip/sdaxcge/pull/292) | `methods/gse218989_cldn4_ici/tables/stats.tsv` |
| PR290 combos | [#290](https://github.com/jinxuanhong1-blip/sdaxcge/pull/290) | `methods/scrna_cldn4_combo/tables/highlighted_combos.tsv` |
| GSE285029 CLDN4 vs IFN | [#239](https://github.com/jinxuanhong1-blip/sdaxcge/pull/239) | `results/w200/A11_GSE285029/correlations.csv` |
| spatial GSE265899 | [#251](https://github.com/jinxuanhong1-blip/sdaxcge/pull/251) | `results/leftover_spatial/tables/geomx_correlations.csv` |
| GSE312098 CLDN4 down | [#254](https://github.com/jinxuanhong1-blip/sdaxcge/pull/254) | `results/c_public_trop2_adc_analogs/tables/key_genes_all_contrasts.tsv` |

---

## Caveats kept on the slides

- GSE50927 is mouse **whole lung**, n=1 vs 1, not a lung-cancer KO.
- GSE218989 ICI response is **null**; the extra is CLDN4 vs CD8A.
- PR290 full-pool mixed grid is weaker (N=145, ρ=−0.137, p=0.129) and is not the claimed cut.
- GSE285029 CLDN4 vs CD8A is null (ρ=+0.067, p=0.31). Epithelial residual drops IFN-compact to partial p≈0.12.
- GSE265899 immune AOIs are NS (ρ=−0.180, n=47, p=0.23). All-AOI mixes compartment.
- GSE312098 is IMMU132 CRC CX-1 2 d, **not SKB264 and not lung**. Junction set is null (MW p=0.99).

Not SKB264. Not a re-cut of any prior slide.
