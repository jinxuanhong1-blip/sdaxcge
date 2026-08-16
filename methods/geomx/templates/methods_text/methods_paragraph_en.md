# Fill-in methods paragraph (English)

Copy into the manuscript and replace every `[…]`. Do not leave a bracket standing. Numbers come from `results/01_qc_*.tsv`, `02_norm_choice.txt`, `03_consensus_correlation.txt`, `04_reporting_block.txt`, `05_reliability.tsv`, and `06_cox_patient.tsv`.

---

## Tissue and GeoMx WTA

Formalin-fixed paraffin-embedded `[resection / biopsy / TMA]` sections from `[n]` patients with `[NSCLC adenocarcinoma / squamous / mixed]` were profiled on the GeoMx Digital Spatial Profiler using the Human Whole Transcriptome Atlas (~18,000 protein-coding targets). Morphology markers were pan-cytokeratin (tumor), CD45 (leukocytes), CD68 (macrophages) and SYTO 13 (nuclei). `[k]` regions of interest (ROIs) were placed per tumor (range `[–]`) and each ROI was segmented into the three compartments, yielding `[N_AOI_in]` areas of illumination (AOIs) before quality control. Public digital count conversion files were obtained from GEO `[GSE271689 / GSE292098]`; no restricted-access counts were used.

## Quality control

AOIs were filtered with GeomxTools using per-segment nuclei and area floors (tumor: `[100]` nuclei / `[5000]` µm²; CD45: `[50]` / `[2000]`; CD68: `[20]` / `[1000]`), a minimum of `[1000]` raw reads, ≥ `[80]`% trimmed / stitched / aligned reads, ≥ `[50]`% sequencing saturation, negative-probe geometric mean ≥ `[10]`, and NTC well counts ≤ `[1000]`. The limit of quantification in AOI *i* was `LOQ_i = max(2, geomean(NegProbe_i) × geoSD(NegProbe_i)²)`. AOIs detecting < `[1]`% of targets above LOQ were removed. Genes were retained if detected above LOQ in ≥ `[5]`% of AOIs **within at least one segment** (union, used for within-compartment analyses); cross-compartment contrasts were restricted to the intersection of the three per-segment gene lists. After QC, `[N_AOI_out]` AOIs and `[N_gene]` genes remained (`[n_tumor]` tumor, `[n_CD45]` CD45, `[n_CD68]` CD68). Counts of AOIs dropped at each step, split by segment, are in Supplementary Table `[ ]`.

## Normalization and batch

`[TMM / Q3 / quantile]` was declared as the primary normalization **before** differential expression was inspected. `[Q3 and quantile]` were run as a sensitivity analysis and compared on relative-log-expression plots and principal components coloured by segment, slide, AOI area and detection rate (Supplementary Figure `[ ]`). Cross-compartment contrasts used jointly normalized counts; within-compartment analyses `[used the same matrix / were re-normalized within segment]`. Slide was `[included as a random effect / shown not to be confounded with cohort]`. Batch-corrected matrices (RUV4 / `removeBatchEffect`) were used for visualization only and were not input to inference.

## Differential expression and mixed models

The independent experimental unit was the patient. Genome-wide differential expression used limma-voom with `duplicateCorrelation` blocked on patient (two rounds when the consensus correlation exceeded 0.5). The consensus correlation was `[ρ]`. *P* values were Benjamini–Hochberg–adjusted **within each contrast**. Co-existing compartments (tumor vs CD45 vs CD68) were additionally modelled with a linear mixed model with a patient random intercept and a random segment slope; Kenward–Roger (or Satterthwaite) degrees of freedom were used. A patient × segment pseudobulk analysis was the prespecified sensitivity analysis; direction disagreements with the voom fit were resolved in favour of the pseudobulk. A *t*-test treating AOIs as independent was not used.

## TACSTD2 and CLDN4

`TACSTD2` and `CLDN4` were prespecified. For each gene we report (i) the fraction of AOIs above LOQ per segment, (ii) the mixed-model tumor–immune log2 difference with 95% confidence interval and the fraction of patients in whom the difference was positive, and (iii) the within-patient ICC and the number of ROIs required for a patient-level reliability ≥ 0.8. Immune-AOI signal was tested against an epithelial spillover index (`EPCAM`, `KRT8`, `KRT18`, `KRT19`, `CDH1`); we report the correlation and the estimate restricted to the lowest spillover tertile. `[SpatialDecon with is_pure_tumor set from PanCK+ AOIs / IHC on the same blocks]` was used as an orthogonal control. The tumor-versus-immune log-fold change is interpreted as a lower bound on compartment specificity.

## Patient-level scores and overall survival

AOI-level expression was aggregated to one value per patient per segment by `[BLUP shrinkage / median / nuclei-weighted mean]`. Four prespecified exposures were tested: tumor-segment expression, immune-segment expression, the tumor–immune difference (specificity), and the within-patient standard deviation across ROIs (heterogeneity). Reliability of the patient mean is reported next to every hazard ratio. The primary survival model was a patient-level Cox model, exposures scaled per standard deviation, adjusted for `[age, sex, stage, histology]` and stratified by `[cohort]`. Proportional hazards were checked with scaled Schoenfeld residuals; functional form was checked with a 3-knot restricted cubic spline. AOI-level Cox models with a patient cluster (robust sandwich variance) or a shared frailty term were sensitivity analyses only; they duplicate the event and do not increase information about it. Cut points `[were not used / were the training-cohort median / tertile, locked in the Yale cohort and applied unchanged to UQ and Athens]`. `[n]` patients and `[e]` deaths entered the primary model (`[e/p]` events per variable). Genome-wide survival screening, if performed, is hypothesis-generating.

## Software

R `[version]`; GeomxTools `[ ]`; standR `[ ]`; edgeR/limma `[ ]`; lme4/lmerTest `[ ]`; survival/coxme `[ ]`. Session information is in Supplementary File `[ ]`.
