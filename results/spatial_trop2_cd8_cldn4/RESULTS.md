# Local TROP2–CD8 after a CLDN4 residual (public lung spatial)

Question: among public lung spatial assays that measure TACSTD2 (TROP2), CLDN4, and CD8 together, does a local CD8 association of TROP2 remain after residualizing CLDN4?

Answer from the matrices read in this run: no stable local association remains, because the total association is already near null. CosMx patient-level partial ρ at 50 µm is -0.012 (Wilcoxon p = 0.3125, n = 5). Visium ring-1 partial correlations in two LUAD series are centered on zero. Public Xenium lung panels opened here do not contain CLDN4, so that residual was not computed.

This is an observational product-method residual, not a causal mediation claim. Same-spot Visium correlation is reported because those whole-transcriptome matrices exist, and it is not the lead. Prior Visium same-spot tests are epithelial-versus-stroma composition (PRs #134, #196, #553).

## Lead: CosMx NSCLC, tumor cell vs nearby CD8 T cells

Dataset: He et al. 2022 CosMx SMI NSCLC FFPE, 8 sections / 5 patients, 960-plex. Mirror: Zenodo 15487520 `cosmx_lung` (counts, coordinates, author `cell_type`). Pixel size 0.18 µm/px, checked by a median nearest-neighbor distance inside 4–30 µm. Tumor cells are author labels beginning with `tumor`. CD8 cells are author `T CD8 naive` and `T CD8 memory`. Outcome is log1p of the CD8 count inside a radius. Exposure and mediator are log1p counts on the index tumor cell. Partial Spearman residualizes ranks on CLDN4 only. OLS mediation z-scores each variable inside the section. Replicate sections are Fisher-averaged (correlations) or mean-averaged (coefficients) to the patient, and the Wilcoxon is on the 5 patients. Cell-level p-values are not used as the claim: neighboring cells are not independent.

50 µm partial ρ(TACSTD2, local CD8 | CLDN4), patient median **-0.012**, Wilcoxon p = **0.3125** (n = 5). Total ρ median -0.024. Standardized total c median -0.027; direct c′ median -0.014 (Wilcoxon p = 0.3125). Sections with partial ρ < 0: 6/8.

TACSTD2 and CLDN4 do co-vary in tumor cells: the standardized a-path median across sections is 0.363 (8/8 sections positive). That coupling does not produce a stable indirect path to local CD8. A proportion (indirect/c) is stored only when |c| ≥ 0.02, which is 3 of 5 patients, and those proportions sit on small total effects. The patient with the largest |total ρ| is Lung12 (total ρ -0.128, partial ρ -0.080; after CLDN4 + KRT8 + library size, partial ρ -0.035). That single patient is not the five-patient result.

100 µm partial ρ median -0.004, Wilcoxon p = 0.4375. Direct c′ median -0.007.

The 50 µm author-CD8 count has median 0 in 7/8 sections (100 µm median ≤ 1 in 7/8 sections; maximum 100 µm median count is 6.0). The outcome is sparse, so a near-null means no detectable gradient, not a precisely estimated zero. With n = 5 patients the two-sided Wilcoxon floor is 0.0625. Section-level signs are not independent (Lung5 has three sections, Lung9 has two).

This does not rewrite the locked CLDN4 exclusion result (cytotoxic neighbor ratios at 50/100 µm). The outcome here is author CD8 T cells, and the exposure is TACSTD2 after CLDN4.

### Patients (50 µm)

| Patient | sections | total ρ | partial ρ \| CLDN4 | partial ρ \| CLDN4+KRT8+nCount | c | c′ | proportion |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Lung12 | 1 | -0.128 | -0.080 | -0.035 | -0.127 | -0.083 | 0.344 |
| Lung13 | 1 | -0.010 | 8.97e-04 | -0.008 | -0.011 | 4.75e-04 | NA |
| Lung5 | 3 | -0.024 | -0.018 | -0.012 | -0.027 | -0.024 | 0.091 |
| Lung6 | 1 | -0.009 | -0.012 | 0.019 | -0.007 | -0.014 | NA |
| Lung9 | 2 | -0.038 | 6.40e-04 | 0.024 | -0.034 | 0.011 | 0.601 |

Per-section counts, detection, and both radii are in `results/spatial_trop2_cd8_cldn4/tables/cosmx_section.tsv`.

## Visium: mediation scored; same-spot is not the lead

Matrices read in this run: GSE189487 (6 early LUAD sections; histology from GEO sample characteristics), GSE273378 (16 stage I LUAD sections), and two 10x CytAssist FFPE demos (LUSC; lung neuroendocrine). Each is whole transcriptome, so TACSTD2, CLDN4, and CD8A are present. Spots are in-tissue with ≥200 counts, log1p(CP10k). Epithelial spots are the top 40% of log KRT8 + log EPCAM inside the section. Local CD8 is the mean log CD8A of other spots inside 1.5× the median nearest-neighbor distance (hex ring 1). E-MTAB-13530 was not re-downloaded: `ftp.ebi.ac.uk` TLS failed in this environment, so no mediation number is reported for it here.

| Series | quantity | n | median | n negative | Wilcoxon p |
| --- | --- | ---: | ---: | ---: | ---: |
| GSE189487 | `local_epi_rho` | 6 | 0.009 | 2/6 | 0.4375 |
| GSE189487 | `local_epi_partial_rho` | 6 | 0.006 | 2/6 | 0.6875 |
| GSE189487 | `local_epi_c_prime` | 6 | -0.002 | 3/6 | 0.6875 |
| GSE189487 | `local_epi_prop` | 2 | 0.025 | 1/2 | NA |
| GSE189487 | `same_epi_partial_rho` | 6 | 0.015 | 2/6 | 1.0000 |
| GSE189487 | `same_epi_prop` | 4 | 0.086 | 1/4 | NA |
| GSE189487 | `same_all_rho` | 6 | 0.002 | 2/6 | 1.0000 |
| GSE189487 | `same_all_partial_rho` | 6 | -0.001 | 3/6 | 0.6875 |
| GSE189487 | `same_all_prop` | 2 | 0.219 | 0/2 | NA |
| GSE273378 | `local_epi_rho` | 16 | -0.007 | 8/16 | 0.9799 |
| GSE273378 | `local_epi_partial_rho` | 16 | -0.007 | 8/16 | 0.9399 |
| GSE273378 | `local_epi_c_prime` | 16 | -0.005 | 8/16 | 0.7436 |
| GSE273378 | `local_epi_prop` | 11 | 0.062 | 4/11 | 1.0000 |
| GSE273378 | `same_epi_partial_rho` | 16 | -0.010 | 10/16 | 0.7436 |
| GSE273378 | `same_epi_prop` | 13 | -0.007 | 7/13 | 0.7354 |
| GSE273378 | `same_all_rho` | 16 | 0.003 | 7/16 | 0.4954 |
| GSE273378 | `same_all_partial_rho` | 16 | 3.49e-04 | 8/16 | 0.9799 |
| GSE273378 | `same_all_prop` | 11 | 0.014 | 5/11 | 0.5771 |

Ring-1 partial ρ is centered on zero (GSE189487 median 0.006, 2/6 negative, p = 0.6875; GSE273378 median -0.007, 8/16 negative, p = 0.9399). Same-spot partial ρ, which is scored only because the spot matrix exists, is also centered on zero. Median neighbor degree on epithelial spots is 6, consistent with one Visium hex ring. Proportions are again undefined when |c| < 0.02, so the proportion rows have a smaller n than the correlation rows.

Column key: `local_epi_partial_rho` is the ring-1 test after CLDN4; `local_epi_prop` is the mediation proportion on that local outcome; `same_all_partial_rho` and `same_epi_partial_rho` are same-spot scores. Same-spot rows are composition-confounded and are not the lead.

10x demos are one section each and are not entered in the Wilcoxon table:

| Section | histology | local partial ρ | local c′ | local proportion | same-spot all-spots partial ρ | same-spot proportion |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 10x_CytAssist_FFPE_LUSC | LUSC | -0.011 | -0.012 | NA | 0.023 | NA |
| 10x_CytAssist_11mm_NEC | lung_neuroendocrine | 0.047 | 0.202 | -0.004 | 0.040 | 0.099 |

## Xenium: no public lung panel here has TACSTD2 + CLDN4 + CD8

CLDN4 is the mediator. Where it is absent, the residual was not computed and no cell-level CD8 distance was invented.

| Dataset | features | TACSTD2 | CLDN4 | CD8A | why not scored |
| --- | ---: | ---: | ---: | ---: | --- |
| 10x Xenium lung preview + add-on (cancer) | 390 | 1 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| 10x Xenium lung preview + add-on (non-diseased) | 390 | 1 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| 10x Xenium hLung cancer multi-tissue preview | 371 | 0 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| 10x Xenium FFPE LUAD multimodal | 371 | 0 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| 10x Xenium FFPE NSCLC IO + add-on | 479 | 0 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| 10x Xenium v1 lung panel post-Xenium LUAD | 289 | 1 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| 10x Xenium Prime 5K post-Xenium LUAD | 5001 | 0 | 0 | 1 | CLDN4 absent; mediator cannot be residualized |
| GSE300007 LUAD TMA2 unimodal/multimodal | 541 | 1 | 0 | 1 | CLDN4 absent on the downloaded features.tsv; LUAD TMA not scored |
| GSE319755 tumor core | 568 | 1 | 0 | 1 | CLDN4 absent; cell coordinates not used |
| GSE319755 adjacent lung | 568 | 1 | 0 | 1 | CLDN4 absent; cell coordinates not used |
| GSE319755 invasive front | 568 | 1 | 0 | 1 | CLDN4 absent; cell coordinates not used |
| GSE311609 NSCLC Prime 5K L1 | 10029 | 0 | 0 | 1 | TACSTD2 and CLDN4 both absent |
| GSE343063 SCLC p2 | 10029 | 0 | 0 | 1 | TACSTD2 and CLDN4 both absent |

GSE311609 and GSE343063 each store 5,001 Gene Expression features plus codeword and control rows (10,029 feature rows in the h5). TACSTD2 and CLDN4 are absent from those gene names. GSE319755 (480 genes plus controls and 27 proteins) and GSE300007 LUAD Xenium have TACSTD2 and CD8A and do not have CLDN4. Official 10x lung / NSCLC `gene_panel.json` files downloaded in this run likewise lack CLDN4; three of them include TACSTD2 and CD8A.

## How to rerun

Raw matrices are not committed. `scripts/download_spatial_trop2_inputs.sh` fetches the Zenodo CosMx zip, the two 10x Visium demos, GSE189487, and GSE273378. Then:

```bash
python3 scripts/spatial_trop2_cd8_after_cldn4.py
```

