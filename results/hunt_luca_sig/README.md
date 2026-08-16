# LuCA/Salcher signature hunt: TACSTD2 vs neutrophils/TRN and LCAM

## Bottom line

For the stated `immune-low / TROP2-high` case, LuCA supports **TACSTD2 as an
epithelial/malignant-state signal, not as a neutrophil, TRN, or LCAM signal**.
It does not support inferring either TRN-high or LCAM-high from TACSTD2-high.
The patient's TRN and LCAM states remain unmeasured because no patient
expression matrix was supplied.

The most defensible interpretation is:

1. **TROP2-high is orthogonal to these immune signatures.** TACSTD2 is absent
   from the final Salcher TAN, NAN, TRN, and generic-neutrophil signatures and
   from both sides of the Leader LCAM signature.
2. **LuCA single-cell expression is compartmentalized.** TACSTD2 was detected
   in 55.7% of malignant cells versus 4.14% of neutrophils. Mean normalized
   expression was 15.6-fold higher in malignant cells. In the
   neutrophil-rich BD Rhapsody assay alone, detection was 87.1% in malignant
   cells and 4.85% in neutrophils.
3. **The residual neutrophil calls are not evidence of a neutrophil program.**
   They are low, assay-dependent, and compatible with ambient RNA/doublets.
   No cell-level significance test is reported because cells are not
   independent replicates and the atlas combines platforms and studies.
4. **TROP2 RNA is not tumor-specific or a protein assay.** TACSTD2 was also
   common in normal epithelial populations (74.8% of alveolar type 1 and
   78.2% of club cells). These RNA data cannot establish tumor-cell surface
   TROP2 abundance or ADC eligibility.

## TRN/neutrophil result

The final Salcher paper defines:

- TAN: 18 genes
- NAN: 20 genes
- TRN: their 38-gene union
- A separate 17-gene major-cell-type `Neutrophils` signature in Table S5

The 38-gene TRN and 17-gene generic-neutrophil signatures overlap at only
`CYP4F3` and `MGAM`. They should not be silently substituted for one another.
The TAN and NAN sets are disjoint by construction.

Salcher et al. associated a high TRN signature with anti-PD-L1
(atezolizumab) treatment failure in retrospective OAK/POPLAR bulk RNA data.
That is cohort-level association, not a validated patient-level predictive
test, and it does not imply that this TROP2-high case is TRN-high.

The LuCA source also uses `TACSTD2` with `AGR2` as a marker for an
`undifferentiated` epithelial/cancer annotation. This is an annotation marker
in analysis code, not evidence that TACSTD2 causes or tracks TRNs.

## LCAM result

LCAM is not a single positive gene list. The original Leader implementation:

1. z-scores each gene across samples;
2. averages genes within each cell-subtype signature;
3. averages the three LCAM-hi subtype scores (IgG plasma, SPP1 macrophage,
   activated T);
4. averages the five LCAM-lo subtype scores (B, alveolar macrophage, cDC2,
   AZU1 macrophage, cDC1); and
5. reports `LCAM-hi - LCAM-lo`.

The subtype rows needed to preserve this weighting are included in
`signatures.tsv`. Averaging all 22 hi genes against all 41 lo genes would not
reproduce the published score.

`Immune-low` is not automatically `LCAM-low`: LCAM was designed to capture an
immune state partly independent of total immune content. Leader et al.
excluded the lowest 10% of immune-content samples from downstream bulk
analysis because the signature was less informative there. Thus, in a truly
immune-poor specimen, an LCAM call may be unstable rather than confidently
low.

The Salcher **preprint** stated that TRN and LCAM scores had low correlation
and interpreted them as independent markers. That LCAM analysis and claim are
absent from the final peer-reviewed Cancer Cell article and final supplement,
so it is retained here only as preprint-level evidence; no correlation
coefficient is claimed.

## What was and was not analyzed

The full LuCA objects were deliberately not downloaded:

- CELLxGENE core H5AD: about 12.6 GB
- CELLxGENE extended H5AD: about 17.6 GB
- Zenodo archives: about 5.0–70.8 GB each

All exceed the requested 2 GB ceiling. Instead, this hunt used:

- final supplemental Table S4 (7.8 MB document) for TAN/NAN/TRN;
- final Table S5 (38 KB workbook) for the generic-neutrophil signature;
- the published Leader LCAM scoring source;
- a gene-restricted CELLxGENE Census 2025-11-08 query for TACSTD2 only
  (1,283,972 cells by one gene; the temporary local H5AD was 35.9 MB and is
  not committed).

The standardized Census view has cell type, donor, assay, disease, and tissue,
but lacks the LuCA custom sample field needed to separate primary tumor from
normal-adjacent samples. A donor pseudobulk correlation of tumor-cell TACSTD2
against TRN or LCAM would therefore merge biological sites and be misleading.
It was not performed. A valid patient-specific comparison requires a tumor
bulk matrix or sample-resolved LuCA object plus sample metadata.

## Files

- `signatures.tsv`: exact published signature genes and LCAM subtype structure
- `signature_overlap.tsv`: set intersections
- `tacstd2_relevant_cell_types.tsv`: overall one-gene Census summaries
- `tacstd2_malignant_vs_neutrophil_by_assay.tsv`: assay-stratified sensitivity check
- `reproduce.py`: fixed-version Census query and table regeneration

Reproduction requires Python with `cellxgene-census`, `pandas`, and `numpy`:

```bash
python results/hunt_luca_sig/reproduce.py
```

## Provenance

- Salcher et al., Cancer Cell 2022:
  [doi:10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008)
- Salcher final [Table S4/document](https://ars.els-cdn.com/content/image/1-s2.0-S1535610822004998-mmc1.pdf)
  and [Table S5](https://ars.els-cdn.com/content/image/1-s2.0-S1535610822004998-mmc3.xlsx)
- LuCA analysis source, revision
  [`a3c0184`](https://github.com/icbi-lab/luca/commit/a3c018451850fbb50b674bb606ea7348182420cb)
- LuCA extended atlas:
  [CELLxGENE collection](https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287),
  dataset `1e6a6ef9-7ec9-4c90-bbfb-2ad3c3165fd1`, Census `2025-11-08`
- Leader et al., Cancer Cell 2021:
  [doi:10.1016/j.ccell.2021.10.009](https://doi.org/10.1016/j.ccell.2021.10.009)
- LCAM scoring source, revision
  [`4a88416`](https://github.com/effiken/Leader_et_al/blob/4a884161d50ed768963603c8a0aea38ea4c9299b/scripts/get_LCAM_scores.R)
- Salcher preprint:
  [doi:10.1101/2022.05.09.491204](https://doi.org/10.1101/2022.05.09.491204)
