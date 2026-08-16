# w200 / A5 — GSE76628 is not KL lung (rejected)

**Do not use GSE76628 for a TROP2-high immune-resistant KL-lung subset.**

GSE76628 (Uhlik et al., Cancer Res 2016, PMID 27197264) is a bulk Affymetrix
Mouse 430 2.0 series of 78 **whole flank-skin** samples from the Ad-VEGF-A164
tumor-surrogate angiogenesis model in **athymic nude mice**, deposited as part
of a gastric-cancer stroma signature paper. It is:

- not lung
- not Kras;Lkb1 (KL)
- not a tumor (VEGF-driven vascular/stromal lesion in normal skin)
- not cell-resolved
- T-cell-deficient (nude)

The original A5 request ("TROP2-high immune-resistant subset") cannot be
answered here. Bulk skin Trop2–immune anti-correlations from an earlier pass
are **tissue-fraction artifacts** and were removed from this folder so they
cannot be cited as KL-lung evidence.

## Where the real analysis is

Public KL lung scRNA-seq that *can* test Tacstd2-high epithelium vs CD8/NK:

**`results/w200/A5_GSE165641/`** — GSE165641, two *KrasLSL-G12D/+; Lkb1fl/fl*
mice, 7,180 cells from lung tumor sections (PMID 34369094).

The original public-only hunt that documented GSE76628 itself remains in
`results/hunt_gse76628/` (PR #41 first commit).
