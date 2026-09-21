# Finding — concordant-4 malignant Palantir + PAGA + Slingshot

**This is not an ICI clock.** Pseudotime orders malignant cell state from an AT2-like pole to a barrier-like pole. It is not immunotherapy exposure, not time-on-treatment, and not MPR/RECIST. GSE205335 is an ICI cohort in the locked four, but RECIST was not a variable in any test here.

ADDITIVE. **CLDN4 only.** Locked concordant-4 is GSE123902 + GSE131907 + GSE205335 + GSE189357. The locked malignant CLDN4 %pos vs T/NK result (n=65, ρ=−0.531, p=1.6×10⁻⁵, I²=0%) is not re-derived and is not re-audited. GSE148071, GSE127465, GSE154826, GSE200563, and E-MTAB-13526 are not added.

Object: **malignant cells only**, same gates as the locked table. GSE131907 uses `Cell_subtype == Malignant cells` (tS1/tS2/tS3 and nLung AT2 are outside that gate and were not added). GSE205335 uses author `Malignant cells` on non-normal tissue. GSE123902 and GSE189357 use marker malignant `(EPCAM|KRT8|KRT18|KRT19)>0` and `PTPRC==0`. Normal GSE123902 libraries were dropped.

## Verdict

Primary test is the patient/sample mean **Palantir probability of the barrier terminal** versus the **locked T/NK fraction**. It is **positive**: n=65, pooled ρ=0.273 (p=0.028); DL ρ=0.343 (p=0.015, I²=10%, k=4). Within cohort: GSE123902 ρ=0.538 (n=13, p=0.058), GSE131907 ρ=0.404 (n=21, p=0.069), GSE189357 ρ=0.600 (n=9, p=0.088), GSE205335 ρ=0.034 (n=22, p=0.879). The ICI cohort is the flat one. This is the opposite sign from the locked malignant CLDN4 %pos vs T/NK association. It does not re-audit that result and it does not replace it. Barrier-fate probability is not a surrogate for CLDN4 %pos.

CLDN4 along the **Slingshot** AT2→barrier lineage, cell bins (descriptive): CLDN4 1.394 → 1.737, barrier score 0.211 → 1.240, AT2 1.646 → -0.635. The rise in CLDN4 is not monotone. Patient-mean CLDN4 vs that lineage's pseudotime is null (n=48, pooled ρ=0.024, p=0.871). n=48 because 17 units have no cells on the barrier lineage, so their mean pseudotime is undefined. Lineage weight vs T/NK on all 65 units is null (ρ=0.077, p=0.543).

**Palantir pseudotime is not that path.** Patient-mean AT2 falls along it (n=65, pooled ρ=-0.644, p=6.97×10⁻⁹), so the clock leaves the AT2 pole. Patient-mean CLDN4 falls (pooled ρ=-0.272, p=0.028). Cell bins: AT2 1.601 → -0.361, but the barrier score peaks in the middle (0.807) and is 0.227 in the last bin; CLDN4 ends at 1.488 after a mid-course drop. The barrier cluster (Leiden 21) does carry high barrier score (1.31), high CLDN4 (2.29), and higher fate probability (mean 0.62). The max-CLDN4 cluster is a different node (Leiden 18, mean CLDN4 2.33) and was excluded from the root by rule. Patient-mean CLDN4 vs barrier-fate probability is negative (pooled ρ=-0.304, p=0.014).

PAGA connects the AT2 root (Leiden 6) directly to the barrier cluster (Leiden 21). Slingshot's MST leaf path is longer: 6 → 13 → 0 → 18 → 20 → 21.

## What was run

- **PAGA** (`scanpy.tl.paga` on Leiden 0.6 after Harmony `batch=dataset`).
- **Palantir** 1.4.5: diffusion maps on Harmony PCs, 400 waypoints, early cell in the AT2-like cluster, two terminals (barrier program; IFN-high alternative). Fate is an absorption probability, not a treatment clock.
- **Slingshot** via `pyslingshot.core getLineages (Street et al. 2018 MST) + lineage-restricted principal curves`. Start cluster = AT2-like Leiden 6. Barrier lineage rule: lineage leaf is the barrier cluster. Lineage: 6 → 13 → 0 → 18 → 20 → 21.

The barrier terminal is the cell nearest the 80th percentile of the barrier/keratin score inside the highest-barrier Leiden cluster. **CLDN4 is not in that score and is not used to pick the cell.** The root cell is inside the highest-AT2 Leiden cluster after excluding the max-CLDN4 cluster, and is at or below that cluster's median CLDN4.

## Honest n

- Cells after QC and cap ≤180/unit: **10824**.
- Units in the object: **65**. Eligible for Spearman (n_cells≥10 and locked T/NK present): **65**.
- By dataset (cells / units / eligible): {'GSE123902': {'cells': 1930, 'units': 13, 'eligible': 13}, 'GSE131907': {'cells': 3516, 'units': 21, 'eligible': 21}, 'GSE189357': {'cells': 1620, 'units': 9, 'eligible': 9}, 'GSE205335': {'cells': 3758, 'units': 22, 'eligible': 22}}.
- PAGA path root → barrier: ['6', '21']. Connected=True.
- Early cell `GSE123902:LX684_69` CLDN4=0.000, AT2=1.706.
- Slingshot lineages: 9. Barrier cluster is a lineage leaf: True.
- SFTPC and SFTPA1 are absent from the GSE123902 public dense matrix, so they are absent from the shared gene set. The AT2 score uses SFTPB, NAPSA, LAMP3, and ABCA3. SFTPC was not imputed.

## Primary and pre-specified tests

Unit = patient (GSE123902 donor, GSE189357, GSE205335) or sample (GSE131907), matching the locked table. Pooled Spearman treats units equally. DL is Fisher-z DerSimonian–Laird across the four cohorts. Cell-level p-values are not the claim.

| Contrast | n | pooled ρ | pooled p | DL ρ | DL p | I² |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| fate_barrier vs T/NK | 65 | 0.273 | 0.028 | 0.343 | 0.015 | 0.10 |
| slingshot barrier PT vs T/NK | 48 | -0.213 | 0.145 | 0.070 | 0.829 | 0.73 |
| slingshot barrier weight vs T/NK | 65 | 0.077 | 0.543 | -0.047 | 0.862 | 0.72 |
| CLDN4 vs Palantir PT | 65 | -0.272 | 0.028 | -0.290 | 0.030 | 0.00 |
| CLDN4 vs Slingshot barrier PT | 48 | 0.024 | 0.871 | -0.190 | 0.287 | 0.14 |
| barrier score vs Palantir PT | 65 | -0.270 | 0.029 | -0.333 | 0.118 | 0.58 |
| AT2 score vs Palantir PT | 65 | -0.644 | 6.97e-09 | -0.735 | 8.18e-12 | 0.00 |
| SFTPC-program AT2 vs fate_barrier | 65 | -0.597 | 1.54e-07 | -0.671 | 2.39e-04 | 0.58 |
| barrier score vs fate_barrier | 65 | -0.185 | 0.139 | -0.157 | 0.249 | 0.00 |
| CLDN4 vs fate_barrier | 65 | -0.304 | 0.014 | -0.354 | 0.007 | 0.00 |

### Within cohort — barrier fate vs locked T/NK

| Dataset | n | ρ | p |
| --- | ---: | ---: | ---: |
| GSE123902 | 13 | 0.538 | 0.058 |
| GSE131907 | 21 | 0.404 | 0.069 |
| GSE205335 | 22 | 0.034 | 0.879 |
| GSE189357 | 9 | 0.600 | 0.088 |

## What this does not claim

- This is not an ICI clock, not a response clock, and not a test of RECIST or MPR.
- The locked CLDN4 %pos vs T/NK ρ=−0.531 was not recomputed. Mean CLDN4 on this capped object is a different score.
- Palantir terminals and Slingshot lineages are phenotypic orderings, not demonstrated developmental lineages.
- Do not write “AT2 differentiates into a barrier tumor because PAGA is connected.”
- Malignant labels are author calls (GSE131907, GSE205335) or the locked marker gate. They are not CNV re-calls.
- No TACSTD2∩CLDN4 dual-high gate. No GSE148071.

## Reproduce

```bash
python3 methods/concordant4_malignant_traj_cldn4/scripts/download.py --out /tmp/c4raw
python3 methods/concordant4_malignant_traj_cldn4/scripts/extract_malignant.py --data /tmp/c4raw --out /tmp/c4_malignant.h5ad
python3 methods/concordant4_malignant_traj_cldn4/scripts/analyze.py --input /tmp/c4_malignant.h5ad
```

