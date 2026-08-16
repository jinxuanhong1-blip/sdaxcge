# Verification report (fable_stk11)

**28 passed / 0 failed**

| Check | Status | Detail |
|---|---|---|
| TCGA-LUAD genotype rows≈566 | PASS | got 566 |
| TCGA-LUAD no duplicate samples | PASS | 0 dups |
| TCGA-LUAD expression rows≈510 | PASS | got 510 |
| TCGA-LUAD TACSTD2 & CLDN4 present | PASS |  |
| TCGA-LUSC genotype rows≈484 | PASS | got 484 |
| TCGA-LUSC no duplicate samples | PASS | 0 dups |
| TCGA-LUSC expression rows≈481 | PASS | got 484 |
| TCGA-LUSC TACSTD2 & CLDN4 present | PASS |  |
| Rizvi2015 genotype rows≈35 | PASS | got 35 |
| Rizvi2015 no duplicate samples | PASS | 0 dups |
| Hellmann2018 genotype rows≈75 | PASS | got 75 |
| Hellmann2018 no duplicate samples | PASS | 0 dups |
| Rizvi2018 genotype rows≈240 | PASS | got 240 |
| Rizvi2018 no duplicate samples | PASS | 0 dups |
| LUAD STK11 freq in [0.08,0.22] | PASS | 0.133 (n=566) |
| LUAD KEAP1 freq in [0.1,0.24] | PASS | 0.180 (n=566) |
| LUAD KRAS freq in [0.25,0.38] | PASS | 0.297 (n=566) |
| LUAD TP53 freq in [0.4,0.6] | PASS | 0.521 (n=566) |
| recompute STK11->TACSTD2 p matches table | PASS | recomputed=3.387e-05 table=3.387e-05 |
| STK11->TACSTD2 direction is lower-in-mutant | PASS | med_mut=12.202 med_wt=12.715 |
| recompute STK11->CLDN4 p matches table | PASS | recomputed=7.725e-06 table=7.725e-06 |
| STK11->CLDN4 direction is lower-in-mutant | PASS | med_mut=12.746 med_wt=13.211 |
| KL co-mutation ->TACSTD2 reproducible | PASS | 9.973e-05 vs 9.973e-05 |
| KL co-mutation ->CLDN4 reproducible | PASS | 5.732e-05 vs 5.732e-05 |
| CD8A Cox HR<1 (protective direction) | PASS | HR=0.775 |
| CD274 Cox HR<1 (protective direction) | PASS | HR=0.686 |
| KRAS-restricted STK11 pooled OR<1 (worse ICI benefit) | PASS | MH-OR=0.359 |
| processed+raw data < 2 GB | PASS | 2.4 MB |
