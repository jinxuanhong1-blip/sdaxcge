# KL > KP Tacstd2 (public; no private 8KL)

Reproduced from the public matrices listed in `tables/kl_kp/`. Test: two-sided exact Mann–Whitney; delta = mean(KL) − mean(KP). Scale matches each deposit (see `focus_contrasts.tsv`).

Locked GSE137244 cell lines: Tacstd2 **+3.238**, Cldn4 **+5.570**, TJ7 **+3.269**, all p = 0.00794.

| accession | setting | n KL vs KP | Tacstd2 Δ | MW p | role |
|---|---|---|---:|---:|---|
| GSE164758 | primary_bulk_tumor | 9 vs 8 | +0.858 | 8.227e-05 | new |
| GSE6135 | mouse_all_histology | 7 vs 5 | +1.795 | 0.002525 | new_tacstd2_same_cohort_as_prior_cldn4 |
| GSE137244 | cultured_cells | 5 vs 5 | +3.238 | 0.007937 | locked_reference |
| GSE6135 | mouse_adeno_only | 4 vs 5 | +0.989 | 0.01587 | sensitivity_adeno_only |
| GSE137396 | in_vivo_nodule | 5 vs 5 | +1.144 | 0.09524 | rescore |
| GSE244452 | syngeneic_bulk_tumor | 3 vs 3 | +7.228 | 0.1 | new |
| GSE274352 | cultured_cells_empty_vector | 3 vs 3 | -1.543 | 0.1 | new |
| GSE274351 | lcm_early_adenoma | 5 vs 5 | -0.860 | 0.5476 | rescore |

Full rows including Cldn4/TJ7 and inventory no-gos: `tables/kl_kp/contrasts.tsv`.
