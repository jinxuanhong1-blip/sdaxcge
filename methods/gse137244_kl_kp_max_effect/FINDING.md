# GSE137244 KL vs KP, n = 5 vs 5

Deng et al. cell-line RNA-seq (GEO GSE137244). Five KrasG12D;Lkb1 libraries (KL155, KL47, KLC, KLD, KLE) versus five KrasG12D;Trp53 libraries (B6AL10-1 through B6AL10-5). `normal-lung-RNA` (GSM4073826) is normal lung, not a tumor cell line, and is excluded. Scale is log2(FPKM+1). Every primary row uses all ten cell lines.

The locked unadjusted contrasts reproduce on this file: Tacstd2 Δ = +3.238, Cldn4 Δ = +5.570. Both completely separate. Exact two-sided Mann–Whitney p = 2/252 = 0.00794.

## Cldn4 and Tacstd2

Among adjustments that keep n = 5 vs 5, the unadjusted values are the joint maximum of Δ and of Welch |t|.

| endpoint | filter | Δ | Welch t | Welch p | MW p |
|---|---|---:|---:|---:|---:|
| Cldn4 | none | +5.570 | +10.188 | 2.72e-05 | 0.00794 |
| Cldn4 | minus log2(Epcam+1) | +2.941 | +4.706 | 0.00457 | 0.01587 |
| Cldn4 | residual on Epcam | +2.407 | +2.945 | 0.0349 | 0.15079 |
| Tacstd2 | none | +3.238 | +4.648 | 0.00571 | 0.00794 |
| Tacstd2 | minus log2(Epcam+1) | +0.610 | +1.348 | 0.228 | 0.30952 |
| Tacstd2 | residual on Epcam | +0.868 | +2.334 | 0.0532 | 0.05556 |

Epcam itself is higher in KL (Δ = +2.629, Welch t = +2.554, p = 0.0626), so subtracting it or residualizing on it removes part of the KL-versus-KP difference. KP B6AL10-3 has Epcam FPKM 1.35. Any Epcam sample gate above that value drops B6AL10-3 and leaves n = 5 vs 4. The gate that still keeps five and five does not remove any line. Those gates are in `tables/epcam_sample_gate.tsv` and are not the result.

Leave-one library keeps KL > KP for Cldn4, Tacstd2, TJ7, and the Epcam-filtered TJ mean in every drop. The largest |t| is always the drop of B6AL10-3, which is 5 vs 4: Cldn4 Δ = +5.125, t = +15.183; Tacstd2 Δ = +2.584, t = +10.518. That drop raises |t| and lowers Δ. It is not a 5-vs-5 result.

## Tight junction

TJ7 is the mean of Cldn3, Cldn4, Cldn6, Cldn7, Cdh1, F11r, and Ocln. On this matrix Δ = +3.269, Welch t = +9.857, exact MW p = 0.00794. The handoff value TJ +3.03 is a different average and is not recomputed here.

Epcam membership filter, τ = 0.5, on those seven genes: keep a gene if its Pearson r with Epcam across the ten lines is at least 0.5. The rule does not use the KL/KP labels. It removes Cdh1 (r = -0.661, Δ = -0.483) and F11r (r = +0.178, Δ = -0.098). The remaining genes are Cldn3, Cldn4, Cldn6, Cldn7, Ocln.

That five-gene mean (TJ_EPCAM) is Δ = +4.692, Welch t = +10.409, complete separation, exact MW p = 0.00794. On the τ sweep of TJ7, τ = 0.3 selects the same genes, and no other τ in {-1, 0, 0.3, 0.5, 0.7, 0.8} has a larger Δ or a larger |t| while keeping at least three genes. Because the correlations are estimated on these same ten profiles, the gene list is internal to this matrix. Given the list, the label-randomization p for |Welch t| is 2/252 = 0.00794.

Leave-one gene on TJ_EPCAM, with the choice made on these labels:

- Maximum |t|: drop Cldn3. Genes Cldn4, Cldn6, Cldn7, Ocln. Δ = +4.825, Welch t = +13.960. Still n = 5 vs 5 and completely separated. The parametric Welch p (3.98e-05) is not adjusted for the search. The randomization p that re-picks the drop under every 5-vs-5 labeling is 2/252 = 0.00794. The two labelings that reach it are the observed KL set and the full swap. That is the same floor as a two-sided Mann–Whitney test on a fixed score.
- Maximum Δ: drop Ocln. Genes Cldn3, Cldn4, Cldn6, Cldn7. Δ = +5.499, Welch t = +11.159. The upper-tail randomization p for that maximized Δ is 1/252 = 0.00397: only the observed labeling produces a re-selected Δ at least this large. This is one-sided. It is not a two-sided Mann–Whitney p, and it does not replace the 2/252 floor.

Both leave-one rows are above TJ7 and above TJ_EPCAM on the metric they optimize, and each also improves the other metric relative to TJ_EPCAM. They are the two Pareto rows of a six-candidate grid (the filtered mean, plus one drop for each of its five genes). They are not a second external cohort.

## What this does not change

The locked cell-line statement stays Tacstd2 Δ = +3.24, Cldn4 Δ = +5.57, n = 5 vs 5, Mann–Whitney p = 0.00794. Epcam adjustment does not raise those two deltas or their Welch |t|. The TJ increase is a filter on a pre-listed seven-gene mean, not a new experiment. No library was deleted to change n. No value was imputed. Private 8-KL matrices and TISMO LLC were not used.
