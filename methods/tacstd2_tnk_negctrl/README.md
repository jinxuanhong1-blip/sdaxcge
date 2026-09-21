# TACSTD2–T/NK attenuation: CLDN4 versus negative controls

Does the patient-level association between TACSTD2 and T/NK move closer to zero when CLDN4 is partialled out than when CLDN3, CLDN7, EPCAM, MUC1, or KRT19 is partialled out?

Result, in `FINDING.md`: CLDN7 attenuates about as much as CLDN4 where that association is negative (five TCGA cohorts, label-swap p = 0.27). The pre-specified seven-cohort pool does separate CLDN4 from CLDN7, but its unadjusted correlation is weak and mixes opposite signs. Concordant-4 has no TACSTD2–T/NK association to attenuate (ρ = −0.112, p = 0.46).

Two public layers, kept separate:

- Concordant-4 scRNA (GSE123902, GSE131907, GSE205335, GSE189357), n = 65. Malignant percent positive is the primary score. Mean log1p is a sensitivity. The malignant gate does not use TACSTD2 or CLDN3.
- TCGA primary tumors, Xena GDC STAR log2(TPM+1). The primary meta-analysis is the locked seven-cohort set (LUAD, BRCA, CESC, KIRC, STAD, BLCA, PAAD). The outcome is the mean of CD3D, CD3E, CD3G, CD8A, NKG7, GNLY, and KLRD1.

Attenuation is partial Spearman ρ minus the unadjusted Spearman ρ. A positive value means a negative association shrank. Percent attenuated is `100 * (ρ_unadj − ρ_partial) / ρ_unadj`. Head-to-head p-values swap the CLDN4 and control labels within cohort (10,000 draws).

```bash
bash methods/tacstd2_tnk_negctrl/download_concordant4.sh /tmp/geo_c4
python3 methods/tacstd2_tnk_negctrl/analyze_concordant4.py
python3 methods/tacstd2_tnk_negctrl/analyze_tcga.py
python3 methods/tacstd2_tnk_negctrl/analyze_tcga_sign_sensitivity.py
python3 methods/tacstd2_tnk_negctrl/plot_summary.py
python3 methods/tacstd2_tnk_negctrl/test_stats.py
```

The sign-sensitivity label-swap (about four minutes) is `analyze_tcga_sign_sensitivity.py --perm`. Without `--perm` it only rebuilds the meta table. Numbers are in `FINDING.md`.
