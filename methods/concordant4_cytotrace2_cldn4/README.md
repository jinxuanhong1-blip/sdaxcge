# Concordant-4 malignant CLDN4 vs CytoTRACE2 and stemness

ADDITIVE. **CLDN4 only.** Public sets only:

GSE123902 + GSE131907 + GSE205335 + GSE189357.

Not GSE148071, GSE127465, GSE154826, GSE200563, or E-MTAB-13526.
No TACSTD2∩CLDN4 dual-high. Does not re-fit the locked T/NK exclusion ρ.

**Question.** Are CLDN4-high malignant cells more differentiated (lower CytoTRACE2 / lower stemness) and more barrier-like?

Patient / donor / sample is the unit. Primary summary is a DerSimonian–Laird meta of the four within-cohort Spearman correlations. CytoTRACE2 is fit inside each dataset.

```bash
python3 -m venv /tmp/ct2venv
/tmp/ct2venv/bin/pip install -r methods/concordant4_cytotrace2_cldn4/requirements.txt
# CPU torch is enough: pip install torch --index-url https://download.pytorch.org/whl/cpu
bash methods/concordant4_cytotrace2_cldn4/scripts/run_all.sh
```

Done when `FINDING.md` and `results/tables/patient_cldn4_vs_potency.tsv` exist.
Raw matrices stay in `/tmp` (not committed).
