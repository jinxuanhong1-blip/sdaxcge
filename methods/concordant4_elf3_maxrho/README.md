# Concordant-4: ELF3 and a 4-gene module vs T/NK

Malignant ELF3, and a four-gene module built on ELF3 + TACSTD2 + CLDN4, against the locked T/NK fraction. Unit = patient / donor / sample, n = 65. The score grid is searched for the largest |DerSimonian–Laird Spearman ρ|. Numbers are computed from the public matrices. The CLDN4 % positive pipeline check has to reproduce ρ = −0.531 before any new score is kept.

Writeup: `FINDING.md`. Methods: `METHODS.md`.

```bash
bash methods/concordant4_elf3_maxrho/scripts/download.sh /tmp/geo_c4
python3 methods/concordant4_elf3_maxrho/scripts/extract_panel.py --geo /tmp/geo_c4
python3 methods/concordant4_elf3_maxrho/scripts/sweep.py --panel /tmp/geo_c4/panel
```
