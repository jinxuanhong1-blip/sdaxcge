# TACSTD2, CLDN4, and immune context

Observational evidence matrix for the claim that a TROP2 (TACSTD2)–immune association depends partly on CLDN4. Platforms: He 2022 CosMx NSCLC and the locked concordant-4 scRNA cohorts.

The one-page matrix is `FINDING.md`. Definitions are in `METHODS.md`.

```bash
python3 methods/trop2_cldn4_mediation/scripts/concordant4_mediation.py
python3 methods/trop2_cldn4_mediation/scripts/write_matrix.py
```

CosMx reads the figshare object (2,755,776,882 bytes):

```bash
curl -L -o /tmp/cosmx_dl/cosmx_human_nsclc_clustered.h5ad \
  https://ndownloader.figshare.com/files/46841842
python3 methods/trop2_cldn4_mediation/scripts/cosmx_mediation.py
python3 methods/trop2_cldn4_mediation/scripts/write_matrix.py
```

The concordant-4 script stops if CLDN4 percent detected versus T/NK fraction is not the locked DerSimonian–Laird ρ of −0.5311678045689989.
