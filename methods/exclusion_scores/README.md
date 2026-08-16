# Exclusion scores

Reproducible research methods for NSCLC immune-checkpoint inhibitor analyses:

- canonical TIDE bridge;
- canonical public IPS implementation;
- sourced IFNG, CTL, CAF, EMT, Bindea and public TIS gene summaries;
- exact Thompson NSCLC EMT/inflammation and GSVA/ssGSEA stromal models;
- TACSTD2-versus-exclusion permutation tests;
- explicitly reconstructed spatial CD8 exclusion scoring;
- Bessede-style survival and treatment-interaction models.

Start with the bilingual [`playbook.md`](playbook.md). Provenance and method
boundaries are in [`SOURCES.md`](SOURCES.md). All repository additions for this
method are contained in this directory.

Quick check:

```bash
python3 -m unittest discover -s methods/exclusion_scores/tests -v
```
