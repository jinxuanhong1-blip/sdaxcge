The Stage-03 `usable_maybe` flag is a text+gene-peek heuristic. Three leftovers
tripped it and were then rejected in Stage 04 / by manual GEO metadata:

- GSE217451 — H1650 si-hMENA cell-line FPKM. Both genes present. Not patients.
- GSE224216 — H2030 si-hMENA cell-line TPM. Both genes present. Not patients.
- GSE224099 — melanoma tumour-infiltrating CD4 T-cell subsets. Not lung
  epithelium and not an ICI-response label.

GSE248378 did *not* trip `usable_maybe` (GEO has no outcome field) but is the
only leftover that can be tested, after joining public Nature Communications
Source Data Fig. 5d via unique ITGAE FPKM values.
