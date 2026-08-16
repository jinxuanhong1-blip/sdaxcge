# sdaxcge

Spatial-transcriptomics hunts on public processed data.

## E-MTAB-13530 Visium — TROP2/CLDN4 vs T cells

**Question.** In human NSCLC 10x Visium (ArrayExpress
[E-MTAB-13530](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13530)),
are TACSTD2/CLDN4-high spots spatially anti-colocalized with T-cell spots?

**Honest answer.** Yes, in all 20 tumor sections / 8 patients
(neighbourhood-enrichment z median −9.3; patient Wilcoxon p = 0.0039).
No, it is **not specific** to TROP2/CLDN4: EPCAM/KRT-high spots show the
same segregation (patient Wilcoxon on TJ − epi: p = 0.47). This is
tumor-nest vs stroma geography, not a TROP2/CLDN4-specific immune-exclusion
phenotype.

Full write-up, tables, and maps: [`results/hunt_visium_tj/`](results/hunt_visium_tj/).

```bash
bash scripts/download_data.sh          # public processed h5 + spatial.tar
python3 scripts/hunt_visium_tj.py      # writes results/hunt_visium_tj/
```
