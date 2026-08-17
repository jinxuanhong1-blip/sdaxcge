# GSE205335 malignant REAL Slingshot/PAGA, CLDN4 only

ADDITIVE. Ahn / Lee ICI (RECIST). Author-malignant cells only.
Root is not a CLDN4-high Leiden cluster. Patient is the unit.
CLDN4 + barrier (CLDN4 excluded) + IFN along pseudotime, and DPT vs RECIST.
No TACSTD2∩CLDN4 dual-high gate. A real lineage plot is required.

Done when `results/tables/lineage_vertices.tsv` and `results/tables/patient_means.tsv` exist.

```bash
pip install -r methods/gse205335_slingshot_real_cldn4/requirements.txt
Rscript -e 'dir.create(Sys.getenv("R_LIBS_USER", "~/R/library"), recursive=TRUE);
  .libPaths(Sys.getenv("R_LIBS_USER", "~/R/library"));
  install.packages("BiocManager", repos="https://cloud.r-project.org");
  BiocManager::install(c("slingshot","SingleCellExperiment"), ask=FALSE, update=FALSE)'
python3 methods/gse205335_slingshot_real_cldn4/scripts/download.py --out /tmp/gse205335_slingshot_real
python3 methods/gse205335_slingshot_real_cldn4/scripts/extract_malignant.py \
  --data /tmp/gse205335_slingshot_real \
  --out /tmp/gse205335_slingshot_real/malignant.h5ad
python3 methods/gse205335_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse205335_slingshot_real/malignant.h5ad \
  --outdir methods/gse205335_slingshot_real_cldn4/results \
  --finding methods/gse205335_slingshot_real_cldn4/FINDING.md
```
