# Concordant-4 barrier ligands: TACSTD2 partial dependence on CLDN4

Public scRNA only (GSE123902, GSE131907, GSE205335, GSE189357).
Malignant senders, T/NK receivers. Barrier ligands F11R, NECTIN2, CDH1, LGALS9.

Question: is barrier-ligand outgoing from TACSTD2-high malignant cells still there
after the ligand is residualized on CLDN4, or after TACSTD2 Q4 and Q1 cells are
matched inside CLDN4 strata?

Rules and the shrink call are in `METHODS.md`. Numbers are in `FINDING.md` after
`scripts/run_partial.py`.

```bash
bash methods/concordant4_barrier_partial_tacstd2/scripts/download.sh /tmp/concordant4_raw
python3 methods/concordant4_barrier_partial_tacstd2/scripts/prepare_131907.py
gcc -O3 -o /tmp/extract_131907 methods/concordant4_barrier_partial_tacstd2/scripts/extract_131907.c -lz
/tmp/extract_131907 \
  /tmp/concordant4_raw/GSE131907/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz \
  /tmp/concordant4_cache/gse131907/genes_request.txt \
  /tmp/concordant4_cache/gse131907/keep_idx.i32 \
  208506 \
  /tmp/concordant4_cache/gse131907
# GSE205335 is a dgCMatrix RDS. Rscript export writes a gene-subset Matrix Market.
Rscript methods/concordant4_barrier_partial_tacstd2/scripts/export_gse205335.R \
  /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz \
  /tmp/concordant4_cache/gse131907/genes_request.txt \
  /tmp/concordant4_cache/gse205335
python3 methods/concordant4_barrier_partial_tacstd2/scripts/run_partial.py
```
