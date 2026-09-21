# MEBOCOST — concordant-4 CLDN4-high malignant → T/NK

Additive metabolite communication layer on the locked concordant-four
scRNA-seq cohorts (GSE123902, GSE131907, GSE205335, GSE189357). Protein
ligand–receptor results stay in the CellChat run. This folder does not
merge private 8KL data and does not add GSE148071.

CINE was not an installable metabolite-CCC tool. MEBOCOST
(https://github.com/kaifuchenlab/MEBOCOST) ran. See `METHODS.md` and
`FINDING.md`.

```bash
bash methods/mebocost_concordant4_cldn4/scripts/install_mebocost.sh
bash methods/mebocost_concordant4_cldn4/scripts/download.sh /tmp/concordant4_raw
# GSE205335 is RDS-only. Peel the double gzip, then:
gzip -dc /tmp/concordant4_raw/GSE205335/GSE205335_Lung_IO_UMI_matrix.rds.gz | gzip -dc > /tmp/gse205335.rds
# gene list is written by the Python loader's database scan; the committed
# subset script expects one symbol per line.
Rscript methods/mebocost_concordant4_cldn4/scripts/subset_gse205335.R \
  /tmp/gse205335.rds /tmp/mebo_wanted_genes.txt /tmp/gse205335_subset
python3 methods/mebocost_concordant4_cldn4/scripts/run_mebocost.py \
  --raw=/tmp/concordant4_raw
```

`/tmp/mebo_wanted_genes.txt` is the enzyme + sensor + lineage gene list.
The Python script rebuilds that list from the cloned MEBOCOST database at
runtime for the other three cohorts.
