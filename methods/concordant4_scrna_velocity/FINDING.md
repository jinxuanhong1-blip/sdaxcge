# Finding — concordant-4 RNA velocity is not estimable

CLDN4-only. Cohorts are the locked four: GSE123902 (13), GSE131907 (21), GSE205335 (22), GSE189357 (9). Catalog n = 65. No added cohort.

**Question.** Do epithelial cells flow toward a CLDN4-high barrier state, in the RNA-velocity sense?

**Answer.** Not estimable. scVelo was not run. None of the four accessions deposits a spliced layer and an unspliced layer. A total-count matrix cannot be split back into those layers.

## Gate

| Accession | Public matrix | Spliced / unspliced file | Public reads |
| --- | --- | --- | --- |
| GSE123902 | SEQC dense gene table (one matrix) | none | SRA PRJNA510251, 17 experiments, 122 public sralite runs, 163261 MB, paired 124 bp. No BAM path |
| GSE131907 | raw UMI table, header is cell barcodes (one matrix); log2TPM companion | none | BioProject PRJNA545296 has 0 SRA runs. Sample supplements are NONE |
| GSE205335 | UMI RDS (R XDR header: byte `X` then newline, after two gzip wraps) plus a cell-identity table | none | BioProject PRJNA844398 has 0 SRA runs. Authors point at EGA EGAD00001008703. Not opened |
| GSE189357 | Cell Ranger 4.0.0 matrix.mtx, one MatrixMarket matrix, feature type `Gene Expression` | none | SRA PRJNA782639, 9 public sralite runs, 163381 MB, paired 187 bp. No BAM path |

GEO supplementary names containing loom, h5ad, spliced, unspliced, velocyto, or intron: **none** (`results/tables/geo_suppl_files.tsv`).

Header checks are prefixes, not full downloads (`results/tables/matrix_header_probe.tsv`). GSE189357 `matrix.mtx` begins `%%MatrixMarket matrix coordinate integer general` and then one dimension line (`33538 15216 19291970` on TD1). That is one count matrix, not spliced plus unspliced plus ambiguous.

## Why the FASTQ was not turned into velocity

GSE123902 and GSE189357 public files are sralite reads, not genome-aligned BAMs and not loom files. velocyto / STARsolo / kallisto-bustools would be a new alignment. The other two cohorts have no public reads, so that alignment would drop 21 + 22 of the 65 units and would no longer be concordant-4.

## Substitutes that were not used

Pseudotime on a neighbor graph, and expression-only vector fields, are different estimators. They are not reported here as RNA velocity. The patient-level CLDN4 versus T/NK result on these 65 units is left as it stands.

## What this does not say

- It does not say cells move toward CLDN4-high, and it does not say they do not.
- It does not rescore exclusion, T/NK fraction, or IFN/MHC.
- Cell-level arrows were not drawn. There is no velocity embedding in this folder.

## Reproduce

```bash
python3 methods/concordant4_scrna_velocity/inventory.py
```

Machine-readable gate: `results/gate.json` (`scvelo_run: false`, `biological_answer: not_estimable`).
