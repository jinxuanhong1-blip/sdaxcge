# Methods

## Question

Do epithelial cells in the concordant-4 pool flow, in the RNA-velocity sense, toward a CLDN4-high barrier state?

## What would have counted as a velocity analysis

scVelo (Bergen et al., Nature Biotechnology 2020), steady-state or dynamical, on cells that already carry spliced and unspliced counts (velocyto / kallisto bus / STARsolo velocity). Labeled new/total counts would also qualify. A single gene-by-cell UMI matrix does not.

The test, had those layers been public:

- Epithelium only. Not a tissue UMAP that mixes immune and stroma.
- Direction is not defined by CLDN4. A CLDN4-high terminal is circular.
- Where an external root exists (GSE123902 uninvolved lung, GSE131907 nLung AT2), that root sets orientation.
- The reported unit would have been the patient or sample (catalog n = 65), not a cell-level p-value.
- Velocity QC (shared spliced and unspliced counts, gene-wise kinetics) can still turn the estimator off. A failed QC is a result.

## What was actually available

`inventory.py` checks three public surfaces.

1. GEO series supplementary file names for loom, h5ad, spliced, unspliced, velocyto, or intron.
2. A short prefix of each deposited matrix: SEQC dense CSV (GSE123902), raw UMI text (GSE131907), Cell Ranger `features.tsv` plus `matrix.mtx` (GSE189357), and the RDS header (GSE205335).
3. SRA experiment counts and runinfo for the linked BioProjects. Aligned BAM in the download path is recorded separately from sralite reads.

GSE205335 raw reads are not in SRA. The BioProject description states EGA accession EGAD00001008703. That controlled archive was not opened.

## What was not substituted

Graph orderings (Slingshot, Palantir, DPT) and expression-only vector fields (scTour) do not use splicing kinetics. They stay out of this folder so a pseudotime trend is not reported as velocity.

Rebuilding a loom from FASTQ is a new quantification. Public sralite reads exist for GSE123902 and GSE189357 only. GSE131907 (PRJNA545296) and GSE205335 (PRJNA844398) have zero public SRA runs. A two-cohort realignment would no longer be concordant-4.

## Locked result that this folder does not touch

Malignant CLDN4-positive fraction versus T/NK on these 65 units is an existing patient-level result. It is not recomputed here and it is not a velocity result.
