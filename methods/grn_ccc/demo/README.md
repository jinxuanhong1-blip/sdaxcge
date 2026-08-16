# Demo — LIANA end-to-end (mechanics + readout)

This demo satisfies the playbook rule *"demo only if a small scRNA <2 GB is
available"* without shipping or downloading a large dataset.

## What it does

`run_liana_demo.py`:

1. **If you provide a small dataset** — set `DEMO_H5AD=/path/to/small.h5ad`
   (with `obs` columns `cell_type`, `epi_state`, `patient_id`) and it runs the
   real LIANA consensus pipeline on it.
2. **Otherwise** — it synthesizes a tiny (~1600 cells, <2 MB, no network)
   labeled AnnData in which TACSTD2/CLDN4-**high** epithelium was *planted* to
   express ~5× fewer T-cell-recruitment chemokines than the **low** state. This
   validates that the LIANA template runs and that the "reduced recruitment
   signal" readout is recovered.

```bash
pip install -r ../env/requirements.txt
python run_liana_demo.py
# or with your own data:
DEMO_H5AD=/data/small_lung.h5ad python run_liana_demo.py
```

## Verified result (synthetic run)

Every recruitment axis (CXCL9/10/11→CXCR3, CCL5→CCR5, CXCL16→CXCR6) is ranked
**weaker** (worse `magnitude_rank`) from the high state than from the low state
— 100% of pairs — matching the planted ground truth. See
`example_output/demo_liana_epi_to_T.csv`.

## Real public demo — GSE207422 (175 MB, <2 GB)

`run_gse207422_demo.py` downloads Hu et al. 2023 (PMID 36869384) processed UMI
matrix + sample metadata, marker-gates epithelium vs T/NK (GEO deposited **no**
cell-type column), splits TACSTD2/CLDN4-high vs -low, then runs the **patient-level
pseudobulk** (primary) and LIANA focus axes (hypothesis). Outputs go to
`gse207422_example/`.

Caveats baked into the script and the playbook §6b: BD Rhapsody not 10x; n=15 /
MPR n=4; ambient not re-estimable; marker gates are not a published annotation.

```bash
python run_gse207422_demo.py
```

## Read this as a teaching artifact, not biology

The synthetic run also pairs the chemokines with receptors like `ADRA2A`,
`GRM7`, `MTNR1A` that are pure background noise in the data. This is exactly the
kind of **false positive** the playbook warns about: LIANA will pair a ligand
with any resource receptor that clears the expression threshold, whether or not
the interaction is biologically plausible. Always filter to the receptors you
care about (here `CXCR3`/`CCR5`/`CXCR6`) and validate — do not read the whole
network as truth.
