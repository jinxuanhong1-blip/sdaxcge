# methods/grn_ccc — GRN & cell-cell communication playbook

A methods package for one biological question on 10x lung ICI (immune-checkpoint
inhibitor) scRNA-seq:

> **Does TACSTD2/CLDN4-high epithelium show a *reduced T-cell recruitment
> signal* (lower CXCL9/10/11, CCL5, CXCL16, …) compared with TACSTD2/CLDN4-low
> epithelium?**

It covers **LIANA / CellChat** (ligand–receptor), **NicheNet** (ligand→target
activity), and **SCENIC / pySCENIC** (TF regulons): how to run them, how they
map onto this question, **where each one overclaims**, and a patient-level
guardrail that is hard to overclaim.

## Contents

| Path | What |
|------|------|
| [`playbook.md`](playbook.md) | The playbook — bilingual (English + 中文). Start here. |
| [`config/gene_sets.yaml`](config/gene_sets.yaml) | Single source of truth for gene panels & parameters. |
| [`config/public_datasets.yaml`](config/public_datasets.yaml) | Public lung ICI scRNA catalog (sizes, ICI status, <2 GB flag). |
| [`scripts/00_preprocess.py`](scripts/00_preprocess.py) | QC, ambient/doublet removal, malignant calling, define `epi_state`. |
| [`scripts/01_liana_lr.py`](scripts/01_liana_lr.py) | LIANA consensus LR (Python), incl. patient-level test. |
| [`scripts/02_cellchat.R`](scripts/02_cellchat.R) | CellChat pathway communication (R). |
| [`scripts/03_nichenet.R`](scripts/03_nichenet.R) | NicheNet ligand-activity (R) — for T-cell *state*, not recruitment. |
| [`scripts/04_pyscenic.sh`](scripts/04_pyscenic.sh) + [`04b_scenic_downstream.py`](scripts/04b_scenic_downstream.py) | SCENIC regulon inference & analysis. |
| [`scripts/05_pseudobulk_guardrail.py`](scripts/05_pseudobulk_guardrail.py) | **Primary** patient-level chemokine test. |
| [`scripts/06_public_lung_loader.py`](scripts/06_public_lung_loader.py) | Load GSE207422 (default) / GSE205335 into the playbook schema. |
| [`scripts/07_scenic_ccc_joint.py`](scripts/07_scenic_ccc_joint.py) | Conservative join of 05 + 01 + 04b; refuses to overclaim. |
| [`demo/`](demo/) | Synthetic LIANA mechanics check + **GSE207422** real public demo (<2 GB). Computed: 8/10 ligands lower in high state, 0/10 FDR<0.05 — do not claim reduced recruitment. |
| [`env/`](env/) | `requirements.txt` (Python) and `environment.yml` (Python + R). |

Scripts are **notebook-style** (`# %%` cells): run top-to-bottom as `.py`/`.R`,
or open the cells in Jupyter / VS Code. They are **templates** — read the `NOTE`
/ `EDIT` comments and adapt paths and labels to your object.

## Quick start

```bash
pip install -r env/requirements.txt      # Python tools (LIANA validated)
python demo/run_liana_demo.py            # mechanics check (synthetic, no download)
python demo/run_gse207422_demo.py        # real public ICI object (downloads 176 MB)
# real analysis on your own object:
python scripts/00_preprocess.py          # -> data/lung_ici_preprocessed.h5ad
# or: python scripts/06_public_lung_loader.py
python scripts/01_liana_lr.py
python scripts/05_pseudobulk_guardrail.py
python scripts/07_scenic_ccc_joint.py
# CellChat / NicheNet / pySCENIC need extra resources (see env/environment.yml)
```

## The one-line methodology

Use LIANA/CellChat/NicheNet/SCENIC to **generate ranked hypotheses**, then let a
**patient-level pseudobulk test** and **orthogonal data** (spatial, bulk/TCGA,
IHC) decide. Never report a cell-level p-value as evidence for a between-patient
claim. See [`playbook.md`](playbook.md) for the full reasoning.
