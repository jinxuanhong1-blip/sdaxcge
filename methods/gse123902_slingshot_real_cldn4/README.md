# GSE123902 real Slingshot/PAGA — CLDN4 only

ADDITIVE. Laughney et al., *Nat Med* 2020 (PMID 32042191) tumor epithelium.
**CLDN4 only.** Donor is the unit. Root is never CLDN4-high.
No TACSTD2∩CLDN4 dual-high gate. Author 36.5 GB H5 is skipped.

Done when `results/tables/lineage.tsv` exists.

```bash
bash methods/gse123902_slingshot_real_cldn4/scripts/install_tools.sh
python3 methods/gse123902_slingshot_real_cldn4/scripts/download.py --out /tmp/gse123902_slingshot
python3 methods/gse123902_slingshot_real_cldn4/scripts/extract_epithelium.py \
  --data /tmp/gse123902_slingshot \
  --out /tmp/gse123902_slingshot/epithelium.h5ad
python3 methods/gse123902_slingshot_real_cldn4/scripts/analyze.py \
  --input /tmp/gse123902_slingshot/epithelium.h5ad \
  --outdir methods/gse123902_slingshot_real_cldn4/results \
  --finding methods/gse123902_slingshot_real_cldn4/FINDING.md
```
