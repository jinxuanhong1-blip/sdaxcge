# scTenifoldKnk CLDN4 virtual knockout (concordant-4 malignant cells)

Additive analysis. The locked patient-level result is unchanged: malignant CLDN4 percent-positive versus T/NK fraction, DerSimonian–Laird ρ = −0.531 on 65 units (PR for the Seurat/Harmony integration). This folder asks a different question: if CLDN4’s outgoing edges are removed from a malignant-cell gene regulatory network, which genes move, and do those genes overlap the observed CLDN4-low versus CLDN4-high contrast.

Cohorts, and only these: GSE123902, GSE131907, GSE205335, GSE189357. Not GSE148071 or GSE127465.

## What the knockout returns

`scTenifoldKnk` 1.1 builds principal-component-regression networks, denoises them by CP tensor decomposition, zeros the knockout gene’s outgoing row, and ranks genes by manifold distance. Columns `distance`, `Z`, `FC`, and `p.value` describe that distance. They are not an expression log fold-change. The signed quantity used here is the denoised outgoing weight `WT[CLDN4, target]`.

## Reproduce

```bash
bash methods/sctenifoldknk_concordant4_cldn4/scripts/download.sh /tmp/geo_c4
python3 methods/sctenifoldknk_concordant4_cldn4/scripts/extract_malignant.py --dataset GSE123902
python3 methods/sctenifoldknk_concordant4_cldn4/scripts/extract_malignant.py --dataset GSE189357
python3 methods/sctenifoldknk_concordant4_cldn4/scripts/extract_malignant.py --dataset GSE131907
Rscript methods/sctenifoldknk_concordant4_cldn4/scripts/extract_gse205335.R /tmp/geo_c4 /tmp/c4_work/GSE205335
python3 methods/sctenifoldknk_concordant4_cldn4/scripts/observed_contrast.py
for ds in GSE123902 GSE131907 GSE205335 GSE189357; do
  python3 methods/sctenifoldknk_concordant4_cldn4/scripts/prepare_grn.py --dataset $ds
  Rscript methods/sctenifoldknk_concordant4_cldn4/scripts/run_knk.R \
    /tmp/c4_work/$ds/knk_input /tmp/c4_work/$ds/knk_CLDN4 CLDN4
done
# control gene is chosen in compare_ko_observed.py; run the same Rscript with that symbol on GSE131907
python3 methods/sctenifoldknk_concordant4_cldn4/scripts/compare_ko_observed.py
```

Raw matrices are not committed. `observed_contrast.py` stops if malignant CLDN4 detection does not match the locked unit table.
