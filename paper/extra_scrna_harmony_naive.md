# Extra atlas n · Harmony joint naive LUAD/NSCLC scRNA

Public treatment-naïve / diagnostic tumor scRNA only. Not an ICI test.

Joint object: GSE131907 tumor sites (Kim 2020) + GSE253013 9 LUAD (Xiang/Sze
2024) + GSE148071 (Wu 2021) + GSE127465 tumor (Zilionis 2019), Harmony-corrected
on a shared lineage panel. Leader/GSE154826 omitted (CD45+ CITE-seq; epithelium
is doublet-gate leak).

Primary unit = tumor donor. Eligible: ≥20 epithelial and ≥20 T/NK cells.

Numbers are written by `methods/scrna_harmony_naive/integrate.py` into
`results/scrna_harmony_naive/association_statistics.tsv`.
