# Extra table: public lung-tumor scRNA, malignant TACSTD2 / CLDN4 before vs after ICI

Zhejiang paired IHC (user A7) is taken as given. This extra is the public single-cell catalog.

Numbers are filled by `scripts/run.sh` into `tables/paper_table_scrna_prepost.tsv` and `tables/combined_stouffer.tsv`.

## Search result

Among public human lung-tumor scRNA series, **one** accession deposits both a pre-ICI/chemo-IO timepoint and a post timepoint with tumor cells: **GSE207422** (3 unmatched pre biopsies, 12 unmatched post resections). GSE337519 claims one paired patient but deposits a single unlabeled 10x library. GSE205335 has no public timepoint labels (multi-sample patients are multi-site, not serial). GSE241934, GSE146100, and GSE291670 are post-only. T-sorted series (GSE179994 and others) have no malignant compartment. Controlled-access paired series (HRA006493) were not used.

## Paper table (sample unit)

See `tables/paper_table_scrna_prepost.tsv` after the run. Primary metric: malignant mean log1p(CP10K), Mann–Whitney U, unmatched.

## Combined estimate

See `tables/combined_stouffer.tsv`. One independent cohort; Stouffer reduces to that study’s signed p.
