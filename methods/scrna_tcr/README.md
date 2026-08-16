# methods/scrna_tcr — public lung ICI scRNA+TCR

Additive high-end scRNA+TCR on **public** files only. Skip EGA.

| Path | What |
| --- | --- |
| [`playbook.md`](playbook.md) | Estimands, files, clone/CXCL13 rules, what is not estimable |
| [`scripts/`](scripts/) | Runnable patient-level analysis |
| [`results/`](results/) | Computed tables and figures (honest n/p) |
| [`data/gse243013_residual_tacstd2_cldn4.tsv`](data/gse243013_residual_tacstd2_cldn4.tsv) | Residual (immune-compartment) TACSTD2/CLDN4 from the public MTX |

Primary cohort: **GSE243013** (TCR + MPR, immune-only). Caushi **GSE176022** GEO processed exists (bulk culture TCR); **GSE176021** VDJ+RDS used; EGA FASTQ not downloaded. Leftover open TCR: **GSE179994**. **GSE185204** listed only (n=3).
