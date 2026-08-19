# Panel-level CLDN4 protein vs CD8 protein

Rule: both proteins must be on the **same** multiplex panel with public cell XY.  
`CD8_only` is the common outcome.

| ID | Tissue | Tech | Public cell XY | CD8 protein | CLDN4 protein | Decision |
|---|---|---|---|---|---|---|
| Sorin2023_LUAD_IMC | lung LUAD | IMC 35-plex | yes (Zenodo 7760826) | yes | no (panCK) | reject |
| Sorin2025_LUAD_driver_IMC | lung LUAD | IMC 35-plex | yes | yes | no | reject |
| Cords2025_LUAD_IMC | lung LUAD | IMC | restricted Zenodo 14827075 | ? | ? | skip restricted |
| CosMx_IO64_protein | multi incl. lung | CosMx protein | if deposited | yes | no | reject |
| Xenium_protein_27 | multi | Xenium protein | if deposited | yes (CD8A) | no | reject |
| GeoMx_IO_protein | multi | DSP ROI | ROI not cell | yes | no | reject |
| PhenoCode_IO60 | multi | PhenoCycler | if deposited | yes | no | reject |
| Hickey2023_intestine_CODEX | intestine (not tumor) | CODEX 57 | yes Zenodo 7311360 | yes | **no** (header) | reject |
| Schurch2020_CRC_CODEX | CRC | CODEX 56 | yes | yes | no | reject |
| Jackson2020_breast_IMC | breast | IMC | yes | yes | no | reject |
| Keren2018_MIBI_TNBC | breast | MIBI | yes | yes | no | reject |
| IMMUcan_panel1_IMC | mixed incl NSCLC/CRC | IMC | example yes | yes (CD8a) | no (panel.csv) | reject |
| OMAP2_intestine | intestine | CODEX | HuBMAP | yes | no | reject |
| OMAP34_esophagus | esophagus | CODEX 31 | images, not dual-protein table | yes | no | reject |
| OMAP_pancreas_12plex | pancreas | CODEX 12 | images | unknown | not listed | reject |
| Karlsen2024_breast_CAF_IMC | breast | IMC 42 | paper | yes | no | reject |
| Scheuermann2024_MACSima | HCC/solid | cyclic IF 118 | code only Zenodo 10057717 | yes | no (Table 1) | reject |
| Ferrara2024_PDAC | pancreas | Visium + IHC + CODEX 51 | Visium + CODEX source | CODEX immune | IHC/RNA only | reject (not same panel) |
| Bolen_CosMx_colon_protein | colon | CosMx protein | Zenodo 14851272 | kit CD8 | kit no CLDN4 | reject |
| GSE271689 | NSCLC | spatial multi-omics seq | GEO | not protein cell table | no | reject |
| GEO CLDN4∩IMC/CODEX (9 IDs) | mixed | RNA/ATAC/blood CyTOF | GEO | n/a | n/a | reject false positives |

**Qualified for metrics: 0.**
