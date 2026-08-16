# Sources

## TISMO

- Site: https://tismo.pku-genomics.org/ (formerly http://tismo.cistrome.org/)
- Paper: Zeng et al., *Nucleic Acids Research* 2022. Tumor Immune Syngeneic MOuse (TISMO) database.
- Endpoint used: `POST https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn`
  - `gene=Tacstd2`
  - `type=3` (table)
  - `icbList` = all in-vivo ICB treatments from `/tismo/gene/getVivoTreatment`
    (antiCTLA4, antiCTLA4&antiPD1, antiCTLA4&antiPDL1, antiPD1, antiPDL1, antiPDL2)
  - `tumorList` = all in-vivo cohorts from `/tismo/gene/getVivoCohort`
- Cancer-type labels: official `cancerType` from `/tismo/metaData/cellLineMeta`.
  Mammary carcinoma / adenocarcinoma / NOS are grouped as **Mammary** for the split.

## Local files

- `tacstd2_vivo_icb.csv` — raw TISMO table (619 rows, 64 slice labels, 605 unique SRX IDs).
- `sample_table.tsv` — same table as TSV.

No GEO FASTQ or count matrices were re-quantified. Values are TISMO's precomputed
in-vivo expression (`value` column) for mouse `Tacstd2`.
