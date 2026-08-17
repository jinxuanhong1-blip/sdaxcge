# FINDING — GSE136961 leftover: CLDN4 vs ICI response / CD8 / CD274

**Verdict: empty.** GSE136961 is public (21 pre-anti-PD-1 NSCLC tumors; author 9 DCB / 12 NDB), but **CLDN4 is not measured**. The deposited matrix is the Oncomine Immune Response Research Assay (GPL24014; 393 TPM rows / 395 unique raw genes). No CLDN family gene is on the panel. Therefore CLDN4 vs response, CLDN4 vs CD8, and CLDN4 vs CD274 are **n = 0**. No values were invented.

CD8A, CD8B, and CD274 *are* on the panel. Those are not CLDN4 tests and are not substituted here.

---

## Honest n

| Item | n | Note |
| --- | --- | --- |
| Public GEO series | **1** | [GSE136961](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136961), public 2020-02-03, PMID 31959763 |
| Tumor samples with deposited TPM | **21** | sample titles D01–D09 (9) and N01–N12 (12) |
| Author ICI outcome | **9 DCB / 12 NDB** | single-agent anti-PD-1; DCB = PR/CR or SD >24 weeks (RECIST 1.1) |
| Genes in `GSE136961_TPM.tsv.gz` | **393** | targeted immune panel, not transcriptome |
| Unique genes in `GSE136961_raw_count.tsv.gz` | **395** | 398 amplicon targets |
| **CLDN4 measured** | **0** | also absent: TACSTD2, TROP2, CLDN3, CLDN7, CLDN18 |
| CLDN4 vs response | **0** | empty |
| CLDN4 vs CD8 (CD8A/CD8B) | **0** | empty |
| CLDN4 vs CD274 | **0** | empty |
| CD8A / CD8B / CD274 on panel | 21 / 21 / 21 | present (`CD8A_14271531`, `CD8B_403510`, `CD274_461569`) |

## Why empty

This leftover is a real public lung ICI cohort with a per-patient outcome definition, but it cannot test CLDN4. The assay is a ~395-gene immune panel. Claudins and TACSTD2 are not among the deposited symbols (exact match plus aliases CLAUDIN4 / CPE-R / WBSCR8 / TROP2 / GA733-1). Re-probe: `python3 methods/gse136961_cldn4/probe.py`.

| Query | TPM | raw | IDs |
| --- | --- | --- | --- |
| CLDN4 | no | no | — |
| TACSTD2 / TROP2 | no | no | — |
| CD274 | yes | yes | `CD274_461569` |
| CD8A | yes | yes | `CD8A_14271531` |
| CD8B | yes | yes | `CD8B_403510` |
| PDCD1 | yes | yes | `PDCD1_143244` |
| PDCD1LG2 | yes | yes | `PDCD1LG2_315423` |

Files: `presence_table.tsv`, `gene_presence.json`, `panel_symbols.tsv`.

## What this is not

- Not a private or missing series. FTP TPM and raw counts downloaded and opened.
- Not a CLDN4 vs DCB/NDB, PFS, CD8, or PD-L1 result. Those tests have no CLDN4 vector.
- Not a CD8A or CD274 vs response substitute. Those genes exist; they are not the leftover question.
- Sample-title D/N matches the author 9/12 split; GEO also deposits PFS/OS. Labels were not needed because CLDN4 n = 0.

## 结论

GSE136961 公开、有抗 PD-1 结局（n=21，9 DCB / 12 NDB），但 Oncomine 免疫 panel **不含 CLDN4**。CLDN4 对疗效 / CD8 / CD274 的检验为 **空，n=0**。未编造表达值。
