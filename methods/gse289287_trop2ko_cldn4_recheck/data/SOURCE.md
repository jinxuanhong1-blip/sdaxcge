# Source

Public GEO series [GSE289287](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289287)
(*Trop-2 governs anti-metastatic desmosomal integrity*, Vacek / Souček lab; public 10 Feb 2026).

File used:

`GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz`

downloaded from

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE289nnn/GSE289287/suppl/GSE289287_DESeq2-Trop2KO_tumors_vs_WT.tsv.gz`

The table is the author DESeq2 contrast of T-47D Trop-2 KO xenografts versus WT xenografts.
It includes the Wald statistics and the per-sample normalized and raw counts that this recheck uses.
No FASTQ and no SRA recount.

FTP listing of `suppl/` (checked 21 Sep 2026) has seven DESeq2 tables. The only Trop-2 table is
`DESeq2-Trop2KO_tumors_vs_WT`. There is no `Trop2KO_cells` table. The 23 GSM records contain
Trop-2 KO samples only as xenografts (GSM8788422–GSM8788425). In vitro samples are WT
(GSM8788410–GSM8788412) and DSG2 KO (GSM8788413–GSM8788418).

Gene sets in `gene_sets.json` are the Hallmark IFN-γ, Hallmark IFN-α, and 21-gene MHC-I / APM
lists already frozen for the GSE289287 prerank GSEA (MSigDB Hallmark symbols; custom APM list).
They are copied here so this folder runs without that branch.
