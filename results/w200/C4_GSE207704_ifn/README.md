# C4 GSE207704 CLDN4 KO — IFN / MHC-I / APM genes (honest)

This folder did not exist in the checkout. Tables here are regenerated from the public GEO processed file.

## What this dataset actually is

- **GSE207704** (PMID 37059993): CRISPR **CLDN4−/− vs WT** in **T47D and MCF7 breast cancer** lines, not lung.
- Deposit is **group-collapsed cufflinks FPKM** (`GSE207704_CLDN4_RNAseq.txt`). Replicates (n=2/group) are already pooled.
- **No gene-level p-value or FDR is valid** from this file. Directions are descriptive.
- Contrast: `log2((KO + 0.5) / (WT + 0.5))`. Positive = higher after CLDN4 loss.
- UP/DOWN calls use `|log2FC| > 0.25` and max FPKM ≥ 1. Consensus requires the two lines not to go opposite ways.

## QC (must pass before any IFN claim)

| Gene | MCF7 WT→KO | T47D WT→KO | mean log2FC | call |
|---|---|---|---|---|
| CLDN4 | 85.85 → 50.85 | 42.51 → 20.44 | −0.89 | DOWN both, but **residual 48–59% mRNA** |
| TACSTD2 | 59.92 → 29.54 | 59.61 → 38.52 | −0.82 | DOWN both |

CLDN4 KO is real at mRNA but incomplete. TACSTD2 falls with it.

## C4 priority panel (IFI27, OAS2, IFIT1, MX1, ISG15, HLA-A)

| Gene | MCF7 | T47D | consensus |
|---|---|---|---|
| IFI27 | **absent** | **absent** | cannot score |
| OAS2 | **absent** | **absent** | cannot score |
| MX1 | **absent** | **absent** | cannot score |
| HLA-A | **absent** | **absent** | cannot score |
| IFIT1 | DOWN (−1.09) | FLAT (−0.03) | DOWN |
| ISG15 | FLAT (+0.11) | DOWN (−1.91) | DOWN (T47D-driven) |

**4/6 priority genes are missing from the deposit** (symbol and Entrez). That is an annotation gap, not a measured zero. The two genes that are present are **not up**.

## IFN / MHC-I / APM: up vs down

Among 82 panel genes excluding QC: **38 measured, 44 absent**.

### UP in KO (consensus)

Only three genes, and only **HERC5** is up in both lines.

- **HERC5** — both lines up (mean log2FC +1.52). The only clear IFN-related induction.
- **IFIT5** — MCF7 up, T47D flat. Weak / line-limited.
- **PSMB10** — MCF7 +0.36, T47D flat (−0.04). Not an MHC-I opening.

### DOWN in KO (consensus)

**IFIT1, ISG15, OAS1, OAS3, USP18, IFI6, IFI35, IFI44, BST2, IFITM1, IFITM3, STAT2, IRF1, IFNAR1, SOCS3, ERAP1, PSME1, PSME2, GBP2**

T47D is the cleaner down set (ISG15, BST2, IFI44, IFI6, OAS3). MCF7 is mixed but not an IFN-high state.

### Discordant (one line up, the other down)

**IFITM2, DDX58, IFIH1, JAK2, HLA-C**

HLA-C is the only classical MHC-I gene in the table: MCF7 up, T47D down.

### Flat

**IFI44L, EIF2AK2, PLSCR1, STAT1, JAK1, IFNAR2, IFNGR1, CALR, CANX, PDIA3, SEC61A1**

STAT1 is flat. Housekeeping APM (CALR/CANX/PDIA3/SEC61A1) is flat.

### Absent from the deposit (cannot evaluate)

**IFI27, OAS2, MX1, HLA-A, HLA-B, HLA-E/F/G, B2M, TAP1, TAP2, TAPBP, PSMB8, PSMB9, NLRC5, IRF7, IRF9, RSAD2, OASL, IFIT2, IFIT3, CXCL9/10/11, CIITA, IDO1**, and others listed in `up_vs_down.tsv`.

You cannot claim MHC-I/APM is unchanged or induced when **HLA-A/B, B2M, TAP1/2, PSMB8/9, NLRC5** are not in the file.

## Verdict

**GSE207704 does not support “CLDN4 KO opens IFN / MHC-I / APM.”**

- Measured ISGs are **mostly down**, especially in T47D.
- There is **no coordinated MHC-I/APM up** among genes that are actually present.
- The C4 six-gene signature is **mostly missing**, and the two present genes go **down**.
- Line discordance is large; a mean-of-lines story hides that.
- This is breast epithelium with residual CLDN4 mRNA, not a lung IFN-high KO.

Mouse lung **GSE50927** is a different experiment and is not evidence for this accession.

## Files

- `gene_stats.tsv` — per-gene FPKM and log2FC
- `up_vs_down.tsv` — the lists above
- `group_summary.tsv` — counts by set
- `analyze.py` lives at `scripts/w200/C4_GSE207704_ifn/analyze.py` (re-downloads GEO if raw is missing)
