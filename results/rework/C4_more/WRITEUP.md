# C4_more write-up

## Answer

**None.** No additional public CLDN4/Cldn4 KD/KO transcriptome was found
outside GSE50927, GSE207704, GSE22493, and GSE245459.

The C4 public analog therefore still rests on the same three genuine
CLDN4-loss expression series: mouse lung GSE50927 (IFN opens) versus the
two cancer-line sets GSE207704 and GSE22493 (IFN closes). GSE245459 is
not a CLDN4 experiment.

## What would have counted

A new hit required a public accession (GEO GSE, ArrayExpress/BioStudies
E-MTAB/E-GEOD unique study, or SRA study with a usable expression matrix)
in which CLDN4/Cldn4 itself is knocked down or knocked out and a
genome-wide transcriptome was deposited. No such extra accession was
located on 2026-08-16.

## What was searched

Live NCBI E-utilities dumps are in `live_search.json` (29 unique GSE
series after union of CLDN4 / Cldn4 / "claudin 4" / claudin-4 / KD-KO
terms, plus explicit lookups of GSE22421, GSE274940, GSE245459, and
GSE334497). Every series title, summary, assay type, and sample-title
list was read. SRA text searches for CLDN4/Cldn4 KD/KO RNA-seq returned
zero runs. BioStudies `E-MTAB AND CLDN4` returned zero hits. OmicsDI
returned only the known GEO mirrors.

The full accept/reject table is `screened.tsv`. Closest non-qualifying
records are `near_misses.tsv`.

## Closest records that still do not qualify

**GSE274940** is the only new claudin-loss RNA-seq that looks tempting.
It compares EpH4 WT versus a complete Cldn-null line (all claudin family
members removed). That is not a CLDN4-specific KD/KO and cannot be used
as a C4 analog.

**GSE22421** is a SKOV-3 microarray after C-CPE treatment. C-CPE binds
CLDN3/4 and can lower CLDN4 protein; it is not genetic knockdown or
knockout.

**GSE245459** (user-excluded) is TACSTD2/TROP2 shRNA in SKOV-3.
**GSE334497** is 4T1 Trop2 KO RNA-seq. Neither perturbs CLDN4.

Two 2025 papers report CLDN4-loss RNA-seq but deposit no public matrix:

- PMID 41016339 — H1688 CRISPR CLDN4 KO RNA-seq.
- PMID 40892111 — CLDN4 shRNA RNA-seq in a pancreatitis model.

Those cannot be used as public analogs.

## Implication for C4

There is still no second public non-cancer epithelial CLDN4-loss
transcriptome to test whether the IFN-open direction of GSE50927
replicates. Cancer-line CLDN4 loss remains represented only by
GSE207704 and GSE22493. This hunt does not add a fourth public
CLDN4 KD/KO transcriptome.
