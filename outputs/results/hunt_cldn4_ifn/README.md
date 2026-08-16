# Public CLDN4-loss → IFN / MHC-I / APM hunt

Search frozen: **2026-08-16**. Direction throughout is **CLDN4 KD/KO/siRNA
relative to its stated control**.

## Bottom line

I found **three accessioned public expression datasets with direct CLDN4/Cldn4
loss** and downloadable processed expression data. Only one contrast is a
credible directional match:

1. **GSE50927, uninjured whole lung, Cldn4 KO vs WT — partial IFN/APM match
   (top hit).** Five of six priority genes/proxies go up. Deposited edgeR
   results support **Oas2 +1.49 log2FC (FDR 0.0015), Isg15 +1.08 (0.0031),
   and Ifit1 +0.58 (0.028)**. The IFI27-family proxy Ifi27l2a (+0.62, FDR
   0.220) and Mx1 (+0.63, 0.648) agree in direction but are not significant.
   The HLA-A proxy H2-K1 is slightly down (-0.07, FDR 1.0). APM support is
   partial: B2m (+0.66, 0.034), Calr (+0.60, 0.0014), and Psmb9 (+0.90,
   0.037) rise significantly, while Tap1/2 and H2-K1 do not.
2. **GSE50927, VILI-high KO vs VILI WT — weak and confounded.** Five of six
   priority genes/proxies are directionally up, but none has FDR < 0.05 and
   these KO lungs were selected for greater injury. This cannot distinguish a
   CLDN4-loss response from a consequence of more severe injury.
3. **All other direct public contrasts fail the requested pattern.**
   GSE50927 VILI-low has only 2/6 priority genes up; GSE22493 SKOV-3 siRNA has
   only 2/5 measured priority genes up and IFIT1 is consistently down;
   GSE207704 MCF7/T47D CRISPR KO has only 2/6 priority genes in its deposited
   table, with no coherent increase.

Therefore the honest conclusion is **one partial public match, not a robust
cross-dataset replication and not evidence that MHC-I broadly rises after
CLDN4 loss**.

## Priority result: GSE50927 baseline

| Human target | Mouse measurement | log2FC | FDR | Call |
|---|---:|---:|---:|---|
| IFI27 | Ifi27l2a (family proxy) | +0.624 | 0.220 | up, not significant |
| OAS2 | Oas2 | +1.489 | 0.0015 | up |
| IFIT1 | Ifit1 | +0.576 | 0.028 | up |
| MX1 | Mx1 | +0.627 | 0.648 | up, not significant |
| ISG15 | Isg15 | +1.075 | 0.0031 | up |
| HLA-A | H2-K1 (mouse proxy) | -0.071 | 1.0 | not up |

The IFI27 mapping is not one-to-one. The deposited table also contains
Ifi27l1 (+0.152, FDR 0.907) and Ifi27l2b (-0.652, FDR 1.0); the main table
uses Ifi27l2a as the closest practical family proxy and labels it explicitly.

## Dataset-by-dataset assessment

### GSE50927 — mouse whole lung RNA-seq

- Direct germline Cldn4 KO; baseline and 2-hour ventilator-induced lung injury
  (VILI) contrasts.
- The baseline contrast is the least confounded and is the only qualifying
  result.
- Important limitation: the paper says two lungs/mice per condition, while
  GEO exposes one GSM per condition and only precomputed edgeR result tables,
  not replicate columns or a raw count matrix. The deposited p-values/FDRs
  therefore cannot be independently audited here. Whole-lung cell-composition
  changes are also possible.
- The VILI-high comparison is outcome-selected (greater injury in KO), so it is
  not a clean causal CLDN4 contrast.

### GSE22493 — SKOV-3 two-color spotted array

- Three paired arrays compare lentiviral CLDN4 siRNA against a CLDN4-high
  control (the metadata calls channel 1 “CLDN4 overexpression (control)”).
- The CLDN4 probe is down in its two non-missing replicates, supporting the
  sign convention.
- Priority median deposited log-ratios: IFI27 -2.32, OAS2 +0.11, IFIT1 -2.12,
  MX1 +0.36; ISG15 is absent; the multiple HLA-A probes are inconsistent.
- Old probe annotations, missing values, strong replicate/probe disagreement,
  and no deposited inferential statistics make this weak negative evidence.

### GSE207704 — MCF7 and T47D CRISPR KO RNA-seq

- Two WT and two KO samples are listed per cell line, but the only processed
  file collapses them to one FPKM value per condition and provides no
  p-values.
- Four priority genes (IFI27, OAS2, MX1, HLA-A) are absent from the deposited
  table. In MCF7, ISG15 is slightly up (+0.11 pseudocount-stabilized log2
  ratio) while IFIT1 is down (-1.60); both are down/non-up in T47D.
- This dataset does not support the requested signature.

## Scope and completeness

“Every public set” cannot be proven globally because repository metadata and
paper data-availability statements can be incomplete. Within the explicit
search below, these are all direct, accessioned CLDN4-loss transcriptome
datasets I could verify and analyze:

- GEO Entrez queries for `CLDN4`, `Cldn4`, `"claudin 4"`, and `"claudin-4"`
  restricted to GSE records (25 unique hits manually screened).
- GEO/SRA supplementary files and metadata; ArrayExpress/BioStudies and
  OmicsDI cross-indexes.
- Web/PubMed searches combining CLDN4/Cldn4 with knockdown, knockout, siRNA,
  shRNA, transcriptome, microarray, RNA-seq, accession, and dataset.
- LINCS L1000 releases GSE92742 and GSE70138: their public signature metadata
  contain no CLDN4 perturbation. The L1000KD2 curated file likewise contains
  zero CLDN4 knockdown signatures; generic portal text suggesting otherwise
  was not accepted as evidence.

Recent H1688 SCLC CLDN4-KO RNA-seq (PMID 41016339) and pancreatic-acinar
CLDN4-KD RNA-seq (PMID 40892111) were reported in papers, but no public
accession or downloadable expression matrix was located. They are leads, not
publicly analyzable sets, and are not counted as matches.

## Files

- `contrast_summary.tsv` — compact ranking and verdicts.
- `priority_signature.tsv` — all six requested genes/proxies in every direct
  contrast, including missingness and caveats.
- `apm_signature.tsv` — predefined MHC-I/APM panel.
- `GSE22493_probe_values.tsv` — raw deposited probe-level ratios used for the
  array assessment.
- `dataset_inventory.tsv` — included and excluded leads with reasons.
- `search_log.md` — queries, URLs, definitions, and analytical decisions.
- `analyze_public_sets.py` — dependency-free reproducible downloader/parser.

Run:

```bash
python3 outputs/results/hunt_cldn4_ifn/analyze_public_sets.py
```

The `_cache/` directory is disposable and git-ignored. Generated TSVs are kept
in version control.

## Primary sources

- [GSE50927](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50927);
  [paper, PMID 25106430](https://pubmed.ncbi.nlm.nih.gov/25106430/)
- [GSE22493](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22493)
- [GSE207704](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207704);
  [paper DOI 10.1186/s13058-023-01646-z](https://doi.org/10.1186/s13058-023-01646-z)
- [LINCS GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742)
  and [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)
- H1688 unpublished-data lead:
  [PMID 41016339](https://pubmed.ncbi.nlm.nih.gov/41016339/)
- Pancreatic-acinar unpublished-data lead:
  [PMID 40892111](https://pubmed.ncbi.nlm.nih.gov/40892111/)
