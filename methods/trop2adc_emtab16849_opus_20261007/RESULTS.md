# TROP2-ADC tumour-cell immune enhancement: E-MTAB-16849 and E-MTAB-16433

Run `trop2adc_emtab16849_opus_20261007`. Scored two public TROP2-ADC single-cell
datasets for the hypothesis that a TROP2-ADC raises tumour-cell antigen-presentation
(APM) and type-I/II interferon (ISG) programs and lowers CLDN4/tight-junction genes.

**Headline: neither dataset supports the APM-up / IFN-up claim.** E-MTAB-16849 is
null and its pre-dose negative control fails, so it cannot be scored as evidence
either way. E-MTAB-16433 is also null for APM and IFN (point estimates slightly
negative), but does show a nominally significant **CLDN4 decrease**, with the
caveat below. Reported exactly as computed; nothing was retuned or sign-flipped.

## Datasets scored

| Accession | Design | Control arm | Constraint-1 status | Unit of inference |
|---|---|---|---|---|
| E-MTAB-16849 | SG time course, human CRC organoid **liver metastases** in mice (HD4246); 0/6/24/48/72/120 h | **IgG1-SN-38 non-targeting ADC**, same experiment | Compliant (targeted ADC vs same-payload non-targeting control) | 9 paired SG-vs-control contrasts |
| E-MTAB-16433 | SG (Trodelvy) vs vehicle, CRC **PDOX** subcutaneous tumours (HD4246), 28 d dosing | vehicle | Compliant (perturbation is a TROP2-ADC) | 4 vs 4 independent mice |

Both are human colorectal cancer models, **not lung**, and are tumour-cell-only
(human-aligned xenograft material). They are cross-indication mechanistic support
at best, not lung-cancer evidence.

### Method

log2(x+1) on CPM pseudobulk; log2FC = treated − control. Core APM panel
(B2M, HLA-A/B/C/E/F, TAP1/2, PSMB8/9, NLRC5, TAPBP), core IFN panel
(STAT1, IRF1, CXCL9/10/11, GBP1/2/4, IRF7, ISG15, MX1, OAS1, IFIT1/3, IDO1),
STING axis (STING1/TMEM173, CGAS, IRF3), plus CLDN4 and TACSTD2 as sanity genes.
Gene-level Wilcoxon and a per-sample panel-score test. E-MTAB-16849 is a paired
design (signed-rank, paired t); E-MTAB-16433 is unpaired 4v4 (rank-sum, Welch t),
with a batch-paired signed-rank sensitivity analysis.

Both annotations use `TMEM173`/`CGAS`; the aliases `STING1`/`MB21D1` are absent.
`CXCL9` is absent from E-MTAB-16433, so its IFN panel has 14 genes.

No raw matrix was downloaded. `stream_extract.sh` pulls each counts matrix through
parallel HTTP range requests and keeps only panel-gene rows plus exact per-cell
library sizes, so peak disk was 17 MB (E-MTAB-16849, 4.96 GB source) and 2.2 MB
(E-MTAB-16433, 749 MB source). Line accounting was verified exactly: 27,594/27,594
and 25,470/25,470 gene rows processed exactly once, constant column count, no
malformed records.

---

## E-MTAB-16849 — SG vs IgG1-SN-38 non-targeting control

### Panel scores (primary, all 9 paired contrasts)

| Panel | Δ score (SG − control) | paired t p | signed-rank p | pairs positive |
|---|---|---|---|---|
| core_APM | **+0.066** | 0.362 | 0.250 | 6/9 |
| core_IFN | **+0.015** | 0.830 | 0.820 | 5/9 |
| STING_axis | +0.048 | 0.070 | 0.074 | 7/9 |
| sanity (CLDN4+TACSTD2) | +0.164 | 0.197 | 0.203 | 7/9 |

Post-dose only (n=8, dropping the 0 h pre-dose pair): core_APM +0.074 (p=0.363),
core_IFN **+0.0005** (p=0.995), STING_axis +0.037 (p=0.149).

### Gene level (all 9 pairs)

Only one gene reaches gene-level significance: **TMEM173/STING1 +0.239, 9/9 pairs
positive, signed-rank p=0.0039** (post-dose: +0.222, 8/8, p=0.0078). This is
fragile — TMEM173 averages 1.55 CPM and is detected in only 3.8% (SG) vs 2.9%
(control) of cells, so the effect rests on a handful of counts.

No APM gene is significant (best: HLA-B +0.151, p=0.164; HLA-A +0.118, p=0.250).
In the IFN panel only IRF1 trends (+0.163, p=0.055).

**Opposite-sign results, reported as found:** several canonical ISGs move *down*
under SG — ISG15 −0.083, MX1 −0.106, IFIT1 −0.109, IFIT3 −0.100, GBP2 −0.115,
STAT1 −0.018 — as do NLRC5 (−0.110) and IRF3 (−0.035).

**CLDN4 moves the wrong way for the hypothesis: +0.107 (up, 7/9 pairs, p=0.164)**,
alongside TACSTD2 +0.221 (p=0.164). Neither is significant, but the CLDN4 sign is
opposite to the predicted tight-junction decrease and is not flipped here.

### Why this dataset cannot be scored as evidence: the pre-dose control fails

0 h tumours were harvested before dosing, so any SG-vs-control difference at 0 h is
pure arm-to-arm technical/compositional offset. At cell level the 0 h contrast is
**larger in magnitude than the mean post-dose contrast for every panel**:

| Panel | Δ at 0 h (pre-dose) | Cohen's d at 0 h | mean Δ post-dose | post-dose range |
|---|---|---|---|---|
| core_APM | +0.308 | 0.186 | +0.087 | −0.488 … +0.415 |
| core_IFN | +0.234 | 0.187 | −0.059 | −0.587 … +0.236 |
| STING_axis | +0.349 | 0.221 | −0.058 | −0.272 … +0.054 |
| sanity | +0.618 | 0.261 | +0.127 | −0.243 … +0.468 |

Consistently, the per-pair core_APM delta correlates strongly with the
CLDN4/TACSTD2 offset across the 9 pairs (Pearson r=0.787, p=0.012): panel deltas
track a shared per-sample offset rather than a specific immune program. Two further
limits: each timepoint has only 1–2 replicates, and the chemokine/IDO1 arm of the
IFN panel is uninformative in this immune-free xenograft (detection: CXCL9 0.1%,
CXCL10 0.5%, CXCL11 1.3%, IDO1 0.5%). The pair `48h_repl1` is a large negative
outlier and also the lowest-depth sample (1,502 cells, 21.8 M counts vs ~55 M
typical); it was deliberately **retained** — no pair was dropped to move a result.

**Verdict (E-MTAB-16849): null / uninterpretable.** core_APM +0.066 (p=0.36) and
core_IFN +0.015 (p=0.83) are indistinguishable from zero, a pre-dose negative
control shows baseline arm offsets larger than any post-dose effect, and CLDN4 goes
up rather than down. Do not cite as support.

### Relationship to E-MTAB-16843

E-MTAB-16843 (reported by the user as core_APM +0.275, 8/8, p=0.0078; core_IFN
+0.268, 8/8 — **not recomputed in this run**) is the *in vitro* companion: per
BioStudies it is the same HD4246 CRC organoid line treated in culture with SG vs
IgG1-SN-38 across 0/3/6/9/12 h plus three post-washout timepoints. E-MTAB-16849 is
the *in vivo* liver-metastasis version of that contrast. The most parsimonious
reading is that the in vitro system delivers uniform, controlled drug exposure,
whereas the in vivo arms carry host-level exposure variability and arm-level batch
offsets large enough to swamp an effect of the size seen in vitro — exactly what the
failing 0 h control shows. So E-MTAB-16849 does **not** replicate E-MTAB-16843, and
it is also not a clean refutation: it is underpowered against its own noise floor.

---

## E-MTAB-16433 — SG (Trodelvy) vs vehicle, CRC PDOX, 4 vs 4 mice

Better design than E-MTAB-16849: four independent tumour-bearing mice per arm,
28 days of dosing, batches alternating between arms. With 4v4 the smallest
attainable two-sided rank-sum p is 0.0286.

### Panel scores (primary, unpaired 4v4)

| Panel | Δ score (SG − vehicle) | Welch t p | rank-sum p |
|---|---|---|---|
| core_APM | **−0.065** | 0.753 | 0.886 |
| core_IFN | **−0.039** | 0.858 | 0.686 |
| STING_axis | +0.121 | 0.147 | 0.114 |
| sanity (CLDN4+TACSTD2) | **−0.385** | 0.046 | 0.0286 |

Batch-paired sensitivity (pairing batches 0-1, 2-3, 4-5, 6-7) agrees: core_APM
−0.065 (3/4 positive, signed-rank p=0.875), core_IFN −0.039 (2/4, p=0.625),
STING_axis +0.121 (3/4, p=0.250), sanity −0.385 (**0/4 positive**, paired t p=0.040).

### Gene level

Significant at rank-sum p=0.0286 (the floor for 4v4):
- **CLDN4 −0.345** (Welch p=0.063) — *down, the predicted direction*
- **TACSTD2 −0.424** (Welch p=0.034) — down
- **GBP2 −0.565** (Welch p=0.072) — *down, opposite to the IFN-up prediction*

Also notable: **GBP1 −0.791** (rank-sum p=0.057, Welch p=0.040), down. TMEM173
trends up (+0.415, rank-sum p=0.114, Welch p=0.072), echoing the one positive in
E-MTAB-16849. IFIT3 (+0.471) and IFIT1 (+0.259) trend up but are far from
significant (p≈0.46–0.89), while HLA-E (−0.342), NLRC5 (−0.205), MX1 (−0.288) and
TAPBP (−0.145) trend down. No APM gene approaches significance.

**Important caveat on CLDN4.** The CLDN4 decrease is real in these data and in the
predicted direction, but TACSTD2 falls by a similar amount (−0.424) and the combined
sanity panel is the single strongest signal in the dataset (−0.385, 0/4 blocks
positive). A TROP2-ADC is expected to deplete TROP2-high cells and/or downregulate
its own target, so a concurrent epithelial/TROP2-program contraction is the simplest
explanation. The CLDN4 drop therefore cannot be cleanly separated from a shift in
tumour-cell composition or overall epithelial program, and should not be presented
as isolated tight-junction remodelling without an internal normaliser.

**Verdict (E-MTAB-16433): APM and IFN null (both slightly negative); CLDN4 down at
the 4v4 significance floor, confounded by a concordant TACSTD2/epithelial decline.**
Usable as modest support for the CLDN4-down arm only, with that caveat stated; it is
evidence *against* tumour-cell APM/ISG induction by SG in this model.

---

## Combined honest read

- Two independent in vivo TROP2-ADC datasets (one with a same-payload non-targeting
  control, one vehicle-controlled with true biological replicates) show **no
  tumour-cell APM induction and no ISG induction**. Point estimates are ≈0 to
  slightly negative, and several canonical ISGs (GBP1/GBP2, MX1, ISG15, IFIT1) move
  down.
- The only reproducible positive across both is **TMEM173/STING1** (+0.239, 9/9,
  p=0.0039 in E-MTAB-16849; +0.415, p=0.11 in E-MTAB-16433). Both rest on a
  low-abundance transcript (1.6 and 7.0 mean CPM), so this is a hypothesis, not a
  result.
- CLDN4 is **inconsistent across datasets**: up in E-MTAB-16849 (+0.107, n.s.), down
  in E-MTAB-16433 (−0.345, p=0.029 but confounded by TACSTD2).
- The APM-up / IFN-up signal reported for E-MTAB-16843 does not reproduce in either
  in vivo dataset. Any claim built on E-MTAB-16843 should be stated as an in vitro
  observation that currently lacks in vivo public replication.
- All of this is CRC, not lung. None of it substitutes for lung-cancer evidence for
  a SKB264+ICI paper.

## Files

| File | Contents |
|---|---|
| `stream_extract.sh` | Parallel range-request streamer; extracts panel rows + exact library sizes without storing the matrix |
| `assemble.py`, `assemble_16433.py` | Build per-cell panel tables and pseudobulk CPM |
| `analyze.py` | E-MTAB-16849 paired scoring (all-pairs and post-dose) |
| `secondary_percell.py` | E-MTAB-16849 cell-level effect sizes per timepoint (descriptive) |
| `diagnostics.py` | Pre-dose negative control and offset-correlation diagnostics |
| `analyze_16433.py` | E-MTAB-16433 4v4 scoring + batch-paired sensitivity |
| `data/*_panel_per_cell.csv.gz` | Per-cell panel counts, library size, sample metadata |
| `data/*_pseudobulk_cpm.csv` | Pseudobulk counts and CPM per group/mouse |
| `results/emtab16849_*.csv` | Gene-level, panel-score, per-pair, cell-level and diagnostic tables |
| `results/emtab16433_*.csv` | Gene-level, panel-score, batch-paired sensitivity tables |
| `results/*_stdout.txt` | Full console output of each analysis |

### Reproduce

```bash
pip install pandas scipy numpy
NW=24 OUT=/tmp/ae/stream bash stream_extract.sh                 # E-MTAB-16849 (~3 min)
python3 assemble.py && python3 analyze.py
python3 secondary_percell.py && python3 diagnostics.py
URL=".../E-MTAB-16433/Files/treatment_HD4246_SG_counts.txt" SIZE=748822660 \
  NW=16 OUT=/tmp/ae/stream433 bash stream_extract.sh            # E-MTAB-16433 (~40 s)
python3 assemble_16433.py && python3 analyze_16433.py
```

Small metadata inputs (`*_metadata.*`, `*_features.*`) come from the same
BioStudies `Files/` directories and are fetched with plain `curl`.
