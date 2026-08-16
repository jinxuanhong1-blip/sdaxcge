# Hunt: public paired pre/post ICB sets where TACSTD2/TROP2 rises

**Verdict:** No public paired ICB lung tumor RNA-seq set confirms a TACSTD2 rise. The only clean large human paired ICB tumor set (GSE91061 melanoma) is 24/43 up, Wilcoxon p=0.079. The largest paired tumor set that includes ICB (NeoTRIP GSE319641, n=150) shows TACSTD2 down. TISMO is 43/61 arm-level mouse contrasts up, not 49/64, and is not within-animal pairing.

User-cited starting points, treated as claims to recompute:

- TISMO “49/64”
- Zhejiang IHC “94 → 121”

## TISMO (mouse syngeneic) — not within-animal pairs

TISMO in vivo ICB samples are treated vs control **arms** of the same cell line in the
same study. Every ICB-treated sample has `Baseline=0`. This is after-vs-without at the
model level, not a longitudinal biopsy.

From `tismo/tismo_summary.json` (Tacstd2 present; 1,491/1,518 samples mapped):

| definition | n | up | frac | sign p |
|---|---:|---:|---:|---:|
| study × line × treatment × timepoint, matched control | 61 | 43 | 0.705 | 0.0019 |
| same, Tacstd2 above floor (max arm mean ≥ 0.5 log2) | 38 | 24 | 0.632 | 0.14 |
| study × line × treatment (timepoints pooled) | 48 | 33 | 0.688 | 0.013 |
| study × line × ICB class | 33 | 24 | 0.727 | 0.014 |
| study × line (all ICB pooled) | 31 | 24 | 0.774 | 0.0033 |

**Does not reproduce 49/64.** Closest raw numbers: 64 contrasts attempted, 3 unmatched,
43/61 up. Median Δ +0.07 log2. Only 7/61 nominally p<0.05 up. Lung: **one** contrast
(GSE155972 LLC anti-PD1+anti-CTLA4, Δ +0.46, MWU p=0.037).

## Human GEO — paired tumor RNA-seq with TACSTD2 actually scored

| GSE | cancer | treatment | n paired | up | Wilcoxon p | notes |
|---|---|---|---:|---:|---:|---|
| GSE91061 | melanoma | nivo | 43 | 24 | 0.079 | only clean large ICB-only pair |
| GSE115821 | melanoma | PD-1 ± CTLA-4 | 6 | 3 | 0.84 | MGH serial subset |
| GSE319641 | TNBC | NACT ± atezo | 150 | 38 | 2.5e-13 | **down**; arm not in GEO |
| GSE179351 | MSS CRC/PDAC | nivo+ipi+RT | 11 | 7 | 0.21 | radiation-confounded |

GSE91061 FPKM pairs: `geo/GSE91061_fpkm_pairs.csv` (rld sensitivity: 22/43 up, p=0.33).
NeoTRIP pairs: `geo/GSE319641_pairs.csv`.

## Lung-specific public sets that looked paired and were not

- **GSE207422** (Hu 2023 NSCLC neoadjuvant PD-1 + chemo): bulk log2TPM is 24 **pre-only**
  biopsies. scRNA metadata is 3 pre + 12 post from **different** patients (P01–P15 once each).
- **GSE248378** (neoadjuvant durvalumab ± RT NSCLC): resected tumors, no pre biopsy in GEO.
- **GSE135222** (Jung NSCLC anti-PD-1/PD-L1): baseline-only; TACSTD2 present.
- **GSE260770** (sintilimab GGO): blood exosomal RNA, not tumor pairs.

## Other seed-series outcomes (honest failures)

E-utilities union: 718 GSE (436 human, 260 mouse) in `geo/geo_candidate_series.json`.
Curated seed of 55 published ICB series walked by `scripts/geo_extract.py`:

- 2 generic “OK” hits: GSE91061 (real) and GSE318645 (PBMC counts; TACSTD2 mostly zero;
  tumor file is single-timepoint GBM resection, not paired tumor).
- The rest failed because the series matrix is empty (RNA-seq), TACSTD2 is absent,
  samples are baseline-only, or pre/post exist without a patient key
  (GSE227666 NeoPembrOV: Pre and Post titles, no patient ID).

Failure table: `geo/geo_paired_tacstd2_failures.csv`.

## Zhejiang IHC (“94 → 121”)

**Not recomputed.** No public patient-level paired IHC table was found. Closest published
numbers are different claims (e.g. Inomata et al. 2025: 94/110 stage III/IV TROP2-positive;
other papers report TROP2 largely stable after mixed anti-cancer therapy). Until a deposit
or supplement with paired H-scores appears, this stays a literature citation.

## Scripts

- `scripts/fetch_tismo.py` / `scripts/analyze_tismo.py`
- `scripts/geo_discover.py` / `scripts/geo_extract.py`
- `scripts/score_known_paired.py` / `scripts/write_hunt_report.py`
