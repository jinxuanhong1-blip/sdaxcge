# CLINICAL_HOOK — Opening paragraph evidence pack

**Scope.** Honest public citations for an Intro **sentence 1** about SKB264 / sacituzumab tirumotecan (sac-TMT; also MK-2870) plus ICI, and class-level TROP2-ADC + PD-1.  
**Rule.** Only numbers that appear in the cited peer-reviewed paper or named abstract. No invented NCTs, HRs, ORRs, or “data on file.”  
**Checked.** ClinicalTrials.gov API + PubMed/EuropePMC/Crossref (access date 2026-09-22). Private CD34 HIS / SKB264 lab work is **out of scope** (not on GitHub).

---

## Intro sentence 1 (use this; published only)

> In PD-L1–positive advanced NSCLC without targetable genomic alterations, first-line sacituzumab tirumotecan (sac-TMT; SKB264/MK-2870) plus pembrolizumab prolonged progression-free survival versus pembrolizumab alone in the randomized phase 3 OptiTROP-Lung05 trial (NCT06448312; median PFS not reached vs 5.7 months; stratified HR 0.35, 95% CI 0.26–0.47).

**Why this sentence.** It states a **published randomized clinical fact** (Lancet interim analysis). It does **not** claim that the combo works because of CLDN4/TROP2 immune exclusion, ADC-driven barrier opening, or PD-L1 induction — those are **not** established by OptiTROP-Lung05.

**Optional one-clause follow-on (still published; separate sentence).** Early-phase OptiTROP-Lung01 (NCT05351788) reported first-line sac-TMT plus the PD-L1 antibody tagitanlimab (KL-A167) with confirmed ORR 40.0% (cohort 1A, Q3W) and 66.7% (cohort 1B, Q2W) in Nature Medicine — note **PD-L1**, not PD-1.

---

## Story gap (what clinical data do *not* fill)

| Layer | Status for Intro |
|---|---|
| Combo **clinical activity** (sac-TMT + PD-1 in 1L PD-L1+ NSCLC) | **Published** — OptiTROP-Lung05 |
| Combo **clinical activity** (sac-TMT + PD-L1 in 1L NSCLC) | **Published** — OptiTROP-Lung01 (single-arm cohorts) |
| Class signal that TROP2-ADC + PD-1 can be active | **Published** for other ADCs (EVOKE-02; TROPION-Lung02), with important phase-3 caveats for SG |
| Mechanism: TROP2/CLDN4 **exclusion → ADC opens cold TME → ICI synergy** | **Speculative** relative to these trials — not a reported primary/secondary biomarker proof in OptiTROP-Lung05/01 abstracts |
| Mechanism: ADC raises PD-L1 / IFN after CLDN4 loss | **Not supported** by public omics catalogs in this repo (see prior C3 / public-analog PRs); do not cite as clinical fact |
| Global confirmatory sac-TMT + pembrolizumab (PD-L1 TPS ≥50%) | **Ongoing, no PubMed results** — NCT06170788 (TroFuse-007 / MK-2870-007) |

Keep CosMx / concordant-4 / TISMO biology in later Intro sentences as **independent public molecular evidence**, not as “shown in OptiTROP.”

---

## Tier A — Published SKB264 / sac-TMT ± ICI (cite for clinical hook)

### A1. OptiTROP-Lung05 — sac-TMT + pembrolizumab vs pembrolizumab (phase 3)

| Field | Verified value |
|---|---|
| Agent | sac-TMT (TROP2 ADC) + pembrolizumab (PD-1) vs pembrolizumab |
| Design | Randomized, open-label, phase 3; China, 68 hospitals |
| Population | Locally advanced/metastatic NSCLC; **no** targetable alterations; PD-L1 TPS ≥1% |
| N randomized | 413 (208 combo vs 205 pembro) |
| Primary endpoint | BICR PFS (ITT) |
| Interim result | Median follow-up 10.5 mo; median PFS **not reached** vs **5.7** mo; stratified **HR 0.35** (95% CI 0.26–0.47) |
| NCT | **NCT06448312** (ClinicalTrials.gov: ACTIVE_NOT_RECRUITING) |
| Citation | Xiong et al. (interim analysis). *Lancet*. 2026. PMID **42214392**. DOI **10.1016/S0140-6736(26)00968-2** |
| ASCO concurrent abstract | JCO 2026;44(16_suppl):8506. DOI **10.1200/JCO.2026.44.16_suppl.8506** (same NCT; reports BICR ORR 70.2% vs 42.0% and immature OS HR 0.55 — use only if citing that abstract explicitly) |

**Safety (Lancet abstract).** Grade ≥3 TEAEs: 115/208 (55%) combo vs 64/204 (31%) pembro.

**Do not over-read.** China-only phase 3; OS immature in the ASCO report; final analysis pending per Lancet. Does not prove TROP2 expression–driven selection or CLDN4 biology.

### A2. OptiTROP-Lung01 — sac-TMT + tagitanlimab / KL-A167 (phase 2)

| Field | Verified value |
|---|---|
| Agent | sac-TMT + **tagitanlimab (KL-A167, anti–PD-L1)** — not pembrolizumab |
| Design | Nonrandomized phase 2 cohorts 1A/1B; not powered for formal hypothesis testing |
| Population | 1L advanced/metastatic NSCLC **without** actionable genomic alterations |
| N | 40 (1A) / 63 (1B) |
| Confirmed ORR | **40.0%** (16/40) Q3W; **66.7%** (42/63) Q2W |
| DCR | 85.0% / 92.1% |
| Median PFS | **15.4** mo (95% CI 6.7–17.9) / **not reached** |
| NCT | **NCT05351788** |
| Citation | Hong et al. *Nat Med*. 2025. PMID **40830660**. DOI **10.1038/s41591-025-03883-5** |

**Labeling rule.** Call this **ADC + PD-L1**, not “SKB264 + PD-1,” unless the sentence explicitly names tagitanlimab/KL-A167.

### A3. sac-TMT monotherapy phase 3 (context only — not +ICI)

Useful so Intro does not imply that only the combo exists, but **do not** use as the +ICI hook.

| Trial | Setting | Key published result | NCT | Citation |
|---|---|---|---|---|
| OptiTROP-Lung04 | EGFRm nonsquamous NSCLC after EGFR-TKI; sac-TMT vs pemetrexed+platinum | Median PFS 8.3 vs 4.3 mo; HR 0.49 (0.39–0.62); OS HR 0.60 (0.44–0.82) | **NCT05870319** | *N Engl J Med* 2026; PMID **41124220**; DOI **10.1056/NEJMoa2512071** |
| Fang et al. BMJ | EGFRm NSCLC after EGFR-TKI **and** platinum; sac-TMT vs docetaxel | BIRC ORR 45% (41/91) vs 16% (7/45) | **NCT05631262** | *BMJ* 2025; PMID **40473437**; DOI **10.1136/bmj-2025-085680** |

Pretreated NSCLC monotherapy expansion (**KL264-01**, **NCT04152499**) has an AACR 2024 abstract (Cancer Res 84(7_Suppl):CT247; DOI 10.1158/1538-7445.AM2024-CT247) with ORR/PFS figures — cite as **abstract-level**, not as a full paper equivalent to A1–A3 journals above.

---

## Tier B — Published class-level TROP2-ADC + PD-1 (not SKB264)

Use only if Intro needs a class clause. Keep agent names explicit.

### B1. Sacituzumab govitecan (SG) + pembrolizumab — EVOKE-02 (phase 2)

| Field | Verified value |
|---|---|
| Design | Multicohort phase 2; SG + pembrolizumab 1L mNSCLC, no AGAs |
| Cohorts | A: PD-L1 TPS ≥50% (n=30); B: TPS <50% (n=62) |
| ORR (IRC) | **66.7%** (A); **29.0%** (B) |
| Median PFS | **13.1** mo (A); **7.0** mo (B) |
| Biomarker note | Trop-2 expression **did not correlate** with greater ORR/PFS in this report |
| NCT | **NCT05186974** |
| Citation | *J Thorac Oncol* 2026; PMID **41173143**; DOI **10.1016/j.jtho.2025.10.016** |

SG + pembrolizumab ± carboplatin cohorts also appear in *Clin Cancer Res* 2026 (PMID **41961582**; same NCT).

### B2. Datopotamab deruxtecan (Dato-DXd) + pembrolizumab — TROPION-Lung02 (phase Ib)

| Field | Verified value |
|---|---|
| Design | Phase Ib; Dato-DXd + pembro ± platinum; no AGAs |
| 1L doublet (treatment-naive) | Confirmed ORR **54.8%**; median PFS **11.2** mo (n=42) |
| 1L triplet | Confirmed ORR **55.6%**; median PFS **6.8** mo (n=54) |
| NCT | **NCT04526691** |
| Citation | *J Thorac Oncol* 2026; PMID **41871716**; DOI **10.1016/j.jtho.2026.103688** |

### B3. Class caveat (honest)

KEYNOTE-D46 / EVOKE-03 (**NCT05609968**; SG + pembrolizumab vs pembrolizumab, PD-L1 TPS ≥50%) was **discontinued** after sponsor/DMC communication that PFS improvement was not statistically significant and OS was unlikely to win — **full peer-reviewed primary paper not identified in PubMed under this NCT as of 2026-09-22**. Do **not** invent HRs. If mentioned, label as **sponsor-discontinued phase 3**, not as a journal result.

---

## Tier C — Speculative / not for sentence 1

Do **not** put these in Intro sentence 1 as established clinical facts:

1. **“SKB264 + PD-1 works by reversing CLDN4/TROP2-mediated T-cell exclusion.”**  
   OptiTROP-Lung05/01 report efficacy/safety; they do not publish a CLDN4 spatial or Tacstd2-exclusion mechanism package.

2. **“TROP2-ADC raises tumor PD-L1 / IFN and thereby sensitizes to PD-1.”**  
   Public RNA catalogs in this project did not support a clean PD-L1-up-after-CLDN4-loss / TROP2-ADC story (C3 and analog hunts). Clinical papers above are not a substitute.

3. **Global 1L sac-TMT + pembrolizumab superiority outside OptiTROP-Lung05.**  
   **NCT06170788** (TroFuse-007 / MK-2870-007; PD-L1 TPS ≥50%) is registered and recruiting; **PubMed hits for NCT06170788 / TroFuse-007 / MK-2870-007 = 0** on 2026-09-22. Cite as **ongoing**, never as completed efficacy.

4. **Any NCT, ORR, HR, or n not listed above.**  
   If it is not in the table, do not invent it.

5. **Private SKB264 / CD34 HIS / KD co-culture numbers.**  
   Per public handoff: not in this repository; another machine’s private data.

---

## Citation ledger (copy-paste; verified IDs only)

| Claim type | Short cite | PMID | DOI | NCT |
|---|---|---|---|---|
| sac-TMT + **PD-1** phase 3 PFS | Lancet OptiTROP-Lung05 | 42214392 | 10.1016/S0140-6736(26)00968-2 | NCT06448312 |
| sac-TMT + **PD-L1** phase 2 | Nat Med OptiTROP-Lung01 | 40830660 | 10.1038/s41591-025-03883-5 | NCT05351788 |
| sac-TMT mono vs chemo (EGFRm, post-TKI) | NEJM OptiTROP-Lung04 | 41124220 | 10.1056/NEJMoa2512071 | NCT05870319 |
| sac-TMT mono vs docetaxel (EGFRm, post-TKI+platinum) | BMJ | 40473437 | 10.1136/bmj-2025-085680 | NCT05631262 |
| SG + pembro phase 2 | JTO EVOKE-02 | 41173143 | 10.1016/j.jtho.2025.10.016 | NCT05186974 |
| Dato-DXd + pembro phase Ib | JTO TROPION-Lung02 | 41871716 | 10.1016/j.jtho.2026.103688 | NCT04526691 |
| Ongoing global sac-TMT + pembro (no results paper) | CT.gov only | — | — | NCT06170788 |

---

## One-paragraph “published vs speculative” blurb (methods/cover letter tone)

**Published:** First-line sac-TMT plus pembrolizumab improved BICR PFS versus pembrolizumab in PD-L1–positive advanced NSCLC without targetable alterations (OptiTROP-Lung05, NCT06448312; Xiong et al., *Lancet*). Earlier OptiTROP-Lung01 (NCT05351788; Hong et al., *Nat Med*) showed first-line sac-TMT plus PD-L1 blockade (tagitanlimab). Class-level TROP2-ADC + PD-1 activity is also reported for SG (EVOKE-02) and Dato-DXd (TROPION-Lung02), while SG’s PD-L1-high phase 3 companion was sponsor-discontinued without a PubMed primary paper under NCT05609968 at last check. **Speculative for Intro:** any claim that these clinical gains prove a CLDN4/TROP2 exclusion-breakage or PD-L1-induction mechanism; any results for NCT06170788; any private SKB264 experiment.

---

## Provenance

- NCTs resolved via `clinicaltrials.gov/api/v2/studies/{NCT}`.
- PMIDs/DOIs resolved via EuropePMC + Crossref (DOIs above resolve to the titled articles).
- No trial identifier in this file was guessed.
)
