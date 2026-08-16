# Results — extra mouse-lung ICI RNA (not TISMO)

**TISMO A4 is taken as given:** Tacstd2 up after ICB in 49/64 models, p=5.8e-5. Not re-audited. GSE155972 excluded.

GSE179500 is the known KL/LKB1 series (Tacstd2 higher when LKB1 is off). Other series below are additive.

Full numbers: `primary_contrasts.tsv`. Per-sample values: `per_sample.tsv`.

## 1. Response labels

No deposited mouse **lung** ICI RNA matrix in this slice has per-mouse responder vs non-responder labels. Closest proxies: GSE297630 (aPD-1-tolerant remaining LLC cells) and GSE246922 (ICB-relapsed CD45− cells). Those are not RECIST-style R/NR.

## 2. ICB-treated tumor RNA (extra series)

### Monotherapy PD-1 / PD-L1 vs control (no TLR / RT / chemo partner)

| Series | Model | n | Tacstd2 log2FC | Welch p | Cldn4 log2FC | Welch p |
|---|---|---|---|---|---|---|
| GSE114601 | KP GEMM nodules, aPD-1 vs vehicle | 2 / 2 | −1.34 | 0.49 | −1.14 | 0.23 |
| GSE157880 | HKP1 lung, aPD-1 0 Gy vs IgG 0 Gy | 3 / 2 | −0.71 | 0.46 | +0.46 | 0.40 |
| GSE262305 | LLC, aPD-L1 vs isotype (FPKM) | 3 / 3 | +0.03 | 0.65 | −0.07 | 0.51 |
| GSE274960 | LL/2, aPD-1 vs IgG (counts) | 3 / 3 | −0.39 | 0.60 | +0.93 | 0.32 |

None of these monotherapy contrasts reach p<0.05. GSE262305 Tacstd2 sits on the FPKM floor (~0.16).

### Combo / confounded ICB (still extra, still public)

| Series | Design | n | Tacstd2 | Cldn4 |
|---|---|---|---|---|
| **GSE239485** | LLC, Poly I:C + aPD-1 vs vehicle | 8 / 8 | **+2.09, Welch p=6.3e-4, MWU p=6.2e-4** | −0.77, p=0.24 |
| GSE239485 | LLC, Poly I:C + aPD-1 + C5aR1i vs vehicle | 8 / 8 | +2.58, p=9.7e-4 | +0.59, p=0.36 |
| GSE169194 | KPM total tumor, A2V+aPD-1 vs IgG | 3 / 3 | +0.05, p=0.82 | −0.24, p=0.44 |
| GSE157880 | HKP1, aPD-1+4 Gy vs IgG 0 Gy | 3 / 3 | +0.02, p=0.98 | +0.26, p=0.70 |
| GSE260596 | 344SQ, aPD-1+aLAIR1 vs aPD-1+aKLH | 3 / 4 | −0.78, p=0.61 | −0.01, p=0.99 |
| GSE197260 | EGFR NSCLC, gef→aPD-1 vs gef→veh | 1 / 1 | −5.55 | −6.51 |

GSE239485 is the only extra ICB tumor RNA with p<0.05 for Tacstd2. There is **no aPD-1 monotherapy arm**; Poly I:C is on every treated sample. Cldn4 does not move with the same contrast.

Radiation-only in GSE157880 (IgG 4 Gy vs IgG 0 Gy) is Tacstd2 +0.59, p=0.67 — the 4 Gy + aPD-1 rise is not ICB-specific.

### ICB-adjacent cell states (not pre vs post bulk tumor)

| Series | Design | n | Tacstd2 | Cldn4 |
|---|---|---|---|---|
| **GSE297630** | LLC aPD-1-tolerant remaining cells vs untreated | 3 / 3 | **+2.05 log2 signal, Welch p=4.2e-5** (MWU p=0.10; n=3 ties MWU) | −0.19, Welch p=0.016 |
| GSE246922 | LLC1 CD45− ICB-relapse vs parental (VST) | 3 / 3 | −0.014, p=0.12 | −0.029, p=0.12 |
| GSE246922 | KP CD45− ICB-relapse vs parental (VST) | 3 / 9 | −0.034, p=0.29 | +0.002, p=0.97 |
| GSE330941 | LLC Ago2KO (ICI-sensitized) vs WT, **no ICI on RNA** | 4 / 4 | +1.25, p=0.056 | +1.37, p=0.053 |

GSE297630 is the strongest extra ICB-linked Tacstd2 call. It is surviving cells after aPD-1, not a labeled R vs NR cohort.

## 3. KL / LKB1 genotype (cite GSE179500; add others)

Not ICB. These test whether Tacstd2 tracks LKB1-off / KL.

| Series | Design | n | Tacstd2 | Cldn4 |
|---|---|---|---|---|
| **GSE179500** (known) | Lkb1 non-restored vs restored, FACS tumor cells | 17 / 11 | sample-level log2FC +0.38, Welch p=0.042, MWU p=0.011. Author DESeq2 Restored vs Non-Restored: log2FC **−0.46**, p=0.0016, **padj=0.032** | sample +0.70, p=0.16; author log2FC −0.77, p=0.018, padj=0.16 |
| **GSE244452** | author DEG KP vs KL subcutaneous tumors | 3 / 3 | KL vs KP **+7.44**, author p=1.0e-69, padj=4.6e-67 (KP means ~9, KL ~1589) | KL vs KP **+8.61**, p=1.3e-25 |
| GSE137396 | KL vs KP GEMM nodules (log, median-centered) | 5 / 5 | +1.02, p=0.081 | +0.43, p=0.58 |
| GSE274351 | KL vs KP LCM adenomas (TPM) | 5 / 5 | −0.35, p=0.73 (sparse zeros) | −0.66, p=0.52 |
| GSE338923 | Lacun3 STK11 KO vs WT, vitro voom | 4 / 4 | **Tacstd2 not in table** | −2.78, p=0.016 |

GSE179500 direction matches the known statement: Tacstd2 is higher when LKB1 is off. GSE244452 is a much larger KL vs KP gap for both genes (subcutaneous KP/KL lines, no ICB on the RNA). GSE137396 is the same direction, n=5, p=0.081. LCM (GSE274351) and in-vitro STK11 KO (GSE338923) do not repeat a Tacstd2-up / Cldn4-up KL pattern.

## 4. Series that exist but were not scored as ICB RNA

| Accession | Why not in the ICB Tacstd2 table |
|---|---|
| GSE182228 | LKB1-deficient LUAD, vehicle / palbo / aPD-1 / combo, n=3. **No processed gene matrix** on GEO (metadata-only series matrix). |
| GSE298051 | KL/KP TNG260 ± aPD-1; deposited files are H3K9ac ChIP, not RNA. |
| GSE114300 | Kras GEMM aPD-1; CD4/CD8 TIL only. |
| GSE165517 | LLC aPD-L1; myeloid Nanostring panel. |
| GSE317011 | LLC-cocultured BMDM ± β-elemene; no ICB on the RNA. |
| GSE236258 | CMT-167 / KLA ± trametinib; no ICB RNA. |
| GSE241978 | CMT-167 AhR KO; no ICB RNA. |
| GSE180963 | KL vs K scRNA, n=2 sections; no ICB. |
| E-MTAB-15883 | LLC + CTX ± aPD-1; PD-1 fate-mapped immune cells only. |

## 5. How this sits next to TISMO

TISMO’s 49/64 Tacstd2-up call is an all-cancer, in-vivo ICB sign test and is left untouched. The extra **lung** ICB RNA that can be scored from public processed tables does not add a second monotherapy aPD-1/aPD-L1 lung series with p<0.05 for Tacstd2. What it does add:

1. **GSE239485** — Tacstd2 up after Poly I:C + aPD-1 in LLC (n=8, p=6.3e-4), with a TLR3 confounder.
2. **GSE297630** — Tacstd2 up in aPD-1-tolerant LLC cells (n=3, Welch p=4.2e-5).
3. **GSE179500 + GSE244452 + GSE137396** — Tacstd2 (and in GSE244452, Cldn4) higher in LKB1-off / KL than in LKB1-on / KP, which is genotype, not ICB.

Cldn4 does not track Tacstd2 on the extra ICB contrasts.
