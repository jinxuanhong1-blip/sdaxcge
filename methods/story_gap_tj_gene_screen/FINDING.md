# STORY GAP FILL: why CLDN4 wins the TJ → immune-cold screen

**Compile-only.** Numbers are quoted from named PRs. **Never fabricated.**
No new GEO / TCGA / CosMx re-fit on this page.

## Slide question

After Tacstd2-high DEG shows TJ / adhesion up, **which TJ/junction gene** tracks fewer T/NK/CD8 across **concordant-4 + TCGA + CosMx**?

## Funnel entry (locked priors)

| Step | Number | Source |
|---|---|---|
| GSE137244 KL>KP | Tacstd2 Δ=+3.24; Cldn4 Δ=+5.57; TJ7 Δ=+3.27; MW p=0.00794 | PR #685 / handoff |
| Concordant-4 Tacstd2 Q4 vs Q1 DEG | TJ family logFC=+0.223 (FDR 0.0042); CLDN4 logFC=+1.632; OCLN +0.989; CDH1 +1.032; TJP1 +0.593; CLDN1 +1.603 | PR #741 |
| Broad TJ score vs T/NK | DL ρ=**+0.137**, p=0.605, I²=70.5% — **null** | PR #741 |
| Locked CLDN4 %pos vs T/NK | ρ=**−0.531**, p=1.65×10⁻⁵, I²=0%, n=65 | PR #644 / #503 / #712 |

**Point of the gap:** broad TJ does not reproduce immune-cold. The member screen does — and CLDN4 leads it on concordant-4.

## Why CLDN4 wins the screen (PPT table)

Wins = most negative keratin-partial malignant **%pos vs patient T/NK** among genes that were actually scored on concordant-4.

| Gene | Tacstd2-DEG logFC (#741) | C4 KRT-partial ρ vs T/NK | C4 rank | TCGA KRT-partial ρ vs CD8 | TCGA rank | CosMx Q4/Q1 CD8+NK @50µm | CosMx |
|---|---:|---:|---:|---:|---:|---:|---|
| **CLDN4** | **+1.632** | **−0.478** (p=2.7e−4, I²=0%) | **#1** | −0.081 | **#5 / 6** | 0.800 (4/8, 3/5) | on panel; **#2 / 2** |
| CLDN7 | +0.713 | −0.352 (p=0.10, I²=56%) | #2 | −0.109 | #3 | — | **OFF 960** |
| CLDN3 | +0.652 | −0.346 (p=0.012, I²=0%) | #3 | −0.060 | #6 | — | **OFF 960** |
| CDH1 | +1.032 | −0.277 (p=0.073) | #4 | −0.126 | #2 | **0.786 (8/8, 5/5)** | on panel; **#1 / 2** |
| OCLN | +0.989 | −0.157 (p=0.27) | #5 | −0.101 | #4 | — | **OFF 960** |
| F11R | +0.447 | −0.070 (p=0.63) | #6 | **−0.131** | **#1** | — | **OFF 960** |
| CLDN1 | +1.603 | **+0.225** (opposite) | loses | not scored solo | GAP | — | **OFF 960** |
| CLDN5 | not in #741 focal table | **not scored** | GAP | **not scored** | GAP | — | **OFF 960** |
| TJP1 | +0.593 | **not scored solo** | GAP | TJ_15 module ρ=−0.291 (LUAD; #69) — not solo TJP1 | module only | — | **OFF 960** |

Sources: C4/TCGA/CosMx ranks for CLDN4/3/7/OCLN/F11R/CDH1 = PR **#748**. CLDN1 C4 = PR **#644**. CosMx panel status = PR **#643** / **#748**.

## What “wins” means (honest)

1. **Concordant-4 membership screen:** CLDN4 is the most negative KRT-partial ρ (−0.478). Unadjusted locked ρ = −0.531. CLDN1 has the **opposite** sign and is separable (Δ vs CLDN4, perm p < 0.0001; #644). Label-swap does **not** separate CLDN4 from CLDN7 (p=0.35; #748).
2. **TCGA bulk:** CLDN4 is **#5 of 6** on the same six-gene panel. F11R leads. Public bulk does **not** pin CLDN4.
3. **CosMx:** Only **CLDN4** among classical claudins is on the 960 panel. Same-cut Q4/Q1 ranks **CDH1** colder than CLDN4. Do **not** invent CosMx numbers for CLDN1/3/5/7/OCLN/F11R/TJP1. Locked CLDN4 cytotoxic ratios **0.36 / 0.52** (PR #698) are a different contrast and are **not** replaced.

## PPT one-liner

> After Tacstd2-high DEG lifts TJ, the **concordant-4 gene screen** picks **CLDN4** (ρ=−0.478; locked −0.531). CosMx/TCGA do **not** uniquely pick CLDN4 over the full TJ list (CDH1 / F11R lead those limbs; most claudins are off the CosMx panel).

## Do not say

- Do not invent CosMx ρ/ratios for off-panel genes.
- Do not invent C4 or TCGA ranks for CLDN5 or solo TJP1.
- Do not write “CLDN4 wins TCGA” or “CLDN4 wins CosMx same-cut.”
- Do not replace locked CosMx 0.36/0.52 with the Q4/Q1 0.80 row.
- Do not claim private KD co-culture on this public page.

## Files

- `ppt/PPT_SLIDE.md` — paste block
- `results/tables/screen_matrix.tsv` — full matrix
- `results/figures/why_cldn4_wins_screen.png` — PPT figure
- `provenance.json` — PR map
