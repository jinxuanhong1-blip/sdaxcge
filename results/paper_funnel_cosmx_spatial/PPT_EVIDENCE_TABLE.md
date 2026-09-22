# PAPER FUNNEL — CosMx spatial corroboration (TACSTD2 cold + CLDN4/TJ program cold)

**One slide table + figure strip.** Numbers are copied from existing CosMx public-analysis PRs only. Nothing was re-fit on the h5ad for this file. **Locked CLDN4 cytotoxic ratios 0.36 / 0.52 are not recomputed and not replaced.**

**Headline (honest):** On He et al. 2022 CosMx NSCLC (figshare 25976224; 8 sections / 5 patients), TACSTD2-high tumor cells are immune-colder at **10–20 µm**, and a CLDN4-containing epithelial adhesion/keratin **program** (not a classical multi-claudin TJ module — only CLDN4 is on the 960-plex) is colder for CD8 neighbors at **50 µm**. The pre-locked CLDN4 cytotoxic exclusion at **50/100 µm** remains **0.36 / 0.52**.

---

## PPT paste table

| # | Limb | Endpoint | Number to show | Supports cold? | One-line caveat |
|---|---|---|---|---|---|
| **L** | **LOCKED CLDN4** | cytotoxic neighbors @ 50 / 100 µm | **ratio 0.36 / 0.52; 8/8 & 5/5; sign P=0.031** | **YES (lock)** | Do not replace. Exclusion, not muzzling |
| L2 | LOCKED companion | GZMB/PRF1/NKG7/IFNG in nearby effectors | hi/lo **1.11–1.22; 0/8 decreased** | muzzling **NO** | Write *exclusion, not muzzling* |
| 1 | **TACSTD2-high cold** | CD8+NK neighbor count, section-median TACSTD2 | **10 µm 0.455; 20 µm 0.670; 8/8 & 5/5; P=0.031** | **YES** | Short-range; 50/100 µm not 8/8 |
| 2 | TACSTD2-high cold | immune-neighbor fraction, same split | **10 µm 0.519; 20 µm 0.643; 8/8 & 5/5** | YES | Companion to row 1 |
| 3 | TACSTD2 vs CLDN4 | immune fraction inside CLDN4-low stratum | 10 µm **0.507**; 20 µm **0.587**; still **8/8 & 5/5** | YES (not just CLDN4) | Attenuation of TACSTD2 after CLDN4 ≈0.12 |
| 4 | **CLDN4/TJ program** | Hotspot-smoothed 10-gene program vs CD8 @ 50 µm | patient-mean **ρ=−0.228; 5/5 & 8/8; shift p=0.005** | **YES** | Classical TJ genes absent except CLDN4 |
| 5 | Program cold fold | outer-quartile CD8 neighbor count | high/low **≈0.61×** (patient-mean tertile ratio **0.623**) | YES | Spec searched in 2,754-score grid |
| 6 | Keratin-free adhesion | CLDN4+CDH1+EPCAM+TACSTD2, same geometry | patient-mean **ρ=−0.213; 5/5 & 8/8** | YES | Keratin residual still −0.081, 5/5 |
| 7 | CLDN4 contact | any CD8 contact <20 µm | pooled OR **0.584**; MH OR **0.560**; **8/8 & 5/5** | YES | Complementary to radius counts |
| 8 | Squidpy calibration | median-split CLDN4 CD8+NK @ 50 µm | median ratio **0.883; 7/8** (not 5/5 donors) | WEAKER | Does **not** replace 0.36/0.52 |

---

## Figure strip (paste into slide)

| File | Use |
|---|---|
| `figures/fig0_ppt_combined_strip.png` | **Primary one-slide strip** (locked + TACSTD2 + program) |
| `figures/fig1_locked_cldn4_cytotoxic_ratios.png` | Lock callout |
| `figures/fig2_tacstd2_cold_neighbors.png` | TACSTD2 10–20 µm |
| `figures/fig3_cldn4_tj_program_cold.png` | Program / contact / lock |

Figures replot **copied** PR numbers only. They are not a new analysis of the CosMx object.

---

## Verdict line for the slide footer

> CosMx: locked CLDN4 cytotoxic exclusion **0.36 / 0.52** (50/100 µm). TACSTD2-high is colder at **10–20 µm** (CD8+NK **0.455 / 0.670**). CLDN4 sits in an adhesion/keratin program with patient-mean **ρ=−0.228** vs CD8@50 µm — classical TJ partners are off-panel. Exclusion, not muzzling.

---

## Do not say

- Do not replace **0.36 / 0.52** with TACSTD2 ratios or with Squidpy 0.883.
- Do not invent CLDN1/CLDN7/F11R/NECTIN2 on the CosMx 960-plex (absent).
- Do not write Visium same-spot correlation as spatial exclusion.
- Do not claim nearby effectors are muzzled (GZMB/PRF1/NKG7/IFNG not down).
- Do not merge private 8KL with this public CosMx object.
