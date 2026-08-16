# A4 TISMO breast / mammary models — Tacstd2 after ICB

Machine-readable headline is in `results/w200/A4_breast/summary.json`.
Do not pool the 29 TISMO boxes into one breast-wide ICB effect.

- Sample rows: 317 (unique samples 303)
- Cohorts: 29; latest-timepoint contrasts: 15
- Honest calls (all cohorts): {'null': 29}
- TISMO p-values reused across distinct cohorts: ['0.022387195755', '0.482356121608']
- Cohorts with both responders and non-responders: none

## Latest-timepoint contrasts

- 4T1_GSE130472_old_antiCTLA4: Δ=-0.413 MWU p=0.054 q=0.414 TISMO p=0.255 reused=no call=null (Non-responders, n=7+8)
- 4T1_GSE130472_old_antiPDL1: Δ=-0.596 MWU p=0.072 q=0.414 TISMO p=0.364 reused=no call=null (Non-responders, n=7+8)
- 4T1_GSE130472_young_antiCTLA4: Δ=0.194 MWU p=0.536 q=0.964 TISMO p=0.482 reused=yes call=null (Responders, n=7+12)
- 4T1_GSE130472_young_antiPDL1: Δ=-0.314 MWU p=0.270 q=0.700 TISMO p=0.482 reused=yes call=null (Responders, n=7+10)
- 4T1_GSE132529_antiPD1: Δ=0.012 MWU p=1.000 q=1.000 TISMO p=1.000 reused=no call=null (Non-responders, n=3+3)
- 4T1_GSE137818_Brca1_KO_antiPD1: Δ=-0.129 MWU p=1.000 q=1.000 TISMO p=0.843 reused=no call=null (Non-responders, n=5+5)
- 4T1_GSE137818_Brca2_KO_antiPD1: Δ=0.058 MWU p=0.764 q=0.964 TISMO p=0.991 reused=no call=null (Responders, n=5+3)
- E0771_GSE174053_high.fat.diet_antiPD1: Δ=0.010 MWU p=0.161 q=0.497 TISMO p=0.459 reused=no call=null (Non-responders, n=8+8)
- EMT6_GSE107801_antiPDL1: Δ=0.045 MWU p=0.290 q=0.700 TISMO p=0.022 reused=yes call=null (Responders, n=10+10)
- EMT6_GSE107801_antiTGFb_trap_antiPDL1: Δ=0.647 MWU p=0.054 q=0.414 TISMO p=0.022 reused=yes call=null (Responders, n=10+10)
- KPB25L_GSE124821_UV_end_antiCTLA4&antiPD1: Δ=-1.342 MWU p=0.171 q=0.497 TISMO p=0.636 reused=no call=null (Responders, n=6+4)
- T11_GSE124821_Apobec_end_antiCTLA4&antiPD1: Δ=0.252 MWU p=0.114 q=0.414 TISMO p=0.999 reused=no call=null (Responders, n=4+4)
- T11_GSE124821_end_antiCTLA4&antiPD1: Δ=-0.007 MWU p=1.000 q=1.000 TISMO p=0.996 reused=no call=null (Non-responders, n=5+5)
- p53-2225L_GSE124821_end_antiCTLA4&antiPD1: Δ=-0.094 MWU p=0.714 q=0.964 TISMO p=0.998 reused=no call=null (Non-responders, n=6+3)
- p53-2336R_GSE124821_end_antiCTLA4&antiPD1: Δ=2.018 MWU p=0.045 q=0.414 TISMO p=2.32e-04 reused=no call=null (Non-responders, n=8+5)
