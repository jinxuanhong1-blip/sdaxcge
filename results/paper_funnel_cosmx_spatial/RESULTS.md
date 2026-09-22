# CosMx spatial corroboration: TACSTD2-high cold neighbors + CLDN4/TJ program cold

## Purpose

PAPER FUNNEL pack for PPT slides. Public CosMx He2022 only (figshare 25976224). Two limbs plus the locked CLDN4 call:

| Limb | Question | Used |
|---|---|---|
| **Locked CLDN4** | cytotoxic neighbor ratio @ 50 / 100 µm | stated lock **0.36 / 0.52** (not recomputed) |
| **TACSTD2-high cold** | fewer immune / CD8+NK neighbors | yes (PR #726) |
| **CLDN4/TJ program cold** | epithelial program vs CD8 neighborhood | yes (PR #643; contact #567) |

No new CosMx download in this PR. No fabricated statistics. Machine table: `evidence_table.tsv`. Slide paste: `PPT_EVIDENCE_TABLE.md`. Figures: `figures/fig0_ppt_combined_strip.png` (and fig1–3).

## Overall verdict

| Limb | Supports immune-cold neighborhood? |
|---|---|
| Locked CLDN4 cytotoxic @ 50/100 µm | **Yes (lock)** — 0.36 / 0.52; 8/8 & 5/5 |
| Nearby effector muzzling | **No** — GZMB/PRF1/NKG7/IFNG hi/lo 1.11–1.22; 0/8↓ |
| TACSTD2-high @ 10–20 µm | **Yes** — CD8+NK 0.455 / 0.670; immune frac 0.519 / 0.643 |
| TACSTD2 @ 50–100 µm | **Partial / no** — not 8/8 |
| CLDN4 epithelial program @ 50 µm | **Yes** — ρ=−0.228; 5/5 & 8/8; shift p=0.005 |
| Classical multi-claudin TJ on CosMx | **Not testable** — only CLDN4 on 960-plex |

## What “TJ program” means here

On this panel CLDN4 is the only claudin. Hotspot partners are adhesion + keratins (CDH1, EPCAM, TACSTD2, EZR, KRT7/8/18/19, CEACAM6). That **program** is the cold CD8 neighborhood feature (PR #643). Do not invent off-panel TJ genes.

## Provenance

See `provenance.json`. Numeric cells cite PR paths that already exist on GitHub. Locked 0.36/0.52 is carried as the public handoff lock (also restated unchanged in PR #698, #629, #726, #643).

## Reproduce display only

```bash
# No analysis rerun. Figures were drawn from copied numbers:
python3 -c "print(open('results/paper_funnel_cosmx_spatial/evidence_table.tsv').readline())"
```
