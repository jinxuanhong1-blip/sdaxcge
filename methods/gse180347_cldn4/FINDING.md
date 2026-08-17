# GSE180347 — CLDN4 / TACSTD2 vs LUAD hot/cold × PD-L1

**Verdict: panel missing. Stop.**

Neither **CLDN4** nor **TACSTD2** is on the NanoString nCounter PanCancer Immune Profiling Panel used in GSE180347 ([GPL29738](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL29738)). No scores vs hot/cold or PD-L1 are computed.

Kümpers et al., *Cancers* 2021 ([PMID 34572789](https://pubmed.ncbi.nlm.nih.gov/34572789/); GEO [GSE180347](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180347)). FFPE LUAD, IHC PD-L1 (TPS >50% vs <1%) × H&E immune infiltrate (hot ≥150 lymphocytes/HPF vs cold <150), then nCounter Immune Profiling.

## Honest n

| Item | n | Source |
|---|---:|---|
| GEO NanoString arrays | **144** | GSM5461051–GSM5461194 |
| Paper / GEO abstract cases | **138** | series summary; PMID 34572789 |
| Design-text groups | PH 23 + PC 7 + NH 53 + NC 55 = **138** | GEO `!Series_overall_design` |
| Array group labels | PH 26 + PC 7 + NH 58 + NC 53 = **144** | sample titles and `group:` characteristics |

144 is the number of public arrays. 138 is the number the paper and the GEO design text claim. The two group tallies do not match; the 6-array gap is not explained in the public metadata. Patient IDs on GEO recycle 01–12 within each cartridge run, so they are lane labels, not 144 unique patient IDs beyond the array count.

Do not treat n as 144 or 138 without saying which count.

## Panel check

Platform table: 784 probes (730 endogenous, 40 housekeeping, 8 negative, 6 positive). Series matrix: 730 quantified genes.

| Query | On GPL29738 | In GSE180347 matrix | Aliases / accessions searched |
|---|---|---|---|
| CLDN4 | **no** | **no** | CLDN4, NM_001305 |
| TACSTD2 | **no** | **no** | TACSTD2, TROP2, TROP-2, EGP1, EGP-1, GA733-1, M1S1, NM_002353 |

No claudin family member is on this panel. `NM_002354` is present and is **EPCAM (TACSTD1)**, not TACSTD2.

Epithelial / IO genes that *are* on the panel (not substitutes): EPCAM, CDH1, CD274.

## Why stop

The assigned test is CLDN4 or TACSTD2 vs hot/cold and PD-L1. Both targets are absent, so the contrast cannot be scored from this series. Nearby panel genes were not used as stand-ins.

## Reproduce

```bash
python3 methods/gse180347_cldn4/check_panel.py
```

Writes `methods/gse180347_cldn4/panel_check.json`.
