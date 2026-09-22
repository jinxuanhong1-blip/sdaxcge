# FINDING — PAPER FUNNEL: TJ / CLDN4 corroboration after Tacstd2 DEG

**Redraw pack only.** This folder does **not** re-download matrices, re-estimate
Spearmans, or invent p-values. It takes already-locked public results and
assembles a PPT-ready corroboration pack that prioritizes **CLDN4 within TJ**.

Locked numbers: `locked_numbers.tsv`. Figures: `figures/`.
Reproduce: `python3 methods/paper_funnel_tj_cldn4_corroboration/make_figures.py`.

---

## Funnel (one sentence)

After the public KL>KP Tacstd2 DEG (GSE137244), the immune-inverse story that
survives across patient scRNA, CosMx space, LR barrier, LSCC protein, and a
CLDN4-high bulk signature is **CLDN4**, not a broad TJ mean.

---

## Locked limbs (quoted)

| Limb | Locked number | Source |
|---|---|---|
| Entry DEG | Tacstd2 **+3.238**, Cldn4 **+5.570**, TJ_TISMO **+3.269**, MW p=**0.00794** (5 vs 5) | [PR #685](https://github.com/jinxuanhong1-blip/sdaxcge/pull/685) |
| Concordant-4 | malignant CLDN4 %pos vs T/NK DL meta **ρ = −0.531** (p=1.65e−5, I²=0%, n=65) | [PR #539](https://github.com/jinxuanhong1-blip/sdaxcge/pull/539) |
| CosMx exclusion | cytotoxic neighbor ratio **0.36** (50 µm) / **0.52** (100 µm); 8/8, 5/5; sign P=0.031 | [PR #698](https://github.com/jinxuanhong1-blip/sdaxcge/pull/698) (locked, not replaced) |
| CosMx not muzzling | nearby GZMB/PRF1/NKG7/IFNG hi/lo **1.11–1.22**, **0/8** decline | public CosMx handoff lock |
| LIANA barrier | T/NK and myeloid barrier family mean Δ > 0 on all three methods; BH q from **1.3e−9** to **4.9e−8** | [PR #579](https://github.com/jinxuanhong1-blip/sdaxcge/pull/579) |
| CPTAC LSCC | CLDN4 protein vs ImmuneScore **ρ = −0.432** (n=78); TJ-15 vs ImmuneScore **ρ = −0.082** (n=108, p=0.40) | [PR #289](https://github.com/jinxuanhong1-blip/sdaxcge/pull/289) |
| Signature #693 | ssGSEA α=0.75, size 163 vs CD8A meta **ρ = −0.533** (I²=0%, n sum=420) | [PR #693](https://github.com/jinxuanhong1-blip/sdaxcge/pull/693) |

Handoff TJ **+3.03** is a different average (stated in PR 685). This pack quotes
the TISMO 7-gene module **+3.269** from that PR.

---

## Why CLDN4 inside TJ (public)

The within-TJ priority argument that is public and protein-level is CPTAC LSCC:

- CLDN4 protein is inverse with ImmuneScore / GEP18 / CD8A (ρ ≈ −0.43 to −0.46, n=78).
- The structural TJ-15 protein score on the same freeze is null vs ImmuneScore / GEP18 / CD8A.
- Partial | WES keeps CLDN4–ImmuneScore negative (ρ=−0.298, p=0.0084).

Surface-molecule ranking that cannot nail CLDN4, and private KD co-culture that
can, are **out of this pack** (handoff rule).

---

## PPT slides (9 PNG / PDF pairs)

| File | Slide job |
|---|---|
| `fig01_funnel_overview` | Six-stage funnel map |
| `fig02_gse137244_entry` | Tacstd2 / Cldn4 / TJ entry DEG |
| `fig03_concordant4_rho` | Locked ρ = −0.531 forest |
| `fig04_cosmx_exclusion` | Exclusion ratios + not-muzzling |
| `fig05_liana_barrier` | Barrier family Δ bars |
| `fig06_cptac_cldn4_vs_tj15` | **Prioritize CLDN4 within TJ** |
| `fig07_signature_693` | Signature #693 meta −0.533 |
| `fig08_master_forest` | Cross-layer Spearman forest |
| `fig09_multipanel` | One appendix multipanel |

---

## Explicit non-claims

- Not a re-audit. Not “claim failed.”
- Visium same-spot correlation is not spatial exclusion.
- LIANA is an expression LR contrast, not CosMx contact.
- Signature #693 meta p is the p of a maximized |ρ| on those cohorts.
- CPTAC LSCC is treatment-naive resection protein, not ICI.
- Private 8-KL / CD34 HIS / CLDN KD co-culture / PDX are not used.
