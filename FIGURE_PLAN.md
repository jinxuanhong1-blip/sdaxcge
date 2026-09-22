# FIGURE_PLAN — Graphical abstract (7 boxes = discovery beats)

**Compile-only.** No new analysis. Numbers are copied from public locks / cited PRs.
Beats match the Intro discovery narrative (sibling `INTRO_NARRATIVE.md` beat list).

**Layout (left → right):** clinical open → TROP2 cold niche (observation) → Tacstd2 DEG → barrier → broad TJ null → CLDN4 screen → wet-lab pin → clinical close.

```
[1 clinical] → [2 TROP2 cold] → [3 DEG→barrier] → [4 TJ null] → [5 CLDN4 wins] → [6 KD pin] → [7 SKB264 close]
```

---

## Box map (which real number goes where)

| Box | Beat | Number to put in the box | Source lock | Do **not** put |
|---:|---|---|---|---|
| **1** | SKB264 + ICI unmet mechanism | **No public analysis lock.** Caption only: TROP2-ADC (SKB264) ± PD-1 — mechanism of cold niche / combo rationale is open. Optional side note (not the box headline): TISMO Tacstd2 **49/64** up after ICB (p=5.8×10⁻⁵) — resistance/up-after-ICB limb, **not** the immune-negative limb. | PUBLIC_HANDOFF; PR #542 / #543 / #587 | Invented ORR/PFS; “SKB264 clears cold niche” as fact |
| **2** | TROP2+ marks immune-neg niches (observation, not causation) | **CosMx** CD8+NK hi/lo **0.455** (10 µm) / **0.670** (20 µm); **8/8 & 5/5**; sign P=**0.031**. Companion: immune-frac **0.519 / 0.643**. **OncoSG** TACSTD2–CD8A purity-partial **ρ=−0.309** (p=4.7×10⁻⁵, n=169). | PR #726 / #738 / #739; PR #736 | Concordant-4 TACSTD2 vs T/NK (null ρ=−0.112); CosMx CLDN4 **0.36/0.52**; mediation |
| **3** | Tacstd2-high malignant DEG → junction / keratin / barrier | **Human concordant-4** (#741): ORA **HALLMARK_APICAL_JUNCTION rank #1** (enrichment 5.89, FDR 1.28×10⁻⁵); GSEA **GOBP_KERATINIZATION rank #2** among positive NES (NES +2.825, FDR 0.00146). Family scores: KERATIN logFC **+0.929**, EPITHELIAL_ADHESION **+0.734**, CLAUDIN_PANEL **+0.426**, TJ **+0.223**. Focal: CLDN4 logFC **+1.632**. | PR #741 | Public mouse GSEA TJ/adhesion/Claudin in top 3 (#744: `any_cohort_highlight_in_top3=false`). Private 8KL TJ#1 only if labeled **PRIVATE** |
| **4** | Broad TJ vs T/NK is null | Broad TJ (206 genes) vs T/NK DL meta **ρ=+0.137** (p=0.605, I²=70.5%, n=64). TJ core (7 genes) **ρ=−0.191** (p=0.406). Honest call: **null**. | PR #741 | “TJ-high tumors exclude T/NK”; mouse KEGG TJ as #1 |
| **5** | Gene-level TJ screen → CLDN4 best | **Primary:** malignant CLDN4 %pos vs T/NK **ρ=−0.531** (p=1.65×10⁻⁵, I²=0%, n=65); joint max of the #712 grid. **CosMx lock:** cytotoxic hi/lo **0.36 / 0.52** at 50/100 µm; **8/8 & 5/5**; sign P=**0.031**; companion effectors **1.11–1.22**, **0/8** decline → *exclusion, not muzzling*. **CPTAC LSCC protein:** CLDN4 vs ImmuneScore **ρ=−0.432** (n=78); TJ-15 **ρ=−0.082** (n=108, p=0.40) — CLDN4 within TJ, not whole TJ. Optional bridge (not mediation): Gygi lung TROP2–CLDN4 protein **ρ=0.693** (n=45). | PR #539 / #712; PUBLIC_HANDOFF / #698; PR #289 / #737; PR #740 | TACSTD2 CosMx 0.455/0.670 as if CLDN4; “CLDN4 mediates TROP2→cold” (#718/#726 null) |
| **6** | Functional pin CLDN1 / 4 / 7 KD | **No public number.** Box label: **planned wet-lab** (private CLDN1/4/7 KD co-culture). Do not invent KD fold-changes or IFN NES from public deposits (lung CLDN4 KD omics largely undeposited — see hunt PRs). | PUBLIC_HANDOFF “private”; planned experiment | Fabricated KD stats; public GSE207704 as if it proves the pin |
| **7** | Close: SKB264 clearing TROP2+ CLDN4-high cold zone | **Hypothesis caption**, not a measured endpoint. Visual: ADC hits TROP2+ / CLDN4-high cold niche → opens space for ICI. Reuse Box 2 + Box 5 numbers as the niche being targeted; do not invent clearance rates. | Narrative close only | Trial response attributed to this mechanism |

---

## ASCII schematic (copy into illustrator brief)

```
┌──────────────┐   ┌──────────────────┐   ┌─────────────────────┐
│ 1  SKB264+ICI│→ │ 2  TROP2+ cold   │→ │ 3  Tacstd2-high DEG │
│  unmet mech. │   │  CosMx 0.455/    │   │  ORA apical junc #1 │
│  (no lock #) │   │  0.670; OncoSG   │   │  keratin GSEA #2    │
│              │   │  ρ=−0.309        │   │  (#741)             │
└──────────────┘   └──────────────────┘   └──────────┬──────────┘
                                                     │
                                                     ▼
┌──────────────┐   ┌──────────────────┐   ┌─────────────────────┐
│ 7  SKB264    │← │ 6  CLDN1/4/7 KD  │← │ 5  CLDN4 wins screen│←┐
│  clears cold │   │  planned wet-lab │   │  ρ=−0.531; CosMx    │ │
│  zone (hyp.) │   │  (PRIVATE / no #)│   │  0.36/0.52; CPTAC   │ │
└──────────────┘   └──────────────────┘   │  CLDN4 −0.432 vs    │ │
                                          │  TJ-15 null          │ │
                                          └──────────▲───────────┘ │
                                                     │             │
                                          ┌──────────┴───────────┐ │
                                          │ 4  Broad TJ vs T/NK  │─┘
                                          │  ρ=+0.137  NULL      │
                                          │  (#741)              │
                                          └──────────────────────┘
```

Flow note: Box 4 is the **negative control / fork** after Box 3 — broad TJ fails, so the arrow into Box 5 is a **gene-level screen among TJ members**, not “TJ module = cold.”

---

## One-line captions (graphical-abstract text budget)

| Box | ≤12-word caption |
|---:|---|
| 1 | SKB264 ± ICI works; cold-niche mechanism unknown |
| 2 | TROP2-high niches are T/NK-poor (CosMx; OncoSG) |
| 3 | Tacstd2-high malignant cells up junction / keratin programs |
| 4 | Whole TJ score does **not** track low T/NK |
| 5 | Among TJ genes, CLDN4 best predicts exclusion |
| 6 | Pin with CLDN1/4/7 KD (wet-lab; private) |
| 7 | Hypothesis: SKB264 opens TROP2+ CLDN4-high cold zones for ICI |

---

## Provenance cheat-sheet (paste into figure legend footnotes)

| Claim | Exact lock string | PR / handoff |
|---|---|---|
| CosMx TACSTD2 short-range | CD8+NK **0.455 / 0.670** @ 10/20 µm; 8/8, 5/5; P=0.031 | #726, #738, #739 |
| OncoSG TACSTD2–CD8 | purity-partial **ρ=−0.309** | #736, #739 |
| Human DEG barrier | ORA apical junction **#1**; keratin GSEA **#2** | #741 |
| Mouse TJ GSEA top3 | **false** in all three cohorts | #744 |
| Broad TJ vs T/NK | **ρ=+0.137**, p=0.605 | #741 |
| CLDN4 vs T/NK | **ρ=−0.531**, I²=0%, n=65 | #539, #712 |
| CosMx CLDN4 exclusion | **0.36 / 0.52** @ 50/100 µm; 8/8, 5/5; P=0.031 | handoff / #698 |
| Not muzzling | effector hi/lo **1.11–1.22**; **0/8**↓ | handoff |
| CPTAC within-TJ | CLDN4 **ρ=−0.432** vs ImmuneScore; TJ-15 **−0.082** | #289, #737 |
| TROP2–CLDN4 protein bridge | Gygi lung **ρ=0.693** (n=45) | #740 |
| TISMO after ICB | Tacstd2 **49/64** up, p=5.8×10⁻⁵ | handoff / #542 |
| KL>KP entry (optional disease context, not a required GA box) | Tacstd2 **+3.24**, Cldn4 **+5.57**, MW p=0.00794 | handoff / #685 |

---

## Forbidden on this figure

- Mediation arrows (TROP2 → CLDN4 → cold) as causal — public mediation **null** (#718, #726).
- Public mouse GSEA “TJ #1” (#744).
- Concordant-4 TACSTD2 %pos vs T/NK as if it were −0.53 (that is **CLDN4**).
- Replacing CosMx CLDN4 **0.36/0.52** with TACSTD2 **0.455/0.670**.
- Merging private 8KL with public CosMx / mouse Harmony.
- Visium same-spot correlation as spatial exclusion.
- Fabricated clinical response or KD effect sizes in boxes 1, 6, or 7.
