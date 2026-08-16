#!/usr/bin/env python3
"""Assemble the hunt inventory + README from computed scores."""
import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "hunt_paired_up")

tismo = json.load(open(os.path.join(OUT, "tismo", "tismo_summary.json")))
known = json.load(open(os.path.join(OUT, "geo", "known_paired_scores.json")))
extract = json.load(open(os.path.join(OUT, "geo", "geo_extract_summary.json")))

# NeoTRIP computed separately
neo = pd.read_csv(os.path.join(OUT, "geo", "GSE319641_pairs.csv"))
neo_rec = {
    "gse": "GSE319641",
    "status": "ok",
    "n_paired_subjects": int(len(neo)),
    "n_up": int((neo.delta > 0).sum()),
    "n_down": int((neo.delta < 0).sum()),
    "frac_up": round(float((neo.delta > 0).mean()), 3),
    "median_delta": float(neo.delta.median()),
    "mean_pre": float(neo.pre.mean()),
    "mean_post": float(neo.post.mean()),
    "wilcoxon_p": 2.5118417913369247e-13,
    "sign_test_p": 1.1582626814448836e-09,
    "unit": "ComBat TPM",
    "title": "NeoTRIP TNBC: pre vs day-1-cycle-2 tumor RNA-seq (Gianni / Bianchini)",
    "cancer": "TNBC",
    "treatment": "NACT ± atezolizumab (arm not in GEO characteristics)",
    "pairing": "same patient, Baseline vs D1C2 tumor biopsy",
    "note": "Largest public paired tumor set touching ICB. TACSTD2 falls. Cannot subset atezo-only from deposited metadata.",
}

rows = []
for r in known["primary"]:
    if r.get("gse") == "GSE115821" and r.get("n_paired_subjects", 0) == 0:
        continue
    rows.append(r)
rows.append(neo_rec)

# inventory of seed outcomes
fail_rows = []
for f in extract.get("failures", []):
    fail_rows.append(
        {
            "gse": f.get("gse"),
            "status": "fail",
            "reason": f.get("reason"),
            "title": (f.get("title") or "")[:160],
            "n_paired_subjects": f.get("n_paired_subjects"),
            "gene_row": f.get("gene_row"),
        }
    )

inv = pd.DataFrame(rows)
inv.to_csv(os.path.join(OUT, "paired_tacstd2_confirmatory_scores.csv"), index=False)
pd.DataFrame(fail_rows).to_csv(os.path.join(OUT, "geo", "geo_seed_failures_classified.csv"), index=False)

summary = {
    "question": "Do public paired pre/post ICB lung or pan-cancer sets show TACSTD2/Tacstd2 rise after ICB?",
    "user_claims": {
        "TISMO": "49/64",
        "Zhejiang_IHC": "94 → 121",
    },
    "tismo": {
        "what_it_is": "Mouse syngeneic treated vs control ARMS, not within-animal pairs. All ICB samples have Baseline=0.",
        "attempted_contrasts": tismo["contrasts_attempted"],
        "with_matched_control": 61,
        "n_up": 43,
        "frac_up": 0.705,
        "sign_test_p": 0.001868101030975251,
        "median_delta_log2": 0.0707,
        "n_up_p_lt_0.05": 7,
        "expressed_models_up": "24/38 (sign p=0.14)",
        "lung_contrasts": "1 (GSE155972 LLC, Δ+0.46, MWU p=0.037)",
        "reproduces_49_of_64": False,
    },
    "human_paired_tumor_rnaseq_with_tacstd2": [
        {
            "gse": "GSE91061",
            "cancer": "melanoma",
            "n": 43,
            "n_up": 24,
            "wilcoxon_p": 0.079,
            "direction": "weak up, not significant",
            "clean_icb_only": True,
        },
        {
            "gse": "GSE115821",
            "cancer": "melanoma",
            "n": 6,
            "n_up": 3,
            "wilcoxon_p": 0.84,
            "direction": "null",
            "clean_icb_only": True,
        },
        {
            "gse": "GSE319641",
            "cancer": "TNBC",
            "n": 150,
            "n_up": 38,
            "wilcoxon_p": 2.5e-13,
            "direction": "DOWN",
            "clean_icb_only": False,
            "note": "NACT ± atezolizumab; arm not deposited",
        },
        {
            "gse": "GSE179351",
            "cancer": "MSS CRC / PDAC",
            "n": 11,
            "n_up": 7,
            "wilcoxon_p": 0.21,
            "direction": "weak up, not significant",
            "clean_icb_only": False,
            "note": "nivo+ipi+radiation",
        },
    ],
    "key_negatives": [
        "GSE207422 NSCLC neoadjuvant PD-1+chemo: bulk TPM is pre-only; scRNA pre/post are different patients.",
        "GSE227666 NeoPembrOV ovarian NACT±pembro: Pre and Post samples exist but no patient ID to pair them.",
        "GSE248378 neoadjuvant durvalumab NSCLC: resected (post) tumors only.",
        "GSE135222 / GSE78220 / GSE145996 / GSE173839 / GSE194040 / GSE176307: TACSTD2 present, baseline-only ICB cohorts.",
        "GSE318645 GBM window trial: tumor is single-timepoint resection; PBMC pairing is not tumor TROP2.",
    ],
    "zhejiang_ihc": {
        "recomputed": False,
        "reason": "No public patient-level paired IHC table. Closest published numbers (Inomata 2025: 94/110 stage III/IV TROP2-positive; other papers report TROP2 largely stable after mixed therapy) do not match a 94→121 paired H-score.",
    },
    "verdict": (
        "No public paired ICB lung tumor RNA-seq set confirms a TACSTD2 rise. "
        "The only clean large human paired ICB tumor set (GSE91061 melanoma) is 24/43 up, Wilcoxon p=0.079. "
        "The largest paired tumor set that includes ICB (NeoTRIP GSE319641, n=150) shows TACSTD2 down. "
        "TISMO is 43/61 arm-level mouse contrasts up, not 49/64, and is not within-animal pairing."
    ),
}
with open(os.path.join(OUT, "hunt_summary.json"), "w") as fh:
    json.dump(summary, fh, indent=2)

readme = f"""# Hunt: public paired pre/post ICB sets where TACSTD2/TROP2 rises

**Verdict:** {summary['verdict']}

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
"""
open(os.path.join(OUT, "README.md"), "w").write(readme)
print(json.dumps(summary, indent=2)[:3000])
print("wrote", os.path.join(OUT, "README.md"))


if __name__ == "__main__":
    pass
