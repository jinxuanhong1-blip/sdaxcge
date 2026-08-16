#!/usr/bin/env python3
"""Screen harvested records for claim C3.

C3 has two arms:
  A. PD-L1/CD274 rises AFTER CLDN4 loss (KD/KO/si/sh/CRISPR/depletion)
  B. PD-L1/CD274 rises AFTER TROP2 ADC exposure (SG, Dato-DXd, SKB264, etc.)

This is a title+abstract screen. It does not treat a hit as proven.
Lookalike claims (claudin-low, CLDN18.2, combo ICI+ADC trials without a
PD-L1 induction readout) are tagged so they are not counted as C3.

Outputs
-------
catalog/screened.tsv          every record with flags and a bucket
catalog/priority.tsv          records that survive the first filter
catalog/screen_summary.tsv    counts by bucket
"""

from __future__ import annotations

import csv
import os
import re
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
INP = os.path.join(ROOT, "catalog", "records.tsv")
OUT = os.path.join(ROOT, "catalog")

# --- tokens ----------------------------------------------------------------
RE_CLDN4 = re.compile(
    r"\b(cldn4|cldn-4|claudin[- ]?4|claudin4)\b", re.I
)
RE_CLDN_FAMILY = re.compile(
    r"\b(claudin[- ]?(low|18(?:\.2)?|3|6|7|9)|cldn18(?:\.2)?|cldn3|cldn6|cldn7|cldn9)\b",
    re.I,
)
RE_PDL1 = re.compile(
    r"\b(pd[- ]?l1|pdl1|cd274|b7[- ]?h1|programmed death[- ]ligand 1)\b", re.I
)
RE_PD1 = re.compile(r"\b(pd[- ]?1|pd1|cd279|pembrolizumab|nivolumab|atezolizumab|"
                    r"durvalumab|cemiplimab|tislelizumab|camrelizumab)\b", re.I)
RE_LOSS = re.compile(
    r"\b(knock[- ]?down|knock[- ]?out|kd\b|ko\b|silenc(?:e|ing)|crispr|"
    r"sirna|shrna|deplet(?:e|ion)|loss[- ]of[- ]function|lof\b|"
    r"null\b|ablat(?:e|ion)|deficient|deficiency)\b",
    re.I,
)
RE_RISE = re.compile(
    r"\b(up[- ]?regulat\w*|induc\w+|increas\w+|elevat\w+|rais\w+|"
    r"augment\w+|enhanc\w+|overexpress\w+|higher|upshift)\b", re.I
)
RE_TROP2 = re.compile(r"\b(trop[- ]?2|tacstd2|ega[- ]?1|ga733[- ]?1)\b", re.I)
RE_TROP2_ADC = re.compile(
    r"\b(sacituzumab|govitecan|trodelvy|datopotamab|dato[- ]?dxd|ds[- ]?1062|"
    r"sacituzumab tirumotecan|skb264|mk[- ]?2870|sac[- ]?tmt|"
    r"trop2[- ]?(adc|antibody[- ]drug)|anti[- ]trop[- ]?2)\b",
    re.I,
)
RE_ADC = re.compile(r"\b(adc|antibody[- ]drug conjugate|govitecan|deruxtecan|"
                    r"dxd|sn[- ]?38|payload)\b", re.I)
RE_TOP1 = re.compile(
    r"\b(sn[- ]?38|irinotecan|topotecan|camptothecin|deruxtecan|dxd|"
    r"topoisomerase[- ]?(i|1|i inhibitor))\b", re.I
)
RE_ICD = re.compile(r"\b(immunogenic cell death|cgas|sting|ifn[- ]?[abg]|"
                    r"type[- ]i interferon|calreticulin|hmgb1)\b", re.I)
RE_LOOKALIKE = re.compile(
    r"\b(claudin[- ]?low|cldn18(?:\.2)?|claudin[- ]?18(?:\.2)?|"
    r"zolbetuximab|osemiTAMAB|osemitamab)\b", re.I
)
RE_COMBO_TRIAL = re.compile(
    r"\b(phase [i123]|phase [123]|randomi[sz]ed|trial|cohort|"
    r"pembrolizumab|nivolumab|atezolizumab|durvalumab|"
    r"combination|plus pembro|plus nivo)\b", re.I
)
RE_REVIEW = re.compile(r"\b(review|editorial|commentary|perspective|consensus|"
                       r"guideline|overview)\b", re.I)
RE_ASSOC = re.compile(
    r"\b(correlat\w+|associat\w+|prognos\w+|co[- ]express\w+|"
    r"signature|biomarker)\b", re.I
)

# Sentence-level "X after/by Y" patterns. Cheap but useful as a flag.
RE_A_CAUSAL = re.compile(
    r"(cldn4|claudin[- ]?4).{0,80}(knock|silenc|crispr|sirna|shrna|deplet|loss|null|deficient)"
    r".{0,80}(pd[- ]?l1|pdl1|cd274)|(pd[- ]?l1|pdl1|cd274).{0,80}"
    r"(cldn4|claudin[- ]?4).{0,80}(knock|silenc|crispr|sirna|shrna|deplet|loss)",
    re.I | re.S,
)
RE_B_CAUSAL = re.compile(
    r"(sacituzumab|datopotamab|dato[- ]?dxd|skb264|trop[- ]?2 adc|sn[- ]?38|deruxtecan)"
    r".{0,100}(pd[- ]?l1|pdl1|cd274).{0,40}(up|induc|increas|elevat|rais)|"
    r"(pd[- ]?l1|pdl1|cd274).{0,80}(up|induc|increas|elevat).{0,80}"
    r"(sacituzumab|datopotamab|dato[- ]?dxd|skb264|trop[- ]?2|sn[- ]?38|deruxtecan)",
    re.I | re.S,
)


def flags(title: str, abstract: str) -> dict:
    text = f"{title}\n{abstract}"
    f = {
        "has_cldn4": bool(RE_CLDN4.search(text)),
        "has_cldn_family_other": bool(RE_CLDN_FAMILY.search(text)) and not bool(RE_CLDN4.search(text)),
        "has_pdl1": bool(RE_PDL1.search(text)),
        "has_pd1_axis": bool(RE_PD1.search(text)),
        "has_loss": bool(RE_LOSS.search(text)),
        "has_rise": bool(RE_RISE.search(text)),
        "has_trop2": bool(RE_TROP2.search(text)),
        "has_trop2_adc": bool(RE_TROP2_ADC.search(text)),
        "has_adc": bool(RE_ADC.search(text)),
        "has_top1": bool(RE_TOP1.search(text)),
        "has_icd": bool(RE_ICD.search(text)),
        "has_lookalike": bool(RE_LOOKALIKE.search(text)),
        "has_combo_trial": bool(RE_COMBO_TRIAL.search(text)),
        "is_reviewish": bool(RE_REVIEW.search(title)),
        "is_assoc_only_language": bool(RE_ASSOC.search(text)) and not bool(RE_LOSS.search(text)),
        "sent_a_causal": bool(RE_A_CAUSAL.search(text)),
        "sent_b_causal": bool(RE_B_CAUSAL.search(text)),
    }
    return f


def bucket(f: dict, arms: str) -> str:
    """Assign one mutually exclusive review bucket.

    Priority order is deliberate: lookalikes and reviews are peeled off
    before anything is allowed to look like a C3 hit.
    """
    if f["has_lookalike"] and not (f["has_cldn4"] and f["has_pdl1"] and f["has_loss"]):
        return "LOOKALIKE_not_C3"
    if f["is_reviewish"]:
        if f["has_trop2_adc"] and f["has_pdl1"]:
            return "REVIEW_trop2_adc_pdl1"
        if f["has_cldn4"] and f["has_pdl1"]:
            return "REVIEW_cldn4_pdl1"
        return "REVIEW_other"
    # Arm A: CLDN4 + PD-L1 + loss language
    if f["has_cldn4"] and f["has_pdl1"] and f["has_loss"]:
        if f["sent_a_causal"] or f["has_rise"]:
            return "PRIORITY_A_cldn4_loss_pdl1"
        return "WATCH_A_cldn4_pdl1_loss_no_rise_verb"
    if f["has_cldn4"] and f["has_pdl1"]:
        return "ASSOC_A_cldn4_and_pdl1_co_mention"
    # Arm B: TROP2 ADC + PD-L1
    if f["has_trop2_adc"] and f["has_pdl1"]:
        if f["sent_b_causal"] or (f["has_rise"] and not f["has_combo_trial"]):
            return "PRIORITY_B_trop2adc_pdl1_induction"
        if f["has_combo_trial"]:
            return "TRIAL_B_trop2adc_plus_ICI_no_induction_claim"
        return "WATCH_B_trop2adc_pdl1_co_mention"
    if f["has_top1"] and f["has_pdl1"] and f["has_rise"]:
        return "MECH_top1_payload_pdl1"  # prior, not TROP2-ADC-specific
    if f["has_trop2"] and f["has_pdl1"]:
        return "ASSOC_B_trop2_and_pdl1_co_mention"
    if f["has_cldn4"] and f["has_loss"]:
        return "CLDN4_loss_no_PDL1"
    if "M" in arms.split("|") and f["has_pdl1"] and (f["has_icd"] or f["has_top1"]):
        return "MECH_bridge_not_C3"
    return "IRRELEVANT_or_weak"


def main() -> None:
    rows = list(csv.DictReader(open(INP), delimiter="\t"))
    out_rows = []
    counts = Counter()
    for r in rows:
        f = flags(r["title"], r["abstract"])
        b = bucket(f, r["arms"])
        counts[b] += 1
        rec = dict(r)
        rec.update({k: ("Y" if v else "N") for k, v in f.items()})
        rec["bucket"] = b
        out_rows.append(rec)

    fieldnames = list(rows[0].keys()) + list(next(iter(out_rows)).keys() - rows[0].keys())
    # keep a stable order
    extra = [
        "has_cldn4", "has_pdl1", "has_loss", "has_rise", "has_trop2",
        "has_trop2_adc", "has_top1", "has_icd", "has_lookalike",
        "has_combo_trial", "is_reviewish", "sent_a_causal", "sent_b_causal",
        "bucket",
    ]
    fieldnames = [c for c in list(rows[0].keys()) + extra if c in (list(rows[0].keys()) + extra)]

    with open(os.path.join(OUT, "screened.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    pri_keep = {
        "PRIORITY_A_cldn4_loss_pdl1",
        "WATCH_A_cldn4_pdl1_loss_no_rise_verb",
        "ASSOC_A_cldn4_and_pdl1_co_mention",
        "PRIORITY_B_trop2adc_pdl1_induction",
        "WATCH_B_trop2adc_pdl1_co_mention",
        "MECH_top1_payload_pdl1",
    }
    pri = [r for r in out_rows if r["bucket"] in pri_keep]
    with open(os.path.join(OUT, "priority.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(pri)

    with open(os.path.join(OUT, "screen_summary.tsv"), "w") as fh:
        fh.write("bucket\tn\n")
        for k, n in counts.most_common():
            fh.write(f"{k}\t{n}\n")

    print("screened", len(out_rows))
    print("priority", len(pri))
    for k, n in counts.most_common():
        print(f"  {n:5d}  {k}")


if __name__ == "__main__":
    main()
