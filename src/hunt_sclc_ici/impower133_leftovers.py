"""Published IMpower133 leftovers — numbers from papers, NOT patient-level data.

No EGA. No fabricated expression matrix. This module writes a citation table
of the only IMpower133 quantities that are actually public (Gay et al. 2021
Cancer Cell text/figure legends; Horn et al. 2018 NEJM trial report).

It is deliberately a static transcription, not an analysis.
"""
from __future__ import annotations

import json
from pathlib import Path

from . import bulk_george as bg

TABLES = bg.TABLES

# All values below are transcribed from the cited papers. They are cohort-level
# leftovers. They do NOT include TACSTD2 or CLDN4 (those genes were never
# reported for IMpower133 in a public table).

LEFTOVERS = [
    {
        "item": "trial_design",
        "value": "ES-SCLC, 1:1 carboplatin/etoposide +/- atezolizumab (anti-PD-L1)",
        "source": "Horn et al. 2018 NEJM (IMpower133 primary)",
        "usable_for_TACSTD2": "no",
    },
    {
        "item": "n_RNAseq_subtyped",
        "value": "276 treatment-naive tumors (Gay text); figure legend also cites 132 atezo + 139 placebo for OS",
        "source": "Gay et al. 2021 Cancer Cell Fig 1F / Fig 3 legend",
        "usable_for_TACSTD2": "no",
    },
    {
        "item": "subtype_distribution",
        "value": "SCLC-A 51%, SCLC-N 23%, SCLC-I 18%, SCLC-P 7%",
        "source": "Gay et al. 2021 Cancer Cell (IMpower133 paragraph)",
        "usable_for_TACSTD2": "no — no gene-level values",
    },
    {
        "item": "Ayers_GEP_in_SCLC-I",
        "value": "18-gene IFN-g T-cell GEP high specifically in SCLC-I (qualitative; heatmap Fig 3G)",
        "source": "Gay et al. 2021 Cancer Cell Fig 3G",
        "usable_for_TACSTD2": "no — heatmap pixels are not a data table",
    },
    {
        "item": "SCLC-I_median_OS_atezo",
        "value": ">18 months (EP+atezo arm)",
        "source": "Gay et al. 2021 Cancer Cell text (not a precise Kaplan-Meier table)",
        "usable_for_TACSTD2": "no",
    },
    {
        "item": "SCLC-I_median_OS_placebo",
        "value": "just over 10 months (EP+placebo arm); comparable to SCLC-A/N on placebo",
        "source": "Gay et al. 2021 Cancer Cell text",
        "usable_for_TACSTD2": "no",
    },
    {
        "item": "SCLC-I_vs_other_OS_HR_atezo_arm",
        "value": "HR=0.566 (95% CI 0.321-0.998)",
        "source": "Gay et al. 2021 Cancer Cell Supplementary Figure 4E-F (quoted in text)",
        "usable_for_TACSTD2": "no",
    },
    {
        "item": "SCLC-I_vs_other_OS_HR_placebo_arm",
        "value": "not significant (SCLC-I not merely prognostic)",
        "source": "Gay et al. 2021 Cancer Cell Supplementary Figure 4E-F",
        "usable_for_TACSTD2": "no",
    },
    {
        "item": "per_patient_RNAseq",
        "value": "EGA EGAS00001004888 / EGAD00001006926-8; DAC-controlled; embargoed full transcriptome",
        "source": "EGA; Gay 2021 data-availability statement",
        "usable_for_TACSTD2": "NO — not public; this hunt does not use it",
    },
    {
        "item": "TACSTD2_or_CLDN4_in_IMpower133",
        "value": "never reported in any public table, figure legend, or supplement we could access",
        "source": "this hunt (negative finding about the public record)",
        "usable_for_TACSTD2": "no — the gene was not published for this trial",
    },
]


def run():
    import pandas as pd
    df = pd.DataFrame(LEFTOVERS)
    df.to_csv(TABLES / "impower133_public_leftovers.tsv", sep="\t", index=False)
    verdict = {
        "can_compute_TACSTD2_vs_SCLC-I_in_IMpower133": False,
        "reason": "Per-patient RNA-seq is EGA-controlled; no public leftover contains TACSTD2 or CLDN4.",
        "what_we_did_instead": [
            "George 2015 bulk (Gay discovery cohort; public; n=81)",
            "Chan 2021 scRNA-seq atlas (public; SCLC donors)",
            "Transcribe Gay/Horn published IMpower133 cohort-level numbers only",
        ],
        "n_leftover_items": len(LEFTOVERS),
    }
    (TABLES / "impower133_verdict.json").write_text(json.dumps(verdict, indent=2))
    print(json.dumps(verdict, indent=2))
    return df


if __name__ == "__main__":
    run()
