#!/usr/bin/env python3
"""Compare galectin-family expression in TACSTD2-high vs -low TCGA lung tumors."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
from scipy.stats import mannwhitneyu, spearmanr


API = "https://www.cbioportal.org/api"
STUDIES = {
    "LUAD": "luad_tcga_pan_can_atlas_2018",
    "LUSC": "lusc_tcga_pan_can_atlas_2018",
}
GENES = {
    "TACSTD2": 4070,
    "LGALS1": 3956,
    "LGALS2": 3957,
    "LGALS3": 3958,
    "LGALS4": 3960,
    "LGALS7": 3963,
    "LGALS7B": 653499,
    "LGALS8": 3964,
    "LGALS9": 3965,
    "LGALS9B": 284194,
    "LGALS9C": 654346,
    "LGALS12": 85329,
    "LGALS13": 29124,
    "LGALS14": 56891,
    "LGALS16": 148003,
}
OUT = Path(__file__).resolve().parents[3] / "results" / "w200" / "A11_galectin"


def get_json(session: requests.Session, url: str) -> object:
    response = session.get(url, timeout=60)
    response.raise_for_status()
    return response.json()


def fetch_study(session: requests.Session, study: str) -> tuple[list[str], dict[str, dict[str, float]]]:
    sample_list = f"{study}_rna_seq_v2_mrna"
    profile = f"{study}_rna_seq_v2_mrna"
    sample_ids = get_json(session, f"{API}/sample-lists/{sample_list}")["sampleIds"]
    response = session.post(
        f"{API}/molecular-profiles/{profile}/molecular-data/fetch",
        params={"projection": "DETAILED"},
        json={"entrezGeneIds": list(GENES.values()), "sampleListId": sample_list},
        timeout=120,
    )
    response.raise_for_status()
    by_sample: dict[str, dict[str, float]] = {sample: {} for sample in sample_ids}
    entrez_to_symbol = {entrez: symbol for symbol, entrez in GENES.items()}
    for row in response.json():
        symbol = entrez_to_symbol.get(row["entrezGeneId"])
        if symbol and row.get("value") is not None:
            by_sample[row["sampleId"]][symbol] = float(row["value"])
    return sample_ids, by_sample


def bh_adjust(rows: list[dict[str, object]], p_key: str, q_key: str) -> None:
    valid = [(i, float(row[p_key])) for i, row in enumerate(rows) if row[p_key] != ""]
    valid.sort(key=lambda item: item[1])
    m = len(valid)
    adjusted = [1.0] * m
    running = 1.0
    for rank_index in range(m - 1, -1, -1):
        rank = rank_index + 1
        running = min(running, valid[rank_index][1] * m / rank)
        adjusted[rank_index] = running
    for (row_index, _), q_value in zip(valid, adjusted):
        rows[row_index][q_key] = q_value


def fmt(value: object) -> str:
    if value == "":
        return ""
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "Content-Type": "application/json"})
    stats: list[dict[str, object]] = []
    expression_rows: list[dict[str, object]] = []
    cohort_counts: dict[str, int] = {}

    for cohort, study in STUDIES.items():
        sample_ids, by_sample = fetch_study(session, study)
        complete = [
            sample
            for sample in sample_ids
            if "TACSTD2" in by_sample[sample]
        ]
        complete.sort(key=lambda sample: (by_sample[sample]["TACSTD2"], sample))
        cohort_counts[cohort] = len(complete)
        split = len(complete) // 2
        low_samples = complete[:split]
        high_samples = complete[-split:]
        group_by_sample = {sample: "low" for sample in low_samples}
        group_by_sample.update({sample: "high" for sample in high_samples})

        for sample in complete:
            row: dict[str, object] = {
                "cohort": cohort,
                "study_id": study,
                "sample_id": sample,
                "tacstd2_group": group_by_sample.get(sample, "median_unassigned"),
            }
            row.update({gene: by_sample[sample].get(gene, "") for gene in GENES})
            expression_rows.append(row)

        tac_samples = low_samples + high_samples
        tac_values = np.asarray([by_sample[s]["TACSTD2"] for s in tac_samples])
        for gene in list(GENES)[1:]:
            paired_samples = [s for s in tac_samples if gene in by_sample[s]]
            x = np.asarray([by_sample[s]["TACSTD2"] for s in paired_samples])
            y = np.asarray([by_sample[s][gene] for s in paired_samples])
            low = np.asarray([by_sample[s][gene] for s in low_samples if gene in by_sample[s]])
            high = np.asarray([by_sample[s][gene] for s in high_samples if gene in by_sample[s]])
            if not len(y) or np.all(y == y[0]):
                stats.append(
                    {
                        "cohort": cohort,
                        "gene": gene,
                        "n": len(y),
                        "n_low": len(low),
                        "n_high": len(high),
                        "median_low": float(np.median(low)) if len(low) else "",
                        "median_high": float(np.median(high)) if len(high) else "",
                        "log2_median_ratio": "",
                        "spearman_rho": "",
                        "spearman_p": "",
                        "mannwhitney_p": "",
                        "spearman_q": "",
                        "mannwhitney_q": "",
                    }
                )
                continue
            rho, rho_p = spearmanr(x, y)
            mw = mannwhitneyu(high, low, alternative="two-sided", method="asymptotic")
            stats.append(
                {
                    "cohort": cohort,
                    "gene": gene,
                    "n": len(y),
                    "n_low": len(low),
                    "n_high": len(high),
                    "median_low": float(np.median(low)),
                    "median_high": float(np.median(high)),
                    "log2_median_ratio": math.log2(
                        (float(np.median(high)) + 1) / (float(np.median(low)) + 1)
                    ),
                    "spearman_rho": float(rho),
                    "spearman_p": float(rho_p),
                    "mannwhitney_p": float(mw.pvalue),
                    "spearman_q": "",
                    "mannwhitney_q": "",
                }
            )

        # The TACSTD2 vector is retained here as a sanity check that both groups are balanced.
        assert len(low_samples) == len(high_samples)
        assert np.median(tac_values[:split]) <= np.median(tac_values[-split:])

    # Correct across both histologies and all tested galectin genes.
    bh_adjust(stats, "spearman_p", "spearman_q")
    bh_adjust(stats, "mannwhitney_p", "mannwhitney_q")

    stats_fields = [
        "cohort",
        "gene",
        "n",
        "n_low",
        "n_high",
        "median_low",
        "median_high",
        "log2_median_ratio",
        "spearman_rho",
        "spearman_p",
        "spearman_q",
        "mannwhitney_p",
        "mannwhitney_q",
    ]
    with (OUT / "galectin_stats.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=stats_fields, delimiter="\t")
        writer.writeheader()
        for row in stats:
            writer.writerow({field: fmt(row[field]) for field in stats_fields})

    expression_fields = ["cohort", "study_id", "sample_id", "tacstd2_group", *GENES]
    with (OUT / "expression_values.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=expression_fields, delimiter="\t")
        writer.writeheader()
        for row in expression_rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in expression_fields})

    by_gene = {
        gene: [row for row in stats if row["gene"] == gene]
        for gene in list(GENES)[1:]
    }
    concordant = []
    discordant = []
    for gene, rows in by_gene.items():
        if len(rows) != 2 or any(row["spearman_rho"] == "" for row in rows):
            continue
        directions = [math.copysign(1, float(row["spearman_rho"])) for row in rows]
        if directions[0] == directions[1] and all(float(row["spearman_q"]) < 0.05 for row in rows):
            concordant.append(gene)
        elif directions[0] != directions[1] and any(float(row["spearman_q"]) < 0.05 for row in rows):
            discordant.append(gene)

    def values_for(gene: str) -> str:
        rows = by_gene[gene]
        return "; ".join(
            f"{row['cohort']} ρ={float(row['spearman_rho']):.2f}, "
            f"high/low log2 ratio={float(row['log2_median_ratio']):+.2f}, "
            f"q={float(row['spearman_q']):.2g}"
            for row in rows
        )

    ranked = sorted(
        concordant,
        key=lambda gene: sum(abs(float(row["spearman_rho"])) for row in by_gene[gene]),
        reverse=True,
    )
    findings = "\n".join(f"- **{gene}:** {values_for(gene)}" for gene in ranked)
    if not findings:
        findings = "- No galectin was significant in the same direction in both histologies."

    summary = f"""# A11 — galectin genes vs TACSTD2-high public lung

## Honest result (≤200 words)
Public TCGA primary tumors were analyzed separately: LUAD n={cohort_counts['LUAD']} and LUSC n={cohort_counts['LUSC']}. TACSTD2-high means the upper within-cohort half; it is not a validated biological subtype.

{findings}

Discordant significant genes: {", ".join(discordant) or "none"}. These are associations, not evidence that TROP2 regulates galectins. Bulk RNA mixes malignant, stromal, and immune cells; TACSTD2 is epithelial, whereas several galectins are also abundant in non-malignant compartments. Histology-stratified agreement is therefore the strongest claim supported here. Protein abundance, spatial co-expression, treatment response, and single-cell malignant-cell effects were not tested.

## Methods
Values are cBioPortal PanCancer Atlas batch-normalized RNA-seq RSEM. Spearman correlation is primary; high/low median ratios and Mann–Whitney tests are descriptive. BH correction covers 28 gene-by-cohort tests. Exact values and the auditable sample-level extract are provided.

## Sources
- [LUAD study](https://www.cbioportal.org/study/summary?id={STUDIES['LUAD']})
- [LUSC study](https://www.cbioportal.org/study/summary?id={STUDIES['LUSC']})
- [cBioPortal API](https://www.cbioportal.org/api/swagger-ui/index.html)
"""
    (OUT / "README.md").write_text(summary)

    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "api": API,
        "studies": STUDIES,
        "molecular_profile_suffix": "_rna_seq_v2_mrna",
        "measurement": "Batch-normalized RNA Seq V2 RSEM",
        "genes": GENES,
        "group_definition": "Within-cohort rank split; equal lower and upper halves",
        "multiple_testing": "Benjamini-Hochberg over 28 gene-by-cohort tests per test family",
        "software": {
            "numpy": np.__version__,
            "requests": requests.__version__,
        },
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
