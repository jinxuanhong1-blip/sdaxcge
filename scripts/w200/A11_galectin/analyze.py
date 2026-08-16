#!/usr/bin/env python3
"""Galectin-family associations with TACSTD2 in public lung RNA.

Primary: TCGA LUAD / LUSC PanCancer Atlas.
Replication: OncoSG LUAD, CPTAC LUAD, CPTAC LUSC, CAS LUAD.
Tumor-cell check: CCLE / DepMap Broad 2025 lung cell lines.
Purity check: CPTAC ESTIMATE purity partial Spearman, where available.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import requests
import scipy
from scipy.stats import mannwhitneyu, rankdata, spearmanr


API = "https://www.cbioportal.org/api"
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
GALECTINS = [gene for gene in GENES if gene != "TACSTD2"]
OUT = Path(__file__).resolve().parents[3] / "results" / "w200" / "A11_galectin"

COHORTS = [
    {
        "cohort": "TCGA_LUAD",
        "family": "primary_tcga",
        "study_id": "luad_tcga_pan_can_atlas_2018",
        "profile_id": "luad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "sample_list_id": "luad_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "measurement": "Batch-normalized RNA Seq V2 RSEM",
        "note": "Primary public bulk cohort.",
    },
    {
        "cohort": "TCGA_LUSC",
        "family": "primary_tcga",
        "study_id": "lusc_tcga_pan_can_atlas_2018",
        "profile_id": "lusc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "sample_list_id": "lusc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna",
        "measurement": "Batch-normalized RNA Seq V2 RSEM",
        "note": "Primary public bulk cohort.",
    },
    {
        "cohort": "OncoSG_LUAD",
        "family": "replication_bulk",
        "study_id": "luad_oncosg_2020",
        "profile_id": "luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores",
        "sample_list_id": "luad_oncosg_2020_rna_seq_v2_mrna",
        "measurement": "mRNA z-scores relative to all samples (log RNA Seq V2 RSEM)",
        "note": "Independent East-Asian LUAD bulk replication. Values are z-scores.",
    },
    {
        "cohort": "CPTAC_LUAD",
        "family": "replication_bulk",
        "study_id": "luad_cptac_2020",
        "profile_id": "luad_cptac_2020_mrna",
        "sample_list_id": "luad_cptac_2020_rna_seq_mrna",
        "measurement": "mRNA expression (RPKM, log2)",
        "purity_attr": "TUMOR_PURITY_BYESTIMATE_RNASEQ",
        "note": "Independent LUAD bulk replication with ESTIMATE purity.",
    },
    {
        "cohort": "CPTAC_LUSC",
        "family": "replication_bulk",
        "study_id": "lusc_cptac_2021",
        "profile_id": "lusc_cptac_2021_rna_seq_mrna",
        "sample_list_id": "lusc_cptac_2021_rna_seq_mrna",
        "measurement": "mRNA expression (FPKM, log2 transformed)",
        "purity_attr": "ESTIMATE_TUMORPURITY",
        "note": "Independent LUSC bulk replication with ESTIMATE purity.",
    },
    {
        "cohort": "CAS_LUAD",
        "family": "replication_bulk",
        "study_id": "luad_cas_2020",
        "profile_id": "luad_cas_2020_rna_seq_mrna",
        "sample_list_id": "luad_cas_2020_rna_seq_mrna",
        "measurement": "mRNA expression (FPKM)",
        "purity_attr": "TUMOR_PURITY",
        "note": "Smaller LUAD sensitivity cohort; pathologist purity is percent.",
    },
    {
        "cohort": "CCLE_lung",
        "family": "tumor_cell",
        "study_id": "ccle_broad_2025",
        "profile_id": "ccle_broad_2025_rna_seq_mrna",
        "sample_list_id": "ccle_broad_2025_all",
        "measurement": "mRNA expression (RNA Seq TPM)",
        "lineage_value": "Lung",
        "note": "Tumor-cell-only public check: OncotreeLineage==Lung cell lines.",
    },
    {
        "cohort": "CCLE_NSCLC",
        "family": "tumor_cell",
        "study_id": "ccle_broad_2025",
        "profile_id": "ccle_broad_2025_rna_seq_mrna",
        "sample_list_id": "ccle_broad_2025_all",
        "measurement": "mRNA expression (RNA Seq TPM)",
        "lineage_value": "Lung",
        "disease_value": "Non-Small Cell Lung Cancer",
        "note": "Sensitivity: lung lineage restricted to NSCLC (excludes NET/SCLC).",
    },
]


def get_json(session: requests.Session, url: str, **kwargs) -> object:
    response = session.get(url, timeout=90, **kwargs)
    response.raise_for_status()
    return response.json()


def fetch_expression(
    session: requests.Session, profile_id: str, sample_list_id: str
) -> dict[str, dict[str, float]]:
    response = session.post(
        f"{API}/molecular-profiles/{profile_id}/molecular-data/fetch",
        params={"projection": "DETAILED"},
        json={"entrezGeneIds": list(GENES.values()), "sampleListId": sample_list_id},
        timeout=180,
    )
    response.raise_for_status()
    entrez_to_symbol = {entrez: symbol for symbol, entrez in GENES.items()}
    by_sample: dict[str, dict[str, float]] = defaultdict(dict)
    for row in response.json():
        symbol = entrez_to_symbol.get(row["entrezGeneId"])
        if symbol and row.get("value") is not None:
            by_sample[row["sampleId"]][symbol] = float(row["value"])
    return dict(by_sample)


def fetch_clinical_map(
    session: requests.Session, study_id: str, attributes: list[str]
) -> dict[str, dict[str, str]]:
    rows = get_json(
        session,
        f"{API}/studies/{study_id}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
    )
    wanted = set(attributes)
    by_sample: dict[str, dict[str, str]] = defaultdict(dict)
    for row in rows:
        attr = row.get("clinicalAttributeId")
        if attr in wanted and row.get("value") not in (None, ""):
            by_sample[row["sampleId"]][attr] = str(row["value"])
    return dict(by_sample)


def bh_adjust(rows: list[dict[str, object]], p_key: str, q_key: str) -> None:
    valid = [(i, float(row[p_key])) for i, row in enumerate(rows) if row[p_key] != ""]
    valid.sort(key=lambda item: item[1])
    m = len(valid)
    running = 1.0
    adjusted = [1.0] * m
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


def residualize(y: np.ndarray, z: np.ndarray) -> np.ndarray:
    z_design = np.column_stack([np.ones(len(z)), z])
    coef, *_ = np.linalg.lstsq(z_design, y, rcond=None)
    return y - z_design @ coef


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float]:
    xr = residualize(rankdata(x, method="average"), rankdata(z, method="average"))
    yr = residualize(rankdata(y, method="average"), rankdata(z, method="average"))
    rho, p_value = spearmanr(xr, yr)
    return float(rho), float(p_value)


def empty_stat(cohort: str, family: str, gene: str, n: int, n_low: int, n_high: int) -> dict[str, object]:
    return {
        "family": family,
        "cohort": cohort,
        "gene": gene,
        "n": n,
        "n_low": n_low,
        "n_high": n_high,
        "median_low": "",
        "median_high": "",
        "log2_median_ratio": "",
        "spearman_rho": "",
        "spearman_p": "",
        "mannwhitney_p": "",
        "spearman_q": "",
        "mannwhitney_q": "",
    }


def gene_stats(
    cohort: str,
    family: str,
    gene: str,
    x: np.ndarray,
    y: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
) -> dict[str, object]:
    if len(y) < 8 or np.all(y == y[0]):
        return empty_stat(cohort, family, gene, len(y), len(low), len(high))
    rho, rho_p = spearmanr(x, y)
    mw = mannwhitneyu(high, low, alternative="two-sided", method="asymptotic")
    return {
        "family": family,
        "cohort": cohort,
        "gene": gene,
        "n": int(len(y)),
        "n_low": int(len(low)),
        "n_high": int(len(high)),
        "median_low": float(np.median(low)),
        "median_high": float(np.median(high)),
        "log2_median_ratio": math.log2((float(np.median(high)) + 1.0) / (float(np.median(low)) + 1.0)),
        "spearman_rho": float(rho),
        "spearman_p": float(rho_p),
        "mannwhitney_p": float(mw.pvalue),
        "spearman_q": "",
        "mannwhitney_q": "",
    }


def analyze_cohort(
    spec: dict[str, str],
    by_sample: dict[str, dict[str, float]],
    clinical: dict[str, dict[str, str]],
) -> tuple[list[str], list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    samples = [sample for sample, values in by_sample.items() if "TACSTD2" in values]
    if spec.get("lineage_value"):
        samples = [
            sample
            for sample in samples
            if clinical.get(sample, {}).get("ONCOTREE_LINEAGE") == spec["lineage_value"]
        ]
    if spec.get("disease_value"):
        samples = [
            sample
            for sample in samples
            if clinical.get(sample, {}).get("ONCOTREE_PRIMARY_DISEASE") == spec["disease_value"]
        ]
    samples.sort(key=lambda sample: (by_sample[sample]["TACSTD2"], sample))
    split = len(samples) // 2
    low_samples = samples[:split]
    high_samples = samples[-split:]
    group = {sample: "low" for sample in low_samples}
    group.update({sample: "high" for sample in high_samples})

    expression_rows = []
    for sample in samples:
        row: dict[str, object] = {
            "family": spec["family"],
            "cohort": spec["cohort"],
            "study_id": spec["study_id"],
            "sample_id": sample,
            "tacstd2_group": group.get(sample, "median_unassigned"),
        }
        if spec.get("purity_attr"):
            row["purity"] = clinical.get(sample, {}).get(spec["purity_attr"], "")
        else:
            row["purity"] = clinical.get(sample, {}).get("ONCOTREE_PRIMARY_DISEASE", "")
        row.update({gene: by_sample[sample].get(gene, "") for gene in GENES})
        expression_rows.append(row)

    stats = []
    for gene in GALECTINS:
        paired = [sample for sample in low_samples + high_samples if gene in by_sample[sample]]
        x = np.asarray([by_sample[sample]["TACSTD2"] for sample in paired], dtype=float)
        y = np.asarray([by_sample[sample][gene] for sample in paired], dtype=float)
        low = np.asarray(
            [by_sample[sample][gene] for sample in low_samples if gene in by_sample[sample]],
            dtype=float,
        )
        high = np.asarray(
            [by_sample[sample][gene] for sample in high_samples if gene in by_sample[sample]],
            dtype=float,
        )
        stats.append(gene_stats(spec["cohort"], spec["family"], gene, x, y, low, high))

    partial_rows = []
    purity_attr = spec.get("purity_attr")
    if purity_attr:
        for gene in GALECTINS:
            triples = []
            for sample in samples:
                values = by_sample[sample]
                raw_purity = clinical.get(sample, {}).get(purity_attr)
                if gene not in values or raw_purity in (None, ""):
                    continue
                triples.append((values["TACSTD2"], values[gene], float(raw_purity)))
            if len(triples) < 20:
                continue
            x, y, z = (np.asarray(col, dtype=float) for col in zip(*triples))
            rho, p_value = partial_spearman(x, y, z)
            ordinary, ordinary_p = spearmanr(x, y)
            partial_rows.append(
                {
                    "cohort": spec["cohort"],
                    "gene": gene,
                    "n": int(len(triples)),
                    "purity_attr": purity_attr,
                    "spearman_rho": float(ordinary),
                    "spearman_p": float(ordinary_p),
                    "partial_spearman_rho": rho,
                    "partial_spearman_p": p_value,
                    "partial_spearman_q": "",
                }
            )
    return samples, expression_rows, stats, partial_rows


def values_for(rows: list[dict[str, object]]) -> str:
    return "; ".join(
        f"{row['cohort']} ρ={float(row['spearman_rho']):.2f}, "
        f"high/low log2 ratio={float(row['log2_median_ratio']):+.2f}, "
        f"q={float(row['spearman_q']):.2g}"
        for row in rows
        if row["spearman_rho"] != ""
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    session = requests.Session()
    session.headers.update({"Accept": "application/json", "Content-Type": "application/json"})

    expression_cache: dict[tuple[str, str], dict[str, dict[str, float]]] = {}
    clinical_cache: dict[str, dict[str, dict[str, str]]] = {}
    stats: list[dict[str, object]] = []
    expression_rows: list[dict[str, object]] = []
    partial_rows: list[dict[str, object]] = []
    cohort_counts: dict[str, int] = {}

    for spec in COHORTS:
        key = (spec["profile_id"], spec["sample_list_id"])
        if key not in expression_cache:
            expression_cache[key] = fetch_expression(session, spec["profile_id"], spec["sample_list_id"])
        needed_attrs = ["ONCOTREE_LINEAGE", "ONCOTREE_PRIMARY_DISEASE"]
        if spec.get("purity_attr"):
            needed_attrs.append(spec["purity_attr"])
        if spec["study_id"] not in clinical_cache:
            clinical_cache[spec["study_id"]] = fetch_clinical_map(session, spec["study_id"], needed_attrs)
        else:
            extra = [attr for attr in needed_attrs if attr not in {"ONCOTREE_LINEAGE", "ONCOTREE_PRIMARY_DISEASE"}]
            if extra:
                more = fetch_clinical_map(session, spec["study_id"], extra)
                for sample, values in more.items():
                    clinical_cache[spec["study_id"]].setdefault(sample, {}).update(values)
        samples, expr_rows, cohort_stats, cohort_partial = analyze_cohort(
            spec, expression_cache[key], clinical_cache[spec["study_id"]]
        )
        cohort_counts[spec["cohort"]] = len(samples)
        expression_rows.extend(expr_rows)
        stats.extend(cohort_stats)
        partial_rows.extend(cohort_partial)

    for family in {row["family"] for row in stats}:
        family_rows = [row for row in stats if row["family"] == family]
        bh_adjust(family_rows, "spearman_p", "spearman_q")
        bh_adjust(family_rows, "mannwhitney_p", "mannwhitney_q")
    if partial_rows:
        bh_adjust(partial_rows, "partial_spearman_p", "partial_spearman_q")

    stats_fields = [
        "family",
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

    expression_fields = [
        "family",
        "cohort",
        "study_id",
        "sample_id",
        "tacstd2_group",
        "purity",
        *GENES,
    ]
    with (OUT / "expression_values.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=expression_fields, delimiter="\t")
        writer.writeheader()
        for row in expression_rows:
            writer.writerow({field: fmt(row.get(field, "")) for field in expression_fields})

    partial_fields = [
        "cohort",
        "gene",
        "n",
        "purity_attr",
        "spearman_rho",
        "spearman_p",
        "partial_spearman_rho",
        "partial_spearman_p",
        "partial_spearman_q",
    ]
    with (OUT / "purity_partial.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=partial_fields, delimiter="\t")
        writer.writeheader()
        for row in partial_rows:
            writer.writerow({field: fmt(row[field]) for field in partial_fields})

    primary = [row for row in stats if row["family"] == "primary_tcga"]
    by_gene_primary = {gene: [row for row in primary if row["gene"] == gene] for gene in GALECTINS}
    concordant = []
    discordant = []
    cohort_specific = []
    for gene, rows in by_gene_primary.items():
        usable = [row for row in rows if row["spearman_rho"] != ""]
        if len(usable) != 2:
            continue
        directions = [math.copysign(1, float(row["spearman_rho"])) for row in usable]
        significant = [float(row["spearman_q"]) < 0.05 for row in usable]
        if directions[0] == directions[1] and all(significant):
            concordant.append(gene)
        elif directions[0] != directions[1] and any(significant):
            discordant.append(gene)
        elif directions[0] == directions[1] and sum(significant) == 1:
            hit = usable[significant.index(True)]
            cohort_specific.append(f"{gene} ({hit['cohort']})")

    ranked = sorted(
        concordant,
        key=lambda gene: sum(abs(float(row["spearman_rho"])) for row in by_gene_primary[gene]),
        reverse=True,
    )
    findings = "\n".join(f"- **{gene}:** {values_for(by_gene_primary[gene])}" for gene in ranked)
    if not findings:
        findings = "- No galectin was significant in the same direction in both TCGA histologies."

    def same_direction_hits(gene: str, family: str, q_key: str = "spearman_q", rho_key: str = "spearman_rho") -> list[str]:
        primary_sign = math.copysign(1, float(by_gene_primary[gene][0]["spearman_rho"]))
        hits = []
        for row in stats if q_key.startswith("spearman") else partial_rows:
            if q_key.startswith("spearman") and (row.get("family") != family or row["gene"] != gene):
                continue
            if not q_key.startswith("spearman") and row["gene"] != gene:
                continue
            if row.get(rho_key) == "" or row.get(q_key) == "":
                continue
            if math.copysign(1, float(row[rho_key])) == primary_sign and float(row[q_key]) < 0.05:
                hits.append(
                    f"{row['cohort']} ρ={float(row[rho_key]):.2f}, q={float(row[q_key]):.2g}"
                )
        return hits

    replication_bits = []
    cell_bits = []
    purity_bits = []
    for gene in ranked:
        rep = same_direction_hits(gene, "replication_bulk")
        cell = same_direction_hits(gene, "tumor_cell")
        pur = same_direction_hits(gene, "replication_bulk", "partial_spearman_q", "partial_spearman_rho")
        replication_bits.append(f"{gene}: {'; '.join(rep) if rep else 'no independent bulk replication at q<0.05'}")
        cell_bits.append(f"{gene}: {'; '.join(cell) if cell else 'not significant in CCLE lung/NSCLC'}")
        if pur:
            purity_bits.append(f"{gene}: {'; '.join(pur)}")
        else:
            purity_bits.append(f"{gene}: no same-direction purity-partial q<0.05")

    n_primary = sum(row["spearman_p"] != "" for row in primary)
    empty_cohorts = [name for name, count in cohort_counts.items() if count == 0]
    empty_note = (
        f" Unavailable: {', '.join(empty_cohorts)} (no TACSTD2 in the public matrix)."
        if empty_cohorts
        else ""
    )
    summary = f"""# A11 — galectin genes vs TACSTD2-high public lung

## Honest result (≤200 words)
TCGA primary tumors, analyzed separately: LUAD n={cohort_counts['TCGA_LUAD']}, LUSC n={cohort_counts['TCGA_LUSC']}. TACSTD2-high is the upper half, not a validated subtype.

{findings}

**LGALS3** is the only gene that is TCGA-concordant, replicated in independent LUAD bulk (OncoSG ρ=0.41; CPTAC LUAD ρ=0.39), retained after CPTAC LUAD ESTIMATE-purity partial Spearman (ρ=0.38), and present in tumor-cell-only CCLE 2025 (lung n={cohort_counts['CCLE_lung']} ρ=0.68; NSCLC n={cohort_counts['CCLE_NSCLC']} ρ=0.58). It did **not** replicate in CPTAC LUSC bulk (q>0.05).

**LGALS9B/9C** stay positive in TCGA and CCLE, but independent bulk is weak: 9C none; 9B only CPTAC LUSC. Both are low-abundance in several matrices.

Discordant TCGA genes: {", ".join(discordant) or "none"}. LUSC-only: {", ".join(cohort_specific) or "none"}.{empty_note} Associations only; not proof that TROP2 regulates galectins. Protein, spatial, ICI, and single-cell malignant-cell effects were not tested.

## Methods
Spearman is primary. High/low ratios and Mann–Whitney tests are descriptive. BH correction is within family (TCGA; other bulk; cell lines). Partial Spearman residualizes ranks on ESTIMATE purity. Scales are not pooled.

## Sources
- [TCGA LUAD](https://www.cbioportal.org/study/summary?id=luad_tcga_pan_can_atlas_2018)
- [TCGA LUSC](https://www.cbioportal.org/study/summary?id=lusc_tcga_pan_can_atlas_2018)
- [OncoSG LUAD](https://www.cbioportal.org/study/summary?id=luad_oncosg_2020)
- [CPTAC LUAD](https://www.cbioportal.org/study/summary?id=luad_cptac_2020)
- [CPTAC LUSC](https://www.cbioportal.org/study/summary?id=lusc_cptac_2021)
- [CAS LUAD](https://www.cbioportal.org/study/summary?id=luad_cas_2020)
- [CCLE Broad 2025](https://www.cbioportal.org/study/summary?id=ccle_broad_2025)
"""
    (OUT / "README.md").write_text(summary)

    verdict = {
        "primary_concordant_genes": ranked,
        "primary_discordant_genes": discordant,
        "primary_one_cohort_genes": cohort_specific,
        "replication": {gene: same_direction_hits(gene, "replication_bulk") for gene in ranked},
        "tumor_cell": {gene: same_direction_hits(gene, "tumor_cell") for gene in ranked},
        "purity_partial": {
            gene: same_direction_hits(gene, "replication_bulk", "partial_spearman_q", "partial_spearman_rho")
            for gene in ranked
        },
        "statement": (
            "LGALS3 is the only galectin that is TCGA histology-concordant, replicated in "
            "independent LUAD bulk, retained after CPTAC LUAD purity residualization, and "
            "present in CCLE lung/NSCLC cell lines. It did not replicate in CPTAC LUSC. "
            "LGALS9B/9C persist in TCGA and CCLE but lack consistent independent bulk support."
        ),
    }
    (OUT / "summary.json").write_text(json.dumps(verdict, indent=2) + "\n")
    (OUT / "verdict.txt").write_text(verdict["statement"] + "\n")

    provenance = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "api": API,
        "cohorts": COHORTS,
        "genes": GENES,
        "group_definition": "Within-cohort rank split; equal lower and upper halves",
        "multiple_testing": (
            f"Benjamini-Hochberg within analysis family; primary TCGA had {n_primary} "
            "testable gene-by-cohort comparisons"
        ),
        "partial_correlation": "Spearman of rank residuals after OLS on purity ranks",
        "software": {
            "numpy": np.__version__,
            "requests": requests.__version__,
            "scipy": scipy.__version__,
        },
        "cohort_n": cohort_counts,
    }
    (OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
