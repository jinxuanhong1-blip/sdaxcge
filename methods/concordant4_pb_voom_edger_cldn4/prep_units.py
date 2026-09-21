#!/usr/bin/env python3
"""Lock concordant-4 units and gene families for the pseudo-bulk DE.

Quartiles are the PR #503 within-cohort malignant CLDN4 %pos split
(rank, then four bins). They are not re-cut on the DE subset.
P4001 stays Q1 on that 22-patient vector and is absent from the
malignant UMI-sum, so it is not a DE unit.
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"

COHORTS = ["GSE123902", "GSE131907", "GSE205335", "GSE189357"]
TUMOR_ORIGINS = {"tLung", "tL/B", "mLN", "PE", "mBrain"}

# PR #503 n_honest.tsv Q1/Q4 labels. Fail if the quantile rule drifts.
LOCKED = {
    "GSE123902": {
        "Q1": ["LX675", "LX682", "LX699", "LX701"],
        "Q4": ["LX653", "LX680", "LX684"],
    },
    "GSE131907": {
        "Q1": ["EBUS_13", "EBUS_15", "EBUS_49", "NS_02", "NS_06", "NS_16"],
        "Q4": ["EBUS_19", "EBUS_28", "NS_03", "NS_04", "NS_07"],
    },
    "GSE205335": {
        "Q1": ["P1015", "P1062", "P1063", "P1090", "P1119", "P4001"],
        "Q4": ["P1016", "P1025", "P1037", "P1084", "P1089", "P1115"],
    },
    "GSE189357": {
        "Q1": ["TD2", "TD4", "TD7"],
        "Q4": ["TD6", "TD9"],
    },
}

CHEMOKINE = [
    "CXCL9", "CXCL10", "CXCL11", "CXCL13", "CXCL16", "CXCL8", "CX3CL1",
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL19", "CCL21", "CCL22",
    "XCL1", "XCL2", "IL15", "IL2", "IL7", "IL21",
    "CXCR3", "CXCR4", "CXCR5", "CCR5", "CCR7",
]


def assign_quartiles(values: pd.Series) -> pd.Series:
    ranks = values.astype(float).rank(method="average")
    qs = pd.qcut(ranks, 4, labels=["Q1", "Q2", "Q3", "Q4"], duplicates="drop")
    return pd.Series(qs.astype(str), index=values.index)


def _pct_percent(pct: pd.Series) -> pd.Series:
    """GSE123902 and GSE189357 store a fraction; the other two store percent."""
    x = pct.astype(float)
    if float(x.max()) <= 1.5:
        return x * 100.0
    return x


def matrix_ids(cohort: str) -> set[str]:
    cols = pd.read_csv(
        DATA / f"{cohort}_malignant_counts.tsv.gz", sep="\t", nrows=0
    ).columns.astype(str)
    return set(cols.tolist()) - {""}


def load_units() -> pd.DataFrame:
    rows = []

    d = pd.read_csv(DATA / "GSE123902_marker_units.tsv", sep="\t")
    tumor = d[d["tissue"].isin(["PRIMARY", "METASTASIS"])].copy()
    tumor = tumor.sort_values(["patient", "tissue"]).drop_duplicates("patient", keep="first")
    el = tumor[tumor["eligible"].astype(str).str.lower() == "true"].copy()
    meta = pd.read_csv(DATA / "GSE123902_malignant_meta.tsv", sep="\t")
    nmap = meta.set_index("patient")["n_malignant_summed"].to_dict()
    q = assign_quartiles(pd.Series(_pct_percent(el["mal_CLDN4_pct"]).to_numpy(), index=el["patient"].astype(str)))
    for r in el.itertuples(index=False):
        pid = str(r.patient)
        rows.append(
            {
                "patient": pid,
                "cohort": "GSE123902",
                "unit": "donor",
                "malig_def": "marker_malig",
                "cldn4_pct": float(_pct_percent(pd.Series([r.mal_CLDN4_pct])).iloc[0]),
                "quartile": q.loc[pid],
                "n_malignant": int(nmap.get(pid, r.n_malignant)),
                "tissue": r.tissue,
                "origin": "",
                "cancer_subtype": "",
                "recist": "",
            }
        )

    d = pd.read_csv(DATA / "GSE131907_samples.tsv", sep="\t")
    mal = d[d["origin"].isin(TUMOR_ORIGINS) & (d["n_malignant"] >= 20)].copy()
    q = assign_quartiles(pd.Series(mal["mal_CLDN4_pct"].to_numpy(float), index=mal["sample"].astype(str)))
    for r in mal.itertuples(index=False):
        pid = str(r.sample)
        rows.append(
            {
                "patient": pid,
                "cohort": "GSE131907",
                "unit": "sample",
                "malig_def": "author_malig",
                "cldn4_pct": float(r.mal_CLDN4_pct),
                "quartile": q.loc[pid],
                "n_malignant": int(r.n_malignant),
                "tissue": "",
                "origin": r.origin,
                "cancer_subtype": "LUAD",
                "recist": "",
            }
        )

    d = pd.read_csv(DATA / "GSE205335_patients.tsv", sep="\t")
    q = assign_quartiles(
        pd.Series(d["mal_CLDN4_pct_pos"].to_numpy(float), index=d["patient"].astype(str))
    )
    for r in d.itertuples(index=False):
        pid = str(r.patient)
        rows.append(
            {
                "patient": pid,
                "cohort": "GSE205335",
                "unit": "patient",
                "malig_def": "author_malig",
                "cldn4_pct": float(r.mal_CLDN4_pct_pos),
                "quartile": q.loc[pid],
                "n_malignant": int(r.n_malignant),
                "tissue": r.tissue,
                "origin": "",
                "cancer_subtype": r.cancer_subtype,
                "recist": r.recist,
            }
        )

    d = pd.read_csv(DATA / "GSE189357_marker_units.tsv", sep="\t")
    el = d[d["eligible"].astype(str).str.lower() == "true"].copy()
    meta = pd.read_csv(DATA / "GSE189357_malignant_meta.tsv", sep="\t")
    nmap = meta.set_index("patient")["n_malignant_summed"].to_dict()
    q = assign_quartiles(pd.Series(_pct_percent(el["mal_CLDN4_pct"]).to_numpy(), index=el["patient"].astype(str)))
    for r in el.itertuples(index=False):
        pid = str(r.patient)
        rows.append(
            {
                "patient": pid,
                "cohort": "GSE189357",
                "unit": "patient",
                "malig_def": "marker_malig",
                "cldn4_pct": float(_pct_percent(pd.Series([r.mal_CLDN4_pct])).iloc[0]),
                "quartile": q.loc[pid],
                "n_malignant": int(nmap.get(pid, r.n_malignant)),
                "tissue": r.tissue,
                "origin": "",
                "cancer_subtype": "AIS-IAC",
                "recist": "",
            }
        )

    out = pd.DataFrame(rows)
    in_mat = {c: matrix_ids(c) for c in COHORTS}
    out["in_count_matrix"] = [
        pid in in_mat[c] for pid, c in zip(out["patient"], out["cohort"])
    ]
    return out


def check_locked(units: pd.DataFrame) -> None:
    for cohort, arms in LOCKED.items():
        sub = units[units["cohort"] == cohort]
        for arm, ids in arms.items():
            got = sorted(sub.loc[sub["quartile"] == arm, "patient"])
            exp = sorted(ids)
            if got != exp:
                raise SystemExit(f"{cohort} {arm} drifted:\n got {got}\n exp {exp}")
    # P4001 is the only locked Q-tail unit missing from a count matrix.
    tails = units[units["quartile"].isin(["Q1", "Q4"]) & ~units["in_count_matrix"]]
    missing = sorted(tails["patient"])
    if missing != ["P4001"]:
        raise SystemExit(f"unexpected units missing from count matrices: {missing}")


def write_families() -> None:
    a8 = json.loads((DATA / "a8_sets.json").read_text())
    sets = a8["sets"]
    krt = set(sets["KRT_EPITHELIAL"])
    tj = set(sets["KEGG_TIGHT_JUNCTION"]) | set(sets["GOBP_TIGHT_JUNCTION_ORGANIZATION"])
    tj |= {g for g in a8["focal_genes"] if g not in krt}
    tj.discard("CLDN4")
    fam = {
        "IFN": set(sets["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
        | set(sets["HALLMARK_INTERFERON_ALPHA_RESPONSE"]),
        "MHC-I/APM": set(sets["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"]),
        "chemokine": set(CHEMOKINE),
        "TJ": tj,
        "keratin": set(krt),
    }
    if "CLDN4" in fam["TJ"]:
        raise SystemExit("CLDN4 leaked into TJ")
    rows = []
    for name, genes in fam.items():
        for g in sorted(genes):
            rows.append({"family": name, "gene": g})
    pd.DataFrame(rows).to_csv(DATA / "families.tsv", sep="\t", index=False)
    print({k: len(v) for k, v in fam.items()})


def main() -> None:
    units = load_units()
    check_locked(units)
    units.to_csv(DATA / "units.tsv", sep="\t", index=False)
    write_families()
    de = units[units["quartile"].isin(["Q1", "Q4"]) & units["in_count_matrix"]]
    print("quartile-vector n", units.groupby("cohort").size().to_dict(), "total", len(units))
    print("DE n", de.groupby(["cohort", "quartile"]).size().to_dict(), "total", len(de))
    sub = units[units["cohort"] == "GSE205335"]
    tails = sub[sub["quartile"].isin(["Q1", "Q4"])]
    print("GSE205335 Q tails by subtype")
    print(pd.crosstab(tails["quartile"], tails["cancer_subtype"]))
    print("GSE205335 Q tails in the count matrix")
    inn = tails[tails["in_count_matrix"]]
    print(pd.crosstab(inn["quartile"], inn["cancer_subtype"]))


if __name__ == "__main__":
    main()
