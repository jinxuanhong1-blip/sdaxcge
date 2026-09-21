#!/usr/bin/env python3
"""Download TCPA RPPA500 disease-specific level-4 abundances.

Source: TCPA patient-cohort CouchDB (RPPA500, release 5), public read.
Disease L4 is the single-histology matrix (TCPA FAQ). The pan-cancer zip is
not the analysis matrix.

LUAD and LUSC are kept as separate datasets. No RNA stand-in.
"""
from __future__ import annotations

import argparse
import json
import urllib.parse
import urllib.request
from pathlib import Path

BASE = "https://tcpa.drbioright.org/rppa500/"

COHORTS = {
    "LUAD": "TCGA-LUAD-L4",
    "LUSC": "TCGA-LUSC-L4",
}

# Antibodies to pull when the protein_id exists. Absence is recorded, not filled.
FETCH_IDS = [
    "CLAUDIN7",
    "ECADHERIN",
    "EMA",
    "KU80",
    "HLADQA1",
    "NCADHERIN",
    "PCADHERIN",
    "ALPHACATENIN",
    "P63",
]


def get_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode("utf-8"))


def protein_ids(dataset_id: str) -> list[str]:
    start = urllib.parse.quote(f'"{dataset_id}-"')
    end = urllib.parse.quote(f'"{dataset_id}-\ufff0"')
    url = f"{BASE}_all_docs?startkey={start}&endkey={end}&limit=2000"
    doc = get_json(url)
    prefix = dataset_id + "-"
    ids = []
    for row in doc["rows"]:
        i = row["id"]
        if i.startswith(prefix):
            ids.append(i[len(prefix) :])
    return ids


def fetch_protein(dataset_id: str, protein_id: str) -> dict:
    url = BASE + urllib.parse.quote(f"{dataset_id}-{protein_id}")
    return get_json(url)


def pancan_missingness(tsv: Path, case_ids: dict[str, set[str]], proteins: list[str]) -> list[dict]:
    """Count non-empty cells in the pan-cancer RPPA500 matrix for the disease-L4 cases.

    Used only to confirm that an all-NA disease-L4 antibody is also empty in the
    pan-cancer file for the same case IDs. Not an analysis matrix.
    """
    with tsv.open() as handle:
        header = handle.readline().rstrip("\n").split("\t")
    idx = {name: i for i, name in enumerate(header)}
    missing = [p for p in proteins if p not in idx]
    if missing:
        raise SystemExit(f"pan-can header missing {missing}")
    counts = {cohort: {p: [0, 0] for p in proteins} for cohort in list(case_ids) + ["other"]}
    with tsv.open() as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            sid = parts[0][:12]
            cohort = "other"
            for name, ids in case_ids.items():
                if sid in ids:
                    cohort = name
                    break
            for protein in proteins:
                raw = parts[idx[protein]]
                empty = raw in ("", "NA", "NaN", "null", "None")
                counts[cohort][protein][1 if empty else 0] += 1
    rows = []
    for cohort, by_protein in counts.items():
        for protein, (n_obs, n_empty) in by_protein.items():
            rows.append(
                {
                    "cohort": cohort,
                    "protein_id": protein,
                    "n_observed": n_obs,
                    "n_empty": n_empty,
                }
            )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/tcpa_rppa_epithelial_ku_hla")
    ap.add_argument(
        "--pancan-tsv",
        default="",
        help="Optional TCPA_TCGA_RPPA500.tsv for a missingness cross-check",
    )
    args = ap.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    inventory = {}
    rows = []
    stored = []
    meta = {
        "source": BASE,
        "release_note": "TCPA RPPA500 patient cohorts, tcpa_release field on each protein doc",
        "normalization": "disease-specific level 4, replicates-based normalization (RBN)",
        "cohorts": {},
    }

    for cohort, dataset_id in COHORTS.items():
        ids = protein_ids(dataset_id)
        inventory[cohort] = ids
        (out / f"protein_ids_{cohort}.txt").write_text("\n".join(ids) + "\n")
        doc_meta = get_json(BASE + urllib.parse.quote(f"DOC-{dataset_id}"))
        documentation = doc_meta.get("documentation", {})
        meta["cohorts"][cohort] = {
            "dataset_id": dataset_id,
            "n_protein_docs": len(ids),
            "documentation": documentation,
        }
        present = [p for p in FETCH_IDS if p in ids]
        for protein_id in present:
            doc = fetch_protein(dataset_id, protein_id)
            abundances = doc.get("abundances") or []
            n_null = sum(1 for a in abundances if a.get("abundance") is None)
            meta["cohorts"][cohort].setdefault("fetched", {})[protein_id] = {
                "n_abundance_rows": len(abundances),
                "n_null": n_null,
                "tcpa_release": doc.get("tcpa_release"),
                "data_set": doc.get("data_set"),
            }
            for a in abundances:
                rows.append(
                    {
                        "cohort": cohort,
                        "dataset_id": dataset_id,
                        "sample_id": a.get("sample_id"),
                        "protein_id": protein_id,
                        "abundance": a.get("abundance"),
                    }
                )
            for item in (doc.get("correlations") or {}).get("spearman") or []:
                other = item.get("protein_id")
                if other in FETCH_IDS and other != protein_id:
                    stored.append(
                        {
                            "cohort": cohort,
                            "protein_id": protein_id,
                            "other_protein_id": other,
                            "tcpa_stored_rho": item.get("correlation"),
                            "tcpa_stored_p": item.get("p_value"),
                        }
                    )

    # long table
    header = ["cohort", "dataset_id", "sample_id", "protein_id", "abundance"]
    lines = ["\t".join(header)]
    for rec in rows:
        val = "" if rec["abundance"] is None else str(rec["abundance"])
        lines.append(
            "\t".join(
                [
                    rec["cohort"],
                    rec["dataset_id"],
                    str(rec["sample_id"]),
                    rec["protein_id"],
                    val,
                ]
            )
        )
    (out / "abundances_long.tsv").write_text("\n".join(lines) + "\n")

    sh = ["cohort\tprotein_id\tother_protein_id\ttcpa_stored_rho\ttcpa_stored_p"]
    for rec in stored:
        sh.append(
            "\t".join(
                [
                    rec["cohort"],
                    rec["protein_id"],
                    rec["other_protein_id"],
                    "" if rec["tcpa_stored_rho"] is None else str(rec["tcpa_stored_rho"]),
                    "" if rec["tcpa_stored_p"] is None else str(rec["tcpa_stored_p"]),
                ]
            )
        )
    (out / "tcpa_stored_spearman.tsv").write_text("\n".join(sh) + "\n")
    (out / "download_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    if args.pancan_tsv:
        case_ids = {}
        for cohort in COHORTS:
            case_ids[cohort] = {
                rec["sample_id"] for rec in rows if rec["cohort"] == cohort and rec["sample_id"]
            }
        check_proteins = ["CLAUDIN7", "ECADHERIN", "EMA", "KU80", "HLADQA1"]
        pan_rows = pancan_missingness(Path(args.pancan_tsv), case_ids, check_proteins)
        ph = ["cohort\tprotein_id\tn_observed\tn_empty"]
        for rec in pan_rows:
            ph.append(
                f"{rec['cohort']}\t{rec['protein_id']}\t{rec['n_observed']}\t{rec['n_empty']}"
            )
        (out / "pancan_lung_missingness.tsv").write_text("\n".join(ph) + "\n")
        print(f"wrote pancan missingness from {args.pancan_tsv}")
    print(f"wrote {out}  abundance rows={len(rows)}  stored corr={len(stored)}")
    for cohort, info in meta["cohorts"].items():
        print(cohort, "proteins", info["n_protein_docs"], "fetched", list(info.get("fetched", )))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
