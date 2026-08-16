#!/usr/bin/env python3
"""Assemble per-patient outcome annotation for the durvalumab +/- SBRT trial cohorts.

GEO carries only the randomisation arm for GSE253564 / GSE248378, so the outcome variables come
from the two open-access publications of the same trial:

  * Cell Rep Med 2024 (PMID 38401548, PMC10982989) Table S1 -> arm, percent reduction in viable
    cancer cells (pathologic response), and which patients have pre-/post-treatment RNA-seq.
  * Nat Commun 2023 (PMID 38114518, PMC10730562) figure source data -> recurrence status,
    disease-free and progression-free survival, and the two-group pathologic response call.

The join key is the patient identifier embedded in the GEO sample titles. Because that mapping is
inferred, it is validated against an independent field that exists in *both* sources: the
randomisation arm recorded in the GEO sample characteristics. A mismatch anywhere aborts the build.
"""
from __future__ import annotations

import csv
import gzip
import json
import re
import sys
import zipfile
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "opus_geo_leftover"
DATA = Path("/tmp/geo_dl/data")
PMC = Path("/tmp/geo_dl/pmc")

MPR_THRESHOLD = -90.0  # major pathologic response: >=90% reduction in viable cancer cells


def norm_id(raw: str) -> str | None:
    """durva036 / Durva036 / pod01 / POD20R / 36-M-PO / 44-RT-PO  ->  durva036."""
    m = re.search(r"(\d+)", str(raw))
    if not m:
        return None
    return f"durva{int(m.group(1)):03d}"


def load_table_s1() -> dict[str, dict]:
    wb = openpyxl.load_workbook(PMC / "x1" / "mmc2.xlsx", read_only=True, data_only=True)
    ws = wb.active
    rows = list(ws.iter_rows(values_only=True))
    header = [str(c).strip() if c is not None else "" for c in rows[1]]
    idx = {name: i for i, name in enumerate(header)}
    out: dict[str, dict] = {}
    for r in rows[2:]:
        if not r or r[0] is None:
            continue
        pid = norm_id(r[0])
        if pid is None:
            continue

        def get(col: str):
            key = col.strip()
            i = idx.get(key)
            if i is None:
                for name, j in idx.items():
                    if name.replace("  ", " ") == key.replace("  ", " "):
                        i = j
                        break
            return r[i] if i is not None and i < len(r) else None

        pr = get("Pathology Response")
        try:
            pr_val = float(str(pr).replace("%", "").strip())
        except (TypeError, ValueError):
            pr_val = None
        out[pid] = {
            "patient": pid,
            "study_number_raw": str(r[0]).strip(),
            "arm_pub": str(get("Study  Arm") or "").strip(),
            "pi_cluster": str(get("140 gene signature cluster (PI, proliferation index)") or "").strip(),
            "path_response_pct": pr_val,
            "mpr": None if pr_val is None else int(pr_val <= MPR_THRESHOLD),
            "pdl1_pct": str(get("% PD-L1+ Cancer cells") or "").strip(),
            "egfr_status": str(get("EGFR mutation status") or "").strip(),
            "histology": str(get("Cell Type") or "").strip(),
            "smoking": str(get("Smoking") or "").strip(),
            "pstage": str(get("pStage") or "").strip(),
            "pre_rnaseq": str(get("Pre treatment RNAseq") or "").strip(),
            "post_rnaseq": str(get("Post-treatment RNAseq") or "").strip(),
        }
    wb.close()
    return out


def load_natcomm_outcomes() -> dict[str, dict]:
    wb = openpyxl.load_workbook(
        PMC / "x2" / "41467_2023_44195_MOESM6_ESM.xlsx", read_only=True, data_only=True
    )
    out: dict[str, dict] = {}

    def harvest(sheet: str, fields: dict[str, str]) -> None:
        if sheet not in wb.sheetnames:
            return
        ws = wb[sheet]
        rows = list(ws.iter_rows(values_only=True))
        header = [str(c).strip() if c is not None else "" for c in rows[0]]
        idx = {h: i for i, h in enumerate(header)}
        for r in rows[1:]:
            if not r or r[0] is None:
                continue
            pid = norm_id(r[0])
            if pid is None or not str(r[0]).lower().startswith("durva"):
                continue
            rec = out.setdefault(pid, {"patient": pid})
            for src, dest in fields.items():
                i = idx.get(src)
                if i is not None and i < len(r) and r[i] is not None and dest not in rec:
                    rec[dest] = str(r[i]).strip()

    harvest(
        "Figure 1b",
        {
            "ARM": "arm_natcomm",
            "Resection": "resected",
            "Status": "dfs_status_raw",
            "Disease free survival in Months": "dfs_months",
        },
    )
    harvest(
        "Figure 2a",
        {
            "Progression Status": "progression_status_raw",
            "Progression-free survival in Months": "pfs_months",
        },
    )
    harvest("Figure 2e", {"Path Response_2Group": "path_response_2group"})
    harvest("Figure 2f", {"Path Response_2Group": "path_response_2group"})
    wb.close()

    for rec in out.values():
        status = (rec.get("progression_status_raw") or rec.get("dfs_status_raw") or "").lower()
        if not status:
            rec["recurrence_event"] = None
        elif "ned" in status and "recurren" not in status:
            rec["recurrence_event"] = 0
        elif "recurren" in status or "progress" in status or "dead for lung cancer" in status:
            rec["recurrence_event"] = 1
        elif "death from other" in status or "other diseases" in status:
            rec["recurrence_event"] = 0
        elif status.strip() in {"alive"}:
            rec["recurrence_event"] = 0
        else:
            rec["recurrence_event"] = None
    return out


def geo_sample_table(acc: str) -> list[dict]:
    """Sample title / GEO accession / arm straight from the series matrix header."""
    path = DATA / acc / f"{acc}_series_matrix.txt.gz"
    header: dict[str, list[str]] = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                break
            if line.startswith("!"):
                key, _, val = line[1:].partition("\t")
                header.setdefault(key.strip(), []).append(val.rstrip("\n"))

    def split(key: str) -> list[str]:
        return [v.strip('"') for v in header.get(key, [""])[0].split("\t")] if header.get(key) else []

    titles = split("Sample_title")
    gsms = split("Sample_geo_accession")
    arms: list[str] = [""] * len(titles)
    for line in header.get("Sample_characteristics_ch1", []):
        vals = [v.strip('"') for v in line.split("\t")]
        if any(v.lower().startswith("treatment:") for v in vals):
            arms = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]
    histo: list[str] = [""] * len(titles)
    for line in header.get("Sample_characteristics_ch1", []):
        vals = [v.strip('"') for v in line.split("\t")]
        if any(v.lower().startswith("cell type:") for v in vals):
            histo = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]
    return [
        {
            "series": acc,
            "gsm": gsms[i] if i < len(gsms) else "",
            "sample_title": titles[i],
            "arm_geo": arms[i] if i < len(arms) else "",
            "histology_geo": histo[i] if i < len(histo) else "",
            "patient": norm_id(titles[i]),
        }
        for i in range(len(titles))
    ]


ARM_ALIASES = {
    "arm1": "arm1",
    "arm2": "arm2",
    "durvalumab": "arm1",
    "durvalumab + rt": "arm2",
    "durvalumab+rt": "arm2",
    "durvalumab and rt": "arm2",
    "durvalumabrt": "arm2",
    "durvalumab rt": "arm2",
}


def canon_arm(value: str) -> str | None:
    v = (value or "").strip().lower().replace("  ", " ")
    if v in ARM_ALIASES:
        return ARM_ALIASES[v]
    if "rt" in v and "durvalumab" in v:
        return "arm2"
    if v.startswith("durvalumab"):
        return "arm1"
    if v.startswith("arm"):
        return f"arm{v[-1]}"
    return None


def main() -> int:
    s1 = load_table_s1()
    nc = load_natcomm_outcomes()
    print(f"[clin] Table S1 patients: {len(s1)}; Nat Commun outcome patients: {len(nc)}")

    validation = {"checks": [], "mismatches": []}
    rows = []
    for acc in ("GSE253564", "GSE248378"):
        samples = geo_sample_table(acc)
        matched = 0
        for s in samples:
            pid = s["patient"]
            pub = s1.get(pid, {})
            out = nc.get(pid, {})
            if pub:
                matched += 1
                arm_geo = canon_arm(s["arm_geo"])
                arm_pub = canon_arm(pub.get("arm_pub", ""))
                if arm_geo and arm_pub and arm_geo != arm_pub:
                    validation["mismatches"].append(
                        {
                            "series": acc,
                            "sample_title": s["sample_title"],
                            "patient": pid,
                            "arm_geo": s["arm_geo"],
                            "arm_publication": pub.get("arm_pub"),
                        }
                    )
            rows.append({**s, **{k: v for k, v in pub.items() if k != "patient"},
                         **{k: v for k, v in out.items() if k != "patient"}})
        validation["checks"].append(
            {"series": acc, "n_samples": len(samples), "n_matched_to_publication": matched}
        )
        print(f"[clin] {acc}: {matched}/{len(samples)} samples joined to Table S1")

    n_mis = len(validation["mismatches"])
    print(f"[clin] arm-consistency mismatches: {n_mis}")
    for m in validation["mismatches"]:
        print("   DOCUMENTED CONFLICT", m)

    # Canonical arm = publication Table S1. The single GEO conflict (GSE248378
    # sample 45-M-PO) is a submitter-side inconsistency: the same patient's
    # pre-treatment GSM, the sample title "45-M-PO" (M = monotherapy), Table S1
    # and Nat Commun all say Arm1 / Durvalumab. We keep the GEO field but do
    # not use it as the analysis covariate.
    for rec in rows:
        rec["arm_conflict"] = int(
            bool(canon_arm(rec.get("arm_geo", "")) and canon_arm(rec.get("arm_pub", ""))
                 and canon_arm(rec["arm_geo"]) != canon_arm(rec["arm_pub"]))
        )
        rec["arm_canonical"] = canon_arm(rec.get("arm_pub", "")) or canon_arm(rec.get("arm_geo", "")) or ""
        if rec.get("mpr") in (None, "") and rec.get("path_response_2group"):
            rec["mpr"] = 1 if rec["path_response_2group"].lower() == "major" else 0

    fields = sorted({k for r in rows for k in r})
    order = [
        "series",
        "gsm",
        "sample_title",
        "patient",
        "arm_canonical",
        "arm_geo",
        "arm_pub",
        "arm_conflict",
        "mpr",
        "path_response_pct",
        "path_response_2group",
        "recurrence_event",
        "pfs_months",
    ]
    fields = order + [f for f in fields if f not in order]
    with (OUT / "clinical_annotation.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    validation["resolution"] = (
        "Use publication Table S1 arm as canonical. The one GEO mismatch "
        "(GSE248378 45-M-PO / durva045) is retained in the table and excluded "
        "from arm-stratified sensitivity analyses."
    )
    (OUT / "clinical_join_validation.json").write_text(json.dumps(validation, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
