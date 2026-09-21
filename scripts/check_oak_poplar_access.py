#!/usr/bin/env python3
"""Record whether OAK/POPLAR linear RNA is publicly reachable.

Queries the public EGA metadata API only. Does not download controlled
files and does not compute ORR, DCR, or OS.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

META = "https://metadata.ega-archive.org"
OUT = Path("results/oak_poplar_access/access_status.json")

DATASETS = {
    "EGAD00001008390": "POPLAR log2(TPM+1)",
    "EGAD00001008391": "OAK log2(TPM+1)",
    "EGAD00001008548": "POPLAR clinical",
    "EGAD00001008549": "OAK clinical",
    "EGAD00001008550": "OAK biomarker",
    "EGAD00001008628": "OAK counts",
    "EGAD00001008629": "OAK CPM",
    "EGAD00001008630": "POPLAR counts",
    "EGAD00001008631": "POPLAR CPM",
    "EGAD00001007703": "OAK+POPLAR FASTQ",
}


def get_json(url: str) -> tuple[int, object]:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=40) as resp:
            return resp.status, json.load(resp)
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", "replace")[:500]
        return exc.code, {"error": body}


def main() -> None:
    rows = []
    for acc, role in DATASETS.items():
        status, meta = get_json(f"{META}/datasets/{acc}")
        fstatus, files = get_json(f"{META}/datasets/{acc}/files")
        file_rows = []
        n_files = len(files) if isinstance(files, list) else None
        # The FASTQ dataset is hundreds of controlled files. Keep the catalog
        # summary only; the TPM and clinical matrices are the ones to request.
        keep_files = isinstance(files, list) and len(files) <= 12
        if keep_files:
            for item in files:
                file_rows.append(
                    {
                        "accession_id": item.get("accession_id"),
                        "extension": item.get("extension"),
                        "filesize": item.get("filesize"),
                        "checksum_type": item.get("unencrypted_checksum_type"),
                        "checksum": item.get("unencrypted_checksum"),
                    }
                )
        rows.append(
            {
                "dataset": acc,
                "role": role,
                "http_status": status,
                "access_type": meta.get("access_type") if isinstance(meta, dict) else None,
                "title": meta.get("title") if isinstance(meta, dict) else None,
                "num_samples": meta.get("num_samples") if isinstance(meta, dict) else None,
                "policy_accession_id": meta.get("policy_accession_id") if isinstance(meta, dict) else None,
                "released_date": meta.get("released_date") if isinstance(meta, dict) else None,
                "files_http_status": fstatus,
                "n_files": n_files,
                "files": file_rows if keep_files else "omitted; see n_files",
            }
        )

    _, dac = get_json(f"{META}/dacs/EGAC00001002120")
    _, policy = get_json(f"{META}/policies/EGAP00001002123")
    contacts = dac.get("contacts") if isinstance(dac, dict) else None
    payload = {
        "checked_at_utc": datetime.now(timezone.utc).isoformat(),
        "linear_tpm_public": False,
        "association_tests_run": False,
        "reason": (
            "EGA metadata access_type is controlled for the POPLAR and OAK "
            "TPM and clinical datasets. No sample-level CLDN4 or TACSTD2 "
            "values were downloaded, and no ORR/DCR/OS contrast was computed."
        ),
        "dac": {
            "accession": "EGAC00001002120",
            "title": dac.get("title") if isinstance(dac, dict) else None,
            "contacts": contacts,
        },
        "policy": {
            "accession": "EGAP00001002123",
            "title": policy.get("title") if isinstance(policy, dict) else None,
            "dac_accession_id": policy.get("dac_accession_id") if isinstance(policy, dict) else None,
        },
        "datasets": rows,
    }
    if any(row["access_type"] != "controlled" for row in rows):
        raise SystemExit("Unexpected public dataset; refusing to continue without review.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"wrote {OUT}")
    print("linear_tpm_public", payload["linear_tpm_public"])
    print("association_tests_run", payload["association_tests_run"])


if __name__ == "__main__":
    main()
