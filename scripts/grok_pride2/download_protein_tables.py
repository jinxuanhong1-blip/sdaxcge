#!/usr/bin/env python3
"""Download protein tables only (never raw MS) for kept PX datasets."""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.request
from pathlib import Path

UA = "grok-pride2/1.0 (research; protein-tables-only; no-raw-ms)"
ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "grok_pride2"
TABLES = RESULTS / "protein_tables"
TABLES.mkdir(parents=True, exist_ok=True)

MAX_BYTES = 80 * 1024 * 1024  # 80 MB hard cap — protein tables, not search dumps
RAW_RE = re.compile(
    r"\.(raw|wiff|wiff\.scan|mzml|mzxml|mgf|d|raw\.gz|mzml\.gz|thermo|bruker)(\.|$)",
    re.I,
)
PROTEIN_RE = re.compile(
    r"(proteinGroups|protein[_\- ]?group|proteins\.txt|protein\.txt|"
    r"pg_matrix|protein_matrix|diann|maxquant|proteome.?discover|"
    r"pd_proteins|mzTab|proteinquant|quantified.?protein|"
    r"differential.?protein|dep[_\-]?table|supplement.*\.(csv|tsv|txt|xlsx)|"
    r"search\.zip|txt\.zip|evidence\.txt|peptides\.txt|"
    r"data_tf_response|dat_cs_all|NSClibrary|"
    r"protein.*\.(csv|tsv|txt|xlsx|xls|zip))",
    re.I,
)
SKIP_NAME_RE = re.compile(
    r"(checksum|readme|\.raw|\.wiff|peaklist|mgf|mzml|mzxml|fasta|"
    r"spectrum|spectra|sdrf$)",
    re.I,
)


def http_download(url: str, dest: Path, max_bytes: int = MAX_BYTES) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    written = 0
    try:
        with urllib.request.urlopen(req, timeout=180) as resp, tmp.open("wb") as fh:
            cl = resp.headers.get("Content-Length")
            if cl and int(cl) > max_bytes:
                return {
                    "ok": False,
                    "reason": f"content-length {cl} exceeds cap {max_bytes}",
                    "bytes": 0,
                }
            while True:
                chunk = resp.read(1024 * 256)
                if not chunk:
                    break
                written += len(chunk)
                if written > max_bytes:
                    fh.close()
                    tmp.unlink(missing_ok=True)
                    return {
                        "ok": False,
                        "reason": f"stream exceeded cap after {written} bytes",
                        "bytes": written,
                    }
                fh.write(chunk)
        tmp.replace(dest)
        return {"ok": True, "reason": "downloaded", "bytes": written}
    except Exception as exc:  # noqa: BLE001
        tmp.unlink(missing_ok=True)
        return {"ok": False, "reason": str(exc), "bytes": written}


def main() -> None:
    inv_path = RESULTS / "px_file_inventory.csv"
    keep_path = RESULTS / "px_keep.csv"
    if not inv_path.exists() or not keep_path.exists():
        raise SystemExit("run search_px.py first")

    keep = {r["accession"] for r in csv.DictReader(keep_path.open())}
    log_rows = []
    for row in csv.DictReader(inv_path.open()):
        acc = row["accession"]
        if acc not in keep:
            continue
        name = row["fileName"]
        url = row.get("download_http") or ""
        size = int(row.get("size_bytes") or 0)
        is_raw = str(row.get("is_raw")).lower() == "true"
        candidate = str(row.get("protein_table_candidate")).lower() == "true"
        dest = TABLES / acc / name
        decision = "skip"
        reason = ""
        if is_raw or RAW_RE.search(name):
            decision, reason = "skip_raw", "raw/peak MS file"
        elif SKIP_NAME_RE.search(name) and not candidate:
            decision, reason = "skip_name", "not a protein table"
        elif not url:
            decision, reason = "skip_nourl", "no HTTP location"
        elif size and size > MAX_BYTES:
            decision, reason = "skip_large", f"size {size} > {MAX_BYTES}"
        elif candidate or PROTEIN_RE.search(name):
            if dest.exists() and dest.stat().st_size > 0:
                decision, reason = "exists", "already downloaded"
            else:
                print(f"GET {acc} {name} ({size})")
                result = http_download(url, dest)
                decision = "ok" if result["ok"] else "fail"
                reason = result["reason"]
                size = result["bytes"] or size
                time.sleep(0.2)
        else:
            decision, reason = "skip_not_protein", "filename not protein-table-like"
        log_rows.append(
            {
                "accession": acc,
                "fileName": name,
                "decision": decision,
                "reason": reason,
                "size_bytes": size,
                "url": url,
                "local_path": str(dest.relative_to(ROOT)) if dest.exists() else "",
            }
        )

    out = RESULTS / "download_log.csv"
    with out.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(log_rows[0].keys()) if log_rows else ["accession"])
        w.writeheader()
        w.writerows(log_rows)
    summary = {
        "n_considered": len(log_rows),
        "by_decision": {},
    }
    from collections import Counter

    summary["by_decision"] = dict(Counter(r["decision"] for r in log_rows))
    (RESULTS / "download_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
