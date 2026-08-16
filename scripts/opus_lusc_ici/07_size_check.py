#!/usr/bin/env python3
"""Confirm processed artefacts stay under the 2 GB slice budget."""

from __future__ import annotations

import json
import sys

from config import (LOGS_DIR, MANIFEST_DIR, MAX_PROCESSED_BYTES,
                    MAX_SINGLE_FILE_BYTES, RESULTS_DIR)


def main() -> int:
    files = [p for p in RESULTS_DIR.rglob("*") if p.is_file()]
    rows = []
    total = 0
    over = []
    for path in sorted(files):
        n = path.stat().st_size
        total += n
        rel = str(path.relative_to(RESULTS_DIR))
        rows.append({"file": rel, "bytes": n})
        if n > MAX_SINGLE_FILE_BYTES:
            over.append(rel)
    report = {
        "n_files": len(files),
        "total_bytes": total,
        "total_mb": round(total / 1e6, 2),
        "budget_bytes": MAX_PROCESSED_BYTES,
        "under_2gb": total < MAX_PROCESSED_BYTES,
        "files_over_90mb": over,
    }
    (MANIFEST_DIR / "processed_size.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    if total >= MAX_PROCESSED_BYTES or over:
        return 1
    (LOGS_DIR / "07_size_check.done").write_text("ok\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
