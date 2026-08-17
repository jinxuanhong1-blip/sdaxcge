#!/usr/bin/env python3
"""Download public TISCH2 GSE148071 expression.h5 + CellMetainfo."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
TISCH_BASE = "https://tisch.compbio.cn/static/data/NSCLC_GSE148071"
FILES = [
    {
        "url": f"{TISCH_BASE}/NSCLC_GSE148071_CellMetainfo_table.tsv",
        "dest": "NSCLC_GSE148071_CellMetainfo_table.tsv",
        "role": "cell_metainfo",
    },
    {
        "url": f"{TISCH_BASE}/NSCLC_GSE148071_expression.h5",
        "dest": "NSCLC_GSE148071_expression.h5",
        "role": "tisch2_expression_h5",
    },
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"HAVE {dest} ({dest.stat().st_size} bytes)", flush=True)
        return
    tmp = dest.with_suffix(dest.suffix + ".partial")
    print(f"GET  {url}", flush=True)
    cmd = ["curl", "-fL", "--retry", "4", "--retry-delay", "8", "-o", str(tmp), url]
    subprocess.run(cmd, check=True)
    tmp.replace(dest)
    print(f"OK   {dest} ({dest.stat().st_size} bytes)", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", type=Path, default=HERE / "data")
    args = ap.parse_args()
    args.datadir.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in FILES:
        dest = args.datadir / spec["dest"]
        fetch(spec["url"], dest)
        rows.append(
            {
                "file": spec["dest"],
                "role": spec["role"],
                "bytes": dest.stat().st_size,
                "sha256": sha256(dest),
                "url": spec["url"],
            }
        )
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "accession": "GSE148071",
        "tisch_id": "NSCLC_GSE148071",
        "paper": "Wu et al. 2021 Nat Commun PMID 33953163",
        "files": rows,
    }
    out = args.datadir / "file_manifest.json"
    out.write_text(json.dumps(manifest, indent=2))
    print(f"wrote {out}", flush=True)


if __name__ == "__main__":
    main()
