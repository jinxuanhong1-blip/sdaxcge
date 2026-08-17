#!/usr/bin/env python3
"""Download leftover TISCH2 CellMetainfo (always) and expression.h5 when needed.

Metainfo first. Skip h5 if:
  - file would exceed 2 GB
  - metainfo has no epithelial/malignant cells (CLDN4 cannot be scored)
  - metainfo has neither epithelial/malignant nor T/NK
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import DATASETS, EPITHELIAL_NONMALIGNANT, MALIGNANT, MAX_H5_BYTES, TISCH_BASE, TNK_LINEAGES


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def curl_head_bytes(url: str) -> int | None:
    cmd = ["curl", "-sI", "-L", "--max-time", "30", url]
    try:
        out = subprocess.check_output(cmd, text=True, errors="replace")
    except subprocess.CalledProcessError:
        return None
    for line in out.splitlines():
        if line.lower().startswith("content-length:"):
            try:
                return int(line.split(":", 1)[1].strip())
            except ValueError:
                return None
    return None


def curl_download(url: str, dest: Path, max_bytes: int | None = None, retries: int = 4) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        if max_bytes is not None and dest.stat().st_size > max_bytes:
            dest.unlink()
            return {"url": url, "path": str(dest), "ok": False, "error": f"existing file >{max_bytes} bytes"}
        return {"url": url, "path": str(dest), "bytes": dest.stat().st_size, "skipped": True, "ok": True}
    delay = 4
    last_err = None
    for attempt in range(1, retries + 1):
        cmd = ["curl", "-fL", "--retry", "3", "--retry-delay", "4", "-o", str(dest), url]
        if max_bytes is not None:
            cmd[1:1] = ["--max-filesize", str(max_bytes)]
        try:
            subprocess.run(cmd, check=True)
            if dest.exists() and dest.stat().st_size > 0:
                if max_bytes is not None and dest.stat().st_size > max_bytes:
                    dest.unlink()
                    return {"url": url, "ok": False, "error": f"downloaded file >{max_bytes} bytes"}
                return {
                    "url": url,
                    "path": str(dest),
                    "bytes": dest.stat().st_size,
                    "skipped": False,
                    "ok": True,
                    "attempt": attempt,
                }
            last_err = "empty file"
        except subprocess.CalledProcessError as e:
            last_err = str(e)
            if dest.exists():
                dest.unlink()
        time.sleep(delay)
        delay *= 2
    return {"url": url, "path": str(dest), "ok": False, "error": last_err}


def _norm(x) -> str:
    return str(x).strip()


def metainfo_gate(meta_path: Path) -> dict:
    import pandas as pd

    meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
    lineage_col = next((c for c in meta.columns if "major-lineage" in c.lower()), None)
    rec = {
        "n_cells": int(len(meta)),
        "lineage_col": lineage_col,
        "n_malignant": 0,
        "n_epithelial_nonmalig": 0,
        "n_tnk": 0,
        "has_cldn4_compartment": False,
        "has_tnk": False,
        "download_h5": False,
        "skip_h5_reason": None,
    }
    if lineage_col is None:
        rec["skip_h5_reason"] = "no major-lineage column"
        return rec
    lin = meta[lineage_col].map(_norm)
    rec["n_malignant"] = int(lin.isin(MALIGNANT).sum())
    rec["n_epithelial_nonmalig"] = int(lin.isin(EPITHELIAL_NONMALIGNANT).sum())
    rec["n_tnk"] = int(lin.isin(TNK_LINEAGES).sum())
    rec["has_cldn4_compartment"] = (rec["n_malignant"] + rec["n_epithelial_nonmalig"]) > 0
    rec["has_tnk"] = rec["n_tnk"] > 0
    if not rec["has_cldn4_compartment"] and not rec["has_tnk"]:
        rec["skip_h5_reason"] = "both CLDN4 compartment and T/NK absent"
        return rec
    if not rec["has_cldn4_compartment"]:
        rec["skip_h5_reason"] = "no epithelial/malignant cells; CLDN4 cannot be scored"
        return rec
    rec["download_h5"] = True
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-dir", default="data/tisch_leftover_merge_cldn4")
    ap.add_argument("--out-dir", default="methods/tisch_leftover_merge_cldn4")
    args = ap.parse_args()
    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

    provenance: dict = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "base": TISCH_BASE,
        "max_h5_bytes": MAX_H5_BYTES,
        "files": {},
        "gates": {},
    }

    for ds in DATASETS:
        meta_url = f"{TISCH_BASE}/{ds}/{ds}_CellMetainfo_table.tsv"
        meta_dest = data_dir / f"{ds}_CellMetainfo_table.tsv"
        print(f">> {ds} CellMetainfo", flush=True)
        rec = curl_download(meta_url, meta_dest)
        if rec.get("ok") and meta_dest.exists():
            rec["sha256"] = sha256_file(meta_dest)
        provenance["files"][f"{ds}_CellMetainfo_table.tsv"] = rec
        print(f"   ok={rec.get('ok')} bytes={rec.get('bytes')} err={rec.get('error')}", flush=True)
        if not rec.get("ok"):
            provenance["gates"][ds] = {"download_h5": False, "skip_h5_reason": "CellMetainfo missing"}
            continue

        gate = metainfo_gate(meta_dest)
        h5_url = f"{TISCH_BASE}/{ds}/{ds}_expression.h5"
        head_bytes = curl_head_bytes(h5_url)
        gate["h5_content_length"] = head_bytes
        if head_bytes is not None and head_bytes > MAX_H5_BYTES:
            gate["download_h5"] = False
            gate["skip_h5_reason"] = f"h5 Content-Length {head_bytes} > 2GB"
        provenance["gates"][ds] = gate
        print(
            f"   gate epi={gate['n_malignant']}+{gate['n_epithelial_nonmalig']} "
            f"tnk={gate['n_tnk']} download_h5={gate['download_h5']} "
            f"reason={gate['skip_h5_reason']}",
            flush=True,
        )
        if not gate["download_h5"]:
            continue

        h5_dest = data_dir / f"{ds}_expression.h5"
        print(f">> {ds} expression.h5", flush=True)
        h5_rec = curl_download(h5_url, h5_dest, max_bytes=MAX_H5_BYTES)
        if h5_rec.get("ok") and h5_dest.exists():
            h5_rec["sha256"] = sha256_file(h5_dest)
        provenance["files"][f"{ds}_expression.h5"] = h5_rec
        print(f"   ok={h5_rec.get('ok')} bytes={h5_rec.get('bytes')} err={h5_rec.get('error')}", flush=True)
        if not h5_rec.get("ok"):
            gate["download_h5"] = False
            gate["skip_h5_reason"] = h5_rec.get("error") or "h5 download failed"
            provenance["gates"][ds] = gate

    dest = out_dir / "download_provenance.json"
    dest.write_text(json.dumps(provenance, indent=2) + "\n")
    print(f"wrote {dest}")
    n_fail = sum(1 for r in provenance["files"].values() if not r.get("ok"))
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
