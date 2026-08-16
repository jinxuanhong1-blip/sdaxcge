#!/usr/bin/env python3
"""Download public matrices for the DepMap IFN-treated conceptual rework.

Sources
-------
- DepMap Public 24Q4 (Figshare+ 10.25452/figshare.plus.27993248.v1):
  Model.csv, expression, CRISPR Chronos gene effect.
- PRISM Repurposing Public 23Q2 (Figshare 10.6084/m9.figshare.23600310.v4):
  extended primary LFC matrix + compound / cell-line metadata.
- CCLE proteomics, Nusinow et al. Cell 2020 (Gygi lab spreadsheet).

There is no public DepMap IFN-α/γ *treated transcriptome*. PRISM is a
small-molecule viability screen; recombinant IFN biologics are not in the
library. This download is the closest public stack (IFN-pathway drugs,
CRISPR IFN-gene effects, protein MHC-I / ISGs).
"""
from __future__ import annotations

import argparse
import hashlib
import sys
import urllib.request
from pathlib import Path

FILES = {
    "Model.csv": {
        "url": "https://ndownloader.figshare.com/files/51065297",
        "md5": "675210d17675f3517b0ce39a3c274f16",
        "min_bytes": 100_000,
    },
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": {
        "url": "https://ndownloader.figshare.com/files/51065489",
        "md5": "71794802b750ce77c422dad0720a40af",
        "min_bytes": 400_000_000,
    },
    "CRISPRGeneEffect.csv": {
        "url": "https://ndownloader.figshare.com/files/51064667",
        "md5": None,
        "min_bytes": 300_000_000,
    },
    "Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv": {
        "url": "https://ndownloader.figshare.com/files/41419821",
        "md5": "fc3c45140a4141650d92ad7e23fb86bd",
        "min_bytes": 50_000_000,
    },
    "Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv": {
        "url": "https://ndownloader.figshare.com/files/41419893",
        "md5": "60dc22963aa20e1cf4758800bcaadfd3",
        "min_bytes": 100_000,
    },
    "Repurposing_Public_23Q2_Cell_Line_Meta_Data.csv": {
        "url": "https://ndownloader.figshare.com/files/41419815",
        "md5": "f2735ba9867fd32ffca83ce29c98a2b3",
        "min_bytes": 10_000,
    },
    "protein_quant_current_normalized.csv.gz": {
        "url": "https://gygi.hms.harvard.edu/data/ccle/protein_quant_current_normalized.csv.gz",
        "md5": None,
        "min_bytes": 5_000_000,
        "fallbacks": [
            "https://gygi.med.harvard.edu/sites/gygi.med.harvard.edu/files/documents/protein_quant_current_normalized.csv.gz",
        ],
    },
}


def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"GET {url} -> {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": "DepMap_IFN_treated/1.0"})
    with urllib.request.urlopen(req, timeout=600) as r, tmp.open("wb") as out:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def fetch(meta: dict, dest: Path, skip_md5: bool) -> None:
    urls = [meta["url"]] + list(meta.get("fallbacks") or [])
    last_err = None
    for url in urls:
        try:
            download(url, dest)
            if dest.stat().st_size >= meta["min_bytes"]:
                return
            last_err = f"too small: {dest.stat().st_size}"
        except Exception as e:  # noqa: BLE001 — try next mirror
            last_err = str(e)
            print(f"WARN {url}: {e}", flush=True)
    raise RuntimeError(f"failed {dest.name}: {last_err}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="data/rework/DepMap_IFN_treated")
    ap.add_argument("--skip-md5", action="store_true")
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name, meta in FILES.items():
        dest = outdir / name
        if dest.exists() and dest.stat().st_size >= meta["min_bytes"]:
            if meta["md5"] and not args.skip_md5:
                got = md5sum(dest)
                if got != meta["md5"]:
                    print(f"MD5 mismatch {name}: {got} != {meta['md5']}; re-downloading", flush=True)
                    dest.unlink()
                else:
                    print(f"OK exists {name} ({dest.stat().st_size} bytes)", flush=True)
                    continue
            else:
                print(f"OK exists {name} ({dest.stat().st_size} bytes)", flush=True)
                continue
        try:
            fetch(meta, dest, args.skip_md5)
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {name}: {e}", file=sys.stderr)
            return 1
        if dest.stat().st_size < meta["min_bytes"]:
            print(f"ERROR {name} too small: {dest.stat().st_size}", file=sys.stderr)
            return 1
        if meta["md5"] and not args.skip_md5:
            got = md5sum(dest)
            if got != meta["md5"]:
                print(f"ERROR MD5 {name}: {got} != {meta['md5']}", file=sys.stderr)
                return 1
        print(f"OK {name} ({dest.stat().st_size} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
