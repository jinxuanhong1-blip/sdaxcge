#!/usr/bin/env python3
"""Download official TISMO in-vivo ICB expression tables (Cldn4 + Tacstd2)."""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

TISMO_DOWN = "https://tismo.pku-genomics.org/rtismo/gene/downVivoExprn"
TISMO_VIVO_META = "https://tismo.pku-genomics.org/tismo/metaData/vivoMeta"
OUT = Path("results/w200/A4_Cldn4/data")


def post_form(url: str, fields: dict[str, str], timeout: int = 180) -> bytes:
    body = urllib.parse.urlencode(fields).encode()
    req = urllib.request.Request(url, data=body, method="POST")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read()


def download_gene(gene: str, dest: Path) -> None:
    raw = post_form(
        TISMO_DOWN,
        {
            "filename": "genetreatment_vivo.csv",
            "type": "3",
            "gene": gene,
            "icbList": '["All"]',
            "tumorList": '["All"]',
        },
    )
    if raw[:1] != b"S" and b"Samples" not in raw[:80]:
        raise RuntimeError(f"{gene}: unexpected payload ({raw[:120]!r})")
    dest.write_bytes(raw)
    print(f"wrote {dest} ({len(raw)} bytes)")


def download_vivo_meta(dest: Path) -> None:
    raw = post_form(TISMO_VIVO_META, {}, timeout=60)
    dest.write_bytes(raw)
    payload = json.loads(raw.decode())
    n = len(payload.get("data") or [])
    print(f"wrote {dest} (n_groups={n})")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    download_gene("Cldn4", OUT / "cldn4_tismo_vivo.csv")
    download_gene("Tacstd2", OUT / "tacstd2_tismo_vivo.csv")
    download_vivo_meta(OUT / "tismo_vivo_meta.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
