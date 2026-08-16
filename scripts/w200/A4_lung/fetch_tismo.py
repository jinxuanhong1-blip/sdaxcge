#!/usr/bin/env python3
"""Download TISMO in-vivo metadata and ICB gene tables used by A4_lung."""

from __future__ import annotations

import argparse
import json
import ssl
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

PORT = "https://tismo.pku-genomics.org/tismo"
RPORT = "https://tismo.pku-genomics.org/rtismo"
ICB_TREATMENTS = [
    "antiCTLA4",
    "antiCTLA4&antiPD1",
    "antiCTLA4&antiPDL1",
    "antiPD1",
    "antiPDL1",
    "antiPDL2",
]
TARGET_GENES = ["Tacstd2", "Cldn4"]
CONTROL_GENES = ["Actb", "Cd8a", "Ifng", "Gzmb", "Cd274", "Epcam"]
CTX = ssl.create_default_context()


def post(url: str, data: dict[str, str] | None = None, json_body: dict | None = None, timeout: int = 180) -> bytes:
    if json_body is not None:
        payload = json.dumps(json_body).encode()
        headers = {"Content-Type": "application/json"}
    else:
        payload = urllib.parse.urlencode(data or {}).encode()
        headers = {"Content-Type": "application/x-www-form-urlencoded"}
    req = urllib.request.Request(url, data=payload, headers=headers, method="POST")
    last_err: Exception | None = None
    for attempt in range(1, 5):
        try:
            with urllib.request.urlopen(req, timeout=timeout, context=CTX) as resp:
                return resp.read()
        except (urllib.error.URLError, TimeoutError) as err:
            last_err = err
            time.sleep(2**attempt)
    raise RuntimeError(f"POST {url} failed after retries: {last_err}")


def write_bytes(path: Path, blob: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(blob)


def fetch_vivo_meta(raw_dir: Path) -> None:
    blob = post(f"{PORT}/metaData/vivoMeta", json_body={})
    payload = json.loads(blob)
    if payload.get("status") != 600200:
        raise RuntimeError(f"vivoMeta unexpected status: {payload.get('status')}")
    dest = raw_dir / "vivo_meta.json"
    dest.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {dest}  n={len(payload['data'])}")


def fetch_gene_icb(raw_dir: Path, gene: str, tumor: str, dest_name: str) -> None:
    icb_list = json.dumps(ICB_TREATMENTS)
    tumor_list = json.dumps([tumor])
    blob = post(
        f"{RPORT}/gene/downVivoExprn",
        data={
            "filename": dest_name,
            "type": "3",
            "gene": gene,
            "icbList": icb_list,
            "tumorList": tumor_list,
        },
    )
    if blob[:1] not in (b"S", b"s") and b"Samples" not in blob[:80]:
        raise RuntimeError(f"{gene} download did not look like a CSV ({blob[:80]!r})")
    dest = raw_dir / dest_name
    write_bytes(dest, blob)
    n = blob.count(b"\n") - 1
    print(f"wrote {dest}  rows≈{n}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("results/w200/A4_lung/raw"),
        help="Directory for raw TISMO exports",
    )
    args = parser.parse_args()
    raw_dir = args.out
    raw_dir.mkdir(parents=True, exist_ok=True)

    fetch_vivo_meta(raw_dir)
    for gene in TARGET_GENES:
        fetch_gene_icb(raw_dir, gene, "All", f"{gene.lower()}_icb_all.csv")
    for gene in CONTROL_GENES:
        fetch_gene_icb(raw_dir, gene, "LLC", f"{gene.lower()}_llc.csv")

    manifest = {
        "source": "https://tismo.pku-genomics.org/",
        "paper": "Zeng et al. Nucleic Acids Res 2022; PMID 34534350",
        "icb_treatments": ICB_TREATMENTS,
        "target_genes": TARGET_GENES,
        "control_genes": CONTROL_GENES,
        "note": (
            "type=3 is TISMO's per-sample table export from /rtismo/gene/downVivoExprn. "
            "Querying tumorList=['All'] is required for Tacstd2/Cldn4 so that "
            "SRX8918393 (present in the All-models table, dropped by the LLC-only "
            "Tacstd2 query) is retained."
        ),
    }
    (raw_dir / "fetch_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print("done")


if __name__ == "__main__":
    main()
