#!/usr/bin/env python3
"""Download TISMO in-vivo ICB Tacstd2 expression for all listed tumors."""

from __future__ import annotations

import json
from pathlib import Path

import requests

PORT = "https://tismo.pku-genomics.org/tismo"
RPORT = "https://tismo.pku-genomics.org/rtismo"
OUT = Path("results/w200/A4_all64/tacstd2_vivo_icb.csv")

HEADERS = {
    "User-Agent": "sdaxcge-A4-recompute/1.0",
    "Origin": "https://tismo.pku-genomics.org",
    "Referer": "https://tismo.pku-genomics.org/",
}


def session() -> requests.Session:
    s = requests.Session()
    s.headers.update(HEADERS)
    return s


def names(endpoint: str, extra: dict | None = None) -> list[str]:
    r = session().post(PORT + endpoint, data=extra or {}, timeout=60)
    r.raise_for_status()
    payload = r.json()
    if payload.get("status") != 600200:
        raise RuntimeError(f"{endpoint} failed: {payload}")
    return [row["name"] for row in payload["data"] if row["name"] != "All"]


def main() -> None:
    treatments = names("/gene/getVivoTreatment")
    tumors = names("/gene/getVivoCohort", {"treatment": ""})
    print("treatments", treatments)
    print("tumors", tumors)

    data = {
        "filename": "vivo.csv",
        "type": "3",
        "gene": "Tacstd2",
        "icbList": json.dumps(treatments),
        "tumorList": json.dumps(tumors),
    }
    r = session().post(RPORT + "/gene/downVivoExprn", data=data, timeout=180)
    r.raise_for_status()
    if r.content[:7] != b"Samples":
        raise RuntimeError(f"unexpected payload ({len(r.content)} bytes): {r.content[:200]!r}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(r.content)
    print(f"wrote {OUT} ({len(r.content)} bytes)")


if __name__ == "__main__":
    main()
