#!/usr/bin/env python3
"""Download TISMO in-vivo ICB + in-vitro IFNg tables for Tacstd2 / Cldn4 / TJ genes."""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "tismo"

VIVO_GENES = ["Tacstd2", "Cldn4", "Cldn3", "Cldn6", "Cldn7", "Cdh1", "F11r", "Ocln"]
VITRO_GENES = VIVO_GENES[:]
CYTOKINES = ["IFNg"]


def write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    print(f"  wrote {path} ({len(text):,} bytes)")


def main() -> int:
    DATA.mkdir(parents=True, exist_ok=True)

    print("[1] vocabularies")
    vivo_tx = tc.names("/gene/getVivoTreatment")
    vivo_models = tc.names("/gene/getVivoCohort", {"treatment": ""})
    vitro_tx = tc.names("/gene/getVitroTreatment")
    vitro_models = tc.names("/gene/getVitroCohort", {"treatment": '"IFNg"'})
    write(
        DATA / "vocabularies.json",
        json.dumps(
            {
                "vivo_treatments": vivo_tx,
                "vivo_models": vivo_models,
                "vitro_treatments": vitro_tx,
                "vitro_ifng_models": vitro_models,
                "api_meta_base": tc.META_BASE,
                "api_r_base": tc.R_BASE,
            },
            indent=2,
        )
        + "\n",
    )

    print("[2] metadata")
    for kind in ("vivoMeta", "vitroMeta", "cellLineMeta"):
        rows = tc.get_metadata(kind)
        write(DATA / f"{kind}.json", json.dumps(rows, indent=2) + "\n")

    print("[3] in vivo ICB expression")
    for gene in VIVO_GENES:
        print(f"  {gene}")
        text = tc.download_vivo_expression(gene, vivo_tx, vivo_models)
        write(DATA / "vivo" / f"{gene}.csv", text)

    print("[4] in vitro IFNg expression")
    for gene in VITRO_GENES:
        print(f"  {gene}")
        text = tc.download_vitro_expression(gene, CYTOKINES, vitro_models)
        write(DATA / "vitro" / f"{gene}.csv", text)

    print("done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
