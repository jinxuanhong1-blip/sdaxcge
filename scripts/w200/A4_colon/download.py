#!/usr/bin/env python3
"""Download the TISMO public tables needed for the colon-only Tacstd2 ICB test.

Writes unmodified payloads to results/w200/A4_colon/data/.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tismo_client as tc

GENE = "Tacstd2"


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    print(f"  wrote {path.relative_to(tc.REPO_ROOT)} ({len(text):,} bytes)")


def main() -> int:
    tc.DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/4] vocabularies")
    treatments = tc.get_vivo_treatments()
    models = tc.get_vivo_cohorts(treatments)
    _write(
        tc.DATA_DIR / "vocabularies.json",
        json.dumps(
            {
                "icb_treatments": treatments,
                "tumor_models": models,
                "api_meta_base": tc.META_BASE,
                "api_r_base": tc.R_BASE,
            },
            indent=2,
        ),
    )

    print("[2/4] metadata")
    for kind in ("vivoMeta", "cellLineMeta"):
        rows = tc.get_metadata(kind)
        _write(tc.DATA_DIR / f"{kind}.json", json.dumps(rows, indent=2))

    cancer_types = tc.load_cell_line_cancer_types()
    colon = tc.colon_models(cancer_types)
    print(f"  colorectal models in cellLineMeta: {colon}")
    print(f"  of those in the ICB Gene-module model list: "
          f"{[m for m in colon if m in models]}")

    print(f"[3/4] {GENE} in vivo expression (all ICB Gene-module models)")
    # Request the full ICB model list. Colon filtering is an analysis step so
    # the denominator (how many models TISMO actually serves) stays visible.
    text = tc.download_vivo_expression(GENE, treatments, models)
    _write(tc.DATA_DIR / f"vivo_expression_{GENE}.csv", text)

    print("[4/4] done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
