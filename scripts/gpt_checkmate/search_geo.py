#!/usr/bin/env python3
"""Record exact-name searches of NCBI GEO DataSets for requested trials."""

from __future__ import annotations

import json
import sys
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "results" / "gpt_checkmate" / "geo_search.json"
API = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
TRIALS = ("017", "057", "227", "9LA", "816", "153")


def main() -> int:
    searches: dict[str, object] = {}
    for trial in TRIALS:
        query = f'"CheckMate {trial}"'
        parameters = {
            "db": "gds",
            "term": query,
            "retmode": "json",
            "retmax": 100,
        }
        url = f"{API}?{urllib.parse.urlencode(parameters)}"
        request = urllib.request.Request(
            url, headers={"User-Agent": "gpt-checkmate-geo-search/1.0"}
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            result = json.load(response)["esearchresult"]
        searches[trial] = {
            "query": query,
            "count": int(result["count"]),
            "ids": result["idlist"],
            "query_translation": result["querytranslation"],
        }

    payload = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "database": "NCBI GEO DataSets (gds)",
        "api": API,
        "scope_note": (
            "Exact trial-name discovery only; a zero is not proof that no indirectly "
            "described or newly deposited dataset exists."
        ),
        "searches": searches,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(f"Wrote GEO search record: {OUTPUT}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, KeyError, TypeError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        raise SystemExit(1)
