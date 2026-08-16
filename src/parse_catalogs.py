"""Parse NGDC OMIX release-list HTML and GSA-Human finished.json into TSVs."""

from __future__ import annotations

import argparse
import csv
import html
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from ngdc_client import fetch, setup_logging  # noqa: E402

OMIX_URL = "https://ngdc.cncb.ac.cn/omix/releaseList"
HRA_URL = "https://ngdc.cncb.ac.cn/gsa-human/json/finished.json"
HRA_REFERER = "https://ngdc.cncb.ac.cn/gsa-human/browse"

ROW_RE = re.compile(
    r"<tr>\s*<td>(OMIX\d+)</td>\s*<td>(.*?)</td>\s*<td[^>]*>(.*?)</td>"
    r"\s*<td[^>]*>(.*?)</td>\s*<td>(.*?)</td>\s*<td>(.*?)</td>",
    re.S,
)


def _text(s: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", s)).strip()


def parse_omix_html(raw: str) -> list[dict]:
    rows = []
    for omix, prj, title, org, access, rel in ROW_RE.findall(raw):
        rows.append(
            {
                "accession": omix.strip(),
                "bioproject": _text(prj),
                "title": _text(title),
                "organism": _text(org),
                "access": access.strip(),
                "release_date": rel.strip(),
                "url": f"https://ngdc.cncb.ac.cn/omix/release/{omix.strip()}",
            }
        )
    return rows


def parse_hra_json(raw: str) -> list[dict]:
    data = json.loads(raw)
    rows = []
    for x in data:
        rows.append(
            {
                "accession": x.get("accession") or "",
                "bioproject": x.get("bioproject") or "",
                "title": x.get("title") or "",
                "organization": x.get("uorganization") or x.get("organization") or "",
                "is_controlled": str(x.get("isControlledAccess")),
                "access": "Controlled" if str(x.get("isControlledAccess")) == "1" else "Open",
                "openaccess_data_type": str(x.get("openaccessDataType")),
                "release_time_string": x.get("lastTimeString") or "",
                "url": f"https://ngdc.cncb.ac.cn/gsa-human/browse/{x.get('accession')}",
            }
        )
    return rows


def write_tsv(path: str, rows: list[dict], fields: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.join(ROOT, "data", "catalogs"))
    args = ap.parse_args()
    setup_logging()

    print("Fetching OMIX releaseList (large HTML)...")
    omix = fetch(OMIX_URL, use_cache=True)
    if not omix["ok"]:
        raise SystemExit(f"OMIX releaseList failed: {omix}")
    omix_rows = parse_omix_html(omix["text"])
    print(f"  OMIX rows: {len(omix_rows)}")

    print("Fetching GSA-Human finished.json...")
    hra = fetch(HRA_URL, referer=HRA_REFERER, use_cache=True)
    if not hra["ok"]:
        raise SystemExit(f"HRA finished.json failed: {hra}")
    hra_rows = parse_hra_json(hra["text"])
    print(f"  HRA rows: {len(hra_rows)}")

    write_tsv(
        os.path.join(args.outdir, "omix_release.tsv"),
        omix_rows,
        ["accession", "bioproject", "title", "organism", "access", "release_date", "url"],
    )
    write_tsv(
        os.path.join(args.outdir, "hra_finished.tsv"),
        hra_rows,
        [
            "accession",
            "bioproject",
            "title",
            "organization",
            "is_controlled",
            "access",
            "openaccess_data_type",
            "release_time_string",
            "url",
        ],
    )
    print("Wrote", args.outdir)


if __name__ == "__main__":
    main()
