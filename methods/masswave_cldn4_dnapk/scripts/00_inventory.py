#!/usr/bin/env python3
"""Inventory PRIDE and MassIVE for CLDN4 knockdown or claudin-4 interactome datasets.

Protein tables only. RAW spectra are not downloaded.
MassIVE keyword search on the ProXI endpoint ignores the keyword and returns an
unfiltered page, so this script downloads the public datasets_json catalog and
filters titles and descriptions locally.
"""

from __future__ import annotations

import csv
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

UA = "cldn4-masswave/1.0 (public-catalog; protein-tables; no-raw)"
ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
NOTES = ROOT / "notes"
RESULTS.mkdir(parents=True, exist_ok=True)
NOTES.mkdir(parents=True, exist_ok=True)

PRIDE_KEYWORDS = [
    "CLDN4",
    "claudin-4",
    "Cldn4",
    "claudin",
    "claudin knockdown",
    "claudin interactome",
    "BioID claudin",
]
CLDN4_RE = re.compile(r"cldn\s*-?\s*4\b|claudin\s*-?\s*4\b", re.I)
CLAUDIN_RE = re.compile(r"claudin|\bcldn\s*-?\s*\d+\b", re.I)
KD_RE = re.compile(
    r"knock-?down|knock-?out|\bKO\b|\bKD\b|silencing|siRNA|shRNA|CRISPR|deficient",
    re.I,
)
INTERACT_RE = re.compile(
    r"interactome|co-?immunoprecipitat|CoIP|AP-MS|BioID|Bio-ID|proximity|APEX|TurboID|pull-?down",
    re.I,
)


def log(msg: str) -> None:
    print(msg, flush=True)


def http_bytes(url: str, timeout: int = 120, retries: int = 3) -> bytes:
    last = None
    for i in range(retries):
        req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(0.8 * (i + 1))
    raise RuntimeError(f"GET failed {url}: {last}")


def http_json(url: str, timeout: int = 90):
    raw = http_bytes(url, timeout=timeout)
    return json.loads(raw.decode("utf-8", errors="replace"))


def pride_keyword(keyword: str, page_size: int = 100) -> list[dict]:
    hits = []
    page = 0
    while page < 10:
        qs = urllib.parse.urlencode({"keyword": keyword, "pageSize": page_size, "page": page})
        data = http_json(f"https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?{qs}") or []
        if not isinstance(data, list) or not data:
            break
        hits.extend(data)
        if len(data) < page_size:
            break
        page += 1
    return hits


def pride_project(accession: str) -> dict:
    return http_json(f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}") or {}


def pride_filenames(accession: str) -> list[str]:
    names = []
    page = 0
    while page < 12:
        url = (
            "https://www.ebi.ac.uk/pride/ws/archive/v2/projects/"
            f"{accession}/files?pageSize=100&page={page}"
        )
        data = http_json(url) or []
        if not isinstance(data, list) or not data:
            break
        names.extend(f.get("fileName") or "" for f in data)
        if len(data) < 100:
            break
        page += 1
    return names


def classify(text: str) -> str:
    cldn4 = bool(CLDN4_RE.search(text))
    claudin = bool(CLAUDIN_RE.search(text))
    kd = bool(KD_RE.search(text))
    interact = bool(INTERACT_RE.search(text))
    if cldn4 and interact:
        return "cldn4_interactome"
    if cldn4 and kd:
        return "cldn4_knockdown"
    if cldn4:
        return "cldn4_other"
    if claudin and kd:
        return "other_claudin_knockdown"
    if claudin and interact:
        return "other_claudin_interactome"
    if claudin:
        return "other_claudin_mention"
    return "keyword_false_positive"


def main() -> None:
    stamp = datetime.now(timezone.utc).isoformat()
    query_log = []
    seen: dict[str, dict] = {}
    for kw in PRIDE_KEYWORDS:
        try:
            hits = pride_keyword(kw)
            err = ""
        except Exception as exc:  # noqa: BLE001
            hits = []
            err = str(exc)
        query_log.append({"api": "pride", "query": kw, "n": len(hits), "error": err})
        log(f"PRIDE {kw!r}: {len(hits)} {err}")
        for rec in hits:
            acc = rec.get("accession")
            if not acc:
                continue
            item = seen.setdefault(acc, {"queries": set(), "title": rec.get("title") or ""})
            item["queries"].add(kw)
            item["title"] = item["title"] or rec.get("title") or ""

    pride_rows = []
    for acc, item in sorted(seen.items()):
        try:
            proj = pride_project(acc)
        except Exception as exc:  # noqa: BLE001
            proj = {}
            log(f"  project fail {acc}: {exc}")
        title = proj.get("title") or item["title"]
        desc = proj.get("projectDescription") or ""
        blob = f"{title}\n{desc}"
        organisms = []
        for org in proj.get("organisms") or []:
            organisms.append(org.get("name") if isinstance(org, dict) else str(org))
        decision = classify(blob)
        cldn4_files = ""
        # Project text often says "claudin family" and not "CLDN4". The CLDN4
        # CoIP is a deposited SEARCH filename.
        if decision in {"other_claudin_interactome", "cldn4_other", "cldn4_interactome"}:
            try:
                names = pride_filenames(acc)
            except Exception as exc:  # noqa: BLE001
                names = []
                log(f"  files fail {acc}: {exc}")
            cldn_names = [n for n in names if re.search(r"cldn\s*-?\s*4", n, re.I)]
            cldn4_files = ";".join(cldn_names)
            if cldn_names and decision != "cldn4_knockdown":
                decision = "cldn4_interactome"
        pride_rows.append(
            {
                "accession": acc,
                "decision": decision,
                "in_scope": decision in {"cldn4_interactome", "cldn4_knockdown"},
                "cldn4_search_files": cldn4_files,
                "title": title,
                "organisms": "; ".join(organisms),
                "queries": "; ".join(sorted(item["queries"])),
                "description": desc.replace("\n", " ")[:800],
                "pride_url": f"https://www.ebi.ac.uk/pride/archive/projects/{acc}",
            }
        )
        time.sleep(0.05)

    pride_path = RESULTS / "inventory_pride.tsv"
    with pride_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(pride_rows[0].keys()) if pride_rows else ["accession"])
        writer.writeheader()
        writer.writerows(pride_rows)
    log(f"wrote {pride_path} n={len(pride_rows)}")

    massive_url = "https://massive.ucsd.edu/ProteoSAFe/datasets_json.jsp"
    log("downloading MassIVE catalog")
    raw = http_bytes(massive_url, timeout=180).decode("utf-8", errors="replace")
    data = json.loads(raw[raw.find("{") :])
    datasets = data.get("datasets") or []
    query_log.append({"api": "massive_catalog", "query": massive_url, "n": len(datasets), "error": ""})
    massive_rows = []
    for rec in datasets:
        title = rec.get("title") or ""
        desc = rec.get("description") or ""
        blob = f"{title}\n{desc}"
        if not CLAUDIN_RE.search(blob):
            continue
        decision = classify(blob)
        massive_rows.append(
            {
                "accession": rec.get("dataset") or "",
                "decision": decision,
                "in_scope": decision in {"cldn4_interactome", "cldn4_knockdown"},
                "mentions_cldn4": bool(CLDN4_RE.search(blob)),
                "title": title.replace("\n", " ")[:300],
                "description": desc.replace("\n", " ")[:500],
                "massive_url": f"https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession={rec.get('dataset')}",
            }
        )
    massive_path = RESULTS / "inventory_massive.tsv"
    fields = [
        "accession",
        "decision",
        "in_scope",
        "mentions_cldn4",
        "title",
        "description",
        "massive_url",
    ]
    with massive_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(massive_rows)
    log(f"wrote {massive_path} n={len(massive_rows)} from catalog {len(datasets)}")

    # ProXI keyword probe: documents that the parameter is ignored.
    proxi_note = {}
    try:
        probe = http_json(
            "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets?keywords=CLDN4&resultType=compact&pageSize=1"
        )
        first = ""
        if isinstance(probe, list) and probe:
            for acc in probe[0].get("accession") or []:
                if acc.get("value", "").startswith("MSV"):
                    first = acc.get("value")
        proxi_note = {
            "url": "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets?keywords=CLDN4&resultType=compact",
            "first_accession": first,
            "keyword_ignored": first not in {"", None} and "CLDN4" not in json.dumps(probe)[:500],
        }
    except Exception as exc:  # noqa: BLE001
        proxi_note = {"error": str(exc)}

    summary = {
        "generated_at": stamp,
        "pride_projects": len(pride_rows),
        "pride_in_scope": [r["accession"] for r in pride_rows if r["in_scope"]],
        "pride_decisions": {},
        "massive_catalog_n": len(datasets),
        "massive_claudin_mentions": len(massive_rows),
        "massive_in_scope": [r["accession"] for r in massive_rows if r["in_scope"]],
        "massive_cldn4_mentions": [r["accession"] for r in massive_rows if r["mentions_cldn4"]],
        "proxi_probe": proxi_note,
    }
    for row in pride_rows:
        summary["pride_decisions"][row["decision"]] = summary["pride_decisions"].get(row["decision"], 0) + 1
    (RESULTS / "inventory_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    (RESULTS / "search_log.json").write_text(
        json.dumps({"generated_at": stamp, "queries": query_log}, indent=2) + "\n"
    )
    log(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
