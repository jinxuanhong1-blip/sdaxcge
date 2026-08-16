#!/usr/bin/env python3
"""Inspect leftover series matrices: platforms, gene mapping, full characteristics.

GSE141479 and GSE99995 have ~40k-row embedded tables (typical microarray probe IDs).
GSE154286 has a 201-row table (targeted panel). This script:
  * dumps every characteristic key/value
  * maps probe IDs via the GEO platform table when needed
  * records whether TACSTD2 / CLDN4 are actually measured
"""
import gzip
import json
import re
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
RES = ROOT / "results" / "w200" / "GEO_2020"
DL = RES / "downloads"
CLIN = RES / "clinical"
CLIN.mkdir(parents=True, exist_ok=True)

TARGETS = {
    "TACSTD2": {"TACSTD2", "TROP2", "TROP-2", "GA733-1", "M1S1", "EGP-1"},
    "CLDN4": {"CLDN4", "CLAUDIN4", "CLAUDIN-4", "CPE-R", "CPER"},
}


def gunzip_text(path):
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        return gzip.decompress(raw).decode("utf-8", errors="replace")
    return raw.decode("utf-8", errors="replace")


def parse_matrix(path):
    text = gunzip_text(path)
    meta = defaultdict(list)
    table_header = None
    table_ids = []
    in_table = False
    platform = None
    for line in text.splitlines():
        if line.startswith("!Series_platform_id"):
            platform = line.split("\t", 1)[-1].strip().strip('"')
        if line.startswith("!") and not in_table:
            key, *vals = line.split("\t")
            meta[key].append([v.strip().strip('"') for v in vals])
        if line.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if line.startswith("!series_matrix_table_end"):
            break
        if in_table:
            if table_header is None:
                table_header = [c.strip().strip('"') for c in line.split("\t")]
                continue
            table_ids.append(line.split("\t", 1)[0].strip().strip('"'))
    return meta, platform, table_header, table_ids


def fetch(url, dest, timeout=90):
    if dest.exists() and dest.stat().st_size > 100:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "w200-GEO-2020/1.0"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                dest.write_bytes(r.read())
            return dest
        except Exception as e:  # noqa: BLE001
            print(f"  retry {attempt} {url} ({e})", file=sys.stderr)
            time.sleep(2 ** attempt)
    return None


def load_platform_map(gpl):
    """Download GPL annotation and return probe -> gene-symbol list."""
    # GPL files live under geo/platforms/GPLxxxxnnn/GPL####/annot/
    num = int(gpl.replace("GPL", ""))
    if num < 1000:
        prefix = "GPL" + str(num).zfill(3)[:1] + "nnn"
    elif num < 10000:
        prefix = "GPL" + str(num)[:2] + "nnn"
    else:
        prefix = "GPL" + str(num)[:3] + "nnn"
    # try annot table then family SOFT
    urls = [
        f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl}/annot/{gpl}.annot.gz",
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}&targ=self&form=text&view=full",
    ]
    dest = DL / "platforms" / f"{gpl}.annot.gz"
    p = fetch(urls[0], dest)
    mapping = {}
    if p and p.stat().st_size > 200:
        text = gunzip_text(p)
        # annot format: ID, Gene symbol, ...
        header = None
        in_table = False
        for line in text.splitlines():
            if line.startswith("#") or line.startswith("^") or line.startswith("!"):
                if line.startswith("!platform_table_begin"):
                    in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            cols = line.split("\t")
            if header is None:
                header = [c.strip() for c in cols]
                continue
            rec = dict(zip(header, cols))
            pid = rec.get("ID") or rec.get("ID_REF") or cols[0]
            # common symbol fields
            symbols = []
            for k in ("Gene symbol", "Gene Symbol", "Symbol", "GENE_SYMBOL",
                      "gene_assignment", "ILMN_Gene", "GeneName"):
                if rec.get(k):
                    symbols.append(rec[k])
            # Clariom D (GPL23126) stores symbols only in gene_assignment
            # ("NM_002353 // TACSTD2 // ..."). A first-pass that looks only
            # at a Gene-symbol column will falsely report zero hits.
            mapping[pid] = ";".join(symbols)
        if mapping:
            return mapping
    # fallback: SOFT via acc.cgi
    dest2 = DL / "platforms" / f"{gpl}.soft.txt"
    req = urllib.request.Request(urls[1], headers={"User-Agent": "w200-GEO-2020/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            dest2.write_bytes(r.read())
        text = dest2.read_text(errors="replace")
        header = None
        in_table = False
        for line in text.splitlines():
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            cols = line.split("\t")
            if header is None:
                header = cols
                continue
            rec = dict(zip(header, cols))
            pid = rec.get("ID") or cols[0]
            symbols = []
            for k, v in rec.items():
                if re.search(r"symbol|gene.?name|ILMN_Gene", k, re.I) and v:
                    symbols.append(v)
            mapping[pid] = ";".join(symbols)
    except Exception as e:  # noqa: BLE001
        print("  platform SOFT fail", gpl, e)
    return mapping


def find_targets(ids, mapping=None):
    hits = {g: [] for g in TARGETS}
    for pid in ids:
        fields = [pid]
        if mapping and pid in mapping:
            fields.append(mapping[pid])
        blob = " ".join(fields).upper()
        for gene, aliases in TARGETS.items():
            if any(re.search(rf"\b{re.escape(a)}\b", blob) for a in aliases):
                hits[gene].append({"probe": pid, "annotation": (mapping or {}).get(pid, "")})
    return hits


def characteristics_table(meta):
    """Build a sample x characteristic table from !Sample_* rows."""
    accs = None
    for key, blocks in meta.items():
        if key == "!Sample_geo_accession":
            accs = blocks[0]
            break
    if not accs:
        return None
    rows = [{"geo_accession": a} for a in accs]
    for key, blocks in meta.items():
        if not key.startswith("!Sample_"):
            continue
        short = key.replace("!Sample_", "")
        for i, block in enumerate(blocks):
            col = short if i == 0 else f"{short}_{i+1}"
            # characteristics often "key: value"
            if short.startswith("characteristics"):
                parsed = []
                keys_seen = []
                for cell in block:
                    if ": " in cell:
                        k, v = cell.split(": ", 1)
                        keys_seen.append(k)
                        parsed.append(v)
                    else:
                        parsed.append(cell)
                if keys_seen and len(set(keys_seen)) == 1:
                    col = keys_seen[0].replace(" ", "_")
                    block = parsed
                elif keys_seen:
                    # mixed keys in one row — keep raw
                    pass
                else:
                    block = parsed
            if len(block) != len(rows):
                continue
            for rec, val in zip(rows, block):
                rec[col] = val
    return rows


def inspect(acc):
    mpath = next((p for p in (DL / acc).glob("*series_matrix.txt.gz")), None)
    if not mpath:
        print(acc, "NO MATRIX")
        return None
    meta, platform, header, ids = parse_matrix(mpath)
    print(f"\n==== {acc} platform={platform} n_probes={len(ids)} n_samples={0 if not header else len(header)-1}")
    # series summary
    for k in ("!Series_title", "!Series_summary", "!Series_overall_design", "!Series_type", "!Series_sample_id"):
        if k in meta:
            val = meta[k][0][0] if meta[k][0] else ""
            print(f"  {k}: {val[:240]}")
    clin = characteristics_table(meta)
    if clin:
        import csv
        out = CLIN / f"{acc}_geo_characteristics.csv"
        cols = list(clin[0].keys())
        with out.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(clin)
        print(f"  wrote {out} ({len(clin)} samples, cols={cols})")
        # unique values for likely outcome columns
        for col in cols:
            vals = sorted({(r.get(col) or "")[:60] for r in clin})
            if 1 < len(vals) <= 20:
                print(f"    {col}: {vals}")

    mapping = None
    # If IDs look like gene symbols already
    direct = find_targets(ids, None)
    n_direct = sum(len(v) for v in direct.values())
    if n_direct == 0 and platform:
        print(f"  downloading platform map {platform} ...")
        mapping = load_platform_map(platform)
        print(f"  platform map size={len(mapping)}")
        mapped = find_targets(ids, mapping)
    else:
        mapped = direct
    for g, hits in mapped.items():
        print(f"  {g}: {len(hits)} probes {hits[:5]}")

    rec = {
        "accession": acc,
        "platform": platform,
        "n_probes": len(ids),
        "n_samples": (len(header) - 1) if header else 0,
        "target_hits": {g: hits for g, hits in mapped.items()},
        "characteristic_columns": list(clin[0].keys()) if clin else [],
    }
    (DL / acc / "inspect.json").write_text(json.dumps(rec, indent=2))
    return rec


def main():
    recs = []
    for acc in ["GSE141479", "GSE154286", "GSE150972", "GSE99995", "GSE124885"]:
        recs.append(inspect(acc))
        time.sleep(0.2)
    (RES / "tables" / "matrix_inspect.json").write_text(json.dumps(recs, indent=2))
    print("\nWrote matrix_inspect.json")


if __name__ == "__main__":
    main()
