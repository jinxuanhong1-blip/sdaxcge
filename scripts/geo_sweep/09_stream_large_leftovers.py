#!/usr/bin/env python3
"""
Stream-extract TACSTD2 / CLDN4 rows from leftover ICI series whose processed
tables exceed the in-memory 300 MB cap but are still < 2 GB.

Does not invent values: only the matching gene row(s) are kept.
Updates notes/geo_sweep/analysis_supp.json in place.
"""
import gzip
import io
import json
import re
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

NOTES = Path("notes/geo_sweep")
SUPP = Path("results/geo_sweep/supp")
SUPP.mkdir(parents=True, exist_ok=True)

GENE_RE = re.compile(
    r"(^|[,;\t\" ])(TACSTD2|TROP2|CLDN4|Tacstd2|Cldn4|ENSMUSG00000051397|"
    r"ENSMUSG00000047501|ENSG00000184292|ENSG00000189143|"
    r"ENSMUST00000058178|ENSMUST00000051401)([,;\t\" \.]|$)",
)

TARGETS = {
    "GSE169688": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169688/suppl/"
        "GSE169688_LLC1.3Days.ExpressionMatrix.csv.gz",
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE169nnn/GSE169688/suppl/"
        "GSE169688_LLC1.1.5Days.Metadata.csv.gz",
    ],
    "GSE235122": [
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE235nnn/GSE235122/suppl/"
        "GSE235122_Matrix_processed_magisterh_melanoma_ICI.csv.gz",
    ],
}


def stream_hits(url, dest_hits, max_header_chars=2_000_000):
    """Download url, keep header + matching gene lines, discard the rest."""
    req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
    n_hit = 0
    n_line = 0
    header = None
    hits = []
    t0 = time.time()
    with urllib.request.urlopen(req, timeout=600) as r:
        with gzip.GzipFile(fileobj=r) as gz:
            for raw in gz:
                n_line += 1
                try:
                    ln = raw.decode("utf-8", "replace")
                except Exception:
                    continue
                if header is None:
                    header = ln
                    continue
                if GENE_RE.search(ln[:400]):
                    hits.append(ln)
                    n_hit += 1
                if n_line % 200000 == 0:
                    print(f"    scanned {n_line} lines, hits={n_hit} "
                          f"({time.time()-t0:.0f}s)", flush=True)
    dest_hits.write_text((header or "") + "".join(hits))
    print(f"  {url.rsplit('/',1)[-1]}: {n_line} lines, {n_hit} hits")
    return header, hits, n_line


def parse_hits(header, hits):
    if not header or not hits:
        return None
    text = header + "".join(hits)
    sample = header
    if sample.count("\t") >= sample.count(",") and "\t" in sample:
        sep = "\t"
    elif sample.count(";") > sample.count(","):
        sep = ";"
    else:
        sep = ","
    df = pd.read_csv(io.StringIO(text), sep=sep)
    return df


def main():
    listing = {o["accession"]: o for o in json.loads((NOTES / "analysis_supp.json").read_text())}
    for acc, urls in TARGETS.items():
        print(f"Streaming {acc} ...", flush=True)
        rec = listing.get(acc, {"accession": acc, "is_ici": "Y", "is_kd_ko": "N",
                                "files_tried": [], "markers": {}})
        rec.setdefault("files_tried", [])
        rec.setdefault("markers", {})
        for url in urls:
            name = url.rsplit("/", 1)[-1]
            if "Metadata" in name:
                dest = SUPP / name
                if not dest.exists():
                    req = urllib.request.Request(url, headers={"User-Agent": "geo-sweep/1.0"})
                    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as f:
                        f.write(r.read())
                rec["files_tried"].append({"name": name, "status": "metadata_saved"})
                continue
            dest = SUPP / (name.replace(".csv.gz", "") + ".marker_rows.csv")
            try:
                header, hits, n_line = stream_hits(url, dest)
            except Exception as e:  # noqa: BLE001
                rec["files_tried"].append({"name": name, "status": f"stream_error: {e}"})
                print("  ERROR", e)
                continue
            rec["files_tried"].append({
                "name": name, "status": "streamed",
                "n_lines_scanned": n_line, "n_hits": len(hits),
            })
            df = parse_hits(header, hits)
            if df is None or df.empty:
                continue
            # identify gene column / values
            for gene, aliases in {
                "TACSTD2": {"TACSTD2", "TROP2", "Tacstd2", "ENSMUSG00000051397",
                            "ENSMUST00000058178", "ENSG00000184292"},
                "CLDN4": {"CLDN4", "Cldn4", "ENSMUSG00000047501",
                          "ENSMUST00000051401", "ENSG00000189143"},
            }.items():
                mask = False
                for c in df.columns[:3]:
                    mask = mask | df[c].astype(str).str.replace(r"\.\d+$", "", regex=True).isin(aliases)
                # also first-col contains
                if not isinstance(mask, pd.Series):
                    mask = df.iloc[:, 0].astype(str).str.contains(
                        "|".join(aliases), case=False, regex=True)
                sub = df.loc[mask] if isinstance(mask, pd.Series) else df
                if sub.empty:
                    continue
                row = sub.iloc[0].to_dict()
                num = {k: float(v) for k, v in row.items()
                       if isinstance(v, (int, float, np.integer, np.floating)) and pd.notna(v)}
                rec["markers"][gene] = {
                    "file": name,
                    "is_de_table": False,
                    "streamed": True,
                    "n_rows_matched": int(sub.shape[0]),
                    "rows": [{k: (None if pd.isna(v) else
                                  (float(v) if isinstance(v, (int, float, np.integer, np.floating))
                                   else str(v)))
                              for k, v in row.items()}],
                    "per_sample_summary": {
                        "n_samples": len(num),
                        "n_values": len(num),
                        "mean": float(np.mean(list(num.values()))) if num else None,
                        "median": float(np.median(list(num.values()))) if num else None,
                        "min": float(np.min(list(num.values()))) if num else None,
                        "max": float(np.max(list(num.values()))) if num else None,
                    } if num else None,
                }
        rec["status"] = "markers_found" if rec["markers"] else rec.get("status", "markers_absent")
        listing[acc] = rec
        print(f"  {acc} -> {rec['status']} ({','.join(rec['markers']) or '-'})")

    out = list(listing.values())
    (NOTES / "analysis_supp.json").write_text(json.dumps(out, indent=2))
    print("Updated analysis_supp.json")


if __name__ == "__main__":
    main()
