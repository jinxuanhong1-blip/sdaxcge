#!/usr/bin/env python3
"""Gate: can scVelo be run on the concordant-4 public objects?

Concordant-4 is GSE123902 + GSE131907 + GSE205335 + GSE189357
(patient/sample n = 65). scVelo needs a spliced layer and an unspliced
layer (or labeled new/total counts). This script only inventories public
GEO supplements and SRA. It does not download count matrices, align reads,
or estimate velocity.

Outputs under results/:
  tables/geo_suppl_files.tsv
  tables/matrix_header_probe.tsv
  tables/sra_public_summary.tsv
  gate.json
"""

from __future__ import annotations

import csv
import datetime
import gzip
import io
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results"
TABLES = OUT / "tables"

SERIES = {
    "GSE123902": {
        "units": 13,
        "suppl": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123902/suppl/",
        "bioproject": "PRJNA510251",
        "sra_study": "SRP173552",
    },
    "GSE131907": {
        "units": 21,
        "suppl": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/",
        "bioproject": "PRJNA545296",
        "sra_study": None,
    },
    "GSE205335": {
        "units": 22,
        "suppl": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/",
        "bioproject": "PRJNA844398",
        "sra_study": None,
        "controlled_raw": "EGAD00001008703",
    },
    "GSE189357": {
        "units": 9,
        "suppl": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189357/suppl/",
        "bioproject": "PRJNA782639",
        "sra_study": "SRP347260",
    },
}

PROBES = [
    {
        "accession": "GSE123902",
        "what": "SEQC dense total-count CSV",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3516nnn/GSM3516666/suppl/GSM3516666_MSK_LX675_NORMAL_dense.csv.gz",
        "kind": "gzip_text",
        "n_lines": 2,
    },
    {
        "accession": "GSE131907",
        "what": "raw UMI matrix (genes x cells)",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131907/suppl/GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz",
        "kind": "gzip_text",
        "n_lines": 2,
    },
    {
        "accession": "GSE189357",
        "what": "Cell Ranger features.tsv",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5699nnn/GSM5699777/suppl/GSM5699777_TD1_features.tsv.gz",
        "kind": "gzip_text",
        "n_lines": 3,
    },
    {
        "accession": "GSE189357",
        "what": "Cell Ranger matrix.mtx",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5699nnn/GSM5699777/suppl/GSM5699777_TD1_matrix.mtx.gz",
        "kind": "gzip_text",
        "n_lines": 4,
    },
    {
        "accession": "GSE205335",
        "what": "cell identity table",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_CellIdentity.txt.gz",
        "kind": "gzip_text",
        "n_lines": 1,
    },
    {
        "accession": "GSE205335",
        "what": "UMI matrix RDS magic",
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE205nnn/GSE205335/suppl/GSE205335_Lung_IO_UMI_matrix.rds.gz",
        "kind": "gzip_magic",
        "n_bytes": 32,
    },
]

VELOCITY_NAME = re.compile(r"loom|h5ad|spliced|unsplic|velocy|intron", re.I)


def fetch(url: str, timeout: int = 120) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": "concordant4-velocity-gate/1.0"})
    last = None
    for i in range(4):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"GET failed {url}: {last}")


def esearch_count(db: str, term: str) -> int:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db="
        + db
        + "&term="
        + urllib.parse.quote(term)
        + "&retmode=json&retmax=0"
    )
    payload = json.loads(fetch(url).decode())
    return int(payload["esearchresult"]["count"])


def esearch_ids(db: str, term: str, retmax: int = 200) -> list[str]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db="
        + db
        + "&term="
        + urllib.parse.quote(term)
        + "&retmode=json&retmax="
        + str(retmax)
    )
    payload = json.loads(fetch(url).decode())
    res = payload["esearchresult"]
    if int(res["count"]) > len(res["idlist"]):
        raise RuntimeError(f"truncated id list for {term}: count={res['count']}")
    return res["idlist"]


def list_suppl(url: str) -> list[dict]:
    html = fetch(url).decode("utf-8", "replace")
    rows = []
    # NCBI FTP index is a <pre> of Apache links, not an HTML table.
    pat = re.compile(
        r'<a href="([^"]+)">[^<]*</a>\s+(\d{4}-\d{2}-\d{2} \d{2}:\d{2})\s+(\S+)'
    )
    for name, modified, size in pat.findall(html):
        if name in ("../", "/") or name.startswith("/") or name.startswith("?"):
            continue
        rows.append(
            {
                "file": urllib.parse.unquote(name),
                "modified": modified.strip(),
                "size": size.strip(),
                "velocity_name": bool(VELOCITY_NAME.search(name)),
            }
        )
    return rows


def clip(text: str, n: int = 180) -> str:
    text = " ".join(text.split())
    return text[:n]


def probe_streaming(spec: dict, max_bytes: int = 250_000) -> dict:
    """Read only the start of a remote file. Multi-GB matrices stay on GEO."""
    req = urllib.request.Request(spec["url"], headers={"User-Agent": "concordant4-velocity-gate/1.0"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        chunk = resp.read(max_bytes)
    out = {
        "accession": spec["accession"],
        "what": spec["what"],
        "url": spec["url"],
        "download_bytes": len(chunk),
        "gzip": chunk[:2] == b"\x1f\x8b",
        "hdf5_magic": chunk[:4] == b"\x89HDF",
        "excerpt": "",
        "single_matrix": "",
        "note": "prefix_only",
    }
    raw_txt = chunk
    wraps = 0
    try:
        while raw_txt[:2] == b"\x1f\x8b" and wraps < 3:
            raw_txt = gzip.GzipFile(fileobj=io.BytesIO(raw_txt)).read(8000)
            wraps += 1
    except EOFError:
        out["note"] = "prefix_only_truncated_gzip"
    out["gzip_wraps"] = wraps
    if spec["kind"] == "gzip_magic":
        out["excerpt"] = repr(raw_txt[:32])
        out["note"] = "rds_xdr_header" if raw_txt[:2] in (b"X\n", b"A\n", b"B\n") else "not_rds_header"
        out["single_matrix"] = "R_serialized_object_not_a_loom"
        return out
    text = raw_txt.decode("utf-8", "replace")
    lines = text.splitlines()[: spec["n_lines"]] or [text[:200]]
    out["excerpt"] = " || ".join(clip(ln) for ln in lines if ln)
    if "MatrixMarket matrix coordinate" in text:
        out["single_matrix"] = "one_matrixmarket_matrix"
        out["note"] = "cellranger_gene_expression_not_three_velocity_layers"
    elif text.startswith(","):
        out["single_matrix"] = "one_dense_gene_table"
    elif text.startswith("Index\t"):
        out["single_matrix"] = "one_umi_table_header_is_cell_barcodes"
    elif "Gene Expression" in text:
        out["single_matrix"] = "cellranger_feature_type_Gene_Expression"
    elif text.startswith("barcode\t"):
        out["single_matrix"] = "cell_metadata_not_counts"
    return out


def sra_summary(bioproject: str) -> dict:
    term = f"{bioproject}[BioProject]"
    count = esearch_count("sra", term)
    row = {
        "bioproject": bioproject,
        "sra_experiment_hits": count,
        "runinfo_rows": 0,
        "public_runs": 0,
        "size_mb": 0,
        "layouts": "",
        "strategies": "",
        "models": "",
        "avg_lengths": "",
        "path_kind": "",
        "aligned_bam_in_path": False,
    }
    if count == 0:
        return row
    ids = esearch_ids("sra", term, retmax=max(count, 1))
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id="
        + ",".join(ids)
        + "&rettype=runinfo&retmode=text"
    )
    text = fetch(url, timeout=180).decode()
    runs = list(csv.DictReader(io.StringIO(text)))
    row["runinfo_rows"] = len(runs)
    row["public_runs"] = sum(1 for r in runs if r.get("Consent") == "public")
    row["size_mb"] = sum(int(r.get("size_MB") or 0) for r in runs)
    row["layouts"] = ",".join(sorted({r.get("LibraryLayout", "") for r in runs}))
    row["strategies"] = ",".join(sorted({r.get("LibraryStrategy", "") for r in runs}))
    row["models"] = ",".join(sorted({r.get("Model", "") for r in runs}))
    row["avg_lengths"] = ",".join(sorted({r.get("avgLength", "") for r in runs}))
    joined = " ".join((r.get("download_path") or "") for r in runs).lower()
    row["aligned_bam_in_path"] = ".bam" in joined
    if "sralite" in joined or ".lite" in joined or joined.endswith(".sra"):
        row["path_kind"] = "sralite_reads_not_bam"
    elif ".fastq" in joined or ".fq" in joined:
        row["path_kind"] = "fastq"
    elif ".bam" in joined:
        row["path_kind"] = "bam"
    else:
        row["path_kind"] = "unspecified"
    return row


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    TABLES.mkdir(parents=True, exist_ok=True)
    suppl_rows = []
    for acc, meta in SERIES.items():
        files = list_suppl(meta["suppl"])
        if not files:
            raise RuntimeError(f"no supplementary files parsed for {acc}")
        for f in files:
            suppl_rows.append({"accession": acc, "units_in_concordant4": meta["units"], **f})
    write_tsv(TABLES / "geo_suppl_files.tsv", suppl_rows)

    probe_rows = []
    for spec in PROBES:
        # huge matrices: stream a prefix. small files also fine this way.
        probe_rows.append(probe_streaming(spec))
        time.sleep(0.3)
    write_tsv(TABLES / "matrix_header_probe.tsv", probe_rows)

    sra_rows = []
    for acc, meta in SERIES.items():
        summary = sra_summary(meta["bioproject"])
        summary["accession"] = acc
        summary["sra_study"] = meta.get("sra_study") or ""
        summary["controlled_raw"] = meta.get("controlled_raw") or ""
        sra_rows.append(summary)
        time.sleep(0.4)
    # stable column order
    sra_fields = [
        "accession",
        "bioproject",
        "sra_study",
        "controlled_raw",
        "sra_experiment_hits",
        "runinfo_rows",
        "public_runs",
        "size_mb",
        "layouts",
        "strategies",
        "models",
        "avg_lengths",
        "path_kind",
        "aligned_bam_in_path",
    ]
    sra_rows = [{k: row[k] for k in sra_fields} for row in sra_rows]
    write_tsv(TABLES / "sra_public_summary.tsv", sra_rows)

    velocity_named = [r for r in suppl_rows if r["velocity_name"]]
    any_bam = any(r["aligned_bam_in_path"] for r in sra_rows)
    public_reads = {r["accession"]: r["sra_experiment_hits"] > 0 for r in sra_rows}
    gate = {
        "question": "Do concordant-4 epithelial cells have an RNA-velocity flow toward a CLDN4-high barrier state?",
        "method": "scVelo (steady-state or dynamical) if spliced and unspliced counts are public",
        "scvelo_run": False,
        "reason": "no_public_spliced_unspliced_layers",
        "velocity_named_supplement_files": [r["file"] for r in velocity_named],
        "aligned_bam_in_public_sra_paths": any_bam,
        "public_fastq_bioproject": public_reads,
        "substitutes_not_used": [
            "graph pseudotime is not RNA velocity",
            "expression-only vector fields (scTour) are not splicing kinetics",
            "realigning public FASTQ would drop the two accessions with no public reads and would not be an available layer",
        ],
        "biological_answer": "not_estimable",
        "locked_tnk_result_retouched": False,
        "retrieved_at_utc": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
    }
    (OUT / "gate.json").write_text(json.dumps(gate, indent=2) + "\n")
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
