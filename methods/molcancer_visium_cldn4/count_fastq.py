#!/usr/bin/env python3
"""Stream-count official Visium FFPE probes from PRJNA1139087 ENA FASTQ.

No Space Ranger image step. Spatial barcodes are mapped with the Visium v1
whitelist + 1-edit correction. R2 (probe insert before polyA) is matched to
the official 10x Human Transcriptome Probe Set v1.0 sequences for a CLDN4 /
CD8 / epithelial panel.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
REFS = ROOT / "data" / "refs"
OUT = ROOT / "data" / "PRJNA1139087" / "counts"

TARGET_GENES = [
    "CLDN4",
    "CD8A",
    "CD8B",
    "KRT8",
    "KRT18",
    "KRT19",
    "EPCAM",
    "CDH1",
    "KRT7",
    "CD3D",
    "CD3E",
    "CD3G",
    "NKG7",
    "PTPRC",
]

SAMPLES = [
    {
        "patient": "PA10",
        "alias": "Tumor-2",
        "srr": "SRR29925400",
        "response": "pCR",
        "timepoint": "post_chemoIO",
        "r1": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/000/SRR29925400/SRR29925400_1.fastq.gz",
        "r2": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/000/SRR29925400/SRR29925400_2.fastq.gz",
    },
    {
        "patient": "PA08",
        "alias": "Tumor-1",
        "srr": "SRR29925401",
        "response": "NMPR",
        "timepoint": "post_chemoIO",
        "r1": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/001/SRR29925401/SRR29925401_1.fastq.gz",
        "r2": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/001/SRR29925401/SRR29925401_2.fastq.gz",
    },
    {
        "patient": "PA09",
        "alias": "Tumor-4",
        "srr": "SRR29925398",
        "response": "NMPR",
        "timepoint": "post_chemoIO",
        "r1": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/098/SRR29925398/SRR29925398_1.fastq.gz",
        "r2": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/098/SRR29925398/SRR29925398_2.fastq.gz",
    },
    {
        "patient": "PA12",
        "alias": "Tumor-3",
        "srr": "SRR29925399",
        "response": "NMPR_responsive_like",
        "timepoint": "post_chemoIO",
        "r1": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/099/SRR29925399/SRR29925399_1.fastq.gz",
        "r2": "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR299/099/SRR29925399/SRR29925399_2.fastq.gz",
    },
]

PROBE_URL = (
    "https://cf.10xgenomics.com/supp/spatial-exp/probeset/"
    "Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv"
)
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36"


def load_probes(path: Path) -> tuple[dict[bytes, str], dict[bytes, str]]:
    exact: dict[bytes, str] = {}
    pref: dict[bytes, str] = {}
    pref_amb: set[bytes] = set()
    want = set(TARGET_GENES)
    with path.open() as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("gene_id"):
                continue
            parts = line.rstrip("\n").split(",")
            if len(parts) < 4:
                continue
            gene_id, seq, probe_id, included = parts[:4]
            if included.upper() != "TRUE":
                continue
            gene = probe_id.split("|")[1] if "|" in probe_id else gene_id
            if gene not in want:
                continue
            b = seq.encode("ascii")
            exact[b] = gene
            p = b[:25]
            if p in pref and pref[p] != gene:
                pref_amb.add(p)
            else:
                pref[p] = gene
    for p in pref_amb:
        pref.pop(p, None)
    return exact, pref


def load_whitelist(coord_path: Path) -> tuple[set[str], dict[str, tuple[int, int]], dict[str, str]]:
    coords: dict[str, tuple[int, int]] = {}
    with coord_path.open() as fh:
        for line in fh:
            bc, c, r = line.split()
            coords[bc] = (int(r), int(c))  # array_row, array_col
    wl = set(coords)
    fix: dict[str, str] = {}
    for bc in wl:
        for i in range(len(bc)):
            for alt in "ACGT":
                if alt == bc[i]:
                    continue
                mut = bc[:i] + alt + bc[i + 1 :]
                if mut in wl:
                    continue
                if mut in fix and fix[mut] != bc:
                    fix[mut] = ""  # ambiguous
                else:
                    fix[mut] = bc
    fix = {k: v for k, v in fix.items() if v}
    return wl, coords, fix


def resolve_barcode(raw: str, wl: set[str], fix: dict[str, str]) -> str | None:
    if raw in wl:
        return raw
    if raw in fix:
        return fix[raw]
    if "N" in raw:
        idxs = [i for i, x in enumerate(raw) if x == "N"]
        if len(idxs) == 1:
            i = idxs[0]
            hits = []
            for alt in "ACGT":
                t = raw[:i] + alt + raw[i + 1 :]
                if t in wl:
                    hits.append(t)
            if len(hits) == 1:
                return hits[0]
    return None


def open_fastq_stream(url: str) -> subprocess.Popen:
    if Path(url).exists():
        cmd = f"gzip -dc '{url}'"
    else:
        cmd = f"curl -fsSL --retry 5 --retry-delay 4 '{url}' | gzip -dc"
    return subprocess.Popen(
        ["bash", "-lc", cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        bufsize=1024 * 1024,
    )


def aria2_get(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1_000_000:
        return dest
    cmd = [
        "aria2c",
        "-x",
        "16",
        "-s",
        "16",
        "-k",
        "1M",
        "--file-allocation=none",
        "--allow-overwrite=true",
        "-d",
        str(dest.parent),
        "-o",
        dest.name,
        url,
    ]
    subprocess.check_call(cmd)
    return dest


def iter_pairs(p1: subprocess.Popen, p2: subprocess.Popen):
    r1 = p1.stdout
    r2 = p2.stdout
    while True:
        h1 = r1.readline()
        if not h1:
            break
        s1 = r1.readline()
        r1.readline()
        r1.readline()
        h2 = r2.readline()
        s2 = r2.readline()
        r2.readline()
        r2.readline()
        if not s1 or not s2:
            break
        yield s1, s2


def probe_insert(seq: bytes) -> bytes:
    i = seq.find(b"AAAAAAAA")
    if i >= 20:
        return seq[:i]
    return seq[:50]


def match_gene(ins: bytes, exact: dict[bytes, str], pref: dict[bytes, str]) -> str | None:
    if ins in exact:
        return exact[ins]
    if len(ins) >= 50 and ins[:50] in exact:
        return exact[ins[:50]]
    if len(ins) >= 25:
        return pref.get(ins[:25])
    return None


def count_sample(sample: dict, exact, pref, wl, coords, fix, out_csv: Path) -> dict:
    print(f"[count] {sample['patient']} {sample['srr']}", flush=True)
    p1 = open_fastq_stream(sample["r1"])
    p2 = open_fastq_stream(sample["r2"])
    n_read = n_bc = n_gene = 0
    reads = defaultdict(int)
    umi_seen: set[tuple[str, str, str]] = set()
    gene_umi = defaultdict(lambda: defaultdict(int))
    try:
        for s1, s2 in iter_pairs(p1, p2):
            n_read += 1
            if n_read % 20_000_000 == 0:
                print(
                    f"  {sample['srr']} reads={n_read:,} bc={n_bc:,} gene={n_gene:,}",
                    flush=True,
                )
            raw = s1[:16].decode("ascii", "ignore")
            bc = resolve_barcode(raw, wl, fix)
            if bc is None:
                continue
            n_bc += 1
            reads[bc] += 1
            umi = s1[16:28].decode("ascii", "ignore")
            gene = match_gene(probe_insert(s2.rstrip(b"\n")), exact, pref)
            if gene is None:
                continue
            key = (bc, umi, gene)
            if key in umi_seen:
                continue
            umi_seen.add(key)
            gene_umi[bc][gene] += 1
            n_gene += 1
    finally:
        for p in (p1, p2):
            if p.poll() is None:
                p.terminate()
        p1.wait()
        p2.wait()
        if p1.returncode not in (0, None, -15) or p2.returncode not in (0, None, -15):
            err = (p1.stderr.read() if p1.stderr else b"") + (p2.stderr.read() if p2.stderr else b"")
            if p1.returncode not in (0, None, -15) and p2.returncode not in (0, None, -15):
                raise RuntimeError(f"fastq stream failed rc={p1.returncode},{p2.returncode} {err[:400]!r}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    genes = TARGET_GENES
    with out_csv.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["barcode", "array_row", "array_col", "n_reads", *genes]
        )
        for bc, n in reads.items():
            r, c = coords[bc]
            g = gene_umi.get(bc, {})
            w.writerow([bc, r, c, n, *[g.get(x, 0) for x in genes]])
    summary = {
        "patient": sample["patient"],
        "srr": sample["srr"],
        "n_fastq_reads": n_read,
        "n_whitelist_reads": n_bc,
        "n_target_umis": n_gene,
        "n_barcodes": len(reads),
        "out": str(out_csv),
    }
    print(summary, flush=True)
    return summary


def ensure_probe_csv() -> Path:
    dest = REFS / "Visium_Human_Transcriptome_Probe_Set_v1.0_GRCh38-2020-A.csv"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    subprocess.check_call(
        ["curl", "-fsSL", "-A", UA, "-o", str(dest), PROBE_URL]
    )
    return dest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--patient", action="append", default=None)
    ap.add_argument(
        "--aria2",
        action="store_true",
        help="Download FASTQ with aria2 (16 connections) then count locally and delete.",
    )
    args = ap.parse_args()
    probe = ensure_probe_csv()
    exact, pref = load_probes(probe)
    print(f"probes exact={len(exact)} prefix25={len(pref)}", flush=True)
    if "CLDN4" not in exact.values():
        raise SystemExit("CLDN4 probe missing")
    wl, coords, fix = load_whitelist(REFS / "visium-v1_coordinates.txt")
    wanted = set(args.patient) if args.patient else {s["patient"] for s in SAMPLES}
    summaries = []
    fq_dir = ROOT / "data" / "PRJNA1139087" / "fastq"
    for sample in SAMPLES:
        if sample["patient"] not in wanted:
            continue
        local = dict(sample)
        r1 = r2 = None
        if args.aria2:
            r1 = fq_dir / f"{sample['srr']}_1.fastq.gz"
            r2 = fq_dir / f"{sample['srr']}_2.fastq.gz"
            print(f"[aria2] {sample['srr']} R1", flush=True)
            aria2_get(sample["r1"], r1)
            print(f"[aria2] {sample['srr']} R2", flush=True)
            aria2_get(sample["r2"], r2)
            local["r1"] = str(r1)
            local["r2"] = str(r2)
        out = OUT / f"{sample['patient']}_spot_counts.csv"
        summaries.append(count_sample(local, exact, pref, wl, coords, fix, out))
        if args.aria2:
            for p in (r1, r2):
                if p is not None and p.exists():
                    p.unlink()
    import json

    (OUT / "count_summary.json").write_text(json.dumps(summaries, indent=2) + "\n")


if __name__ == "__main__":
    main()
