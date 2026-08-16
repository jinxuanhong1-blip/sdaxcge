#!/usr/bin/env python3
"""Catalog public ELF3 ChIP and test locus overlap at TACSTD2 / CLDN4.

Sources:
  * ENCODE experiment search (target.label=ELF3)
  * ChIP-Atlas hg38 ELF3 target-gene matrix (5 kb) + MACS2 bed05 peaks
  * GEO esearch: ELF3 AND ChIP-seq AND Homo sapiens (lung filter is a
    post-hoc read of titles/summaries, not a second search used to hide hits)

There is no public lung / NSCLC / A549 ELF3 ChIP-seq in those catalogs.
Locus overlap is therefore reported only for the available non-lung
biosamples, and is not a lung-binding claim.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

UA = (
    "sdaxcge-w200-A10-ELF3/1.0 "
    "(reproducible public ChIP catalog; "
    "+https://github.com/jinxuanhong1-blip/sdaxcge)"
)

# Ensembl GRCh38 / hg38 (looked up 2026-08-16).
LOCI = {
    "TACSTD2": {
        "chrom": "chr1",
        "gene_start": 58575433,
        "gene_end": 58577252,
        "strand": -1,
        "ensembl": "ENSG00000184292",
        # Canonical protein-coding: ENST00000371225, single exon, minus.
        "coding_tss": 58577252,
        "notes": "Single-exon protein-coding gene; TSS = gene end.",
    },
    "CLDN4": {
        "chrom": "chr7",
        "gene_start": 73799542,
        "gene_end": 73832690,
        "strand": 1,
        "ensembl": "ENSG00000189143",
        # Canonical protein-coding: ENST00000340958 TSS=73830996.
        # Gene span starts at ENST00000466411 (protein_coding_CDS_not_defined).
        "coding_tss": 73830996,
        "alt_tss": 73799542,
        "notes": (
            "Gene span includes a ~31 kb upstream non-coding/retained isoform "
            "(CLDN4-204). Coding TSS is ENST00000340958."
        ),
    },
}

# ChIP-Atlas ELF3 columns (hg38 target matrix, 5 kb). Biosample notes are
# from GEO sample records, not inferred from the cell-line name alone.
CHIP_ATLAS_EXPS = [
    {
        "srx": "SRX1389378",
        "gsm": "GSM1919984",
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC-1 KLF5-KO ELF3 ChIP (GSE64557, PDAC).",
    },
    {
        "srx": "SRX1389380",
        "gsm": "GSM1919986",
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC-1 ELF3 ChIP (GSE64557, PDAC).",
    },
    {
        "srx": "SRX21154913",
        "gsm": None,
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC-1 ELF3 ChIP (ChIP-Atlas; PDAC).",
    },
    {
        "srx": "SRX21154914",
        "gsm": None,
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC-1 ELF3 ChIP (ChIP-Atlas; PDAC).",
    },
    {
        "srx": "SRX21154915",
        "gsm": None,
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC-1 ELF3 ChIP (ChIP-Atlas; PDAC).",
    },
    {
        "srx": "SRX21154916",
        "gsm": None,
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC-1 ELF3 ChIP (ChIP-Atlas; PDAC).",
    },
    {
        "srx": "SRX825394",
        "gsm": "GSM1574273",
        "cell": "CFPAC-1",
        "tissue": "pancreas",
        "lung": False,
        "note": "CFPAC1.ELF3 (GSE64557 / related PDAC ELF3 ChIP).",
    },
    {
        "srx": "SRX6060975",
        "gsm": "GSM3886298",
        "cell": "ESO-26",
        "tissue": "esophagus",
        "lung": False,
        "note": "Eso26 ELF3-Flag (GSE132680, esophageal adenocarcinoma).",
    },
    {
        "srx": "SRX8939702",
        "gsm": "GSM4725867",
        "cell": "HBDEC2",
        "tissue": "biliary_epithelium",
        "lung": False,
        "note": (
            "HBDEC2 is an immortalized biliary epithelial line "
            "(GSE156165), not bronchus/lung."
        ),
    },
    {
        "srx": "SRX2636285",
        "gsm": None,
        "cell": "HepG2",
        "tissue": "liver",
        "lung": False,
        "note": "HepG2 ELF3 (ChIP-Atlas; ENCODE-related).",
    },
    {
        "srx": "SRX2636286",
        "gsm": None,
        "cell": "HepG2",
        "tissue": "liver",
        "lung": False,
        "note": "HepG2 ELF3 (ChIP-Atlas; ENCODE-related).",
    },
]

ENCODE_EXP = {
    "accession": "ENCSR770AOR",
    "peak_file": "ENCFF080FAU",
    "peak_url": "https://www.encodeproject.org/files/ENCFF080FAU/@@download/ENCFF080FAU.bed.gz",
    "output_type": "optimal IDR thresholded peaks",
    "assembly": "GRCh38",
    "cell": "HepG2",
    "tissue": "liver",
    "lung": False,
    "note": "Only released ENCODE ELF3 TF ChIP-seq experiment (2026-08-16 search).",
}

CHIP_ATLAS_TARGETS = (
    "https://chip-atlas.dbcls.jp/data/hg38/target/ELF3.5.tsv"
)
CHIP_ATLAS_BED = (
    "https://chip-atlas.dbcls.jp/data/hg38/eachData/bed05/{srx}.05.bed"
)

WINDOWS_BP = {
    "gene_body": 0,
    "gene_pm10kb": 10_000,
    "coding_tss_pm2kb": 2_000,
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(url: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    rec = {"url": url, "path": str(dest)}
    if dest.exists() and dest.stat().st_size > 0:
        rec["status"] = "cached"
        rec["bytes"] = dest.stat().st_size
        rec["sha256"] = sha256_file(dest)
        print(f"[skip] {dest.name} ({rec['bytes']:,} bytes)", flush=True)
        return rec
    print(f"[get ] {url}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=180) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)
    rec["status"] = "downloaded"
    rec["bytes"] = dest.stat().st_size
    rec["sha256"] = sha256_file(dest)
    print(f"[done] {dest.name}  size={rec['bytes']:,}", flush=True)
    return rec


def window_interval(gene: str, kind: str) -> tuple[str, int, int]:
    loc = LOCI[gene]
    chrom = loc["chrom"]
    if kind == "gene_body":
        return chrom, loc["gene_start"], loc["gene_end"]
    if kind == "gene_pm10kb":
        return chrom, loc["gene_start"] - 10_000, loc["gene_end"] + 10_000
    if kind == "coding_tss_pm2kb":
        tss = loc["coding_tss"]
        return chrom, tss - 2_000, tss + 2_000
    if kind == "alt_tss_pm2kb":
        tss = loc["alt_tss"]
        return chrom, tss - 2_000, tss + 2_000
    raise KeyError(kind)


def overlaps(chrom: str, start: int, end: int, pchrom: str, pstart: int, pend: int) -> bool:
    return chrom == pchrom and pend > start and pstart < end


def read_bed(path: Path):
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            p = line.rstrip("\n").split("\t")
            yield p[0], int(p[1]), int(p[2]), p


def extract_chip_atlas_scores(tsv: Path, genes: list[str]) -> tuple[list[str], dict[str, list[str]]]:
    with tsv.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        wanted = {}
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if cols and cols[0] in genes:
                wanted[cols[0]] = cols
    return header, wanted


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--cache-dir", default="/tmp/a10_elf3_data")
    p.add_argument("--out-dir", default="results/w200/A10_ELF3")
    args = p.parse_args()
    cache = Path(args.cache_dir)
    out = Path(args.out_dir)
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    beds = cache / "chip_beds"
    beds.mkdir(exist_ok=True)

    manifest: dict = {
        "downloaded_at_utc": datetime.now(timezone.utc).isoformat(),
        "user_agent": UA,
        "files": [],
        "loci_hg38": LOCI,
    }

    # --- ENCODE peaks ---
    enc_path = cache / "ELF3_HepG2_ENCFF080FAU_optIDR.bed.gz"
    rec = fetch(ENCODE_EXP["peak_url"], enc_path)
    rec["name"] = enc_path.name
    manifest["files"].append(rec)

    # --- ChIP-Atlas target matrix ---
    tgt_path = cache / "ELF3_ChIPAtlas_targets_5kb.tsv"
    rec = fetch(CHIP_ATLAS_TARGETS, tgt_path)
    rec["name"] = tgt_path.name
    manifest["files"].append(rec)

    header, scores = extract_chip_atlas_scores(tgt_path, ["TACSTD2", "CLDN4"])
    score_rows = []
    # columns: Target_genes, ELF3|Average, SRX...|cell, ..., STRING
    exp_cols = [c for c in header[2:] if c != "STRING" and "|" in c]
    for gene in ("TACSTD2", "CLDN4"):
        cols = scores.get(gene)
        if not cols:
            continue
        by_hdr = dict(zip(header, cols))
        for col in exp_cols:
            srx, cell = col.split("|", 1)
            score_rows.append(
                {
                    "gene": gene,
                    "srx": srx,
                    "chip_atlas_cell": cell,
                    "macs2_score_or_zero": by_hdr.get(col, ""),
                    "chip_atlas_average": by_hdr.get("ELF3|Average", ""),
                    "window": "ChIP-Atlas official 5kb TSS/gene target call",
                }
            )
    with (out / "chip_atlas_gene_scores.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(score_rows[0].keys()))
        w.writeheader()
        w.writerows(score_rows)

    # --- ChIP-Atlas peak beds ---
    bed_paths: dict[str, Path] = {}
    for exp in CHIP_ATLAS_EXPS:
        dest = beds / f"{exp['srx']}.05.bed"
        rec = fetch(CHIP_ATLAS_BED.format(srx=exp["srx"]), dest)
        rec["name"] = dest.name
        manifest["files"].append(rec)
        bed_paths[exp["srx"]] = dest

    # --- Overlaps ---
    window_kinds = ["gene_body", "gene_pm10kb", "coding_tss_pm2kb"]
    overlap_rows = []
    encode_near = []

    def record_peaks(source: str, meta: dict, path: Path) -> None:
        peaks = list(read_bed(path))
        n_peaks = len(peaks)
        for gene, loc in LOCI.items():
            kinds = list(window_kinds)
            if "alt_tss" in loc:
                kinds.append("alt_tss_pm2kb")
            for kind in kinds:
                chrom, a, b = window_interval(gene, kind)
                hits = [pk for pk in peaks if overlaps(chrom, a, b, pk[0], pk[1], pk[2])]
                overlap_rows.append(
                    {
                        "source": source,
                        "srx_or_file": meta.get("srx") or meta.get("peak_file"),
                        "cell": meta["cell"],
                        "tissue": meta["tissue"],
                        "is_lung": meta["lung"],
                        "n_peaks_in_file": n_peaks,
                        "gene": gene,
                        "window": kind,
                        "window_chrom": chrom,
                        "window_start": a,
                        "window_end": b,
                        "n_overlapping_peaks": len(hits),
                        "peak_intervals": ";".join(f"{c}:{s}-{e}" for c, s, e, _ in hits),
                    }
                )
                if source.startswith("ENCODE") and hits:
                    for c, s, e, raw in hits:
                        encode_near.append(
                            raw
                            + [
                                gene,
                                kind,
                                f"{a}-{b}",
                            ]
                        )

    if enc_path.exists():
        record_peaks("ENCODE_HepG2_optIDR_ENCFF080FAU", ENCODE_EXP, enc_path)
    for exp in CHIP_ATLAS_EXPS:
        path = bed_paths[exp["srx"]]
        if path.exists() and path.stat().st_size > 0:
            record_peaks("ChIP-Atlas_bed05", exp, path)

    with (out / "chip_locus_overlap.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(overlap_rows[0].keys()))
        w.writeheader()
        w.writerows(overlap_rows)

    if encode_near:
        with (out / "encode_hepg2_peaks_near_loci.bed").open("w") as f:
            for raw in encode_near:
                f.write("\t".join(raw) + "\n")
    else:
        (out / "encode_hepg2_peaks_near_loci.bed").write_text(
            "# no ENCODE optimal-IDR ELF3 peaks in the requested windows\n"
        )

    # --- Experiment catalog ---
    exp_rows = [
        {
            "catalog": "ENCODE",
            "id": ENCODE_EXP["accession"],
            "peak_file": ENCODE_EXP["peak_file"],
            "cell": ENCODE_EXP["cell"],
            "tissue": ENCODE_EXP["tissue"],
            "is_lung": False,
            "note": ENCODE_EXP["note"],
        }
    ]
    for exp in CHIP_ATLAS_EXPS:
        exp_rows.append(
            {
                "catalog": "ChIP-Atlas",
                "id": exp["srx"],
                "peak_file": f"{exp['srx']}.05.bed",
                "cell": exp["cell"],
                "tissue": exp["tissue"],
                "is_lung": False,
                "note": exp["note"],
            }
        )
    # GEO survey (2026-08-16): ELF3 AND ChIP-seq AND Homo sapiens returned
    # 16 GDS/GEO records. None is a lung/NSCLC/A549 ELF3 ChIP. Extra
    # non-lung hit not in the ChIP-Atlas ELF3.5 matrix: HEK293 (GSE280165).
    exp_rows.append(
        {
            "catalog": "GEO_survey",
            "id": "GSE280165",
            "peak_file": "",
            "cell": "HEK293",
            "tissue": "kidney_cell_line",
            "is_lung": False,
            "note": (
                "2024 HEK293 ELF3 ChIP (genomic dark-matter TF screen). "
                "Not lung. Not used for locus overlap (no hg38 peak bed pulled)."
            ),
        }
    )
    exp_rows.append(
        {
            "catalog": "GEO_survey",
            "id": "GSE57177",
            "peak_file": "",
            "cell": "mouse_lung_tumor",
            "tissue": "lung_mouse",
            "is_lung": False,
            "note": (
                "GEO text hit because the series is lung ChIP-seq, but the "
                "target is not ELF3 (Smad4-loss mouse model). Not an ELF3 ChIP."
            ),
        }
    )
    with (out / "chip_experiments.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(exp_rows[0].keys()))
        w.writeheader()
        w.writerows(exp_rows)

    n_lung = sum(1 for r in exp_rows if r["is_lung"])
    # Collapse overlap to a readable call.
    def any_hit(source_prefix: str, gene: str, window: str, tissue: str | None = None) -> bool:
        for r in overlap_rows:
            if not r["source"].startswith(source_prefix):
                continue
            if r["gene"] != gene or r["window"] != window:
                continue
            if tissue and r["tissue"] != tissue:
                continue
            if r["n_overlapping_peaks"] > 0:
                return True
        return False

    chip_verdict = {
        "lung_elf3_chip_public": False,
        "n_lung_experiments_in_catalog": n_lung,
        "statement": (
            "No public human lung / NSCLC ELF3 ChIP-seq was found in ENCODE, "
            "ChIP-Atlas (ELF3.5 target matrix), or a GEO ELF3+ChIP-seq survey. "
            "ELF3 binding at TACSTD2/CLDN4 therefore cannot be tested in lung "
            "from public ChIP. Available non-lung ChIP (CFPAC-1 pancreas, "
            "HepG2 liver, ESO-26 esophagus, HBDEC2 biliary) is mixed and is "
            "not a lung result."
        ),
        "non_lung_calls": {
            "TACSTD2_coding_TSS_pm2kb_CFPAC1": any_hit(
                "ChIP-Atlas", "TACSTD2", "coding_tss_pm2kb", "pancreas"
            ),
            "TACSTD2_coding_TSS_pm2kb_HepG2_ENCODE_IDR": any_hit(
                "ENCODE", "TACSTD2", "coding_tss_pm2kb", "liver"
            ),
            "TACSTD2_coding_TSS_pm2kb_ESO26_or_HBDEC2": any_hit(
                "ChIP-Atlas", "TACSTD2", "coding_tss_pm2kb", "esophagus"
            )
            or any_hit("ChIP-Atlas", "TACSTD2", "coding_tss_pm2kb", "biliary_epithelium"),
            "CLDN4_coding_TSS_pm2kb_CFPAC1": any_hit(
                "ChIP-Atlas", "CLDN4", "coding_tss_pm2kb", "pancreas"
            ),
            "CLDN4_alt_TSS_pm2kb_CFPAC1": any_hit(
                "ChIP-Atlas", "CLDN4", "alt_tss_pm2kb", "pancreas"
            ),
            "CLDN4_alt_TSS_pm2kb_HepG2_ENCODE_IDR": any_hit(
                "ENCODE", "CLDN4", "alt_tss_pm2kb", "liver"
            ),
            "CLDN4_coding_TSS_pm2kb_HepG2_ENCODE_IDR": any_hit(
                "ENCODE", "CLDN4", "coding_tss_pm2kb", "liver"
            ),
        },
    }
    manifest["chip_verdict"] = chip_verdict
    (out / "chip_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "chip_verdict.json").write_text(json.dumps(chip_verdict, indent=2) + "\n")
    print(json.dumps(chip_verdict, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
