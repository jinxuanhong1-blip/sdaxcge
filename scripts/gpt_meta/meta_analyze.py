#!/usr/bin/env python3
"""Reproduce the TACSTD2/CLDN4 ICI-response meta-analysis.

Only Python's standard library is required. GEO files are downloaded into
results/gpt_meta/raw when absent. All generated files remain in the three
gpt_meta directories requested for this analysis.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
import statistics
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "results/gpt_meta/raw"
OUT = ROOT / "results/gpt_meta"
GENES = ("TACSTD2", "CLDN4")
ENSEMBL = {"TACSTD2": "ENSG00000184292", "CLDN4": "ENSG00000189143"}

URLS = {
    "GSE126044_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/soft/GSE126044_family.soft.gz",
    "GSE126044_counts.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126044/suppl/GSE126044_counts.txt.gz",
    "GSE135222_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/soft/GSE135222_family.soft.gz",
    "GSE135222_expression.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/suppl/GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
    "GSE136961_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/soft/GSE136961_family.soft.gz",
    "GSE136961_TPM.tsv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136961/suppl/GSE136961_TPM.tsv.gz",
    "GSE166449_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/soft/GSE166449_family.soft.gz",
    "GSE166449_TPM.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166449/suppl/GSE166449_Raw_gene_TPM_matrix.txt.gz",
}


def download() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        path = RAW / name
        if not path.exists():
            print(f"Downloading {name}")
            req = urllib.request.Request(url, headers={"User-Agent": "gpt-meta-repro/1.0"})
            with urllib.request.urlopen(req, timeout=120) as response, path.open("wb") as out:
                while chunk := response.read(1024 * 1024):
                    out.write(chunk)


def parse_soft(accession: str) -> list[dict[str, object]]:
    samples: list[dict[str, object]] = []
    current: dict[str, object] | None = None
    with gzip.open(RAW / f"{accession}_family.soft.gz", "rt", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                current = {"gsm": line.split(" = ", 1)[1], "characteristics": {}}
            elif current and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif current and line.startswith("!Sample_description = "):
                current["description"] = line.split(" = ", 1)[1]
            elif current and line.startswith("!Sample_characteristics_ch1 = "):
                value = line.split(" = ", 1)[1]
                if ": " in value:
                    key, val = value.split(": ", 1)
                    cast = current["characteristics"]
                    assert isinstance(cast, dict)
                    cast[key.strip().lower()] = val.strip()
        if current:
            samples.append(current)
    return samples


def read_matrix(name: str) -> tuple[list[str], dict[str, list[float]]]:
    rows: dict[str, list[float]] = {}
    with gzip.open(RAW / name, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        samples = [x.strip('"') for x in header[1:]]
        for row in reader:
            if not row:
                continue
            rows[row[0].strip('"')] = [float(x) for x in row[1:]]
    return samples, rows


def labels_126044() -> dict[str, str]:
    result = {}
    for sample in parse_soft("GSE126044"):
        title = str(sample["title"]).removeprefix("RNA-seq_")
        chars = sample["characteristics"]
        assert isinstance(chars, dict)
        label = chars["patient response"]
        result[title] = "R" if label == "responder" else "NR"
    return result


def labels_135222() -> dict[str, str]:
    """Classify durable benefit from GEO PFS using the published 6-month rule."""
    result = {}
    for sample in parse_soft("GSE135222"):
        title = str(sample["title"]).replace(" ", "")
        chars = sample["characteristics"]
        assert isinstance(chars, dict)
        pfs_days = float(chars["pfs.time"])
        result[title] = "R" if pfs_days >= 180 else "NR"
    return result


def labels_166449() -> dict[str, str]:
    result = {}
    for sample in parse_soft("GSE166449"):
        title = str(sample["title"]).lower()
        description = str(sample["description"])
        result[description] = "NR" if "nonresponder" in title else "R"
    return result


def sample_rows() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    long_rows: list[dict[str, object]] = []
    availability: list[dict[str, object]] = []

    # GSE126044: raw counts -> log2(CPM + 1).
    samples, matrix = read_matrix("GSE126044_counts.txt.gz")
    labels = labels_126044()
    library_sizes = [sum(values[i] for values in matrix.values()) for i in range(len(samples))]
    for gene in GENES:
        present = gene in matrix
        availability.append({"cohort": "GSE126044", "gene": gene, "available": present})
        if present:
            for i, sample in enumerate(samples):
                value = math.log2(matrix[gene][i] / library_sizes[i] * 1_000_000 + 1)
                long_rows.append(row("GSE126044", sample, labels[sample], gene, value, "log2(CPM+1)"))

    # GSE135222: supplied TPM; Ensembl IDs are explicitly mapped above.
    samples, matrix = read_matrix("GSE135222_expression.tsv.gz")
    labels = labels_135222()
    stripped = {key.split(".", 1)[0]: values for key, values in matrix.items()}
    for gene in GENES:
        gene_id = ENSEMBL[gene]
        present = gene_id in stripped
        availability.append({"cohort": "GSE135222", "gene": gene, "available": present})
        if present:
            for i, sample in enumerate(samples):
                value = math.log2(stripped[gene_id][i] + 1)
                long_rows.append(row("GSE135222", sample, labels[sample], gene, value, "log2(TPM+1)"))

    # GSE136961 is a targeted 395-gene panel. Record non-availability.
    _, matrix = read_matrix("GSE136961_TPM.tsv.gz")
    symbols = {key.rsplit("_", 1)[0] for key in matrix}
    for gene in GENES:
        availability.append({"cohort": "GSE136961", "gene": gene, "available": gene in symbols})

    # GSE166449 values are already log-scale in the deposited matrix.
    samples, matrix = read_matrix("GSE166449_TPM.txt.gz")
    labels = labels_166449()
    for gene in GENES:
        present = gene in matrix
        availability.append({"cohort": "GSE166449", "gene": gene, "available": present})
        if present:
            for i, sample in enumerate(samples):
                long_rows.append(row("GSE166449", sample, labels[sample], gene, matrix[gene][i], "deposited log-scale TPM"))

    return long_rows, availability


def row(cohort: str, sample: str, response: str, gene: str, value: float, scale: str) -> dict[str, object]:
    return {
        "cohort": cohort,
        "sample": sample,
        "response": response,
        "gene": gene,
        "expression": value,
        "analysis_scale": scale,
    }


def hedges_g(responder: list[float], nonresponder: list[float]) -> dict[str, float]:
    n_r, n_nr = len(responder), len(nonresponder)
    mean_r, mean_nr = statistics.mean(responder), statistics.mean(nonresponder)
    sd_r, sd_nr = statistics.stdev(responder), statistics.stdev(nonresponder)
    df = n_r + n_nr - 2
    pooled = math.sqrt(((n_r - 1) * sd_r**2 + (n_nr - 1) * sd_nr**2) / df)
    d = (mean_r - mean_nr) / pooled
    correction = 1 - 3 / (4 * df - 1)
    g = correction * d
    variance = (n_r + n_nr) / (n_r * n_nr) + g**2 / (2 * df)
    se = math.sqrt(variance)
    z = g / se
    return {
        "n_responder": n_r,
        "n_nonresponder": n_nr,
        "mean_responder": mean_r,
        "mean_nonresponder": mean_nr,
        "hedges_g": g,
        "se": se,
        "variance": variance,
        "ci_low": g - 1.96 * se,
        "ci_high": g + 1.96 * se,
        "p_value": math.erfc(abs(z) / math.sqrt(2)),
    }


def cohort_effects(long_rows: list[dict[str, object]], availability: list[dict[str, object]]) -> list[dict[str, object]]:
    effects: list[dict[str, object]] = []
    for item in availability:
        cohort, gene = str(item["cohort"]), str(item["gene"])
        matching = [x for x in long_rows if x["cohort"] == cohort and x["gene"] == gene]
        if not matching:
            effects.append({
                "cohort": cohort,
                "gene": gene,
                "status": "not estimable: gene absent from deposited matrix",
            })
            continue
        responders = [float(x["expression"]) for x in matching if x["response"] == "R"]
        nonresponders = [float(x["expression"]) for x in matching if x["response"] == "NR"]
        stats = hedges_g(responders, nonresponders)
        effects.append({"cohort": cohort, "gene": gene, "status": "estimated", **stats})
    return effects


def pool(items: list[dict[str, object]], random: bool) -> dict[str, float]:
    effects = [float(x["hedges_g"]) for x in items]
    variances = [float(x["variance"]) for x in items]
    fixed_weights = [1 / v for v in variances]
    fixed_mean = sum(w * y for w, y in zip(fixed_weights, effects)) / sum(fixed_weights)
    q = sum(w * (y - fixed_mean) ** 2 for w, y in zip(fixed_weights, effects))
    df = len(items) - 1
    c = sum(fixed_weights) - sum(w**2 for w in fixed_weights) / sum(fixed_weights)
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    weights = [1 / (v + tau2) for v in variances] if random else fixed_weights
    estimate = sum(w * y for w, y in zip(weights, effects)) / sum(weights)
    se = math.sqrt(1 / sum(weights))
    z = estimate / se
    return {
        "k": len(items),
        "estimate": estimate,
        "se": se,
        "ci_low": estimate - 1.96 * se,
        "ci_high": estimate + 1.96 * se,
        "p_value": math.erfc(abs(z) / math.sqrt(2)),
        "Q": q,
        "Q_df": df,
        "I2_percent": max(0.0, (q - df) / q * 100) if q > 0 else 0.0,
        "tau2_DL": tau2,
    }


def meta_summary(effects: list[dict[str, object]]) -> list[dict[str, object]]:
    result = []
    for gene in GENES:
        items = [x for x in effects if x["gene"] == gene and x["status"] == "estimated"]
        for model, is_random in (("fixed", False), ("random_DL", True)):
            result.append({"gene": gene, "model": model, **pool(items, is_random)})
    return result


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for record in rows:
            writer.writerow(record)


def write_manifest() -> None:
    rows = []
    for name, url in URLS.items():
        content = (RAW / name).read_bytes()
        rows.append({
            "file": name,
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
            "url": url,
        })
    write_tsv(OUT / "source_manifest.tsv", rows, ["file", "bytes", "sha256", "url"])


def forest_svg(effects: list[dict[str, object]], summary: list[dict[str, object]]) -> None:
    estimated = [x for x in effects if x["status"] == "estimated"]
    pooled = [x for x in summary if x["model"] == "random_DL"]
    rows: list[tuple[str, str, float, float, float, bool]] = []
    for gene in GENES:
        for item in estimated:
            if item["gene"] == gene:
                rows.append((gene, str(item["cohort"]), float(item["hedges_g"]), float(item["ci_low"]), float(item["ci_high"]), False))
        item = next(x for x in pooled if x["gene"] == gene)
        rows.append((gene, "Random-effects", float(item["estimate"]), float(item["ci_low"]), float(item["ci_high"]), True))
    width, height = 900, 90 + 42 * len(rows)
    xmin = min(-2.5, min(x[3] for x in rows) - 0.2)
    xmax = max(2.5, max(x[4] for x in rows) + 0.2)
    left, right = 250, 820
    scale = lambda value: left + (value - xmin) / (xmax - xmin) * (right - left)
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<style>text{font-family:Arial,sans-serif;font-size:14px}.gene{font-weight:bold}.pool{font-weight:bold}</style>',
        '<rect width="100%" height="100%" fill="white"/>',
        '<text x="20" y="28" font-size="18" font-weight="bold">ICI response: expression in responders vs nonresponders</text>',
        f'<line x1="{scale(0):.1f}" y1="48" x2="{scale(0):.1f}" y2="{height-35}" stroke="#888" stroke-dasharray="4 4"/>',
    ]
    last_gene = None
    for i, (gene, cohort, est, low, high, is_pool) in enumerate(rows):
        y = 70 + i * 42
        if gene != last_gene:
            svg.append(f'<text class="gene" x="20" y="{y}">{gene}</text>')
            last_gene = gene
        css = ' class="pool"' if is_pool else ""
        svg.append(f'<text{css} x="100" y="{y}">{cohort}</text>')
        svg.append(f'<line x1="{scale(low):.1f}" y1="{y-5}" x2="{scale(high):.1f}" y2="{y-5}" stroke="#2457a6" stroke-width="2"/>')
        if is_pool:
            x = scale(est)
            svg.append(f'<polygon points="{x-7:.1f},{y-5:.1f} {x:.1f},{y-12:.1f} {x+7:.1f},{y-5:.1f} {x:.1f},{y+2:.1f}" fill="#c43d3d"/>')
        else:
            svg.append(f'<circle cx="{scale(est):.1f}" cy="{y-5}" r="5" fill="#2457a6"/>')
        svg.append(f'<text{css} x="835" y="{y}" text-anchor="end">{est:.2f}</text>')
    for tick in range(math.ceil(xmin), math.floor(xmax) + 1):
        x = scale(tick)
        svg.append(f'<line x1="{x:.1f}" y1="{height-30}" x2="{x:.1f}" y2="{height-25}" stroke="black"/>')
        svg.append(f'<text x="{x:.1f}" y="{height-8}" text-anchor="middle">{tick}</text>')
    svg.append(f'<text x="{(left+right)/2:.1f}" y="{height-8}" text-anchor="middle" transform="translate(0,-20)">Hedges g (positive = higher in responders)</text>')
    svg.append("</svg>")
    (OUT / "forest_plot.svg").write_text("\n".join(svg), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    download()
    long_rows, availability = sample_rows()
    effects = cohort_effects(long_rows, availability)
    summary = meta_summary(effects)
    write_tsv(
        OUT / "sample_expression.tsv",
        long_rows,
        ["cohort", "sample", "response", "gene", "expression", "analysis_scale"],
    )
    write_tsv(
        OUT / "cohort_effects.tsv",
        effects,
        ["cohort", "gene", "status", "n_responder", "n_nonresponder", "mean_responder",
         "mean_nonresponder", "hedges_g", "se", "ci_low", "ci_high", "p_value"],
    )
    write_tsv(
        OUT / "meta_summary.tsv",
        summary,
        ["gene", "model", "k", "estimate", "se", "ci_low", "ci_high", "p_value",
         "Q", "Q_df", "I2_percent", "tau2_DL"],
    )
    write_manifest()
    forest_svg(effects, summary)
    print(f"Wrote {len(long_rows)} sample-gene rows, {len(effects)} cohort rows, and {len(summary)} pooled rows.")


if __name__ == "__main__":
    main()
