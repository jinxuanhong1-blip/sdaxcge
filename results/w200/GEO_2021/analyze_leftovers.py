#!/usr/bin/env python3
"""Analyze the two outcome-bearing cohorts among 2021 GEO leftovers."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test


BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE190nnn"
FILES = {
    "GSE190265_TPM_France3.csv.gz":
        f"{BASE}/GSE190265/suppl/GSE190265_TPM_France3.csv.gz",
    "GSE190265_samples_info_France3.csv.gz":
        f"{BASE}/GSE190265/suppl/GSE190265_samples_info_France3.csv.gz",
    "GSE190266_TPM_France4.csv.gz":
        f"{BASE}/GSE190266/suppl/GSE190266_TPM_France4.csv.gz",
    "GSE190266_series_matrix.txt.gz":
        f"{BASE}/GSE190266/matrix/GSE190266_series_matrix.txt.gz",
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("/tmp/GEO_2021"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("survival_results.csv"),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path(__file__).with_name("input_manifest.csv"),
    )
    parser.add_argument("--no-download", action="store_true")
    return parser.parse_args()


def obtain_files(data_dir: Path, allow_download: bool) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    for filename, url in FILES.items():
        destination = data_dir / filename
        if destination.is_file() and destination.stat().st_size:
            continue
        if not allow_download:
            raise FileNotFoundError(destination)
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, destination)


def parse_france4_metadata(path: Path) -> pd.DataFrame:
    sample_names: list[str] | None = None
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt") as source:
        for line in source:
            fields = [value.strip('"') for value in line.rstrip("\n").split("\t")]
            if fields[0] == "!Sample_title":
                sample_names = fields[1:]
            elif fields[0] == "!Sample_characteristics_ch1":
                characteristics.append(fields[1:])
    if sample_names is None:
        raise ValueError("GSE190266 sample titles not found")

    def characteristic(prefix: str, cast: type) -> list[object]:
        row = next(values for values in characteristics if values[0].startswith(prefix))
        return [cast(value.split(": ", 1)[1]) for value in row]

    return pd.DataFrame(
        {
            "time": characteristic("pfs_time", float),
            "event": characteristic("pfs_evt", int),
        },
        index=sample_names,
    )


def load_cohorts(data_dir: Path) -> list[tuple[str, pd.DataFrame, pd.DataFrame, list[str]]]:
    france3_expression = pd.read_csv(
        data_dir / "GSE190265_TPM_France3.csv.gz", sep=";", index_col=0
    )
    france3_metadata = pd.read_csv(
        data_dir / "GSE190265_samples_info_France3.csv.gz", sep=";"
    ).set_index("sample")
    france3_metadata = france3_metadata.rename(
        columns={"time_PFS": "time", "evtPFS": "event"}
    )

    france4_expression = pd.read_csv(
        data_dir / "GSE190266_TPM_France4.csv.gz",
        sep=";",
        decimal=",",
        index_col=0,
    )
    france4_metadata = parse_france4_metadata(
        data_dir / "GSE190266_series_matrix.txt.gz"
    )
    return [
        (
            "GSE190265",
            france3_expression,
            france3_metadata,
            ["TACSTD2", "CLDN4"],
        ),
        (
            "GSE190266",
            france4_expression,
            france4_metadata,
            ["CLDN4"],  # TACSTD2 is absent from the deposited processed matrix.
        ),
    ]


def analyze(
    dataset: str,
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    genes: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for gene in genes:
        frame = metadata[["time", "event"]].join(expression[[gene]], how="inner")
        frame = frame.dropna()
        frame["log2_tpm_plus_1"] = np.log2(frame[gene] + 1)
        model = CoxPHFitter().fit(
            frame[["time", "event", "log2_tpm_plus_1"]],
            duration_col="time",
            event_col="event",
        )
        coefficient = model.summary.loc["log2_tpm_plus_1"]
        median = frame["log2_tpm_plus_1"].median()
        high = frame["log2_tpm_plus_1"] > median
        logrank = logrank_test(
            frame.loc[high, "time"],
            frame.loc[~high, "time"],
            frame.loc[high, "event"],
            frame.loc[~high, "event"],
        )
        rows.append(
            {
                "dataset": dataset,
                "gene": gene,
                "n": len(frame),
                "events": int(frame["event"].sum()),
                "cox_hr_per_log2_tpm_plus_1": coefficient["exp(coef)"],
                "cox_ci95_low": coefficient["exp(coef) lower 95%"],
                "cox_ci95_high": coefficient["exp(coef) upper 95%"],
                "cox_p": coefficient["p"],
                "median_split_logrank_p": logrank.p_value,
                "high_n": int(high.sum()),
                "low_n": int((~high).sum()),
            }
        )
    return rows


def add_bh_q_values(rows: list[dict[str, object]]) -> None:
    p_values = np.array([float(row["cox_p"]) for row in rows])
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1)
    q_values = np.empty_like(adjusted)
    q_values[order] = adjusted
    for row, q_value in zip(rows, q_values):
        row["cox_bh_q_across_3_tests"] = q_value


def write_manifest(data_dir: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", newline="") as output:
        writer = csv.writer(output)
        writer.writerow(["file", "bytes", "sha256", "url"])
        for filename, url in FILES.items():
            path = data_dir / filename
            writer.writerow(
                [filename, path.stat().st_size, hashlib.sha256(path.read_bytes()).hexdigest(), url]
            )


def main() -> None:
    args = arguments()
    obtain_files(args.data_dir, allow_download=not args.no_download)
    rows: list[dict[str, object]] = []
    for cohort in load_cohorts(args.data_dir):
        rows.extend(analyze(*cohort))
    add_bh_q_values(rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output, index=False)
    write_manifest(args.data_dir, args.manifest)
    print(f"Wrote {len(rows)} tests to {args.output}")


if __name__ == "__main__":
    main()
