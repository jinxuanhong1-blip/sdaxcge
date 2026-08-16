#!/usr/bin/env python3
"""IFN-treated lung-line GEO add-on for the DepMap TACSTD2 analysis.

DepMap Public 24Q4 has no IFN-treated transcriptomes. This script tests
whether recombinant IFN changes TACSTD2 in public processed GEO lung-line
datasets, with ISG genes as treatment-quality controls.
"""

from __future__ import annotations

import argparse
import gzip
import io
import json
import re
import tarfile
import urllib.request
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


CONTROL_GENES = ["TACSTD2", "MX1", "ISG15", "STAT1", "CXCL10"]
ENSG = {
    "TACSTD2": "ENSG00000184292",
    "MX1": "ENSG00000157601",
    "ISG15": "ENSG00000187608",
    "STAT1": "ENSG00000115415",
    "CXCL10": "ENSG00000169245",
}
FILES = {
    "GSE5542_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5542/matrix/GSE5542_series_matrix.txt.gz",
    "GSE60572_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE60nnn/GSE60572/matrix/GSE60572_series_matrix.txt.gz",
    "GPL96.annot.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz",
    "GPL10558.annot.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/annot/GPL10558.annot.gz",
    "GSE156295_count_hg38.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156295/suppl/GSE156295_count_hg38.txt.gz",
    "GSE156295_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE156nnn/GSE156295/matrix/GSE156295_series_matrix.txt.gz",
    "GSE178640_RAW.tar":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178640/suppl/GSE178640_RAW.tar",
    "GSE178640_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178640/matrix/GSE178640_series_matrix.txt.gz",
    "GSE215771_Expression_Values.csv.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/GSE215771/suppl/GSE215771_Expression_Values.csv.gz",
    "GSE215771_filelist.txt":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/GSE215771/suppl/filelist.txt",
    "GSE215771_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE215nnn/GSE215771/matrix/GSE215771_series_matrix.txt.gz",
    "GSE109720_COUNTS_genes_LUNGCANCER_10.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109720/suppl/GSE109720_COUNTS_genes_LUNGCANCER_10.txt.gz",
    "GSE109720_series_matrix.txt.gz":
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109720/matrix/GSE109720_series_matrix.txt.gz",
}


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as response, dest.open("wb") as handle:
        while True:
            chunk = response.read(1 << 20)
            if not chunk:
                break
            handle.write(chunk)


def parse_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    meta_rows: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if in_table:
                table_lines.append(line)
                continue
            if line.startswith("!Sample_"):
                key, *values = line.rstrip("\n").split("\t")
                values = [v.strip('"') for v in values]
                if key in meta_rows and len(meta_rows[key]) == len(values):
                    meta_rows[key] = [
                        f"{old}; {new}" for old, new in zip(meta_rows[key], values)
                    ]
                else:
                    meta_rows[key] = values
    meta = pd.DataFrame(meta_rows)
    if "!Sample_geo_accession" in meta:
        meta = meta.set_index("!Sample_geo_accession")
    expr = pd.DataFrame()
    if table_lines:
        expr = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
        expr.index = expr.index.astype(str)
    return meta, expr


def annot_map(path: Path, genes: list[str]) -> dict[str, list[str]]:
    mapping = {g: [] for g in genes}
    with gzip.open(path, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("ID"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            probe, symbol = fields[0], fields[2]
            if symbol in mapping:
                mapping[symbol].append(probe)
    return mapping


def log2p1(values: np.ndarray) -> np.ndarray:
    return np.log2(np.asarray(values, dtype=float) + 1.0)


def contrast_stats(treated: np.ndarray, control: np.ndarray) -> dict[str, float]:
    treated = np.asarray(treated, dtype=float)
    control = np.asarray(control, dtype=float)
    treated = treated[np.isfinite(treated)]
    control = control[np.isfinite(control)]
    result = {
        "n_treated": int(len(treated)),
        "n_control": int(len(control)),
        "mean_treated": float(np.mean(treated)) if len(treated) else np.nan,
        "mean_control": float(np.mean(control)) if len(control) else np.nan,
        "log2fc": np.nan,
        "p": np.nan,
    }
    if len(treated) and len(control):
        result["log2fc"] = result["mean_treated"] - result["mean_control"]
    if len(treated) >= 2 and len(control) >= 2 and (
        np.unique(treated).size > 1 or np.unique(control).size > 1
    ):
        result["p"] = float(stats.ttest_ind(treated, control, equal_var=False).pvalue)
    return result


def collapse_probes(expr: pd.DataFrame, probes: list[str]) -> pd.Series:
    present = [p for p in probes if p in expr.index]
    if not present:
        return pd.Series(np.nan, index=expr.columns)
    return expr.loc[present].apply(pd.to_numeric, errors="coerce").median(axis=0)


def cpm_log(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    return np.log2((counts.div(lib, axis=1) * 1e6) + 1.0)


def collect_rows(rows: list[dict], **kwargs) -> None:
    rows.append(kwargs)


def gse5542(geo: Path, rows: list[dict], values: list[dict]) -> None:
    meta, expr = parse_series_matrix(geo / "GSE5542_series_matrix.txt.gz")
    probes = annot_map(geo / "GPL96.annot.gz", CONTROL_GENES)
    log_expr = np.log2(expr.apply(pd.to_numeric, errors="coerce") + 1.0)
    titles = meta["!Sample_title"]
    parsed = titles.str.extract(
        r"(?P<treatment>Infergen\+IFN_Gamma|IFN_Gamma|Infergen|untreated)(?:_(?P<time>6hr|24hr))?",
        expand=True,
    )
    parsed["time"] = parsed["time"].fillna(
        titles.str.extract(r"(6hr|24hr)", expand=False)
    )
    parsed.loc[parsed["treatment"].eq("untreated") & parsed["time"].isna(), "time"] = "24hr"
    gene_series = {g: collapse_probes(log_expr, probes[g]) for g in CONTROL_GENES}
    for time in ["6hr", "24hr"]:
        control = parsed["treatment"].eq("untreated") & parsed["time"].eq(time)
        for treatment, ifn_type in [("Infergen", "IFNA"), ("IFN_Gamma", "IFNG")]:
            treated = parsed["treatment"].eq(treatment) & parsed["time"].eq(time)
            add_gene_contrasts(
                rows, values, "GSE5542", "A549", ifn_type, time, "array_log2",
                treated.index[treated], treated.index[control], gene_series, primary=True,
            )
        combo = parsed["treatment"].eq("Infergen+IFN_Gamma") & parsed["time"].eq(time)
        add_gene_contrasts(
            rows, values, "GSE5542", "A549", "IFNA+IFNG", time, "array_log2",
            list(combo.index[combo]), list(parsed.index[control]),
            gene_series, primary=False,
        )


def gse60572(geo: Path, rows: list[dict], values: list[dict]) -> None:
    meta, expr = parse_series_matrix(geo / "GSE60572_series_matrix.txt.gz")
    probes = annot_map(geo / "GPL10558.annot.gz", CONTROL_GENES)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    chars = meta["!Sample_characteristics_ch1"]
    treated_with = chars.str.extract(r"treated with: ([^;]+)", expand=False).str.strip()
    time = chars.str.extract(r"time point: ([^;]+)", expand=False).str.replace(" hrs", "h")
    gene_series = {g: collapse_probes(expr, probes[g]) for g in CONTROL_GENES}
    for hours in ["6h", "24h"]:
        control = (
            treated_with.str.contains("none|non_treated|non-treated|untreated", case=False, na=False)
            & ~treated_with.str.contains("IFN|ppp5", case=False, na=False)
            & time.eq(hours)
        )
        if not control.any():
            control = treated_with.str.contains("RNAiMax", case=False, na=False) & time.eq(hours)
        for dose in ["IFN_2b_100", "IFN_2b_1000"]:
            treated = treated_with.eq(dose) & time.eq(hours)
            add_gene_contrasts(
                rows, values, "GSE60572", "A549", "IFNB", f"{hours}_{dose}",
                "array_log2", list(treated.index[treated]), list(control.index[control]),
                gene_series, primary=bool(control.any()),
                note="" if control.any() else "no matched 24h untreated samples deposited",
            )


def gse156295(geo: Path, rows: list[dict], values: list[dict]) -> None:
    counts = pd.read_csv(geo / "GSE156295_count_hg38.txt.gz", sep="\t").set_index("Gene")
    log = cpm_log(counts)
    gene_series = {g: log.loc[g] if g in log.index else pd.Series(np.nan, index=log.columns)
                   for g in CONTROL_GENES}
    for line in ["A549", "HTBE"]:
        treated = [c for c in log.columns if c.startswith(f"{line}-IFN")]
        control = [c for c in log.columns if c.startswith(f"{line}-Mock")]
        add_gene_contrasts(
            rows, values, "GSE156295", line, "IFNI", "unspecified", "log2CPM",
            treated, control, gene_series, primary=line == "A549",
        )


def gse178640(geo: Path, rows: list[dict], values: list[dict]) -> None:
    counts = {}
    with tarfile.open(geo / "GSE178640_RAW.tar") as archive:
        for member in archive.getmembers():
            handle = archive.extractfile(member)
            if handle is None:
                continue
            raw = gzip.decompress(handle.read()).decode()
            frame = pd.read_csv(io.StringIO(raw), sep="\t")
            sample = Path(member.name).name.split("_", 1)[0]
            series = frame.set_index("Identifier")["matchCounts"]
            series.index = series.index.str.replace(r"\.\d+$", "", regex=True)
            counts[sample] = series
    counts = pd.DataFrame(counts).fillna(0)
    log = cpm_log(counts)
    meta, _ = parse_series_matrix(geo / "GSE178640_series_matrix.txt.gz")
    titles = meta["!Sample_title"]
    gene_series = {}
    for gene, ensembl in ENSG.items():
        gene_series[gene] = log.loc[ensembl] if ensembl in log.index else pd.Series(np.nan, index=log.columns)
    treated = titles.index[titles.str.startswith("DMSO_IFN")]
    control = titles.index[titles.str.startswith("DMSO_no")]
    add_gene_contrasts(
        rows, values, "GSE178640", "A549", "IFNA2", "DMSO_pretreated", "log2CPM",
        list(treated), list(control), gene_series, primary=True,
    )


def gse215771(geo: Path, rows: list[dict], values: list[dict]) -> None:
    filelist = (geo / "GSE215771_filelist.txt").read_text()
    barcode_to_gsm = {}
    for line in filelist.splitlines():
        match = re.search(r"(GSM\d+)_([a-z])\.([A-Z-]+)\.", line)
        if match:
            barcode_to_gsm[match.group(3)] = match.group(1)
    meta, _ = parse_series_matrix(geo / "GSE215771_series_matrix.txt.gz")
    expr = pd.read_csv(geo / "GSE215771_Expression_Values.csv.gz", low_memory=False)
    tpm_cols = [c for c in expr.columns if c.endswith("TPM")]
    mapped = {}
    for col in tpm_cols:
        match = re.search(r"([A-Z]{10}-[A-Z]{10})", col)
        if match and match.group(1) in barcode_to_gsm:
            mapped[barcode_to_gsm[match.group(1)]] = col
    gene_series = {}
    for gene in CONTROL_GENES:
        row = expr.loc[expr["Name"].eq(gene)]
        series = pd.Series(np.nan, index=mapped)
        if not row.empty:
            series = pd.Series({gsm: row.iloc[0][col] for gsm, col in mapped.items()}, dtype=float)
        gene_series[gene] = np.log2(series + 1.0)
    genotype = meta["!Sample_characteristics_ch1"].str.extract(r"genotype: ([^;]+)", expand=False)
    treatment = meta["!Sample_characteristics_ch1"].str.extract(r"treatment: ([^;]+)", expand=False)
    for geno, primary in [("WT", True), ("IRF1KO", False), ("NF2KO", False)]:
        treated = genotype.eq(geno) & treatment.str.contains("IFN", na=False)
        control = genotype.eq(geno) & treatment.eq("Control")
        add_gene_contrasts(
            rows, values, "GSE215771", f"A549_{geno}", "IFNG", "24h", "log2TPM",
            list(treated.index[treated]), list(control.index[control]),
            gene_series, primary=primary,
        )


def gse109720(geo: Path, rows: list[dict], values: list[dict]) -> None:
    with gzip.open(geo / "GSE109720_COUNTS_genes_LUNGCANCER_10.txt.gz", "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        samples = header[1:]
        records = []
        for line in handle:
            left, *vals = line.rstrip("\n").split("\t")
            ensembl, gene, _ = left.split(",", 2)
            if gene in CONTROL_GENES:
                records.append((gene, [float(v) for v in vals]))
    counts = pd.DataFrame({gene: vals for gene, vals in records}, index=samples).T
    # Library sizes are unavailable from the gene subset; use log2(count+1).
    log = np.log2(counts + 1.0)
    meta, _ = parse_series_matrix(geo / "GSE109720_series_matrix.txt.gz")
    titles = list(meta["!Sample_title"])
    if len(titles) != len(samples):
        raise ValueError("GSE109720 sample-title count does not match count columns")
    rename = dict(zip(samples, meta.index))
    log = log.rename(columns=rename)
    gene_series = {g: log.loc[g] if g in log.index else pd.Series(np.nan, index=log.columns)
                   for g in CONTROL_GENES}
    title_by_gsm = pd.Series(titles, index=meta.index)
    for line, control_pat, treated_pat in [
        ("EBC-1", r"^EBC1R\d+$", r"^EBC1IFN\d+$"),
        ("NCI-H1573", r"^H1573B\d+$", r"^H1573IFN\d+$"),
        ("NCI-H1993", r"^H1993B\d+$", r"^H1993IFN\d+$"),
        ("NCI-H596", r"^H596B\d+$", r"^H596IFN\d+$"),
    ]:
        control = title_by_gsm.index[title_by_gsm.str.match(control_pat)]
        treated = title_by_gsm.index[title_by_gsm.str.match(treated_pat)]
        add_gene_contrasts(
            rows, values, "GSE109720", line, "IFNG", "unspecified", "log2count",
            list(treated), list(control), gene_series, primary=True,
            note="column-to-sample order assumed from series matrix; ISG QC required",
        )


def add_gene_contrasts(
    rows: list[dict],
    values: list[dict],
    dataset: str,
    line: str,
    ifn_type: str,
    timepoint: str,
    scale: str,
    treated_ids: list[str],
    control_ids: list[str],
    gene_series: dict[str, pd.Series],
    primary: bool,
    note: str = "",
) -> None:
    isg_log2fc = {}
    for gene in CONTROL_GENES:
        series = gene_series[gene]
        treated = series.reindex(treated_ids).to_numpy()
        control = series.reindex(control_ids).to_numpy()
        stats_row = contrast_stats(treated, control)
        if gene != "TACSTD2":
            isg_log2fc[gene] = stats_row["log2fc"]
        collect_rows(
            rows,
            dataset=dataset,
            cell_line=line,
            ifn_type=ifn_type,
            timepoint=timepoint,
            scale=scale,
            gene=gene,
            primary=primary,
            note=note,
            **stats_row,
        )
        for sample_id, value, arm in (
            [(s, series.get(s, np.nan), "treated") for s in treated_ids]
            + [(s, series.get(s, np.nan), "control") for s in control_ids]
        ):
            values.append({
                "dataset": dataset, "cell_line": line, "ifn_type": ifn_type,
                "timepoint": timepoint, "gene": gene, "sample_id": sample_id,
                "arm": arm, "value": value, "scale": scale,
            })
    mx1 = isg_log2fc.get("MX1", np.nan)
    cxcl10 = isg_log2fc.get("CXCL10", np.nan)
    isg15 = isg_log2fc.get("ISG15", np.nan)
    passed = bool(
        pd.notna(mx1) and mx1 >= 1
        or (pd.notna(isg15) and isg15 >= 1 and (pd.isna(mx1) or mx1 > -0.5))
        or (pd.notna(cxcl10) and cxcl10 >= 1 and (pd.isna(mx1) or mx1 > -0.5))
    )
    if pd.notna(mx1) and mx1 <= -1:
        passed = False
    for row in rows[::-1]:
        if row["dataset"] == dataset and row["cell_line"] == line and row["ifn_type"] == ifn_type and row["timepoint"] == timepoint:
            row["isg_qc_pass"] = passed
            if row["gene"] == "TACSTD2":
                break


def forest_plot(contrasts: pd.DataFrame, output: Path) -> None:
    data = contrasts[
        contrasts["gene"].eq("TACSTD2") & contrasts["primary"] & contrasts["isg_qc_pass"]
    ].copy()
    if data.empty:
        data = contrasts[contrasts["gene"].eq("TACSTD2") & contrasts["primary"]].copy()
    data["label"] = data.apply(
        lambda r: f"{r.dataset} {r.cell_line} {r.ifn_type} {r.timepoint}", axis=1
    )
    data = data.sort_values("log2fc")
    sns.set_theme(style="whitegrid", context="talk")
    fig, axis = plt.subplots(figsize=(11, max(3.5, 0.45 * len(data))), constrained_layout=True)
    colors = np.where(data["log2fc"] >= 0, "#c23b23", "#2c6e9e")
    axis.barh(data["label"], data["log2fc"], color=colors, alpha=0.85)
    axis.axvline(0, color="black", linewidth=1)
    axis.set_xlabel("TACSTD2 log2 fold-change (IFN minus control)")
    axis.set_title("IFN-treated lung lines: TACSTD2")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def qc_plot(contrasts: pd.DataFrame, output: Path) -> None:
    data = contrasts[contrasts["primary"] & contrasts["gene"].isin(["TACSTD2", "MX1", "ISG15", "CXCL10"])].copy()
    data["label"] = data.apply(
        lambda r: f"{r.dataset} {r.cell_line} {r.ifn_type} {r.timepoint}", axis=1
    )
    sns.set_theme(style="whitegrid", context="talk")
    fig, axis = plt.subplots(figsize=(12, 6.5), constrained_layout=True)
    sns.barplot(data=data, x="label", y="log2fc", hue="gene", ax=axis)
    axis.axhline(0, color="black", linewidth=1)
    axis.tick_params(axis="x", rotation=75, labelsize=8)
    axis.set_ylabel("log2 fold-change")
    axis.set_xlabel("")
    axis.set_title("IFN treatment QC and TACSTD2")
    fig.savefig(output, dpi=200)
    plt.close(fig)


def update_report(report: Path, contrasts: pd.DataFrame, inventory: pd.DataFrame) -> None:
    tac = contrasts[contrasts["gene"].eq("TACSTD2")].copy()
    primary = tac[tac["primary"]]
    passed = primary[primary["isg_qc_pass"]]
    failed = primary[~primary["isg_qc_pass"]]
    n_sig = int(((passed["q"] < 0.05) & passed["log2fc"].notna()).sum()) if "q" in passed else 0
    induced = passed[passed["q"] < 0.05] if "q" in passed else passed.iloc[0:0]
    direction = "no BH-significant TACSTD2 induction or repression"
    if n_sig:
        pos = int((induced["log2fc"] > 0).sum())
        neg = int((induced["log2fc"] < 0).sum())
        direction = f"{pos} induced and {neg} repressed TACSTD2 contrasts at q<0.05 among ISG-passing tests"
    def md(frame: pd.DataFrame) -> str:
        shown = frame[[
            "dataset", "cell_line", "ifn_type", "timepoint", "n_treated", "n_control",
            "log2fc", "p", "q", "isg_qc_pass",
        ]].copy()
        shown["log2fc"] = shown["log2fc"].map(lambda x: "NA" if pd.isna(x) else f"{x:.3f}")
        shown["p"] = shown["p"].map(lambda x: "NA" if pd.isna(x) else (f"{x:.2e}" if x < 0.001 else f"{x:.3f}"))
        shown["q"] = shown["q"].map(lambda x: "NA" if pd.isna(x) else (f"{x:.2e}" if x < 0.001 else f"{x:.3f}"))
        return shown.to_markdown(index=False)
    section = [
        "## IFN-treated lung lines (GEO; not DepMap)",
        "",
        "DepMap Public 24Q4 has **no IFN-treated RNA-seq or microarray profiles**. "
        "PRISM also lacks recombinant IFN biologics. The treatment test therefore uses "
        "public GEO processed matrices of lung lines treated with IFN protein.",
        "",
        f"- Primary TACSTD2 contrasts with ISG QC pass: {len(passed)}. {direction}.",
        f"- Primary contrasts that failed ISG QC (MX1/ISG15/CXCL10 log2FC all < 1): {len(failed)}.",
        "- Most usable public IFN time courses are A549. In DepMap 24Q4 A549 TACSTD2 is "
        "very low (log2(TPM+1)=0.20), so these experiments test induction from a low basal state, "
        "not regulation in TROP2-high NSCLC.",
        "- A549-heavy GEO data cannot confirm or refute the basal DepMap TACSTD2–IFN "
        "correlation. They only ask whether IFN itself moves TACSTD2.",
        "",
        "Primary TACSTD2 contrasts (BH across TACSTD2 primary tests):",
        "",
        md(primary.sort_values(["isg_qc_pass", "q", "p"], ascending=[False, True, True])),
        "",
        "ISG QC requires MX1 log2FC >= 1, or ISG15/CXCL10 log2FC >= 1 without MX1 collapse. "
        "A contrast with MX1 log2FC <= -1 is treated as a failed IFN response. "
        "GSE109720 sample columns were aligned to series-matrix order; failed-QC lines "
        "are not interpreted as IFN-responsive. BH q-values are only among ISG-passing "
        "primary TACSTD2 tests.",
        "",
        "Inventory:",
        "",
        inventory.to_markdown(index=False),
        "",
        "",
    ]
    text = report.read_text() if report.exists() else ""
    start = text.find("## IFN-treated lung lines")
    end = text.find("## Robustness and limitations")
    block = "\n".join(section)
    if start != -1 and end != -1:
        text = text[:start] + block + text[end:]
    elif end != -1:
        text = text[:end] + block + text[end:]
    else:
        text = text.rstrip() + "\n\n" + block
    # bottom-line bullet
    bullet = (
        f"- DepMap has no IFN-treated profiles. In public GEO lung-line IFN experiments, "
        f"{direction}."
    )
    marker = "- The analysis is observational and cross-sectional."
    if bullet not in text and marker in text:
        text = text.replace(marker, bullet + "\n" + marker)
    report.write_text(text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo-dir", type=Path, default=Path(__file__).parent / "geo")
    parser.add_argument("--out-dir", type=Path, default=Path(__file__).parent)
    args = parser.parse_args()
    args.geo_dir.mkdir(parents=True, exist_ok=True)
    tables = args.out_dir / "tables"
    figures = args.out_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    for name, url in FILES.items():
        download(url, args.geo_dir / name)

    rows: list[dict] = []
    values: list[dict] = []
    gse5542(args.geo_dir, rows, values)
    gse60572(args.geo_dir, rows, values)
    gse156295(args.geo_dir, rows, values)
    gse178640(args.geo_dir, rows, values)
    gse215771(args.geo_dir, rows, values)
    gse109720(args.geo_dir, rows, values)

    contrasts = pd.DataFrame(rows)
    tac_primary = contrasts["gene"].eq("TACSTD2") & contrasts["primary"]
    contrasts["q"] = np.nan
    valid = tac_primary & contrasts["isg_qc_pass"] & contrasts["p"].notna()
    if valid.any():
        contrasts.loc[valid, "q"] = multipletests(contrasts.loc[valid, "p"], method="fdr_bh")[1]
    contrasts = contrasts.sort_values(["gene", "primary", "dataset", "cell_line", "ifn_type", "timepoint"])
    values_frame = pd.DataFrame(values)

    inventory = pd.DataFrame([
        {
            "dataset": "DepMap 24Q4", "included": False,
            "reason": "No IFN-treated transcriptome; basal RNA/CRISPR only",
        },
        {
            "dataset": "GSE5542", "included": True,
            "reason": "A549 IFNA-con1 / IFNG 6h and 24h, n=4, processed series matrix",
        },
        {
            "dataset": "GSE60572", "included": True,
            "reason": "A549 IFNB 100/1000 IU; only 6h has matched untreated samples",
        },
        {
            "dataset": "GSE156295", "included": True,
            "reason": "A549 and HTBE type I IFN vs mock, n=2, processed counts",
        },
        {
            "dataset": "GSE178640", "included": True,
            "reason": "A549 IFN-alpha2 after DMSO, n=3, processed Ensembl counts",
        },
        {
            "dataset": "GSE215771", "included": True,
            "reason": "A549 IFNG 24h WT/IRF1KO/NF2KO, n=3, processed TPM",
        },
        {
            "dataset": "GSE109720", "included": True,
            "reason": "EBC-1/H1573/H1993/H596 IFNG; ISG QC applied because some lines do not induce ISGs",
        },
        {
            "dataset": "GSE241977", "included": False,
            "reason": "A549 AhR KO vs control; no IFN-versus-mock contrast in the deposited samples",
        },
        {
            "dataset": "GSE261920", "included": False,
            "reason": "lncRNA-focused A549 IFN series; not used as a coding-gene TACSTD2 test",
        },
    ])

    contrasts.to_csv(tables / "ifn_treated_contrasts.csv", index=False)
    values_frame.to_csv(tables / "ifn_treated_sample_values.csv", index=False)
    inventory.to_csv(tables / "ifn_treated_inventory.csv", index=False)
    forest_plot(contrasts, figures / "ifn_treated_tacstd2_forest.png")
    qc_plot(contrasts, figures / "ifn_treated_isg_qc.png")
    update_report(args.out_dir / "REPORT.md", contrasts, inventory)

    summary = {
        "depmap_ifn_treated_profiles": False,
        "primary_tacstd2_contrasts": int(tac_primary.sum()),
        "isg_qc_pass_primary": int((tac_primary & contrasts["isg_qc_pass"]).sum()),
        "bh_significant_isg_pass": int((
            tac_primary & contrasts["isg_qc_pass"] & contrasts["q"].lt(0.05)
        ).sum()),
    }
    (args.out_dir / "ifn_treated_qc.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
