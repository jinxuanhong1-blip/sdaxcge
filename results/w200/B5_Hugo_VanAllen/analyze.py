#!/usr/bin/env python3
"""Hugo 2016 and Van Allen 2015: TACSTD2/CLDN4 versus melanoma ICI response."""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import urllib.request
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
SEED = 2015
N_BOOT = 20_000

URLS = {
    DATA / "GSE78220_PatientFPKM.xlsx": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/suppl/"
        "GSE78220_PatientFPKM.xlsx"
    ),
    DATA / "GSE78220_series_matrix.txt.gz": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE78nnn/GSE78220/matrix/"
        "GSE78220_series_matrix.txt.gz"
    ),
    DATA / "TPM_RSEM_VAScience2015.txt": (
        "https://raw.githubusercontent.com/vanallenlab/"
        "VanAllen_CTLA4_Science_RNASeq_TPM/master/TPM_RSEM_VAScience2015.txt"
    ),
    DATA / "skcm_dfci_2015_clinical_patient.txt": (
        "https://media.githubusercontent.com/media/cBioPortal/datahub/master/"
        "public/skcm_dfci_2015/data_clinical_patient.txt"
    ),
}

HUGO_RESPONSE = {
    "Complete Response": "CR",
    "Partial Response": "PR",
    "Progressive Disease": "PD",
}
VA_GENES = {
    "TACSTD2": "ENSG00000184292",
    "CLDN4": "ENSG00000189143",
}


def download_inputs() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for path, url in URLS.items():
        if path.exists() and path.stat().st_size > 1000:
            continue
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, path)
        if path.stat().st_size < 1000:
            raise RuntimeError(f"Download too small: {path} from {url}")


def bootstrap_auc_ci(
    responders: np.ndarray, nonresponders: np.ndarray, draws: int = N_BOOT
) -> tuple[float, float]:
    rng = np.random.default_rng(SEED)
    estimates = np.empty(draws)
    n_r, n_nr = len(responders), len(nonresponders)
    for i in range(draws):
        a = rng.choice(responders, n_r, replace=True)
        b = rng.choice(nonresponders, n_nr, replace=True)
        estimates[i] = stats.mannwhitneyu(a, b).statistic / (n_r * n_nr)
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def summarize(responders: np.ndarray, nonresponders: np.ndarray) -> dict:
    test = stats.mannwhitneyu(
        responders, nonresponders, alternative="two-sided", method="asymptotic"
    )
    n_r, n_nr = len(responders), len(nonresponders)
    auc = float(test.statistic / (n_r * n_nr))
    auc_lo, auc_hi = bootstrap_auc_ci(responders, nonresponders)
    detected_r = int(np.sum(responders > 0))
    detected_nr = int(np.sum(nonresponders > 0))
    fisher = stats.fisher_exact(
        [
            [detected_r, n_r - detected_r],
            [detected_nr, n_nr - detected_nr],
        ],
        alternative="two-sided",
    )
    return {
        "n_responder": n_r,
        "n_nonresponder": n_nr,
        "median_responder": float(np.median(responders)),
        "median_nonresponder": float(np.median(nonresponders)),
        "mann_whitney_u": float(test.statistic),
        "mann_whitney_p": float(test.pvalue),
        "auc_responder_higher": auc,
        "auc_bootstrap_95ci_low": float(auc_lo),
        "auc_bootstrap_95ci_high": float(auc_hi),
        "detected_responder": detected_r,
        "detected_nonresponder": detected_nr,
        "detection_fisher_odds_ratio": float(fisher.statistic),
        "detection_fisher_p": float(fisher.pvalue),
    }


def parse_hugo_metadata() -> pd.DataFrame:
    matrix = DATA / "GSE78220_series_matrix.txt.gz"
    sample_rows: list[list[str]] = []
    with gzip.open(matrix, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_"):
                sample_rows.append(next(csv.reader([line], delimiter="\t")))

    titles = next(row[1:] for row in sample_rows if row[0] == "!Sample_title")
    accessions = next(
        row[1:] for row in sample_rows if row[0] == "!Sample_geo_accession"
    )
    if len(titles) != 28 or len(accessions) != 28:
        raise ValueError(f"Expected 28 Hugo samples; found {len(titles)}")

    char_rows = [row[1:] for row in sample_rows if row[0] == "!Sample_characteristics_ch1"]
    records = []
    for i, (title, gsm) in enumerate(zip(titles, accessions, strict=True)):
        fields: dict[str, str] = {}
        for row in char_rows:
            value = row[i].strip()
            if not value or ": " not in value:
                continue
            key, val = value.split(": ", 1)
            fields[key] = val
        recist = fields.get("anti-pd-1 response")
        if recist not in HUGO_RESPONSE:
            raise ValueError(f"Unexpected Hugo response for {title}: {recist}")
        patient = title.rstrip("AB")
        records.append(
            {
                "cohort": "Hugo_2016",
                "sample": title,
                "patient": patient,
                "gsm": gsm,
                "recist": HUGO_RESPONSE[recist],
                "recist_raw": recist,
                "os_days": fields.get("overall survival (days)"),
                "vital_status": fields.get("vital status"),
                "previous_mapki": fields.get("previous mapki"),
                "study_site": fields.get("study site"),
                "gender": fields.get("gender"),
                "age": fields.get("age (yrs)"),
            }
        )
    return pd.DataFrame(records)


def load_hugo() -> tuple[pd.DataFrame, pd.DataFrame]:
    meta = parse_hugo_metadata()
    expr = pd.read_excel(DATA / "GSE78220_PatientFPKM.xlsx")
    if expr.shape != (25268, 29) or expr.columns[0] != "Gene":
        raise ValueError(f"Unexpected Hugo FPKM shape {expr.shape}")
    missing = [g for g in ("TACSTD2", "CLDN4") if g not in set(expr["Gene"])]
    if missing:
        raise ValueError(f"Missing Hugo genes: {missing}")

    sample_rows = []
    for _, row in meta.iterrows():
        title = row["sample"]
        matches = [c for c in expr.columns[1:] if c.startswith(f"{title}.")]
        if len(matches) != 1:
            raise ValueError(f"Could not map Hugo sample {title} to FPKM column")
        col = matches[0]
        visit = "on_treatment" if col.endswith(".OnTx") else "pretreatment"
        recist = row["recist"]
        sample_rows.append(
            {
                **row.to_dict(),
                "expr_column": col,
                "visit": visit,
                "tacstd2_raw": float(expr.loc[expr["Gene"] == "TACSTD2", col].iloc[0]),
                "cldn4_raw": float(expr.loc[expr["Gene"] == "CLDN4", col].iloc[0]),
                "binary_recist": (
                    "Responder" if recist in {"CR", "PR"} else "Nonresponder"
                ),
            }
        )
    samples = pd.DataFrame(sample_rows)
    samples["tacstd2_log2"] = np.log2(samples["tacstd2_raw"] + 1)
    samples["cldn4_log2"] = np.log2(samples["cldn4_raw"] + 1)

    pre = samples[samples["visit"] == "pretreatment"].copy()
    patient = pre.groupby("patient", as_index=False).agg(
        cohort=("cohort", "first"),
        recist=("recist", "first"),
        binary_recist=("binary_recist", "first"),
        previous_mapki=("previous_mapki", "first"),
        os_days=("os_days", "first"),
        vital_status=("vital_status", "first"),
        n_biopsies=("sample", "size"),
        tacstd2_log2=("tacstd2_log2", "mean"),
        cldn4_log2=("cldn4_log2", "mean"),
        tacstd2_raw=("tacstd2_raw", "mean"),
        cldn4_raw=("cldn4_raw", "mean"),
    )
    if set(patient["recist"]) - {"CR", "PR", "PD"}:
        raise ValueError("Hugo patient-level RECIST has unexpected labels")
    return samples, patient


def load_vanallen() -> pd.DataFrame:
    tpm = pd.read_csv(DATA / "TPM_RSEM_VAScience2015.txt", sep="\t")
    if tpm.shape[1] != 43:
        raise ValueError(f"Expected 42 Van Allen RNA samples; found {tpm.shape[1] - 1}")
    gene_col = tpm.columns[0]
    tpm[gene_col] = tpm[gene_col].astype(str).str.replace(r"\.\d+$", "", regex=True)
    for symbol, ensg in VA_GENES.items():
        n = int((tpm[gene_col] == ensg).sum())
        if n != 1:
            raise ValueError(f"Expected one {symbol} ({ensg}) row; found {n}")

    clin = pd.read_csv(
        DATA / "skcm_dfci_2015_clinical_patient.txt", sep="\t", comment="#"
    )
    clin = clin.set_index("PATIENT_ID")

    rows = []
    for col in tpm.columns[1:]:
        if not col.startswith("MEL-IPI_"):
            raise ValueError(f"Unexpected Van Allen sample name: {col}")
        patient = col.replace("MEL-IPI_", "")
        if patient not in clin.index:
            recist = ""
            os_months = float("nan")
            os_status = ""
        else:
            recist = clin.loc[patient, "DURABLE_CLINICAL_BENEFIT"]
            os_months = clin.loc[patient, "OS_MONTHS"]
            os_status = clin.loc[patient, "OS_STATUS"]
        recist = "" if pd.isna(recist) else str(recist)
        if recist in {"CR", "PR"}:
            binary_recist = "Responder"
        elif recist == "PD":
            binary_recist = "Nonresponder"
        else:
            binary_recist = ""
        if recist in {"CR", "PR"} or (
            recist == "SD" and pd.notna(os_months) and float(os_months) > 12
        ):
            paper_cb = "Clinical_benefit"
        elif recist == "PD" or (
            recist == "SD" and pd.notna(os_months) and float(os_months) <= 12
        ):
            paper_cb = "No_clinical_benefit"
        else:
            paper_cb = ""
        expr = {
            symbol: float(tpm.loc[tpm[gene_col] == ensg, col].iloc[0])
            for symbol, ensg in VA_GENES.items()
        }
        rows.append(
            {
                "cohort": "VanAllen_2015",
                "sample": col,
                "patient": patient,
                "visit": "pretreatment",
                "recist": recist,
                "os_months": os_months,
                "os_status": os_status,
                "binary_recist": binary_recist,
                "paper_cb": paper_cb,
                "tacstd2_raw": expr["TACSTD2"],
                "cldn4_raw": expr["CLDN4"],
                "tacstd2_log2": math.log2(expr["TACSTD2"] + 1),
                "cldn4_log2": math.log2(expr["CLDN4"] + 1),
            }
        )
    out = pd.DataFrame(rows)
    if len(out) != 42:
        raise ValueError(f"Expected 42 Van Allen patients; found {len(out)}")
    return out


def groups(df: pd.DataFrame, gene_col: str, label_col: str, pos: str, neg: str):
    a = df.loc[df[label_col] == pos, gene_col].to_numpy(dtype=float)
    b = df.loc[df[label_col] == neg, gene_col].to_numpy(dtype=float)
    return a, b


def spearman_pair(df: pd.DataFrame) -> dict:
    rho, p = stats.spearmanr(df["tacstd2_log2"], df["cldn4_log2"])
    return {"n": int(len(df)), "spearman_rho": float(rho), "spearman_p": float(p)}


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_figure(panels: list[dict]) -> None:
    plt.style.use("seaborn-v0_8-whitegrid")
    fig, axes = plt.subplots(2, 2, figsize=(9.2, 7.4))
    rng = np.random.default_rng(SEED)
    colors = ["#6b7280", "#2563eb"]
    for ax, panel in zip(axes.ravel(), panels, strict=True):
        ordered = [panel["nr"], panel["r"]]
        bp = ax.boxplot(
            ordered,
            positions=[1, 2],
            widths=0.55,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for box, color in zip(bp["boxes"], colors, strict=True):
            box.set_facecolor(color)
            box.set_alpha(0.35)
        for pos, values, color in zip([1, 2], ordered, colors, strict=True):
            jitter = rng.uniform(-0.12, 0.12, len(values))
            ax.scatter(
                pos + jitter,
                values,
                s=26,
                color=color,
                alpha=0.85,
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
        ax.set_xticks([1, 2], [f"NR\nn={len(ordered[0])}", f"R\nn={len(ordered[1])}"])
        ax.set_title(panel["title"], fontsize=10, weight="bold")
        ax.set_ylabel(panel["ylabel"])
        ax.text(
            0.5,
            0.98,
            f"MW p={panel['p']:.3g}   AUC={panel['auc']:.2f}",
            transform=ax.transAxes,
            ha="center",
            va="top",
            fontsize=8.5,
        )
    fig.suptitle(
        "Melanoma ICI RNA: TACSTD2 and CLDN4 versus RECIST response",
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "tacstd2_cldn4_vs_response.png", dpi=220, bbox_inches="tight")
    fig.savefig(ROOT / "tacstd2_cldn4_vs_response.svg", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    download_inputs()
    hugo_samples, hugo_patients = load_hugo()
    vanallen = load_vanallen()

    analyses = []
    figure_panels = []

    def add_analysis(
        *,
        analysis: str,
        cohort: str,
        gene: str,
        endpoint: str,
        scale: str,
        r: np.ndarray,
        nr: np.ndarray,
        panel: bool = False,
        ylabel: str | None = None,
        title: str | None = None,
    ) -> dict:
        stats_row = summarize(r, nr)
        row = {
            "analysis": analysis,
            "cohort": cohort,
            "gene": gene,
            "endpoint": endpoint,
            "scale": scale,
            **stats_row,
        }
        analyses.append(row)
        if panel:
            figure_panels.append(
                {
                    "title": title,
                    "ylabel": ylabel,
                    "r": r,
                    "nr": nr,
                    "p": stats_row["mann_whitney_p"],
                    "auc": stats_row["auc_responder_higher"],
                }
            )
        return row

    for gene, col in (("TACSTD2", "tacstd2_log2"), ("CLDN4", "cldn4_log2")):
        r, nr = groups(hugo_patients, col, "binary_recist", "Responder", "Nonresponder")
        add_analysis(
            analysis=f"hugo_pre_patient_{gene.lower()}_recist",
            cohort="Hugo_2016",
            gene=gene,
            endpoint="RECIST_CRPR_vs_PD",
            scale="log2(FPKM+1)",
            r=r,
            nr=nr,
            panel=True,
            ylabel="log2(FPKM + 1)",
            title=f"Hugo 2016 anti–PD-1  {gene}",
        )

    hugo_no_pt28 = hugo_patients[hugo_patients["patient"] != "Pt28"]
    for gene, col in (("TACSTD2", "tacstd2_log2"), ("CLDN4", "cldn4_log2")):
        r, nr = groups(hugo_no_pt28, col, "binary_recist", "Responder", "Nonresponder")
        add_analysis(
            analysis=f"hugo_pre_patient_{gene.lower()}_recist_drop_Pt28",
            cohort="Hugo_2016",
            gene=gene,
            endpoint="RECIST_CRPR_vs_PD_drop_Pt28_outlier",
            scale="log2(FPKM+1)",
            r=r,
            nr=nr,
        )

    hugo_pre_samples = hugo_samples[hugo_samples["visit"] == "pretreatment"]
    for gene, col in (("TACSTD2", "tacstd2_log2"), ("CLDN4", "cldn4_log2")):
        r, nr = groups(hugo_pre_samples, col, "binary_recist", "Responder", "Nonresponder")
        add_analysis(
            analysis=f"hugo_pre_sample_{gene.lower()}_recist",
            cohort="Hugo_2016",
            gene=gene,
            endpoint="RECIST_CRPR_vs_PD_sample_level",
            scale="log2(FPKM+1)",
            r=r,
            nr=nr,
        )

    for gene, col in (("TACSTD2", "tacstd2_log2"), ("CLDN4", "cldn4_log2")):
        r, nr = groups(vanallen, col, "binary_recist", "Responder", "Nonresponder")
        add_analysis(
            analysis=f"vanallen_pre_{gene.lower()}_recist",
            cohort="VanAllen_2015",
            gene=gene,
            endpoint="RECIST_CRPR_vs_PD",
            scale="log2(TPM+1)",
            r=r,
            nr=nr,
            panel=True,
            ylabel="log2(TPM + 1)",
            title=f"Van Allen 2015 anti–CTLA-4  {gene}",
        )
        r, nr = groups(
            vanallen, col, "paper_cb", "Clinical_benefit", "No_clinical_benefit"
        )
        add_analysis(
            analysis=f"vanallen_pre_{gene.lower()}_paperCB",
            cohort="VanAllen_2015",
            gene=gene,
            endpoint="paper_clinical_benefit",
            scale="log2(TPM+1)",
            r=r,
            nr=nr,
        )

    make_figure(figure_panels)
    write_csv(ROOT / "summary.csv", analyses)

    hugo_out = hugo_samples.copy()
    hugo_out["paper_cb"] = ""
    hugo_out["os_months"] = ""
    va_out = vanallen.copy()
    va_out["gsm"] = ""
    va_out["expr_column"] = va_out["sample"]
    va_out["recist_raw"] = va_out["recist"]
    va_out["previous_mapki"] = ""
    va_out["study_site"] = ""
    va_out["gender"] = ""
    va_out["age"] = ""
    va_out["os_days"] = ""
    va_out["vital_status"] = ""
    cols = [
        "cohort",
        "sample",
        "patient",
        "gsm",
        "expr_column",
        "visit",
        "recist",
        "recist_raw",
        "binary_recist",
        "paper_cb",
        "tacstd2_raw",
        "cldn4_raw",
        "tacstd2_log2",
        "cldn4_log2",
        "previous_mapki",
        "os_days",
        "os_months",
        "os_status" if "os_status" in va_out.columns else "vital_status",
        "vital_status",
        "study_site",
        "gender",
        "age",
    ]
    hugo_out["os_status"] = hugo_out["vital_status"]
    sample_level = pd.concat([hugo_out[cols], va_out[cols]], ignore_index=True)
    sample_level.to_csv(ROOT / "sample_level.csv", index=False)

    hugo_patients.to_csv(ROOT / "hugo_patient_level.csv", index=False)

    corr_rows = [
        {"cohort": "Hugo_2016_pretreatment_patients", **spearman_pair(hugo_patients)},
        {
            "cohort": "VanAllen_2015_pretreatment",
            **spearman_pair(vanallen),
        },
    ]
    write_csv(ROOT / "gene_correlation.csv", corr_rows)

    with (ROOT / "provenance.tsv").open("w") as handle:
        handle.write("file\tsha256\turl\n")
        for path, url in URLS.items():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            handle.write(f"{path.name}\t{digest}\t{url}\n")

    print("Hugo sample RECIST:", dict(Counter(hugo_samples["recist"])))
    print("Hugo sample visit:", dict(Counter(hugo_samples["visit"])))
    print("Hugo patient RECIST:", dict(Counter(hugo_patients["recist"])))
    print("Van Allen RECIST:", dict(Counter(vanallen["recist"])))
    print("Van Allen paper CB:", dict(Counter(vanallen["paper_cb"])))
    for row in analyses:
        print(
            f"{row['analysis']}: R={row['n_responder']} NR={row['n_nonresponder']} "
            f"AUC={row['auc_responder_higher']:.3f} "
            f"p={row['mann_whitney_p']:.4g}"
        )
    for row in corr_rows:
        print(
            f"{row['cohort']}: TACSTD2–CLDN4 Spearman ρ={row['spearman_rho']:.3f} "
            f"p={row['spearman_p']:.4g}"
        )


if __name__ == "__main__":
    main()
