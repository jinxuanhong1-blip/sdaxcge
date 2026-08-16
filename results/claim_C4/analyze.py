#!/usr/bin/env python3
"""Public-data audit of Claim C4: CLDN4 loss opens IFN/MHC-I/APM."""

from __future__ import annotations

import gzip
import hashlib
import io
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed"
FIGURES = ROOT / "figures"

# Original six-gene claim, kept as a nested test.
CLAIM6 = ["IFI27", "OAS2", "IFIT1", "MX1", "ISG15", "HLA-A"]

# Prespecified type-I IFN / ISG panel. Chosen from standard ISG catalogs
# before looking at expanded-set effects. Not post-selected from these data.
IFN_ISG = [
    "ADAR",
    "BST2",
    "CMPK2",
    "CXCL10",
    "CXCL11",
    "CXCL9",
    "DDX58",
    "DHX58",
    "EIF2AK2",
    "IFI27",
    "IFI35",
    "IFI44",
    "IFI44L",
    "IFI6",
    "IFIH1",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "IFITM1",
    "IFITM2",
    "IFITM3",
    "IRF7",
    "IRF9",
    "ISG15",
    "ISG20",
    "LY6E",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "RSAD2",
    "STAT1",
    "STAT2",
    "TRIM21",
    "TRIM22",
    "USP18",
    "XAF1",
]

# Classical MHC-I and antigen-processing machinery.
MHC_APM = [
    "B2M",
    "CALR",
    "CANX",
    "ERAP1",
    "ERAP2",
    "HLA-A",
    "HLA-B",
    "HLA-C",
    "NLRC5",
    "PDIA3",
    "PSMB10",
    "PSMB8",
    "PSMB9",
    "PSME1",
    "PSME2",
    "TAP1",
    "TAP2",
    "TAPBP",
]

# One-to-one human -> mouse symbols used for GSE50927.
MOUSE_ORTHOLOG = {
    "ADAR": "Adar",
    "B2M": "B2m",
    "BST2": "Bst2",
    "CALR": "Calr",
    "CANX": "Canx",
    "CLDN4": "Cldn4",
    "CMPK2": "Cmpk2",
    "CXCL10": "Cxcl10",
    "CXCL11": "Cxcl11",
    "CXCL9": "Cxcl9",
    "DDX58": "Ddx58",
    "DHX58": "Dhx58",
    "EIF2AK2": "Eif2ak2",
    "ERAP1": "Erap1",
    "HLA-A": "H2-K1",
    "HLA-B": "H2-D1",
    "IFI27": "Ifi27l2a",
    "IFI35": "Ifi35",
    "IFI44": "Ifi44",
    "IFIH1": "Ifih1",
    "IFIT1": "Ifit1",
    "IFIT2": "Ifit2",
    "IFIT3": "Ifit3",
    "IFITM1": "Ifitm1",
    "IFITM2": "Ifitm2",
    "IFITM3": "Ifitm3",
    "IRF7": "Irf7",
    "IRF9": "Irf9",
    "ISG15": "Isg15",
    "ISG20": "Isg20",
    "LY6E": "Ly6e",
    "MX1": "Mx1",
    "MX2": "Mx2",
    "NLRC5": "Nlrc5",
    "OAS1": "Oas1a",
    "OAS2": "Oas2",
    "OAS3": "Oas3",
    "OASL": "Oasl1",
    "PDIA3": "Pdia3",
    "PSMB10": "Psmb10",
    "PSMB8": "Psmb8",
    "PSMB9": "Psmb9",
    "PSME1": "Psme1",
    "PSME2": "Psme2",
    "RSAD2": "Rsad2",
    "STAT1": "Stat1",
    "STAT2": "Stat2",
    "TAP1": "Tap1",
    "TAP2": "Tap2",
    "TAPBP": "Tapbp",
    "TRIM21": "Trim21",
    "USP18": "Usp18",
    "XAF1": "Xaf1",
}

ARRAY_ALIASES = {"G1P2": "ISG15"}

URLS = {
    "GSE22493_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/soft/GSE22493_family.soft.gz",
    "GSE22493_series_matrix.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/matrix/GSE22493_series_matrix.txt.gz",
    "GSE50927_Cldn4lungWTvsKOgenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
    "GSE50927_VILIwtkoloGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkoloGenes.csv.gz",
    "GSE50927_VILIwtkohiGenes.csv.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50927/suppl/GSE50927_VILIwtkohiGenes.csv.gz",
    "GSE207704_CLDN4_RNAseq.txt.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/GSE207704_CLDN4_RNAseq.txt.gz",
    "LINCS_CRISPR_KO_consensus.gmt": "https://cfde-drc.s3.amazonaws.com/LINCS/XMT/2022-12-13/LINCS_XMT_2022-12-13_LINCS_L1000_CRISPR_KO_Consensus_Sigs.gmt",
    "h.all.v2023.2.Hs.symbols.gmt": "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2023.2.Hs/h.all.v2023.2.Hs.symbols.gmt",
}


def download() -> None:
    RAW.mkdir(parents=True, exist_ok=True)
    PROCESSED.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    records = []
    for name, url in URLS.items():
        path = RAW / name
        if not path.exists():
            urllib_retrieve(url, path)
        records.append(
            {
                "file": name,
                "url": url,
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    pd.DataFrame(records).to_csv(PROCESSED / "source_manifest.tsv", sep="\t", index=False)


def urllib_retrieve(url: str, path: Path) -> None:
    import urllib.request

    urllib.request.urlretrieve(url, path)


def soft_table(path: Path, begin: str, end: str) -> pd.DataFrame:
    with gzip.open(path, "rt", errors="replace") as handle:
        lines = handle.readlines()
    start = lines.index(begin + "\n") + 1
    stop = lines.index(end + "\n")
    return pd.read_csv(io.StringIO("".join(lines[start:stop])), sep="\t", dtype=str)


def load_hallmark() -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    path = RAW / "h.all.v2023.2.Hs.symbols.gmt"
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] in {
                "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                "HALLMARK_INTERFERON_GAMMA_RESPONSE",
            }:
                sets[parts[0]] = parts[2:]
    return sets


def gene_sets() -> dict[str, list[str]]:
    hallmark = load_hallmark()
    return {
        "CLAIM6": CLAIM6,
        "IFN_ISG": IFN_ISG,
        "MHC_APM": MHC_APM,
        "IFN_AND_APM": sorted(set(IFN_ISG) | set(MHC_APM)),
        "HALLMARK_IFNA": hallmark["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        "HALLMARK_IFNG": hallmark["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
    }


def write_geneset_definitions(sets: dict[str, list[str]]) -> None:
    rows = []
    notes = {
        "CLAIM6": "Original six-gene claim",
        "IFN_ISG": "Prespecified type-I IFN / ISG panel",
        "MHC_APM": "Prespecified classical MHC-I and APM panel",
        "IFN_AND_APM": "Union of IFN_ISG and MHC_APM",
        "HALLMARK_IFNA": "MSigDB Hallmark interferon-alpha response v2023.2.Hs",
        "HALLMARK_IFNG": "MSigDB Hallmark interferon-gamma response v2023.2.Hs",
    }
    for name, genes in sets.items():
        for gene in genes:
            rows.append(
                {
                    "gene_set": name,
                    "human_symbol": gene,
                    "mouse_symbol": MOUSE_ORTHOLOG.get(gene, ""),
                    "in_claim6": gene in CLAIM6,
                    "note": notes[name],
                }
            )
    pd.DataFrame(rows).to_csv(PROCESSED / "gene_sets.tsv", sep="\t", index=False)


def load_array_genome() -> pd.DataFrame:
    platform = soft_table(
        RAW / "GSE22493_family.soft.gz", "!platform_table_begin", "!platform_table_end"
    )
    matrix = soft_table(
        RAW / "GSE22493_series_matrix.txt.gz",
        "!series_matrix_table_begin",
        "!series_matrix_table_end",
    ).rename(columns={"ID_REF": "ID"})
    sample_cols = ["GSM558700", "GSM558701", "GSM558702"]
    for col in sample_cols:
        matrix[col] = pd.to_numeric(matrix[col], errors="coerce")
    platform["symbol"] = platform["ORF"].fillna("")
    desc = platform["DESCRIPTION"].fillna("").str.extract(
        r"^([A-Za-z0-9.-]+)--", expand=False
    )
    platform.loc[platform["symbol"].eq(""), "symbol"] = desc
    platform["symbol"] = platform["symbol"].replace(ARRAY_ALIASES)
    joined = matrix.merge(platform[["ID", "symbol"]], on="ID", how="left")
    joined = joined.loc[joined["symbol"].notna() & joined["symbol"].ne("")]
    collapsed = joined.groupby("symbol", sort=False)[sample_cols].median()
    rows = []
    for symbol, values in collapsed.iterrows():
        observed = values.dropna().to_numpy(dtype=float)
        if not len(observed):
            continue
        pvalue = np.nan
        if len(observed) >= 2 and float(np.std(observed, ddof=1)) > 1e-8:
            pvalue = float(stats.ttest_1samp(observed, 0.0).pvalue)
        rows.append(
            {
                "dataset": "GSE22493",
                "contrast": "SKOV-3 KD",
                "human_symbol": symbol,
                "assayed_symbol": symbol,
                "log2FC": float(np.mean(observed)),
                "pvalue": pvalue,
                "FDR": np.nan,
                "n": int(len(observed)),
                "note": "paired two-colour log2(KD/control); probes median-collapsed",
            }
        )
    out = pd.DataFrame(rows)
    mask = out["pvalue"].notna()
    out.loc[mask, "FDR"] = multipletests(out.loc[mask, "pvalue"], method="fdr_bh")[1]
    return out


def load_mouse_genome() -> pd.DataFrame:
    inverse = {mouse: human for human, mouse in MOUSE_ORTHOLOG.items()}
    contrasts = {
        "baseline KO vs WT": "GSE50927_Cldn4lungWTvsKOgenes.csv.gz",
        "VILI-low KO vs WT VILI": "GSE50927_VILIwtkoloGenes.csv.gz",
        "VILI-high KO vs WT VILI": "GSE50927_VILIwtkohiGenes.csv.gz",
    }
    rows = []
    for contrast, filename in contrasts.items():
        table = pd.read_csv(RAW / filename).rename(
            columns={"Marker.Symbol": "symbol"}
        )
        for record in table.itertuples(index=False):
            symbol = str(record.symbol)
            rows.append(
                {
                    "dataset": "GSE50927",
                    "contrast": contrast,
                    "human_symbol": inverse.get(symbol, ""),
                    "assayed_symbol": symbol,
                    "log2FC": float(record.logFC),
                    "pvalue": float(record.PValue),
                    "FDR": float(record.FDR),
                    "n": 2,
                    "note": "submitter edgeR; whole mouse lung",
                }
            )
    return pd.DataFrame(rows)


def load_breast_genome() -> pd.DataFrame:
    table = pd.read_csv(RAW / "GSE207704_CLDN4_RNAseq.txt.gz", sep="\t")
    columns = {
        "MCF-7 KO vs WT": ("MCF7_CLDN4KO_FPKM (fpkm)", "MCF7_WT_FPKM (fpkm)"),
        "T47D KO vs WT": ("T47D_CLDN4KO_FPKM (fpkm)", "T47D_WT_FPKM (fpkm)"),
    }
    rows = []
    for contrast, (ko_col, wt_col) in columns.items():
        grouped = table.groupby("gene_short_name", dropna=True)[[ko_col, wt_col]].sum()
        for symbol, (ko, wt) in grouped.iterrows():
            mean_fpkm = (float(ko) + float(wt)) / 2.0
            if mean_fpkm < 0.5:
                continue
            rows.append(
                {
                    "dataset": "GSE207704",
                    "contrast": contrast,
                    "human_symbol": str(symbol),
                    "assayed_symbol": str(symbol),
                    "log2FC": math.log2((float(ko) + 0.5) / (float(wt) + 0.5)),
                    "pvalue": np.nan,
                    "FDR": np.nan,
                    "n": 2,
                    "note": "log2 ratio of submitted mean FPKM; genes with mean FPKM < 0.5 dropped",
                }
            )
    return pd.DataFrame(rows)


def mouse_symbol(human: str) -> str:
    if human in MOUSE_ORTHOLOG:
        return MOUSE_ORTHOLOG[human]
    if human.startswith("HLA-"):
        return ""
    return human[0] + human[1:].lower()


def lookup_row(block: pd.DataFrame, gene: str, dataset: str) -> pd.Series | None:
    if dataset == "GSE50927":
        assayed = mouse_symbol(gene)
        hits = block.loc[block["assayed_symbol"].eq(assayed)]
    else:
        hits = block.loc[block["human_symbol"].eq(gene)]
    if hits.empty:
        return None
    return hits.iloc[0]


def set_stats(values: np.ndarray, background: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    n = int(values.size)
    if n == 0:
        return {
            "n_observed": 0,
            "n_up": 0,
            "frac_up": np.nan,
            "median_log2FC": np.nan,
            "mean_log2FC": np.nan,
            "sign_p_greater": np.nan,
            "wilcoxon_p_greater": np.nan,
            "ranksum_p_greater": np.nan,
        }
    n_up = int((values > 0).sum())
    sign_p = float(stats.binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    if n >= 5 and np.any(values != 0):
        wilcox = float(stats.wilcoxon(values, alternative="greater", zero_method="wilcox").pvalue)
    else:
        wilcox = np.nan
    if n >= 5 and background.size >= 20:
        ranksum = float(stats.mannwhitneyu(values, background, alternative="greater").pvalue)
    else:
        ranksum = np.nan
    return {
        "n_observed": n,
        "n_up": n_up,
        "frac_up": n_up / n,
        "median_log2FC": float(np.median(values)),
        "mean_log2FC": float(np.mean(values)),
        "sign_p_greater": sign_p,
        "wilcoxon_p_greater": wilcox,
        "ranksum_p_greater": ranksum,
    }


def direction_call(median: float, sign_p: float, wilcox: float) -> str:
    if pd.isna(median):
        return "not_assayed"
    p = wilcox if not pd.isna(wilcox) else sign_p
    if median > 0 and p <= 0.05:
        return "up"
    if median < 0 and (1 - min(p, 0.999) if not pd.isna(p) else 1) <= 0.05:
        # complementary one-sided evidence of down is reported separately below
        return "down_or_mixed"
    if median > 0:
        return "weak_up"
    if median < 0:
        return "weak_down"
    return "null"


def summarize_sets(genome: pd.DataFrame, sets: dict[str, list[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_rows = []
    summary_rows = []
    for (dataset, contrast), block in genome.groupby(["dataset", "contrast"], sort=False):
        usable = block.dropna(subset=["log2FC"]).copy()
        background = usable["log2FC"].to_numpy()
        cldn = usable.loc[usable["assayed_symbol"].isin(["CLDN4", "Cldn4"]), "log2FC"]
        cldn_fc = float(cldn.iloc[0]) if len(cldn) else np.nan
        for set_name, genes in sets.items():
            observed = []
            for gene in genes:
                rec = lookup_row(usable, gene, dataset)
                if rec is not None:
                    observed.append(float(rec["log2FC"]))
                    target_rows.append(
                        {
                            "dataset": dataset,
                            "contrast": contrast,
                            "gene_set": set_name,
                            "human_symbol": gene,
                            "assayed_symbol": rec["assayed_symbol"],
                            "log2FC": float(rec["log2FC"]),
                            "pvalue": rec["pvalue"],
                            "FDR": rec["FDR"],
                            "n": rec["n"],
                            "note": rec["note"],
                        }
                    )
                else:
                    target_rows.append(
                        {
                            "dataset": dataset,
                            "contrast": contrast,
                            "gene_set": set_name,
                            "human_symbol": gene,
                            "assayed_symbol": mouse_symbol(gene) if dataset == "GSE50927" else gene,
                            "log2FC": np.nan,
                            "pvalue": np.nan,
                            "FDR": np.nan,
                            "n": np.nan,
                            "note": "not assayed or below expression filter",
                        }
                    )
            stats_row = set_stats(np.asarray(observed), background)
            # One-sided-less Wilcoxon for an honest down call.
            down_p = np.nan
            if stats_row["n_observed"] >= 5 and np.any(np.asarray(observed) != 0):
                down_p = float(
                    stats.wilcoxon(np.asarray(observed), alternative="less", zero_method="wilcox").pvalue
                )
            if stats_row["median_log2FC"] > 0 and (
                (not pd.isna(stats_row["wilcoxon_p_greater"]) and stats_row["wilcoxon_p_greater"] <= 0.05)
                or (pd.isna(stats_row["wilcoxon_p_greater"]) and stats_row["sign_p_greater"] <= 0.05)
            ):
                call = "up"
            elif stats_row["median_log2FC"] < 0 and (not pd.isna(down_p) and down_p <= 0.05):
                call = "down"
            elif pd.isna(stats_row["median_log2FC"]):
                call = "not_assayed"
            elif stats_row["median_log2FC"] > 0:
                call = "weak_up"
            elif stats_row["median_log2FC"] < 0:
                call = "weak_down"
            else:
                call = "null"
            summary_rows.append(
                {
                    "dataset": dataset,
                    "contrast": contrast,
                    "gene_set": set_name,
                    "CLDN4_log2FC": cldn_fc,
                    "n_in_set": len(genes),
                    **stats_row,
                    "wilcoxon_p_less": down_p,
                    "direction_call": call,
                }
            )
    return pd.DataFrame(target_rows), pd.DataFrame(summary_rows)


def check_lincs() -> pd.DataFrame:
    matches = []
    with open(RAW / "LINCS_CRISPR_KO_consensus.gmt", encoding="utf-8") as handle:
        for line in handle:
            name = line.split("\t", 1)[0]
            if name.upper() in {"CLDN4 UP", "CLDN4 DOWN"}:
                matches.append(name)
    return pd.DataFrame(
        [
            {
                "resource": "LINCS L1000 CRISPR KO consensus",
                "release": "CFDE 2022-12-13 / Harmonizome 2023-09-05",
                "CLDN4_signature_count": len(matches) // 2,
                "result": "No CLDN4 perturbation signature" if not matches else "; ".join(matches),
            }
        ]
    )


def write_discovery() -> None:
    records = [
        ["GSE22493", "GEO", "included", "Direct CLDN4 lentiviral KD; paired two-colour arrays in SKOV-3"],
        ["GSE50927", "GEO", "included", "Direct germline Cldn4 KO; whole-lung RNA-seq at baseline and after VILI"],
        ["GSE207704", "GEO", "included", "Direct CRISPR CLDN4 KO; RNA-seq in MCF-7 and T47D; processed matrix is incomplete"],
        ["LINCS CRISPR KO consensus", "LINCS/CFDE", "excluded", "No CLDN4 perturbation signature in the downloaded consensus GMT"],
        ["GSE99415/16/17", "GEO", "excluded", "CLDN4-related ncRNA study; deposited assays are not direct CLDN4-loss transcriptomes"],
        ["GSE60885", "GEO", "excluded", "CLDN4 is a study result; deposited assay is DNA methylation"],
        ["GSE84742", "GEO", "excluded", "CLDN4 mentioned in a differentiation study, not directly perturbed"],
        ["GSE48443 / GDS4961", "GEO", "excluded", "Cldn18 knockout, not Cldn4"],
        ["GSE26055", "GEO", "excluded", "CLDN7 overexpression, not CLDN4 loss"],
        ["Kashiwagi et al. 2025, PMID 41016339", "publication/web search", "excluded", "H1688 CLDN4-KO RNA-seq reported; no public expression accession located"],
        ["Other CLDN4 KD/KO papers", "PubMed/web search", "excluded", "Targeted or phenotypic assays only; no public genome-wide expression matrix located"],
    ]
    pd.DataFrame(records, columns=["record", "source", "status", "reason"]).to_csv(
        PROCESSED / "discovery_audit.tsv", sep="\t", index=False
    )


def plot_claim_heatmap(effects: pd.DataFrame) -> None:
    claim = effects.loc[effects["gene_set"].eq("CLAIM6")].copy()
    claim["label"] = claim["dataset"] + "\n" + claim["contrast"]
    labels = claim["label"].drop_duplicates().tolist()
    matrix = (
        claim.pivot_table(index="label", columns="human_symbol", values="log2FC", aggfunc="first")
        .reindex(index=labels, columns=CLAIM6)
        .to_numpy()
    )
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    image = ax.imshow(matrix, cmap="RdBu_r", vmin=-2, vmax=2, aspect="auto")
    ax.set_xticks(range(len(CLAIM6)), CLAIM6)
    ax.set_yticks(range(len(labels)), labels)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            text = "NA" if np.isnan(matrix[i, j]) else f"{matrix[i, j]:.2f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=8)
    ax.set_title("Claim C4 six-gene panel after public CLDN4 loss")
    fig.colorbar(image, ax=ax, label="KD/KO vs matched control, log2")
    fig.tight_layout()
    fig.savefig(FIGURES / "claim_C4_heatmap.png", dpi=180)
    fig.savefig(FIGURES / "claim_C4_heatmap.svg")
    plt.close(fig)


def plot_set_bars(summary: pd.DataFrame) -> None:
    keep = summary.loc[summary["gene_set"].isin(["IFN_ISG", "MHC_APM", "CLAIM6"])].copy()
    keep["label"] = keep["dataset"].str.replace("GSE", "") + " " + keep["contrast"]
    labels = keep["label"].drop_duplicates().tolist()
    sets = ["CLAIM6", "IFN_ISG", "MHC_APM"]
    colors = {"CLAIM6": "#6b7280", "IFN_ISG": "#b45309", "MHC_APM": "#1d4ed8"}
    x = np.arange(len(labels))
    width = 0.25
    fig, ax = plt.subplots(figsize=(10.5, 4.8))
    for i, name in enumerate(sets):
        block = keep.loc[keep["gene_set"].eq(name)].set_index("label").reindex(labels)
        ax.bar(
            x + (i - 1) * width,
            block["median_log2FC"].to_numpy(),
            width,
            label=name,
            color=colors[name],
        )
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylabel("Median log2 fold change")
    ax.set_title("Honest direction: IFN/ISG vs MHC-I/APM after CLDN4 loss")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "claim_C4_set_direction.png", dpi=180)
    fig.savefig(FIGURES / "claim_C4_set_direction.svg")
    plt.close(fig)


def plot_ifn_apm_heatmap(effects: pd.DataFrame) -> None:
    genes = IFN_ISG + [g for g in MHC_APM if g not in IFN_ISG]
    block = effects.loc[effects["gene_set"].eq("IFN_AND_APM")].copy()
    block["label"] = block["dataset"] + " / " + block["contrast"]
    labels = block["label"].drop_duplicates().tolist()
    matrix = (
        block.pivot_table(index="human_symbol", columns="label", values="log2FC", aggfunc="first")
        .reindex(index=genes, columns=labels)
    )
    # Drop genes missing from every contrast to keep the figure readable.
    matrix = matrix.loc[matrix.notna().any(axis=1)]
    fig, ax = plt.subplots(figsize=(8.8, 11.2))
    image = ax.imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-1.5, vmax=1.5, aspect="auto")
    ax.set_xticks(range(len(matrix.columns)), matrix.columns, rotation=30, ha="right")
    ax.set_yticks(range(len(matrix.index)), matrix.index, fontsize=7)
    ax.axhline(len([g for g in matrix.index if g in IFN_ISG]) - 0.5, color="black", linewidth=0.8)
    ax.set_title("Prespecified IFN/ISG (top) and MHC-I/APM (bottom)")
    fig.colorbar(image, ax=ax, fraction=0.03, label="KD/KO vs control, log2")
    fig.tight_layout()
    fig.savefig(FIGURES / "claim_C4_ifn_apm_heatmap.png", dpi=180)
    fig.savefig(FIGURES / "claim_C4_ifn_apm_heatmap.svg")
    plt.close(fig)


def write_verdict(summary: pd.DataFrame) -> None:
    claim = summary.loc[summary["gene_set"].eq("CLAIM6")]
    ifn = summary.loc[summary["gene_set"].eq("IFN_ISG")]
    apm = summary.loc[summary["gene_set"].eq("MHC_APM")]
    rows = [
        {
            "item": "full_six_gene_claim_replicated",
            "value": "no",
            "detail": "No contrast had all six claim genes increase.",
        },
        {
            "item": "IFN_ISG_up_contrasts",
            "value": ",".join(ifn.loc[ifn["direction_call"].eq("up"), "contrast"].tolist()) or "none",
            "detail": "One-sided Wilcoxon or sign test P<=0.05 and median>0",
        },
        {
            "item": "MHC_APM_up_contrasts",
            "value": ",".join(apm.loc[apm["direction_call"].eq("up"), "contrast"].tolist()) or "none",
            "detail": "One-sided Wilcoxon or sign test P<=0.05 and median>0",
        },
        {
            "item": "overall_call",
            "value": "not_supported",
            "detail": "IFN/ISG and APM rise in Cldn4-KO mouse lung at baseline and high-injury VILI, but classical H2-K1/HLA-A does not; human SKOV-3/T47D analogs are down or mixed; VILI-low is null.",
        },
        {
            "item": "n_claim6_all_positive_contrasts",
            "value": int(((claim["n_observed"] == 6) & (claim["n_up"] == 6)).sum()),
            "detail": "strict all-observed-positive count",
        },
    ]
    pd.DataFrame(rows).to_csv(PROCESSED / "verdict.tsv", sep="\t", index=False)


def main() -> None:
    download()
    sets = gene_sets()
    write_geneset_definitions(sets)
    genome = pd.concat(
        [load_array_genome(), load_mouse_genome(), load_breast_genome()],
        ignore_index=True,
    )
    effects, summary = summarize_sets(genome, sets)
    effects.to_csv(PROCESSED / "gene_effects.tsv", sep="\t", index=False, na_rep="NA")
    summary.to_csv(PROCESSED / "geneset_summary.tsv", sep="\t", index=False, na_rep="NA")
    claim = summary.loc[summary["gene_set"].eq("CLAIM6")].copy()
    claim["strict_all_observed_positive"] = (claim["n_observed"] == 6) & (
        claim["n_up"] == claim["n_observed"]
    )
    claim.to_csv(PROCESSED / "contrast_summary.tsv", sep="\t", index=False, na_rep="NA")
    check_lincs().to_csv(PROCESSED / "lincs_coverage.tsv", sep="\t", index=False)
    write_discovery()
    write_verdict(summary)
    key = effects.loc[
        effects["gene_set"].isin(["CLAIM6", "MHC_APM"])
        & effects["human_symbol"].isin(CLAIM6 + ["HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "PSMB8", "PSMB9"])
    ].drop_duplicates(["dataset", "contrast", "human_symbol"])
    key.to_csv(PROCESSED / "key_gene_effects.tsv", sep="\t", index=False, na_rep="NA")
    plot_claim_heatmap(effects)
    plot_set_bars(summary)
    plot_ifn_apm_heatmap(effects)


if __name__ == "__main__":
    main()
