#!/usr/bin/env python3
"""Signed NHEJ gene set from Yamamoto et al., Mol Cancer Ther 2022, scored in
public CLDN4-loss transcriptomes GSE207704 and GSE22493.

The 2022 paper did not deposit a knockdown transcriptome or the RPPA matrix.
Table S1 is a TCGA ovarian Firehose reanalysis (CLDN4-high vs CLDN4-low).
Signs for the NHEJ set come from that table (q < 0.05). TP53BP1 and XRCC1,
the only proteins the paper states fall after CLDN4 knockdown, are scored
as a separate pre-specified effector pair.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import urllib.request
import warnings
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "input"
CACHE = ROOT / "cache"
OUT = ROOT / "output"
TABLE_S1 = INPUT / "yamamoto2022_table_s1.tsv"

# Classical NHEJ (KEGG 2021 Human "Non-homologous end-joining", hsa03450)
# plus 53BP1-axis factors that Reactome R-HSA-5693571 places in NHEJ and that
# are not histone peptides or the BRCA1-A complex (those antagonize NHEJ).
# XRCC1 is the paper's second RPPA effector; it is SSBR/BER, not classical NHEJ.
NHEJ_MEMBERSHIP = {
    "DCLRE1C": "KEGG_NHEJ",
    "DNTT": "KEGG_NHEJ",
    "FEN1": "KEGG_NHEJ_also_altNHEJ",
    "LIG4": "KEGG_NHEJ",
    "MRE11": "KEGG_NHEJ",
    "NHEJ1": "KEGG_NHEJ",
    "POLL": "KEGG_NHEJ",
    "POLM": "KEGG_NHEJ",
    "PRKDC": "KEGG_NHEJ",
    "RAD50": "KEGG_NHEJ",
    "XRCC4": "KEGG_NHEJ",
    "XRCC5": "KEGG_NHEJ",
    "XRCC6": "KEGG_NHEJ",
    "TP53BP1": "Reactome_NHEJ_53BP1",
    "RIF1": "Reactome_NHEJ_53BP1",
    "RNF8": "Reactome_NHEJ_DSB_response",
    "RNF168": "Reactome_NHEJ_DSB_response",
    "H2AX": "Reactome_NHEJ_DSB_response",
    "MDC1": "Reactome_NHEJ_DSB_response",
    "PAXIP1": "Reactome_NHEJ_DSB_response",
    "ATM": "Reactome_NHEJ_DSB_response",
    "NBN": "Reactome_NHEJ_MRN",
    "TDP1": "Reactome_NHEJ",
    "TDP2": "Reactome_NHEJ",
    "MAD2L2": "shieldin",
    "SHLD1": "shieldin",
    "SHLD2": "shieldin",
    "SHLD3": "shieldin",
    "XRCC1": "paper_RPPA_SSBR_not_classical_NHEJ",
}

# Symbols collapsed onto the Table S1 / HGNC symbol used in NHEJ_MEMBERSHIP.
ALIAS = {
    "H2AFX": "H2AX",
    "MRE11A": "MRE11",
    "NBS1": "NBN",
    "G22P1": "XRCC6",
    "TTRAP": "TDP2",
    "WHSC1": "NSD2",
    "FAM175A": "ABRAXAS1",
    "C20ORF196": "SHLD1",
    "FAM35A": "SHLD2",
    "RINN1": "SHLD3",
    "REV7": "MAD2L2",
    "TMEM173": "STING1",
    "53BP1": "TP53BP1",
}

PAPER_NAMED_REPAIR = ["POLR2J", "BRIP1", "NUP160", "PRKDC", "TP53BP1", "XRCC1"]
RPPA_DOWN = {"TP53BP1", "XRCC1"}  # text: protein reduced after CLDN4 shRNA
Q_CUTOFF = 0.05
N_PERM = 10000
RNG_SEED = 20220401
FPKM_FLOOR = 1.0
FPKM_PSEUDO = 0.1

GSE207704_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/"
    "GSE207704_CLDN4_RNAseq.txt.gz"
)
GSE22493_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE22nnn/GSE22493/matrix/"
    "GSE22493_series_matrix.txt.gz"
)
GSM558700_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM558nnn/GSM558700/suppl/"
    "GSM558700.txt.gz"
)


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 1000:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as resp:
        dest.write_bytes(resp.read())


def canon(symbol: str) -> str:
    s = symbol.strip().upper()
    return ALIAS.get(s, s)


def load_table_s1(path: Path) -> dict[str, dict]:
    out = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene = canon(row["gene"])
            rec = {
                "gene": gene,
                "log2_ratio_high_over_low": float(row["log2_ratio_high_over_low"]),
                "p_value": float(row["p_value"]),
                "q_value": float(row["q_value"]),
                "higher_expression_in": row["higher_expression_in"].strip(),
            }
            # Keep the more significant row if an alias collides.
            prev = out.get(gene)
            if prev is None or rec["q_value"] < prev["q_value"]:
                out[gene] = rec
    return out


def sign_from_table(rec: dict) -> int:
    higher = rec["higher_expression_in"]
    if higher == "CLDN4 High":
        return 1
    if higher == "CLDN4 Low":
        return -1
    ratio = rec["log2_ratio_high_over_low"]
    if ratio > 0:
        return 1
    if ratio < 0:
        return -1
    return 0


def load_gse207704(path: Path) -> dict[str, dict[str, float]]:
    """Sum cufflinks FPKM by gene. Columns are condition-level, not replicates."""
    sums: dict[str, dict[str, float]] = {}
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        cols = {
            "MCF7_KO": idx["MCF7_CLDN4KO_FPKM (fpkm)"],
            "MCF7_WT": idx["MCF7_WT_FPKM (fpkm)"],
            "T47D_KO": idx["T47D_CLDN4KO_FPKM (fpkm)"],
            "T47D_WT": idx["T47D_WT_FPKM (fpkm)"],
        }
        name_i = idx["gene_short_name"]
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= name_i:
                continue
            gene = canon(parts[name_i])
            if not gene or gene in {"-", "NA"}:
                continue
            bucket = sums.setdefault(gene, {k: 0.0 for k in cols})
            for key, col in cols.items():
                try:
                    val = float(parts[col])
                except ValueError:
                    val = 0.0
                if math.isfinite(val):
                    bucket[key] += val
    return sums


def log2fc_fpkm(ko: float, wt: float) -> float | None:
    if max(ko, wt) < FPKM_FLOOR:
        return None
    return math.log2((ko + FPKM_PSEUDO) / (wt + FPKM_PSEUDO))


def symbol_from_operon_name(name: str) -> str | None:
    name = name.strip().strip('"')
    pref = name.split("--", 1)[0].strip() if "--" in name else name.strip()
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9\-]{1,20}", pref) and " " not in pref and ";" not in pref:
        return canon(pref)
    return None


def load_operon_symbols(path: Path) -> dict[str, str]:
    mapping = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        header = None
        for line in fh:
            if line.startswith("Index\t"):
                header = line.rstrip("\n").split("\t")
                continue
            if header is None:
                continue
            if line.startswith("END DATA"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            rec = dict(zip(header, parts))
            sym = symbol_from_operon_name(rec.get("Name", ""))
            if sym:
                mapping[rec["Index"]] = sym
    return mapping


def load_gse22493_matrix(path: Path) -> dict[str, list[float]]:
    mat = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        started = False
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                started = True
                continue
            if not started:
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if line.startswith('"ID_REF"') or line.startswith("ID_REF"):
                continue
            parts = line.rstrip("\n").split("\t")
            idx = parts[0].strip().strip('"')
            vals = []
            for item in parts[1:4]:
                item = item.strip().strip('"')
                if item == "" or item.lower() == "null":
                    vals.append(math.nan)
                else:
                    try:
                        vals.append(float(item))
                    except ValueError:
                        vals.append(math.nan)
            while len(vals) < 3:
                vals.append(math.nan)
            mat[idx] = vals[:3]
    return mat


def collapse_probes(mat: dict[str, list[float]], symbols: dict[str, str]) -> dict[str, dict]:
    """Median of probes within each replicate, then mean across replicates."""
    by_gene: dict[str, list[list[float]]] = {}
    for idx, vals in mat.items():
        gene = symbols.get(idx)
        if gene is None:
            continue
        by_gene.setdefault(gene, []).append(vals)
    out = {}
    for gene, probes in by_gene.items():
        arr = np.array(probes, dtype=float)
        if not np.isfinite(arr).any():
            continue
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            rep = np.nanmedian(arr, axis=0)
        finite = rep[np.isfinite(rep)]
        if finite.size < 2:
            continue
        out[gene] = {
            "log2fc": float(np.mean(finite)),
            "n_reps": int(finite.size),
            "n_probes": int(arr.shape[0]),
            "rep_log2fc": [None if not math.isfinite(x) else float(x) for x in rep],
        }
    return out


def entrez_count(db: str, term: str) -> int:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db={db}&retmode=json&retmax=0&term=" + urllib.request.quote(term)
    )
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.loads(resp.read().decode())
    return int(payload["esearchresult"]["count"])


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def concordance_stat(signs: np.ndarray, logfc: np.ndarray) -> float:
    """Mean of sign_with_CLDN4 * (-log2FC). Positive means the KD moves genes
    the way the reported CLDN4 association predicts."""
    return float(np.mean(signs * (-logfc)))


def perm_p(obs: float, null: np.ndarray) -> float:
    return float((np.sum(np.abs(null) >= abs(obs)) + 1) / (null.size + 1))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)

    table = load_table_s1(TABLE_S1)
    n_sig = sum(1 for r in table.values() if r["q_value"] < Q_CUTOFF)

    fpkm_path = CACHE / "GSE207704_CLDN4_RNAseq.txt.gz"
    matrix_path = CACHE / "GSE22493_series_matrix.txt.gz"
    gsm_path = CACHE / "GSM558700.txt.gz"
    download(GSE207704_URL, fpkm_path)
    download(GSE22493_URL, matrix_path)
    download(GSM558700_URL, gsm_path)

    fpkm = load_gse207704(fpkm_path)
    symbols = load_operon_symbols(gsm_path)
    matrix = load_gse22493_matrix(matrix_path)
    skov = collapse_probes(matrix, symbols)

    contrasts = {}
    # GSE207704 condition-level FPKM
    for label, ko_key, wt_key in [
        ("GSE207704_T47D_CLDN4ko", "T47D_KO", "T47D_WT"),
        ("GSE207704_MCF7_CLDN4ko", "MCF7_KO", "MCF7_WT"),
    ]:
        logfc = {}
        for gene, vals in fpkm.items():
            d = log2fc_fpkm(vals[ko_key], vals[wt_key])
            if d is not None:
                logfc[gene] = d
        contrasts[label] = logfc

    # Equal-weight mean of the two cell-line log2FC values, genes measured in both.
    both = set(contrasts["GSE207704_T47D_CLDN4ko"]) & set(contrasts["GSE207704_MCF7_CLDN4ko"])
    contrasts["GSE207704_mean_T47D_MCF7"] = {
        g: 0.5
        * (
            contrasts["GSE207704_T47D_CLDN4ko"][g]
            + contrasts["GSE207704_MCF7_CLDN4ko"][g]
        )
        for g in both
    }
    contrasts["GSE22493_SKOV3_CLDN4kd"] = {g: rec["log2fc"] for g, rec in skov.items()}

    # Signed NHEJ set: membership AND Table S1 q<0.05. XRCC1 stays in the
    # membership table but is excluded from the NHEJ score (SSBR/BER).
    signed_rows = []
    for gene, source in sorted(NHEJ_MEMBERSHIP.items()):
        rec = table.get(gene)
        row = {
            "gene": gene,
            "membership": source,
            "in_table_s1": int(rec is not None),
            "log2_ratio_high_over_low": "" if rec is None else f"{rec['log2_ratio_high_over_low']:.6g}",
            "p_value": "" if rec is None else f"{rec['p_value']:.6g}",
            "q_value": "" if rec is None else f"{rec['q_value']:.6g}",
            "higher_expression_in": "" if rec is None else rec["higher_expression_in"],
            "q_lt_0.05": int(rec is not None and rec["q_value"] < Q_CUTOFF),
            "sign_with_cldn4": "" if rec is None else sign_from_table(rec),
            "rppa_protein_down_in_cldn4_kd": int(gene in RPPA_DOWN),
            "in_signed_nhej_set": int(
                source != "paper_RPPA_SSBR_not_classical_NHEJ"
                and rec is not None
                and rec["q_value"] < Q_CUTOFF
            ),
        }
        signed_rows.append(row)
    write_tsv(
        OUT / "nhej_signed_geneset.tsv",
        signed_rows,
        list(signed_rows[0].keys()),
    )

    signed_genes = [r["gene"] for r in signed_rows if r["in_signed_nhej_set"] == 1]
    signed_sign = {
        r["gene"]: int(r["sign_with_cldn4"])
        for r in signed_rows
        if r["in_signed_nhej_set"] == 1
    }

    # Background for permutation: Table S1 q<0.05 genes, excluding XRCC1 so the
    # null matches the NHEJ filter's transcriptome definition.
    sig_genes = {
        g
        for g, rec in table.items()
        if rec["q_value"] < Q_CUTOFF and g != "XRCC1"
    }

    rng = np.random.default_rng(RNG_SEED)
    summary_rows = []
    gene_score_rows = []

    for contrast, logfc_map in contrasts.items():
        measured_sig = sorted(g for g in sig_genes if g in logfc_map)
        set_genes = [g for g in signed_genes if g in logfc_map]
        missing = [g for g in signed_genes if g not in logfc_map]
        signs = np.array([signed_sign[g] for g in set_genes], dtype=float)
        vals = np.array([logfc_map[g] for g in set_genes], dtype=float)
        obs = concordance_stat(signs, vals) if set_genes else float("nan")
        n_down = int(np.sum(vals < 0)) if set_genes else 0
        n_concordant = int(np.sum(np.sign(-vals) == signs)) if set_genes else 0
        mean_logfc = float(np.mean(vals)) if set_genes else float("nan")

        bg_genes = measured_sig
        bg_signs = np.array([sign_from_table(table[g]) for g in bg_genes], dtype=float)
        bg_vals = np.array([logfc_map[g] for g in bg_genes], dtype=float)
        k = len(set_genes)
        null = np.empty(N_PERM, dtype=float)
        null_mean = np.empty(N_PERM, dtype=float)
        if k >= 2 and len(bg_genes) > k:
            for i in range(N_PERM):
                draw = rng.choice(len(bg_genes), size=k, replace=False)
                null[i] = concordance_stat(bg_signs[draw], bg_vals[draw])
                null_mean[i] = float(np.mean(bg_vals[draw]))
            p_conc = perm_p(obs, null)
            p_mean = perm_p(mean_logfc, null_mean)
        else:
            null[:] = np.nan
            null_mean[:] = np.nan
            p_conc = float("nan")
            p_mean = float("nan")

        # Unsigned core NHEJ panel (all detected classical/53BP1-axis genes,
        # not only q<0.05). Tests the directional claim NHEJ transcript down.
        panel = [
            g
            for g, src in NHEJ_MEMBERSHIP.items()
            if src != "paper_RPPA_SSBR_not_classical_NHEJ" and g in logfc_map
        ]
        panel_vals = np.array([logfc_map[g] for g in panel], dtype=float)
        panel_mean = float(np.mean(panel_vals)) if panel else float("nan")
        # Permute among all scored genes in the contrast.
        universe = np.array([logfc_map[g] for g in sorted(logfc_map)], dtype=float)
        null_panel = np.empty(N_PERM, dtype=float)
        if len(panel) >= 2 and universe.size > len(panel):
            for i in range(N_PERM):
                draw = rng.choice(universe.size, size=len(panel), replace=False)
                null_panel[i] = float(np.mean(universe[draw]))
            p_panel = perm_p(panel_mean, null_panel)
        else:
            p_panel = float("nan")

        summary_rows.append(
            {
                "contrast": contrast,
                "n_signed_set_measured": k,
                "n_signed_set_missing": len(missing),
                "signed_set_missing": ",".join(missing),
                "n_sign_plus": int(np.sum(signs > 0)) if set_genes else 0,
                "n_sign_minus": int(np.sum(signs < 0)) if set_genes else 0,
                "n_concordant": n_concordant,
                "concordant_fraction": (n_concordant / k) if k else "",
                "mean_signed_concordance": obs,
                "perm_p_concordance": p_conc,
                "mean_log2fc_signed_set": mean_logfc,
                "n_log2fc_negative": n_down,
                "perm_p_mean_log2fc": p_mean,
                "n_background_q05": len(bg_genes),
                "n_perm": N_PERM,
                "n_core_nhej_detected": len(panel),
                "mean_log2fc_core_nhej": panel_mean,
                "frac_core_nhej_down": (float(np.mean(panel_vals < 0)) if panel else ""),
                "perm_p_core_nhej_mean": p_panel,
                "CLDN4_log2fc": logfc_map.get("CLDN4", ""),
                "TP53BP1_log2fc": logfc_map.get("TP53BP1", ""),
                "XRCC1_log2fc": logfc_map.get("XRCC1", ""),
                "STING1_log2fc": logfc_map.get("STING1", ""),
            }
        )

        for gene in sorted(set(signed_genes) | RPPA_DOWN | {"CLDN4", "STING1"} | set(PAPER_NAMED_REPAIR)):
            rec = table.get(gene)
            gene_score_rows.append(
                {
                    "contrast": contrast,
                    "gene": gene,
                    "in_signed_nhej_set": int(gene in signed_sign),
                    "sign_with_cldn4": signed_sign.get(gene, ""),
                    "table_s1_log2_ratio": "" if rec is None else rec["log2_ratio_high_over_low"],
                    "table_s1_q": "" if rec is None else rec["q_value"],
                    "log2fc_kd_vs_control": logfc_map.get(gene, ""),
                    "measured": int(gene in logfc_map),
                }
            )

    # GSE22493 replicate-level t-tests for the effector pair and CLDN4.
    ttest_rows = []
    seen = []
    for gene in ["CLDN4", "TP53BP1", "XRCC1", "STING1", "PRKDC", *signed_genes]:
        if gene in seen:
            continue
        seen.append(gene)
        rec = skov.get(gene)
        if rec is None:
            ttest_rows.append(
                {
                    "gene": gene,
                    "n_reps": 0,
                    "mean_log2fc": "",
                    "t_stat": "",
                    "ttest_p": "",
                    "rep_log2fc": "",
                }
            )
            continue
        reps = np.array([x for x in rec["rep_log2fc"] if x is not None], dtype=float)
        if reps.size >= 2 and np.std(reps, ddof=1) > 0:
            t_stat, p = stats.ttest_1samp(reps, 0.0)
        else:
            t_stat, p = float("nan"), float("nan")
        ttest_rows.append(
            {
                "gene": gene,
                "n_reps": int(reps.size),
                "n_probes": rec["n_probes"],
                "mean_log2fc": float(np.mean(reps)),
                "t_stat": t_stat,
                "ttest_p": p,
                "rep_log2fc": ",".join(f"{x:.4g}" for x in reps),
            }
        )

    write_tsv(OUT / "set_scores.tsv", summary_rows, list(summary_rows[0].keys()))
    write_tsv(OUT / "gene_scores.tsv", gene_score_rows, list(gene_score_rows[0].keys()))
    write_tsv(OUT / "gse22493_gene_ttests.tsv", ttest_rows, list(ttest_rows[0].keys()))

    # Data-availability queries. Failures are recorded, not invented.
    queries = [
        ("gds", "Yamamoto[All Fields] AND Bitler[All Fields] AND CLDN4[All Fields]"),
        ("gds", "35373300[PMID]"),
        ("gds", "MCT-21-0827"),
        ("sra", "Yamamoto[All Fields] AND Bitler[All Fields] AND CLDN4[All Fields]"),
        ("bioproject", "Yamamoto[All Fields] AND Bitler[All Fields] AND CLDN4[All Fields]"),
        ("gds", "GSE207704[Accession]"),
        ("gds", "GSE22493[Accession]"),
    ]
    avail_rows = []
    for db, term in queries:
        try:
            count = entrez_count(db, term)
            note = ""
        except Exception as exc:  # noqa: BLE001
            count = ""
            note = f"query_failed: {exc}"
        avail_rows.append({"database": db, "term": term, "count": count, "note": note})
    write_tsv(OUT / "data_availability_queries.tsv", avail_rows, ["database", "term", "count", "note"])

    # Named-gene snapshot from Table S1.
    named_rows = []
    named_seen = set()
    for gene in PAPER_NAMED_REPAIR + ["STING1", "TMEM173", "CLDN4"]:
        gene = canon(gene)
        if gene in named_seen:
            continue
        named_seen.add(gene)
        rec = table.get(gene)
        if rec is None:
            named_rows.append({"gene": gene, "in_table_s1": 0})
            continue
        named_rows.append(
            {
                "gene": rec["gene"],
                "in_table_s1": 1,
                "log2_ratio_high_over_low": rec["log2_ratio_high_over_low"],
                "p_value": rec["p_value"],
                "q_value": rec["q_value"],
                "higher_expression_in": rec["higher_expression_in"],
                "q_lt_0.05": int(rec["q_value"] < Q_CUTOFF),
            }
        )
    write_tsv(OUT / "paper_named_genes_table_s1.tsv", named_rows, list(named_rows[0].keys()))

    meta = {
        "table_s1_genes": len(table),
        "table_s1_q_lt_0.05": n_sig,
        "signed_nhej_genes": signed_genes,
        "gse207704_genes": len(fpkm),
        "gse22493_genes_ge2reps": len(skov),
        "gse22493_probes": len(matrix),
        "gse22493_probes_with_symbol": len(symbols),
        "fpkm_floor": FPKM_FLOOR,
        "fpkm_pseudo": FPKM_PSEUDO,
        "n_perm": N_PERM,
        "seed": RNG_SEED,
    }
    (OUT / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(json.dumps(meta, indent=2))
    print("--- set scores ---")
    for row in summary_rows:
        print(
            row["contrast"],
            "k",
            row["n_signed_set_measured"],
            "conc",
            row["mean_signed_concordance"],
            "p",
            row["perm_p_concordance"],
            "meanLFC",
            row["mean_log2fc_signed_set"],
            "panel",
            row["mean_log2fc_core_nhej"],
            "p_panel",
            row["perm_p_core_nhej_mean"],
            "CLDN4",
            row["CLDN4_log2fc"],
            "53BP1",
            row["TP53BP1_log2fc"],
            "XRCC1",
            row["XRCC1_log2fc"],
        )


if __name__ == "__main__":
    main()
