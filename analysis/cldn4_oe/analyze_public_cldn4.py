#!/usr/bin/env python3
"""Score NHEJ and interferon on open CLDN4-high versus CLDN4-low GEO profiles.

No public GEO/SRA series is a CLDN4 cDNA overexpression RNA-seq experiment.
The open contrasts scored here are:

* GSE207704: parental versus CLDN4-knockout RNA-seq in MCF7 and T47D
  (Murakami-Chiba et al., Breast Cancer Research 2023). Replicates are
  quantified from the deposited SRA runs with Salmon. The GEO FPKM table
  is also scored, and it is replicate-collapsed and missing many genes.
* GSE22493: two-color microarray of parental SKOV3 cells, labeled
  "CLDN4 overexpression (control)", versus CLDN4 siRNA (Gao, 2010).
  This is a knockdown array, not an overexpression RNA-seq.

The CLDN4-high direction is parental/control over CLDN4-low. That is the
direction in which an overexpression experiment would be read: NHEJ up
and interferon down, if that pattern holds.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent
GENESETS = ROOT / "genesets"
TABLES = ROOT / "tables"

HISTONE_RE = re.compile(r"^H2[AB]C|^H3C|^H3-|^H4C")
LIGAND_RE = re.compile(r"^IFNA\d+$|^IFNB1$|^IFNK$|^IFNW1$")

PRIMARY = (
    "KEGG_NHEJ",
    "Hallmark_IFN_alpha",
    "Hallmark_IFN_gamma",
)


def load_list(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def load_sets() -> dict[str, list[str]]:
    kegg = load_list(GENESETS / "Non-homologous_end-joining.txt")
    ifna = load_list(GENESETS / "Interferon_Alpha_Response.txt")
    ifng = load_list(GENESETS / "Interferon_Gamma_Response.txt")
    react_nhej = load_list(GENESETS / "Nonhomologous_End-Joining_(NHEJ).txt")
    react_ifn = load_list(GENESETS / "Interferon_Alpha_Beta_Signaling.txt")
    return {
        "KEGG_NHEJ": kegg,
        "Hallmark_IFN_alpha": ifna,
        "Hallmark_IFN_gamma": ifng,
        "Reactome_NHEJ_no_histone": [g for g in react_nhej if not HISTONE_RE.match(g)],
        "Reactome_IFN_ab_no_ligand": [g for g in react_ifn if not LIGAND_RE.match(g)],
    }


def load_alias_groups() -> dict[str, set[str]]:
    """Map every symbol and synonym in the frozen table onto its alias group."""
    groups: dict[str, set[str]] = {}
    text = (GENESETS / "hgnc_aliases.tsv").read_text().splitlines()[1:]
    for line in text:
        query, official, syn = line.split("\t")
        names = {query, official}
        if syn:
            names.update(s for s in syn.split(",") if s)
        bucket: set[str] = set()
        for name in names:
            bucket |= groups.get(name, set())
        bucket |= names
        for name in bucket:
            groups[name] = bucket
    return groups


def resolve(symbol: str, available: set[str], groups: dict[str, set[str]]) -> str | None:
    if symbol in available:
        return symbol
    # Set order is hash-randomized; sort so a synonym collision is stable.
    for alt in sorted(groups.get(symbol, ())):
        if alt in available:
            return alt
    return None


def bh(pvals: list[float]) -> list[float]:
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty(n, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        val = ranked[i] * n / (i + 1) if np.isfinite(ranked[i]) else 1.0
        prev = min(prev, val)
        q[i] = prev
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out.tolist()


def mannwhitney(set_fc: np.ndarray, bg_fc: np.ndarray) -> tuple[float, float]:
    if len(set_fc) < 5 or len(bg_fc) < 50:
        return float("nan"), float("nan")
    _u, p = stats.mannwhitneyu(set_fc, bg_fc, alternative="two-sided")
    return float(_u), float(p)


def write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")


def score_logfc(logfc: dict[str, float], sets: dict[str, list[str]], groups, label: str) -> tuple[list[dict], list[dict]]:
    available = set(logfc)
    set_rows = []
    gene_rows = []
    bg_genes = list(logfc)
    for name, genes in sets.items():
        mapped = []
        for g in genes:
            hit = resolve(g, available, groups)
            if hit is not None and hit not in {m[1] for m in mapped}:
                mapped.append((g, hit))
        values = np.array([logfc[hit] for _, hit in mapped], dtype=float)
        bg = np.array([logfc[g] for g in bg_genes if g not in {hit for _, hit in mapped}], dtype=float)
        _u, p = mannwhitney(values, bg)
        set_rows.append(
            {
                "contrast": label,
                "set": name,
                "role": "primary" if name in PRIMARY else "sensitivity",
                "n_set_genes": len(genes),
                "n_detected": int(len(values)),
                "mean_log2fc": f"{float(np.mean(values)):.4f}" if len(values) else "",
                "median_log2fc": f"{float(np.median(values)):.4f}" if len(values) else "",
                "frac_positive": f"{float(np.mean(values > 0)):.3f}" if len(values) else "",
                "rank_p": f"{p:.6g}" if np.isfinite(p) else "",
                "q_bh_primary": "",
            }
        )
        for query, hit in mapped:
            gene_rows.append(
                {
                    "contrast": label,
                    "set": name,
                    "query_symbol": query,
                    "matrix_symbol": hit,
                    "log2fc": f"{logfc[hit]:.4f}",
                }
            )
    primary_idx = [i for i, r in enumerate(set_rows) if r["role"] == "primary" and r["rank_p"] != ""]
    qvals = bh([float(set_rows[i]["rank_p"]) for i in primary_idx])
    for i, q in zip(primary_idx, qvals):
        set_rows[i]["q_bh_primary"] = f"{q:.6g}"
    return set_rows, gene_rows


def load_fpkm(path: Path) -> dict[str, dict[str, float]]:
    best: dict[str, tuple[float, dict[str, float]]] = {}
    cols = {
        "MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)",
        "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
        "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)",
        "T47D_WT": "T47D_WT_FPKM (fpkm)",
    }
    with gzip.open(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: header.index(name) for name in ["gene_short_name", *cols.values()]}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            symbol = parts[idx["gene_short_name"]].strip()
            if not symbol or symbol == "-":
                continue
            rec = {}
            ok = True
            for key, col in cols.items():
                try:
                    rec[key] = float(parts[idx[col]])
                except ValueError:
                    ok = False
            if not ok:
                continue
            mean_v = sum(rec.values()) / 4
            if symbol not in best or mean_v > best[symbol][0]:
                best[symbol] = (mean_v, rec)
    return {symbol: rec for symbol, (_mean, rec) in best.items()}


def fpkm_logfc(expr: dict[str, dict[str, float]], high: str, low: str, min_fpkm: float = 1.0) -> dict[str, float]:
    out = {}
    for symbol, rec in expr.items():
        if max(rec[high], rec[low]) < min_fpkm:
            continue
        out[symbol] = math.log2((rec[high] + 0.1) / (rec[low] + 0.1))
    return out


def load_tx2gene(fasta: Path) -> dict[str, str]:
    tx2gene = {}
    opener = gzip.open if str(fasta).endswith(".gz") else open
    with opener(fasta, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            token = header.split()[0]
            symbol = ""
            for field in header.split():
                if field.startswith("gene_symbol:"):
                    symbol = field.split(":", 1)[1]
            if symbol:
                tx2gene[token] = symbol
                tx2gene[token.split(".")[0]] = symbol
    return tx2gene


def load_salmon_counts(quant_dir: Path, tx2gene: dict[str, str]) -> dict[str, dict[str, float]]:
    samples = sorted(p.parent.name for p in quant_dir.glob("*/quant.sf"))
    counts: dict[str, dict[str, float]] = {}
    for sample in samples:
        path = quant_dir / sample / "quant.sf"
        with path.open() as fh:
            header = fh.readline().rstrip("\n").split("\t")
            name_i = header.index("Name")
            read_i = header.index("NumReads")
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                name = parts[name_i].split("|")[0]
                symbol = tx2gene.get(name) or tx2gene.get(name.split(".")[0])
                if not symbol:
                    continue
                rec = counts.setdefault(symbol, {s: 0.0 for s in samples})
                rec[sample] += float(parts[read_i])
    # genes missing from a sample stay 0 because setdefault fills every sample
    return counts


def counts_to_cpm(counts: dict[str, dict[str, float]]) -> dict[str, dict[str, float]]:
    samples = next(iter(counts.values())).keys()
    lib = {sample: sum(rec[sample] for rec in counts.values()) for sample in samples}
    cpm = {}
    for symbol, rec in counts.items():
        cpm[symbol] = {sample: (rec[sample] / lib[sample]) * 1e6 if lib[sample] else 0.0 for sample in samples}
    return cpm


def mean_cpm_logfc(cpm: dict[str, dict[str, float]], high_samples: list[str], low_samples: list[str], min_cpm: float = 1.0) -> dict[str, float]:
    out = {}
    for symbol, rec in cpm.items():
        high = float(np.mean([rec[s] for s in high_samples]))
        low = float(np.mean([rec[s] for s in low_samples]))
        if max(high, low) < min_cpm:
            continue
        out[symbol] = math.log2((high + 0.1) / (low + 0.1))
    return out


def replicate_set_deltas(cpm, sets, groups, high_samples, low_samples, min_cpm=1.0) -> list[dict]:
    """Mean log2(CPM+0.1) of each set in each sample, then WT minus KO."""
    detected = {g for g, rec in cpm.items() if max(rec[s] for s in high_samples + low_samples) >= min_cpm}
    rows = []
    for name, genes in sets.items():
        mapped = []
        for g in genes:
            hit = resolve(g, detected, groups)
            if hit is not None and hit not in mapped:
                mapped.append(hit)
        if len(mapped) < 5:
            continue
        sample_means = {}
        for sample in high_samples + low_samples:
            vals = [math.log2(cpm[g][sample] + 0.1) for g in mapped]
            sample_means[sample] = float(np.mean(vals))
        high_mean = float(np.mean([sample_means[s] for s in high_samples]))
        low_mean = float(np.mean([sample_means[s] for s in low_samples]))
        row = {
            "set": name,
            "role": "primary" if name in PRIMARY else "sensitivity",
            "n_detected": len(mapped),
            "mean_log2cpm_high": f"{high_mean:.4f}",
            "mean_log2cpm_low": f"{low_mean:.4f}",
            "delta_high_minus_low": f"{(high_mean - low_mean):.4f}",
        }
        for sample, value in sample_means.items():
            row[sample] = f"{value:.4f}"
        rows.append(row)
    return rows


def load_platform(path: Path) -> dict[str, str]:
    mapping = {}
    with gzip.open(path, "rt", errors="replace") as fh:
        in_table = False
        for line in fh:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table or line.startswith("ID\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[2].strip():
                mapping[parts[0]] = parts[2].strip()
    return mapping


def load_array_ratios(matrix_path: Path, probe_to_gene: dict[str, str]) -> dict[str, np.ndarray]:
    """Deposited series-matrix VALUE is log2(knockdown / parental control)."""
    buckets: dict[str, list[list[float]]] = {}
    with gzip.open(matrix_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"):
                next(fh)
                break
        for line in fh:
            if line.startswith("!series_matrix_table_end"):
                break
            parts = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]
            if len(parts) < 4:
                continue
            gene = probe_to_gene.get(parts[0])
            if not gene:
                continue
            slot = buckets.setdefault(gene, [[], [], []])
            for i, raw in enumerate(parts[1:4]):
                if raw and raw.upper() != "NA":
                    slot[i].append(float(raw))
    out = {}
    for gene, arrs in buckets.items():
        vals = []
        for xs in arrs:
            vals.append(float(np.mean(xs)) if xs else float("nan"))
        if sum(not math.isnan(v) for v in vals) >= 2:
            out[gene] = np.array(vals, dtype=float)
    return out


def array_set_tests(ratios: dict[str, np.ndarray], sets, groups, sign: float, label: str) -> tuple[list[dict], list[dict]]:
    gene_mean = {}
    for gene, arr in ratios.items():
        m = float(np.nanmean(arr))
        if not math.isnan(m):
            gene_mean[gene] = sign * m
    set_rows, gene_rows = score_logfc(gene_mean, sets, groups, label)
    # Replace rank p as the headline for this dataset: the replicate unit is the array.
    available = set(gene_mean)
    by_set = {r["set"]: r for r in set_rows}
    for name, genes in sets.items():
        mapped = []
        for g in genes:
            hit = resolve(g, available, groups)
            if hit is not None and hit not in mapped:
                mapped.append(hit)
        scores = []
        for i in range(3):
            xs = []
            for gene in mapped:
                v = ratios[gene][i]
                if not math.isnan(v):
                    xs.append(sign * v)
            if xs:
                scores.append(float(np.mean(xs)))
        scores_a = np.array(scores, dtype=float)
        if len(scores_a) == 3 and np.std(scores_a) > 0:
            _t, p = stats.ttest_1samp(scores_a, 0.0)
        else:
            p = float("nan")
        row = by_set[name]
        row["array_means"] = ",".join(f"{v:.4f}" for v in scores_a)
        row["array_mean"] = f"{float(np.mean(scores_a)):.4f}" if len(scores_a) else ""
        row["array_t_p"] = f"{float(p):.6g}" if np.isfinite(p) else ""
        row["array_q_bh_primary"] = ""
    primary = [r for r in set_rows if r["role"] == "primary" and r.get("array_t_p")]
    qvals = bh([float(r["array_t_p"]) for r in primary])
    for row, q in zip(primary, qvals):
        row["array_q_bh_primary"] = f"{q:.6g}"
    return set_rows, gene_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geo", type=Path, default=Path("/tmp/geo"))
    args = parser.parse_args()
    geo = args.geo
    sets = load_sets()
    groups = load_alias_groups()
    TABLES.mkdir(parents=True, exist_ok=True)

    set_rows: list[dict] = []
    gene_rows: list[dict] = []
    rep_rows: list[dict] = []

    fpkm_path = geo / "GSE207704_CLDN4_RNAseq.txt.gz"
    if fpkm_path.exists():
        expr = load_fpkm(fpkm_path)
        cldn4 = expr.get("CLDN4", {})
        (TABLES / "gse207704_cldn4_fpkm.tsv").write_text(
            "sample\tfpkm\n" + "".join(f"{k}\t{v}\n" for k, v in cldn4.items())
        )
        for line, high, low in (("MCF7", "MCF7_WT", "MCF7_KO"), ("T47D", "T47D_WT", "T47D_KO")):
            logfc = fpkm_logfc(expr, high, low)
            srows, grows = score_logfc(logfc, sets, groups, f"GSE207704_deposited_FPKM_{line}_WT_over_KO")
            set_rows.extend(srows)
            gene_rows.extend(grows)

    quant_dir = geo / "quant"
    fasta = geo / "ref" / "pc.fa"
    if quant_dir.exists() and any(quant_dir.glob("*/quant.sf")) and fasta.exists():
        tx2gene = load_tx2gene(fasta)
        counts = load_salmon_counts(quant_dir, tx2gene)
        cpm = counts_to_cpm(counts)
        map_rows = []
        for meta_path in sorted(quant_dir.glob("*/aux_info/meta_info.json")):
            meta = json.loads(meta_path.read_text())
            map_rows.append(
                {
                    "sample": meta_path.parent.parent.name,
                    "percent_mapped": f"{float(meta.get('percent_mapped', float('nan'))):.2f}",
                    "num_mapped": meta.get("num_mapped", ""),
                    "num_processed": meta.get("num_processed", ""),
                }
            )
        if map_rows:
            write_tsv(
                TABLES / "gse207704_salmon_mapping.tsv",
                map_rows,
                ["sample", "percent_mapped", "num_mapped", "num_processed"],
            )
        cldn4_cpm = cpm.get("CLDN4", {})
        with (TABLES / "gse207704_cldn4_salmon_cpm.tsv").open("w") as fh:
            fh.write("sample\tcpm\tnumreads\n")
            for sample, value in cldn4_cpm.items():
                fh.write(f"{sample}\t{value:.4f}\t{counts['CLDN4'][sample]:.2f}\n")
        designs = {
            "MCF7": (["MCF7_WT_rep1", "MCF7_WT_rep2"], ["MCF7_KO_rep1", "MCF7_KO_rep2"]),
            "T47D": (["T47D_WT_rep1", "T47D_WT_rep2"], ["T47D_KO_rep1", "T47D_KO_rep2"]),
        }
        for line, (high_s, low_s) in designs.items():
            if any(s not in next(iter(cpm.values())) for s in high_s + low_s):
                continue
            logfc = mean_cpm_logfc(cpm, high_s, low_s)
            srows, grows = score_logfc(logfc, sets, groups, f"GSE207704_salmon_{line}_WT_over_KO")
            set_rows.extend(srows)
            gene_rows.extend(grows)
            for row in replicate_set_deltas(cpm, sets, groups, high_s, low_s):
                row["cell_line"] = line
                rep_rows.append(row)

    platform = geo / "GPL10555_family.soft.gz"
    matrix = geo / "GSE22493_series_matrix.txt.gz"
    array_rows: list[dict] = []
    if platform.exists() and matrix.exists():
        ratios = load_array_ratios(matrix, load_platform(platform))
        cldn4 = ratios.get("CLDN4")
        if cldn4 is not None:
            (TABLES / "gse22493_cldn4_log2_kd_over_control.tsv").write_text(
                "array\tGSM\tlog2_KD_over_control\n"
                + "".join(
                    f"{i+1}\t{gsm}\t{'' if math.isnan(v) else f'{v:.4f}'}\n"
                    for i, (gsm, v) in enumerate(zip(("GSM558700", "GSM558701", "GSM558702"), cldn4))
                )
            )
        srows, grows = array_set_tests(
            ratios, sets, groups, sign=-1.0, label="GSE22493_parental_over_siRNA"
        )
        array_rows.extend(srows)
        gene_rows.extend(grows)

    set_cols = [
        "contrast", "set", "role", "n_set_genes", "n_detected",
        "mean_log2fc", "median_log2fc", "frac_positive", "rank_p", "q_bh_primary",
    ]
    write_tsv(TABLES / "set_scores.tsv", set_rows, set_cols)
    array_cols = set_cols + ["array_means", "array_mean", "array_t_p", "array_q_bh_primary"]
    write_tsv(TABLES / "gse22493_set_scores.tsv", array_rows, array_cols)
    write_tsv(
        TABLES / "gene_log2fc.tsv",
        gene_rows,
        ["contrast", "set", "query_symbol", "matrix_symbol", "log2fc"],
    )
    if rep_rows:
        sample_cols = sorted({k for row in rep_rows for k in row if k.endswith("_rep1") or k.endswith("_rep2")})
        write_tsv(
            TABLES / "gse207704_salmon_replicate_set_means.tsv",
            rep_rows,
            ["cell_line", "set", "role", "n_detected", "mean_log2cpm_high", "mean_log2cpm_low", "delta_high_minus_low", *sample_cols],
        )
    print(f"wrote tables under {TABLES}")
    for row in set_rows:
        if row["role"] == "primary":
            print(row["contrast"], row["set"], "n", row["n_detected"], "mean", row["mean_log2fc"], "p", row["rank_p"], "q", row["q_bh_primary"])
    for row in array_rows:
        if row["role"] == "primary":
            print(row["contrast"], row["set"], "arrays", row.get("array_means"), "t_p", row.get("array_t_p"), "q", row.get("array_q_bh_primary"))


if __name__ == "__main__":
    main()
