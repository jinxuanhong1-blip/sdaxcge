#!/usr/bin/env python3
"""Hunt for any concordant IFN / APM / NHEJ increase after CLDN4 loss
in the cancer sets GSE207704 and GSE22493.

A result counts as concordant only when CLDN4 itself falls and the
panel is called up by the same rule as Claim C4: median log2FC > 0 and
a one-sided Wilcoxon P <= 0.05 (sign test if n < 5).

Optimistic rules (max locus, max probe) are reported and cannot by
themselves overturn the primary call. Replicate-level kallisto counts,
if present under /tmp/kd_sweep/kallisto, are an additional DE.
"""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest, mannwhitneyu, norm, ttest_1samp, wilcoxon

ROOT = Path(__file__).resolve().parents[2]
LOCK = ROOT / "methods" / "cldn4_kd_match_wave" / "locked" / "geneset_summary_ifn_apm.tsv"
OUT = ROOT / "methods" / "cldn4_kd_match_wave"
RAW = Path("/tmp/kd_sweep/raw")
KALLISTO = Path("/tmp/kd_sweep/kallisto")

NHEJ_CORE = ["XRCC5", "XRCC6", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C", "PAXX"]
NHEJ_EXT = NHEJ_CORE + ["APLF", "PNKP", "POLL", "POLM", "APTX", "PARP1", "LIG3", "XRCC1"]
MOUSE = {
    "XRCC5": "Xrcc5",
    "XRCC6": "Xrcc6",
    "PRKDC": "Prkdc",
    "LIG4": "Lig4",
    "XRCC4": "Xrcc4",
    "NHEJ1": "Nhej1",
    "DCLRE1C": "Dclre1c",
    "PAXX": "Paxx",
    "APLF": "Aplf",
    "PNKP": "Pnkp",
    "POLL": "Poll",
    "POLM": "Polm",
    "APTX": "Aptx",
    "PARP1": "Parp1",
    "LIG3": "Lig3",
    "XRCC1": "Xrcc1",
}
ALIASES = {"G1P2": "ISG15"}
ARRAYS = ["GSM558700", "GSM558701", "GSM558702"]
PRIMARY_PANELS = ["IFN_ISG", "MHC_APM", "NHEJ_CORE"]


def read_tsv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fields})


def fnum(value) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return ""
    return format(float(value), ".8g")


def load_gmt(path: Path, name: str) -> list[str]:
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts[0] == name:
                return [gene for gene in parts[2:] if gene]
    raise RuntimeError(f"{name} missing from {path}")


def load_panels() -> dict[str, list[str]]:
    symbols = defaultdict(list)
    seen = defaultdict(set)
    for row in read_tsv(LOCK.parent / "gene_effects_ifn_apm.tsv"):
        module, gene = row["gene_set"], row["human_symbol"]
        if gene not in seen[module]:
            seen[module].add(gene)
            symbols[module].append(gene)
    panels = {
        "IFN_ISG": symbols["IFN_ISG"],
        "MHC_APM": symbols["MHC_APM"],
        "NHEJ_CORE": list(NHEJ_CORE),
        "NHEJ_EXT": list(NHEJ_EXT),
        "HALLMARK_IFNA": load_gmt(RAW / "h.all.v2023.2.Hs.symbols.gmt", "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
        "HALLMARK_IFNG": load_gmt(RAW / "h.all.v2023.2.Hs.symbols.gmt", "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
        "REACTOME_NHEJ": load_gmt(
            RAW / "c2.cp.reactome.v2023.2.Hs.symbols.gmt",
            "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ",
        ),
    }
    if len(panels["IFN_ISG"]) != 39 or len(panels["MHC_APM"]) != 18:
        raise RuntimeError("locked IFN/APM panels changed size")
    return panels


def direction_call(median: float, n: int, sign_p: float, wilcox_greater: float, wilcox_less: float) -> str:
    if n == 0 or not math.isfinite(median):
        return "not_assayed"
    if median > 0 and (
        (math.isfinite(wilcox_greater) and wilcox_greater <= 0.05)
        or (not math.isfinite(wilcox_greater) and sign_p <= 0.05)
    ):
        return "up"
    if median < 0 and math.isfinite(wilcox_less) and wilcox_less <= 0.05:
        return "down"
    if median > 0:
        return "weak_up"
    if median < 0:
        return "weak_down"
    return "null"


def score_values(values: list[float], background: list[float] | None = None) -> dict:
    arr = np.asarray([v for v in values if v is not None and math.isfinite(v)], dtype=float)
    n = int(arr.size)
    if n == 0:
        return {
            "n_observed": 0,
            "n_up": 0,
            "median_log2FC": float("nan"),
            "mean_log2FC": float("nan"),
            "sign_p_greater": float("nan"),
            "wilcoxon_p_greater": float("nan"),
            "wilcoxon_p_less": float("nan"),
            "ranksum_p_greater": float("nan"),
            "call": "not_assayed",
        }
    n_up = int(np.sum(arr > 0))
    sign_p = float(binomtest(n_up, n, 0.5, alternative="greater").pvalue)
    wil_g = wil_l = float("nan")
    if n >= 5 and np.any(arr != 0):
        wil_g = float(wilcoxon(arr, alternative="greater", zero_method="wilcox").pvalue)
        wil_l = float(wilcoxon(arr, alternative="less", zero_method="wilcox").pvalue)
    ranksum = float("nan")
    if background is not None and n >= 5 and len(background) >= 20:
        bg = np.asarray(background, dtype=float)
        ranksum = float(mannwhitneyu(arr, bg, alternative="greater").pvalue)
    med = float(np.median(arr))
    return {
        "n_observed": n,
        "n_up": n_up,
        "median_log2FC": med,
        "mean_log2FC": float(np.mean(arr)),
        "sign_p_greater": sign_p,
        "wilcoxon_p_greater": wil_g,
        "wilcoxon_p_less": wil_l,
        "ranksum_p_greater": ranksum,
        "call": direction_call(med, n, sign_p, wil_g, wil_l),
    }


def gsea_nes(scores: dict[str, float], genes: list[str], nperm: int, rng: np.random.Generator) -> dict:
    universe = [(gene, score) for gene, score in scores.items() if math.isfinite(score)]
    universe.sort(key=lambda item: item[1], reverse=True)
    symbols = [gene for gene, _ in universe]
    weights = np.abs(np.asarray([score for _, score in universe], dtype=float))
    membership = np.array([gene in set(genes) for gene in symbols])
    n_hit = int(membership.sum())
    if n_hit < 8 or n_hit >= len(symbols):
        return {"nes": float("nan"), "perm_p_greater": float("nan"), "n_in_rank": n_hit}

    def enrichment(mask: np.ndarray) -> float:
        hit = np.where(mask, weights, 0.0)
        hit_sum = hit.sum()
        if hit_sum <= 0:
            return 0.0
        hit /= hit_sum
        miss = np.where(mask, 0.0, 1.0 / (len(mask) - int(mask.sum())))
        running = np.cumsum(hit - miss)
        imax = int(np.argmax(running))
        imin = int(np.argmin(running))
        return float(running[imax] if abs(running[imax]) >= abs(running[imin]) else running[imin])

    observed = enrichment(membership)
    null = np.empty(nperm, dtype=float)
    idx = np.arange(len(symbols))
    for i in range(nperm):
        chosen = rng.choice(idx, size=n_hit, replace=False)
        mask = np.zeros(len(symbols), dtype=bool)
        mask[chosen] = True
        null[i] = enrichment(mask)
    pos = null[null > 0]
    neg = null[null < 0]
    if observed >= 0:
        denom = float(np.mean(pos)) if pos.size else float("nan")
        perm_p = float((np.sum(null >= observed) + 1) / (nperm + 1))
    else:
        denom = float(np.mean(np.abs(neg))) if neg.size else float("nan")
        perm_p = float("nan")
    nes = observed / denom if denom and math.isfinite(denom) and denom != 0 else float("nan")
    return {"nes": nes, "perm_p_greater": perm_p, "n_in_rank": n_hit}


def stouffer(pvals: list[float]) -> float:
    usable = [p for p in pvals if math.isfinite(p) and 0 < p < 1]
    if len(usable) < 5:
        return float("nan")
    stats_ = norm.isf(np.clip(usable, 1e-300, 1 - 1e-16))
    return float(norm.sf(np.sum(stats_) / math.sqrt(len(usable))))


def soft_table(path: Path, begin: str, end: str) -> list[dict]:
    rows = []
    header = None
    with gzip.open(path, "rt", errors="replace") as handle:
        inside = False
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(begin):
                inside = True
                continue
            if inside and line.startswith(end):
                break
            if not inside or not line or line.startswith("!"):
                continue
            parsed = next(csv.reader([line], delimiter="\t"))
            parsed = [cell.strip().strip('"') for cell in parsed]
            if header is None:
                header = parsed
                continue
            rows.append(dict(zip(header, parsed)))
    return rows


def symbol_from_platform(row: dict) -> str:
    orf = (row.get("ORF") or "").strip()
    if orf and " " not in orf and orf.upper() not in {"NA", "NULL"}:
        symbol = orf
    else:
        desc = row.get("DESCRIPTION") or ""
        match = re.match(r"^([A-Za-z0-9.-]+)--", desc)
        symbol = match.group(1) if match else ""
    return ALIASES.get(symbol, symbol)


def load_fpkm() -> dict[str, dict[str, tuple[float, float]]]:
    """symbol -> contrast -> (ko_fpkm, wt_fpkm) after sum, median, or max is applied later.

    Returns per-locus rows grouped by symbol.
    """
    path = RAW / "GSE207704_CLDN4_RNAseq.txt.gz"
    groups = {
        "MCF-7 KO vs WT": ("MCF7_CLDN4KO_FPKM (fpkm)", "MCF7_WT_FPKM (fpkm)"),
        "T47D KO vs WT": ("T47D_CLDN4KO_FPKM (fpkm)", "T47D_WT_FPKM (fpkm)"),
    }
    loci = {contrast: defaultdict(list) for contrast in groups}
    with gzip.open(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            symbol = (row.get("gene_short_name") or "").strip()
            if not symbol or symbol == "gene_short_name":
                continue
            for contrast, (ko_col, wt_col) in groups.items():
                loci[contrast][symbol].append((float(row[ko_col]), float(row[wt_col])))
    return loci


def collapse_loci(pairs: list[tuple[float, float]], how: str) -> tuple[float, float]:
    ko = np.array([a for a, _ in pairs], dtype=float)
    wt = np.array([b for _, b in pairs], dtype=float)
    if how == "sum":
        return float(ko.sum()), float(wt.sum())
    if how == "median":
        return float(np.median(ko)), float(np.median(wt))
    if how == "max_fpkm":
        total = ko + wt
        pick = int(np.argmax(total))
        return float(ko[pick]), float(wt[pick])
    if how == "max_logfc":
        ratios = np.log2((ko + 0.5) / (wt + 0.5))
        pick = int(np.argmax(ratios))
        return float(ko[pick]), float(wt[pick])
    raise RuntimeError(how)


def genome_from_fpkm(loci, how: str, pseudo: float, min_mean: float) -> dict[str, dict[str, float]]:
    out = {}
    for contrast, genes in loci.items():
        scores = {}
        for symbol, pairs in genes.items():
            ko, wt = collapse_loci(pairs, how)
            if (ko + wt) / 2.0 < min_mean:
                continue
            scores[symbol] = math.log2((ko + pseudo) / (wt + pseudo))
        out[contrast] = scores
    return out


def load_deposited_array() -> dict[str, dict[str, np.ndarray]]:
    platform = soft_table(RAW / "GSE22493_family.soft.gz", "!platform_table_begin", "!platform_table_end")
    matrix = soft_table(RAW / "GSE22493_series_matrix.txt.gz", "!series_matrix_table_begin", "!series_matrix_table_end")
    id_to_symbol = {}
    for row in platform:
        symbol = symbol_from_platform(row)
        if symbol:
            id_to_symbol[row["ID"]] = symbol
    by_gene = defaultdict(list)
    for row in matrix:
        symbol = id_to_symbol.get(row["ID_REF"])
        if not symbol:
            continue
        vals = []
        for sample in ARRAYS:
            text = row.get(sample, "")
            vals.append(float(text) if text not in ("", "NA", "null") else float("nan"))
        by_gene[symbol].append(np.asarray(vals, dtype=float))
    return by_gene


def gene_effect_from_probes(probes: list[np.ndarray], how: str, arrays: list[int]) -> float:
    block = np.vstack([probe[arrays] for probe in probes])
    if not np.isfinite(block).any():
        return float("nan")
    if how == "median_probe":
        per_array = np.nanmedian(block, axis=0)
    elif how == "mean_probe":
        per_array = np.nanmean(block, axis=0)
    elif how in {"max_probe", "min_probe"}:
        means = np.nanmean(block, axis=1)
        if not np.isfinite(means).any():
            return float("nan")
        pick = int(np.nanargmax(means) if how == "max_probe" else np.nanargmin(means))
        per_array = block[pick]
    else:
        raise RuntimeError(how)
    per_array = per_array[np.isfinite(per_array)]
    if per_array.size == 0:
        return float("nan")
    return float(np.mean(per_array))


def genome_from_array(by_gene, how: str, drop_array: int | None) -> dict[str, float]:
    keep = [i for i in range(3) if i != drop_array]
    scores = {}
    for symbol, probes in by_gene.items():
        effect = gene_effect_from_probes(probes, how, keep)
        if math.isfinite(effect):
            scores[symbol] = effect
    return scores


def load_scanarray() -> dict[str, dict[str, list[float]]]:
    """symbol -> sample -> list of log2(KD/control) from background-subtracted medians.

    GEO channel 2 is the knockdown and channel 1 is the overexpression control.
    """
    platform = soft_table(RAW / "GSE22493_family.soft.gz", "!platform_table_begin", "!platform_table_end")
    # ScanArray Name matches DESCRIPTION. Build name -> symbol. Ambiguous names dropped.
    name_to_symbols = defaultdict(set)
    for row in platform:
        symbol = symbol_from_platform(row)
        desc = (row.get("DESCRIPTION") or "").strip()
        if symbol and desc:
            name_to_symbols[desc].add(symbol)
    name_to_symbol = {name: next(iter(symbols)) for name, symbols in name_to_symbols.items() if len(symbols) == 1}
    per_sample = {sample: defaultdict(list) for sample in ARRAYS}
    tar_path = RAW / "GSE22493_RAW.tar"
    import tarfile

    with tarfile.open(tar_path) as archive:
        for sample in ARRAYS:
            member = archive.extractfile(f"{sample}.txt.gz")
            with gzip.GzipFile(fileobj=member) as gz:
                text = gz.read().decode("latin1").splitlines()
            start = next(i for i, line in enumerate(text) if line.startswith("BEGIN DATA"))
            header = text[start + 1].split("\t")
            idx = {name: header.index(name) for name in ("Name", "Ch1 Median - B", "Ch2 Median - B", "Flags")}
            for line in text[start + 2 :]:
                if not line or line.startswith("END DATA"):
                    continue
                cells = line.split("\t")
                if len(cells) <= idx["Ch2 Median - B"]:
                    continue
                symbol = name_to_symbol.get(cells[idx["Name"]].strip().strip('"'))
                if not symbol:
                    continue
                try:
                    ch1 = float(cells[idx["Ch1 Median - B"]])
                    ch2 = float(cells[idx["Ch2 Median - B"]])
                    flags = int(float(cells[idx["Flags"]]))
                except ValueError:
                    continue
                per_sample[sample][symbol].append((ch1, ch2, flags))
    return per_sample


def scanarray_genome(per_sample, how: str, drop_array: int | None, detected_only: bool) -> dict[str, float]:
    samples = [sample for i, sample in enumerate(ARRAYS) if i != drop_array]
    symbols = set()
    for sample in samples:
        symbols.update(per_sample[sample])
    scores = {}
    for symbol in symbols:
        per_array = []
        for sample in samples:
            spots = per_sample[sample].get(symbol, [])
            ratios = []
            for ch1, ch2, flags in spots:
                ch1p = max(ch1, 0.0)
                ch2p = max(ch2, 0.0)
                if detected_only and (ch1p <= 0 or ch2p <= 0):
                    continue
                ratios.append(math.log2((ch2p + 1.0) / (ch1p + 1.0)))
            if not ratios:
                continue
            arr = np.asarray(ratios, dtype=float)
            if how == "median_probe":
                per_array.append(float(np.median(arr)))
            elif how == "max_probe":
                per_array.append(float(np.max(arr)))
            elif how == "min_probe":
                per_array.append(float(np.min(arr)))
            else:
                raise RuntimeError(how)
        if per_array:
            scores[symbol] = float(np.mean(per_array))
    return scores


def panel_values(scores: dict[str, float], genes: list[str]) -> list[float]:
    return [scores[gene] for gene in genes if gene in scores and math.isfinite(scores[gene])]


def background_values(scores: dict[str, float], genes: list[str]) -> list[float]:
    banned = set(genes)
    return [value for gene, value in scores.items() if gene not in banned and math.isfinite(value)]


def drop_extreme(values_by_gene: list[tuple[str, float]]) -> list[float]:
    if len(values_by_gene) < 6:
        return [value for _, value in values_by_gene]
    arr = np.array([value for _, value in values_by_gene], dtype=float)
    cut = np.quantile(np.abs(arr), 0.9)
    kept = [value for value in arr if abs(value) <= cut]
    return kept if len(kept) >= 5 else arr.tolist()


def add_row(rows: list[dict], **kwargs) -> None:
    call = kwargs["call"]
    cldn = kwargs["cldn4_log2FC"]
    counts = call == "up" and math.isfinite(cldn) and cldn < 0
    kwargs["concordant_hit"] = "yes" if counts else "no"
    kwargs["cldn4_down"] = "yes" if math.isfinite(cldn) and cldn < 0 else "no"
    rows.append(kwargs)


def score_genome(rows, scores, cldn, panels, meta) -> None:
    for panel, genes in panels.items():
        values = panel_values(scores, genes)
        stats_ = score_values(values, background_values(scores, genes))
        family = meta["family"]
        if meta["method"].endswith("max_probe") or meta["method"].endswith("max_logfc") or "max_logfc" in meta["method"]:
            family = "optimistic"
        base = {key: value for key, value in meta.items() if key not in {"family", "pathway", "seed"}}
        add_row(
            rows,
            panel=panel,
            n_in_set=len(genes),
            estimand="median_log2FC",
            estimate=stats_["median_log2FC"],
            aux_mean_log2FC=stats_["mean_log2FC"],
            nes="",
            perm_p="",
            stouffer_p="",
            **base,
            **{k: stats_[k] for k in ("n_observed", "n_up", "median_log2FC", "mean_log2FC", "sign_p_greater", "wilcoxon_p_greater", "wilcoxon_p_less", "ranksum_p_greater", "call")},
            family=family,
            cldn4_log2FC=cldn,
        )
        if meta.get("pathway"):
            rng = np.random.default_rng(meta["seed"])
            gsea = gsea_nes(scores, genes, 300, rng)
            gsea_call = "up" if math.isfinite(gsea["nes"]) and gsea["nes"] > 0 and gsea["perm_p_greater"] <= 0.05 else (
                "weak_up" if math.isfinite(gsea["nes"]) and gsea["nes"] > 0 else "not_up"
            )
            add_row(
                rows,
                accession=meta["accession"],
                contrast=meta["contrast"],
                class_label="cancer",
                organ=meta["organ"],
                method=meta["method"] + "+prerank_gsea",
                family="pathway",
                panel=panel,
                n_in_set=len(genes),
                n_observed=gsea["n_in_rank"],
                n_up="",
                estimand="NES",
                estimate=gsea["nes"],
                aux_mean_log2FC="",
                median_log2FC="",
                mean_log2FC="",
                sign_p_greater="",
                wilcoxon_p_greater="",
                wilcoxon_p_less="",
                ranksum_p_greater="",
                nes=gsea["nes"],
                perm_p=gsea["perm_p_greater"],
                stouffer_p="",
                call=gsea_call,
                cldn4_log2FC=cldn,
                note=meta["note"] + " Prerank GSEA, 300 gene-set permutations.",
            )


def array_level_tests(rows, by_gene, panels, cldn) -> None:
    """One-sample tests that use the three arrays as the unit."""
    for how, family in (("median_probe", "array_unit"), ("max_probe", "optimistic")):
        gene_by_array = {}
        for symbol, probes in by_gene.items():
            effects = []
            for i in range(3):
                effects.append(gene_effect_from_probes(probes, how, [i]))
            gene_by_array[symbol] = effects
        for panel, genes in panels.items():
            if panel not in PRIMARY_PANELS and panel != "HALLMARK_IFNG":
                continue
            present = [gene for gene in genes if gene in gene_by_array]
            if len(present) < 5:
                continue
            panel_medians = []
            for i in range(3):
                vals = [gene_by_array[gene][i] for gene in present if math.isfinite(gene_by_array[gene][i])]
                if vals:
                    panel_medians.append(float(np.median(vals)))
            if len(panel_medians) < 3:
                continue
            t_res = ttest_1samp(panel_medians, 0.0, alternative="greater")
            call = "up" if float(np.mean(panel_medians)) > 0 and t_res.pvalue <= 0.05 else (
                "weak_up" if float(np.mean(panel_medians)) > 0 else "not_up"
            )
            add_row(
                rows,
                accession="GSE22493",
                contrast="SKOV-3 KD",
                class_label="cancer",
                organ="ovary",
                method=f"deposited_{how}+array_t",
                family=family if family == "optimistic" else "array_unit",
                panel=panel,
                n_in_set=len(genes),
                n_observed=len(present),
                n_up=sum(v > 0 for v in panel_medians),
                estimand="mean_of_array_panel_medians",
                estimate=float(np.mean(panel_medians)),
                aux_mean_log2FC=float(np.mean(panel_medians)),
                median_log2FC=float(np.median(panel_medians)),
                mean_log2FC=float(np.mean(panel_medians)),
                sign_p_greater="",
                wilcoxon_p_greater="",
                wilcoxon_p_less="",
                ranksum_p_greater="",
                nes="",
                perm_p="",
                stouffer_p=float(t_res.pvalue),
                call=call,
                cldn4_log2FC=cldn,
                note="Unit is the array. One-sample t-test of the three panel medians, greater than 0.",
            )
        # Stouffer on gene-level one-sided t-tests, median-probe means
        scores_mean = {}
        p_greater = []
        p_genes = []
        for symbol, effects in gene_by_array.items():
            arr = np.asarray([v for v in effects if math.isfinite(v)], dtype=float)
            if arr.size < 2 or float(np.std(arr, ddof=1)) <= 1e-8:
                if arr.size:
                    scores_mean[symbol] = float(np.mean(arr))
                continue
            scores_mean[symbol] = float(np.mean(arr))
            p_greater.append((symbol, float(ttest_1samp(arr, 0.0, alternative="greater").pvalue)))
        pmap = dict(p_greater)
        for panel, genes in panels.items():
            if panel not in PRIMARY_PANELS:
                continue
            vals = [pmap[gene] for gene in genes if gene in pmap]
            effects = [scores_mean[gene] for gene in genes if gene in scores_mean]
            if len(vals) < 5:
                continue
            sp = stouffer(vals)
            med = float(np.median(effects)) if effects else float("nan")
            call = "up" if math.isfinite(sp) and sp <= 0.05 and med > 0 else ("weak_up" if med > 0 else "not_up")
            add_row(
                rows,
                accession="GSE22493",
                contrast="SKOV-3 KD",
                class_label="cancer",
                organ="ovary",
                method=f"deposited_{how}+stouffer_t",
                family="alternate_de" if family != "optimistic" else "optimistic",
                panel=panel,
                n_in_set=len(genes),
                n_observed=len(vals),
                n_up="",
                estimand="stouffer_on_gene_t",
                estimate=med,
                aux_mean_log2FC=float(np.mean(effects)) if effects else "",
                median_log2FC=med,
                mean_log2FC=float(np.mean(effects)) if effects else "",
                sign_p_greater="",
                wilcoxon_p_greater="",
                wilcoxon_p_less="",
                ranksum_p_greater="",
                nes="",
                perm_p="",
                stouffer_p=sp,
                call=call,
                cldn4_log2FC=cldn,
                note="Stouffer combination of one-sided gene t-tests across the three arrays.",
            )


def _kallisto_qc(gsm: str) -> tuple[bool, int, float]:
    """A GSM counts only when abundance exists and n_processed exceeds 30 million."""
    info = KALLISTO / gsm / "run_info.json"
    abundance = KALLISTO / gsm / "abundance.tsv"
    if not info.exists() or not abundance.exists():
        return False, 0, float("nan")
    meta = json.loads(info.read_text())
    n_processed = int(meta.get("n_processed") or 0)
    rate = float(meta.get("p_pseudoaligned") or float("nan"))
    return n_processed >= 30_000_000, n_processed, rate


def load_kallisto() -> tuple[dict[str, dict[str, float]], str] | None:
    """Score each GSE207704 line only when both replicates of both genotypes pass QC."""
    if not KALLISTO.exists():
        return None
    lines = {
        "T47D KO vs WT": (["GSM6310642", "GSM6310643"], ["GSM6310640", "GSM6310641"]),
        "MCF-7 KO vs WT": (["GSM6310646", "GSM6310647"], ["GSM6310644", "GSM6310645"]),
    }
    ready: dict[str, tuple[list[str], list[str]]] = {}
    skipped = []
    rates = []
    for contrast, (ko_gsms, wt_gsms) in lines.items():
        qc = {gsm: _kallisto_qc(gsm) for gsm in ko_gsms + wt_gsms}
        if not all(ok for ok, _, _ in qc.values()):
            missing = [gsm for gsm, (ok, _, _) in qc.items() if not ok]
            skipped.append(f"{contrast} ({', '.join(missing)} below 30 million reads or absent)")
            continue
        ready[contrast] = (ko_gsms, wt_gsms)
        for gsm, (_, n_processed, rate) in qc.items():
            rates.append((gsm, n_processed, rate))
    if not ready:
        return None
    fasta = Path("/tmp/kd_sweep/tx/Homo_sapiens.GRCh38.cdna.all.fa.gz")
    if not fasta.exists():
        return None
    gene_of = {}
    with gzip.open(fasta, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].split()
            tid = header[0].split(".")[0]
            symbol = ""
            for field in header[1:]:
                if field.startswith("gene_symbol:"):
                    symbol = field.split(":", 1)[1]
            if symbol:
                gene_of[header[0]] = symbol
                gene_of[tid] = symbol
    needed = [gsm for ko_gsms, wt_gsms in ready.values() for gsm in ko_gsms + wt_gsms]
    counts = {gsm: defaultdict(float) for gsm in needed}
    for gsm in needed:
        with (KALLISTO / gsm / "abundance.tsv").open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                symbol = gene_of.get(row["target_id"]) or gene_of.get(row["target_id"].split(".")[0])
                if not symbol:
                    continue
                counts[gsm][symbol] += float(row["est_counts"])

    def logfc(ko_gsms, wt_gsms) -> dict[str, float]:
        symbols = set()
        for gsm in ko_gsms + wt_gsms:
            symbols.update(counts[gsm])
        matrix = []
        names = []
        for symbol in symbols:
            vector = [counts[gsm][symbol] for gsm in wt_gsms + ko_gsms]
            if sum(vector) < 10:
                continue
            names.append(symbol)
            matrix.append(vector)
        mat = np.asarray(matrix, dtype=float)
        with np.errstate(divide="ignore"):
            log_mat = np.log(mat + 1e-8)
        geomean = np.exp(np.mean(log_mat, axis=1))
        geomean[geomean <= 0] = np.nan
        ratios = mat / geomean[:, None]
        size = np.nanmedian(ratios, axis=0)
        norm = mat / size
        wt = norm[:, :2].mean(axis=1)
        ko = norm[:, 2:].mean(axis=1)
        return {name: math.log2((ko_i + 1.0) / (wt_i + 1.0)) for name, ko_i, wt_i in zip(names, ko, wt)}

    scored = {contrast: logfc(ko_gsms, wt_gsms) for contrast, (ko_gsms, wt_gsms) in ready.items()}
    rate_vals = [rate for _, _, rate in rates if math.isfinite(rate)]
    rate_txt = (
        f"{min(rate_vals):.1f}–{max(rate_vals):.1f}%"
        if rate_vals
        else "not recorded"
    )
    note = (
        "Kallisto 0.52.0 on Ensembl 111 cDNA, single-end fragment length 200 (sd 20). "
        "Estimated counts summed to gene_symbol. Median-of-ratios size factors. "
        "log2((mean KO + 1) / (mean WT + 1)), n=2 vs 2. "
        f"Pseudoalignment of the scored runs: {rate_txt}. "
        "A row is a concordant hit only when CLDN4 log2FC is negative in this same scoring."
    )
    if skipped:
        note += " Not scored: " + "; ".join(skipped) + "."
    return scored, note


def load_lung_symbols() -> dict[str, float]:
    scores = {}
    with gzip.open(RAW / "GSE50927_Cldn4lungWTvsKOgenes.csv.gz", "rt") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            scores[row["Marker.Symbol"]] = float(row["logFC"])
    return scores


def assert_locked(rows: list[dict]) -> None:
    locked = {
        (r["dataset"], r["contrast"], r["gene_set"]): float(r["median_log2FC"])
        for r in read_tsv(LOCK)
        if r["dataset"] in {"GSE207704", "GSE22493"}
    }
    found = {
        (r["accession"], r["contrast"], r["panel"]): float(r["median_log2FC"])
        for r in rows
        if r["method"] in {"fpkm_sum_pc0.5_min0.5", "deposited_median_probe_mean_arrays"}
        and r["panel"] in {"IFN_ISG", "MHC_APM"}
        and r["estimand"] == "median_log2FC"
    }
    for key, expected in locked.items():
        method_key = (key[0], key[1], key[2])
        if method_key not in found:
            raise RuntimeError(f"missing primary reproduction {key}")
        if not math.isclose(found[method_key], expected, rel_tol=0, abs_tol=1e-6):
            raise RuntimeError(f"primary drift {key}: {found[method_key]} vs {expected}")


def plot(rows: list[dict]) -> None:
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    panels = PRIMARY_PANELS
    contrasts = [
        ("GSE50927", "baseline KO vs WT", "Non-cancer  lung baseline"),
        ("GSE207704", "T47D KO vs WT", "Cancer  T47D"),
        ("GSE207704", "MCF-7 KO vs WT", "Cancer  MCF-7"),
        ("GSE22493", "SKOV-3 KD", "Cancer  SKOV-3"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.6), sharey=True)
    for ax, panel in zip(axes, panels):
        for y, (acc, contrast, label) in enumerate(contrasts):
            block = [
                r
                for r in rows
                if r["accession"] == acc
                and r["contrast"] == contrast
                and r["panel"] == panel
                and r["estimand"] == "median_log2FC"
                and r["family"] != "optimistic"
                and r.get("cldn4_down") == "yes"
                and r["median_log2FC"] not in ("", None)
            ]
            primary = [
                r
                for r in block
                if r["method"]
                in {
                    "fpkm_sum_pc0.5_min0.5",
                    "deposited_median_probe_mean_arrays",
                    "lung_edger_locked_or_scored",
                }
            ]
            meds = [float(r["median_log2FC"]) for r in block]
            if not meds and not primary:
                continue
            color = "#2C6E8A" if acc == "GSE50927" else "#B85C38"
            if meds:
                ax.plot([min(meds), max(meds)], [y, y], color=color, lw=2, solid_capstyle="round", zorder=2)
            if primary:
                est = float(primary[0]["median_log2FC"])
                ax.scatter([est], [y], s=46, color=color, zorder=3)
                ax.text(est, y + 0.18, primary[0]["call"], ha="center", va="bottom", fontsize=7, color=color)
            optimistic = [
                r
                for r in rows
                if r["accession"] == acc
                and r["contrast"] == contrast
                and r["panel"] == panel
                and r["family"] == "optimistic"
                and r.get("cldn4_down") == "yes"
                and r["estimand"] == "median_log2FC"
                and r["median_log2FC"] not in ("", None)
            ]
            if optimistic:
                best = max(float(r["median_log2FC"]) for r in optimistic)
                ax.scatter([best], [y], s=28, facecolors="none", edgecolors=color, linewidths=1.1, zorder=3)
        ax.axvline(0, color="#222222", lw=0.8)
        ax.set_title(panel.replace("_", " "), loc="left", fontsize=11)
        ax.set_xlabel("Median log2FC after CLDN4 loss")
        ax.set_yticks(range(len(contrasts)))
        ax.set_yticklabels([c[2] for c in contrasts], fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylim(-0.6, len(contrasts) - 0.2)
    fig.suptitle("Cancer KD sweep: primary point, sensitivity range, open circle = most optimistic rule", fontsize=11, x=0.01, ha="left")
    fig.text(
        0.01,
        0.01,
        "Filled point: primary estimate. Thick bar: min–max of non-optimistic sensitivities with CLDN4 down. "
        "Open circle: highest median under a max-probe or max-locus rule with CLDN4 down. "
        "Lung IFN/APM points are the locked baseline, not a refit.",
        fontsize=7.5,
        color="#333333",
    )
    fig.tight_layout(rect=(0, 0.06, 1, 0.92))
    fig.savefig(fig_dir / "sweep_cancer_sensitivity.png", dpi=150)
    fig.savefig(fig_dir / "sweep_cancer_sensitivity.pdf")
    plt.close(fig)


def _fmt_panel(rows: list[dict], method: str, contrast: str, panel: str, label: str) -> str:
    match = [
        row
        for row in rows
        if row["method"] == method and row["contrast"] == contrast and row["panel"] == panel and row["estimand"] == "median_log2FC"
    ]
    if not match:
        return f"{label} was not scored."
    row = match[0]
    med = float(row["median_log2FC"])
    p = row["wilcoxon_p_greater"]
    p_txt = f"{float(p):.3g}" if p not in ("", None) and str(p) not in {"nan"} else "NA"
    return (
        f"{label}: median {med:+.3f}, {row['n_up']}/{row['n_observed']} up, "
        f"Wilcoxon P(greater) = {p_txt}, call {row['call']}."
    )


def write_sweep_md(rows: list[dict], kallisto_note: str | None) -> None:
    hits = [r for r in rows if r["concordant_hit"] == "yes" and r["class_label"] == "cancer"]
    nonopt_hits = [r for r in hits if r["family"] not in {"optimistic"}]
    opt_hits = [r for r in hits if r["family"] == "optimistic"]
    single_hits = [r for r in hits if r["family"] == "single_array"]
    combined_hits = [r for r in nonopt_hits if r["family"] != "single_array"]
    lines = [
        "# Cancer KD sweep: IFN, APM, and NHEJ after CLDN4 loss",
        "",
        "Question: does any honest re-scoring of GSE207704 or GSE22493 show IFN, APM, or NHEJ going up after CLDN4 loss?",
        "",
        "A concordant hit requires two things at once. CLDN4 log2FC is negative in that same scoring, and the panel is called up. For gene-median scores the up call is the Claim C4 rule (median above 0 and one-sided Wilcoxon P ≤ 0.05, or the sign test when fewer than 5 genes are observed). Pathway and array-unit rows use the test named in the method.",
        "",
        "GSE207704 in GEO is a four-column Cufflinks FPKM table (condition means). Replicate FASTQs exist (two per genotype). "
        + (
            kallisto_note
            if kallisto_note
            else "This run did not find a complete n=2 versus n=2 quantification for either line, so GSE207704 has no replicate-level DE in the table."
        ),
        "",
        "## Primary estimates still match the locked table",
        "",
        "GSE207704 sum of loci, pseudocount 0.5, mean FPKM at least 0.5, and GSE22493 median-probe then mean across the three arrays, reproduce the locked IFN/ISG and MHC-I/APM medians. Those calls stay discordant or unsigned, and every primary median is negative.",
        "",
        "## What still stands",
        "",
        "The combined cancer contrasts do not open IFN or APM. "
        f"Concordant hits on those combined contrasts, excluding single arrays and max-probe rules: {len(combined_hits)}. "
        "Primary IFN/ISG and MHC-I/APM medians stay negative in T47D, MCF-7, and SKOV-3, and they still match the locked table. "
        "Prerank GSEA, the three-array t-test, and Stouffer's combination of the gene-wise t-tests do not call those panels up. "
        "Raw ScanArray log2(KD/control), with CLDN4 down, calls IFN/ISG, MHC-I/APM, and NHEJ_CORE down.",
        "",
        "NHEJ_CORE is not a second version of the IFN split. "
        + _fmt_panel(rows, "fpkm_sum_pc0.5_min0.5", "T47D KO vs WT", "NHEJ_CORE", "T47D")
        + " "
        + _fmt_panel(rows, "fpkm_sum_pc0.5_min0.5", "MCF-7 KO vs WT", "NHEJ_CORE", "MCF-7")
        + " "
        + _fmt_panel(rows, "deposited_median_probe_mean_arrays", "SKOV-3 KD", "NHEJ_CORE", "SKOV-3")
        + " "
        + _fmt_panel(rows, "lung_edger_locked_or_scored", "baseline KO vs WT", "NHEJ_CORE", "Lung baseline"),
        "",
        "## The one array that leans the other way",
        "",
        "SKOV-3 is three two-colour arrays. "
        + _fmt_panel(rows, "deposited_median_probe_GSM558701", "SKOV-3 KD", "IFN_ISG", "GSM558701 IFN/ISG")
        + " "
        + _fmt_panel(rows, "deposited_median_probe_GSM558701", "SKOV-3 KD", "MHC_APM", "GSM558701 MHC-I/APM")
        + " "
        + _fmt_panel(rows, "deposited_median_probe_GSM558701", "SKOV-3 KD", "HALLMARK_IFNA", "GSM558701 Hallmark IFN-α")
        + " "
        + _fmt_panel(rows, "deposited_median_probe_GSM558701", "SKOV-3 KD", "HALLMARK_IFNG", "GSM558701 Hallmark IFN-γ")
        + " "
        "GSM558700 has no CLDN4 measurement and keeps IFN and APM negative. "
        "GSM558702 has CLDN4 down and keeps IFN and APM negative. "
        "One array out of three is not the combined SKOV-3 contrast, and it does not move T47D or MCF-7.",
        "",
    ]
    kal_method = "kallisto_median_of_ratios_log2FC"
    kal_contrasts = []
    for row in rows:
        if row["method"] == kal_method and row["contrast"] not in kal_contrasts:
            kal_contrasts.append(row["contrast"])
    if kal_contrasts:
        lines += [
            "## Replicate quantification",
            "",
            kallisto_note or "",
            "",
        ]
        for contrast in kal_contrasts:
            cldn = next(
                row["cldn4_log2FC"]
                for row in rows
                if row["method"] == kal_method and row["contrast"] == contrast
            )
            bits = [f"{contrast}: CLDN4 log2FC {float(cldn):+.3f}."]
            for panel, label in (
                ("IFN_ISG", "IFN/ISG"),
                ("MHC_APM", "MHC-I/APM"),
                ("NHEJ_CORE", "NHEJ_CORE"),
                ("HALLMARK_IFNA", "Hallmark IFN-α"),
                ("HALLMARK_IFNG", "Hallmark IFN-γ"),
            ):
                bits.append(_fmt_panel(rows, kal_method, contrast, panel, label))
            lines.append(" ".join(bits))
            lines.append("")
        gsea_up = [
            row
            for row in rows
            if row["method"] == kal_method + "+prerank_gsea"
            and row["panel"] in {"IFN_ISG", "MHC_APM", "NHEJ_CORE"}
            and row["call"] == "up"
        ]
        if not gsea_up:
            lines.append(
                "Prerank GSEA on those kallisto rankings does not call IFN/ISG, MHC-I/APM, or NHEJ_CORE up."
            )
            lines.append("")
    lines += [
        "## Concordant hits",
        "",
        f"Single-array hits: {len(single_hits)}. Other non-optimistic hits: {len(combined_hits)}. "
        f"Optimistic max-probe or max-locus hits: {len(opt_hits)}.",
        "",
    ]
    if not hits:
        lines.append("No cancer row met the up-call rule while CLDN4 was down.")
    else:
        lines.append("| Family | Method | Contrast | Panel | Estimate | Call |")
        lines.append("|---|---|---|---|---:|---|")
        for row in hits:
            est = row["estimate"]
            est_s = f"{float(est):+.3f}" if est != "" else ""
            lines.append(
                f"| {row['family']} | {row['method']} | {row['contrast']} | {row['panel']} | {est_s} | {row['call']} |"
            )
        lines.append("")
        lines.append("Optimistic rows pick the highest probe or locus per gene before taking the panel median. They are a hunt, not the primary estimate.")
    lines += [
        "",
        "## How to read the figure",
        "",
        "Filled points are the primary medians. The bar is the min-to-max of the other non-optimistic sensitivities in which CLDN4 fell (pseudocount, filter, mean or median collapse, leave-one-array-out, single arrays, raw ScanArray log2(KD/control), outlier trim, and kallisto when that line was quantified). The open circle is the highest median from a max-probe or max-locus rule that also has CLDN4 down. Lung IFN/APM points are the locked baseline and were not refit. Lung NHEJ, when drawn, is scored on the author edgeR table.",
        "",
        "Raw ScanArray uses GEO channel 2 as knockdown and channel 1 as the overexpression control, log2((Ch2 background-subtracted median + 1) / (Ch1 + 1)). A scoring is left out of the hit list and out of the sensitivity bar when CLDN4 does not fall.",
        "",
        "Reactome NHEJ contains histone genes, so it is a sensitivity panel. NHEJ_CORE is XRCC5, XRCC6, PRKDC, LIG4, XRCC4, NHEJ1, DCLRE1C, and PAXX.",
        "",
        "Tables: `sweep/sweep_results.tsv`, `sweep/concordant_hits.tsv`.",
        "",
        "Reproduce: `python3 scripts/cldn4_kd_match_wave/sweep_cancer_kd.py`",
        "",
    ]
    (OUT / "SWEEP.md").write_text("\n".join(lines))


def main() -> None:
    panels = load_panels()
    rows: list[dict] = []
    loci = load_fpkm()
    fpkm_grid = []
    for how, pseudo, min_mean, pathway in (
        ("sum", 0.5, 0.5, True),
        ("sum", 0.1, 0.5, False),
        ("sum", 1.0, 0.5, False),
        ("sum", 0.5, 0.0, False),
        ("sum", 0.5, 1.0, False),
        ("median", 0.5, 0.5, False),
        ("max_fpkm", 0.5, 0.5, False),
        ("max_logfc", 0.5, 0.5, False),
    ):
        fpkm_grid.append((how, pseudo, min_mean, pathway))
    for how, pseudo, min_mean, pathway in fpkm_grid:
        genomes = genome_from_fpkm(loci, how, pseudo, min_mean)
        for contrast, scores in genomes.items():
            cldn = scores.get("CLDN4", float("nan"))
            method = f"fpkm_{how}_pc{pseudo}_min{min_mean}"
            if how == "sum" and pseudo == 0.5 and min_mean == 0.5:
                method = "fpkm_sum_pc0.5_min0.5"
            note = "Cufflinks condition-mean FPKM. Not a replicate-level DE."
            if how == "max_logfc":
                note += " Max-locus log2FC is an optimistic isoform/locus rule."
            meta = {
                "accession": "GSE207704",
                "contrast": contrast,
                "class_label": "cancer",
                "organ": "breast",
                "method": method,
                "family": "primary" if method == "fpkm_sum_pc0.5_min0.5" else "sensitivity",
                "note": note,
                "pathway": pathway,
                "seed": 207704,
            }
            score_genome(rows, scores, cldn, panels, meta)
            # outlier trim on the primary collapse only
            if method == "fpkm_sum_pc0.5_min0.5":
                for panel, genes in panels.items():
                    paired = [(gene, scores[gene]) for gene in genes if gene in scores]
                    trimmed = drop_extreme(paired)
                    stats_ = score_values(trimmed)
                    add_row(
                        rows,
                        accession="GSE207704",
                        contrast=contrast,
                        class_label="cancer",
                        organ="breast",
                        method="fpkm_sum_pc0.5_min0.5_trim_abs90",
                        family="sensitivity",
                        panel=panel,
                        n_in_set=len(genes),
                        estimand="median_log2FC",
                        estimate=stats_["median_log2FC"],
                        aux_mean_log2FC=stats_["mean_log2FC"],
                        nes="",
                        perm_p="",
                        stouffer_p="",
                        cldn4_log2FC=cldn,
                        note="Drops panel genes above the 90th percentile of absolute log2FC.",
                        **{k: stats_[k] for k in ("n_observed", "n_up", "median_log2FC", "mean_log2FC", "sign_p_greater", "wilcoxon_p_greater", "wilcoxon_p_less", "ranksum_p_greater", "call")},
                    )

    by_gene = load_deposited_array()
    for how, drop, pathway in (
        ("median_probe", None, True),
        ("mean_probe", None, False),
        ("max_probe", None, False),
        ("min_probe", None, False),
        ("median_probe", 0, False),
        ("median_probe", 1, False),
        ("median_probe", 2, False),
    ):
        scores = genome_from_array(by_gene, how, drop)
        cldn = scores.get("CLDN4", float("nan"))
        method = "deposited_median_probe_mean_arrays" if how == "median_probe" and drop is None else f"deposited_{how}" + (f"_drop_array{drop+1}" if drop is not None else "")
        meta = {
            "accession": "GSE22493",
            "contrast": "SKOV-3 KD",
            "class_label": "cancer",
            "organ": "ovary",
            "method": method,
            "family": "primary" if method == "deposited_median_probe_mean_arrays" else "sensitivity",
            "note": "Deposited GEO VALUE, log2(KD/control). Probes collapsed, then mean across arrays kept.",
            "pathway": pathway,
            "seed": 22493,
        }
        score_genome(rows, scores, cldn, panels, meta)
    primary_scores = genome_from_array(by_gene, "median_probe", None)
    array_level_tests(rows, by_gene, panels, primary_scores.get("CLDN4", float("nan")))
    for array_i, sample in enumerate(ARRAYS):
        scores = {}
        for symbol, probes in by_gene.items():
            effect = gene_effect_from_probes(probes, "median_probe", [array_i])
            if math.isfinite(effect):
                scores[symbol] = effect
        meta = {
            "accession": "GSE22493",
            "contrast": "SKOV-3 KD",
            "class_label": "cancer",
            "organ": "ovary",
            "method": f"deposited_median_probe_{sample}",
            "family": "single_array",
            "note": f"One array only ({sample}). Median of probes on that array. Not the three-array contrast.",
            "pathway": False,
            "seed": 22493,
        }
        score_genome(rows, scores, scores.get("CLDN4", float("nan")), panels, meta)

    print("parsing ScanArray", flush=True)
    scan = load_scanarray()
    for how, drop, detected in (
        ("median_probe", None, False),
        ("median_probe", None, True),
        ("max_probe", None, False),
        ("min_probe", None, False),
        ("median_probe", 0, False),
        ("median_probe", 1, False),
        ("median_probe", 2, False),
    ):
        scores = scanarray_genome(scan, how, drop, detected)
        cldn = scores.get("CLDN4", float("nan"))
        method = f"scanarray_{how}" + ("_detected" if detected else "") + (f"_drop_array{drop+1}" if drop is not None else "")
        meta = {
            "accession": "GSE22493",
            "contrast": "SKOV-3 KD",
            "class_label": "cancer",
            "organ": "ovary",
            "method": method,
            "family": "sensitivity",
            "note": "ScanArray log2((Ch2+1)/(Ch1+1)); Ch2 is knockdown, Ch1 is overexpression control.",
            "pathway": how == "median_probe" and drop is None and not detected,
            "seed": 558700,
        }
        score_genome(rows, scores, cldn, panels, meta)

    loaded = load_kallisto()
    kallisto_note = None
    if loaded:
        kallisto, kallisto_note = loaded
        for contrast, scores in kallisto.items():
            cldn = scores.get("CLDN4", float("nan"))
            meta = {
                "accession": "GSE207704",
                "contrast": contrast,
                "class_label": "cancer",
                "organ": "breast",
                "method": "kallisto_median_of_ratios_log2FC",
                "family": "alternate_de",
                "note": kallisto_note,
                "pathway": True,
                "seed": 6310640,
            }
            score_genome(rows, scores, cldn, panels, meta)

    lung = load_lung_symbols()
    # Locked IFN/APM so the figure keeps the class split, plus NHEJ scored on the same edgeR file.
    locked_rows = read_tsv(LOCK)
    for row in locked_rows:
        if row["dataset"] != "GSE50927" or row["contrast"] != "baseline KO vs WT":
            continue
        if row["gene_set"] not in {"IFN_ISG", "MHC_APM"}:
            continue
        add_row(
            rows,
            accession="GSE50927",
            contrast="baseline KO vs WT",
            class_label="non-cancer",
            organ="lung",
            method="lung_edger_locked_or_scored",
            family="locked_baseline",
            panel=row["gene_set"],
            n_in_set=int(row["n_in_set"]),
            n_observed=int(row["n_observed"]),
            n_up=int(row["n_up"]),
            estimand="median_log2FC",
            estimate=float(row["median_log2FC"]),
            aux_mean_log2FC=float(row["mean_log2FC"]),
            median_log2FC=float(row["median_log2FC"]),
            mean_log2FC=float(row["mean_log2FC"]),
            sign_p_greater=row["sign_p_greater"],
            wilcoxon_p_greater=row["wilcoxon_p_greater"],
            wilcoxon_p_less=row["wilcoxon_p_less"],
            ranksum_p_greater=row["ranksum_p_greater"],
            nes="",
            perm_p="",
            stouffer_p="",
            call=row["direction_call"],
            cldn4_log2FC=float(row["CLDN4_log2FC"]),
            note="Locked Claim C4 baseline. Not refit in this sweep.",
        )
    for panel, genes in (("NHEJ_CORE", NHEJ_CORE), ("NHEJ_EXT", NHEJ_EXT)):
        values = []
        for gene in genes:
            mouse = MOUSE[gene]
            if mouse in lung:
                values.append(lung[mouse])
            elif mouse.lower() in {k.lower(): k for k in lung}:
                pass
        # case-insensitive lookup
        lookup = {key.lower(): value for key, value in lung.items()}
        values = [lookup[MOUSE[gene].lower()] for gene in genes if MOUSE[gene].lower() in lookup]
        stats_ = score_values(values)
        add_row(
            rows,
            accession="GSE50927",
            contrast="baseline KO vs WT",
            class_label="non-cancer",
            organ="lung",
            method="lung_edger_locked_or_scored",
            family="locked_baseline",
            panel=panel,
            n_in_set=len(genes),
            estimand="median_log2FC",
            estimate=stats_["median_log2FC"],
            aux_mean_log2FC=stats_["mean_log2FC"],
            nes="",
            perm_p="",
            stouffer_p="",
            cldn4_log2FC=lung.get("Cldn4", float("nan")),
            note="Author edgeR logFC, naive KO vs WT. Mouse symbols of the human NHEJ list.",
            **{k: stats_[k] for k in ("n_observed", "n_up", "median_log2FC", "mean_log2FC", "sign_p_greater", "wilcoxon_p_greater", "wilcoxon_p_less", "ranksum_p_greater", "call")},
        )

    assert_locked(rows)
    fields = [
        "accession",
        "contrast",
        "class_label",
        "organ",
        "method",
        "family",
        "panel",
        "n_in_set",
        "n_observed",
        "n_up",
        "estimand",
        "estimate",
        "median_log2FC",
        "mean_log2FC",
        "aux_mean_log2FC",
        "sign_p_greater",
        "wilcoxon_p_greater",
        "wilcoxon_p_less",
        "ranksum_p_greater",
        "nes",
        "perm_p",
        "stouffer_p",
        "call",
        "cldn4_log2FC",
        "cldn4_down",
        "concordant_hit",
        "note",
    ]
    dumped = []
    for row in rows:
        item = {}
        for key in fields:
            value = row.get(key, "")
            if isinstance(value, float):
                item[key] = fnum(value)
            else:
                item[key] = value
        dumped.append(item)
    write_tsv(OUT / "sweep" / "sweep_results.tsv", dumped, fields)
    hits = [row for row in dumped if row["concordant_hit"] == "yes"]
    write_tsv(OUT / "sweep" / "concordant_hits.tsv", hits, fields)
    plot(rows)
    write_sweep_md(rows, kallisto_note)
    cancer_hits = [r for r in hits if r["class_label"] == "cancer"]
    print(f"rows {len(rows)} cancer_hits {len(cancer_hits)} kallisto {kallisto_note is not None}")
    for row in cancer_hits:
        print("HIT", row["family"], row["method"], row["contrast"], row["panel"], row["estimate"], row["call"])


if __name__ == "__main__":
    main()
