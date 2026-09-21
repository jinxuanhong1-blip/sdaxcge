#!/usr/bin/env python3
"""Method sweep: CLDN4 down plus IFN, APM, and NHEJ up.

Pre-specified primary row, fixed before ranking:
  transform log2(x+1); for counts, x is CPM
  set summary = median log2FC
  IFN = MSigDB HALLMARK_INTERFERON_GAMMA_RESPONSE
  APM = Reactome antigen presentation, folding, assembly and peptide loading of class I MHC
  NHEJ = Reactome nonhomologous end joining
  contrasts = ADC monotherapy versus control

The grid also varies pseudocount, mean versus median, alternate gene sets,
and the combination arms. signed_score = (-CLDN4 log2FC) + IFN + APM + NHEJ
using that row's summary. A sign hit requires CLDN4 log2FC < 0 and all three
set summaries > 0. Thresholds are applied to those rows; they are not a
second hidden search.

E-MTAB-16433 uses the mouse pseudobulk counts deposited with PR 294.
ftp.ebi.ac.uk TLS failed here, so the cell matrix was not re-downloaded.
E-MTAB-16843 and E-MTAB-16849 were not scored for the same reason.
"""

from __future__ import annotations

import csv
import gzip
import math
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"

REACTOME_APM = "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC"
REACTOME_NHEJ = "REACTOME_NONHOMOLOGOUS_END_JOINING_NHEJ"
REACTOME_IFNG = "REACTOME_INTERFERON_GAMMA_SIGNALING"
HALLMARK_IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"
HALLMARK_IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE"

CANONICAL_APM = {
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP",
    "PSMB8", "PSMB9", "PSMB10", "PSME1", "PSME2", "ERAP1", "ERAP2",
    "CALR", "PDIA3", "CANX", "NLRC5",
}
CLASSICAL_NHEJ = {"XRCC6", "XRCC5", "PRKDC", "LIG4", "XRCC4", "NHEJ1", "DCLRE1C", "PAXX"}

SET_BUNDLES = {
    "primary_hallmarkG_reactomeAPM_reactomeNHEJ": {
        "ifn_name": HALLMARK_IFNG,
        "apm_name": REACTOME_APM,
        "nhej_name": REACTOME_NHEJ,
    },
    "alt_hallmarkA_canonicalAPM_coreNHEJ": {
        "ifn_name": HALLMARK_IFNA,
        "apm_name": "CANONICAL_APM_18",
        "nhej_name": "CLASSICAL_NHEJ_8",
    },
    "alt_reactomeIFNG_reactomeAPM_reactomeNHEJ": {
        "ifn_name": REACTOME_IFNG,
        "apm_name": REACTOME_APM,
        "nhej_name": REACTOME_NHEJ,
    },
}


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 research-hunt"})
    with urllib.request.urlopen(req, timeout=180) as resp, dest.open("wb") as out:
        out.write(resp.read())


def load_gmt(path: Path) -> dict[str, set[str]]:
    out = {}
    for line in path.read_text().splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 3:
            out[parts[0]] = set(parts[2:])
    return out


def read_table(path: Path):
    raw = gzip.open(path, "rb").read()
    if raw.startswith(b"\xff\xfe"):
        text = raw.decode("utf-16")
    else:
        text = raw.decode("utf-8")
    return list(csv.DictReader(text.splitlines(), delimiter="\t"))


def collapse_symbols(rows, value_cols, biotype_filter: bool):
    best = {}
    for row in rows:
        if biotype_filter and row.get("gene_biotype") not in (None, "", "protein_coding"):
            if row.get("gene_biotype") != "protein_coding":
                continue
        symbol = (row.get("gene_name") or row.get("barcode") or "").strip()
        if not symbol or symbol == "-":
            continue
        try:
            vals = np.asarray([float(row[c]) for c in value_cols], dtype=float)
        except (KeyError, ValueError):
            continue
        mean = float(vals.mean())
        prev = best.get(symbol)
        if prev is None or mean > prev[0]:
            best[symbol] = (mean, vals)
    symbols = sorted(best)
    mat = np.vstack([best[s][1] for s in symbols]) if symbols else np.zeros((0, len(value_cols)))
    return symbols, mat


def log2fc(treat: np.ndarray, ctrl: np.ndarray, pseudo: float) -> np.ndarray:
    return np.log2(treat + pseudo).mean(axis=1) - np.log2(ctrl + pseudo).mean(axis=1)


def welch_gene(treat: np.ndarray, ctrl: np.ndarray, pseudo: float, gene_index: int):
    a = np.log2(treat[gene_index] + pseudo)
    b = np.log2(ctrl[gene_index] + pseudo)
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and np.allclose(a[0], b[0]):
        return float("nan")
    res = stats.ttest_ind(a, b, equal_var=False, alternative="two-sided")
    return float(res.pvalue)


def paired_gene(treat: np.ndarray, ctrl: np.ndarray, pseudo: float, gene_index: int):
    d = np.log2(treat[gene_index] + pseudo) - np.log2(ctrl[gene_index] + pseudo)
    if np.allclose(d, 0):
        return float("nan")
    res = stats.ttest_rel(np.log2(treat[gene_index] + pseudo), np.log2(ctrl[gene_index] + pseudo))
    return float(res.pvalue)


def set_sample_score(mat_arm: np.ndarray, idx: np.ndarray, pseudo: float) -> np.ndarray:
    """Mean log2(expr+pseudo) across genes, one number per sample."""
    if idx.size == 0:
        return np.full(mat_arm.shape[1], np.nan)
    return np.log2(mat_arm[idx] + pseudo).mean(axis=0)


def welch_vec(a: np.ndarray, b: np.ndarray) -> float:
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return float("nan")
    if np.allclose(a, a[0]) and np.allclose(b, b[0]) and np.allclose(a[0], b[0]):
        return float("nan")
    return float(stats.ttest_ind(a, b, equal_var=False).pvalue)


def fmt(x):
    if x is None or (isinstance(x, float) and not math.isfinite(x)):
        return "NA"
    if isinstance(x, float):
        return f"{x:.6g}"
    return str(x)


def write_tsv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def prepare_sets(hallmark, reactome):
    sets = {}
    for name in (HALLMARK_IFNG, HALLMARK_IFNA):
        sets[name] = hallmark[name]
    for name in (REACTOME_APM, REACTOME_NHEJ, REACTOME_IFNG):
        sets[name] = reactome[name]
    sets["CANONICAL_APM_18"] = set(CANONICAL_APM)
    sets["CLASSICAL_NHEJ_8"] = set(CLASSICAL_NHEJ)
    return sets


def evaluate_contrast(name, tissue, adc_only, symbols, mat, treat_idx, ctrl_idx, paired, sets, transforms):
    rows = []
    sym_index = {s: i for i, s in enumerate(symbols)}
    if "CLDN4" not in sym_index:
        raise SystemExit(f"CLDN4 missing from {name}")
    gi = sym_index["CLDN4"]
    treat = mat[:, treat_idx]
    ctrl = mat[:, ctrl_idx]
    for transform_name, pseudo, kind in transforms:
        fc = log2fc(treat, ctrl, pseudo)
        cldn4 = float(fc[gi])
        cldn4_p = welch_gene(treat, ctrl, pseudo, gi)
        cldn4_paired_p = paired_gene(treat, ctrl, pseudo, gi) if paired else float("nan")
        for bundle_name, bundle in SET_BUNDLES.items():
            members = {
                "ifn": sets[bundle["ifn_name"]],
                "apm": sets[bundle["apm_name"]],
                "nhej": sets[bundle["nhej_name"]],
            }
            present = {}
            for key, genes in members.items():
                present[key] = np.array([sym_index[g] for g in genes if g in sym_index], dtype=int)
            for summary_name in ("median", "mean"):
                stats_out = {}
                sample_p = {}
                for key, idx in present.items():
                    vals = fc[idx]
                    stats_out[key] = float(np.median(vals) if summary_name == "median" else np.mean(vals))
                    st = set_sample_score(treat, idx, pseudo)
                    sc = set_sample_score(ctrl, idx, pseudo)
                    sample_p[key] = welch_vec(st, sc)
                signed = (-cldn4) + stats_out["ifn"] + stats_out["apm"] + stats_out["nhej"]
                sign_hit = cldn4 < 0 and stats_out["ifn"] > 0 and stats_out["apm"] > 0 and stats_out["nhej"] > 0
                primary = (
                    bundle_name.startswith("primary_")
                    and summary_name == "median"
                    and transform_name in {"log2(x+1)", "log2(CPM+1)"}
                    and adc_only
                )
                rows.append({
                    "contrast": name,
                    "tissue": tissue,
                    "adc_monotherapy": "yes" if adc_only else "no",
                    "n_treat": int(treat.shape[1]),
                    "n_control": int(ctrl.shape[1]),
                    "paired": "yes" if paired else "no",
                    "transform": transform_name,
                    "expr_kind": kind,
                    "gene_sets": bundle_name,
                    "ifn_set": bundle["ifn_name"],
                    "apm_set": bundle["apm_name"],
                    "nhej_set": bundle["nhej_name"],
                    "summary": summary_name,
                    "n_ifn_genes": int(present["ifn"].size),
                    "n_apm_genes": int(present["apm"].size),
                    "n_nhej_genes": int(present["nhej"].size),
                    "CLDN4_log2FC": cldn4,
                    "CLDN4_welch_p": cldn4_p,
                    "CLDN4_paired_p": cldn4_paired_p,
                    "IFN_summary": stats_out["ifn"],
                    "APM_summary": stats_out["apm"],
                    "NHEJ_summary": stats_out["nhej"],
                    "IFN_sample_welch_p": sample_p["ifn"],
                    "APM_sample_welch_p": sample_p["apm"],
                    "NHEJ_sample_welch_p": sample_p["nhej"],
                    "signed_score": signed,
                    "sign_hit_CLDN4down_sets_up": "yes" if sign_hit else "no",
                    "primary_row": "yes" if primary else "no",
                    "pass_sign": "yes" if sign_hit else "no",
                    "pass_moderate_CLDN4lt-0.25_sets_gt0": "yes" if (cldn4 < -0.25 and sign_hit) else "no",
                    "pass_strict_CLDN4lt-0.5_sets_gt0.25": "yes" if (
                        cldn4 < -0.5 and stats_out["ifn"] > 0.25 and stats_out["apm"] > 0.25 and stats_out["nhej"] > 0.25
                    ) else "no",
                })
    return rows


def numeric_rows(rows):
    out = []
    for row in rows:
        r = dict(row)
        for k, v in list(r.items()):
            if isinstance(v, float):
                r[k] = fmt(v)
        out.append(r)
    return out


def plot_primary(rows):
    primary = [r for r in rows if r["primary_row"] == "yes"]
    # stable display order
    order = [
        "GSE312098_IMMU132_vs_control",
        "GSE311016_IMMU132_vs_control",
        "GSE304294_IMMU132_vs_control",
        "E-MTAB-16433_SG_vs_vehicle",
    ]
    primary = sorted(primary, key=lambda r: order.index(r["contrast"]) if r["contrast"] in order else 99)
    labels = [r["contrast"].replace("_", "\n") for r in primary]
    components = [
        ("-CLDN4 log2FC", [-r["CLDN4_log2FC"] for r in primary], "#b45309"),
        ("IFN median", [r["IFN_summary"] for r in primary], "#1d4e89"),
        ("APM median", [r["APM_summary"] for r in primary], "#0f766e"),
        ("NHEJ median", [r["NHEJ_summary"] for r in primary], "#6b21a8"),
    ]
    x = np.arange(len(primary))
    width = 0.18
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    for i, (lab, vals, color) in enumerate(components):
        ax.bar(x + (i - 1.5) * width, vals, width, label=lab, color=color)
    ax.axhline(0, color="#333", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel("Primary-method effect (log2)")
    ax.set_title("ADC monotherapy vs control\nlog2(x+1) or log2(CPM+1), median gene-set log2FC")
    ax.legend(fontsize=8, frameon=False, ncol=4, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGS / "sweep_primary_adc_only.png", dpi=160)
    fig.savefig(FIGS / "sweep_primary_adc_only.pdf")
    plt.close(fig)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312098/suppl/GSE312098_gene_fpkm.txt.gz",
        DATA / "GSE312098_gene_fpkm.txt.gz",
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE311nnn/GSE311016/suppl/GSE311016_gene_fpkm.txt.gz",
        DATA / "GSE311016_gene_fpkm.txt.gz",
    )
    download(
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304294/suppl/GSE304294_gene_fpkm.txt.gz",
        DATA / "GSE304294_gene_fpkm.txt.gz",
    )
    download(
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt",
        DATA / "h.all.v2024.1.Hs.symbols.gmt",
    )
    download(
        "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt",
        DATA / "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
    )
    download(
        "https://raw.githubusercontent.com/jinxuanhong1-blip/sdaxcge/cursor/trop2-adc-cldn4-extra-41de/methods/trop2_adc_cldn4_extra/tables/EMTAB16433_mouse_pseudobulk_counts.tsv.gz",
        DATA / "EMTAB16433_mouse_pseudobulk_counts.tsv.gz",
    )
    hallmark = load_gmt(DATA / "h.all.v2024.1.Hs.symbols.gmt")
    reactome = load_gmt(DATA / "c2.cp.reactome.v2024.1.Hs.symbols.gmt")
    sets = prepare_sets(hallmark, reactome)
    fpkm_transforms = [("log2(x+1)", 1.0, "FPKM"), ("log2(x+0.1)", 0.1, "FPKM")]
    count_transforms = [("log2(CPM+1)", 1.0, "CPM"), ("log2(CPM+0.1)", 0.1, "CPM"), ("log2(count+1)", 1.0, "count")]

    all_rows = []

    # GSE312098
    rows = read_table(DATA / "GSE312098_gene_fpkm.txt.gz")
    cols = [f"X_{i}" for i in range(1, 13)]
    symbols, mat = collapse_symbols(rows, cols, True)
    col = {c: i for i, c in enumerate(cols)}
    ctrl = [col[c] for c in ("X_1", "X_2", "X_3")]
    immu = [col[c] for c in ("X_4", "X_5", "X_6")]
    combo = [col[c] for c in ("X_10", "X_11", "X_12")]
    all_rows += evaluate_contrast(
        "GSE312098_IMMU132_vs_control", "CRC CX-1 cell line, 48 h", True,
        symbols, mat, immu, ctrl, False, sets, fpkm_transforms,
    )
    all_rows += evaluate_contrast(
        "GSE312098_combination_vs_control", "CRC CX-1, IMMU132+GSK2606414, 48 h", False,
        symbols, mat, combo, ctrl, False, sets, fpkm_transforms,
    )

    # GSE311016 paired PDX
    rows = read_table(DATA / "GSE311016_gene_fpkm.txt.gz")
    models = ["36", "82", "83", "114", "196"]
    cols = []
    for m in models:
        cols += [f"C_{m}", f"T_{m}"]
    symbols, mat = collapse_symbols(rows, cols, True)
    # rebuild in this column order
    col = {c: i for i, c in enumerate(cols)}
    ctrl = [col[f"C_{m}"] for m in models]
    immu = [col[f"T_{m}"] for m in models]
    all_rows += evaluate_contrast(
        "GSE311016_IMMU132_vs_control", "CRC PDX, day 29, paired by model", True,
        symbols, mat, immu, ctrl, True, sets, fpkm_transforms,
    )

    # GSE304294
    rows = read_table(DATA / "GSE304294_gene_fpkm.txt.gz")
    cols = ["OX1_1", "OX1_2", "OX1_3", "OX2_1", "OX2_2", "OX4_1", "OX4_2", "OX4_3"]
    symbols, mat = collapse_symbols(rows, cols, True)
    col = {c: i for i, c in enumerate(cols)}
    ctrl = [col[c] for c in ("OX1_1", "OX1_2", "OX1_3")]
    immu = [col[c] for c in ("OX2_1", "OX2_2")]
    combo = [col[c] for c in ("OX4_1", "OX4_2", "OX4_3")]
    all_rows += evaluate_contrast(
        "GSE304294_IMMU132_vs_control", "ESCC KYSE30, 1 day, n=2 vs 3", True,
        symbols, mat, immu, ctrl, False, sets, fpkm_transforms,
    )
    all_rows += evaluate_contrast(
        "GSE304294_combination_vs_control", "ESCC KYSE30, IMMU132+IACS010759, 1 day", False,
        symbols, mat, combo, ctrl, False, sets, fpkm_transforms,
    )

    # E-MTAB-16433 pseudobulk counts
    rows = read_table(DATA / "EMTAB16433_mouse_pseudobulk_counts.tsv.gz")
    # first column header is "barcode" but values are gene symbols
    treat_cols = [c for c in rows[0].keys() if c.startswith("trodelvy_")]
    ctrl_cols = [c for c in rows[0].keys() if c.startswith("vehicle_")]
    cols = treat_cols + ctrl_cols
    symbols, mat = collapse_symbols(rows, cols, False)
    # CPM
    lib = mat.sum(axis=0)
    cpm = mat / lib * 1e6
    treat_i = list(range(len(treat_cols)))
    ctrl_i = list(range(len(treat_cols), len(cols)))
    # log library sizes for the provenance note
    lib_note = {
        "trodelvy_total_counts": [float(x) for x in lib[treat_i]],
        "vehicle_total_counts": [float(x) for x in lib[ctrl_i]],
    }
    all_rows += evaluate_contrast(
        "E-MTAB-16433_SG_vs_vehicle", "CRC PDOX HD42466, 28 d, mouse pseudobulk", True,
        symbols, cpm, treat_i, ctrl_i, False, sets,
        [("log2(CPM+1)", 1.0, "CPM"), ("log2(CPM+0.1)", 0.1, "CPM")],
    )
    all_rows += evaluate_contrast(
        "E-MTAB-16433_SG_vs_vehicle", "CRC PDOX HD42466, 28 d, mouse pseudobulk", True,
        symbols, mat, treat_i, ctrl_i, False, sets,
        [("log2(count+1)", 1.0, "count")],
    )

    # stringify
    raw_rows = all_rows
    out_rows = numeric_rows(all_rows)
    write_tsv(TABLES / "sweep_all_rows.tsv", out_rows)
    write_tsv(TABLES / "sweep_primary.tsv", [r for r in out_rows if r["primary_row"] == "yes"])
    sign_hits = [r for r in out_rows if r["sign_hit_CLDN4down_sets_up"] == "yes"]
    sign_hits = sorted(sign_hits, key=lambda r: -float(r["signed_score"]))
    write_tsv(TABLES / "sweep_sign_hits.tsv", sign_hits)

    # rank among ADC-only primary-transform rows and among all sign hits
    def score_of(r):
        return r["signed_score"]

    adc_sign = [r for r in raw_rows if r["adc_monotherapy"] == "yes" and r["sign_hit_CLDN4down_sets_up"] == "yes"]
    winner = max(adc_sign, key=score_of) if adc_sign else None
    primary = [r for r in raw_rows if r["primary_row"] == "yes"]
    primary_winner = max(primary, key=score_of) if primary else None
    strict = [r for r in raw_rows if r["pass_strict_CLDN4lt-0.5_sets_gt0.25"] == "yes"]
    moderate = [r for r in raw_rows if r["pass_moderate_CLDN4lt-0.25_sets_gt0"] == "yes" and r["adc_monotherapy"] == "yes"]

    summary = []
    summary.append({"item": "n_grid_rows", "value": str(len(raw_rows))})
    summary.append({"item": "n_sign_hits_all_rows", "value": str(sum(r["sign_hit_CLDN4down_sets_up"] == "yes" for r in raw_rows))})
    summary.append({"item": "n_adc_only_sign_hits", "value": str(len(adc_sign))})
    summary.append({"item": "n_adc_only_moderate_threshold", "value": str(len(moderate))})
    summary.append({"item": "n_strict_threshold_any_row", "value": str(len(strict))})
    summary.append({"item": "emtab16433_library_totals", "value": str(lib_note)})
    if winner:
        summary.append({"item": "strongest_adc_only_sign_hit", "value": winner["contrast"]})
        summary.append({"item": "strongest_transform", "value": winner["transform"]})
        summary.append({"item": "strongest_sets", "value": winner["gene_sets"]})
        summary.append({"item": "strongest_summary", "value": winner["summary"]})
        summary.append({"item": "strongest_signed_score", "value": fmt(winner["signed_score"])})
        summary.append({"item": "strongest_CLDN4", "value": fmt(winner["CLDN4_log2FC"])})
        summary.append({"item": "strongest_IFN", "value": fmt(winner["IFN_summary"])})
        summary.append({"item": "strongest_APM", "value": fmt(winner["APM_summary"])})
        summary.append({"item": "strongest_NHEJ", "value": fmt(winner["NHEJ_summary"])})
    if primary_winner:
        summary.append({"item": "strongest_primary_contrast", "value": primary_winner["contrast"]})
        summary.append({"item": "strongest_primary_score", "value": fmt(primary_winner["signed_score"])})
        summary.append({"item": "strongest_primary_sign_hit", "value": primary_winner["sign_hit_CLDN4down_sets_up"]})
    write_tsv(TABLES / "sweep_summary.tsv", summary)
    plot_primary(raw_rows)

    print("rows", len(raw_rows), "sign hits", len(sign_hits), "strict", len(strict))
    print("PRIMARY")
    for r in primary:
        print(r["contrast"], "CLDN4", fmt(r["CLDN4_log2FC"]), "p", fmt(r["CLDN4_welch_p"]),
              "IFN", fmt(r["IFN_summary"]), "APM", fmt(r["APM_summary"]), "NHEJ", fmt(r["NHEJ_summary"]),
              "score", fmt(r["signed_score"]), "hit", r["sign_hit_CLDN4down_sets_up"],
              "n", r["n_ifn_genes"], r["n_apm_genes"], r["n_nhej_genes"])
    if winner:
        print("WINNER", winner["contrast"], winner["transform"], winner["gene_sets"], winner["summary"], fmt(winner["signed_score"]))
    print("library", lib_note)


if __name__ == "__main__":
    main()
