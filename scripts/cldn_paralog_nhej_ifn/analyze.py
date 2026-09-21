#!/usr/bin/env python3
"""Score NHEJ and IFN modules on open CLDN1/3/7 knockdown RNA profiles vs CLDN4.

Effect size is the difference in mean log2 expression (KD minus control).
The replicate-level test is a Welch t-test on the per-sample module mean.
A joint NHEJ-down / IFN-up call requires the mean and the median gene effect
to share that sign, the on-target claudin to fall (log2 difference < -0.5),
and at least 3 samples in each arm. Nominal p < 0.05 on the two pre-specified
modules (NHEJ_CORE, HALLMARK interferon-alpha) is recorded; it is not a
genome-wide DESeq2 FDR.

Unreplicated matrices (group-mean FPKM, author edgeR with one GEO library per
arm) contribute direction only.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[2]
DATA = Path("/tmp/cldn_kd")
OUT = ROOT / "results" / "cldn_paralog_nhej_ifn"
SETS_PATH = Path(__file__).resolve().parent / "gene_sets.tsv"

PRIMARY = ("NHEJ_CORE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")
MODULES = (
    "NHEJ_CORE",
    "NHEJ_53BP1_SHIELDIN",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "HALLMARK_DNA_REPAIR",
)
SHORT = {
    "NHEJ_CORE": "NHEJ_CORE",
    "NHEJ_53BP1_SHIELDIN": "NHEJ_53BP1",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE": "IFN_ALPHA",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": "IFN_GAMMA",
    "HALLMARK_DNA_REPAIR": "DNA_REPAIR",
}


def load_sets() -> dict[str, set[str]]:
    df = pd.read_csv(SETS_PATH, sep="\t")
    out = {}
    for name, sub in df.groupby("set"):
        out[name] = set(sub["gene"].str.upper())
    missing = [m for m in MODULES if m not in out]
    if missing:
        raise SystemExit(f"gene sets missing: {missing}")
    return out


def collapse_upper(df: pd.DataFrame) -> pd.DataFrame:
    """Genes x samples. Index uppercased; duplicate symbols summed."""
    x = df.copy()
    x.index = x.index.astype(str).str.upper()
    x = x[~x.index.isin(["", "NAN", "NONE", "-"])]
    x = x.groupby(level=0).sum(numeric_only=True)
    return x


def read_series_matrix(path: Path) -> tuple[pd.DataFrame, list[str]]:
    titles = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("!Sample_title"):
                titles = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!series_matrix_table_begin"):
                break
        df = pd.read_csv(fh, sep="\t", index_col=0)
    df = df.loc[~df.index.astype(str).str.startswith("!")]
    df.index = df.index.astype(str).str.strip('"')
    df.columns = [c.strip('"') for c in df.columns]
    return df.apply(pd.to_numeric, errors="coerce"), titles


def log2cpm(counts: pd.DataFrame) -> pd.DataFrame:
    lib = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(lib, axis=1) * 1e6
    return np.log2(cpm + 1.0)


def welch(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2:
        return np.nan
    if np.nanstd(a) == 0 and np.nanstd(b) == 0:
        return 1.0 if np.allclose(a.mean(), b.mean()) else 0.0
    return float(stats.ttest_ind(a, b, equal_var=False, alternative="two-sided").pvalue)


def score_modules(
    log_expr: pd.DataFrame,
    kd_cols: list[str],
    ctrl_cols: list[str],
    sets: dict[str, set[str]],
    *,
    ratio_mode: bool = False,
) -> tuple[pd.DataFrame, dict]:
    """log_expr is genes x samples, already log2 (or log2 ratio if ratio_mode)."""
    use = log_expr.copy()
    rows = []
    gene_delta = None
    if ratio_mode:
        # each column is already a KD/reference log2 ratio
        gene_delta = use[kd_cols].mean(axis=1)
        n_kd, n_ctrl = len(kd_cols), 0
    else:
        gene_delta = use[kd_cols].mean(axis=1) - use[ctrl_cols].mean(axis=1)
        n_kd, n_ctrl = len(kd_cols), len(ctrl_cols)
    background = gene_delta.replace([np.inf, -np.inf], np.nan).dropna()
    for name in MODULES:
        members = [g for g in use.index if g in sets[name]]
        sub = gene_delta.reindex(members).replace([np.inf, -np.inf], np.nan).dropna()
        if ratio_mode:
            sample_scores = use.reindex(members).mean(axis=0).reindex(kd_cols)
            p_sample = (
                float(stats.ttest_1samp(sample_scores.dropna(), 0.0).pvalue)
                if sample_scores.dropna().shape[0] >= 2
                else np.nan
            )
            delta = float(sample_scores.mean()) if len(sample_scores) else np.nan
        else:
            sample_scores = use.reindex(members).mean(axis=0)
            delta = float(sample_scores[kd_cols].mean() - sample_scores[ctrl_cols].mean())
            p_sample = welch(sample_scores[kd_cols].to_numpy(), sample_scores[ctrl_cols].to_numpy())
        rest = background.drop(index=sub.index, errors="ignore")
        if len(sub) >= 5 and len(rest) >= 20:
            p_mw = float(stats.mannwhitneyu(sub, rest, alternative="two-sided").pvalue)
        else:
            p_mw = np.nan
        if len(sub) >= 5 and np.any(sub.to_numpy() != 0):
            try:
                p_wil = float(stats.wilcoxon(sub.to_numpy(), alternative="two-sided").pvalue)
            except ValueError:
                p_wil = np.nan
        else:
            p_wil = np.nan
        rows.append(
            {
                "module": name,
                "module_short": SHORT[name],
                "n_genes": int(len(sub)),
                "delta_mean": delta,
                "delta_median": float(sub.median()) if len(sub) else np.nan,
                "frac_genes_down": float((sub < 0).mean()) if len(sub) else np.nan,
                "frac_genes_up": float((sub > 0).mean()) if len(sub) else np.nan,
                "p_sample": p_sample,
                "p_gene_mw": p_mw,
                "p_gene_wilcoxon": p_wil,
                "n_kd": n_kd,
                "n_ctrl": n_ctrl,
            }
        )
    return pd.DataFrame(rows), {"gene_delta": gene_delta}


def target_stats(log_expr, kd_cols, ctrl_cols, target: str, ratio_mode=False) -> dict:
    target = target.upper()
    if target not in log_expr.index:
        return {"target": target, "target_log2fc": np.nan, "target_p": np.nan, "target_present": False}
    if ratio_mode:
        vals = log_expr.loc[target, kd_cols].astype(float)
        p = float(stats.ttest_1samp(vals, 0.0).pvalue) if vals.shape[0] >= 2 else np.nan
        return {
            "target": target,
            "target_log2fc": float(vals.mean()),
            "target_p": p,
            "target_present": True,
        }
    kd = log_expr.loc[target, kd_cols].astype(float)
    ct = log_expr.loc[target, ctrl_cols].astype(float)
    return {
        "target": target,
        "target_log2fc": float(kd.mean() - ct.mean()),
        "target_p": welch(kd.to_numpy(), ct.to_numpy()),
        "target_present": True,
    }


def arm_down(mean: float, median: float) -> bool:
    """A module is down only when the mean and the median both move by a real amount."""
    return bool(np.isfinite(mean) and np.isfinite(median) and mean <= -0.10 and median <= -0.05)


def arm_up(mean: float, median: float) -> bool:
    return bool(np.isfinite(mean) and np.isfinite(median) and mean >= 0.10 and median >= 0.05)


def annotate_call(mod: pd.DataFrame, tgt: dict, *, allow_nominal: bool) -> str:
    if not tgt.get("target_present") or not np.isfinite(tgt.get("target_log2fc", np.nan)):
        return "invalid_target_absent"
    if tgt["target_log2fc"] >= -0.5:
        return "invalid_target_not_down"
    nhej = mod.loc[mod["module"] == "NHEJ_CORE"].iloc[0]
    ifn = mod.loc[mod["module"] == "HALLMARK_INTERFERON_ALPHA_RESPONSE"].iloc[0]
    nhej_down = arm_down(nhej["delta_mean"], nhej["delta_median"])
    ifn_up = arm_up(ifn["delta_mean"], ifn["delta_median"])
    if not (nhej_down and ifn_up):
        return "joint_pattern_absent"
    if not allow_nominal:
        return "direction_only"
    if (
        nhej["n_kd"] >= 3
        and nhej["n_ctrl"] >= 3
        and nhej["p_sample"] < 0.05
        and ifn["p_sample"] < 0.05
    ):
        return "joint_nominal_p"
    return "direction_only"


def bh_within(mod: pd.DataFrame) -> pd.DataFrame:
    out = mod.copy()
    p = out["p_sample"].to_numpy(dtype=float)
    q = np.full(len(out), np.nan)
    ok = np.isfinite(p)
    if ok.sum() >= 1:
        q[ok] = multipletests(p[ok], method="fdr_bh")[1]
    out["q_sample_bh"] = q
    return out


def human_symbol_map_from_counts(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", usecols=["gene_id", "gene_name"])
    df["gene_id"] = df["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    df = df.dropna().drop_duplicates("gene_id")
    return dict(zip(df["gene_id"], df["gene_name"].astype(str)))


def load_hucc_counts(path: Path, kd_prefix: str, ctrl_prefix: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    df = pd.read_csv(path, sep="\t")
    meta_cols = [c for c in df.columns if not str(c).startswith(("A1_", "CDX_", "PDM"))]
    # keep count columns by numeric dtype after setting index
    df = df.set_index("gene_name")
    count_cols = [c for c in df.columns if c not in meta_cols and c != "gene_id"]
    # meta columns that are not samples
    sample_cols = [c for c in count_cols if pd.api.types.is_numeric_dtype(df[c])]
    counts = collapse_upper(df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0))
    kd = [c for c in counts.columns if kd_prefix in c]
    ctrl = [c for c in counts.columns if ctrl_prefix in c]
    return counts, kd, ctrl


def gpl6104_map() -> dict[str, str]:
    path = DATA / "annot" / "GPL6104.annot.gz"
    mapping = {}
    with gzip.open(path, "rt") as fh:
        header = None
        for line in fh:
            if line.startswith("ID\t"):
                header = line.rstrip("\n").split("\t")
                break
        assert header is not None
        id_i = header.index("ID")
        sy_i = header.index("Gene symbol")
        for line in fh:
            if not line or line.startswith(("!", "#", "^")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(id_i, sy_i):
                continue
            sym = parts[sy_i].strip()
            if sym and sym != "NA":
                mapping[parts[id_i].strip()] = sym
    return mapping


def gpl10555_map() -> dict[str, str]:
    path = DATA / "annot" / "GPL10555_family.soft.gz"
    mapping = {}
    with gzip.open(path, "rt") as fh:
        started = False
        for line in fh:
            if line.startswith("!platform_table_begin"):
                started = True
                header = next(fh).rstrip("\n").split("\t")
                id_i = header.index("ID")
                orf_i = header.index("ORF")
                continue
            if not started:
                continue
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= orf_i:
                continue
            sym = parts[orf_i].strip()
            if sym:
                mapping[parts[id_i].strip()] = sym
    return mapping


def ensure_gpl10787() -> None:
    """Stream the GEO platform table and keep probe id plus gene symbol."""
    dest = DATA / "annot" / "GPL10787_symbols.tsv"
    if dest.exists() and dest.stat().st_size > 0:
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    import urllib.request

    url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10787&targ=gpl&form=text&view=data"
    req = urllib.request.Request(url, headers={"User-Agent": "research/1.0"})
    with urllib.request.urlopen(req, timeout=600) as resp, dest.open("w") as out:
        out.write("id\tsymbol\n")
        started = False
        for raw in resp:
            line = raw.decode("utf-8", "replace").rstrip("\n")
            if line.startswith("!platform_table_begin"):
                started = True
                continue
            if not started:
                continue
            if line.startswith("!platform_table_end"):
                break
            if line.startswith("ID\t"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            sym = parts[6].strip()
            if sym:
                out.write(parts[0].strip() + "\t" + sym + "\n")


def gpl10787_map() -> dict[str, str]:
    ensure_gpl10787()
    df = pd.read_csv(DATA / "annot" / "GPL10787_symbols.tsv", sep="\t", dtype=str)
    df = df.dropna()
    df = df[df["symbol"].str.len() > 0]
    return dict(zip(df["id"], df["symbol"]))


def collapse_probes(probe_df: pd.DataFrame, mapping: dict[str, str], log: bool) -> pd.DataFrame:
    x = probe_df.copy()
    x["symbol"] = x.index.map(lambda i: mapping.get(str(i), ""))
    x = x[x["symbol"].astype(str).str.len() > 0]
    val_cols = [c for c in x.columns if c != "symbol"]
    if log:
        x[val_cols] = np.log2(x[val_cols].clip(lower=0) + 1.0)
    # median across probes
    g = x.groupby(x["symbol"].str.upper())[val_cols].median()
    return g


def mouse_ensembl_map() -> dict[str, str]:
    df = pd.read_csv(DATA / "GSE274940_raw_counts.csv.gz")
    df = df.dropna(subset=["ENSEMBL", "ALIAS"])
    df["ENSEMBL"] = df["ENSEMBL"].astype(str)
    df["ALIAS"] = df["ALIAS"].astype(str)
    return dict(zip(df["ENSEMBL"], df["ALIAS"]))


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    sets = load_sets()
    module_rows = []
    gene_rows = []
    summary_rows = []

    def add(contrast_id, gene, perturbation, accession, species, tissue, assay, contrast_class, allow_nominal, log_expr, kd, ctrl, ratio_mode=False, note=""):
        mod, extra = score_modules(log_expr, kd, ctrl, sets, ratio_mode=ratio_mode)
        mod = bh_within(mod)
        tgt = target_stats(log_expr, kd, ctrl, gene, ratio_mode=ratio_mode)
        call = annotate_call(mod, tgt, allow_nominal=allow_nominal and not ratio_mode)
        mod.insert(0, "contrast_id", contrast_id)
        mod.insert(1, "accession", accession)
        mod.insert(2, "gene", gene)
        module_rows.append(mod)
        # focus genes
        gd = extra["gene_delta"]
        focus = sorted(sets["NHEJ_CORE"] | sets["NHEJ_53BP1_SHIELDIN"] | {gene.upper()})
        for g in focus:
            if g in gd.index and np.isfinite(gd.loc[g]):
                gene_rows.append(
                    {
                        "contrast_id": contrast_id,
                        "gene_symbol": g,
                        "log2_diff": float(gd.loc[g]),
                        "in_nhej_core": g in sets["NHEJ_CORE"],
                        "in_shieldin": g in sets["NHEJ_53BP1_SHIELDIN"],
                        "is_target": g == gene.upper(),
                    }
                )
        nhej = mod.loc[mod.module == "NHEJ_CORE"].iloc[0]
        ifn = mod.loc[mod.module == "HALLMARK_INTERFERON_ALPHA_RESPONSE"].iloc[0]
        gfn = mod.loc[mod.module == "HALLMARK_INTERFERON_GAMMA_RESPONSE"].iloc[0]
        dna = mod.loc[mod.module == "HALLMARK_DNA_REPAIR"].iloc[0]
        summary_rows.append(
            {
                "contrast_id": contrast_id,
                "accession": accession,
                "gene": gene,
                "perturbation": perturbation,
                "species": species,
                "tissue_or_model": tissue,
                "assay": assay,
                "contrast_class": contrast_class,
                "n_kd": int(nhej["n_kd"]),
                "n_ctrl": int(nhej["n_ctrl"]),
                "target_log2fc": tgt["target_log2fc"],
                "target_p_sample": tgt["target_p"],
                "nhej_delta": nhej["delta_mean"],
                "nhej_median": nhej["delta_median"],
                "nhej_p": nhej["p_sample"],
                "nhej_q": nhej["q_sample_bh"],
                "nhej_n_genes": int(nhej["n_genes"]),
                "ifn_alpha_delta": ifn["delta_mean"],
                "ifn_alpha_median": ifn["delta_median"],
                "ifn_alpha_p": ifn["p_sample"],
                "ifn_alpha_q": ifn["q_sample_bh"],
                "ifn_alpha_n_genes": int(ifn["n_genes"]),
                "ifn_gamma_delta": gfn["delta_mean"],
                "ifn_gamma_median": gfn["delta_median"],
                "ifn_gamma_p": gfn["p_sample"],
                "dna_repair_delta": dna["delta_mean"],
                "dna_repair_median": dna["delta_median"],
                "dna_repair_p": dna["p_sample"],
                "pattern": call,
                "note": note,
            }
        )
        print(
            f"{contrast_id}: target {tgt['target_log2fc']:.3f}  "
            f"NHEJ {nhej['delta_mean']:.3f} (p={nhej['p_sample']})  "
            f"IFNa {ifn['delta_mean']:.3f} (p={ifn['p_sample']})  {call}",
            flush=True,
        )

    # ----- CLDN1 RNA-seq counts -----
    specs = [
        (
            "GSE312708_HuCCA1_CLDN1KO_vs_ctrl",
            DATA / "GSE312708_gene_count.tsv.gz",
            "KO",
            "CTRL",
            "HuCC-A1 CLDN1 KO vs control, in vitro",
            "genetic_ko_rnaseq",
            True,
            "CLDN1 raw counts stay high (control 58k–71k, KO 43k–68k). mRNA loss is not confirmed.",
        ),
        (
            "GSE312711_HuCCA1_CDX_CLDN1KO_vs_ctrl",
            DATA / "GSE312711_gene_count.tsv.gz",
            "KO",
            "CT",
            "HuCC-A1 CLDN1 KO vs control, cell-line xenograft",
            "genetic_ko_rnaseq",
            True,
            "CLDN1 counts remain high in the knockout tumors (about 51k–198k versus 102k–264k in controls).",
        ),
        (
            "GSE312713_PDM217_shCLDN1_vs_ctrl",
            DATA / "GSE312713_gene_count.tsv.gz",
            "sh",
            "CTL",
            "CCA organoid PDM217 shCLDN1 vs control",
            "genetic_kd_rnaseq",
            True,
            "3 vs 3",
        ),
    ]
    for cid, path, kd_key, ctrl_key, tissue, klass, allow, note in specs:
        counts, kd, ctrl = load_hucc_counts(path, kd_key, ctrl_key)
        print(cid, "samples", kd, ctrl, flush=True)
        if not kd or not ctrl:
            raise SystemExit(f"sample parse failed for {cid}")
        add(cid, "CLDN1", "KO" if "KO" in cid else "shRNA", cid.split("_")[0], "human", tissue, "RNA-seq counts", klass, allow, log2cpm(counts), kd, ctrl, note=note)

    # ----- GSE234513 TPM, column order checked against CLDN1 -----
    sym = human_symbol_map_from_counts(DATA / "GSE312708_gene_count.tsv.gz")
    tpm = pd.read_csv(DATA / "GSE234513_tpms.csv.gz")
    tpm["gene_id"] = tpm["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)
    tpm["symbol"] = tpm["gene_id"].map(sym)
    tpm = tpm.dropna(subset=["symbol"])
    val_cols = [c for c in tpm.columns if c.startswith("Yasu_")]
    # natural order Yasu_1 .. Yasu_12
    val_cols = sorted(val_cols, key=lambda s: int(s.split("_")[1]))
    mat = tpm.groupby(tpm["symbol"].str.upper())[val_cols].max()
    logtpm = np.log2(mat + 1.0)
    # GEO GSM order is the only column key deposited (Yasu_1..12).
    # Validate: CLDN1 should be lower in the shCLDN1 blocks than in shNTC.
    cldn1 = logtpm.loc["CLDN1", val_cols]
    print("GSE234513 CLDN1 log2(TPM+1):", {c: round(float(cldn1[c]), 3) for c in val_cols}, flush=True)
    # Expected GSM order: 1-2 DLD NTC, 3-6 DLD sh, 7-8 LoVo NTC, 9-12 LoVo sh
    dld_ntc, dld_sh = val_cols[0:2], val_cols[2:6]
    lovo_ntc, lovo_sh = val_cols[6:8], val_cols[8:12]
    # GEO sample order is DLD-1 shNTC, shCLDN1#1, shCLDN1#2, then the same for LoVo,
    # two libraries each. CLDN1 falls in the first hairpin block of each line and
    # does not fall in the second, so the hairpins are scored separately.
    blocks = {
        "GSE234513_DLD1_sh1_vs_shNTC": (val_cols[2:4], val_cols[0:2], "DLD-1"),
        "GSE234513_DLD1_sh2_vs_shNTC": (val_cols[4:6], val_cols[0:2], "DLD-1"),
        "GSE234513_LoVo_sh1_vs_shNTC": (val_cols[8:10], val_cols[6:8], "LoVo"),
        "GSE234513_LoVo_sh2_vs_shNTC": (val_cols[10:12], val_cols[6:8], "LoVo"),
    }
    for cid, (kd_cols, ctrl_cols, line) in blocks.items():
        add(
            cid,
            "CLDN1",
            "shRNA liver tumors, one hairpin",
            "GSE234513",
            "human",
            f"{line} LIN28B-high liver tumors after portal-vein injection",
            "RNA-seq TPM",
            "genetic_kd_rnaseq",
            False,
            logtpm,
            kd_cols,
            ctrl_cols,
            note=(
                "Yasu_1..12 assigned in GEO sample order (shNTC, sh1, sh2 within each line). "
                "n=2 vs 2, so the call cannot be joint_nominal. "
                "A hairpin is interpretable only when CLDN1 log2(TPM+1) falls by more than 0.5."
            ),
        )

    # ----- GSE296175 siCLDN1 and PDS mAb -----
    gist = pd.read_excel(DATA / "GSE296175.xlsx", sheet_name="counts")
    gist["symbol"] = gist["Transcript"].astype(str).str.replace(r"-\d+$", "", regex=True)
    gcols = [c for c in gist.columns if c not in ("Transcript", "symbol")]
    gcounts = gist.groupby("symbol")[gcols].sum()
    gcounts = collapse_upper(gcounts)
    si = [c for c in gcounts.columns if c.startswith("siCLDN1")]
    unt = [c for c in gcounts.columns if c.startswith("GIST-T1R")]
    pds = [c for c in gcounts.columns if c.startswith("PDS")]
    print("GIST columns", si, unt, pds, flush=True)
    add(
        "GSE296175_siCLDN1_vs_untreated",
        "CLDN1",
        "siRNA",
        "GSE296175",
        "human",
        "imatinib-resistant GIST cells (deposited titles mix GIST-T1R and GIST 430)",
        "RNA-seq counts",
        "genetic_kd_rnaseq",
        True,
        log2cpm(gcounts),
        si,
        unt,
        note="3 vs 3. PDS-0330 arm excluded from this contrast.",
    )
    add(
        "GSE296175_PDS0330_vs_untreated",
        "CLDN1",
        "CLDN1 monoclonal antibody PDS-0330, not genetic KD",
        "GSE296175",
        "human",
        "same GIST cultures, antibody arm",
        "RNA-seq counts",
        "mab_not_genetic",
        False,
        log2cpm(gcounts),
        pds,
        unt,
        note="Scored only as a non-genetic comparator. Not a knockdown.",
    )

    # ----- CLDN3 germline liver -----
    mmap = mouse_ensembl_map()
    liver = pd.read_csv(DATA / "GSE159914_data.csv.gz", index_col=0)
    liver.index = liver.index.map(lambda i: mmap.get(str(i), ""))
    liver = liver[liver.index.astype(str).str.len() > 0]
    liver = collapse_upper(liver.apply(pd.to_numeric, errors="coerce").fillna(0))
    print("GSE159914 mapped genes", liver.shape[0], "Cldn3" in liver.index or "CLDN3" in liver.index, flush=True)
    base_ko = ["CTRL_1", "CTRL_2", "CTRL_3"]
    base_wt = ["CTRL_5", "CTRL_6", "CTRL_7"]
    ph_ko = ["PH_126", "PH_127", "PH_128"]
    ph_wt = ["PH_123", "PH_124", "PH_131"]
    aged_ko = ["Aged_KO_1", "Aged_KO_2", "Aged_KO_3"]
    aged_wt = ["Aged_WT_1", "Aged_WT_2", "Aged_WT_3"]
    add(
        "GSE159914_Cldn3KO_vs_WT_baseline_liver",
        "CLDN3",
        "germline knockout",
        "GSE159914",
        "mouse",
        "adult female liver, no surgery",
        "RNA-seq counts",
        "genetic_ko_rnaseq",
        True,
        log2cpm(liver),
        base_ko,
        base_wt,
        note="CTRL_1/2/3 are CLDN3-/-; CTRL_5/6/7 are CLDN3+/+. Technical replicate CTRL5_Tec2 excluded.",
    )
    add(
        "GSE159914_Cldn3KO_vs_WT_48h_hepatectomy",
        "CLDN3",
        "germline knockout, regenerating liver",
        "GSE159914",
        "mouse",
        "liver 48 h after two-thirds hepatectomy",
        "RNA-seq counts",
        "genetic_ko_rnaseq",
        True,
        log2cpm(liver),
        ph_ko,
        ph_wt,
        note="Regeneration is a second factor. Not a cancer knockdown.",
    )
    add(
        "GSE159914_Cldn3KO_vs_WT_aged_liver",
        "CLDN3",
        "germline knockout, aged",
        "GSE159914",
        "mouse",
        "aged female liver, baseline",
        "RNA-seq counts",
        "genetic_ko_rnaseq",
        True,
        log2cpm(liver),
        aged_ko,
        aged_wt,
        note="Separate aged cohort in the same series.",
    )

    # ----- all-claudin null, not paralog-specific -----
    null = pd.read_csv(DATA / "GSE274940_raw_counts.csv.gz")
    null = null.set_index("ALIAS")
    null_cols = ["WT1", "WT2", "WT3", "KO1", "KO2", "KO3"]
    null_counts = collapse_upper(null[null_cols].apply(pd.to_numeric, errors="coerce").fillna(0))
    # no single target; use CLDN3 as a representative that should fall, but call rule needs one gene.
    # Score with a synthetic check: report CLDN4 as the tracked gene AND override pattern if we want.
    add(
        "GSE274940_EpH4_allCldnNull_vs_WT",
        "CLDN4",
        "deletion of the entire claudin family, not CLDN4 alone",
        "GSE274940",
        "mouse",
        "EpH4 mammary epithelial cells",
        "RNA-seq counts",
        "all_claudin_null",
        False,
        log2cpm(null_counts),
        ["KO1", "KO2", "KO3"],
        ["WT1", "WT2", "WT3"],
        note="Comparator only. Every claudin is gone, so a CLDN4-specific claim cannot be read from it.",
    )

    # ----- CLDN7 shRNA organoids -----
    fc = pd.read_csv(DATA / "GSE273512_fc.txt.gz", sep="\t")
    fc = fc.dropna(subset=["GeneID"])
    cols = [c for c in fc.columns if c not in ("EnsmblID", "GeneID")]
    fc_counts = collapse_upper(fc.set_index("GeneID")[cols].apply(pd.to_numeric, errors="coerce").fillna(0))
    ctrl = [c for c in fc_counts.columns if c.startswith("Control")]
    sh1 = [c for c in fc_counts.columns if "sh1" in c]
    sh2 = [c for c in fc_counts.columns if "sh2" in c]
    sh_all = sh1 + sh2
    print("Cldn7 samples", ctrl, sh1, sh2, flush=True)
    log_fc = log2cpm(fc_counts)
    add(
        "GSE273512_shCldn7_pooled_vs_shControl",
        "CLDN7",
        "shRNA, two hairpins pooled",
        "GSE273512",
        "mouse",
        "MMTV-Neu mammary organoids in collagen",
        "RNA-seq counts",
        "genetic_kd_rnaseq",
        True,
        log_fc,
        sh_all,
        ctrl,
        note="Pooled hairpins. sh1 and sh2 are also scored separately for concordance.",
    )
    add(
        "GSE273512_sh1Cldn7_vs_shControl",
        "CLDN7",
        "shRNA hairpin 1",
        "GSE273512",
        "mouse",
        "MMTV-Neu mammary organoids in collagen",
        "RNA-seq counts",
        "genetic_kd_rnaseq",
        True,
        log_fc,
        sh1,
        ctrl,
        note="4 vs 4",
    )
    add(
        "GSE273512_sh2Cldn7_vs_shControl",
        "CLDN7",
        "shRNA hairpin 2",
        "GSE273512",
        "mouse",
        "MMTV-Neu mammary organoids in collagen",
        "RNA-seq counts",
        "genetic_kd_rnaseq",
        True,
        log_fc,
        sh2,
        ctrl,
        note="4 vs 4",
    )

    # ----- CLDN7 siRNA array -----
    raw = pd.read_csv(DATA / "GSE26055_non-normalized.txt.gz", sep="\t")
    raw = raw.set_index("ID")
    sig_cols = [c for c in raw.columns if c.startswith("AVG_Signal")]
    probe = raw[sig_cols].apply(pd.to_numeric, errors="coerce")
    probe.columns = [c.replace("AVG_Signal-", "") for c in probe.columns]
    mapped = collapse_probes(probe, gpl6104_map(), log=True)
    print("GSE26055 genes", mapped.shape, "cols", list(mapped.columns), flush=True)
    # column names from the file
    def cols_containing(df, *bits):
        hit = []
        for c in df.columns:
            if all(b.lower() in c.lower() for b in bits):
                hit.append(c)
        return hit

    for line, token in (("OVCA420", "420"), ("OVCAR2", "R2")):
        kd = cols_containing(mapped, token, "cldn7")
        ct = cols_containing(mapped, token, "cont")
        print(line, kd, ct, flush=True)
        add(
            f"GSE26055_{line}_siCLDN7_vs_ctrl",
            "CLDN7",
            "siRNA",
            "GSE26055",
            "human",
            f"{line} ovarian cancer cells",
            "Illumina HumanRef-8 array",
            "genetic_kd_array",
            False,
            mapped,
            kd,
            ct,
            note="n=2 vs 2. Series matrix values are empty; non-normalized AVG_Signal was log2-transformed. Call capped at direction_only.",
        )

    # ----- CLDN7 germline intestine array -----
    sm, _titles = read_series_matrix(DATA / "GSE256329_series_matrix.txt.gz")
    intest = collapse_probes(sm, gpl10787_map(), log=True)
    # GSM order from the series matrix header
    # WT SI, KO SI, WT LI, KO LI, WT SI, KO SI, WT LI, KO LI
    gsm = list(intest.columns)
    print("GSE256329 gsm", gsm, "genes", intest.shape[0], "CLDN7" in intest.index, flush=True)
    add(
        "GSE256329_Cldn7KO_vs_WT_small_intestine",
        "CLDN7",
        "germline knockout",
        "GSE256329",
        "mouse",
        "postnatal day-3 small intestine",
        "Agilent SurePrint G3 mouse array",
        "genetic_ko_array",
        False,
        intest,
        [gsm[1], gsm[5]],
        [gsm[0], gsm[4]],
        note="n=2 vs 2 whole-tissue array. Not a cancer knockdown.",
    )
    add(
        "GSE256329_Cldn7KO_vs_WT_large_intestine",
        "CLDN7",
        "germline knockout",
        "GSE256329",
        "mouse",
        "postnatal day-3 large intestine",
        "Agilent SurePrint G3 mouse array",
        "genetic_ko_array",
        False,
        intest,
        [gsm[3], gsm[7]],
        [gsm[2], gsm[6]],
        note="n=2 vs 2 whole-tissue array. Not a cancer knockdown.",
    )

    # ----- CLDN4 public contrast -----
    fpkm = pd.read_csv(DATA / "GSE207704_CLDN4_RNAseq.txt.gz", sep="\t")
    fpkm = fpkm.dropna(subset=["gene_short_name"])
    fcols = [c for c in fpkm.columns if "FPKM" in c]
    fmat = fpkm.groupby(fpkm["gene_short_name"].str.upper())[fcols].max()
    flog = np.log2(fmat + 0.25)
    rename = {
        "MCF7_CLDN4KO_FPKM (fpkm)": "MCF7_KO",
        "MCF7_WT_FPKM (fpkm)": "MCF7_WT",
        "T47D_CLDN4KO_FPKM (fpkm)": "T47D_KO",
        "T47D_WT_FPKM (fpkm)": "T47D_WT",
    }
    flog = flog.rename(columns=rename)
    for line in ("T47D", "MCF7"):
        add(
            f"GSE207704_{line}_CLDN4KO_vs_WT",
            "CLDN4",
            "CRISPR knockout, group-mean FPKM",
            "GSE207704",
            "human",
            f"{line} breast cancer cells",
            "RNA-seq group-mean FPKM",
            "genetic_ko_author_table",
            False,
            flog,
            [f"{line}_KO"],
            [f"{line}_WT"],
            note="One deposited FPKM per genotype. No replicate-level p. Gene-set p in the module table is across genes, not mice or cultures.",
        )

    edge = pd.read_csv(DATA / "GSE50927_naive.csv.gz")
    edge["symbol"] = edge["Marker.Symbol"].astype(str).str.upper()
    edge = edge.dropna(subset=["symbol"]).drop_duplicates("symbol")
    # author logFC is already KO minus WT (Cldn4 is negative in the prior audit)
    elog = edge.set_index("symbol")[["logFC"]].rename(columns={"logFC": "KO_minus_WT"})
    # represent as a one-column ratio matrix
    add(
        "GSE50927_naiveLung_Cldn4KO_vs_WT",
        "CLDN4",
        "germline knockout, author edgeR",
        "GSE50927",
        "mouse",
        "naive whole lung",
        "author edgeR logFC",
        "genetic_ko_author_table",
        False,
        elog,
        ["KO_minus_WT"],
        [],
        ratio_mode=True,
        note="GEO lists one library per condition. logFC is the author table (KO vs WT). One-sample gene tests are not biological-replicate tests. IFN-up on this series was already reported in the mouse Cldn4 KO audit.",
    )
    nhej_syms = sets["NHEJ_CORE"]
    edge_n = edge[edge["symbol"].isin(nhej_syms)][
        ["symbol", "logFC", "PValue", "FDR", "logCPM"]
    ].sort_values("symbol")
    edge_n.to_csv(OUT / "GSE50927_NHEJ_CORE_author_edgeR.tsv", sep="\t", index=False)

    # two-color KD vs CLDN4 overexpression
    arr, _ = read_series_matrix(DATA / "GSE22493_series_matrix.txt.gz")
    amap = gpl10555_map()
    amat = collapse_probes(arr, amap, log=False)  # already a log ratio
    print("GSE22493 CLDN4 ratios", amat.loc["CLDN4"].to_dict() if "CLDN4" in amat.index else "missing", flush=True)
    add(
        "GSE22493_CLDN4KD_vs_CLDN4OE",
        "CLDN4",
        "siRNA lentivirus versus CLDN4 overexpression, two-color array",
        "GSE22493",
        "human",
        "SKOV-3-IP-Luc",
        "two-color microarray log ratio",
        "confounded_kd_vs_overexpression",
        False,
        amat,
        list(amat.columns),
        [],
        ratio_mode=True,
        note="Channel 2 is CLDN4 knockdown and channel 1 is CLDN4 overexpression. This is not knockdown versus parental cells.",
    )

    modules = pd.concat(module_rows, ignore_index=True)
    genes = pd.DataFrame(gene_rows)
    summary = pd.DataFrame(summary_rows)
    modules.to_csv(OUT / "module_stats.tsv", sep="\t", index=False)
    genes.to_csv(OUT / "nhej_and_target_genes.tsv", sep="\t", index=False)
    summary.to_csv(OUT / "contrast_summary.tsv", sep="\t", index=False)
    plot_summary(summary)
    print("wrote", OUT, flush=True)


def plot_summary(summary: pd.DataFrame) -> None:
    df = summary.copy()
    # Short labels keep the comparison readable.
    df["label"] = df["gene"] + "  " + df["contrast_id"].str.replace("_", " ", regex=False)
    fig, ax = plt.subplots(figsize=(11.2, 9.0))
    y = np.arange(len(df))[::-1]
    ax.axvline(0, color="#888888", lw=0.8)
    ax.scatter(df["nhej_delta"], y + 0.12, color="#1f4e79", s=28, label="NHEJ core", zorder=3)
    ax.scatter(df["ifn_alpha_delta"], y - 0.12, color="#c45911", s=28, label="IFN-alpha hallmark", zorder=3)
    for yi, (_, r) in zip(y, df.iterrows()):
        ax.plot(
            [r["nhej_delta"], r["ifn_alpha_delta"]],
            [yi + 0.12, yi - 0.12],
            color="#e6e6e6",
            lw=1,
            zorder=1,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(df["label"], fontsize=6.5)
    ax.set_xlabel("Mean log2 difference (perturbed minus control)")
    ax.set_title("Open claudin-loss profiles: NHEJ core and type I IFN")
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / "fig_nhej_ifn_deltas.png", dpi=160)
    fig.savefig(OUT / "fig_nhej_ifn_deltas.pdf")
    plt.close()


if __name__ == "__main__":
    run()
