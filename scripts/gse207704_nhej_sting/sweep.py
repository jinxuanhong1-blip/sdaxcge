#!/usr/bin/env python3
"""Exhaustive direction sweep for GSE207704.

Thesis slice (pre-set, same quantification and same contrast):
  a c-NHEJ set with mean log2FC <= -0.10 and more genes down than up
  AND a STING or IFN set with mean log2FC >= +0.10 and more genes up than down.
  A set qualifies only if at least 3 member genes are measured.
  |mean log2FC| < 0.10 is recorded as a sign-only near miss, not a slice.

Positive log2FC = higher after CLDN4 knockout.
"""
from __future__ import annotations

import json
import math
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "gse207704_nhej_sting"
DATA = OUT / "data"
TAB = OUT / "tables"
RAW_GZ = DATA / "GSE207704_CLDN4_RNAseq.txt.gz"
SET_DIR = Path(os.environ.get("GSE207704_SETS", "/tmp/gse207704_sra/sets"))
QUANT = Path(os.environ.get("GSE207704_QUANT", "/tmp/gse207704_sra/quant"))
FASTA = Path(
    os.environ.get(
        "GSE207704_FASTA",
        "/tmp/gse207704_sra/ref/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    )
)

PSEUDO = 0.5
BAR = 0.10
ORA_CUT = 0.25
MIN_GENES = 3

CNHEJ7 = ["PRKDC", "XRCC4", "LIG4", "RIF1", "TP53BP1", "XRCC5", "XRCC6"]
ENZYMATIC5 = ["PRKDC", "XRCC4", "LIG4", "XRCC5", "XRCC6"]
KU = ["XRCC5", "XRCC6"]
LIGASE = ["XRCC4", "LIG4", "NHEJ1"]
STING5 = ["CGAS", "STING1", "TBK1", "IRF3", "STAT1"]
STING_NO_STAT1 = ["CGAS", "STING1", "TBK1", "IRF3"]
# Genes shared by Reactome STING and c-NHEJ. Removed in the pure STING slice.
STING_NHEJ_OVERLAP = ["PRKDC", "XRCC5", "XRCC6"]

SAMPLES = [
    ("T47D_WT_rep1", "T47D", "WT"),
    ("T47D_WT_rep2", "T47D", "WT"),
    ("T47D_KO_rep1", "T47D", "KO"),
    ("T47D_KO_rep2", "T47D", "KO"),
    ("MCF7_WT_rep1", "MCF7", "WT"),
    ("MCF7_WT_rep2", "MCF7", "WT"),
    ("MCF7_KO_rep1", "MCF7", "KO"),
    ("MCF7_KO_rep2", "MCF7", "KO"),
]


def load_gmt(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_enrichr(path: Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > 2:
                sets[parts[0]] = [g for g in parts[2:] if g]
    return sets


def load_list(path: Path) -> list[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]


def build_hgnc(path: Path) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """approved symbol -> all names; every name -> approved symbols."""
    frame = pd.read_csv(path, sep="\t", dtype=str, low_memory=False).fillna("")
    symbol_names: dict[str, set[str]] = {}
    name_to_approved: dict[str, set[str]] = {}
    for row in frame.itertuples(index=False):
        approved = str(row.symbol).strip()
        if not approved:
            continue
        names = {approved}
        for field in (row.alias_symbol, row.prev_symbol):
            for piece in str(field).split("|"):
                piece = piece.strip()
                if piece:
                    names.add(piece)
        symbol_names[approved] = names
        for name in names:
            name_to_approved.setdefault(name.upper(), set()).add(approved)
    return symbol_names, name_to_approved


def resolve_gene(
    symbol: str,
    present: set[str],
    symbol_names: dict[str, set[str]],
    name_to_approved: dict[str, set[str]],
) -> str | None:
    """Map a set symbol onto the matrix.

    An old alias may point at its current approved symbol. A current symbol
    may point at a previous symbol still used in Ensembl 90. A name that is
    now the approved symbol of a different gene is not used (PAN2 is both an
    old NLRP4 alias and a different gene).
    """
    if symbol in present:
        return symbol
    upper = symbol.upper()
    if upper in present:
        return upper
    candidates: list[str] = []
    for approved in sorted(name_to_approved.get(upper, ())):
        if approved in present and approved not in candidates:
            candidates.append(approved)
        for name in sorted(symbol_names.get(approved, ())):
            if name not in present or name in candidates:
                continue
            owners = name_to_approved.get(name.upper(), set())
            if owners <= {approved}:
                candidates.append(name)
    if len(candidates) == 1:
        return candidates[0]
    return None


def load_cufflinks() -> pd.DataFrame:
    raw = pd.read_csv(RAW_GZ, sep="\t", low_memory=False)
    raw = raw.rename(columns={"gene_short_name": "symbol"})
    cols = {
        "MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)",
        "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
        "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)",
        "T47D_WT": "T47D_WT_FPKM (fpkm)",
    }
    raw = raw[raw["symbol"].notna() & (raw["symbol"].astype(str).str.len() > 0)].copy()
    raw["symbol"] = raw["symbol"].astype(str).str.strip()
    for col in cols.values():
        raw[col] = pd.to_numeric(raw[col], errors="coerce").fillna(0.0)
    raw["_mean"] = raw[list(cols.values())].mean(axis=1)
    primary = (
        raw.sort_values("_mean", ascending=False)
        .drop_duplicates("symbol")
        .set_index("symbol")
        .rename(columns={v: k for k, v in cols.items()})
    )
    lg = np.log2(primary[list(cols)] + PSEUDO)
    ranks = pd.DataFrame(
        {
            "MCF7": lg["MCF7_KO"] - lg["MCF7_WT"],
            "T47D": lg["T47D_KO"] - lg["T47D_WT"],
        }
    )
    ranks["mean"] = ranks.mean(axis=1)
    return ranks


def tx_to_gene(fasta_gz: Path) -> dict[str, str]:
    import gzip

    mapping: dict[str, str] = {}
    opener = gzip.open if str(fasta_gz).endswith(".gz") else open
    with opener(fasta_gz, "rt") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].split()
            tx = header[0].split(".")[0]
            symbol = ""
            for field in header:
                if field.startswith("gene_symbol:"):
                    symbol = field.split(":", 1)[1]
            if symbol:
                mapping[header[0]] = symbol
                mapping[tx] = symbol
    return mapping


def gene_counts_from_kallisto(tx2gene: dict[str, str]) -> pd.DataFrame | None:
    frames = []
    for label, _line, _geno in SAMPLES:
        path = QUANT / label / "abundance.tsv"
        if not path.exists():
            return None
        abundance = pd.read_csv(path, sep="\t")
        abundance["gene"] = abundance["target_id"].map(
            lambda tx: tx2gene.get(tx, tx2gene.get(tx.split(".")[0], ""))
        )
        abundance = abundance[abundance["gene"] != ""]
        summed = abundance.groupby("gene")["est_counts"].sum()
        frames.append(summed.rename(label))
    counts = pd.concat(frames, axis=1).fillna(0.0)
    return counts


def size_factor_log2fc(counts: pd.DataFrame, ko: list[str], wt: list[str]) -> pd.Series:
    positive = counts.clip(lower=0)
    logged = np.log(positive.replace(0, np.nan))
    geo = np.exp(logged.mean(axis=1))
    use = geo.notna() & (geo > 0)
    ratios = positive.loc[use].div(geo.loc[use], axis=0)
    size = ratios.median(axis=0).replace(0, np.nan)
    norm = positive.div(size, axis=1)
    ko_mean = norm[ko].mean(axis=1)
    wt_mean = norm[wt].mean(axis=1)
    return np.log2((ko_mean + 1.0) / (wt_mean + 1.0))


def both_sides(counts: pd.DataFrame, ko: list[str], wt: list[str]) -> pd.Series:
    """1 if both KO replicates are below both WT replicates on normalized counts."""
    positive = counts.clip(lower=0)
    logged = np.log(positive.replace(0, np.nan))
    geo = np.exp(logged.mean(axis=1))
    use = geo.notna() & (geo > 0)
    ratios = positive.loc[use].div(geo.loc[use], axis=0)
    size = ratios.median(axis=0).replace(0, np.nan)
    norm = positive.div(size, axis=1)
    ko_max = norm[ko].max(axis=1)
    wt_min = norm[wt].min(axis=1)
    ko_min = norm[ko].min(axis=1)
    wt_max = norm[wt].max(axis=1)
    down = ko_max < wt_min
    up = ko_min > wt_max
    out = pd.Series("overlap", index=norm.index)
    out[down] = "both_KO_below_both_WT"
    out[up] = "both_KO_above_both_WT"
    return out


def deseq_log2fc(counts: pd.DataFrame, line: str) -> pd.DataFrame | None:
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as exc:
        print("pydeseq2 unavailable", exc)
        return None
    labels = [lab for lab, ln, _ in SAMPLES if ln == line]
    sub = counts[labels].copy()
    sub = sub.loc[sub.sum(axis=1) >= 10]
    rounded = sub.round().astype(int)
    meta = pd.DataFrame(
        {"condition": ["WT" if "WT" in lab else "KO" for lab in labels]},
        index=labels,
    )
    # samples x genes
    dds = DeseqDataSet(
        counts=rounded.T,
        metadata=meta,
        design="~condition",
        ref_level=["condition", "WT"],
        refit_cooks=False,
        quiet=True,
        n_cpus=2,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=["condition", "KO", "WT"],
        cooks_filter=False,
        independent_filter=True,
        quiet=True,
        n_cpus=2,
    )
    stats.summary()
    result = stats.results_df.copy()
    result.index.name = "gene"
    return result


def collect_sets() -> dict[str, dict]:
    """name -> {genes, family, source}."""
    sets: dict[str, dict] = {}

    def add(name: str, genes: list[str], family: str, source: str) -> None:
        cleaned = []
        seen = set()
        for gene in genes:
            gene = gene.strip()
            if gene and gene not in seen:
                seen.add(gene)
                cleaned.append(gene)
        sets[name] = {"genes": cleaned, "family": family, "source": source}

    add("user_cNHEJ_7", CNHEJ7, "nhej", "user panel")
    add("enzymatic_cNHEJ_5", ENZYMATIC5, "nhej", "user panel without RIF1/TP53BP1")
    add("ku_XRCC5_XRCC6", KU, "nhej", "Ku heterodimer")
    add("ligase_XRCC4_LIG4_NHEJ1", LIGASE, "nhej", "ligase complex")
    add("STING_core_5", STING5, "sting_ifn", "user STING core")
    add("STING_core_no_STAT1", STING_NO_STAT1, "sting_ifn", "STING core without STAT1")

    apm = json.loads((DATA / "ifn_apm_sets.json").read_text())["sets"]
    add("APM_21_repo", apm["CUSTOM_MHC_I_ANTIGEN_PRESENTATION"], "apm", "repo MHC-I/APM list")
    add(
        "hallmark_IFNa_repo_a8",
        apm["HALLMARK_INTERFERON_ALPHA_RESPONSE"],
        "sting_ifn",
        "a8 Hallmark 2020 IFN-α",
    )
    add(
        "hallmark_IFNg_repo_a8",
        apm["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
        "sting_ifn",
        "a8 Hallmark 2020 IFN-γ",
    )

    for version, filename in [
        ("2023.1", "h.all.v2023.1.Hs.symbols.gmt"),
        ("2023.2", "h.all.v2023.2.Hs.symbols.gmt"),
        ("2024.1", "h.all.v2024.1.Hs.symbols.gmt"),
        ("2025.1", "h.all.v2025.1.Hs.symbols.gmt"),
    ]:
        gmt = load_gmt(SET_DIR / filename)
        add(f"hallmark_IFNa_{version}", gmt["HALLMARK_INTERFERON_ALPHA_RESPONSE"], "sting_ifn", f"MSigDB {version}")
        add(f"hallmark_IFNg_{version}", gmt["HALLMARK_INTERFERON_GAMMA_RESPONSE"], "sting_ifn", f"MSigDB {version}")
        add(f"hallmark_DNA_REPAIR_{version}", gmt["HALLMARK_DNA_REPAIR"], "damage", f"MSigDB {version}")

    reactome = load_enrichr(SET_DIR / "enrichr_reactome.txt")
    wanted = {
        "Nonhomologous End-Joining (NHEJ) R-HSA-5693571": ("reactome_NHEJ", "nhej"),
        "HDR Thru MMEJ (alt-NHEJ) R-HSA-5685939": ("reactome_altNHEJ_MMEJ", "damage"),
        "STING Mediated Induction Of Host Immune Responses R-HSA-1834941": ("reactome_STING", "sting_ifn"),
        "Regulation Of Innate Immune Responses To Cytosolic DNA R-HSA-3134975": (
            "reactome_cytosolic_DNA",
            "sting_ifn",
        ),
        "Interferon Alpha/Beta Signaling R-HSA-909733": ("reactome_IFNab", "sting_ifn"),
        "Interferon Gamma Signaling R-HSA-877300": ("reactome_IFNg", "sting_ifn"),
        "Interferon Signaling R-HSA-913531": ("reactome_IFN_signaling", "sting_ifn"),
        "DNA Double Strand Break Response R-HSA-5693606": ("reactome_DSB_response", "damage"),
    }
    for raw_name, (name, family) in wanted.items():
        if raw_name not in reactome:
            raise SystemExit(f"missing reactome set {raw_name}")
        add(name, reactome[raw_name], family, "Enrichr Reactome_2022")
    sting_pure = [g for g in sets["reactome_STING"]["genes"] if g not in STING_NHEJ_OVERLAP]
    add("reactome_STING_drop_NHEJ_genes", sting_pure, "sting_ifn", "Reactome STING minus PRKDC/XRCC5/XRCC6")

    add("GO_NHEJ_0006303", load_list(SET_DIR / "go_nhej.txt"), "nhej", "QuickGO GO:0006303")
    add(
        "GO_type1_IFN_0060337",
        load_list(SET_DIR / "go_type1_ifn_signaling.txt"),
        "sting_ifn",
        "QuickGO GO:0060337",
    )
    add(
        "GO_response_IFNB_0035456",
        load_list(SET_DIR / "go_response_to_ifnb.txt"),
        "sting_ifn",
        "QuickGO GO:0035456",
    )
    return sets


def score_one(rank: pd.Series, genes: list[str], present: set[str], symbol_names, name_to_approved) -> dict:
    mapped = []
    missing = []
    used = set()
    for gene in genes:
        hit = resolve_gene(gene, present, symbol_names, name_to_approved)
        if hit is None or hit in used:
            if hit is None:
                missing.append(gene)
            continue
        used.add(hit)
        mapped.append(hit)
    values = rank.reindex(mapped).dropna()
    n = int(len(values))
    n_pos = int((values > 0).sum()) if n else 0
    n_neg = int((values < 0).sum()) if n else 0
    if n_pos + n_neg == 0:
        sign_p = np.nan
    else:
        sign_p = float(binomtest(n_pos, n_pos + n_neg, 0.5, alternative="two-sided").pvalue)
    if n >= 3:
        bg = rank.drop(index=list(used), errors="ignore")
        _u, mw_p = mannwhitneyu(values, bg, alternative="two-sided")
        mw_p = float(mw_p)
    else:
        mw_p = np.nan
    universe = rank.dropna()
    down = set(universe.index[universe < -ORA_CUT])
    up = set(universe.index[universe > ORA_CUT])
    measured = set(values.index)

    def ora(target: set[str]) -> float:
        if n < MIN_GENES or not target:
            return np.nan
        a = len(measured & target)
        b = len(measured - target)
        c = len(target - measured)
        d = len(set(universe.index) - measured - target)
        _odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        return float(p)

    mean = float(values.mean()) if n else np.nan
    median = float(values.median()) if n else np.nan
    if n < MIN_GENES or not np.isfinite(mean):
        direction = "too_few"
    elif mean <= -BAR and n_neg > n_pos:
        direction = "down"
    elif mean >= BAR and n_pos > n_neg:
        direction = "up"
    elif mean < 0 and n_neg > n_pos:
        direction = "sign_down_below_bar"
    elif mean > 0 and n_pos > n_neg:
        direction = "sign_up_below_bar"
    else:
        direction = "mixed_or_flat"
    return {
        "n_set": len(genes),
        "n_measured": n,
        "n_missing": len(missing),
        "mean_log2FC": mean,
        "median_log2FC": median,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "sign_p": sign_p,
        "mw_p": mw_p,
        "ora_up_p": ora(up),
        "ora_down_p": ora(down),
        "direction": direction,
        "mapped_genes": ",".join(mapped),
        "missing_genes": ",".join(missing[:40]),
    }


def gsea_table(rank: pd.Series, sets: dict[str, dict]) -> dict[str, dict]:
    eligible = {}
    for name, spec in sets.items():
        members = spec["genes"]
        # GSEA uses symbols as written; alias expansion is in the mean test.
        # Re-resolve into the rank index so version/alias differences are scored.
        eligible[name] = members
    # Build resolved membership outside; caller passes resolved lists.
    scored = gsea_prerank(rank.dropna().sort_values(ascending=False), eligible, nperm=NPERM, seed=SEED)
    if scored.empty:
        return {}
    scored["fdr"] = bh_fdr(scored["nom_p"])
    out = {}
    for row in scored.itertuples(index=False):
        out[row.term] = {"nes": float(row.nes), "gsea_p": float(row.nom_p), "gsea_fdr": float(row.fdr)}
    return out


def rank_sources(symbol_names, name_to_approved) -> dict[str, pd.DataFrame]:
    sources = {"cufflinks_group_mean_FPKM": load_cufflinks()}
    if FASTA.exists() and all((QUANT / lab / "abundance.tsv").exists() for lab, _, _ in SAMPLES):
        print("loading kallisto counts")
        tx2gene = tx_to_gene(FASTA)
        counts = gene_counts_from_kallisto(tx2gene)
        if counts is not None:
            t47_cols = ["T47D_KO_rep1", "T47D_KO_rep2", "T47D_WT_rep1", "T47D_WT_rep2"]
            mcf_cols = ["MCF7_KO_rep1", "MCF7_KO_rep2", "MCF7_WT_rep1", "MCF7_WT_rep2"]
            t47 = size_factor_log2fc(counts, ["T47D_KO_rep1", "T47D_KO_rep2"], ["T47D_WT_rep1", "T47D_WT_rep2"])
            mcf = size_factor_log2fc(counts, ["MCF7_KO_rep1", "MCF7_KO_rep2"], ["MCF7_WT_rep1", "MCF7_WT_rep2"])
            # Counts below 10 are not a direction. log2((4+1)/(1+1)) on TMEM173 is not induction.
            t47 = t47.where(counts[t47_cols].mean(axis=1) >= 10)
            mcf = mcf.where(counts[mcf_cols].mean(axis=1) >= 10)
            ranks = pd.DataFrame({"T47D": t47, "MCF7": mcf})
            ranks["mean"] = ranks.mean(axis=1, skipna=False)
            sources["kallisto_ensembl90_norm_log2FC"] = ranks
            counts.to_csv(TAB / "kallisto_gene_counts.tsv", sep="\t")
            side_rows = []
            for line, ko, wt in [
                ("T47D", ["T47D_KO_rep1", "T47D_KO_rep2"], ["T47D_WT_rep1", "T47D_WT_rep2"]),
                ("MCF7", ["MCF7_KO_rep1", "MCF7_KO_rep2"], ["MCF7_WT_rep1", "MCF7_WT_rep2"]),
            ]:
                side = both_sides(counts, ko, wt)
                side_rows.append(side.rename(line))
                try:
                    de = deseq_log2fc(counts, line)
                except Exception as exc:
                    print(f"DESeq2 failed for {line}: {exc}")
                    de = None
                if de is not None:
                    de.to_csv(TAB / f"deseq2_{line}.tsv", sep="\t")
            pd.concat(side_rows, axis=1).to_csv(TAB / "kallisto_replicate_separation.tsv", sep="\t")
    return sources


def alias_report(present_by_source: dict[str, set[str]], symbol_names, name_to_approved) -> pd.DataFrame:
    rows = []
    queries = {
        "STING1": ["STING1", "TMEM173", "STING", "ERIS", "MPYS", "MITA", "NET23", "FLJ38577"],
        "CGAS": ["CGAS", "MB21D1", "C6orf150"],
    }
    # official names from HGNC
    for approved in ("STING1", "CGAS"):
        for name in symbol_names.get(approved, ()):
            queries[approved].append(name)
    for source, present in present_by_source.items():
        for approved, names in queries.items():
            seen = []
            for name in names:
                if name in seen:
                    continue
                seen.append(name)
                hit = name if name in present else ""
                rows.append(
                    {
                        "source": source,
                        "query_gene": approved,
                        "alias": name,
                        "found_as": hit,
                        "official_hgnc": "|".join(sorted(symbol_names.get(approved, []))),
                    }
                )
    return pd.DataFrame(rows)


def fmt(value, digits=3) -> str:
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return "NA"
    return f"{value:.{digits}f}"


def write_sweep_md(scores: pd.DataFrame, joints: pd.DataFrame, near: pd.DataFrame, alias: pd.DataFrame, sources: list[str]) -> None:
    fastq_done = "kallisto_ensembl90_norm_log2FC" in sources
    if joints is None or joints.empty or "kind" not in getattr(joints, "columns", []):
        full_hits = pd.DataFrame()
    else:
        full_hits = joints[joints["kind"] == "thesis_slice"]
    n_full = len(full_hits)
    if fastq_done and n_full == 0:
        status = "FINAL discordant cancer-line KD"
    elif fastq_done and n_full > 0:
        status = "directional slice found — see joint table; not a both-line FDR result unless the table says so"
    else:
        status = "cufflinks sweep only; FASTQ quantification not in this file yet"

    def mini(frame: pd.DataFrame, cols: list[str], limit: int = 40) -> str:
        if frame is None or frame.empty:
            return "_none_"
        view = frame[cols].head(limit)
        header = "| " + " | ".join(cols) + " |"
        sep = "|" + "|".join(["---"] * len(cols)) + "|"
        lines = [header, sep]
        for row in view.itertuples(index=False):
            cells = []
            for value in row:
                if isinstance(value, float):
                    cells.append(fmt(value))
                else:
                    cells.append(str(value))
            lines.append("| " + " | ".join(cells) + " |")
        return "\n".join(lines)

    focus_names = [
        "user_cNHEJ_7",
        "enzymatic_cNHEJ_5",
        "reactome_NHEJ",
        "GO_NHEJ_0006303",
        "STING_core_5",
        "reactome_STING",
        "reactome_STING_drop_NHEJ_genes",
        "reactome_cytosolic_DNA",
        "hallmark_IFNa_2025.1",
        "hallmark_IFNg_2025.1",
        "hallmark_IFNa_2023.1",
        "APM_21_repo",
        "hallmark_DNA_REPAIR_2025.1",
        "reactome_altNHEJ_MMEJ",
    ]
    focus = scores[scores["set"].isin(focus_names) & scores["contrast"].isin(["T47D", "MCF7"])]
    focus = focus.sort_values(["source", "set", "contrast"])
    sting_alias = alias[(alias["query_gene"] == "STING1") & (alias["found_as"] != "")]
    text = f"""# Sweep — GSE207704 c-NHEJ / STING / IFN

Status: **{status}**

Positive log2FC means higher after CLDN4 knockout. A thesis slice, fixed before reading these set scores, is a c-NHEJ set and a STING/IFN set in the **same** quantification and the **same** contrast, each with at least 3 measured genes, mean log2FC at least 0.10 in the thesis direction, and a majority of genes on that side. Sign-only results inside ±0.10 are near misses.

## What was tried

- Cufflinks group-mean FPKM (the GEO supplementary table), T47D, MCF7, and the mean of the two lines.
- Open SRA runs SRR20029118–SRR20029125 (public S3 `sra-pub-run-odp`). ENA HTTPS FASTQ failed an SSL handshake. `fasterq-dump` segfaulted; `fastq-dump` in 8 million-read chunks did not. Spots are 51 bp single-end, not the 100 bp length named in the paper. {"Kallisto 0.51.1 on Ensembl 90 cDNA, single-end fragment prior `-l 200 -s 20`, transcript counts summed to gene symbols. Set tests use size-factor log2FC and drop genes with mean estimated count < 10 in that contrast. PyDESeq2 Wald is per line." if fastq_done else "Kallisto quantification was not finished when this file was written."}
- HGNC aliases for STING1 (approved symbol, previous symbol TMEM173, aliases STING, ERIS, MPYS, MITA, NET23, FLJ38577) and CGAS (MB21D1, C6orf150), matched to symbols actually in each matrix. The cufflinks locus window chr5:139475534–139482790 (Ensembl 90 TMEM173) has no row.
- Hallmark IFN-α, IFN-γ, and DNA repair from MSigDB 2023.1, 2023.2, 2024.1, and 2025.1, plus the repo’s Hallmark 2020 lists. IFN-α membership is identical across 2023.1–2025.1. IFN-γ differs by METTL7B (2023.1) versus TMT1B (2024.1/2025.1).
- Reactome 2022 via Enrichr: NHEJ, alt-NHEJ/MMEJ, STING, cytosolic DNA, IFN-α/β, IFN-γ, IFN signaling, DSB response. Reactome STING includes PRKDC, XRCC5, and XRCC6; those three were also dropped in a separate set.
- QuickGO GO:0006303 (NHEJ), GO:0060337 (type I IFN signaling), GO:0035456 (response to IFNβ).
- QuickGO GO:0031040 (micronucleus) and GO:0032125 (micronucleus organization): the terms exist and returned **zero** annotations, so no micronucleus gene set was scored.
- Methods other than GSEA: mean, median, sign test, Mann–Whitney against the rest of the rank, and Fisher exact overlap with genes at log2FC < −0.25 or > +0.25. Prerank GSEA (p=1, 1000 perms, seed 42) only when at least 8 genes are in the rank.
- Line emphasis is the T47D column versus the MCF7 column, not a pooled test.

## STING1 aliases found

{mini(sting_alias, ["source", "alias", "found_as"]) if len(sting_alias) else "_No STING1 alias is present in the scored matrices._"}

## Focused set scores

{mini(focus, ["source", "set", "contrast", "n_measured", "mean_log2FC", "n_pos", "n_neg", "direction", "nes"], limit=80)}

## Thesis slices (mean |log2FC| ≥ 0.10, majority of genes, ≥3 genes, same source and contrast)

{mini(full_hits, ["source", "contrast", "nhej_set", "nhej_mean", "ifn_set", "ifn_mean"]) if n_full else "_No pair met the bar._"}

## Sign-only near misses (right signs, |mean| < 0.10)

{mini(near, ["source", "contrast", "nhej_set", "nhej_mean", "ifn_set", "ifn_mean"], limit=30)}

Full rows: `tables/sweep_set_scores.tsv`, `tables/sweep_joint.tsv`, `tables/sweep_aliases.tsv`.
"""
    (OUT / "SWEEP.md").write_text(text)


def main() -> None:
    TAB.mkdir(parents=True, exist_ok=True)
    symbol_names, name_to_approved = build_hgnc(SET_DIR / "hgnc_complete_set.txt")
    sets = collect_sets()
    # freeze the exact sets used
    frozen = {name: spec["genes"] for name, spec in sets.items()}
    (DATA / "sweep_sets.json").write_text(json.dumps({"sets": frozen, "bar": BAR}, indent=2))
    sources = rank_sources(symbol_names, name_to_approved)
    rows = []
    present_by_source = {}
    for source, ranks in sources.items():
        present = set(ranks.index)
        present_by_source[source] = present
        for contrast in ["T47D", "MCF7", "mean"]:
            rank = ranks[contrast].dropna()
            resolved_for_gsea = {}
            scored_rows = []
            for name, spec in sets.items():
                stats = score_one(rank, spec["genes"], present, symbol_names, name_to_approved)
                mapped = [g for g in stats["mapped_genes"].split(",") if g]
                resolved_for_gsea[name] = mapped
                stats.update(
                    {
                        "source": source,
                        "contrast": contrast,
                        "set": name,
                        "family": spec["family"],
                        "set_source": spec["source"],
                    }
                )
                scored_rows.append(stats)
            gsea = {}
            eligible = {k: v for k, v in resolved_for_gsea.items() if len(v) >= 8}
            if eligible:
                scored = gsea_prerank(rank.sort_values(ascending=False), eligible, nperm=NPERM, seed=SEED)
                if not scored.empty:
                    scored["fdr"] = bh_fdr(scored["nom_p"])
                    for rec in scored.to_dict(orient="records"):
                        gsea[rec["term"]] = rec
            for stats in scored_rows:
                extra = gsea.get(stats["set"], {})
                stats["nes"] = extra.get("nes", np.nan)
                stats["gsea_p"] = extra.get("nom_p", np.nan)
                stats["gsea_fdr"] = extra.get("fdr", np.nan)
                rows.append(stats)
    scores = pd.DataFrame(rows)
    scores.to_csv(TAB / "sweep_set_scores.tsv", sep="\t", index=False, float_format="%.6g")
    alias = alias_report(present_by_source, symbol_names, name_to_approved)
    alias.to_csv(TAB / "sweep_aliases.tsv", sep="\t", index=False)

    joint_rows = []
    near_rows = []
    for (source, contrast), block in scores.groupby(["source", "contrast"]):
        nhej = block[(block["family"] == "nhej") & (block["n_measured"] >= MIN_GENES)]
        ifn = block[(block["family"] == "sting_ifn") & (block["n_measured"] >= MIN_GENES)]
        for nrec in nhej.itertuples(index=False):
            for irec in ifn.itertuples(index=False):
                n_down = nrec.direction == "down"
                i_up = irec.direction == "up"
                n_sign = nrec.direction == "sign_down_below_bar"
                i_sign = irec.direction == "sign_up_below_bar"
                if n_down and i_up:
                    joint_rows.append(
                        {
                            "kind": "thesis_slice",
                            "source": source,
                            "contrast": contrast,
                            "nhej_set": nrec.set,
                            "nhej_mean": nrec.mean_log2FC,
                            "nhej_n": nrec.n_measured,
                            "ifn_set": irec.set,
                            "ifn_mean": irec.mean_log2FC,
                            "ifn_n": irec.n_measured,
                        }
                    )
                elif (n_down or n_sign) and (i_up or i_sign) and not (n_down and i_up):
                    near_rows.append(
                        {
                            "source": source,
                            "contrast": contrast,
                            "nhej_set": nrec.set,
                            "nhej_mean": nrec.mean_log2FC,
                            "nhej_direction": nrec.direction,
                            "ifn_set": irec.set,
                            "ifn_mean": irec.mean_log2FC,
                            "ifn_direction": irec.direction,
                        }
                    )
    joints = pd.DataFrame(joint_rows)
    near = pd.DataFrame(near_rows)
    joints.to_csv(TAB / "sweep_joint.tsv", sep="\t", index=False, float_format="%.6g")
    near.to_csv(TAB / "sweep_near_miss.tsv", sep="\t", index=False, float_format="%.6g")
    write_sweep_md(scores, joints, near, alias, list(sources))
    print(f"sources {list(sources)}")
    print(f"thesis slices {0 if joints.empty else len(joints)}")
    print(f"near misses {0 if near.empty else len(near)}")
    focus = scores[scores["set"].isin(["user_cNHEJ_7", "STING_core_5", "hallmark_IFNa_2025.1", "hallmark_IFNg_2025.1", "reactome_NHEJ", "reactome_STING"])]
    print(focus[["source", "contrast", "set", "n_measured", "mean_log2FC", "direction", "nes"]].to_string(index=False))


if __name__ == "__main__":
    main()
