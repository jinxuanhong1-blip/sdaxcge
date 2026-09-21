#!/usr/bin/env python3
"""GSE207704 CLDN4 knockout breast lines: NHEJ, STING, IFN/APM scores vs WT.

Public GEO processed file only (GSE207704_CLDN4_RNAseq.txt.gz). The series
has 2 biological replicates per genotype, but the deposited table is already
collapsed to one cufflinks FPKM per group. WT is the parental "No treatment"
arm (the negative-control comparator). There is no siRNA NC arm.

Gene score = unweighted mean of log2(FPKM + 0.5) across mapped genes.
delta = score(KO) - score(WT). Positive = higher after CLDN4 loss.

GSEA is preranked (weighted KS p=1, 1000 gene-set permutations, seed=42),
same engine as scripts/gse207704_nhej_sting_ifn/gsea_core.py. Positive NES =
enriched at the CLDN4-loss end of the rank.
"""
from __future__ import annotations

import math
import sys
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gsea_core import NPERM, SEED, bh_fdr, gsea_prerank  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "methods" / "gse207704_nhej_sting_ifn"
DATA = OUT / "data" / "GSE207704_CLDN4_RNAseq.txt.gz"
GMT = OUT / "reactome_sets.gmt"
GEO_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207704/suppl/"
    "GSE207704_CLDN4_RNAseq.txt.gz"
)
PSEUDO = 0.5
LOW_FPKM = 1.0
WEAK = 0.25

# Deposited symbols are Ensembl/GRCh38 v90 (2017). Map current HGNC names.
ALIAS = {
    "CGAS": "MB21D1",
    "STING1": "TMEM173",
    "H2AX": "H2AFX",
    "RIGI": "DDX58",
}

NHEJ_CORE = ["PRKDC", "LIG4", "XRCC4", "XRCC5", "XRCC6", "NHEJ1"]
STING_AXIS = ["CGAS", "STING1", "TBK1", "IRF3"]
IFN_GENES = [
    "IFI27", "OAS2", "IFIT1", "MX1", "ISG15",
    "IFIT2", "IFIT3", "IFIT5", "MX2", "OAS1", "OAS3", "OASL", "RSAD2",
    "USP18", "IFI6", "IFI35", "IFI44", "IFI44L", "IFI16", "BST2",
    "IFITM1", "IFITM2", "IFITM3", "DDX58", "IFIH1", "DDX60", "HERC5",
    "XAF1", "EIF2AK2", "LY6E", "CMPK2", "EPSTI1", "SAMD9L", "TRIM22",
    "ZBP1", "PLSCR1",
    "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "JAK1", "JAK2",
    "IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "SOCS1", "SOCS3",
    "GBP1", "GBP2", "GBP4", "GBP5", "CXCL9", "CXCL10", "CXCL11",
    "CIITA", "IDO1", "IFNG", "IFNB1",
]
APM_GENES = [
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M", "NLRC5",
    "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "PSMB10", "ERAP1", "ERAP2",
    "CALR", "CANX", "PDIA3", "SEC61A1", "PSME1", "PSME2",
]
SCORE_SETS = {
    "NHEJ_CORE": NHEJ_CORE,
    "STING_AXIS": STING_AXIS,
    "IFN": IFN_GENES,
    "APM": APM_GENES,
}
# Pathway GSEA uses the Reactome GMT (min size 8). The two named panels are
# scored separately and also preranked with a lower floor, because 6 and 4
# genes sit under the usual GSEA size cutoff.
SMALL_GSEA = ("NHEJ_CORE", "STING_AXIS")
LINES = ("T47D", "MCF7")


def log2fc(ko: float, wt: float) -> float:
    return math.log2((ko + PSEUDO) / (wt + PSEUDO))


def direction(lfc: float, max_fpkm: float) -> str:
    if max_fpkm < LOW_FPKM:
        return "LOW"
    if lfc > WEAK:
        return "UP"
    if lfc < -WEAK:
        return "DOWN"
    return "FLAT"


def ensure_matrix() -> None:
    if DATA.exists() and DATA.stat().st_size > 1000:
        return
    DATA.parent.mkdir(parents=True, exist_ok=True)
    print(f"downloading {GEO_URL}")
    urllib.request.urlretrieve(GEO_URL, DATA)


def load_fpkm() -> pd.DataFrame:
    raw = pd.read_csv(DATA, sep="\t", low_memory=False)
    cols = {
        "MCF7_KO": "MCF7_CLDN4KO_FPKM (fpkm)",
        "MCF7_WT": "MCF7_WT_FPKM (fpkm)",
        "T47D_KO": "T47D_CLDN4KO_FPKM (fpkm)",
        "T47D_WT": "T47D_WT_FPKM (fpkm)",
    }
    df = raw.rename(columns={"gene_short_name": "symbol"})
    df = df[df["symbol"].notna()].copy()
    df["symbol"] = df["symbol"].astype(str).str.strip()
    df = df[df["symbol"].ne("") & df["symbol"].ne("nan")]
    for c in cols.values():
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df["_m"] = df[list(cols.values())].mean(axis=1)
    # Same locus pick as the earlier GSE207704 prerank: highest mean FPKM.
    fpkm = (
        df.sort_values("_m", ascending=False)
        .drop_duplicates("symbol")
        .set_index("symbol")[list(cols.values())]
        .rename(columns={v: k for k, v in cols.items()})
    )
    return fpkm


def resolve(symbol: str, index: pd.Index) -> str | None:
    if symbol in index:
        return symbol
    alt = ALIAS.get(symbol)
    if alt and alt in index:
        return alt
    return None


def gene_table(fpkm: pd.DataFrame) -> pd.DataFrame:
    rows = []
    catalog = [("QC", "CLDN4"), ("QC", "TACSTD2")]
    for set_name, genes in SCORE_SETS.items():
        for g in genes:
            catalog.append((set_name, g))
    for set_name, symbol in catalog:
        hit = resolve(symbol, fpkm.index)
        rec = {
            "set": set_name,
            "symbol": symbol,
            "deposited_symbol": hit or "",
            "alias_used": bool(hit and hit != symbol),
            "mapped": hit is not None,
        }
        if hit is None:
            rec.update(
                {
                    "note": "absent from cufflinks table (annotation gap, not a measured zero)",
                    "dir_T47D": "ABSENT",
                    "dir_MCF7": "ABSENT",
                    "consensus": "ABSENT",
                }
            )
        else:
            r = fpkm.loc[hit]
            for line in LINES:
                ko = float(r[f"{line}_KO"])
                wt = float(r[f"{line}_WT"])
                lfc = log2fc(ko, wt)
                rec[f"{line}_WT"] = wt
                rec[f"{line}_KO"] = ko
                rec[f"log2FC_{line}"] = lfc
                rec[f"dir_{line}"] = direction(lfc, max(ko, wt))
            rec["log2FC_mean"] = 0.5 * (rec["log2FC_T47D"] + rec["log2FC_MCF7"])
            dirs = {rec["dir_T47D"], rec["dir_MCF7"]}
            if "UP" in dirs and "DOWN" in dirs:
                rec["consensus"] = "DISCORDANT"
            elif dirs <= {"UP"} or (dirs <= {"UP", "FLAT", "LOW"} and "UP" in dirs):
                rec["consensus"] = "UP" if "UP" in dirs else "FLAT"
            elif dirs <= {"DOWN"} or (dirs <= {"DOWN", "FLAT", "LOW"} and "DOWN" in dirs):
                rec["consensus"] = "DOWN" if "DOWN" in dirs else "FLAT"
            elif dirs <= {"LOW"}:
                rec["consensus"] = "LOW"
            else:
                rec["consensus"] = "FLAT"
            rec["note"] = "descriptive; group-mean FPKM, no gene-level FDR"
        rows.append(rec)
    return pd.DataFrame(rows)


def score_table(genes: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for set_name, members in SCORE_SETS.items():
        sub = genes[genes["set"].eq(set_name) & genes["symbol"].isin(members)]
        # one row per requested symbol
        sub = sub.drop_duplicates("symbol")
        absent = sub.loc[~sub["mapped"], "symbol"].tolist()
        mapped = sub[sub["mapped"]]
        for line in LINES:
            expressed = mapped[mapped[f"dir_{line}"].ne("LOW")] if len(mapped) else mapped
            low = mapped[mapped[f"dir_{line}"].eq("LOW")] if len(mapped) else mapped

            def _delta(frame: pd.DataFrame) -> tuple[float, float, float]:
                if frame.empty:
                    return (np.nan, np.nan, np.nan)
                wt = np.log2(frame[f"{line}_WT"].astype(float) + PSEUDO)
                ko = np.log2(frame[f"{line}_KO"].astype(float) + PSEUDO)
                return (float(wt.mean()), float(ko.mean()), float(ko.mean() - wt.mean()))

            wt_e, ko_e, d_e = _delta(expressed)
            _, _, d_all = _delta(mapped)
            rows.append(
                {
                    "set": set_name,
                    "line": line,
                    "contrast": f"{line} CLDN4-/- vs WT(NC)",
                    "n_requested": len(members),
                    "n_mapped": int(len(mapped)),
                    "n_expressed": int(len(expressed)),
                    "n_low": int(len(low)),
                    "n_absent": int(len(absent)),
                    "absent_genes": ",".join(absent) if absent else ".",
                    "score_WT_NC": wt_e,
                    "score_KO": ko_e,
                    "delta_KO_minus_WT": d_e,
                    "delta_all_mapped": d_all,
                    "n_up": int((expressed[f"dir_{line}"] == "UP").sum()) if len(expressed) else 0,
                    "n_down": int((expressed[f"dir_{line}"] == "DOWN").sum()) if len(expressed) else 0,
                    "n_flat": int((expressed[f"dir_{line}"] == "FLAT").sum()) if len(expressed) else 0,
                    "up_genes": ",".join(expressed.loc[expressed[f"dir_{line}"].eq("UP"), "symbol"]) or ".",
                    "down_genes": ",".join(expressed.loc[expressed[f"dir_{line}"].eq("DOWN"), "symbol"]) or ".",
                    "score_definition": (
                        "mean log2(FPKM+0.5) of mapped genes with max(WT,KO) FPKM>=1; "
                        "delta = KO minus parental WT"
                    ),
                }
            )
    out = pd.DataFrame(rows)
    # Descriptive both-line label. Not a p-value.
    labels = {}
    for set_name, sub in out.groupby("set"):
        d = {r.line: r.delta_KO_minus_WT for r in sub.itertuples()}
        a, b = d["T47D"], d["MCF7"]
        if not np.isfinite(a) or not np.isfinite(b):
            lab = "not_scored"
        elif a <= -0.10 and b <= -0.10:
            lab = "both_lines_down"
        elif a >= 0.10 and b >= 0.10:
            lab = "both_lines_up"
        elif (a > 0.10 and b < -0.10) or (a < -0.10 and b > 0.10):
            lab = "lines_discordant"
        else:
            lab = "near_zero"
        labels[set_name] = lab
    out["both_lines"] = out["set"].map(labels)
    return out


def read_gmt(path: Path) -> dict[str, list[str]]:
    sets = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        sets[parts[0]] = parts[2:]
    return sets


def ranks_from(fpkm: pd.DataFrame) -> dict[str, pd.Series]:
    lg = np.log2(fpkm.astype(float) + PSEUDO)
    ranks = {
        "T47D": (lg["T47D_KO"] - lg["T47D_WT"]).dropna().sort_values(ascending=False),
        "MCF7": (lg["MCF7_KO"] - lg["MCF7_WT"]).dropna().sort_values(ascending=False),
    }
    both = pd.concat([ranks["T47D"], ranks["MCF7"]], axis=1).mean(axis=1)
    ranks["mean_both_lines"] = both.dropna().sort_values(ascending=False)
    return ranks


def run_gsea(fpkm: pd.DataFrame) -> pd.DataFrame:
    library = read_gmt(GMT)
    library.update({k: list(v) for k, v in SCORE_SETS.items()})
    # Resolve aliases onto deposited symbols, drop unresolved, keep unique.
    resolved = {}
    for name, genes in library.items():
        hits = []
        for g in genes:
            hit = resolve(g, fpkm.index)
            if hit and hit not in hits:
                hits.append(hit)
        resolved[name] = hits

    rank_map = ranks_from(fpkm)
    rows = []
    for contrast, rank in rank_map.items():
        pathway_names = [n for n in resolved if n not in SMALL_GSEA]
        small_names = [n for n in SMALL_GSEA if n in resolved]
        pathway = gsea_prerank(
            rank, {n: resolved[n] for n in pathway_names}, nperm=NPERM, seed=SEED, min_size=8
        )
        small = gsea_prerank(
            rank, {n: resolved[n] for n in small_names}, nperm=NPERM, seed=SEED, min_size=3, max_size=15
        )
        if len(pathway):
            pathway["fdr"] = bh_fdr(pathway["nom_p"])
            pathway["fdr_family"] = "pathway_sets_min8"
            pathway["size_note"] = "standard floor"
        if len(small):
            small["fdr"] = np.nan
            small["fdr_family"] = "not_applied_set_below_8"
            small["size_note"] = "named panel below the 8-gene GSEA floor; nominal p only"
        for frame, kind in ((pathway, "pathway"), (small, "named_panel")):
            if not len(frame):
                continue
            frame = frame.copy()
            frame["contrast"] = contrast
            frame["kind"] = kind
            frame["n_genes_ranked"] = int(rank.shape[0])
            frame["n_perm"] = NPERM
            frame["seed"] = SEED
            rows.append(frame)
        # Record requested sets that fell under the floor.
        tested = set()
        for frame in (pathway, small):
            if len(frame):
                tested.update(frame["term"])
        for name, hits in resolved.items():
            if name in tested:
                continue
            floor = 3 if name in SMALL_GSEA else 8
            if len(hits) < floor:
                rows.append(
                    pd.DataFrame(
                        [
                            {
                                "term": name,
                                "es": np.nan,
                                "nes": np.nan,
                                "nom_p": np.nan,
                                "n_set_in_rank": len(hits),
                                "mean_stat": np.nan,
                                "lead_genes": "",
                                "n_lead": 0,
                                "fdr": np.nan,
                                "fdr_family": "skipped",
                                "size_note": f"skipped: {len(hits)} genes in rank < floor {floor}",
                                "contrast": contrast,
                                "kind": "skipped",
                                "n_genes_ranked": int(rank.shape[0]),
                                "n_perm": NPERM,
                                "seed": SEED,
                            }
                        ]
                    )
                )
    return pd.concat(rows, ignore_index=True)


def _fmt(x, nd=3) -> str:
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    return f"{float(x):.{nd}f}"


def write_finding(fpkm: pd.DataFrame, genes: pd.DataFrame, scores: pd.DataFrame, gsea: pd.DataFrame) -> None:
    def lfc(symbol: str, line: str) -> float:
        hit = resolve(symbol, fpkm.index)
        r = fpkm.loc[hit]
        return log2fc(float(r[f"{line}_KO"]), float(r[f"{line}_WT"]))

    def gene_line(symbol: str) -> str:
        sub = genes[genes["symbol"].eq(symbol)].iloc[0]
        if not sub["mapped"]:
            return f"| {symbol} | absent | absent | absent | ABSENT |"
        dep = sub["deposited_symbol"]
        label = symbol if dep == symbol else f"{symbol} ({dep})"
        return (
            f"| {label} | {_fmt(sub['log2FC_T47D'])} {sub['dir_T47D']} | "
            f"{_fmt(sub['log2FC_MCF7'])} {sub['dir_MCF7']} | {_fmt(sub['log2FC_mean'])} | {sub['consensus']} |"
        )

    def score_line(set_name: str) -> str:
        sub = scores[scores["set"].eq(set_name)]
        bits = []
        for r in sub.itertuples():
            bits.append(
                f"{r.line} Δ {_fmt(r.delta_KO_minus_WT)} "
                f"(n={r.n_expressed}/{r.n_requested} expressed, up {r.n_up}, down {r.n_down}, flat {r.n_flat})"
            )
        label = sub["both_lines"].iloc[0]
        return f"{set_name}: {label}. " + "; ".join(bits)

    def gsea_md(kind: str) -> str:
        sub = gsea[gsea["kind"].eq(kind)].copy()
        if sub.empty:
            return "_none_\n"
        lines = [
            "| Contrast | Set | NES | nominal p | FDR | n in rank | mean log2FC |",
            "|---|---|---:|---:|---:|---:|---:|",
        ]
        order = ["T47D", "MCF7", "mean_both_lines"]
        sub["_o"] = sub["contrast"].map({k: i for i, k in enumerate(order)})
        sub = sub.sort_values(["_o", "term"])
        for r in sub.itertuples():
            fdr = "NA" if r.fdr_family.startswith("not_applied") or not np.isfinite(r.fdr) else _fmt(r.fdr, 3)
            lines.append(
                f"| {r.contrast} | {r.term} | {_fmt(r.nes)} | {_fmt(r.nom_p)} | {fdr} | "
                f"{int(r.n_set_in_rank)} | {_fmt(r.mean_stat)} |"
            )
        return "\n".join(lines) + "\n"

    nhej_rows = "\n".join(gene_line(g) for g in NHEJ_CORE)
    sting_rows = "\n".join(gene_line(g) for g in STING_AXIS)
    cldn4_t = lfc("CLDN4", "T47D")
    cldn4_m = lfc("CLDN4", "MCF7")
    tac_t = lfc("TACSTD2", "T47D")
    tac_m = lfc("TACSTD2", "MCF7")

    text = f"""# FINDING — GSE207704 CLDN4 knockout: NHEJ, STING, IFN/APM scores vs WT

Additive public evidence. This file does not replace the CosMx, concordant-4, GSE137244, TCGA keratin, or TISMO results. It is not a lung dataset and it is not SKB264.

## Design

GSE207704 (PMID 37059993) is CRISPR **CLDN4−/− vs parental WT** in **T47D and MCF7** breast cancer lines. GEO lists 2 biological replicates per genotype (GSM6310640–GSM6310647). The only processed matrix, `GSE207704_CLDN4_RNAseq.txt.gz`, is **group-collapsed cufflinks FPKM** (one column per genotype per line). Sample-level *t* tests and gene-level FDR are not available from this file.

The parental WT arm is described as genotype WT, “No treatment”. That is the negative-control comparator used here. The series has **no siRNA non-targeting arm**. Calling the experiment “KD vs NC” overstates the design: it is a knockout versus parental WT.

Rank and scores use log2((FPKM + 0.5) / (FPKM_WT + 0.5)) after keeping the highest-mean-FPKM locus per symbol. Positive = higher after CLDN4 loss. Per-gene UP/DOWN uses |log2FC| > 0.25 and max FPKM ≥ 1. The gene score is the unweighted mean of log2(FPKM + 0.5) on genes that clear that expression floor. Δ = score(KO) − score(WT).

## QC

| Line | CLDN4 log2FC (KO vs WT) | TACSTD2 log2FC |
|---|---:|---:|
| T47D | {_fmt(cldn4_t)} | {_fmt(tac_t)} |
| MCF7 | {_fmt(cldn4_m)} | {_fmt(tac_m)} |

CLDN4 mRNA falls in both lines but is not abolished (residual FPKM remains). TACSTD2 falls with it. Symbols in the rank: {fpkm.shape[0]}.

## NHEJ core (PRKDC, LIG4, XRCC4, XRCC5, XRCC6, NHEJ1)

All six genes are in the deposit. A consensus of UP means one line clears |log2FC| > 0.25 and the other line is not DOWN. No core NHEJ gene is UP in both lines.

| Gene | T47D log2FC | MCF7 log2FC | mean | consensus |
|---|---|---|---:|---|
{nhej_rows}

{score_line("NHEJ_CORE")}

The 6-gene panel is under the 8-gene GSEA floor. Its NES is reported with a nominal permutation p and **no BH-FDR**.

## STING axis (CGAS, STING1, TBK1, IRF3)

CGAS is deposited as **MB21D1** (GRCh38 v90 symbol). STING1 (TMEM173, Entrez 340061) is **not in the file**, so the 4-gene axis cannot be scored complete. The three measured genes sit near zero.

| Gene | T47D log2FC | MCF7 log2FC | mean | consensus |
|---|---|---|---:|---|
{sting_rows}

{score_line("STING_AXIS")}

## IFN and APM gene scores

Same gene lists as the earlier GSE207704 IFN/MHC-I panel. HLA-A is in the APM score, not the IFN score. Genes missing from the cufflinks table are excluded from the mean, not entered as zero. Many classical APM genes (HLA-A/B, B2M, TAP1/2, PSMB8/9, NLRC5) are absent, so the APM score is only the genes the file actually contains.

{score_line("IFN")}

{score_line("APM")}

## GSEA

Engine: weighted KS *p* = 1, {NPERM} gene-set permutations, seed {SEED}. Positive NES = enriched among genes that rise after CLDN4 loss. BH-FDR is computed **inside the pathway sets** (Reactome plus the IFN and APM panels) for each contrast. Named NHEJ and STING panels are not in that FDR family.

### Pathway sets (size floor 8)

{gsea_md("pathway")}

### Named panels (size floor 3, FDR not applied)

{gsea_md("named_panel")}

Reactome STING (R-HSA-1834941) also contains PRKDC, XRCC5, and XRCC6, so that pathway is not a pure cGAS–STING1–TBK1–IRF3 test. The STING_AXIS score above is the four-gene readout.

## What this accession supports

- **NHEJ core is not a shared program after CLDN4 loss.** T47D Δ = +0.232 is carried by LIG4 (+1.038) and PRKDC (+0.324); the other four T47D genes are FLAT. MCF7 Δ = +0.018 and every MCF7 core gene is FLAT (PRKDC −0.249 sits just under the 0.25 cutoff). XRCC5 is slightly higher in both lines (+0.164, +0.168) and still FLAT. Named-panel NES is +1.253 (nominal p 0.108) in T47D and −0.451 (nominal p 0.445) in MCF7. Reactome NHEJ (29 genes in rank) is not enriched (FDR 0.397 and 0.526).
- **STING axis does not move, and STING1 is missing.** CGAS is MB21D1 in this table and is FLAT (+0.138 T47D, +0.197 MCF7). TBK1 and IRF3 are FLAT. STING1/TMEM173 is absent, so 3 of 4 genes are scored (T47D Δ +0.014, MCF7 Δ +0.125). The 3-gene NES is nominal only. Reactome STING has 8 of 16 genes in the rank, mixes in PRKDC/XRCC5/XRCC6, and is not enriched (FDR 0.448 and 0.526).
- **IFN gene score is lower in both lines.** T47D Δ −0.469 (12 down, 2 up among 28 expressed genes). MCF7 Δ −0.167 (12 down, 5 up among 25). Genes down in both lines include OAS3, USP18, IFI6, and IFI44. HERC5 is the IFN gene up in both lines (+1.108 T47D, +1.924 MCF7). GSEA on the mapped IFN panel: T47D NES −1.842 FDR 0.010; MCF7 NES −1.524 FDR 0.052; mean of lines NES −1.893 FDR 0.007. Reactome interferon alpha/beta signaling is down in both lines (T47D NES −1.794 FDR 0.010; MCF7 NES −1.540 FDR 0.035). The GSEA mean log2FC is a little more negative than the expressed-gene score because low-FPKM genes stay in the rank.
- **APM is not opened.** Nine of 22 requested genes are present. HLA-A/B, B2M, TAP1/2, PSMB8/9, and NLRC5 are absent. Among genes that are present, ERAP1 is down in both lines, HLA-C is discordant (T47D down, MCF7 up), PSMB10 is up only in MCF7, and CALR/CANX/PDIA3/SEC61A1 are flat. APM GSEA FDR is 0.347 (T47D) and 0.526 (MCF7).

Permutation p-values test where the gene set sits on this collapsed rank. They are not tests of the n=2 biological replicates. Replicate-level uncertainty is not estimable from the deposited file.
"""
    (OUT / "FINDING.md").write_text(text)


def plot(genes: pd.DataFrame, scores: pd.DataFrame, gsea: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(13.2, 5.6), constrained_layout=True)

    heat_symbols = ["CLDN4", *NHEJ_CORE, *STING_AXIS]
    mat = []
    ylabels = []
    for symbol in heat_symbols:
        sub = genes[genes["symbol"].eq(symbol)].iloc[0]
        if symbol == "CLDN4":
            ylabels.append("CLDN4 (QC)")
        elif sub["mapped"] and sub["deposited_symbol"] != symbol:
            ylabels.append(f"{symbol}\n({sub['deposited_symbol']})")
        else:
            ylabels.append(symbol + ("" if sub["mapped"] else "\n(absent)"))
        if sub["mapped"]:
            mat.append([sub["log2FC_T47D"], sub["log2FC_MCF7"]])
        else:
            mat.append([np.nan, np.nan])
    mat = np.array(mat, dtype=float)
    ax = axes[0]
    vmax = max(1.0, np.nanmax(np.abs(mat)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks([0, 1], ["T47D", "MCF7"])
    ax.set_yticks(range(len(ylabels)), ylabels, fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            txt = "absent" if not np.isfinite(val) else f"{val:+.2f}"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7, color="black")
    ax.set_title("log2FC, KO vs WT")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="log2FC")

    ax = axes[1]
    set_order = ["NHEJ_CORE", "STING_AXIS", "IFN", "APM"]
    x = np.arange(len(set_order))
    width = 0.36
    for i, line in enumerate(LINES):
        vals = [
            float(scores[(scores["set"].eq(s)) & (scores["line"].eq(line))]["delta_KO_minus_WT"].iloc[0])
            for s in set_order
        ]
        ax.bar(x + (i - 0.5) * width, vals, width=width, label=line)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x, ["NHEJ\n6 genes", "STING\n3 of 4", "IFN", "APM"], fontsize=8)
    ax.set_ylabel("Δ gene score (KO − WT)")
    ax.set_title("Mean log2 score vs parental WT")
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    labeled = []
    for i, line in enumerate(LINES):
        for j, s in enumerate(set_order):
            val = float(scores[(scores["set"].eq(s)) & (scores["line"].eq(line))]["delta_KO_minus_WT"].iloc[0])
            labeled.append(val)
            ax.text(
                x[j] + (i - 0.5) * width,
                val + (0.025 if val >= 0 else -0.025),
                f"{val:+.2f}",
                ha="center",
                va="bottom" if val >= 0 else "top",
                fontsize=6,
                color="#333333",
            )
    ymin = min(labeled) - 0.12
    ymax = max(labeled) + 0.12
    ax.set_ylim(ymin, ymax)

    ax = axes[2]
    path = gsea[gsea["kind"].eq("pathway") & gsea["contrast"].isin(LINES)].copy()
    terms = [
        "REACTOME_NONHOMOLOGOUS_END_JOINING",
        "REACTOME_STING_MEDIATED_INDUCTION",
        "IFN",
        "APM",
        "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
        "REACTOME_INTERFERON_GAMMA_SIGNALING",
        "REACTOME_CLASS_I_MHC_PEPTIDE_LOADING",
    ]
    short = {
        "REACTOME_NONHOMOLOGOUS_END_JOINING": "Reactome NHEJ",
        "REACTOME_STING_MEDIATED_INDUCTION": "Reactome STING",
        "IFN": "IFN panel",
        "APM": "APM panel",
        "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING": "IFN-α/β",
        "REACTOME_INTERFERON_GAMMA_SIGNALING": "IFN-γ",
        "REACTOME_CLASS_I_MHC_PEPTIDE_LOADING": "MHC-I load",
    }
    y = np.arange(len(terms))
    for i, line in enumerate(LINES):
        vals = []
        for term in terms:
            hit = path[(path["term"].eq(term)) & (path["contrast"].eq(line))]
            vals.append(float(hit["nes"].iloc[0]) if len(hit) else np.nan)
        ax.scatter(vals, y + (i - 0.5) * 0.18, label=line, s=28)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_yticks(y, [short[t] for t in terms], fontsize=8)
    ax.set_xlabel("NES (positive = up after KO)")
    ax.set_title("Prerank GSEA")
    ax.legend(frameon=False, fontsize=8)
    fig.suptitle("GSE207704 CLDN4−/− vs parental WT (collapsed FPKM)", fontsize=12)
    fig_dir = OUT / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / "fig_nhej_sting_ifn_scores.png", dpi=160)
    fig.savefig(fig_dir / "fig_nhej_sting_ifn_scores.pdf")
    plt.close(fig)


def main() -> None:
    ensure_matrix()
    fpkm = load_fpkm()
    genes = gene_table(fpkm)
    scores = score_table(genes)
    gsea = run_gsea(fpkm)

    # Locked QC from the same matrix and pseudocount used in the earlier prerank.
    cldn4 = genes[genes["symbol"].eq("CLDN4")].iloc[0]
    assert abs(cldn4["log2FC_T47D"] - (-1.039)) < 0.01, cldn4["log2FC_T47D"]
    assert abs(cldn4["log2FC_MCF7"] - (-0.750)) < 0.01, cldn4["log2FC_MCF7"]
    assert fpkm.shape[0] == 12635, fpkm.shape[0]
    for g in NHEJ_CORE:
        assert bool(genes.loc[genes["symbol"].eq(g), "mapped"].iloc[0]), g

    table_dir = OUT / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    genes.to_csv(table_dir / "gene_stats.tsv", sep="\t", index=False, float_format="%.6g")
    scores.to_csv(table_dir / "gene_scores.tsv", sep="\t", index=False, float_format="%.6g")
    gsea.to_csv(table_dir / "gsea_prerank.tsv", sep="\t", index=False, float_format="%.6g")
    write_finding(fpkm, genes, scores, gsea)
    plot(genes, scores, gsea)
    print(scores[["set", "line", "n_expressed", "n_requested", "delta_KO_minus_WT", "both_lines"]].to_string(index=False))
    print(gsea[["contrast", "term", "nes", "nom_p", "fdr", "n_set_in_rank", "kind"]].to_string(index=False))
    print("wrote", OUT)


if __name__ == "__main__":
    main()
