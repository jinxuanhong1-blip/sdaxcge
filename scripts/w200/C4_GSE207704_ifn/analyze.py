#!/usr/bin/env python3
"""C4 GSE207704 CLDN4 KO: IFN / MHC-I / APM gene directions (honest).

Source: GEO processed cufflinks FPKM (GSE207704_CLDN4_RNAseq.txt).
Design: T47D and MCF7, CLDN4 CRISPR KO vs WT. Deposit is group-collapsed
(n=2 pooled per genotype per line). No gene-level p-value/FDR is valid.

log2FC = log2((KO + 0.5) / (WT + 0.5)). Positive = higher in KO.
"""
from __future__ import annotations

import csv
import gzip
import io
import math
import os
import urllib.request
from collections import defaultdict

GEO_URL = (
    "https://www.ncbi.nlm.nih.gov/geo/download/"
    "?acc=GSE207704&format=file&file=GSE207704_CLDN4_RNAseq.txt.gz"
)
PSEUDO = 0.5
LOW_FPKM = 1.0
WEAK = 0.25

# Entrez used only as a second lookup; many C4 genes are missing from the deposit.
PANEL = [
    # QC
    ("CLDN4", "QC", "1364"),
    ("TACSTD2", "QC", "4070"),
    # C4 priority IFN / MHC-I
    ("IFI27", "priority_ifn_mhci", "3429"),
    ("OAS2", "priority_ifn_mhci", "4939"),
    ("IFIT1", "priority_ifn_mhci", "3434"),
    ("MX1", "priority_ifn_mhci", "4599"),
    ("ISG15", "priority_ifn_mhci", "9636"),
    ("HLA-A", "priority_ifn_mhci", "3105"),
    # Type I ISGs
    ("IFIT2", "ifn_type1", "3433"),
    ("IFIT3", "ifn_type1", "3437"),
    ("IFIT5", "ifn_type1", "24138"),
    ("MX2", "ifn_type1", "4600"),
    ("OAS1", "ifn_type1", "4938"),
    ("OAS3", "ifn_type1", "4940"),
    ("OASL", "ifn_type1", "8638"),
    ("RSAD2", "ifn_type1", "91543"),
    ("USP18", "ifn_type1", "11274"),
    ("IFI6", "ifn_type1", "2537"),
    ("IFI35", "ifn_type1", "3430"),
    ("IFI44", "ifn_type1", "10561"),
    ("IFI44L", "ifn_type1", "10964"),
    ("IFI16", "ifn_type1", "3428"),
    ("BST2", "ifn_type1", "684"),
    ("IFITM1", "ifn_type1", "8519"),
    ("IFITM2", "ifn_type1", "10581"),
    ("IFITM3", "ifn_type1", "10410"),
    ("DDX58", "ifn_type1", "23586"),
    ("IFIH1", "ifn_type1", "64135"),
    ("DDX60", "ifn_type1", "55601"),
    ("HERC5", "ifn_type1", "51191"),
    ("XAF1", "ifn_type1", "54739"),
    ("EIF2AK2", "ifn_type1", "5610"),
    ("LY6E", "ifn_type1", "4061"),
    ("CMPK2", "ifn_type1", "129607"),
    ("EPSTI1", "ifn_type1", "94240"),
    ("SAMD9L", "ifn_type1", "219285"),
    ("TRIM22", "ifn_type1", "10346"),
    ("ZBP1", "ifn_type1", "81030"),
    ("PLSCR1", "ifn_type1", "5359"),
    # Signaling
    ("STAT1", "ifn_signaling", "6772"),
    ("STAT2", "ifn_signaling", "6773"),
    ("IRF1", "ifn_signaling", "3659"),
    ("IRF7", "ifn_signaling", "3665"),
    ("IRF9", "ifn_signaling", "10379"),
    ("JAK1", "ifn_signaling", "3716"),
    ("JAK2", "ifn_signaling", "3717"),
    ("IFNAR1", "ifn_signaling", "3454"),
    ("IFNAR2", "ifn_signaling", "3455"),
    ("IFNGR1", "ifn_signaling", "3459"),
    ("IFNGR2", "ifn_signaling", "3460"),
    ("SOCS1", "ifn_signaling", "8651"),
    ("SOCS3", "ifn_signaling", "9021"),
    # MHC-I / APM
    ("HLA-B", "mhc1_apm", "3106"),
    ("HLA-C", "mhc1_apm", "3107"),
    ("HLA-E", "mhc1_apm", "3133"),
    ("HLA-F", "mhc1_apm", "3134"),
    ("HLA-G", "mhc1_apm", "3135"),
    ("B2M", "mhc1_apm", "567"),
    ("NLRC5", "mhc1_apm", "84166"),
    ("TAP1", "mhc1_apm", "6890"),
    ("TAP2", "mhc1_apm", "6891"),
    ("TAPBP", "mhc1_apm", "6892"),
    ("PSMB8", "mhc1_apm", "5696"),
    ("PSMB9", "mhc1_apm", "5698"),
    ("PSMB10", "mhc1_apm", "5699"),
    ("ERAP1", "mhc1_apm", "51752"),
    ("ERAP2", "mhc1_apm", "64167"),
    ("CALR", "mhc1_apm", "811"),
    ("CANX", "mhc1_apm", "821"),
    ("PDIA3", "mhc1_apm", "2923"),
    ("SEC61A1", "mhc1_apm", "29927"),
    ("PSME1", "mhc1_apm", "5720"),
    ("PSME2", "mhc1_apm", "5721"),
    # Type II / chemokines
    ("GBP1", "ifn_type2", "2633"),
    ("GBP2", "ifn_type2", "2634"),
    ("GBP4", "ifn_type2", "115361"),
    ("GBP5", "ifn_type2", "115362"),
    ("CXCL9", "ifn_type2", "4283"),
    ("CXCL10", "ifn_type2", "3627"),
    ("CXCL11", "ifn_type2", "6373"),
    ("CIITA", "ifn_type2", "4261"),
    ("IDO1", "ifn_type2", "3620"),
    ("IFNG", "ifn_type2", "3458"),
    ("IFNB1", "ifn_type2", "3456"),
]

COLS = {
    "mcf7_ko": "MCF7_CLDN4KO_FPKM (fpkm)",
    "mcf7_wt": "MCF7_WT_FPKM (fpkm)",
    "t47d_ko": "T47D_CLDN4KO_FPKM (fpkm)",
    "t47d_wt": "T47D_WT_FPKM (fpkm)",
}


def log2fc(ko: float, wt: float, pseudo: float = PSEUDO) -> float:
    return math.log2((ko + pseudo) / (wt + pseudo))


def direction(lfc: float, max_fpkm: float) -> str:
    if max_fpkm < LOW_FPKM:
        return "LOW"
    if lfc > WEAK:
        return "UP"
    if lfc < -WEAK:
        return "DOWN"
    return "FLAT"


def consensus(d_mcf7: str, d_t47d: str, mean_lfc: float) -> str:
    if d_mcf7 == "ABSENT" or d_t47d == "ABSENT":
        return "ABSENT"
    dirs = {d_mcf7, d_t47d}
    if dirs <= {"LOW"}:
        return "LOW"
    if "UP" in dirs and "DOWN" in dirs:
        return "DISCORDANT"
    if dirs <= {"UP", "FLAT", "LOW"} and (d_mcf7 == "UP" or d_t47d == "UP" or mean_lfc > WEAK):
        if "UP" in dirs:
            return "UP"
    if dirs <= {"DOWN", "FLAT", "LOW"} and (d_mcf7 == "DOWN" or d_t47d == "DOWN" or mean_lfc < -WEAK):
        if "DOWN" in dirs:
            return "DOWN"
    return "FLAT"


def load_matrix(path: str) -> list[dict]:
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        # Duplicate trailing gene_id column: keep Entrez from the last field.
        raw = fh.read()
    reader = csv.reader(io.StringIO(raw), delimiter="\t")
    header = next(reader)
    rows = []
    for parts in reader:
        if len(parts) < 10:
            continue
        rows.append(
            {
                "tracking_id": parts[0],
                "cuff_id": parts[1],
                "gene_short_name": parts[2],
                "mcf7_ko": float(parts[5] or 0),
                "mcf7_wt": float(parts[6] or 0),
                "t47d_ko": float(parts[7] or 0),
                "t47d_wt": float(parts[8] or 0),
                "refseq": parts[9] if len(parts) > 9 else "",
                "entrez": parts[10] if len(parts) > 10 else "",
            }
        )
    return rows


def index_rows(rows: list[dict]):
    by_name = defaultdict(list)
    by_entrez = defaultdict(list)
    for r in rows:
        name = (r["gene_short_name"] or "").strip().upper()
        if name:
            by_name[name].append(r)
        entrez = (r["entrez"] or "").strip()
        if entrez.isdigit():
            by_entrez[entrez].append(r)
    return by_name, by_entrez


def collapse(hits: list[dict]) -> dict | None:
    if not hits:
        return None
    keys = ["mcf7_ko", "mcf7_wt", "t47d_ko", "t47d_wt"]
    # Prefer the cufflinks locus with the highest total FPKM (canonical transcript).
    best = max(hits, key=lambda h: sum(h[k] for k in keys))
    out = {k: best[k] for k in keys}
    out["n_loci"] = len(hits)
    out["loci"] = best["tracking_id"]
    return out


def analyze(rows: list[dict]) -> list[dict]:
    by_name, by_entrez = index_rows(rows)
    out = []
    for symbol, group, entrez in PANEL:
        hits = by_name.get(symbol.upper(), [])
        mapped = "symbol"
        if not hits:
            hits = by_entrez.get(entrez, [])
            mapped = "entrez" if hits else "absent"
        collapsed = collapse(hits)
        rec = {
            "symbol": symbol,
            "group": group,
            "entrez": entrez,
            "mapping": mapped,
            "n_loci": 0,
            "MCF7_WT": "",
            "MCF7_KO": "",
            "T47D_WT": "",
            "T47D_KO": "",
            "log2FC_MCF7": "",
            "log2FC_T47D": "",
            "log2FC_mean": "",
            "dir_MCF7": "ABSENT",
            "dir_T47D": "ABSENT",
            "consensus": "ABSENT",
            "note": "not in deposited cufflinks table (annotation gap, not a measured zero)",
        }
        if collapsed:
            m_wt, m_ko = collapsed["mcf7_wt"], collapsed["mcf7_ko"]
            t_wt, t_ko = collapsed["t47d_wt"], collapsed["t47d_ko"]
            lfc_m = log2fc(m_ko, m_wt)
            lfc_t = log2fc(t_ko, t_wt)
            lfc_mean = 0.5 * (lfc_m + lfc_t)
            rec.update(
                {
                    "n_loci": collapsed["n_loci"],
                    "MCF7_WT": f"{m_wt:.4g}",
                    "MCF7_KO": f"{m_ko:.4g}",
                    "T47D_WT": f"{t_wt:.4g}",
                    "T47D_KO": f"{t_ko:.4g}",
                    "log2FC_MCF7": f"{lfc_m:.4f}",
                    "log2FC_T47D": f"{lfc_t:.4f}",
                    "log2FC_mean": f"{lfc_mean:.4f}",
                    "dir_MCF7": direction(lfc_m, max(m_wt, m_ko)),
                    "dir_T47D": direction(lfc_t, max(t_wt, t_ko)),
                    "note": "descriptive only; collapsed FPKM, no gene-level FDR",
                }
            )
            rec["consensus"] = consensus(rec["dir_MCF7"], rec["dir_T47D"], lfc_mean)
        out.append(rec)
    return out


def write_tsv(path: str, rows: list[dict], fields: list[str]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def summarize(rows: list[dict]) -> list[dict]:
    groups = []
    seen = []
    for r in rows:
        if r["group"] not in seen:
            seen.append(r["group"])
    for g in seen:
        sub = [r for r in rows if r["group"] == g]
        measured = [r for r in sub if r["consensus"] != "ABSENT"]
        groups.append(
            {
                "group": g,
                "n_panel": len(sub),
                "n_measured": len(measured),
                "n_absent": sum(r["consensus"] == "ABSENT" for r in sub),
                "n_up": sum(r["consensus"] == "UP" for r in sub),
                "n_down": sum(r["consensus"] == "DOWN" for r in sub),
                "n_discordant": sum(r["consensus"] == "DISCORDANT" for r in sub),
                "n_flat": sum(r["consensus"] == "FLAT" for r in sub),
                "n_low": sum(r["consensus"] == "LOW" for r in sub),
                "up_genes": ",".join(r["symbol"] for r in sub if r["consensus"] == "UP") or ".",
                "down_genes": ",".join(r["symbol"] for r in sub if r["consensus"] == "DOWN") or ".",
                "discordant_genes": ",".join(r["symbol"] for r in sub if r["consensus"] == "DISCORDANT") or ".",
                "absent_genes": ",".join(r["symbol"] for r in sub if r["consensus"] == "ABSENT") or ".",
            }
        )
    return groups


def main() -> None:
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
    outdir = os.path.join(root, "results", "w200", "C4_GSE207704_ifn")
    raw_dir = os.path.join(outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)
    raw_gz = os.path.join(raw_dir, "GSE207704_CLDN4_RNAseq.txt.gz")
    if not os.path.exists(raw_gz):
        print(f"downloading {GEO_URL}")
        urllib.request.urlretrieve(GEO_URL, raw_gz)
    rows = load_matrix(raw_gz)
    stats = analyze(rows)
    fields = [
        "symbol",
        "group",
        "entrez",
        "mapping",
        "n_loci",
        "MCF7_WT",
        "MCF7_KO",
        "T47D_WT",
        "T47D_KO",
        "log2FC_MCF7",
        "log2FC_T47D",
        "log2FC_mean",
        "dir_MCF7",
        "dir_T47D",
        "consensus",
        "note",
    ]
    write_tsv(os.path.join(outdir, "gene_stats.tsv"), stats, fields)
    write_tsv(
        os.path.join(outdir, "group_summary.tsv"),
        summarize(stats),
        [
            "group",
            "n_panel",
            "n_measured",
            "n_absent",
            "n_up",
            "n_down",
            "n_discordant",
            "n_flat",
            "n_low",
            "up_genes",
            "down_genes",
            "discordant_genes",
            "absent_genes",
        ],
    )
    # Simple up vs down lists (IFN/MHC-I/APM only; exclude QC)
    immuno = [r for r in stats if r["group"] != "QC"]
    lists = [
        {
            "list": "UP_in_KO_consensus",
            "definition": "same-direction UP (log2FC>+0.25, max FPKM>=1) or UP+FLAT; both lines not opposite",
            "genes": ",".join(r["symbol"] for r in immuno if r["consensus"] == "UP") or ".",
        },
        {
            "list": "DOWN_in_KO_consensus",
            "definition": "same-direction DOWN (log2FC<-0.25, max FPKM>=1) or DOWN+FLAT; both lines not opposite",
            "genes": ",".join(r["symbol"] for r in immuno if r["consensus"] == "DOWN") or ".",
        },
        {
            "list": "DISCORDANT",
            "definition": "one line UP and the other DOWN",
            "genes": ",".join(r["symbol"] for r in immuno if r["consensus"] == "DISCORDANT") or ".",
        },
        {
            "list": "FLAT",
            "definition": "|log2FC|<=0.25 in both lines or no qualifying UP/DOWN",
            "genes": ",".join(r["symbol"] for r in immuno if r["consensus"] == "FLAT") or ".",
        },
        {
            "list": "ABSENT_from_deposit",
            "definition": "symbol and Entrez missing from GSE207704_CLDN4_RNAseq.txt",
            "genes": ",".join(r["symbol"] for r in immuno if r["consensus"] == "ABSENT") or ".",
        },
        {
            "list": "LOW_expression",
            "definition": "max FPKM < 1 in the line(s) that would otherwise be called",
            "genes": ",".join(r["symbol"] for r in immuno if r["consensus"] == "LOW") or ".",
        },
    ]
    write_tsv(
        os.path.join(outdir, "up_vs_down.tsv"),
        lists,
        ["list", "definition", "genes"],
    )
    print(f"wrote {outdir}")
    for rec in lists:
        print(f"{rec['list']}: {rec['genes']}")


if __name__ == "__main__":
    main()
