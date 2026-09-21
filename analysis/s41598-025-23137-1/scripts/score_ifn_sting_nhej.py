#!/usr/bin/env python3
"""Score IFN, STING, and NHEJ genes on open matrices related to
Sci Rep 2025 s41598-025-23137-1.

The 2025 CRISPRi CLDN4 experiment has no deposited RNA-seq or proteome.
This script scores two open resources that the search did return:

1. TCGA OV PanCancer Atlas (cBioPortal study ov_tcga_pan_can_atlas_2018),
   the public cohort cited by the paper. RNA-seq is RSEM; protein is the
   CPTAC mass-spec profile in the same study. Spearman is against CLDN4.
2. GSE22493, the only GEO series whose design is CLDN4 silencing in an
   ovarian line. It is a 2010 SKOV-3 two-color array (Gao), not the 2025
   OVCAR3/OVCA429 CRISPRi experiment.

Raw matrices are downloaded to a cache directory and are not written next
to the result tables.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import time
import urllib.request
from pathlib import Path

from scipy import stats

STUDY = "ov_tcga_pan_can_atlas_2018"
RNA_PROFILE = f"{STUDY}_rna_seq_v2_mrna"
PROT_PROFILE = f"{STUDY}_protein_quantification"
CBIO = "https://www.cbioportal.org/api"

# Aliases are queried separately. TCGA PanCancer RNA uses the 2018 symbol
# (TMEM173, MB21D1); a missing modern symbol is recorded as not in the profile.
PANELS = {
    "IFN_ISG": [
        "IFNB1",
        "IFNA1",
        "ISG15",
        "MX1",
        "MX2",
        "OAS1",
        "OAS2",
        "OAS3",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "EIF2AK2",
        "IRF7",
        "IRF9",
        "STAT1",
        "STAT2",
        "RSAD2",
        "BST2",
        "IFI44",
        "IFI44L",
        "IFI6",
        "IFITM1",
        "CXCL10",
        "CCL5",
        "ISG20",
        "USP18",
        "OASL",
        "RIGI",  # cBioPortal symbol for DDX58
        "DDX58",
        "IFIH1",
        "MAVS",
        "SOCS1",
    ],
    "STING": [
        "CGAS",
        "MB21D1",
        "STING1",
        "TMEM173",
        "TBK1",
        "IKBKE",
        "IRF3",
        "TREX1",
        "ENPP1",
        "RAB7A",
        "BECN1",
        "MAP1LC3B",
        "ATG5",
        "ULK1",
    ],
    "NHEJ": [
        "XRCC4",
        "XRCC5",
        "XRCC6",
        "LIG4",
        "PRKDC",
        "NHEJ1",
        "DCLRE1C",
        "PAXX",
        "C9orf142",
        "TP53BP1",
        "POLL",
        "POLM",
        "PNKP",
        "APTX",
        "RIF1",
        "MAD2L2",
        "SHLD1",
        "SHLD2",
        "FAM35A",
        "SHLD3",
    ],
}

PAPER_ISG = [
    "ISG15",
    "MX1",
    "OAS1",
    "IFIT1",
    "EIF2AK2",
    "IRF7",
    "STAT1",
    "RSAD2",
    "BST2",
    "IFI44",
]


def bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    order = sorted(range(m), key=lambda i: (math.inf if math.isnan(pvals[i]) else pvals[i], i))
    q = [math.nan] * m
    running = 1.0
    for rank_from_end, i in enumerate(reversed(order)):
        p = pvals[i]
        if math.isnan(p):
            continue
        rank = m - rank_from_end
        val = min(running, p * m / rank)
        running = val
        q[i] = min(1.0, val)
    return q


def spearman(x: list[float], y: list[float]) -> tuple[float, float]:
    if len(x) < 3:
        return math.nan, math.nan
    rho, p = stats.spearmanr(x, y)
    return float(rho), float(p)


def http_json(url: str, payload: dict | list | None = None, attempts: int = 4):
    data = None
    headers = {"Accept": "application/json", "User-Agent": "cldn4-inventory/1.0"}
    if payload is not None:
        data = json.dumps(payload).encode()
        headers["Content-Type"] = "application/json"
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(url, data=data, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode())
        except Exception as exc:  # noqa: BLE001 — retry transient API errors
            last = exc
            time.sleep(1.5 * (i + 1))
    raise RuntimeError(f"GET/POST failed {url}: {last}")


def fetch_entrez(symbols: list[str]) -> dict[str, int]:
    found: dict[str, int] = {}
    # cBioPortal accepts HUGO symbols in chunks.
    for start in range(0, len(symbols), 40):
        chunk = symbols[start : start + 40]
        rows = http_json(
            f"{CBIO}/genes/fetch?geneIdType=HUGO_GENE_SYMBOL",
            chunk,
        )
        for row in rows:
            found[row["hugoGeneSymbol"]] = int(row["entrezGeneId"])
    return found


def fetch_profile(profile: str, entrez: int, cache: Path) -> dict[str, float]:
    dest = cache / f"{profile}_{entrez}.json"
    if dest.exists():
        rows = json.loads(dest.read_text())
    else:
        url = (
            f"{CBIO}/molecular-profiles/{profile}/molecular-data"
            f"?entrezGeneId={entrez}&sampleListId={STUDY}_all"
            f"&projection=SUMMARY&pageSize=100000&pageNumber=0"
        )
        rows = http_json(url)
        dest.write_text(json.dumps(rows))
        time.sleep(0.12)
    out = {}
    for row in rows:
        val = row.get("value")
        if val is None:
            continue
        out[row["sampleId"]] = float(val)
    return out


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(cols)]
    for row in rows:
        cells = []
        for c in cols:
            v = row.get(c, "")
            if isinstance(v, float):
                if math.isnan(v):
                    cells.append("")
                else:
                    cells.append(f"{v:.6g}")
            else:
                cells.append(str(v))
        lines.append("\t".join(cells))
    path.write_text("\n".join(lines) + "\n")


def score_profile(
    profile: str,
    anchor_symbol: str,
    symbols: list[tuple[str, str]],
    entrez: dict[str, int],
    cache: Path,
    transform,
) -> tuple[list[dict], dict[str, dict[str, float]]]:
    anchor_id = entrez[anchor_symbol]
    anchor = {s: transform(v) for s, v in fetch_profile(profile, anchor_id, cache).items()}
    stored = {anchor_symbol: anchor}
    rows = []
    for panel, symbol in symbols:
        eid = entrez.get(symbol)
        if eid is None:
            rows.append(
                {
                    "profile": profile,
                    "panel": panel,
                    "gene": symbol,
                    "entrez": "",
                    "n": 0,
                    "rho_vs_CLDN4": math.nan,
                    "p": math.nan,
                    "status": "symbol_not_in_cbioportal",
                }
            )
            continue
        values = fetch_profile(profile, eid, cache)
        stored[symbol] = {s: transform(v) for s, v in values.items()}
        paired_x = []
        paired_y = []
        for sample, y in stored[symbol].items():
            if sample in anchor:
                paired_x.append(anchor[sample])
                paired_y.append(y)
        if symbol == anchor_symbol:
            rho, p = 1.0, 0.0
            status = "anchor"
        elif len(paired_x) < 3:
            rho, p = math.nan, math.nan
            status = "too_few_paired_samples"
        else:
            rho, p = spearman(paired_x, paired_y)
            status = "ok"
        rows.append(
            {
                "profile": profile,
                "panel": panel,
                "gene": symbol,
                "entrez": eid,
                "n": len(paired_x),
                "rho_vs_CLDN4": rho,
                "p": p,
                "status": status,
            }
        )
    # Modern and 2018 symbols can resolve to one Entrez id. Keep the first
    # symbol in panel order for the FDR and drop the alias from the signature.
    seen_entrez: set[int] = set()
    for r in rows:
        if r["status"] != "ok":
            continue
        eid = int(r["entrez"])
        if eid in seen_entrez:
            r["status"] = "alias_same_entrez"
            r["p"] = math.nan
            stored.pop(r["gene"], None)
        else:
            seen_entrez.add(eid)
    testable = [r for r in rows if r["status"] == "ok"]
    qs = bh([r["p"] for r in testable])
    for r, q in zip(testable, qs):
        r["q_bh"] = q
    for r in rows:
        r.setdefault("q_bh", math.nan)
    return rows, stored


def zscore(vals: dict[str, float]) -> dict[str, float]:
    xs = list(vals.values())
    n = len(xs)
    if n < 3:
        return {}
    mu = sum(xs) / n
    var = sum((x - mu) ** 2 for x in xs) / (n - 1)
    if var <= 0:
        return {}
    sd = math.sqrt(var)
    return {s: (v - mu) / sd for s, v in vals.items()}


def panel_signatures(stored: dict[str, dict[str, float]], panel_of: dict[str, str], anchor: str) -> list[dict]:
    anchor_vals = stored[anchor]
    out = []
    panels = sorted(set(panel_of.values()))
    for panel in panels:
        genes = [g for g, p in panel_of.items() if p == panel and g in stored and g != anchor]
        # Drop alias rows that are empty.
        genes = [g for g in genes if stored[g]]
        zmaps = {g: zscore(stored[g]) for g in genes}
        samples = []
        sig = []
        cl = []
        for sample, cval in anchor_vals.items():
            zs = [zmaps[g][sample] for g in genes if sample in zmaps[g]]
            if len(zs) < max(3, len(genes) // 2):
                continue
            samples.append(sample)
            sig.append(sum(zs) / len(zs))
            cl.append(cval)
        rho, p = spearman(cl, sig)
        out.append(
            {
                "panel": panel,
                "n_genes_in_mean": len(genes),
                "n_samples": len(samples),
                "rho_signature_vs_CLDN4": rho,
                "p": p,
            }
        )
    if out:
        qs = bh([r["p"] for r in out])
        for r, q in zip(out, qs):
            r["q_bh"] = q
    return out


def load_platform(path: Path) -> dict[str, str]:
    mapping = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!platform_table_begin"):
                in_table = True
                continue
            if line.startswith("!platform_table_end"):
                break
            if not in_table:
                continue
            parts = line.rstrip("\n").split("\t")
            if not parts or parts[0] in {"ID", "ID_REF"}:
                continue
            probe = parts[0]
            desc = parts[1] if len(parts) > 1 else ""
            orf = parts[2].strip() if len(parts) > 2 else ""
            symbol = (orf or desc.split("--")[0]).strip()
            if not symbol:
                continue
            symbol = symbol.split()[0].strip().upper()
            if symbol:
                mapping[probe] = symbol
    return mapping


def load_series_matrix(path: Path) -> dict[str, list[float]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    probes: dict[str, list[float]] = {}
    with opener(path, "rt", errors="replace") as handle:
        in_table = False
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                in_table = True
                next(handle)
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            if not in_table:
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            if len(parts) < 4:
                continue
            vals = []
            for cell in parts[1:4]:
                if cell in {"", "null", "NA"}:
                    vals.append(math.nan)
                else:
                    vals.append(float(cell))
            probes[parts[0]] = vals
    return probes


def score_gse22493(matrix_path: Path, platform_path: Path) -> list[dict]:
    platform = load_platform(platform_path)
    matrix = load_series_matrix(matrix_path)
    wanted = {"CLDN4": "anchor"}
    for panel, genes in PANELS.items():
        for gene in genes:
            wanted.setdefault(gene, panel)
    by_gene: dict[str, list[list[float]]] = {}
    for probe, vals in matrix.items():
        symbol = platform.get(probe, "").upper()
        if symbol in wanted:
            by_gene.setdefault(symbol, []).append(vals)
    rows = []
    for symbol, panel in wanted.items():
        probes = by_gene.get(symbol, [])
        if not probes:
            rows.append(
                {
                    "dataset": "GSE22493",
                    "panel": panel,
                    "gene": symbol,
                    "n_probes": 0,
                    "mean_log_ratio_kd_over_oe": math.nan,
                    "sd": math.nan,
                    "t_p_n3": math.nan,
                    "status": "not_on_array",
                }
            )
            continue
        # Gene-level value: median across probes of each replicate, then mean.
        per_rep = []
        for rep in range(3):
            xs = [p[rep] for p in probes if not math.isnan(p[rep])]
            if not xs:
                per_rep.append(math.nan)
                continue
            xs.sort()
            mid = len(xs) // 2
            med = xs[mid] if len(xs) % 2 else 0.5 * (xs[mid - 1] + xs[mid])
            per_rep.append(med)
        finite = [v for v in per_rep if not math.isnan(v)]
        mean = sum(finite) / len(finite)
        if len(finite) >= 2:
            var = sum((v - mean) ** 2 for v in finite) / (len(finite) - 1)
            sd = math.sqrt(var)
        else:
            sd = math.nan
        if len(finite) == 3 and sd > 0:
            tstat = mean / (sd / math.sqrt(3))
            p = float(2 * stats.t.sf(abs(tstat), 2))
        elif len(finite) == 3 and sd == 0:
            p = 0.0 if mean != 0 else 1.0
        else:
            p = math.nan
        rows.append(
            {
                "dataset": "GSE22493",
                "panel": panel,
                "gene": symbol,
                "n_probes": len(probes),
                "rep1": per_rep[0],
                "rep2": per_rep[1],
                "rep3": per_rep[2],
                "mean_log_ratio_kd_over_oe": mean,
                "sd": sd,
                "t_p_n3": p,
                "status": "ok",
            }
        )
    testable = [r for r in rows if r["status"] == "ok" and r["panel"] != "anchor"]
    qs = bh([r["t_p_n3"] for r in testable])
    for r, q in zip(testable, qs):
        r["q_bh"] = q
    for r in rows:
        r.setdefault("q_bh", math.nan)
        r.setdefault("rep1", math.nan)
        r.setdefault("rep2", math.nan)
        r.setdefault("rep3", math.nan)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cache",
        type=Path,
        default=Path("/tmp/nhej-sting/cbio_cache"),
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "tables",
    )
    parser.add_argument(
        "--gse-matrix",
        type=Path,
        default=Path("/tmp/nhej-sting/geo/GSE22493_series_matrix.txt.gz"),
    )
    parser.add_argument(
        "--gse-platform",
        type=Path,
        default=Path("/tmp/nhej-sting/geo/GPL10555_family.soft.gz"),
    )
    args = parser.parse_args()
    args.cache.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    symbols = ["CLDN4"]
    panel_of = {}
    pairs = []
    for panel, genes in PANELS.items():
        for gene in genes:
            symbols.append(gene)
            panel_of[gene] = panel
            pairs.append((panel, gene))
    # unique, keep order
    seen = set()
    uniq = []
    for s in symbols:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    entrez = fetch_entrez(uniq)
    (args.outdir / "cbioportal_entrez.json").write_text(json.dumps(entrez, indent=2) + "\n")

    rna_rows, rna_stored = score_profile(
        RNA_PROFILE,
        "CLDN4",
        pairs,
        entrez,
        args.cache,
        transform=lambda v: math.log2(v + 1.0),
    )
    prot_rows, prot_stored = score_profile(
        PROT_PROFILE,
        "CLDN4",
        pairs,
        entrez,
        args.cache,
        transform=lambda v: v,
    )
    cols = ["profile", "panel", "gene", "entrez", "n", "rho_vs_CLDN4", "p", "q_bh", "status"]
    write_tsv(args.outdir / "tcga_ov_rnaseq_cldn4_spearman.tsv", rna_rows, cols)
    write_tsv(args.outdir / "tcga_ov_cptac_protein_cldn4_spearman.tsv", prot_rows, cols)

    sig_rows = panel_signatures(rna_stored, panel_of, "CLDN4")
    for row in sig_rows:
        row["profile"] = RNA_PROFILE
        row["transform"] = "mean of per-gene z-scores of log2(RSEM+1)"
    sig_prot = panel_signatures(prot_stored, panel_of, "CLDN4")
    for row in sig_prot:
        row["profile"] = PROT_PROFILE
        row["transform"] = "mean of per-gene z-scores of CPTAC protein values"
    write_tsv(
        args.outdir / "panel_signature_vs_cldn4.tsv",
        sig_rows + sig_prot,
        ["profile", "panel", "n_genes_in_mean", "n_samples", "rho_signature_vs_CLDN4", "p", "q_bh", "transform"],
    )

    if args.gse_matrix.exists() and args.gse_platform.exists():
        gse_rows = score_gse22493(args.gse_matrix, args.gse_platform)
        write_tsv(
            args.outdir / "gse22493_kd_over_oe_logratio.tsv",
            gse_rows,
            [
                "dataset",
                "panel",
                "gene",
                "n_probes",
                "rep1",
                "rep2",
                "rep3",
                "mean_log_ratio_kd_over_oe",
                "sd",
                "t_p_n3",
                "q_bh",
                "status",
            ],
        )
    else:
        print("GSE22493 files missing; skipped", args.gse_matrix, args.gse_platform)

    print("wrote", args.outdir)


if __name__ == "__main__":
    main()
