#!/usr/bin/env python3
"""Download the three extra public LUAD RNA cohorts for the A9 TJ supplement.

All sources are open. Nothing here is dbGaP / EGA / controlled-access.

  OncoSG LUAD 2020 (Chen et al.) via cBioPortal
    molecular profile: luad_oncosg_2020_rna_seq_v2_mrna_median_all_sample_Zscores
    sample list:       luad_oncosg_2020_rna_seq_v2_mrna
    unit: z-score of log RNA-seq V2 RSEM (Spearman is rank-invariant)

  CPTAC LUAD RNA, LinkedOmics / CPTAC pancancer freeze v1.2
    LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt
    https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/data_freeze_v1.2_reorganized/LUAD/

  GSE31210 Okayama et al. 2012, Japanese stage I–II LUAD, Affymetrix U133 Plus 2.0
    series matrix + GPL570.annot
    Tumors only (`tissue: primary lung tumor`). Normals are dropped.

GSE72094 is Moffitt (US) LUAD, not East-Asian, and is not downloaded.

Raw matrices are written under results/.../data/raw/ and are gitignored.
Compact TJ+TACSTD2 extracts are written next to them for the analyzer.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from genesets import CURRENT_SYMBOL, TARGET, TJ_UNIVERSE

HGNC_URL = (
    "https://www.genenames.org/cgi-bin/download/custom?"
    "col=gd_app_sym&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id"
    "&status=Approved&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"
)

CBIO = "https://www.cbioportal.org/api"
ONCOSG_STUDY = "luad_oncosg_2020"
ONCOSG_PROFILE = f"{ONCOSG_STUDY}_rna_seq_v2_mrna_median_all_sample_Zscores"
ONCOSG_SAMPLES = f"{ONCOSG_STUDY}_rna_seq_v2_mrna"

CPTAC_URL = (
    "https://cptac-pancancer-data.s3.us-west-2.amazonaws.com/"
    "data_freeze_v1.2_reorganized/LUAD/"
    "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
)

GSE31210_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31210/"
    "matrix/GSE31210_series_matrix.txt.gz"
)
GPL570_ANNOT = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

UA = (
    "sdaxcge-w200-A9-extra-tj/1.0 "
    "(reproducible public LUAD extract; +https://github.com/jinxuanhong1-blip/sdaxcge)"
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def http_json(method: str, url: str, **kw):
    import requests

    last = None
    for attempt in range(5):
        try:
            r = requests.request(method, url, timeout=180, **kw)
            r.raise_for_status()
            return r.json()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2**attempt)
    raise last


def download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    print(f"GET {url}\n -> {dest}", flush=True)
    with urllib.request.urlopen(req, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    tmp.replace(dest)


def catalog_query_symbols() -> list[str]:
    """Symbols to send to APIs, including current names for retired catalog symbols."""
    out = [TARGET]
    for g in TJ_UNIVERSE:
        out.append(g)
        if g in CURRENT_SYMBOL:
            out.append(CURRENT_SYMBOL[g])
    return list(dict.fromkeys(out))


def to_catalog_symbol(symbol: str) -> str:
    inv = {v: k for k, v in CURRENT_SYMBOL.items()}
    return inv.get(symbol, symbol)


def fetch_oncosg(out: Path) -> dict:
    genes = catalog_query_symbols()
    g = http_json(
        "POST",
        f"{CBIO}/genes/fetch",
        params={"geneIdType": "HUGO_GENE_SYMBOL"},
        json=genes,
    )
    entrez2hugo = {x["entrezGeneId"]: to_catalog_symbol(x["hugoGeneSymbol"]) for x in g}
    looked_up = {to_catalog_symbol(x["hugoGeneSymbol"]) for x in g}
    missing_lookup = sorted(set([TARGET, *TJ_UNIVERSE]) - looked_up)
    data = http_json(
        "POST",
        f"{CBIO}/molecular-profiles/{ONCOSG_PROFILE}/molecular-data/fetch",
        params={"projection": "SUMMARY"},
        json={"entrezGeneIds": list(entrez2hugo.keys()), "sampleListId": ONCOSG_SAMPLES},
    )
    recs = [
        (entrez2hugo[d["entrezGeneId"]], d["sampleId"], float(d["value"]))
        for d in data
        if d.get("value") is not None and d["entrezGeneId"] in entrez2hugo
    ]
    df = pd.DataFrame(recs, columns=["gene", "sample", "value"])
    expr = df.pivot_table(index="gene", columns="sample", values="value")
    expr.to_csv(out / "oncosg_tj_expr.tsv", sep="\t")
    clin = http_json(
        "GET",
        f"{CBIO}/studies/{ONCOSG_STUDY}/clinical-data",
        params={"clinicalDataType": "SAMPLE", "projection": "SUMMARY"},
    )
    purity = pd.Series(
        {
            d["sampleId"]: float(d["value"])
            for d in clin
            if d.get("clinicalAttributeId") == "PURITY" and d.get("value") not in (None, "", "NA")
        },
        name="PURITY",
    )
    purity.to_csv(out / "oncosg_purity.tsv", sep="\t", header=True)
    present = sorted(expr.index.astype(str))
    return {
        "study": ONCOSG_STUDY,
        "profile": ONCOSG_PROFILE,
        "sample_list": ONCOSG_SAMPLES,
        "unit": "z-score of log RNA-seq V2 RSEM vs all samples",
        "n_samples": int(expr.shape[1]),
        "n_genes_present": int(expr.shape[0]),
        "genes_present": present,
        "genes_missing_from_api_lookup": missing_lookup,
        "genes_requested_absent_in_matrix": sorted(set([TARGET, *TJ_UNIVERSE]) - set(present)),
        "n_purity": int(purity.notna().sum()),
        "portal": "https://www.cbioportal.org/study/summary?id=luad_oncosg_2020",
    }


def load_hgnc(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df.columns = ["symbol", "previous", "alias", "ensembl"]
    return df


def hgnc_symbol_to_ensembl(hgnc: pd.DataFrame) -> dict[str, str]:
    """Map catalog / current / previous / alias symbols to ENSG (no version)."""
    m: dict[str, str] = {}
    for _, row in hgnc.iterrows():
        ensg = str(row["ensembl"]).strip()
        if not ensg.startswith("ENSG"):
            continue
        names = [row["symbol"]]
        names.extend([x.strip() for x in str(row["previous"]).split(",") if x.strip()])
        names.extend([x.strip() for x in str(row["alias"]).split(",") if x.strip()])
        for name in names:
            m.setdefault(name, ensg)
    return m


def extract_cptac(raw_path: Path, hgnc_path: Path, out: Path) -> dict:
    df = pd.read_csv(raw_path, sep="\t", index_col=0)
    df.index = [str(x).split(".")[0] for x in df.index]
    if df.index.duplicated().any():
        df = df.groupby(level=0).mean()
    hgnc = load_hgnc(hgnc_path)
    sym2ensg = hgnc_symbol_to_ensembl(hgnc)
    mapped = {}
    unmapped = []
    for g in [TARGET, *TJ_UNIVERSE]:
        query = CURRENT_SYMBOL.get(g, g)
        ensg = sym2ensg.get(g) or sym2ensg.get(query)
        if ensg and ensg in df.index:
            mapped[g] = ensg
        else:
            unmapped.append(g)
    sub = df.loc[list(mapped.values())].astype(float)
    sub.index = list(mapped.keys())
    if sub.index.duplicated().any():
        sub = sub.groupby(level=0).mean()
    sub.to_csv(out / "cptac_luad_tj_expr.tsv", sep="\t")
    pd.Series(mapped, name="ensembl").to_csv(out / "cptac_symbol_to_ensembl.tsv", sep="\t")
    return {
        "source_file": raw_path.name,
        "url": CPTAC_URL,
        "unit": "log2(RSEM coding UQ 1500) tumor",
        "gene_id_in_source": "Ensembl gene ID with version; version stripped; mapped via HGNC",
        "hgnc": HGNC_URL,
        "n_samples": int(sub.shape[1]),
        "n_genes_present": int(sub.shape[0]),
        "genes_present": sorted(sub.index.astype(str)),
        "genes_requested_absent": unmapped,
        "n_genes_in_full_matrix": int(df.shape[0]),
        "linkedomics": "https://www.linkedomics.org/data_download/CPTAC-pancan-LUAD/",
        "freeze": "data_freeze_v1.2_reorganized",
    }


def parse_gpl570(annot_path: Path) -> pd.DataFrame:
    rows = []
    with gzip.open(annot_path, "rt", encoding="latin-1", errors="replace") as fh:
        in_table = False
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            if not in_table:
                if line.startswith("ID\t") or line.startswith("ID "):
                    header = line.split("\t")
                    in_table = True
                continue
            if line.startswith("!dataset_table_end") or line.startswith("!platform_table_end"):
                break
            parts = line.split("\t")
            rec = {header[i]: parts[i] if i < len(parts) else "" for i in range(len(header))}
            rows.append(rec)
    plat = pd.DataFrame(rows)
    # GEO .annot uses "Gene symbol"; some dumps use "Gene Symbol"
    sym_col = None
    for c in plat.columns:
        if c.lower().replace(" ", "") in {"genesymbol", "genesymbols"}:
            sym_col = c
            break
    if sym_col is None:
        raise SystemExit(f"GPL570 annot has no Gene symbol column: {list(plat.columns)}")
    plat = plat.rename(columns={"ID": "probe_id", sym_col: "gene_symbol"})
    plat["gene_symbol"] = plat["gene_symbol"].astype(str).str.strip()
    # drop multi-mappers
    plat["n_symbols"] = plat["gene_symbol"].apply(lambda s: 0 if s in {"", "nan", "NA"} else s.count("///") + 1)
    plat.loc[plat["gene_symbol"].isin({"", "nan", "NA", "---"}), "n_symbols"] = 0
    return plat[["probe_id", "gene_symbol", "n_symbols"]]


def parse_gse31210(matrix_path: Path, plat: pd.DataFrame, out: Path) -> dict:
    with gzip.open(matrix_path, "rt", encoding="latin-1", errors="replace") as fh:
        meta_rows = []
        expr = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("!Sample_"):
                parts = line.split("\t")
                key = parts[0][1:]  # drop leading !
                meta_rows.append((key, parts[1:]))
            elif line.startswith("!series_matrix_table_begin"):
                expr = pd.read_csv(fh, sep="\t", index_col=0)
                break
    if expr is None:
        raise SystemExit("GSE31210 series matrix table not found")
    expr = expr.apply(pd.to_numeric, errors="coerce")
    # last row can be !series_matrix_table_end
    expr = expr[~expr.index.astype(str).str.startswith("!")]
    expr.index = expr.index.astype(str)
    expr.columns = [str(c).strip('"') for c in expr.columns]
    # GEO series matrix is MAS5 linear intensity; tests use log2(MAS5+1).
    expr = np.log2(expr.clip(lower=0) + 1)

    # rebuild sample metadata
    acc_row = next(v for k, v in meta_rows if k == "Sample_geo_accession")
    acc = [str(x).strip('"') for x in acc_row]
    meta = pd.DataFrame({"geo_accession": acc})
    # characteristics are repeated keys; collect as ch1_0, ch1_1, ...
    ch_i = 0
    tissue = None
    for key, vals in meta_rows:
        vals = [str(x).strip('"') for x in vals]
        if key == "Sample_characteristics_ch1":
            col = f"ch1_{ch_i}"
            meta[col] = vals
            if vals and vals[0].lower().startswith("tissue:"):
                tissue = [v.split(":", 1)[-1].strip() for v in vals]
            ch_i += 1
        elif key in {"Sample_title", "Sample_source_name_ch1"}:
            meta[key] = vals
    if tissue is None:
        raise SystemExit("GSE31210 tissue characteristic not found")
    meta["tissue"] = tissue
    meta["is_tumor"] = meta["tissue"].str.lower().str.contains("primary lung tumor")
    meta.to_csv(out / "gse31210_sample_meta.tsv", sep="\t", index=False)

    tumor_ids = meta.loc[meta["is_tumor"], "geo_accession"].tolist()
    missing_cols = [s for s in tumor_ids if s not in expr.columns]
    if missing_cols:
        raise SystemExit(f"GSE31210 tumor accessions missing from matrix: {missing_cols[:5]}")
    expr_t = expr.loc[:, tumor_ids]

    plat = plat.set_index("probe_id")
    want = set([TARGET, *TJ_UNIVERSE])
    unique = plat[(plat["n_symbols"] == 1) & (plat["gene_symbol"].isin(want))].copy()
    collapsed = {}
    probe_used = {}
    for gene, block in unique.groupby("gene_symbol"):
        probes = [p for p in block.index if p in expr_t.index]
        if not probes:
            continue
        means = expr_t.loc[probes].mean(axis=1)
        best = str(means.idxmax())
        collapsed[gene] = expr_t.loc[best]
        probe_used[gene] = best
    sub = pd.DataFrame(collapsed).T
    sub.to_csv(out / "gse31210_tj_expr.tsv", sep="\t")
    pd.Series(probe_used, name="probe_id").to_csv(out / "gse31210_probe_used.tsv", sep="\t")
    return {
        "accession": "GSE31210",
        "platform": "GPL570",
        "citation": "Okayama et al. Cancer Res 2012; Japanese stage I-II LUAD",
        "unit": "log2(MAS5+1); series matrix is MAS5 linear intensity",
        "n_samples_in_series": int(meta.shape[0]),
        "n_tumors": int(sub.shape[1]),
        "n_normals_excluded": int((~meta["is_tumor"]).sum()),
        "tumor_rule": "tissue == 'primary lung tumor'",
        "n_genes_present": int(sub.shape[0]),
        "genes_present": sorted(sub.index.astype(str)),
        "genes_requested_absent": sorted(want - set(sub.index)),
        "probe_collapse": "max-mean unique-mapped probe; multi-mapped probes dropped",
        "geo": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31210",
        "why_not_gse72094": (
            "GSE72094 is Moffitt (US) LUAD (Schabath et al.), not an East-Asian series. "
            "GSE31210 is the East-Asian LUAD GEO used here."
        ),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--out-dir", default="results/w200/A9_extra_tj_intersection/data")
    p.add_argument("--force", action="store_true")
    args = p.parse_args()

    out = Path(args.out_dir)
    raw = out / "raw"
    out.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)

    manifest = {
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "cohorts": {},
        "raw_files": {},
    }

    print("=== OncoSG via cBioPortal ===", flush=True)
    manifest["cohorts"]["OncoSG_LUAD"] = fetch_oncosg(out)

    hgnc_dest = raw / "hgnc_symbol_ensembl.tsv"
    if args.force or not hgnc_dest.exists():
        download(HGNC_URL, hgnc_dest)
    else:
        print(f"exists, skipping: {hgnc_dest}", flush=True)
    manifest["raw_files"][hgnc_dest.name] = {
        "url": HGNC_URL,
        "bytes": hgnc_dest.stat().st_size,
        "sha256": sha256_file(hgnc_dest),
    }

    cptac_dest = raw / "LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt"
    if args.force or not cptac_dest.exists():
        download(CPTAC_URL, cptac_dest)
    else:
        print(f"exists, skipping: {cptac_dest}", flush=True)
    manifest["raw_files"][cptac_dest.name] = {
        "url": CPTAC_URL,
        "bytes": cptac_dest.stat().st_size,
        "sha256": sha256_file(cptac_dest),
    }
    print("=== CPTAC LUAD RNA extract ===", flush=True)
    manifest["cohorts"]["CPTAC_LUAD_RNA"] = extract_cptac(cptac_dest, hgnc_dest, out)

    gse_dest = raw / "GSE31210_series_matrix.txt.gz"
    gpl_dest = raw / "GPL570.annot.gz"
    if args.force or not gse_dest.exists():
        download(GSE31210_MATRIX, gse_dest)
    else:
        print(f"exists, skipping: {gse_dest}", flush=True)
    if args.force or not gpl_dest.exists():
        download(GPL570_ANNOT, gpl_dest)
    else:
        print(f"exists, skipping: {gpl_dest}", flush=True)
    manifest["raw_files"][gse_dest.name] = {
        "url": GSE31210_MATRIX,
        "bytes": gse_dest.stat().st_size,
        "sha256": sha256_file(gse_dest),
    }
    manifest["raw_files"][gpl_dest.name] = {
        "url": GPL570_ANNOT,
        "bytes": gpl_dest.stat().st_size,
        "sha256": sha256_file(gpl_dest),
    }
    print("=== GSE31210 extract ===", flush=True)
    plat = parse_gpl570(gpl_dest)
    manifest["cohorts"]["GSE31210_LUAD"] = parse_gse31210(gse_dest, plat, out)

    (out / "SOURCES.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: {ck: cv for ck, cv in v.items() if ck != "genes_present"} if isinstance(v, dict) else v for k, v in manifest["cohorts"].items()}, indent=2))
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
