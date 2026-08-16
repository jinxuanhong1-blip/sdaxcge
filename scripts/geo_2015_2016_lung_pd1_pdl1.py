#!/usr/bin/env python3
"""GEO 2015-2016 lung PD-1/PD-L1 leftover series: TACSTD2 / CLDN4.

Pipeline (reproducible; every accession comes from NCBI, not a hand list):
  1. esearch db=gds for GSE, PDAT 2015/01/01-2016/12/31, lung AND PD-1/PD-L1
     synonyms (PD-1, PD1, PD-L1, PDL1, PDCD1, CD274, B7-H1, "programmed death",
     "programmed cell death").
  2. esummary metadata -> series_metadata.tsv.
  3. Ensembl REST xrefs/symbol for TACSTD2/CLDN4 (human + mouse) ->
     ensembl_id_lookup.tsv. IDs are whatever Ensembl returns that day.
  4. Extract TACSTD2/CLDN4 from processed GEO files only:
       - array: series matrix + GPL annotation (ID / Gene symbol / Alias)
       - RNA-seq: supplementary count/FPKM/RPKM tables, matching gene symbol
         or Ensembl ID (version / _chr suffix stripped)
       - SuperSeries: follow Series_relation SubSeries; do not invent GSEs
  5. Write per-series tables, tacstd2_cldn4_summary.tsv, processing_status.tsv.

Honesty: missing genes, missing tables, and non-lung NCBI hits are recorded
as such. Nothing is imputed.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import re
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from statistics import mean, median

OUTDIR = Path(__file__).resolve().parent.parent / "results" / "w200" / "GEO_2015_2016"
CACHE = Path("/tmp/geo_cache")
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
GEO_FTP = "https://ftp.ncbi.nlm.nih.gov/geo"
ENSEMBL_REST = "https://rest.ensembl.org"

SEARCH_TERM = (
    'lung[All Fields] AND ("PD-1"[All Fields] OR "PD1"[All Fields] OR '
    '"PD-L1"[All Fields] OR "PDL1"[All Fields] OR "PDCD1"[All Fields] OR '
    '"CD274"[All Fields] OR "B7-H1"[All Fields] OR '
    '"programmed death"[All Fields] OR "programmed cell death"[All Fields]) '
    "AND gse[Entry Type] AND 2015/01/01:2016/12/31[PDAT]"
)

# Seed aliases from NCBI Gene otheraliases (human TACSTD2/CLDN4, mouse Tacstd2
# fetched 2016-08-16). Extra aliases from Ensembl xrefs are merged at runtime.
TARGETS: dict[str, set[str]] = {
    "TACSTD2": {
        "TACSTD2", "TROP2", "TROP-2", "M1S1", "GA733-1", "GA7331",
        "EGP-1", "EGP1", "GP50", "LY97",
    },
    "CLDN4": {
        "CLDN4", "CPE-R", "CPER", "CPETR", "CPETR1", "WBSCR8", "HCPE-R",
    },
}
DESC_PHRASES = {
    "TACSTD2": (
        "tumor-associated calcium signal transducer 2",
        "cell surface glycoprotein trop-2",
        "cell surface glycoprotein trop2",
    ),
    "CLDN4": (
        "claudin 4",
        "claudin-4",
        "clostridium perfringens enterotoxin receptor 1",
    ),
}
# Do not use a bare "^ID$" alone on every GPL (some IDs are probe numbers).
# NanoString GPL19965 documents "#ID = Official Symbol".
SYMBOL_COL_RE = re.compile(
    r"(?i)gene.?symbol|^symbol$|^gene$|gene_assignment|^orf$|^GeneSymbol$|"
    r"^Alias$|Official Full Name|Official Symbol"
)
TOKEN_SPLIT_RE = re.compile(r"\s*(?:///|//|[;,|])\s*")
LUNG_TERMS = ("lung", "nsclc", "sclc", "pulmonary")
PD_TERMS = (
    "pd-1", "pd1", "pd-l1", "pdl1", "pdcd1", "cd274", "b7-h1",
    "anti-pd", "checkpoint", "nivolumab", "pembrolizumab",
    "programmed death", "programmed cell death",
)

ENSEMBL_TO_GENE: dict[str, str] = {}


def http_get(url: str, binary: bool = False, retries: int = 4, headers: dict | None = None):
    last = None
    req_headers = {"User-Agent": "geo-2015-2016-leftover/1.0"}
    if headers:
        req_headers.update(headers)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=req_headers)
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
            return data if binary else data.decode("utf-8", errors="replace")
        except urllib.error.HTTPError as e:
            if e.code in (404, 410):
                raise
            last = e
        except Exception as e:
            last = e
        time.sleep(2 ** (attempt + 1))
        print(f"  retry {attempt + 1} for {url}: {last}", file=sys.stderr)
    raise last


def cached_download(url: str, fname: str) -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / fname
    if not path.exists():
        print(f"  downloading {url}")
        path.write_bytes(http_get(url, binary=True))
    return path


def gse_dir(acc: str) -> str:
    return f"{GEO_FTP}/series/{acc[:-3]}nnn/{acc}"


def ftp_listing_hrefs(url: str) -> list[str]:
    try:
        html = http_get(url)
    except urllib.error.HTTPError as e:
        if e.code in (404, 410):
            return []
        raise
    return sorted(set(re.findall(r'href="([^"]+)"', html)))


# ---------------------------------------------------------------- step 1+2
def search_and_summarize():
    url = (EUTILS + "esearch.fcgi?db=gds&retmax=500&retmode=json&term="
           + urllib.parse.quote(SEARCH_TERM))
    res = json.loads(http_get(url))["esearchresult"]
    ids = res["idlist"]
    print(f"esearch: {res['count']} GSE series found")
    if not ids:
        sys.exit("No series returned by NCBI; aborting (nothing to invent).")
    url = EUTILS + "esummary.fcgi?db=gds&retmode=json&id=" + ",".join(ids)
    summ = json.loads(http_get(url))["result"]
    series = []
    for uid in summ["uids"]:
        r = summ[uid]
        series.append({
            "accession": r["accession"],
            "pdat": r["pdat"],
            "gdstype": r["gdstype"],
            "n_samples": r.get("n_samples", ""),
            "taxon": r["taxon"],
            "gpl": ";".join("GPL" + g for g in r.get("gpl", "").split(";") if g),
            "title": r["title"],
            "summary": re.sub(r"\s+", " ", r.get("summary", "")).strip(),
        })
    series.sort(key=lambda s: s["accession"])
    cols = ["accession", "pdat", "gdstype", "n_samples", "taxon", "gpl", "title", "summary"]
    with open(OUTDIR / "series_metadata.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for s in series:
            f.write("\t".join(str(s[c]).replace("\t", " ") for c in cols) + "\n")
    write_lung_pd_flags(series)
    return series


def write_lung_pd_flags(series):
    """Record whether title+summary literally contain lung / PD terms."""
    with open(OUTDIR / "lung_pd_flags.tsv", "w") as f:
        f.write("accession\tcontains_lung_term\tcontains_pd_term\tnote\n")
        for s in series:
            text = f"{s['title']} {s['summary']}".lower()
            lung = any(t in text for t in LUNG_TERMS)
            pd = any(t in text for t in PD_TERMS)
            if lung and pd:
                note = "title/summary mention both lung and PD-1/PD-L1-class terms"
            elif lung:
                note = "title/summary mention lung but not PD-1/PD-L1-class terms; NCBI still returned this GSE"
            elif pd:
                note = "title/summary mention PD-1/PD-L1-class terms but not lung; NCBI still returned this GSE"
            else:
                note = "neither lung nor PD terms in title/summary; NCBI still returned this GSE"
            f.write(f"{s['accession']}\t{str(lung).lower()}\t{str(pd).lower()}\t{note}\n")


# ---------------------------------------------------------- Ensembl lookup
def lookup_ensembl_ids():
    """Official Ensembl gene IDs via REST xrefs/symbol. Written to TSV."""
    queries = [
        ("homo_sapiens", "TACSTD2", "TACSTD2"),
        ("homo_sapiens", "CLDN4", "CLDN4"),
        ("mus_musculus", "Tacstd2", "TACSTD2"),
        ("mus_musculus", "Cldn4", "CLDN4"),
    ]
    rows = []
    for species, symbol, gene in queries:
        url = (f"{ENSEMBL_REST}/xrefs/symbol/{species}/{urllib.parse.quote(symbol)}"
               "?content-type=application/json")
        time.sleep(0.2)
        data = json.loads(http_get(url, headers={"Content-Type": "application/json"}))
        ids = [d["id"] for d in data if d.get("type") == "gene"]
        if not ids:
            print(f"  Ensembl returned no gene id for {species} {symbol}", file=sys.stderr)
        for eid in ids:
            ENSEMBL_TO_GENE[eid] = gene
            rows.append((species, symbol, gene, eid, "xrefs/symbol"))
        # pull extra display aliases from the first gene id
        if ids:
            xurl = (f"{ENSEMBL_REST}/xrefs/id/{ids[0]}"
                    "?content-type=application/json")
            time.sleep(0.2)
            try:
                xrefs = json.loads(http_get(xurl, headers={"Content-Type": "application/json"}))
            except Exception as e:
                print(f"  Ensembl xrefs/id failed for {ids[0]}: {e}", file=sys.stderr)
                continue
            for x in xrefs:
                disp = (x.get("display_id") or "").strip()
                if disp:
                    TARGETS[gene].add(disp.upper())
    with open(OUTDIR / "ensembl_id_lookup.tsv", "w") as f:
        f.write("species\tquery_symbol\ttarget_gene\tensembl_id\tsource\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"Ensembl IDs: {ENSEMBL_TO_GENE}")


def strip_ensembl(token: str) -> str | None:
    m = re.match(r"(ENS[A-Z]*G\d+)", token.strip(), re.I)
    return m.group(1).upper().replace("ENSG", "ENSG").replace("ENSMUSG", "ENSMUSG") if m else None


def match_target(field_value: str) -> str | None:
    if not field_value:
        return None
    for token in TOKEN_SPLIT_RE.split(field_value.strip()):
        t = token.strip().upper()
        for gene, aliases in TARGETS.items():
            if t in aliases:
                return gene
        eid = strip_ensembl(token)
        if eid and eid in ENSEMBL_TO_GENE:
            return ENSEMBL_TO_GENE[eid]
    low = field_value.lower()
    # require phrase match; do not treat "claudin 15" as CLDN4
    for gene, phrases in DESC_PHRASES.items():
        for p in phrases:
            if p in low:
                # reject "claudin 4" hitting "claudin 41" etc.
                if gene == "CLDN4" and re.search(r"claudin[- ]4\d", low):
                    continue
                return gene
    return None


# ------------------------------------------------------------ array handling
def platform_symbol_map(gpl: str) -> dict[str, str]:
    stub = f"{GEO_FTP}/platforms/{gpl[:-3]}nnn/{gpl}"
    text = None
    try:
        path = cached_download(f"{stub}/annot/{gpl}.annot.gz", f"{gpl}.annot.gz")
        text = gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
    except urllib.error.HTTPError:
        url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gpl}"
               "&targ=self&form=text&view=data")
        path = cached_download(url, f"{gpl}.table.txt")
        text = path.read_text(errors="replace")

    lines = text.splitlines()
    try:
        start = next(i for i, l in enumerate(lines)
                     if l.lower().startswith("!platform_table_begin"))
    except StopIteration:
        raise RuntimeError(f"{gpl}: no platform table found")

    # GPL19965: "#ID = Official Symbol"
    id_is_symbol = any(
        re.search(r"(?i)^#ID\s*=\s*Official Symbol", l) for l in lines[:start]
    )
    header = lines[start + 1].split("\t")
    sym_cols = [i for i, h in enumerate(header) if SYMBOL_COL_RE.search(h)]
    if id_is_symbol and 0 not in sym_cols:
        sym_cols.insert(0, 0)
    if not sym_cols:
        raise RuntimeError(f"{gpl}: no gene-symbol-like column in {header}")
    probes = {}
    for line in lines[start + 2:]:
        if line.lower().startswith("!platform_table_end"):
            break
        parts = line.split("\t")
        for ci in sym_cols:
            if ci < len(parts):
                gene = match_target(parts[ci])
                if gene:
                    probes[parts[0]] = gene
                    break
    return probes


def series_matrices(acc: str) -> list[str]:
    hrefs = ftp_listing_hrefs(f"{gse_dir(acc)}/matrix/")
    return [h for h in hrefs if h.endswith("series_matrix.txt.gz")]


def parse_series_matrix(path: Path):
    text = gzip.decompress(path.read_bytes()).decode("utf-8", errors="replace")
    lines = text.splitlines()
    platform = next((l.split("\t")[1].strip('"') for l in lines
                     if l.startswith("!Series_platform_id")), None)
    gsm_titles, gsm_ids = [], []
    rows = {}
    in_table = False
    for l in lines:
        if l.startswith("!series_matrix_table_begin"):
            in_table = True
            continue
        if l.startswith("!series_matrix_table_end"):
            break
        if l.startswith("!Sample_title"):
            gsm_titles = [x.strip('"') for x in l.split("\t")[1:]]
        if not in_table:
            continue
        parts = [x.strip('"') for x in l.split("\t")]
        if parts[0] == "ID_REF":
            gsm_ids = parts[1:]
        else:
            rows[parts[0]] = parts[1:]
    return platform, gsm_ids, gsm_titles, rows


def process_array_series(acc: str):
    matrices = series_matrices(acc)
    if not matrices:
        raise RuntimeError("no series matrix on GEO FTP")
    out_rows = []
    for m in matrices:
        path = cached_download(f"{gse_dir(acc)}/matrix/{m}", m)
        platform, gsm_ids, gsm_titles, rows = parse_series_matrix(path)
        if not rows:
            raise RuntimeError(f"{m}: series matrix contains no data table")
        probes = platform_symbol_map(platform)
        hit_probes = {p: g for p, g in probes.items() if p in rows}
        titles = gsm_titles if len(gsm_titles) == len(gsm_ids) else [""] * len(gsm_ids)
        for probe, gene in sorted(hit_probes.items()):
            for gsm, title, val in zip(gsm_ids, titles, rows[probe]):
                out_rows.append((gene, platform, probe, gsm, title, val))
    return out_rows


# ----------------------------------------------------------- RNA-seq handlers
def clean_bam_sample(name: str) -> str:
    s = name.strip().strip("./")
    s = re.sub(r"/accepted_hits\.bam$", "", s)
    return s.split("/")[-1]


def process_gse84789(acc: str):
    """Gencode v19 featureCounts; Geneid is versioned ENSG."""
    fname = f"{acc}_fcount_primary_output.txt.gz"
    path = cached_download(f"{gse_dir(acc)}/suppl/{fname}", fname)
    out = []
    with gzip.open(path, "rt", errors="replace") as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if header is None:
                header = parts
                continue
            gene = match_target(parts[0])
            if not gene:
                continue
            samples = header[6:] if len(header) > 6 else header[1:]
            values = parts[6:] if len(parts) > 6 else parts[1:]
            for s, v in zip(samples, values):
                out.append((gene, "RNA-seq featureCounts",
                            f"{fname}:{parts[0]}", clean_bam_sample(s), "", v))
    return out


def process_gse78220(acc: str):
    import pandas as pd
    fname = f"{acc}_PatientFPKM.xlsx"
    path = cached_download(f"{gse_dir(acc)}/suppl/{fname}", fname)
    df = pd.read_excel(path)
    gene_col = df.columns[0]
    out = []
    for _, row in df.iterrows():
        gene = match_target(str(row[gene_col]))
        if not gene:
            continue
        for col in df.columns[1:]:
            out.append((gene, "RNA-seq FPKM", f"{fname}:{row[gene_col]}",
                        str(col), "", str(row[col])))
    return out


def process_gse76356(acc: str):
    """Per-sample tables: Gene is ENSMUSG*_chr; value = Unique_RPKM (col 4)."""
    fname = f"{acc}_RAW.tar"
    path = cached_download(f"{gse_dir(acc)}/suppl/{fname}", fname)
    out = []
    with tarfile.open(path) as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            raw = tar.extractfile(member).read()
            if member.name.endswith(".gz"):
                raw = gzip.decompress(raw)
            text = raw.decode("utf-8", errors="replace")
            sample = member.name.split("_")[0]
            lines = text.splitlines()
            if not lines:
                continue
            header = lines[0].split("\t")
            try:
                rpkm_i = header.index("Unique_RPKM")
            except ValueError:
                rpkm_i = 3
            for line in lines[1:]:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= rpkm_i:
                    continue
                # Gene id + Additional_Annotation (tabs inside annotation)
                gene = match_target(parts[0]) or match_target("\t".join(parts[6:]))
                if gene:
                    out.append((gene, "RNA-seq Unique_RPKM",
                                f"{member.name}:{parts[0]}", sample, "", parts[rpkm_i]))
    return out


def series_soft_relations(acc: str) -> list[str]:
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
           "&targ=self&form=text&view=brief")
    text = http_get(url)
    subs = []
    for line in text.splitlines():
        if line.startswith("!Series_relation"):
            m = re.search(r"(GSE\d+)", line)
            if m and ("SuperSeries of" in line or "SubSeries of" in line):
                subs.append(m.group(1))
    return subs


def process_gse81257_deseq():
    """RNA-seq SubSeries of GSE81258: DESeq contrast table, not per-sample counts."""
    acc = "GSE81257"
    fname = f"{acc}_gene_expression_cqnnorm_deseq.csv.gz"
    path = cached_download(f"{gse_dir(acc)}/suppl/{fname}", fname)
    with gzip.open(path, "rt", errors="replace") as f:
        rows = list(csv.reader(f))
    if len(rows) < 4:
        raise RuntimeError(f"{fname}: unexpected DESeq table shape")
    groups, metrics = rows[0], rows[1]
    out = []
    for row in rows[3:]:
        if not row:
            continue
        gene = match_target(row[0])
        if not gene:
            continue
        for i, val in enumerate(row[1:], start=1):
            grp = groups[i] if i < len(groups) else ""
            met = metrics[i] if i < len(metrics) else f"col{i}"
            out.append((gene, "GSE81257 DESeq (cqn-norm)",
                        f"{fname}:{row[0]}:{grp}:{met}",
                        f"{grp}:{met}", "", val))
    return out


def process_gse81258(acc: str):
    """SuperSeries: no series-level expression. Follow Series_relation only."""
    rel = series_soft_relations(acc)
    print(f"  {acc} Series_relation accessions: {rel}")
    if "GSE81257" not in rel:
        raise RuntimeError(
            f"{acc} SuperSeries relations {rel} do not include GSE81257; "
            "not inventing a SubSeries"
        )
    return process_gse81257_deseq()


RNASEQ_HANDLERS = {
    "GSE84789": process_gse84789,
    "GSE78220": process_gse78220,
    "GSE76356": process_gse76356,
    "GSE81258": process_gse81258,
}


# -------------------------------------------------------------------- driver
def summarize(acc: str, rows):
    # GSE81258/GSE81257 is a DESeq contrast table. Do not average log2FC
    # with p-values; numeric summary uses baseMean only.
    if acc == "GSE81258":
        rows = [r for r in rows if str(r[3]).endswith(":baseMean")]
    stats = []
    for gene in TARGETS:
        vals = []
        sources = set()
        for g, _plat, probe, _gsm, _t, v in rows:
            if g != gene:
                continue
            sources.add(probe)
            try:
                if v != "":
                    vals.append(float(v))
            except ValueError:
                pass
        rec = {"accession": acc, "gene": gene, "n_probes_or_rows": len(sources),
               "n_values": len(vals),
               "mean": f"{mean(vals):.4g}" if vals else "NA",
               "median": f"{median(vals):.4g}" if vals else "NA",
               "min": f"{min(vals):.4g}" if vals else "NA",
               "max": f"{max(vals):.4g}" if vals else "NA"}
        stats.append(rec)
    return stats


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)
    series = search_and_summarize()
    lookup_ensembl_ids()

    status, all_stats = [], []
    for s in series:
        acc = s["accession"]
        print(f"== {acc} ({s['gdstype']}) ==")
        rows, note = None, ""
        try:
            if acc in RNASEQ_HANDLERS:
                rows = RNASEQ_HANDLERS[acc](acc)
                if acc == "GSE81258":
                    note = ("SuperSeries has no series-level expression table; "
                            "processed GSE81257 DESeq table named in Series_relation")
                else:
                    note = "processed from supplementary RNA-seq expression table"
            elif "Expression profiling by array" in s["gdstype"]:
                rows = process_array_series(acc)
                note = "processed from series matrix + platform annotation"
            else:
                note = ("not applicable: no processed expression table "
                        f"(series type: {s['gdstype']})")
        except Exception as e:
            note = f"processing failed: {e}"
        if rows is not None:
            out = OUTDIR / f"{acc}_TACSTD2_CLDN4_expression.tsv"
            with open(out, "w") as f:
                f.write("gene\tplatform_or_source\tprobe_or_row\tsample\tsample_title\tvalue\n")
                for r in rows:
                    f.write("\t".join(r) + "\n")
            all_stats.extend(summarize(acc, rows))
            n_genes = len({r[0] for r in rows})
            if rows:
                note += f"; {len(rows)} values, {n_genes}/2 target genes found"
            else:
                note += ("; 0 values, 0/2 target genes found in the processed "
                         "table (gene absent from platform/file, not imputed)")
        status.append({"accession": acc, "taxon": s["taxon"],
                       "gdstype": s["gdstype"], "status": note})
        print(f"  {note}")

    with open(OUTDIR / "processing_status.tsv", "w") as f:
        f.write("accession\ttaxon\tgdstype\tstatus\n")
        for st in status:
            f.write("\t".join(st[c] for c in ("accession", "taxon", "gdstype", "status")) + "\n")

    cols = ["accession", "gene", "n_probes_or_rows", "n_values", "mean", "median", "min", "max"]
    with open(OUTDIR / "tacstd2_cldn4_summary.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for st in all_stats:
            f.write("\t".join(str(st[c]) for c in cols) + "\n")
    print("done")


if __name__ == "__main__":
    main()
