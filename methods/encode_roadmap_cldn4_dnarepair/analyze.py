#!/usr/bin/env python3
"""Exploratory ENCODE / Roadmap test: is CLDN4 linked to DNA-repair genes
in epithelial Hi-C, ChIA-PET, or ChIP-derived chromatin state?

Definitions are fixed below before the intersection is interpreted.
Nothing here is a regulatory or clinical claim.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.join(ROOT, "results", "encode_roadmap_cldn4_dnarepair")
CACHE = "/tmp/encode_cldn4_cache"
UA = {"User-Agent": "sdaxcge-cldn4-dnarepair/0.1", "Accept": "application/json"}

# Pre-specified windows. Gene coordinates from Ensembl are 1-based inclusive
# and are converted to 0-based half-open BED intervals before overlap tests.
LOOP_FLANK = 5000
TSS_PROMOTER = 2000
DOMAIN_MAX = 2_000_000
NEAR_BP = 1_000_000
LOCAL_PEAK_MAX = 50_000

KEGG_MAPS = {
    "hsa03410": "base excision repair",
    "hsa03420": "nucleotide excision repair",
    "hsa03430": "mismatch repair",
    "hsa03440": "homologous recombination",
    "hsa03450": "non-homologous end-joining",
    "hsa03460": "Fanconi anemia",
}

# Official Roadmap STD_NAME rows (EID_metadata.tab). Primary = epithelial
# cell or carcinoma line. Mucosa = epithelial-rich tissue, not a pure isolate.
# E096 bulk lung is a comparator, not an epithelial call.
ROADMAP_PRIMARY = [
    "E027", "E028", "E057", "E058", "E114", "E117", "E118", "E119", "E127",
]
ROADMAP_MUCOSA = [
    "E075", "E079", "E094", "E101", "E102", "E106", "E109", "E110",
]
ROADMAP_LUNG_BULK = ["E096"]

EXCLUDE_BIOSAMPLE = (
    "endothelial",
    "mesothelial",
    "sk-mel",
    "rpmi7951",
    "sk-n-mc",
    "hek293",
)

ASSAYS_3D = ["intact Hi-C", "in situ Hi-C", "dilution Hi-C", "ChIA-PET"]
CHIP_BIOSAMPLES = [
    "A549",
    "NCI-H460",
    "HCT116",
    "HepG2",
    "MCF-7",
    "MCF 10A",
    "HeLa-S3",
    "keratinocyte",
    "mammary epithelial cell",
    "Panc1",
]
CHIP_TARGETS = [
    ("Histone ChIP-seq", "H3K27ac"),
    ("Histone ChIP-seq", "H3K4me3"),
    ("Histone ChIP-seq", "H3K4me1"),
    ("Histone ChIP-seq", "H3K27me3"),
    ("TF ChIP-seq", "CTCF"),
    ("TF ChIP-seq", "POLR2A"),
    ("TF ChIP-seq", "RAD21"),
]
HISTONE_PEAK_PREFERENCE = [
    "replicated peaks",
    "pseudoreplicated peaks",
    "optimal IDR thresholded peaks",
    "IDR thresholded peaks",
]
TF_PEAK_PREFERENCE = [
    "IDR thresholded peaks",
    "optimal IDR thresholded peaks",
    "conservative IDR thresholded peaks",
    "replicated peaks",
]

os.makedirs(OUT, exist_ok=True)
os.makedirs(CACHE, exist_ok=True)


def http_json(url, data=None, headers=None, timeout=180, attempts=4):
    hdr = dict(UA)
    if headers:
        hdr.update(headers)
    payload = None if data is None else data
    last = None
    for i in range(attempts):
        req = urllib.request.Request(url, data=payload, headers=hdr)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return json.load(resp)
        except urllib.error.HTTPError as exc:
            last = exc
            body = exc.read().decode("utf-8", "replace")[:300]
            if exc.code == 404:
                raise
            if exc.code in (403, 429, 500, 502, 503, 504):
                time.sleep(2 ** i)
                continue
            raise RuntimeError(f"HTTP {exc.code} {url[:140]} {body}") from exc
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"failed {url[:160]}: {last}")


def http_bytes(url, timeout=300, attempts=4):
    hdr = {"User-Agent": UA["User-Agent"]}
    last = None
    for i in range(attempts):
        req = urllib.request.Request(url, headers=hdr)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except Exception as exc:  # noqa: BLE001
            last = exc
            time.sleep(2 ** i)
    raise RuntimeError(f"download failed {url[:160]}: {last}")


def encode_search(params):
    q = dict(params)
    q["format"] = "json"
    q.setdefault("limit", "all")
    url = "https://www.encodeproject.org/search/?" + urllib.parse.urlencode(q)
    try:
        return http_json(url)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return {"total": 0, "@graph": []}
        raise


def open_maybe_gzip(raw, name):
    if raw[:2] == b"\x1f\x8b" or name.endswith(".gz"):
        return gzip.GzipFile(fileobj=io.BytesIO(raw))
    return io.BytesIO(raw)


def overlaps(a1, a2, b1, b2):
    return a1 < b2 and b1 < a2


def gap(a1, a2, b1, b2):
    if overlaps(a1, a2, b1, b2):
        return 0
    if a2 <= b1:
        return b1 - a2
    return a1 - b2


def to_bed(start_1based, end_inclusive):
    return int(start_1based) - 1, int(end_inclusive)


def stratum_of(description):
    """Split parental maps from tagged-untreated degron lines and from acute treatment.

    Dexamethasone 0-hour labels still contain the word "treated". Those are
    vehicle/time-zero controls and stay in the baseline stratum.
    """
    dl = (description or "").lower()
    if re.search(r"deletion|knock-?down", dl):
        return "treated"
    if "dexamethasone" in dl:
        if re.search(r"\b0\s*hours?\b", dl):
            return "baseline"
        return "treated"
    tagged = bool(re.search(r"degron|auxin|ostir|maid|iaa", dl))
    if tagged:
        if "untreated" in dl or re.search(r"\bcontrol\b", dl):
            return "tagged_untreated"
        return "treated"
    if re.search(r"\btreated\b", dl) and "untreated" not in dl:
        return "treated"
    return "baseline"


def epithelial_call(term, organ_slims):
    term_l = (term or "").lower()
    if any(x in term_l for x in EXCLUDE_BIOSAMPLE):
        return None
    slims = organ_slims or []
    if isinstance(slims, str):
        slims = [slims]
    if "epithelium" in slims:
        return "encode_epithelium_slim"
    if term == "A549":
        return "lung_carcinoma_line"
    return None


def fetch_kegg_repair_genes():
    print("KEGG DNA repair gene set", flush=True)
    symbol_pathways = defaultdict(set)
    entrez_for_symbol = {}
    # hsa id -> pathways, then symbols
    hsa_path = defaultdict(set)
    for map_id, name in KEGG_MAPS.items():
        raw = http_bytes(f"https://rest.kegg.jp/link/hsa/{map_id}").decode()
        ids = []
        for line in raw.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                ids.append(parts[1])
                hsa_path[parts[1]].add(name)
        print(f"  {map_id} {name}: {len(ids)}", flush=True)
    # KEGG hsa:N is the NCBI Gene id for human.
    id_list = sorted(hsa_path)
    id_to_symbol = {}
    for i in range(0, len(id_list), 40):
        chunk = [x.split(":")[1] for x in id_list[i : i + 40]]
        body = json.dumps({"ids": chunk, "fields": "symbol,name"}).encode()
        hits = http_json(
            "https://mygene.info/v3/gene",
            data=body,
            headers={"Content-Type": "application/json"},
        )
        for hit in hits:
            if hit.get("notfound") or "symbol" not in hit:
                continue
            id_to_symbol[hit["query"]] = hit["symbol"]
    missing = []
    for hsa, paths in hsa_path.items():
        entrez = hsa.split(":")[1]
        sym = id_to_symbol.get(entrez)
        if not sym:
            missing.append(entrez)
            continue
        symbol_pathways[sym].update(paths)
        entrez_for_symbol[sym] = entrez
    symbols = sorted(symbol_pathways)
    print(f"  symbols {len(symbols)} unmapped entrez {len(missing)}", flush=True)
    coords = {}
    for i in range(0, len(symbols), 40):
        chunk = symbols[i : i + 40]
        body = json.dumps({"symbols": chunk}).encode()
        data = http_json(
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
            data=body,
            headers={"Content-Type": "application/json", "Accept": "application/json"},
        )
        for sym, rec in data.items():
            if not isinstance(rec, dict):
                continue
            chrom = str(rec.get("seq_region_name"))
            if not chrom.isdigit() and chrom not in ("X", "Y", "MT"):
                continue
            start, end = to_bed(rec["start"], rec["end"])
            strand = int(rec.get("strand") or 0)
            tss = start if strand == 1 else end - 1
            coords[sym] = {
                "symbol": sym,
                "ensembl": rec.get("id"),
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "tss": tss,
                "entrez": entrez_for_symbol.get(sym),
                "pathways": sorted(symbol_pathways[sym]),
            }
    print(f"  GRCh38 coords {len(coords)}", flush=True)
    return coords, sorted(missing)


def fetch_cldn4():
    body = json.dumps({"symbols": ["CLDN4"]}).encode()
    data = http_json(
        "https://rest.ensembl.org/lookup/symbol/homo_sapiens",
        data=body,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    rec = data["CLDN4"]
    start, end = to_bed(rec["start"], rec["end"])
    strand = int(rec["strand"])
    tss = start if strand == 1 else end - 1
    return {
        "symbol": "CLDN4",
        "ensembl": rec["id"],
        "chrom": str(rec["seq_region_name"]),
        "start": start,
        "end": end,
        "strand": strand,
        "tss": tss,
        "description": rec.get("description"),
    }


def neighborhood_genes(chrom, center_start, center_end, pad=1_500_000):
    """Protein-coding genes on a tiled Ensembl overlap around CLDN4."""
    left = max(1, center_start - pad)
    right = center_end + pad
    genes = []
    seen = set()
    pos = left
    while pos < right:
        # Ensembl overlap rejects windows much above 5 Mb.
        end = min(right, pos + 4_000_000)
        region = f"{chrom}:{pos}-{end}"
        url = (
            "https://rest.ensembl.org/overlap/region/human/"
            + region
            + "?feature=gene&content-type=application/json"
        )
        rows = http_json(url, headers={"Accept": "application/json", "Content-Type": "application/json"})
        for g in rows:
            if g.get("id") in seen:
                continue
            seen.add(g.get("id"))
            if g.get("biotype") != "protein_coding":
                continue
            s, e = to_bed(g["start"], g["end"])
            strand = int(g.get("strand") or 0)
            genes.append(
                {
                    "symbol": g.get("external_name"),
                    "ensembl": g.get("id"),
                    "chrom": chrom,
                    "start": s,
                    "end": e,
                    "strand": strand,
                    "tss": s if strand == 1 else e - 1,
                    "description": (g.get("description") or "").split("[")[0].strip(),
                }
            )
        pos = end
    return genes


def distance_table(cldn, repair):
    rows = []
    for sym, g in repair.items():
        chrom_same = g["chrom"] == cldn["chrom"]
        if chrom_same:
            g_bp = gap(cldn["start"], cldn["end"], g["start"], g["end"])
            tss_bp = abs(g["tss"] - cldn["tss"])
        else:
            g_bp = None
            tss_bp = None
        rows.append(
            {
                "symbol": sym,
                "entrez": g["entrez"],
                "chrom": g["chrom"],
                "start0": g["start"],
                "end0": g["end"],
                "tss0": g["tss"],
                "strand": g["strand"],
                "same_chrom_as_CLDN4": chrom_same,
                "gap_bp": g_bp if g_bp is not None else "",
                "tss_distance_bp": tss_bp if tss_bp is not None else "",
                "within_1Mb_gap": bool(chrom_same and g_bp <= NEAR_BP),
                "within_2Mb_gap": bool(chrom_same and g_bp <= 2 * NEAR_BP),
                "pathways": ";".join(g["pathways"]),
            }
        )
    rows.sort(
        key=lambda r: (
            0 if r["same_chrom_as_CLDN4"] else 1,
            r["gap_bp"] if r["gap_bp"] != "" else 10**18,
            r["symbol"],
        )
    )
    return rows


def write_tsv(path, rows, fieldnames):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)


def load_experiments():
    print("ENCODE experiment inventory", flush=True)
    found = {}
    for assay in ASSAYS_3D:
        data = encode_search(
            {
                "type": "Experiment",
                "assay_title": assay,
                "status": "released",
                "biosample_ontology.organ_slims": "epithelium",
            }
        )
        for e in data.get("@graph") or []:
            found[e["accession"]] = assay
        print(f"  epithelium slim {assay}: {data.get('total')}", flush=True)
        data = encode_search(
            {
                "type": "Experiment",
                "assay_title": assay,
                "status": "released",
                "biosample_ontology.term_name": "A549",
            }
        )
        for e in data.get("@graph") or []:
            found[e["accession"]] = assay
        print(f"  A549 {assay}: {data.get('total')}", flush=True)
    print(f"  unique experiments {len(found)}", flush=True)
    return found


def fetch_experiment(accession):
    cache = os.path.join(CACHE, f"{accession}.json")
    if os.path.exists(cache) and os.path.getsize(cache) > 100:
        with open(cache) as fh:
            return json.load(fh)
    data = http_json(f"https://www.encodeproject.org/experiments/{accession}/?format=json")
    with open(cache, "w") as fh:
        json.dump(data, fh)
    return data


def classify_experiment(exp, assay_hint):
    bios = exp.get("biosample_ontology") or {}
    term = bios.get("term_name")
    slims = bios.get("organ_slims") or []
    call = epithelial_call(term, slims)
    desc = exp.get("description") or ""
    target = exp.get("target") or {}
    target_label = target.get("label") if isinstance(target, dict) else None
    files = []
    for f in exp.get("files") or []:
        if not isinstance(f, dict):
            continue
        if f.get("status") != "released":
            continue
        if f.get("assembly") != "GRCh38":
            continue
        ot = f.get("output_type")
        ff = f.get("file_format")
        if ot == "loops" and ff in ("bedpe", "tsv", "bed"):
            kind = "loops"
        elif ot == "contact domains" and ff in ("bed", "bedpe"):
            kind = "domains"
        else:
            continue
        files.append(
            {
                "accession": f.get("accession"),
                "href": "https://www.encodeproject.org" + f.get("href", ""),
                "file_format": ff,
                "file_size": f.get("file_size"),
                "output_type": ot,
                "kind": kind,
                "biological_replicates": f.get("biological_replicates") or [],
            }
        )
    return {
        "experiment": exp.get("accession"),
        "assay": exp.get("assay_title") or assay_hint,
        "biosample": term,
        "classification": bios.get("classification"),
        "organ_slims": slims,
        "epithelial_call": call,
        "description": desc,
        "stratum": stratum_of(desc),
        "target": target_label,
        "files": files,
    }


def select_pooled(files):
    """Keep the replicate-union file(s) when both single-rep and pooled calls exist."""
    if not files:
        return []
    width = [len(f["biological_replicates"] or []) for f in files]
    best = max(width)
    # A file with an empty replicate list is treated as pooled-unknown, kept
    # only when nothing has a replicate annotation.
    if best == 0:
        return files
    return [f for f, w in zip(files, width) if w == best]


def parse_intervals_from_line(line, file_format, kind):
    if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 3:
        return None
    # header
    if parts[1] in ("start", "x1", "start1") or not parts[1][:1].isdigit() and parts[1][:1] != "-":
        try:
            int(parts[1])
        except ValueError:
            return None
    try:
        if file_format == "bed" or (kind == "domains" and file_format == "bed"):
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            return [(chrom, start, end)]
        if len(parts) < 6:
            return None
        c1, s1, e1 = parts[0], int(parts[1]), int(parts[2])
        c2, s2, e2 = parts[3], int(parts[4]), int(parts[5])
        if kind == "domains":
            if c1 != c2:
                return [(c1, s1, e1), (c2, s2, e2)]
            return [(c1, min(s1, s2), max(e1, e2))]
        return [(c1, s1, e1), (c2, s2, e2)]
    except ValueError:
        return None


def stream_hits(file_rec, cldn, repair_by_chrom):
    """Return loop links and domain co-occupancy for one file."""
    raw = http_bytes(file_rec["href"])
    name = file_rec["href"].rsplit("/", 1)[-1]
    fh = open_maybe_gzip(raw, name)
    c_chrom = "chr" + cldn["chrom"]
    c1, c2 = cldn["start"] - LOOP_FLANK, cldn["end"] + LOOP_FLANK
    text = io.TextIOWrapper(fh, encoding="utf-8", errors="replace")
    loop_hits = []
    domain_hits = []
    n_loops_cldn = 0
    n_domains_cldn = 0
    n_records = 0
    for line in text:
        recs = parse_intervals_from_line(line, file_rec["file_format"], file_rec["kind"])
        if not recs:
            continue
        n_records += 1
        if file_rec["kind"] == "domains":
            for chrom, start, end in recs:
                if chrom == c_chrom and overlaps(start, end, cldn["start"], cldn["end"]):
                    n_domains_cldn += 1
                    length = end - start
                    for g in repair_by_chrom.get(cldn["chrom"], []):
                        if overlaps(start, end, g["start"], g["end"]):
                            domain_hits.append(
                                {
                                    "repair_gene": g["symbol"],
                                    "domain_chrom": chrom,
                                    "domain_start0": start,
                                    "domain_end0": end,
                                    "domain_bp": length,
                                    "within_2Mb_cap": length <= DOMAIN_MAX,
                                    "pathways": ";".join(g["pathways"]),
                                }
                            )
            continue
        # loops: two anchors
        if len(recs) != 2:
            continue
        (a_c, a_s, a_e), (b_c, b_s, b_e) = recs
        a_hit = a_c == c_chrom and overlaps(a_s, a_e, c1, c2)
        b_hit = b_c == c_chrom and overlaps(b_s, b_e, c1, c2)
        if not a_hit and not b_hit:
            continue
        n_loops_cldn += 1
        pairs = []
        if a_hit:
            pairs.append((b_c, b_s, b_e, a_s, a_e))
        if b_hit:
            pairs.append((a_c, a_s, a_e, b_s, b_e))
        for p_c, p_s, p_e, s_s, s_e in pairs:
            chrom_name = p_c[3:] if p_c.startswith("chr") else p_c
            for g in repair_by_chrom.get(chrom_name, []):
                g1, g2 = g["start"] - LOOP_FLANK, g["end"] + LOOP_FLANK
                tss1, tss2 = g["tss"] - TSS_PROMOTER, g["tss"] + TSS_PROMOTER
                if overlaps(p_s, p_e, g1, g2):
                    loop_hits.append(
                        {
                            "repair_gene": g["symbol"],
                            "partner_chrom": p_c,
                            "partner_start0": p_s,
                            "partner_end0": p_e,
                            "cldn_anchor_start0": s_s,
                            "cldn_anchor_end0": s_e,
                            "gene_body_flank": True,
                            "promoter_both": overlaps(s_s, s_e, cldn["tss"] - TSS_PROMOTER, cldn["tss"] + TSS_PROMOTER)
                            and overlaps(p_s, p_e, tss1, tss2),
                            "pathways": ";".join(g["pathways"]),
                        }
                    )
    return {
        "n_records": n_records,
        "n_cldn_loops": n_loops_cldn,
        "n_cldn_domains": n_domains_cldn,
        "loop_hits": loop_hits,
        "domain_hits": domain_hits,
    }


def analyze_3d(cldn, repair):
    accession_assay = load_experiments()
    print(f"Fetching {len(accession_assay)} experiment objects", flush=True)
    records = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        futs = {pool.submit(fetch_experiment, acc): (acc, assay) for acc, assay in accession_assay.items()}
        done = 0
        for fut in as_completed(futs):
            acc, assay = futs[fut]
            done += 1
            try:
                exp = fut.result()
            except Exception as exc:  # noqa: BLE001
                print(f"  experiment failed {acc}: {exc}", flush=True)
                continue
            records.append(classify_experiment(exp, assay))
            if done % 25 == 0:
                print(f"  experiments {done}/{len(futs)}", flush=True)
    repair_by_chrom = defaultdict(list)
    for g in repair.values():
        repair_by_chrom[g["chrom"]].append(g)

    manifest = []
    loop_rows = []
    domain_rows = []
    # Download selected files. Primary cohort is epithelial baseline.
    jobs = []
    for rec in records:
        if not rec["epithelial_call"]:
            for f in rec["files"]:
                manifest.append({**_file_manifest(rec, f), "used": "no", "reason": "not_epithelial_call"})
            continue
        if rec["stratum"] == "treated":
            for f in rec["files"]:
                manifest.append({**_file_manifest(rec, f), "used": "no", "reason": "treated"})
            continue
        by_kind = defaultdict(list)
        for f in rec["files"]:
            by_kind[f["kind"]].append(f)
        for kind, files in by_kind.items():
            chosen = select_pooled(files)
            chosen_ids = {f["accession"] for f in chosen}
            for f in files:
                if f["accession"] not in chosen_ids:
                    manifest.append({**_file_manifest(rec, f), "used": "no", "reason": "single_replicate_when_pooled_exists"})
            for f in chosen:
                jobs.append((rec, f))

    print(f"Downloading {len(jobs)} loop/domain files", flush=True)

    def _run(job):
        rec, f = job
        try:
            hit = stream_hits(f, cldn, repair_by_chrom)
            err = ""
        except Exception as exc:  # noqa: BLE001
            hit = None
            err = str(exc)[:300]
        return rec, f, hit, err

    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = [pool.submit(_run, job) for job in jobs]
        done = 0
        for fut in as_completed(futs):
            rec, f, hit, err = fut.result()
            done += 1
            if done % 10 == 0 or done == len(futs):
                print(f"  files {done}/{len(futs)}", flush=True)
            base = _file_manifest(rec, f)
            if err:
                manifest.append({**base, "used": "error", "reason": err, "n_records": "", "n_cldn_loops": "", "n_cldn_domains": ""})
                continue
            manifest.append(
                {
                    **base,
                    "used": "yes",
                    "reason": rec["stratum"],
                    "n_records": hit["n_records"],
                    "n_cldn_loops": hit["n_cldn_loops"],
                    "n_cldn_domains": hit["n_cldn_domains"],
                }
            )
            for h in hit["loop_hits"]:
                loop_rows.append({**base, **h})
            for h in hit["domain_hits"]:
                domain_rows.append({**base, **h})
    return records, manifest, loop_rows, domain_rows


def _file_manifest(rec, f):
    return {
        "experiment": rec["experiment"],
        "assay": rec["assay"],
        "biosample": rec["biosample"],
        "epithelial_call": rec["epithelial_call"],
        "stratum": rec["stratum"],
        "target": rec["target"] or "",
        "description": (rec["description"] or "").replace("\t", " ")[:240],
        "file": f["accession"],
        "kind": f["kind"],
        "file_format": f["file_format"],
        "file_size": f["file_size"],
        "biological_replicates": ",".join(str(x) for x in (f["biological_replicates"] or [])),
    }


def experiment_level(manifest, loop_rows, domain_rows):
    """One row per experiment x kind that was actually scanned."""
    used = [m for m in manifest if m.get("used") == "yes"]
    loops_by = defaultdict(list)
    domains_by = defaultdict(list)
    for row in loop_rows:
        loops_by[row["experiment"]].append(row)
    for row in domain_rows:
        domains_by[row["experiment"]].append(row)
    out = []
    seen = set()
    for m in used:
        key = (m["experiment"], m["kind"])
        if key in seen:
            continue
        seen.add(key)
        if m["kind"] == "loops":
            hits = loops_by.get(m["experiment"], [])
            genes = sorted({h["repair_gene"] for h in hits})
            promo = sorted({h["repair_gene"] for h in hits if str(h.get("promoter_both")) in ("True", "true", True)})
            out.append(
                {
                    **{k: m[k] for k in ("experiment", "assay", "biosample", "epithelial_call", "stratum", "target", "description")},
                    "kind": "loops",
                    "files_scanned": sum(1 for x in used if x["experiment"] == m["experiment"] and x["kind"] == "loops"),
                    "cldn_touching_records": sum(int(x["n_cldn_loops"] or 0) for x in used if x["experiment"] == m["experiment"] and x["kind"] == "loops"),
                    "repair_gene_hit": bool(genes),
                    "repair_genes": ";".join(genes),
                    "promoter_repair_genes": ";".join(promo),
                    "domain_le_2Mb": "",
                    "repair_genes_le_2Mb": "",
                }
            )
        else:
            hits = domains_by.get(m["experiment"], [])
            genes = sorted({h["repair_gene"] for h in hits})
            short = sorted({h["repair_gene"] for h in hits if str(h.get("within_2Mb_cap")) in ("True", "true", True)})
            out.append(
                {
                    **{k: m[k] for k in ("experiment", "assay", "biosample", "epithelial_call", "stratum", "target", "description")},
                    "kind": "domains",
                    "files_scanned": sum(1 for x in used if x["experiment"] == m["experiment"] and x["kind"] == "domains"),
                    "cldn_touching_records": sum(int(x["n_cldn_domains"] or 0) for x in used if x["experiment"] == m["experiment"] and x["kind"] == "domains"),
                    "repair_gene_hit": bool(genes),
                    "repair_genes": ";".join(genes),
                    "promoter_repair_genes": "",
                    "domain_le_2Mb": bool(short),
                    "repair_genes_le_2Mb": ";".join(short),
                }
            )
    return out


def roadmap_states(cldn, repair_near):
    print("Roadmap chromHMM", flush=True)
    meta_raw = http_bytes("https://egg2.wustl.edu/roadmap/data/byFileType/metadata/EID_metadata.tab").decode()
    meta = {}
    header = None
    for line in meta_raw.splitlines():
        parts = line.split("\t")
        if header is None:
            header = parts
            continue
        rec = dict(zip(header, parts))
        meta[rec["EID"]] = rec
    points = [("CLDN4_TSS", cldn["tss"])]
    for g in repair_near:
        points.append((g["symbol"] + "_TSS", g["tss"]))
    # also gene-body coverage is summarized by the TSS state only; record the
    # segment that contains each point.
    wanted = ROADMAP_PRIMARY + ROADMAP_MUCOSA + ROADMAP_LUNG_BULK
    rows = []
    base = (
        "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/"
        "ChmmModels/coreMarks/jointModel/final/"
    )
    for eid in wanted:
        info = meta.get(eid, {})
        if eid in ROADMAP_PRIMARY:
            panel = "epithelial_cell_or_line"
        elif eid in ROADMAP_MUCOSA:
            panel = "mucosal_tissue"
        else:
            panel = "bulk_lung_comparator"
        url = base + f"{eid}_15_coreMarks_hg38lift_mnemonics.bed.gz"
        try:
            raw = http_bytes(url)
        except Exception as exc:  # noqa: BLE001
            print(f"  {eid} download failed {exc}", flush=True)
            for label, tss in points:
                rows.append(
                    {
                        "eid": eid,
                        "panel": panel,
                        "std_name": info.get("STD_NAME", ""),
                        "group": info.get("GROUP", ""),
                        "mnemonic": info.get("MNEMONIC", ""),
                        "anatomy": info.get("ANATOMY", ""),
                        "query": label,
                        "state": "",
                        "segment_start0": "",
                        "segment_end0": "",
                        "error": str(exc)[:200],
                    }
                )
            continue
        fh = gzip.GzipFile(fileobj=io.BytesIO(raw))
        # Collect chr7 segments once, then assign points.
        segs = []
        for line in io.TextIOWrapper(fh, encoding="utf-8", errors="replace"):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] not in ("chr7", "7"):
                continue
            try:
                s, e = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            segs.append((s, e, parts[3]))
        print(f"  {eid} {info.get('STD_NAME','')} chr7 segments {len(segs)}", flush=True)
        for label, tss in points:
            state, ss, ee = "", "", ""
            for s, e, name in segs:
                if s <= tss < e:
                    state, ss, ee = name, s, e
                    break
            rows.append(
                {
                    "eid": eid,
                    "panel": panel,
                    "std_name": info.get("STD_NAME", ""),
                    "group": info.get("GROUP", ""),
                    "mnemonic": info.get("MNEMONIC", ""),
                    "anatomy": info.get("ANATOMY", ""),
                    "query": label,
                    "state": state,
                    "segment_start0": ss,
                    "segment_end0": ee,
                    "error": "",
                }
            )
    return rows


def _empty_chip_row(q, picked, error, windows):
    row = {
        "biosample": q["biosample"],
        "assay": q["assay"],
        "target": q["target"],
        "n_released": q["n_released"],
        "n_baseline": q["n_baseline"],
        "picked": picked,
        "file": "",
        "output_type": "",
        "file_size": "",
        "error": error,
        "peaks_overlapping_both_promoters": "",
        "peaks_overlapping_both_promoters_le_50kb": "",
    }
    for name in windows:
        row[name] = ""
    return row


def chip_peak_survey(cldn, rfc2):
    """One baseline peak file per epithelial biosample x target, TSS ± 2 kb."""
    print("ENCODE ChIP peak survey", flush=True)
    windows = {
        "CLDN4_promoter": (cldn["tss"] - TSS_PROMOTER, cldn["tss"] + TSS_PROMOTER),
        "RFC2_promoter": (rfc2["tss"] - TSS_PROMOTER, rfc2["tss"] + TSS_PROMOTER),
        "CLDN4_body": (cldn["start"], cldn["end"]),
        "RFC2_body": (rfc2["start"], rfc2["end"]),
    }
    queries = []
    for biosample in CHIP_BIOSAMPLES:
        for assay, target in CHIP_TARGETS:
            data = encode_search(
                {
                    "type": "Experiment",
                    "assay_title": assay,
                    "status": "released",
                    "biosample_ontology.term_name": biosample,
                    "target.label": target,
                }
            )
            graph = data.get("@graph") or []
            baseline = []
            for e in graph:
                desc = e.get("description") or ""
                if stratum_of(desc) != "baseline":
                    continue
                baseline.append(e["accession"])
            queries.append(
                {
                    "biosample": biosample,
                    "assay": assay,
                    "target": target,
                    "n_released": data.get("total") or 0,
                    "n_baseline": len(baseline),
                    "candidates": baseline[:4],
                }
            )
            print(f"  inventory {biosample} {target}: released {data.get('total') or 0} baseline {len(baseline)}", flush=True)
    rows = []
    for q in queries:
        if not q["candidates"]:
            rows.append(_empty_chip_row(q, "", "no_baseline_experiment", windows))
            continue
        chosen = None
        chosen_exp = ""
        last_err = "no_GRCh38_bed_peaks"
        pref = HISTONE_PEAK_PREFERENCE if q["assay"].startswith("Histone") else TF_PEAK_PREFERENCE
        for acc in q["candidates"]:
            try:
                exp = fetch_experiment(acc)
            except Exception as exc:  # noqa: BLE001
                last_err = str(exc)[:200]
                continue
            candidates = []
            for f in exp.get("files") or []:
                if not isinstance(f, dict):
                    continue
                if f.get("status") != "released" or f.get("assembly") != "GRCh38":
                    continue
                if f.get("file_format") != "bed":
                    continue
                if f.get("output_type") not in pref:
                    continue
                if (f.get("file_size") or 0) > 30_000_000:
                    continue
                candidates.append(f)
            if not candidates:
                chosen_exp = acc
                continue
            candidates.sort(key=lambda f: (pref.index(f["output_type"]), -len(f.get("biological_replicates") or [])))
            chosen = candidates[0]
            chosen_exp = acc
            break
        if chosen is None:
            rows.append(_empty_chip_row(q, chosen_exp, last_err, windows))
            continue
        href = "https://www.encodeproject.org" + chosen["href"]
        try:
            raw = http_bytes(href)
            fh = open_maybe_gzip(raw, chosen["href"])
            counts = {k: 0 for k in windows}
            both_promoter = 0
            both_local = 0
            for line in io.TextIOWrapper(fh, encoding="utf-8", errors="replace"):
                parts = line.split("\t")
                if len(parts) < 3 or parts[0] not in ("chr7", "7"):
                    continue
                try:
                    s, e = int(parts[1]), int(parts[2])
                except ValueError:
                    continue
                hit_names = []
                for name, (ws, we) in windows.items():
                    if overlaps(s, e, ws, we):
                        counts[name] += 1
                        hit_names.append(name)
                if "CLDN4_promoter" in hit_names and "RFC2_promoter" in hit_names:
                    both_promoter += 1
                    if (e - s) <= LOCAL_PEAK_MAX:
                        both_local += 1
            row = _empty_chip_row(q, chosen_exp, "", windows)
            row.update(counts)
            row["file"] = chosen["accession"]
            row["output_type"] = chosen["output_type"]
            row["file_size"] = chosen.get("file_size") or ""
            row["peaks_overlapping_both_promoters"] = both_promoter
            row["peaks_overlapping_both_promoters_le_50kb"] = both_local
            rows.append(row)
            print(
                f"  {q['biosample']} {q['target']} {chosen['accession']} "
                f"CLDN4p {counts['CLDN4_promoter']} RFC2p {counts['RFC2_promoter']}",
                flush=True,
            )
        except Exception as exc:  # noqa: BLE001
            rows.append(_empty_chip_row(q, chosen_exp, str(exc)[:200], windows))
    return rows


def _state_short(state):
    if not state:
        return ""
    return state.split("_", 1)[-1]


def make_figures(cldn, rfc2, level_rows, chrom_rows, neighborhood):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    window_path = os.path.join(OUT, "domains_in_window.tsv")
    window_rows = []
    if os.path.exists(window_path):
        with open(window_path) as fh:
            window_rows = list(csv.DictReader(fh, delimiter="\t"))

    # Fixed tracks so the one spanning HCT116 call is visible next to the
    # fine-scale calls that keep CLDN4 and RFC2 apart.
    track_plan = [
        ("A549 parental", "ENCSR444WCZ"),
        ("A549 0 h dex", "ENCSR662QKG"),
        ("HCT116 intact", "ENCSR477GZK"),
        ("HCT116 in situ", "ENCSR637QCS"),
        ("HMEC intact", "ENCSR707XVJ"),
        ("MCF-7 intact", "ENCSR660LPJ"),
        ("HepG2 intact", "ENCSR888DEJ"),
        ("MCF10A intact", "ENCSR370TFL"),
        ("keratinocyte", "ENCSR869CSI"),
        ("esophagus mucosa", "ENCSR015CVF"),
        ("breast epithelium", "ENCSR326YHP"),
        ("Panc1 intact", "ENCSR584RBV"),
    ]
    fig, ax = plt.subplots(figsize=(11.4, 7.4))
    x0 = 73_200_000
    x1 = 74_600_000
    label_genes = {"CLDN3", "CLDN4", "ELN", "RFC2", "LIMK1"}
    genes = [g for g in neighborhood if g["symbol"] and overlaps(g["start"], g["end"], x0, x1)]
    for g in genes:
        color = "#b2182b" if g["symbol"] == "CLDN4" else ("#2166ac" if g["symbol"] == "RFC2" else "#737373")
        ax.add_patch(Rectangle((g["start"], 0.2), max(g["end"] - g["start"], 1), 0.45, color=color, lw=0))
        if g["symbol"] in label_genes:
            ax.text((g["start"] + g["end"]) / 2, 0.78, g["symbol"], ha="center", va="bottom", fontsize=8, color=color)
    tracks = []
    for label, exp in track_plan:
        rows = [r for r in window_rows if r["experiment"] == exp and r["stratum"] == "baseline"]
        if rows:
            tracks.append((label, exp, rows))
    for i, (label, exp, rows) in enumerate(tracks):
        y = 2.2 + i * 0.85
        ax.text(x0 - 30000, y + 0.2, label, ha="right", va="center", fontsize=8)
        seen_iv = set()
        for r in rows:
            s, e = int(r["domain_start0"]), int(r["domain_end0"])
            key = (s, e)
            if key in seen_iv:
                continue
            seen_iv.add(key)
            covers = r["covers_CLDN4"] == "True" and r["covers_RFC2"] == "True" and int(r["domain_bp"]) <= DOMAIN_MAX
            color = "#2166ac" if covers else "#bdbdbd"
            ax.add_patch(
                Rectangle((max(s, x0), y), max(min(e, x1) - max(s, x0), 1), 0.5, color=color, lw=0)
            )
    ax.axvline(cldn["tss"], color="#b2182b", lw=0.7, ls="--")
    ax.axvline(rfc2["tss"], color="#2166ac", lw=0.7, ls="--")
    ax.set_xlim(x0, x1)
    ax.set_ylim(0, 2.2 + len(tracks) * 0.85 + 0.4)
    ax.set_yticks([])
    ax.set_xlabel("chr7, GRCh38")
    ax.set_title("Blue domain overlaps both CLDN4 and RFC2 (length ≤ 2 Mb)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_locus_domains.png"), dpi=160)
    plt.close()
    _ = level_rows

    # chromHMM heatmap: primary epithelial samples, two TSS columns
    primary = [r for r in chrom_rows if r["panel"] == "epithelial_cell_or_line" and r["state"]]
    eids = []
    for r in primary:
        if r["eid"] not in eids:
            eids.append(r["eid"])
    queries = []
    for r in primary:
        if r["query"] not in queries:
            queries.append(r["query"])
    if not eids or not queries:
        return
    state_colors = {
        "TssA": "#ff0000",
        "TssAFlnk": "#ff4500",
        "TxFlnk": "#32cd32",
        "Tx": "#008000",
        "TxWk": "#006400",
        "EnhG": "#c2e105",
        "Enh": "#ffff00",
        "ZNF/Rpts": "#66cdaa",
        "Het": "#8a91d0",
        "TssBiv": "#cd5c5c",
        "BivFlnk": "#e9967a",
        "EnhBiv": "#bdb76b",
        "ReprPC": "#808080",
        "ReprPCWk": "#c0c0c0",
        "Quies": "#f5f5f5",
    }
    lookup = {(r["eid"], r["query"]): _state_short(r["state"]) for r in primary}
    names = {r["eid"]: r["std_name"] for r in primary}
    fig, ax = plt.subplots(figsize=(8.4, 0.46 * len(eids) + 1.6))
    for i, eid in enumerate(eids):
        for j, q in enumerate(queries):
            st = lookup.get((eid, q), "")
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=state_colors.get(st, "#ffffff"), edgecolor="white", lw=1.2))
            ax.text(j + 0.5, i + 0.5, st, ha="center", va="center", fontsize=7)
        ax.text(-0.1, i + 0.5, f"{eid}  {names.get(eid, '')}", ha="right", va="center", fontsize=7)
    ax.set_xlim(0, len(queries))
    ax.set_ylim(len(eids), 0)
    ax.set_xticks([j + 0.5 for j in range(len(queries))])
    ax.set_xticklabels(queries, rotation=30, ha="right")
    ax.set_yticks([])
    ax.set_title("Roadmap 15-state chromHMM at the TSS (hg38 lift)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_chromhmm_tss.png"), dpi=160)
    plt.close()


def collect_window_domains(manifest, cldn, rfc2):
    """Second pass is avoided: slice domains from files already listed as used.

    Domain files are small. Re-read baseline epithelial domain files and keep
    intervals overlapping the CLDN4–RFC2 window for the locus figure.
    """
    x0, x1 = 73_200_000, 74_600_000
    rows = []
    files = [
        m
        for m in manifest
        if m.get("used") == "yes" and m.get("kind") == "domains" and m.get("stratum") == "baseline"
    ]
    print(f"Window slice from {len(files)} baseline domain files", flush=True)
    for m in files:
        href = None
        # href not stored; reconstruct ENCODE download path from format
        ext = "bedpe.gz" if m["file_format"] == "bedpe" else "bed.gz"
        # tsv domains are not expected. bed files are gzipped on the download URL.
        href = f"https://www.encodeproject.org/files/{m['file']}/@@download/{m['file']}.{ext}"
        try:
            raw = http_bytes(href)
        except Exception:
            alt = f"https://www.encodeproject.org/files/{m['file']}/@@download/{m['file']}.{m['file_format']}"
            try:
                raw = http_bytes(alt)
                href = alt
            except Exception as exc:  # noqa: BLE001
                print(f"  window skip {m['file']} {exc}", flush=True)
                continue
        fh = open_maybe_gzip(raw, href)
        for line in io.TextIOWrapper(fh, encoding="utf-8", errors="replace"):
            recs = parse_intervals_from_line(line, m["file_format"], "domains")
            if not recs:
                continue
            for chrom, start, end in recs:
                if chrom not in ("chr7", "7"):
                    continue
                if not overlaps(start, end, x0, x1):
                    continue
                rows.append(
                    {
                        "experiment": m["experiment"],
                        "assay": m["assay"],
                        "biosample": m["biosample"],
                        "stratum": m["stratum"],
                        "file": m["file"],
                        "domain_start0": start,
                        "domain_end0": end,
                        "domain_bp": end - start,
                        "covers_CLDN4": overlaps(start, end, cldn["start"], cldn["end"]),
                        "covers_RFC2": overlaps(start, end, rfc2["start"], rfc2["end"]),
                    }
                )
    return rows


def write_hic_loop_distances(manifest, cldn, rfc2):
    """For baseline Hi-C loops that touch CLDN4, record the partner gap to RFC2.

    This is a near-miss check. It does not change the pre-specified hit rule
    (partner must overlap the repair gene ± 5 kb).
    """
    if not rfc2:
        return []
    files = [
        r
        for r in manifest
        if r.get("used") == "yes"
        and r.get("kind") == "loops"
        and r.get("stratum") == "baseline"
        and r.get("assay") != "ChIA-PET"
        and int(r.get("n_cldn_loops") or 0) > 0
    ]
    c1, c2 = cldn["start"] - LOOP_FLANK, cldn["end"] + LOOP_FLANK
    rows = []
    print(f"Hi-C near-miss scan on {len(files)} CLDN4-touching loop files", flush=True)
    for m in files:
        ext = "bedpe.gz" if m["file_format"] == "bedpe" else ("tsv" if m["file_format"] == "tsv" else "bed.gz")
        href = f"https://www.encodeproject.org/files/{m['file']}/@@download/{m['file']}.{ext}"
        try:
            raw = http_bytes(href)
        except Exception:
            href = f"https://www.encodeproject.org/files/{m['file']}/@@download/{m['file']}.{m['file_format']}"
            raw = http_bytes(href)
        fh = open_maybe_gzip(raw, href)
        min_gap = None
        n = 0
        closest = None
        n_trans = 0
        for line in io.TextIOWrapper(fh, encoding="utf-8", errors="replace"):
            recs = parse_intervals_from_line(line, m["file_format"], "loops")
            if not recs or len(recs) != 2:
                continue
            (a_c, a_s, a_e), (b_c, b_s, b_e) = recs
            a_hit = a_c == "chr7" and overlaps(a_s, a_e, c1, c2)
            b_hit = b_c == "chr7" and overlaps(b_s, b_e, c1, c2)
            if not a_hit and not b_hit:
                continue
            n += 1
            partners = []
            if a_hit:
                partners.append((b_c, b_s, b_e))
            if b_hit:
                partners.append((a_c, a_s, a_e))
            for pc, ps, pe in partners:
                if pc != "chr7":
                    n_trans += 1
                    continue
                g = gap(ps, pe, rfc2["start"], rfc2["end"])
                if min_gap is None or g < min_gap:
                    min_gap = g
                    closest = (pc, ps, pe)
        rows.append(
            {
                "experiment": m["experiment"],
                "assay": m["assay"],
                "biosample": m["biosample"],
                "file": m["file"],
                "n_cldn_loops": n,
                "n_trans_partners": n_trans,
                "min_gap_to_RFC2_bp": "" if min_gap is None else min_gap,
                "closest_partner": "" if closest is None else f"{closest[0]}:{closest[1]}-{closest[2]}",
            }
        )
    if rows:
        write_tsv(os.path.join(OUT, "hic_cldn4_loop_distance_to_rfc2.tsv"), rows, list(rows[0].keys()))
    return rows


def summarize(cldn, dist_rows, level_rows, chrom_rows, chip_rows, neighborhood, kegg_n):
    def _near(rows, bp):
        return [r for r in rows if r["gap_bp"] != "" and int(r["gap_bp"]) <= bp and r["symbol"] != "CLDN4"]

    near_1 = [r for r in dist_rows if r["within_1Mb_gap"] in (True, "True")]
    # CLDN4 itself is not in the repair set.
    def baseline_hi(kind):
        return [
            r
            for r in level_rows
            if r["kind"] == kind and r["stratum"] == "baseline" and r["assay"] != "ChIA-PET"
        ]

    def baseline_chia(kind):
        return [r for r in level_rows if r["kind"] == kind and r["stratum"] == "baseline" and r["assay"] == "ChIA-PET"]

    hi_domains = baseline_hi("domains")
    hi_loops = baseline_hi("loops")
    chia_loops = baseline_chia("loops")

    def count_hit(rows, field):
        return sum(1 for r in rows if str(r.get(field)) in ("True", "true", True))

    def genes_of(rows, field):
        genes = set()
        for r in rows:
            if str(r.get(field)) in ("True", "true", True):
                for g in (r.get("repair_genes_le_2Mb") or r.get("repair_genes") or "").split(";"):
                    if g:
                        genes.add(g)
        return sorted(genes)

    # lung lines
    lung_bios = {"A549", "NCI-H460"}

    def lung(rows):
        return [r for r in rows if r["biosample"] in lung_bios]

    chrom_primary = [r for r in chrom_rows if r["panel"] == "epithelial_cell_or_line"]
    def state_counts(query):
        c = defaultdict(int)
        for r in chrom_primary:
            if r["query"] == query and r["state"]:
                c[_state_short(r["state"])] += 1
        return dict(c)

    chip_ok = [r for r in chip_rows if r.get("file")]
    both = sum(int(r.get("peaks_overlapping_both_promoters_le_50kb") or 0) for r in chip_ok)

    local_genes = []
    for g in neighborhood:
        if not g["symbol"]:
            continue
        g_bp = gap(cldn["start"], cldn["end"], g["start"], g["end"])
        if g_bp <= NEAR_BP:
            local_genes.append((g_bp, g["symbol"], g["description"]))
    local_genes.sort()

    summary = {
        "cldn4": cldn,
        "kegg_genes_with_coords": kegg_n,
        "nearest_repair": dist_rows[:8],
        "n_repair_within_1Mb": len(near_1),
        "n_repair_within_2Mb": sum(1 for r in dist_rows if r["within_2Mb_gap"] in (True, "True")),
        "baseline_hic_domain_experiments": len(hi_domains),
        "baseline_hic_domain_hit_any": count_hit(hi_domains, "repair_gene_hit"),
        "baseline_hic_domain_hit_le2mb": count_hit(hi_domains, "domain_le_2Mb"),
        "baseline_hic_loop_experiments": len(hi_loops),
        "baseline_hic_loop_gene_hit": count_hit(hi_loops, "repair_gene_hit"),
        "baseline_chia_loop_experiments": len(chia_loops),
        "baseline_chia_loop_gene_hit": count_hit(chia_loops, "repair_gene_hit"),
        "lung_domain_experiments": len(lung(hi_domains)),
        "lung_domain_le2mb": count_hit(lung(hi_domains), "domain_le_2Mb"),
        "lung_loop_experiments": len(lung(hi_loops)),
        "lung_loop_hit": count_hit(lung(hi_loops), "repair_gene_hit"),
        "chromhmm_primary_n": len({r["eid"] for r in chrom_primary if r["state"]}),
        "chromhmm_CLDN4": state_counts("CLDN4_TSS"),
        "chromhmm_nearest": {},
        "chip_files_scored": len(chip_ok),
        "chip_both_promoter_peaks_le_50kb": both,
        "local_genes_1Mb": [{"gap_bp": g, "symbol": s, "description": d} for g, s, d in local_genes],
    }
    # nearest repair symbol for chromhmm key
    if near_1:
        sym = near_1[0]["symbol"]
        summary["nearest_symbol"] = sym
        summary["nearest_gap_bp"] = int(near_1[0]["gap_bp"])
        summary["nearest_tss_bp"] = int(near_1[0]["tss_distance_bp"])
        summary["nearest_pathways"] = near_1[0]["pathways"]
        summary["chromhmm_nearest"] = state_counts(sym + "_TSS")
    else:
        summary["nearest_symbol"] = None
    # per-biosample domain hit list for the write-up
    summary["domain_le2mb_by_biosample"] = sorted(
        {
            (r["biosample"], r["assay"], r["repair_genes_le_2Mb"])
            for r in hi_domains
            if str(r.get("domain_le_2Mb")) in ("True", "true", True)
        }
    )
    summary["loop_hit_by_biosample"] = sorted(
        {
            (r["biosample"], r["assay"], r["target"] or "", r["repair_genes"])
            for r in hi_loops + chia_loops
            if str(r.get("repair_gene_hit")) in ("True", "true", True)
        }
    )
    summary["domain_miss_biosamples"] = sorted({r["biosample"] for r in hi_domains if str(r.get("domain_le_2Mb")) not in ("True", "true", True)})
    summary["hic_loop_biosamples_scanned"] = sorted({r["biosample"] for r in hi_loops})
    summary["hic_domain_biosamples_scanned"] = sorted({r["biosample"] for r in hi_domains})
    return summary


def main():
    cldn = fetch_cldn4()
    repair, unmapped = fetch_kegg_repair_genes()
    print(
        f"CLDN4 chr{cldn['chrom']}:{cldn['start']}-{cldn['end']} (0-based) TSS {cldn['tss']}",
        flush=True,
    )
    dist_rows = distance_table(cldn, repair)
    dist_fields = list(dist_rows[0].keys())
    write_tsv(os.path.join(OUT, "kegg_dna_repair_distance_to_cldn4.tsv"), dist_rows, dist_fields)
    near = [r for r in dist_rows if r["within_1Mb_gap"] in (True, "True")]
    print("Repair genes with gap <= 1 Mb:", [(r["symbol"], r["gap_bp"]) for r in near], flush=True)
    neighborhood = neighborhood_genes(cldn["chrom"], cldn["start"], cldn["end"], pad=1_500_000)
    write_tsv(
        os.path.join(OUT, "protein_coding_within_1p5Mb.tsv"),
        [
            {
                **g,
                "gap_to_CLDN4_bp": gap(cldn["start"], cldn["end"], g["start"], g["end"]),
                "in_kegg_dna_repair": g["symbol"] in repair,
            }
            for g in sorted(neighborhood, key=lambda g: g["start"])
        ],
        ["symbol", "ensembl", "chrom", "start", "end", "strand", "tss", "description", "gap_to_CLDN4_bp", "in_kegg_dna_repair"],
    )
    records, manifest, loop_rows, domain_rows = analyze_3d(cldn, repair)
    write_tsv(
        os.path.join(OUT, "encode_file_manifest.tsv"),
        manifest,
        [
            "experiment",
            "assay",
            "biosample",
            "epithelial_call",
            "stratum",
            "target",
            "description",
            "file",
            "kind",
            "file_format",
            "file_size",
            "biological_replicates",
            "used",
            "reason",
            "n_records",
            "n_cldn_loops",
            "n_cldn_domains",
        ],
    )
    if loop_rows:
        write_tsv(os.path.join(OUT, "loop_links_to_dna_repair.tsv"), loop_rows, list(loop_rows[0].keys()))
    else:
        write_tsv(
            os.path.join(OUT, "loop_links_to_dna_repair.tsv"),
            [],
            ["experiment", "assay", "biosample", "repair_gene"],
        )
    if domain_rows:
        write_tsv(os.path.join(OUT, "domain_overlaps_dna_repair.tsv"), domain_rows, list(domain_rows[0].keys()))
    else:
        write_tsv(
            os.path.join(OUT, "domain_overlaps_dna_repair.tsv"),
            [],
            ["experiment", "assay", "biosample", "repair_gene"],
        )
    level = experiment_level(manifest, loop_rows, domain_rows)
    if level:
        write_tsv(os.path.join(OUT, "experiment_summary.tsv"), level, list(level[0].keys()))
    # inventory of experiments that were not epithelial, for the audit
    inv = []
    for rec in records:
        inv.append(
            {
                "experiment": rec["experiment"],
                "assay": rec["assay"],
                "biosample": rec["biosample"],
                "classification": rec["classification"],
                "organ_slims": ";".join(rec["organ_slims"]) if isinstance(rec["organ_slims"], list) else rec["organ_slims"],
                "epithelial_call": rec["epithelial_call"] or "",
                "stratum": rec["stratum"],
                "target": rec["target"] or "",
                "n_loop_or_domain_files": len(rec["files"]),
                "description": (rec["description"] or "").replace("\t", " ")[:240],
            }
        )
    write_tsv(
        os.path.join(OUT, "encode_experiment_inventory.tsv"),
        inv,
        list(inv[0].keys()) if inv else ["experiment"],
    )
    repair_near_genes = [repair[r["symbol"]] for r in near if r["symbol"] in repair]
    chrom_rows = roadmap_states(cldn, repair_near_genes)
    write_tsv(
        os.path.join(OUT, "roadmap_chromhmm_tss.tsv"),
        chrom_rows,
        list(chrom_rows[0].keys()) if chrom_rows else ["eid"],
    )
    rfc2 = repair.get("RFC2")
    chip_rows = []
    if rfc2:
        chip_rows = chip_peak_survey(cldn, rfc2)
        write_tsv(os.path.join(OUT, "encode_chip_tss_peaks.tsv"), chip_rows, list(chip_rows[0].keys()))
    if rfc2:
        write_hic_loop_distances(manifest, cldn, rfc2)
    window_rows = collect_window_domains(manifest, cldn, rfc2) if rfc2 else []
    if window_rows:
        write_tsv(os.path.join(OUT, "domains_in_window.tsv"), window_rows, list(window_rows[0].keys()))
    summary = summarize(cldn, dist_rows, level, chrom_rows, chip_rows, neighborhood, len(repair))
    summary["kegg_unmapped_entrez"] = unmapped
    summary["n_experiments_seen"] = len(records)
    summary["n_files_used"] = sum(1 for m in manifest if m.get("used") == "yes")
    # JSON-friendly cldn
    with open(os.path.join(OUT, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)
    try:
        make_figures(cldn, rfc2, level, chrom_rows, neighborhood)
    except Exception as exc:  # noqa: BLE001
        print("figure failed", exc, flush=True)
    print(json.dumps({k: summary[k] for k in summary if k not in ("nearest_repair", "local_genes_1Mb", "cldn4")}, indent=2))


if __name__ == "__main__":
    main()
