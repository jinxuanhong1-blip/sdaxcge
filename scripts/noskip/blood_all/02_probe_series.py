#!/usr/bin/env python3
"""Probe GEO series: blood vs tumor, ICI, response labels, supplements, gene presence hints."""
from __future__ import annotations

import gzip
import json
import re
import time
import urllib.request
from pathlib import Path

OUT = Path("results/noskip/blood_all")
CACHE = Path("/tmp/geo_blood/soft")
CACHE.mkdir(parents=True, exist_ok=True)

BLOOD_RE = re.compile(
    r"\b(blood|pbmc|pbmcs|plasma|serum|paxgene|platelet|tep|circulating|"
    r"peripheral blood|whole blood|buffy)\b",
    re.I,
)
TUMOR_RE = re.compile(r"\b(tumor|tumour|biopsy|resect|FFPE|tissue|malignant)\b", re.I)
ICI_RE = re.compile(
    r"(PD-?1|PD-?L1|CTLA-?4|nivolumab|pembrolizumab|atezolizumab|durvalumab|"
    r"avelumab|ipilimumab|cemiplimab|tislelizumab|camrelizumab|sintilimab|"
    r"toripalimab|immune checkpoint|immunotherapy|\bICI\b|\bICB\b|anti-PD)",
    re.I,
)
LUNG_RE = re.compile(
    r"(NSCLC|SCLC|lung cancer|lung adenocarcinoma|lung squamous|"
    r"non-small.cell lung|small.cell lung|\bLUAD\b|\bLUSC\b)",
    re.I,
)
RESP_RE = re.compile(
    r"(responder|non-?responder|RECIST|\bCR\b|\bPR\b|\bSD\b|\bPD\b|"
    r"durable|DCB|NDB|MPR|pCR|PFS|progress)",
    re.I,
)
GENE_EXPR_RE = re.compile(
    r"(Expression profiling|RNA-seq|RNA seq|transcriptom|microarray|scRNA|single-cell RNA)",
    re.I,
)
NOT_GENE_RE = re.compile(
    r"(microRNA|miRNA|non-coding RNA|methylation|ATAC|ChIP|TCR-seq|"
    r"flow cytometr|mass cytometr|CyTOF|CUT&Tag|5-hydroxymethyl)",
    re.I,
)


def fetch(url: str, dest: Path, retries: int = 4) -> Path:
    if dest.exists() and dest.stat().st_size > 200:
        return dest
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "noskip-blood-ici/1.0"})
            with urllib.request.urlopen(req, timeout=90) as r:
                dest.write_bytes(r.read())
            return dest
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"{url}: {last}")


def parse_miniml(acc: str) -> dict:
    """Use GEO accession brief XML via esearch/efetch gds, plus series matrix header if possible."""
    # NCBI GEO Miniml
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=self&form=xml&view=quick"
    dest = CACHE / f"{acc}.xml"
    try:
        fetch(url, dest)
        text = dest.read_text(errors="replace")
    except Exception as e:
        return {"accession": acc, "error": str(e)}

    def grab(tag: str) -> str:
        m = re.search(rf"<{tag}[^>]*>(.*?)</{tag}>", text, re.S | re.I)
        return re.sub(r"<[^>]+>", "", m.group(1)).strip() if m else ""

    title = grab("Title") or grab("title")
    summary = grab("Summary") or grab("summary")
    # characteristics from sample snippets
    chars = re.findall(r"<Characteristic[^>]*>(.*?)</Characteristic>", text, re.S | re.I)
    chars = [re.sub(r"<[^>]+>", " ", c).strip() for c in chars]
    return {
        "accession": acc,
        "xml_title": title,
        "xml_summary": summary[:2000],
        "xml_n_characteristic_fields": len(chars),
        "xml_characteristic_preview": " | ".join(chars[:30])[:1500],
        "xml_bytes": dest.stat().st_size,
    }


def list_suppl(acc: str) -> list[dict]:
    nnn = acc[3:]
    prefix = f"GSE{nnn[:-3]}nnn" if len(nnn) > 3 else f"GSE{nnn}nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}/suppl/"
    dest = CACHE / f"{acc}_suppl.html"
    try:
        fetch(url, dest)
        html = dest.read_text(errors="replace")
    except Exception as e:
        return [{"error": str(e)}]
    files = []
    for m in re.finditer(r'href="([^"]+)"', html):
        href = m.group(1)
        if href in ("../", "./") or href.endswith("/"):
            continue
        files.append({"file": href.split("/")[-1], "url": urllib_join(url, href)})
    return files


def urllib_join(base: str, href: str) -> str:
    if href.startswith("http"):
        return href
    return base.rstrip("/") + "/" + href.lstrip("/")


def classify(title: str, summary: str, gdstype: str) -> dict:
    blob = f"{title}\n{summary}\n{gdstype}"
    is_blood = bool(BLOOD_RE.search(blob))
    is_tumor = bool(TUMOR_RE.search(blob))
    is_ici = bool(ICI_RE.search(blob))
    is_lung = bool(LUNG_RE.search(blob))
    has_resp_mention = bool(RESP_RE.search(blob))
    gene_expr = bool(GENE_EXPR_RE.search(blob))
    not_gene = bool(NOT_GENE_RE.search(blob)) and not gene_expr
    # if both RNA and miRNA mentioned, keep gene_expr True
    if GENE_EXPR_RE.search(blob) and NOT_GENE_RE.search(blob):
        not_gene = False
        gene_expr = True
    human = "Homo sapiens" in blob or "human" in blob.lower()
    return {
        "is_human": human,
        "is_lung": is_lung,
        "is_ici": is_ici,
        "is_blood_compartment": is_blood,
        "mentions_tumor": is_tumor,
        "mentions_response": has_resp_mention,
        "looks_gene_expression": gene_expr,
        "looks_non_gene_assay": not_gene,
    }


def main() -> None:
    raw = json.loads((OUT / "geo_search_raw.json").read_text())
    rows = raw["rows"]
    # unique GSE
    by_acc = {}
    for r in rows:
        acc = r.get("accession") or ""
        if acc.startswith("GSE") and acc not in by_acc:
            by_acc[acc] = r

    catalog = []
    for i, (acc, r) in enumerate(sorted(by_acc.items()), 1):
        print(f"[{i}/{len(by_acc)}] probe {acc}", flush=True)
        cls = classify(r.get("title", ""), r.get("summary", ""), r.get("gdstype", ""))
        try:
            suppl = list_suppl(acc)
        except Exception as e:
            suppl = [{"error": str(e)}]
        time.sleep(0.2)
        rec = {
            **{k: r.get(k) for k in ("accession", "title", "taxon", "n_samples", "gdstype", "pdat", "suppfile", "ftp", "summary")},
            **cls,
            "suppl_files": [s.get("file") for s in suppl if s.get("file")],
            "suppl_n": sum(1 for s in suppl if s.get("file")),
            "suppl_error": next((s.get("error") for s in suppl if s.get("error")), ""),
        }
        catalog.append(rec)

    (OUT / "probe_catalog.json").write_text(json.dumps(catalog, indent=2))
    keys = [
        "accession", "title", "taxon", "n_samples", "gdstype", "pdat",
        "is_human", "is_lung", "is_ici", "is_blood_compartment", "mentions_tumor",
        "mentions_response", "looks_gene_expression", "looks_non_gene_assay",
        "suppl_n", "suppl_error",
    ]
    lines = ["\t".join(keys)]
    for rec in catalog:
        lines.append("\t".join(str(rec.get(k, "")).replace("\t", " ").replace("\n", " ")[:400] for k in keys))
    (OUT / "probe_catalog.tsv").write_text("\n".join(lines) + "\n")
    print(f"wrote {len(catalog)} series")


if __name__ == "__main__":
    main()
