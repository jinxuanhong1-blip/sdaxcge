"""Verify EVERY candidate accession against authoritative GEO records.

For each GSE from 01_search.py this fetches:
  * the SOFT "brief" record (status/public date, submission date, last update,
    series type, pubmed id, sample ids) -> proves the accession is REAL and
    records its true metadata (no invented IDs);
  * the suppl/ and matrix/ FTP directory listings -> enumerates processed files
    and their sizes so we can apply the "open processed <2GB" download rule.

It then classifies each series by organism (human?), disease (lung?), and
immune-checkpoint relevance (ICI / PD-1 / PD-L1 / CTLA-4?), and flags likely
treatment-response labels.

Outputs:
  results/fable_geo_2015_2018/verified_series.tsv
  results/fable_geo_2015_2018/verified_series.json
  results/fable_geo_2015_2018/suppl_files.tsv
"""
import json
import os
import re
import sys
import html

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eutils import _get  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "fable_geo_2015_2018")

SIZE_RE = re.compile(r'<a href="([^"/][^"]*)">[^<]*</a>\s+[\d-]+\s+[\d:]+\s+([\d.]+[KMGT]?)', re.I)


def parse_size(s):
    s = s.strip()
    if not s or s == "-":
        return 0
    mult = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    if s[-1] in mult:
        return int(float(s[:-1]) * mult[s[-1]])
    return int(float(s))


def human(n):
    for unit in ["B", "K", "M", "G", "T"]:
        if n < 1024 or unit == "T":
            return f"{n:.1f}{unit}"
        n /= 1024


def ftp_listing(acc):
    """Return list of (filename, size_bytes) for suppl/ + matrix/ of a series."""
    prefix = acc[:-3] + "nnn" if len(acc) > 6 else acc[:3] + "nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{acc}"
    files = []
    for sub in ("suppl", "matrix"):
        url = f"{base}/{sub}/"
        try:
            page = _get(url)
        except Exception:
            continue
        for m in SIZE_RE.finditer(page):
            name, size = m.group(1), m.group(2)
            if name in ("Parent Directory",):
                continue
            files.append((sub, html.unescape(name), parse_size(size), size))
    return files


def soft_brief(acc):
    url = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
           f"&targ=self&form=text&view=brief")
    txt = _get(url)
    d = {"n_samples_soft": 0, "pubmed": ""}
    for line in txt.splitlines():
        if line.startswith("!Series_status"):
            d["status"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_submission_date"):
            d["submission_date"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_last_update_date"):
            d["last_update_date"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_type"):
            d.setdefault("types", []).append(line.split("=", 1)[1].strip())
        elif line.startswith("!Series_pubmed_id"):
            d["pubmed"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_sample_id"):
            d["n_samples_soft"] += 1
        elif line.startswith("!Series_title"):
            d["title_soft"] = line.split("=", 1)[1].strip()
    return d


LUNG_RE = re.compile(r"\b(lung|nsclc|pulmonary|adenocarcinoma of the lung|lung adenocarcinoma|"
                     r"squamous cell lung|lung carcinoma|lung cancer|lung tumor|lung neoplasm)\b", re.I)
ICI_RE = re.compile(r"(pd-?1|pdcd1|pd-?l1|cd274|ctla-?4|nivolumab|pembrolizumab|atezolizumab|"
                    r"durvalumab|avelumab|ipilimumab|tremelimumab|immune checkpoint|"
                    r"checkpoint blockade|checkpoint inhibitor|anti-pd|immunotherap)", re.I)
# 'PD 0332991' (palbociclib) is NOT an anti-PD-1 agent; guard against that false hit.
FALSE_PD_RE = re.compile(r"pd[\s-]?0332991|palbociclib", re.I)
TREAT_RE = re.compile(r"treated with|treatment|response|responder|non-?responder|"
                      r"clinical outcome|survival|blockade|anti-pd|nivolumab|pembrolizumab|"
                      r"atezolizumab|durvalumab|ipilimumab", re.I)


def classify(taxon, title, summary):
    text = f"{title} {summary}"
    is_human = "Homo sapiens" in taxon
    is_lung = bool(LUNG_RE.search(text))
    ici_hit = bool(ICI_RE.search(text))
    only_false_pd = bool(FALSE_PD_RE.search(text)) and not re.search(
        r"pd-?l1|pdcd1|ctla|nivolumab|pembrolizumab|atezolizumab|durvalumab|"
        r"avelumab|ipilimumab|tremelimumab|immune checkpoint|checkpoint blockade|"
        r"checkpoint inhibitor|anti-pd-?1|anti-pd-?l1|immunotherap", text, re.I)
    is_ici = ici_hit and not only_false_pd
    treat = bool(TREAT_RE.search(text))
    return is_human, is_lung, is_ici, treat


def main():
    raw = json.load(open(os.path.join(OUT, "search_raw_summaries.json")))
    summaries = raw["summaries"]
    series = {uid: s for uid, s in summaries.items() if s.get("entrytype") == "GSE"}

    records = []
    supp_rows = []
    for uid, s in sorted(series.items(), key=lambda kv: kv[1].get("accession", "")):
        acc = s["accession"]
        taxon = s.get("taxon", "")
        title = s.get("title", "")
        summary = s.get("summary", "")
        print(f"verifying {acc} ...", flush=True)
        sb = soft_brief(acc)
        files = ftp_listing(acc)
        total = sum(f[2] for f in files)
        is_human, is_lung, is_ici, treat = classify(taxon, title, summary)
        rec = {
            "accession": acc,
            "uid": uid,
            "taxon": taxon,
            "is_human": is_human,
            "is_lung": is_lung,
            "is_ici": is_ici,
            "likely_response_labels": treat,
            "status": sb.get("status", ""),
            "submission_date": sb.get("submission_date", ""),
            "last_update_date": sb.get("last_update_date", ""),
            "types": "; ".join(sb.get("types", [])),
            "pubmed": sb.get("pubmed", ""),
            "n_samples": s.get("n_samples", sb.get("n_samples_soft", "")),
            "gpl": s.get("gpl", ""),
            "n_suppl_files": len(files),
            "suppl_total_size": human(total),
            "suppl_total_bytes": total,
            "title": title,
        }
        records.append(rec)
        for sub, name, size_b, size_h in files:
            supp_rows.append({
                "accession": acc, "location": sub, "filename": name,
                "size": size_h, "size_bytes": size_b,
                "under_2gb": size_b < 2 * 1024**3,
            })

    with open(os.path.join(OUT, "verified_series.json"), "w") as f:
        json.dump(records, f, indent=1)

    cols = ["accession", "taxon", "is_human", "is_lung", "is_ici",
            "likely_response_labels", "status", "submission_date",
            "last_update_date", "types", "pubmed", "n_samples", "gpl",
            "n_suppl_files", "suppl_total_size", "suppl_total_bytes", "title"]
    with open(os.path.join(OUT, "verified_series.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in records:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    scols = ["accession", "location", "filename", "size", "size_bytes", "under_2gb"]
    with open(os.path.join(OUT, "suppl_files.tsv"), "w") as f:
        f.write("\t".join(scols) + "\n")
        for r in supp_rows:
            f.write("\t".join(str(r[c]) for c in scols) + "\n")

    n_h = sum(r["is_human"] for r in records)
    n_hl = sum(r["is_human"] and r["is_lung"] for r in records)
    n_hli = sum(r["is_human"] and r["is_lung"] and r["is_ici"] for r in records)
    print(f"\nverified {len(records)} series | human={n_h} | human+lung={n_hl} | "
          f"human+lung+ICI={n_hli}")
    print("-> results/fable_geo_2015_2018/verified_series.tsv")


if __name__ == "__main__":
    main()
