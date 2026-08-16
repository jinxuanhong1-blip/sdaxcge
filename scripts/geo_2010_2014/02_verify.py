"""Verify EVERY candidate accession against authoritative GEO records.

For each GSE from 01_search.py this fetches:
  * the SOFT brief record (status, dates, type, pubmed, sample ids)
  * the suppl/ and matrix/ FTP listings (processed files + sizes)

Then classifies organism / lung / checkpoint-ICI / broader immunotherapy,
and flags leftover status (this window has no prior GEO ICI PR to subtract).

Outputs:
  verified_series.tsv / .json
  suppl_files.tsv
  leftover_shortlist.tsv
"""
import html
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from eutils import _get  # noqa: E402

ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "w200", "GEO_2010_2014")

SIZE_RE = re.compile(
    r'<a href="([^"/][^"]*)">[^<]*</a>\s+[\d-]+\s+[\d:]+\s+([\d.]+[KMGT]?)',
    re.I,
)

# Already-mined GEO ICI windows in this repo (later years). None of those
# accessions can appear in 2010-2014; listed so leftover is explicit.
ALREADY_COVERED_WINDOWS = [
    "results/fable_geo_2015_2018/",
    "results/fable_geo_2019_2021/",
    "results/fable_geo_2022_2023/",
    "results/fable_geo_2026/",
]


def parse_size(s):
    s = s.strip()
    if not s or s == "-":
        return 0
    mult = {"K": 1024, "M": 1024**2, "G": 1024**3, "T": 1024**4}
    if s[-1] in mult:
        return int(float(s[:-1]) * mult[s[-1]])
    return int(float(s))


def human(n):
    n = float(n)
    for unit in ["B", "K", "M", "G", "T"]:
        if n < 1024 or unit == "T":
            return f"{n:.1f}{unit}"
        n /= 1024


def ftp_listing(acc):
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
    url = (
        f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}"
        f"&targ=self&form=text&view=brief"
    )
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
        elif line.startswith("!Series_platform_id"):
            d.setdefault("platforms_soft", []).append(line.split("=", 1)[1].strip())
    return d


LUNG_RE = re.compile(
    r"\b(lung|nsclc|pulmonary|adenocarcinoma of the lung|lung adenocarcinoma|"
    r"squamous cell lung|lung carcinoma|lung cancer|lung tumor|lung neoplasm)\b",
    re.I,
)
# Broader immunotherapy (cancer vaccine, HLA epitopes, DNA vaccine, etc.).
IO_RE = re.compile(
    r"(immunotherap|immune checkpoint|checkpoint blockade|checkpoint inhibitor|"
    r"anti-pd|nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|"
    r"ipilimumab|tremelimumab|mage-?a3|cancer vaccine|t[- ]?cell epitope)",
    re.I,
)
# Strict PD-1 / PD-L1 / CTLA-4 checkpoint ICI (drugs or blockade, not just the gene).
CHECKPOINT_RE = re.compile(
    r"(nivolumab|pembrolizumab|atezolizumab|durvalumab|avelumab|ipilimumab|"
    r"tremelimumab|anti-pd-?1|anti-pd-?l1|anti-ctla|checkpoint blockade|"
    r"checkpoint inhibitor|immune checkpoint|pd-?1 blockade|pd-?l1 blockade)",
    re.I,
)
# Gene-level PD-L1 mention without treatment (e.g. "elevated Pdl1 expression").
PDL1_GENE_RE = re.compile(r"\b(pd-?l1|pdl1|cd274|pd-?1|pdcd1|ctla-?4|ctla4)\b", re.I)
FALSE_PD_RE = re.compile(r"pd[\s-]?0332991|palbociclib", re.I)
TREAT_RE = re.compile(
    r"treated with|treatment|response|responder|non-?responder|"
    r"clinical outcome|survival|blockade|anti-pd|nivolumab|pembrolizumab|"
    r"atezolizumab|durvalumab|ipilimumab|mage-?a3",
    re.I,
)


def classify(taxon, title, summary):
    text = f"{title} {summary}"
    is_human = "Homo sapiens" in taxon
    is_lung = bool(LUNG_RE.search(text))
    only_false_pd = bool(FALSE_PD_RE.search(text)) and not CHECKPOINT_RE.search(text)
    is_io = bool(IO_RE.search(text)) and not only_false_pd
    is_checkpoint_ici = bool(CHECKPOINT_RE.search(text)) and not only_false_pd
    is_pdl1_gene = bool(PDL1_GENE_RE.search(text))
    treat = bool(TREAT_RE.search(text))
    return is_human, is_lung, is_io, is_checkpoint_ici, is_pdl1_gene, treat


def fallout_reason(is_human, is_lung, is_io, is_checkpoint_ici):
    if is_human and is_lung and is_checkpoint_ici:
        return "in_scope_checkpoint_ici"
    if is_human and is_lung and is_io:
        return "leftover_human_lung_immunotherapy_not_checkpoint_ici"
    if (not is_human) and is_lung and is_pdl1_gene:
        return "mouse_lung_pdl1_expression_not_ici_treatment"
    if (not is_human) and is_lung:
        return "mouse_lung_not_ici"
    if is_human and not is_lung:
        return "human_not_lung"
    if not is_human:
        return "non_human"
    return "human_lung_not_io"


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
        is_human, is_lung, is_io, is_checkpoint_ici, is_pdl1_gene, treat = classify(
            taxon, title, summary
        )
        reason = fallout_reason(is_human, is_lung, is_io, is_checkpoint_ici)
        rec = {
            "accession": acc,
            "uid": uid,
            "taxon": taxon,
            "is_human": is_human,
            "is_lung": is_lung,
            "is_immunotherapy": is_io,
            "is_checkpoint_ici": is_checkpoint_ici,
            "is_pdl1_gene_mention": is_pdl1_gene,
            "likely_response_labels": treat,
            "leftover_reason": reason,
            "already_covered_by_later_window": False,
            "status": sb.get("status", ""),
            "submission_date": sb.get("submission_date", ""),
            "last_update_date": sb.get("last_update_date", ""),
            "types": "; ".join(sb.get("types", [])),
            "pubmed": sb.get("pubmed", "") or ";".join(str(p) for p in s.get("pubmedids", [])),
            "n_samples": s.get("n_samples", sb.get("n_samples_soft", "")),
            "gpl": s.get("gpl", ""),
            "pdat": s.get("pdat", ""),
            "n_suppl_files": len(files),
            "suppl_total_size": human(total),
            "suppl_total_bytes": total,
            "title": title,
        }
        records.append(rec)
        for sub, name, size_b, size_h in files:
            supp_rows.append(
                {
                    "accession": acc,
                    "location": sub,
                    "filename": name,
                    "size": size_h,
                    "size_bytes": size_b,
                    "under_2gb": size_b < 2 * 1024**3,
                }
            )

    with open(os.path.join(OUT, "verified_series.json"), "w") as f:
        json.dump(records, f, indent=1)

    cols = [
        "accession",
        "taxon",
        "is_human",
        "is_lung",
        "is_immunotherapy",
        "is_checkpoint_ici",
        "is_pdl1_gene_mention",
        "likely_response_labels",
        "leftover_reason",
        "status",
        "submission_date",
        "last_update_date",
        "types",
        "pubmed",
        "n_samples",
        "gpl",
        "pdat",
        "n_suppl_files",
        "suppl_total_size",
        "suppl_total_bytes",
        "title",
    ]
    with open(os.path.join(OUT, "verified_series.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in records:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    scols = ["accession", "location", "filename", "size", "size_bytes", "under_2gb"]
    with open(os.path.join(OUT, "suppl_files.tsv"), "w") as f:
        f.write("\t".join(scols) + "\n")
        for r in supp_rows:
            f.write("\t".join(str(r[c]) for c in scols) + "\n")

    # Leftover shortlist: everything that is not a later-window duplicate
    # (none are) plus an honest scope tag.
    short_cols = [
        "accession",
        "leftover_reason",
        "is_human",
        "is_lung",
        "is_immunotherapy",
        "is_checkpoint_ici",
        "n_samples",
        "gpl",
        "pdat",
        "pubmed",
        "title",
    ]
    with open(os.path.join(OUT, "leftover_shortlist.tsv"), "w") as f:
        f.write("\t".join(short_cols) + "\n")
        for r in records:
            f.write("\t".join(str(r[c]) for c in short_cols) + "\n")

    n = len(records)
    n_h = sum(r["is_human"] for r in records)
    n_hl = sum(r["is_human"] and r["is_lung"] for r in records)
    n_hlio = sum(r["is_human"] and r["is_lung"] and r["is_immunotherapy"] for r in records)
    n_hli = sum(r["is_human"] and r["is_lung"] and r["is_checkpoint_ici"] for r in records)
    print(
        f"\nverified {n} series | human={n_h} | human+lung={n_hl} | "
        f"human+lung+IO={n_hlio} | human+lung+checkpoint_ICI={n_hli}"
    )
    print("-> results/w200/GEO_2010_2014/verified_series.tsv")


if __name__ == "__main__":
    main()
