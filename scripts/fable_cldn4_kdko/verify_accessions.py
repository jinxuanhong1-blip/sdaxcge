"""Verify each candidate GEO accession by pulling sample-level metadata.

For every GSE we fetch the quick GSM view (titles, source, organism,
characteristics) and the series self record (summary, platform, supp
files). We then classify whether the series contains a *true CLDN4
perturbation* (KD/KO/shRNA/siRNA/CRISPR), print a triage table, and
save everything to results/ for the writeup.
"""
import re
import sys
import time
import json
import urllib.request

CANDIDATES = [
    "GSE50927", "GSE22493", "GSE207704", "GSE99415", "GSE99416",
    "GSE99417", "GSE84742", "GSE65107", "GSE11759", "GSE244820",
    "GSE26055", "GSE48443", "GSE60885",
]

GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def fetch(acc, targ):
    url = f"{GEO}?acc={acc}&targ={targ}&form=text&view=quick"
    for i in range(5):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "fable_cldn4_kdko"})
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:  # noqa
            time.sleep(2 * (i + 1))
            last = e
    raise RuntimeError(f"fetch failed {acc} {targ}: {last}")


CLDN4_RE = re.compile(r"cldn4|claudin[- ]?4", re.I)
PERTURB_RE = re.compile(r"knock[- ]?down|knock[- ]?out|shRNA|siRNA|CRISPR|sgRNA|silenc|\bKO\b|\bKD\b|deplet|deficient|-/-|null", re.I)


def parse_series(txt):
    d = {"title": "", "summary": "", "supp": [], "platform": [], "type": ""}
    for line in txt.splitlines():
        if line.startswith("!Series_title"):
            d["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Series_summary"):
            d["summary"] += " " + line.split("=", 1)[1].strip()
        elif line.startswith("!Series_type"):
            d["type"] += line.split("=", 1)[1].strip() + "; "
        elif line.startswith("!Series_supplementary_file"):
            d["supp"].append(line.split("=", 1)[1].strip())
        elif line.startswith("!Series_platform_id"):
            d["platform"].append(line.split("=", 1)[1].strip())
    return d


def parse_samples(txt):
    samples = []
    cur = None
    for line in txt.splitlines():
        if line.startswith("^SAMPLE"):
            if cur:
                samples.append(cur)
            cur = {"gsm": line.split("=", 1)[1].strip(), "title": "",
                   "source": "", "organism": "", "char": []}
        elif cur is None:
            continue
        elif line.startswith("!Sample_title"):
            cur["title"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_source_name_ch1"):
            cur["source"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_organism_ch1"):
            cur["organism"] = line.split("=", 1)[1].strip()
        elif line.startswith("!Sample_characteristics_ch1"):
            cur["char"].append(line.split("=", 1)[1].strip())
    if cur:
        samples.append(cur)
    return samples


def main():
    out = {}
    for acc in CANDIDATES:
        try:
            series = parse_series(fetch(acc, "self"))
            samples = parse_samples(fetch(acc, "gsm"))
        except Exception as e:  # noqa
            print(f"!! {acc}: {e}")
            continue
        # Does any sample mention CLDN4 + perturbation?
        sample_blob = " ".join(
            s["title"] + " " + s["source"] + " " + " ".join(s["char"]) for s in samples)
        series_blob = series["title"] + " " + series["summary"]
        has_cldn4_pert_sample = bool(
            CLDN4_RE.search(sample_blob) and PERTURB_RE.search(sample_blob))
        out[acc] = {"series": series, "samples": samples,
                    "n_samples": len(samples),
                    "has_cldn4_perturb_sample": has_cldn4_pert_sample}
        organisms = sorted({s["organism"] for s in samples if s["organism"]})
        print(f"\n===== {acc} ({len(samples)} samples) {organisms} =====")
        print("TYPE:", series["type"])
        print("TITLE:", series["title"])
        print("SUPP :", series["supp"])
        for s in samples:
            print(f"  {s['gsm']:12} | {s['title']:38} | {s['source'][:30]:30} | {'; '.join(s['char'])[:80]}")

    with open("../../results/fable_cldn4_kdko/geo_verified.json", "w") as f:
        json.dump(out, f, indent=1)
    print("\nSaved geo_verified.json")


if __name__ == "__main__":
    main()
