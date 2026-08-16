"""Parse per-sample clinical metadata from the downloaded series_matrix files for
every in-scope (human + lung + ICI) series, and flag which series actually carry
treatment / response / survival labels.

Outputs:
  results/fable_geo_2015_2018/clinical/<ACC>_clinical.tsv     (tidy per-sample)
  results/fable_geo_2015_2018/clinical_label_index.tsv        (label availability)
"""
import glob
import gzip
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "fable_geo_2015_2018")
DATA = os.path.join(OUT, "data")
CLIN = os.path.join(OUT, "clinical")
os.makedirs(CLIN, exist_ok=True)

RESP_KEY = re.compile(r"resp|recist|benefit|outcome|nivol|pembro|atezo|durva|"
                      r"ipili|treatment|therapy|drug|pfs|surv( |_)|vital|"
                      r"os_|overall survival|progression", re.I)


def read_matrix_meta(path):
    op = gzip.open if path.endswith(".gz") else open
    samples = []
    char_rows = []
    title = src = None
    with op(path, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!Sample_geo_accession"):
                samples = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_title"):
                title = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_source_name_ch1"):
                src = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                char_rows.append(vals)
            elif line.startswith("!series_matrix_table_begin"):
                break
    return samples, title, src, char_rows


def main():
    accs = sorted({os.path.basename(os.path.dirname(p))
                   for p in glob.glob(os.path.join(DATA, "*", "*series_matrix.txt.gz"))})
    index = []
    for acc in accs:
        # merge all platform matrices for this series
        rows_by_gsm = {}
        keys = set()
        for path in sorted(glob.glob(os.path.join(DATA, acc, "*series_matrix.txt.gz"))):
            samples, title, src, char_rows = read_matrix_meta(path)
            if not samples:
                continue
            for i, gsm in enumerate(samples):
                rec = rows_by_gsm.setdefault(gsm, {"gsm": gsm})
                if title and i < len(title):
                    rec["title"] = title[i]
                if src and i < len(src):
                    rec["source"] = src[i]
                for row in char_rows:
                    if i < len(row) and row[i]:
                        cell = row[i]
                        if ":" in cell:
                            k, v = cell.split(":", 1)
                            k = k.strip().lower().replace(" ", "_")
                            rec[k] = v.strip()
                            keys.add(k)
        if not rows_by_gsm:
            continue
        keys = ["title", "source"] + sorted(keys)
        cols = ["gsm"] + keys
        outp = os.path.join(CLIN, f"{acc}_clinical.tsv")
        with open(outp, "w") as f:
            f.write("\t".join(cols) + "\n")
            for gsm, rec in sorted(rows_by_gsm.items()):
                f.write("\t".join(str(rec.get(c, "")) for c in cols) + "\n")
        resp_keys = [k for k in keys if RESP_KEY.search(k)]
        index.append({
            "accession": acc,
            "n_samples": len(rows_by_gsm),
            "n_char_keys": len(keys),
            "response_like_keys": ";".join(resp_keys) if resp_keys else "",
            "has_response_labels": bool(resp_keys),
        })
        print(f"{acc}: {len(rows_by_gsm)} samples, response-like keys: {resp_keys}")

    icols = ["accession", "n_samples", "n_char_keys", "has_response_labels",
             "response_like_keys"]
    with open(os.path.join(OUT, "clinical_label_index.tsv"), "w") as f:
        f.write("\t".join(icols) + "\n")
        for r in index:
            f.write("\t".join(str(r[c]) for c in icols) + "\n")
    print(f"\nwrote clinical tables for {len(index)} series -> {CLIN}")


if __name__ == "__main__":
    main()
