"""Parse per-sample characteristics from downloaded series matrices.

Flags likely response / survival / treatment labels. Does not invent labels.
Outputs:
  clinical/<ACC>_clinical.tsv
  clinical_label_index.tsv
"""
import gzip
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "w200", "GEO_2010_2014")
DATA = os.path.join(OUT, "data")
CLIN = os.path.join(OUT, "clinical")
os.makedirs(CLIN, exist_ok=True)

CHAR_RE = re.compile(r"^!Sample_characteristics_ch1\t(.+)$")
TITLE_RE = re.compile(r"^!Sample_title\t(.+)$")
ACC_RE = re.compile(r"^!Sample_geo_accession\t(.+)$")
SRC_RE = re.compile(r"^!Sample_source_name_ch1\t(.+)$")


def _unq(cell):
    return cell.strip().strip('"')


def parse_matrix(path):
    opener = gzip.open if path.endswith(".gz") else open
    samples = []
    titles = []
    sources = []
    chars = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            if ln.startswith("!series_matrix_table_begin"):
                break
            m = ACC_RE.match(ln.rstrip("\n"))
            if m:
                samples = [_unq(c) for c in m.group(1).split("\t")]
                continue
            m = TITLE_RE.match(ln.rstrip("\n"))
            if m:
                titles = [_unq(c) for c in m.group(1).split("\t")]
                continue
            m = SRC_RE.match(ln.rstrip("\n"))
            if m:
                sources = [_unq(c) for c in m.group(1).split("\t")]
                continue
            m = CHAR_RE.match(ln.rstrip("\n"))
            if m:
                chars.append([_unq(c) for c in m.group(1).split("\t")])
    n = len(samples)
    rows = []
    for i in range(n):
        rec = {
            "geo_accession": samples[i] if i < n else "",
            "title": titles[i] if i < len(titles) else "",
            "source_name": sources[i] if i < len(sources) else "",
        }
        extra = []
        for ch in chars:
            val = ch[i] if i < len(ch) else ""
            if ":" in val:
                k, v = val.split(":", 1)
                key = re.sub(r"[^A-Za-z0-9_]+", "_", k.strip().lower()).strip("_")
                if key and key not in rec:
                    rec[key] = v.strip()
                else:
                    extra.append(val)
            elif val:
                extra.append(val)
        rec["other_characteristics"] = " | ".join(extra)
        rows.append(rec)
    return rows


def flag_labels(rows):
    keys = set()
    for r in rows:
        keys.update(r.keys())
    blob = " ".join(keys).lower()
    has_response = bool(
        re.search(r"resp|recist|cr_pr|responder|clinical_outcome|mage", blob)
    )
    has_survival = bool(re.search(r"surviv|vital|os_|pfs|dead|censor", blob))
    has_treatment = bool(re.search(r"treat|drug|therapy|vaccine|immun", blob))
    return has_response, has_survival, has_treatment, sorted(keys)


def main():
    index = []
    if not os.path.isdir(DATA):
        print("no data/ yet")
        return
    for acc in sorted(os.listdir(DATA)):
        d = os.path.join(DATA, acc)
        if not os.path.isdir(d):
            continue
        matrices = [
            os.path.join(d, n)
            for n in os.listdir(d)
            if "series_matrix" in n.lower()
        ]
        if not matrices:
            index.append(
                {
                    "accession": acc,
                    "n_samples": 0,
                    "has_response_label": False,
                    "has_survival_label": False,
                    "has_treatment_label": False,
                    "n_char_keys": 0,
                    "char_keys": "",
                    "note": "no_series_matrix",
                }
            )
            continue
        # Prefer the first matrix; multi-platform series get one clinical table
        # with a platform column if present in the filename.
        all_rows = []
        for path in sorted(matrices):
            rows = parse_matrix(path)
            plat = ""
            m = re.search(r"-(GPL\d+)_series_matrix", os.path.basename(path))
            if m:
                plat = m.group(1)
            for r in rows:
                r["matrix_file"] = os.path.basename(path)
                r["platform_from_filename"] = plat
            all_rows.extend(rows)
        if not all_rows:
            continue
        keys = []
        for r in all_rows:
            for k in r:
                if k not in keys:
                    keys.append(k)
        outp = os.path.join(CLIN, f"{acc}_clinical.tsv")
        with open(outp, "w") as f:
            f.write("\t".join(keys) + "\n")
            for r in all_rows:
                f.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")
        has_r, has_s, has_t, _ = flag_labels(all_rows)
        index.append(
            {
                "accession": acc,
                "n_samples": len(all_rows),
                "has_response_label": has_r,
                "has_survival_label": has_s,
                "has_treatment_label": has_t,
                "n_char_keys": len(keys),
                "char_keys": ";".join(keys),
                "note": "",
            }
        )
        print(f"{acc}: n={len(all_rows)} response={has_r} survival={has_s} treat={has_t}")

    cols = [
        "accession",
        "n_samples",
        "has_response_label",
        "has_survival_label",
        "has_treatment_label",
        "n_char_keys",
        "char_keys",
        "note",
    ]
    with open(os.path.join(OUT, "clinical_label_index.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in index:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    print("-> results/w200/GEO_2010_2014/clinical_label_index.tsv")


if __name__ == "__main__":
    main()
