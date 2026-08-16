"""Build the final curated master table joining:
  * authoritative verification (verified_series.tsv / .json)
  * scope classification (human AND lung AND ICI)
  * clinical label availability (clinical_label_index.tsv)
  * target-gene (TACSTD2/CLDN4/CD274) measurability (gene_availability_audit.tsv)
  * download status (download_manifest.tsv)

Output: results/fable_geo_2015_2018/master_summary.tsv
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "fable_geo_2015_2018")


def load_tsv(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in f]


def main():
    verified = json.load(open(os.path.join(OUT, "verified_series.json")))
    labels = {r["accession"]: r for r in load_tsv(os.path.join(OUT, "clinical_label_index.tsv"))}
    audit = load_tsv(os.path.join(OUT, "analysis", "gene_availability_audit.tsv"))
    gene_by_acc = {}
    for r in audit:
        acc = r["accession"]
        g = gene_by_acc.setdefault(acc, {"TACSTD2": False, "CLDN4": False, "CD274": False})
        for k in ("TACSTD2", "CLDN4", "CD274"):
            g[k] = g[k] or (r[k] == "True")
    manifest = load_tsv(os.path.join(OUT, "download_manifest.tsv"))
    dl_by_acc = {}
    for r in manifest:
        acc = r["accession"]
        d = dl_by_acc.setdefault(acc, {"n": 0, "bytes": 0, "skipped": 0})
        if r["status"] in ("downloaded", "already_present"):
            d["n"] += 1
            d["bytes"] += int(r["stored_bytes"] or 0)
        elif r["status"] == "skipped_ge_2gb":
            d["skipped"] += 1

    rows = []
    for v in verified:
        acc = v["accession"]
        in_scope = v["is_human"] and v["is_lung"] and v["is_ici"]
        lab = labels.get(acc, {})
        gene = gene_by_acc.get(acc, {})
        dl = dl_by_acc.get(acc, {})
        rows.append({
            "accession": acc,
            "organism": v["taxon"],
            "in_scope_human_lung_ici": in_scope,
            "is_human": v["is_human"], "is_lung": v["is_lung"], "is_ici": v["is_ici"],
            "public_status": v["status"],
            "submission_date": v["submission_date"],
            "series_type": v["types"],
            "n_samples": v["n_samples"],
            "pubmed": v["pubmed"],
            "has_response_or_survival_labels": lab.get("has_response_labels", ""),
            "label_keys": lab.get("response_like_keys", ""),
            "TACSTD2_measured": gene.get("TACSTD2", ""),
            "CLDN4_measured": gene.get("CLDN4", ""),
            "CD274_measured": gene.get("CD274", ""),
            "downloaded_files": dl.get("n", 0),
            "downloaded_MB": round(dl.get("bytes", 0) / 1024**2, 1),
            "files_skipped_ge_2gb": dl.get("skipped", 0),
            "title": v["title"],
        })
    rows.sort(key=lambda r: (not r["in_scope_human_lung_ici"], r["accession"]))
    cols = ["accession", "organism", "in_scope_human_lung_ici", "is_human",
            "is_lung", "is_ici", "public_status", "submission_date", "series_type",
            "n_samples", "pubmed", "has_response_or_survival_labels", "label_keys",
            "TACSTD2_measured", "CLDN4_measured", "CD274_measured",
            "downloaded_files", "downloaded_MB", "files_skipped_ge_2gb", "title"]
    with open(os.path.join(OUT, "master_summary.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    n_scope = sum(r["in_scope_human_lung_ici"] for r in rows)
    print(f"master_summary.tsv: {len(rows)} verified series, {n_scope} in-scope "
          f"(human+lung+ICI)")


if __name__ == "__main__":
    main()
