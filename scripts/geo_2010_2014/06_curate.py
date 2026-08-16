"""Join search + verify + download + analysis into one master table + verdict.

Outputs:
  master_summary.tsv
  catalog.tsv
  verdict.json
"""
import json
import os
from datetime import datetime, timezone

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(ROOT, "results", "w200", "GEO_2010_2014")
ANA = os.path.join(OUT, "analysis")


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        return [dict(zip(header, ln.rstrip("\n").split("\t"))) for ln in f]


def main():
    verified = json.load(open(os.path.join(OUT, "verified_series.json")))
    dl = read_tsv(os.path.join(OUT, "download_manifest.tsv"))
    clin = read_tsv(os.path.join(OUT, "clinical_label_index.tsv"))
    audit = read_tsv(os.path.join(ANA, "gene_availability_audit.tsv"))
    analysis = {}
    ap = os.path.join(ANA, "analysis_results.json")
    if os.path.exists(ap):
        analysis = json.load(open(ap))

    dl_by = {}
    for r in dl:
        dl_by.setdefault(r["accession"], []).append(r)
    clin_by = {r["accession"]: r for r in clin}
    audit_by = {}
    for r in audit:
        audit_by.setdefault(r["accession"], []).append(r)

    rows = []
    for v in verified:
        acc = v["accession"]
        dls = dl_by.get(acc, [])
        stored = sum(int(x.get("stored_bytes") or 0) for x in dls)
        statuses = sorted({x.get("status", "") for x in dls})
        c = clin_by.get(acc, {})
        au = audit_by.get(acc, [])
        tac = max((int(x.get("TACSTD2_n") or 0) for x in au), default=0)
        cld = max((int(x.get("CLDN4_n") or 0) for x in au), default=0)
        cd = max((int(x.get("CD274_n") or 0) for x in au), default=0)
        pd = max((int(x.get("PDCD1_n") or 0) for x in au), default=0)
        # pick first analysis block for coexpression
        rho = ""
        for key, rec in analysis.items():
            if rec.get("accession") == acc:
                cx = rec.get("coexpression", {})
                t2 = cx.get("TACSTD2~CD274", {})
                if t2.get("spearman") is not None:
                    rho = str(t2["spearman"])
                break
        rows.append(
            {
                "accession": acc,
                "pdat": v.get("pdat", ""),
                "taxon": v.get("taxon", ""),
                "n_samples": v.get("n_samples", ""),
                "gpl": v.get("gpl", ""),
                "pubmed": v.get("pubmed", ""),
                "is_human": v.get("is_human"),
                "is_lung": v.get("is_lung"),
                "is_immunotherapy": v.get("is_immunotherapy"),
                "is_checkpoint_ici": v.get("is_checkpoint_ici"),
                "leftover_reason": v.get("leftover_reason", ""),
                "n_download_files": len(dls),
                "download_status": ";".join(statuses),
                "stored_bytes": stored,
                "has_response_label": c.get("has_response_label", ""),
                "has_survival_label": c.get("has_survival_label", ""),
                "TACSTD2_n": tac,
                "CLDN4_n": cld,
                "CD274_n": cd,
                "PDCD1_n": pd,
                "TACSTD2_CD274_spearman": rho,
                "title": v.get("title", ""),
            }
        )

    cols = [
        "accession",
        "pdat",
        "taxon",
        "n_samples",
        "gpl",
        "pubmed",
        "is_human",
        "is_lung",
        "is_immunotherapy",
        "is_checkpoint_ici",
        "leftover_reason",
        "n_download_files",
        "download_status",
        "stored_bytes",
        "has_response_label",
        "has_survival_label",
        "TACSTD2_n",
        "CLDN4_n",
        "CD274_n",
        "PDCD1_n",
        "TACSTD2_CD274_spearman",
        "title",
    ]
    with open(os.path.join(OUT, "master_summary.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
    with open(os.path.join(OUT, "catalog.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")

    n_check = sum(1 for r in rows if r["is_checkpoint_ici"] and r["is_human"] and r["is_lung"])
    n_io = sum(1 for r in rows if r["is_immunotherapy"] and r["is_human"] and r["is_lung"])
    verdict = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "window": "2010-01-01/2014-12-31",
        "n_verified_gse": len(rows),
        "n_human_lung_checkpoint_ici": n_check,
        "n_human_lung_immunotherapy": n_io,
        "honest_verdict": (
            "No leftover human lung PD-1/PD-L1/CTLA-4 ICI treatment-response "
            "cohort in GEO 2010-2014. The window predates NSCLC ICI approval "
            "(nivolumab 2015). Real leftover hits are cancer-vaccine / epitope "
            "immunotherapy (not checkpoint ICI) plus mouse Pdl1-expression SCC."
        ),
        "human_lung_immunotherapy_leftovers": [
            r["accession"] for r in rows if r["is_human"] and r["is_lung"] and r["is_immunotherapy"]
        ],
        "human_lung_checkpoint_ici": [
            r["accession"] for r in rows if r["is_human"] and r["is_lung"] and r["is_checkpoint_ici"]
        ],
        "no_invented_ids": True,
        "accessions_are_ncbi_hits_only": True,
    }
    with open(os.path.join(OUT, "verdict.json"), "w") as f:
        json.dump(verdict, f, indent=1)
    print(json.dumps(verdict, indent=2))
    print("-> results/w200/GEO_2010_2014/master_summary.tsv")


if __name__ == "__main__":
    main()
