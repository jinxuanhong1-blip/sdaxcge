#!/usr/bin/env python3
"""Confirm TACSTD2/CLDN4 presence and summarize lung ICI proteomes (processed tables only)."""

from __future__ import annotations

import csv
import gzip
import json
import math
import re
import statistics
from collections import defaultdict
from pathlib import Path

ROOT = Path("results/grok_pride")
PROC = ROOT / "processed"
OUT = ROOT
OUT.mkdir(parents=True, exist_ok=True)

# Official UniProt + aliases
HUMAN = {
    "TACSTD2": {
        "uniprot": {"P09758"},
        "genes": {"TACSTD2", "TACD2", "TROP2", "EGP1", "GA7331", "M1S1", "EGP-1", "GA733-1"},
        "names": {"tumor-associated calcium signal transducer 2", "tumor associated calcium signal transducer 2"},
    },
    "CLDN4": {
        "uniprot": {"O14493"},
        "genes": {"CLDN4", "CLD4", "CPE-R", "CPER", "WBSCR8", "CPETR1"},
        "names": {"claudin-4", "claudin 4"},
    },
}
MOUSE = {
    "Tacstd2": {
        "uniprot": {"Q8BGV3"},
        "genes": {"TACSTD2", "TACD2", "TROP2", "EGP1"},
        "names": {"tumor-associated calcium signal transducer 2"},
    },
    "Cldn4": {
        "uniprot": {"O35054"},
        "genes": {"CLDN4", "CLD4", "CPE-R"},
        "names": {"claudin-4", "claudin 4"},
    },
}

TOKEN_RE = re.compile(r"[A-Z0-9][A-Z0-9\-_]+", re.I)


def hit_human(text: str) -> list[str]:
    t = text.upper()
    found = []
    for gene, spec in HUMAN.items():
        if any(u in t for u in spec["uniprot"]):
            found.append(gene)
            continue
        if any(re.search(rf"(?<![A-Z0-9]){re.escape(g)}(?![A-Z0-9])", t) for g in spec["genes"]):
            found.append(gene)
    return found


def hit_mouse(text: str) -> list[str]:
    t = text.upper()
    found = []
    for gene, spec in MOUSE.items():
        if any(u in t for u in spec["uniprot"]):
            found.append(gene)
            continue
        if any(re.search(rf"(?<![A-Z0-9]){re.escape(g)}(?![A-Z0-9])", t) for g in spec["genes"]):
            found.append(gene)
    return found


def parse_pxd042091_all() -> list[dict]:
    p = PROC / "PXD042091" / "2_dat_cs_all.txt"
    rows = []
    with p.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        # first col empty index name; last is accessions
        samples = [h for h in header if h and h != "accessions"]
        acc_idx = header.index("accessions") if "accessions" in header else -1
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc = parts[acc_idx] if acc_idx >= 0 and acc_idx < len(parts) else parts[-1]
            vals = []
            for i, h in enumerate(header):
                if not h or h == "accessions":
                    continue
                if i < len(parts) and parts[i] not in ("", "NA", "NaN"):
                    try:
                        vals.append(float(parts[i]))
                    except ValueError:
                        pass
            rows.append({"accession": acc, "n_quant": len(vals), "median": statistics.median(vals) if vals else None, "mean": statistics.mean(vals) if vals else None})
    return rows


def parse_pxd042091_tf() -> list[dict]:
    p = PROC / "PXD042091" / "6_data_tf_response.txt"
    rows = []
    with p.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        # protein, R*, NR*, p_val, foldchange, log_pval
        r_cols = [i for i, h in enumerate(header) if h.startswith("R") and not h.startswith("NR")]
        nr_cols = [i for i, h in enumerate(header) if h.startswith("NR")]
        p_idx = header.index("p_val") if "p_val" in header else None
        fc_idx = header.index("foldchange") if "foldchange" in header else None
        prot_idx = header.index("protein") if "protein" in header else 1
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= prot_idx:
                continue
            prot = parts[prot_idx]
            def nums(idxs):
                out = []
                for i in idxs:
                    if i < len(parts) and parts[i] not in ("", "NA", "NaN"):
                        try:
                            out.append(float(parts[i]))
                        except ValueError:
                            pass
                return out
            r = nums(r_cols)
            nr = nums(nr_cols)
            rec = {
                "protein": prot,
                "n_R": len(r),
                "n_NR": len(nr),
                "median_R": statistics.median(r) if r else None,
                "median_NR": statistics.median(nr) if nr else None,
                "p_val": float(parts[p_idx]) if p_idx is not None and p_idx < len(parts) and parts[p_idx] not in ("", "NA") else None,
                "foldchange": float(parts[fc_idx]) if fc_idx is not None and fc_idx < len(parts) and parts[fc_idx] not in ("", "NA") else None,
            }
            rows.append(rec)
    return rows


def parse_mztab_proteins(path: Path) -> list[dict]:
    opener = gzip.open if path.suffix == ".gz" or path.name.endswith(".gz") else open
    proteins = []
    with opener(path, "rt", errors="replace") as f:
        headers = None
        for line in f:
            if line.startswith("PRH"):
                headers = line.rstrip("\n").split("\t")
                continue
            if not line.startswith("PRT") or headers is None:
                continue
            parts = line.rstrip("\n").split("\t")
            rec = {headers[i]: parts[i] if i < len(parts) else "" for i in range(len(headers))}
            proteins.append(rec)
    return proteins


def main() -> None:
    summary = {
        "datasets": {},
        "target_hits": [],
        "notes": [],
    }

    # ---- PXD042091 ----
    all_rows = parse_pxd042091_all()
    tf_rows = parse_pxd042091_tf()
    print(f"PXD042091 all proteins: {len(all_rows)}")
    print(f"PXD042091 TF-response proteins: {len(tf_rows)}")

    hits_all = [r for r in all_rows if hit_human(r["accession"])]
    hits_tf = [r for r in tf_rows if hit_human(r["protein"])]
    print("PXD042091 target hits in all:", hits_all)
    print("PXD042091 target hits in TF:", hits_tf)

    # also dump a few accession examples and gene-token inventory
    gene_tokens = set()
    for r in all_rows:
        acc = r["accession"]
        if "_" in acc:
            gene_tokens.add(acc.split("_", 1)[1])
    print("sample gene tokens containing CLD/TAC/TROP:", sorted(g for g in gene_tokens if any(x in g.upper() for x in ("CLD", "TAC", "TROP", "EGP"))))

    # write protein list
    (OUT / "PXD042091_protein_accessions.tsv").write_text(
        "accession\tn_quant\tmedian\tmean\thits\n"
        + "\n".join(
            f"{r['accession']}\t{r['n_quant']}\t{r['median']}\t{r['mean']}\t{','.join(hit_human(r['accession']))}"
            for r in all_rows
        )
        + "\n"
    )
    (OUT / "PXD042091_tf_response_targets.tsv").write_text(
        "protein\tn_R\tn_NR\tmedian_R\tmedian_NR\tp_val\tfoldchange\thits\n"
        + "\n".join(
            f"{r['protein']}\t{r['n_R']}\t{r['n_NR']}\t{r['median_R']}\t{r['median_NR']}\t{r['p_val']}\t{r['foldchange']}\t{','.join(hit_human(r['protein']))}"
            for r in tf_rows
            if hit_human(r["protein"]) or True  # keep full table? too big. write only if we want full
        )
        + "\n"
    )
    # rewrite TF as full but that's 6k lines - OK
    with (OUT / "PXD042091_tf_stats.tsv").open("w") as f:
        f.write("protein\tn_R\tn_NR\tmedian_R\tmedian_NR\tp_val\tfoldchange\thits\n")
        for r in tf_rows:
            hits = ",".join(hit_human(r["protein"]))
            f.write(f"{r['protein']}\t{r['n_R']}\t{r['n_NR']}\t{r['median_R']}\t{r['median_NR']}\t{r['p_val']}\t{r['foldchange']}\t{hits}\n")

    # signature proteins from paper
    sig = ["ATG9A", "DCDC2", "HPS5", "FIL1L", "LZTL1", "PGTA", "SPTN2", "ATG9", "DCDC", "HPS5", "FLII", "LZTFL1", "RABGGTA", "SPTBN2"]
    sig_hits = [r for r in tf_rows if any(s in r["protein"].upper() for s in sig)]
    print("signature-like proteins:", [r["protein"] for r in sig_hits[:30]])

    # related epithelial / tight-junction proteins
    related = []
    keys = ("CLD", "OCLN", "TJP", "EPCAM", "CDH1", "KRT", "MUC1", "CEACAM", "CD274", "PDCD1", "PDCD1LG2")
    for r in tf_rows:
        prot = r["protein"].upper()
        if any(k in prot for k in keys):
            related.append(r)
    print("related epithelial/IO proteins n=", len(related))
    for r in related:
        print(" ", r["protein"], "p", r["p_val"], "fc", r["foldchange"], "medR", r["median_R"], "medNR", r["median_NR"])

    summary["datasets"]["PXD042091"] = {
        "n_proteins_all": len(all_rows),
        "n_proteins_tf": len(tf_rows),
        "tacstd2_present": any("TACSTD2" in hit_human(r["accession"]) for r in all_rows) or any("TACSTD2" in hit_human(r["protein"]) for r in tf_rows),
        "cldn4_present": any("CLDN4" in hit_human(r["accession"]) for r in all_rows) or any("CLDN4" in hit_human(r["protein"]) for r in tf_rows),
        "hits_all": hits_all,
        "hits_tf": hits_tf,
        "related_n": len(related),
        "related": [
            {"protein": r["protein"], "p_val": r["p_val"], "foldchange": r["foldchange"], "median_R": r["median_R"], "median_NR": r["median_NR"]}
            for r in related
        ],
    }

    # ---- PXD059688 mzTab ----
    mz = PROC / "PXD059688" / "20230904_Tumor_AA.mzTab.gz"
    proteins = parse_mztab_proteins(mz)
    print(f"PXD059688 PRT rows: {len(proteins)}")
    if proteins:
        print("PRT keys", list(proteins[0].keys())[:25])
        print("example", {k: proteins[0][k] for k in list(proteins[0])[:12]})

    mz_hits = []
    for rec in proteins:
        blob = " ".join(str(v) for v in rec.values())
        hs = hit_mouse(blob) or hit_human(blob)
        if hs:
            mz_hits.append({"hits": hs, "accession": rec.get("accession", ""), "description": rec.get("description", "")[:200], "rec": {k: rec[k] for k in rec if "abundance" in k.lower() or k in ("accession", "description", "taxid", "species")}})
    print("PXD059688 target hits:", len(mz_hits))
    for h in mz_hits:
        print(" ", h["hits"], h["accession"], h["description"])

    # write all protein accessions + descriptions
    with (OUT / "PXD059688_proteins.tsv").open("w") as f:
        f.write("accession\tdescription\thits\n")
        for rec in proteins:
            blob = " ".join(str(v) for v in rec.values())
            hs = ",".join(hit_mouse(blob) or hit_human(blob))
            f.write(f"{rec.get('accession','')}\t{rec.get('description','').replace(chr(9),' ')}\t{hs}\n")

    # abundance columns for targets + related
    related_mz = []
    for rec in proteins:
        blob = " ".join(str(v) for v in rec.values()).upper()
        if any(k in blob for k in ("TACSTD", "TACD2", "TROP2", "CLDN4", "CLD4", "CLAUDIN-4", "CLAUDIN 4", "EPCAM", "CDH1", "OCLN", "TJP1", "CD274", "PDCD1")):
            related_mz.append(rec)
    print("PXD059688 related n", len(related_mz))

    summary["datasets"]["PXD059688"] = {
        "n_prt": len(proteins),
        "tacstd2_present": any("Tacstd2" in h["hits"] or "TACSTD2" in h["hits"] for h in mz_hits),
        "cldn4_present": any("Cldn4" in h["hits"] or "CLDN4" in h["hits"] for h in mz_hits),
        "hits": [
            {"hits": h["hits"], "accession": h["accession"], "description": h["description"]}
            for h in mz_hits
        ],
    }

    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2, default=str))
    print("wrote analysis_summary.json")


if __name__ == "__main__":
    main()
