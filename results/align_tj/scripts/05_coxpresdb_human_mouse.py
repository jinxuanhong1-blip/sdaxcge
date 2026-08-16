"""COXPRESdb v8 public coexpression (human Hsa-u.c4-0 and mouse Mmu-u.c4-0).

For each query gene we pull the top 1000 coexpressed genes (mutual rank + logit
score) and record where TACSTD2/CLDN4 and the user TJ/TF lists sit.

This is corpus-wide (all tissues), not lung-restricted — stated as such.

Outputs:
  tables/coxpresdb_query_top1000.tsv     (all returned neighbors)
  tables/coxpresdb_target_hits.tsv       (query x partner rank/MR if in top 1000)
"""
import json
import time
import urllib.request
import pandas as pd
import common as C

# Entrez Gene IDs
HUMAN = {
    "TACSTD2": 4070, "CLDN1": 9076, "CLDN4": 1364, "CLDN7": 1366,
    "F11R": 50848, "PARD3": 56288,
    "ELF3": 1999, "GRHL1": 29841, "KLF4": 9314, "TFAP2A": 7020,
}
MOUSE = {
    "Tacstd2": 56753, "Cldn1": 12737, "Cldn4": 12740, "Cldn7": 53624,
    "F11r": 16456, "Pard3": 93742,
    "Elf3": 13710, "Grhl1": 75077, "Klf4": 16600, "Tfap2a": 21418,
}

QUERIES_H = ["TACSTD2", "CLDN4", "ELF3", "GRHL1", "KLF4", "TFAP2A"]
QUERIES_M = ["Tacstd2", "Cldn4", "Elf3", "Grhl1", "Klf4", "Tfap2a"]
PARTNERS_H = list(HUMAN)
PARTNERS_M = list(MOUSE)


def fetch(gene_id, db, topn=1000, retries=4):
    url = f"https://coxpresdb.jp/api4/?gene={gene_id}&db={db}&topN={topn}"
    last = None
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "align-tj/1.0"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            last = e
            time.sleep(2 ** i)
    raise RuntimeError(f"COXPRESdb failed {url}: {last}")


def harvest(symbol_to_id, queries, partners, db, species):
    id_to_sym = {v: k for k, v in symbol_to_id.items()}
    partner_ids = {symbol_to_id[s] for s in partners}
    long_rows, hit_rows = [], []
    for q in queries:
        qid = symbol_to_id[q]
        print(f"  {species} {q} ({qid}) ...")
        js = fetch(qid, db)
        recs = js["result_set"][0]["results"]
        for rank, rec in enumerate(recs, start=1):
            gid = int(rec["gene"])
            long_rows.append({
                "species": species, "db": db, "query": q, "query_id": qid,
                "rank": rank, "neighbor_id": gid,
                "neighbor": id_to_sym.get(gid, ""),
                "mutual_rank": rec.get("mutual_rank"),
                "logit_score": rec.get("logit_score"),
            })
            if gid in partner_ids:
                hit_rows.append({
                    "species": species, "db": db, "query": q,
                    "partner": id_to_sym[gid],
                    "rank_in_top1000": rank,
                    "mutual_rank": rec.get("mutual_rank"),
                    "logit_score": rec.get("logit_score"),
                    "in_top1000": True,
                })
        found = {r["partner"] for r in hit_rows if r["query"] == q}
        for p in partners:
            if p == q:
                continue
            if p not in found:
                hit_rows.append({
                    "species": species, "db": db, "query": q, "partner": p,
                    "rank_in_top1000": None, "mutual_rank": None,
                    "logit_score": None, "in_top1000": False,
                })
    return long_rows, hit_rows


def main():
    print("COXPRESdb human...")
    lh, hh = harvest(HUMAN, QUERIES_H, PARTNERS_H, "hsa-u", "human")
    print("COXPRESdb mouse...")
    lm, hm = harvest(MOUSE, QUERIES_M, PARTNERS_M, "mmu-u", "mouse")
    long = pd.DataFrame(lh + lm)
    hits = pd.DataFrame(hh + hm)
    long.to_csv(f"{C.TABLES}/coxpresdb_query_top1000.tsv", sep="\t", index=False)
    hits.to_csv(f"{C.TABLES}/coxpresdb_target_hits.tsv", sep="\t", index=False)
    # compact view of the pairs the user asked about
    focus = hits[hits.partner.str.upper().isin(
        ["TACSTD2", "CLDN4"] + [x.upper() for x in C.USER_TFS + C.USER_TJ_GENES]
    )].copy()
    print(focus.sort_values(["species", "query", "rank_in_top1000"],
                            na_position="last").to_string(index=False))


if __name__ == "__main__":
    main()
