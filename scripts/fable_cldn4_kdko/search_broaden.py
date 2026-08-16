"""Second-pass exhaustive search: SRA + broad GEO, to catch CLDN4
perturbation series where the KD/KO signal only appears at sample level.
"""
import json
from collections import OrderedDict
from eutils import esearch, esummary

# ---- GEO gds, looser (no Entry Type filter, catch DataSets/subseries) ----
gds_terms = [
    '(CLDN4[Title]) OR ("claudin-4"[Title]) OR ("claudin 4"[Title])',
    'CLDN4 AND (siRNA OR shRNA OR CRISPR OR knockdown OR knockout)',
    '("claudin 4" OR "claudin-4" OR CLDN4) AND (lung OR pulmonary OR alveolar)',
]
print("==== GEO gds broad ====")
gds_ids = OrderedDict()
for t in gds_terms:
    r = esearch("gds", t, retmax=1000)
    print(f"[{r.get('count')}] {t}")
    for i in r.get("idlist", []):
        gds_ids.setdefault(i, t)
summ = esummary("gds", list(gds_ids.keys()))
gds_rows = []
for uid, m in summ.items():
    if m.get("entryType") != "GSE":
        continue
    gds_rows.append({"accession": m.get("accession"), "taxon": m.get("taxon"),
                     "n": m.get("n_samples"), "type": m.get("gdsType"),
                     "title": m.get("title")})
gds_rows.sort(key=lambda x: x["accession"])
print(f"\nUnique GSE from broad GEO: {len(gds_rows)}")
for r in gds_rows:
    print(f"  {r['accession']:12} {str(r['taxon'])[:20]:20} n={str(r['n']):3} {r['title'][:75]}")

# ---- SRA: sample-level CLDN4 perturbation ----
print("\n==== SRA ====")
sra_terms = [
    'CLDN4 AND (knockdown OR knockout OR shRNA OR siRNA OR CRISPR)',
    '"claudin 4" AND (knockdown OR knockout OR shRNA OR siRNA OR CRISPR)',
]
sra_ids = OrderedDict()
for t in sra_terms:
    r = esearch("sra", t, retmax=300)
    print(f"[{r.get('count')}] {t}")
    for i in r.get("idlist", []):
        sra_ids.setdefault(i, t)
ssum = esummary("sra", list(sra_ids.keys())[:200])
# Extract study accessions (SRP/PRJNA) + titles from expxml blob
import re
studies = {}
for uid, m in ssum.items():
    exp = m.get("expxml", "")
    title = ""
    mt = re.search(r"<Title>(.*?)</Title>", exp)
    if mt:
        title = mt.group(1)
    srp = re.findall(r"(SRP\d+|PRJNA\d+|GSE\d+)", exp)
    key = srp[0] if srp else uid
    studies.setdefault(key, {"title": title, "n": 0, "any_srp": set()})
    studies[key]["n"] += 1
    for s in srp:
        studies[key]["any_srp"].add(s)
print(f"\nSRA studies (grouped): {len(studies)}")
for k, v in sorted(studies.items()):
    print(f"  {k:14} runs~{v['n']:3} {sorted(v['any_srp'])} :: {v['title'][:70]}")

out = {"gds_broad": gds_rows, "sra_studies": {k: {"title": v["title"], "n": v["n"], "acc": sorted(v["any_srp"])} for k, v in studies.items()}}
with open("../../results/fable_cldn4_kdko/search_broaden.json", "w") as f:
    json.dump(out, f, indent=1)
print("\nsaved search_broaden.json")
