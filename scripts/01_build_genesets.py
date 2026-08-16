#!/usr/bin/env python3
"""Build mouse keratin / tight-junction gene sets (GMT) from GO + MGI annotations.

Sources:
- go-basic.obo (GO release, downloaded from purl.obolibrary.org)
- mgi.gaf.gz   (mouse GO annotations, current.geneontology.org)

GO annotations are propagated: a gene annotated to any descendant (via is_a /
part_of) of a target term is included in that target term's set.
Two hand-curated sets (keratin family by symbol pattern; core tight-junction
components) are added because GO annotation of mouse keratins is incomplete.
"""
import gzip
import re
from collections import defaultdict
from pathlib import Path

ANNOT_DIR = Path("/workspace/data/annotations")
OUT_DIR = Path("/workspace/results/w200/A8_mouse/genesets")
OUT_DIR.mkdir(parents=True, exist_ok=True)

TARGET_TERMS = {
    "GO:0045095": "GO_KERATIN_FILAMENT",
    "GO:0031424": "GO_KERATINIZATION",
    "GO:0070268": "GO_CORNIFICATION",
    "GO:0005882": "GO_INTERMEDIATE_FILAMENT",
    "GO:0005923": "GO_BICELLULAR_TIGHT_JUNCTION",
    "GO:0070160": "GO_TIGHT_JUNCTION",  # occluding junction
    "GO:0120193": "GO_TIGHT_JUNCTION_ORGANIZATION",
    "GO:0120192": "GO_TIGHT_JUNCTION_ASSEMBLY",
    "GO:0043296": "GO_APICAL_JUNCTION_COMPLEX",
    # context / specificity controls (other junction types)
    "GO:0030057": "GO_DESMOSOME",
    "GO:0005912": "GO_ADHERENS_JUNCTION",
    "GO:0061436": "GO_ESTABLISHMENT_OF_SKIN_BARRIER",
}

# Compact epithelial-barrier panel centered on Tacstd2 / Cldn4 (A8 claim).
# Kept small so GSEA is not just "any epithelium".
BARRIER_CORE = [
    "Tacstd2",
    "Cldn1", "Cldn3", "Cldn4", "Cldn7", "Cldn8",
    "Ocln", "Marveld2", "Marveld3",
    "Tjp1", "Tjp2", "Tjp3",
    "F11r", "Jam2", "Jam3",
    "Cgn", "Crb3",
    "Krt5", "Krt7", "Krt8", "Krt17", "Krt18", "Krt19",
    "Cdh1", "Epcam",
]

CURATED_TJ_CORE = [
    # claudins
    *[f"Cldn{i}" for i in list(range(1, 28)) + [34]],
    "Cldn34a", "Cldn34b1", "Cldn34b2", "Cldn34b3", "Cldn34b4", "Cldn34c1",
    # occludin / MARVEL family, JAMs, scaffolds, angulins
    "Ocln", "Marveld2", "Marveld3",
    "Tjp1", "Tjp2", "Tjp3",
    "F11r", "Jam2", "Jam3", "Igsf5", "Esam",
    "Cgn", "Cgnl1", "Crb3", "Ildr1", "Ildr2", "Lsr",
    "Pard3", "Pard6a", "Pard6b", "Pard6g", "Mpdz", "Patj", "Inadl",
    "Amot", "Amotl1", "Amotl2", "Mpp5", "Symp",
]


def parse_obo_children(path):
    """Return child -> parents map using is_a and part_of, ignoring obsoletes."""
    parents = defaultdict(set)
    cur_id, obsolete, alt_ids = None, False, []
    alt_map = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur_id, obsolete, alt_ids = None, False, []
            elif line.startswith("id: GO:"):
                cur_id = line[4:]
            elif line.startswith("alt_id: GO:") and cur_id:
                alt_map[line[8:]] = cur_id
            elif line.startswith("is_obsolete: true"):
                obsolete = True
            elif line.startswith("is_a: GO:") and cur_id and not obsolete:
                parents[cur_id].add(line[6:].split(" ")[0])
            elif line.startswith("relationship: part_of GO:") and cur_id and not obsolete:
                parents[cur_id].add(line.split(" ")[2])
    children = defaultdict(set)
    for child, ps in parents.items():
        for p in ps:
            children[p].add(child)
    return children, alt_map


def descendants(term, children):
    seen, stack = {term}, [term]
    while stack:
        t = stack.pop()
        for c in children.get(t, ()):  # BFS/DFS over is_a + part_of
            if c not in seen:
                seen.add(c)
                stack.append(c)
    return seen


def main():
    children, alt_map = parse_obo_children(ANNOT_DIR / "go-basic.obo")

    term_closure = {}  # descendant GO id -> set of target names it rolls up to
    for tid, name in TARGET_TERMS.items():
        for d in descendants(tid, children):
            term_closure.setdefault(d, set()).add(name)

    sets = defaultdict(set)
    all_symbols = set()
    with gzip.open(ANNOT_DIR / "mgi.gaf.gz", "rt") as fh:
        for line in fh:
            if line.startswith("!"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 13:
                continue
            symbol, qualifier, goid, evidence = f[2], f[3], f[4], f[6]
            if "NOT" in qualifier or evidence == "ND":
                continue
            goid = alt_map.get(goid, goid)
            all_symbols.add(symbol)
            for name in term_closure.get(goid, ()):
                sets[name].add(symbol)

    # curated keratin family: canonical Krt gene symbols annotated in MGI
    krt_re = re.compile(r"^Krt\d+[a-z0-9\-]*$")
    sets["CURATED_KERATIN_FAMILY"] = {s for s in all_symbols if krt_re.match(s)}
    sets["CURATED_TIGHT_JUNCTION_CORE"] = {
        s for s in CURATED_TJ_CORE if s in all_symbols
    }
    sets["BARRIER_CORE_TACSTD2_CLDN4"] = {s for s in BARRIER_CORE if s in all_symbols}

    gmt = OUT_DIR / "mouse_keratin_tj.gmt"
    with open(gmt, "w") as out, open(OUT_DIR / "geneset_sizes.tsv", "w") as sz:
        sz.write("gene_set\tn_genes\n")
        for name in sorted(sets):
            genes = sorted(sets[name])
            out.write(f"{name}\tGO/MGI-derived mouse gene set\t" + "\t".join(genes) + "\n")
            sz.write(f"{name}\t{len(genes)}\n")
            print(f"{name}: {len(genes)} genes")


if __name__ == "__main__":
    main()
