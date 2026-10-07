#!/usr/bin/env bash
# Download every public input used by the GSE334497 reanalysis.
# GEO files go to raw/ (committed); large reference files go to cache/ (git-ignored).
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
RAW="$HERE/raw"; CACHE="$HERE/cache"
mkdir -p "$RAW" "$CACHE"

geo=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE334nnn/GSE334497
curl -sSfL -o "$RAW/GSE334497_normalized_counts.csv.gz" "$geo/suppl/GSE334497_normalized_counts.csv.gz"
curl -sSfL -o "$RAW/GSE334497_series_matrix.txt.gz"     "$geo/matrix/GSE334497_series_matrix.txt.gz"
curl -sSfL -o "$RAW/GSE334497_geo_brief.txt" \
  "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334497&targ=self&form=text&view=brief"

# Paper: Wu B ... Ellisen LW, J Immunother Cancer 2026;14(4):e012265, PMID 41932810, PMC13052784
curl -sSfL -o "$CACHE/pubmed_41932810.xml" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=41932810&retmode=xml"
curl -sSfL -o "$CACHE/PMC13052784.xml" \
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC13052784&retmode=xml"
curl -sSfL -o "$CACHE/PMC13052784_supp.zip" \
  "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13052784/supplementaryFiles"

# Annotation used by the submitters (GRCm38) -> Ensembl release 102 gene models.
curl -sSfL -o "$CACHE/Mus_musculus.GRCm38.102.gtf.gz" \
  "https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"

# MSigDB 2024.1 (human collections + mouse-Ensembl -> human-ortholog chip, as used by GSEA desktop).
ms=https://data.broadinstitute.org/gsea-msigdb/msigdb
curl -sSfL -o "$CACHE/msigdb.v2024.1.Hs.symbols.gmt" "$ms/release/2024.1.Hs/msigdb.v2024.1.Hs.symbols.gmt"
curl -sSfL -o "$CACHE/Mouse_Ensembl_Gene_ID_Human_Orthologs_MSigDB.v2024.1.Hs.chip" \
  "$ms/annotations/human/Mouse_Ensembl_Gene_ID_Human_Orthologs_MSigDB.v2024.1.Hs.chip"
for c in h.all c2.cgp c2.cp.biocarta c2.cp.kegg_legacy c2.cp.reactome c2.cp.wikipathways c5.go.bp c5.go.cc c5.go.mf; do
  curl -sSfL -o "$CACHE/$c.v2024.1.Hs.symbols.gmt" "$ms/release/2024.1.Hs/$c.v2024.1.Hs.symbols.gmt"
done

(cd "$CACHE" && python3 - <<'EOF'
import zipfile; zipfile.ZipFile('PMC13052784_supp.zip').extractall('PMC13052784_supp')
EOF
)

zcat "$CACHE/Mus_musculus.GRCm38.102.gtf.gz" | awk -F'\t' '$3=="gene"' | python3 -c '
import sys, re
print("ensembl_gene_id\tgene_name\tgene_biotype\tchrom")
for line in sys.stdin:
    f = line.rstrip("\n").split("\t"); a = f[8]
    gid = re.search(r"gene_id \"([^\"]+)\"", a).group(1)
    m = re.search(r"gene_name \"([^\"]+)\"", a)
    b = re.search(r"gene_biotype \"([^\"]+)\"", a).group(1)
    print(gid, m.group(1) if m else "", b, f[0], sep="\t")
' > "$CACHE/ensembl102_mouse_genes.tsv"

(cd "$RAW" && md5sum GSE334497_normalized_counts.csv.gz GSE334497_series_matrix.txt.gz) > "$RAW/MD5SUMS"
(cd "$CACHE" && md5sum Mus_musculus.GRCm38.102.gtf.gz msigdb.v2024.1.Hs.symbols.gmt \
   Mouse_Ensembl_Gene_ID_Human_Orthologs_MSigDB.v2024.1.Hs.chip) >> "$RAW/MD5SUMS"
cat "$RAW/MD5SUMS"
