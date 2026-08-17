# GSE68465 processed source

Public GEO series matrix + official GPL96 annotation only. No CEL / MAS5 reprocess. No ICI labels.

| File | URL |
|---|---|
| `GSE68465_series_matrix.txt.gz` | https://ftp.ncbi.nlm.nih.gov/geo/series/GSE68nnn/GSE68465/matrix/GSE68465_series_matrix.txt.gz |
| `GPL96.annot.gz` | https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL96/annot/GPL96.annot.gz |

Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68465

Paper: Shedden et al., *Nat Med* 2008; PMID 18641660 (Director's Challenge / jacob-00182). The paper text says 442 LUAD. The deposited matrix has **462** arrays: **443** `disease_state = Lung Adenocarcinoma` and **19** `Normal` (Stratagene). This folder uses the GEO filter (**n = 443**). Do not write n=462 as the tumour n.

Platform: GPL96 Affymetrix HG-U133A. Deposited values are MAS5 intensities (range ~0.06–1.3e5; no negatives). Primary transform here: `log2(MAS5)`.

CLDN4 only. TACSTD2 is not an anchor. Named CLDN4 probe is `201428_at`.
