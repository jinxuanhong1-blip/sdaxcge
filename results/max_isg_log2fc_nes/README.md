# Cell-line ISG maxima

Largest Hallmark interferon-stimulated-gene log2 fold change and NES for:

- HCT116 Ku80-AID and Ku86-flox, GSE294709
- HCT116 DNA-PKcs knockout, GSE285698
- MCF-7 53BP1 knockout, GSE84986

HEK293 Ku70 from GSE294709 is included as a separate cell line in the same series.

See `WRITEUP.md`. Run `python3 scripts/max_isg_log2fc_nes/analyze.py`. Raw GEO matrices are downloaded to `/tmp/max_isg` and are not committed.
