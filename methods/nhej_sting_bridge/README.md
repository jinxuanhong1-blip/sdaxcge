# NHEJ-component loss vs IFN / STING / APM

Public test of the middle step in CLDN4 → NHEJ → IFN: after an NHEJ component
is lost, do interferon, STING-axis, or MHC-I antigen-presentation transcripts
go up?

Series (deposited counts, not a precomputed DEG table):

- **GSE180581** — HEK293T monoallelic knockout plus siRNA of Ku70 (`XRCC6`),
  Ku80 (`XRCC5`), or DNA-PKcs (`PRKDC`), each versus siControl in the same
  heterozygous line (n = 3). PMID 34849385, 35430316.
- **GSE135274** — HeLa `XRCC4(-/-)` versus scrambled gRNA, no mirin (n = 2
  experiments). Read 1 and read 2 are technical and are summed. Every sample
  was transfected with a TALEN NHEJ reporter. PMID 35054780.

```bash
python3 methods/nhej_sting_bridge/analyze.py
```

The script downloads the GEO count files and the HGNC table into
`methods/nhej_sting_bridge/data/` (gitignored) and writes tables and figures
under `results/nhej_sting_bridge/`. Narrative: `results/nhej_sting_bridge/RESULTS.md`.
