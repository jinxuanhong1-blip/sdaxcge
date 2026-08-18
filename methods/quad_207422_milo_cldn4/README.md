# methods/quad_207422_milo_cldn4

Additive **CLDN4-only** Milo-style neighbourhood DA on the public QUAD
merge: GSE207422 + GSE131907 + GSE148071 + GSE205335.

- GSE207422 is **included** by importing PR #323 (null, n=7). That
  207422-only Milo is **not** re-run.
- No TACSTD2 gate. No dual-high.
- Independent unit = patient / sample.
- Graph **per dataset** (GSE131907 also splits tLung / mBrain). Harmony
  was not used.
- miloR is not used. Same documented fallback as the single-dataset
  folders: sample-level Spearman on neighbourhood proportions +
  SpatialFDR k-distance weights.

`FINDING.md` and `tables/nhoods.tsv` are written by `scripts/merge_quad.py`.
