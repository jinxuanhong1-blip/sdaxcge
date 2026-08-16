# Package install attempt (this environment)

| Tool | Result |
|---|---|
| **R / Rscript** | Not on PATH. `which R` failed. |
| **CellChat** | Not installed. Requires R + Bioconductor. **No CellChat output was written.** |
| **liana 1.8.1** | Installed with scanpy 1.12.3 / anndata 0.13.2 via `pip`. |
| **Documented CellPhoneDB score** | Always run (numpy/pandas/scipy). |

If LIANA import or `li.mt.cellphonedb` fails at runtime, `results/summary.json`
records `liana_status` and the documented score remains the primary table.
