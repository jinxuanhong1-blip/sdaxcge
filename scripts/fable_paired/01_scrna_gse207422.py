"""GSE207422 scRNA-seq: TACSTD2/CLDN4 dynamics in the epithelial compartment.

Design of the public data:
  * 3 pre-treatment biopsies + 12 post-treatment surgical resections (15 patients),
    NSCLC treated with neoadjuvant anti-PD-1 + platinum chemotherapy.
  * Pathologic response annotated per sample (MPR vs NMPR/pCR).

Because TACSTD2 (TROP2) and CLDN4 are epithelial genes, we restrict to the
epithelial compartment (EPCAM+ / PTPRC-) and quantify their expression per
sample (pseudobulk) so that changes are not merely driven by shifts in immune
infiltration. We contrast pre vs post treatment and, within post-treatment
tumors, responders (MPR) vs non-responders (NMPR).

The 176 MB UMI matrix is streamed once (awk for library sizes, grep for the
marker/target rows); a compact processed table is cached under results/.
"""
from __future__ import annotations

import io
import subprocess
import shutil

import numpy as np
import pandas as pd

import config as C
from common import compare_groups, strip_box, savefig, PALETTE, pstars
import matplotlib.pyplot as plt


CACHE_PANEL = C.TABLES_DIR / "gse207422_scrna_panel_percell.csv.gz"
CACHE_LIBSIZE = C.TABLES_DIR / "gse207422_scrna_libsize.csv.gz"


def _stream_extract() -> tuple[pd.DataFrame, pd.Series]:
    """Return (panel genes x cells raw UMI) and per-cell library size."""
    matrix = C.raw_path("GSE207422_sc_expr")
    if CACHE_PANEL.exists() and CACHE_LIBSIZE.exists():
        print(f"[cache] reading {CACHE_PANEL.name} / {CACHE_LIBSIZE.name}")
        panel = pd.read_csv(CACHE_PANEL, index_col=0)
        lib = pd.read_csv(CACHE_LIBSIZE, index_col=0)["libsize"]
        return panel, lib

    if not (shutil.which("awk") and shutil.which("grep")):
        raise RuntimeError("awk and grep are required for the streaming extraction")

    # 1) per-cell library size (sum over all genes)
    print("[sc  ] computing per-cell library sizes (awk)...")
    awk = (
        r"""NR==1{for(i=2;i<=NF;i++)h[i]=$i; n=NF; next}"""
        r"""{for(i=2;i<=NF;i++)s[i]+=$i}"""
        r"""END{print "cell\tlibsize"; for(i=2;i<=n;i++) printf "%s\t%d\n", h[i], s[i]}"""
    )
    p1 = subprocess.run(
        f"zcat {matrix} | awk -F'\\t' '{awk}'",
        shell=True, capture_output=True, text=True, check=True,
    )
    lib = pd.read_csv(io.StringIO(p1.stdout), sep="\t", index_col="cell")["libsize"]

    # 2) panel gene rows
    print("[sc  ] extracting marker/target gene rows (grep)...")
    genes_re = "|".join(C.SCRNA_PANEL)
    p2 = subprocess.run(
        f"zcat {matrix} | grep -P '^(Gene|{genes_re})\\t'",
        shell=True, capture_output=True, text=True, check=True,
    )
    panel = pd.read_csv(io.StringIO(p2.stdout), sep="\t", index_col="Gene")
    panel = panel.loc[[g for g in C.SCRNA_PANEL if g in panel.index]]

    # align columns
    panel = panel[lib.index]
    panel.to_csv(CACHE_PANEL)
    lib.to_frame("libsize").to_csv(CACHE_LIBSIZE)
    print(f"[sc  ] cached panel {panel.shape} and libsize ({lib.shape[0]} cells)")
    return panel, lib


def _load_meta() -> pd.DataFrame:
    meta = pd.read_excel(C.raw_path("GSE207422_sc_meta"))
    meta = meta.dropna(subset=["Sample"]).copy()
    meta = meta[meta["Sample"].astype(str).str.startswith("BD_immune")]
    res = meta["Resource"].str.lower()
    meta["timepoint"] = np.where(res.str.contains("pre"), "pre", "post")
    pr = meta["Pathologic Response"].astype(str).str.upper()
    meta["response"] = np.where(pr.str.contains("NMPR"), "non_responder",
                        np.where(pr.str.contains("MPR") | pr.str.contains("PCR"),
                                 "responder", "unknown"))
    return meta.set_index("Sample")


def main() -> None:
    panel, lib = _stream_extract()
    meta = _load_meta()

    cells = panel.columns
    sample = pd.Series(cells, index=cells).str.rsplit("_", n=1).str[0]  # BD_immuneNN
    libsize = lib.reindex(cells).astype(float)

    # normalize: counts per 10k, log1p
    def norm(gene: str) -> pd.Series:
        raw = panel.loc[gene].reindex(cells).astype(float)
        return np.log1p(raw / libsize.replace(0, np.nan) * 1e4)

    epcam = panel.loc["EPCAM"].reindex(cells).astype(float)
    ptprc = panel.loc["PTPRC"].reindex(cells).astype(float)
    is_epi = (epcam > 0) & (ptprc == 0)

    cell_df = pd.DataFrame({
        "cell": cells,
        "sample": sample.values,
        "libsize": libsize.values,
        "is_epithelial": is_epi.values,
        "TACSTD2": norm("TACSTD2").values,
        "CLDN4": norm("CLDN4").values,
        "EPCAM_raw": epcam.values,
        "PTPRC_raw": ptprc.values,
    })
    cell_df = cell_df.join(meta[["timepoint", "response", "Pathologic Response",
                                 "RECIST", "Patient"]], on="sample")
    cell_df.to_csv(C.TABLES_DIR / "gse207422_scrna_cells_annotated.csv.gz", index=False)

    epi = cell_df[cell_df["is_epithelial"]].copy()
    print(f"[sc  ] epithelial cells: {len(epi):,} / {len(cell_df):,} "
          f"({100*len(epi)/len(cell_df):.1f}%)")

    # sample-level pseudobulk within epithelial cells
    grp = epi.groupby("sample")
    pseudo = grp.agg(
        n_epi=("cell", "size"),
        TACSTD2_mean=("TACSTD2", "mean"),
        CLDN4_mean=("CLDN4", "mean"),
        TACSTD2_pctpos=("TACSTD2", lambda s: float((s > 0).mean() * 100)),
        CLDN4_pctpos=("CLDN4", lambda s: float((s > 0).mean() * 100)),
    )
    pseudo = pseudo.join(meta[["timepoint", "response", "Pathologic Response", "RECIST"]])
    # require a minimal epithelial cell count for a stable pseudobulk value
    pseudo_use = pseudo[pseudo["n_epi"] >= 20].copy()
    pseudo.to_csv(C.TABLES_DIR / "gse207422_scrna_pseudobulk_epithelial.csv")
    print(pseudo.sort_values("timepoint").to_string())

    # ---------- statistics ----------
    rows = []
    for gene in C.TARGET_GENES:
        col = f"{gene}_mean"
        pre = pseudo_use.loc[pseudo_use.timepoint == "pre", col].values
        post = pseudo_use.loc[pseudo_use.timepoint == "post", col].values
        rows.append(compare_groups(
            "GSE207422_scRNA_epi", gene, "post_vs_pre(sample pseudobulk)",
            post, pre, "post", "pre").as_row())

        post_df = pseudo_use[pseudo_use.timepoint == "post"]
        r = post_df.loc[post_df.response == "responder", col].values
        nr = post_df.loc[post_df.response == "non_responder", col].values
        rows.append(compare_groups(
            "GSE207422_scRNA_epi", gene, "post: responder_vs_nonresponder",
            r, nr, "responder(MPR)", "non_responder(NMPR)").as_row())

    stats_df = pd.DataFrame(rows)
    stats_df.to_csv(C.TABLES_DIR / "gse207422_scrna_stats.csv", index=False)
    print(stats_df.to_string())

    # ---------- figure ----------
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for j, gene in enumerate(C.TARGET_GENES):
        col = f"{gene}_mean"
        # pre vs post
        groups = {
            "pre": pseudo_use.loc[pseudo_use.timepoint == "pre", col].values,
            "post": pseudo_use.loc[pseudo_use.timepoint == "post", col].values,
        }
        strip_box(axes[0, j], groups,
                  {"pre": PALETTE["pre"], "post": PALETTE["post"]},
                  f"{gene} (log1p CP10k, epi mean)",
                  f"GSE207422 {gene}: pre vs post")
        # responder vs non-responder among post
        post_df = pseudo_use[pseudo_use.timepoint == "post"]
        groups2 = {
            "responder\n(MPR)": post_df.loc[post_df.response == "responder", col].values,
            "non-responder\n(NMPR)": post_df.loc[post_df.response == "non_responder", col].values,
        }
        strip_box(axes[1, j], groups2,
                  {"responder\n(MPR)": PALETTE["responder"],
                   "non-responder\n(NMPR)": PALETTE["non_responder"]},
                  f"{gene} (log1p CP10k, epi mean)",
                  f"GSE207422 {gene}: post-treatment by response")
    fig.suptitle("GSE207422 (NSCLC neoadjuvant anti-PD-1 + chemo) - epithelial compartment",
                 fontsize=11)
    savefig(fig, C.FIGURES_DIR / "fig1_gse207422_scrna_epithelial.png")


if __name__ == "__main__":
    main()
