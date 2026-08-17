#!/usr/bin/env python3
"""Subset GSE127465 tumor epithelium and gate malignant+CLDN4.

If author PatientN-specific malignant cells with CLDN4>0 are missing,
write a stop note (honest n=0) and do not build a trajectory object.
"""
from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import sys

import numpy as np
import pandas as pd
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gene_sets import AUTHOR_LEFTOVER_EPI, AUTHOR_MALIGNANT_SUBSTR  # noqa: E402

N_CELLS = 54773
N_GENES = 41861
PAT_SPEC_RE = re.compile(r"specific", re.I)


def load_genes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        genes = [ln.strip().split("\t")[0] for ln in handle if ln.strip()]
    if len(genes) != N_GENES:
        raise SystemExit(f"genes {len(genes)} != {N_GENES}")
    return genes


def load_meta(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    if len(df) != N_CELLS:
        raise SystemExit(f"metadata rows {len(df)} != {N_CELLS}")
    return df


def is_tumor_epithelium(meta: pd.DataFrame) -> pd.Series:
    tumor = meta["Tissue"].astype(str) == "tumor"
    maj = meta["Major cell type"].astype(str)
    leftover = maj.isin(AUTHOR_LEFTOVER_EPI)
    malignant = maj.str.contains(AUTHOR_MALIGNANT_SUBSTR, case=False, regex=False)
    return tumor & (leftover | malignant)


def stream_mtx_rows(mtx_path: Path, keep_rows_0: np.ndarray, n_genes: int) -> sparse.csr_matrix:
    """Matrix Market 54773 cells x 41861 genes, 1-based. Keep selected cell rows."""
    keep_set = set(int(i) + 1 for i in keep_rows_0)  # 1-based
    row_map = {int(i) + 1: k for k, i in enumerate(keep_rows_0)}
    rs: list[int] = []
    cs: list[int] = []
    vs: list[float] = []
    with gzip.open(mtx_path, "rt") as handle:
        for line in handle:
            if line.startswith("%"):
                continue
            banner = line.split()
            nrows, ncols, nnz = map(int, banner)
            if nrows != N_CELLS or ncols != n_genes:
                raise SystemExit(f"MTX banner {nrows}x{ncols} != {N_CELLS}x{n_genes}")
            print(f"MTX banner {nrows} x {ncols} nnz={nnz}", flush=True)
            break
        for i, line in enumerate(handle, 1):
            parts = line.split()
            if len(parts) < 3:
                continue
            r = int(parts[0])
            if r not in keep_set:
                continue
            cs.append(int(parts[1]) - 1)
            rs.append(row_map[r])
            vs.append(float(parts[2]))
            if i % 10_000_000 == 0:
                print(f"  MTX lines {i:,} kept={len(vs):,}", flush=True)
    print(f"kept nnz={len(vs):,} cells={len(keep_rows_0)}", flush=True)
    return sparse.csr_matrix(
        (np.asarray(vs, dtype=np.float32), (rs, cs)),
        shape=(len(keep_rows_0), n_genes),
    )


def write_stop(outdir: Path, finding: Path, gate: dict) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(exist_ok=True)
    (outdir / "figures").mkdir(exist_ok=True)
    (outdir / "tables" / "STOP_no_malignant_cldn4.json").write_text(
        json.dumps(gate, indent=2)
    )
    finding.parent.mkdir(parents=True, exist_ok=True)
    finding.write_text(
        "\n".join(
            [
                "# Finding — GSE127465 REAL Slingshot/PAGA, CLDN4 only",
                "",
                "**STOP.** Author malignant + CLDN4 is missing. Trajectory was not run.",
                "",
                "ADDITIVE. **CLDN4 only.** Zilionis et al., *Immunity* 2019, PMID 30979687 "
                "(GSE127465 human inDrops). No TACSTD2∩CLDN4 dual-high gate. "
                "Root would not have been CLDN4-high. n=7 patients would have been thin.",
                "",
                "## Honest n",
                "",
                f"- Human patients deposited: **7** (p1–p7).",
                f"- Tumor cells in metadata: **{gate.get('n_tumor')}**.",
                f"- Tumor epithelium (author Type I/II + club + ciliated + Patient*-specific): "
                f"**{gate.get('n_tumor_epithelium')}**.",
                f"- Author malignant (Patient*-specific in tumor): **{gate.get('n_malignant')}**.",
                f"- Malignant with CLDN4 > 0: **{gate.get('n_malignant_cldn4_pos')}**.",
                f"- CLDN4 gene present: **{gate.get('cldn4_present')}**.",
                "- Inferential n for a trajectory test: **n=0** (gate failed).",
                "",
                "Done criterion: stop note with honest n=0. No lineage table.",
                "",
            ]
        )
    )
    print(json.dumps({"ok": False, "stop": True, "gate": gate}, indent=2))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/tmp/gse127465_slingshot")
    ap.add_argument("--out", default="/tmp/gse127465_slingshot/epithelium.h5ad")
    ap.add_argument(
        "--results",
        default="methods/gse127465_slingshot_real_cldn4/results",
    )
    ap.add_argument(
        "--finding",
        default="methods/gse127465_slingshot_real_cldn4/FINDING.md",
    )
    args = ap.parse_args()
    data = Path(args.data)
    meta = load_meta(data / "GSE127465_human_cell_metadata_54773x25.tsv.gz")
    genes = load_genes(data / "GSE127465_gene_names_human_41861.tsv.gz")
    if "CLDN4" not in genes:
        gate = {
            "cldn4_present": False,
            "n_tumor": int((meta["Tissue"].astype(str) == "tumor").sum()),
            "n_tumor_epithelium": int(is_tumor_epithelium(meta).sum()),
            "n_malignant": int(
                (
                    (meta["Tissue"].astype(str) == "tumor")
                    & meta["Major cell type"].astype(str).str.contains("specific", case=False)
                ).sum()
            ),
            "n_malignant_cldn4_pos": 0,
            "reason": "CLDN4 absent from gene table",
        }
        write_stop(Path(args.results), Path(args.finding), gate)
        raise SystemExit(0)

    epi_mask = is_tumor_epithelium(meta)
    keep_idx = np.flatnonzero(epi_mask.to_numpy())
    print(f"tumor epithelium cells={keep_idx.size}", flush=True)
    X = stream_mtx_rows(
        data / "GSE127465_human_counts_normalized_54773x41861.mtx.gz",
        keep_idx,
        len(genes),
    )
    obs = meta.loc[epi_mask].copy()
    obs["cell_id"] = (
        obs["Library"].astype(str) + "_" + obs["Barcode"].astype(str)
    )
    obs = obs.reset_index(drop=True)
    maj = obs["Major cell type"].astype(str)
    obs["is_malignant"] = maj.str.contains("specific", case=False, regex=False)
    obs["is_type2"] = maj.eq("Type II cells")
    obs["lineage"] = np.where(
        obs["is_malignant"],
        "Malignant",
        maj,
    )
    cldn4_i = genes.index("CLDN4")
    cldn4 = np.asarray(X[:, cldn4_i].todense()).ravel()
    obs["CLDN4"] = cldn4
    n_mal = int(obs["is_malignant"].sum())
    n_mal_pos = int(((obs["is_malignant"]) & (obs["CLDN4"] > 0)).sum())
    gate = {
        "cldn4_present": True,
        "n_tumor": int((meta["Tissue"].astype(str) == "tumor").sum()),
        "n_tumor_epithelium": int(len(obs)),
        "n_malignant": n_mal,
        "n_malignant_cldn4_pos": n_mal_pos,
        "n_type2": int(obs["is_type2"].sum()),
        "patients": sorted(obs["Patient"].astype(str).unique().tolist()),
        "lineage_counts": obs["lineage"].value_counts().to_dict(),
        "reason": None if n_mal_pos > 0 else "no CLDN4+ malignant cells",
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.results).mkdir(parents=True, exist_ok=True)
    (Path(args.results) / "tables").mkdir(exist_ok=True)
    (Path(args.results) / "tables" / "gate.json").write_text(json.dumps(gate, indent=2))
    if n_mal == 0 or n_mal_pos == 0:
        write_stop(Path(args.results), Path(args.finding), gate)
        raise SystemExit(0)

    import anndata as ad

    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=pd.Index(genes, name="gene")))
    adata.obs_names = obs["cell_id"].astype(str).tolist()
    adata.obs_names_make_unique()
    adata.layers["normalized"] = adata.X.copy()
    adata.write_h5ad(args.out)
    print(json.dumps({"ok": True, "out": args.out, "gate": gate}, indent=2, default=str))


if __name__ == "__main__":
    main()
