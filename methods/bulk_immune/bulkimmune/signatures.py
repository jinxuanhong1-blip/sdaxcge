"""Signature library: curated published sets + MSigDB + user GMTs.

Design rule: a signature is only hard-coded here when the source paper prints
the full gene list in its text (the IFN-gamma, GEP and TLS sets). Everything
larger comes from a downloaded, version-stamped file. That way a reviewer can
always answer "which version of this gene set did you use".
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

__all__ = ["read_gmt", "write_gmt", "SignatureLibrary", "cytolytic_activity", "signature_zscore_mean"]


def read_gmt(path: str | Path) -> dict[str, list[str]]:
    sets: dict[str, list[str]] = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            sets[fields[0]] = [g for g in fields[2:] if g]
    return sets


def write_gmt(sets: dict[str, Iterable[str]], path: str | Path, description: str = "na") -> None:
    with open(path, "w") as handle:
        for name, genes in sets.items():
            handle.write("\t".join([name, description, *genes]) + "\n")


@dataclass
class SignatureLibrary:
    """Named gene sets plus the citation for each."""

    sets: dict[str, list[str]]
    meta: dict[str, dict]

    @classmethod
    def load(
        cls,
        curated_json: str | Path,
        msigdb_gmt: str | Path | None = None,
        msigdb_keep: Iterable[str] | None = None,
        extra_gmt: Iterable[str | Path] = (),
    ) -> "SignatureLibrary":
        with open(curated_json) as handle:
            payload = json.load(handle)
        sets: dict[str, list[str]] = {}
        meta: dict[str, dict] = {}
        for name, entry in payload["signatures"].items():
            sets[name] = list(entry["genes"])
            meta[name] = {
                "citation": entry.get("citation", ""),
                "note": entry.get("note", ""),
                "family": entry.get("family", ""),
                "source": "curated_signatures.json",
            }

        if msigdb_gmt is not None:
            hallmark = read_gmt(msigdb_gmt)
            wanted = set(msigdb_keep) if msigdb_keep is not None else set(hallmark)
            for name, genes in hallmark.items():
                if name in wanted:
                    sets[name] = genes
                    meta[name] = {
                        "citation": "MSigDB (Liberzon A et al. Cell Syst 2015)",
                        "note": "",
                        "family": "hallmark",
                        "source": str(Path(msigdb_gmt).name),
                    }

        for gmt in extra_gmt:
            for name, genes in read_gmt(gmt).items():
                sets[name] = genes
                meta[name] = {
                    "citation": "user-supplied",
                    "note": "",
                    "family": "custom",
                    "source": str(Path(gmt).name),
                }

        return cls(sets=sets, meta=meta)

    def by_family(self, family: str) -> dict[str, list[str]]:
        return {k: v for k, v in self.sets.items() if self.meta[k].get("family") == family}

    def harmonise(self, hgnc) -> dict[str, dict[str, str]]:
        """Map every signature onto current HGNC symbols; return the renames."""
        renames: dict[str, dict[str, str]] = {}
        for name, genes in self.sets.items():
            mapped, renamed = hgnc.harmonise(genes)
            self.sets[name] = [
                new if new is not None else old for old, new in zip(genes, mapped)
            ]
            if renamed:
                renames[name] = renamed
        return renames

    def table(self) -> pd.DataFrame:
        rows = [
            {
                "signature": name,
                "n_genes": len(genes),
                "family": self.meta[name].get("family", ""),
                "source": self.meta[name].get("source", ""),
                "citation": self.meta[name].get("citation", ""),
            }
            for name, genes in self.sets.items()
        ]
        return pd.DataFrame(rows)


def cytolytic_activity(expr_tpm: pd.DataFrame, offset: float = 1.0) -> pd.Series:
    """CYT = geometric mean of GZMA and PRF1 in TPM (Rooney et al. 2015).

    Deliberately not an ssGSEA score: the published definition is a geometric
    mean on the linear scale, and a two-gene "ssGSEA" score is close to
    meaningless anyway (the rank walk has nothing to integrate over).
    """
    missing = [g for g in ("GZMA", "PRF1") if g not in expr_tpm.index]
    if missing:
        raise KeyError(f"missing {missing} for cytolytic activity")
    if float(np.nanmax(expr_tpm.to_numpy())) < 50:
        raise ValueError("CYT expects linear TPM, not log-transformed values")
    values = expr_tpm.loc[["GZMA", "PRF1"]] + offset
    return np.exp(np.log(values).mean(axis=0)).rename("CYT")


def signature_zscore_mean(
    expr_log: pd.DataFrame,
    gene_sets: dict[str, list[str]],
    min_genes: int = 2,
) -> pd.DataFrame:
    """Mean of per-gene z-scores -- the other common way to score a signature.

    Worth computing alongside ssGSEA: the two answer slightly different
    questions (rank position of the set vs average standardised level) and a
    result that only holds under one of them is fragile. Like ssGSEA's
    normalisation, the z-score is taken across the samples in hand, so this
    score is also cohort-relative.
    """
    z = expr_log.sub(expr_log.mean(axis=1), axis=0).div(
        expr_log.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    out = {}
    for name, genes in gene_sets.items():
        present = [g for g in dict.fromkeys(genes) if g in z.index]
        if len(present) < min_genes:
            continue
        out[name] = z.loc[present].mean(axis=0, skipna=True)
    return pd.DataFrame(out)
