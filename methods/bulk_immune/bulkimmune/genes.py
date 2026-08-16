"""Gene identifier harmonisation.

Every deconvolution method here is keyed on HGNC symbols, and every one of them
silently returns a *number* when the identifier namespace is wrong -- the score
is simply computed on whatever fraction of the signature happened to match. The
functions below make the matching explicit and auditable.

Two failure modes dominate in practice:

* **Ensembl versions.** ``ENSG00000000003.10`` never matches ``ENSG00000000003``.
* **Symbol drift.** Signature files published between 2011 and 2018 still use
  withdrawn symbols (``FYB`` -> ``FYB1``, ``KIAA1524`` -> ``CIP2A``, ...). Mapping
  aliases and previous symbols back onto current ones typically recovers a
  handful of genes per signature, which matters for the 6- and 12-gene sets.

The alias table is the HGNC complete set (downloaded by
``scripts/00_fetch_resources.py``); nothing is hard-coded.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

__all__ = ["HGNC", "strip_ensembl_version", "collapse_duplicates", "match_report"]


def strip_ensembl_version(ids: Iterable[str]) -> list[str]:
    """``ENSG00000000003.10`` -> ``ENSG00000000003`` (leaves other ids alone)."""
    out = []
    for gid in ids:
        s = str(gid)
        if s.startswith("ENS") and "." in s:
            s = s.split(".", 1)[0]
        out.append(s)
    return out


@dataclass
class HGNC:
    """HGNC symbol/Ensembl lookup built from ``hgnc_complete_set.txt``."""

    table: pd.DataFrame
    _ensembl_to_symbol: dict[str, str] = field(default_factory=dict, repr=False)
    _alias_to_symbol: dict[str, str] = field(default_factory=dict, repr=False)
    _approved: set[str] = field(default_factory=set, repr=False)

    @classmethod
    def from_file(cls, path: str | Path) -> "HGNC":
        table = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
        obj = cls(table=table)
        obj._build()
        return obj

    def _build(self) -> None:
        t = self.table
        approved = t.loc[t["status"].eq("Approved")] if "status" in t else t
        self._approved = set(approved["symbol"].dropna())

        ens = approved.dropna(subset=["ensembl_gene_id"])
        self._ensembl_to_symbol = dict(zip(ens["ensembl_gene_id"], ens["symbol"]))

        # Previous symbols are authoritative renames; aliases are weaker. Fill
        # previous symbols last so that they win over ambiguous aliases.
        alias: dict[str, str] = {}
        for column in ("alias_symbol", "prev_symbol"):
            if column not in approved:
                continue
            for symbol, blob in zip(approved["symbol"], approved[column]):
                if not isinstance(blob, str):
                    continue
                for old in blob.split("|"):
                    old = old.strip()
                    # Never let an alias shadow a currently approved symbol.
                    if old and old not in self._approved:
                        alias[old] = symbol
        self._alias_to_symbol = alias

    def ensembl_to_symbol(self, ids: Sequence[str]) -> list[str | None]:
        clean = strip_ensembl_version(ids)
        return [self._ensembl_to_symbol.get(g) for g in clean]

    def canonical_symbol(self, symbol: str) -> str | None:
        """Approved symbol for ``symbol``, resolving previous symbols/aliases."""
        if symbol in self._approved:
            return symbol
        return self._alias_to_symbol.get(symbol)

    def harmonise(self, symbols: Sequence[str]) -> tuple[list[str | None], dict[str, str]]:
        """Map a list of symbols to approved symbols; also return the renames."""
        mapped: list[str | None] = []
        renamed: dict[str, str] = {}
        for s in symbols:
            c = self.canonical_symbol(s)
            if c is not None and c != s:
                renamed[s] = c
            mapped.append(c)
        return mapped, renamed


def collapse_duplicates(
    expr: pd.DataFrame,
    method: str = "max_mean",
) -> pd.DataFrame:
    """Collapse rows sharing a gene identifier.

    ``max_mean`` keeps, for each identifier, the row with the highest mean
    expression (the WGCNA/``collapseRows`` convention used by most deconvolution
    papers). ``sum`` adds them, which is only correct for counts of distinct
    transcripts of the same gene. ``mean`` averages them, which dilutes a real
    probe with a dead one and is almost never what you want -- it is offered
    because some published pipelines use it and reproducing them requires it.
    """
    if not expr.index.has_duplicates:
        return expr.copy()

    if method == "sum":
        out = expr.groupby(level=0, sort=False).sum()
    elif method == "mean":
        out = expr.groupby(level=0, sort=False).mean()
    elif method == "max_mean":
        order = expr.mean(axis=1).to_numpy()
        tmp = expr.assign(__order=order)
        tmp = tmp.sort_values("__order", ascending=False)
        out = tmp[~tmp.index.duplicated(keep="first")].drop(columns="__order")
        out = out.loc[[g for g in dict.fromkeys(expr.index) if g in out.index]]
    else:
        raise ValueError(f"unknown method={method!r}")
    return out


def match_report(
    expr_index: Sequence[str],
    gene_sets: dict[str, Sequence[str]],
) -> pd.DataFrame:
    """How much of each signature is actually present in the matrix."""
    present = set(expr_index)
    rows = []
    for name, members in gene_sets.items():
        members = list(dict.fromkeys(members))
        found = [g for g in members if g in present]
        rows.append(
            {
                "signature": name,
                "n_genes": len(members),
                "n_found": len(found),
                "frac_found": len(found) / len(members) if members else np.nan,
                "missing": ",".join(g for g in members if g not in present),
            }
        )
    return pd.DataFrame(rows).sort_values("frac_found")
