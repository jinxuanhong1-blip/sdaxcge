"""bulkimmune: bulk RNA immune deconvolution and exclusion scoring for ICI cohorts.

The package is deliberately dependency-light (numpy/pandas/scipy/sklearn/statsmodels)
so that the whole pipeline runs without an R installation. Where a method is only
defined by an R package or a web server (xCell spillover data, CIBERSORTx, TIP),
the algorithm is reimplemented here against the *published* parameters and the
original resource files are downloaded by ``scripts/00_fetch_resources.py`` rather
than vendored into the repository.

See ``playbook.md`` for the methodological guidance and ``docs/VALIDATION.md`` for
the numerical agreement with the reference R implementations.
"""

__version__ = "1.0.0"

__all__ = [
    "genes",
    "preprocess",
    "ssgsea",
    "signatures",
    "mcpcounter",
    "estimate",
    "xcell",
    "cibersort",
    "tide",
    "tip",
    "batch",
    "stats",
    "geo",
    "tcga",
    "pipeline",
    "concordance",
]
