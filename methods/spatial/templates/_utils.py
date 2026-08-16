"""Shared helpers for the spatial methods templates.
空间方法模板的共享工具函数。

Nothing here computes or fabricates study results; these are I/O, config, and
gene-panel utilities used by the numbered templates.
此处不计算、不伪造研究结果；仅为编号模板提供 I/O、配置与基因面板工具。
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable

import yaml


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_CONFIG = REPO_ROOT / "methods" / "spatial" / "config" / "config.yaml"


def load_config(path: str | os.PathLike | None = None) -> dict:
    """Load the playbook YAML config. 读取手册 YAML 配置。"""
    cfg_path = Path(path) if path else DEFAULT_CONFIG
    with open(cfg_path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh)
    cfg["_config_path"] = str(cfg_path)
    return cfg


def resolve(path_like: str) -> Path:
    """Resolve a config path relative to the repo root if not absolute."""
    p = Path(path_like)
    return p if p.is_absolute() else (REPO_ROOT / p)


def ensure_dir(path: str | os.PathLike) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def assert_panel(var_names: Iterable[str], cfg: dict, platform: str) -> None:
    """Panel gate. 面板门控。

    Whole-transcriptome platforms (visium, geomx_wta) always pass. Targeted
    imaging panels are only allowed through if BOTH required genes are present.
    全转录组平台恒通过；靶向面板仅当两个必需基因都存在时通过。
    """
    required = set(cfg["panel_gate"]["required_genes"])
    present = set(var_names)
    missing = sorted(required - present)

    wta_like = platform in {"visium", "geomx_wta"}
    if not missing:
        return
    if wta_like:
        # A WTA/Visium object missing these genes is a data problem, warn loudly.
        raise ValueError(
            f"[panel_gate] {platform} object is missing {missing}. "
            "For whole-transcriptome data these genes must be present; check the "
            "gene id namespace (symbol vs Ensembl) before proceeding."
        )
    if cfg["panel_gate"]["abort_if_missing_on_targeted_panel"]:
        raise ValueError(
            f"[panel_gate] Targeted panel is missing {missing}. This playbook only "
            "analyses TACSTD2/CLDN4 niches when BOTH are on the panel. Aborting; "
            "do not substitute proxy genes."
        )


def genes_in(adata_var_names: Iterable[str], genes: Iterable[str]) -> list[str]:
    """Subset a marker list to those actually present. 仅保留实际存在的标记。"""
    present = set(adata_var_names)
    return [g for g in genes if g in present]
