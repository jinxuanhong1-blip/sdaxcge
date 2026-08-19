#!/usr/bin/env python3
"""Download official 10x Xenium human lung / NSCLC gene panels and test CLDN4.

Gate (user instruction): confirm CLDN4 is on the lung panel (base or custom
add-on). If absent, report panel genes and stop — no cell-level neighborhood
analysis, no substitute genes, no GSE300007, no private 8-KL.
"""

from __future__ import annotations

import csv
import json
import urllib.request
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "data" / "cache" / "xenium_panels"
OUT = ROOT / "results" / "xenium_official_lung"

UA = "Mozilla/5.0 (compatible; xenium-lung-panel-check/1.0)"

# Official 10x.com/datasets human lung / NSCLC / FFPE lung Xenium releases.
# gene_panel.json is the authoritative per-run target list (base + add-on).
DATASETS = [
    {
        "dataset_id": "xenium_human_lung_preview_cancer",
        "title": "Xenium Human Lung Preview Data — invasive adenocarcinoma",
        "page": "https://www.10xgenomics.com/datasets/xenium-human-lung-preview-data-1-standard",
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/1.3.0/"
            "Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE/"
            "Xenium_Preview_Human_Lung_Cancer_With_Add_on_2_FFPE_gene_panel.json"
        ),
        "tissue": "FFPE invasive lung adenocarcinoma (Avaden)",
        "role": "lung_preview_addon",
    },
    {
        "dataset_id": "xenium_human_lung_preview_nondiseased",
        "title": "Xenium Human Lung Preview Data — non-diseased lung",
        "page": "https://www.10xgenomics.com/datasets/xenium-human-lung-preview-data-1-standard",
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/1.3.0/"
            "Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE/"
            "Xenium_Preview_Human_Non_diseased_Lung_With_Add_on_FFPE_gene_panel.json"
        ),
        "tissue": "FFPE healthy adult lung (Avaden)",
        "role": "lung_preview_addon",
    },
    {
        "dataset_id": "xenium_hlung_cancer_multitissue_preview",
        "title": "Human Lung Cancer Preview Data (Multi-Tissue and Cancer Panel)",
        "page": (
            "https://www.10xgenomics.com/datasets/"
            "human-lung-cancer-preview-data-xenium-human-multi-tissue-and-cancer-panel-1-standard"
        ),
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/1.5.0/"
            "Xenium_V1_hLung_cancer_section/"
            "Xenium_V1_hLung_cancer_section_gene_panel.json"
        ),
        "tissue": "FFPE lung cancer microstrip (AcePix)",
        "role": "multitissue_377",
    },
    {
        "dataset_id": "xenium_ffpe_luad_multimodal",
        "title": "Preview Data: FFPE Human Lung Cancer with Multimodal Cell Segmentation",
        "page": (
            "https://www.10xgenomics.com/datasets/"
            "preview-data-ffpe-human-lung-cancer-with-xenium-multimodal-cell-segmentation-1-standard"
        ),
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/2.0.0/"
            "Xenium_V1_humanLung_Cancer_FFPE/"
            "Xenium_V1_humanLung_Cancer_FFPE_gene_panel.json"
        ),
        "tissue": "FFPE lung adenocarcinoma (BioIVT)",
        "role": "multitissue_377",
    },
    {
        "dataset_id": "xenium_ffpe_nsclc_io_addon",
        "title": "FFPE Human Lung Cancer with Human Immuno-Oncology Profiling Panel and Custom Add-on",
        "page": (
            "https://www.10xgenomics.com/datasets/"
            "ffpe-human-lung-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard"
        ),
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/2.0.0/"
            "Xenium_V1_Human_Lung_Cancer_Addon_FFPE/"
            "Xenium_V1_Human_Lung_Cancer_Addon_FFPE_gene_panel.json"
        ),
        "tissue": "FFPE NSCLC Stage I-B Grade 2 (Discovery Life Sciences)",
        "role": "io_addon",
    },
    {
        "dataset_id": "xenium_postxenium_v1_lung_panel",
        "title": "Post-Xenium Technical Note — Experiment 1 (Xenium v1 lung panel)",
        "page": "https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote",
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/3.0.0/"
            "Xenium_V1_Human_Lung_Cancer_FFPE/"
            "Xenium_V1_Human_Lung_Cancer_FFPE_gene_panel.json"
        ),
        "tissue": "FFPE lung adenocarcinoma (invasive acinar)",
        "role": "lung_base_official",
    },
    {
        "dataset_id": "xenium_postxenium_prime5k",
        "title": "Post-Xenium Technical Note — Experiment 2 (Xenium Prime 5K)",
        "page": "https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote",
        "panel_url": (
            "https://cf.10xgenomics.com/samples/xenium/3.0.0/"
            "Xenium_Prime_Human_Lung_Cancer_FFPE/"
            "Xenium_Prime_Human_Lung_Cancer_FFPE_gene_panel.json"
        ),
        "tissue": "FFPE lung adenocarcinoma (same donor series as Exp. 1)",
        "role": "prime5k",
    },
]

OFFICIAL_LUNG_PANEL_JSON = (
    "https://cdn.10xgenomics.com/raw/upload/v1697489942/"
    "software-support/Xenium-panels/gene_panel_json_files/"
    "human_lung_gene_expression.json"
)

CLDN4_ENSEMBL = "ENSG00000189143"

# Requested analysis genes + close family / controls.
MARKERS = [
    "CLDN4",
    "CLDN1",
    "CLDN3",
    "CLDN5",
    "CLDN7",
    "CLDN18",
    "CD8A",
    "CD8B",
    "CD3D",
    "CD3E",
    "CD3G",
    "KRT8",
    "KRT18",
    "KRT7",
    "KRT5",
    "EPCAM",
    "TACSTD2",
    "PTPRC",
    "NKX2-1",
]


def download(url: str, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as resp:
        dest.write_bytes(resp.read())
    return dest


def parse_targets(obj: dict) -> list[dict]:
    rows = []
    for t in obj.get("payload", {}).get("targets", []):
        typ = t.get("type") or {}
        data = typ.get("data") or {}
        src = (t.get("source") or {}).get("identity") or {}
        rows.append(
            {
                "name": data.get("name") or "",
                "ensembl": data.get("id") or "",
                "descriptor": typ.get("descriptor") or "",
                "source_name": src.get("name") or "",
                "design_id": src.get("design_id") or "",
                "category": (t.get("source") or {}).get("category") or "",
            }
        )
    return rows


def gene_names(rows: list[dict], category: str | None = None) -> list[str]:
    out = []
    for r in rows:
        if r["descriptor"] != "gene" or not r["name"]:
            continue
        if category and r["category"] != category:
            continue
        out.append(r["name"])
    return sorted(set(out))


def has_cldn4(rows: list[dict]) -> bool:
    for r in rows:
        if r["descriptor"] != "gene":
            continue
        if r["name"] == "CLDN4" or r["ensembl"] == CLDN4_ENSEMBL:
            return True
    return False


def write_gene_list(path: Path, genes: list[str]) -> None:
    path.write_text("\n".join(genes) + ("\n" if genes else ""))


def svg_presence_grid(
    path: Path,
    columns: list[str],
    rows: list[str],
    present: dict[tuple[str, str], bool],
    title: str,
    col_labels: dict[str, str] | None = None,
) -> None:
    """Simple presence heatmap (no matplotlib)."""
    cell_w, cell_h = 46, 24
    left, top = 130, 170
    width = max(720, left + cell_w * len(columns) + 40)
    height = top + cell_h * len(rows) + 40
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}" font-family="Helvetica, Arial, sans-serif">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="16" y="28" font-size="15" font-weight="700" fill="#111">'
        f"{_xml(title)}</text>",
        '<text x="16" y="48" font-size="11" fill="#444">'
        "Green = on panel (gene target). Gray = absent. Official 10x gene_panel.json.</text>",
        '<rect x="16" y="58" width="12" height="12" fill="#1b7f4a" stroke="#0f5132"/>',
        '<text x="32" y="68" font-size="11" fill="#222">present</text>',
        '<rect x="96" y="58" width="12" height="12" fill="#e9ecef" stroke="#adb5bd"/>',
        '<text x="112" y="68" font-size="11" fill="#222">absent</text>',
    ]
    for i, col in enumerate(columns):
        label = (col_labels or {}).get(col, col)
        x = left + i * cell_w + cell_w / 2
        y = top - 8
        parts.append(
            f'<text x="{x:.1f}" y="{y:.1f}" font-size="10" fill="#222" '
            f'text-anchor="end" transform="rotate(-60 {x:.1f} {y:.1f})">{_xml(label)}</text>'
        )
    for j, gene in enumerate(rows):
        y = top + j * cell_h
        parts.append(
            f'<text x="{left - 8}" y="{y + cell_h * 0.7:.1f}" font-size="11" '
            f'text-anchor="end" fill="#111">{_xml(gene)}</text>'
        )
        for i, col in enumerate(columns):
            x = left + i * cell_w
            on = present.get((col, gene), False)
            fill = "#1b7f4a" if on else "#e9ecef"
            stroke = "#0f5132" if on else "#adb5bd"
            parts.append(
                f'<rect x="{x + 2}" y="{y + 2}" width="{cell_w - 4}" height="{cell_h - 4}" '
                f'rx="3" fill="{fill}" stroke="{stroke}"/>'
            )
            mark = "Y" if on else "—"
            tfill = "#fff" if on else "#6c757d"
            parts.append(
                f'<text x="{x + cell_w / 2:.1f}" y="{y + cell_h * 0.68:.1f}" '
                f'font-size="10" text-anchor="middle" fill="{tfill}">{mark}</text>'
            )
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n")


def _xml(s: str) -> str:
    return (
        s.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def main() -> int:
    CACHE.mkdir(parents=True, exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)

    catalog_rows = []
    panel_sets: dict[str, set[str]] = {}
    detailed = []

    # Official predesigned lung panel (support site, not a run).
    official_path = download(OFFICIAL_LUNG_PANEL_JSON, CACHE / "official_human_lung_gene_expression.json")
    official_obj = json.loads(official_path.read_text())
    official_rows = parse_targets(official_obj)
    official_genes = gene_names(official_rows)
    panel_sets["lung_base_hLung_v1.1"] = set(official_genes)
    write_gene_list(OUT / "lung_base_hLung_v1.1_genes.txt", official_genes)
    official_cldn4 = has_cldn4(official_rows)
    panel_meta = official_obj.get("payload", {}).get("panel", {})
    catalog_rows.append(
        {
            "dataset_id": "official_predesigned_lung_panel",
            "title": panel_meta.get("description") or "Xenium Human Lung Gene Expression Panel",
            "page": "https://www.10xgenomics.com/support/software/xenium-panel-designer/latest/tutorials/pre-designed-xenium-v1",
            "tissue": "predesigned panel (not a tissue run)",
            "panel_name": panel_meta.get("identity", {}).get("name", ""),
            "design_id": panel_meta.get("identity", {}).get("design_id", ""),
            "n_gene_targets": len(official_genes),
            "n_cldn4": int(official_cldn4),
            "cldn_family": ",".join(sorted(g for g in official_genes if g.startswith("CLDN"))),
            "has_CD8A": int("CD8A" in official_genes),
            "has_KRT8": int("KRT8" in official_genes),
            "has_EPCAM": int("EPCAM" in official_genes),
            "role": "lung_base_official",
        }
    )

    preview_addon_genes: list[str] | None = None
    io_addon_genes: list[str] | None = None

    for ds in DATASETS:
        dest = CACHE / f"{ds['dataset_id']}_gene_panel.json"
        path = download(ds["panel_url"], dest)
        obj = json.loads(path.read_text())
        rows = parse_targets(obj)
        genes = gene_names(rows)
        panel_sets[ds["dataset_id"]] = set(genes)
        meta = obj.get("payload", {}).get("panel", {})
        cldn4 = has_cldn4(rows)
        cldn_fam = sorted(g for g in genes if g.startswith("CLDN"))
        catalog_rows.append(
            {
                "dataset_id": ds["dataset_id"],
                "title": ds["title"],
                "page": ds["page"],
                "tissue": ds["tissue"],
                "panel_name": meta.get("description") or meta.get("identity", {}).get("name", ""),
                "design_id": meta.get("identity", {}).get("design_id", ""),
                "n_gene_targets": len(genes),
                "n_cldn4": int(cldn4),
                "cldn_family": ",".join(cldn_fam),
                "has_CD8A": int("CD8A" in genes),
                "has_KRT8": int("KRT8" in genes),
                "has_EPCAM": int("EPCAM" in genes),
                "role": ds["role"],
            }
        )
        write_gene_list(OUT / f"{ds['dataset_id']}_genes.txt", genes)

        # Split add-on vs base for the two add-on runs.
        if ds["role"] == "lung_preview_addon" and preview_addon_genes is None:
            preview_addon_genes = gene_names(rows, category="current")
            preview_base = gene_names(rows, category="base")
            write_gene_list(OUT / "lung_preview_addon_PD346_genes.txt", preview_addon_genes)
            write_gene_list(OUT / "lung_preview_base_PD339_genes.txt", preview_base)
            panel_sets["lung_preview_addon_only"] = set(preview_addon_genes)
            panel_sets["lung_preview_base_only"] = set(preview_base)
        if ds["role"] == "io_addon" and io_addon_genes is None:
            io_addon_genes = gene_names(rows, category="current")
            io_base = gene_names(rows, category="base")
            write_gene_list(OUT / "io_addon_A_genes.txt", io_addon_genes)
            write_gene_list(OUT / "io_base_hImmune_v1_genes.txt", io_base)
            panel_sets["io_addon_only"] = set(io_addon_genes)

        src_counts = Counter((r["category"], r["source_name"]) for r in rows if r["descriptor"] == "gene")
        detailed.append(
            {
                "dataset_id": ds["dataset_id"],
                "sources": {f"{a}|{b}": n for (a, b), n in src_counts.items()},
                "cldn4": cldn4,
            }
        )

    # Marker presence table.
    presence_cols = [
        "lung_base_hLung_v1.1",
        "lung_preview_base_only",
        "lung_preview_addon_only",
        "xenium_hlung_cancer_multitissue_preview",
        "xenium_ffpe_nsclc_io_addon",
        "xenium_postxenium_prime5k",
    ]
    col_labels = {
        "lung_base_hLung_v1.1": "Lung base v1.1 (289)",
        "lung_preview_base_only": "Lung preview base (292)",
        "lung_preview_addon_only": "Lung preview +100 add-on",
        "xenium_hlung_cancer_multitissue_preview": "Multi-tissue 377",
        "xenium_ffpe_nsclc_io_addon": "IO 380 + add-on A",
        "xenium_postxenium_prime5k": "Prime 5K (5001)",
    }

    present_map: dict[tuple[str, str], bool] = {}
    with (OUT / "marker_presence.tsv").open("w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["gene"] + [col_labels[c] for c in presence_cols])
        for gene in MARKERS:
            row = [gene]
            for col in presence_cols:
                on = gene in panel_sets.get(col, set())
                present_map[(col, gene)] = on
                row.append("1" if on else "0")
            w.writerow(row)

    svg_presence_grid(
        OUT / "fig_marker_presence.svg",
        presence_cols,
        MARKERS,
        present_map,
        "Requested / control genes on official 10x Xenium lung panels",
        col_labels,
    )

    cldn_rows = ["CLDN1", "CLDN2", "CLDN3", "CLDN4", "CLDN5", "CLDN6", "CLDN7", "CLDN8", "CLDN18"]
    cldn_present: dict[tuple[str, str], bool] = {}
    for gene in cldn_rows:
        for col in presence_cols:
            cldn_present[(col, gene)] = gene in panel_sets.get(col, set())
    svg_presence_grid(
        OUT / "fig_cldn_family.svg",
        presence_cols,
        cldn_rows,
        cldn_present,
        "CLDN family on official 10x Xenium lung panels (CLDN4 is absent)",
        col_labels,
    )

    # Catalog
    with (OUT / "catalog.tsv").open("w", newline="") as f:
        fields = list(catalog_rows[0].keys())
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(catalog_rows)

    lung_base_cldn4 = official_cldn4
    lung_addon_cldn4 = bool(preview_addon_genes) and ("CLDN4" in preview_addon_genes)
    any_run_cldn4 = any(int(r["n_cldn4"]) for r in catalog_rows)

    gate = {
        "cldn4_on_official_lung_base": lung_base_cldn4,
        "cldn4_on_lung_preview_addon": lung_addon_cldn4,
        "cldn4_on_any_official_10x_lung_xenium_run": any_run_cldn4,
        "action": "STOP",
        "reason": (
            "CLDN4 is not a gene target on the Xenium Human Lung Gene Expression "
            "Panel (base) or the official lung preview custom add-on. Per "
            "instructions: report panel genes and stop. No cell-level CLDN4-high "
            "vs CD8 neighborhood analysis was run."
        ),
        "n_official_lung_base_genes": len(official_genes),
        "n_lung_preview_addon_genes": len(preview_addon_genes or []),
        "excluded": ["GSE300007", "private_8KL", "claim-failed framing"],
        "datasets_checked": [r["dataset_id"] for r in catalog_rows],
    }
    (OUT / "gate.json").write_text(json.dumps(gate, indent=2) + "\n")
    (OUT / "parse_notes.json").write_text(json.dumps(detailed, indent=2) + "\n")

    print(json.dumps(gate, indent=2))
    print(f"Wrote outputs under {OUT}")
    if lung_base_cldn4 or lung_addon_cldn4:
        print("CLDN4 PRESENT — downstream cell-level analysis would proceed.")
        return 0
    print("CLDN4 ABSENT on lung base and lung add-on. STOP.")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
