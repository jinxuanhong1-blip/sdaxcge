# Cell-table column dictionary / 细胞表列字典

One row = one cell. Spatial metrics are computed **per `image_id`**.

一行 = 一个细胞。空间指标按 **`image_id`** 计算。

## Identity / 标识

| Canonical name | Aliases (case-insensitive) | Required | Notes |
|---|---|---|---|
| `cell_id` | `cell`, `object_id`, `cellid`, `name` | yes | Unique within `image_id`. |
| `image_id` | `roi_id`, `roi`, `image`, `filename`, `slide_id` | yes | One spatial graph per value. |
| `sample_id` | `patient_id`, `case_id`, `donor` | no | For later clinical join only. |

## Coordinates / 坐标

| Canonical name | Aliases | Required | Notes |
|---|---|---|---|
| `x` | `centroid_x`, `cell: centroid x µm`, `x_centroid`, `x_um` | for spatial | Prefer µm. |
| `y` | `centroid_y`, `cell: centroid y µm`, `y_centroid`, `y_um` | for spatial | Prefer µm. |

If coordinates are pixels, pass `--pixel-size-um` to `analyze_cell_table.py`.

若坐标是像素，运行时传 `--pixel-size-um`。

## Marker intensities / 标记强度

Any of these may be missing. The script uses what exists.

下列列均可缺失；脚本用存在的。

| Canonical name | Aliases (contains, case-insensitive) | Typical QuPath export |
|---|---|---|
| `trop2` | `trop2`, `tacstd2` | `Cell: TROP2 mean` |
| `trop2_nucleus` | `trop2` + `nucleus` / `nuclear` | `Nucleus: TROP2 mean` |
| `trop2_cytoplasm` | `trop2` + `cytoplasm` / `cyto` | `Cytoplasm: TROP2 mean` |
| `trop2_membrane` | `trop2` + `membrane` / `mem` | `Membrane: TROP2 mean` |
| `cldn4` | `cldn4`, `claudin4`, `claudin-4`, `claudin 4` | `Cell: CLDN4 mean` |
| `cldn4_nucleus` | `cldn4`/`claudin` + `nucleus` | `Nucleus: CLDN4 mean` |
| `cldn4_cytoplasm` | `cldn4`/`claudin` + `cytoplasm` | `Cytoplasm: CLDN4 mean` |
| `cldn4_membrane` | `cldn4`/`claudin` + `membrane` | `Membrane: CLDN4 mean` |
| `panck` | `panck`, `pan-ck`, `pan_ck`, `cytokeratin`, `ae1` | `Cell: PanCK mean` |
| `cd8` | `cd8` | `Cell: CD8 mean` |
| `pd_l1` | `pd-l1`, `pdl1`, `cd274` | `Cell: PD-L1 mean` |
| `dapi` | `dapi`, `dna`, `ir193`, `ir191` | `Nucleus: DAPI mean` |

QuPath measurement names often look like `Cell: TROP2 mean`. The mapper treats them as intensities if they contain the marker token and a summary (`mean`, `median`). Max / std-dev columns are ignored for phenotyping.

QuPath 测量名常为 `Cell: TROP2 mean`。若含标记词与 `mean`/`median`，映射为强度。`max` / `std` 不用于表型。

## Phenotype / 表型

| Canonical name | Aliases | Notes |
|---|---|---|
| `phenotype` | `class`, `classification`, `celltype`, `cell_type` | Free text. Tokens `TROP2`, `CLDN4`, `PanCK`/`CK`/`tumor`, `CD8` are parsed if present (QuPath combined classes like `PanCK: TROP2: CLDN4`). |

## QC / 质控

| Canonical name | Aliases | Notes |
|---|---|---|
| `area_um2` | `cell: area`, `area` | Drop debris/clumps. |
| `nucleus_area` | `nucleus: area` | Optional. |

## Header-only example / 仅表头示例

See `cell_table.header.csv`. It is empty of rows on purpose.

该文件故意没有数据行。
