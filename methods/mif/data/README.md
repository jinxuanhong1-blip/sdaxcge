# Cell tables go here / 细胞表放这里

**English.** Drop a user-supplied cell table as `cells.csv` or `cells.parquet`. Nothing in this folder is a public Bessede or TROP2/CLDN4 image dump. Bessede multiplex IF images and cell-level tables are not public (Clin Cancer Res 2024;30:779–785, data-availability statement). This directory is gitignored except this README.

**中文。** 将使用者提供的细胞表存为 `cells.csv` 或 `cells.parquet`。此文件夹不是公开的 Bessede 或 TROP2/CLDN4 图像库。Bessede 多重免疫荧光图像与单细胞表不公开。除本 README 外，目录内容被 gitignore。

Required columns are documented in `../schema/cell_table_columns.md`.

Then:

```bash
python ../python/analyze_cell_table.py
```
