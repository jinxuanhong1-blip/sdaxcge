/**
 * methods/mif/qupath/02_export_cell_table.groovy
 *
 * Export detections (cells) to a TSV that analyze_cell_table.py can read.
 * 将检测（细胞）导出为 analyze_cell_table.py 可读的 TSV。
 *
 * Columns: cell_id, image_id, x, y, phenotype, then all QuPath measurements
 * (nucleus / cytoplasm / membrane / cell means if computed).
 *
 * This does not create data. It only writes what is already in the current image.
 * 不创造数据，只写出当前图像里已有的内容。
 *
 * No images are included in this repository.
 */
import qupath.lib.objects.PathCellObject

def imageData = getCurrentImageData()
if (imageData == null) {
    println "No image open. 未打开图像。"
    return
}

def cells = getDetectionObjects().findAll { it instanceof PathCellObject }
if (cells.isEmpty()) {
    println "No PathCellObject detections. Run 01_stardist_cells.groovy first."
    println "没有细胞检测。请先运行 01_stardist_cells.groovy。"
    return
}

def defaultName = getProjectEntry() != null ?
        getProjectEntry().getImageName() + '.cells.tsv' :
        'cells.tsv'

def outFile = getQuPath().getDialogHelper().promptToSaveFile(
        "Export cell table / 导出细胞表",
        null,
        defaultName,
        "Tab-separated",
        ".tsv"
)
if (outFile == null) {
    println "Export cancelled. 已取消导出。"
    return
}

def mlNames = new LinkedHashSet<String>()
cells.each { c ->
    c.getMeasurementList().getMeasurementNames().each { mlNames.add(it) }
}

def img = imageData.getServer().getShortServerName()
outFile.withWriter { w ->
    def header = ['cell_id', 'image_id', 'x', 'y', 'phenotype'] + mlNames.toList()
    w.writeLine(header.join('\t'))
    cells.each { c ->
        def roi = c.getROI()
        def row = [
                c.getID().toString(),
                img,
                roi.getCentroidX(),
                roi.getCentroidY(),
                c.getPathClass() != null ? c.getPathClass().toString() : ''
        ]
        def ml = c.getMeasurementList()
        mlNames.each { name ->
            def v = ml.getMeasurementValue(name)
            row.add(v == null || Double.isNaN((double) v) ? '' : v)
        }
        w.writeLine(row.join('\t'))
    }
}

println "Wrote ${cells.size()} cells → ${outFile.getAbsolutePath()}"
println "Copy the TSV to methods/mif/data/cells.csv (or pass --table)."
println "将 TSV 复制到 methods/mif/data/cells.csv（或传 --table）。"
