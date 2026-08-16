/**
 * methods/mif/qupath/01_stardist_cells.groovy
 *
 * StarDist nucleus detection + cell expansion for multiplex IF or IMC.
 * Requires the QuPath StarDist extension and a local .pb model.
 *
 * 多重 IF 或 IMC 的 StarDist 核检测 + 细胞膨胀。
 * 需要 QuPath StarDist 扩展和本地 .pb 模型.
 *
 * mIF:  dsb2018_heavy_augment.pb on channel DAPI, pixelSize ~0.5
 * IMC:  same model often works on Ir193 if you set pixelSize ~1.0 and inspect
 *
 * Cite: Schmidt et al. MICCAI 2018; Bankhead et al. Sci Rep 2017.
 *
 * No images are included in this repository.
 * 本仓库不含图像。
 *
 * Usage: select parent annotation(s) (tumor), then Run.
 * 用法：先选中肿瘤父标注，再 Run。
 */
import qupath.ext.stardist.StarDist2D

// --- edit these ---
def pathModel = '/PATH/TO/dsb2018_heavy_augment.pb'
def detectChannel = 'DAPI'   // IMC: 'Ir193'
double pixelSizeUm = 0.5     // IMC: try 1.0
double probThreshold = 0.5
double cellExpansionUm = 3.0
double minNucleusArea = 8.0  // µm^2; drop dust
double maxNucleusArea = 400.0
// ------------------

def imageData = getCurrentImageData()
def parents = getSelectedObjects()
if (parents.isEmpty()) {
    println "Select a parent annotation first (tumor ROI). 请先选择肿瘤父标注。"
    return
}

def stardist = StarDist2D.builder(pathModel)
        .threshold(probThreshold)
        .channels(detectChannel)
        .normalizePercentiles(1, 99)
        .pixelSize(pixelSizeUm)
        .cellExpansion(cellExpansionUm)
        .cellConstrainScale(1.5)
        .measureShape()
        .measureIntensity()
        .includeProbability(true)
        .build()

stardist.detectObjects(imageData, parents)
stardist.close()

// Drop impossible nuclei (area in µm^2 if calibration exists).
def toRemove = getDetectionObjects().findAll { d ->
    def area = d.getMeasurementList().getMeasurementValue('Nucleus: Area')
    if (area == null || Double.isNaN(area))
        return false
    return area < minNucleusArea || area > maxNucleusArea
}
if (!toRemove.isEmpty()) {
    removeObjects(toRemove, true)
    println "Removed ${toRemove.size()} detections outside nucleus area [${minNucleusArea}, ${maxNucleusArea}] µm^2."
}

println "StarDist done. Review membrane TROP2/CLDN4 overlays before any threshold."
println "StarDist 完成。阈值前请目视检查膜 TROP2/CLDN4。"
