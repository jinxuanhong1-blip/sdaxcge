/**
 * methods/mif/qupath/00_set_channels.groovy
 *
 * Rename fluorescence / IMC channels to short marker names.
 * Edit CHANNELS to match YOUR panel. Bessede 5-plex had no CLDN4.
 *
 * 将荧光 / IMC 通道改为短标记名。按你的面板改 CHANNELS。
 * Bessede 五色没有 CLDN4。
 *
 * This repo ships no images. Run only on locally authorized OME-TIFFs.
 * 本仓库不含图像。仅在本地授权 OME-TIFF 上运行。
 *
 * QuPath: Automate → Show script editor → Run.
 */
import qupath.lib.gui.scripting.QPEx

// Example TROP2/CLDN4 panel (not Bessede). Length must match channel count.
// 示例 TROP2/CLDN4 面板（不是 Bessede）。长度必须等于通道数。
def CHANNELS = [
    'DAPI',    // IMC: use 'Ir193' or 'DNA1'
    'TROP2',
    'CLDN4',
    'PanCK',
    'CD8',
    // 'PDL1',
]

def server = getCurrentServer()
int n = server.nChannels()
if (CHANNELS.size() != n) {
    println "Refusing to rename: script has ${CHANNELS.size()} names, image has ${n} channels."
    println "拒绝改名：脚本 ${CHANNELS.size()} 个名字，图像 ${n} 个通道。"
    println "Current names / 当前名:"
    (0..<n).each { i -> println "  ${i}: ${server.getChannel(i).getName()}" }
    return
}

setChannelNames(*CHANNELS as String[])
setImageType('FLUORESCENCE')
println "Channel names set. Image type = Fluorescence (also correct for IMC)."
println "通道名已设置。图像类型 = Fluorescence（IMC 也用这个）。"
