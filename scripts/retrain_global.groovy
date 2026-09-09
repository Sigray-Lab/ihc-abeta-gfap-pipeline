/*
 * retrain_global.groovy — third arm: ONE offset and ONE gain per section.
 *
 *   QuPath script -p <project> -i <any one image> --args "<outName>,<umPerPx>,<annoDir>,<mapFile>,<anchorCsv>,<outDir>"
 *   (-i is required: without it QuPath runs a project script once PER IMAGE, i.e. 121 trainings)
 *
 * Why this exists. Exposure is constant across the endpoint set (1840 ms on all 67
 * positives) yet Cy3 tissue medians span 4.1x, so dividing by exposure corrects nothing.
 * The unnormalised classifier tracks that residual brightness (r = +0.73 in antibody-free
 * controls); z130 removes it but, being a local-contrast filter, cannot represent signal
 * broader than its own neighbourhood. This arm sits between them: rescale each section
 * once, globally, and leave every spatial frequency intact.
 *
 * The map is  I' = (I - glass) / (tissue_p50 - glass), so glass -> 0 and typical
 * non-plaque neuropil -> 1 in every section. Both anchors come from
 * an anchor CSV with columns image, glass_p50, t_p50 (see docs/normalisation.md).
 *
 * Why the tissue anchor is not circular here: plaque burden is ~2% of tissue area, so the
 * median of in-tissue Cy3 is plaque-free by construction and barely moves if burden
 * doubles. That is the difference between this and whole-histogram matching, which forces
 * the high tail -- the endpoint itself -- to agree between sections.
 *
 * Implementation note. QuPath bakes one op chain into a classifier, so a per-image
 * constant cannot live in a shared chain. Training therefore builds a SEPARATE op per
 * training image and concatenates the resulting sample matrices before fitting one forest.
 * Application uses `make_global_classifiers.py`, which writes one small classifier JSON
 * per section with that section's constants prepended. No image is ever rewritten.
 */

import qupath.lib.images.servers.ColorTransforms
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.io.PathIO
import qupath.opencv.ml.OpenCVClassifiers
import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.opencv.ops.ImageOps
import qupath.opencv.tools.MultiscaleFeatures.MultiscaleFeature
import qupath.process.gui.commands.ml.PixelClassifierTraining
import qupath.lib.classifiers.pixel.PixelClassifierMetadata
import org.bytedeco.opencv.opencv_ml.RTrees
import org.bytedeco.opencv.opencv_ml.TrainData
import org.bytedeco.opencv.opencv_core.Mat
import static org.bytedeco.opencv.global.opencv_core.vconcat
import static org.bytedeco.opencv.global.opencv_ml.ROW_SAMPLE
import org.bytedeco.javacpp.indexer.IntIndexer
import static org.bytedeco.opencv.global.opencv_core.CV_32S
import static org.bytedeco.opencv.global.opencv_core.CV_32F
import org.bytedeco.javacpp.indexer.FloatIndexer
import qupath.lib.objects.classes.PathClass

def parts   = (args[0] as String).split(',')
String OUT  = parts[0]
double UM   = parts[1] as double
String ANNO = parts[2]
String MAPF = parts[3]
String ANCH = parts[4]
String ODIR = parts[5]

def MAP = [:]
new File(MAPF).eachLine { l ->
    if (l.startsWith('#') || !l.contains(',')) return
    def b = l.split(','); MAP[b[0].trim()] = b[1].trim()
}
// image -> [glass, gain]; gain is tissue_p50 - glass, the multiplicative term
def ANCHOR = [:]
def hdr = null
new File(ANCH).eachLine { l ->
    def b = l.split(',')
    if (hdr == null) { hdr = b as List; return }
    def m = [:]; hdr.eachWithIndex { h, i -> m[h] = b[i] }
    double g = m['glass_p50'] as double
    ANCHOR[m['image']] = [g, (m['t_p50'] as double) - g]
}
println "  anchors: ${ANCHOR.size()} sections"

def feats = [MultiscaleFeature.GAUSSIAN, MultiscaleFeature.LAPLACIAN,
             MultiscaleFeature.WEIGHTED_STD_DEV, MultiscaleFeature.GRADIENT_MAGNITUDE]
// Feature sigmas are declared in MICROMETRES and converted to pixels at the working
// resolution, so that a resolution sweep changes sampling WITHOUT changing physical
// feature scale. Until 2026-08-19 these were hard-coded pixel sigmas [1, 2, 4]; run that
// way against the config's downsample sweep [1,2,4,8] (0.325-2.60 um/px), pixel size and
// physical scale move together and the sweep answers nothing. External review 2026-08-19.
def FEATURE_SIGMAS_UM = [1.3, 2.6, 5.2]      // == [1, 2, 4] px at the current 1.30 um/px
def scaleOps = FEATURE_SIGMAS_UM.collect { um ->
    def s = um / UM
    ImageOps.Filters.features(feats, s, s)
}
println("  feature sigmas: " + FEATURE_SIGMAS_UM + " um = " +
        FEATURE_SIGMAS_UM.collect { (it / UM).round(3) } + " px at ${UM} um/px")

def opFor = { double glass, double gain ->
    ImageOps.buildImageDataOp(ColorTransforms.createChannelExtractor('Cy3'))
            .appendOps([ImageOps.Core.subtract(glass),
                        ImageOps.Core.divide(gain),
                        ImageOps.Core.splitMerge(scaleOps)] as qupath.opencv.ops.ImageOp[])
}

def cal = new PixelCalibration.Builder().pixelSizeMicrons(UM, UM).build()
def project = getProject()
def byName = [:]; project.getImageList().each { byName[it.getImageName()] = it }

Mat allX = null, allY = null
def labelMap = null
// One canonical class->index map for every image. Order matches the shipped model so
// downstream code reading class indices is unaffected.
def CANON = ['Abeta': 0, 'Ignore*': 1, 'Negative': 2]
def PC_BY_NAME = [:]
int nImg = 0
new File(ANNO).listFiles().sort{ it.name }.each { f ->
    if (!f.name.endsWith('_annotations.geojson')) return
    def src = f.name.replace('_annotations.geojson','')
    def target = MAP[src]
    if (target == null || byName[target] == null) { println "  SKIP ${src}"; return }
    def a = ANCHOR[target]
    if (a == null) { println "  SKIP ${src}: no anchor for ${target}"; return }
    def entry = byName[target]
    def imageData = entry.readImageData()
    imageData.getHierarchy().clearAll()
    imageData.getHierarchy().addObjects(PathIO.readObjects(f))
    def training = new PixelClassifierTraining(opFor(a[0], a[1]))
    training.setResolution(cal)
    def data = training.createTrainingData([imageData])
    // A file containing a single class yields no training data: QuPath cannot form a
    // contrast from one label. The other arms pooled all images before training so their
    // single-class files still contributed Ignore* samples; this arm trains per image, so
    // those two files (34_s02, 49_s01 -- Ignore* only) are necessarily dropped. The
    // 2026-08-18 ablation showed removing them changes the repeat error by less than one
    // standard error, so the arms remain comparable, but it IS a difference and is logged.
    if (data == null) { println "  SKIP ${src}: single-class file, no contrast to learn"; return }
    def td = data.getTrainData()
    // ---- LABEL-MAP FIX (external review 2026-09-09, finding F1) ------------------------
    // createTrainingData() is called PER IMAGE, and QuPath numbers the classes it finds in
    // THAT image from 0. A three-class image yields {Abeta:0, Ignore*:1, Negative:2}; a
    // control image with only Ignore*/Negative yields {Ignore*:0, Negative:1}. The old code
    // kept the first image's map and concatenated raw responses, so a control's Ignore*
    // became Abeta and its Negative became Ignore* -- 1,857 of 4,782 samples (38.8 %) in
    // the shipped v3_global. Every image's responses are now remapped to ONE canonical
    // map before concatenation, and the remap count is logged per image.
    def lm = data.getLabelMap()                       // PathClass -> Integer, this image
    int[] remap = new int[lm.size()]
    lm.each { pc, ix ->
        def name = pc.toString()
        if (!CANON.containsKey(name))
            throw new RuntimeException("unexpected class '${name}' in ${src}; canonical map is ${CANON}")
        remap[ix as Integer] = CANON[name]
        PC_BY_NAME[name] = pc
    }
    def X = td.getSamples(), Y = td.getResponses()
    int changed = 0
    if (Y.depth() == CV_32S) {
        IntIndexer yi = Y.createIndexer()
        for (long i = 0; i < Y.rows(); i++) {
            int v = yi.get(i, 0L); int nv = remap[v]
            if (nv != v) { yi.put(i, 0L, nv); changed++ }
        }
        yi.release()
    } else if (Y.depth() == CV_32F) {
        FloatIndexer yf = Y.createIndexer()
        for (long i = 0; i < Y.rows(); i++) {
            int v = Math.round(yf.get(i, 0L)); int nv = remap[v]
            if (nv != v) { yf.put(i, 0L, (float) nv); changed++ }
        }
        yf.release()
    } else {
        throw new RuntimeException("unexpected response depth ${Y.depth()} in ${src}")
    }
    println "  LABELS ${src}: image map=${lm} -> canonical ${CANON}; remapped ${changed} of ${Y.rows()} samples"
    labelMap = CANON
    if (allX == null) { allX = X; allY = Y }
    else {
        def nx = new Mat(); vconcat(allX, X, nx); allX = nx
        def ny = new Mat(); vconcat(allY, Y, ny); allY = ny
    }
    nImg++
    println "  ${src} -> ${target}  glass=${String.format('%.1f',a[0])} gain=${String.format('%.1f',a[1])}  samples=${X.rows()}"
}
if (allX == null) { println 'NO TRAINING DATA'; return }
println "  pooled: ${allX.rows()} samples from ${nImg} images, ${allX.cols()} features"

def td = TrainData.create(allX, ROW_SAMPLE, allY)
def statModel = OpenCVClassifiers.createStatModel(RTrees.class)
statModel.train(td)

def byLabel = [:]
CANON.each { name, ix -> byLabel[ix as Integer] = PC_BY_NAME[name] ?: PathClass.fromString(name) }
def metadata = new PixelClassifierMetadata.Builder()
        .inputResolution(cal).inputShape(512, 512)
        .setChannelType(qupath.lib.images.servers.ImageServerMetadata.ChannelType.CLASSIFICATION)
        .classificationLabels(byLabel).build()
// reference chain uses identity constants; per-section constants are injected at apply time
def classifier = PixelClassifiers.createClassifier(statModel, opFor(0.0d, 1.0d), metadata, true)
def outPath = java.nio.file.Paths.get(ODIR, OUT + '.json')
PixelClassifiers.writeClassifier(classifier, outPath)
println "WROTE ${outPath}"
