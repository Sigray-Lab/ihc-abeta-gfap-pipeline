/*
 * retrain.groovy — rebuild the Abeta classifier headlessly, with per-image
 * normalisation, from the wet-lab scientist's own training annotations.
 *
 *   QuPath script -p <project.qpproj> --args "<sigma>,<outName>,<umPerPx>" retrain.groovy
 *
 * sigma = local-normalisation Gaussian sigma in PIXELS at the working resolution.
 *         0 disables it, reproducing the delivered classifier as a control.
 *
 * Why local normalisation and not the "Normalization" option in the training
 * dialog: that one fits ONE offset/scale per feature from the pooled training data
 * and applies it to every image identically. The trees are then trained on those
 * normalised features, so their thresholds are already in normalised units --
 * a fixed monotonic per-feature rescale cannot change a decision tree's output at
 * all. It would look like a fix in the JSON and change nothing in the numbers.
 * Local normalisation divides each pixel by variation measured in its own
 * neighbourhood, inside each image, which is what actually removes the
 * section-to-section brightness dependence.
 */

import qupath.lib.images.servers.ColorTransforms
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.io.PathIO
import qupath.lib.objects.classes.PathClass
import qupath.opencv.ml.OpenCVClassifiers
import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.opencv.ops.ImageOps
import qupath.opencv.tools.MultiscaleFeatures.MultiscaleFeature
import qupath.process.gui.commands.ml.PixelClassifierTraining
import qupath.lib.classifiers.pixel.PixelClassifierMetadata
import org.bytedeco.opencv.opencv_ml.RTrees
import static org.bytedeco.opencv.global.opencv_core.setRNGSeed

// args: "<sigmaMean>,<outName>,<umPerPx>,<annotationDir>,<mappingFile>,<outputDir>[,<sigmaVar>[,<stdRadius>,<noiseFloor>[,<seed>]]]"
//
// sigmaMean and sigmaVar are in PIXELS at the working resolution.
//   sigmaVar = 0  -> Gaussian local MEAN SUBTRACTION only. Removes an additive
//                   offset and leaves multiplicative gain untouched: doubling the
//                   input doubles the output. This is what an earlier sweep ran,
//                   and it is why that sweep could not have fixed a gain problem.
//   sigmaVar > 0  -> local mean subtraction AND division by local standard
//                   deviation. Gain-invariant: doubling the input leaves the
//                   output unchanged. Verified directly against
//                   LocalNormalization.gaussianNormalize2D.
def parts = (args[0] as String).split(',')
double SIGMA   = parts[0] as double
String OUT     = parts[1]
double UM      = parts[2] as double
String ANNO_DIR = parts[3]
String MAP_FILE = parts[4]
String OUT_DIR  = parts[5]
double SIGMA_VAR = parts.length > 6 ? (parts[6] as double) : 0.0
// Optional 8th/9th args switch to a NOISE-FLOORED z-score:
//   Z(x) = (I - localMean_sigma) / max(localStd_radius, floor)
// The built-in localNormalization has no floor, so in flat background it divides by
// noise and manufactures structure -- which is exactly how the gain-invariant runs
// lost specificity. clip() supplies the floor.
int    STD_RADIUS = parts.length > 7 ? (parts[7] as int)    : 0
double NOISE_FLOOR = parts.length > 8 ? (parts[8] as double) : 0.0
// MEASURED 2026-08-18: this knob does nothing. RTrees bootstraps samples and
// subsamples features per split, so a seed *should* change the forest -- but training
// the same configuration at seed 1 and seed 2 produced byte-identical classifiers
// (sha256 c4e837287390bf21, 3518448 bytes, both). setRNGSeed runs and is logged, yet
// does not reach the RNG the forest draws from; cv::theRNG() is thread-local and
// training runs on another thread.
//
// Kept because the empirical fact it establishes is worth having: same annotations and
// settings give a bit-for-bit identical classifier every run. Any difference between two
// of our classifiers is therefore caused by the setting under test, not by the draw.
// It also means forest-realisation sensitivity cannot be probed this way -- perturb the
// training set instead (see the *_v1files / *_noIgOnly / *_noCtrl ablations).
// 0 = leave OpenCV's default.
int    SEED = parts.length > 9 ? (parts[9] as int) : 0

// Which annotation file belongs on which image. Read from a mapping file rather
// than hard-coded, because that mapping IS part of the blinding key --seven rows of it
// would put tube -> code pairs into this repository.
//
// Format: one "sourceStem,targetImageName" per line, e.g.  Image_NN_s01,<code>_s01
//
// NOTE the remap is not a rename. A rescan's series _01 is a DIFFERENT physical
// section from the original's _01 (it is physical section 03), so matching by name,
// or by pixel dimensions, silently lands annotations on the wrong tissue -- in this
// cohort, on a negative control. Build the mapping from the manifest, not by eye.
def MAP = [:]
new File(MAP_FILE).eachLine { line ->
    line = line.trim()
    if (!line || line.startsWith('#')) return
    def bits = line.split(',')
    MAP[bits[0].trim()] = bits[1].trim()
}
println "  mapping: ${MAP.size()} annotation file(s)"
println "  normalisation: " + (STD_RADIUS > 0
    ? "FLOORED Z-SCORE  sigmaMean=${SIGMA}px stdRadius=${STD_RADIUS}px floor=${NOISE_FLOOR} (gain-invariant + noise floor)"
    : SIGMA_VAR > 0 ? "sigmaMean=${SIGMA}px sigmaVar=${SIGMA_VAR}px (GAIN-INVARIANT, no floor)"
    : SIGMA > 0 ? "sigmaMean=${SIGMA}px (mean subtraction only)" : "(none)")

def project = getProject()
def byName = [:]
project.getImageList().each{ byName[it.getImageName()] = it }

// ---- 1. load each training image and attach its annotations -----------------
def imageDataList = []
MAP.each { src, target ->
    def entry = byName[target]
    if (entry == null) { println "MISSING image ${target}"; return }
    def f = new File("${ANNO_DIR}/${src}_annotations.geojson")
    if (!f.exists()) { println "MISSING annotations ${f}"; return }
    def imageData = entry.readImageData()
    imageData.getHierarchy().clearAll()
    def objs = PathIO.readObjects(f)
    imageData.getHierarchy().addObjects(objs)
    def counts = objs.groupBy{ it.getPathClass()?.toString() }.collectEntries{ k,v -> [k, v.size()] }
    println "  ${src} -> ${target}: ${objs.size()} objects ${counts}"
    imageDataList << imageData
}
if (imageDataList.isEmpty()) { println 'NO TRAINING DATA'; return }

// ---- 2. feature op: Cy3 -> [local normalisation] -> 12 multiscale features --
def feats = [MultiscaleFeature.GAUSSIAN, MultiscaleFeature.LAPLACIAN,
             MultiscaleFeature.WEIGHTED_STD_DEV, MultiscaleFeature.GRADIENT_MAGNITUDE]
def scaleOps = [1.0, 2.0, 4.0].collect { s -> ImageOps.Filters.features(feats, s, s) }

def ops = []
if (STD_RADIUS > 0) {
    def ident = ImageOps.Core.multiply([1.0d] as double[])
    ops << ImageOps.Core.splitDivide(
        ImageOps.Core.splitSubtract(ident, ImageOps.Filters.gaussianBlur(SIGMA)),
        ImageOps.Core.sequential(ImageOps.Filters.stdDev(STD_RADIUS),
                                 ImageOps.Core.clip(NOISE_FLOOR, 1e9d)))
} else if (SIGMA > 0) {
    ops << ImageOps.Normalize.localNormalization(SIGMA, SIGMA_VAR)
}
ops << ImageOps.Core.splitMerge(scaleOps)

def op = ImageOps.buildImageDataOp(ColorTransforms.createChannelExtractor('Cy3'))
                 .appendOps(ops as qupath.opencv.ops.ImageOp[])

// ---- 3. train ---------------------------------------------------------------
def cal = new PixelCalibration.Builder().pixelSizeMicrons(UM, UM).build()
def training = new PixelClassifierTraining(op)
training.setResolution(cal)

def data = training.createTrainingData(imageDataList)
def trainData = data.getTrainData()
def labels = data.getLabelMap()
println "  training samples: ${trainData.getNSamples()}  features: ${trainData.getNVars()}  labels: ${labels}"

if (SEED != 0) {
    setRNGSeed(SEED)
    println "  RNG seed: ${SEED}"
} else {
    println "  RNG seed: unset (OpenCV default)"
}
def statModel = OpenCVClassifiers.createStatModel(RTrees.class)
statModel.train(trainData)

// ---- 4. assemble and save ---------------------------------------------------
// labelMap is PathClass -> Integer; the metadata wants Integer -> PathClass
def byLabel = [:]
labels.each { pc, ix -> byLabel[ix as Integer] = pc }
def metadata = new PixelClassifierMetadata.Builder()
        .inputResolution(cal)
        .inputShape(512, 512)
        .setChannelType(qupath.lib.images.servers.ImageServerMetadata.ChannelType.CLASSIFICATION)
        .classificationLabels(byLabel)
        .build()

def classifier = PixelClassifiers.createClassifier(statModel, op, metadata, true)
def outPath = java.nio.file.Paths.get(OUT_DIR, OUT + '.json')
PixelClassifiers.writeClassifier(classifier, outPath)
println "WROTE ${outPath}"
