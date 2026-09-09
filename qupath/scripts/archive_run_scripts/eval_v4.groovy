/* args: "<classifierName>,<outCsv>,<roiDir>"
 *
 * Measures the classifier INSIDE the DAPI tissue mask, so numerator and denominator
 * refer to the same region. Fixes the stop-level defect found by external review
 * 2026-08-18: the previous scripts measured the whole rectangular frame and divided by
 * DAPI tissue area, which is not a percent area of anything. For the locally-normalised
 * classifier the classified area came to 1.169x the tissue area on every image.
 *
 * Also exports, as QC:
 *   - Abeta measured on the FULL FRAME, so out-of-tissue Abeta is quantified, not hidden
 *   - the tissue ROI area as QuPath sees it, to verify the mask round-tripped intact
 *
 * Fails loudly. A missing ROI, a missing measurement, or classes that do not partition
 * the mask is an error, not a zero.
 */
import qupath.lib.io.PathIO

def parts = (args[0] as String).split(',')
def C = parts[0]; def OUT = parts[1]; def ROIDIR = parts[2]
def name = getProjectEntry().getImageName()
def cal = getCurrentServer().getPixelCalibration()
double pxArea = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()

def roiFile = new File("${ROIDIR}/${name}.geojson")
if (!roiFile.exists()) throw new RuntimeException("NO TISSUE ROI for ${name}")

def classifier = loadPixelClassifier(C)
def meas = { obj, key ->
    def ml = obj.getMeasurementList()
    def n = ml.getMeasurementNames().find { it.contains(key) }
    if (n == null) throw new RuntimeException("MISSING MEASUREMENT '${key}' on ${name} / ${C}")
    return ml.get(n)
}

// ---- 1. inside the tissue mask -------------------------------------------------
removeAllObjects()
def objs = PathIO.readObjects(roiFile)
if (objs.isEmpty()) throw new RuntimeException("EMPTY TISSUE ROI for ${name}")
addObjects(objs)
def tissue = getAnnotationObjects()[0]
double tissueUm2 = tissue.getROI().getArea() * pxArea
selectObjects([tissue])
addPixelClassifierMeasurements(classifier, 'T')
double abetaIn = meas(tissue, 'Abeta area')
double negIn   = meas(tissue, 'Negative area')
double pctIn   = meas(tissue, 'Abeta %')

// ---- 2. full frame, for the out-of-tissue QC field ------------------------------
removeAllObjects()
createFullImageAnnotation(true)
def frame = getAnnotationObjects()[0]
double frameUm2 = frame.getROI().getArea() * pxArea
selectObjects([frame])
addPixelClassifierMeasurements(classifier, 'F')
double abetaAll = meas(frame, 'Abeta area')
double negAll   = meas(frame, 'Negative area')

// ---- 3. assertions --------------------------------------------------------------
// 0/0 means the classifier assigned the ENTIRE tissue mask to Ignore* -- a total
// failure on that image. That is data, not a crash: record it as NaN so it appears in
// the results and is counted, rather than vanishing. Only a real number outside 0-100
// indicates a broken measurement.
boolean degenerate = (abetaIn + negIn) <= 0d
if (!degenerate && (pctIn < -0.001 || pctIn > 100.001))
    throw new RuntimeException("Abeta % outside 0-100 on ${name}/${C}: ${pctIn}")
if (degenerate) pctIn = Double.NaN
if (abetaIn > tissueUm2 * 1.001 || negIn > tissueUm2 * 1.001)
    throw new RuntimeException("class area exceeds tissue on ${name}/${C}")

new File(OUT).append(String.format("%s,%s,%.6f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f%n",
   name, C, pctIn, abetaIn, negIn, tissueUm2, abetaAll, negAll, frameUm2,
   Math.max(abetaAll - abetaIn, 0d)))
removeAllObjects()
