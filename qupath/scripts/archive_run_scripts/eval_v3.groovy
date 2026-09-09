/* args: "<classifierName>,<outCsv>"
 *
 * Measures BOTH candidate endpoints in one pass, so the "does diffuse haze count?"
 * question can be answered from the same run rather than re-running:
 *
 *   total Abeta area                -> "all immunoreactivity"
 *   object-size histogram           -> any size-based definition, decided later
 *
 * Exporting a histogram rather than every object because a section can contain
 * >200k connected components at this resolution, most of them single noise pixels.
 * MIN_OBJ_UM2 prunes those before they are ever created.
 *
 * Ignore* area is recovered as (full annotation area - Abeta - Negative), because
 * QuPath omits Ignore* from its own measurements by design.
 */
import qupath.lib.objects.PathObject

def parts = (args[0] as String).split(',')
def C = parts[0]
def OUT = parts[1]
double MIN_OBJ_UM2 = 20.0
def EDGES = [20, 80, 500, 2000, 5000, 20000, 100000] as double[]

def name = getProjectEntry().getImageName()
def cal = getCurrentServer().getPixelCalibration()
double pxArea = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()

removeAllObjects()
createFullImageAnnotation(true)
selectAnnotations()
def classifier = loadPixelClassifier(C)
addPixelClassifierMeasurements(classifier, 'M')
def ann = getAnnotationObjects()[0]
def ml = ann.getMeasurementList()
def g = { n -> ml.getMeasurementNames().find{ it.contains(n) }?.with{ ml.get(it) } ?: 0d }
double abeta = g('Abeta area'), negative = g('Negative area'), pct = g('Abeta %')
double totalUm2 = ann.getROI().getArea() * pxArea

// object-level decomposition of the Abeta class
selectAnnotations()
createDetectionsFromPixelClassifier(classifier, MIN_OBJ_UM2, 0.0, "SPLIT")
def dets = getDetectionObjects().findAll { it.getPathClass()?.toString() == 'Abeta' }
long[] counts = new long[EDGES.length + 1]
double[] areas = new double[EDGES.length + 1]
dets.each { d ->
    double a = d.getROI().getArea() * pxArea
    int b = EDGES.length
    for (int i = 0; i < EDGES.length; i++) { if (a < EDGES[i]) { b = i; break } }
    counts[b]++; areas[b] += a
}
def sb = new StringBuilder()
sb.append(String.format("%s,%s,%.6f,%.1f,%.1f,%.1f,%d", name, C, pct, abeta, negative, totalUm2, dets.size()))
for (int i = 0; i < counts.length; i++) sb.append(String.format(",%d,%.1f", counts[i], areas[i]))
new File(OUT).append(sb.toString() + System.lineSeparator())
removeAllObjects()
