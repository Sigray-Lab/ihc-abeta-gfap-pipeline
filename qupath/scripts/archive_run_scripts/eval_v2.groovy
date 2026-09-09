/* args: "<classifierName>,<outCsv>"
 * Captures Abeta and Negative area AND the full annotation area, so Ignore* area
 * can be derived. QuPath omits Ignore* from its own measurements by design, which
 * is exactly why it has to be recovered here: an over-eager Ignore* silently
 * deletes signal from the numerator and nothing else would reveal it.
 */
def parts = (args[0] as String).split(',')
def C = parts[0]
def OUT = parts[1]
def name = getProjectEntry().getImageName()
removeAllObjects()
createFullImageAnnotation(true)
selectAnnotations()
addPixelClassifierMeasurements(loadPixelClassifier(C), 'M')
def ann = getAnnotationObjects()[0]
def ml = ann.getMeasurementList()
if (System.getenv('DUMP_NAMES') == '1')
    println "MEASUREMENT NAMES: ${ml.getMeasurementNames()}"
def g = { n -> ml.getMeasurementNames().find{ it.contains(n) }?.with{ ml.get(it) } ?: 0d }
// full annotation area in um^2, from the ROI itself rather than from any class sum
def cal = getCurrentServer().getPixelCalibration()
double totalUm2 = ann.getROI().getArea() * cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()
new File(OUT).append(String.format("%s,%s,%.6f,%.1f,%.1f,%.1f%n",
   name, C, g('Abeta %'), g('Abeta area'), g('Negative area'), totalUm2))
removeAllObjects()
