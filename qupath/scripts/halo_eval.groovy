/* Blob-recovery experiment, stage 2. Run per synthetic image, no project needed:
 *
 *   QuPath script --image <host>_blobs.ome.tif --args "<classifier.json>,<blobs.json>,<out.csv>,<label>" halo_eval.groovy
 *
 * For every injected blob: Abeta area inside its footprint disc (radius 1.665 sigma) and inside a
 * 2x disc, plus Abeta area over the whole tissue-bearing frame. Recovery = footprint Abeta / footprint area.
 * Also run on the matching _clean image to get the host's own false-positive baseline.
 */
import qupath.lib.classifiers.pixel.PixelClassifiers
import qupath.lib.objects.PathObjects
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import groovy.json.JsonSlurper

def parts = (args[0] as String).split(',')
def CLS = parts[0]; def BLOBS = parts[1]; def OUT = parts[2]; def LABEL = parts[3]
setImageType('FLUORESCENCE')
def classifier = PixelClassifiers.readClassifier(java.nio.file.Paths.get(CLS))
def meta = new JsonSlurper().parse(new File(BLOBS))
def cal = getCurrentServer().getPixelCalibration()
double pxArea = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()
def plane = ImagePlane.getDefaultPlane()
removeAllObjects()
def objs = []
meta.blobs.eachWithIndex { b, i ->
    [1.0d, 2.0d].each { k ->
        double r = (b.foot_px as double) * k
        def roi = ROIs.createEllipseROI((b.cx as double) - r, (b.cy as double) - r, 2 * r, 2 * r, plane)
        def ann = PathObjects.createAnnotationObject(roi); ann.setName("blob${i}_x${k as int}"); objs << ann
    }
}
createFullImageAnnotation(true)
def frame = getAnnotationObjects()[0]; frame.setName('frame')
addObjects(objs)
selectObjects(getAnnotationObjects())
addPixelClassifierMeasurements(classifier, 'H')
def meas = { obj, key -> def ml = obj.getMeasurementList(); def n = ml.getMeasurementNames().find { it.contains(key) }; n == null ? 0d : ml.get(n) }
def f = new File(OUT)
if (!f.exists()) f.append("image,classifier,blob,scale,fwhm_um,amp_factor,local_bg,foot_area_um2,abeta_um2,negative_um2,frame_abeta_um2\n")
double frameAb = meas(frame, 'Abeta area')
getAnnotationObjects().findAll { it.getName()?.startsWith('blob') }.each { ann ->
    def (bi, k) = ann.getName().replace('blob','').split('_x')
    def b = meta.blobs[bi as int]
    f.append(String.format("%s,%s,%s,%s,%s,%s,%.1f,%.1f,%.1f,%.1f,%.1f%n", meta.host + (LABEL ? "_" + LABEL : ""), new File(CLS).name.replace('.json',''),
        bi, k, b.fwhm_um, b.amp_factor, b.local_bg as double, ann.getROI().getArea() * pxArea,
        meas(ann, 'Abeta area'), meas(ann, 'Negative area'), frameAb))
}
println "DONE ${meta.host} ${new File(CLS).name}"
