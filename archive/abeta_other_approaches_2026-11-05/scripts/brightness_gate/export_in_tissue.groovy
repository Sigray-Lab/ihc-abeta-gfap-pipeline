/* args: "<classifier>,<outGeojson>,<roiDir>" -- export the classifier's Abeta objects INSIDE the tissue ROI (for images whose full frame has unreadable tiles) */
import qupath.lib.io.PathIO
def parts=(args[0] as String).split(','); def C=parts[0]; def OUT=parts[1]; def ROIDIR=parts[2]
def name=getProjectEntry().getImageName()
removeAllObjects()
def objs=PathIO.readObjects(new File("${ROIDIR}/${name}.geojson")); addObjects(objs); selectObjects(objs)
createDetectionsFromPixelClassifier(loadPixelClassifier(C), 0.0, 0.0)
def dets=getDetectionObjects().findAll{ it.getPathClass()?.toString()=='Abeta' }
exportObjectsToGeoJson(dets, OUT, "FEATURE_COLLECTION")
println "EXPORTED ${name} ${C}: ${dets.size()} Abeta objects inside tissue ROI"
removeAllObjects()
