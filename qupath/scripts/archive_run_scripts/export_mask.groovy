def C = args[0]
removeAllObjects(); createFullImageAnnotation(true); selectAnnotations()
createDetectionsFromPixelClassifier(loadPixelClassifier(C), 25.0, 0.0)
def dets = getDetectionObjects().findAll { it.getPathClass()?.toString() == 'Abeta' }
def name = getProjectEntry().getImageName()
exportObjectsToGeoJson(dets, "/tmp/masks/${name}__${C}.geojson", "FEATURE_COLLECTION")
println "EXPORTED ${name} ${C}: ${dets.size()} objects"
removeAllObjects()
