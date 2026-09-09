def C = args[0]
removeAllObjects(); createFullImageAnnotation(true); selectAnnotations()
addPixelClassifierMeasurements(loadPixelClassifier(C), 'M')
def ml = getAnnotationObjects()[0].getMeasurementList()
def g = { n -> ml.getMeasurementNames().find{ it.contains(n) }?.with{ ml.get(it) } ?: 0d }
new File('/tmp/full_results.csv').append(String.format("%s,%s,%.6f,%.1f,%.1f%n",
   getProjectEntry().getImageName(), C, g('Abeta %'), g('Abeta area'), g('Negative area')))
removeAllObjects()
