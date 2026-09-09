/* Import the fixed validation boxes for one image into the blinded project.
 *   QuPath script -p <project> -i <image> --save --args "<boxes_dir>" import_reference_boxes.groovy
 * Boxes come from scripts/make_reference_boxes.py (class RefBox). Existing RefBox annotations
 * on the image are replaced; everything else is left alone.
 */
import qupath.lib.io.PathIO
def dir = args[0] as String
def name = getProjectEntry().getImageName()
def f = new File("${dir}/${name}.geojson")
if (!f.exists()) { println "no boxes for ${name}"; return }
def old = getAnnotationObjects().findAll { it.getPathClass()?.toString() == 'RefBox' }
removeObjects(old, true)
def objs = PathIO.readObjects(f)
addObjects(objs)
println "imported ${objs.size()} RefBox annotations into ${name}"
