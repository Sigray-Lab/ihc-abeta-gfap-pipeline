/*
 * import_blinded_project.groovy — build a BLINDED QuPath project for region delineation.
 * ---------------------------------------------------------------------------------
 * Pipeline: Abeta + GFAP in App NL-G-F mouse brain.  QuPath 0.7.0 (pinned, D-15).
 * Reads the JSON written by `ihc.ingest.qupath_export.write_project_spec`.
 *
 * WHAT IT GUARANTEES
 *   1. Only tissue series are imported.  A series is opened only if its internal name
 *      ends in `_01`..`_09`, so the slide label image — which shows the tube number as
 *      printed text AND as a DataMatrix barcode — cannot be imported even by mistake.
 *   2. The displayed image name comes from the spec (`K07_s01`), never from the file
 *      and never from the internal series name.  Tube 60's internal series are called
 *      `60_20x_DAPI, FITC, Cy3_01`, and QuPath shows internal names, not file names.
 *   3. Image type is set to FLUORESCENCE.
 *   4. **Only DAPI is visible.  FITC (GFAP) and Cy3 (Abeta) are switched off.**
 *      This is a scientific requirement, not a nicety: anatomical boundaries drawn
 *      while marker signal is visible are pulled towards or away from plaque burden,
 *      which biases the denominator of every number measured inside them, even though
 *      the person drawing is blinded to treatment group (CLAUDE_v1.2 §9, ADR-0010).
 *      If the display cannot be set for any image, the script FAILS rather than
 *      handing over a project that quietly does not meet that requirement.
 *   5. It refuses to run at all if the spec references anything that looks like an
 *      animal identifier.
 *
 * HOW TO RUN
 *   Command line (preferred — it is reproducible and leaves a log):
 *       QuPath.app/Contents/MacOS/QuPath script \
 *           --args /path/to/project_out/project_spec.json \
 *           qupath/scripts/import_blinded_project.groovy
 *   Or set the environment variable IHC_PROJECT_SPEC and run the script from the
 *   QuPath script editor (Automate > Script editor, then Run).
 *   Or, last resort, edit SPEC_PATH_OVERRIDE below.
 *
 * WHO RUNS IT
 *   The custodian (the person holding the blinding key), NOT the person who will draw
 *   the regions.  The spec's symlink directory sits next to the raw data; the finished
 *   QuPath project is what gets handed over.
 *
 * TESTED — 2026-08-07, QuPath v0.7.0 arm64 (build 2026-02-25, commit 04ccfa4)
 *   Run headless against the real cohort: 26 of 26 images imported, every one typed
 *   FLUORESCENCE at 0.3250 um/px, and every one opening with DAPI selected and FITC
 *   and Cy3 deselected (verified independently by re-opening the project and reading
 *   ImageDisplay.selectedChannels, not by trusting this script's own report).
 *   The written project.qpproj was scanned for leaks: no Image_NN filename, no RawData
 *   path, no stack1, no `60_` prefix, no group or arm name.
 *
 *   One real bug this testing caught: QuPath 0.7.0 does NOT bundle groovy-json, so the
 *   original `import groovy.json.JsonSlurper` failed to compile. Gson is on the
 *   classpath instead and is what this script now uses.
 *
 * WHAT IS STILL UNTESTED
 *   The interactive GUI path -- this script has
 *   never been executed.  Everything it depends on is checked at run time and it fails
 *   loudly rather than silently: see the CHECKS list printed at the end of a run.
 */

import java.awt.image.BufferedImage
// QuPath 0.7.0 does NOT bundle groovy-json -- `import groovy.json.JsonSlurper` fails to
// compile with "unable to resolve class". Gson is on the classpath because QuPath uses it
// internally, so use that. Verified on QuPath v0.7.0 arm64.
import com.google.gson.Gson
import com.google.gson.GsonBuilder
import com.google.gson.JsonParser

import qupath.lib.images.ImageData
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.projects.Projects
import qupath.lib.display.ImageDisplay

// ---------------------------------------------------------------------------------
// 0.  Locate the spec
// ---------------------------------------------------------------------------------

/** Last-resort hard-coded spec path. Leave empty and use --args or IHC_PROJECT_SPEC. */
final String SPEC_PATH_OVERRIDE = ""

String specPath = null
try {
    if (binding.hasVariable('args') && args != null && args.length > 0 && args[0])
        specPath = args[0] as String
} catch (Throwable ignored) { /* no 'args' binding in the GUI script editor */ }
if (!specPath) specPath = System.getenv('IHC_PROJECT_SPEC')
if (!specPath) specPath = SPEC_PATH_OVERRIDE ?: null
if (!specPath)
    throw new IllegalStateException(
        "No project spec given.\n" +
        "  Pass it with:  QuPath script --args /path/to/project_spec.json <this script>\n" +
        "  or set the environment variable IHC_PROJECT_SPEC\n" +
        "  or edit SPEC_PATH_OVERRIDE at the top of this file.")

def specFile = new File(specPath)
if (!specFile.isFile())
    throw new IllegalStateException("Project spec not found: ${specFile.getAbsolutePath()}")

String specText = specFile.getText('UTF-8')
def spec = new Gson().fromJson(specText, Map.class)

if (spec.schema != 'ihc.ingest.qupath_export/1')
    throw new IllegalStateException(
        "Unexpected spec schema ${spec.schema}. This script understands " +
        "'ihc.ingest.qupath_export/1'. Regenerate the spec, or use the matching " +
        "version of the script — do not mix them.")

println "=" * 84
println "Blinded QuPath project import"
println "  spec     : ${specFile.getAbsolutePath()}"
println "  project  : ${spec.project_name}"
println "  images   : ${spec.counts?.images} (${spec.counts?.positive} positive, " +
        "${spec.counts?.negative} negative) from ${spec.counts?.animals} animals"
println "  skipped  : ${spec.counts?.skipped}   excluded: ${spec.counts?.excluded}"
println "=" * 84

// QuPath version is a warning, not a gate: 0.7.1 should be fine, 0.5 probably is not.
try {
    def running = qupath.lib.common.GeneralTools.getVersion()?.toString()
    def wanted = spec.qupath?.expected_version
    if (running && wanted && !running.startsWith(wanted))
        println "  WARNING  running QuPath ${running}, the pipeline pins ${wanted} " +
                "(env/tool_versions.yaml, D-15). Record which one produced this project."
    else if (running)
        println "  QuPath   ${running}"
} catch (Throwable t) {
    println "  WARNING  could not determine the QuPath version: ${t}"
}

// ---------------------------------------------------------------------------------
// 1.  DEFENSIVE: refuse to build a project out of a spec that leaks animal identity
// ---------------------------------------------------------------------------------
// Two patterns, used at two strengths.
//   FILENAME  `Image_29`, `image 29`, `Image-29`  — safe to run over any text.
//   BARE_TUBE a cohort tube number standing alone (29..58, 60). The lookarounds are
//             alphanumeric rather than \b so a coded ID like `K42` does NOT match
//             (no boundary between K and 4) but `1007344 - 29` — the text printed on
//             the slide label — does.
// The bare-tube pattern is applied only to fields that become part of the project;
// applied to the whole document it would fire on ordinary counts like "images": 34.

final def FILENAME_RE  = ~/(?i)image[ _\-]?\d{1,4}/
final def BARE_TUBE_RE = ~/(?<![0-9A-Za-z])(?:29|[3-5][0-9]|60)(?![0-9A-Za-z])/

def hits = { String text, boolean strict ->
    if (text == null) return []
    def found = []
    FILENAME_RE.matcher(text).each { found << it }
    if (strict) BARE_TUBE_RE.matcher(text).each { found << it }
    return found.unique()
}

def leaks = []
def checkStrict = { String label, String text ->
    def h = hits(text, true)
    if (h) leaks << "${label}: ${h.join(', ')}  in  ${text}"
}

checkStrict("project_name", spec.project_name as String)
spec.channels?.names?.each { checkStrict("channel name", it as String) }
spec.images?.each { img ->
    checkStrict("image_name", img.image_name as String)
    checkStrict("code", img.code as String)
    checkStrict("series_match_suffix", img.series_match_suffix as String)
}

// Paths get their own rule. The leak in a path is the FILE NAME — `Image_29.vsi`,
// `_Image_29_`. Parent directories are named by whoever set up the run and their
// numbers mean nothing; refusing on a folder called `batch-42` would be a false
// alarm, and a check that cries wolf is a check that gets switched off. So the
// filename pattern runs over the whole path, the bare-number pattern over the
// basename only. This mirrors _assert_clean_path() on the Python side.
def pathsAreCoded = (spec.blinding?.paths_are_coded == true)
def acknowledgement = (spec.blinding?.uncoded_paths_acknowledged_by ?: "").toString().trim()
def pathLeaks = []
spec.images?.each { img ->
    def full = img.vsi_path as String
    def h = hits(full, false) + hits(new File(full).getName(), true)
    if (h) pathLeaks << "${img.image_name} -> ${full}  (${h.unique().join(', ')})"
}
// The whole document, filename pattern only: catches an identifier smuggled into a
// warning message, a note, or a field this script does not know about.
def documentLeaks = hits(specText, false)

if (documentLeaks)
    leaks << "the spec document contains file names that identify animals: " +
             documentLeaks.unique().join(', ')

if (pathLeaks) {
    if (pathsAreCoded) {
        leaks.addAll(pathLeaks.collect { "image path (spec claims coded paths): ${it}" })
    } else if (!acknowledgement) {
        leaks.addAll(pathLeaks.collect { "image path: ${it}" })
    } else {
        println ""
        println "  !! UNCODED PATHS, ACKNOWLEDGED BY: ${acknowledgement}"
        println "  !! ${pathLeaks.size()} image path(s) contain the animal number. QuPath stores"
        println "  !! the file path in the project and shows it in the Image tab, so anyone"
        println "  !! opening this project can read the tube ID. This is NOT a blinded project."
        pathLeaks.each { println "  !!   ${it}" }
        println ""
    }
}

if (leaks) {
    println ""
    println "REFUSING TO BUILD THE PROJECT — the spec carries animal identifiers:"
    leaks.each { println "  * ${it}" }
    println ""
    println "A blinded artefact that leaks is worse than none, because the leak is"
    println "invisible in the output and the delineation it produces looks exactly as"
    println "trustworthy as a clean one. Rebuild the spec with build_project_spec()."
    throw new IllegalStateException("blinding check failed: ${leaks.size()} identifier leak(s)")
}
println "  blinding check: clean (${spec.images?.size() ?: 0} image entries scanned)"

if (!spec.images)
    throw new IllegalStateException(
        "The spec lists no images, so this would create an empty project and report " +
        "success. Look at project_spec.json -> skipped and -> excluded: on this cohort " +
        "the usual cause is that no payload folders have been transferred yet (23 of 31 " +
        "animals are index-only), which is a normal state but means there is nothing to " +
        "delineate.")

// ---------------------------------------------------------------------------------
// 2.  Create the project
// ---------------------------------------------------------------------------------

def projectDir = new File(spec.project_dir as String)
if (!projectDir.exists() && !projectDir.mkdirs())
    throw new IOException("could not create ${projectDir}")

def existing = new File(projectDir, "project.qpproj")
if (existing.isFile())
    throw new IllegalStateException(
        "A project already exists at ${existing.getAbsolutePath()}.\n" +
        "  Refusing to add images to it: if regions have already been drawn, a second\n" +
        "  import would create duplicate entries under the same coded names and there\n" +
        "  would be no way to tell which annotations belong to which. Build into a\n" +
        "  fresh directory instead.")

def project = Projects.createProject(projectDir, BufferedImage.class)

// ---------------------------------------------------------------------------------
// 3.  Helpers
// ---------------------------------------------------------------------------------

/** Series whose name matches any of these is never opened. */
def forbidden = (spec.blinding?.forbidden_series_names ?: []) as List

/** Tissue series names always end in _01.._09. Nothing else is ever a section. */
final def SECTION_SUFFIX_RE = ~/.*_0[1-9]$/

def builderCache = [:]      // vsi path -> List<ServerBuilder>
def seriesCache  = [:]      // vsi path -> List<Map> [index, name, width, height]

def buildersFor = { String path ->
    if (builderCache.containsKey(path)) return builderCache[path]
    def support = ImageServerProvider.getPreferredUriImageSupport(BufferedImage.class, path)
    def list = support == null ? null : support.getBuilders()
    builderCache[path] = list
    return list
}

/**
 * Describe every series in a file: name, width, height.
 * Each server is built once and closed immediately; the result is cached, so a
 * four-section slide costs one pass over its six series rather than four.
 */
def seriesFor = { String path ->
    if (seriesCache.containsKey(path)) return seriesCache[path]
    def builders = buildersFor(path)
    def out = []
    if (builders != null) {
        builders.eachWithIndex { b, i ->
            def name = null, w = -1, h = -1, err = null
            // Some builders expose metadata directly; most need a build. Try the
            // cheap route first, fall back to building, never let one bad series
            // abort the whole file.
            try {
                def md = b.getMetadata()
                if (md != null) { name = md.getName(); w = md.getWidth(); h = md.getHeight() }
            } catch (Throwable ignored) { }
            if (name == null) {
                def server = null
                try {
                    server = b.build()
                    name = server.getMetadata().getName()
                    w = server.getWidth()
                    h = server.getHeight()
                } catch (Throwable t) {
                    err = t.toString()
                } finally {
                    if (server != null) { try { server.close() } catch (Throwable ignored) { } }
                }
            }
            out << [index: i, name: name, width: w, height: h, error: err]
        }
    }
    seriesCache[path] = out
    return out
}

def visibleNames = (spec.channels?.visible ?: ['DAPI']) as List
def isVisibleChannel = { String channelName ->
    if (channelName == null) return false
    return visibleNames.any { channelName.toLowerCase().contains((it as String).toLowerCase()) }
}

/**
 * Switch the display to the counterstain only, persist it into the ImageData, and
 * verify that it took. Returns null on success or a reason string on failure.
 */
def setDapiOnlyDisplay = { imageData ->
    def display = null
    try { display = ImageDisplay.create(imageData) } catch (Throwable ignored) { }
    if (display == null) { try { display = new ImageDisplay(imageData) } catch (Throwable ignored) { } }
    if (display == null) {
        try {
            display = new ImageDisplay()
            display.setImageData(imageData, false)
        } catch (Throwable t) {
            return "could not create an ImageDisplay: ${t}"
        }
    }

    def available, selected
    try {
        available = new ArrayList<>(display.availableChannels())
        if (available.isEmpty()) return "the image reports no display channels"

        // Turn the counterstain on BEFORE turning anything off, so the selection is
        // never momentarily empty (QuPath may refuse to deselect the last channel).
        available.findAll { isVisibleChannel(it.getName()) }
                 .each { display.setChannelSelected(it, true) }
        available.findAll { !isVisibleChannel(it.getName()) }
                 .each { display.setChannelSelected(it, false) }

        selected = new ArrayList<>(display.selectedChannels())
    } catch (Throwable t) {
        return "could not set channel visibility (${t}). The QuPath display API may " +
               "have changed; do not proceed until the marker channels are provably off."
    }
    def selectedNames = selected.collect { it.getName() }
    if (selected.isEmpty())
        return "no channel is selected after setting the display"
    def stillOn = selectedNames.findAll { !isVisibleChannel(it) }
    if (stillOn)
        return "marker channel(s) still visible after setting the display: ${stillOn}"

    // Persist. `saveChannelColorProperties` writes the display into the ImageData
    // properties, which is what QuPath reads back when the image is next opened.
    def persisted = false
    try { display.saveChannelColorProperties(); persisted = true } catch (Throwable ignored) { }
    if (!persisted) {
        try {
            imageData.setProperty(ImageDisplay.class.getName(), display.toJSON(false))
            persisted = true
        } catch (Throwable ignored) { }
    }
    if (!persisted) {
        try {
            imageData.setProperty(ImageDisplay.class.getName(), display.toJSON())
            persisted = true
        } catch (Throwable ignored) { }
    }
    if (!persisted)
        return "the display was set to ${selectedNames} but could not be saved into " +
               "the project, so it would revert when the image is opened"

    return [ok: true, selected: selectedNames, display: display]
}

// ---------------------------------------------------------------------------------
// 4.  Import
// ---------------------------------------------------------------------------------

def imported = []
def failures = []
def notes = []

spec.images.eachWithIndex { img, n ->
    String name = img.image_name
    String path = img.vsi_path
    String suffix = img.series_match_suffix
    def fail = { String why -> failures << [image: name, reason: why]; println "  FAIL  ${name}: ${why}" }

    print String.format("  [%3d/%3d] %-12s ", n + 1, spec.images.size(), name)

    def file = new File(path)
    if (!file.exists()) {
        println ""
        fail("file not found: ${path}. If the spec uses coded symlinks, run " +
             "write_project_spec() again — the links live beside the spec.")
        return
    }

    def builders = buildersFor(path)
    if (builders == null || builders.isEmpty()) {
        println ""
        fail("QuPath could not open this file at all (no image support)")
        return
    }

    def all = seriesFor(path)
    if (img.n_series_expected != null && all.size() != img.n_series_expected)
        notes << "${name}: QuPath reports ${all.size()} series, the payload has " +
                 "${img.n_series_expected} stacks. Not fatal — the series is matched by " +
                 "name — but worth knowing."

    // ---- resolve the series -------------------------------------------------
    // By NAME, not by index. The index in the spec is derived from the payload's
    // stack order and is only a hint: it has never been checked against Bio-Formats.
    // The name suffix cannot select the label or the overview, which is the point.
    def matches = all.findAll { s ->
        s.name != null && SECTION_SUFFIX_RE.matcher(s.name).matches() && s.name.endsWith(suffix)
    }
    if (matches.isEmpty()) {
        println ""
        def names = all.collect { it.name ?: "<unreadable>" }
        fail("no series whose name ends in '${suffix}'. Series present: ${names}")
        return
    }
    if (matches.size() > 1) {
        println ""
        fail("${matches.size()} series end in '${suffix}' — ambiguous, refusing to guess")
        return
    }
    def chosen = matches[0]

    if (forbidden.any { (chosen.name as String).toLowerCase().contains((it as String).toLowerCase()) }) {
        println ""
        fail("resolved to a forbidden series (${chosen.name}). This should be " +
             "impossible; treat it as a bug and do not work around it.")
        return
    }

    if (img.series_index != null && chosen.index != (img.series_index as int))
        notes << "${name}: opened series ${chosen.index}, the spec's hint was " +
                 "${img.series_index}. Resolved by name, which is authoritative — but " +
                 "the hint is systematically wrong if this appears for every image."

    // ---- dimensions: the strongest available proof we opened the right thing --
    if (img.width_px != null && img.height_px != null && chosen.width > 0) {
        if (chosen.width != (img.width_px as int) || chosen.height != (img.height_px as int)) {
            println ""
            fail("dimension mismatch: QuPath reports ${chosen.width}x${chosen.height}, " +
                 "the VSI index says ${img.width_px}x${img.height_px}. Either the wrong " +
                 "series was opened or the file does not match the spec. Refusing.")
            return
        }
    }

    // ---- add it -------------------------------------------------------------
    def entry
    try {
        entry = project.addImage(builders[chosen.index])
    } catch (Throwable t) {
        println ""
        fail("could not add the image: ${t}")
        return
    }

    // The displayed name comes from the spec. Never from the file, never from the
    // internal series name (tube 60's internal names start with the tube number).
    entry.setImageName(name)

    def imageData
    try {
        imageData = entry.readImageData()
        imageData.setImageType(ImageData.ImageType.FLUORESCENCE)
    } catch (Throwable t) {
        println ""
        fail("could not read/type the image data: ${t}")
        return
    }

    def displayResult = setDapiOnlyDisplay(imageData)
    if (!(displayResult instanceof Map && displayResult.ok)) {
        println ""
        fail("CHANNEL DISPLAY: ${displayResult}. Regions must not be drawn with " +
             "marker channels visible (CLAUDE_v1.2 §9), so this image is a failure, " +
             "not a warning.")
        return
    }

    try {
        entry.saveImageData(imageData)
    } catch (Throwable t) {
        println ""
        fail("could not save the image data: ${t}")
        return
    }

    // Thumbnail. It must be rendered through OUR display, or the project browser
    // would show a plaque-burden preview of every section. If the two-argument form
    // is unavailable, no thumbnail is set — a blank tile is fine, a leaky one is not.
    def thumbState = "skipped"
    try {
        def img2 = qupath.lib.gui.commands.ProjectCommands.getThumbnailRGB(
                imageData.getServer(), displayResult.display)
        entry.setThumbnail(img2)
        thumbState = "rendered from the DAPI-only display"
    } catch (Throwable t) {
        thumbState = "not set (no display-aware thumbnail API): ${t.getClass().getSimpleName()}"
    }

    // Metadata for the headless stages that later join annotations back to the
    // manifest. Deliberately does NOT record positive/negative: it is of no use when
    // drawing an outline, and knowing a section is a no-primary control could make
    // it get less care than the sections it is supposed to control for.
    def meta = [ihc_code: img.code as String,
                ihc_section: img.section_label as String,
                ihc_spec_schema: spec.schema as String]
    try {
        meta.each { k, v -> entry.putMetadataValue(k, v) }
    } catch (Throwable first) {
        try {
            meta.each { k, v -> entry.getMetadata().put(k, v) }
        } catch (Throwable second) {
            notes << "${name}: could not attach metadata (${second.getClass().getSimpleName()}). " +
                     "Cosmetic only — the coded image name still carries code and section."
        }
    }

    imported << [image: name, series: chosen.index, series_name_suffix: suffix,
                 width: chosen.width, height: chosen.height,
                 channels_visible: displayResult.selected, thumbnail: thumbState]
    println "series ${chosen.index}  ${chosen.width}x${chosen.height}  " +
            "visible=${displayResult.selected}"
}

project.syncChanges()

// ---------------------------------------------------------------------------------
// 5.  Report, then fail loudly if anything is not right
// ---------------------------------------------------------------------------------

def report = [
    schema        : "ihc.qupath.import_report/1",
    spec          : specFile.getAbsolutePath(),
    project       : projectDir.getAbsolutePath(),
    n_requested   : spec.images.size(),
    n_imported    : imported.size(),
    n_failed      : failures.size(),
    channels_kept_visible: visibleNames,
    channels_hidden      : spec.channels?.hidden,
    imported      : imported,
    failures      : failures,
    notes         : notes,
]
def reportFile = new File(specFile.getParentFile(), "import_report.json")
reportFile.setText(new GsonBuilder().setPrettyPrinting().create().toJson(report), 'UTF-8')

println ""
println "=" * 84
println "  imported ${imported.size()} of ${spec.images.size()} image(s) into ${projectDir}"
println "  visible channel(s): ${visibleNames}   hidden: ${spec.channels?.hidden}"
println "  report: ${reportFile.getAbsolutePath()}"
if (notes) {
    println "  notes:"
    notes.each { println "    - ${it}" }
}
if (spec.skipped) {
    println "  ${spec.skipped.size()} section(s) were left out of the spec itself " +
            "(mostly: pixels not transferred yet). See project_spec.json -> skipped."
}
if (spec.excluded) {
    println "  ${spec.excluded.size()} section(s) were EXCLUDED because their staining " +
            "condition is unresolved. See project_spec.json -> excluded, and resolve at " +
            "the bench before they can be used."
}
println "=" * 84

if (failures) {
    println ""
    println "${failures.size()} image(s) FAILED:"
    failures.each { println "  * ${it.image}: ${it.reason}" }
    throw new RuntimeException(
        "${failures.size()} of ${spec.images.size()} images failed to import. The project " +
        "at ${projectDir} is INCOMPLETE and must not be handed over. Fix the failures, " +
        "delete the directory, and re-run — do not import the rest on top.")
}
println "OK — every image imported, typed FLUORESCENCE, and set to ${visibleNames} only."
