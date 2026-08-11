/*
 * make_threshold_classifier.groovy — build a simple intensity thresholder as a BASELINE.
 * ------------------------------------------------------------------------------------
 * QuPath 0.7.0. Run headless:
 *
 *   QuPath script --args /path/to/output_dir qupath/scripts/make_threshold_classifier.groovy
 *
 * WHAT THIS IS
 *   A single-channel intensity threshold, wrapped as a QuPath pixel classifier so it can
 *   be loaded and applied exactly like a trained one. It says: "Cy3 above N counts is
 *   Abeta, below is Negative." Nothing else — no texture, no shape, no size filter.
 *
 * WHAT IT IS FOR
 *   A baseline and a cross-check, NOT a replacement for the trained classifier. Its value
 *   is that the threshold was chosen objectively rather than by eye: the 52 secondary-only
 *   negative-control sections give a direct false-positive measurement, and the threshold
 *   is the LOWEST one that holds the 95th-percentile negative below 0.05 % area — lowest,
 *   so that as much real signal as possible is retained.
 *
 * WHAT IT CANNOT DO, AND WHY THE TRAINED CLASSIFIER IS STILL NEEDED
 *   A threshold cannot tell a plaque from any other bright thing. At the calibrated value
 *   one negative control still reports 0.23 % "Abeta" — a dust speck or bright edge, which
 *   a threshold has no way to reject. Separation is therefore good but NOT complete:
 *   the weakest positive (0.11 %) sits below the worst negative (0.23 %). Learning to
 *   reject artefacts is exactly what the `Ignore*` class and the trained features add.
 *
 * CALIBRATION (2026-08-11, 53 images at 1.3 um/px, tubes 51/60 excluded as off-exposure)
 *   threshold   positives (median)   negatives (median)   negatives (95th pct)
 *         600              1.749 %             0.009 %               0.897 %
 *         900              0.888 %             0.002 %               0.046 %   <- chosen
 *        2000              0.146 %             0.000 %               0.001 %
 *   Raising it further keeps buying specificity at the cost of real signal; 900 is where
 *   the worst-case negative first drops below 0.05 %.
 *
 *   Sanity check: median 0.89 % Abeta area inside tissue, spread 0.11–5.65 % across
 *   animals. The earlier whole-image method reported a median of 13.1 %.
 *
 * REMEMBER the background is essentially zero here (5th percentile 0–35 counts across all
 * images), which is why a single raw threshold transfers between animals at all. That
 * would not hold on a stain with a variable background.
 */

import qupath.opencv.ml.pixel.PixelClassifiers
import qupath.lib.images.servers.PixelCalibration
import qupath.lib.objects.classes.PathClass
import qupath.lib.classifiers.pixel.PixelClassifier
import qupath.lib.io.GsonTools

int    CHANNEL   = 2        // DAPI 0, FITC 1, Cy3 2
double THRESHOLD = 900.0    // raw counts, see calibration above
double PIXEL_UM  = 1.3      // resolution the threshold was calibrated at
String NAME      = 'Abeta_threshold_900'

String outDir = '.'
try {
    if (binding.hasVariable('args') && args != null && args.length > 0 && args[0])
        outDir = args[0] as String
} catch (Throwable ignored) { }

def cal = new PixelCalibration.Builder().pixelSizeMicrons(PIXEL_UM, PIXEL_UM).build()
def classifier = PixelClassifiers.createThresholdClassifier(
        cal, CHANNEL, THRESHOLD,
        PathClass.fromString('Negative'), PathClass.fromString('Abeta'))

def gson = GsonTools.getInstance(true)
def dir = new File(outDir)
dir.mkdirs()
def out = new File(dir, NAME + '.json')

// Serialise against the DECLARED interface, not the concrete class: the polymorphic type
// adapter only writes the `pixel_classifier_type` discriminator that way, and without it
// QuPath cannot read the file back.
out.text = gson.toJson(classifier, PixelClassifier.class)

def back = gson.fromJson(out.text, PixelClassifier.class)   // fail loudly, not at load time
println "wrote ${out.absolutePath} (${out.length()} bytes)"
println "round-trip OK: ${back.getClass().simpleName}, " +
        "${back.getMetadata().getInputResolution().getPixelWidthMicrons()} um/px, " +
        "labels ${back.getMetadata().getClassificationLabels()}"
