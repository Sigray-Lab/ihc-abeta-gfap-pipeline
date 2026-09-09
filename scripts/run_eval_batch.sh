#!/bin/bash
# Batch-evaluate a pixel classifier inside per-image ROIs, N QuPath instances in parallel.
#
#   scripts/run_eval_batch.sh <nworkers> <roi_dir> <classifier_spec> <out_csv> [region_tag]
#
# classifier_spec is a literal name (v2_z130_s1) or a template with {image}
# (gx_{image}, rx_Hippocampus_{image}). Images are the *.geojson stems in roi_dir, so a
# region directory evaluates exactly the sections that have that region. Each worker gets
# its own copy of the project (QuPath writes into a project when it opens it) and its own
# output file; outputs are merged at the end with a `region` column when region_tag is given.
# Replaces the /tmp-only run_full_parallel.sh that produced the August measurements
# (archived under qupath/scripts/archive_run_scripts/).
set -u
N=${1:?nworkers}; ROI=${2:?roi_dir}; SPEC=${3:?classifier_spec}; OUT=${4:?out_csv}; TAG=${5:-}
Q=$HOME/Applications/QuPath.app/Contents/MacOS/QuPath-0.7.0-arm64
SRC=$HOME/ihc_work/qupath
HERE=$(cd "$(dirname "$0")/.." && pwd)
GROOVY=$HERE/qupath/scripts/eval_in_tissue.groovy
HDR=$HERE/qupath/scripts/eval_in_tissue_header.txt
WORK=$(mktemp -d "${TMPDIR:-/tmp}/ihc_eval.XXXXXX")
LOG=$OUT.log; : > "$LOG"
say(){ echo "$(date +%T) $*" >> "$LOG"; }
ls "$ROI"/*.geojson | sed 's|.*/||; s|\.geojson$||' | sort > "$WORK/all.txt"
python3 -c "import json,sys; print('\n'.join(sorted(e['image_name'] for e in json.load(open(sys.argv[1]))['images'])))" "$SRC/project_spec.json" | sort > "$WORK/proj.txt"
comm -12 "$WORK/all.txt" "$WORK/proj.txt" > "$WORK/images.txt"     # only ROIs whose image is in the project
TOTAL=$(wc -l < "$WORK/images.txt" | tr -d ' ')
say "$TOTAL images, $N workers, spec=$SPEC roi=$ROI"
for w in $(seq 0 $((N-1))); do
  awk -v w=$w -v n=$N 'NR % n == w' "$WORK/images.txt" > "$WORK/w$w.txt"
  cp -R "$SRC" "$WORK/p$w"
  ( cp "$HDR" "$WORK/w$w.csv"
    while read -r IMG; do
      [ -z "$IMG" ] && continue
      C=${SPEC//\{image\}/$IMG}
      "$Q" script -p "$WORK/p$w/qupath/project.qpproj" -i "$IMG" --args "$C,$WORK/w$w.csv,$ROI" "$GROOVY" >/dev/null 2>&1 \
        || say "WORKER$w FAILED $C $IMG"
    done < "$WORK/w$w.txt"
    say "worker $w done ($(( $(wc -l < "$WORK/w$w.csv") - 1 )) rows)" ) &
done
wait
python3 - "$WORK" "$OUT" "$TAG" <<'PY'
import sys,glob,pandas as pd
w,out,tag=sys.argv[1:4]
d=pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{w}/w*.csv"))],ignore_index=True)
if tag: d.insert(0,"region",tag)
d.to_csv(out,index=False); print(f"{len(d)} rows -> {out}")
PY
say "DONE"; rm -rf "$WORK"
