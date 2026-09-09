#!/bin/bash
# export v4 (gx_<image>) Abeta objects for every section in the tissue-ROI set, N workers, own project copy each
S=/private/tmp/claude-501/-Users-pontusps-Documents/2021b057-59a2-47a1-9a78-0e1252f4a4ae/scratchpad
Q=$HOME/Applications/QuPath.app/Contents/MacOS/QuPath-0.7.0-arm64
N=3; OUT=$S/gate/geojson; LOG=$S/gate/export.log; : > "$LOG"
ls $HOME/ihc_work/tissue_rois/*.geojson | sed 's|.*/||; s|\.geojson$||' | sort > "$S/gate/all.txt"
python3 -c "import json,sys; print('\n'.join(sorted(e['image_name'] for e in json.load(open(sys.argv[1]))['images'])))" "$HOME/ihc_work/qupath/project_spec.json" | sort > "$S/gate/proj.txt"
comm -12 "$S/gate/all.txt" "$S/gate/proj.txt" > "$S/gate/images.txt"
echo "$(date +%T) $(wc -l < "$S/gate/images.txt") images, $N workers" >> "$LOG"
for w in $(seq 0 $((N-1))); do
  awk -v w=$w -v n=$N 'NR % n == w' "$S/gate/images.txt" > "$S/gate/w$w.txt"
  rm -rf "$S/gate/p$w"; cp -R "$HOME/ihc_work/qupath" "$S/gate/p$w"
  ( while read -r IMG; do
      [ -z "$IMG" ] && continue; [ -s "$OUT/$IMG.geojson" ] && continue
      "$Q" script -p "$S/gate/p$w/qupath/project.qpproj" -i "$IMG" --args "gx_$IMG,$OUT/$IMG.geojson,0" "$HOME/ihc_work/auto_negatives/export_fp.groovy" >/dev/null 2>&1 \
        || echo "$(date +%T) WORKER$w FAILED $IMG" >> "$LOG"
    done < "$S/gate/w$w.txt"; echo "$(date +%T) worker $w done" >> "$LOG" ) &
done
wait; echo "$(date +%T) DONE $(ls "$OUT" | wc -l) files" >> "$LOG"; rm -rf "$S"/gate/p[0-9]
