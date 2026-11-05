#!/usr/bin/env python3
"""Summarise the outline ratings of all coders. Writes ~/ihc_work/manual_coding/ratings_summary.csv and ratings_all.csv."""
import pandas as pd
from pathlib import Path
M = Path.home() / "ihc_work/manual_coding"; files = list((M / "coders").glob("*/ratings.csv"))
if not files: raise SystemExit("no ratings yet")
d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True); d.to_csv(M / "ratings_all.csv", index=False)
order = ["too overinclusive", "about right", "too underinclusive", "unsure"]
unknown = set(d.verdict) - set(order)
if unknown: raise SystemExit(f"unknown verdict labels {unknown}: fix the list above rather than drop them")
print(f"{len(d)} ratings from {d.coder.nunique()} coder(s) on {d.groupby(['image','box']).ngroups} boxes\n")
print("== verdicts per method (all coders pooled)"); t = pd.crosstab(d.method, d.verdict).reindex(columns=[c for c in order if c in d.verdict.unique()], fill_value=0); t["about_right_share"] = (t.get("about right", 0) / t.sum(axis=1)).round(2); print(t.to_string())
print("\n== per coder"); print(pd.crosstab([d.coder, d.method], d.verdict).to_string())
d["kind"] = d.label.str.contains("control").map({True: "control", False: "stained"})
print("\n== stained vs control boxes"); print(pd.crosstab([d.kind, d.method], d.verdict).to_string())
if d.coder.nunique() >= 2:
    w = d.pivot_table(index=["image", "box", "method"], columns="coder", values="verdict", aggfunc="first").dropna(); coders = list(w.columns)
    agree = (w[coders[0]] == w[coders[1]]).mean() if len(coders) >= 2 else float("nan")
    print(f"\n== coder agreement on identical (box, method) pairs: {agree:.2f} over {len(w)} pairs"); print(pd.crosstab(w[coders[0]], w[coders[1]]).to_string())
    pm = w.reset_index().groupby("method").apply(lambda g: (g[coders[0]] == g[coders[1]]).mean(), include_groups=False).round(2); print("\nagreement per method:"); print(pm.to_string())
s = pd.crosstab([d.method], d.verdict); s.to_csv(M / "ratings_summary.csv"); print(f"\ncomments:\n" + "\n".join(f"  [{r.coder}] {r.image} {r.box} {r.method}: {r.comment}" for r in d[d.comment.notna() & (d.comment.astype(str).str.strip() != '')].itertuples()))
