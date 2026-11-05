#!/usr/bin/env python3
"""Agreement between two coders, box by box: dots within MATCH_UM of each other count as the same deposit.
Writes ~/ihc_work/manual_coding/comparison.csv and prints the table. Usage: compare_coders.py [coderA] [coderB]"""
import sys, json, math
from pathlib import Path
import pandas as pd
M = Path.home() / "ihc_work/manual_coding"; MATCH_UM = 15.0
coders = sorted(p.name for p in (M / "coders").iterdir() if p.is_dir()); A, B = (sys.argv[1:3] + coders)[:2] if len(coders) >= 2 or len(sys.argv) >= 3 else (None, None)
if not A or not B: sys.exit(f"need two coders with saved boxes; found {coders}")
rows = []
for fa in sorted((M / "coders" / A).glob("*.json")):
    fb = M / "coders" / B / fa.name
    if not fb.exists(): continue
    a, b = json.loads(fa.read_text()), json.loads(fb.read_text()); px = a["px_um"]
    da = [d for d in a["dots"] if d["label"] == "deposit"]; db = [d for d in b["dots"] if d["label"] == "deposit"]; used = set(); m = 0
    for d in da:
        best = None
        for j, o in enumerate(db):
            if j in used: continue
            dist = math.hypot(d["x"] - o["x"], d["y"] - o["y"]) * px
            if dist <= MATCH_UM and (best is None or dist < best[0]): best = (dist, j)
        if best: used.add(best[1]); m += 1
    meta = json.loads((M / "data" / fa.name).read_text())
    rows.append(dict(image=a["image"], box=a["box"], label=meta["label"], **{f"{A}_deposits": len(da), f"{B}_deposits": len(db)}, agreed=m, **{f"only_{A}": len(da) - m, f"only_{B}": len(db) - m},
                     agreement=round(2 * m / max(len(da) + len(db), 1), 2), v4_pct_of_box=meta["v4_pct_of_box"], gate3_pct_of_box=meta["gate3_pct_of_box"]))
t = pd.DataFrame(rows); t.to_csv(M / "comparison.csv", index=False); print(t.to_string(index=False)) if len(t) else print("no boxes saved by both coders yet")
if len(t): print(f"\noverall: agreed {t.agreed.sum()}, only {A} {t[f'only_{A}'].sum()}, only {B} {t[f'only_{B}'].sum()}, Dice-style agreement {2*t.agreed.sum()/max(t[f'{A}_deposits'].sum()+t[f'{B}_deposits'].sum(),1):.2f}")
