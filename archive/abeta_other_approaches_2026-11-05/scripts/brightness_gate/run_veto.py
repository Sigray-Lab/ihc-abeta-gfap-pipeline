import sys, json, time
from pathlib import Path
from multiprocessing import Pool
sys.path.insert(0, str(Path(__file__).parent)); from veto import veto_section
G = Path(__file__).parent; GJ = G / "geojson"; R = G / "results_v6"; R.mkdir(exist_ok=True); LOG = G / "run_veto.log"
def one(p):
    try:
        res = veto_section(p.stem, p); (R / f"{p.stem}.json").write_text(json.dumps(res)); return f"ok {p.stem} v4 {100*res['abeta_v4_um2']/res['tissue_um2']:.2f}% v6 {100*res['abeta_v6_um2']/res['tissue_um2']:.2f}% vetoed {100*res['vetoed_share']:.0f}%"
    except Exception as e: return f"FAIL {p.stem}: {e!r}"
if __name__ == "__main__":
    todo = [p for p in sorted(GJ.glob("*.geojson")) if p.stat().st_size > 0 and not (R / f"{p.stem}.json").exists()]
    with Pool(6) as pool:
        for line in pool.imap_unordered(one, todo): LOG.open("a").write(f"{time.strftime('%H:%M:%S')} {line}\n")
    LOG.open("a").write(f"{time.strftime('%H:%M:%S')} VETO DONE {len(list(R.glob('*.json')))} results\n")
