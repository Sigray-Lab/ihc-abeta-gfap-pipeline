#!/usr/bin/env python3
"""Plaque Coder: a small viewer for manual coding of Abeta deposits inside fixed boxes.

Two coders (Pontus, Aris) each click the deposits they see in the same boxes; the clicks are saved per coder and
box and can be compared afterwards (compare_coders.py, or the Compare mode here). A Practice mode walks through
spots the classifiers flagged, from typical deposits to things that only one method calls amyloid, and reveals
what each method says. Data is prepared by prep_boxes.py into ~/ihc_work/manual_coding/data/.

Keys: 1 Cy3, 2 DAPI, 3 GFAP, 4 merge | z undo | s save | n next box | p previous box | + / - zoom | left click = deposit,
right click = uncertain, click on an existing dot removes it."""
import json, sys, time, math
from pathlib import Path
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
from PIL import Image, ImageTk, ImageDraw
M = Path.home() / "ihc_work/manual_coding"; DATA = M / "data"; CODERS = M / "coders"
COLORS = {"deposit": "#ff3b30", "uncertain": "#ffd60a"}; OTHER = "#30d5ff"; MATCH_UM = 15.0
def stretch(a, lo_p=1.0, hi_p=99.7):
    lo, hi = np.percentile(a, [lo_p, hi_p]); hi = max(hi, lo + 1)
    return np.clip((a.astype(np.float32) - lo) / (hi - lo), 0, 1)
class Coder(tk.Tk):
    def __init__(self):
        super().__init__(); self.title("Plaque Coder"); self.geometry("1380x980")
        self.boxes = sorted(DATA.glob("*.json"), key=lambda p: p.name)
        if not self.boxes: messagebox.showerror("Plaque Coder", f"No boxes found in {DATA}. Run prep_boxes.py first."); self.destroy(); return
        pre = next((sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--coder" and i + 1 < len(sys.argv)), None)
        self.coder = (pre or simpledialog.askstring("Who is coding?", "Your name (Pontus / Aris):", parent=self) or "anon").strip().capitalize()
        (CODERS / self.coder).mkdir(parents=True, exist_ok=True)
        self.mode = tk.StringVar(value="Code"); self.channel = tk.StringVar(value="merge"); self.hi_p = tk.DoubleVar(value=99.7); self.zoom = 1.0
        self.show_v4 = tk.BooleanVar(value=False); self.show_gate3 = tk.BooleanVar(value=False); self.show_veto = tk.BooleanVar(value=False); self.show_other = tk.BooleanVar(value=False)
        self.i = 0; self.dots = []; self.undo = []; self.spot_i = -1; self.revealed = False; self.dirty = False
        self._build(); self.load_box(0)
        self.protocol("WM_DELETE_WINDOW", self.on_close)
    # ---------------- layout ----------------
    def _build(self):
        top = ttk.Frame(self); top.pack(side="top", fill="x", padx=6, pady=4)
        ttk.Label(top, text=f"Coder: {self.coder}", font=("Helvetica", 13, "bold")).pack(side="left", padx=(0, 16))
        for m in ("Code", "Practice", "Compare"): ttk.Radiobutton(top, text=m, value=m, variable=self.mode, command=self.redraw).pack(side="left")
        ttk.Label(top, text="   Channel:").pack(side="left")
        for txt, val, key in (("Cy3 (Aβ)", "cy3", "1"), ("DAPI", "dapi", "2"), ("GFAP", "fitc", "3"), ("Merge", "merge", "4")):
            ttk.Radiobutton(top, text=f"{txt} [{key}]", value=val, variable=self.channel, command=self.render).pack(side="left")
        ttk.Label(top, text="   Contrast:").pack(side="left"); ttk.Scale(top, from_=98.0, to=99.99, variable=self.hi_p, command=lambda e: self.render(), length=140).pack(side="left")
        ttk.Button(top, text="Zoom − [-]", command=lambda: self.set_zoom(self.zoom / 1.5)).pack(side="left", padx=(16, 2)); ttk.Button(top, text="Zoom + [+]", command=lambda: self.set_zoom(self.zoom * 1.5)).pack(side="left")
        self.status = ttk.Label(self, text="", anchor="w", font=("Helvetica", 12)); self.status.pack(side="bottom", fill="x", padx=8, pady=4)
        body = ttk.Frame(self); body.pack(fill="both", expand=True)
        right = ttk.Frame(body, width=330); right.pack(side="right", fill="y", padx=6); right.pack_propagate(False)
        ttk.Label(right, text="Boxes", font=("Helvetica", 12, "bold")).pack(anchor="w")
        self.listbox = tk.Listbox(right, height=16, font=("Menlo", 11), exportselection=False); self.listbox.pack(fill="x"); self.listbox.bind("<<ListboxSelect>>", lambda e: self.load_box(self.listbox.curselection()[0]) if self.listbox.curselection() else None)
        nav = ttk.Frame(right); nav.pack(fill="x", pady=4); ttk.Button(nav, text="◀ Prev [p]", command=lambda: self.load_box(self.i - 1)).pack(side="left", expand=True, fill="x"); ttk.Button(nav, text="Next [n] ▶", command=lambda: self.load_box(self.i + 1)).pack(side="left", expand=True, fill="x")
        self.info = tk.Text(right, height=6, wrap="word", font=("Helvetica", 11), relief="flat", bg=self.cget("bg")); self.info.pack(fill="x", pady=4)
        cf = ttk.LabelFrame(right, text="Your dots"); cf.pack(fill="x", pady=4)
        self.count_lbl = ttk.Label(cf, text="", font=("Helvetica", 12)); self.count_lbl.pack(anchor="w", padx=6)
        bf = ttk.Frame(cf); bf.pack(fill="x", padx=4, pady=4)
        ttk.Button(bf, text="Undo [z]", command=self.do_undo).pack(side="left", expand=True, fill="x"); ttk.Button(bf, text="Clear box", command=self.do_clear).pack(side="left", expand=True, fill="x")
        ttk.Button(cf, text="Save [s]", command=self.save).pack(fill="x", padx=4, pady=(0, 4))
        self.done_var = tk.BooleanVar(value=False); ttk.Checkbutton(cf, text="Mark this box as finished", variable=self.done_var, command=lambda: self._mark_dirty()).pack(anchor="w", padx=6, pady=(0, 4))
        of = ttk.LabelFrame(right, text="Overlays (what the machine says)"); of.pack(fill="x", pady=4)
        ttk.Checkbutton(of, text="v4 classifier outline (red)", variable=self.show_v4, command=self.redraw).pack(anchor="w", padx=6)
        ttk.Checkbutton(of, text="3× brightness gate (yellow fill)", variable=self.show_gate3, command=self.render).pack(anchor="w", padx=6)
        ttk.Checkbutton(of, text="DAPI-veto removals (magenta fill)", variable=self.show_veto, command=self.render).pack(anchor="w", padx=6)
        ttk.Checkbutton(of, text="other coder's dots (cyan)", variable=self.show_other, command=self.redraw).pack(anchor="w", padx=6)
        pf = ttk.LabelFrame(right, text="Practice spots"); pf.pack(fill="x", pady=4); self.pf = pf
        sb = ttk.Frame(pf); sb.pack(fill="x", padx=4, pady=2)
        ttk.Button(sb, text="◀ spot", command=lambda: self.goto_spot(self.spot_i - 1)).pack(side="left", expand=True, fill="x"); ttk.Button(sb, text="spot ▶", command=lambda: self.goto_spot(self.spot_i + 1)).pack(side="left", expand=True, fill="x")
        ttk.Button(pf, text="Reveal what the machine says", command=self.reveal).pack(fill="x", padx=4, pady=2)
        self.spot_txt = tk.Text(pf, height=7, wrap="word", font=("Helvetica", 11), relief="flat", bg=self.cget("bg")); self.spot_txt.pack(fill="x", padx=4)
        self.cmp_txt = tk.Text(right, height=5, wrap="word", font=("Helvetica", 11), relief="flat", bg=self.cget("bg")); self.cmp_txt.pack(fill="x", pady=4)
        cv = ttk.Frame(body); cv.pack(side="left", fill="both", expand=True)
        self.canvas = tk.Canvas(cv, bg="black", highlightthickness=0); hs = ttk.Scrollbar(cv, orient="horizontal", command=self.canvas.xview); vs = ttk.Scrollbar(cv, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=hs.set, yscrollcommand=vs.set); vs.pack(side="right", fill="y"); hs.pack(side="bottom", fill="x"); self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Button-1>", lambda e: self.click(e, "deposit")); self.canvas.bind("<Button-2>", lambda e: self.click(e, "uncertain")); self.canvas.bind("<Button-3>", lambda e: self.click(e, "uncertain")); self.canvas.bind("<Control-Button-1>", lambda e: self.click(e, "uncertain"))
        for k, f in (("1", lambda e: self._ch("cy3")), ("2", lambda e: self._ch("dapi")), ("3", lambda e: self._ch("fitc")), ("4", lambda e: self._ch("merge")), ("z", lambda e: self.do_undo()), ("s", lambda e: self.save()), ("n", lambda e: self.load_box(self.i + 1)), ("p", lambda e: self.load_box(self.i - 1)), ("plus", lambda e: self.set_zoom(self.zoom * 1.5)), ("minus", lambda e: self.set_zoom(self.zoom / 1.5)), ("equal", lambda e: self.set_zoom(self.zoom * 1.5))):
            self.bind(f"<KeyPress-{k}>", f)
        self.bind("<Command-z>", lambda e: self.do_undo()); self.bind("<Command-s>", lambda e: self.save())
        self.refresh_list()
    def _ch(self, v): self.channel.set(v); self.render()
    # ---------------- data ----------------
    def coder_file(self, coder=None): return CODERS / (coder or self.coder) / f"{self.meta['image']}_{self.meta['box']}.json"
    def refresh_list(self):
        self.listbox.delete(0, "end")
        for p in self.boxes:
            m = json.loads(p.read_text()); f = CODERS / self.coder / f"{m['image']}_{m['box']}.json"; mark = " "
            if f.exists(): d = json.loads(f.read_text()); mark = "✓" if d.get("done") else "·"
            self.listbox.insert("end", f"{mark} {m['image']} {m['box']}  {m['side_um']:.0f} µm  {'CONTROL' if m['kind']=='control' else ''}")
    def load_box(self, i):
        if hasattr(self, "meta") and self.dirty: self.save(quiet=True)
        i = max(0, min(i, len(self.boxes) - 1)); self.i = i; self.meta = json.loads(self.boxes[i].read_text()); stem = self.boxes[i].stem
        z = np.load(DATA / f"{stem}.npz"); self.img = {k: z[k] for k in ("cy3", "dapi", "fitc")}
        self.masks = {k: np.array(Image.open(DATA / f"{stem}_{k}.png")) > 0 for k in ("gate3", "veto") if (DATA / f"{stem}_{k}.png").exists()}
        f = self.coder_file(); d = json.loads(f.read_text()) if f.exists() else {}
        self.dots = [dict(x=float(q["x"]), y=float(q["y"]), label=q.get("label", "deposit")) for q in d.get("dots", [])]; self.done_var.set(bool(d.get("done", False))); self.undo = []; self.dirty = False
        self.spot_i = -1; self.revealed = False; self.listbox.selection_clear(0, "end"); self.listbox.selection_set(i); self.listbox.see(i)
        self.info.delete("1.0", "end"); self.info.insert("end", f"{self.meta['image']} {self.meta['box']} — {self.meta['side_um']:.0f} µm box\n{self.meta['label']}\n\nLeft click = deposit, right click = uncertain, click a dot to remove it.")
        self.zoom = min(1.0, 900 / max(self.img['cy3'].shape)); self.render(); self.spot_txt.delete("1.0", "end")
    # ---------------- rendering ----------------
    def base_rgb(self):
        hi = float(self.hi_p.get()); ch = self.channel.get()
        if ch == "merge": r, g, b = stretch(self.img["cy3"], 1, hi), stretch(self.img["fitc"], 1, hi), stretch(self.img["dapi"], 1, hi)
        else: v = stretch(self.img[ch], 1, hi); r, g, b = v, v, v
        rgb = np.dstack([r, g, b])
        if self.show_gate3.get() and "gate3" in self.masks: m = self.masks["gate3"]; rgb[m] = 0.55 * rgb[m] + 0.45 * np.array([1, 1, 0])
        if self.show_veto.get() and "veto" in self.masks: m = self.masks["veto"]; rgb[m] = 0.45 * rgb[m] + 0.55 * np.array([1, 0, 1])
        return Image.fromarray((rgb * 255).astype(np.uint8))
    def render(self):
        im = self.base_rgb(); H, W = im.size[1], im.size[0]; z = self.zoom
        if abs(z - 1) > 1e-6: im = im.resize((max(1, int(W * z)), max(1, int(H * z))), Image.BILINEAR)
        self.photo = ImageTk.PhotoImage(im); self.canvas.delete("all"); self.canvas.create_image(0, 0, anchor="nw", image=self.photo); self.canvas.configure(scrollregion=(0, 0, im.size[0], im.size[1])); self.redraw()
    def redraw(self):
        self.canvas.delete("ov"); z = self.zoom; px = self.meta["px_um"]
        if self.show_v4.get():
            for poly in self.meta.get("v4_polys", []):
                if len(poly) >= 3: self.canvas.create_polygon([c * z for p in poly for c in p], outline="#ff5050", fill="", width=1, tags="ov")
        if self.show_other.get() or self.mode.get() == "Compare":
            for o in self.other_dots():
                r = 7; self.canvas.create_oval(o["x"] * z - r, o["y"] * z - r, o["x"] * z + r, o["y"] * z + r, outline=OTHER, width=2, tags="ov")
        for d in self.dots:
            r = 6; self.canvas.create_oval(d["x"] * z - r, d["y"] * z - r, d["x"] * z + r, d["y"] * z + r, fill=COLORS[d["label"]], outline="white", width=1.5, tags="ov")
        if self.mode.get() == "Practice" and 0 <= self.spot_i < len(self.meta.get("spots", [])):
            s = self.meta["spots"][self.spot_i]; r = s["r"] * z; self.canvas.create_oval(s["x"] * z - r, s["y"] * z - r, s["x"] * z + r, s["y"] * z + r, outline="#00ff9d", width=3, tags="ov")
        n_dep = sum(d["label"] == "deposit" for d in self.dots); n_unc = len(self.dots) - n_dep
        self.count_lbl.configure(text=f"deposits: {n_dep}    uncertain: {n_unc}")
        self.status.configure(text=f"{self.meta['image']} {self.meta['box']}   zoom {self.zoom:.2f}   mode {self.mode.get()}   {'(unsaved)' if self.dirty else 'saved'}   box {self.i+1}/{len(self.boxes)}")
        self.update_compare()
    def set_zoom(self, z): self.zoom = max(0.2, min(4.0, z)); self.render()
    # ---------------- clicks ----------------
    def click(self, e, label):
        if self.mode.get() == "Compare": return
        x, y = self.canvas.canvasx(e.x) / self.zoom, self.canvas.canvasy(e.y) / self.zoom; px = self.meta["px_um"]
        near = [d for d in self.dots if math.hypot(d["x"] - x, d["y"] - y) * px < 6.0]
        if near: d = near[0]; self.dots.remove(d); self.undo.append(("add", d))
        else: d = dict(x=x, y=y, label=label); self.dots.append(d); self.undo.append(("remove", d))
        self._mark_dirty(); self.redraw()
    def do_undo(self):
        if not self.undo: return
        op, d = self.undo.pop()
        if op == "remove" and d in self.dots: self.dots.remove(d)
        elif op == "add": self.dots.append(d)
        self._mark_dirty(); self.redraw()
    def do_clear(self):
        if self.dots and messagebox.askyesno("Clear", f"Remove all {len(self.dots)} dots in this box?"): self.undo.extend(("add", d) for d in self.dots); self.dots = []; self._mark_dirty(); self.redraw()
    def _mark_dirty(self): self.dirty = True; self.status.configure(text=self.status.cget("text").replace("saved", "(unsaved)"))
    def save(self, quiet=False):
        px = self.meta["px_um"]; x0, y0 = self.meta["bounds"][0], self.meta["bounds"][1]
        out = dict(coder=self.coder, image=self.meta["image"], box=self.meta["box"], saved_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), px_um=px, bounds=self.meta["bounds"], done=bool(self.done_var.get()),
                   n_deposit=sum(d["label"] == "deposit" for d in self.dots), n_uncertain=sum(d["label"] == "uncertain" for d in self.dots),
                   dots=[dict(x=round(d["x"], 1), y=round(d["y"], 1), x_fullres=round(x0 + d["x"], 1), y_fullres=round(y0 + d["y"], 1), label=d["label"]) for d in self.dots])
        self.coder_file().write_text(json.dumps(out, indent=1)); self.dirty = False; self.refresh_list(); self.listbox.selection_set(self.i); self.redraw()
        if not quiet: self.status.configure(text=f"saved {self.coder_file().name}  ({out['n_deposit']} deposits, {out['n_uncertain']} uncertain)")
    def on_close(self):
        if self.dirty: self.save(quiet=True)
        self.destroy()
    # ---------------- practice ----------------
    def goto_spot(self, k):
        spots = self.meta.get("spots", [])
        if not spots: self.spot_txt.delete("1.0", "end"); self.spot_txt.insert("end", "No practice spots in this box."); return
        self.mode.set("Practice"); self.spot_i = k % len(spots); self.revealed = False; s = spots[self.spot_i]
        self.spot_txt.delete("1.0", "end"); self.spot_txt.insert("end", f"Spot {self.spot_i+1} of {len(spots)}  ({s['area_um2']} µm²)\nLook at the green circle. Is this a deposit? Decide, then press Reveal.")
        z = self.zoom; self.canvas.xview_moveto(max(0, (s["x"] * z - self.canvas.winfo_width() / 2) / max(1, self.canvas.bbox("all")[2]))); self.canvas.yview_moveto(max(0, (s["y"] * z - self.canvas.winfo_height() / 2) / max(1, self.canvas.bbox("all")[3]))); self.redraw()
    def reveal(self):
        spots = self.meta.get("spots", [])
        if not (0 <= self.spot_i < len(spots)): return
        s = spots[self.spot_i]; self.revealed = True
        self.spot_txt.delete("1.0", "end"); self.spot_txt.insert("end", f"Spot {self.spot_i+1}/{len(spots)} — {s['category'].upper()}\n{s['text']}\n\nnumbers: area {s['area_um2']} µm², Cy3 {s['cy3_contrast']}× background, DAPI {s['dapi_contrast']}×, survives 3× gate: {int(100*s['gate3_survival'])}%")
    # ---------------- compare ----------------
    def other_dots(self):
        out = []
        for c in CODERS.iterdir():
            if c.name == self.coder or not c.is_dir(): continue
            f = c / f"{self.meta['image']}_{self.meta['box']}.json"
            if f.exists(): out.extend(dict(x=float(q["x"]), y=float(q["y"]), label=q.get("label", "deposit"), coder=c.name) for q in json.loads(f.read_text()).get("dots", []))
        return out
    def update_compare(self):
        self.cmp_txt.delete("1.0", "end")
        if self.mode.get() != "Compare" and not self.show_other.get(): return
        other = [o for o in self.other_dots() if o["label"] == "deposit"]; mine = [d for d in self.dots if d["label"] == "deposit"]
        if not other: self.cmp_txt.insert("end", "No other coder has saved this box yet."); return
        px = self.meta["px_um"]; used = set(); matched = 0
        for d in mine:
            best = None
            for j, o in enumerate(other):
                if j in used: continue
                dist = math.hypot(d["x"] - o["x"], d["y"] - o["y"]) * px
                if dist <= MATCH_UM and (best is None or dist < best[0]): best = (dist, j)
            if best: used.add(best[1]); matched += 1
        who = other[0]["coder"]; self.cmp_txt.insert("end", f"Compare with {who} (match within {MATCH_UM:.0f} µm):\n  agreed: {matched}\n  only you: {len(mine) - matched}\n  only {who}: {len(other) - matched}\n  cyan rings = {who}'s dots")
if __name__ == "__main__":
    app = Coder()
    if "--selftest" in sys.argv: app.after(1500, app.destroy)
    app.mainloop()
