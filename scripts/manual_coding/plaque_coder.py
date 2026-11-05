#!/usr/bin/env python3
"""Plaque Coder — box viewer for manual assessment of Abeta classification.

Modes
  Code      click deposits (left = deposit, right = uncertain, click a dot to remove); saved per coder and box.
  Practice  step through spots the classifiers flagged and reveal what each method says (teaching).
  Rate      outline rating: each box shows the candidate outlines one at a time, by name, fixed order; press g (too
            overinclusive), a (about right), t (too underinclusive) or u (unsure), optional comment, auto-advance. Ratings go to coders/<name>/ratings.csv; tabulate_ratings.py unblinds and summarises.
  Compare   the other coder's dots and agreement counts.
Keys: 1 Cy3  2 DAPI  3 GFAP  4 merge | o overlay on/off | z undo | s save | n/p next/prev box | +/- zoom
      Rate: g a t u verdicts, left/right arrows = previous/next outline."""
import json, sys, time, math, csv
from pathlib import Path
import numpy as np, cv2
import tkinter as tk
from tkinter import ttk, messagebox, simpledialog
from PIL import Image, ImageTk
M = Path.home() / "ihc_work/manual_coding"; DATA = M / "data"; CODERS = M / "coders"
COLORS = {"deposit": "#ff3b30", "uncertain": "#ffd60a"}; OTHER = "#30d5ff"; MATCH_UM = 15.0; OUTLINE = "#ff2bd6"
METHODS = ["v4", "gate1.5", "gate2", "gate3", "v6", "v5"]          # candidate outlines; v5 only where exported
VERDICTS = {"g": "too overinclusive", "a": "about right", "t": "too underinclusive", "u": "unsure"}
METHOD_NAMES = {"v4": "v4 GLOBAL (as is)", "gate1.5": "v4 + 1.5× brightness gate", "gate2": "v4 + 2× brightness gate", "gate3": "v4 + 3× brightness gate", "v6": "v6 = v4 minus DAPI veto", "v5": "v5 (retrained with control negatives)"}
def stretch(a, lo_p=1.0, hi_p=99.7):
    lo, hi = np.percentile(a, [lo_p, hi_p]); hi = max(hi, lo + 1)
    return np.clip((a.astype(np.float32) - lo) / (hi - lo), 0, 1)
class Coder(tk.Tk):
    def __init__(self):
        super().__init__(); self.title("Plaque Coder"); self.geometry("1400x1000")
        self.all_boxes = sorted(DATA.glob("*.json"), key=lambda p: (json.loads(p.read_text()).get("set", "practice") != "rating", p.name))
        self.boxes = self.all_boxes
        if not self.boxes: messagebox.showerror("Plaque Coder", f"No boxes found in {DATA}. Run prep_boxes.py first."); self.destroy(); return
        pre = next((sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == "--coder" and i + 1 < len(sys.argv)), None)
        self.coder = (pre or simpledialog.askstring("Who is coding?", "Your name (Pontus / Aris):", parent=self) or "anon").strip().capitalize()
        (CODERS / self.coder).mkdir(parents=True, exist_ok=True)
        self.mode = tk.StringVar(value="Rate"); self.channel = tk.StringVar(value="cy3"); self.hi_p = tk.DoubleVar(value=99.7); self.zoom = 1.0
        self.show_v4 = tk.BooleanVar(value=False); self.show_gate3 = tk.BooleanVar(value=False); self.show_veto = tk.BooleanVar(value=False); self.show_other = tk.BooleanVar(value=False)
        self.overlay_on = True; self.i = 0; self.dots = []; self.undo = []; self.spot_i = -1; self.dirty = False; self.m_i = 0; self.contour_cache = {}
        self._build(); self.load_box(0); self.protocol("WM_DELETE_WINDOW", self.on_close)
    # ---------------- layout ----------------
    def _build(self):
        top = ttk.Frame(self); top.pack(side="top", fill="x", padx=6, pady=4)
        ttk.Label(top, text=f"Coder: {self.coder}", font=("Helvetica", 13, "bold")).pack(side="left", padx=(0, 16))
        for m in ("Rate", "Code", "Practice", "Compare"): ttk.Radiobutton(top, text=m, value=m, variable=self.mode, command=self.mode_changed).pack(side="left")
        ttk.Label(top, text="   Channel:").pack(side="left")
        for txt, val, key in (("Cy3 (Aβ)", "cy3", "1"), ("DAPI", "dapi", "2"), ("GFAP", "fitc", "3"), ("Merge", "merge", "4")):
            ttk.Radiobutton(top, text=f"{txt} [{key}]", value=val, variable=self.channel, command=self.render).pack(side="left")
        ttk.Label(top, text="   Contrast:").pack(side="left"); ttk.Scale(top, from_=98.0, to=99.99, variable=self.hi_p, command=lambda e: self.render(), length=140).pack(side="left")
        ttk.Button(top, text="Zoom − [-]", command=lambda: self.set_zoom(self.zoom / 1.5)).pack(side="left", padx=(16, 2)); ttk.Button(top, text="Zoom + [+]", command=lambda: self.set_zoom(self.zoom * 1.5)).pack(side="left")
        self.status = ttk.Label(self, text="", anchor="w", font=("Helvetica", 12)); self.status.pack(side="bottom", fill="x", padx=8, pady=4)
        body = ttk.Frame(self); body.pack(fill="both", expand=True)
        right = ttk.Frame(body, width=340); right.pack(side="right", fill="y", padx=6); right.pack_propagate(False)
        ttk.Label(right, text="Boxes", font=("Helvetica", 12, "bold")).pack(anchor="w")
        self.listbox = tk.Listbox(right, height=14, font=("Menlo", 11), exportselection=False); self.listbox.pack(fill="x")
        self.listbox.bind("<<ListboxSelect>>", lambda e: self.load_box(self.listbox.curselection()[0]) if self.listbox.curselection() else None)
        nav = ttk.Frame(right); nav.pack(fill="x", pady=4); ttk.Button(nav, text="◀ Prev box [p]", command=lambda: self.load_box(self.i - 1)).pack(side="left", expand=True, fill="x"); ttk.Button(nav, text="Next box [n] ▶", command=lambda: self.load_box(self.i + 1)).pack(side="left", expand=True, fill="x")
        self.info = tk.Text(right, height=4, wrap="word", font=("Helvetica", 11), relief="flat", bg=self.cget("bg")); self.info.pack(fill="x", pady=2)
        # ---- rating panel
        rf = ttk.LabelFrame(right, text="Rate this outline"); rf.pack(fill="x", pady=4); self.rf = rf
        self.rate_lbl = ttk.Label(rf, text="", font=("Helvetica", 13, "bold")); self.rate_lbl.pack(anchor="w", padx=6, pady=(4, 0))
        self.rate_prog = ttk.Label(rf, text="", font=("Helvetica", 11)); self.rate_prog.pack(anchor="w", padx=6)
        rb = ttk.Frame(rf); rb.pack(fill="x", padx=4, pady=4)
        for key, name in VERDICTS.items(): ttk.Button(rb, text=f"{name} [{key}]", command=lambda k=key: self.verdict(k)).pack(side="top", fill="x", pady=1)
        ttk.Label(rf, text="Comment (optional, saved with the verdict):").pack(anchor="w", padx=6)
        self.comment = ttk.Entry(rf); self.comment.pack(fill="x", padx=6, pady=(0, 4))
        rn = ttk.Frame(rf); rn.pack(fill="x", padx=4, pady=(0, 4))
        ttk.Button(rn, text="◀ outline [←]", command=lambda: self.goto_method(self.m_i - 1)).pack(side="left", expand=True, fill="x"); ttk.Button(rn, text="overlay on/off [o]", command=self.toggle_overlay).pack(side="left", expand=True, fill="x"); ttk.Button(rn, text="outline ▶ [→]", command=lambda: self.goto_method(self.m_i + 1)).pack(side="left", expand=True, fill="x")
        self.rated_txt = tk.Text(rf, height=7, wrap="word", font=("Menlo", 10), relief="flat", bg=self.cget("bg")); self.rated_txt.pack(fill="x", padx=4)
        # ---- dots panel
        cf = ttk.LabelFrame(right, text="Your dots (Code mode)"); cf.pack(fill="x", pady=4)
        self.count_lbl = ttk.Label(cf, text="", font=("Helvetica", 12)); self.count_lbl.pack(anchor="w", padx=6)
        bf = ttk.Frame(cf); bf.pack(fill="x", padx=4, pady=4)
        ttk.Button(bf, text="Undo [z]", command=self.do_undo).pack(side="left", expand=True, fill="x"); ttk.Button(bf, text="Clear box", command=self.do_clear).pack(side="left", expand=True, fill="x"); ttk.Button(bf, text="Save [s]", command=self.save).pack(side="left", expand=True, fill="x")
        self.done_var = tk.BooleanVar(value=False); ttk.Checkbutton(cf, text="box finished (dots)", variable=self.done_var, command=self._mark_dirty).pack(anchor="w", padx=6)
        # ---- overlays + practice
        of = ttk.LabelFrame(right, text="Overlays (Code / Practice only)"); of.pack(fill="x", pady=4); self.of = of
        self.cb_v4 = ttk.Checkbutton(of, text="v4 outline (red)", variable=self.show_v4, command=self.redraw); self.cb_v4.pack(anchor="w", padx=6)
        self.cb_g3 = ttk.Checkbutton(of, text="3× gate (yellow fill)", variable=self.show_gate3, command=self.render); self.cb_g3.pack(anchor="w", padx=6)
        self.cb_vt = ttk.Checkbutton(of, text="DAPI-veto removals (magenta fill)", variable=self.show_veto, command=self.render); self.cb_vt.pack(anchor="w", padx=6)
        ttk.Checkbutton(of, text="other coder's dots (cyan)", variable=self.show_other, command=self.redraw).pack(anchor="w", padx=6)
        pf = ttk.LabelFrame(right, text="Practice spots"); pf.pack(fill="x", pady=4)
        sb = ttk.Frame(pf); sb.pack(fill="x", padx=4, pady=2)
        ttk.Button(sb, text="◀ spot", command=lambda: self.goto_spot(self.spot_i - 1)).pack(side="left", expand=True, fill="x"); ttk.Button(sb, text="spot ▶", command=lambda: self.goto_spot(self.spot_i + 1)).pack(side="left", expand=True, fill="x"); ttk.Button(sb, text="Reveal", command=self.reveal).pack(side="left", expand=True, fill="x")
        self.spot_txt = tk.Text(pf, height=5, wrap="word", font=("Helvetica", 11), relief="flat", bg=self.cget("bg")); self.spot_txt.pack(fill="x", padx=4)
        self.cmp_txt = tk.Text(right, height=4, wrap="word", font=("Helvetica", 11), relief="flat", bg=self.cget("bg")); self.cmp_txt.pack(fill="x", pady=2)
        cv = ttk.Frame(body); cv.pack(side="left", fill="both", expand=True)
        self.canvas = tk.Canvas(cv, bg="black", highlightthickness=0); hs = ttk.Scrollbar(cv, orient="horizontal", command=self.canvas.xview); vs = ttk.Scrollbar(cv, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=hs.set, yscrollcommand=vs.set); vs.pack(side="right", fill="y"); hs.pack(side="bottom", fill="x"); self.canvas.pack(fill="both", expand=True)
        self.canvas.bind("<Button-1>", lambda e: self.click(e, "deposit")); self.canvas.bind("<Button-2>", lambda e: self.click(e, "uncertain")); self.canvas.bind("<Button-3>", lambda e: self.click(e, "uncertain")); self.canvas.bind("<Control-Button-1>", lambda e: self.click(e, "uncertain"))
        keys = {"1": lambda e: self._ch("cy3"), "2": lambda e: self._ch("dapi"), "3": lambda e: self._ch("fitc"), "4": lambda e: self._ch("merge"), "z": lambda e: self.do_undo(), "s": lambda e: self.save(), "n": lambda e: self.load_box(self.i + 1), "p": lambda e: self.load_box(self.i - 1),
                "plus": lambda e: self.set_zoom(self.zoom * 1.5), "equal": lambda e: self.set_zoom(self.zoom * 1.5), "minus": lambda e: self.set_zoom(self.zoom / 1.5), "o": lambda e: self.toggle_overlay(), "Left": lambda e: self.goto_method(self.m_i - 1), "Right": lambda e: self.goto_method(self.m_i + 1)}
        for k, v in VERDICTS.items(): keys[k] = (lambda e, kk=k: self.verdict(kk))
        for k, f in keys.items(): self.bind(f"<KeyPress-{k}>", lambda e, f=f: None if isinstance(self.focus_get(), ttk.Entry) and e.keysym not in ("Left", "Right") else f(e))
        self.bind("<Command-z>", lambda e: self.do_undo()); self.bind("<Command-s>", lambda e: self.save()); self.bind("<Return>", lambda e: self.canvas.focus_set())
        self.refresh_list()
    def _ch(self, v): self.channel.set(v); self.render()
    def mode_changed(self):
        want = [p for p in self.all_boxes if json.loads(p.read_text()).get("set") == "rating"] if self.mode.get() == "Rate" else self.all_boxes
        if want and want != self.boxes:
            cur = self.boxes[self.i]; self.boxes = want; self.refresh_list()
            if cur in self.boxes: self.i = self.boxes.index(cur); self.listbox.selection_clear(0, "end"); self.listbox.selection_set(self.i)
            else: self.load_box(0); return
        if self.mode.get() == "Rate":
            for v in (self.show_v4, self.show_gate3, self.show_veto): v.set(False)
            for cb in (self.cb_v4, self.cb_g3, self.cb_vt): cb.state(["disabled"])
        else:
            for cb in (self.cb_v4, self.cb_g3, self.cb_vt): cb.state(["!disabled"])
        self.render()
    # ---------------- data ----------------
    def coder_file(self, coder=None): return CODERS / (coder or self.coder) / f"{self.meta['image']}_{self.meta['box']}.json"
    def ratings_file(self): return CODERS / self.coder / "ratings.csv"
    def load_ratings(self):
        f = self.ratings_file(); out = {}
        if f.exists():
            for r in csv.DictReader(open(f)): out[(r["image"], r["box"], r["method"])] = r
        return out
    def refresh_list(self):
        self.listbox.delete(0, "end"); rt = self.load_ratings()
        for p in self.boxes:
            m = json.loads(p.read_text()); f = CODERS / self.coder / f"{m['image']}_{m['box']}.json"; mark = " "
            if f.exists(): d = json.loads(f.read_text()); mark = "✓" if d.get("done") else "·"
            n_meth = len(self.methods_for(m)); n_rated = sum((m["image"], m["box"], me) in rt for me in self.methods_for(m))
            self.listbox.insert("end", f"{mark} {m['image']:9s} {m['box']}  {m['side_um']:.0f}µm  {'CTRL ' if m['kind']=='control' else ''}{'' if m.get('set')=='rating' else '(practice) '}rated {n_rated}/{n_meth}")
    def methods_for(self, meta): return [m for m in METHODS if m != "v5" or "v5_polys" in meta]   # empty outline = a result to rate, not a missing method
    def load_box(self, i):
        if hasattr(self, "meta") and self.dirty: self.save(quiet=True)
        i = max(0, min(i, len(self.boxes) - 1)); self.i = i; self.meta = json.loads(self.boxes[i].read_text()); stem = self.boxes[i].stem
        z = np.load(DATA / f"{stem}.npz"); self.img = {k: z[k] for k in ("cy3", "dapi", "fitc")}
        self.masks = {k: np.array(Image.open(DATA / f"{stem}_{k}.png")) > 0 for k in ("gate3", "gate2", "gate15", "veto", "v4") if (DATA / f"{stem}_{k}.png").exists()}
        self.contour_cache = {}
        f = self.coder_file(); d = json.loads(f.read_text()) if f.exists() else {}
        self.dots = [dict(x=float(q["x"]), y=float(q["y"]), label=q.get("label", "deposit")) for q in d.get("dots", [])]; self.done_var.set(bool(d.get("done", False))); self.undo = []; self.dirty = False
        self.spot_i = -1; self.listbox.selection_clear(0, "end"); self.listbox.selection_set(i); self.listbox.see(i)
        # blind order of methods for this coder and box
        self.order = self.methods_for(self.meta)
        rt = self.load_ratings(); self.m_i = next((k for k, m in enumerate(self.order) if (self.meta["image"], self.meta["box"], m) not in rt), 0)
        self.info.delete("1.0", "end"); self.info.insert("end", f"{self.meta['image']} {self.meta['box']} — {self.meta['side_um']:.0f} µm box\n{self.meta['label']}")
        self.zoom = min(1.0, 900 / max(self.img['cy3'].shape)); self.overlay_on = True; self.mode_changed(); self.spot_txt.delete("1.0", "end")
    # ---------------- rendering ----------------
    def base_rgb(self):
        hi = float(self.hi_p.get()); ch = self.channel.get()
        if ch == "merge": r, g, b = stretch(self.img["cy3"], 1, hi), stretch(self.img["fitc"], 1, hi), stretch(self.img["dapi"], 1, hi)
        else: v = stretch(self.img[ch], 1, hi); r, g, b = v, v, v
        rgb = np.dstack([r, g, b])
        if self.mode.get() != "Rate":
            if self.show_gate3.get() and "gate3" in self.masks: m = self.masks["gate3"]; rgb[m] = 0.55 * rgb[m] + 0.45 * np.array([1, 1, 0])
            if self.show_veto.get() and "veto" in self.masks: m = self.masks["veto"]; rgb[m] = 0.45 * rgb[m] + 0.55 * np.array([1, 0, 1])
        return Image.fromarray((rgb * 255).astype(np.uint8))
    def render(self):
        im = self.base_rgb(); H, W = im.size[1], im.size[0]; z = self.zoom
        if abs(z - 1) > 1e-6: im = im.resize((max(1, int(W * z)), max(1, int(H * z))), Image.BILINEAR)
        self.photo = ImageTk.PhotoImage(im); self.canvas.delete("all"); self.canvas.create_image(0, 0, anchor="nw", image=self.photo); self.canvas.configure(scrollregion=(0, 0, im.size[0], im.size[1])); self.redraw()
    def method_mask(self, m):
        mk = self.masks
        if m == "v4": return mk.get("v4")
        if m == "gate1.5": return mk.get("gate15")
        if m == "gate2": return mk.get("gate2")
        if m == "gate3": return mk.get("gate3")
        if m == "v6": return (mk["v4"] & ~mk["veto"]) if "v4" in mk and "veto" in mk else None
        return None
    def method_contours(self, m):
        if m in self.contour_cache: return self.contour_cache[m]
        if m == "v5": polys = [np.array(p, dtype=np.float32) for p in self.meta.get("v5_polys", []) if len(p) >= 3]
        else:
            mask = self.method_mask(m)
            if mask is None: polys = []
            else:
                cs, _ = cv2.findContours(mask.astype(np.uint8), cv2.RETR_LIST, cv2.CHAIN_APPROX_SIMPLE); polys = [c[:, 0, :].astype(np.float32) for c in cs if len(c) >= 3]
        self.contour_cache[m] = polys; return polys
    def redraw(self):
        self.canvas.delete("ov"); z = self.zoom; mode = self.mode.get()
        if mode == "Rate":
            if self.overlay_on and self.order:
                m = self.order[self.m_i]
                for poly in self.method_contours(m): self.canvas.create_polygon([c * z for p in poly for c in p], outline=OUTLINE, fill="", width=1.5, tags="ov")
            self.update_rate_panel()
        else:
            if self.show_v4.get():
                for poly in self.meta.get("v4_polys", []):
                    if len(poly) >= 3: self.canvas.create_polygon([c * z for p in poly for c in p], outline="#ff5050", fill="", width=1, tags="ov")
            if self.show_other.get() or mode == "Compare":
                for o in self.other_dots(): r = 7; self.canvas.create_oval(o["x"] * z - r, o["y"] * z - r, o["x"] * z + r, o["y"] * z + r, outline=OTHER, width=2, tags="ov")
            for d in self.dots: r = 6; self.canvas.create_oval(d["x"] * z - r, d["y"] * z - r, d["x"] * z + r, d["y"] * z + r, fill=COLORS[d["label"]], outline="white", width=1.5, tags="ov")
            if mode == "Practice" and 0 <= self.spot_i < len(self.meta.get("spots", [])):
                s = self.meta["spots"][self.spot_i]; r = s["r"] * z; self.canvas.create_oval(s["x"] * z - r, s["y"] * z - r, s["x"] * z + r, s["y"] * z + r, outline="#00ff9d", width=3, tags="ov")
        n_dep = sum(d["label"] == "deposit" for d in self.dots); self.count_lbl.configure(text=f"deposits: {n_dep}    uncertain: {len(self.dots) - n_dep}")
        self.status.configure(text=f"{self.meta['image']} {self.meta['box']}   zoom {self.zoom:.2f}   mode {mode}   box {self.i+1}/{len(self.boxes)}   {'overlay ON' if self.overlay_on else 'overlay OFF'}")
        self.update_compare()
    def set_zoom(self, z): self.zoom = max(0.2, min(4.0, z)); self.render()
    def toggle_overlay(self): self.overlay_on = not self.overlay_on; self.redraw()
    # ---------------- rating ----------------
    def goto_method(self, k):
        if not self.order: return
        self.m_i = k % len(self.order); self.overlay_on = True; self.contour_cache.pop(None, None); self.redraw()
    def update_rate_panel(self):
        rt = self.load_ratings(); cur = self.order[self.m_i]
        self.rate_lbl.configure(text=f"{METHOD_NAMES.get(cur, cur)}  ({self.m_i+1} of {len(self.order)})" + ("\nThis method outlines NOTHING in this box" if not self.method_contours(cur) else ""))
        done = [(m, rt[(self.meta["image"], self.meta["box"], m)]["verdict"]) for m in self.order if (self.meta["image"], self.meta["box"], m) in rt]
        n_boxes_done = sum(all((json.loads(p.read_text())["image"], json.loads(p.read_text())["box"], m) in rt for m in self.methods_for(json.loads(p.read_text()))) for p in self.boxes)
        self.rate_prog.configure(text=f"this box: {len(done)}/{len(self.order)} rated      boxes complete: {n_boxes_done}/{len(self.boxes)}")
        self.rated_txt.delete("1.0", "end"); self.rated_txt.insert("end", "\n".join(f"{l}: {v}" for l, v in done) or "nothing rated yet in this box")
    def verdict(self, key):
        if self.mode.get() != "Rate" or not self.order: return
        m = self.order[self.m_i]; f = self.ratings_file(); new = f.exists() is False
        rows = [r for r in (csv.DictReader(open(f)) if f.exists() else []) if not (r["image"] == self.meta["image"] and r["box"] == self.meta["box"] and r["method"] == m)]
        rows.append(dict(coder=self.coder, image=self.meta["image"], box=self.meta["box"], method=m, verdict=VERDICTS[key], comment=self.comment.get().strip(), saved_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), label=self.meta["label"], side_um=self.meta["side_um"]))
        tmp = f.with_suffix(".tmp")
        with open(tmp, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        tmp.replace(f)
        self.comment.delete(0, "end"); self.canvas.focus_set(); self.refresh_list(); self.listbox.selection_set(self.i)
        rt = self.load_ratings(); nxt = next((k for k, mm in enumerate(self.order) if (self.meta["image"], self.meta["box"], mm) not in rt), None)
        if nxt is None:
            self.update_rate_panel(); self.status.configure(text=f"box {self.meta['image']} {self.meta['box']} fully rated — moving to the next box"); self.after(600, lambda: self.load_box(self.i + 1) if self.i + 1 < len(self.boxes) else None)
        else: self.goto_method(nxt)
    # ---------------- clicks / dots ----------------
    def click(self, e, label):
        if self.mode.get() in ("Compare", "Rate"): self.canvas.focus_set(); return
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
    def _mark_dirty(self): self.dirty = True
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
        self.mode.set("Practice"); self.mode_changed(); self.spot_i = k % len(spots); s = spots[self.spot_i]
        self.spot_txt.delete("1.0", "end"); self.spot_txt.insert("end", f"Spot {self.spot_i+1} of {len(spots)}  ({s['area_um2']} µm²)\nLook at the green circle. Deposit or not? Decide, then Reveal.")
        z = self.zoom; bb = self.canvas.bbox("all") or (0, 0, 1, 1); self.canvas.xview_moveto(max(0, (s["x"] * z - self.canvas.winfo_width() / 2) / max(1, bb[2]))); self.canvas.yview_moveto(max(0, (s["y"] * z - self.canvas.winfo_height() / 2) / max(1, bb[3]))); self.redraw()
    def reveal(self):
        spots = self.meta.get("spots", [])
        if not (0 <= self.spot_i < len(spots)): return
        s = spots[self.spot_i]; self.spot_txt.delete("1.0", "end"); self.spot_txt.insert("end", f"Spot {self.spot_i+1}/{len(spots)} — {s['category'].upper()}\n{s['text']}\nnumbers: area {s['area_um2']} µm², Cy3 {s['cy3_contrast']}× background, DAPI {s['dapi_contrast']}×, survives 3× gate: {int(100*s['gate3_survival'])}%")
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
        if not other: self.cmp_txt.insert("end", "No other coder has saved dots for this box yet."); return
        px = self.meta["px_um"]; used = set(); matched = 0
        for d in mine:
            best = None
            for j, o in enumerate(other):
                if j in used: continue
                dist = math.hypot(d["x"] - o["x"], d["y"] - o["y"]) * px
                if dist <= MATCH_UM and (best is None or dist < best[0]): best = (dist, j)
            if best: used.add(best[1]); matched += 1
        who = other[0]["coder"]; self.cmp_txt.insert("end", f"Compare with {who} (15 µm match): agreed {matched}, only you {len(mine) - matched}, only {who} {len(other) - matched}")
if __name__ == "__main__":
    app = Coder()
    if "--selftest" in sys.argv: app.after(1500, app.destroy)
    app.mainloop()
