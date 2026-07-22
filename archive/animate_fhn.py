"""
animate_fhn.py  —  Stochastic FHN: excitability + full loop animation
No LaTeX required (Text + Unicode throughout).

Acts
----
1. Phase portrait
2. Deterministic excitability: sub-threshold vs supra-threshold perturbation
3. Stochastic full loop with 4 annotated phases
4. v(t) colour-coded inset

Usage:
    cd "/Users/solomonwilliams/Reduced Neuronal Model Project"
    manim -pql animate_fhn.py FHNScene    # fast preview
    manim -pqh animate_fhn.py FHNScene    # high quality
"""

from manim import *
import numpy as np
from scipy.optimize import brentq

# ── Parameters ────────────────────────────────────────────────────────────────
I_EXT = -0.1
A, B  = 0.7, 0.8
EPS   = 0.12
SIGMA = 0.35
DT    = 2e-3
SEED  = 21

def fp_eq(v): return v - v**3/3 - (v + A)/B + I_EXT
V_FP  = brentq(fp_eq, -2.5, -1.0)
W_FP  = (V_FP + A) / B
V_SAD = brentq(lambda v: v - v**3/3 - W_FP + I_EXT, -1.0, 0.5)  # saddle on middle branch


# ── Pre-compute trajectories ──────────────────────────────────────────────────

def det_traj(v0, w0=W_FP, max_steps=80000):
    """Deterministic trajectory, runs until returned to fixed point."""
    v, w = v0, w0
    vs, ws = [v], [w]
    fired = False
    for _ in range(max_steps):
        v = v + (v - v**3/3 - w + I_EXT)*DT
        w = w + EPS*(v + A - B*w)*DT
        vs.append(v); ws.append(w)
        if v >= 1.2:
            fired = True
        if len(vs) > 300 and abs(v - V_FP) < 0.005 and abs(w - W_FP) < 0.005:
            break
    return np.array(vs), np.array(ws)

def stoch_traj(sigma, seed=SEED, max_steps=20000):
    """Stochastic trajectory, runs until full loop complete."""
    rng = np.random.default_rng(seed)
    sqrt_dt = np.sqrt(DT)
    v, w = V_FP, W_FP
    vs, ws = [v], [w]
    fired = False
    for _ in range(max_steps):
        v = v + (v - v**3/3 - w + I_EXT)*DT + sigma*sqrt_dt*rng.standard_normal()
        w = w + EPS*(v + A - B*w)*DT
        vs.append(v); ws.append(w)
        if v >= 1.2:
            fired = True
        if fired and abs(v - V_FP) < 0.005 and abs(w - W_FP) < 0.005:
            break
    return np.array(vs), np.array(ws)

# Sub-threshold: perturb to v=-0.65, spirals back to FP
V_SUB, W_SUB = det_traj(-0.65)

# Supra-threshold: perturb to v=-0.45, fires full loop back to FP
V_SUP, W_SUP = det_traj(-0.45)
SUP_FIRE  = int(np.where(V_SUP >= 1.2)[0][0])
SUP_PEAK  = SUP_FIRE + int(np.argmax(V_SUP[SUP_FIRE:]))
_ap       = np.where((np.arange(len(V_SUP)) > SUP_PEAK) & (V_SUP < 1.0))[0]
SUP_RFOLD = int(_ap[0]) if len(_ap) else len(V_SUP)-1
SUP_END   = len(V_SUP) - 1

# Stochastic full loop
V_S, W_S  = stoch_traj(SIGMA)
FIRE_IDX  = int(np.where(V_S >= 1.2)[0][0])
PEAK_IDX  = FIRE_IDX + int(np.argmax(V_S[FIRE_IDX:]))
_ap2      = np.where((np.arange(len(V_S)) > PEAK_IDX) & (V_S < 1.0))[0]
RFOLD_IDX = int(_ap2[0]) if len(_ap2) else len(V_S)-1


# ── Colours ───────────────────────────────────────────────────────────────────
C_BG    = "#0f1117"
C_AX    = "#dddddd"
C_VNULL = "#f0c040"
C_WNULL = "#60b8f0"
C_MFLD  = "#80e8a0"
C_FOLD  = "#ff9944"
C_SUB   = "#60b8f0"   # sub-threshold: blue
C_SUP   = "#ff9944"   # supra-threshold: orange
C_PHASE = {
    "noise":   "#ff6b6b",
    "jump_r":  "#ff9944",
    "drift_r": "#c084fc",
    "jump_l":  "#60b8f0",
}


# ── Helper: animate a trajectory segment ─────────────────────────────────────

def animate_seg(scene, ax, VS, WS, i0, i1, color, n_chunks=50, speed=0.012):
    chunk = max(1, (i1 - i0) // n_chunks)
    dot = Dot(ax.c2p(float(np.clip(VS[i0], -2.28, 2.28)),
                     float(np.clip(WS[i0], -1.38, 1.38))),
              radius=0.09, color=color)
    scene.add(dot)
    for i in range(i0 + chunk, i1 + 1, chunk):
        vp = float(np.clip(VS[i], -2.28, 2.28))
        wp = float(np.clip(WS[i], -1.38, 1.38))
        seg = VMobject(color=color, stroke_width=2.5, stroke_opacity=0.8)
        pts = [ax.c2p(float(np.clip(VS[j], -2.28, 2.28)),
                      float(np.clip(WS[j], -1.38, 1.38)))
               for j in range(max(i0, i-chunk), i+1)]
        seg.set_points_as_corners(pts if len(pts) > 1 else pts*2)
        scene.play(dot.animate.move_to(ax.c2p(vp, wp)),
                   Create(seg), run_time=speed, rate_func=linear)
    return dot


class FHNScene(Scene):

    def construct(self):
        self.camera.background_color = C_BG
        ax = self._axes()
        self.play(Create(ax), run_time=0.8)
        self._build_portrait(ax)
        self._act2_excitability(ax)
        self._act3_stochastic(ax)
        self._timeseries(ax)
        self.wait(2)

    # ── Axes ──────────────────────────────────────────────────────────────────

    def _axes(self):
        ax = Axes(
            x_range=[-2.3, 2.3, 1], y_range=[-1.4, 1.4, 0.5],
            x_length=9, y_length=6,
            axis_config={"color": C_AX, "stroke_width": 1.5,
                         "include_tip": True, "tip_length": 0.18,
                         "include_numbers": False},
        ).shift(DOWN*0.1)
        ax.add(Text("v  (voltage)",  font_size=20, color=C_AX)
                 .next_to(ax.x_axis.get_end(), RIGHT, buff=0.08))
        ax.add(Text("w  (recovery)", font_size=20, color=C_AX)
                 .next_to(ax.y_axis.get_end(), UP, buff=0.08))
        return ax

    # ── Act 1: Phase portrait ─────────────────────────────────────────────────

    def _build_portrait(self, ax):
        title = Text("FitzHugh–Nagumo  |  Phase Plane",
                     font_size=24, color=C_AX).to_edge(UP, buff=0.2)
        self.play(FadeIn(title))

        vnull = ax.plot(lambda v: v - v**3/3 + I_EXT,
                        x_range=[-2.3, 2.3], color=C_VNULL, stroke_width=2.5)
        vnull_lbl = Text("v̇ = 0", font_size=19, color=C_VNULL)\
                      .next_to(ax.c2p(1.4, 1.4 - 1.4**3/3 + I_EXT), UR, buff=0.08)
        wnull = ax.plot(lambda v: (v + A)/B,
                        x_range=[-2.3, 2.3], color=C_WNULL, stroke_width=2.5)
        wnull_lbl = Text("ẇ = 0", font_size=19, color=C_WNULL)\
                      .next_to(ax.c2p(0.5, (0.5+A)/B), UR, buff=0.08)
        self.play(Create(vnull), FadeIn(vnull_lbl), run_time=1.0)
        self.play(Create(wnull), FadeIn(wnull_lbl), run_time=0.8)

        # Slow manifold
        v_l = np.linspace(-2.28, -1.02, 300)
        mfld = VMobject(color=C_MFLD, stroke_width=5, stroke_opacity=0.55)
        mfld.set_points_as_corners([ax.c2p(v, v - v**3/3 + I_EXT) for v in v_l])
        mfld_lbl = Text("slow manifold", font_size=16, color=C_MFLD)\
                     .next_to(ax.c2p(-1.9, -1.0), DOWN, buff=0.06)
        self.play(Create(mfld), FadeIn(mfld_lbl), run_time=0.7)

        # Fixed point
        fp = Dot(ax.c2p(V_FP, W_FP), radius=0.11, color="#ffffff")
        fp_lbl = Text("stable fixed point", font_size=15, color="#ffffff")\
                   .next_to(fp, UR, buff=0.1)
        self.play(GrowFromCenter(fp), FadeIn(fp_lbl))

        # Saddle point on middle branch
        sad = Dot(ax.c2p(V_SAD, W_FP), radius=0.09, color="#aaaaaa")
        sad_lbl = Text("saddle\n(threshold)", font_size=14, color="#aaaaaa")\
                    .next_to(sad, UP, buff=0.1)
        self.play(GrowFromCenter(sad), FadeIn(sad_lbl), run_time=0.5)

        # Fold markers
        for vf, name in [(-1.0, "left fold"), (1.0, "right fold")]:
            wf = vf - vf**3/3 + I_EXT
            fd = Dot(ax.c2p(vf, wf), radius=0.09, color=C_FOLD)
            fl = Text(name, font_size=14, color=C_FOLD)\
                   .next_to(fd, DOWN if vf < 0 else UP, buff=0.1)
            self.play(GrowFromCenter(fd), FadeIn(fl), run_time=0.4)

        # Param box
        params = VGroup(
            Text(f"ε = {EPS}", font_size=18, color=C_AX),
            Text(f"σ = {SIGMA}", font_size=18, color=C_PHASE["noise"]),
            Text(f"I = {I_EXT}", font_size=18, color=C_AX),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.09).to_corner(DR, buff=0.35)
        box = SurroundingRectangle(params, color=C_AX, fill_color=C_BG,
                                   fill_opacity=0.8, stroke_width=1, buff=0.12)
        self.play(FadeIn(box), FadeIn(params))
        self.wait(0.6)

    # ── Act 2: Deterministic excitability ─────────────────────────────────────

    def _act2_excitability(self, ax):
        hdr = Text("σ = 0  |  Deterministic excitability",
                   font_size=22, color=C_AX).to_corner(UL, buff=0.28)
        sub_hdr = Text("The saddle point is the threshold —\n"
                       "cross it and the system must fire a full spike.",
                       font_size=16, color=GRAY).next_to(hdr, DOWN, buff=0.1,
                                                          aligned_edge=LEFT)
        self.play(FadeIn(hdr), FadeIn(sub_hdr))
        self.wait(0.5)

        # ── Sub-threshold ──
        sub_lbl = Text("below threshold  →  returns to rest",
                       font_size=18, color=C_SUB).to_corner(UL, buff=0.28)
        self.play(Transform(hdr, sub_lbl), FadeOut(sub_hdr), run_time=0.4)

        # Mark starting point
        start_sub = Dot(ax.c2p(-0.65, W_FP), radius=0.10, color=C_SUB)
        arr_sub = Arrow(ax.c2p(-0.65, W_FP + 0.25), ax.c2p(-0.65, W_FP + 0.05),
                        color=C_SUB, stroke_width=2, buff=0.02)
        self.play(GrowFromCenter(start_sub), Create(arr_sub), run_time=0.4)
        self.play(FadeOut(arr_sub), run_time=0.2)

        dot_sub = animate_seg(self, ax, V_SUB, W_SUB, 0, len(V_SUB)-1,
                               C_SUB, n_chunks=60, speed=0.018)
        self.wait(0.3)
        self.play(FadeOut(dot_sub), FadeOut(start_sub))

        # ── Supra-threshold ──
        sup_lbl = Text("above threshold  →  must fire a full spike",
                       font_size=18, color=C_SUP).to_corner(UL, buff=0.28)
        self.play(Transform(hdr, sup_lbl), run_time=0.4)

        start_sup = Dot(ax.c2p(-0.45, W_FP), radius=0.10, color=C_SUP)
        arr_sup = Arrow(ax.c2p(-0.45, W_FP + 0.25), ax.c2p(-0.45, W_FP + 0.05),
                        color=C_SUP, stroke_width=2, buff=0.02)
        self.play(GrowFromCenter(start_sup), Create(arr_sup), run_time=0.4)
        self.play(FadeOut(arr_sup), run_time=0.2)

        # Animate full supra-threshold loop in phases
        phases_sup = [
            (0,           SUP_FIRE,  C_SUP,              0.018, None),
            (SUP_FIRE,    SUP_PEAK,  C_PHASE["jump_r"],  0.010, "SPIKE!"),
            (SUP_PEAK,    SUP_RFOLD, C_PHASE["drift_r"], 0.020, None),
            (SUP_RFOLD,   SUP_END,   C_PHASE["jump_l"],  0.018, None),
        ]
        dot_sup = Dot(ax.c2p(float(np.clip(V_SUP[0],-2.28,2.28)),
                             float(np.clip(W_SUP[0],-1.38,1.38))),
                      radius=0.09, color=C_SUP)
        self.add(dot_sup)

        for (i0, i1, col, spd, flash_lbl) in phases_sup:
            chunk = max(1, (i1 - i0) // 50)
            for i in range(i0 + chunk, i1 + 1, chunk):
                vp = float(np.clip(V_SUP[i], -2.28, 2.28))
                wp = float(np.clip(W_SUP[i], -1.38, 1.38))
                seg = VMobject(color=col, stroke_width=2.5, stroke_opacity=0.8)
                pts = [ax.c2p(float(np.clip(V_SUP[j],-2.28,2.28)),
                              float(np.clip(W_SUP[j],-1.38,1.38)))
                       for j in range(max(i0, i-chunk), i+1)]
                seg.set_points_as_corners(pts if len(pts) > 1 else pts*2)
                self.play(dot_sup.animate.move_to(ax.c2p(vp, wp)),
                          Create(seg), run_time=spd, rate_func=linear)

            if flash_lbl:
                fl = Flash(dot_sup.get_center(), color=col,
                           line_length=0.3, num_lines=12, flash_radius=0.4)
                sl = Text(flash_lbl, font_size=24, color=col, weight=BOLD)\
                       .next_to(dot_sup, UP, buff=0.18)
                self.play(fl, FadeIn(sl), run_time=0.35)
                self.wait(0.25)
                self.play(FadeOut(sl), run_time=0.15)

        rest_lbl = Text("returns to rest  ✓", font_size=18, color=C_MFLD)\
                     .to_corner(UL, buff=0.28)
        self.play(Transform(hdr, rest_lbl), run_time=0.4)
        self.wait(0.7)
        self.play(FadeOut(hdr), FadeOut(dot_sup), FadeOut(start_sup))

    # ── Act 3: Stochastic full loop ───────────────────────────────────────────

    def _act3_stochastic(self, ax):
        phases = [
            (0,          FIRE_IDX,  C_PHASE["noise"],   0.018, "① noise-driven escape",      None),
            (FIRE_IDX,   PEAK_IDX,  C_PHASE["jump_r"],  0.010, "② fast jump → right branch", "SPIKE!"),
            (PEAK_IDX,   RFOLD_IDX, C_PHASE["drift_r"], 0.022, "③ slow drift down",           None),
            (RFOLD_IDX,  len(V_S)-1,C_PHASE["jump_l"],  0.018, "④ fast jump back → rest",     None),
        ]

        lbl = Text(f"σ = {SIGMA}  |  Stochastic: noise crosses the threshold",
                   font_size=20, color=C_PHASE["noise"]).to_corner(UL, buff=0.28)
        self.play(FadeIn(lbl))

        dot = Dot(ax.c2p(float(np.clip(V_S[0],-2.28,2.28)),
                         float(np.clip(W_S[0],-1.38,1.38))),
                  radius=0.10, color=C_PHASE["noise"])
        self.play(GrowFromCenter(dot))

        phase_lbl = Text("", font_size=19).to_corner(UL, buff=0.55)

        for (i0, i1, col, spd, label, flash_str) in phases:
            new_pl = Text(label, font_size=18, color=col)\
                       .next_to(lbl, DOWN, buff=0.08, aligned_edge=LEFT)
            self.play(Transform(phase_lbl, new_pl), run_time=0.25)
            dot.set_color(col)

            chunk = max(1, (i1 - i0) // 55)
            for i in range(i0 + chunk, i1 + 1, chunk):
                vp = float(np.clip(V_S[i], -2.28, 2.28))
                wp = float(np.clip(W_S[i], -1.38, 1.38))
                seg = VMobject(color=col, stroke_width=2.8, stroke_opacity=0.8)
                pts = [ax.c2p(float(np.clip(V_S[j],-2.28,2.28)),
                              float(np.clip(W_S[j],-1.38,1.38)))
                       for j in range(max(i0, i-chunk), i+1)]
                seg.set_points_as_corners(pts if len(pts) > 1 else pts*2)
                self.play(dot.animate.move_to(ax.c2p(vp, wp)),
                          Create(seg), run_time=spd, rate_func=linear)

            if flash_str:
                fl = Flash(dot.get_center(), color=col,
                           line_length=0.3, num_lines=12, flash_radius=0.42)
                sl = Text(flash_str, font_size=26, color=col, weight=BOLD)\
                       .next_to(dot, UP, buff=0.18)
                self.play(fl, FadeIn(sl), run_time=0.35)
                self.wait(0.25)
                self.play(FadeOut(sl), run_time=0.15)

        rest = Text("returns to rest  ✓", font_size=18, color=C_MFLD)\
                 .next_to(lbl, DOWN, buff=0.08, aligned_edge=LEFT)
        self.play(Transform(phase_lbl, rest), run_time=0.35)
        self.wait(0.8)
        self.play(FadeOut(lbl), FadeOut(phase_lbl), FadeOut(dot))

    # ── v(t) time series ──────────────────────────────────────────────────────

    def _timeseries(self, ax):
        n = len(V_S)
        t = np.arange(n) * DT

        ts = Axes(
            x_range=[0, t[-1]*1.05, t[-1]/4],
            y_range=[-2.1, 2.1, 1],
            x_length=3.8, y_length=1.6,
            axis_config={"color": GRAY, "stroke_width": 1.1,
                         "include_tip": False, "include_numbers": False},
        ).to_corner(DR, buff=0.3).shift(UP*0.1)

        bg    = BackgroundRectangle(ts, color=C_BG, fill_opacity=0.9, buff=0.1)
        title = Text("v(t)", font_size=16, color=GRAY).next_to(ts, UP, buff=0.05)

        thresh = Line(ts.c2p(0, 1.0), ts.c2p(t[-1], 1.0),
                      color=C_FOLD, stroke_width=1.2).set_opacity(0.65)
        thresh_lbl = Text("threshold", font_size=11, color=C_FOLD)\
                       .next_to(ts.c2p(0, 1.0), LEFT, buff=0.05)

        self.play(FadeIn(bg), Create(ts), FadeIn(title), run_time=0.5)
        self.play(Create(thresh), FadeIn(thresh_lbl), run_time=0.3)

        for (i0, i1, col) in [
            (0,          FIRE_IDX,  C_PHASE["noise"]),
            (FIRE_IDX,   PEAK_IDX,  C_PHASE["jump_r"]),
            (PEAK_IDX,   RFOLD_IDX, C_PHASE["drift_r"]),
            (RFOLD_IDX,  n-1,       C_PHASE["jump_l"]),
        ]:
            i1 = min(i1, n-1)
            if i1 <= i0: continue
            seg = VMobject(color=col, stroke_width=2.0)
            pts = [ts.c2p(t[i], float(np.clip(V_S[i], -2.1, 2.1)))
                   for i in range(i0, i1+1)]
            seg.set_points_as_corners(pts)
            self.play(Create(seg), run_time=0.5, rate_func=linear)
