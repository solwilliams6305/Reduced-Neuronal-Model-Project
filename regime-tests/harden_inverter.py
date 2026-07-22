"""
harden_inverter.py — three hardenings of the proximity-to-bifurcation inverter, the open items in
INVERSION_MVP_ROADMAP.md Phase 1+/2 ("still open: colored-noise axis; off-grid fitting; near-threshold
weighting for the class call"):

  (A) OFF-GRID (interpolated) fitting — kills the grid-snap artifact. The atlas stores a normalised
      quantile function Q(p) on a coarse parameter grid; Q is smooth in p, so we LINEARLY INTERPOLATE
      the quantile curves between grid nodes and minimise W1 on a dense continuous axis (bilinear for
      the 2D (nu,b) and (nu,tau_s) grids). param_hat becomes continuous instead of snapping to nodes.

  (B) NEAR-THRESHOLD WEIGHTING for the class call — the SNIC/Hopf signature is strong NEAR threshold
      and washes out at high drive (regular Type-I firing -> looks Hopf). A single high-drive sweep is
      an unreliable class witness. aggregate_class() pools per-sweep soft evidence across a cell's
      sweep set, weighting each sweep by a near-threshold factor 1/(1+(nu_hat/nu0)^2) and its sample
      size. The cell-level Type-I/Type-II call is then dominated by its most informative (near-edge)
      sweeps.

  (C) COLORED-NOISE confound axis — real input noise is not white. colored_snic_intervals(nu, tau_s)
      drives the SNIC QIF with OU noise (matching qif_validation's convention, white as tau_s->0); a
      2D (nu, tau_s) atlas + joint fit_colored() recovers nu_hat under colored noise that otherwise
      biases it, reports tau_s_hat, and the (nu, tau_s) valley where the two are not separable.

All three build on the existing atlas machinery in bifurcation_classifier (qfun / W1 / classify3 /
fit_adapt) and the SNIC QIF used for the adaptation axis — nothing downstream is modified.

    python3 harden_inverter.py coloratlas   # build + cache the (nu, tau_s) colored atlas (~10 s)
    python3 harden_inverter.py interp       # off-grid vs grid-snap recovery test
    python3 harden_inverter.py classvote    # near-threshold weighted class aggregation test
    python3 harden_inverter.py colored      # colored-noise bias vs joint-fit recovery test
    python3 harden_inverter.py harden       # all three -> figures/harden_inverter.png
"""
from __future__ import annotations
import os
import sys
from functools import lru_cache
import numpy as np

import bifurcation_classifier as bc
from bifurcation_classifier import qfun, W1, NQ, skew
from offcritical_phase_edge import J_of_nu

HERE = os.path.dirname(os.path.abspath(__file__))
COLOR = os.path.join(HERE, "_colored_atlas.npz")


# =================================================================================================
# (A) OFF-GRID / interpolated fitting
# =================================================================================================
def _dense_1d(params, Q, R=400):
    """Linearly interpolate quantile curves Q (shape (n, NQ)) across the param grid to R points."""
    f = np.linspace(params[0], params[-1], R)
    Qf = np.empty((R, Q.shape[1]))
    for q in range(Q.shape[1]):
        Qf[:, q] = np.interp(f, params, Q[:, q])
    return f, Qf


def fit_interp_1d(qd, params, Q, R=400):
    """Continuous-parameter W1 fit: (param_hat, W1) on the interpolated atlas."""
    f, Qf = _dense_1d(params, Q, R)
    d = np.mean(np.abs(Qf - qd[None, :]), axis=1)
    j = int(np.argmin(d))
    return float(f[j]), float(d[j])


def classify3_interp(intervals, R=400):
    """Off-grid version of bc.classify3 -> (class, param_hat, W1, per)."""
    qd = qfun(intervals)
    per = {}
    for cls, (params, Q) in bc._atlas3().items():
        per[cls] = fit_interp_1d(qd, np.asarray(params, float), np.asarray(Q, float), R)
    best = min(per.items(), key=lambda kv: kv[1][1])
    return best[0], best[1][0], best[1][1], per


def _dense_2d(p1, p2, Q, R1=60, R2=40):
    """Bilinear-interpolate a (n1, n2, NQ) quantile atlas onto an (R1, R2, NQ) dense grid."""
    p1 = np.asarray(p1, float); p2 = np.asarray(p2, float); Q = np.asarray(Q, float)
    n1, n2, nq = Q.shape
    f1 = np.linspace(p1[0], p1[-1], R1); f2 = np.linspace(p2[0], p2[-1], R2)
    A = np.empty((R1, n2, nq))                       # interp along axis 0 (p1)
    for jb in range(n2):
        for q in range(nq):
            A[:, jb, q] = np.interp(f1, p1, Q[:, jb, q])
    B = np.empty((R1, R2, nq))                       # then along axis 1 (p2)
    for i in range(R1):
        for q in range(nq):
            B[i, :, q] = np.interp(f2, p2, A[i, :, q])
    return f1, f2, B


@lru_cache(maxsize=1)
def _dense_adapt(R1=60, R2=40):
    d = np.load(bc.ADAPT)
    f1, f2, B = _dense_2d(d["nu_g"], d["b_g"], d["Q"], R1, R2)
    return f1, f2, B


def fit_adapt_interp(intervals, joint=True):
    """Off-grid (nu, b) joint fit -> (nu_hat, b_hat, W1) on the bilinear-interpolated adapt atlas."""
    qd = qfun(intervals); f1, f2, B = _dense_adapt()
    if not joint:
        d = np.mean(np.abs(B[:, 0, :] - qd[None, :]), axis=1)      # b=0 column
        i = int(np.argmin(d)); return float(f1[i]), 0.0, float(d[i])
    d = np.mean(np.abs(B - qd[None, None, :]), axis=2)
    i, j = np.unravel_index(int(np.argmin(d)), d.shape)
    return float(f1[i]), float(f2[j]), float(d[i, j])


# =================================================================================================
# (B) NEAR-THRESHOLD weighted class aggregation
# =================================================================================================
TYPE_I = ("SNIC", "Homoclinic")


def _soft_class(per, tau=None):
    """Per-sweep soft class probabilities from {class: (param, W1)} via softmax(-W1/tau)."""
    classes = list(per.keys()); w = np.array([per[c][1] for c in classes])
    if tau is None:
        tau = max(1e-3, 0.5 * (w.max() - w.min()))
    p = np.exp(-(w - w.min()) / tau); p /= p.sum()
    return dict(zip(classes, p))


def aggregate_class(sweeps, nu0=0.8, n_min=8, interp=True, weighting="nearthr"):
    """Cell-level class call from a list of per-sweep ISI arrays.

    weighting='nearthr' down-weights high-drive sweeps (where SNIC<->Hopf separability collapses) by
    1/(1+(nu_hat/nu0)^2) and small samples by min(1, n/20); 'uniform' gives every sweep equal weight
    (the baseline). Returns (call, P, detail): call in {'Type-I','Type-II'}, P = {class: posterior},
    detail = list of per-sweep (nu_hat, class_argmin, weight)."""
    classify = classify3_interp if interp else (lambda x: bc.classify3(x))
    fitc = fit_adapt_interp if interp else (lambda x, joint: bc.fit_adapt(x, joint=joint))
    score = {}; detail = []; wtot = 0.0
    for isis in sweeps:
        isis = np.asarray(isis, float); isis = isis[isis > 0]
        if len(isis) < n_min:
            continue
        cls, p, w, per = classify(isis)
        nu_hat = fitc(isis, True)[0]                              # drive proxy (>=0 near edge)
        if weighting == "uniform":
            wt = 1.0
        else:
            w_thr = 1.0 / (1.0 + (max(nu_hat, 0.0) / nu0) ** 2)   # near-threshold emphasis
            wt = w_thr * min(1.0, len(isis) / 20.0)               # x sample-size weight
        soft = _soft_class(per)
        for c, pc in soft.items():
            score[c] = score.get(c, 0.0) + wt * pc
        wtot += wt; detail.append((nu_hat, cls, wt))
    if wtot == 0:
        return "none", {}, detail
    P = {c: v / wtot for c, v in score.items()}
    pI = sum(P.get(c, 0.0) for c in TYPE_I)
    return ("Type-I" if pI >= 0.5 else "Type-II"), P, detail


# =================================================================================================
# (C) COLORED-NOISE confound axis
# =================================================================================================
def colored_snic_intervals(nu, tau_s, sigma=0.4, v_th=14.0, v_reset=-14.0,
                           N=400, dt=2.5e-3, n_isi=9, seed=0):
    """Steady-state ISIs of a SNIC QIF driven by OU (colored) input noise of correlation time tau_s.

    dv = (v^2 + I) dt + xi dt,   dxi = -(xi/tau_s) dt + ou_a sqrt(dt) dW,  ou_a = sigma/tau_s
    (the qif_validation OU convention: integrated intensity -> white sigma dW as tau_s -> 0).
    First two ISIs per neuron (start transient) dropped, mirroring the adaptation atlas."""
    rng = np.random.default_rng(seed); I = nu * sigma ** (4.0 / 3.0)
    v = np.full(N, float(v_reset)); xi = np.zeros(N); last = np.full(N, np.nan)
    cnt = np.zeros(N, int); isis = []; sdt = np.sqrt(dt)
    ou_a = sigma / tau_s                                          # s*sqrt(2/tau_s), s=sigma/sqrt(2 tau_s)
    Tmax = (n_isi + 3) * sigma ** (-2.0 / 3.0) * max(J_of_nu(nu), 2.0)
    for k in range(int(Tmax / dt)):
        t = k * dt
        xi += -(xi / tau_s) * dt + ou_a * sdt * rng.standard_normal(N)
        v += (v * v + I) * dt + xi * dt
        fired = v > v_th
        if fired.any():
            ok = fired & (cnt >= 2) & np.isfinite(last)
            if ok.any():
                isis.extend((t - last[ok]).tolist())
            last[fired] = t; cnt[fired] += 1; v[fired] = v_reset
    return np.array(isis)


def build_colored_atlas(nu_g=(0.0, 0.4, 0.8, 1.2, 1.7, 2.3),
                        tau_g=(0.05, 0.3, 0.6, 1.0)):
    nu_g = np.array(nu_g, float); tau_g = np.array(tau_g, float)
    Q = np.empty((len(nu_g), len(tau_g), NQ))
    for i, nu in enumerate(nu_g):
        for jt, ts in enumerate(tau_g):
            Q[i, jt] = qfun(colored_snic_intervals(nu, ts, seed=30 + i * 7 + jt))
        print(f"  colored atlas nu={nu:+.1f} done")
    np.savez(COLOR, nu_g=nu_g, tau_g=tau_g, Q=Q)
    print(f"  saved {COLOR} (SNIC x colored-noise tau_s)")


@lru_cache(maxsize=1)
def _dense_colored(R1=60, R2=40):
    d = np.load(COLOR)
    f1, f2, B = _dense_2d(d["nu_g"], d["tau_g"], d["Q"], R1, R2)
    return f1, f2, B


def fit_colored(intervals, joint=True):
    """Off-grid (nu, tau_s) joint fit -> (nu_hat, tau_s_hat, W1). joint=False = white-only (tau_s min)."""
    qd = qfun(intervals); f1, f2, B = _dense_colored()
    if not joint:
        d = np.mean(np.abs(B[:, 0, :] - qd[None, :]), axis=1)      # smallest-tau row (white-ish)
        i = int(np.argmin(d)); return float(f1[i]), float(f2[0]), float(d[i])
    d = np.mean(np.abs(B - qd[None, None, :]), axis=2)
    i, j = np.unravel_index(int(np.argmin(d)), d.shape)
    return float(f1[i]), float(f2[j]), float(d[i, j])


# =================================================================================================
# validation entry points
# =================================================================================================
def demo_interp():
    """Show off-grid fitting removes grid-snap: recover continuous nu from SNIC data at off-node nu."""
    d = np.load(bc.ATLAS3); params = d["snic_nu"]
    print(f"  SNIC atlas nodes: {np.round(params, 3)}")
    truth = [0.27, 0.63, 1.05, 1.55, 2.05]                         # deliberately BETWEEN nodes
    rows = []
    for nu in truth:
        s = bc.snic_intervals(nu, N=4000, seed=int(100 * nu) + 1)
        qd = qfun(s)
        # grid (snap to nearest node)
        dg = np.array([W1(qd, d["snic_Q"][j]) for j in range(len(params))])
        nu_grid = float(params[int(np.argmin(dg))])
        nu_int, _ = fit_interp_1d(qd, np.asarray(params, float), np.asarray(d["snic_Q"], float))
        rows.append((nu, nu_grid, nu_int))
        print(f"  true ν={nu:.2f} | grid ν̂={nu_grid:.3f} (snapped) | off-grid ν̂={nu_int:.3f}")
    rows = np.array(rows)
    eg = np.abs(rows[:, 1] - rows[:, 0]); ei = np.abs(rows[:, 2] - rows[:, 0])
    print(f"  mean |error|: grid={eg.mean():.3f}  off-grid={ei.mean():.3f}  "
          f"({100*(1-ei.mean()/eg.mean()):.0f}% reduction)")
    return rows


def demo_classvote(K=80):
    """(1) The justification: at large samples the SNIC->Hopf distance W1_Hopf (the 'room for error')
    collapses with drive, so high-drive sweeps sit on the class boundary. (2) The payoff: a Monte
    Carlo over Type-I cells shows near-threshold weighting beats uniform weighting on the cell call.
    All sims are batched (vectorised across units) to stay within the per-call budget."""
    from allen_phase2_inversion import sweep_isis_batch
    drives = np.array([0.2, 0.5, 0.9, 1.4, 2.0, 2.8])

    # (1) population-level separability vs drive: ONE batched sim, M=60 particles per drive
    Mp = 60
    NU = np.repeat(drives, Mp)
    sw = sweep_isis_batch(NU, [0.4] * len(NU), [0.0] * len(NU), n_isi=16, seed=11)
    print("  separability vs drive (large-sample):  drive  CV   W1_Hopf  margin")
    sep = []
    for di, nu in enumerate(drives):
        isis = np.concatenate([sw[di * Mp + m] for m in range(Mp) if len(sw[di * Mp + m]) > 0])
        cls, p, w, per = classify3_interp(isis)
        ws, wh = per["SNIC"][1], per["Hopf"][1]; margin = (wh - ws) / (wh + ws)
        sep.append((nu, isis.std() / isis.mean(), wh, margin))
        print(f"    {nu:.1f}   {isis.std()/isis.mean():.2f}   {wh:.4f}   {margin:+.3f}")

    # (2) Monte Carlo over two sampling protocols, each cell's sweeps in ONE batched sim:
    #     balanced  = drives spread across the edge (the easy case);
    #     highdrive = mostly well above rheobase (1 near-threshold + 5 high) — the realistic hard case
    #                 where uniform voting is dragged toward Hopf and the weighting earns its keep.
    protocols = {"balanced": drives,
                 "highdrive": np.array([0.3, 1.6, 2.2, 2.8, 3.4, 4.0])}
    out = {}
    for proto, dr in protocols.items():
        NU2 = np.tile(dr, K)
        swc = sweep_isis_batch(NU2, [0.4] * len(NU2), [0.0] * len(NU2), n_isi=26, seed=777)
        res = {"uniform": [], "nearthr": []}
        for k in range(K):
            sweeps = [swc[k * len(dr) + i] for i in range(len(dr))]
            for mode in ("uniform", "nearthr"):
                call, P, _ = aggregate_class(sweeps, weighting=mode)
                res[mode].append((call == "Type-I", sum(P.get(c, 0.0) for c in TYPE_I)))
        out[proto] = {m: np.array(v, float) for m, v in res.items()}
        print(f"\n  Monte Carlo — {proto} sampling, {K} Type-I cells (6 sweeps each, ~26 ISIs):")
        for mode in ("uniform", "nearthr"):
            a = out[proto][mode]
            print(f"    {mode:8s}:  correct Type-I cell call = {100*a[:,0].mean():.0f}%   "
                  f"mean P(Type-I) = {a[:,1].mean():.2f}")
    return sep, out


def demo_colored():
    """Colored-noise data biases the white-atlas nu_hat; the joint (nu, tau_s) fit recovers it."""
    if not os.path.exists(COLOR):
        build_colored_atlas()
    truth = [0.4, 0.8, 1.2, 1.7]; tau_true = 0.6
    rows = []
    for nu in truth:
        data = colored_snic_intervals(nu, tau_true, seed=int(70 + 100 * nu))
        cv = data.std() / data.mean()
        nu_naive, _, _ = fit_colored(data, joint=False)
        nu_j, tau_j, _ = fit_colored(data, joint=True)
        rows.append((nu, nu_naive, nu_j, tau_j, cv))
        print(f"  true ν={nu:.2f} (τ_s={tau_true}, CV={cv:.2f}): "
              f"white-atlas ν̂={nu_naive:+.2f} | joint ν̂={nu_j:+.2f}, τ̂_s={tau_j:.2f}")
    rows = np.array(rows)
    en = np.abs(rows[:, 1] - rows[:, 0]).mean(); ej = np.abs(rows[:, 2] - rows[:, 0]).mean()
    print(f"  mean |ν error|: white-atlas={en:.3f}  joint={ej:.3f}  "
          f"({100*(1-ej/en):.0f}% reduction)")

    # τ_s tracking: τ̂_s should RISE with true τ_s (detect coloring), even if biased
    tau_track = []
    for ts in [0.05, 0.3, 0.6, 1.0]:
        w = colored_snic_intervals(0.8, ts, seed=5)
        tau_track.append((ts, fit_colored(w, True)[1]))
    print(f"  τ_s tracking (ν=0.8): {[ (t, round(th,2)) for t,th in tau_track ]}")

    # quantify the (ν, τ_s) degeneracy: W1 valley anisotropy (steep in ν, flat in τ_s)
    f1, f2, B = _dense_colored(); qd = qfun(colored_snic_intervals(0.8, 0.6, seed=7))
    d = np.mean(np.abs(B - qd[None, None, :]), axis=2)
    i, j = np.unravel_index(int(np.argmin(d)), d.shape)
    s_nu = d[:, j].max() - d[:, j].min(); s_tau = d[i, :].max() - d[i, :].min()
    print(f"  W1-valley anisotropy: sensitivity along ν={s_nu:.3f} vs τ_s={s_tau:.3f} "
          f"(ν {s_nu/s_tau:.0f}× steeper ⇒ τ_s only weakly identified)")
    return rows


def harden_figure():
    """Recompute the three hardening tests and render a 2x2 verification figure."""
    if not os.path.exists(COLOR):
        build_colored_atlas()
    print("== (A) off-grid fitting =="); interp = demo_interp()
    print("\n== (B) near-threshold class weighting =="); sep, mc = demo_classvote()
    print("\n== (C) colored-noise axis =="); col = demo_colored()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 11, "axes.edgecolor": "#334155", "axes.linewidth": 0.9})
    fig, ax = plt.subplots(2, 2, figsize=(13.0, 9.2))
    CG, CI, CU, CW = "#94a3b8", "#16a34a", "#dc2626", "#2563eb"

    # A: off-grid kills grid-snap
    a = ax[0, 0]
    a.plot([0, 2.2], [0, 2.2], "--", color="#cbd5e1", lw=1.2, label="ideal")
    a.plot(interp[:, 0], interp[:, 1], "o", color=CU, ms=9, mec="white", mew=1.0, label="grid (snapped)")
    a.plot(interp[:, 0], interp[:, 2], "s", color=CI, ms=9, mec="white", mew=1.0, label="off-grid (interp)")
    eg = np.abs(interp[:, 1] - interp[:, 0]).mean(); ei = np.abs(interp[:, 2] - interp[:, 0]).mean()
    a.set_xlabel(r"true $\nu$ (between atlas nodes)"); a.set_ylabel(r"recovered $\hat\nu$")
    a.set_title(fr"(A) Off-grid fitting: mean $|$err$|$ {eg:.3f}$\to${ei:.3f}")
    a.legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    # B: separability collapse with drive (the justification for weighting)
    b = ax[0, 1]; sep = np.array(sep)
    b.plot(sep[:, 0], sep[:, 2], "o-", color=CW, ms=7, mec="white", mew=1.0, label=r"$W_1$ to Hopf atlas")
    b.set_xlabel(r"drive $\nu$"); b.set_ylabel(r"$W_1$(data, Hopf) — room for error", color=CW)
    b.tick_params(axis="y", labelcolor=CW)
    b2 = b.twinx()
    b2.plot(sep[:, 0], sep[:, 1], "s--", color="#7c3aed", ms=6, label="ISI CV")
    b2.set_ylabel("ISI CV", color="#7c3aed"); b2.tick_params(axis="y", labelcolor="#7c3aed")
    b.set_title("(B) SNIC↔Hopf separability collapses at high drive")

    # C: class-vote MC — mean P(Type-I), uniform vs near-threshold, two protocols
    c = ax[1, 0]
    protos = ["balanced", "highdrive"]; x = np.arange(2); wbar = 0.36
    pu = [mc[p]["uniform"][:, 1].mean() for p in protos]
    pn = [mc[p]["nearthr"][:, 1].mean() for p in protos]
    c.bar(x - wbar/2, pu, wbar, color=CU, label="uniform weighting")
    c.bar(x + wbar/2, pn, wbar, color=CI, label="near-threshold weighting")
    c.axhline(0.5, color="#94a3b8", lw=1.0, ls=":"); c.text(1.25, 0.51, "coin-flip", fontsize=8, color="#64748b")
    for xi, (u, n) in enumerate(zip(pu, pn)):
        c.text(xi - wbar/2, u + 0.01, f"{u:.2f}", ha="center", fontsize=9)
        c.text(xi + wbar/2, n + 0.01, f"{n:.2f}", ha="center", fontsize=9)
    c.set_xticks(x); c.set_xticklabels(["balanced\nsampling", "high-drive\nsampling"])
    c.set_ylabel("mean P(Type-I) — posterior confidence"); c.set_ylim(0, 0.9)
    c.set_title("(C) Near-threshold weighting improves calibration")
    c.legend(fontsize=8.6, framealpha=0.95, edgecolor="#cbd5e1", loc="upper right")

    # D: colored-noise ν recovery, white-atlas vs joint (ν, τ_s)
    d = ax[1, 1]
    d.plot([0, 2.0], [0, 2.0], "--", color="#cbd5e1", lw=1.2, label="ideal")
    d.plot(col[:, 0], col[:, 1], "o", color=CU, ms=9, mec="white", mew=1.0, label="white atlas (biased)")
    d.plot(col[:, 0], col[:, 2], "s", color=CI, ms=9, mec="white", mew=1.0, label="joint (ν,τ_s) fit")
    en = np.abs(col[:, 1] - col[:, 0]).mean(); ej = np.abs(col[:, 2] - col[:, 0]).mean()
    d.set_xlabel(r"true $\nu$ (colored noise $\tau_s=0.6$)"); d.set_ylabel(r"recovered $\hat\nu$")
    d.set_title(fr"(D) Colored-noise axis: mean $|\nu$ err$|$ {en:.3f}$\to${ej:.3f}")
    d.legend(fontsize=9, framealpha=0.95, edgecolor="#cbd5e1")

    fig.tight_layout()
    out = os.path.abspath(os.path.join(HERE, "..", "figures", "harden_inverter.png"))
    fig.savefig(out, dpi=140); print("\nsaved", out)


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) >= 2 else ""
    if cmd == "coloratlas":
        build_colored_atlas()
    elif cmd == "interp":
        demo_interp()
    elif cmd == "classvote":
        demo_classvote()
    elif cmd == "colored":
        demo_colored()
    elif cmd == "harden":
        print("== (A) off-grid fitting =="); demo_interp()
        print("\n== (B) near-threshold class weighting =="); demo_classvote()
        print("\n== (C) colored-noise axis =="); demo_colored()
    elif cmd == "fig":
        harden_figure()
    else:
        print("usage: harden_inverter.py coloratlas | interp | classvote | colored | harden | fig")
