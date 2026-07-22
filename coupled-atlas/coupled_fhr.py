"""
coupled_fhr.py — Two electrically coupled FitzHugh–Rinzel (FHR) units.
=====================================================================

Single unit (same convention as regime-tests/mmo_fhr_alpha_fine.py):

    v' = v - v^3/3 - w + y + I        (fast)
    w' = eps (v + A - B w)            (slow recovery)
    y' = eps*delta (c - v)           (super-slow)   -> three timescales

Defaults A=0.7, B=0.8, eps=0.08, delta=0.2, I=0.30.  The folded-node MMO
("staircase") band is c in [-0.95, -0.74]; deeper c (toward -1) => more
small-amplitude oscillations (SAOs) per large spike.

Electrical (diffusive) coupling enters the FAST variable only — this is the
physically correct gap-junction term, well-defined here because FHR has a
smooth spike (unlike Izhikevich, whose reset makes I_gap during the spike
ill-defined; see COUPLED_REDUCED_NEURON_LIT_MAP.md §II):

    v_i' += g * Phi( v_j(t - tau) - v_i )

Coupling forms Phi (the "coupling form" novelty axis):
    'linear'      Phi(x) = x                 standard gap junction / diffusive
    'tanh'        Phi(x) = tanh(k x)/k        saturating (large gaps clipped)
    'rectifying'  Phi(x) = max(x, 0)          one-way / rectifying junction

Novelty axes exposed:
    g       coupling strength
    form    coupling form (above)
    tau     coupling/transmission delay (time units)
    sigma   additive white noise on v  (Euler–Maruyama;  matches kernel.py)
    + per-unit heterogeneity (each unit has its own A,B,eps,delta,I,c)

Integration: RK4 for the deterministic case (sigma=0, matches the FHR script),
Euler–Maruyama when sigma>0.  Delay handled by a circular history buffer with
the delayed partner voltage frozen across RK4 substeps (method of steps; valid
for the small delays this project targets).

This module is intentionally dependency-light (numpy only), mirroring kernel.py.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field, replace


# ---------------------------------------------------------------------------
# Single-unit parameters
# ---------------------------------------------------------------------------
@dataclass
class FHR:
    """FitzHugh–Rinzel single-unit parameters (defaults = the repo's working tuning)."""
    A: float = 0.7
    B: float = 0.8
    eps: float = 0.08
    delta: float = 0.2
    I: float = 0.30
    c: float = -0.80          # control parameter; MMO band c in [-0.95, -0.74]


# ---------------------------------------------------------------------------
# Coupling forms
# ---------------------------------------------------------------------------
def _phi_linear(x):     return x
def _phi_tanh(x, k=2.0): return np.tanh(k * x) / k
def _phi_rect(x):       return np.where(x > 0.0, x, 0.0)

_FORMS = {"linear": _phi_linear, "tanh": _phi_tanh, "rectifying": _phi_rect}


# ---------------------------------------------------------------------------
# Coupled-pair configuration
# ---------------------------------------------------------------------------
@dataclass
class CoupledFHR:
    """Two electrically coupled FHR units (gap-junction / diffusive coupling)."""
    unit1: FHR = field(default_factory=FHR)
    unit2: FHR = field(default_factory=FHR)
    g: float = 0.0               # coupling strength
    form: str = "linear"         # 'linear' | 'tanh' | 'rectifying'
    tau: float = 0.0             # coupling delay (time units)
    sigma: float = 0.0           # additive noise amplitude on v

    # ---- helpers ----
    def _arrays(self):
        u = (self.unit1, self.unit2)
        return (np.array([x.A for x in u]),
                np.array([x.B for x in u]),
                np.array([x.eps for x in u]),
                np.array([x.delta for x in u]),
                np.array([x.I for x in u]),
                np.array([x.c for x in u]))

    @staticmethod
    def identical(unit: FHR | None = None, **kw) -> "CoupledFHR":
        """Two identical units with the given single-unit parameters."""
        u = unit if unit is not None else FHR()
        return CoupledFHR(unit1=replace(u), unit2=replace(u), **kw)

    @staticmethod
    def mismatched(base: FHR | None = None, dc: float = 0.0, dI: float = 0.0,
                   **kw) -> "CoupledFHR":
        """Two units split symmetrically in c (by dc) and/or I (by dI)."""
        b = base if base is not None else FHR()
        u1 = replace(b, c=b.c - dc / 2.0, I=b.I - dI / 2.0)
        u2 = replace(b, c=b.c + dc / 2.0, I=b.I + dI / 2.0)
        return CoupledFHR(unit1=u1, unit2=u2, **kw)


# ---------------------------------------------------------------------------
# Integrator
# ---------------------------------------------------------------------------
def integrate(model: CoupledFHR,
              T: float = 7000.0,
              dt: float = 0.04,
              x0=None,
              t_warmup: float = 500.0,
              stride: int = 1,
              rng: np.random.Generator | None = None,
              seed: int | None = None) -> dict:
    """
    Integrate the coupled pair.

    Returns dict with:
        t : (nrec,)        time
        v : (nrec, 2)      fast variables  [unit1, unit2]
        w : (nrec, 2)      recovery
        y : (nrec, 2)      super-slow
        dt, model          echoed back
    """
    A, B, eps, delta, I, c = model._arrays()
    phi = _FORMS[model.form]
    g, sigma, tau = model.g, model.sigma, model.tau
    stochastic = sigma > 0.0

    if rng is None:
        rng = np.random.default_rng(seed)

    # initial state: split the two units in phase so coupling has something to do
    if x0 is None:
        V = np.array([0.5, -0.5]); W = np.zeros(2); Y = np.zeros(2)
    else:
        V, W, Y = (np.array(z, dtype=float) for z in x0)

    lag = int(round(tau / dt))
    buflen = lag + 1
    Vbuf = np.tile(V.copy(), (buflen, 1))      # circular buffer of past v
    head = 0                                    # Vbuf[head] is most-recent v

    sqrt_dt = np.sqrt(dt)

    def drift(V, W, Y, Vd):
        # Vd = delayed partner-voltage source (length-2 array of v at t-tau)
        v_partner = Vd[::-1]                    # partner of 0 is 1, of 1 is 0
        coup = g * phi(v_partner - V)
        dV = V - V**3 / 3.0 - W + Y + I + coup
        dW = eps * (V + A - B * W)
        dY = eps * delta * (c - V)
        return dV, dW, dY

    def delayed_v():
        return Vbuf[(head - lag) % buflen]

    def advance(V, W, Y):
        nonlocal head
        Vd = delayed_v()
        if stochastic:
            dV, dW, dY = drift(V, W, Y, Vd)
            V = V + dV * dt + sigma * sqrt_dt * rng.standard_normal(2)
            W = W + dW * dt
            Y = Y + dY * dt
        else:
            k1 = drift(V, W, Y, Vd)
            k2 = drift(V + .5*dt*k1[0], W + .5*dt*k1[1], Y + .5*dt*k1[2], Vd)
            k3 = drift(V + .5*dt*k2[0], W + .5*dt*k2[1], Y + .5*dt*k2[2], Vd)
            k4 = drift(V + dt*k3[0],    W + dt*k3[1],    Y + dt*k3[2],    Vd)
            V = V + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.0
            W = W + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.0
            Y = Y + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6.0
        head = (head + 1) % buflen
        Vbuf[head] = V
        return V, W, Y

    # warm-up (discard transient)
    for _ in range(int(round(t_warmup / dt))):
        V, W, Y = advance(V, W, Y)

    # record
    n = int(round(T / dt))
    nrec = n // stride
    t = np.empty(nrec); vrec = np.empty((nrec, 2))
    wrec = np.empty((nrec, 2)); yrec = np.empty((nrec, 2))
    ti = 0.0
    ri = 0
    for k in range(n):
        V, W, Y = advance(V, W, Y)
        ti += dt
        if (k % stride) == 0 and ri < nrec:
            t[ri] = ti; vrec[ri] = V; wrec[ri] = W; yrec[ri] = Y
            ri += 1

    return {"t": t[:ri], "v": vrec[:ri], "w": wrec[:ri], "y": yrec[:ri],
            "dt": dt * stride, "model": model}


# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------
def peaks(v: np.ndarray) -> np.ndarray:
    """Indices of local maxima (same rule as the repo's mmo_fhr peak detector)."""
    dv = np.diff(v)
    return np.where((dv[:-1] > 0) & (dv[1:] <= 0))[0] + 1


def mmo_counts(v: np.ndarray, thr: float = 0.0) -> dict:
    """Count large (peak>thr) vs small (peak<=thr) oscillations -> MMO signature."""
    pk = v[peaks(v)]
    nL = int(np.sum(pk > thr))
    nS = int(np.sum(pk <= thr))
    return {"nL": nL, "nS": nS, "ratio_S_per_L": (nS / nL) if nL else np.nan,
            "n_peaks": int(pk.size)}


def sync_chi(v: np.ndarray) -> float:
    """
    Golomb–Hansel synchrony measure chi in [0,1] for two columns of v(t).
    chi^2 = Var_t( mean_i v_i ) / mean_i Var_t( v_i ).
    chi -> 1 fully synchronous, chi -> 0 incoherent.
    """
    vbar = v.mean(axis=1)
    num = np.var(vbar)
    den = np.mean(np.var(v, axis=0))
    return float(np.sqrt(num / den)) if den > 0 else np.nan


def large_spike_times(t: np.ndarray, v: np.ndarray, thr: float = 0.0) -> np.ndarray:
    """Times of large-amplitude spikes (peak voltage > thr) for one unit."""
    idx = peaks(v)
    idx = idx[v[idx] > thr]
    return t[idx]


def phase_difference(t: np.ndarray, v: np.ndarray, thr: float = 0.0) -> float:
    """
    Mean relative phase between the two units' large spikes, in [0, 0.5]
    (0 = in-phase, 0.5 = anti-phase), using nearest-spike time differences
    normalised by unit-1's mean period.  NaN if either train is too short.
    """
    t1 = large_spike_times(t, v[:, 0], thr)
    t2 = large_spike_times(t, v[:, 1], thr)
    if t1.size < 3 or t2.size < 2:
        return np.nan
    period = np.median(np.diff(t1))
    if not np.isfinite(period) or period <= 0:
        return np.nan
    # nearest unit-2 spike to each interior unit-1 spike
    d = []
    for tt in t1[1:-1]:
        j = np.argmin(np.abs(t2 - tt))
        dphi = ((t2[j] - tt) / period) % 1.0
        if dphi > 0.5:
            dphi = 1.0 - dphi          # fold to [0, 0.5]
        d.append(dphi)
    return float(np.mean(d)) if d else np.nan


# ---------------------------------------------------------------------------
# quick self-test
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # single-unit MMO sanity check at g=0
    m = CoupledFHR.identical(FHR(c=-0.85), g=0.0)
    out = integrate(m, T=4000.0, t_warmup=1500.0, dt=0.04)
    print("unit1 MMO counts:", mmo_counts(out["v"][:, 0]))
    print("unit2 MMO counts:", mmo_counts(out["v"][:, 1]))
    print("chi (g=0):", round(sync_chi(out["v"]), 3))
