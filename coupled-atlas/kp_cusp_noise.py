"""
kp_cusp_noise.py — noisy peel-off law of the Kristiansen–Pedersen coupled-FHN cusp.
===================================================================================

At the located cusp/MMO regime (c≈0.99, g≈-0.12, repulsive), add additive noise and
measure the peel-off level distribution — the value of the slow variable w_1 at the
escape from the upper slow manifold (the fold at v=1). Compare the coupled-cusp pair
(g=-0.12) to the single-unit baseline (g=0, a plain relaxation fold).

Question: does the cusp escape law differ in SHAPE from the single-fold escape?
Per Kristiansen–Pedersen the cusp reduces to the WEBER (parabolic-cylinder) equation —
a folded-node-type 2nd-order rung — NOT the literal 3rd-order Pearcey. So we expect the
cusp escape to be a non-trivial, structured (canard-funnel) law, distinct from the bare
relaxation fold.

Honest caveat: the full-model peel-off carries OUTER (Gaussian) corrections on top of the
inner edge law, so this is a SHAPE comparison, not a clean inner-law identification; the
definitive rung assignment is K-P's Weber reduction.

Output: figures/kp_cusp_noise.png + printed summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)


def run(g, c=0.99, eps=0.015, sigma=0.02, M=700, T=2600.0, dt=0.02, warmup=700.0, seed=0):
    rng = np.random.default_rng(seed); sdt = np.sqrt(dt)
    v = np.tile([0.2, -0.2], (M, 1)).astype(float); w = np.zeros((M, 2))

    def step(v, w, noise):
        coup = g * (v[:, ::-1] - v)
        vn = v + (-v**3 + 3*v - w + coup) * dt + noise
        wn = w + eps * (v - c) * dt
        return vn, wn

    for _ in range(int(warmup/dt)):
        v, w = step(v, w, 0.0)

    on_up = np.zeros(M, bool); wmax = np.full(M, -9.0); peels = []
    n = int(T/dt)
    for _ in range(n):
        noise = sigma * sdt * rng.standard_normal((M, 2))
        v, w = step(v, w, noise)
        v1, w1 = v[:, 0], w[:, 0]
        on_up |= (v1 > 1.2)
        wmax = np.where(on_up, np.maximum(wmax, w1), wmax)
        esc = on_up & (v1 < 0.5)
        if esc.any():
            peels.append(wmax[esc].copy())
        on_up = np.where(esc, False, on_up)
        wmax = np.where(esc, -9.0, wmax)
    return np.concatenate(peels) if peels else np.array([])


def moments(x):
    m = x.mean(); d = x - m; var = np.mean(d**2)
    return m, np.sqrt(var), np.mean(d**3)/var**1.5, np.mean(d**4)/var**2 - 3.0, x.size


def main():
    t0 = time.time()
    print("=" * 70)
    print("KP coupled-FHN cusp — noisy peel-off law (c=0.99, ε=0.015, σ=0.02)")
    print("=" * 70)
    res = {}
    for label, g in [("single fold (g=0)", 0.0), ("coupled cusp (g=-0.12)", -0.12)]:
        x = run(g, seed=7)
        res[label] = x
        if x.size > 20:
            m, s, sk, ku, n = moments(x)
            print(f"  {label:>26}: mean {m:.3f}  std {s:.4f}  skew {sk:+.3f}  kurt {ku:+.3f}  (n={n})")
        else:
            print(f"  {label:>26}: too few escapes ({x.size})")

    a = res["single fold (g=0)"]; b = res["coupled cusp (g=-0.12)"]
    if a.size > 20 and b.size > 20:
        se = np.sqrt(6.0 / min(a.size, b.size))
        dsk = moments(b)[2] - moments(a)[2]
        print(f"\n  skew(cusp) − skew(fold) = {dsk:+.3f}  (s.e.≈{se:.3f}):  "
              f"{'DIFFERENT shape' if abs(dsk) > 2*se else 'not resolved'}")
    print("\n  Rung assignment (from Kristiansen–Pedersen): the slow-fast cusp reduces to the")
    print("  WEBER equation — a folded-node-type 2nd-order rung — NOT the 3rd-order Pearcey.")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for label, col in [("single fold (g=0)", "#888780"), ("coupled cusp (g=-0.12)", "#7a3b8f")]:
        x = res[label]
        if x.size > 20:
            z = (x - x.mean()) / x.std()
            ax[0].hist(z, bins=55, range=(-4, 4), density=True, histtype="step", lw=1.8,
                       color=col, label=f"{label} (skew {moments(x)[2]:+.2f})")
    gg = np.linspace(-4, 4, 200)
    ax[0].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[0].set_xlabel("standardized peel-off w₁"); ax[0].set_ylabel("density")
    ax[0].set_title("(A) peel-off law: single fold vs coupled cusp")
    ax[0].legend(fontsize=8.5, frameon=False)

    # ladder context (skew comparison)
    labels = ["fold/TW\n(n=2)", "cusp/Pearcey\n(n=3, ODE ladder)", "Weber/k≈2\n(multicritical)",
              "KP cusp\n(this, full FHN)"]
    vals = [0.22, 0.10, 0.60, moments(b)[2] if b.size > 20 else np.nan]
    cols = ["#b3402b", "#1f3b73", "#2c7d59", "#7a3b8f"]
    ax[1].bar(range(4), vals, color=cols)
    ax[1].set_xticks(range(4)); ax[1].set_xticklabels(labels, fontsize=8)
    ax[1].axhline(0.22, color="#b3402b", lw=0.8, ls="--")
    ax[1].set_ylabel("peel-off skewness")
    ax[1].set_title("(B) where the KP cusp lands among the ladders")

    fig.suptitle("Kristiansen–Pedersen coupled-FHN cusp — noisy peel-off (Weber rung, per K-P)",
                 fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "kp_cusp_noise.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
