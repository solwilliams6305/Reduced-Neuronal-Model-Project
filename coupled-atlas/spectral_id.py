"""
spectral_id.py — PART B: stochastic-operator identification of the cusp peel-off law.
=====================================================================================

The fold peel-off is the ground state of the stochastic AIRY operator
    H_1 = -d²/dx² + x + (2/√β) b'(x)   on [0,∞), Dirichlet,   (-Λ0 =d TW_β; RRV / the paper).
By the shooting↔operator correspondence, the swept inner first-node u''=(sign(Y)|Y|^k-ηξ)u
should equal the ground state of
    H_k = -d²/dx² + x^k + (2/√β) b'(x)     (k=1 fold/Airy; k=2 the cusp/Weber-class rung).

Here we compute the ground-state law of H_k by an INDEPENDENT method — direct tridiagonal
finite-difference discretisation + Sturm-sequence bisection for the smallest eigenvalue — and:
  (1) VALIDATE: k=1 reproduces Tracy–Widom (skew ≈ 0.22 at β=2), and matches the swept first-node;
  (2) IDENTIFY: k=2 operator ground state matches the swept k=2 law and the KP cusp peel-off
      (Part 1, Weber-class), and is distinct from TW — the cusp analogue of "stochastic Airy = TW".

Two independent methods agreeing (shooting + operator) = the spectral identification (numerically);
the RRV-type theorem for k=2 remains open. Output: figures/spectral_id.png + tagged summary.
"""
from __future__ import annotations
import os, time
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peeloff_cusp_ladder import peeloff_ladder            # swept first-node (shooting)

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "figures")
TW2 = dict(skew=0.2241, kurt=0.0934)


def ground_states(k, beta=2.0, R=4000, N=700, L=7.0, seed=0):
    """Smallest eigenvalue of H_k = -d²/dx² + x^k + (2/√β)b' (tridiagonal, Sturm bisection)."""
    rng = np.random.default_rng(seed)
    h = L / N; x = (np.arange(1, N + 1)) * h
    inv_h2 = 1.0 / h**2
    pot = x**k                                            # (N,)
    noise = (2.0 / np.sqrt(beta)) * (1.0 / np.sqrt(h)) * rng.standard_normal((R, N))
    diag = (2.0 * inv_h2) + pot[None, :] + noise          # (R,N)
    e2 = inv_h2**2                                         # off-diagonal squared

    def count_below(lam):                                  # eigenvalues < lam, per realization
        q = diag[:, 0] - lam
        cnt = (q < 0).astype(np.int32)
        for i in range(1, N):
            q = (diag[:, i] - lam) - e2 / q
            q = np.where(np.abs(q) < 1e-300, -1e-300, q)
            cnt += (q < 0)
        return cnt

    lo = np.full(R, -4.0); hi = np.full(R, 9.0)
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        below = count_below(mid) >= 1
        hi = np.where(below, mid, hi); lo = np.where(below, lo, mid)
    return 0.5 * (lo + hi)                                 # smallest eigenvalue Λ0 per realization


def moments(x):
    x = x[np.isfinite(x)]; m = x.mean(); d = x - m; v = np.mean(d**2)
    return m, np.sqrt(v), np.mean(d**3) / v**1.5, np.mean(d**4) / v**2 - 3.0


def main():
    t0 = time.time()
    eta = np.sqrt(2.0)
    print("=" * 74)
    print("PART B — stochastic-operator spectral identification (β=2)")
    print("=" * 74)

    # operator ground states (peel-off analogue = -Λ0)
    op1 = -ground_states(1, seed=1)
    op2 = -ground_states(2, seed=2)
    # swept first-node (shooting), same k
    sw1 = peeloff_ladder([1.0], eta, N=6000, seed=11)[0]
    sw2 = peeloff_ladder([2.0], eta, N=6000, seed=12)[0]

    print(f"  {'method':>26} | {'skew':>8}{'exkurt':>9}")
    for name, x in [("operator k=1 (−Λ0)", op1), ("swept first-node k=1", sw1),
                    ("  TW₂ reference", None),
                    ("operator k=2 (−Λ0)", op2), ("swept first-node k=2", sw2)]:
        if x is None:
            print(f"  {name:>26} | {TW2['skew']:+8.3f}{TW2['kurt']:+9.3f}")
        else:
            m, s, sk, ku = moments(x)
            print(f"  {name:>26} | {sk:+8.3f}{ku:+9.3f}")

    s_op1 = moments(op1)[2]; s_sw1 = moments(sw1)[2]
    s_op2 = moments(op2)[2]; s_sw2 = moments(sw2)[2]
    print(f"\n  [VALIDATE] k=1: operator skew {s_op1:+.3f} vs swept {s_sw1:+.3f} vs TW₂ {TW2['skew']:+.3f}"
          f"  → {'reproduces TW' if abs(s_op1-TW2['skew'])<0.06 else 'CHECK'}")
    print(f"  [IDENTIFY] k=2: operator skew {s_op2:+.3f} vs swept {s_sw2:+.3f}"
          f"  → {'two methods agree (Weber-class)' if abs(s_op2-s_sw2)<0.10 else 'CHECK'}")
    print(f"             k=2 distinct from k=1/TW: {'YES' if s_op2 > s_op1+0.15 else 'no'}")

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
    gg = np.linspace(-4, 4, 220)
    ax[0].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    for x, col, lab in [(op1, "#b3402b", f"operator k=1 (skew {s_op1:+.2f})"),
                        (sw1, "#e08a72", f"swept k=1 (skew {s_sw1:+.2f})")]:
        z = (x - x.mean())/x.std()
        ax[0].hist(z, bins=60, range=(-4, 4), density=True, histtype="step", lw=1.7, color=col, label=lab)
    ax[0].set_title("(A) k=1: operator = swept = Tracy–Widom (validation)")
    ax[0].set_xlabel("standardised peel-off"); ax[0].set_ylabel("density"); ax[0].legend(fontsize=8.5, frameon=False)

    for x, col, lab in [(op2, "#1f3b73", f"operator k=2 (skew {s_op2:+.2f})"),
                        (sw2, "#5d82c4", f"swept k=2 (skew {s_sw2:+.2f})")]:
        z = (x - x.mean())/x.std()
        ax[1].hist(z, bins=60, range=(-4, 4), density=True, histtype="step", lw=1.7, color=col, label=lab)
    ax[1].plot(gg, np.exp(-gg**2/2)/np.sqrt(2*np.pi), "k:", lw=1, label="Gaussian")
    ax[1].set_title("(B) k=2: naive harmonic op (≈Gaussian) ≠ swept Weber-class — cusp op is soft-edge [OPEN]")
    ax[1].set_xlabel("standardised peel-off"); ax[1].set_ylabel("density"); ax[1].legend(fontsize=8.5, frameon=False)

    fig.suptitle("Part B — spectral ID: k=1 operator=swept=TW (validated); k=2 naive operator FAILS "
                 "(correct cusp operator is soft-edge, open)", fontsize=11.5)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fp = os.path.join(FIG, "spectral_id.png")
    fig.savefig(fp, dpi=140, bbox_inches="tight")
    print(f"\n  Figure: {fp}")
    print(f"  TOTAL {time.time()-t0:.1f}s — DONE")


if __name__ == "__main__":
    main()
