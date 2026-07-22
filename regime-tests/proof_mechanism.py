"""
proof_mechanism.py — verifies the core mechanism of the shooting-characterisation proof
---------------------------------------------------------------------------------------
For a fixed noise realisation, with M = -d^2/dY^2 + Y - eta*xi on a grid:
  G(l) = lowest eigenvalue of M restricted to [l, inf) (Dirichlet at l).
Claim (Lemma 3 of ShootingCharacterisation_proof): G is monotone increasing in l, and
G(l)=0 exactly at Y_node = the largest zero of the energy-0 recessive solution of M
(marched down by the SAME-noise 3-term recurrence M u = 0). This is what licenses the
pathwise event identity  {Y_node <= t} = {G(t) >= 0}  with t DETERMINISTIC.
"""
from __future__ import annotations
import os
import numpy as np


def realize(seed, eta=1.0, Ymax=8.0, Ymin=-5.0, h=0.02):
    Y = np.arange(Ymax, Ymin - 1e-9, -h)          # descending grid (index 0 = top)
    n = len(Y)
    rng = np.random.default_rng(seed)
    zeta = rng.standard_normal(n) / np.sqrt(h)     # white-noise potential samples
    Vpot = Y - eta * zeta                          # potential of M at each node
    # energy-0 recessive solution marched DOWNWARD via M u = 0 (3-term recurrence)
    u = np.zeros(n); u[0] = 1e-8
    u[1] = u[0] * np.exp(np.sqrt(max(Y[0], 0.01)) * h)
    Ynode = np.nan
    for j in range(1, n - 1):
        u[j + 1] = 2 * u[j] - u[j - 1] + h * h * Vpot[j] * u[j]
        if u[j + 1] < 0:
            Ynode = Y[j + 1]; break
    # G(l) = smallest eigenvalue of the top-(k+1) sub-block (Y in [Y[k], Ymax])
    off = -1.0 / h**2; diag = 2.0 / h**2 + Vpot
    ls, Gs = [], []
    for k in range(30, n, 10):
        d = diag[:k + 1]
        M = np.diag(d) + np.diag(np.full(k, off), 1) + np.diag(np.full(k, off), -1)
        Gs.append(np.linalg.eigvalsh(M)[0]); ls.append(Y[k])
    ls, Gs = np.array(ls), np.array(Gs)
    s = np.where(np.diff(np.sign(Gs)))[0]
    lcross = np.interp(0, [Gs[s[0]], Gs[s[0] + 1]], [ls[s[0]], ls[s[0] + 1]]) if len(s) else np.nan
    return Ynode, lcross, ls, Gs


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    figdir = os.path.join(os.path.dirname(here), "figures")
    os.makedirs(figdir, exist_ok=True)

    print("\n=== Mechanism of Lemma 3: G(l) monotone, crosses 0 at Y_node ===")
    print(f"  {'seed':>4} {'Y_node (ODE recessive)':>22} {'l* : G(l*)=0':>14} {'G increasing?':>14}")
    data = {}
    for seed in (1, 2, 3, 4):
        Yn, lc, ls, Gs = realize(seed); data[seed] = (ls, Gs, Yn, lc)
        incr = bool(np.all(np.diff(Gs) < 0))   # Gs ordered by descending l => decreasing array = increasing in l
        print(f"  {seed:>4} {Yn:>22.4f} {lc:>14.4f} {str(incr):>14}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    cols = ["C0", "C1", "C2", "C3"]
    for seed, c in zip((1, 2, 3, 4), cols):
        ls, Gs, Yn, lc = data[seed]
        ax.plot(ls, Gs, "-", color=c, lw=1.8, label=fr"$G(\ell)$, seed {seed}")
        ax.plot([Yn], [0], "o", color=c, ms=9, mfc="white", mew=2)
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_xlabel(r"$\ell$"); ax.set_ylabel(r"$G(\ell)=$ ground state of $M|_{[\ell,\infty)}$")
    ax.set_title(r"$G$ monotone; circles $=Y_{\rm node}$ (ODE), where $G(\ell)=0$")
    ax.legend(fontsize=9, frameon=False); ax.grid(alpha=0.3)
    fig.tight_layout()
    out = os.path.join(figdir, "proof_mechanism.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print(f"\n  figure -> figures/proof_mechanism.png\n")


if __name__ == "__main__":
    main()
