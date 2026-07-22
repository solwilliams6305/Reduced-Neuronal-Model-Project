"""
tw_breakdown.py — How large can the noise be before the Tracy-Widom (single committed
first-passage) reduction of the noisy folded limit cycle breaks down?

Everything is gauged in the EFFECTIVE noise  eta = sigma / sqrt(eps)  (eps = ramp rate
through the fold = the role of eps_2).  Raw sigma is not universal; eta is.  We therefore
report eta* and CONFIRM it by collapsing two different ramp rates onto the same curve.

We compare, at matched eta:

  (TW reference)  canonical inner Riccati   dR = (R^2 - Y) dT + eta dB,  Y = Y0 - T,
                  whose first-blow-up level Y_node is TW_beta (beta = 4/eta^2) by Thm 1.
                  Has NO recrossing by construction (R -> +inf is one-way).

  (full system)   Bautin / fold-of-limit-cycles normal form (subcritical Hopf + SNPO):
                  dz = [ mu z + i*w*z + z(|z|^2 - |z|^4) ] dt + sigma dW,   (a = 1)
                  radial drift dr/dt = r(mu + r^2 - r^4); stable cycle, unstable cycle
                  (threshold) and stable origin coexist for mu in (mu_fold, 0),
                  mu_fold = -1/4 (fold of cycles), mu = 0 (subcritical Hopf).
                  mu(t) = mu0 - eps t is ramped slowly down through the fold.
                  The STABLE ORIGIN is the second attractor that lets the trajectory
                  shrink-and-pop (recross) -- exactly the effect that kills the TW law.

Breakdown diagnostics vs eta (all free of the unknown blow-down constants):
  * P_recross : fraction of passages that, after first dropping toward rest, climb back
                onto the cycle while still inside the bistable window.  TW assumes 0.
  * KS_std    : two-sample Kolmogorov-Smirnov distance between the STANDARDISED full-system
                peel-off and the STANDARDISED canonical Y_node at the same eta.  Standardising
                (subtract mean, divide by sd) removes the blow-down constants -> pure shape
                test of TW-ness.  Compared against the Monte-Carlo floor KS_floor.
  * div       : mean gap mu_first - mu_last between the first-passage and last-passage
                peel-off (the consistency gap of the previous discussion); ~0 when the
                escape is a single committed event, grows once recrossing sets in.
"""
from __future__ import annotations
import numpy as np

MU_FOLD = -0.25
OMEGA = 1.0


# ----------------------------- canonical Riccati (TW) -----------------------------
def canonical_Ynode(eta, N=1200, Y0=6.0, dT=4e-3, Resc=8.0, seed=0):
    rng = np.random.default_rng(seed)
    R = np.full(N, -np.sqrt(Y0))          # start on the attracting canard branch R ~ -sqrt(Y)
    Tnode = np.full(N, np.nan); done = np.zeros(N, bool)
    sdt = np.sqrt(dT); steps = int((Y0 + 6.0) / dT)
    for k in range(steps):
        if done.all():
            break
        Y = Y0 - k * dT
        R = R + (R * R - Y) * dT + eta * sdt * rng.standard_normal(N)
        R = np.minimum(R, 40.0)
        cr = (~done) & (R >= Resc); Tnode[cr] = k * dT; done |= cr
    Yn = Y0 - Tnode
    return Yn[np.isfinite(Yn)]


# --------------------- Bautin fold-of-cycles full slow passage --------------------
def bautin_passage(sigma, eps, N=400, mu0=0.20, mu_end=-0.50, dt=0.05, seed=0):
    rng = np.random.default_rng(seed)
    rs = np.sqrt((1 + np.sqrt(1 + 4 * mu0)) / 2)        # stable-cycle radius at mu0
    ph = rng.uniform(0, 2 * np.pi, N)
    x = rs * np.cos(ph); y = rs * np.sin(ph)
    r_hi, r_mid, r_lo = 0.85, 0.55, 0.30
    nsteps = int((mu0 - mu_end) / (eps * dt)); sdt = np.sqrt(dt)
    mu_first = np.full(N, np.nan); mu_last = np.full(N, mu0)
    everLow = np.zeros(N, bool); recross = np.zeros(N, bool); gotfirst = np.zeros(N, bool)
    for k in range(nsteps):
        mu = mu0 - eps * dt * k
        s = x * x + y * y; g = mu + s - s * s
        nx = rng.standard_normal(N); ny = rng.standard_normal(N)
        xn = x + (mu * x - OMEGA * y + x * g) * dt + sigma * sdt * nx
        yn = y + (mu * y + OMEGA * x + y * g) * dt + sigma * sdt * ny
        x = np.clip(xn, -5, 5); y = np.clip(yn, -5, 5)
        r = np.sqrt(x * x + y * y)
        nf = (~gotfirst) & (r < r_lo); mu_first[nf] = mu; gotfirst |= nf
        hi = r > r_hi; mu_last[hi] = mu
        recross |= everLow & hi & (mu > MU_FOLD + 1e-3)
        everLow |= r < r_mid
    return mu_first, mu_last, recross


def bautin_trace(sigma, eps, N=6, mu0=0.20, mu_end=-0.50, dt=0.05, sub=25, seed=0):
    """Record r(mu) for a few individual walkers, to visualise clean escape vs shrink-and-pop."""
    rng = np.random.default_rng(seed)
    rs = np.sqrt((1 + np.sqrt(1 + 4 * mu0)) / 2)
    ph = rng.uniform(0, 2 * np.pi, N); x = rs * np.cos(ph); y = rs * np.sin(ph)
    nsteps = int((mu0 - mu_end) / (eps * dt)); sdt = np.sqrt(dt)
    mus = []; R = []
    for k in range(nsteps):
        mu = mu0 - eps * dt * k
        s = x * x + y * y; g = mu + s - s * s
        xn = x + (mu * x - OMEGA * y + x * g) * dt + sigma * sdt * rng.standard_normal(N)
        yn = y + (mu * y + OMEGA * x + y * g) * dt + sigma * sdt * rng.standard_normal(N)
        x = np.clip(xn, -5, 5); y = np.clip(yn, -5, 5)
        if k % sub == 0:
            mus.append(mu); R.append(np.sqrt(x * x + y * y).copy())
    return np.array(mus), np.array(R)


# ---------------------------------- statistics -----------------------------------
def ks_std(a, b):
    a = (a - a.mean()) / a.std(); b = (b - b.mean()) / b.std()
    sa = np.sort(a); sb = np.sort(b); grid = np.concatenate([sa, sb])
    Fa = np.searchsorted(sa, grid, side="right") / len(sa)
    Fb = np.searchsorted(sb, grid, side="right") / len(sb)
    return float(np.max(np.abs(Fa - Fb)))


def skew(v):
    v = v - v.mean(); return float((v ** 3).mean() / (v.var() ** 1.5))


# ------------------------------------- main --------------------------------------
def run(etas, eps, canon, seed0=1, N_b=1500, keep=()):
    """canon: dict eta -> (Yn, Yn2). keep: etas at which to stash standardised samples."""
    rows = []; samples = {}
    for i, eta in enumerate(etas):
        sigma = eta * np.sqrt(eps)
        Yn, Yn2 = canon[eta]
        mf, ml, rc = bautin_passage(sigma, eps, N=N_b, seed=seed0 + 200 + i)
        ok = np.isfinite(mf); mfx = mf[ok]; mlx = ml[ok]
        ks_floor = ks_std(Yn, Yn2)
        ks_full = ks_std(mfx, Yn) if mfx.size > 20 else np.nan
        div = float(np.mean(mfx - mlx)) if mfx.size > 0 else np.nan
        rows.append(dict(eta=eta, sigma=sigma, P_recross=float(rc.mean()),
                         KS_full=ks_full, KS_floor=ks_floor, div=div, n_esc=int(mfx.size)))
        if eta in keep and mfx.size > 20:
            samples[eta] = ((mfx - mfx.mean()) / mfx.std(), (Yn - Yn.mean()) / Yn.std())
    return rows, samples


def show(tag, rows):
    print(f"\n=== {tag} ===")
    print(f"{'eta':>5} {'sigma':>7} {'P_recross':>10} {'KS_full':>8} {'KS_floor':>9} {'div':>7} {'nesc':>5}")
    for r in rows:
        print(f"{r['eta']:5.2f} {r['sigma']:7.4f} {r['P_recross']:10.3f} {r['KS_full']:8.3f} "
              f"{r['KS_floor']:9.3f} {r['div']:7.3f} {r['n_esc']:5d}")


def crossing(etas, vals, level):
    e = np.asarray(etas, float); v = np.asarray(vals, float)
    for j in range(1, len(e)):
        if v[j - 1] < level <= v[j]:
            return float(e[j - 1] + (level - v[j - 1]) * (e[j] - e[j - 1]) / (v[j] - v[j - 1]))
    return float("nan")


if __name__ == "__main__":
    etas = [0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4]
    canon = {eta: (canonical_Ynode(eta, N=2500, seed=1 + i),
                   canonical_Ynode(eta, N=2500, seed=501 + i)) for i, eta in enumerate(etas)}
    r1, samp = run(etas, 0.004, canon, seed0=1, N_b=1500, keep=(0.8, 2.0))
    show("eps = 0.004", r1)
    r2, _ = run(etas, 0.002, canon, seed0=7, N_b=1500)
    show("eps = 0.002 (collapse check)", r2)

    P1 = [r["P_recross"] for r in r1]; P2 = [r["P_recross"] for r in r2]
    eta_star_05 = np.nanmean([crossing(etas, P1, 0.05), crossing(etas, P2, 0.05)])
    eta_star_01 = np.nanmean([crossing(etas, P1, 0.01), crossing(etas, P2, 0.01)])
    print(f"\n  eta* (P_recross=1%)  ~ {eta_star_01:.2f}")
    print(f"  eta* (P_recross=5%)  ~ {eta_star_05:.2f}")

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].plot(etas, P1, "o-", color="#10b981", label=r"$\epsilon=0.004$")
    ax[0].plot(etas, P2, "s--", color="#4ea1ff", label=r"$\epsilon=0.002$")
    ax[0].axhline(0.05, color="#94a3b8", lw=1, ls=":")
    ax[0].axvline(eta_star_05, color="#f87272", lw=1.2, label=fr"$\eta^*\approx{eta_star_05:.1f}$")
    ax[0].set_xlabel(r"$\eta=\sigma/\sqrt{\epsilon}$"); ax[0].set_ylabel(r"$P_{\rm recross}$")
    ax[0].set_title("Recrossing onset (collapses in $\\eta$)"); ax[0].legend(frameon=False)

    KSf = [r["KS_full"] for r in r1]; KSfl = [r["KS_floor"] for r in r1]
    ax[1].plot(etas, KSf, "o-", color="#10b981", label="full vs TW")
    ax[1].plot(etas, KSfl, "-", color="#94a3b8", label="MC floor")
    ax[1].fill_between(etas, 0, KSfl, color="#94a3b8", alpha=0.25)
    ax[1].axvline(eta_star_05, color="#f87272", lw=1.2)
    ax[1].set_xlabel(r"$\eta$"); ax[1].set_ylabel("standardised KS distance")
    ax[1].set_title("Peel-off shape departs from TW"); ax[1].legend(frameon=False)

    for eta, col, sd in [(0.7, "#10b981", 11), (2.0, "#f87272", 12)]:
        mus, R = bautin_trace(eta * np.sqrt(0.004), 0.004, N=6, seed=sd)
        for j in range(R.shape[1]):
            ax[2].plot(mus, R[:, j], color=col, lw=0.9, alpha=0.7)
    ax[2].axvline(MU_FOLD, color="#94a3b8", lw=1.2, ls="--")
    ax[2].text(MU_FOLD, 1.28, "fold", color="#64748b", fontsize=8, ha="center")
    ax[2].axvline(0.0, color="#94a3b8", lw=1, ls=":")
    ax[2].text(0.0, 1.28, "Hopf", color="#64748b", fontsize=8, ha="center")
    ax[2].set_xlabel(r"$\mu$ (ramped down $\rightarrow$)"); ax[2].set_ylabel(r"radius $r$")
    ax[2].invert_xaxis()
    ax[2].set_title(r"green $\eta$=0.7: clean escape / red $\eta$=2: shrink-and-pop")
    fig.tight_layout()
    out = "/sessions/trusting-fervent-archimedes/mnt/Reduced Neuronal Model Project/figures/tw_breakdown.png"
    fig.savefig(out, dpi=130); print("saved", out)
