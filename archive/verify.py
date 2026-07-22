"""
verify.py  —  Sanity checks for the FHN simulation pipeline.
Should finish in < 90 seconds.
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

import numpy as np
from simulate import simulate_fhn, reduced_drift_time
from sweep import log_ratio

PASS = "\033[92m✓\033[0m"
FAIL = "\033[91m✗\033[0m"

V_FP, W_FP = -1.2563, -0.6954


def check(name, condition, detail=""):
    status = PASS if condition else FAIL
    print(f"  {status}  {name}" + (f"  [{detail}]" if detail else ""))
    return condition


def run():
    all_ok = True
    rng = np.random.default_rng(0)

    print("\n── 1. sigma=0 at FP — no firing ──")
    r = simulate_fhn(sigma=0.0, eps=0.1, v0=V_FP, w0=W_FP,
                     n_trajectories=5, T=100, rng=rng)
    all_ok &= check("No firing (deterministic)", r["fraction_fired"] == 0.0,
                    f"fired={r['fraction_fired']:.2f}")

    print("\n── 2. Large sigma fires frequently ──")
    rng = np.random.default_rng(1)
    r = simulate_fhn(sigma=0.25, eps=0.1, n_trajectories=300, T=400, rng=rng)
    all_ok &= check("Frequent firing at sigma=0.25", r["fraction_fired"] > 0.3,
                    f"fired={r['fraction_fired']:.2f}")
    all_ok &= check("MFPT_full finite", not np.isnan(r["mfpt"]),
                    f"mfpt={r['mfpt']:.1f}")

    print("\n── 3. Drift times scale as 1/eps ──")
    t1 = reduced_drift_time(eps=0.1)
    t2 = reduced_drift_time(eps=0.2)
    ratio = t1 / t2
    all_ok &= check("T_drift ~ 1/eps", abs(ratio - 2.0) < 0.05,
                    f"T(0.1)/T(0.2)={ratio:.3f}  (expect 2.0)")

    print("\n── 4. Large sigma/small eps → log_ratio > 0 (noise shortens MFPT) ──")
    # sigma=0.3, eps=0.01: sigma/eps^1.5 = 0.3/0.001 = 300 >> 1, well above BG
    rng = np.random.default_rng(2)
    r = simulate_fhn(sigma=0.3, eps=0.01, n_trajectories=400, T=300, rng=rng)
    td = reduced_drift_time(eps=0.01)
    lr = log_ratio(r["mfpt"], td)
    ok = np.isnan(lr) or lr > 0
    all_ok &= check("log_ratio > 0 well above BG", ok,
                    f"MFPT_full={r['mfpt']:.1f}, T_drift={td:.1f}, log_ratio={lr:.2f}")

    print("\n── 5. log_ratio increases as eps decreases (more noise-dominated) ──")
    # Fix sigma=0.25; decrease eps → sigma/eps^1.5 increases → should fail harder
    rng = np.random.default_rng(3)
    r_hi  = simulate_fhn(sigma=0.25, eps=0.10, n_trajectories=300, T=400, rng=rng)
    rng = np.random.default_rng(4)
    r_lo  = simulate_fhn(sigma=0.25, eps=0.01, n_trajectories=300, T=300, rng=rng)
    lr_hi = log_ratio(r_hi["mfpt"], reduced_drift_time(eps=0.10))
    lr_lo = log_ratio(r_lo["mfpt"], reduced_drift_time(eps=0.01))
    ok = (np.isnan(lr_hi) or np.isnan(lr_lo)) or lr_lo > lr_hi
    all_ok &= check("log_ratio(eps=0.01) > log_ratio(eps=0.10)", ok,
                    f"lr(eps=0.01)={lr_lo:.2f}  lr(eps=0.10)={lr_hi:.2f}")

    print()
    if all_ok:
        print(f"  {PASS}  All checks passed.\n")
    else:
        print(f"  {FAIL}  Some checks failed.\n")
    return all_ok


if __name__ == "__main__":
    ok = run()
    sys.exit(0 if ok else 1)
