"""
verify_kernel.py — sanity checks for kernel.py
Should finish in < 2 minutes.
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
from kernel import FHN2D, simulate_kernel, reduced_prediction, compute_cv, log_ratio

PASS = "\033[92m✓\033[0m"
FAIL = "\033[91m✗\033[0m"

def check(name, condition, detail=""):
    s = PASS if condition else FAIL
    print(f"  {s}  {name}" + (f"  [{detail}]" if detail else ""))
    return condition

all_ok = True
rng = np.random.default_rng(0)

# 1. Model geometry
print("\n── 1. FHN2D geometry ──")
m = FHN2D(I=-0.1, a=0.7, b=0.8)
print(f"  {m}")
all_ok &= check("V_FP < -1 (left branch)", m.V_FP < -1.0, f"V_FP={m.V_FP:.4f}")
all_ok &= check("ΔU > 0", m.delta_U > 0, f"ΔU={m.delta_U:.5f}")
all_ok &= check("regime=excitable", m.regime == "excitable", m.regime)
all_ok &= check("w_fold_left = I + 2/3", abs(m.w_fold_left - (m.I + 2/3)) < 1e-10)

# 2. No firing at sigma=0
print("\n── 2. sigma=0 — no firing ──")
res = simulate_kernel(m, sigma=0.0, eps=0.1, mode="first_passage",
                      n_trajectories=5, T=100, rng=np.random.default_rng(1))
all_ok &= check("No firing (deterministic)", res["fraction_fired"] == 0.0,
                f"fired={res['fraction_fired']:.2f}")

# 3. Large sigma fires frequently
print("\n── 3. Large sigma fires ──")
res = simulate_kernel(m, sigma=0.25, eps=0.1, mode="first_passage",
                      n_trajectories=300, T=400, rng=np.random.default_rng(2))
all_ok &= check("Frequent firing", res["fraction_fired"] > 0.3,
                f"fired={res['fraction_fired']:.2f}")
all_ok &= check("MFPT finite", not np.isnan(res["mfpt"]), f"mfpt={res['mfpt']:.1f}")

# 4. T_drift scales as 1/eps
print("\n── 4. T_drift ~ 1/eps ──")
t1 = reduced_prediction(m, eps=0.1, mode="first_passage")
t2 = reduced_prediction(m, eps=0.2, mode="first_passage")
ratio = t1 / t2
all_ok &= check("T_drift ~ 1/eps", abs(ratio - 2.0) < 0.05,
                f"T(0.1)/T(0.2)={ratio:.3f}  (expect 2.0)")

# 5. log_ratio > 0 above BG curve
print("\n── 5. log_ratio > 0 well above BG ──")
res = simulate_kernel(m, sigma=0.3, eps=0.01, mode="first_passage",
                      n_trajectories=300, T=300, rng=np.random.default_rng(3))
td  = reduced_prediction(m, eps=0.01)
lr  = log_ratio(res["mfpt"], td)
all_ok &= check("log_ratio > 0", np.isnan(lr) or lr > 0,
                f"MFPT={res['mfpt']:.1f}, T_drift={td:.1f}, lr={lr:.2f}")

# 6. ISI mode collects sequences
print("\n── 6. ISI sequence mode ──")
m_tonic = FHN2D(I=0.5, a=0.7, b=0.8)   # well past SNIC
print(f"  {m_tonic}")
res = simulate_kernel(m_tonic, sigma=0.05, eps=0.1, mode="isi_sequence",
                      n_trajectories=50, n_isi=5, T=2000,
                      rng=np.random.default_rng(4))
all_ok &= check("ISI mean finite", not np.isnan(res["isi_mean"]),
                f"isi_mean={res['isi_mean']:.2f}")
all_ok &= check("CV finite", not np.isnan(res["cv"]),
                f"cv={res['cv']:.3f}")
all_ok &= check("CV < 1 (regular tonic)", np.isnan(res["cv"]) or res["cv"] < 1.0,
                f"cv={res['cv']:.3f}")

# 7. T_cycle for tonic mode
print("\n── 7. T_cycle (tonic reduced prediction) ──")
tc = reduced_prediction(m_tonic, eps=0.1, mode="isi_sequence")
all_ok &= check("T_cycle finite", not np.isnan(tc), f"T_cycle={tc:.2f}")

# 8. compute_cv
print("\n── 8. compute_cv ──")
isis = np.array([1.0, 1.1, 0.9, 1.05, 0.95])
cv   = compute_cv(isis)
all_ok &= check("CV ≈ 0.07 for near-regular ISIs", abs(cv - np.std(isis)/np.mean(isis)) < 1e-10,
                f"cv={cv:.4f}")

print()
if all_ok:
    print(f"  {PASS}  All checks passed.\n")
else:
    print(f"  {FAIL}  Some checks failed.\n")
sys.exit(0 if all_ok else 1)
