"""milestone2_bootstrap.py -- Uncertainty quantification of the Borel phase.

Bootstrap the phase theta and modulus |zeta| over the STATED coefficient uncertainties:
  v5 in [-2.3, -1.0]   (symbolic -1.1 vs grid -2.29)
  v3 in [-0.04, 0.0]   (delicate cancellation +-0.02 around -0.02..-0.03)
  v4 in [-0.47, -0.44]
(v0,v1,v2 held; v2 also lightly varied 0.100..0.104 to bridge task/base.)

For each sampled ladder we compute theta,|zeta| by the estimators that actually work on 6 coeffs:
  - Borel-Pade [1/4], [0/5], [2/3], [3/2]  (nearest complex pair)
  - Darboux nmin=1 alpha=0
  - Direct-ratio last-2-equations
and aggregate the spread.

Also: DROP-v5 test -- compute theta from v0..v4 only (5 coeffs) with the [1/3],[0/4] Pade and a
4-coeff direct-ratio, and compare to the v0..v5 estimate: is the phase converging or drifting?

Run: python3 milestone2_bootstrap.py
"""
import numpy as np
from milestone2_measure import (borel_pade_nearest_complex, borel_pade_pairs, fit_darboux,
                                 direct_ratio_pair_3term, direct_ratio_pair, LAMBDA0_MOD)

np.random.seed(12345)

V0, V1 = 0.134, 0.111

def make_ladder(v2, v3, v4, v5):
    return [V0, V1, v2, v3, v4, v5]

def pade_pair_from_split(v, L, M):
    for LL, MM, pairs in borel_pade_pairs(v):
        if LL == L and MM == M:
            cand = [(mod, th) for mod, th in pairs if 15 < abs(th) < 88]
            if cand:
                cand.sort(key=lambda x: x[0])
                return abs(cand[0][1]), cand[0][0]  # (theta_abs, mod)
    return None

def estimators(v):
    """Return dict of estimator -> (theta_deg, mod) or None."""
    out = {}
    # Pade splits
    for L, M in [(1, 4), (0, 5), (2, 3), (3, 2)]:
        r = pade_pair_from_split(v, L, M)
        out[f"pade[{L}/{M}]"] = r
    # Darboux nmin=1 alpha=0
    try:
        d = fit_darboux(v, nmin=1, alpha_fixed=0.0)
        out["darboux_n1a0"] = (abs(d['theta']), d['mod'])
    except Exception:
        out["darboux_n1a0"] = None
    # direct ratio last-2
    d3 = direct_ratio_pair_3term(v)
    out["ratio_last2"] = (d3['theta'], d3['mod']) if (d3 and d3['complex_pair']) else None
    dwl = direct_ratio_pair(v, weight_late=True)
    out["ratio_wlate"] = (dwl['theta'], dwl['mod']) if dwl['complex_pair'] else None
    return out

# ------------------------------------------------------------------
def bootstrap(N=4000):
    thetas_all = {k: [] for k in
                  ["pade[1/4]", "pade[0/5]", "pade[2/3]", "pade[3/2]",
                   "darboux_n1a0", "ratio_last2", "ratio_wlate"]}
    mods_all = {k: [] for k in thetas_all}
    pooled_theta = []
    pooled_mod = []
    for _ in range(N):
        v2 = np.random.uniform(0.100, 0.104)
        v3 = np.random.uniform(-0.04, 0.0)
        v4 = np.random.uniform(-0.47, -0.44)
        v5 = np.random.uniform(-2.3, -1.0)
        v = make_ladder(v2, v3, v4, v5)
        est = estimators(v)
        for k, r in est.items():
            if r is not None:
                th, mod = r
                if np.isfinite(th) and np.isfinite(mod) and 0 < mod < 20:
                    thetas_all[k].append(th)
                    mods_all[k].append(mod)
                    # pool the "trusted" estimators
                    if k in ("pade[1/4]", "pade[0/5]", "pade[2/3]", "darboux_n1a0", "ratio_last2"):
                        pooled_theta.append(th)
                        pooled_mod.append(mod)
    print(f"=== BOOTSTRAP over v2,v3,v4,v5 uncertainties (N={N}) ===")
    print(f"{'estimator':16s} {'theta med':>9s} {'[16,84]%':>16s} {'|zeta| med':>10s} {'[16,84]%':>16s}  n")
    for k in thetas_all:
        t = np.array(thetas_all[k]); m = np.array(mods_all[k])
        if len(t) < 10:
            print(f"{k:16s}  (insufficient complex-pair detections: n={len(t)})")
            continue
        tlo, tmed, thi = np.percentile(t, [16, 50, 84])
        mlo, mmed, mhi = np.percentile(m, [16, 50, 84])
        print(f"{k:16s} {tmed:9.2f} [{tlo:6.2f},{thi:6.2f}] {mmed:10.3f} [{mlo:5.2f},{mhi:5.2f}]  {len(t)}")
    pt = np.array(pooled_theta); pm = np.array(pooled_mod)
    tlo, tmed, thi = np.percentile(pt, [16, 50, 84])
    t025, t975 = np.percentile(pt, [2.5, 97.5])
    mlo, mmed, mhi = np.percentile(pm, [16, 50, 84])
    print(f"\nPOOLED (trusted estimators): theta = {tmed:.1f} deg, "
          f"68% CI [{tlo:.1f},{thi:.1f}], 95% CI [{t025:.1f},{t975:.1f}]")
    print(f"                              |zeta| = {mmed:.3f}, 68% CI [{mlo:.3f},{mhi:.3f}]")
    print(f"  P(theta > 45 deg) = {np.mean(pt > 45)*100:.1f}%   P(theta > 50) = {np.mean(pt>50)*100:.1f}%")
    print(f"  P(|zeta| within 10% of |lambda0|={LAMBDA0_MOD}) = "
          f"{np.mean(np.abs(pm-LAMBDA0_MOD)/LAMBDA0_MOD < 0.10)*100:.1f}%")
    print(f"  P(|zeta| within 20% of |lambda0|) = "
          f"{np.mean(np.abs(pm-LAMBDA0_MOD)/LAMBDA0_MOD < 0.20)*100:.1f}%")
    return pt, pm

# ------------------------------------------------------------------
def drop_v5_test():
    """Compare v0..v4 (drop v5) vs v0..v5 (include v5) phase estimates. Convergence vs drift."""
    print("\n=== DROP-v5 test: does theta converge or drift when v5 is added? ===")
    v2, v3, v4, v5 = 0.100, -0.02, -0.45, -1.1
    v4list = [V0, V1, v2, v3, v4]              # v0..v4 (5 coeffs)
    v5list = [V0, V1, v2, v3, v4, v5]          # v0..v5 (6 coeffs)

    def pade_all(v):
        res = []
        for L, M, pairs in borel_pade_pairs(v):
            cand = [(mod, th) for mod, th in pairs if 15 < abs(th) < 88]
            if cand:
                cand.sort(key=lambda x: x[0])
                res.append((L, M, abs(cand[0][1]), cand[0][0]))
        return res

    print("  v0..v4 (drop v5) Borel-Pade complex pairs:")
    for L, M, th, mod in pade_all(v4list):
        print(f"     [{L}/{M}]: theta={th:6.2f}  |zeta|={mod:.3f}")
    d4 = direct_ratio_pair_3term(v4list)
    if d4 and d4['complex_pair']:
        print(f"     ratio_last2: theta={d4['theta']:6.2f} |zeta|={d4['mod']:.3f}")
    try:
        dar4 = fit_darboux(v4list, nmin=1, alpha_fixed=0.0)
        print(f"     darboux_n1a0: theta={abs(dar4['theta']):6.2f} |zeta|={dar4['mod']:.3f}")
    except Exception as e:
        print(f"     darboux failed: {e}")

    print("  v0..v5 (include v5) Borel-Pade complex pairs:")
    for L, M, th, mod in pade_all(v5list):
        print(f"     [{L}/{M}]: theta={th:6.2f}  |zeta|={mod:.3f}")
    d5 = direct_ratio_pair_3term(v5list)
    if d5 and d5['complex_pair']:
        print(f"     ratio_last2: theta={d5['theta']:6.2f} |zeta|={d5['mod']:.3f}")
    dar5 = fit_darboux(v5list, nmin=1, alpha_fixed=0.0)
    print(f"     darboux_n1a0: theta={abs(dar5['theta']):6.2f} |zeta|={dar5['mod']:.3f}")

    # v5 sweep: how theta moves as v5 sweeps its stated range
    print("  v5 sweep (Pade [1/4] and [3/2], theta_deg / |zeta|):")
    for v5v in [-1.0, -1.1, -1.5, -1.9, -2.29]:
        v = [V0, V1, v2, v3, v4, v5v]
        row = f"     v5={v5v:+.2f}: "
        for L, M, th, mod in pade_all(v):
            row += f"[{L}/{M}]={th:.1f}/{mod:.2f}  "
        print(row)


if __name__ == "__main__":
    drop_v5_test()
    bootstrap(N=4000)
