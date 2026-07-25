"""
Swallowtail (q=3) trans-series analysis: sign pattern / factorial growth, Domb-Sykes, Borel-Pade
with the real instanton as positive-axis anchor, complex pair, median-Borel vs the beta=2 truth.

Mirrors the cusp paper's Section "Validation: the resummed trans-series reconstructs W":
    x = eta^2 = 4/beta,  f(x) = Var/eta^2 = sum_n v_n x^n,
    B(u) = sum_n (v_n/n!) u^n  ~ diagonal Pade,
    f_pm(x) = int_0^{inf e^{+-i phi}} e^{-t} B(x t) dt,   median = (f_+ + f_-)/2 = Re f_+.
Ground truth: Var(beta=2) = 0.328  =>  f(2) = 0.164.

ANCHOR CORRECTION (instanton_action_q.py, this session): the left-tail FW action is
    I(s) = s^{2q+1} / [2(2q+1)]   ->  q=2: s^5/10 (matches the cusp paper + instanton_action.py),
                                      q=3: s^7/14.
SWALLOWTAIL_TRANSSERIES_NOTES.md quoted s^7/28 (i.e. 1/[4(2q+1)]), which is a factor 2 too small:
at q=2 that formula gives s^5/20, contradicting the numerically verified s^5/10.  Verified here at
both q by the same BVP machinery (ratio-to-notes -> 1.86 at q=2, 1.96 at q=3).

Coefficients are merged from swtl_lowrungs.json (v0..v2) and swtl_results.json (v3..v6), BOTH
produced by swtl_production.py so the ladder is on one normalization (v0 = 0.04969 vs the exact
0.04953187 is the 0.3% grid check).

Run:  python3 swtl_borel.py
"""
import json, os, sys
import numpy as np
from math import factorial, pi
from numpy.polynomial import polynomial as P

HERE = os.path.dirname(os.path.abspath(__file__))
# Ground truth Var(Y*) at beta=2, REPRODUCED from the FP-PDE by swtl_groundtruth.py (0.3324/0.3309
# at dp=0.02/0.01; q=2 regression gate passes).  The notes' asserted 0.328 is confirmed to 1%.
VAR_TRUTH_BETA2 = 0.331
Q = 3
INSTANTON_C = 1.0 / (2 * (2 * Q + 1))   # I(s) = C * s^{2q+1};  q=3 -> 1/14
V0_EXACT = 0.04953187

# ----------------------------------------------------------------- coefficients
def robust_extrap(ns, vs):
    """Return (best, lo, hi, ests, note) for v(n->inf) from several estimators.

    Two failure modes of the campaign's quadratic-in-1/n Richardson are screened here:

      (a) NON-MONOTONE grid sequence -> the 1/n^2 coefficient is fitting noise
          (v4: 0.0480, 0.0456, 0.0468 -> campaign reported +0.0858, ~2x every raw value);
      (b) WRONG-SIDE extrapolation -> the extrapolant sits on the opposite side of the last
          data point from the direction the sequence is actually moving
          (v5: -1.794 -> -1.603 -> -1.538 is RISING, yet the campaign reported -1.962,
          below all three points).

    A credible n->inf limit must lie on the far side of the last grid point, in the
    direction of travel.  Estimators violating that are dropped; the surviving spread is
    reported as the honest band.
    """
    ns = np.asarray(ns, float); vs = np.asarray(vs, float)
    ests = {'raw(last)': vs[-1]}
    resid = {}                                  # per-estimator max |fit - data|
    if len(ns) >= 2:
        A = np.vstack([np.ones(2), 1 / ns[-2:]]).T
        ests['lin 1/n (last 2)'] = float(np.linalg.solve(A, vs[-2:])[0])
    if len(ns) >= 3:
        A = np.vstack([np.ones(3), 1 / ns[-3:]]).T
        c, *_ = np.linalg.lstsq(A, vs[-3:], rcond=None)
        ests['lin 1/n (tail 3)'] = float(c[0])
        resid['lin 1/n (tail 3)'] = float(np.max(np.abs(A @ c - vs[-3:])))
        A = np.vstack([np.ones_like(ns), 1 / ns]).T
        c, *_ = np.linalg.lstsq(A, vs, rcond=None)
        ests['lin 1/n (all)'] = float(c[0]); resid['lin 1/n (all)'] = float(np.max(np.abs(A @ c - vs)))
        A = np.vstack([np.ones_like(ns), 1 / ns, 1 / ns ** 2]).T
        c, *_ = np.linalg.lstsq(A, vs, rcond=None)
        ests['quad (campaign)'] = float(c[0]); resid['quad (campaign)'] = float(np.max(np.abs(A @ c - vs)))

    # PLATEAU / NOISE-FLOOR test.  Under exact 1/n convergence each gap shrinks by a predictable
    # factor; if the observed last gap collapses far below that, the sequence has left the 1/n
    # regime and further extrapolation is fitting the noise floor, not the tail.
    # (v5: gaps +0.0363, +0.0025, +0.0030 over n=24..36 -- 1/n predicts ~+0.027 for the 28->32 step,
    #  10x the observed, and the last two gaps are non-monotone => converged at the ~0.003 level.)
    plateau = False
    if len(ns) >= 4:
        gaps = np.diff(vs)
        h = 1.0 / ns
        # Scan only the TAIL (final four points).  The collapse can happen one step before the
        # end -- v5 collapses at 24->28->32 (0.0363 -> 0.0025 against a 1/n prediction of 0.0272,
        # i.e. 0.09x) and its last two gaps are then both ~noise -- but scanning the whole grid
        # would also fire on the early exit from the pre-asymptotic regime, which is NOT
        # convergence (v6's 10->12->14 drops 9.51 -> 1.18 while its tail gaps, 0.57 and 0.38, are
        # still tapering normally at ~0.8x per step).
        for i in range(max(2, len(vs) - 2), len(vs)):
            denom = h[i - 2] - h[i - 1]
            if denom <= 0:
                continue
            pred = abs(gaps[i - 2]) * ((h[i - 1] - h[i]) / denom)
            if pred > 0 and abs(gaps[i - 1]) < 0.3 * pred:
                plateau = True
                break

    monotone = len(vs) < 3 or (np.all(np.diff(vs) > 0) or np.all(np.diff(vs) < 0))
    d = np.sign(vs[-1] - vs[-2]) if len(vs) >= 2 else 0.0   # direction of travel in n
    travel = abs(vs[-1] - vs[0])                            # total motion actually observed
    flags = []
    keep = {}
    for name, val in ests.items():
        if name != 'raw(last)' and d != 0 and np.sign(val - vs[-1]) == -d:
            flags.append(f'{name} on WRONG SIDE ({val:+.4f} vs last {vs[-1]:+.4f}, trend {"up" if d>0 else "down"})')
            continue
        # (c) BAD FIT: the estimator's own model does not describe the grid data, so its
        # intercept is meaningless.  Screening on *residual* rather than on how far the
        # extrapolant reaches: for genuine 1/n convergence the remaining motion is expected
        # to be n_first/(n_last-n_first) times the observed travel (2x over 24..36, 2.5x
        # over 10..14), so a reach-based cut would wrongly reject the converged low rungs.
        # (v6 'lin 1/n (all)' = +17.13 comes from a line fitted to strongly curved data.)
        if name in resid and travel > 0 and resid[name] > 0.25 * travel:
            flags.append(f'{name} BAD FIT ({val:+.4f}: max resid {resid[name]:.4f} '
                         f'vs travel {travel:.4f})')
            continue
        # (d) SIGN FLIP: every grid value shares a sign, but the extrapolant crosses zero.
        # (v6: all three raw values lie in [-22.5, -11.8], yet 'lin 1/n (all)' lands at +17.1.)
        if name != 'raw(last)' and np.all(vs > 0) != np.all(vs < 0) and np.sign(val) == -np.sign(vs[-1]):
            flags.append(f'{name} SIGN FLIP ({val:+.4f}; all grid values {"positive" if vs[-1]>0 else "negative"})')
            continue
        keep[name] = val
    if not keep:
        keep = {'raw(last)': vs[-1]}
    if plateau:
        # Converged to the noise floor: the raw last value IS the answer, and the band only needs
        # to admit the mild residual drift the tail fit still suggests.
        cands = [vs[-1]] + [keep[k] for k in ('lin 1/n (tail 3)', 'lin 1/n (last 2)') if k in keep]
        best = float(vs[-1])
        note = '; '.join(['PLATEAU (converged to noise floor; extrapolation would fit noise)'] + flags)
        return best, min(cands), max(cands), ests, note
    vals = list(keep.values())
    if monotone and 'quad (campaign)' in keep:
        best = keep['quad (campaign)']
    elif 'lin 1/n (tail 3)' in keep:
        best = keep['lin 1/n (tail 3)']
    elif 'lin 1/n (last 2)' in keep:
        best = keep['lin 1/n (last 2)']
    else:
        best = float(np.median(vals))
    note = '; '.join(([] if monotone else ['NON-MONOTONE grid']) + flags)
    return best, min(vals), max(vals), ests, note


def _max_idx(k):
    """Highest Wiener-chaos index Y_idx that v_k's assembly needs."""
    try:
        from _v6_driver import vk_terms
        return max(p for pr in vk_terms(k) for p in pr[:2])
    except Exception:
        return 0


def load_ladder(verbose=True):
    """Assemble the ladder, DISCARDING grid points with n < max chaos index.

    A grid of n points cannot resolve the functional Y_idx when idx > n, and such points are not
    merely noisy -- they are systematically wrong.  The evidence is sharp: v6 needs idx 13 and its
    n=10,12 points swing 74% (-22.45, -12.94) while n=14..20 taper smoothly to -9.18; v5 (idx 11)
    and v4 (idx 9) were computed entirely above their thresholds and converge cleanly.  Keeping the
    under-resolved points corrupts every global fit through them.
    """
    v, band, notes = {}, {}, {}
    for fn in ('swtl_lowrungs.json', 'swtl_results.json'):
        p = os.path.join(HERE, fn)
        if not os.path.exists(p):
            continue
        with open(p) as f:
            d = json.load(f)
        for k, grids in d.items():
            got = sorted((int(n), val) for n, val in grids.items() if n != 'extrap_ninf')
            need = _max_idx(int(k))
            dropped = [n for n, _ in got if n < need]
            got = [(n, val) for n, val in got if n >= need]
            if dropped and verbose:
                print(f"  [v{k}] dropped under-resolved grids n={dropped} (need n >= max chaos idx {need})")
            if not got:
                if verbose:
                    print(f"  [v{k}] NO usable grid points (all below idx {need}) -- excluded from the ladder")
                continue
            best, lo, hi, ests, note = robust_extrap([g[0] for g in got], [g[1] for g in got])
            v[int(k)] = best; band[int(k)] = (lo, hi); notes[int(k)] = (got, ests, note)
    if verbose:
        print("=" * 78)
        print("LADDER  (production engine, q=3)   Var/eta^2 = sum v_k (eta^2)^k")
        print("=" * 78)
        for k in sorted(v):
            got, ests, note = notes[k]
            print(f"  v{k} = {v[k]:+.5f}   band [{band[k][0]:+.5f}, {band[k][1]:+.5f}]"
                  f"   grid {[(n, round(x, 6)) for n, x in got]}")
            print(f"        estimators: " + ", ".join(f"{a}={b:+.5f}" for a, b in ests.items())
                  + (f"   <<< {note}" if note else ""))
        if 0 in v:
            print(f"\n  NORMALIZATION CHECK: v0 = {v[0]:+.6f} vs exact {V0_EXACT:.8f}"
                  f"  ratio {v[0]/V0_EXACT:.4f}  "
                  f"{'OK (grid-level)' if abs(v[0]/V0_EXACT - 1) < 0.02 else '*** MISMATCH ***'}")
    ks = sorted(v)
    contiguous = [v[k] for k in range(0, max(ks) + 1) if k in v] if ks == list(range(len(ks))) else None
    return v, band, contiguous


# ----------------------------------------------------------------- (1) signs / growth
def sign_and_growth(v):
    print("\n" + "=" * 78)
    print("(1) SIGN PATTERN AND FACTORIAL GROWTH")
    print("=" * 78)
    print("  signs: " + " ".join('+' if x > 0 else '-' for x in v))
    print(f"  {'k':>2} {'v_k':>12} {'v_k/v_{k-1}':>12} {'ratio/k':>10} {'b_k=v_k/k!':>12} {'b_k/b_{k-1}':>12}")
    b = [v[n] / factorial(n) for n in range(len(v))]
    for k in range(len(v)):
        r = v[k] / v[k - 1] if k and v[k - 1] != 0 else np.nan
        rb = b[k] / b[k - 1] if k and b[k - 1] != 0 else np.nan
        print(f"  {k:>2} {v[k]:>+12.5f} {r:>12.4f} {r/k if k else np.nan:>10.4f} "
              f"{b[k]:>+12.5f} {rb:>+12.4f}")
    print("\n  Reading: factorial (Gevrey-1) growth => v_k/v_{k-1} ~ (k+alpha-1)/zeta grows LINEARLY in k,")
    print("  equivalently b_k/b_{k-1} -> 1/zeta (a constant).  A sign flip / non-monotone ratio")
    print("  signals a COMPLEX pair (oscillation cos(k*theta - phi)) rather than a real singularity.")


# ----------------------------------------------------------------- (2) Domb-Sykes
def domb_sykes(v):
    print("\n" + "=" * 78)
    print("(2) DOMB-SYKES")
    print("=" * 78)
    b = np.array([v[n] / factorial(n) for n in range(len(v))])
    ks, rs = [], []
    for k in range(1, len(b)):
        if b[k - 1] != 0:
            ks.append(k); rs.append(b[k] / b[k - 1])
    ks = np.array(ks, float); rs = np.array(rs)
    print("  On Borel coefficients b_k = v_k/k!  (r_k -> 1/|zeta| with slope giving the exponent):")
    for k, r in zip(ks, rs):
        print(f"    k={k:.0f}  r_k={r:+.5f}   1/k={1/k:.4f}   implied |zeta|={abs(1/r) if r else np.inf:.4f}")
    good = np.isfinite(rs) & (np.abs(rs) > 1e-14)
    if good.sum() >= 2:
        sl, ic = np.polyfit(1 / ks[good], rs[good], 1)
        print(f"\n  Linear fit r_k = {ic:+.5f} + {sl:+.5f}/k  =>  1/zeta = {ic:+.5f}"
              f"  =>  zeta = {1/ic if ic else np.inf:+.4f}")
        print(f"  exponent:  slope/intercept = -(1+g) => g = {-(sl/ic) - 1 if ic else np.nan:+.4f}")
        if ic < 0:
            print("  NOTE: negative intercept => nearest singularity on the NEGATIVE axis or a")
            print("        complex pair with |theta|>90deg; Domb-Sykes assumes a real positive one,")
            print("        so treat this as evidence for the complex pair, not a literal zeta.")
    print("\n  Raw-coefficient ratios v_k/v_{k-1} vs k (linear growth <=> factorial divergence):")
    kk = np.arange(1, len(v), dtype=float)
    rr = np.array([v[k] / v[k - 1] if v[k - 1] else np.nan for k in range(1, len(v))])
    ok = np.isfinite(rr)
    if ok.sum() >= 2:
        sl2, ic2 = np.polyfit(kk[ok], rr[ok], 1)
        print(f"    fit v_k/v_{{k-1}} = {sl2:+.4f} k + {ic2:+.4f}  (slope 1/zeta_eff = {sl2:+.4f}"
              f" => |zeta_eff| = {abs(1/sl2) if sl2 else np.inf:.3f})")


# ----------------------------------------------------------------- Pade / Borel
def pade(c, L, M):
    c = np.asarray(c, float)
    A = np.zeros((M, M)); rhs = np.zeros(M)
    for i in range(1, M + 1):
        for j in range(1, M + 1):
            k = L + i - j
            A[i - 1, j - 1] = c[k] if k >= 0 else 0.0
        rhs[i - 1] = -c[L + i]
    q = np.concatenate([[1.0], np.linalg.solve(A, rhs)])
    p = np.array([sum(c[k - j] * q[j] for j in range(0, min(k, M) + 1)) for k in range(L + 1)])
    return p, q


def borel_pade_poles(v, L, M):
    b = [v[n] / factorial(n) for n in range(len(v))]
    _, q = pade(b, L, M)
    return sorted(P.polyroots(q), key=abs)


def report_poles(v, s_anchor_range=(1.0, 3.0)):
    print("\n" + "=" * 78)
    print("(3) BOREL-PADE POLES  +  COMPLEX PAIR  +  REAL INSTANTON ANCHOR")
    print("=" * 78)
    A_of_s = lambda s: INSTANTON_C * s ** (2 * Q + 1)
    s_of_A = lambda A: (A / INSTANTON_C) ** (1.0 / (2 * Q + 1))
    print(f"  Real FW instanton (CORRECTED): I(s) = s^{2*Q+1}/{1/INSTANTON_C:.0f}"
          f"   [notes said /{4*(2*Q+1)} -- factor 2 too small, see instanton_action_q.py]")
    print(f"  A positive-axis Borel pole at u = zeta_r corresponds to escape depth"
          f" s = (14 zeta_r)^(1/7):")
    for A in (0.5, 1.0, 1.5, 2.0, 3.0):
        print(f"    zeta_r = {A:.2f}  <->  s = {s_of_A(A):.3f}")
    print(f"  and conversely s in [{s_anchor_range[0]}, {s_anchor_range[1]}] "
          f"<-> zeta_r in [{A_of_s(s_anchor_range[0]):.3f}, {A_of_s(s_anchor_range[1]):.3f}]")
    K = len(v) - 1
    print(f"\n  Poles of [L/M] Pade to B(u) = sum v_n u^n/n!   ({K+1} coefficients):")
    found_pair, found_real = [], []
    for M in range(1, K + 1):
        L = K - M
        if L < 0:
            continue
        try:
            roots = borel_pade_poles(v, L, M)
        except np.linalg.LinAlgError:
            print(f"    [{L}/{M}]: singular"); continue
        desc = ", ".join(f"{abs(r):.3f}@{np.degrees(np.angle(r)):+.1f}deg" for r in roots)
        print(f"    [{L}/{M}]: {desc}")
        for r in roots:
            th = np.degrees(np.angle(r))
            if abs(th) < 12 and r.real > 0:
                found_real.append(abs(r))
            elif 12 <= abs(th) <= 168:
                found_pair.append((abs(r), abs(th)))
    if found_real:
        med = float(np.median(found_real))
        print(f"\n  REAL POSITIVE poles: {[round(x,3) for x in sorted(found_real)]}")
        print(f"    median |zeta_r| = {med:.3f}  ->  implied instanton depth s = {s_of_A(med):.3f}"
              f"   (I(s) = {med:.3f})")
    else:
        print("\n  No real positive pole isolated at this order.")
    if found_pair:
        mods = np.array([p[0] for p in found_pair]); ths = np.array([p[1] for p in found_pair])
        print(f"  COMPLEX PAIR candidates: |zeta| median {np.median(mods):.3f} "
              f"(range {mods.min():.3f}-{mods.max():.3f}), "
              f"theta median {np.median(ths):.1f}deg (range {ths.min():.1f}-{ths.max():.1f})")
    else:
        print("  No complex pair isolated at this order.")
    return found_real, found_pair


# ----------------------------------------------------------------- (4) median Borel
def median_borel(v, x, L=None, M=None, phi_deg=35.0, Rmax=60.0, npts=6000):
    """f(x) = int_0^{inf e^{i phi}} e^{-t} B_pade(x t) dt ; median = Re f_+ (real coeffs)."""
    K = len(v) - 1
    if M is None:
        M = K // 2
    if L is None:
        L = K - M
    b = [v[n] / factorial(n) for n in range(len(v))]
    p, q = pade(b, L, M)
    phi = np.radians(phi_deg)
    r = np.linspace(1e-9, Rmax, npts)
    t = r * np.exp(1j * phi)
    u = x * t
    num = np.polyval(p[::-1], u); den = np.polyval(q[::-1], u)
    integ = np.exp(-t) * (num / den) * np.exp(1j * phi)
    return float(np.real(np.trapz(integ, r))), (L, M)


def resummation(v):
    print("\n" + "=" * 78)
    print("(4) MEDIAN-BOREL RESUMMATION vs GROUND TRUTH")
    print("=" * 78)
    f_truth = VAR_TRUTH_BETA2 / 2.0
    print(f"  Ground truth: Var(beta=2) = {VAR_TRUTH_BETA2}  =>  f(2) = Var/eta^2 = {f_truth:.4f}")
    K = len(v) - 1
    # HEADLINE = the W paper's actual protocol: the DIAGONAL Pade (for 7 coefficients, [3/3]),
    # with "median" meaning the median of the two LATERAL Borel sums, (f_+ + f_-)/2 = Re f_+.
    # It does NOT mean a median over Pade orders -- that is a different, much worse statistic
    # (it returns 0.136 on the cusp at x=2, against the published 0.235 and truth 0.237).
    print(f"\n  DIAGONAL [{K//2}/{K - K//2}] Pade, median of lateral Borel sums (paper protocol):")
    print(f"  {'x':>5} {'beta':>6} {'naive':>12} {'diagonal':>10}   (phi-spread)")
    for x, beta in ((1.0, 4.0), (1.5, 8 / 3), (2.0, 2.0)):
        naive = sum(v[n] * x ** n for n in range(len(v)))
        vals = []
        for phi in (20.0, 25.0, 30.0, 35.0, 45.0):
            try:
                f, _ = median_borel(v, x, L=K - K // 2, M=K // 2, phi_deg=phi)
                if np.isfinite(f):
                    vals.append(f)
            except (np.linalg.LinAlgError, ValueError):
                pass
        if vals:
            a = np.array(vals)
            mark = f"   <-- truth {f_truth:.4f}" if x == 2.0 else ""
            print(f"  {x:>5.2f} {beta:>6.2f} {naive:>12.3f} {np.median(a):>10.4f}   "
                  f"[{a.min():+.4f},{a.max():+.4f}]{mark}")
    fdiag, _ = median_borel(v, 2.0, L=K - K // 2, M=K // 2, phi_deg=35.0)
    print(f"\n  >>> DIAGONAL f(2) = {fdiag:+.4f}  vs truth {f_truth:.4f}  "
          f"=> Var = {2*fdiag:+.4f} vs {VAR_TRUTH_BETA2}  (rel err {abs(2*fdiag/VAR_TRUTH_BETA2-1)*100:.1f}%)")

    print(f"\n  --- secondary: median over ALL [L/M] and phi (a robustness check, NOT the protocol) ---")
    print(f"  {'x':>5} {'beta':>6} {'naive':>12} {'median resum (various [L/M], phi)':>12}")
    for x, beta in ((1.0, 4.0), (1.5, 8/3), (2.0, 2.0)):
        naive = sum(v[n] * x ** n for n in range(len(v)))
        vals = []
        for M in range(1, K + 1):
            L = K - M
            if L < 0 or M < 1:
                continue
            for phi in (25.0, 35.0, 45.0, 60.0):
                try:
                    f, _ = median_borel(v, x, L=L, M=M, phi_deg=phi)
                    if np.isfinite(f) and abs(f) < 1e4:
                        vals.append((f, L, M, phi))
                except (np.linalg.LinAlgError, ValueError):
                    pass
        if vals:
            arr = np.array([a[0] for a in vals])
            med = float(np.median(arr))
            iqr = (float(np.percentile(arr, 25)), float(np.percentile(arr, 75)))
            print(f"  {x:>5.2f} {beta:>6.2f} {naive:>12.3f}   median over {len(vals)} fits = {med:+.4f}"
                  f"   IQR [{iqr[0]:+.4f},{iqr[1]:+.4f}]   full range [{arr.min():+.3f},{arr.max():+.3f}]")
            if x == 2.0:
                print(f"\n  >>> f(2) resummed = {med:+.4f}   vs ground truth {f_truth:.4f}"
                      f"   =>  Var = {2*med:+.4f} vs {VAR_TRUTH_BETA2}"
                      f"   (rel err {abs(2*med/VAR_TRUTH_BETA2 - 1)*100:.1f}%)")
                spread = (iqr[1] - iqr[0]) / max(abs(med), 1e-9)
                print(f"      IQR spread / |median| = {spread*100:.0f}%  ->  "
                      f"{'STABLE' if spread < 0.25 else 'SCATTERED: not a reconstruction'}")
        else:
            print(f"  {x:>5.2f} {beta:>6.2f} {naive:>12.3f}   (no usable Pade)")


# ----------------------------------------------------------------- main
def anchored_darboux(v, label="q=3", zr_grid=None, alphas=(-0.5, 0.0, 0.5), verbose=True):
    """Instanton-ANCHORED Darboux/Dingle fit: complex pair + a real term whose position is FIXED.

    Raw Borel-Pade poles mix the complex pair with the real Freidlin-Wentzell singularity, and that
    contamination is exactly what biased the cusp's early (retracted) theta estimates.  Here the
    real term sits at a PROFILED z_r (scanned, not fitted), so the pair is extracted at each
    assumed instanton position:

        v_n ~ 2C |zeta|^-(n+1) Gamma(n+1+a) cos((n+1)theta - phi)  +  C_r z_r^-(n+1) Gamma(n+1+a)

    Five free parameters (C, |zeta|, theta, phi, C_r) on 7 points, with a profiled instead of
    fitted -- over-fit protection, as the cusp paper describes.  theta is trustworthy only if it is
    flat across the (z_r, alpha) profile.
    """
    from scipy.optimize import least_squares
    from scipy.special import gamma as Gamma
    ns = np.arange(len(v), dtype=float)
    y = np.asarray(v, float)
    scale = np.abs(y) + 1e-3
    if zr_grid is None:
        zr_grid = [0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 3.0]

    def model(p, zr, a):
        C, zeta, theta, phi, Cr = p
        pair = 2 * C * zeta ** (-(ns + 1)) * Gamma(ns + 1 + a) * np.cos((ns + 1) * theta - phi)
        real = Cr * zr ** (-(ns + 1)) * Gamma(ns + 1 + a)
        return pair + real

    rows = []
    for zr in zr_grid:
        for a in alphas:
            best = None
            for th0 in np.radians([25, 35, 45, 55, 65]):
                for z0 in (0.8, 1.2, 1.8, 2.5):
                    p0 = [0.05, z0, th0, 0.5, 0.0]
                    lo = [1e-8, 0.3, np.radians(5), -2 * pi, -50.0]
                    hi = [50.0, 6.0, np.radians(120), 2 * pi, 50.0]
                    try:
                        r = least_squares(lambda p: (model(p, zr, a) - y) / scale, p0,
                                          bounds=(lo, hi), xtol=1e-14, ftol=1e-14)
                    except ValueError:
                        continue
                    if best is None or r.cost < best.cost:
                        best = r
            if best is None:
                continue
            C, zeta, theta, phi, Cr = best.x
            rows.append((zr, a, zeta, np.degrees(theta), Cr, float(np.max(np.abs(best.fun)))))
    if verbose:
        print("\n" + "=" * 78)
        print(f"(7) INSTANTON-ANCHORED DARBOUX FIT ({label})")
        print("=" * 78)
        print(f"  real term pinned at z_r (profiled); s = (14 z_r)^(1/7) is the implied escape depth")
        print(f"  {'z_r':>6} {'s':>6} {'alpha':>6} {'|zeta|':>8} {'theta':>8} {'C_r':>9} {'maxresid':>9}")
        for zr, a, zeta, th, Cr, res in rows:
            s = (zr / INSTANTON_C) ** (1.0 / (2 * Q + 1))
            print(f"  {zr:>6.2f} {s:>6.3f} {a:>6.1f} {zeta:>8.3f} {th:>8.1f} {Cr:>+9.3f} {res:>9.3f}")
    if rows:
        ths = np.array([r[3] for r in rows]); zs = np.array([r[2] for r in rows])
        # theta is near-degenerate with alpha but almost independent of z_r, so quoting one
        # pooled spread hides the structure: report the z_r-spread AT FIXED alpha (that is what
        # the anchor actually constrains) and carry the alpha dependence as a systematic.
        if verbose:
            print(f"\n  theta vs the ANCHOR, at fixed alpha (this is what z_r pins):")
        for a in alphas:
            sel = [r for r in rows if r[1] == a and r[5] < 0.5]
            if not sel:
                continue
            t = np.array([r[3] for r in sel]); z = np.array([r[2] for r in sel])
            if verbose:
                print(f"    alpha={a:+.1f}: theta = {np.median(t):5.1f} deg  "
                      f"(z_r-spread {t.max()-t.min():.1f} deg)   |zeta| = {np.median(z):.2f}")
        good = np.array([r[5] for r in rows]) < 0.5
        use_t, use_z = (ths[good], zs[good]) if good.sum() >= 3 else (ths, zs)
        a0 = [r for r in rows if r[1] == 0.0 and r[5] < 0.5]
        t0 = np.array([r[3] for r in a0]) if a0 else use_t
        if verbose:
            print(f"\n  HEADLINE (alpha=0): theta = {np.median(t0):.1f} deg, "
                  f"z_r-spread {t0.max()-t0.min():.1f} deg"
                  f"  ->  {'PINNED by the anchor' if t0.max()-t0.min() < 5 else 'NOT pinned'}")
            print(f"  systematic from alpha in {alphas}: theta ranges {use_t.min():.1f}-{use_t.max():.1f} deg")
        return np.median(t0), use_t.min(), use_t.max()
    return np.nan, np.nan, np.nan


def sensitivity(v, band):
    """Push the per-coefficient extrapolation bands through the pair location and f(2).

    The top rungs (v4, v5, v6) are the least converged, and they are exactly the ones that
    control the large-order behaviour -- so a central-value-only answer would overstate what
    the ladder supports.  Corner-sample the bands and report the induced spread.
    """
    print("\n" + "=" * 78)
    print("(5) SENSITIVITY TO THE EXTRAPOLATION BANDS")
    print("=" * 78)
    K = len(v) - 1
    top = [k for k in range(max(0, K - 2), K + 1)]
    print(f"  Corner-sampling v{top[0]}..v{top[-1]} over their bands "
          f"({2**len(top)} corners); lower rungs held at central values.")
    import itertools
    f_truth = VAR_TRUTH_BETA2 / 2.0
    rows = []
    for corner in itertools.product(*[[band[k][0], band[k][1]] for k in top]):
        vv = list(v)
        for k, val in zip(top, corner):
            vv[k] = val
        # pair location
        mods, ths = [], []
        for M in range(1, K + 1):
            L = K - M
            if L < 0:
                continue
            try:
                for r in borel_pade_poles(vv, L, M):
                    th = abs(np.degrees(np.angle(r)))
                    if 12 <= th <= 168:
                        mods.append(abs(r)); ths.append(th)
            except np.linalg.LinAlgError:
                pass
        # f(2)
        fs = []
        for M in range(1, K + 1):
            L = K - M
            if L < 0:
                continue
            for phi in (25.0, 35.0, 45.0, 60.0):
                try:
                    f, _ = median_borel(vv, 2.0, L=L, M=M, phi_deg=phi)
                    if np.isfinite(f) and abs(f) < 1e4:
                        fs.append(f)
                except (np.linalg.LinAlgError, ValueError):
                    pass
        rows.append((corner,
                     np.median(mods) if mods else np.nan,
                     np.median(ths) if ths else np.nan,
                     np.median(fs) if fs else np.nan))
    print(f"  {'corner (v' + ','.join('v'+str(k) for k in top) + ')':>34} {'|zeta|':>8} {'theta':>8} {'f(2)':>9}")
    for corner, m, t, f in rows:
        cs = ",".join(f"{c:+.3f}" for c in corner)
        print(f"  {cs:>34} {m:>8.3f} {t:>8.1f} {f:>9.4f}")
    arr_m = np.array([r[1] for r in rows]); arr_t = np.array([r[2] for r in rows])
    arr_f = np.array([r[3] for r in rows])
    ok = np.isfinite(arr_t)
    if ok.any():
        print(f"\n  ACROSS BANDS:  |zeta| = {np.nanmin(arr_m):.2f}-{np.nanmax(arr_m):.2f}"
              f"   theta = {np.nanmin(arr_t):.1f}-{np.nanmax(arr_t):.1f} deg"
              f"   f(2) = {np.nanmin(arr_f):.3f}-{np.nanmax(arr_f):.3f}  (truth {f_truth:.3f})")
        print(f"  => theta is {'insensitive' if np.nanmax(arr_t)-np.nanmin(arr_t) < 8 else 'sensitive'}"
              f" to the extrapolation bands; f(2) is "
              f"{'insensitive' if np.nanmax(arr_f)-np.nanmin(arr_f) < 0.05 else 'sensitive'} to them.")
        print("  CAUTION: insensitivity to the BANDS is not correctness. It says the grid")
        print("  extrapolation is not the limiting error -- it says nothing about truncation")
        print("  error in K (see section 6) or about agreement with the ground truth.")


# Cusp (q=2) ladder, W paper eq:coeffs -- used to validate this file's Borel machinery against a
# published answer before trusting it at q=3.  Paper: median resummation gives f(2)=0.235 vs the
# FP-PDE ground truth 0.237.
CUSP_V = [0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90]
CUSP_TRUTH = {1.0: 0.232, 1.5: 0.239, 2.0: 0.237, 2.25: 0.234}
CUSP_PAPER = {1.0: 0.224, 1.5: 0.234, 2.0: 0.235, 2.25: 0.234}


def validate_against_cusp():
    print("=" * 78)
    print("VALIDATION: same machinery on the CUSP (q=2) ladder, where the answer is published")
    print("=" * 78)
    print(f"  cusp v0..v6 = {CUSP_V}")
    K = len(CUSP_V) - 1
    print("\n  Borel-Pade poles (paper: complex pair at theta~50deg, |zeta|~1.9):")
    for M in range(1, K + 1):
        L = K - M
        try:
            roots = borel_pade_poles(CUSP_V, L, M)
        except np.linalg.LinAlgError:
            continue
        print(f"    [{L}/{M}]: " + ", ".join(
            f"{abs(r):.3f}@{np.degrees(np.angle(r)):+.1f}deg" for r in roots))
    print(f"\n  {'x':>5} {'diag [3/3] phi=35':>18} {'paper':>8} {'truth':>8}")
    ok = True
    for x in (1.0, 1.5, 2.0, 2.25):
        f, LM = median_borel(CUSP_V, x, L=3, M=3, phi_deg=35.0)
        print(f"  {x:>5.2f} {f:>18.4f} {CUSP_PAPER[x]:>8.3f} {CUSP_TRUTH[x]:>8.3f}")
        if abs(f - CUSP_PAPER[x]) > 0.03:
            ok = False
    print(f"\n  => machinery {'REPRODUCES' if ok else 'DOES NOT reproduce'} the published cusp"
          f" resummation; q=3 results {'are' if ok else 'are NOT'} trustworthy on code grounds.")
    # the anchored fit must also recover the cusp's published theta ~ 50 deg
    med, lo, hi = anchored_darboux(CUSP_V, label="CUSP q=2 -- paper says theta = 50 +- 2 deg",
                                   zr_grid=[1.0, 1.5, 1.9, 2.5, 3.8], alphas=(-0.5, 0.0, 0.5))
    print(f"  => anchored fit on the cusp gives theta = {med:.1f} deg (paper: 50 +- 2); "
          f"{'CONSISTENT' if abs(med - 50) < 8 else 'INCONSISTENT -- fit not trustworthy'}")
    return ok


def order_stability(v):
    """theta and f(2) as a function of how many coefficients are used.

    This is the diagnostic the cusp paper needed and lacked: its 5/6-coefficient reading
    'theta ~ 54-63deg, robustly != 45deg' was RETRACTED once the 7th coefficient landed.
    Stability against the extrapolation bands says nothing about stability in K.
    """
    print("\n" + "=" * 78)
    print("(6) ORDER STABILITY (theta, f(2) vs number of coefficients)")
    print("=" * 78)
    print(f"  {'K (uses v0..vK)':>16} {'pair |zeta|':>12} {'pair theta':>12} {'f(2) median':>12}")
    for K in range(3, len(v)):
        vv = v[:K + 1]
        mods, ths, fs = [], [], []
        for M in range(1, K + 1):
            L = K - M
            if L < 0:
                continue
            try:
                for r in borel_pade_poles(vv, L, M):
                    th = abs(np.degrees(np.angle(r)))
                    if 12 <= th <= 168:
                        mods.append(abs(r)); ths.append(th)
            except np.linalg.LinAlgError:
                pass
            for phi in (25.0, 35.0, 45.0, 60.0):
                try:
                    f, _ = median_borel(vv, 2.0, L=L, M=M, phi_deg=phi)
                    if np.isfinite(f) and abs(f) < 1e4:
                        fs.append(f)
                except (np.linalg.LinAlgError, ValueError):
                    pass
        print(f"  {K:>16} {np.median(mods) if mods else np.nan:>12.3f} "
              f"{np.median(ths) if ths else np.nan:>12.1f} {np.median(fs) if fs else np.nan:>12.4f}")
    print("  If theta wanders by >10deg as K increases, it is NOT determined by this ladder.")


if __name__ == "__main__":
    if '--validate' in sys.argv:
        validate_against_cusp()
        print()
    v, band, contiguous = load_ladder()
    if contiguous is None or len(contiguous) < 3:
        print("\nNot enough contiguous coefficients yet; rerun when the campaign has advanced.")
        sys.exit(0)
    print(f"\n  -> using contiguous ladder v0..v{len(contiguous)-1}")
    sign_and_growth(contiguous)
    domb_sykes(contiguous)
    report_poles(contiguous)
    resummation(contiguous)
    anchored_darboux(contiguous)
    sensitivity(contiguous, band)
    order_stability(contiguous)
