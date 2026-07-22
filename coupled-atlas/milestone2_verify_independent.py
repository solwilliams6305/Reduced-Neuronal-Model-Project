"""Independent adversarial re-check of the Borel phase, written from scratch.
Does NOT import the agent's modules. Uses mpmath for high precision Pade to rule out
conditioning artifacts, and an independent bootstrap.
"""
import numpy as np
from math import factorial
import mpmath as mp
mp.mp.dps = 40

V_TASK = [0.134, 0.111, 0.100, -0.02, -0.45, -1.1]
V_BASE = [0.134, 0.111, 0.104, -0.030, -0.449, -1.1]

def borel(v):
    return [mp.mpf(str(v[n]))/mp.factorial(n) for n in range(len(v))]

def pade_poles_mp(v, L, M):
    """High-precision Pade [L/M] of the Borel series; return complex singularities zeta=1/root(q)."""
    b = borel(v)
    # solve for denominator q (q0=1): sum_{j=1..M} q_j b_{L+i-j} = -b_{L+i}, i=1..M
    A = mp.matrix(M, M); rhs = mp.matrix(M, 1)
    for i in range(1, M+1):
        for j in range(1, M+1):
            k = L + i - j
            A[i-1, j-1] = b[k] if k >= 0 else mp.mpf(0)
        rhs[i-1] = -b[L+i]
    qtail = mp.lu_solve(A, rhs)
    qcoef = [mp.mpf(1)] + [qtail[i] for i in range(M)]
    # polynomial q(t) = sum qcoef[k] t^k ; its roots t are Borel singularities (=zeta)
    roots = mp.polyroots(qcoef[::-1], maxsteps=200, extraprec=100) if M>0 else []
    out = []
    for r in roots:
        rc = complex(r)
        if abs(rc) > 1e-9 and abs(rc) < 100:
            out.append(rc)
    return sorted(out, key=abs)

def nearest_complex(v, L, M):
    try:
        roots = pade_poles_mp(v, L, M)
    except Exception:
        return None
    cand = [r for r in roots if 15 < abs(np.degrees(np.angle(r))) < 88]
    if not cand:
        return None
    cand.sort(key=abs)
    r = cand[0]
    return abs(np.degrees(np.angle(r))), abs(r)

print("=== Independent high-precision Borel-Pade (mpmath dps=40) ===")
for name, v in [("TASK", V_TASK), ("BASE", V_BASE)]:
    print(f" {name} ladder:")
    K = len(v)-1
    for M in range(1, K+1):
        L = K - M
        if L < 0: continue
        nc = nearest_complex(v, L, M)
        if nc:
            print(f"   [{L}/{M}] complex pair: theta={nc[0]:.2f}  |zeta|={nc[1]:.3f}")
        else:
            print(f"   [{L}/{M}] no complex pair in 15-88deg")

# Independent bootstrap using ONLY high-precision Pade[1/4],[2/3],[0/5] (drop [3/2] which is 2-pole only)
print("\n=== Independent bootstrap (mpmath Pade, N=1500) ===")
rng = np.random.default_rng(2024)
V0, V1 = 0.134, 0.111
def sample_theta(ranges, N):
    T, Mo = [], []
    for _ in range(N):
        v2 = rng.uniform(*ranges['v2']); v3 = rng.uniform(*ranges['v3'])
        v4 = rng.uniform(*ranges['v4']); v5 = rng.uniform(*ranges['v5'])
        v = [V0, V1, v2, v3, v4, v5]
        for L, M in [(1,4),(2,3),(0,5)]:
            nc = nearest_complex(v, L, M)
            if nc and np.isfinite(nc[0]) and np.isfinite(nc[1]) and 0<nc[1]<20:
                T.append(nc[0]); Mo.append(nc[1])
    return np.array(T), np.array(Mo)

FULL = dict(v2=(0.100,0.104), v3=(-0.04,0.0), v4=(-0.47,-0.44), v5=(-2.3,-1.0))
COND = dict(FULL); COND['v5'] = (-1.3,-1.0)
for tag, r in [("FULL v5 range", FULL), ("trust symbolic v5", COND)]:
    T, Mo = sample_theta(r, 1500)
    lo,med,hi = np.percentile(T,[16,50,84]); a,b=np.percentile(T,[2.5,97.5])
    print(f" {tag:22s}: theta={med:.1f} 68%[{lo:.1f},{hi:.1f}] 95%[{a:.1f},{b:.1f}] "
          f"P(>45)={np.mean(T>45)*100:.0f}% |zeta|med={np.median(Mo):.3f} "
          f"P(|z|in10%of1.258)={np.mean(np.abs(Mo-1.258)/1.258<0.1)*100:.1f}%")
