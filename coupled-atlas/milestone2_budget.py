"""milestone2_budget.py -- error-budget decomposition: which coefficient dominates the phase spread.

Fast bootstrap using ONLY the Pade[1/4],[2/3],[0/5] + ratio_last2 estimators (no slow Darboux).
Vary one coefficient at a time (others at nominal) to attribute the phase variance.
Also give the CONDITIONAL estimate under the 'symbolic v5' prior (v5 in [-1.3,-1.0], trusting -1.1).
"""
import numpy as np
from milestone2_measure import borel_pade_pairs, direct_ratio_pair_3term, LAMBDA0_MOD

V0, V1 = 0.134, 0.111
NOM = dict(v2=0.100, v3=-0.02, v4=-0.45, v5=-1.1)

def ladder(d): return [V0, V1, d['v2'], d['v3'], d['v4'], d['v5']]

def pade_split(v, L, M):
    for LL, MM, pairs in borel_pade_pairs(v):
        if LL == L and MM == M:
            cand = [(mod, th) for mod, th in pairs if 15 < abs(th) < 88]
            if cand:
                cand.sort(key=lambda x: x[0]); return abs(cand[0][1]), cand[0][0]
    return None

def trusted(v):
    ths, mods = [], []
    for L, M in [(1, 4), (2, 3), (0, 5)]:
        r = pade_split(v, L, M)
        if r: ths.append(r[0]); mods.append(r[1])
    d = direct_ratio_pair_3term(v)
    if d and d['complex_pair'] and 15 < d['theta'] < 88:
        ths.append(d['theta']); mods.append(d['mod'])
    return ths, mods

def sample(ranges, N=6000, seed=1):
    rng = np.random.default_rng(seed)
    T, M = [], []
    for _ in range(N):
        d = dict(NOM)
        for k, (lo, hi) in ranges.items():
            d[k] = rng.uniform(lo, hi)
        ths, mods = trusted(ladder(d))
        T += ths; M += mods
    return np.array(T), np.array(M)

FULL = dict(v2=(0.100, 0.104), v3=(-0.04, 0.0), v4=(-0.47, -0.44), v5=(-2.3, -1.0))
ONLY = {k: {k: FULL[k]} for k in FULL}

if __name__ == "__main__":
    print("=== One-at-a-time phase spread (others nominal) ===")
    print(f"{'vary':6s} {'theta med':>9s} {'68% CI':>16s} {'sd(theta)':>10s} {'|zeta| med':>10s}")
    for k in ["v2", "v3", "v4", "v5"]:
        T, M = sample(ONLY[k], N=4000)
        lo, med, hi = np.percentile(T, [16, 50, 84])
        print(f"{k:6s} {med:9.2f} [{lo:6.2f},{hi:6.2f}] {T.std():10.2f} {np.median(M):10.3f}")

    print("\n=== FULL joint bootstrap (all four) ===")
    T, M = sample(FULL, N=8000)
    lo, med, hi = np.percentile(T, [16, 50, 84]); a, b = np.percentile(T, [2.5, 97.5])
    ml, mm, mh = np.percentile(M, [16, 50, 84])
    print(f"theta = {med:.1f}  68%[{lo:.1f},{hi:.1f}]  95%[{a:.1f},{b:.1f}]  sd={T.std():.1f}")
    print(f"|zeta|= {mm:.3f} 68%[{ml:.3f},{mh:.3f}]")
    print(f"P(theta>45)={np.mean(T>45)*100:.1f}%  P(theta>50)={np.mean(T>50)*100:.1f}%")
    print(f"P(|zeta| in 10% of {LAMBDA0_MOD})={np.mean(np.abs(M-LAMBDA0_MOD)/LAMBDA0_MOD<0.1)*100:.1f}%")

    print("\n=== CONDITIONAL: trust symbolic v5 (v5 in [-1.3,-1.0]), others full ===")
    cond = dict(FULL); cond['v5'] = (-1.3, -1.0)
    T, M = sample(cond, N=8000)
    lo, med, hi = np.percentile(T, [16, 50, 84]); a, b = np.percentile(T, [2.5, 97.5])
    ml, mm, mh = np.percentile(M, [16, 50, 84])
    print(f"theta = {med:.1f}  68%[{lo:.1f},{hi:.1f}]  95%[{a:.1f},{b:.1f}]  sd={T.std():.1f}")
    print(f"|zeta|= {mm:.3f} 68%[{ml:.3f},{mh:.3f}]")
    print(f"P(theta>45)={np.mean(T>45)*100:.1f}%  P(theta>50)={np.mean(T>50)*100:.1f}%")

    print("\n=== Nominal point (all at stated central values) ===")
    ths, mods = trusted(ladder(NOM))
    print(f"trusted thetas: {[round(x,1) for x in ths]}  -> mean {np.mean(ths):.1f} sd {np.std(ths):.1f}")
    print(f"trusted |zeta|: {[round(x,3) for x in mods]}  -> mean {np.mean(mods):.3f}")
