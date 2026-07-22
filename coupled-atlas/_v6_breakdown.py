"""Compute v6 and separately the contribution of boundary+order>5 moments (the convention-caveat
set), so v6 can be reported as (solid pure-bulk + dense-boundary) plus a quantified correction."""
import sys, time
import chaos_diagram as CD
import chaos_transfer as CT

def run(n):
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
    u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
    Us, s, Vs, Ye = CD.load_Y('_yexprs_13.txt')
    Y = {k: CD.Ybase(k, Us, s, Vs, Ye, u0ps, Vst) for k in range(1, 14)}
    pairs = [(7, 7, 1)] + [(j, 14 - j, 2) for j in range(1, 7)]

    full = {}   # cache moment
    bnd_hi = 0.0     # accumulated contribution of boundary+maxk>5 moments (with their coef)
    def mom(atoms):
        return CT.moment_hybrid(atoms)

    def is_bnd_hi(key):
        if not any(a[0] == 's' for a in key):
            return False
        mk = max((a[1] for a in key if a[0] == 'U'), default=0)
        return mk > CT.DENSE_MAXORD

    # accumulate v6 and the bnd-hi part; track mean() pieces too
    tot = 0.0; tot_bhi = 0.0
    t0 = time.time()
    for (a, b, coef) in pairs:
        A = Y[a]; B = Y[b]
        # prod_moment - mean*mean
        pm = 0.0; pm_bhi = 0.0
        for ca, aa in A:
            for cb, ab in B:
                if ca == 0.0 or cb == 0.0:
                    continue
                atoms = aa + ab
                key = tuple(sorted(atoms))
                m = mom(atoms)
                pm += ca * cb * m
                if is_bnd_hi(key):
                    pm_bhi = pm_bhi + ca * cb * m
        # means (a==b gives Var = prod - mean^2; cross Cov = prod - meanA*meanB)
        mA = sum(c * mom(at) for c, at in A)
        mB = mA if a == b else sum(c * mom(at) for c, at in B)
        val = pm - mA * mB
        tot += coef * val
        tot_bhi += coef * pm_bhi
        print(f"   ({a},{b})x{coef}: {val:+.5f}   bnd_hi_part={pm_bhi:+.5f}  [{time.time()-t0:.0f}s]", flush=True)
    print(f"v6(n={n}) = {tot:+.5f}   | boundary+order>5 part = {tot_bhi:+.5f}  ({100*tot_bhi/tot:.1f}% of v6)", flush=True)
    return tot, tot_bhi

if __name__ == "__main__":
    for n in [int(x) for x in sys.argv[1:]] or [14]:
        run(n)
