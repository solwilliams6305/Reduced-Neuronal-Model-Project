"""Driver: compute v_k via the transfer engine (chaos_transfer.moment_transfer) and
cross-check against the dense engine where feasible."""
import sys, time
import chaos_diagram as CD
import chaos_transfer as CT

def prod_moment(A, B, mom):
    tot = 0.0
    for ca, aa in A:
        for cb, ab in B:
            if ca == 0.0 or cb == 0.0:
                continue
            tot += ca * cb * mom(aa + ab)
    return tot

def mean(A, mom):
    return sum(c * mom(at) for c, at in A)

def vk_terms(k):
    """v_k = Var(Y_{k+1}) + 2 sum_{j=1..k} Cov(Y_j, Y_{2k+2-j}).  Returns list of (a,b) pairs
    where the coefficient is 1 for (k+1,k+1) [the Var] and 2 for the crosses."""
    m = k + 1
    pairs = [(m, m, 1)]
    for j in range(1, k + 1):
        a, b = j, 2 * m - j
        if a < b:
            pairs.append((a, b, 2))
    return pairs

def compute_vk(k, n, mom, maxU=None, verbose=False):
    setup_ord = min(k + 1, CT.DENSE_MAXORD)
    CD.setup(n_grid=n, MAXORD=setup_ord)
    CT.clear_all()
    u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
    Us, s, Vs, Ye = CD.load_Y(fname=_yfile(2 * k + 1))
    need = set()
    for a, b, _ in vk_terms(k):
        need.add(a); need.add(b)
    for a, b, _ in vk_terms(k):
        pass
    Y = {}
    for idx in sorted(set([p for pr in vk_terms(k) for p in pr[:2]])):
        Y[idx] = CD.Ybase(idx, Us, s, Vs, Ye, u0ps, Vst)
    tot = 0.0
    t0 = time.time()
    for a, b, coef in vk_terms(k):
        if a == b:
            mm = mean(Y[a], mom)
            val = prod_moment(Y[a], Y[a], mom) - mm * mm
        else:
            val = prod_moment(Y[a], Y[b], mom) - mean(Y[a], mom) * mean(Y[b], mom)
        tot += coef * val
        if verbose:
            print(f"    ({a},{b}) x{coef}: {val:+.5f}  [{time.time()-t0:.0f}s]", flush=True)
    return tot

def _yfile(order):
    import os
    for avail in [7, 9, 11, 13, 15]:
        if avail >= order and os.path.exists(f'_yexprs_{avail}.txt'):
            return f'_yexprs_{avail}.txt'
    return '_yexprs_13.txt' if os.path.exists('_yexprs_13.txt') else '_yexprs_11.txt'

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else 'v2'
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    k = int(which[1:])
    eng = sys.argv[3] if len(sys.argv) > 3 else 'hybrid'
    mom = {'transfer': CT.moment_transfer, 'dense': CD.moment}.get(eng, CT.moment_hybrid)
    t0 = time.time()
    v = compute_vk(k, n, mom, verbose=True)
    print(f"[{eng}] v{k}(n={n}) = {v:+.5f}   [{time.time()-t0:.0f}s]", flush=True)
