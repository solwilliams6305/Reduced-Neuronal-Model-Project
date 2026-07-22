"""Multi-core v7 (and general v_k) assembly.

The assembly  v_k = sum_(a,b) coef*[<Y_a Y_b> - <Y_a><Y_b>]  reduces to a large set
of INDEPENDENT moments <sorted(atoms)> that are memoized.  We (1) enumerate every
distinct moment key, (2) compute them across worker processes, (3) assemble from the
key->value table.  Reproduces the serial _v7_n10_validate value exactly.

Usage:  python3 _v7_parallel.py <n> [nproc] [k]
        n     grid size
        nproc worker processes (default: cpu_count)
        k     which v_k (default 7)
"""
import os, sys, time, math
os.environ.setdefault('CT_NUMBA', '1')
import multiprocessing as mp
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import vk_terms, _yfile

# ---- worker: setup GEO once, then compute a chunk of moment keys ----
_WN = [None]
def _init(n, k):
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    CT.clear_all()
    _WN[0] = n

def _work(keys):
    out = []
    for key in keys:
        out.append((key, CT.moment_hybrid(list(key))))
    return out

def _chunks(lst, m):
    """split list into m roughly-equal chunks (round-robin to balance cost variance)."""
    out = [[] for _ in range(m)]
    for i, x in enumerate(lst):
        out[i % m].append(x)
    return [c for c in out if c]

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else os.cpu_count()
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    t0 = time.time()

    # ---- build Y (main process) ----
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
    Us, s, Vs, Ye = CD.load_Y(fname=_yfile(2 * k + 1))
    pairs = vk_terms(k)
    need = sorted(set(p for pr in pairs for p in pr[:2]))
    Y = {idx: CD.Ybase(idx, Us, s, Vs, Ye, u0ps, Vst) for idx in need}
    print(f"[{time.time()-t0:.0f}s] Y built, indices {need}", flush=True)

    # ---- Phase 1: enumerate distinct moment keys ----
    def mkey(atoms):
        return tuple(sorted(atoms))
    keyset = set()
    # product moments <Y_a Y_b>
    for a, b, _ in pairs:
        for ca, aa in Y[a]:
            if ca == 0.0:
                continue
            for cb, ab in Y[b]:
                if cb == 0.0:
                    continue
                keyset.add(mkey(aa + ab))
    # mean moments <Y_a>
    for idx in need:
        for c, at in Y[idx]:
            if c != 0.0:
                keyset.add(mkey(at))
    keys = list(keyset)
    print(f"[{time.time()-t0:.0f}s] Phase 1: {len(keys)} distinct moments to compute "
          f"(nproc={nproc})", flush=True)

    # ---- Phase 2: compute in parallel ----
    M = {}
    if nproc <= 1:
        _init(n, k)
        for key in keys:
            M[key] = CT.moment_hybrid(list(key))
    else:
        chunks = _chunks(keys, nproc * 6)     # 6x oversubscription for load balance
        done = 0
        with mp.Pool(nproc, initializer=_init, initargs=(n, k)) as pool:
            for res in pool.imap_unordered(_work, chunks):
                for key, val in res:
                    M[key] = val
                done += len(res)
                print(f"    [{time.time()-t0:.0f}s] {done}/{len(keys)} moments", flush=True)
    print(f"[{time.time()-t0:.0f}s] Phase 2 done ({len(M)} moments)", flush=True)

    # ---- Phase 3: assemble ----
    def prod(A, B):
        tot = 0.0
        for ca, aa in A:
            if ca == 0.0:
                continue
            for cb, ab in B:
                if cb == 0.0:
                    continue
                tot += ca * cb * M[mkey(aa + ab)]
        return tot
    def mean(A):
        return sum(c * M[mkey(at)] for c, at in A if c != 0.0)

    tot = 0.0
    for a, b, coef in pairs:
        if a == b:
            mm = mean(Y[a]); val = prod(Y[a], Y[a]) - mm * mm
        else:
            val = prod(Y[a], Y[b]) - mean(Y[a]) * mean(Y[b])
        tot += coef * val
        print(f"  ({a},{b}) x{coef}: {val:+.6f}", flush=True)

    print(f"\nv{k}(n={n}) = {tot:+.6f}   [{time.time()-t0:.0f}s total]", flush=True)
    if k == 7 and n == 10:
        anchor = -178.60840
        print(f"anchor {anchor:+.5f}  |diff|={abs(tot-anchor):.2e}  "
              f"{'MATCH' if abs(tot-anchor) < 1e-3 else 'MISMATCH'}", flush=True)
    return tot

if __name__ == "__main__":
    main()
