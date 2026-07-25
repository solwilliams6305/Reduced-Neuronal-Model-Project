"""
SWALLOWTAIL (q=3) third-cumulant (skew) ladder -- the non-circular route to the Borel data.

    kappa_3(Y*) = eta^4 ( t_0 + t_1 eta^2 + t_2 eta^4 + ... ),   t_j = sum_{a+b+c = 4+2j} cum3.

Why this and not more variance rungs.  The conformal map and the direct non-perturbative extraction
both failed for the SAME circular reason: each needs the Borel data (A, theta) as input.  The skew
ladder does not.  Its late-order behaviour is controlled by the SAME Borel singularities -- same
instanton actions, same theta -- with different amplitudes, so t_0..t_J are extra INDEPENDENT
constraints on the same unknowns rather than a restatement of them.

It is also cheaper per coefficient.  t_j needs chaos index only up to 2+2j, whereas v_k needs 2k+1:

    t_5 -> idx 12      v_6 -> idx 13      v_7 -> idx 15 (10+ days)

and by the resolution rule (a grid of n points cannot resolve Y_idx for idx > n) t_j only needs
n >= 2+2j.  Every _ybase_q3_*.pkl up to idx 15 is already built, so no symbolic work is required.

Implementation.  Reuses the q=3 monkeypatches from swtl_production (the ones that passed the q=2
term-by-term gate) by importing it, rather than duplicating them -- so the backbone, the node
V-derivative table and the q-tagged Ybase cache cannot drift out of sync with the variance ladder.
Assembly follows _w_kappa3.py (the cusp's kappa_3 driver).

NOTE on multiprocessing: use 'spawn', NOT 'fork'.  numba's threading layer is not fork-safe and the
parent has already initialized it by the time the pool is created -- forking deadlocks (observed:
workers sit at 0% CPU indefinitely).  'spawn' is also correct for the patches: it re-imports THIS
module in the child, which re-runs the module-level `import swtl_production` above and so re-applies
the q=3 monkeypatches.  _init() asserts the node is the q=3 one, so a patch failure is caught loudly
instead of silently returning q=2 moments.

Usage:  python3 swtl_kappa3.py <n> [jmax]        # jmax => t_0..t_jmax
Writes swtl_kappa3_n{n}.json
"""
import os, sys, time, itertools, json

os.environ.setdefault('CT_NUMBA', '1')
import multiprocessing as mp

# --- q=3 patches: import the validated driver (Q comes from argv[1] at its import time) ---
_saved_argv = sys.argv
sys.argv = ['swtl_production', '3', '0', '12']
import swtl_production as _SP          # noqa: F401  (import applies CE.V, CD._get_sub, CD.Ybase)
sys.argv = _saved_argv

import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import _yfile

Q = 3


NODE_Q3 = -2.046707          # first Bessel-1/5 peel-off level; the q=2 node is elsewhere


def _init(n):
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    CT.clear_all()
    # guard: under 'spawn' the child must have re-applied the q=3 patches.  If it did not, the
    # backbone is the cusp's and every moment would be silently wrong.
    node = CD.GEO['s'][0]
    if abs(node - NODE_Q3) > 1e-3:
        raise RuntimeError(f"worker has the WRONG backbone: node={node:.6f}, expected "
                           f"{NODE_Q3:.6f} -- q=3 monkeypatches did not reach this process")


def _work(keys):
    return [(k, CT.moment_hybrid(list(k))) for k in keys]


def chunks(lst, m):
    out = [[] for _ in range(m)]
    for i, x in enumerate(lst):
        out[i % m].append(x)
    return [c for c in out if c]


def triples(m):
    """multisets 1<=a<=b<=c with a+b+c=m, and their symmetric multiplicity."""
    res = []
    for a in range(1, m // 3 + 1):
        for b in range(a, (m - a) // 2 + 1):
            c = m - a - b
            if c >= b:
                mu = 1 if a == b == c else (3 if (a == b or b == c) else 6)
                res.append((a, b, c, mu))
    return res


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    jmax = int(sys.argv[2]) if len(sys.argv) > 2 else 5
    nproc = max(1, (os.cpu_count() or 4) - 2)
    ms = [4 + 2 * j for j in range(jmax + 1)]
    max_idx = max(ms) - 2                      # largest single index appearing in any triple
    t0 = time.time()

    if n < max_idx:
        print(f"REFUSING: n={n} < max chaos idx {max_idx} needed for t_{jmax}. "
              f"Under-resolved grids are systematically wrong (see the resolution rule); "
              f"use n >= {max_idx} or lower jmax.", flush=True)
        return

    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
    print(f"[q={Q}] n={n} jmax={jmax}  node Yst={CD.GEO['s'][0]:.6f}  u0ps={u0ps:.4g}  "
          f"max idx {max_idx}", flush=True)
    Us, s, Vs, Ye = CD.load_Y(fname=_yfile(max_idx))

    tri_list = {m: triples(m) for m in ms}
    need = sorted({i for m in ms for (a, b, c, _) in tri_list[m] for i in (a, b, c)})
    Y = {k: CD.Ybase(k, Us, s, Vs, Ye, u0ps, Vst) for k in need}
    print(f"[{time.time()-t0:.0f}s] Y built, indices {need} "
          f"(sizes {[len(Y[k]) for k in need]})", flush=True)

    def mkey(atoms):
        return tuple(sorted(atoms))

    keyset = set()

    def add_prod(*idxs):
        for combo in itertools.product(*[Y[i] for i in idxs]):
            if all(c != 0.0 for c, _ in combo):
                atoms = []
                for _, at in combo:
                    atoms += at
                keyset.add(mkey(atoms))

    seen_pair = set()
    for m in ms:
        for a, b, c, mu in tri_list[m]:
            add_prod(a, b, c)
            for (x, y) in ((a, b), (a, c), (b, c)):
                if (x, y) not in seen_pair:
                    seen_pair.add((x, y)); add_prod(x, y)
            for x in (a, b, c):
                add_prod(x)
    keys = list(keyset)
    print(f"[{time.time()-t0:.0f}s] {len(keys)} distinct moments, nproc={nproc}", flush=True)

    M = {}
    with mp.Pool(nproc, initializer=_init, initargs=(n,)) as pool:
        for res in pool.imap_unordered(_work, chunks(keys, nproc * 8)):
            for k, val in res:
                M[k] = val
    print(f"[{time.time()-t0:.0f}s] {len(M)} moments computed", flush=True)

    def prod(*idxs):
        tot = 0.0
        for combo in itertools.product(*[Y[i] for i in idxs]):
            cf = 1.0; atoms = []; ok = True
            for c, at in combo:
                if c == 0.0:
                    ok = False; break
                cf *= c; atoms += at
            if ok:
                tot += cf * M[mkey(atoms)]
        return tot

    meanc = {a: sum(c * M[mkey(at)] for c, at in Y[a] if c != 0.0) for a in need}
    paircache = {}

    def pair(a, b):
        key = (min(a, b), max(a, b))
        if key not in paircache:
            paircache[key] = prod(a, b)
        return paircache[key]

    def cum3(a, b, c):
        return (prod(a, b, c) - meanc[a] * pair(b, c) - meanc[b] * pair(a, c)
                - meanc[c] * pair(a, b) + 2 * meanc[a] * meanc[b] * meanc[c])

    print(f"\n kappa_3 ladder (q=3) at n={n}:", flush=True)
    ladder = []
    for j, m in enumerate(ms):
        tj = sum(mu * cum3(a, b, c) for a, b, c, mu in tri_list[m])
        ladder.append(tj)
        print(f"  t_{j} (m={m}, idx<={m-2}) = {tj:+.6f}   [{time.time()-t0:.0f}s]", flush=True)

    V0 = 0.04953187          # exact q=3 one-loop variance coefficient
    print(f"\n leading skew coefficient t_0/v0^1.5 = {ladder[0]/V0**1.5:+.4f}"
          f"   (cusp analogue is ~1.16-1.20)", flush=True)
    json.dump({'q': Q, 'n': n, 'ms': ms, 't': ladder},
              open(f'swtl_kappa3_n{n}.json', 'w'), indent=2)
    print(f" wrote swtl_kappa3_n{n}.json", flush=True)


if __name__ == "__main__":
    try:
        mp.set_start_method('spawn')    # fork deadlocks with numba; spawn re-applies the patches
    except RuntimeError:
        pass
    main()
