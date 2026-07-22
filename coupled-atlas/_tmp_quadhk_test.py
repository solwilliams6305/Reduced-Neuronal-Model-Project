"""Self-test / benchmark for the Held-Karp subset-DP rewrite of _moment_quad.

Compares NEW _moment_quad vs:
  (a) the OLD factorial-enumeration implementation (_moment_quad_ref, inlined below), and
  (b) the generic python DP  (_moment_bulk_gen with QUAD_ENABLE=False, CT_NUMBA=0).

Cases: 2..6 chains x 0..6 legs (mixed head types), pure-chain sets, plus
7-chain and 8-chain + 4-leg cases, and an 8-chain+4-leg n=12 timing benchmark
(old enumeration vs new subset-DP).
"""
import os, time
os.environ['CT_NUMBA'] = '0'           # force pure python DP for the (b) baseline

import numpy as np
import chaos_diagram as CD
import chaos_transfer as ct
from itertools import combinations as _combinations, permutations as _permutations


# ---------------------------------------------------------------------------
# OLD implementation, inlined verbatim (the pre-rewrite factorial enumeration).
# ---------------------------------------------------------------------------
def _moment_quad_ref(gatoms):
    G = ct._geo(); n = G['n']; w = G['w']; U = G['U']; Gm = G['G']
    if sum(o for o, _ in gatoms) % 2 != 0:
        return 0.0
    chain_heads = [h for (o, h) in gatoms if o == 2]
    legvecs = [h * U for (o, h) in gatoms if o == 1]
    tkey = {}; btypes = []; bcount = []
    for b in legvecs:
        k = b.tobytes()
        if k not in tkey:
            tkey[k] = len(btypes); btypes.append(b); bcount.append(0)
        bcount[tkey[k]] += 1
    NC = len(chain_heads)
    triu = np.triu(np.ones((n, n)), 1)
    Ts = []
    for h in chain_heads:
        A = (h[:, None] * Gm * U[None, :]) * triu
        S = 0.5 * (A + A.T)
        Ts.append(S * w[None, :])
    bw = [b * w for b in btypes]

    cyc_memo = {}
    def cyc_val(seq):
        v = cyc_memo.get(seq)
        if v is None:
            M = Ts[seq[0]]
            for a in seq[1:]:
                M = M @ Ts[a]
            v = float(np.trace(M))
            cyc_memo[seq] = v
        return v

    path_memo = {}
    def path_val(u, seq, vt):
        key = (u, seq, vt)
        v = path_memo.get(key)
        if v is None:
            x = bw[u]
            for a in seq:
                x = x @ Ts[a]
            v = float(np.dot(x, btypes[vt]))
            path_memo[key] = v
        return v

    all_chains = frozenset(range(NC))
    asm_memo = {}
    def assemble(counts, chains):
        if not any(counts) and not chains:
            return 1.0
        key = (counts, chains)
        v = asm_memo.get(key)
        if v is not None:
            return v
        tot = 0.0
        rem = [c for c in chains]
        if any(counts):
            u = next(i for i, c in enumerate(counts) if c > 0)
            c1 = list(counts); c1[u] -= 1
            for vt in range(len(counts)):
                nv = c1[vt]
                if nv == 0:
                    continue
                c2 = list(c1); c2[vt] -= 1; c2t = tuple(c2)
                mult = nv
                for r in range(len(rem) + 1):
                    for sub in _combinations(rem, r):
                        rest = chains - frozenset(sub)
                        for seq in _permutations(sub):
                            tot += mult * (2.0 ** r) * path_val(u, seq, vt) * assemble(c2t, rest)
        else:
            c0 = min(chains)
            others = sorted(chains - {c0})
            for r in range(1, len(others) + 1):
                for sub in _combinations(others, r):
                    rest = chains - {c0} - frozenset(sub)
                    wgt = 2.0 ** (r + 1) if r >= 2 else 2.0
                    seen = set()
                    for perm in _permutations(sub):
                        if perm[::-1] in seen:
                            continue
                        seen.add(perm)
                        tot += wgt * cyc_val((c0,) + perm) * assemble(counts, rest)
        asm_memo[key] = tot
        return tot

    return assemble(tuple(bcount), all_chains)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def relerr(a, b):
    d = abs(a - b)
    s = max(abs(a), abs(b))
    return d / s if s > 0 else d


def make_heads(n, seed):
    """A small library of distinct head vectors so 'mixed head types' is real."""
    rng = np.random.default_rng(seed)
    lib = [
        np.ones(n),
        np.linspace(0.3, 1.7, n),
        np.cos(np.linspace(0, 3.0, n)) + 1.5,
        rng.uniform(0.5, 1.5, n),
        np.sin(np.linspace(0.2, 2.0, n)) + 1.2,
    ]
    return lib


def build_gatoms(n, n_chain, leg_spec, seed=0):
    """leg_spec: list of head-library indices for the order-1 legs."""
    lib = make_heads(n, seed)
    ga = []
    for i in range(n_chain):
        ga.append((2, lib[i % len(lib)]))
    for li in leg_spec:
        ga.append((1, lib[li % len(lib)]))
    return ga


def gen_dp(gatoms):
    """Generic python DP baseline (QUAD_ENABLE off, numba off)."""
    old = ct.QUAD_ENABLE
    ct.QUAD_ENABLE = False
    try:
        return ct._moment_bulk_gen(gatoms)
    finally:
        ct.QUAD_ENABLE = old


# ---------------------------------------------------------------------------
# test driver
# ---------------------------------------------------------------------------
def run_grid(n, do_dp=True, dp_max_total=None):
    print(f"\n===== n_grid = {n} =====", flush=True)
    CD.setup(n_grid=n, MAXORD=2)
    max_new_ref = 0.0
    max_new_dp = 0.0
    results = []

    # 2..6 chains x 0..6 legs, mixed head types
    for nc in range(2, 7):
        for nleg in range(0, 7):
            # parity: total order must be even -> 2*nc + nleg even -> nleg even
            leg_spec = [ (j + nc) for j in range(nleg) ]   # mix head types
            ga = build_gatoms(n, nc, leg_spec, seed=nc * 100 + nleg)
            total = 2 * nc + nleg
            new = ct._moment_quad(ga)
            ref = _moment_quad_ref(ga)
            e_ref = relerr(new, ref)
            max_new_ref = max(max_new_ref, e_ref)
            tag = f"nc={nc} legs={nleg}"
            dp_str = ""
            if do_dp and (dp_max_total is None or total <= dp_max_total):
                dp = gen_dp(ga)
                e_dp = relerr(new, dp)
                max_new_dp = max(max_new_dp, e_dp)
                dp_str = f"  rel(new,dp)={e_dp:.2e}"
            results.append((tag, new, e_ref))
            print(f"  {tag:16s} val={new: .6e}  rel(new,ref)={e_ref:.2e}{dp_str}", flush=True)

    # pure-chain sets (no legs), 2..6 chains -> even totals only (all even)
    print("  --- pure-chain sets ---", flush=True)
    for nc in range(2, 7):
        ga = build_gatoms(n, nc, [], seed=nc * 7)
        new = ct._moment_quad(ga)
        ref = _moment_quad_ref(ga)
        e_ref = relerr(new, ref)
        max_new_ref = max(max_new_ref, e_ref)
        dp_str = ""
        if do_dp and (dp_max_total is None or 2 * nc <= dp_max_total):
            dp = gen_dp(ga)
            e_dp = relerr(new, dp)
            max_new_dp = max(max_new_dp, e_dp)
            dp_str = f"  rel(new,dp)={e_dp:.2e}"
        print(f"  pure nc={nc:2d}       val={new: .6e}  rel(new,ref)={e_ref:.2e}{dp_str}", flush=True)

    print(f"  >>> MAX rel(new,ref) [n={n}] = {max_new_ref:.3e}")
    if do_dp:
        print(f"  >>> MAX rel(new,dp)  [n={n}] = {max_new_dp:.3e}", flush=True)
    return max_new_ref, max_new_dp


def run_big_cases(n):
    """7-chain and 8-chain + 4-leg. new vs ref (exact, slow) always;
    vs generic DP only where affordable."""
    print(f"\n===== big cases at n_grid = {n} =====", flush=True)
    CD.setup(n_grid=n, MAXORD=2)
    out = []
    for nc in (7, 8):
        leg_spec = [j + nc for j in range(4)]        # 4 mixed legs
        ga = build_gatoms(n, nc, leg_spec, seed=nc * 1000)

        t0 = time.time(); new = ct._moment_quad(ga); t_new = time.time() - t0
        t0 = time.time(); ref = _moment_quad_ref(ga); t_ref = time.time() - t0
        e_ref = relerr(new, ref)
        print(f"  nc={nc}+4legs  new={new: .6e}  ref={ref: .6e}  "
              f"rel(new,ref)={e_ref:.2e}  t_new={t_new:.3f}s t_ref={t_ref:.3f}s", flush=True)
        # generic DP for 7/8 chains is n^{nc} states -> infeasible even at n=8;
        # the exact factorial reference above is the validation baseline here
        # (it is proven equal to the DP on all smaller cases to machine eps).
        out.append((nc, e_ref))
    return out


def run_benchmark():
    """8 chains (order 2, mixed heads) + 4 legs at n=12: old vs new timing."""
    print("\n===== BENCHMARK: 8 chains (order2, mixed) + 4 legs, n_grid=12 =====", flush=True)
    CD.setup(n_grid=12, MAXORD=2)
    nc = 8
    leg_spec = [j + nc for j in range(4)]
    ga = build_gatoms(12, nc, leg_spec, seed=424242)

    t0 = time.time(); new = ct._moment_quad(ga); t_new = time.time() - t0
    t0 = time.time(); ref = _moment_quad_ref(ga); t_ref = time.time() - t0
    e = relerr(new, ref)
    print(f"  new (subset-DP)   = {new: .8e}   t = {t_new:.4f} s", flush=True)
    print(f"  old (enumeration) = {ref: .8e}   t = {t_ref:.4f} s")
    print(f"  rel(new,old)      = {e:.3e}", flush=True)
    if t_new > 0:
        print(f"  speedup           = {t_ref / t_new:.1f}x", flush=True)
    return t_new, t_ref, e


if __name__ == '__main__':
    # DP baseline is O(n^{#chains}); cap it to affordable chain counts.
    # n=8 : run DP up to total order 10 (<=5 chains, or 4 chains+2 legs, etc.)
    # n=12: DP only for small totals (<=8) so it stays in the seconds range.
    m8_ref, m8_dp = run_grid(8, do_dp=True, dp_max_total=10)
    m12_ref, m12_dp = run_grid(12, do_dp=True, dp_max_total=8)
    big8 = run_big_cases(8)
    bench = run_benchmark()

    print("\n===== SUMMARY =====", flush=True)
    print(f"n=8  : max rel(new,ref)={m8_ref:.3e}  max rel(new,dp)={m8_dp:.3e}")
    print(f"n=12 : max rel(new,ref)={m12_ref:.3e}  max rel(new,dp)={m12_dp:.3e}", flush=True)
    worst = max(m8_ref, m12_ref, m8_dp, m12_dp, *[e for _, e in big8])
    print(f"WORST rel error overall = {worst:.3e}  (target < 1e-12)", flush=True)
    print(f"benchmark n=12 8chain+4leg: t_new={bench[0]:.4f}s t_old={bench[1]:.4f}s rel={bench[2]:.2e}")
    print("PASS" if worst < 1e-12 else "FAIL", flush=True)
