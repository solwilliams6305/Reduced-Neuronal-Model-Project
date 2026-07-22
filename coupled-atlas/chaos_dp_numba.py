"""
Numba (nopython) reimplementation of the pure-python dict-DP in
chaos_transfer._moment_bulk_gen.  A single general kernel reproduces BOTH the
pure-chain path (no order-1 legs) and the chains+legs path, because with an
empty leg bath the general recursion collapses to the pure-chain one (odd-sized
chain subsets give odd participant counts m and are skipped, so only even |A|
survive -- exactly the pure-chain rule).

DP semantics reproduced exactly (see _moment_bulk_gen docstring / TASK spec):
  * sweep grid time t = 0..n-1
  * state = per-chain (slots_done c_a, last_time L_a)  (L = -1 initially; a chain
    is 'done' when c_a == korder_a, stored as (korder_a, -1))  +  per-leg-type
    remaining counts rem_i.
  * at each t, any subset A of not-done chains advances one slot each, and any
    leg-type counts (l_1..l_NL) with l_i <= rem_i activate; m = |A| + sum(l) must
    be even (m == 0 is the carry / do-nothing transition).
  * transition factor = amp * prod_{a in A} advfactor(a) * dfact(m) * w[t]^(m/2)
                        * prod_i C(rem_i, l_i) * lval_i[t]^l_i,
    where advfactor(a) = heads[a][t] if c_a == 0 else Gm[L_a, t], times U[t] iff
    this is the chain's LAST slot (c_a + 1 == korder_a); new L_a = t (or done).
  * prune: skip any state where a not-done chain needs more remaining slots than
    remaining grid times.
  * answer = amplitude of the all-done, all-legs-consumed state after the sweep.
  * dfact(m) = (m-1)!!.

State <-> int64 encoding (mixed radix):
  per chain a:  sc_a = c_a*(n+1) + (L_a + 1),  radix Rc_a = (korder_a+1)*(n+1)
  chains combined positionally (place[a] = prod_{b<a} Rc_b)
  per leg type i appended in its own radix (lcount_i + 1) above the chain block.

Everything float64; tiny summation-order fp differences vs the python DP are
expected and acceptable (target rel agreement ~1e-12).

This module is imported lazily/guarded by chaos_transfer so machines without
numba keep working (routing is opt-in via CT_NUMBA=1 -> NUMBA_DP).
"""
import numpy as np
import chaos_diagram as CD

from numba import njit
from numba.typed import Dict
from numba.core import types

_INT64 = types.int64
_FLOAT64 = types.float64


def _dfact_int(m):
    """(m-1)!! for even m>=0; m==0 -> 1.  (matches chaos_transfer._dfact)."""
    r = 1
    k = m - 1
    while k > 1:
        r *= k
        k -= 2
    return r


@njit(cache=True)
def _dp_kernel(n, w, U, Gm, korder, heads, lcount, lval, dfact, comb2d):
    NC = korder.shape[0]
    NL = lcount.shape[0]
    np1 = n + 1

    # ---- mixed-radix places ----
    Rc = np.empty(NC, np.int64)
    chain_place = np.empty(NC, np.int64)
    p = 1
    for a in range(NC):
        Rc[a] = (korder[a] + 1) * np1
        chain_place[a] = p
        p *= Rc[a]
    leg_place = np.empty(NL, np.int64)
    leg_radix = np.empty(NL, np.int64)
    for i in range(NL):
        leg_place[i] = p
        leg_radix[i] = lcount[i] + 1
        p *= leg_radix[i]

    # ---- start / done keys ----
    start = np.int64(0)
    for i in range(NL):
        start += lcount[i] * leg_place[i]        # rem_i = lcount_i, all chains c=0,L=-1 -> sc=0
    done = np.int64(0)
    for a in range(NC):
        done += (korder[a] * np1) * chain_place[a]   # c=korder,L=-1 -> sc=korder*np1

    states = Dict.empty(_INT64, _FLOAT64)
    states[start] = 1.0

    # ---- scratch (reused across states) ----
    cc = np.empty(NC, np.int64)
    LL = np.empty(NC, np.int64)
    rem = np.empty(NL, np.int64)
    notdone = np.empty(NC, np.int64)
    advf = np.empty(NC, np.float64)
    adv_delta = np.empty(NC, np.int64)
    legs = np.empty(NL, np.int64)

    for t in range(n):
        wt = w[t]
        Ut = U[t]
        remaining_t = n - t
        new = Dict.empty(_INT64, _FLOAT64)

        for key in states:
            amp = states[key]

            # decode chains
            nnd = 0
            prune = False
            for a in range(NC):
                sub = (key // chain_place[a]) % Rc[a]
                c = sub // np1
                L = (sub % np1) - 1
                cc[a] = c
                LL[a] = L
                if c < korder[a]:
                    if korder[a] - c > remaining_t:
                        prune = True
                        break
                    notdone[nnd] = a
                    nnd += 1
            if prune:
                continue

            # decode legs
            for i in range(NL):
                rem[i] = (key // leg_place[i]) % leg_radix[i]

            # per-not-done-chain advance factor + key delta
            for b in range(nnd):
                a = notdone[b]
                c = cc[a]
                L = LL[a]
                if c == 0:
                    f = heads[a, t]
                else:
                    f = Gm[L, t]
                nc = c + 1
                if nc == korder[a]:
                    f = f * Ut
                    newsc = korder[a] * np1                 # done: (korder,-1)
                else:
                    newsc = nc * np1 + (t + 1)              # (nc, t)
                advf[b] = f
                oldsc = c * np1 + (L + 1)
                adv_delta[b] = (newsc - oldsc) * chain_place[a]

            # number of leg combos (odometer size)
            ncombo = 1
            for i in range(NL):
                ncombo *= (rem[i] + 1)

            # enumerate chain subsets as bitmasks over the not-done list
            nmask = 1 << nnd
            for mask in range(nmask):
                r = 0
                cf = 1.0
                cdelta = np.int64(0)
                mm = mask
                b = 0
                while mm:
                    if mm & 1:
                        r += 1
                        cf *= advf[b]
                        cdelta += adv_delta[b]
                    mm >>= 1
                    b += 1
                basekey = key + cdelta

                # leg odometer
                for idx in range(ncombo):
                    q = idx
                    Lsum = 0
                    legdelta = np.int64(0)
                    for i in range(NL):
                        li = q % (rem[i] + 1)
                        q //= (rem[i] + 1)
                        legs[i] = li
                        Lsum += li
                        legdelta += li * leg_place[i]
                    m = r + Lsum
                    if m == 0:
                        # carry (mask==0, all legs zero) -> same state
                        if key in new:
                            new[key] += amp
                        else:
                            new[key] = amp
                        continue
                    if (m & 1) == 1:
                        continue
                    fac = amp * cf * dfact[m] * (wt ** (m // 2))
                    for i in range(NL):
                        li = legs[i]
                        if li:
                            fac *= comb2d[rem[i], li] * (lval[i, t] ** li)
                    if fac == 0.0:
                        continue
                    newkey = basekey - legdelta
                    if newkey in new:
                        new[newkey] += fac
                    else:
                        new[newkey] = fac

        states = new

    return states.get(done, 0.0)


def moment_bulk_numba(gatoms):
    """Thin wrapper: converts [(order, headvec)] atoms to arrays and calls the
    numba DP kernel.  Reproduces both DP paths of _moment_bulk_gen.

    Returns the moment as a float, or None if the mixed-radix state space would
    overflow int64 (caller falls back to the python DP)."""
    G = CD.GEO
    if 'G' not in G:
        raise RuntimeError("call chaos_diagram.setup(...) first")
    n = int(G['n'])
    w = np.ascontiguousarray(G['w'], dtype=np.float64)
    U = np.ascontiguousarray(G['U'], dtype=np.float64)
    Gm = np.ascontiguousarray(G['G'], dtype=np.float64)

    if not gatoms:
        return 1.0
    if sum(o for o, _ in gatoms) % 2 != 0:
        return 0.0

    chains = [(o, h) for (o, h) in gatoms if o >= 2]
    legvecs = [h for (o, h) in gatoms if o == 1]

    # group order-1 legs by head-type (byte-identical), exactly like the python DP
    tkey = {}
    ltypes = []
    lcount = []
    for h in legvecs:
        kb = h.tobytes()
        if kb not in tkey:
            tkey[kb] = len(ltypes)
            ltypes.append(h)
            lcount.append(0)
        lcount[tkey[kb]] += 1
    lval = [h * U for h in ltypes]

    NC = len(chains)
    NL = len(ltypes)
    korder = np.array([o for o, _ in chains], dtype=np.int64)
    if NC:
        heads = np.ascontiguousarray(np.array([h for _, h in chains], dtype=np.float64))
    else:
        heads = np.empty((0, n), dtype=np.float64)
    lcount_a = np.array(lcount, dtype=np.int64)
    if NL:
        lval_a = np.ascontiguousarray(np.array(lval, dtype=np.float64))
    else:
        lval_a = np.empty((0, n), dtype=np.float64)

    # overflow guard: mixed-radix state space must fit comfortably in int64
    p = 1
    for o in korder:
        p *= (int(o) + 1) * (n + 1)
    for lc in lcount:
        p *= (int(lc) + 1)
    if p > (1 << 62):
        return None

    maxm = NC + int(sum(lcount))
    dfact = np.array([_dfact_int(m) for m in range(maxm + 1)], dtype=np.float64)
    maxleg = max(lcount) if lcount else 0
    comb2d = np.zeros((maxleg + 1, maxleg + 1), dtype=np.float64)
    for i in range(maxleg + 1):
        c = 1.0
        for j in range(i + 1):
            comb2d[i, j] = c
            c = c * (i - j) / (j + 1)

    val = _dp_kernel(n, w, U, Gm, korder, heads, lcount_a, lval_a, dfact, comb2d)
    return float(val)
