"""
Pure-numpy (no numba, no scipy) vectorization of the pure-python dict-DP in
chaos_transfer._moment_bulk_gen.  Drop-in replacement for chaos_dp_numba.moment_bulk_numba:
reproduces BOTH DP paths (pure-chain and chains+legs) to machine precision, but vectorizes
the per-timestep transition over the whole state set with numpy instead of a Python loop.

Why this exists:  numba is unavailable in this environment (cannot be pip-installed), so the
CT_NUMBA fast path in chaos_transfer silently fell back to the Python dict-DP, whose n^{#chains}
Python-level state loop makes v7's Var(Y8) many-chain moments take minutes each (~22 h for the
grid).  This module keeps the EXACT same, validated DP semantics (state = per-chain (slots_done,
last_time) + per-leg-type remaining counts; even total participants per time; (m-1)!! w^{m/2}
weight; bond factor G[L,t] kept EXPLICIT so it is numerically stable -- unlike the rank-2
separable DP) and just executes it with array ops.

DP semantics (identical to chaos_transfer._moment_bulk_gen / chaos_dp_numba._dp_kernel):
  * sweep grid time t = 0..n-1
  * state = per-chain (c_a, L_a): c_a slots done, L_a last grid-time (-1 before first slot;
    a done chain is (korder_a, -1)); + per-leg-type remaining rem_i.
  * at each t, a subset A of not-done chains advances one slot each, and leg-type counts
    l_i <= rem_i activate; m = |A| + sum(l) must be even (m == 0 = carry).
  * factor = amp * prod_{a in A} advf(a) * dfact(m) * w[t]^(m/2) * prod_i C(rem_i,l_i) lval_i[t]^l_i,
    advf(a) = heads[a,t] if c_a==0 else Gm[L_a,t], times U[t] iff c_a+1==korder_a.
  * prune: drop any state where a not-done chain needs more slots than remaining grid times.
  * answer = amplitude of the all-done / all-legs-consumed state after the sweep.

State <-> int64 mixed-radix encoding (same as chaos_dp_numba):
  per chain a:  sc_a = c_a*(n+1) + (L_a+1),  radix Rc_a=(korder_a+1)*(n+1),  place = prod_{b<a} Rc_b
  per leg type i appended above the chain block in radix (lcount_i+1).

Returns the moment as float, or None if the mixed-radix state space would overflow int64
(caller falls back to the python DP), matching the numba wrapper's contract.
"""
import numpy as np
from itertools import product as _iproduct
import chaos_diagram as CD

# cap on the number of (newkey, contribution) pairs buffered before we compact via
# unique+bincount.  Keeps peak memory bounded on large grids (n=14 many-chain moments).
_FLUSH = 1_000_000
# hard cap on the live state count: if a moment's reachable state set exceeds this, bail out
# with None so the caller routes it elsewhere (dense) instead of OOM-ing the box.  ~5M int64+
# float keys/vals plus the bounded flush buffer stays comfortably under a few GB.  In practice
# the router only sends maxk>=6 moments here (an order>=6 chain forces heavy pruning), so the
# reachable set is small; this cap is just a safety net.
_MAXSTATES = 5_000_000
# diagnostic: peak live-state count of the most recent call (for tuning / probes)
LAST_PEAK = 0


def _dfact_int(m):
    r = 1
    k = m - 1
    while k > 1:
        r *= k
        k -= 2
    return r


def _compact(keys_list, vals_list):
    """concat buffered (keys, vals) and sum duplicates -> (ukeys, uvals)."""
    if not keys_list:
        return (np.empty(0, np.int64), np.empty(0, np.float64))
    keys = np.concatenate(keys_list)
    vals = np.concatenate(vals_list)
    if keys.size == 0:
        return (np.empty(0, np.int64), np.empty(0, np.float64))
    uk, inv = np.unique(keys, return_inverse=True)
    uv = np.zeros(uk.shape[0], np.float64)
    np.add.at(uv, inv, vals)
    return uk, uv


def moment_bulk_numpy(gatoms):
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
    korder = [int(o) for o, _ in chains]
    heads = [np.ascontiguousarray(h, dtype=np.float64) for _, h in chains]
    np1 = n + 1

    # ---- mixed-radix places (int64) ----
    Rc = [(korder[a] + 1) * np1 for a in range(NC)]
    chain_place = [0] * NC
    p = 1
    for a in range(NC):
        chain_place[a] = p
        p *= Rc[a]
    leg_place = [0] * NL
    for i in range(NL):
        leg_place[i] = p
        p *= (lcount[i] + 1)
    if p > (1 << 62):
        return None                     # overflow -> caller falls back to python DP

    # ---- start / done keys ----
    start = 0
    for i in range(NL):
        start += lcount[i] * leg_place[i]
    done = 0
    for a in range(NC):
        done += (korder[a] * np1) * chain_place[a]

    # dfact / comb tables
    maxm = NC + int(sum(lcount))
    dfact = [float(_dfact_int(m)) for m in range(maxm + 1)]
    maxleg = max(lcount) if lcount else 0
    comb2d = np.zeros((maxleg + 1, maxleg + 1), dtype=np.float64)
    for i in range(maxleg + 1):
        c = 1.0
        for j in range(i + 1):
            comb2d[i, j] = c
            c = c * (i - j) / (j + 1)

    # precompute the chain-subset list (bitmask -> tuple of chain indices)
    subsets = []
    for mask in range(1 << NC):
        A = tuple(a for a in range(NC) if (mask >> a) & 1)
        subsets.append(A)
    # precompute leg combos (tuple of l_i)
    legcombos = list(_iproduct(*[range(lcount[i] + 1) for i in range(NL)])) if NL else [()]

    # state arrays
    keys = np.array([start], dtype=np.int64)
    amps = np.array([1.0], dtype=np.float64)
    global LAST_PEAK
    LAST_PEAK = 1

    for t in range(n):
        if keys.size == 0:
            break
        if keys.size > _MAXSTATES:
            return None                     # too big -> caller routes to dense instead
        if keys.size > LAST_PEAK:
            LAST_PEAK = keys.size
        wt = w[t]
        Ut = U[t]
        remaining_t = n - t
        Gm_t = Gm[:, t]                      # bond factor to time t (indexed by last time L)
        heads_t = [heads[a][t] for a in range(NC)]
        lval_t = [lval[i][t] for i in range(NL)]

        # decode chains: c_a, L_a  (vectorized over states)
        cc = []
        LL = []
        alive = np.ones(keys.shape[0], dtype=bool)
        for a in range(NC):
            sub = (keys // chain_place[a]) % Rc[a]
            c = sub // np1
            L = (sub % np1) - 1
            cc.append(c)
            LL.append(L)
            notdone_a = c < korder[a]
            # prune: a not-done chain needing more slots than remaining times kills the state
            need = korder[a] - c
            alive &= ~(notdone_a & (need > remaining_t))
        if not alive.all():
            keys = keys[alive]
            amps = amps[alive]
            if keys.size == 0:
                break
            cc = [c[alive] for c in cc]
            LL = [L[alive] for L in LL]
        # decode legs: rem_i
        rem = []
        for i in range(NL):
            rem.append((keys // leg_place[i]) % (lcount[i] + 1))

        # per-chain advance factor and key delta (vectorized over surviving states)
        advf = []
        cdelta = []
        notdone = []
        for a in range(NC):
            c = cc[a]
            L = LL[a]
            nd = c < korder[a]
            notdone.append(nd)
            # bond/head factor; use safe index for L (clip -1 -> 0, masked by c==0 anyway)
            Lsafe = np.where(L >= 0, L, 0)
            f = np.where(c == 0, heads_t[a], Gm_t[Lsafe])
            last = (c + 1) == korder[a]
            f = np.where(last, f * Ut, f)
            advf.append(f)
            oldsc = c * np1 + (L + 1)
            newsc = np.where(last, korder[a] * np1, (c + 1) * np1 + (t + 1))
            cdelta.append((newsc - oldsc) * chain_place[a])

        buf_k = []
        buf_v = []
        buf_n = 0

        for A in subsets:
            r = len(A)
            # subset validity: every chain in A must be not-done
            if A:
                mA = notdone[A[0]].copy()
                for a in A[1:]:
                    mA &= notdone[a]
            else:
                mA = None                     # empty subset -> all states (carry handled below)
            # advance factor & key delta for this subset
            if A:
                cf = advf[A[0]].copy()
                dk = cdelta[A[0]].copy()
                for a in A[1:]:
                    cf = cf * advf[a]
                    dk = dk + cdelta[a]
            else:
                cf = None
                dk = None

            for lc in legcombos:
                # flush BEFORE growing the buffer further, so a moment with many leg-types
                # cannot accumulate #legcombos * #states rows before compacting (that was the
                # OOM: the check used to sit after the whole legcombo loop).
                if buf_n >= _FLUSH:
                    uk, uv = _compact(buf_k, buf_v)
                    if uk.size > _MAXSTATES:
                        return None
                    buf_k = [uk]; buf_v = [uv]; buf_n = uk.size
                Lsum = int(sum(lc))
                m = r + Lsum
                if m == 0:
                    # carry: same state, contribution = amp  (only the empty-subset,zero-legs term)
                    buf_k.append(keys)
                    buf_v.append(amps)
                    buf_n += keys.size
                    continue
                if m % 2 == 1:
                    continue
                # participation mask: legs available
                if A:
                    mask = mA.copy()
                else:
                    mask = np.ones(keys.shape[0], dtype=bool)
                for i in range(NL):
                    if lc[i] > 0:
                        mask &= (rem[i] >= lc[i])
                if not mask.any():
                    continue
                idx = np.nonzero(mask)[0]
                fac = amps[idx] * (dfact[m] * (wt ** (m // 2)))
                if A:
                    fac = fac * cf[idx]
                newk = keys[idx].copy()
                if A:
                    newk = newk + dk[idx]
                for i in range(NL):
                    li = lc[i]
                    if li > 0:
                        fac = fac * (comb2d[rem[i][idx], li] * (lval_t[i] ** li))
                        newk = newk - li * leg_place[i]
                # drop exact zeros to keep buffers small
                nz = fac != 0.0
                if not nz.all():
                    newk = newk[nz]
                    fac = fac[nz]
                if newk.size:
                    buf_k.append(newk)
                    buf_v.append(fac)
                    buf_n += newk.size

        keys, amps = _compact(buf_k, buf_v)
        if keys.size > _MAXSTATES:
            return None                     # post-transition blow-up -> route to dense

    # read off the done state
    if keys.size == 0:
        return 0.0
    hit = np.nonzero(keys == done)[0]
    if hit.size == 0:
        return 0.0
    return float(amps[hit[0]])
