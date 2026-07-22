"""
Transfer-operator (grid-time sweep) Wiener-chaos moment engine.

Drop-in replacement for chaos_diagram.moment(atoms): computes <prod atoms> EXACTLY,
but never materializes a symmetric n^k tensor (the k! symmetrization wall that blocks v6).

Math.  Each 'U' atom = I_k(chain), chain_a(t_1<...<t_k)=head_a(t_1) prod G(t_p,t_{p+1}) u0(t_k),
supported on the ordered simplex.  For a product of such multiple Wiener-Ito integrals, the
diagram/Wick formula gives the moment as a sum over complete contractions (matchings) of legs
across DIFFERENT atoms; a contraction between two legs is nonzero only when they sit at the SAME
grid time t and then contributes the noise covariance w[t].  Hence:

    <prod_a I_{k_a}(chain_a)> = sum over per-atom strictly-increasing time-tuples {tau^a}
        [ prod_a chain_a(tau^a) ] * prod_t  M(occ_t) * w[t]^{occ_t/2}

where occ_t = # atoms with a slot at time t (must be even) and M(occ)=(occ-1)!! counts the
pairings at t.  Equivalently, at each time we pick a partial matching of the atoms that place a
slot there, each pair contributing adv_a(t) adv_b(t) w[t].  This is evaluated by a left-to-right
transfer over grid time; the only state is each active atom's last slot time (for its next bond),
so a self-contraction <U_k^2> costs O(n^2 k) and a general coupling stays tractable.

Boundary legs s_j = xi^(j)(node): a boundary leg pairs with one bulk leg and localizes it to the
node with j derivatives (the validated `beval` rule); boundary-boundary pairs are dropped (renorm).

Validated to machine precision against chaos_diagram.moment on all boundary-free and boundary
moments through v5 before use (see _selftest / test_transfer.py).
"""
import os
import numpy as np
from math import factorial
import chaos_diagram as CD

# Numba DP routing for _moment_bulk_gen: default ON (gated 2026-07-15: battery max
# rel 4.3e-15, ladder v3/v4 <= 1.3e-12, v6(n=10)=+80.94197 vs anchor in 390s vs ~2.1h,
# adversarial radix/overflow/boundary hunts all pass).  Escape hatch: CT_NUMBA=0.
# Import of chaos_dp_numba stays lazy/guarded (see below), so machines without a
# working numba silently keep the pure-python DP.
NUMBA_DP = os.environ.get('CT_NUMBA', '1') != '0'
_numba_bulk = None            # cached moment_bulk_numba callable (None until first use)
_numba_failed = False         # set True if the numba import/compile ever fails

def _get_numba_bulk():
    """Lazily import chaos_dp_numba.moment_bulk_numba; cache it.  Returns None
    (and disables further attempts) if numba is unavailable."""
    global _numba_bulk, _numba_failed
    if _numba_bulk is not None:
        return _numba_bulk
    if _numba_failed:
        return None
    try:
        from chaos_dp_numba import moment_bulk_numba
    except Exception:
        _numba_failed = True
        return None
    _numba_bulk = moment_bulk_numba
    return _numba_bulk

# Pure-numpy vectorized DP: the numba-free acceleration of the same dict-DP.  This is the
# workhorse when numba is unavailable (it cannot be installed in this environment).  It
# reproduces the pure-python DP to machine precision (validated: Var(Y8) maxk>=6 battery, rel
# err <= 1e-14) but vectorizes the per-timestep transition with numpy -> 20-60x faster, so
# v7's Var(Y8) many-chain moments take <1 s each instead of minutes.  Returns None if the
# reachable state set exceeds its safety cap (then we fall back to the python DP).
NUMPY_DP = os.environ.get('CT_NUMPY', '1') != '0'
_numpy_bulk = None
_numpy_failed = False

def _get_numpy_bulk():
    """Lazily import chaos_dp_numpy.moment_bulk_numpy; cache it. Returns None if unavailable."""
    global _numpy_bulk, _numpy_failed
    if _numpy_bulk is not None:
        return _numpy_bulk
    if _numpy_failed:
        return None
    try:
        from chaos_dp_numpy import moment_bulk_numpy
    except Exception:
        _numpy_failed = True
        return None
    _numpy_bulk = moment_bulk_numpy
    return _numpy_bulk

# ---------- geometry accessors (share chaos_diagram.GEO) ----------
def _geo():
    G = CD.GEO
    if 'G' not in G:
        raise RuntimeError("call chaos_diagram.setup(...) first")
    return G

def _dfact(m):  # (m-1)!! for even m>=0 ; m==0 -> 1
    r = 1
    k = m - 1
    while k > 1:
        r *= k
        k -= 2
    return r

# ---------- atom heads ----------
def _head_vec(a):
    G = _geo()
    _, k, m = a
    return G['Gnode'] if m == 0 else G['H']

# ---------- partial matchings of a list of items ----------
def _partial_matchings(items):
    """yield every set (as list) of disjoint unordered pairs drawn from items, including []."""
    if len(items) < 2:
        yield []
        return
    first = items[0]
    rest = items[1:]
    # option 1: first is unmatched at this time
    for m in _partial_matchings(rest):
        yield m
    # option 2: first pairs with some j in rest
    for idx in range(len(rest)):
        j = rest[idx]
        remain = rest[:idx] + rest[idx+1:]
        for m in _partial_matchings(remain):
            yield [(first, j)] + m

# ---------- vectorized 2-chain fast path (dominant Var terms) ----------
def _bond_sq_masked():
    G = _geo()
    key = '_Gsq_lower'
    if key not in G:
        Gm = G['G']
        n = G['n']
        Gsq = Gm * Gm
        mask = np.triu(np.ones((n, n)), 1)    # t' < t : Gsq[t',t] with row t' < col t
        G[key] = Gsq * mask
    return G[key]

def _moment_two(order, headA, headB):
    """<I_k(chainA) I_k(chainB)> both order k, forward transfer sweep O(n^2 k)."""
    G = _geo(); n = G['n']; w = G['w']; U = G['U']
    Gsq = _bond_sq_masked()          # Gsq[t', t] nonzero for t'<t
    f = headA * headB * w            # f_1[t]
    for _ in range(order - 1):
        # f_m[t] = w[t] * sum_{t'<t} f_{m-1}[t'] Gsq[t',t]
        f = w * (f @ Gsq)
    return float(np.dot(f, U * U))

# ---------- generalized bulk atom: (order, headvec) chain; tail always u0, bonds G ----------
from math import comb as _comb
from itertools import combinations as _combinations, permutations as _permutations

QUAD_ENABLE = True   # route all-order<=2 bulk moments through _moment_quad (set False to force the DP)

def _moment_quad(gatoms, extra_legs=None):
    """All atoms order<=2 (I2 quadratic + I1 linear Gaussian forms): connected-diagram expansion.

    extra_legs: optional list of pre-computed order-1 leg vectors given by their lval (= head*u0
    analogue), e.g. beval-reduced boundary legs from _moment_dense_bnd.  They join the leg bath
    directly (no head*U step).

    Moment = E[prod_a xi^T S_a xi * prod_j b_j^T xi], xi_t independent N(0,w[t]).  Connected pieces:
      cycle over chains (a1..ak), k>=2:  2^{k-1} tr(S_{a1} W S_{a2} W ... S_{ak} W)
      path  u -(a1..am)- v (legs u,v):   2^m  b_u^T W S_{a1} W ... S_{am} W b_v
    (weights verified against the transfer DP to machine precision; single-chain cycles vanish
    because the strict-simplex kernel has zero diagonal).  Assembled over all partitions by a
    counted matching recursion -- O(poly(#atoms) * n^2), replacing the exponential-state DP for
    the dominant slow class of v6 (few chains + many legs)."""
    G = _geo(); n = G['n']; w = G['w']; U = G['U']; Gm = G['G']
    if (sum(o for o, _ in gatoms) + (len(extra_legs) if extra_legs else 0)) % 2 != 0:
        return 0.0
    chain_heads = [h for (o, h) in gatoms if o == 2]
    legvecs = [h * U for (o, h) in gatoms if o == 1] + (list(extra_legs) if extra_legs else [])
    # dedupe legs by content -> counted types
    tkey = {}; btypes = []; bcount = []
    for b in legvecs:
        k = b.tobytes()
        if k not in tkey:
            tkey[k] = len(btypes); btypes.append(b); bcount.append(0)
        bcount[tkey[k]] += 1
    NC = len(chain_heads)
    # chain matrices: A[t1,t2] = head(t1) G(t1,t2) u0(t2) on t1<t2; S = (A+A^T)/2; T = S @ W
    triu = np.triu(np.ones((n, n)), 1)
    Ts = []
    for h in chain_heads:
        A = (h[:, None] * Gm * U[None, :]) * triu
        S = 0.5 * (A + A.T)
        Ts.append(S * w[None, :])          # T = S W
    bw = [b * w for b in btypes]           # b^T W

    # ------------------------------------------------------------------
    # Per-SUBSET aggregated piece values via Held-Karp subset DP over the
    # (small) chain set.  Replaces the factorial _permutations() sums in
    # assemble() with resummed aggregates -- SAME pieces, SAME weights.
    #
    #   paths:  Pagg(u, S, v) = SUM over ALL orderings seq of S of
    #                           bw[u] @ T_{seq} @ btypes[v]
    #           (weight 2^|S| still applied per piece in assemble)
    #   cycles: Cycagg(C)     = SUM over ALL orderings pi of C\{min C} of
    #                           tr(T_{min C} @ T_pi)
    #           closed contribution in assemble is 2^(|C|-1) * Cycagg(C),
    #           which EXACTLY equals the old
    #             sum_{distinct arrangement} wgt * cyc_val
    #           because (i) for |C|>=3 reflection is a fixed-point-free
    #           involution on orderings so sum_all = 2*sum_distinct and
    #           2^|C| * sum_distinct = 2^(|C|-1) * sum_all, and (ii) for
    #           |C|=2 there is a single ordering and 2^1 * it == old weight 2.
    #
    # Subsets are bitmasks over range(NC).  All aggregates depend only on the
    # chain CONTENT of a subset, not on 'rest', so they are precomputed once.
    # ------------------------------------------------------------------
    NL = len(bcount)
    full = (1 << NC) - 1

    # ---- path aggregates: for each leg-type u, psum_u[S] = length-n vector
    #      = sum over ALL orderings of S of bw[u] @ T_{seq}.  Held-Karp on
    #      xend[S][a] = same restricted to orderings ending at chain a. ----
    # index of set bits
    bits = [[a for a in range(NC) if (S >> a) & 1] for S in range(1 << NC)]
    psum = [None] * NL   # psum[u][S] -> np.ndarray(n)
    for u in range(NL):
        xend = [None] * (1 << NC)                 # xend[S] -> np.ndarray(NC, n)
        ps = [None] * (1 << NC)
        ps[0] = bw[u].copy()                      # empty interior -> bw[u]
        # singletons
        for a in range(NC):
            X = np.zeros((NC, n))
            X[a] = bw[u] @ Ts[a]
            xend[1 << a] = X
        # by increasing popcount
        order = sorted(range(1, 1 << NC), key=lambda S: bin(S).count('1'))
        for S in order:
            sb = bits[S]
            if len(sb) == 1:
                ps[S] = xend[S][sb[0]]
                continue
            X = np.zeros((NC, n))
            for a in sb:
                Sm = S ^ (1 << a)
                # sum over end-chains b of xend[Sm][b], then @ T_a
                prev = xend[Sm]
                acc = np.zeros(n)
                for b in bits[Sm]:
                    acc = acc + prev[b]
                X[a] = acc @ Ts[a]
            xend[S] = X
            ps[S] = X.sum(axis=0)
        psum[u] = ps

    def Pagg(u, S, v):
        return float(np.dot(psum[u][S], btypes[v]))

    # ---- cycle aggregates: Cycagg[C] = sum over ALL orderings of C\{c0}
    #      of tr(T_c0 @ T_ordering), c0 = min(C).  Held-Karp on matrices
    #      Qend[S][a] = sum over orderings of S ending at a of T_{seq}. ----
    Cycagg = {}
    for C in range(1, 1 << NC):
        cb = bits[C]
        if len(cb) < 2:
            continue                              # single-chain cycle vanishes
        c0 = cb[0]
        others = C ^ (1 << c0)
        ob = bits[others]
        # Held-Karp over 'others' only, Qend[S][a] = n x n matrix
        Qend = {}
        Qend[0] = None                            # empty product handled specially
        for a in ob:
            Qend[(1 << a)] = {a: Ts[a]}
        osub = sorted([S for S in range(1, 1 << NC) if (S & ~others) == 0 and S],
                      key=lambda S: bin(S).count('1'))
        for S in osub:
            sb = bits[S]
            if len(sb) == 1:
                continue                          # already seeded
            d = {}
            for a in sb:
                Sm = S ^ (1 << a)
                prevd = Qend[Sm]
                acc = None
                for b in bits[Sm]:
                    M = prevd[b]
                    acc = M if acc is None else (acc + M)
                d[a] = acc @ Ts[a]
            Qend[S] = d
        # close: tr(T_c0 @ sum_a Qend[others][a])
        prevd = Qend[others]
        acc = None
        for b in ob:
            M = prevd[b]
            acc = M if acc is None else (acc + M)
        Cycagg[C] = float(np.trace(Ts[c0] @ acc))

    all_chains = frozenset(range(NC))
    def _mask(chains):
        m = 0
        for c in chains:
            m |= (1 << c)
        return m

    asm_memo = {}
    def assemble(counts, chains):
        if not any(counts) and not chains:
            return 1.0
        key = (counts, chains)
        v = asm_memo.get(key)
        if v is not None:
            return v
        tot = 0.0
        cmask = _mask(chains)
        if any(counts):
            u = next(i for i, c in enumerate(counts) if c > 0)
            c1 = list(counts); c1[u] -= 1
            for vt in range(len(counts)):
                nv = c1[vt]
                if nv == 0:
                    continue
                c2 = list(c1); c2[vt] -= 1; c2t = tuple(c2)
                mult = nv
                # ordered chain sequence between the two legs, resummed per subset
                sub = cmask
                while True:                       # iterate all submasks of cmask
                    rest = frozenset(chains - {c for c in range(NC) if (sub >> c) & 1}) \
                        if sub else chains
                    r = bin(sub).count('1')
                    val = Pagg(u, sub, vt)
                    if val != 0.0:
                        tot += mult * (2.0 ** r) * val * assemble(c2t, rest)
                    if sub == 0:
                        break
                    sub = (sub - 1) & cmask
        else:
            c0 = min(chains)
            c0bit = 1 << c0
            othersmask = cmask ^ c0bit
            # cycle set C = {c0} union (nonempty submask of others)
            sub = othersmask
            while True:
                if sub:                           # cycle size >= 2
                    C = c0bit | sub
                    rest = frozenset(chains - {c0} -
                                     {c for c in range(NC) if (sub >> c) & 1})
                    ca = Cycagg.get(C)
                    if ca:
                        r = bin(sub).count('1')   # cycle size = r+1
                        tot += (2.0 ** (r + 1 - 1)) * ca * assemble(counts, rest)
                if sub == 0:
                    break
                sub = (sub - 1) & othersmask
        asm_memo[key] = tot
        return tot

    return assemble(tuple(bcount), all_chains)

def _chain_tensor(order, head):
    """Weight-absorbed SYMMETRIC order-tensor of a chain with the given head vector
    (mirrors chaos_engine.chain / chaos_diagram.atom_kernel, but for an ARBITRARY head, so it
    also serves boundary-REDUCED chains whose head is node-pinned by reduce_head).  order<=5
    only (the symmetrization is order! transposes of an n^order tensor)."""
    G = _geo(); Gm = G['G']; U = G['U']; sw = G['sw']; n = G['n']
    if order == 1:
        return (head * U) * sw                      # 1-tensor, weights absorbed
    triu = np.triu(np.ones((n, n)), 1)
    Gmask = Gm * triu
    P = np.asarray(head, dtype=float).copy()
    for _ in range(2, order + 1):
        P = P[..., None] * Gmask                     # extend ordered chain: bake i<j bonds
    T = P * U                                        # tail u0 on last axis
    s = np.zeros_like(T)
    for perm in _permutations(range(order)):
        s += np.transpose(T, perm)
    T = s / factorial(order)                         # symmetrize
    wt = sw.copy()
    for _ in range(order - 1):
        wt = np.multiply.outer(wt, sw)               # absorb sw^{⊗order}
    return T * wt

_DENSE_LETTERS = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

def _tn_tensors(tensors, M):
    """Contract explicit weight-absorbed tensors per multigraph M (bulk-bulk edges only) --
    the boundary-free core of chaos_diagram._tn on an arbitrary tensor list."""
    N = len(tensors)
    letters = iter(_DENSE_LETTERS)
    leg_labels = {i: [] for i in range(N)}
    for (i, j), m in M.items():
        for _ in range(m):
            L = next(letters); leg_labels[i].append(L); leg_labels[j].append(L)
    ops = []; subs = []
    for i in range(N):
        T = tensors[i]
        if T.ndim != len(leg_labels[i]):
            return 0.0
        ops.append(T); subs.append(''.join(leg_labels[i]))
    if not ops:
        return 1.0
    if len(ops) == 1:
        return float(ops[0].sum() if subs[0] == '' else np.einsum(subs[0] + '->', ops[0]))
    return float(np.einsum(','.join(subs) + '->', *ops, optimize='greedy'))

def _moment_bulk_dense(gatoms):
    """Dense multigraph contraction on (order, head) chains -- the SAME algorithm as
    chaos_diagram.moment (validated), but on the transfer engine's explicit-head representation
    so it also handles boundary-REDUCED chains.  Fast for MANY-CHAIN, LOW-order moments (max
    order <= DENSE_MAXORD), the exact regime where the DP's n^{#chains} state set explodes and
    the numpy DP would bail.  Never materializes a tensor above the largest single chain."""
    orders = [o for o, _ in gatoms]
    if sum(orders) % 2 != 0:
        return 0.0
    tensors = [_chain_tensor(o, h) for o, h in gatoms]
    N = len(gatoms)
    tot = 0.0
    for M in CD._enum_multigraphs(orders, [False] * N):
        mult = 1.0
        for o in orders:
            mult *= factorial(o)
        for e, m in M.items():
            mult /= factorial(m)
        if mult != 0.0:
            tot += mult * _tn_tensors(tensors, M)
    return tot

def _moment_bulk_gen(gatoms):
    """gatoms: list of (order:int>=1, headvec:np.ndarray len n).  chain = head(t1) prod G u0(t_last).

    Transfer DP over grid time.  At each time the contractions there pair up the participants
    (advancing chain-slots + activated order-1 legs); every perfect matching of m participants gives
    the same value, so we weight by (m-1)!! instead of enumerating matchings.  Order-1 atoms are a
    counted 'leg bath' grouped by head-type, so <chain . leg^p> costs O(n^2 * k * p), not 2^p."""
    G = _geo()
    n = G['n']; w = G['w']; U = G['U']; Gm = G['G']
    if not gatoms:
        return 1.0
    if sum(o for o, _ in gatoms) % 2 != 0:
        return 0.0
    if QUAD_ENABLE and all(o <= 2 for o, _ in gatoms):
        return _moment_quad(gatoms)      # Gaussian quadratic/linear-form fast path (validated vs DP)
    # many-chain, low-order: the DP's n^{#chains} state set explodes (7-chain (2,2,2,2,2,3,3):
    # transfer >25 s vs dense 0.3 s), so use the dense multigraph contraction instead.  Applies
    # here (not just in moment_hybrid) so it ALSO catches boundary-REDUCED moments coming through
    # _mt_recurse, whose modified heads dense-proper (CD.moment) cannot represent.
    _maxk = max((o for o, _ in gatoms), default=0)
    _nchain = sum(1 for o, _ in gatoms if o >= 2)
    if DENSE_ROUTE and 3 <= _maxk <= DENSE_MAXORD:
        _nlegtypes = len(set(h.tobytes() for o, h in gatoms if o == 1))
        if _nchain >= DENSE_NCHAIN or _nlegtypes >= DENSE_LEGTYPES:
            return _moment_bulk_dense(gatoms)
    chains = [(o, h) for (o, h) in gatoms if o >= 2]
    legvecs = [h for (o, h) in gatoms if o == 1]
    # group order-1 legs by head-type
    tkey = {}; ltypes = []; lcount = []
    for h in legvecs:
        k = h.tobytes()
        if k not in tkey:
            tkey[k] = len(ltypes); ltypes.append(h); lcount.append(0)
        lcount[tkey[k]] += 1
    lval = [h * U for h in ltypes]           # leg factor per type: h(t) u0(t)
    NL = len(ltypes)
    lcount0 = tuple(lcount)

    # fast path: exactly two chains, no legs
    if len(chains) == 2 and NL == 0:
        (o0, h0), (o1, h1) = chains
        return _moment_two(o0, h0, h1) if o0 == o1 else 0.0
    # opt-in numba DP: reproduces BOTH the pure-chain and chains+legs python DP
    # below.  Guarded/lazy import; returns None on state-space overflow (falls
    # through to the python DP).
    if NUMBA_DP and not _numba_failed:
        fn = _get_numba_bulk()
        if fn is not None:
            _r = fn(gatoms)
            if _r is not None:
                return _r
    # pure-numpy vectorized DP (numba-free): same DP, vectorized with numpy.  This is the
    # fast path when numba is absent.  Returns None if the state set exceeds its cap, in
    # which case we fall through to the pure-python DP (correct but slow).
    if NUMPY_DP and not _numpy_failed:
        fn = _get_numpy_bulk()
        if fn is not None:
            _r = fn(gatoms)
            if _r is not None:
                return _r
    if SEP_ENABLE:
        return _moment_sep(gatoms)       # rank-2 separable DP (n-independent state; validated vs DP)
    if not chains:
        # legs only: even total handled below via DP with no chains
        pass

    korder = [o for o, _ in chains]
    heads = [h for _, h in chains]
    NC = len(chains)
    korder_t = tuple(korder)

    # ---- pure-chain fast path (no order-1 legs): only subsets of chains advance, |A| even ----
    if NL == 0:
        start = tuple((0, -1) for _ in range(NC))
        states = {start: 1.0}
        rng = range(NC)
        for t in range(n):
            wt = w[t]; Gcol = Gm[:, t]; Ut = U[t]
            remaining_t = n - t
            new = {}
            new_get = new.get
            for phase, amp in states.items():
                notdone = [a for a in rng if phase[a][0] < korder_t[a]]
                if any(korder_t[a] - phase[a][0] > remaining_t for a in notdone):
                    continue
                # carry (empty A)
                new[phase] = new_get(phase, 0.0) + amp
                # per-chain advance factor at this t
                advf = {}
                for a in notdone:
                    c, Lst = phase[a]
                    f = heads[a][t] if c == 0 else Gcol[Lst]
                    if c + 1 == korder_t[a]:
                        f *= Ut
                    advf[a] = f
                nd = len(notdone)
                for r in range(2, nd + 1, 2):        # |A| must be even
                    dfac_w = _dfact(r) * (wt ** (r // 2))
                    for A in _combinations(notdone, r):
                        cf = amp * dfac_w
                        for a in A:
                            cf *= advf[a]
                        if cf == 0.0:
                            continue
                        newphase = list(phase)
                        for a in A:
                            c = phase[a][0] + 1
                            newphase[a] = (korder_t[a], -1) if c == korder_t[a] else (c, t)
                        key = tuple(newphase)
                        new[key] = new_get(key, 0.0) + cf
                        new_get = new.get
            states = new
        done = tuple((korder_t[a], -1) for a in range(NC))
        return float(states.get(done, 0.0))

    # ---- general path: chains + order-1 leg bath ----
    def adv_factor(a, phase, t):
        c, L = phase
        f = heads[a][t] if c == 0 else Gm[L, t]
        if c + 1 == korder[a]:
            f = f * U[t]
        return f

    _lc_cache = {}
    def leg_choices(rem):
        opts = _lc_cache.get(rem)
        if opts is None:
            opts = list(_iproduct(*[range(r + 1) for r in rem]))
            _lc_cache[rem] = opts
        return opts

    start = (tuple((0, -1) for _ in range(NC)), lcount0)
    states = {start: 1.0}
    for t in range(n):
        wt = w[t]
        remaining_t = n - t
        lvt = [lv[t] for lv in lval]
        new = {}
        for (phase, rem), amp in states.items():
            notdone = [a for a in range(NC) if phase[a][0] < korder[a]]
            if any(korder[a] - phase[a][0] > remaining_t for a in notdone):
                continue
            for r in range(len(notdone) + 1):
                for A in _combinations(notdone, r):
                    cf = 1.0
                    for a in A:
                        cf *= adv_factor(a, phase[a], t)
                    for legs in leg_choices(rem):
                        L = sum(legs)
                        m = len(A) + L
                        if m == 0:
                            new[(phase, rem)] = new.get((phase, rem), 0.0) + amp
                            continue
                        if m % 2 != 0:
                            continue
                        fac = amp * cf * _dfact(m) * (wt ** (m // 2))
                        for ti in range(NL):
                            a_tau = legs[ti]
                            if a_tau:
                                fac *= _comb(rem[ti], a_tau) * (lvt[ti] ** a_tau)
                        if fac == 0.0:
                            continue
                        newphase = list(phase)
                        for a in A:
                            c, _2 = phase[a]; nc = c + 1
                            newphase[a] = (korder[a], -1) if nc == korder[a] else (nc, t)
                        newrem = tuple(rem[ti] - legs[ti] for ti in range(NL))
                        key = (tuple(newphase), newrem)
                        new[key] = new.get(key, 0.0) + fac
        states = new
    done = (tuple((korder[a], -1) for a in range(NC)), tuple(0 for _ in range(NL)))
    return float(states.get(done, 0.0))

SEP_ENABLE = False  # rank-2 separable DP: numerically DEAD in double precision -- U spans ~16 orders
                    # of magnitude over the interval, so the branch products U*PSI exceed the true bond
                    # G by ~1e4 and the per-bond cancellation compounds (measured: rel error 1e5 even for
                    # single-bond chains, 1e20 for the slow class).  Kept for the record; do not enable.

def _moment_sep(gatoms):
    """Rank-2 separable transfer DP.  G(t,t') = U(t)PSI(t') - PSI(t)U(t') is rank 2, so expanding
    every bond makes each chain kernel a sum of PER-SLOT-FACTOR products on the ordered simplex:
    the DP state per chain is (slots consumed, pending B-factor eps) -- 2k states, NOT the last grid
    time.  Total state space is n-independent (vs n^{#chains} for the generic sweep), which is what
    makes >=3-chain moments and large grids affordable.  Wick combinatorics at each time step is
    identical to the generic DP ((m-1)!! w^{m/2} for m participants).  Validated vs the DP."""
    G = _geo()
    n = G['n']; w = G['w']; U = G['U']; PSI = G['PSI']
    if sum(o for o, _ in gatoms) % 2 != 0:
        return 0.0
    chains = [(o, h) for (o, h) in gatoms if o >= 2]
    legvecs = [h for (o, h) in gatoms if o == 1]
    tkey = {}; ltypes = []; lcount = []
    for h in legvecs:
        kb = h.tobytes()
        if kb not in tkey:
            tkey[kb] = len(ltypes); ltypes.append(h); lcount.append(0)
        lcount[tkey[kb]] += 1
    lval = [h * U for h in ltypes]
    NL = len(ltypes); lcount0 = tuple(lcount)
    NC = len(chains)
    korder = [o for o, _ in chains]
    A2 = (U, PSI); B2 = (PSI, U); SGN = (1.0, -1.0)   # G(t,t') = sum_eps SGN A2[eps](t) B2[eps](t')

    # per-chain advance options: opts[a][(c,pend)] = [(vec over t, newstate), ...]; newstate None=done
    opts = []
    for a in range(NC):
        k = korder[a]; h = chains[a][1]; d = {}
        if k == 1:
            raise ValueError("order-1 atoms belong to the leg bath")
        # first slot (c=0): head * A_eps (sign on A side); k==2 -> next is last
        d[(0, -1)] = [(h * A2[e] * SGN[e], (1, e)) for e in range(2)]
        for c in range(1, k - 1):
            for p in range(2):
                d[(c, p)] = [(B2[p] * A2[e] * SGN[e], (c + 1, e)) for e in range(2)]
        for p in range(2):
            d[(k - 1, p)] = [(B2[p] * U, None)]
        opts.append(d)

    DONE = (10**9, 0)
    _lc_cache = {}
    def leg_choices(rem):
        v = _lc_cache.get(rem)
        if v is None:
            v = list(_iproduct(*[range(r + 1) for r in rem]))
            _lc_cache[rem] = v
        return v

    start = (tuple((0, -1) for _ in range(NC)), lcount0)
    states = {start: 1.0}
    for t in range(n):
        wt = w[t]
        remaining_t = n - t
        lvt = [lv[t] for lv in lval]
        new = {}
        new_get = new.get
        for (phase, rem), amp in states.items():
            notdone = [a for a in range(NC) if phase[a] != DONE]
            if any(korder[a] - phase[a][0] > remaining_t for a in notdone):
                continue
            for r in range(len(notdone) + 1):
                for A in _combinations(notdone, r):
                    # branch over each advancing chain's options
                    for branch in _iproduct(*[opts[a][phase[a]] for a in A]):
                        cf = 1.0
                        for (vec, _ns) in branch:
                            cf *= vec[t]
                        if cf == 0.0 and r > 0:
                            continue
                        for legs in leg_choices(rem):
                            L = sum(legs)
                            m = r + L
                            if m == 0:
                                key = (phase, rem)
                                new[key] = new_get(key, 0.0) + amp
                                continue
                            if m % 2 != 0:
                                continue
                            fac = amp * cf * _dfact(m) * (wt ** (m // 2))
                            for ti in range(NL):
                                a_tau = legs[ti]
                                if a_tau:
                                    fac *= _comb(rem[ti], a_tau) * (lvt[ti] ** a_tau)
                            if fac == 0.0:
                                continue
                            newphase = list(phase)
                            for idx, a in enumerate(A):
                                ns = branch[idx][1]
                                newphase[a] = DONE if ns is None else ns
                            newrem = tuple(rem[ti] - legs[ti] for ti in range(NL))
                            key = (tuple(newphase), newrem)
                            new[key] = new_get(key, 0.0) + fac
        states = new
    done = (tuple(DONE for _ in range(NC)), tuple(0 for _ in range(NL)))
    return float(states.get(done, 0.0))

# ---------- boundary reduction: pin a chain's first slot to the node with d^j, x bfac ----------
def _dj_at_node(vec_x_times_something):
    """apply d^j/dx^j at node (grid index 0) along axis 0; returns the axis-0=0 slice."""
    raise NotImplementedError

def _reduce_head(head, j):
    """A boundary leg s_j pins the FIRST slot (at node) of a chain with this head.
    Returns the new head vector over the NEXT slot tau:  bfac * d^j_x[ head(x) G(x,tau) ]|_{x=node}.
    (head(x)=G(node,x) for U*0, or dG/dY(node,x) for U*1; the pinned factor is head(x)*G(x,tau)."""
    G = _geo(); h = G['h']; Gm = G['G']; bfac = G.get('bfac', 0.5)
    M = head[:, None] * Gm            # M[x,tau] = head(x) G(x,tau)
    for _ in range(j):
        M = np.gradient(M, h, axis=0)
    return bfac * M[0]                # evaluate x=node (index 0) -> vector over tau

def _reduce_scalar(head, j):
    """Order-1 chain (head(t) u0(t)) with its single slot pinned to node by s_j -> scalar."""
    G = _geo(); h = G['h']; U = G['U']; bfac = G.get('bfac', 0.5)
    v = head * U                      # head(x) u0(x)
    for _ in range(j):
        v = np.gradient(v, h)
    return bfac * v[0]

# ---------- FAST dense (beval) boundary path for all-order<=2-chain moments ----------
# The certified v0..v6 anchors use the dense/beval boundary convention.  CD.moment computes it
# correctly but its multigraph enumeration blows up when a boundary moment has many order-1 legs
# (the slow class).  For moments whose chains are all order<=2, we instead reduce each boundary leg
# with the validated CD.beval (numerically, on a tiny order-2 tensor) and contract the remaining
# order<=2 bulk with the fast quad Gaussian path -- machine-precision match to CD.moment (validated
# in _test_fastbnd), but quad-fast.  Moments containing an order>=3 chain fall back to CD.moment
# (beval on them yields a tensor, not a leg; and they carry fewer legs so CD.moment is affordable).
def _chain2_tensor(head):
    """weight-absorbed symmetric order-2 tensor for a chain with this head (matches atom_kernel)."""
    G = _geo(); Gm = G['G']; U = G['U']; sw = G['sw']; n = G['n']
    triu = np.triu(np.ones((n, n)), 1)
    A = (head[:, None] * Gm * U[None, :]) * triu
    T = 0.5 * (A + A.T)
    return T * sw[:, None] * sw[None, :]

def _beval_leg2(head, j):
    """order-2 chain reduced by boundary s_j (dense/beval) -> order-1 'beval-leg' lval."""
    G = _geo(); sw = G['sw']
    vA = CD.beval(_chain2_tensor(head), 0, j)     # weight-absorbed vector over the free slot
    return vA / sw                                 # un-absorb -> leg lval

def _reduce_scalar_lval(lval, j):
    """boundary s_j reduces an order-1 leg given by its lval (= head*u0): bfac * d^j[lval]|node."""
    G = _geo(); h = G['h']; bfac = G.get('bfac', 0.5)
    v = lval.copy()
    for _ in range(j):
        v = np.gradient(v, h)
    return bfac * v[0]

def _dense_bnd_rec(chains, leglvals, bnds):
    """recursive beval boundary reduction (dense convention); base case -> quad."""
    if not bnds:
        return _moment_quad([(2, h) for h in chains], extra_legs=leglvals)
    j = bnds[0]; rest = bnds[1:]; tot = 0.0
    for i in range(len(leglvals)):
        scal = _reduce_scalar_lval(leglvals[i], j)
        if scal != 0.0:
            tot += scal * _dense_bnd_rec(chains, leglvals[:i] + leglvals[i + 1:], rest)
    for c in range(len(chains)):
        bl = _beval_leg2(chains[c], j)             # order-2 chain -> beval-leg; factor 2 = order!
        tot += 2.0 * _dense_bnd_rec(chains[:c] + chains[c + 1:], leglvals + [bl], rest)
    return tot

def _moment_dense_bnd(atoms):
    """Fast dense (beval) boundary moment for all-order<=2-chain moments; CD.moment fallback else."""
    if any(a[0] == 'U' and a[1] >= 3 for a in atoms):
        return CD.moment(atoms)                    # order>=3 chain -> beval yields a tensor
    G = _geo(); U = G['U']
    chains = [_head_vec(a) for a in atoms if a[0] == 'U' and a[1] == 2]
    leglvals = [_head_vec(a) * U for a in atoms if a[0] == 'U' and a[1] == 1]
    bnds = [a[1] for a in atoms if a[0] == 's']
    return _dense_bnd_rec(chains, leglvals, bnds)

# ---------- top-level moment (bulk + boundary), matches chaos_diagram.moment ----------
from itertools import product as _iproduct

_TMEMO = {}
def clear():
    _TMEMO.clear()
    _RMEMO.clear()

# Hybrid: dense engine is fast+validated for max bulk order <= DENSE_MAXORD; the transfer sweep
# is used only for genuinely high-order moments (>= 6-chaos) that the dense n^k symmetrization
# cannot reach.  Both are memoized here; setup with MAXORD=DENSE_MAXORD (cheap) suffices because
# the transfer sweep only needs G/Gnode/H/U (all O(n^2), always built).
DENSE_MAXORD = 5
# pure-bulk router threshold: a moment with >= DENSE_NCHAIN order>=2 chains and max order in
# [3, DENSE_MAXORD] is sent to the dense engine (fast for many low-order chains) instead of the
# transfer DP (whose n^{#chains} state set explodes there).  5 is safe: the numpy DP handles
# <=4-chain moments (incl. many-leg ones) comfortably, and dense handles >=5-chain low-order ones.
DENSE_NCHAIN = 5          # >=5 order>=2 chains at low order -> dense: only here does the numpy DP's
                          # n^{#chains} state set genuinely explode (7-chain (2,2,2,2,2,3,3): dense
                          # 0.3 s).  nchain<=4 stays on the numpy DP: its COUNTED leg-bath is far
                          # faster than dense's multigraph enum when a moment has many order-1 legs
                          # (e.g. (1^7,2,2,2,3): dense 8 s vs numpy DP <0.5 s).
DENSE_LEGTYPES = 9        # >= this many DISTINCT order-1 leg heads -> dense (safety net only): the
                          # numpy DP leg odometer is prod(count_i+1), fine for the <=8 distinct-head
                          # boundary-reduced moments; dense is a fallback for pathological cases.
DENSE_ROUTE = True        # route many-chain low-order bulk moments to the dense contractor (set
                          # False to force the DP everywhere -- used only for cross-validation)
BND_ENGINE = 'transfer'   # convention for boundary moments: 'transfer' (uniform, all orders) or 'dense'
_HMEMO = {}

def _hafnian(C, idxs, memo):
    """hafnian of the covariance submatrix on index-tuple idxs (sorted). C: dict (i,j)->cov, i<j."""
    if not idxs:
        return 1.0
    v = memo.get(idxs)
    if v is not None:
        return v
    a = idxs[0]; rest = idxs[1:]; s = 0.0
    for k in range(len(rest)):
        b = rest[k]
        c = C[(a, b)]
        if c != 0.0:
            s += c * _hafnian(C, rest[:k] + rest[k+1:], memo)
    memo[idxs] = s
    return s

def _degree1_moment(atoms):
    """all atoms are degree-1 (order-1 'U' or boundary 's'): moment = hafnian of pairwise covariances
    (boundary-boundary pair = 0 by renorm).  Fast (O(2^m) not O((m-1)!!) with per-matching work)."""
    m = len(atoms)
    if m % 2 != 0:
        return 0.0
    C = {}
    for i in range(m):
        for j in range(i + 1, m):
            C[(i, j)] = CD.moment([atoms[i], atoms[j]])
    return _hafnian(C, tuple(range(m)), {})

def moment_hybrid(atoms):
    if not atoms:
        return 1.0
    key = tuple(sorted(atoms))
    v = _HMEMO.get(key)
    if v is not None:
        return v
    if _connectivity_zero(atoms):
        v = 0.0
    elif all(a[0] == 's' or a[1] == 1 for a in atoms):
        v = _degree1_moment(atoms)           # hafnian, dense-convention pairwise covariances
    else:
        maxk = max((a[1] for a in atoms if a[0] == 'U'), default=0)
        has_bnd = any(a[0] == 's' for a in atoms)
        if has_bnd:
            # Boundary moments are individually grid-divergent (delta(0) self-contractions); they are
            # finite only after cancellation in the assembled v_k, which needs ONE consistent
            # convention.  BND_ENGINE selects it: 'transfer' (reduce_head recursive Wick -- one
            # convention for ALL orders, incl. >5 which dense cannot reach) or 'dense' (beval; the
            # convention the certified v0..v6 anchors were built with, capped at 5-chaos).
            #
            # FAST dense-boundary path: reduce_head (transfer) and beval (dense) are IDENTICAL for
            # boundary legs s_j with j<=1 (verified to ~1e-16); they differ only for j>=2 (a
            # discretization branch).  So in 'dense' mode we only pay the slow CD.moment multigraph
            # for moments that actually contain a j>=2 leg; all-j<=1 boundary moments go through the
            # fast transfer engine and give the identical dense value.  This cuts the slow dense set
            # by ~2/3 (measured on v6: 4316 maxk<=5 boundary moments -> only 1358 need CD.moment).
            if BND_ENGINE == 'dense' and maxk <= DENSE_MAXORD:
                maxj = max(a[1] for a in atoms if a[0] == 's')
                if maxj >= 2:
                    v = _moment_dense_bnd(atoms)         # beval convention; fast for order<=2 chains,
                                                         # CD.moment fallback for order>=3 chains
                else:
                    v = moment_transfer(atoms)           # j<=1: transfer == dense, but FAST
            else:
                v = moment_transfer(atoms)               # order>5 boundary: transfer (dense can't build)
        else:
            # pure bulk: transfer == dense to machine precision (no boundary convention), so we route
            # each moment to whichever engine is fast for its SHAPE.  The two engines are complementary
            # (measured on Var(Y8)):
            #   * dense (CD.moment) builds tensors only up to the max single-atom order, so it is fast
            #     for MANY-CHAIN, LOW-order moments -- e.g. 7-chain (2,2,2,2,2,3,3): dense 0.3 s, while
            #     the transfer DP's n^{#chains} state set explodes (>25 s).  But dense cannot build an
            #     order>=6 atom (setup MAXORD=5) and its multigraph enum blows up with many order-1 legs.
            #   * transfer's _moment_bulk_gen routes internally to the quad Gaussian path (all order<=2,
            #     kills the many-leg class), the 2-chain sweep, and the numpy vectorized DP (order>=3;
            #     an order>=6 chain forces <=6 chains + heavy pruning, so the DP stays small and fast).
            # Rule: many low-order chains -> dense; everything else -> transfer.  (The old 2026-07-17
            # "route everything to transfer" fixed a dense many-leg hang but REINTRODUCED the many-chain
            # explosion that stalled v7 -- this shape-aware split fixes both.)
            maxk = max((a[1] for a in atoms if a[0] == 'U'), default=0)
            nchain = sum(1 for a in atoms if a[0] == 'U' and a[1] >= 2)
            if 3 <= maxk <= DENSE_MAXORD and nchain >= DENSE_NCHAIN:
                v = CD.moment(atoms)         # many-chain, low-order: dense is fast; transfer explodes
            else:
                v = moment_transfer(atoms)   # quad / 2-chain / numpy-DP
    _HMEMO[key] = v
    return v

def clear_all():
    _TMEMO.clear(); _HMEMO.clear(); _RMEMO.clear(); CD.clear_cache()

def moment_transfer(atoms):
    if not atoms:
        return 1.0
    key = tuple(sorted(atoms))
    v = _TMEMO.get(key)
    if v is not None:
        return v
    v = _moment_transfer_impl(atoms)
    _TMEMO[key] = v
    return v

def _connectivity_zero(atoms):
    """True if some atom must have an unpaired leg (moment = 0): a bulk atom of order k with
    k > (sum of all other bulk orders) + (# boundary legs).  Boundary legs pair only with bulk."""
    bulkorders = [a[1] for a in atoms if a[0] == 'U']
    nb = sum(1 for a in atoms if a[0] == 's')
    if not bulkorders:
        return nb > 0            # only boundary legs -> can't pair (bnd-bnd dropped) -> 0
    tot = sum(bulkorders)
    for k in bulkorders:
        if k > (tot - k) + nb:
            return True
    if (tot + nb) % 2 != 0:
        return True
    return False

def _moment_transfer_impl(atoms):
    if _connectivity_zero(atoms):
        return 0.0
    # generalized items: ('s', j) boundary, or (order, headvec) bulk chain
    items = []
    for a in atoms:
        items.append(('s', a[1]) if a[0] == 's' else (a[1], _head_vec(a)))
    return _mt_recurse(items)

def _items_key(items):
    parts = []
    for it in items:
        if it[0] == 's':
            parts.append(('s', it[1], b''))
        else:
            parts.append(('u', it[0], it[1].tobytes()))
    return tuple(sorted(parts))

_RMEMO = {}
def _mt_recurse(items):
    """Recursive Wick over boundary legs: pull one boundary s_j, sum over the bulk leg it contracts
    with (contracting two boundary legs is dropped -- the renormalization), recurse.  Pulling one at
    a time (not reducing a chain by several boundaries at once) is what makes multiple distinct
    boundary legs correct.  When no boundary remains, contract the pure-bulk chains via the sweep.
    Memoized (branches reach identical sub-item-sets) -- essential for multi-boundary moments."""
    key = _items_key(items)
    cached = _RMEMO.get(key)
    if cached is not None:
        return cached
    bidx = next((i for i, it in enumerate(items) if it[0] == 's'), None)
    if bidx is None:
        v = _moment_bulk_gen(items)
        _RMEMO[key] = v
        return v
    j = items[bidx][1]
    rest = items[:bidx] + items[bidx + 1:]
    tot = 0.0
    for i, it in enumerate(rest):
        if it[0] == 's':
            continue                                  # drop boundary-boundary (renorm)
        order, head = it
        others = rest[:i] + rest[i + 1:]
        if order == 1:
            scal = _reduce_scalar(head, j)
            if scal != 0.0:
                tot += scal * _mt_recurse(others)
        else:
            newhead = _reduce_head(head, j)
            tot += _mt_recurse(others + [(order - 1, newhead)])
    _RMEMO[key] = tot
    return tot

