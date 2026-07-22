"""
Hybrid HARD/SOFT Wiener-chaos moment engine (REFERENCE, clarity over speed).

Computes  M = < prod_a I_{k_a}(chain_a) . prod_j leg_j >  for one Gaussian grid
field xi_t ~ N(0, w[t]) independent, EXACTLY, with DP state ~ n^(#order>=3 chains)
instead of n^(#chains).

-------------------------------------------------------------------------------
DERIVATION (exact hybrid factorization)
-------------------------------------------------------------------------------
Atoms split into three sectors (all in the SAME field xi):

  HARD  : chains of order k>=3.  Chain a is a multilinear form on the STRICT
          simplex t1<...<tk with kernel
              K_a(t1..tk) = head_a(t1) * prod_p Gm[t_p,t_{p+1}] * U[t_k].
          Its k "ports" are xi_{t1..tk}.  We NEVER materialize this k-tensor; we
          sweep its slots left-to-right (state = (slots_done, last_time)), and at
          each slot the port must be CONTRACTED (Wick-paired) with something.

  SOFT2 : chains of order 2.  Chain a  ==  xi^T S_a xi  with the zero-diagonal
          symmetric matrix  S_a = (A_a + A_a^T)/2,  A_a[t,t'] = head_a(t) Gm[t,t'] U(t')
          on t<t'  (identical to chaos_transfer._moment_quad's S).

  LEGS  : chains of order 1.  Leg j == b_j^T xi  with  b_j = head_j * U.

The moment is the Gaussian expectation
      M = E[ HARD(xi) * prod_a (xi^T S_a xi) * prod_j (b_j^T xi) ].
Wick's theorem = sum over complete matchings of ALL slots (hard ports + the 2
slots of each S_a + the 1 slot of each b_j), a pair at times t,t' contributing
the covariance  W[t,t'] = w[t] delta_{t,t'}.  Organize each matching by the graph
on the ATOM set.  The connected pieces are exactly those of the all-order<=2 case
(_moment_quad) PLUS pieces that touch hard ports:

  * closed soft pieces (no hard port), each consuming a disjoint set of soft
    atoms, multiply in as a scalar:
        - CYCLE over chains a1..am (m>=2):  2^{m-1} tr(S_{a1} W ... S_{am} W)
        - PATH  legu -(a1..am)- legv    :  2^m  b_u^T W S_{a1} W ... S_{am} W b_v
    (single-chain cycles vanish: zero diagonal.)  These are EXACTLY the pieces
    _moment_quad assembles; we reuse that machinery on the unused soft atoms.

  * BRIDGE piece: connects two hard ports at times t,t' through an ORDERED
    subset B={a1..am} of soft chains:
        Sigma_B[t,t'] = 2^{|B|} * sum_{orderings of B} (W S_{a1} W ... S_{am} W)[t,t']
    (|B|=0  ->  W[t,t'] : a bare hard-hard pairing.)

  * SOURCE piece: connects one hard port at time t to a leg-type v through an
    ORDERED subset B of soft chains:
        sigma_{B,v}[t] = 2^{|B|} * sum_{orderings of B} (W S_{a1} W ... S_{am} W b_v)[t]
    (|B|=0  ->  (W b_v)[t].)

All four weight rules (bridge 0S/1S/2S, source 0S/1S, closed cycle/path) are the
_moment_quad rules and are verified against the reference DP to ratio 1.0 in
_tmp_hybrid_proof.py.

The full moment is then, EXACTLY:

  M = sum over PARTITIONS of the soft atoms {S_a} and legs {b_j} into
        (a set of BRIDGE pieces) U (a set of SOURCE pieces) U (a closed remainder)
      of   [ Z_closed(remainder) ]
         * [ contraction of the hard ports against the labeled BAG of bridge/source
             pieces produced by this partition ].

The "contraction of hard ports against a bag" is a generalized Wick-with-mean
matching: every hard port is matched either to the OTHER end of a BRIDGE (pairing
two ports, kernel Sigma_B) or to a SOURCE (a single port, vector sigma_{B,v}),
each bag item used exactly once, and this respects each hard chain's internal
slot time-ordering.  It is evaluated by the hard simplex-DP whose only extra state
(beyond per-chain (slots_done,last_time)) is a bitmask of REMAINING bag items.

Because the number of soft atoms per moment is small (<=4 chains + a few legs in
the target workload), the partition sum and the bag are small; the hard sweep
keeps state ~ n^(#hard chains).  This module is the plain-numpy REFERENCE; speed
comes later.

-------------------------------------------------------------------------------
Public entry point:  moment_hybrid_core(gatoms)
  gatoms : list of (order:int>=1, headvec:np.ndarray length n)     [same as
           chaos_transfer._moment_bulk_gen].  Returns the moment as a float.
-------------------------------------------------------------------------------
"""
import numpy as np
from itertools import combinations, permutations, product as iproduct
import chaos_diagram as CD
import chaos_transfer as CT


def _geo():
    G = CD.GEO
    if 'G' not in G:
        raise RuntimeError("call chaos_diagram.setup(...) first")
    return G


# --------------------------------------------------------------------------
# soft-sector linear-algebra primitives (built once from the grid geometry)
# --------------------------------------------------------------------------
def _build_soft(soft_heads, leg_vecs):
    """Return (Ws, Ts, W, bws, btypes) where
         S_a  = (A+A^T)/2, A[t,t']=h(t)Gm[t,t']U(t') on t<t'    (zero diagonal)
         T_a  = W S_a   (used for bridges/sources as prefactor rows)
       and Sigma/sigma below are assembled from  W S_{a1} W ... .
    We store  S list, and W as a diagonal vector.  btypes = distinct leg vectors."""
    G = _geo()
    n = G['n']; w = G['w']; U = G['U']; Gm = G['G']
    triu = np.triu(np.ones((n, n)), 1)
    Ss = []
    for h in soft_heads:
        A = (h[:, None] * Gm * U[None, :]) * triu
        Ss.append(0.5 * (A + A.T))
    W = np.diag(w)
    return Ss, W, leg_vecs


# --------------------------------------------------------------------------
# soft connected-piece values via Held-Karp over soft-chain subsets
# (identical weights to chaos_transfer._moment_quad, reused here directly)
# --------------------------------------------------------------------------
class _SoftKernels:
    """Precompute, for EVERY subset B of soft chains, the summed ordered products
        P_B = 2^{|B|} * sum_{orderings of B}  W S_{a1} W S_{a2} W ... W S_{a|B|} W
    as an n x n matrix (B=empty -> W).  Then:
       bridge kernel for subset B between ports  =  P_B                (n x n)
       source vector for subset B into leg v     =  P_B @ b_v          (n,)
    Also the closed-soft scalar for any subset of chains + legs via the
    connected cycle/path assembly (delegated to _moment_quad-style logic)."""

    def __init__(self, Ss, W, btypes):
        self.Ss = Ss
        self.W = W
        self.n = W.shape[0]
        self.NC = len(Ss)
        self.btypes = btypes            # list of distinct leg vectors
        self._P = {}                    # subset-mask -> P_B  (n x n)
        self._build_P()

    def _build_P(self):
        W = self.W; Ss = self.Ss; NC = self.NC
        # P_empty = W
        self._P[0] = W.copy()
        # WS_a  building blocks
        WS = [W @ S for S in Ss]        # (W S_a)
        # Held-Karp: for subset B, sum over orderings of (2 W S_{a1} W ... W S_{am} W).
        # Represent Q_end[B][a] = sum over orderings of B ENDING at chain a of
        #    2^{|B|} * W S_{a1} W ... W S_{a} W    (already closed with trailing W).
        # Recurrence: Q_end[B][a] = ( sum_{b in B\a} Q_end[B\a][b] ) @ (S_a W) * 2 ... careful:
        # Instead accumulate the ordered product OP_end[B][a] = sum orderings of B ending at a of
        #    W S_{a1} W ... W S_{a} W    (NO 2^ factor); apply 2^{|B|} at the end.
        OP_end = {}
        for a in range(NC):
            OP_end[(1 << a) , a] = W @ Ss[a] @ W
        masks_by_pop = sorted([m for m in range(1, 1 << NC)], key=lambda m: bin(m).count('1'))
        for B in masks_by_pop:
            bits = [a for a in range(NC) if (B >> a) & 1]
            if len(bits) == 1:
                a = bits[0]
                # already seeded
                pass
            else:
                for a in bits:
                    Bm = B ^ (1 << a)
                    acc = None
                    for b in range(NC):
                        if (Bm >> b) & 1:
                            M = OP_end[(Bm, b)]
                            acc = M if acc is None else (acc + M)
                    # extend: ... W S_b W  -> ... W S_b W S_a W  ==  acc @ Ss[a] @ W
                    OP_end[(B, a)] = acc @ Ss[a] @ self.W
            # P_B = 2^{|B|} * sum_a OP_end[B][a]
            tot = None
            for a in bits:
                M = OP_end[(B, a)]
                tot = M if tot is None else (tot + M)
            self._P[B] = (2.0 ** len(bits)) * tot

    def bridge(self, B):
        return self._P[B]

    def source(self, B, v):
        return self._P[B] @ self.btypes[v]

    # ---- closed soft scalar on a subset of chains (mask) + a leg multiset ----
    def closed(self, chain_mask, leg_counts):
        """Z_closed = < prod_{a in chain_mask}(xi^T S_a xi) prod (b_v^T xi)^{leg_counts[v]} >
        assembled from cycles (over chains) and paths (leg -...- leg through chains),
        EXACTLY chaos_transfer._moment_quad's assembly restricted to these atoms.
        Uses the same 2^{...} weights.  Returns 0.0 if the total slot count is odd."""
        chains = [a for a in range(self.NC) if (chain_mask >> a) & 1]
        nlegs = sum(leg_counts)
        # build the gatom-style list and call _moment_quad directly (all order<=2)
        G = _geo()
        gatoms = []
        # reconstruct order-2 head is not needed: _moment_quad rebuilds S from head;
        # but we already have S.  Simpler: reuse the SAME connected assembly here.
        return self._closed_assemble(chain_mask, tuple(leg_counts))

    def _closed_assemble(self, chain_mask, leg_counts):
        # Delegate to a local connected-diagram assembler identical to _moment_quad,
        # but taking S-matrices and leg vectors directly.
        return _closed_soft_value(self.Ss, self.W, self.btypes, chain_mask, leg_counts)


# --------------------------------------------------------------------------
# closed soft value (cycles + leg-leg paths) -- mirrors _moment_quad exactly
# --------------------------------------------------------------------------
def _closed_soft_value(Ss, W, btypes, chain_mask, leg_counts):
    """< prod_{a in chain_mask} xi^T S_a xi  *  prod_v (b_v^T xi)^{leg_counts[v]} >.
    Connected pieces: cycles over chain-subsets, paths between two legs through a
    chain-subset.  Assembled by the same counted matching recursion as _moment_quad.
    """
    n = W.shape[0]
    chains = [a for a in range(len(Ss)) if (chain_mask >> a) & 1]
    NC = len(chains)
    remap = {a: i for i, a in enumerate(chains)}
    Ts = [W @ Ss[a] for a in chains]        # T = W S  (note: _moment_quad uses S W; symmetric-trace
                                            # invariant.  We keep W S for consistency with sources.)
    bw = [W @ b for b in btypes]            # W b_v
    NL = len(btypes)
    total_slots = 2 * NC + sum(leg_counts)
    if total_slots % 2 != 0:
        return 0.0
    if NC == 0 and sum(leg_counts) == 0:
        return 1.0

    bits = [[a for a in range(NC) if (S >> a) & 1] for S in range(1 << NC)]

    # ---- path aggregates: psum[u][S] = sum over ALL orderings of chain-subset S of
    #      bw[u]^T (S_{a1} W)(S_{a2} W)... ending vector, then dot b_v ----
    # We compute Pagg(u,S,v) = sum_orderings  b_u^T W S_{a1} W S_{a2} W ... S_{am} W b_v.
    # Held-Karp on end-vectors.
    Tchain = [Ss[a] @ W for a in chains]    # S_a W  (right-multiplication chain)
    psum = [None] * NL
    for u in range(NL):
        xend = [None] * (1 << NC)
        ps = [None] * (1 << NC)
        ps[0] = bw[u].copy()
        for a in range(NC):
            X = np.zeros((NC, n))
            X[a] = bw[u] @ Tchain[a]
            xend[1 << a] = X
        order = sorted(range(1, 1 << NC), key=lambda S: bin(S).count('1'))
        for S in order:
            sb = bits[S]
            if len(sb) == 1:
                ps[S] = xend[S][sb[0]]
                continue
            X = np.zeros((NC, n))
            for a in sb:
                Sm = S ^ (1 << a)
                prev = xend[Sm]
                acc = np.zeros(n)
                for b in bits[Sm]:
                    acc = acc + prev[b]
                X[a] = acc @ Tchain[a]
            xend[S] = X
            ps[S] = X.sum(axis=0)
        psum[u] = ps

    def Pagg(u, S, v):
        return float(np.dot(psum[u][S], btypes[v]))

    # ---- cycle aggregates: Cycagg[C] = sum over orderings of C\{min} of
    #      tr(T_{min} T_ordering), T = S W ----
    TSW = [Ss[a] @ W for a in chains]
    Cycagg = {}
    for C in range(1, 1 << NC):
        cb = bits[C]
        if len(cb) < 2:
            continue
        c0 = cb[0]
        others = C ^ (1 << c0)
        ob = bits[others]
        Qend = {}
        for a in ob:
            Qend[(1 << a)] = {a: TSW[a]}
        osub = sorted([S for S in range(1, 1 << NC) if (S & ~others) == 0 and S],
                      key=lambda S: bin(S).count('1'))
        for S in osub:
            sb = bits[S]
            if len(sb) == 1:
                continue
            d = {}
            for a in sb:
                Sm = S ^ (1 << a)
                prevd = Qend[Sm]
                acc = None
                for b in bits[Sm]:
                    M = prevd[b]
                    acc = M if acc is None else (acc + M)
                d[a] = acc @ TSW[a]
            Qend[S] = d
        prevd = Qend[others]
        acc = None
        for b in ob:
            M = prevd[b]
            acc = M if acc is None else (acc + M)
        Cycagg[C] = float(np.trace(TSW[c0] @ acc))

    all_chains = frozenset(range(NC))

    def _mask(chs):
        m = 0
        for c in chs:
            m |= (1 << c)
        return m

    memo = {}

    def assemble(counts, chset):
        if not any(counts) and not chset:
            return 1.0
        key = (counts, chset)
        v = memo.get(key)
        if v is not None:
            return v
        tot = 0.0
        cmask = _mask(chset)
        if any(counts):
            u = next(i for i, c in enumerate(counts) if c > 0)
            c1 = list(counts); c1[u] -= 1
            for vt in range(len(counts)):
                nv = c1[vt]
                if nv == 0:
                    continue
                c2 = list(c1); c2[vt] -= 1; c2t = tuple(c2)
                mult = nv
                sub = cmask
                while True:
                    rest = frozenset(chset - {c for c in range(NC) if (sub >> c) & 1}) if sub else chset
                    r = bin(sub).count('1')
                    val = Pagg(u, sub, vt)
                    if val != 0.0:
                        tot += mult * (2.0 ** r) * val * assemble(c2t, rest)
                    if sub == 0:
                        break
                    sub = (sub - 1) & cmask
        else:
            c0 = min(chset)
            c0bit = 1 << c0
            othersmask = cmask ^ c0bit
            sub = othersmask
            while True:
                if sub:
                    C = c0bit | sub
                    rest = frozenset(chset - {c0} - {c for c in range(NC) if (sub >> c) & 1})
                    ca = Cycagg.get(C)
                    if ca:
                        r = bin(sub).count('1')
                        tot += (2.0 ** (r + 1 - 1)) * ca * assemble(counts, rest)
                if sub == 0:
                    break
                sub = (sub - 1) & othersmask
        memo[key] = tot
        return tot

    return assemble(tuple(leg_counts), all_chains)


# --------------------------------------------------------------------------
# enumerate soft partitions into { bridges, sources, closed remainder }
# --------------------------------------------------------------------------
def _iter_soft_partitions(NC, leg_counts, nports):
    """Yield (bag, closed_chain_mask, closed_leg_counts).

    `bag` is a list of pieces, each:
        ('bridge', chain_subset_mask)              -> arity 2 (pairs two hard ports)
        ('source', chain_subset_mask, legtype v)   -> arity 1 (one hard port -> leg v)
    The bag must have exactly  nb bridges and ns sources with  2*nb + ns == nports
    (so the hard ports can be matched).  Bridge/source chain subsets are DISJOINT;
    each source consumes one leg of its type; leftover chains -> closed_chain_mask,
    leftover legs -> closed_leg_counts.

    Chain subsets 'B' are UNORDERED: the sum over internal orderings is already
    inside Sigma_B / sigma_{B,v}.  Empty subsets ARE allowed (bare hard-hard bond
    W, or bare port->leg Wb).  To avoid double-counting interchangeable bridges we
    generate bridge chain-subsets in a canonical (non-decreasing bitmask) order and
    divide out nothing (canonical order already yields each multiset once).  For
    sources we enumerate the (legtype, chain-subset) content; two sources with the
    same content are also generated canonically.
    """
    NL = len(leg_counts)
    total_legs = sum(leg_counts)
    all_chain_mask = (1 << NC) - 1

    # choose number of sources ns and bridges nb
    for ns in range(0, nports + 1):
        rem = nports - ns
        if rem % 2 != 0:
            continue
        nb = rem // 2
        if ns > total_legs:
            continue
        # choose which leg types the ns sources terminate on (a multiset over types
        # with per-type count <= leg_counts[v]); ns sources total.
        for leg_choice in _multisets_leg(NL, leg_counts, ns):
            # leg_choice: tuple length NL, sum == ns, entry <= leg_counts[v]
            # Now distribute soft chains among: nb bridges + ns sources + closed.
            # Sources are distinguishable by their leg type; bridges are
            # interchangeable.  We assign each chain to one "slot": one of the nb
            # bridge slots, one of the ns source slots, or closed.  To avoid
            # over/under counting from interchangeable slots we canonicalize AFTER
            # building.  Simplest exact-and-clear: enumerate all assignments of each
            # chain to a labeled slot, build the bag, and DEDUPE by canonical key.
            src_slots = []
            for v in range(NL):
                for _ in range(leg_choice[v]):
                    src_slots.append(v)
            # labeled slots: 0..nb-1 bridges, nb..nb+ns-1 sources(with leg type)
            nslots = nb + ns
            seen = set()
            for assign in _iter_chain_assign(NC, nslots + 1):
                # assign[a] in 0..nslots ; nslots == closed
                bridges = [0] * nb
                sources = [0] * ns
                closed_chain_mask = 0
                ok = True
                for a in range(NC):
                    s = assign[a]
                    if s < nb:
                        bridges[s] |= (1 << a)
                    elif s < nslots:
                        sources[s - nb] |= (1 << a)
                    else:
                        closed_chain_mask |= (1 << a)
                bag = []
                for m in bridges:
                    bag.append(('bridge', m))
                for i in range(ns):
                    bag.append(('source', sources[i], src_slots[i]))
                # canonical key: multiset of bridge masks + multiset of (legtype,mask)
                bkey = tuple(sorted(bridges))
                skey = tuple(sorted((src_slots[i], sources[i]) for i in range(ns)))
                key = (bkey, skey)
                if key in seen:
                    continue
                seen.add(key)
                closed_legs = tuple(leg_counts[v] - leg_choice[v] for v in range(NL))
                yield bag, closed_chain_mask, closed_legs


def _multisets_leg(NL, leg_counts, ns):
    """all tuples c length NL, 0<=c[v]<=leg_counts[v], sum(c)==ns."""
    def rec(v, left):
        if v == NL:
            if left == 0:
                yield ()
            return
        for k in range(min(leg_counts[v], left) + 1):
            for tail in rec(v + 1, left - k):
                yield (k,) + tail
    yield from rec(0, ns)


def _iter_chain_assign(NC, nvals):
    """all tuples length NC with entries in range(nvals)."""
    if NC == 0:
        yield ()
        return
    for combo in iproduct(range(nvals), repeat=NC):
        yield combo


# --------------------------------------------------------------------------
# hard-port contraction against a labeled bag (bridges + sources), via the
# simplex-DP over hard chains.  State = per-chain (slots_done,last_time) +
# remaining-bag bitmask.  Each port matches ONE remaining bag item:
#   - a SOURCE item (arity 1): consume it, weight sigma_item[t].
#   - a BRIDGE item (arity 2): a port "opens" it (records the open end's time as a
#     pending half-edge keyed by item), a later port "closes" it with kernel
#     Sigma_item[t_open,t].  To keep it a forward DP with small state, we treat a
#     bridge as: the FIRST port to hit it stores its time in the state; the SECOND
#     port closes it.  We therefore also track, per still-open bridge, the open
#     time.  Since bridges are few, we carry a tuple of (item_idx, open_time).
# --------------------------------------------------------------------------
HYBRID_SWEEP = True   # use the vectorized tensor-return contractor (kills open-time blowup)


def _contract_hard_sweep(hard, bag, SK):
    """FAST rewrite of _contract_hard (chaos_hybrid.HYBRID_SWEEP).

    SAME recursion as _contract_hard (chain-major, slot-in-chain; each port either
    opens a bridge, closes an open bridge, or consumes a source), but instead of
    memoizing every half-open bridge's SCALAR open-time in the state key (which
    reintroduced n^(#open bridges) and walled v7), `solve` RETURNS A TENSOR over
    the open bridges' open-times and memoizes only on the DISCRETE state
    (ci, si, lastt, remain, open-bridge indices).  The open-time dependence lives
    on numpy axes and is summed/contracted in C, not enumerated in Python.

    Axis convention: solve(...) returns an ndarray of shape (n,)*len(opens); axis i
    is the open-time of bridge  opens[i]  (opens is a sorted tuple of bag indices)."""
    G = _geo()
    n = G['n']; U = G['U']; Gm = G['G']
    NC = len(hard)
    korder = tuple(o for o, _ in hard)
    heads = [h for _, h in hard]
    nports = sum(korder)

    nbr = sum(1 for p in bag if p[0] == 'bridge')
    nsr = sum(1 for p in bag if p[0] == 'source')
    if 2 * nbr + nsr != nports:
        return 0.0
    if NC == 0:
        return 1.0

    # bag kernels/vectors, indexed by bag position
    bridge_ker = {}      # bag_idx -> n x n
    source_vec = {}      # bag_idx -> n
    is_bridge = {}
    for bidx, p in enumerate(bag):
        if p[0] == 'bridge':
            bridge_ker[bidx] = SK.bridge(p[1]); is_bridge[bidx] = True
        else:
            source_vec[bidx] = SK.source(p[1], p[2]); is_bridge[bidx] = False
    NB = len(bag)

    memo = {}

    def solve(ci, si, lastt, remain, opens):
        # advance past finished chains
        while ci < NC and si == korder[ci]:
            ci += 1; si = 0; lastt = -1
        if ci == NC:
            if not remain and not opens:
                return np.array(1.0)               # 0-dim
            return None                             # dead branch (0)
        key = (ci, si, lastt, remain, opens)
        v = memo.get(key)
        if v is not None:
            return v if v is not False else None
        L = len(opens)
        R = np.zeros((n,) * L)
        head = heads[ci]
        is_last = (si + 1 == korder[ci])
        rem_slots_here = korder[ci] - si          # incl this one
        for t in range(lastt + 1, n - (rem_slots_here - 1)):
            f = head[t] if si == 0 else Gm[lastt, t]
            if is_last:
                f = f * U[t]
            if f == 0.0:
                continue
            # (a) close an open bridge at position p
            for p in range(L):
                oi = opens[p]
                col = bridge_ker[oi][:, t]         # over oi's open-time
                child = opens[:p] + opens[p + 1:]
                sub = solve(ci, si + 1, t, remain, child)   # axes = child (L-1)
                if sub is None:
                    continue
                col_r = col.reshape([n if a == p else 1 for a in range(L)])
                sub_e = np.expand_dims(sub, p) if L > 1 else col_r * 0 + sub  # insert axis p
                R += f * col_r * sub_e
            # (b) open a new bridge (from remaining bridge items)
            for bidx in remain:
                if not is_bridge[bidx]:
                    continue
                child = tuple(sorted(opens + (bidx,)))
                pos = child.index(bidx)
                sub = solve(ci, si + 1, t, remain - {bidx}, child)  # axes = child (L+1)
                if sub is None:
                    continue
                R += f * np.take(sub, t, axis=pos)  # this bridge opened at t -> index its axis
            # (c) consume a source
            for bidx in remain:
                if is_bridge[bidx]:
                    continue
                sv = source_vec[bidx][t]
                if sv == 0.0:
                    continue
                sub = solve(ci, si + 1, t, remain - {bidx}, opens)  # axes = opens (L)
                if sub is None:
                    continue
                R += (f * sv) * sub
        memo[key] = R if R.any() else False
        return R if R.any() else None

    out = solve(0, 0, -1, frozenset(range(NB)), ())
    return 0.0 if out is None else float(out)


def _contract_hard(hard, bag, SK):
    """hard: list of (order, headvec) with order>=3.  bag: list of pieces
    (('bridge',mask) or ('source',mask,v)).  SK: _SoftKernels.
    Returns the contraction of all hard ports against exactly-this bag (every bag
    item used exactly once, every hard port matched)."""
    if HYBRID_SWEEP:
        return _contract_hard_sweep(hard, bag, SK)
    G = _geo()
    n = G['n']; U = G['U']; Gm = G['G']
    NC = len(hard)
    korder = [o for o, _ in hard]
    heads = [h for _, h in hard]
    nports = sum(korder)

    # arity check: #ports must equal 2*#bridges + #sources
    nbr = sum(1 for p in bag if p[0] == 'bridge')
    nsr = sum(1 for p in bag if p[0] == 'source')
    if 2 * nbr + nsr != nports:
        return 0.0
    if NC == 0:
        # no hard ports: bag must be empty (handled by arity check) -> factor 1
        return 1.0

    # Precompute bag kernels/vectors.
    bridge_ker = []      # list of n x n
    source_vec = []      # list of n
    bag_kind = []        # 'B' or 'S' with index into bridge_ker/source_vec
    for p in bag:
        if p[0] == 'bridge':
            bridge_ker.append(SK.bridge(p[1]))
            bag_kind.append(('B', len(bridge_ker) - 1))
        else:
            source_vec.append(SK.source(p[1], p[2]))
            bag_kind.append(('S', len(source_vec) - 1))
    NB = len(bag)

    # DP state:
    #   phase   : tuple over chains of (slots_done, last_time)   (done -> (korder,-1))
    #   remain  : frozenset of bag-item indices not yet consumed
    #   opens   : tuple sorted of (bag_item_idx, open_time) for bridges half-open
    # We advance ONE hard slot per micro-step? No: the simplex-DP advances subsets
    # of chains at each grid time t, but here EACH port must contract individually
    # (it pairs to a bag item), and a bridge pairs two ports possibly on different
    # chains / different times.  So we cannot use the (m-1)!! bulk weight; we must
    # enumerate port-to-bag matchings explicitly.  We therefore advance the DP
    # slot-by-slot in a fixed port order (chain-major, slot within chain), and at
    # each port choose its bag partner.  The chain's simplex ordering is enforced
    # by requiring the port's grid time > last_time of that chain, and the port's
    # kernel factor uses Gm[last_time, t] (or head at first slot) * U at last slot.
    #
    # Port order: we process chains in order; within a chain, slots in order; but a
    # bridge may pair a port of chain c1 with a LATER-processed port of chain c2.
    # We handle this by the 'opens' set: opening a bridge stores the time; closing
    # it later applies Sigma[open_time, t].  A source is applied immediately.

    # Enumerate ports in a canonical order and recurse.  To bound the state we do a
    # straightforward recursive DP with memoization over (chain progress, remain,
    # opens).  Grid times are summed inside.

    from functools import lru_cache
    heads_t = heads
    korder_t = tuple(korder)

    # We recurse chain-by-chain, slot-by-slot, carrying:
    #   ci        : current chain index
    #   si        : current slot index within chain ci
    #   lastt     : last grid time used in chain ci (-1 if si==0)
    #   remain    : tuple(sorted) remaining bag indices
    #   opens     : tuple(sorted) of (bag_idx, open_time) currently half-open
    # We sum over the grid time t for this port and its bag choice.
    # Amplitude accumulates the chain kernel factors.

    memo = {}

    def solve(ci, si, lastt, remain, opens):
        # advance to next chain if this one is done
        if si == korder_t[ci]:
            ci += 1; si = 0; lastt = -1
        if ci == NC:
            # all ports placed; must have consumed all bag items and no open bridges
            if not remain and not opens:
                return 1.0
            return 0.0
        key = (ci, si, lastt, remain, opens)
        v = memo.get(key)
        if v is not None:
            return v
        tot = 0.0
        head = heads_t[ci]
        is_last = (si + 1 == korder_t[ci])
        # choose grid time t for this port (must exceed lastt for simplex ordering;
        # first slot t in [0,n)); also must leave room for remaining slots of chain.
        rem_slots_here = korder_t[ci] - si  # including this one
        for t in range(lastt + 1, n - (rem_slots_here - 1)):
            # chain kernel factor for this slot
            if si == 0:
                f = head[t]
            else:
                f = Gm[lastt, t]
            if is_last:
                f = f * U[t]
            if f == 0.0:
                continue
            # bag choice for this port:
            # (a) close an OPEN bridge:
            for oi, (bidx, ot) in enumerate(opens):
                kind, kk = bag_kind[bidx]
                # only bridges are in opens
                ker = bridge_ker[kk][ot, t]
                if ker == 0.0:
                    continue
                newopens = opens[:oi] + opens[oi + 1:]
                sub = solve(ci, si + 1, t, remain, newopens)
                if sub != 0.0:
                    tot += f * ker * sub
            # (b) OPEN a new bridge (from remaining bridge items):
            for bidx in remain:
                kind, kk = bag_kind[bidx]
                if kind != 'B':
                    continue
                newremain = tuple(x for x in remain if x != bidx)
                newopens = tuple(sorted(opens + ((bidx, t),)))
                sub = solve(ci, si + 1, t, newremain, newopens)
                if sub != 0.0:
                    tot += f * sub   # kernel applied at close time
            # (c) consume a SOURCE:
            for bidx in remain:
                kind, kk = bag_kind[bidx]
                if kind != 'S':
                    continue
                sv = source_vec[kk][t]
                if sv == 0.0:
                    continue
                newremain = tuple(x for x in remain if x != bidx)
                sub = solve(ci, si + 1, t, newremain, opens)
                if sub != 0.0:
                    tot += f * sv * sub
        memo[key] = tot
        return tot

    remain0 = tuple(range(NB))
    return solve(0, 0, -1, remain0, tuple())


# --------------------------------------------------------------------------
# public entry
# --------------------------------------------------------------------------
def moment_hybrid_core(gatoms):
    """Exact moment of  prod atoms, atoms = (order, headvec).  Same semantics as
    chaos_transfer._moment_bulk_gen.  Folds the soft sector (order-2 chains + legs)
    analytically and sweeps only the order>=3 hard chains."""
    G = _geo()
    n = G['n']; U = G['U']
    if not gatoms:
        return 1.0
    if sum(o for o, _ in gatoms) % 2 != 0:
        return 0.0

    hard = [(o, h) for (o, h) in gatoms if o >= 3]
    soft_heads = [h for (o, h) in gatoms if o == 2]
    legvecs_raw = [h * U for (o, h) in gatoms if o == 1]

    # dedupe legs by content -> counted types
    tkey = {}; btypes = []; leg_counts = []
    for b in legvecs_raw:
        k = b.tobytes()
        if k not in tkey:
            tkey[k] = len(btypes); btypes.append(b); leg_counts.append(0)
        leg_counts[tkey[k]] += 1
    leg_counts = tuple(leg_counts)

    # no hard chains -> pure soft: delegate to the validated quad path
    if not hard:
        return CT._moment_quad(gatoms)

    Ss, W, _ = _build_soft(soft_heads, btypes)
    NC = len(Ss)
    SK = _SoftKernels(Ss, W, btypes)

    nports = sum(o for o, _ in hard)
    total = 0.0
    for bag, closed_chain_mask, closed_legs in _iter_soft_partitions(NC, leg_counts, nports):
        zc = SK.closed(closed_chain_mask, closed_legs)
        if zc == 0.0:
            continue
        # combinatorial weight of this bag:
        #  - the port-DP consumes bag pieces as DISTINCT labeled items, overcounting
        #    interchangeable identical pieces by prod (multiplicity)!  -> divide.
        #  - legs of the same type are identical ATOMS; a source terminating on type
        #    v consumes one physical leg.  For used_v sources on type v out of
        #    leg_counts[v] identical legs, the number of distinct physical-leg
        #    assignments to the (labeled) source slots is  P(leg_counts[v], used_v)
        #    = leg_counts[v]! / (leg_counts[v]-used_v)!  -> multiply.
        from math import factorial
        # multiplicity division (identical pieces by content+type)
        cnt = {}
        for p in bag:
            cnt[p] = cnt.get(p, 0) + 1
        sym = 1
        for c in cnt.values():
            sym *= factorial(c)
        # leg permutation multiply
        used = [0] * len(leg_counts)
        for p in bag:
            if p[0] == 'source':
                used[p[2]] += 1
        legperm = 1
        for v in range(len(leg_counts)):
            lc = leg_counts[v]; uv = used[v]
            legperm *= factorial(lc) // factorial(lc - uv)
        weight = legperm / sym
        hc = _contract_hard(hard, bag, SK)
        if hc == 0.0:
            continue
        total += zc * hc * weight
    return float(total)
