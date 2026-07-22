"""Prototype + validate a FAST dense-boundary path for all-order<=2 boundary moments:
reduce each boundary leg with the validated CD.beval (numerically, on small tensors), representing
an order-2 chain reduced by a boundary leg as an order-1 'beval-leg' (its un-absorbed lval), then
contract the remaining order<=2 bulk via the quad Gaussian path. Compare to CD.moment (dense)."""
import os, sys, time, itertools
os.environ.setdefault('CT_NUMBA', '0')
import numpy as np
import chaos_diagram as CD
import chaos_transfer as CT

n = 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
G = CD.GEO; U = G['U']; Gm = G['G']; sw = G['sw']; hgrid = G['h']; bfac = G.get('bfac', 0.5)
triu = np.triu(np.ones((n, n)), 1)

def chain2_tensor(head):
    """weight-absorbed symmetric order-2 tensor for a chain with this head."""
    A = (head[:, None] * Gm * U[None, :]) * triu     # head(s1)G(s1,s2)u0(s2), s1<s2
    T = 0.5 * (A + A.T)                                # symmetrize
    return T * sw[:, None] * sw[None, :]              # absorb sw on both axes

def beval_leg2(head, j):
    """order-2 chain reduced by s_j -> order-1 'beval-leg' lval (un-absorbed)."""
    T = chain2_tensor(head)
    vA = CD.beval(T, 0, j)                             # weight-absorbed vector over s2
    return vA / sw                                     # un-absorb -> leg lval (= head'*u0 analogue)

def reduce_scalar_lval(lval, j):
    """boundary s_j reduces an order-1 leg given by its lval (=head*u0): bfac * d^j[lval]|node."""
    v = lval.copy()
    for _ in range(j):
        v = np.gradient(v, hgrid)
    return bfac * v[0]

def quad_ext(order2_heads, leg_lvals):
    """quad Gaussian moment of (order-2 chains + order-1 legs-by-lval).  Mirrors _moment_quad but
    takes leg lvals directly (so beval-legs need no head)."""
    w = G['w']; NC = len(order2_heads)
    # chain matrices S=(A+A^T)/2, T=S W
    Ts = []
    for h in order2_heads:
        A = (h[:, None] * Gm * U[None, :]) * triu
        S = 0.5 * (A + A.T)
        Ts.append(S * w[None, :])
    # dedupe legs
    tkey = {}; btypes = []; bcount = []
    for b in leg_lvals:
        kb = b.tobytes()
        if kb not in tkey:
            tkey[kb] = len(btypes); btypes.append(b); bcount.append(0)
        bcount[tkey[kb]] += 1
    if sum(len(order2_heads) * [2]) % 2 != (sum(bcount)) % 2:
        pass
    # total legs parity
    if (2 * NC + sum(bcount)) % 2 != 0:
        return 0.0
    bw = [b * w for b in btypes]
    # reuse the exact assembler from _moment_quad by calling it via a shim: build gatoms with
    # standard order-2 chains and represent beval-legs as order-1 atoms whose h*U == lval.
    # Simpler: replicate _moment_quad's assemble on these Ts and btypes.
    from math import comb as _comb
    NL = len(bcount); full = (1 << NC) - 1
    bits = [[a for a in range(NC) if (S >> a) & 1] for S in range(1 << NC)]
    psum = [None] * NL
    for u in range(NL):
        xend = [None] * (1 << NC); ps = [None] * (1 << NC); ps[0] = bw[u].copy()
        for a in range(NC):
            X = np.zeros((NC, n)); X[a] = bw[u] @ Ts[a]; xend[1 << a] = X
        order = sorted(range(1, 1 << NC), key=lambda S: bin(S).count('1'))
        for S in order:
            sb = bits[S]
            if len(sb) == 1:
                ps[S] = xend[S][sb[0]]; continue
            X = np.zeros((NC, n))
            for a in sb:
                Sm = S ^ (1 << a); prev = xend[Sm]; acc = np.zeros(n)
                for b in bits[Sm]: acc = acc + prev[b]
                X[a] = acc @ Ts[a]
            xend[S] = X; ps[S] = X.sum(axis=0)
        psum[u] = ps
    def Pagg(u, S, v): return float(np.dot(psum[u][S], btypes[v]))
    Cycagg = {}
    for C in range(1, 1 << NC):
        cb = bits[C]
        if len(cb) < 2: continue
        c0 = cb[0]; others = C ^ (1 << c0); ob = bits[others]
        Qend = {}
        for a in ob: Qend[(1 << a)] = {a: Ts[a]}
        osub = sorted([S for S in range(1, 1 << NC) if (S & ~others) == 0 and S], key=lambda S: bin(S).count('1'))
        for S in osub:
            sb = bits[S]
            if len(sb) == 1: continue
            d = {}
            for a in sb:
                Sm = S ^ (1 << a); prevd = Qend[Sm]; acc = None
                for b in bits[Sm]:
                    M = prevd[b]; acc = M if acc is None else (acc + M)
                d[a] = acc @ Ts[a]
            Qend[S] = d
        prevd = Qend[others]; acc = None
        for b in ob:
            M = prevd[b]; acc = M if acc is None else (acc + M)
        Cycagg[C] = float(np.trace(Ts[c0] @ acc))
    all_chains = frozenset(range(NC))
    def _mask(chains):
        m = 0
        for c in chains: m |= (1 << c)
        return m
    memo = {}
    def assemble(counts, chains):
        if not any(counts) and not chains: return 1.0
        key = (counts, chains); v = memo.get(key)
        if v is not None: return v
        tot = 0.0; cmask = _mask(chains)
        if any(counts):
            u = next(i for i, c in enumerate(counts) if c > 0)
            c1 = list(counts); c1[u] -= 1
            for vt in range(len(counts)):
                nv = c1[vt]
                if nv == 0: continue
                c2 = list(c1); c2[vt] -= 1; c2t = tuple(c2); mult = nv
                sub = cmask
                while True:
                    rest = frozenset(chains - {c for c in range(NC) if (sub >> c) & 1}) if sub else chains
                    r = bin(sub).count('1'); val = Pagg(u, sub, vt)
                    if val != 0.0: tot += mult * (2.0 ** r) * val * assemble(c2t, rest)
                    if sub == 0: break
                    sub = (sub - 1) & cmask
        else:
            c0 = min(chains); c0bit = 1 << c0; othersmask = cmask ^ c0bit
            sub = othersmask
            while True:
                if sub:
                    C = c0bit | sub
                    rest = frozenset(chains - {c0} - {c for c in range(NC) if (sub >> c) & 1})
                    ca = Cycagg.get(C)
                    if ca:
                        r = bin(sub).count('1'); tot += (2.0 ** (r + 1 - 1)) * ca * assemble(counts, rest)
                if sub == 0: break
                sub = (sub - 1) & othersmask
        memo[key] = tot; return tot
    return assemble(tuple(bcount), all_chains)

def dense_bnd_fast(chains, leg_lvals, bnds):
    """recursive boundary reduction (dense/beval convention) for all-order<=2 moments.
    chains: list of (2, head); leg_lvals: list of order-1 lval vectors; bnds: list of j."""
    if not bnds:
        return quad_ext([h for _, h in chains], leg_lvals)
    j = bnds[0]; rest_b = bnds[1:]; tot = 0.0
    # reduce onto each leg
    for i in range(len(leg_lvals)):
        scal = reduce_scalar_lval(leg_lvals[i], j)
        if scal != 0.0:
            tot += scal * dense_bnd_fast(chains, leg_lvals[:i] + leg_lvals[i+1:], rest_b)
    # reduce onto each order-2 chain -> becomes a beval-leg.  Factor 2 = the chain's order! that
    # CD.moment carries for an order-2 atom (the boundary can hit either of its 2 legs).
    for c in range(len(chains)):
        bl = beval_leg2(chains[c][1], j)
        tot += 2.0 * dense_bnd_fast(chains[:c] + chains[c+1:], leg_lvals + [bl], rest_b)
    return tot

# ---- validate vs CD.moment on all-order<=2 boundary moments ----
def moment_via_fast(atoms):
    chains = [(2, CT._head_vec(a)) for a in atoms if a[0] == 'U' and a[1] == 2]
    legs = [CT._head_vec(a) * U for a in atoms if a[0] == 'U' and a[1] == 1]
    bnds = [a[1] for a in atoms if a[0] == 's']
    return dense_bnd_fast(chains, legs, bnds)

tests = [
    [('U',2,0),('U',2,0),('U',1,0),('s',2)],
    [('U',2,0),('U',2,0),('s',2),('s',3)],
    [('U',2,0),('U',2,1),('U',1,0),('U',1,1),('s',2)],
    [('U',2,0),('U',2,0),('U',2,0),('U',1,0),('s',3)],
    [('U',2,0)]*3 + [('U',1,0)]*3 + [('s',2)],
    [('U',2,0)]*2 + [('U',1,0)]*5 + [('s',4)],
    [('U',2,1),('U',2,0),('U',1,0),('U',1,0),('s',2),('s',3)],
]
maxre = 0.0
for atoms in tests:
    if sum(a[1] if a[0]=='U' else 1 for a in atoms) % 2 != 0:
        continue
    CT.BND_ENGINE='dense'; CT.clear_all()
    vd = CD.moment(atoms)
    vf = moment_via_fast(atoms)
    re = abs(vd-vf)/max(abs(vd),abs(vf),1e-300); maxre=max(maxre,re)
    desc = [(a[1] if a[0]=='U' else ('s',a[1])) for a in atoms]
    print(f'  {desc}: dense={vd:+.5e} fast={vf:+.5e} re={re:.1e}', flush=True)
print(f'max relerr {maxre:.2e}  ' + ('PASS' if maxre < 1e-9 else 'FAIL'), flush=True)
