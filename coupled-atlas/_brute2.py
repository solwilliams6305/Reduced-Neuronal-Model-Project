"""Definitive independent moment: leg-level Wick with boundary-boundary DROPPED (renormalization).

Enumerate every perfect matching of the legs (bulk atom I_k -> k legs; boundary s_j -> 1 leg),
forbidding intra-atom pairs and boundary-boundary pairs.  For each matching, contract:
  value = sum over one grid index per edge of  prod_edges w[x_e]
          * prod_bulk h_a(its legs' indices) * prod_bnd c^{(j)}(its leg's index),
with raw symmetric bulk kernels h_a and boundary vector c^{(j)}[x] = bfac * D^j[node,x]/w[x].
Independent of dense's multigraph grouping and of the transfer DP.
"""
import numpy as np, itertools
import chaos_diagram as CD
from chaos_engine import _sym

def _bulk_raw(a, G):
    n = G['n']; U = G['U']; Gm = G['G']
    _, k, m = a
    head = G['Gnode'] if m == 0 else G['H']
    if k == 1:
        return (head * U).copy()
    # ordered chain tensor
    triu = np.triu(np.ones((n, n)), 1)
    Gmask = Gm * triu
    P = head.copy()
    for _ in range(2, k + 1):
        P = P[..., None] * Gmask
    T = P * U
    return _sym(T)   # raw symmetric, no weights

def _bnd_vec(a, G):
    n = G['n']; w = G['w']; h = G['h']; bfac = G.get('bfac', 0.5)
    j = a[1]
    D = np.eye(n)
    for _ in range(j):
        D = np.gradient(D, h, axis=0)
    return bfac * D[0] / w    # c[x]

def _matchings(pool, legs):
    """pool: tuple of original leg-indices remaining. yield perfect matchings as lists of (i,j) pairs
    (original indices), no same atom, no boundary-boundary."""
    if not pool:
        yield []
        return
    first = pool[0]
    a0 = legs[first]
    for kk in range(1, len(pool)):
        j = pool[kk]
        aj = legs[j]
        if a0[0] == aj[0]:
            continue                      # same atom
        if a0[1] and aj[1]:
            continue                      # bnd-bnd
        rest = tuple(pool[t] for t in range(1, len(pool)) if t != kk)
        for sub in _matchings(rest, legs):
            yield [(first, j)] + sub

def brute2_moment(atoms):
    G = CD.GEO; n = G['n']; w = G['w']
    # build legs and kernels
    kernels = []; legs = []; isbnd = []
    for aid, a in enumerate(atoms):
        if a[0] == 's':
            kernels.append(('b', _bnd_vec(a, G)))
            legs.append((aid, True))
        else:
            k = a[1]
            kernels.append(('u', _bulk_raw(a, G), k))
            for _ in range(k):
                legs.append((aid, False))
    L = len(legs)
    if L % 2 != 0:
        return 0.0
    tot = 0.0
    letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQR'
    for match in _matchings(tuple(range(L)), legs):
        # assign an edge letter to each pair; each leg gets its edge's letter
        leg_letter = [None] * L
        edge_letters = []
        for ei, (i, j) in enumerate(match):
            leg_letter[i] = letters[ei]; leg_letter[j] = letters[ei]
            edge_letters.append(letters[ei])
        # build einsum: per atom, gather its legs' letters
        ops = []; subs = []
        # map atom -> its leg positions in order
        atom_legpos = {}
        for pos, (aid, b) in enumerate(legs):
            atom_legpos.setdefault(aid, []).append(pos)
        wvec_needed = edge_letters
        for aid, a in enumerate(atoms):
            kk = kernels[aid]
            positions = atom_legpos[aid]
            lett = ''.join(leg_letter[p] for p in positions)
            if kk[0] == 'b':
                ops.append(kk[1]); subs.append(lett)
            else:
                ops.append(kk[1]); subs.append(lett)
        # weight per edge
        for el in edge_letters:
            ops.append(w); subs.append(el)
        val = np.einsum(','.join(subs) + '->', *ops, optimize='greedy')
        tot += float(val)
    return tot
