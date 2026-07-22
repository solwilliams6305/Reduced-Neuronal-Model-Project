"""Independent ground-truth moment via discrete-Gaussian Isserlis (no multigraph, no transfer DP).

White noise on the grid: g[i], <g_i g_j> = w_i delta_ij.
Bulk atom U_{k,m} = sum_{i1<...<ik} chain(i1..ik) g_{i1}..g_{ik}   (ordered, multilinear).
Boundary s_j defined so <s_j , I_1(f)> = bfac * d^j f(node): s_j = sum_i (bfac * D^j[node,i]/w_i) g_i,
where D^j is the j-th finite-difference (np.gradient) operator matrix.  Moment of a product of atoms
= Gaussian expectation of the product polynomial (Wick/Isserlis over the g's).
"""
import numpy as np, itertools
import chaos_diagram as CD

def _grad_matrix(n, h, j):
    M = np.eye(n)
    for _ in range(j):
        M = np.gradient(M, h, axis=0)
    return M   # M[a,i] = (d^j/dx^j applied to e_i)[a]

def _atom_poly(a, G):
    """return dict {sorted tuple of g-indices : coeff} for one atom."""
    n = G['n']; U = G['U']; Gm = G['G']; w = G['w']; h = G['h']
    if a[0] == 's':
        j = a[1]; bfac = G.get('bfac', 0.5)
        D = _grad_matrix(n, h, j)
        poly = {}
        for i in range(n):
            c = bfac * D[0, i] / w[i]
            if c != 0.0:
                poly[(i,)] = poly.get((i,), 0.0) + c
        return poly
    # bulk U_{k,m}
    _, k, m = a
    head = G['Gnode'] if m == 0 else G['H']
    poly = {}
    if k == 1:
        for i in range(n):
            poly[(i,)] = head[i] * U[i]
        return poly
    for idx in itertools.combinations(range(n), k):
        val = head[idx[0]]
        for p in range(k - 1):
            val *= Gm[idx[p], idx[p + 1]]
        val *= U[idx[-1]]
        poly[idx] = val
    return poly

def _poly_mul(p, q):
    out = {}
    for ka, va in p.items():
        for kb, vb in q.items():
            key = tuple(sorted(ka + kb))
            out[key] = out.get(key, 0.0) + va * vb
    return out

def _gauss_expect(monomial, w):
    """<prod g_{i}> for a multiset of indices (tuple) with <g_i g_j>=w_i delta_ij. Wick."""
    idx = list(monomial)
    if len(idx) % 2 != 0:
        return 0.0
    # group counts
    from collections import Counter
    cnt = Counter(idx)
    # <prod g_i^{c_i}> = prod over i of ( (c_i-1)!! w_i^{c_i/2} ) if all c_i even else 0
    val = 1.0
    for i, c in cnt.items():
        if c % 2 != 0:
            return 0.0
        # (c-1)!!
        d = 1
        kk = c - 1
        while kk > 1:
            d *= kk; kk -= 2
        val *= d * (w[i] ** (c // 2))
    return val

def brute_moment(atoms):
    G = CD.GEO
    if not atoms:
        return 1.0
    w = G['w']
    poly = _atom_poly(atoms[0], G)
    for a in atoms[1:]:
        poly = _poly_mul(poly, _atom_poly(a, G))
    tot = 0.0
    for mono, coef in poly.items():
        if coef == 0.0:
            continue
        tot += coef * _gauss_expect(mono, w)
    return tot
