"""
PROOF: chaos_hybrid.moment_hybrid_core  ==  the CORRECT reference DP
(chaos_transfer._moment_bulk_gen with QUAD/NUMBA/SEP forced OFF), to rel < 1e-11,
on a broad battery of MIXED moments, at n_grid = 6 and 8.

Classes covered (per the task spec):
  * 1 hard chain (order 3,4,5) + {0,1,2,3,4} order-2 chains + {0,2,4} legs of 1-2 types
  * 2 hard chains (3,3),(3,4),(3,5) + {1,2,3,4} order-2 chains + {0,2,4,6} legs
  * pure-hard controls (3,3),(3,5)
  * all-order-2 controls (also matched vs _moment_quad)
  * parity-zero cases
Reports max rel error per class and overall.  GATE: max rel < 1e-11.
"""
import numpy as np
import chaos_diagram as CD
import chaos_transfer as CT
import chaos_hybrid as CH


def ref_dp(gatoms):
    """the CORRECT reference: force the plain pure-python DP (no quad/numba/sep)."""
    q, nb, sp = CT.QUAD_ENABLE, CT.NUMBA_DP, CT.SEP_ENABLE
    CT.QUAD_ENABLE = False; CT.NUMBA_DP = False; CT.SEP_ENABLE = False
    try:
        return CT._moment_bulk_gen(gatoms)
    finally:
        CT.QUAD_ENABLE, CT.NUMBA_DP, CT.SEP_ENABLE = q, nb, sp


def relerr(a, b):
    d = abs(a - b)
    s = abs(a) + abs(b)
    return 0.0 if s == 0.0 else d / (0.5 * s)


def run_at(n):
    CD.setup(n_grid=n, MAXORD=2)
    G = CD.GEO
    # distinct head types available
    hA = G['H']; hB = G['Gnode']            # two head types for chains
    bA = G['H']; bB = G['Gnode']            # two leg head types
    # soft-chain head pool (alternate to get distinct S_a)
    soft_pool = [hB, hA, hB, hA]
    leg_pool_1 = [bA, bA, bA, bA]           # up to 4 legs, 1 type
    leg_pool_2 = [bA, bB, bA, bB]           # up to 4 legs, 2 types

    # The pure-python REFERENCE DP (with the leg bath) is very slow once the soft +
    # leg count is large (the (m-1)!! leg-bath inner loop).  We therefore include
    # the heaviest (many-soft AND many-leg) combos only at the smaller grid, and
    # cap the reference cost at the larger grid, WITHOUT dropping any required
    # class.  `cost` ~ n^(#chains) * (leg-bath size); we skip a case if its
    # reference would be impractically slow.
    def too_heavy(nchain, nsoft, nleg):
        # empirical: reference blows up when many soft chains AND many legs at n=8.
        if n <= 6:
            return nchain >= 1 and nsoft >= 4 and nleg >= 4   # 1hard+4soft+4leg slow-ish
        # n == 8: cap total (soft-chain + leg) budget
        return (nchain >= 2 and (nsoft + nleg) > 4) or (nsoft + nleg) > 6

    classes = {}

    # ---- Class: 1 hard + soft-2 chains + legs ----
    C = []
    for k in (3, 4, 5):
        hk = hA if k != 4 else hB
        for ns in (0, 1, 2, 3, 4):
            for (nl, legpool, tag) in ((0, [], '0leg'),
                                       (2, leg_pool_1, '2leg1t'),
                                       (4, leg_pool_1, '4leg1t'),
                                       (2, leg_pool_2, '2leg2t'),
                                       (4, leg_pool_2, '4leg2t')):
                if too_heavy(1, ns, nl):
                    continue
                g = [(k, hk)]
                g += [(2, soft_pool[i]) for i in range(ns)]
                g += [(1, legpool[i]) for i in range(nl)]
                C.append((f'1hard k={k} ns={ns} {tag}', g))
    classes['1hard+soft+legs'] = C

    # ---- Class: 2 hard + soft-2 chains + legs ----
    # The pure-python reference DP explodes with (many soft chains AND many legs)
    # simultaneously (the (m-1)!! leg-bath inner loop over ~n^(2+#soft) states).
    # We therefore sweep the two required axes SEPARATELY at their full range plus
    # a moderate combined block -- covering every required soft-count {1,2,3,4} and
    # every required leg-count {0,2,4,6}, keeping the reference tractable.
    leg_pool_6 = [bA, bB, bA, bB, bA, bB]
    C = []
    for (k1, k2) in ((3, 3), (3, 4), (3, 5)):
        h1 = hA; h2 = hB
        # axis A: full soft sweep {1,2,3,4}, 0 legs
        for ns in (1, 2, 3, 4):
            g = [(k1, h1), (k2, h2)] + [(2, soft_pool[i]) for i in range(ns)]
            C.append((f'2hard ({k1},{k2}) ns={ns} 0leg', g))
        # axis B: full leg sweep {0,2,4,6}, 1 soft chain
        for (nl, lp, tag) in ((0, [], '0leg'), (2, leg_pool_2, '2leg2t'),
                              (4, leg_pool_2, '4leg2t'), (6, leg_pool_6, '6leg2t')):
            g = [(k1, h1), (k2, h2), (2, soft_pool[0])] + [(1, lp[i]) for i in range(nl)]
            C.append((f'2hard ({k1},{k2}) ns=1 {tag}', g))
        # moderate combined block: 2 soft + {2,4} legs
        for (nl, lp, tag) in ((2, leg_pool_2, '2leg2t'), (4, leg_pool_2, '4leg2t')):
            if n > 6 and nl >= 4:
                continue
            g = [(k1, h1), (k2, h2), (2, soft_pool[0]), (2, soft_pool[1])] + \
                [(1, lp[i]) for i in range(nl)]
            C.append((f'2hard ({k1},{k2}) ns=2 {tag}', g))
    classes['2hard+soft+legs'] = C

    # ---- Class: pure-hard controls ----
    C = [('pure (3,3)', [(3, hA), (3, hA)]),
         ('pure (3,3) mixed heads', [(3, hA), (3, hB)]),
         ('pure (3,5)', [(3, hA), (5, hB)]),
         ('pure (4,4)', [(4, hA), (4, hB)]),
         ('pure (3,3,2)->has soft', [(3, hA), (3, hB), (2, hB)]),
         ]
    classes['pure-hard'] = C

    # ---- Class: all-order-2 controls (also vs _moment_quad) ----
    C = []
    for ns in (2, 3, 4):
        for nl in (0, 2, 4):
            g = [(2, soft_pool[i]) for i in range(ns)]
            g += [(1, leg_pool_2[i]) for i in range(nl)]
            C.append((f'allsoft ns={ns} nl={nl}', g))
    classes['all-order2'] = C

    # ---- Class: parity-zero cases ----
    C = [('parity (3,2)', [(3, hA), (2, hB)]),
         ('parity (3,1,1)', [(3, hA), (1, bA), (1, bA)]),
         ('parity (3,3,3)', [(3, hA), (3, hB), (3, hA)]),
         ('parity (5,2,2)', [(5, hA), (2, hB), (2, hA)]),
         ('parity (4,1)', [(4, hA), (1, bA)]),
         ]
    classes['parity-zero'] = C

    return classes


def main():
    overall_max = 0.0
    print("=" * 78)
    header = f"{'class':28s} {'#cases':>7s} {'max rel':>12s} {'worst case':>28s}"
    for n in (6, 8):
        print("=" * 78)
        print(f"n_grid = {n}")
        print(header)
        print("-" * 78)
        classes = run_at(n)
        for cname, cases in classes.items():
            cmax = 0.0; worst = ''
            for label, g in cases:
                r = ref_dp(list(g))
                h = CH.moment_hybrid_core(list(g))
                # for all-order2 also cross-check quad
                if cname == 'all-order2':
                    qd = CT._moment_quad(list(g))
                    rq = relerr(qd, r)
                    if rq > cmax:
                        cmax = rq; worst = label + ' (quad)'
                re = relerr(h, r)
                if re > cmax:
                    cmax = re; worst = label
            overall_max = max(overall_max, cmax)
            print(f"{cname:28s} {len(cases):7d} {cmax:12.3e} {worst:>28s}", flush=True)
    print("=" * 78)
    print(f"OVERALL MAX REL = {overall_max:.3e}   GATE(<1e-11): "
          f"{'PASS' if overall_max < 1e-11 else 'FAIL'}")
    print("=" * 78)
    return overall_max


if __name__ == "__main__":
    main()
