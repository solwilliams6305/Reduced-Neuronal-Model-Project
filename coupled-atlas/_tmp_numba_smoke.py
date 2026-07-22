"""Smoke test for chaos_dp_numba: validate the numba DP kernel against the
pure-python dict-DP in chaos_transfer._moment_bulk_gen, and time the slow class.

Run:  python _tmp_numba_smoke.py
"""
import time
import numpy as np

import chaos_diagram as CD
import chaos_transfer as CT
from chaos_dp_numba import moment_bulk_numba


def py_dp(gatoms):
    """Force the python dict-DP: disable QUAD/SEP so order<=2 cases also route
    through the DP for an apples-to-apples comparison.  (Two-chain-no-leg cases
    still hit the validated _moment_two fast path, which is a fine reference.)"""
    old_q, old_s, old_nb = CT.QUAD_ENABLE, CT.SEP_ENABLE, CT.NUMBA_DP
    CT.QUAD_ENABLE = False
    CT.SEP_ENABLE = False
    CT.NUMBA_DP = False
    try:
        return CT._moment_bulk_gen(gatoms)
    finally:
        CT.QUAD_ENABLE, CT.SEP_ENABLE, CT.NUMBA_DP = old_q, old_s, old_nb


def relerr(a, b):
    d = abs(a - b)
    s = max(abs(a), abs(b))
    return d / s if s > 0 else d


def main():
    CD.setup(n_grid=16, MAXORD=5)
    G = CD.GEO
    n = G['n']
    Gnode = G['Gnode']
    H = G['H']
    # a couple of extra distinct head vectors to make >2 leg types / distinct chains
    rng = np.random.default_rng(0)
    Hb = G['PSI'].copy()
    Hc = (G['U'] * 0.3 + G['PSI'] * 0.7).copy()

    def chain(o, which):
        h = {'g': Gnode, 'h': H, 'b': Hb, 'c': Hc}[which]
        return (o, h)

    def leg(which):
        h = {'g': Gnode, 'h': H, 'b': Hb, 'c': Hc}[which]
        return (1, h)

    cases = []
    # (1) 3 chains mixed orders (3,3,2) + 4 legs of 2 types
    cases.append(("3 chains (3,3,2) + 4 legs/2 types",
                  [chain(3, 'g'), chain(3, 'h'), chain(2, 'b'),
                   leg('g'), leg('g'), leg('h'), leg('h')]))
    # (2) (3,2,2,2) + legs
    cases.append(("4 chains (3,2,2,2) + 3 legs/2 types",
                  [chain(3, 'g'), chain(2, 'h'), chain(2, 'b'), chain(2, 'c'),
                   leg('g'), leg('h'), leg('h')]))
    # (3) pure chains (3,3,2)  -> pure-chain DP path
    cases.append(("pure chains (3,3,2)",
                  [chain(3, 'g'), chain(3, 'h'), chain(2, 'b')]))
    # (4) single chain + legs
    cases.append(("single chain (3) + 3 legs/2 types",
                  [chain(3, 'g'), leg('h'), leg('h'), leg('b')]))
    # (5) order-7 pair + extras (2 chains order 7 + a chain order2 + 2 legs)
    cases.append(("order-7 pair + chain(2) + 2 legs",
                  [chain(7, 'g'), chain(7, 'h'), chain(2, 'b'), leg('g'), leg('h')]))
    # (6) parity-zero case (odd total): wrapper must return 0.0
    cases.append(("parity-zero (3,2) total odd",
                  [chain(3, 'g'), chain(2, 'h')]))
    # (7) pure chains (4,4) equal order -> two-chain fast path reference
    cases.append(("pure chains (4,4)",
                  [chain(4, 'g'), chain(4, 'h')]))
    # (8) legs-only bath (no chains): 4 legs / 2 types
    cases.append(("legs only: 4 legs / 2 types",
                  [leg('g'), leg('g'), leg('h'), leg('h')]))
    # (9) 3 chains (3,3,2) + 6 legs / 3 types  (stress leg-type count)
    cases.append(("3 chains (3,3,2) + 6 legs/3 types",
                  [chain(3, 'g'), chain(3, 'h'), chain(2, 'b'),
                   leg('g'), leg('g'), leg('h'), leg('h'), leg('b'), leg('b')]))

    print(f"grid n={n}\n")
    print(f"{'case':44s} {'py_dp':>16s} {'numba':>16s} {'rel_err':>10s}")
    worst = 0.0
    for name, at in cases:
        vp = py_dp(at)
        vn = moment_bulk_numba(at)
        re = relerr(vp, vn)
        worst = max(worst, re)
        print(f"{name:44s} {vp:16.9e} {vn:16.9e} {re:10.2e}")
    print(f"\nworst rel err over spot-checks: {worst:.3e}")

    # ---- timing: (3,3,2) + 8 legs at n=16 ----
    timing_atoms = [chain(3, 'g'), chain(3, 'h'), chain(2, 'b')] + \
                   [leg('g')] * 4 + [leg('h')] * 4
    # correctness of the timed class first
    vp = py_dp(timing_atoms)
    vn = moment_bulk_numba(timing_atoms)
    print(f"\n(3,3,2)+8 legs/2 types  py={vp:.9e}  numba={vn:.9e}  rel={relerr(vp,vn):.2e}")

    # warm the numba compile (already compiled from spot-checks, but be safe)
    moment_bulk_numba(timing_atoms)

    REP = 5
    t0 = time.perf_counter()
    for _ in range(REP):
        py_dp(timing_atoms)
    t_py = (time.perf_counter() - t0) / REP
    t0 = time.perf_counter()
    for _ in range(REP):
        moment_bulk_numba(timing_atoms)
    t_nb = (time.perf_counter() - t0) / REP
    print(f"\ntiming (3,3,2)+8 legs, n=16, avg of {REP}:")
    print(f"  python DP : {t_py*1e3:8.2f} ms")
    print(f"  numba DP  : {t_nb*1e3:8.2f} ms   speedup x{t_py/t_nb:.1f}")

    # ---- routing check: CT_NUMBA unset -> behavior unchanged ----
    assert CT.NUMBA_DP is False, "NUMBA_DP should default False when CT_NUMBA unset/0"
    print("\nrouting: NUMBA_DP defaults to", CT.NUMBA_DP, "(CT_NUMBA unset) -> python DP unchanged")


if __name__ == "__main__":
    main()
