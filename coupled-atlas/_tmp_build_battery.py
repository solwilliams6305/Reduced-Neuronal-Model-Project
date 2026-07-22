"""Build a REFERENCE BATTERY from the trusted moment engine (chaos_transfer, pure-python DP;
NO numba). Computes/pickles reference moment values spanning the DP cost classes, plus ladder-level
compute_vk references, plus known full-engine v6/v7 anchors as metadata.

Output pickle: /Users/solomonwilliams/Reduced Neuronal Model Project/coupled-atlas/_tmp_battery.pkl
  { 'battery' : [ (n_grid, atoms_tuple, value), ... ],   # atoms in engine format
    'ladders' : { ('v3',10): val, ('v3',12): val, ('v4',10): val, ('v4',12): val },
    'anchors' : { ('v6',10): +80.941972, ('v6',12): -1.73703, ('v7',10): -178.60840 },
    'meta'    : { ... } }
"""
import os, sys, time, pickle, platform

# Guarantee we do NOT engage any numba path.
os.environ.pop('CT_NUMBA', None)

import chaos_diagram as CD
import chaos_transfer as CT
import _v6_driver as drv

OUT = "/Users/solomonwilliams/Reduced Neuronal Model Project/coupled-atlas/_tmp_battery.pkl"

# ---- helpers for building leg lists ----
def legs(n0, n1=0):
    return [('U', 1, 0)] * n0 + [('U', 1, 1)] * n1

# ---------------------------------------------------------------------------
# MOMENT SPECS.  Each entry: (atoms, grids) where grids is the list of n_grid
# to evaluate at.  Default grids = both 8 and 12.
# ---------------------------------------------------------------------------
BOTH = [8, 12]

PURE = [
    [('U', 2, 0), ('U', 2, 0)],                                  # 2 chains ord2, m=0
    [('U', 2, 1), ('U', 2, 1)],                                  # 2 chains ord2, m=1
    [('U', 2, 0), ('U', 2, 1)],                                  # mixed heads
    [('U', 3, 0), ('U', 3, 0)],                                  # 2 chains ord3
    [('U', 3, 0), ('U', 3, 1)],                                  # 2 chains ord3 mixed
    [('U', 4, 0), ('U', 4, 0)],                                  # 2 chains ord4
    [('U', 4, 0), ('U', 2, 0)],                                  # ord4 x ord2
    [('U', 5, 0), ('U', 5, 0)],                                  # 2 chains ord5
    [('U', 5, 0), ('U', 3, 0)],                                  # ord5 x ord3
    [('U', 6, 0), ('U', 6, 0)],                                  # 2 chains ord6
    [('U', 6, 1), ('U', 6, 0)],                                  # ord6 mixed heads
    [('U', 7, 0), ('U', 7, 0)],                                  # order-7 pair (v7)
    [('U', 8, 0), ('U', 8, 0)],                                  # order-8 pair (v7)
    [('U', 2, 0), ('U', 2, 0), ('U', 2, 0)],                     # 3 chains ord2
    [('U', 3, 0), ('U', 3, 0), ('U', 2, 0)],                     # 3 chains (3,3,2)
    [('U', 3, 1), ('U', 3, 0), ('U', 2, 1)],                     # 3 chains mixed heads
    [('U', 4, 0), ('U', 3, 0), ('U', 3, 0)],                     # 3 chains (4,3,3)
    [('U', 2, 0), ('U', 2, 0), ('U', 2, 0), ('U', 2, 0)],        # 4 chains ord2
    [('U', 3, 0), ('U', 2, 0), ('U', 2, 0), ('U', 2, 0)],        # 4 chains (3,2,2,2)
    [('U', 2, 1), ('U', 2, 0), ('U', 2, 1), ('U', 2, 0)],        # 4 chains mixed heads
]

CHAINS_LEGS = [
    [('U', 2, 0)] + legs(2),                                     # 1 chain + 2 legs
    [('U', 2, 0)] + legs(1, 1),                                  # 1 chain + 2 legs mixed
    [('U', 3, 0)] + legs(3),                                     # 1 chain + 3 legs
    [('U', 4, 0)] + legs(1, 1),                                  # 1 chain ord4 + 2 legs
    [('U', 4, 0)] + legs(2, 2),                                  # 1 chain ord4 + 4 legs
    [('U', 3, 0)] + legs(3, 2),                                  # 1 chain + 5 legs mixed
    [('U', 2, 0), ('U', 2, 0)] + legs(2),                        # 2 chains + 2 legs
    [('U', 2, 0), ('U', 2, 1)] + legs(2, 2),                     # 2 chains + 4 legs mixed
    [('U', 3, 0), ('U', 2, 0)] + legs(3),                        # 2 chains + 3 legs
    [('U', 3, 0), ('U', 3, 1)] + legs(3, 3),                     # 2 chains ord3 + 6 legs mixed
    [('U', 4, 0), ('U', 3, 0)] + legs(3),                        # 2 chains + 3 legs (tot 10)
    [('U', 2, 0), ('U', 2, 0), ('U', 2, 0)] + legs(2),           # 3 chains + 2 legs
    [('U', 2, 0)] * 3 + legs(4, 4),                              # 3 chains ord2 + 8 legs (v6 class)
]

# heavy DP-class entries (give explicit grids so the 42s n=12 case runs once)
HEAVY = [
    ([('U', 3, 0), ('U', 3, 0), ('U', 2, 0)] + legs(4, 4), BOTH),   # v6 bottleneck (3,3,2)+8 legs
    ([('U', 3, 0), ('U', 2, 0), ('U', 2, 0), ('U', 2, 0)] + legs(3, 2), BOTH),  # 4-chain (3,2,2,2)+5 legs (nonzero)
]

BOUNDARY = [
    [('U', 3, 0), ('U', 2, 0), ('s', 0)],                        # 2 chains + 1 s
    [('U', 2, 0), ('U', 1, 0), ('s', 0)],                        # 1 chain + 1 leg + 1 s
    [('U', 2, 0), ('U', 2, 0), ('U', 1, 0), ('s', 0)],           # 2 chains + 1 leg + 1 s
    [('U', 2, 0), ('U', 2, 0), ('s', 0), ('s', 1)],              # 2 chains + 2 s
    [('U', 3, 0), ('U', 3, 0), ('s', 0), ('s', 1)],              # 2 chains ord3 + 2 s
    [('U', 2, 0), ('U', 2, 0)] + legs(1, 1) + [('s', 0), ('s', 1)],  # 2 chains + 2 legs + 2 s
    [('U', 3, 0), ('U', 2, 0), ('U', 2, 1), ('s', 0), ('U', 1, 0), ('U', 1, 1)],  # test_sep, 1 s
    [('U', 4, 0), ('U', 3, 1), ('U', 2, 0), ('s', 1), ('s', 0), ('U', 1, 0)],     # test_sep, 2 s
]

# parity-zero / connectivity-zero (expect exactly 0.0)
ZERO = [
    [('U', 2, 0)],                                               # lone chain: connectivity zero
    [('U', 3, 0), ('U', 2, 0)],                                  # bulk tot 5 odd: parity zero
    [('U', 2, 0), ('U', 1, 0)],                                  # tot 3 odd: parity zero
    [('U', 5, 0), ('U', 2, 0)],                                  # order 5 > 2: connectivity zero
    [('U', 4, 0), ('U', 1, 0)],                                  # order 4 > 1: connectivity zero
    [('s', 0)],                                                  # boundary only: zero
    [('s', 0), ('s', 1)],                                        # boundary only: zero
    [('U', 3, 0), ('U', 2, 0), ('U', 2, 0), ('U', 2, 0)] + legs(3, 3),  # 4-chain (3,2,2,2)+6 legs: tot 15 odd -> 0
]

# assemble (atoms, grids) work items
work = []
for a in PURE + CHAINS_LEGS + BOUNDARY + ZERO:
    work.append((a, BOTH))
for a, g in HEAVY:
    work.append((a, g))

# ---------------------------------------------------------------------------
# compute moments
# ---------------------------------------------------------------------------
CT.QUAD_ENABLE = True      # trusted default fast path (validated to <1e-10 vs the DP)
CT.SEP_ENABLE = False      # trusted default

battery = []
timings = []   # (dt, n, atoms_tuple)
cur_n = None
t_start = time.time()
print("=== moment battery ===", flush=True)
for atoms, grids in work:
    for n in grids:
        if n != cur_n:
            CD.setup(n_grid=n, MAXORD=2)
            cur_n = n
        CT.clear_all()
        t0 = time.time()
        v = CT.moment_transfer(atoms)
        dt = time.time() - t0
        battery.append((n, tuple(atoms), v))
        timings.append((dt, n, tuple(atoms)))
        flag = "  <<slow" if dt > 3.0 else ""
        print(f"  n={n:2d} [{dt:7.2f}s] {v:+.8e}  {atoms}{flag}", flush=True)

print(f"\n[moments done: {len(battery)} entries in {time.time()-t_start:.0f}s]", flush=True)

# ---------------------------------------------------------------------------
# ladder-level references (full compute_vk with the trusted transfer engine)
# ---------------------------------------------------------------------------
print("\n=== ladder references (compute_vk, trusted engine) ===", flush=True)
ladders = {}
for k, n in [(3, 10), (3, 12), (4, 10), (4, 12)]:
    t0 = time.time()
    v = drv.compute_vk(k, n, CT.moment_transfer)
    dt = time.time() - t0
    ladders[(f"v{k}", n)] = v
    print(f"  v{k}(n={n}) = {v:+.6f}   [{dt:.1f}s]", flush=True)

# known full-engine anchors (recorded, not recomputed here)
anchors = {
    ("v6", 10): +80.941972,
    ("v6", 12): -1.73703,
    ("v7", 10): -178.60840,
}

meta = {
    "built": time.strftime("%Y-%m-%d %H:%M:%S"),
    "engine": "chaos_transfer (pure-python DP, QUAD_ENABLE=True, SEP_ENABLE=False, no numba)",
    "python": platform.python_version(),
    "CT_NUMBA_env": os.environ.get("CT_NUMBA"),
    "DENSE_MAXORD": CT.DENSE_MAXORD,
    "note": "battery values are trusted-engine references for validating a faster DP; "
            "atoms format ('U',k,m) bulk / ('s',j) boundary; anchors are known full-engine values.",
}

out = {"battery": battery, "ladders": ladders, "anchors": anchors, "meta": meta}
with open(OUT, "wb") as f:
    pickle.dump(out, f)

# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------
print(f"\n=== wrote {OUT} ===", flush=True)
print(f"battery entries : {len(battery)}", flush=True)
print(f"ladders         : {ladders}", flush=True)
print(f"anchors (meta)  : {anchors}", flush=True)
timings.sort(reverse=True)
print("\nslowest 6 moment entries:", flush=True)
for dt, n, atoms in timings[:6]:
    print(f"  [{dt:7.2f}s] n={n} {atoms}", flush=True)
print(f"\ntotal wall time: {time.time()-t_start:.0f}s", flush=True)
