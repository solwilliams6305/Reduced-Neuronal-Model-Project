"""Diagnostic: profile the Var(Y8) moment set at n=10. Time every distinct moment via the
production engine; log the slow ones by structural class (chain-orders, #legs, boundary), so we
find whatever class the numba/quad fixes did NOT cover. Prints running totals so a hang is visible."""
import os, sys, time, collections
os.environ.setdefault('CT_NUMBA', '1')
import chaos_diagram as CD
import chaos_transfer as CT
from _v6_driver import _yfile

n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
CT.clear_all()
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us, s, Vs, Ye = CD.load_Y(fname=_yfile(15))
Y8 = CD.Ybase(8, Us, s, Vs, Ye, u0ps, Vst)
print(f"Y8: {len(Y8)} base monomials at n={n}", flush=True)

def cls(atoms):
    chains = tuple(sorted(a[1] for a in atoms if a[0] == 'U' and a[1] >= 2))
    legs = sum(1 for a in atoms if a[0] == 'U' and a[1] == 1)
    bnd = sum(1 for a in atoms if a[0] == 's')
    return (chains, legs, bnd)

# enumerate distinct moments in <Y8 Y8>
seen = {}
for ca, aa in Y8:
    for cb, ab in Y8:
        if ca == 0.0 or cb == 0.0:
            continue
        key = tuple(sorted(aa + ab))
        if key not in seen:
            seen[key] = key
print(f"distinct moments in <Y8 Y8>: {len(seen)}", flush=True)

slow = []
class_time = collections.defaultdict(float)
class_count = collections.defaultdict(int)
t0 = time.time(); done = 0
for key in seen:
    atoms = list(key)
    t1 = time.time()
    CT.moment_hybrid(atoms)
    dt = time.time() - t1
    c = cls(atoms)
    class_time[c] += dt; class_count[c] += 1
    done += 1
    if dt > 1.0:
        slow.append((dt, c))
        print(f"  SLOW {dt:6.2f}s  chains={c[0]} legs={c[1]} bnd={c[2]}", flush=True)
    if done % 2000 == 0:
        print(f"  ... {done}/{len(seen)} moments, {time.time()-t0:.0f}s elapsed", flush=True)

print(f"\nTOTAL {done} moments in {time.time()-t0:.0f}s", flush=True)
print("Top classes by total time:", flush=True)
for c, tt in sorted(class_time.items(), key=lambda x: -x[1])[:15]:
    print(f"  {tt:8.1f}s  x{class_count[c]:5d}  chains={c[0]} legs={c[1]} bnd={c[2]}", flush=True)
