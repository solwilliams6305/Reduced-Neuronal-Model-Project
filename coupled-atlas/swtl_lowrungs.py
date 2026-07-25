"""
Swallowtail q=3 LOW rungs (v0,v1,v2) on the PRODUCTION engine.

Purpose: the notes' v0..v2 came from the weak-noise engine (weaknoise_q.py, delta->0), while the
campaign's v3..v6 come from the Wick/transfer production engine. Domb-Sykes / Borel-Pade on a
ladder spliced from two engines is meaningless unless the normalizations agree. This recomputes
v0,v1,v2 with the SAME driver as v3..v6 on the same grids, so the ladder is homogeneous -- and
v0 doubles as an absolute check, since v0 is known exactly (0.04953187).

Cheap: v0 needs Ybase idx {1}, v1 {1..3}, v2 {1..5} -- all already cached by the campaign.
Writes swtl_lowrungs.json (same shape as swtl_results.json) so the analysis can merge them.
"""
import subprocess, json, re, os, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'swtl_lowrungs.json')
SCHEDULE = [(0, [24, 30, 36]), (1, [24, 30, 36]), (2, [24, 30, 36])]


def load():
    if os.path.exists(RESULTS):
        with open(RESULTS) as f:
            return json.load(f)
    return {}


def save(d):
    with open(RESULTS, 'w') as f:
        json.dump(d, f, indent=2)


def richardson(ns, vs):
    import numpy as np
    ns = np.array(ns, float); vs = np.array(vs, float)
    A = (np.vstack([np.ones_like(ns), 1 / ns, 1 / ns ** 2]).T if len(ns) >= 3
         else np.vstack([np.ones_like(ns), 1 / ns]).T)
    c, *_ = np.linalg.lstsq(A, vs, rcond=None)
    return float(c[0])


def main():
    res = load()
    print(f"=== swallowtail v0..v2 on the production engine (-> {RESULTS}) ===", flush=True)
    for k, grids in SCHEDULE:
        kk = str(k); res.setdefault(kk, {})
        for n in grids:
            if str(n) in res[kk]:
                print(f"  v{k}(n={n}) cached = {res[kk][str(n)]:+.6f}", flush=True); continue
            t0 = time.time()
            p = subprocess.run([sys.executable, os.path.join(HERE, 'swtl_production.py'), '3', str(k), str(n)],
                               capture_output=True, text=True, cwd=HERE)
            out = p.stdout + p.stderr
            m = re.search(rf'v{k}\(q=3, n={n}\)\s*=\s*([-+0-9.eE]+)', out)
            if not m:
                print(f"  v{k}(n={n}) FAILED [{time.time()-t0:.0f}s]:\n" + "\n".join(out.splitlines()[-6:]), flush=True)
                continue
            res[kk][str(n)] = float(m.group(1)); save(res)
            print(f"  v{k}(n={n}) = {float(m.group(1)):+.6f}  [{time.time()-t0:.0f}s]", flush=True)
        got = sorted((int(n), v) for n, v in res[kk].items() if n != 'extrap_ninf')
        if len(got) >= 2:
            vinf = richardson([g[0] for g in got], [g[1] for g in got])
            res[kk]['extrap_ninf'] = vinf; save(res)
            print(f"  ==> v{k}(n->inf) = {vinf:+.6f}   from {got}", flush=True)
    print("\n=== LOW RUNGS COMPLETE ===", flush=True)
    if '0' in res and 'extrap_ninf' in res['0']:
        v0 = res['0']['extrap_ninf']; exact = 0.04953187
        print(f"  v0 = {v0:+.6f}   exact = {exact:.8f}   ratio = {v0/exact:.6f}"
              f"   {'NORMALIZATION MATCHES' if abs(v0/exact-1)<0.02 else '*** MISMATCH ***'}", flush=True)


if __name__ == "__main__":
    main()
