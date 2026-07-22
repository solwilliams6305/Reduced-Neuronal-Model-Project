"""Term-wise grid extrapolation for v_k: parse per-term driver logs, fit each (a,b) term's
n-sequence, assemble v_k(inf) = sum of extrapolated terms with a fit-spread error bar.

The assembled v6(n) residue oscillates (terms ~+-40 cancel to ~+-2) so fitting the residue is
ill-conditioned; the individual terms converge smoothly and are the right objects to extrapolate.

Models per term (fit on the smooth branch, default n>=12):
  R2:  v + a/n + b/n^2         (algebraic boundary error)
  GEO: v + A*q^n               (geometric transient), fitted by log-diff regression
Report both; the spread is the systematic error estimate.

Usage: python3 _extrap_terms.py <log1> [<log2> ...]   (logs concatenated; each '--- vK n=N ---'
block contributes that grid's per-term values; single-grid logs from _v6_driver runs also parse
via the filename convention v6q_nNN.log)
"""
import re, sys
import numpy as np

def parse_logs(paths):
    data = {}   # (a,b) -> {n: val};  'total' -> {n: val}
    for p in paths:
        txt = open(p).read()
        n_from_name = None
        m = re.search(r'_n(\d+)\.log', p)
        if m:
            n_from_name = int(m.group(1))
        cur_n = n_from_name
        for line in txt.splitlines():
            mh = re.match(r'--- v\d+ n=(\d+) ---', line)
            if mh:
                cur_n = int(mh.group(1)); continue
            mt = re.match(r'\s*\((\d+),(\d+)\) x(\d+): ([+-][\d.]+)', line)
            if mt and cur_n is not None:
                a, b, coef, val = int(mt.group(1)), int(mt.group(2)), int(mt.group(3)), float(mt.group(4))
                data.setdefault((a, b, coef), {})[cur_n] = val
            mv = re.match(r'.*v(\d+)\(n=(\d+)\) = ([+-][\d.]+)', line)
            if mv:
                data.setdefault('total', {})[int(mv.group(2))] = float(mv.group(3))
    return data

def fit_r2(ns, ys):
    A = np.vstack([np.ones_like(ns), 1.0/ns, 1.0/ns**2]).T
    c, *_ = np.linalg.lstsq(A, ys, rcond=None)
    resid = A @ c - ys
    return c[0], np.max(np.abs(resid))

def fit_geo(ns, ys):
    """v + A q^n via diffs: d_i = y_{i+1}-y_i = A q^{n_i}(q^h - 1) with uniform h -> q from ratio."""
    d = np.diff(ys)
    if len(d) < 2 or np.any(d[:-1] == 0) or np.any(np.sign(d[1:]) != np.sign(d[:-1])):
        return None, None
    r = d[1:] / d[:-1]
    if np.any(r <= 0) or np.any(r >= 1.05):
        return None, None
    h = ns[1] - ns[0]
    q = np.median(r) ** (1.0 / h)
    tail = d[-1] * (r.mean() / (1 - r.mean())) if r.mean() < 1 else 0.0
    return ys[-1] + tail, abs(d[-1])

def report(data, nmin=12, label=''):
    tot_r2 = 0.0; tot_geo = 0.0; ok_geo = True
    print(f"--- term-wise extrapolation {label} (fit n>={nmin}) ---")
    for key in sorted(k for k in data if k != 'total'):
        a, b, coef = key
        pts = sorted((n, v) for n, v in data[key].items() if n >= nmin)
        ns = np.array([p[0] for p in pts], float); ys = np.array([p[1] for p in pts], float)
        if len(ns) < 3:
            print(f"  ({a},{b}) x{coef}: <3 pts, skipped"); continue
        vr, rr = fit_r2(ns, ys)
        vg, dg = fit_geo(ns, ys)
        tot_r2 += coef * vr
        if vg is None:
            ok_geo = False
            print(f"  ({a},{b}) x{coef}: R2={vr:+.4f} (resid {rr:.3f})  GEO=n/a   last={ys[-1]:+.4f}")
        else:
            tot_geo += coef * vg
            print(f"  ({a},{b}) x{coef}: R2={vr:+.4f} (resid {rr:.3f})  GEO={vg:+.4f}  last={ys[-1]:+.4f}")
    print(f"  assembled: R2 -> {tot_r2:+.4f}" + (f" | GEO -> {tot_geo:+.4f}" if ok_geo else " | GEO incomplete"))
    if 'total' in data:
        pts = sorted(data['total'].items())
        print("  raw totals: " + "  ".join(f"n={n}:{v:+.4f}" for n, v in pts))
    return tot_r2, (tot_geo if ok_geo else None)

if __name__ == "__main__":
    data = parse_logs(sys.argv[1:])
    report(data)
