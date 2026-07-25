"""
Finer grids for the swallowtail v5 and v6 -- the two coefficients that limit the Borel analysis.

Motivation (SWALLOWTAIL_TRANSSERIES_NOTES.md, "CAMPAIGN RESULT"): with grids n=16,20,24 (v5) and
n=10,12,14 (v6), v6 is pinned only to a factor ~2.5 and v5 to ~25%.  Both sequences are still
rising steeply, so the n->inf extrapolation -- not the number of rungs -- is what blocks the
complex-pair phase and the median-Borel resummation.  This extends both grids.

Parallel, unlike swtl_campaign.py: every _ybase_q3_*.pkl (idx 1..13) is already built, so the jobs
share no writable symbolic state, and the per-job dense caches are keyed _ckpt_{k}_{n}_*.pkl.  The
two tracks use DISJOINT n so no two concurrent jobs can target the same cache file.  Each job is
single-threaded (~100% of one core, ~0.5 GB), so 6 at once fits comfortably on 10 cores / 17 GB.

Cost estimate from the observed n^~4.8 scaling (v6: 297s at n=12 -> 627s at n=14):
    v5  n=28 ~14min, n=32 ~27min, n=36 ~47min
    v6  n=16 ~20min, n=18 ~35min, n=20 ~58min
=> ~1h for the whole fan-out in parallel, vs ~3.5h sequential.

Results are merged into swtl_results.json by THIS process only (children never write it), so the
existing checkpoints cannot be corrupted by a race.  Log: swtl_fine.log

Run:  python3 swtl_fine.py
"""
import subprocess, json, re, os, sys, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, 'swtl_results.json')

JOBS = [(5, 28), (5, 32), (5, 36), (6, 16), (6, 18), (6, 20)]
MAXPAR = 6


def run_job(job):
    k, n = job
    t0 = time.time()
    p = subprocess.run([sys.executable, os.path.join(HERE, 'swtl_production.py'), '3', str(k), str(n)],
                       capture_output=True, text=True, cwd=HERE)
    out = p.stdout + p.stderr
    m = re.search(rf'v{k}\(q=3, n={n}\)\s*=\s*([-+0-9.eE]+)', out)
    dt = time.time() - t0
    if not m:
        print(f"  v{k}(n={n}) FAILED [{dt/60:.1f}min]\n" + "\n".join(out.splitlines()[-6:]), flush=True)
        return k, n, None, dt
    val = float(m.group(1))
    print(f"  v{k}(n={n}) = {val:+.6f}  [{dt/60:.1f}min]", flush=True)
    return k, n, val, dt


def main():
    print(f"=== swallowtail FINE grids: {JOBS} ({MAXPAR} in parallel) ===", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=MAXPAR) as ex:
        results = list(ex.map(run_job, JOBS))
    # single-writer merge
    res = json.load(open(RESULTS)) if os.path.exists(RESULTS) else {}
    for k, n, val, dt in results:
        if val is None:
            continue
        res.setdefault(str(k), {})[str(n)] = val
    json.dump(res, open(RESULTS, 'w'), indent=2)
    print(f"\n=== FINE GRIDS DONE in {(time.time()-t0)/60:.1f} min; merged into swtl_results.json ===",
          flush=True)
    for k in ('5', '6'):
        if k in res:
            got = sorted((int(n), v) for n, v in res[k].items() if n != 'extrap_ninf')
            print(f"  v{k}: {[(n, round(v,6)) for n, v in got]}", flush=True)
    print("\nNow rerun:  python3 swtl_borel.py", flush=True)


if __name__ == "__main__":
    main()
