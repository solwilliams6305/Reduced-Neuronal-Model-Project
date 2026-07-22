"""Mean ladder  m_j = <Y_{2j}>  for  <Y*> = -2.188 + sum_{j>=1} m_j eta^{2j}.
These are SINGLE moments (no products) -> cheap.  Computes m_1..m_6 = <Y_2>,<Y_4>,
<Y_6>,<Y_8>,<Y_10>,<Y_12> at grid n (arg), parallel over the distinct moments.
Usage: python3 _w_mean_ladder.py <n> [nproc]
Anchor check: m_1 = <Y_2> should be ~ +0.212.
"""
import os, sys, time
os.environ.setdefault('CT_NUMBA','1')
import multiprocessing as mp
import chaos_diagram as CD, chaos_transfer as CT
from _v6_driver import _yfile

def _init(n):
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
def _work(keys):
    return [(k, CT.moment_hybrid(list(k))) for k in keys]

def main():
    n = int(sys.argv[1]) if len(sys.argv)>1 else 16
    nproc = int(sys.argv[2]) if len(sys.argv)>2 else os.cpu_count()
    t0=time.time()
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    u0ps=CD.GEO['u0ps']; Vst=CD.GEO['Vst']
    Us,s,Vs,Ye=CD.load_Y(fname=_yfile(15))
    ks=[2,4,6,8,10,12]
    Y={k: CD.Ybase(k,Us,s,Vs,Ye,u0ps,Vst) for k in ks}
    print(f"[{time.time()-t0:.0f}s] Y built, sizes "+", ".join(f"|Y{k}|={len(Y[k])}" for k in ks), flush=True)

    def mkey(atoms): return tuple(sorted(atoms))
    keyset=set()
    for k in ks:
        for c,at in Y[k]:
            if c!=0.0: keyset.add(mkey(at))
    keys=list(keyset)
    print(f"[{time.time()-t0:.0f}s] {len(keys)} distinct moments, nproc={nproc}", flush=True)

    M={}
    def chunks(lst,m):
        out=[[] for _ in range(m)]
        for i,x in enumerate(lst): out[i%m].append(x)
        return [c for c in out if c]
    with mp.Pool(nproc, initializer=_init, initargs=(n,)) as pool:
        for res in pool.imap_unordered(_work, chunks(keys, nproc*6)):
            for k,v in res: M[k]=v
    print(f"[{time.time()-t0:.0f}s] moments done", flush=True)

    print(f"\n mean ladder at n={n}:")
    for j,k in enumerate(ks, start=1):
        mj=sum(c*M[mkey(at)] for c,at in Y[k] if c!=0.0)
        print(f"  m_{j} = <Y_{k}> = {mj:+.5f}", flush=True)

if __name__=="__main__":
    main()
