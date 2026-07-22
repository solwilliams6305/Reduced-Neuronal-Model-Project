"""Third-cumulant weak-noise ladder for Y*:
     kappa_3(Y*) = sum_{a,b,c>=1} eta^{a+b+c} cum3(Y_a,Y_b,Y_c)
                 = eta^4 (t_0 + t_1 eta^2 + t_2 eta^4 + ...),   t_j = sum_{a+b+c=4+2j} cum3.
cum3(a,b,c) = <Ya Yb Yc> - <Ya><Yb Yc> - <Yb><Ya Yc> - <Yc><Ya Yb> + 2<Ya><Yb><Yc>.
Parallel over the distinct moments (means, pairs, triples).  Grid-extrapolated.

Usage: python3 _w_kappa3.py <n> [jmax] [nproc]     jmax => t_0..t_jmax (m up to 4+2*jmax)
Leading check: t_0 ~ 0.057  (skew ~ 1.16 eta,  t_0 = 1.16 * v0^1.5, v0=0.134).
"""
import os, sys, time, itertools
os.environ.setdefault('CT_NUMBA','1')
import multiprocessing as mp
import chaos_diagram as CD, chaos_transfer as CT
from _v6_driver import _yfile

def _init(n):
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD); CT.clear_all()
def _work(keys):
    return [(k, CT.moment_hybrid(list(k))) for k in keys]
def chunks(lst,m):
    out=[[] for _ in range(m)]
    for i,x in enumerate(lst): out[i%m].append(x)
    return [c for c in out if c]

def triples(m):
    """multisets a<=b<=c, a+b+c=m, a>=1, with symmetric multiplicity."""
    res=[]
    for a in range(1,m//3+1):
        for b in range(a,(m-a)//2+1):
            c=m-a-b
            if c>=b:
                mult = {1:{ (a==b==c):6}}  # placeholder
                if a==b==c: mu=1
                elif a==b or b==c: mu=3
                else: mu=6
                res.append((a,b,c,mu))
    return res

def main():
    n=int(sys.argv[1]) if len(sys.argv)>1 else 16
    jmax=int(sys.argv[2]) if len(sys.argv)>2 else 3
    nproc=int(sys.argv[3]) if len(sys.argv)>3 else os.cpu_count()
    ms=[4+2*j for j in range(jmax+1)]
    t0=time.time()
    CD.setup(n_grid=n, MAXORD=CT.DENSE_MAXORD)
    u0ps=CD.GEO['u0ps']; Vst=CD.GEO['Vst']
    Us,s,Vs,Ye=CD.load_Y(fname=_yfile(15))
    need=set()
    tri_list={}
    for m in ms:
        tri_list[m]=triples(m)
        for a,b,c,mu in tri_list[m]:
            need.update([a,b,c])
    need=sorted(need)
    Y={k: CD.Ybase(k,Us,s,Vs,Ye,u0ps,Vst) for k in need}
    print(f"[{time.time()-t0:.0f}s] Y built, indices {need}", flush=True)

    def mkey(atoms): return tuple(sorted(atoms))
    # enumerate distinct moment keys: means <Ya>, pairs <Ya Yb>, triples <Ya Yb Yc>
    keyset=set()
    def add_prod(*idxs):
        pools=[Y[i] for i in idxs]
        for combo in itertools.product(*pools):
            if all(c!=0.0 for c,_ in combo):
                atoms=[]
                for _,at in combo: atoms+=at
                keyset.add(mkey(atoms))
    seen_pair=set(); seen_tri=set()
    for m in ms:
        for a,b,c,mu in tri_list[m]:
            add_prod(a,b,c)
            for (x,y) in [(a,b),(a,c),(b,c)]:
                if (x,y) not in seen_pair: seen_pair.add((x,y)); add_prod(x,y)
            for x in (a,b,c):
                add_prod(x)
    keys=list(keyset)
    print(f"[{time.time()-t0:.0f}s] {len(keys)} distinct moments, nproc={nproc}", flush=True)

    M={}
    with mp.Pool(nproc, initializer=_init, initargs=(n,)) as pool:
        done=0
        for res in pool.imap_unordered(_work, chunks(keys, nproc*8)):
            for k,v in res: M[k]=v
            done+=len(res)
    print(f"[{time.time()-t0:.0f}s] {len(M)} moments computed", flush=True)

    # assemble
    def prod(*idxs):
        pools=[Y[i] for i in idxs]; tot=0.0
        for combo in itertools.product(*pools):
            cf=1.0; atoms=[]
            ok=True
            for c,at in combo:
                if c==0.0: ok=False; break
                cf*=c; atoms+=at
            if ok: tot+=cf*M[mkey(atoms)]
        return tot
    from functools import lru_cache
    mean=lambda a: sum(c*M[mkey(at)] for c,at in Y[a] if c!=0.0)
    meanc={a:mean(a) for a in need}
    paircache={}
    def pair(a,b):
        key=(min(a,b),max(a,b))
        if key not in paircache: paircache[key]=prod(a,b)
        return paircache[key]
    def cum3(a,b,c):
        return (prod(a,b,c) - meanc[a]*pair(b,c) - meanc[b]*pair(a,c)
                - meanc[c]*pair(a,b) + 2*meanc[a]*meanc[b]*meanc[c])

    print(f"\n kappa_3 ladder at n={n}:")
    ladder=[]
    for j,m in enumerate(ms):
        tj=0.0
        for a,b,c,mu in tri_list[m]:
            tj += mu*cum3(a,b,c)
        ladder.append(tj)
        print(f"  t_{j} (m={m}) = {tj:+.6f}", flush=True)
    # leading skew check
    v0=0.134
    print(f"\n leading skew coeff = t_0/v0^1.5 = {ladder[0]/v0**1.5:+.4f}   (expected ~1.16)", flush=True)
    import json; json.dump({'n':n,'ms':ms,'t':ladder}, open(f'_w_kappa3_n{n}.json','w'))

if __name__=="__main__":
    main()
