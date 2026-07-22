"""One-time build of _ybase_14.pkl and _ybase_15.pkl (base-atom decompositions of Y14, Y15).
Grid-independent; any small grid works for the setup constants."""
import time
import chaos_diagram as CD

t0 = time.time()
CD.setup(n_grid=10, MAXORD=2)
u0ps = CD.GEO['u0ps']; Vst = CD.GEO['Vst']
Us, s, Vs, Ye = CD.load_Y('_yexprs_15.txt')
for idx in [14, 15]:
    terms = CD.Ybase(idx, Us, s, Vs, Ye, u0ps, Vst)
    print(f'Y{idx}: {len(terms)} base monomials [{time.time()-t0:.0f}s]', flush=True)
