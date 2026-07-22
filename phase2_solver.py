"""
Phase 2 — uniform-grid quasipotential solver (Freidlin-Wentzell geometric action),
an ordered-upwind / Dijkstra member of the OLIM family, used to test hypothesis H1
(resolution breakdown of a uniform mesh at a slow-fast fold).

Quasipotential convention: SDE dx=b dt+sigma dW, isotropic noise, rate~exp(-V/sigma^2),
H(x,p)=1/2|p|^2 + b.p.  Geometric (Heymann-Vanden-Eijnden) action of a path:
   S = int ( |b||x'| - b.x' ) ds        (>=0 by Cauchy-Schwarz; = V at the minimiser)
On a uniform grid this is a shortest-path (Dijkstra) problem with nonnegative edge
weights w(X->Y) = |b_mid| |Y-X| - b_mid.(Y-X),  b_mid = (b(X)+b(Y))/2.

Gradient sanity: b=-gradU  =>  S = 2 int|gradU| = 2 Delta U = V  (matches V=2U).
"""
import numpy as np, heapq, json, sys, time

# 16-neighbour stencil (king + knight): better isotropy than 8-neighbour.
OFF = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1),
       (1,2),(1,-2),(-1,2),(-1,-2),(2,1),(2,-1),(-2,1),(-2,-1)]

def dijkstra_quasipotential(BV, BW, h, src_ij):
    """BV,BW: drift components on the grid (shape (Nv,Nw)). h: mesh spacing (square).
       src_ij: (i,j) source node (attractor, V=0). Returns V array."""
    Nv, Nw = BV.shape
    V = np.full(Nv*Nw, np.inf)
    done = np.zeros(Nv*Nw, bool)
    def idx(i,j): return i*Nw + j
    s = idx(*src_ij); V[s] = 0.0
    pq = [(0.0, s)]
    offs = [(di,dj, np.hypot(di*h, dj*h), di*h, dj*h) for (di,dj) in OFF]
    while pq:
        d,u = heapq.heappop(pq)
        if done[u]: continue
        done[u] = True
        i,j = divmod(u, Nw)
        bxu, byu = BV[i,j], BW[i,j]
        for di,dj,L,dxp,dyp in offs:
            ni,nj = i+di, j+dj
            if ni<0 or ni>=Nv or nj<0 or nj>=Nw: continue
            v = idx(ni,nj)
            if done[v]: continue
            bmx = 0.5*(bxu + BV[ni,nj]); bmy = 0.5*(byu + BW[ni,nj])
            w = L*np.hypot(bmx,bmy) - (bmx*dxp + bmy*dyp)
            if w < 0: w = 0.0
            nd = d + w
            if nd < V[v]:
                V[v] = nd
                heapq.heappush(pq,(nd,v))
    return V.reshape(Nv,Nw)

def make_grid(vlo,vhi,wlo,whi,h):
    vs = np.arange(vlo, vhi+1e-12, h); ws = np.arange(wlo, whi+1e-12, h)
    return vs, ws

# ---------------------------------------------------------------- validation
def lyap_sigma1(A):
    """Solve A S + S A^T = I (2x2) via Kronecker; exact stationary cov scale."""
    n = A.shape[0]; I = np.eye(n)
    K = np.kron(I, A) + np.kron(A, I)
    S = np.linalg.solve(K, I.reshape(-1)).reshape(n,n)
    return 0.5*(S+S.T)

def validate_linear(beta, hs):
    """Linear non-gradient SDE dx=-Ax dt+sigma dW, A=[[1,beta],[0,1]].
       Exact quasipotential V(x)=0.5 x^T Sigma1^{-1} x,  A Sigma1 + Sigma1 A^T = I."""
    A = np.array([[1.0, beta],[0.0,1.0]])
    Sig1 = lyap_sigma1(A); P = np.linalg.inv(Sig1)   # V = 0.5 x^T P x
    out = []
    for h in hs:
        vs, ws = make_grid(-2,2,-2,2,h)
        VV, WW = np.meshgrid(vs, ws, indexing='ij')
        BV = -(A[0,0]*VV + A[0,1]*WW); BW = -(A[1,0]*VV + A[1,1]*WW)
        si = int(np.argmin(np.abs(vs))); sj = int(np.argmin(np.abs(ws)))
        Vnum = dijkstra_quasipotential(BV, BW, h, (si,sj))
        Vex = 0.5*(P[0,0]*VV**2 + 2*P[0,1]*VV*WW + P[1,1]*WW**2)
        # compare on an annulus away from the source singularity and away from walls
        R = np.hypot(VV-vs[si], WW-ws[sj])
        m = (R>0.5) & (R<1.5)
        rel = np.abs(Vnum[m]-Vex[m])/np.maximum(Vex[m],1e-9)
        out.append((h, float(np.median(rel)), float(np.percentile(rel,90)), vs.size, ws.size))
        print(f"  beta={beta}  h={h:.4f} grid={vs.size}x{ws.size}  "
              f"median rel.err={np.median(rel):.3%}  p90={np.percentile(rel,90):.3%}")
    return A.tolist(), out

# ---------------------------------------------------------------- FHN fold
def fhn_drift(VV, WW, eps, a):
    BV = VV - VV**3/3.0 - WW
    BW = eps*(VV - a)
    return BV, BW

def fold_geometry(a):
    """Rest state v*=a, w*=a-a^3/3 (on left branch, a<-1). Lower fold v=-1, w_f=-2/3.
       delta = w* - w_f (leading-order fold parameter). Saddle at v=-1+sqrt(delta)."""
    wstar = a - a**3/3.0; wf = -2.0/3.0; delta = wstar - wf
    v_saddle = -1.0 + np.sqrt(delta)
    dV_frozen = (8.0/3.0)*delta**1.5      # quasipotential barrier = 2*DeltaU
    return dict(a=a, vstar=a, wstar=wstar, wf=wf, delta=delta,
                v_saddle=v_saddle, dV_frozen=dV_frozen)

def fold_experiment(eps, hs, a, half_w=0.30, vlo=-2.0, vhi=-0.35):
    g = fold_geometry(a); wstar = g['wstar']; vsad = g['v_saddle']
    res = []
    for h in hs:
        t0 = time.time()
        vs, ws = make_grid(vlo, vhi, wstar-half_w, wstar+half_w, h)
        VV, WW = np.meshgrid(vs, ws, indexing='ij')
        BV, BW = fhn_drift(VV, WW, eps, a)
        si = int(np.argmin(np.abs(vs-g['vstar']))); sj = int(np.argmin(np.abs(ws-wstar)))
        Vnum = dijkstra_quasipotential(BV, BW, h, (si,sj))
        # barrier = min of V over the saddle column v=v_saddle (the threshold ridge)
        ci = int(np.argmin(np.abs(vs - vsad)))
        col = Vnum[ci,:]; B = float(np.nanmin(col[np.isfinite(col)]))
        # also min over fold line v=-1 (cross-check)
        cf = int(np.argmin(np.abs(vs - (-1.0))))
        colf = Vnum[cf,:]; Bf = float(np.nanmin(colf[np.isfinite(colf)]))
        res.append(dict(h=h, B=B, B_foldline=Bf, nv=int(vs.size), nw=int(ws.size),
                        t=round(time.time()-t0,2)))
        print(f"  eps={eps:g} h={h:.4f} grid={vs.size}x{ws.size} "
              f"B={B:.5f} (frozen {g['dV_frozen']:.5f})  {time.time()-t0:.1f}s")
    return dict(eps=eps, a=a, geom=g, runs=res)

# ---------------------------------------------------------------- CLI
if __name__ == "__main__":
    mode = sys.argv[1]
    OUT = "/sessions/confident-pensive-heisenberg/mnt/outputs"
    if mode == "validate":
        allres = {}
        for beta in [0.0, 1.0]:
            A, out = validate_linear(beta, [0.08,0.05,0.03,0.02,0.013])
            allres[str(beta)] = out
        json.dump(allres, open(f"{OUT}/phase2_validation.json","w"), indent=2)
        print("saved validation.")
    elif mode == "fold":
        eps = float(sys.argv[2]); a = float(sys.argv[3]); hs = [float(x) for x in sys.argv[4:]]
        r = fold_experiment(eps, hs, a)
        fn = f"{OUT}/fold_eps_{eps:g}.json"
        json.dump(r, open(fn,"w"), indent=2, default=float)
        print("saved", fn)
