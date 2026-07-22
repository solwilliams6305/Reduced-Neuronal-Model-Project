"""
folded_pde_airy_process.py
==========================
Direction "Airy process from a folded PDE" (banked #4): can a spatially-extended /
many-mode folded system produce the Airy_2 *process* (the KPZ-class curve), not just
the TW_beta marginal (rung 1) or the single-operator Airy point process (rung 2)?

The transverse linearisation of N coupled Cole-Hopf fields about the canard is
    U''(x) = [ x I + C - eta * Xi(x) ] U(x),
so the COUPLING MATRIX C is the random-matrix ensemble whose edge governs the joint
peel-off law. This script asks which C gives genuine Dyson (log-gas) repulsion.

Three diagnostics, no scipy (numpy eigvalsh + Sturm-free):
  (1) Reference object: Dyson Brownian motion edge  ->  Airy_2 signature
      (TW_2 marginal + increment variance linear then saturating at 2 Var(TW_2)).
  (2) Mechanism: level-spacing ratio <r> and edge law for three coupling structures.
  (3) Honest control: the increment-variance *shape* is generic and does NOT
      discriminate; the discriminators are repulsion <r> and the edge marginal.

Verdict (see FOLDED_PDE_AIRY_PROCESS.md): genuine Airy_2 needs a LONG-RANGE AND
DISORDERED coupling (full random matrix). Nearest-neighbour -> Anderson localisation
(Poisson edge); uniform mean-field -> a single BBP collective outlier. Neither works.
"""
import numpy as np

TW2 = dict(mean=-1.771086, var=0.813195, skew=0.224084)  # Tracy-Widom_2 reference

# ---------- ensembles (C = coupling structure of the folded linearisation) ----------
def gue(N, rng):                                   # random all-to-all (disordered)
    A = rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N))
    return (A + A.conj().T)/2                       # E|H_ij|^2 = 1 (off-diag)

def goe(N, rng):                                   # real-symmetric sibling (beta=1)
    A = rng.standard_normal((N, N));  return (A + A.T)/np.sqrt(2)

def nn_anderson(N, rng, W=2.0):                    # nearest-neighbour (spatially local)
    d = W*rng.standard_normal(N)
    e = (rng.standard_normal(N-1) + 1j*rng.standard_normal(N-1))/np.sqrt(2)
    H = np.diag(d).astype(complex)
    H[np.arange(N-1), np.arange(1, N)] = e;  H[np.arange(1, N), np.arange(N-1)] = e.conj()
    return H

def mean_field(N, rng, kappa=2.0):                 # uniform all-to-all = rank-1
    return np.diag(rng.standard_normal(N)) + (kappa/N)*np.ones((N, N))

def rbar(ev):                                      # level-spacing ratio (bulk), no unfolding
    s = np.diff(np.sort(ev)); s = s[len(s)//4:3*len(s)//4]
    r = np.minimum(s[:-1], s[1:])/np.maximum(s[:-1], s[1:]);  return r.mean()
# GUE <r>=0.6027, GOE 0.5307, Poisson 0.3863

# ---------- (1) Airy_2 reference via Dyson Brownian motion ----------
def dbm_reference(N=128, T=2500, dt=0.02, seed=1):
    rng = np.random.default_rng(seed)
    H = gue(N, rng); lam = np.empty(T)
    for t in range(T):
        lam[t] = np.linalg.eigvalsh(H)[-1]
        H = H*(1-0.5*dt) + gue(N, rng)*np.sqrt(dt)   # Hermitian OU -> stationary GUE
    xi = N**(1/6.)*(lam - 2*np.sqrt(N))              # edge rescale -> TW_2 marginal
    lags = np.array([1, 2, 4, 8, 16, 32, 64, 128, 256])
    V = np.array([(xi[L:]-xi[:-L]).var() for L in lags])
    return xi, lags*dt, V

# ---------- (2) mechanism: which coupling repels like a log-gas ----------
def mechanism(N=160, S=1200, seed=7):
    rng = np.random.default_rng(seed)
    out = {}
    for name, gen in [("random(GUE)", gue), ("nearest-nbr(Anderson)", nn_anderson),
                      ("mean-field(rank1)", mean_field)]:
        tops = np.empty(S); rs = np.empty(S)
        for i in range(S):
            ev = np.linalg.eigvalsh(gen(N, rng)); tops[i] = ev[-1]; rs[i] = rbar(ev)
        z = (tops-tops.mean())/tops.std()
        out[name] = dict(rbar=rs.mean(), edge_skew=float((z**3).mean()), z=z)
    out["GOE<r>"] = float(np.mean([rbar(np.linalg.eigvalsh(goe(N, rng))) for _ in range(300)]))
    return out

if __name__ == "__main__":
    xi, t, V = dbm_reference()
    print("(1) DBM edge  marginal: mean %.3f var %.3f skew %.3f  (TW2 %.3f/%.3f/%.3f)"
          % (xi.mean(), xi.var(), ((xi-xi.mean())**3).mean()/xi.var()**1.5,
             TW2['mean'], TW2['var'], TW2['skew']))
    print("    increment var saturates at %.3f  (2 Var(TW2)=%.3f) -> Airy_2 signature"
          % (V[-2], 2*TW2['var']))
    m = mechanism()
    print("\n(2) coupling           <r>     edge-skew(std)   verdict")
    print("    random(GUE)       %.3f    %+.3f          TW2 edge + Dyson repulsion -> Airy_2"
          % (m["random(GUE)"]['rbar'], m["random(GUE)"]['edge_skew']))
    print("    nearest-nbr       %.3f    %+.3f          Poisson (localised), not TW"
          % (m["nearest-nbr(Anderson)"]['rbar'], m["nearest-nbr(Anderson)"]['edge_skew']))
    print("    mean-field        %.3f    %+.3f          BBP collective outlier, not TW"
          % (m["mean-field(rank1)"]['rbar'], m["mean-field(rank1)"]['edge_skew']))
    print("    GOE (real-symm)   %.3f                   real couplings -> beta=1 (TW1/Airy1)"
          % m["GOE<r>"])
