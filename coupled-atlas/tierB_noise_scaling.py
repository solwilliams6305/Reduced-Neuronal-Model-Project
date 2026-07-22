import numpy as np
# Robust collapse test: evolve swept pitchfork to a FIXED inner time s (mu=sqrt(eps)*s), measure Var(u)/sqrt(eps)
# = Var(u_bar). Exact rescaling => depends on (sigma,eps) ONLY via eta=sigma_u/sqrt(eps) (beta=8eps/sigma^2).
def var_ubar_at_s(sigma, eps, s_final, N, M=1.0, seed=0):
    sig_u=sigma/np.sqrt(2.0); dt=0.05*np.sqrt(eps)
    rng=np.random.default_rng(seed)
    t_final=M/eps + s_final/np.sqrt(eps); n=int(t_final/dt); sq=np.sqrt(dt)
    u=np.zeros(N)
    for i in range(n):
        mu=-M+eps*(i*dt); u+=(mu*u-u**3)*dt+sig_u*sq*rng.standard_normal(N)
        np.clip(u,-5,5,out=u)
    return np.var(u)/np.sqrt(eps), np.mean(np.abs(u))/eps**0.25
for s in [-0.5, 0.0, 0.5]:
    print(f"--- inner time s={s} (mu={s}*sqrt(eps)) ---")
    for name,(sig,eps) in {'A1 e=.04 s=.4':(0.4,0.04),'A2 e=.01 s=.2':(0.2,0.01),'B e=.04 s=.2 (b=8)':(0.2,0.04)}.items():
        beta=8*eps/sig**2; v,ma=var_ubar_at_s(sig,eps,s,8000,seed=1)
        print(f"   {name:18} beta={beta:.1f}  Var(u_bar)={v:.4f}  mean|u_bar|={ma:.4f}")
print("=> A1,A2 (same beta=2) should match Var(u_bar); B (beta=8, weaker noise) smaller")
