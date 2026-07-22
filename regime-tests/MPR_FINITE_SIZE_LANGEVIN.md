# Finite-size Langevin for the QIF mean-field, and its edge reduction

*Theory track. The rigorous object under the collective two-edge result: the $O(N^{-1/2})$ system-size
Langevin correction to the Montbrió–Pazó–Roxin (MPR) mean-field, its reduction at a macroscopic
bifurcation, and the honest identification of which edge the studied network actually realises. Code:
`mpr_finite_size.py` (+ figure). Tags **[R]** proved / **[N]** numerical / **[H]** heuristic.*

## 1. Setup and the exact $N\to\infty$ closure  [R]

$N$ all-to-all QIF/theta neurons,
$$\dot\theta_j=(1-\cos\theta_j)+(1+\cos\theta_j)\,\big(\eta_j+Jr(t)-a(t)\big),\qquad
\eta_j\sim g(\eta)=\frac{\Delta/\pi}{(\eta-\bar\eta)^2+\Delta^2},$$
with population rate $r$ (flux through $\theta=\pi$) and adaptation $\tau_a\dot a=-a+\alpha r$. With the
Kuramoto order parameter $Z=\langle e^{i\theta}\rangle$ and Montbrió's conformal map
$$W\equiv \pi r+iv=\frac{1-\bar Z}{1+\bar Z},$$
the Lorentzian + Ott–Antonsen ansatz closes the dynamics **exactly** as $N\to\infty$ on $(r,v)$:
$$\boxed{\;\dot r=\Delta/\pi+2rv,\qquad \dot v=v^2+(\bar\eta-a)+Jr-\pi^2r^2\;}\tag{MPR}$$
(Montbrió, Pazó, Roxin, *PRX* 2015). This is the deterministic backbone; the adaptation makes it a
slow–fast system whose fast $(r,v)$ subsystem is swept by $a$.

## 2. The system-size expansion  [R structure / H covariance]

For finite $N$ the empirical measure fluctuates about the OA manifold. Writing the Kuramoto–Daido
moments $z_m=\langle e^{im\theta}\rangle$ (OA manifold: $z_m=z_1^m$), the discreteness of $N$
oscillators is a fluctuating force on the moment hierarchy. The van Kampen / system-size expansion
(Buice–Cowan; Bressloff, for the empirical-measure dynamics) adapted through the OA/circular-cumulant
closure (Goldobin et al.) gives, to leading order, a **closed Langevin equation for the order
parameter**,
$$\dot z=G(z)+N^{-1/2}\,\sqrt{Q(z)}\,\xi(t),\qquad \langle\xi(t)\bar\xi(t')\rangle=\delta(t-t'),$$
with $Q(z)$ the fluctuation covariance fixed by the instantaneous one-oscillator density (the second
circular cumulant, the leading departure from OA, feeds the closure). Pushed through the conformal map,
$$\boxed{\;d\begin{pmatrix}r\\v\end{pmatrix}=F_{\rm MPR}(r,v)\,dt
+N^{-1/2}\,B(r,v)\,dW_t\;,\qquad BB^{\mathsf T}=\Sigma(r,v).}\tag{Langevin}$$
**Robust facts** (framework-level): the macroscopic noise is (i) $O(N^{-1/2})$, (ii) multiplicative
(state-dependent), (iii) carried by the mean-field variables. The **explicit** covariance $\Sigma(r,v)$
requires the finite-size fluctuation covariance of the QIF order parameter — the circular-cumulant
correction beyond OA — which is technical and is the honest **[H]** step here (a definite computation,
left to the program below).

## 3. Direct validation: the noise is $O(N^{-1/2})$  [N]

We test the central scaling against the **direct finite-$N$ theta network** (no mean-field assumed; $r$
read from $Z$ via the conformal map), at a stable monostable fixed point ($\bar\eta=-8$, $\alpha=0$),
measuring the macroscopic rate variance versus $N$:

| $N$ | 250 | 500 | 1000 | 2000 | 4000 |
|---|---|---|---|---|---|
| $\mathrm{Var}(r)$ | $1.15\times10^{-4}$ | $5.78\times10^{-5}$ | $2.87\times10^{-5}$ | $1.45\times10^{-5}$ | $7.68\times10^{-6}$ |
| $\mathrm{Var}(r)\cdot N$ | 0.0287 | 0.0289 | 0.0287 | 0.0289 | 0.0307 |

$\mathrm{Var}(r)\propto 1/N$ (**log–log slope $-0.98$**, prediction $-1$); $\mathrm{Var}(r)\cdot N$ is
constant across a $16\times$ range. The $O(N^{-1/2})$ Langevin noise amplitude is confirmed directly.
(A run at $\bar\eta=-5$ jumped to the up-state — an independent sighting of the coexisting rest state of
§5.) Figure `figures/mpr_finite_size.png`, panel A.

## 4. Reduction at a macroscopic bifurcation  [H]

At a codim-1 bifurcation of $F_{\rm MPR}$ with critical point $x_*=(r_*,v_*)$, zero eigenvalue and
left/right null vectors $\ell,e$, project (Langevin) onto the slow direction $R=\ell\!\cdot\!(x-x_*)$:
$$\dot R=\underbrace{\alpha(\bar\eta-\bar\eta^*)}_{\mu}+\tfrac{\kappa}{2}R^2
+N^{-1/2}\,s\,\xi_R(t),\qquad s^2=\ell\,\Sigma(x_*)\,\ell^{\mathsf T},\ \ \kappa=\ell\!\cdot\! D^2F[e,e].$$
Rescaling carries this to the **canonical noisy saddle-node** $dR=(\mu+R^2)\,d\tau+\sigma\,dB$ with
effective noise $\sigma_{\rm eff}\propto N^{-1/2}s$ — *exactly the inner equation of the paper's phase
channel*. **If the macroscopic bifurcation is a SNIC** (saddle-node on the rate cycle), the phase-channel
results then transfer to the collective variable, with concrete size scalings:
$$\omega_{\rm coll}\propto\sigma_{\rm eff}^{2/3}\propto N^{-1/3},\qquad
\nu=\mu/\sigma_{\rm eff}^{4/3}\ \Rightarrow\ \text{critical window in }\bar\eta\ \propto N^{-2/3},\qquad
\mathrm{CV}\to0.57,$$
and the collective inter-burst law is the universal quartic-FPT shape. This is the would-be *theorem*:
the finite-size collective phase fluctuations are the phase-edge universality, with $1/\sqrt N$ as the
noise intensity.

## 5. Which edge does the studied network realise? — an honest correction  [N]

The reduction in §4 is conditional on a **SNIC**. At the adaptation operating point
$(J,\Delta,\tau_a,\alpha)=(15,1,15,5)$ the macroscopic bifurcation is **not** a SNIC. The deterministic
diagnostic (figure panel B): as $\bar\eta\to-2.312$ the burst **period stays finite** ($\sim16\to17$,
no divergence) and the **amplitude stays steady** ($\sim2.1$), then the cycle **vanishes abruptly** —
the signature of a **fold of limit cycles** (the *amplitude* edge), with a **coexisting rest state**.
Consistently, the collective $\sigma^{2/3}$ test fails here (mean inter-burst interval vs $\sigma$ gives
slope $-0.18$, not $-2/3$; $\mathrm{CV}\approx0.78$, not $\to0.57$): the noise-induced bursts are
**rate-governed bistable hopping** (rest $\leftrightarrow$ cycle), the paper's "coexisting rest state,
two-way, rate-governed" regime — **not** a clean edge universality.

**This corrects the earlier `NEURAL_MASS_EDGE.md` claim** of a collective *phase* edge: that classifier
call ("SNIC") was a misfire on skewed noise-activated escape intervals. The honest picture at this
operating point is the amplitude edge (fold of cycles) plus a coexisting rest, giving rate-governed
hopping rather than the $\sigma^{2/3}$ phase law.

## 6. Status and program

- **[R]** OA/Lorentzian closure $\Rightarrow$ MPR (exact, cited).
- **[R structure]** system-size expansion $\Rightarrow$ an $O(N^{-1/2})$ multiplicative macroscopic
  Langevin; **[H]** the explicit covariance $\Sigma(r,v)$ (circular-cumulant correction) is not computed
  here.
- **[N]** $\mathrm{Var}(r)\propto 1/N$ (slope $-0.98$) — the $N^{-1/2}$ noise amplitude, validated
  directly on the finite-$N$ network.
- **[N]** the operating-point macroscopic edge is a **fold of cycles** (amplitude edge) + coexisting
  rest, not a SNIC; the collective $\sigma^{2/3}$ phase law is therefore **not** realised here.
- **[H]** the saddle-node reduction $\Rightarrow$ canonical noisy saddle-node and the
  $N^{-1/3}$ / $N^{-2/3}$ phase-edge predictions **await a macroscopic SNIC**.

**The program.** (1) *Locate a macroscopic SNIC* in MPR+adaptation parameter space (an onset where the
rate period diverges) and test the $N^{-1/3}$ collective-timescale law — the clean collective phase
edge. (2) Compute the explicit $\Sigma(r,v)$ (the circular-cumulant fluctuation covariance) to turn §4
into a theorem. (3) The **collective Tracy–Widom amplitude statistics at the present fold of cycles** —
the amplitude edge the network *does* have — is the natural companion result. Each is a concrete,
citable contribution to the next-generation neural-mass programme, which is otherwise deterministic.
```
python3 mpr_finite_size.py scaling   # Var(r) ~ 1/N  (the N^{-1/2} validation)
python3 mpr_finite_size.py diag       # the fold-of-cycles diagnostic
python3 mpr_finite_size.py fig        # -> figures/mpr_finite_size.png
```
