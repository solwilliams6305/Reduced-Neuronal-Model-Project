# Closing the two gaps: the finite-size covariance Σ, and the β∝N test

*The two steps that would turn the canard reduction (`MPR_CANARD_REDUCTION.md`) into a theorem with
constants. Both now done as far as they go: **Step 1 succeeds up to an O(1) coupling renormalisation;
Step 2 corroborates the TW magnitude but not yet the β∝N scaling.** Tags **[R]/[N]/[H]**.*

## Step 1 — the finite-size covariance Σ  [R scaling + leading constant / H renormalisation]

**Derivation.** The Kuramoto order parameter is the empirical mean $Z=\frac1N\sum_j e^{i\theta_j}$. Under
molecular chaos (leading order), its finite-$N$ covariance is just the variance of an empirical mean of
$N$ points on the circle with one-point density $\rho$:
$$\langle|\delta Z|^2\rangle=\frac{1-|Z|^2}{N},\qquad
\langle\delta Z^2\rangle=\frac{z_2-z_1^2}{N}=0\ \ (\text{OA: }z_2=Z^2).$$
Montbrió's conformal map $e^{i\theta}=(1+iv)/(1-iv)$ and the Lorentzian (OA) density give, by residues,
$Z=\dfrac{(1-\pi r)+iv}{(1+\pi r)-iv}$, whose **inverse is explicit**,
$$r=-\frac1\pi\,\mathrm{Re}\!\frac{Z-1}{Z+1},\qquad
\frac{\partial r}{\partial Z}=-\frac{1}{\pi(Z+1)^2}.$$
Hence the macroscopic rate variance is **closed-form**:
$$\boxed{\ \mathrm{Var}(r)\,N=2\,|\partial_Z r|^2(1-|Z|^2)=\frac{2(1-|Z|^2)}{\pi^2\,|Z+1|^4}\ }$$
and likewise $\Sigma_{(r,v)}=\tfrac1N M(r,v)$ in closed form — the explicit $\Sigma$ the Langevin needed.

**Test.** At the fixed point $\bar\eta=-8$ ($r^*=0.0596$, $|Z|=0.955$): predicted $\mathrm{Var}(r)N=0.081$
vs **measured $0.029$** (finite-$N$ theta network). The **scaling is exact** ($\propto1/N$, slope $-0.98$),
the **order of magnitude right**, but the constant is **$\approx2.8\times$ too large**.

**Why — two candidates ruled out, mechanism identified.**
- *Not the linear-response Jacobian.* The $(r,v)$ Jacobian at the fixed point gives eigenvalues
  $\{-4.06,-6.63\}$ (coupled) vs $-5.35\pm0.37i$ (no field feedback) — coupling *slows* relaxation
  (ratio $0.76$), and the 2×2 Lyapunov balance $AC+CA^{\mathsf T}=-D$ returns $\mathrm{Var}(r)N=0.084$,
  i.e. **essentially the bare value** — the feedback does not close the gap.
- *Not the sampling protocol.* Re-running with **random i.i.d.** Lorentzian $\eta$ (not deterministic
  quantiles) gives $\mathrm{Var}(r)N=0.035$ — still $\sim2.3\times$ below $0.081$.
- *The mechanism is sub-Poissonian phase correlations.* The measured equal-time order-parameter variance
  ($\sim0.03$) is **below** the i.i.d. (Poissonian) value $(1-|Z|^2)/N=0.081$: the deterministic coupled
  flow **anti-correlates** the phases (self-organisation, more uniform than random), suppressing the
  fluctuations $\sim2.4\times$. This is a genuine many-body effect — the **pair-correlation / second
  circular cumulant** — that the leading (uncorrelated) closure omits.

So: **$\beta(N)=cN$, scaling exact, $c$ to leading order (within $\sim2.4\times$)**; the residual **[H]**
is now sharply identified as the sub-Poissonian pair-correlation constant (a circular-cumulant
computation), not the field-feedback dressing I first guessed.

## Step 2 — the β∝N scaling of the peel-off  [N]

The reduction predicts the peel-off level $\sim\mathrm{TW}_\beta$ with $\beta=2/\eta^2$ and
$\eta\propto\sigma_{\rm eff}\propto N^{-1/2}$, i.e. $\beta\propto N$. Lowering the (macroscopic) noise
$\sigma$ should therefore **raise** $\beta$ and **lower** the peel-off skew toward $0$ (TW$_\infty$); the
width should scale as the edge power. Measured (peel-off $=\max a$, fold regime, $300$ trials):

| $\sigma\ (\propto N^{-1/2})$ | 0.06 | 0.12 | 0.22 |
|---|---|---|---|
| skew | $+0.18$ | $+0.20$ | $+0.19$ |
| std | $0.018$ | $0.031$ | $0.050$ |

**Read honestly:** the skew sits at **$+0.19$ — TW-signed and TW-magnitude** (between $\mathrm{TW}_2$
$=0.22$ and $\mathrm{TW}_4=0.17$) — but it is **roughly constant in $\sigma$**, not growing as
$\beta\propto\sigma^{-2}$ would require, and the width scales **$\sim\sigma^{0.8}$**, not the edge
$\sigma^{2/3}$. So the measurement is **not in the asymptotic edge regime** at these parameters: the
TW-magnitude skew is real corroboration, but the **β∝N scaling is not confirmed**.

**Direct finite-$N$ test (the proper version).** Bursting **finite-$N$ theta network** with a
per-burst peel-off estimator, $N=400$ vs $1600$:

| | $N=400$ | $N=1600$ |
|---|---|---|
| peel-off skew | $+0.21$ | $+0.21$ |
| peel-off std | $0.51$ | $0.51$ |

The skew is **$+0.21$, TW-magnitude** (between $\mathrm{TW}_2$/$\mathrm{TW}_4$) — and remarkably
**consistent**: the peel-off carries a positive, TW-magnitude skew across *every* setting tried
($+0.19$ to $+0.29$). **But it is $N$-independent**, and the width does not shrink with $N$: at these
parameters the peel-off spread is dominated by **deterministic** burst-to-burst variability (the
slow–fast canard sensitivity), with the finite-size noise sub-dominant. So **β∝N is not observable
here** — not because the noise isn't $N^{-1/2}$ (it is, validated), but because it is swamped by the
deterministic spread in this observable.

**The fix (now precise).** A clean β∝N test needs a regime where the finite-size noise *dominates* the
peel-off: a **sharp canard** ($\varepsilon=1/\tau_a\ll1$) where the deterministic burst-to-burst spread
collapses and the peel-off is set by the noise, then the per-passage overshoot $\mu_f-\min\mu$ scales as
$N^{-1/3}$ with $\mathrm{TW}_\beta$, $\beta\propto N$. That is a singular-perturbation regime
($\tau_a\to\infty$) the present integration cannot reach economically — a dedicated, stiff,
sharp-canard computation.

## Net status of the theory track

- **[R/N]** Finite-size Langevin: derived; $O(N^{-1/2})$ noise validated (slope $-0.98$).
- **[R structure]** Canard reduction at the critical-manifold fold $\Rightarrow$ stochastic Airy
  $\Rightarrow$ **collective $\mathrm{TW}_\beta$**, $\beta\propto N$.
- **[R + leading constant]** $\Sigma$ in **closed form** $\mathrm{Var}(r)N=2(1-|Z|^2)/\pi^2|Z+1|^4$;
  scaling exact, constant within $\sim2.8\times$ (coupling renormalisation = the one residual **[H]**).
- **[N]** Peel-off skew TW-signed & TW-magnitude ($+0.19$); **β∝N not yet confirmed** (not asymptotic
  regime).

**Two clean, finite tasks remain, both now sharply posed:** (1) the single O(1) self-consistent
renormalisation of $\Sigma$ (pins $c$ in $\beta=cN$); (2) the β∝N / full-shape $\mathrm{TW}_\beta$ fit at
a sharp-canard, finite-$N$ operating point with the per-passage estimator. The structural claim —
**the QIF network inherits the paper's Tracy–Widom edge, with $1/\sqrt N$ as the noise intensity** —
stands; what remains is constants and a cleaner numerical confirmation, not new structure.
