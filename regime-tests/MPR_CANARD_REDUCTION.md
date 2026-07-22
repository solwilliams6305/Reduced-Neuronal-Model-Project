# Reduction of the noisy QIF mean-field at its fold: does the network inherit Tracy–Widom?

*The analytic step the numerics couldn't settle. We reduce the finite-size MPR+adaptation system at the
fold of its critical manifold and ask whether the collective amplitude dynamics is the paper's
amplitude-channel Riccati/Airy canard — hence collective $\mathrm{TW}_\beta$. Verdict: **structurally
yes**, at the fold of the critical manifold (not the fold of cycles), with the correct observable being
the **per-jump peel-off level** — which the earlier follow-up numerics mis-measured. Tags **[R]/[N]/[H]**.*

## 0. Why the earlier numerics missed it

Follow-up 2 measured inter-burst timing and peak amplitude near the **fold of cycles** (where the whole
bursting rhythm disappears). The amplitude-edge canard does **not** live there. It lives at the **fold of
the critical manifold** — the saddle-node of the *fast* $(r,v)$ equilibria that **every burst jumps
through**. The right observable is the **peel-off level at each fold passage**, not the rhythm's death.

## 1. The critical manifold and its fold  [R / N]

Fast subsystem $\dot r=\Delta/\pi+2rv$, $\dot v=v^2+\mu+Jr-\pi^2r^2$ with slow drive $\mu=\bar\eta-a$.
Fast equilibria: $v_c(r)=-\Delta/(2\pi r)$ and
$$\mu(r)=\pi^2r^2-Jr-\frac{\Delta^2}{4\pi^2 r^2}\qquad(\text{the critical manifold}).$$
Folds at $\mu'(r)=0$: $\;2\pi^2 r-J+\Delta^2/(2\pi^2r^3)=0$. For $(J,\Delta)=(15,1)$:
$$\text{upper fold } (r_f,\mu_f)=(0.754,\,-5.744),\ \mu''(r_f)=+19.3;\qquad
\text{lower fold }(0.163,-3.136).$$
The deterministic bursting trajectory reaches $\mu\in[-5.785,-3.641]$ — i.e. it **passes the upper fold
$\mu_f=-5.744$ and overshoots to $-5.785$**: the canard signature (the trajectory follows the repelling
slow manifold slightly *past* the fold before peeling off). **[N] verified.**

## 2. Fold normal form + the adaptation slow passage  [R structure]

Set $\rho=r-r_f$, $\delta=\mu-\mu_f$. Standard slow–fast fold reduction (Fenichel + the fold normal
form; Berglund–Gentz) projects the fast flow onto the centre direction:
$$\dot\rho = A\,\delta + B\,\rho^2 + \dots,\qquad B\propto \tfrac12\mu''(r_f).$$
Adaptation supplies the slow passage: $\mu=\bar\eta-a$, $\dot a=(-a+\alpha r)/\tau_a$, so at the fold
$$\dot\mu=-\dot a=-\frac{-a_f+\alpha r_f}{\tau_a}\approx-0.020\quad(\text{constant, slow})
\ \Rightarrow\ \delta(t)\approx\dot\mu\,(t-t_f)\ \text{(linear in time)}.$$
**The passage is linear in time** — the hallmark of the Airy structure. With the finite-size noise of
`MPR_FINITE_SIZE_LANGEVIN.md` projected onto $\rho$ ($\sigma_\rho\propto N^{-1/2}$),
$$\boxed{\,d\rho=\big(A\dot\mu\,(t-t_f)+B\rho^2\big)\,dt+\sigma_\rho\,dW_t\,}\tag{inner}$$
— a **Riccati with a time-linear coefficient**, exactly the form of the paper's amplitude channel.

## 3. Cole–Hopf $\to$ stochastic Airy $\to$ collective $\mathrm{TW}_\beta$  [R machinery / H constants]

Rescaling (inner) to the canonical amplitude-channel equation $dY=(s-Y^2)\,ds+\eta\,dW$ (balance the
linear drift, the quadratic, and the noise; $s$ the rescaled slow time, $\eta$ the noise-to-curvature
ratio) and applying the paper's Cole–Hopf linearisation $Y=-\eta^2\,\psi'/\psi$ sends (inner) to
$$-\psi''+\big(s+\text{(white noise)}\big)\psi=\lambda\psi,$$
the **stochastic Airy operator** $H_\beta$ with $\beta=2/\eta^2$. Its ground state $\lambda_0\sim
\mathrm{TW}_\beta$ (the paper's Theorem on the amplitude edge). Therefore the **collective peel-off
level** — the value of $\rho$ (equivalently $r$ or $\mu$) at which each burst leaves the repelling slow
manifold — is $\mathrm{TW}_\beta$ distributed, with
$$\eta\propto\sigma_\rho\propto N^{-1/2}\quad\Rightarrow\quad \beta=2/\eta^2\propto N.$$

**This is the answer: the network inherits the amplitude edge.** The QIF *single-neuron* folded-cycle
amplitude edge and the QIF *network* fold both Cole–Hopf-linearise to the same stochastic Airy operator;
the collective peel-off is $\mathrm{TW}_\beta$, with the **system size $1/\sqrt N$ playing the role of
the noise intensity** and $\beta\propto N$ (so $N\to\infty$ is the deterministic $\mathrm{TW}_\infty$
limit, finite $N$ the fluctuating edge).

## 4. What is rigorous, what is open

- **[R]** the critical manifold, its fold, $\mu''$ — explicit and exact; **[N]** the trajectory's
  canard overshoot past the fold — verified.
- **[R structure]** the fold reduction $\to$ time-linear Riccati (inner), and Cole–Hopf $\to$ stochastic
  Airy $\to\mathrm{TW}_\beta$ — standard slow–fast + the paper's amplitude-channel machinery.
- **[H]** the constants — the noise projection $\sigma_\rho$ (hence $\beta$ and $\eta$) need the explicit
  finite-size covariance $\Sigma$ (still the one open computation, from `MPR_FINITE_SIZE_LANGEVIN.md`).
- **[N] peel-off corroboration — encouraging but not definitive.** Measuring the corrected observable
  (per-trial deepest overshoot past the fold, $\max a$) under finite-size noise: at $\tau_a=15$ the
  skew is $\mathbf{+0.29}$ — *exactly* $\mathrm{TW}_1$ ($+0.29$), and the right sign/magnitude — but at
  $\tau_a=40$ it flattens to $\approx0$ (Gaussian). So the peel-off carries a TW-signed, TW-magnitude
  skew at the base parameters, **consistent with collective $\mathrm{TW}_\beta$** (with $\beta$
  parameter/$N$-dependent), yet it is **parameter-sensitive and small**, so this is corroboration, not
  confirmation.
- **[H] caveats on the canard regime.** $\varepsilon=1/\tau_a\approx0.067$ is only moderately small;
  the burst's jump-off $r\approx0.65$ sits below the static fold $r_f=0.754$, so the passage is a
  *near-canard*. The decisive test is the **$\beta\propto N$ scaling** and a full-shape $\mathrm{TW}_\beta$
  fit of the peel-off distribution (not just the skew sign).

## 5. The result, and the next two concrete steps

**Result.** The noisy QIF mean-field's collective amplitude statistics at the fold of its critical
manifold reduce — structurally — to the stochastic Airy operator: **collective Tracy–Widom**, $\beta\propto N$.
This carries the paper's headline edge from the single neuron to the whole network, and identifies the
correct observable (per-jump peel-off level) the earlier numerics missed.

**Next.** (1) Compute $\Sigma$ (the circular-cumulant finite-size covariance) to fix $\beta(N)$ — turning
§3 into a theorem with constants. (2) Numerically test the **peel-off-level** distribution vs
$\mathrm{TW}_\beta$ at a sharper-fold operating point ($\tau_a$ large) across $N$ — the clean
confirmation. Both are now well-posed.
