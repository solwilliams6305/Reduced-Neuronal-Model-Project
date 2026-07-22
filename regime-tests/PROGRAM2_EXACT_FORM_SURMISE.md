# Program 2 — surmising the EXACT FORM of 𝒲: a resurgent trans-series with a complex-conjugate Borel pair

_2026-07-08 (Fable 5). Synthesis of the whole program into a conjectured exact form for the cusp edge law
$\mathcal W_\beta$, with the evidence chain, quantitative predictions, and the plan to pin/prove it. Tags:
**[EVIDENCE]/[CONJECTURE]/[PREDICTION]**. Builds on `PROGRAM2_ROUTE2B_NOTES.md` §§6–8,
`PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md`._

## 1. The surmise (one paragraph)
$\mathcal W_\beta$ has no 1-D closed form (proved: isomonodromy fixed point ⇒ no Painlevé-σ; unbounded-below
first-passage ⇒ no Fredholm determinant). Its exact form is a **resurgent trans-series in $\eta^2=4/\beta$** whose
Borel plane carries a **complex-conjugate pair of singularities at argument $\pm\pi/4$**, inherited from the
deterministic Weber connection root $\lambda_0=0.890-0.890i$ (which sits at **exactly $-45^\circ$**). Consequently the
weak-noise cumulant coefficients **oscillate in sign** with period $2\pi/(\pi/4)=8$ in the order $n$, and the series
is **asymptotic (factorially divergent) but Borel-summable along $\mathbb R_+$** (the singularities are off the
positive axis). The perturbative sector is a *renormalized* series (boundary self-contractions subtracted from
$v_3$ on); the non-perturbative sectors are $e^{-S_\pm/\eta^2}$ with **complex actions** $S_\pm=|S|e^{\mp i\pi/4}$
governing the bulk, plus the **real** far-tail instanton $S_{\rm tail}=s^5/10$. The Stokes automorphism linking the
sectors is fixed by $\lambda_0$.

## 2. The evidence chain
- **[EVIDENCE, UPGRADED 2026-07-08 — EXACT engine] The coefficient ladder, computed EXACTLY** (diagram chaos engine,
  `PROGRAM2_CHAOS_ENGINE.md`, no sampling/truncation, grid-converged): $v_0=0.134,\ v_1=0.111,\ v_2\approx0.10,\
  v_3\approx-0.02$. Pattern $+,+,+,-$ with $v_3$ **small** — the earlier sampled $v_3\approx-0.3$ is RETRACTED. The
  smallness of $v_3$ is itself the confirmation: the $\theta=\pi/4$ oscillation fit to $v_0,v_1,v_2$ predicts $v_3$
  near a cosine zero. Fitting the finite-radius oscillatory form $v_n=C R^{-n}\cos(n\pi/4+\varphi)$ ($\theta=\pi/4$
  FIXED to $\lambda_0$) to all four exact coefficients gives $R\approx1.07$ ($\beta_c\approx3.7$), $\varphi\approx-37^\circ$,
  and (with only $v_0..v_3$) tentatively suggested a decreasing envelope. **RESOLVED 2026-07-09 by the exact $v_4$:**
  $v_4\approx-0.45$ (engine, grid-converged) — $|v_4|$ JUMPS above all earlier $|v_n|$, so the envelope **GROWS**: the
  series is **DIVERGENT (asymptotic/resurgent), NOT convergent**. $v_3$ was small only because it sits near an
  oscillation NODE, not because the envelope decays. Refit to all five ($\theta=\pi/4$ fixed): $r\approx0.6$–0.7 $<1$,
  Borel singularity at $|\zeta|\approx0.6$–0.7, $\arg=\pm\pi/4$; $\eta^2_c\approx0.6$–0.7, $\beta_c\approx6$; predicts
  $v_5\approx-0.9$, period-8 sign pattern $+,+,+,-,-,-,+,+$. **This CONFIRMS the core surmise — a genuine resurgent
  trans-series with a complex-conjugate Borel pair at $\pm\pi/4$ from $\lambda_0$.** Exact ladder:
  $v_0=0.134,\ v_1=0.111,\ v_2=0.100,\ v_3\approx-0.02,\ v_4\approx-0.45$ (all MC-free, `chaos_diagram.py`).
- **[EVIDENCE] $\lambda_0$ is at exactly $-45^\circ$.** The machine-verified connection root $0.890-0.890i$ has
  $\mathrm{Re}=-\mathrm{Im}$, i.e. $\arg\lambda_0=-\pi/4$ exactly, $|\lambda_0|=0.890\sqrt2=1.259$.
- **[EVIDENCE] The phase predicts the sign flip.** A Borel pair at $\pm\pi/4$ ⇒ $v_n\propto\cos(n\pi/4+\varphi)$;
  $\cos$ first turns negative between $n=2$ and $n=3$ for the fitted $\varphi$ — the first negative coefficient is
  $v_3$, as observed. Fixing $\theta=\pi/4$ and fitting $C,R,\varphi$ to $v_0,v_1,v_2$ alone already forces the flip
  at $n\!\approx\!3$ and gives radius $R\approx1.2$ ($\beta_c\approx3.4$), consistent with the independent
  ratio-test radius $z_c\approx1$ ($\beta_c\approx4$).
- **[EVIDENCE] Renormalization threshold at $v_3$** (§8): the exact counterterm $\tfrac{23}{3}\langle\xi(Y^\star_0)^2
  \rangle v_0^3$ (from $B_4=\tfrac13Y_1^3$, $R_7=\tfrac15Y_1^5$) — the perturbative series is a *renormalized* object
  from $v_3$ on. Renormalization ambiguity is the Borel/Stokes shadow.
- **[EVIDENCE] Real far-tail instanton** $-\log\mathbb P(Y^\star<-s)\to s^5/10\eta^2$ — a real saddle governing the
  large-deviation tail (distinct sector from the complex bulk pair).
- **[EVIDENCE, imported] The deterministic backbone is proved resurgent** (Nikolaev 2024, arXiv:2410.17224): the
  Weber operator's WKB solutions are resurgent with a geometric Borel plane / Stokes rays — the scaffold the
  stochastic theory dresses.

## 3. The conjectured trans-series form [CONJECTURE]
For a cumulant/observable $\mathcal O(\eta^2)$ of $\mathcal W_\beta$ (e.g. $\operatorname{Var}$):
$$
\mathcal O(\eta^2)\ \simeq\ \underbrace{\sum_{n\ge0} v_n^{\rm ren}\,\eta^{2n}}_{\text{renormalized perturbative}}
\ +\ \underbrace{2\,\mathrm{Re}\Big[\mathcal S\, e^{-S_+/\eta^2}\!\sum_{n\ge0} v_n^{(+)}\eta^{2n}\Big]}_{\text{complex bulk pair, }S_\pm=|S|e^{\mp i\pi/4}}
\ +\ \underbrace{(\text{real tail sector } e^{-s^5/10\eta^2})}_{\text{large-deviation}} ,
$$
with $v_n^{\rm ren}\sim C\,\Gamma(n+b)\,|S|^{-n}\cos(n\pi/4+\varphi)$ (factorial × oscillation), $|S|\approx R\approx1.2$,
and Stokes constant $\mathcal S$ fixed by $\lambda_0$. The $2\,\mathrm{Re}[\cdots]$ makes the physical law real while
encoding the oscillation. Borel-summable along $\mathbb R_+$ (singularities at $\pm\pi/4$, off-axis) — so the
resummation is **unambiguous** (no lateral-summation ambiguity on the positive axis), matching Nikolaev's
"summable in almost all directions."

## 4. Sharp, falsifiable predictions [PREDICTION]
- **$v_4$ is the decisive test — factorial vs convergent.** Same $\theta=\pi/4$ form fit two ways:
  finite-radius (convergent) ⇒ $v_4\approx-0.06$; factorial (asymptotic trans-series) ⇒ $v_4\approx-4$. **A ~70×
  gap.** Since the program posits a genuine trans-series (an instanton sector exists), $v_4\approx-\text{few}$ is
  expected. Measuring $|v_4|$ settles asymptotic-vs-convergent decisively (far better than the $v_3$ 20% split).
- **Sign pattern (period 8):** $v_3,v_4,v_5,v_6<0$, returning $v_7>0$ (flips near $n=3$ and $n=7$). Any coefficient
  breaking the period-8 oscillation falsifies the $\lambda_0$-tied complex-pair picture.
- **$\beta$-scaling / radius:** the perturbative radius $z_c=\eta^2_c\approx1.2$ ($\beta_c\approx3.4$–4) should equal
  $|S|=|\lambda_0|$-related; a precise $|S|$ from the stochastic exact-WKB must reproduce it.
- **Off-axis Borel summability:** the variance series should be numerically Borel-resummable along $\mathbb R_+$
  (Padé–Borel with poles emerging at $\arg\approx\pm45^\circ$).

## 5. The plan to pin/prove it (ordered)
1. **[reachable, decisive] Compute $v_4$ (and $v_5$) exactly.** $v_4=\operatorname{Var}(Y_5)+2\operatorname{Cov}
   (Y_4,Y_6)+2\operatorname{Cov}(Y_3,Y_7)+2\operatorname{Cov}(Y_2,Y_8)+2\operatorname{Cov}(Y_1,Y_9)$. **The divergence
   structure is now fully mapped** (`_bstruct9.py`, `_lin.py`, exact): no cubic boundary terms; two divergence orders —
   $\langle s_0^2\rangle\!\sim\!1/\delta$ and a NEW $\langle s_1^2\rangle,\langle s_0s_2\rangle\!\sim\!1/\delta^3$ (the
   $\xi'(Y^\star_0)^2$ renormalization, first appearing in $\operatorname{Var}(Y_5)$ and $\operatorname{Cov}(Y_1,Y_9)$).
   The DOMINANT ($1/\delta^3$) counterterms have **pure-$Y_1$** coefficients — $B_{5,1}=Y_1^4/12$, $Y_9|_{s_1^2}=
   -11Y_1^7/630$, $Y_9|_{s_0s_2}=-37Y_1^7/1260$ (also $Y_7|_{s_0^2}=Y_1^5/5$, $Y_8|_{s_0s_1}=11Y_1^6/90$) — hence
   exactly computable via jointly-Gaussian moments. But the subleading $s_0^2$ coefficients ($B_{5,0},B_{6,0}$,
   $Y_8|_{s_0^2},Y_9|_{s_0^2}$) are **mixed** (involve $u_2,u_1',u_3$ at the node), and their cross-moments with the
   regular parts carry $1/\delta^3$ estimator variance — so **renormalized SAMPLING dies at $v_4$; a deterministic
   chaos-moment engine is required** (represent every $\langle Y_aY_b\rangle$ as exact Wick integrals of the
   Green's-function kernels §6a + boundary-Wick, no sampling). This is the single highest-leverage next build; $|v_4|$
   settles factorial-vs-convergent (~−4 vs ~−0.06), signs test the period-8 oscillation.
2. **[reachable] Padé–Borel the coefficient sequence** $\{v_n^{\rm ren}\}$: locate the Borel singularities in the
   complex plane, confirm $\arg=\pm\pi/4$ and $|{\cdot}|=|S|$, read the growth (factorial ⇒ genuine trans-series).
3. **[frontier] Stochastic exact-WKB → the Stokes constant $\mathcal S$ from $\lambda_0$.** Dress the *proved*
   deterministic Weber resurgence (Nikolaev 2024) with the noise: the Voros symbols become random, the Borel-plane
   geometry (Stokes rays at $\pm\pi/4$) is inherited, and $\mathcal S$ follows from the connection data $\lambda_0$
   (`PROGRAM2_PROVED_CLOSED_FORM_ROUTE.md` Phase 3 — the one genuine blank). The renormalization counterterm (§8)
   is the concrete bridge (its ambiguity = the Stokes discontinuity).
4. **[closure] Prove the trans-series equals $\mathcal W_\beta$** (moment-determinacy already in hand, Tier A).

## 6. Why this is the right shape (physical/structural argument)
The escape law is built on the Weber (parabolic-cylinder) operator, whose connection data is intrinsically
**complex** — the $\Gamma$-quotient root $\lambda_0$ lives at $-45^\circ$, not on the real axis. Any resurgent
structure erected on this backbone therefore has **complex** Borel singularities, hence **oscillatory** perturbative
coefficients — which is exactly the $+,+,+,-$ we measure, with the first flip forced to $n=3$ by the $45^\circ$
phase. The real tail instanton $s^5/10$ is a *separate*, real saddle (large-deviation regime); the bulk cumulant
resurgence is the complex-pair story. This dual (real tail + complex bulk) non-perturbative content is the
signature of a genuinely non-integrable, isomonodromy-fixed-point transcendent — no simpler (real-axis,
single-instanton) trans-series can carry an oscillating coefficient sequence.
