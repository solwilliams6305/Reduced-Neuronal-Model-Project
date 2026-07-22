# Cracking the Weber-TW block — a strategy

_June 2026. A worked-out plan of attack for the one genuinely open object: a canonical / RRV-type
characterization of the cusp (Weber) escape law. Supersedes the soft-edge framing in
`DYSON_WEBER_NOTES.md` for the marginal. Probe: `coupled-atlas/ladder_kurtosis_probe.py`._

---

## 1. The block, stated as constraints

We want the cusp analogue of "stochastic Airy operator = Tracy–Widom." Any correct answer must
reproduce **every** signature we have measured:

- **marginal moments** skew ≈ **+0.60**, excess kurtosis ≈ **−0.23** (β=2) — *sub-Gaussian, light
  tails*;
- it is the first-node law of the **swept 2nd-order parabolic-cylinder (Weber)** equation
  $u''=(\mathrm{sign}(Y)Y^2-\eta\xi)u$;
- it is **not** the ground state of any static self-adjoint operator (skew saturates, kurtosis sign
  wrong) and **not** the 1st-order pitchfork escape (kurtosis sign wrong);
- it sits in a **negative-kurtosis valley** along the catastrophe ladder $V=\mathrm{sign}(Y)|Y|^q$:
  exkurt $= +1.6,\,+0.79,\,+0.04,\,-0.22,\,-0.23,\,-0.13,\,+0.14$ for $q=0.5\dots3$ — i.e. the fold
  ($q=1$) sits at the **zero-crossing** (consistent with TW's small +kurtosis) and the cusp ($q=2$)
  near the **minimum**. Skew rises monotonically with $q$.
- the point process has genuine consecutive-level rigidity (the Weber analogue of Airy's ≈0.50).

The kurtosis valley is the key new clue: **the turning order $q$ controls the kurtosis sign**, so the
object is intrinsically a *ladder* phenomenon, with the cusp as a distinguished (most sub-Gaussian)
point — not a perturbation of the soft edge.

---

## 2. The reframe — RRV is really a *Riccati-explosion* theorem; generalise THAT

The mistake was hunting for a static self-adjoint operator. RRV's own proof does not use one
essentially — it uses the **Riccati diffusion**. Apply Cole–Hopf $p=u'/u$ to the swept inner
equation:

$$\frac{dp}{dY}= \big(W(Y)-p^2\big)-\eta\,\xi(Y),\qquad W(Y)=\mathrm{sign}(Y)|Y|^{q}.$$

As $Y$ sweeps down, $p\to-\infty$ at each node of $u$; **successive explosions = successive peel-offs
= the point process**, and the **first explosion = the escape marginal**. So:

> **Reframe.** The cusp law is the first-explosion law of the swept-Weber Riccati diffusion
> $dp=(\mathrm{sign}(Y)Y^2-p^2)\,dY+\eta\,dW$. Its generator
> $\;\mathcal L_Y=\tfrac{\eta^2}{2}\partial_p^2+(\mathrm{sign}(Y)Y^2-p^2)\,\partial_p\;$
> is a genuine (non-self-adjoint) diffusion generator, time-inhomogeneous through $W(Y)$.

This is exactly the "directional, non-self-adjoint" object the numerics forced on us — now it is a
*standard probabilistic object* (a 1-D diffusion with explosion), not a mystery. For $q=1$ it is RRV's
Riccati diffusion, whose explosion law is TW via Painlevé II. The block is: **solve the $q=2$ case.**

---

## 3. Why sub-Gaussian — the mechanism, now grounded

Past the turning point the WKB phase $\Theta(Y)=\int\sqrt{|W|}\,dY\sim |Y|^{(q+2)/2}$ accumulates, and
the first node lands at $\Theta\approx\pi/2$. The node window width is $\delta Y\sim 1/\Theta'(Y)$:
the **faster the accumulation (larger $q$), the tighter and more bounded the window ⇒ lighter tails ⇒
exkurt $<0$**; the slow (sub-linear, small $q$) regime is broad-windowed and heavy-tailed (exkurt
$\gg0$, e.g. $+1.6$ at $q=0.5$). The fold ($q=1$, linear phase) is the marginal case sitting at
exkurt $\approx0$. This *predicts the whole valley* and explains why no static/1st-order reduction can
work: they discard the oscillatory phase that bounds the window.

---

## 4. Four attack routes, ranked

### Route A — the explosion-probability PDE → Painlevé IV  *(primary, most likely to give the exact law)*
Let $h(p,Y)=\mathbb P(\text{no explosion above }Y\mid p)$. It solves the backward Kolmogorov equation
$$\partial_Y h+(\mathrm{sign}(Y)Y^2-p^2)\,\partial_p h+\tfrac{\eta^2}{2}\partial_p^2 h=0,\qquad
h(-\infty)=0,\;h(+\infty)=1,$$
and the escape CDF is read off as $Y\to-\infty$. For $q=1$ this BK problem is equivalent to the
Hastings–McLeod Painlevé II representation of TW.
**Conjecture (the lead):** since the deterministic skeleton is the **parabolic-cylinder (Weber)**
equation, and parabolic-cylinder functions $D_\nu$ are the classical special-function solutions of
**Painlevé IV** (exactly as Airy is for Painlevé II), the Weber-TW law is governed by a **Painlevé IV
transcendent** — the cusp analogue of Hastings–McLeod.
**First steps:** (i) derive the BK PDE precisely with the correct Stratonovich/sweep signs; (ii) solve
it numerically (1-D in $p$, march in $Y$) and confirm it reproduces the simulated escape law + the
$q$-valley — this *validates the reframe* cheaply; (iii) look for a self-similar / Bäcklund reduction
to PIV; check the PIV special-parameter family against the moments.

### Route B — the parabolic-cylinder ("Weber") Fredholm kernel  *(the exact-law route, parallel to A)*
TW $=\det(1-K_{\mathrm{Ai}})$ with $K_{\mathrm{Ai}}(x,y)=\int_0^\infty \mathrm{Ai}(x+s)\mathrm{Ai}(y+s)\,ds$.
Build the analogue from parabolic-cylinder functions,
$K_{W}(x,y)=\int_0^\infty D_\nu(x+s)D_\nu(y+s)\,ds$ (or the PIV integrable kernel), and compute the gap
probability $\det(1-K_W)|_{(t,\infty)}$.
**First step:** numerically assemble $K_W$ on a grid (Nyström), get the largest-eigenvalue CDF, extract
skew/exkurt and compare to $(+0.60,-0.23)$ and the valley. **Caveat & cross-check:** the *Pearcey*
kernel (two merging edges) is the other cusp object; we found "Weber, not Pearcey" at the ODE-order
level, so test BOTH kernels' marginals — only the one with **negative** excess kurtosis can be ours
(see §5).

### Route C — the bounded-window canonical reduced model  *(fast, explanatory; pins the mechanism)*
Formalise §3: model the first node as the first passage of a noisy accumulating phase
$d\Theta=\sqrt{|W(Y)|}\,dY+(\text{phase noise})$ to threshold, with $W=|Y|^q$. Derive
$\mathrm{exkurt}(q)$ and test against the measured valley **including the $q=3$ up-turn** (which a
correct model must reproduce). Success gives a 1-parameter canonical object and the cleanest proof of
the mechanism. **First step:** simulate the reduced phase model; fit the valley.

### Route D — literature anchor: the noisy folded node  *(may shortcut the whole thing)*
The deterministic folded-node SAO count is *already* a Weber-equation result (Wechselberger;
Desroches–Krauskopf–Osinga). The **stochastic** version — Berglund, Gentz & Kuehn, *Hunting French
ducks in a noisy environment* and *From random Poincaré maps to stochastic mixed-mode oscillation
patterns* — derives the distribution of the SAO count / escape via a **random Poincaré-map / Markov
kernel** built on the Weber dynamics. This is the closest prior art.
**First step:** read those two papers; check whether their Poincaré kernel, in our scaling, **is** the
Weber-TW object (its spectrum/contraction giving our law), and position our edge-statistics framing +
the kurtosis valley relative to their SAO-count results. *(Verify the exact statements against the
papers — these are recalled from memory.)*

---

## 5. Moment-fingerprint cross-check
$(\text{skew},\text{exkurt})=(+0.60,-0.23)$ is **incompatible** with every soft-edge / extreme-value
family — TW ($+0.22,+0.09$), Gamma, inverse-Gaussian, Gumbel all have *positive* excess kurtosis. It
*is* compatible with **bounded / Beta-like** laws (which are routinely platykurtic). This corroborates
the bounded-window picture and is a quick filter: fit a Beta (and a reflected/Generalised-Gumbel) to
the standardised cusp law by KS — a phenomenological identification that constrains Routes A–C.

---

## 6. Recommended order of attack
1. **Route A step (ii)** — ✅ **DONE** (`RICCATI_EXPLOSION_NOTES.md`). The Riccati-explosion SDE
   reproduces the swept-Weber law across the ladder (cusp +0.59/−0.29 vs +0.64/−0.23; valley
   reproduced) — **reframe confirmed**. The deterministic BK-PDE is formulated; an *implicit*
   (Crank–Nicolson) solver is the remaining coding step.
1b. **Tail asymptotics** — ✅ **DONE** (`EXACT_LAW_TAILS_NOTES.md`). Freidlin–Wentzell gives
   $\alpha_{\rm L}=2q+1$ ($F\sim e^{-|s|^{2q+1}/(4(2q+1))}$) and $\alpha_{\rm R}=\tfrac32 q$
   ($S\sim e^{-(4/3)|s|^{3q/2}}$), **exact at $q=1$ = TW (constants and all)**. Cusp: left
   $e^{-|s|^5/20}$, right $e^{-(4/3)|s|^3}$ — the $|s|^5$ parabolic-cylinder/PIV signature; MC slopes
   2.75 ($q{=}1$), 5.23 ($q{=}2$) confirm. These are necessary conditions the PIV identification must
   satisfy — and they match.
2. **Route C** — the bounded-window model; nail the mechanism and the valley shape. *(Cheap, high
   insight.)*
3. **Route B** — assemble the parabolic-cylinder kernel; test its gap-probability moments (and the
   Pearcey marginal as the foil). *(Heavier; the exact-law payoff.)*
4. **Route A step (iii) + Route D** — the PIV reduction and the Berglund–Gentz–Kuehn anchor, in
   parallel. *(The theorem / the citation home.)*

## 7. What counts as "solved"
Either (a) an **exact law** — a Fredholm determinant of a named (parabolic-cylinder / PIV) kernel, or a
Painlevé-IV transcendent representation — reproducing skew/exkurt$(q)$ and the point-process rigidity;
or (b) a **rigorous statement** of the form *"the first explosion of the swept-Weber Riccati diffusion
is distributed as [named law]"* — the cusp analogue of RRV. Route A is the most direct path to (b),
Route B to (a); Route C explains *why* and Route D may already contain pieces of both.
