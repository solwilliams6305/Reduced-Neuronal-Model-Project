# Brainstorm: new angles on `a_min` (#3) — cast wide, try things, fail freely

**Target restated.** `a_min(c)` = smallest SAO amplitude in the folded-node
funnel. Empirically `a_min ~ exp(−C·μ^{−p})`, `p ≈ 0.7` (exponentially small),
all data in the transition regime `μ < ε^{1/2}`, and the realised value is set by
the global return ✕ the transition-regime inner solution. It's the keystone of
the whole MMO noise/staircase side (#4, #5, #6). Standard matched asymptotics give
only the (loose, algebraic) upper bound. So this is a **beyond-all-orders,
transition-regime, partly-global** object — exactly the kind of thing that yields
to ideas borrowed from elsewhere.

**The single highest-value first move (frames every theory below):** *pin the
exponent.* Is `p` really `0.7`, or is it a clean `2/3` or `3/4` masked by fixed
`ε`? The measured `0.7` comes from one `ε=0.08`; `a_min` almost certainly obeys
`exp(−c/(μ^a ε^b))`, which at fixed `ε` collapses to `exp(−c'/μ^a)`. **Sweep `ε`**
to separate the two, and the clean exponent — the thing all the analytic routes
must reproduce — will emerge. Cheap, decisive, and it kills or confirms each idea
below by its predicted exponent.

---

## Tier 1 — most promising (the exponential-asymptotics family)

**1. Complex-time singularities / exact WKB (Voros–Écalle).** `a_min` is
exponentially small, so it is governed by the *nearest singularity of the solution
in complex time* — the "instanton" of the deterministic flow. The canard fold is
Airy (`V=−u'/u`, `u=Ai(−ξ)`; see `TONIC_CMID_AIRY.md`) and `V` blows up at the
complex zeros of `Ai`; the folded node is the Weber (parabolic-cylinder) analogue.
The exponent is `Re/Im ∮√Q dτ` around the nearest complex turning point, and the
prefactor is the **Stokes constant**. Exact WKB (Delabaere–Dillinger–Pham, Voros)
is *built* for Schrödinger-type equations — and the Weber inner equation
`ψ''+(c₀−τ²/4)ψ=0` is exactly Schrödinger with a parabolic potential. This is the
most direct, most rigorous attack on the exponent `p` and the constant `C`. **Try:
compute the Borel/Stokes structure of the transition-regime inner equation; the
softening `p<1` is a turning-point coalescence (μ crossing ε^{1/2}).**

**2. Resurgence from the divergent algebraic series (cheap, do this now).** The
sector-scaling series for the canard amplitude *diverges* (factorially). By the
resurgence dictionary, its large-order growth `a_n ~ A·Γ(n+β)/C^n` encodes the
exponentially small `a_min`: `C` is the exponent base, `A` the Stokes
prefactor. So you can read `a_min`'s non-perturbative parameters off the
*perturbative* series — compute ~6–10 terms of the algebraic expansion, fit the
large-order growth, done. Semi-numerical, no connection problem, and it pins `p`
and `C` independently of every other route here. **This is the lowest-effort
high-information experiment.**

**3. Uniform large-order parabolic-cylinder asymptotics (Dunster/Olver).** The
"execute KW2010 but with modern tools" route. `a_min` = recessive/dominant ratio
of the Weber solution of order `a ~ 1/μ`. Dunster's *uniform* asymptotics for
parabolic cylinder functions of large order cover the turning-point region in one
expansion — which is precisely the transition regime that breaks the standard
sector scaling. The exponentially small ratio falls out of the uniform connection
formula. **Try: identify the order `a(μ,ε)` in the transition regime, then read
`a_min` from Dunster's large-`a` PCF connection.**

## Tier 2 — genuinely creative cross-overs

**4. Painlevé II / Riemann–Hilbert.** The canard ↔ Airy is the *linear* shadow of
a *nonlinear* transition; the natural nonlinear generalisation that interpolates
Airy is **Painlevé II** (Hastings–McLeod connects to Airy at one end). If the
folded-node inner problem with the cubic retained is PII-type, then `a_min` is a
**PII connection constant** — and PII's connection problem is *solved* exactly by
Riemann–Hilbert / isomonodromy (Its–Kapaev). That would give `a_min` in closed
form including the exponentially small piece. **Try: check whether the
transition-regime inner equation reduces to PII (or PXXXIV, the Airy-adjacent
one); if so, the answer is in the RH literature.** High payoff, real chance the
reduction holds because the geometry is exactly "fold + slow drift."

**5. Arnold tongues / circle-map renormalisation.** The funnel return map *is* a
circle map; the SAO sectors *are* Arnold tongues; `a_min` is the trajectory's
depth in the **narrowest accessible tongue**. Tongue widths have universal
scaling (Feigenbaum–Kadanoff–Shenker renormalisation at criticality; exponentially
thin tongues for high rotation number away from it). This reframes `a_min` as a
*tongue-width* problem and ties it self-consistently to the staircase exponent `α`
(#6). **Try: extract the circle map from the Poincaré return, identify the deepest
rotation number selected by the global return, apply tongue-width asymptotics.**
Bonus: this is the natural home for the "global return" piece — the return *is*
the circle map.

**6. Nonlinear Landau–Zener.** The funnel passage is a slow sweep through an
avoided crossing (weak vs strong eigendirection). Standard LZ gives an
exponentially small transition `exp(−πΔ²/v)` (`p=1`). But the fold is *nonlinear*
(quadratic), and **nonlinear Landau–Zener** (as in mean-field/BEC sweeps) produces
*anomalous, fractional* exponents — a natural candidate for the observed `p≈0.7`.
**Try: linearise-about-weak-canard into a 2-level LZ with a quadratic nonlinearity;
the nonlinear-LZ exponent should reproduce `p`.** If it lands, `a_min`'s exponent
has a one-line physical origin.

**7. Renormalisation of the funnel.** The per-turn geometric growth `ln R = κμ` is
a *self-similarity*: rescaling by `R` each turn is an RG step. `a_min` is where the
RG flow terminates (inner ↔ global matching). A funnel renormalisation operator
would give `p` as a universal scaling dimension and `C` from the fixed point.
Speculative but the self-similar cascade is manifestly there.

**8. Splitting of separatrices (Lazutkin/Gelfreich/Treschev).** `a_min` = closest
approach to the weak canard = the *exponentially small splitting* between strong
and weak invariant manifolds. The splitting-of-separatrices machinery computes
exactly such `exp(−C/μ)` gaps with prefactors via complex-time singularities — the
rigorous version of idea 1, specialised to manifold splitting. Mature theory,
rarely pointed at folded nodes.

## Tier 3 — wild long-shots (try without worrying)

**9. Multifractal / thermodynamic formalism.** The devil's staircase has a
multifractal spectrum; `a_min` ↔ the most-singular point of the staircase measure.
The pressure function + Legendre transform of the staircase could deliver
`a_min`'s scaling as a thermodynamic quantity. Real link to #6.

**10. Modular / hyperbolic geometry.** Mode-locking's Farey/Stern–Brocot tree is
`PSL(2,ℤ)`; the staircase is the Minkowski `?`-function. `a_min` at the deepest
sector ↔ a cusp depth / geodesic excursion in the modular surface, and
`exp(−C/μ)` ↔ `exp(−hyperbolic depth)`. Wild, but the number theory is genuinely
there.

**11. `q`-deformation / quantum dilogarithm.** `R = e^{κμ}` is a `q` (`q=e^{κμ}`);
the SAO amplitudes form a `q`-series and `a_min ~ q^{s_max}` is a `q`-exponential.
The quantum-dilogarithm/`q`-special-function calculus might sum the funnel exactly.

**12. Tracy–Widom / Airy-kernel edge.** Our Airy structure + "smallest amplitude"
smells like a soft-edge extreme-value problem (Airy kernel, TW left tail
`~exp(−c|s|^3)`). If `a_min` maps to a soft-edge fluctuation, its tail is a TW
object. Long shot but the Airy connection is not nothing.

**13. p-adic / ultrametric tree.** The nested sectors form a tree (Stern–Brocot);
`a_min` lives at a deep node; an ultrametric analysis of tree depth ↔ `a_min`.

**14. Data-driven (PINN / symbolic regression / sparse identification).** Train a
neural-ODE / PINN on the transition-regime inner equation, or run symbolic
regression on `a_min(μ,ε)` over an `ε`-`μ` grid, to *discover* the functional form
(and settle `p`). Pragmatic, and it feeds idea-by-idea falsification.

**15. Computer-assisted proof (CAPD / interval arithmetic).** Not a closed form,
but rigorous `a_min(c)` bounds — enough to *prove* the exponent `p` and validate
whichever analytic route wins.

---

## "Tried some stuff" — three I'd actually start

- **(A) Resurgence on the divergent series (idea 2)** — cheapest, pins `p` and `C`
  from the perturbative coefficients alone. Do this first; it referees the rest.
- **(B) Exact-WKB / complex-time on the Weber inner equation (idea 1/8)** — the
  rigorous exponent: the Stokes constant of `ψ''+(c₀−τ²/4)ψ=0` in the transition
  regime, with the `p<1` softening coming from turning-point coalescence at
  `μ~ε^{1/2}`. This is the principled closed-ish route.
- **(C) Painlevé-II check (idea 4)** — highest payoff per hour: *if* the inner
  equation is PII/PXXXIV, `a_min` is a known RH connection constant and you're
  done. One afternoon to confirm or kill the reduction.

And the unifying realisation: every route splits `a_min` into the same two
factors the audit named — a **transition-regime inner exponential**
(ideas 1–4, 6, 8: exact-WKB / PCF / Painlevé / LZ all attack *this*) ✕ a
**global-return prefactor** (ideas 5, 7, 9: circle maps / RG / multifractal attack
*this*). The cleanest possible outcome would be: inner exponent from exact-WKB or
Painlevé, global prefactor from the circle-map tongue width — and the keystone
`a_min` falls, taking #4, #5, #6 with it.

**First concrete step, today, no theory required:** the `ε`-sweep + resurgence
fit to pin the exponent. Everything above is graded by whether it predicts that
number.
