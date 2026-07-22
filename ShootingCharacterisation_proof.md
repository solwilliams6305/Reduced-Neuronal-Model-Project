# A complete proof of the shooting characterisation: the noisy-Airy inner exit measure is Tracy–Widom

Solomon Williams, University of Edinburgh

## Abstract

We prove that the inner exit measure of the noisy folded limit cycle — the law of the
first-node level `Y_node` of the swept noisy Airy equation — coincides with the law of
minus the ground-state eigenvalue of the stochastic Airy operator, and hence with
Tracy–Widom: `Y_node =_d −Λ₀(β) =_d TW_β`, `η=2/√β`. The proof is elementary given the
random-Schrödinger framework of Ramírez–Rider–Virág (RRV): the swept solution is the
**energy-zero** solution of a random Schrödinger operator `M`, its first zero is the
level at which the truncated `M` has ground state `0`, and that level is read off through
a **monotone spectral function** whose defining event involves only a **deterministic**
translation — so the stationarity and reflection symmetry of white noise close the
identity. The single nontrivial random shift, which obstructed earlier arguments, is
bypassed by Dirichlet domain monotonicity.

## 1. Setup and standing facts

**The two objects.** Fix `η>0`, `β=4/η²`. Let `ξ̃` be a white noise on `ℝ`.

*Object 1 (the inner exit level).* In the slow variable `Y` the noisy canard solves
```
u_YY = (Y − η ξ̃) u .                                                          (1)
```
(Equation (1) holds because `Y=Y₀−T` gives `∂_Y² = ∂_T²`, and the time white noise
becomes a white noise in `Y`.) Let `ψ` be the solution of (1) recessive
(square-integrable) as `Y→+∞`; its largest zero is the **inner exit level** `Y_node`.

*Object 2 (the SAO ground state).* The stochastic Airy operator
```
H_β = −d²/dx² + x + (2/√β) b′(x)    on [0,∞), Dirichlet at 0,                 (2)
```
(`b′` white noise, `2/√β=η`) has a.s. a simple lowest eigenvalue `Λ₀(β)`.

**The operator `M`.** Write (1) as `Mu=0` with the random Schrödinger operator
```
M = −d²/dY² + Y − η ξ̃(Y) ,                                                    (3)
```
so `ψ` is the energy-`0` recessive solution of `M`. For `ℓ∈ℝ` let `M_ℓ` be `M` on
`L²([ℓ,∞))` with Dirichlet at `ℓ`.

We use four standing facts, each classical in its area.

- **(F1) White-noise invariances.** For every *deterministic* `c∈ℝ`, `ξ̃(c+·) =_d ξ̃`
  (stationarity) and `−ξ̃ =_d ξ̃` (reflection). Hence `−η ξ̃(c+·) =_d (2/√β) b′`.
- **(F2) Random Schrödinger operator.** For a.e. realisation and every `ℓ`, `M_ℓ` is
  self-adjoint, bounded below, with discrete spectrum near its bottom and a *simple*
  lowest eigenvalue whose eigenfunction is *strictly positive* on `(ℓ,∞)`; likewise for
  `H_β`. (The RRV construction [RRV]; confining potential `Y→+∞` plus white-noise term.
  See also [FN] for 1D Schrödinger with white-noise potential.)
- **(F3) Weyl limit point at `+∞`.** Since the potential `→+∞`, (1) is limit-point at
  `+∞`: for each energy there is a unique-up-to-scale solution square-integrable near
  `+∞` [CL]. In particular `ψ` exists, unique up to scale.
- **(F4) RRV edge.** `−Λ₀(β) =_d TW_β` [RRV].

## 2. The theorem

> **Theorem (shooting characterisation).** `Y_node =_d −Λ₀(β) =_d TW_β`, `η=2/√β`.
> Moreover, if `u^{Y₀}` is the solution of (1) started at `Y=Y₀` on the canard branch
> (`u_Y/u = −√Y₀`) and `Y_node(Y₀)` its largest zero, then `Y_node(Y₀) → Y_node` a.s. as
> `Y₀→∞`; hence `Y_node(Y₀) →_d −Λ₀(β)`.

The proof is four lemmas; the crux is Lemma 3.

**Lemma 1 (the exit level is well defined).** *A.s. `ψ` exists (F3), is real and not
identically zero, is positive for all large `Y`, and has a largest zero `Y_node∈ℝ`.*

*Proof.* Existence/uniqueness up to scale is F3; take `ψ` real. For `Y` large the
potential `Y−ηξ̃>0`, so by (1) `ψ` is convex where positive; the recessive branch is
positive for all large `Y`. As `Y→−∞` the potential `→−∞` and `ψ` oscillates with
infinitely many isolated zeros (a nontrivial 2nd-order linear ODE solution has only
simple zeros). The zero set is bounded above and discrete, so it has a maximum
`Y_node`. ∎

**Lemma 2 (forgetting: the sweep limit).** *A.s. `Y_node(Y₀) → Y_node` as `Y₀→∞`,
exponentially fast.*

*Proof.* Use the Riccati variable `P=u_Y/u`, solving `P_Y=(Y−ηξ̃)−P²`. On `{Y>0}` the
WKB branches are `P≈±√Y`; linearising about `P=−√Y` gives `∂_P[(Y−ηξ̃)−P²]=−2P≈+2√Y>0`,
so along *decreasing* `Y` the branch `P=−√Y` is attracting. For two solutions
`δ=P^{(1)}−P^{(2)}` obeys `δ_Y=−(P^{(1)}+P^{(2)})δ`, so for `Y` decreasing from `Y₀`,
```
|δ(Y)| = |δ(Y₀)| exp(−∫_Y^{Y₀}(P^{(1)}+P^{(2)})dY′) ≤ |δ(Y₀)| e^{−2∫_Y^{Y₀}√Y′(1+o(1))dY′} .
```
The canard datum `P^{Y₀}(Y₀)=−√Y₀` differs from the recessive value by the WKB error
`δ(Y₀)=O(Y₀^{−1})` (F3); integrating down to any fixed `Y_*` contracts this to
`O(Y₀^{−1} e^{−(4/3)Y₀^{3/2}})`, uniformly. Hence `P^{Y₀}→P_ψ` locally uniformly on
`{Y≤Y_*}`, and the largest zero `Y_node(Y₀)` (where `P→−∞` transversally) converges to
`Y_node`. ∎

**Lemma 3 (monotone spectral function).** *Let `G(ℓ)` be the lowest eigenvalue of `M_ℓ`
(F2). A.s.:*
- *(i) `G` is strictly increasing on `ℝ`;*
- *(ii) `G(Y_node)=0`;*
- *(iii) for every `t∈ℝ`, `{Y_node ≤ t} = {G(t) ≥ 0}`.*

*Proof.* **(i)** For `ℓ₁<ℓ₂`, the form domain of `M_{ℓ₂}` embeds into that of `M_{ℓ₁}`
by extension by zero, with equal forms on the image; by min–max (Dirichlet domain
monotonicity [RS]) `G(ℓ₁)≤G(ℓ₂)`. Strictly: the ground state `φ₁` of `M_{ℓ₁}` is
strictly positive on `(ℓ₁,∞)` (F2), so `φ₁(ℓ₂)>0` and `φ₁` is inadmissible for `M_{ℓ₂}`
(which needs vanishing at `ℓ₂`), so the variational minimum strictly increases.
**(ii)** Restrict `ψ` to `[Y_node,∞)`: there `ψ(Y_node)=0` (Dirichlet), `ψ∈L²`
(recessive decay), `Mψ=0`, and `ψ` has no zero in `(Y_node,∞)` since `Y_node` is the
*largest* zero (Lemma 1). A nodeless `L²` eigenfunction is the ground state (F2), so `0`
is the lowest eigenvalue of `M_{Y_node}`: `G(Y_node)=0`. **(iii)** By (i),(ii),
`G(t)≥0=G(Y_node) ⟺ t≥Y_node`. ∎

**Lemma 4 (deterministic translation).** *For each fixed `t∈ℝ`, `G(t) =_d Λ₀(β)+t`.*

*Proof.* Fix `t`; substitute the *deterministic* shift `Y=t+x`, `x≥0`:
```
M_t = −d²/dx² + (t+x) − η ξ̃(t+x) = [ −d²/dx² + x − η ξ̃(t+x) ] + t =: H′_t + t ,
```
an operator on `[0,∞)`, Dirichlet at `0`. By F1, `ξ̃(t+·) =_d ξ̃` and
`−η ξ̃(t+·) =_d (2/√β) b′`, so `H′_t =_d H_β`; in particular
`Λ₀(H′_t) =_d Λ₀(β)`. Since `M_t = H′_t + t`, `G(t)=Λ₀(H′_t)+t`, whence
`G(t) =_d Λ₀(β)+t`. ∎

**Proof of the Theorem.** Fix `t∈ℝ`. By Lemma 3(iii) (a *pathwise* event identity) and
then Lemma 4 (the *law* of `G(t)` at the fixed level `t`),
```
P(Y_node ≤ t) = P(G(t) ≥ 0) = P(Λ₀(β)+t ≥ 0) = P(Λ₀(β) ≥ −t) = P(−Λ₀(β) ≤ t).
```
As `t` is arbitrary, `Y_node =_d −Λ₀(β)`; by F4, `=_d TW_β`. The sweep statement is
Lemma 2 plus the continuous-mapping theorem. ∎

## 3. Discussion

**What is elementary and what is cited.** Given F1–F4 the proof is elementary: Lemmas 1,
3, 4 use only Sturm theory, Dirichlet domain monotonicity and the substitution `Y=t+x`;
Lemma 2 is a Gronwall estimate on the Riccati. The only heavy input is F2 — the rigorous
meaning of a Schrödinger operator with a white-noise potential and the
simplicity/positivity of its ground state. That is precisely the RRV framework we are
matching to, applied to `M`, which is an SAO *in law* (Lemma 4 at `t=0`); the proof
imports no machinery beyond what the statement already presupposes.

**Why earlier arguments stalled.** A direct "reflection" argument translates by the
*random* level `Y_node`, producing a noise `ξ̃(Y_node+·)` shifted by a noise-dependent
amount — not a clean white noise. Lemma 3 removes this: the event `{Y_node ≤ t}` is
recast as `{G(t) ≥ 0}` with `t` *deterministic*, after which F1 applies verbatim.

**Numerical confirmation.** (i) `Y_node =_d −Λ₀(β)` is confirmed by computing `−Λ₀` as
the smallest eigenvalue of the tridiagonal discretisation of `H_β` and `Y_node` by ODE
shooting — two methods with no shared code path — agreeing in mean/variance/skewness at
`β=1,2,4` and matching the exact `TW_β` densities. (ii) The mechanism of Lemma 3 is
confirmed directly: for fixed realisations, `G(ℓ)` (smallest eigenvalue of the
discretised `M_ℓ`) is monotone increasing and crosses `0` at the recessive solution's
largest zero `Y_node` computed from the same noise.

**Consequence.** The inner exit measure of the noisy folded limit cycle is `TW_β` with
`β=4/η²`, exactly; the early-escape probability is the tail mass `1−F_β(0)`, and the
folded-limit-cycle canard escape lies in the Tracy–Widom / KPZ edge universality class.

## References

- [RRV] J. Ramírez, B. Rider, B. Virág, *Beta ensembles, stochastic Airy spectrum, and a
  diffusion*, J. Amer. Math. Soc. **24** (2011) 919–944.
- [FN] M. Fukushima, S. Nakao, *On spectra of the Schrödinger operator with a white
  Gaussian noise potential*, Z. Wahrsch. **37** (1977) 267–274.
- [CL] E. Coddington, N. Levinson, *Theory of Ordinary Differential Equations*,
  McGraw–Hill (1955), Ch. 9 (limit-point/limit-circle).
- [RS] M. Reed, B. Simon, *Methods of Modern Mathematical Physics IV*, Academic Press
  (1978), §XIII (min–max, Dirichlet–Neumann bracketing).
- [TW] C. Tracy, H. Widom, *Level-spacing distributions and the Airy kernel*, Comm. Math.
  Phys. **159** (1994) 151–174.
```
```
*Companion (mechanism check): `proof_mechanism.py` (G(ℓ) monotone, crosses 0 at Y_node);
dual-method confirmation: `folded_cycle_shooting_duality.py`.*
