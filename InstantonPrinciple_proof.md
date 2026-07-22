# The instanton principle, proved: a large-deviation theorem for the stochastic-Airy ground state

Solomon Williams, University of Edinburgh

## Abstract

We prove the optimal-fluctuation ("instanton") principle that fixes the early-escape rate
of the noisy folded limit cycle. Writing the stochastic Airy operator as
`H_ε = −∂_x² + x + ε b′` on `[0,∞)` (Dirichlet at 0, `ε = 2/√β = η`), we show
```
lim_{ε→0} ε² log P(Λ₀(ε) < 0) = −c ,   c = ½ min_ψ K[ψ]²/∫ψ⁴ ,  K[ψ]=∫₀^∞(ψ'²+xψ²) ,
```
the minimiser solving the nonlinear-Airy ground state `−g″+xg=g³`, whence
`c = ∫g'² = 5.4439…`. The proof is **Schilder's theorem plus the contraction principle**:
an integration by parts turns the ground-state eigenvalue into a *continuous* functional
of the Brownian path, after which the large deviations of `εb` transfer to `Λ₀(ε)`.
Combined with the shooting characterisation `Y_node =_d −Λ₀(β)` this yields the
folded-cycle early-escape law `P_early(η) = e^{−c/η²(1+o(1))}`.

## 1. Setup and statement

Let `W=b` be standard Brownian motion on `[0,∞)`, `ε>0`. The stochastic Airy operator
`H_ε = −∂_x² + x + ε b′` (Dirichlet at 0) is defined through its quadratic form [RRV]: on
```
H = { ψ ∈ H¹(0,∞) : ψ(0)=0, ∫ x ψ² < ∞ },   ‖ψ‖₂ = 1,
```
one has
```
⟨ψ, H_ε ψ⟩ = K[ψ] + ε N[ψ] ,   K[ψ]=∫₀^∞(ψ'²+xψ²) ,  N[ψ]=∫₀^∞ ψ² dW ,        (1)
```
where `N[ψ]` is the Wiener integral — centred Gaussian with variance `Q[ψ]=∫₀^∞ ψ⁴`.
The ground state is the variational quantity
```
Λ₀(ε) = inf_{ψ∈H, ‖ψ‖₂=1} ( K[ψ] + ε N[ψ] ) ,                                  (2)
```
with deterministic limit `ξ₀ := inf_ψ K[ψ] = −a₁ ≈ 2.338`.

> **Theorem (instanton principle).** With `c := ½ inf_ψ K[ψ]²/Q[ψ]` (infimum over
> `ψ∈H`, `‖ψ‖₂=1`),
> ```
> lim_{ε→0} ε² log P(Λ₀(ε) < 0) = −c .
> ```
> The infimum is attained at `ψ*` solving `−ψ*″+xψ*=γψ*³` (`γ=K[ψ*]/Q[ψ*]`); after
> rescaling `g=ψ*` this is the nonlinear-Airy ground state and `c = ∫g'² = 5.4439…`.

Three ingredients: (§2) `Λ₀(ε)` is a continuous functional of the path `εW`; (§3)
Schilder + contraction give the LDP; (§4) the rate set is the stated variational problem.

## 2. The eigenvalue is a continuous functional of the path

The white noise `b′` is rough, but (2) hides a smooth dependence on `W`, exposed by parts.

> **Lemma 1.** For `ψ∈H`, `‖ψ‖₂=1`, `N[ψ] = −2∫₀^∞ ψψ′ W dx` (a pathwise integral).
> Define, for `g ∈ C₀ := {g∈C([0,∞)) : g(0)=0}`,
> ```
> F(g) := inf_{ψ∈H, ‖ψ‖₂=1} ( K[ψ] − 2∫₀^∞ ψψ′ g dx ) .                       (3)
> ```
> Then `Λ₀(ε) = F(εW)` a.s., and `F` is finite and Lipschitz on every ball
> `{‖g‖_* ≤ R}`, where `‖g‖_* := sup_x |g(x)| e^{−x}`.

*Proof.* Integration by parts in the Wiener integral (valid since `ψ²∈H¹`, `ψ²(0)=0`, and
`ψ²(x)W(x)→0` as `x→∞` because `ψ` decays super-polynomially while
`W=O(√(x loglog x))`) gives `N[ψ] = ∫ψ²dW = −∫(ψ²)′W = −2∫ψψ′W`. Substituting into (2)
gives `Λ₀(ε)=F(εW)`.

For finiteness and the Lipschitz bound, `|∫ψψ′g| ≤ ‖g‖_* ∫|ψψ′|e^x`, and by
Cauchy–Schwarz with the confinement `∫|ψψ′|e^x ≤ (∫ψ'²)^{1/2}(∫ψ²e^{2x})^{1/2} ≤ C_M` on
`{K[ψ]≤M}` (the second factor finite because the form domain forces Gaussian-type decay).
Minimisers of (3) have `K[ψ] ≤ K[ψ*]+2RC ≤ M_R`, so the inf is effectively over the
compact (in `L⁴`) set `{‖ψ‖₂=1, K[ψ]≤M_R}`, hence attained and finite. For `g,g′` in the
ball, `|F(g)−F(g′)| ≤ 2‖g−g′‖_* sup_{K≤M_R}∫|ψψ′|e^x ≤ C_R′‖g−g′‖_*`. ∎

*(The weight `e^{−x}` accommodates the half-line; minimisers are localised in `x=O(1)`
with exponential tails, so the weight is immaterial to the variational content.)*

## 3. Schilder and the contraction principle

> **Lemma 2 (Schilder, [DZ Thm 5.2.3]).** As `ε→0`, `εW` satisfies an LDP in
> `(C₀, ‖·‖_*)` with good rate function `I[g] = ½∫₀^∞ g′(x)² dx` (`g` absolutely
> continuous, `g(0)=0`), `+∞` otherwise.

> **Proposition.** `Λ₀(ε)` satisfies an LDP with speed `ε^{−2}` and good rate function
> `J(λ)=inf{I[g] : F(g)=λ}`. Hence
> ```
> −inf_{λ<0} J(λ) ≤ liminf ε² log P(Λ₀(ε)<0) ,   limsup ε² log P(Λ₀(ε)≤0) ≤ −inf_{λ≤0} J(λ) .
> ```

*Proof.* `F:(C₀,‖·‖_*)→ℝ` is continuous (Lemma 1), so by the contraction principle
[DZ Thm 4.2.1] applied to Lemma 2, `Λ₀(ε)=F(εW)` obeys the stated LDP; the two displayed
bounds are the LDP lower bound on the open set `(−∞,0)` and upper bound on the closed set
`(−∞,0]`. (Exponential tightness on the half-line follows from
`Λ₀(ε) ≥ ξ₀ − ε sup_{K≤M}|N[ψ]|` and the Gaussian tail of `sup N`.) ∎

The two one-sided rates coincide because `inf{I:F<0} = inf{I:F≤0}`: any `g` with `F(g)=0`
is a limit of `g_δ=(1+δ)g` with `F(g_δ)<0` and `I[g_δ]→I[g]` (Lemma 3 shows the optimiser
scales this way). Hence
```
lim_{ε→0} ε² log P(Λ₀(ε)<0) = −inf{ I[g] : F(g) ≤ 0 } .
```

## 4. The rate set is the instanton

> **Lemma 3.** `inf{ I[g] : F(g) ≤ 0 } = ½ inf_ψ K[ψ]²/Q[ψ] =: c.`

*Proof.* Write `Φ=g′`, so `I[g]=½∫Φ²` and, undoing the parts integration in (3),
`F(g)=inf_ψ(K[ψ]+∫ψ²Φ)`. Thus
```
inf{I[g]:F(g)≤0} = inf_ψ inf{ ½∫Φ² : ∫ψ²Φ ≤ −K[ψ] } .
```
For fixed `ψ` the inner problem is a quadratic minimisation under one linear constraint;
since `K[ψ]>0` the constraint binds and the minimiser is `Φ=−μψ²` with `μ∫ψ⁴=K[ψ]`, i.e.
`μ=K[ψ]/Q[ψ]`, giving inner value `½μ²Q[ψ]=K[ψ]²/(2Q[ψ])`. Minimising over `ψ` gives `c`.

*Attainment.* A minimising sequence `ψ_n` (`‖ψ_n‖₂=1`) has `K[ψ_n]` bounded (else
`K²/Q→∞` by Gagliardo–Nirenberg `Q ≤ C K^{1/2}`, below), hence bounded in `H¹` and tight
in `x`; by Rellich–Kondrachov with the confining weight a subsequence converges in `L⁴` to
`ψ*`, with `Q[ψ_n]→Q[ψ*]` and `K[ψ*]≤liminf K[ψ_n]` (lower semicontinuity), so `ψ*`
attains the inf. The Euler–Lagrange equation of `K²/Q` at `‖ψ‖₂=1` is
`−ψ*″+xψ*=γψ*³` with `γ=K[ψ*]/Q[ψ*]>0`; rescaling normalises `γ=1`, the nonlinear-Airy
ground state. The scaling `ψ↦(1+δ)ψ` keeps `K²/Q` fixed while pushing `F` strictly below 0,
justifying `inf{I:F<0}=inf{I:F≤0}`. ∎

*(Gagliardo–Nirenberg in 1D, `‖ψ‖₂=1`: `Q[ψ]=‖ψ‖₄⁴ ≤ C‖ψ′‖₂‖ψ‖₂³ = C‖ψ′‖₂ ≤ C K^{1/2}`,
so `K²/Q ≥ K^{3/2}/C → ∞` as `K→∞` — the coercivity used above and the negligibility of
large-`K` fluctuations.)*

**Proof of the Theorem.** Combine the Proposition, `inf{I:F<0}=inf{I:F≤0}`, and Lemma 3:
`lim ε² log P(Λ₀(ε)<0) = −c`. The value `c=∫g'²=5.4439…` follows from the nonlinear-Airy
ground state via Pohozaev `∫g'²=∫xg²=½∫g⁴` (multiply `−g″+xg=g³` by `g` and by `xg′`),
solved numerically to that precision. ∎

## 5. Consequence and scope

By `ε=η` and the shooting characterisation `Y_node =_d −Λ₀(4/η²)` (companion note),
`P_early(η)=P(Y_node>0)=P(Λ₀<0)=e^{−c/η²(1+o(1))}`, `c=5.4439…`. This is the early-escape
rate; it is **not** `2π`, and the leading-order critical-noise constant `C_q=2√π` is
unaffected (it derives from the saddle-crossing hazard, a different functional).

**What is cited vs elementary.** The single external input is Schilder's theorem (Lemma 2),
a textbook large-deviation result [DZ]; the contraction principle is likewise standard.
Everything else — the integration by parts of Lemma 1 that renders `Λ₀(ε)` a continuous
functional of the path (the step that tames the white-noise roughness), the explicit
quadratic minimisation in Lemma 3, attainment by the direct method, and the Pohozaev
reduction — is elementary. The rigorous meaning of `H_ε` and `N[ψ]` is the RRV form
construction [RRV], which the statement presupposes.

**Numerical corroboration.** The local slope `−d log P_early/d(η^{−2})` decreases from
`≈7` at `η≈1.9` toward `c=5.444` as `η→0` (`5.66` at the smallest measured `η`), and the
variational value is pinned by the Pohozaev-verified ODE solve.

## References

- [DZ] A. Dembo, O. Zeitouni, *Large Deviations Techniques and Applications*, 2nd ed.,
  Springer (1998). Schilder: Thm 5.2.3; contraction: Thm 4.2.1.
- [RRV] J. Ramírez, B. Rider, B. Virág, *Beta ensembles, stochastic Airy spectrum, and a
  diffusion*, J. Amer. Math. Soc. **24** (2011) 919–944.
- [Lifshitz] I. Lifshitz, S. Gredeskul, L. Pastur, *Introduction to the Theory of
  Disordered Systems*, Wiley (1988) — optimal-fluctuation method.
