# Composition error propagation for the noisy folded limit cycle

**The composed chart map concentrates, θ-uniformly: errors add and the fold renews**

Solomon Williams, University of Edinburgh

## Abstract

We prove the propagation step of Path A: the composed noisy transition map
`Π = Π₃∘κ₂₃∘Π₂∘κ₁₂∘Π₁` across the blow-up charts `K₁` (entry) → `K₂` (inner/fold) →
`K₃` (exit) concentrates around the deterministic JKK passage with errors that **add**
rather than compound, uniformly in the angle `θ∈S¹`. The mechanism is a telescoping
(discrete Gronwall) bound for the transverse deviation together with the **renewal**
supplied by the fold: the entry-to-peel-off map is an exponential contraction, so any
entry deviation is forgotten and cannot accumulate across `K₂`. Consequently the exit law
is the inner exit measure (Tracy–Widom, established separately) transported by the bounded
exit chart, plus sub-dominant Gaussian corrections, and the critical-noise law
`σ_*^A = 2√π √ε₂ √(ac/b)` survives the composition. The argument reduces the full
composition to per-chart inputs that are either standard Berglund–Gentz concentration
(`K₁,K₃`) or already proved (`K₂` inner exit measure, forgetting, θ-uniformity).

## 1. Framework

After blow-up the passage is a finite chain of slow–fast SDE segments. We track the
**transverse coordinate** `x` (the deviation of the noisy trajectory from the
deterministic canard `γ`) across sections `Σ₀ → Σ₁ → Σ₂ → Σ₃`: `Σ₀` the entry to `K₁`,
`Σ₁=Σ₁₂` the `K₁/K₂` overlap, `Σ₂=Σ₂₃` the `K₂/K₃` overlap, `Σ₃` the exit. The
chart-overlap maps `κ_{i,i+1}` (smooth blow-down Jacobians) are absorbed into the
transition maps. Writing `x_i` for the state on `Σ_i` and `x̄_i` for the deterministic
state, each stage is a random map
```
x_i = f_i(x_{i-1}) + ζ_i ,    x̄_i = f_i(x̄_{i-1}) ,                            (1)
```
with `f_i` the deterministic transition (JKK) and `ζ_i` the noise increment of the `i`th
segment. By the strong Markov property the `ζ_i` are conditionally independent given the
`x_{i-1}`. Set `δ_i := x_i − x̄_i` (the accumulated transverse error). We use two
structural inputs proved elsewhere.

> **Input 1 (θ-uniform deterministic bounds).** Each `f_i` is Lipschitz with constant
> `L_i`, and the `L_i` (with the contraction and quasipotential constants) are bounded
> above/below by θ-independent constants depending only on min/max of `a,b,c` on `S¹`.
> *(Proved: θ-uniformity lemma; the maps and derivatives are JKK Theorem 3.2.)*

> **Input 2 (fold renewal / forgetting).** The entry-to-peel-off map (`K₁` then the
> attracting phase of `K₂` up to the fold) is an exponential contraction with Lipschitz
> constant `Λ_c ≤ e^{−(4/3)Y_*^{3/2}}` (contraction rate `2√Y` integrated to the entry
> level `Y_*`), so `Λ_c → 0`; equivalently `L₂ ≤ Λ_c` in (1). *(Proved: the canard is
> attracting — the "forgetting" Lemma R1 of the shooting characterisation; confirmed
> numerically — the peel-off law is invariant to entry offset/noise within the basin.)*

## 2. The propagation lemma

> **Lemma 1 (errors add).** For the chain (1),
> ```
> |δ_m| ≤ (∏_{i=1}^m L_i)|δ₀| + Σ_{j=1}^m (∏_{i=j+1}^m L_i)|ζ_j| .
> ```
> Moreover, if `P(|ζ_j|>h_j | x_{j-1}) ≤ p_j(h_j)`, then for any radii `h_j`,
> ```
> P( |δ_m| > (∏_i L_i)|δ₀| + Σ_j (∏_{i>j} L_i) h_j ) ≤ Σ_{j=1}^m p_j(h_j) ,
> ```
> all constants θ-uniform by Input 1.

*Proof.* From (1), `δ_i = ζ_i + (f_i(x_{i-1}) − f_i(x̄_{i-1}))`, so
`|δ_i| ≤ |ζ_i| + L_i|δ_{i-1}|`. Iterating this scalar recurrence (discrete Gronwall) gives
the deterministic estimate. The event `{|δ_m| > …}` is contained in `∪_j {|ζ_j|>h_j}` (if
every `|ζ_j|≤h_j` the deterministic estimate gives the bound); the union bound and the
conditional tails `p_j`, with the tower property over the Markov chain, give the
probabilistic bound. θ-uniformity is inherited from the θ-uniform `L_i, p_j`. ∎

**Why errors do not compound.** The only blow-up route is `∏L_i ≫ 1`. Two facts forbid
it: the products are θ-uniformly bounded (Input 1; JKK's composed map is a diffeomorphism
with controlled derivatives), and — decisively — the fold factor obeys `L₂ ≤ Λ_c → 0`
(Input 2). Hence every product containing the fold stage is exponentially small. This is
the **renewal**: the fold resets the transverse memory, so deviations entering `K₂` are
annihilated, not amplified.

## 3. Per-chart inputs

> **Input 3 (`K₁`, attracting Fenichel tube; Berglund–Gentz).** On `K₁` the deviation
> from the attracting slow manifold is an OU-type process with θ-uniform contraction rate
> `κ₁`; `P(sup|ζ₁|>h) ≤ C₁(t₁/η²) e^{−κ₁h²/2η²}`, and `f₁` contracts onto `γ`. [BG, Ch. 5]

> **Input 4 (`K₂`, inner exit measure).** The peel-off level on `Σ₂` has the inner exit
> measure `μ_η`: in the inner scaling `(Y,R)=(η^{4/3}y, η^{2/3}r)` it is the law of
> `−Λ₀(β) =_d TW_β`, `β=4/η²`, **independent of the entry** (Input 2); the early-escape
> probability is `1−F_β(0)`. Its tail is the `TW_β` tail, `P(|ζ₂|>h) ≤ p₂(h)`. *(Proved:
> reduction theorem and shooting characterisation.)*

> **Input 5 (`K₃`, exit tube; Berglund–Gentz).** On `K₃` the trajectory, having left the
> canard, follows the deterministic exit solution with a θ-uniform Gaussian tube,
> `P(sup|ζ₃|>h) ≤ C₃(t₃/η²) e^{−κ₃h²/2η²}`, and the exit map `f₃` has θ-uniform Lipschitz
> constant `L₃ ≤ L̄`. [BG Ch. 5; BGK]

## 4. The composition theorem

> **Theorem (composition error propagation).** Under Inputs 1–5, with `δ₀=0`:
> - **(a) concentration.** For radii `h₁,h₂,h₃`,
>   `P(|δ₃| > L̄h₂ + h₃ + Λ_c h₁) ≤ p₁(h₁)+p₂(h₂)+p₃(h₃)`, θ-uniformly; the entry term
>   `Λ_c h₁` is exponentially small.
> - **(b) exit law.** As `η→0`, the exit transverse coordinate converges in law to
>   `f₃(μ_η)` — the inner exit measure transported by the exit chart — uniformly in `θ`;
>   the sub-dominant correction is the `O(η)` Gaussian `ζ₃`.
> - **(c) critical noise.** The escape probability of the whole passage equals the inner
>   escape probability `1−F_β(0)` up to the sub-dominant corrections of (a); hence escape
>   becomes likely at `σ_*^A = 2√π √ε₂ √(ac/b)`, uniformly in `θ`.

*Proof.* Apply Lemma 1 with `m=3` and the per-chart tails (Inputs 3–5). The Lipschitz
products are `∏_{i=1}^3 L_i = L₃L₂L₁ ≤ L̄ Λ_c L₁`, and in the sum
`∏_{i>1}L_i = L₃L₂ ≤ L̄Λ_c` (multiplies `|ζ₁|`), `∏_{i>2}L_i = L₃ ≤ L̄` (multiplies
`|ζ₂|`), empty product `1` (multiplies `|ζ₃|`). With `δ₀=0`,
`|δ₃| ≤ L̄Λ_c|ζ₁| + L̄|ζ₂| + |ζ₃|`, and the union bound gives (a).

For (b): take `h₁,h₃ = O(η)` (the BG tube radii), so `L̄Λ_c h₁` and `h₃` are `O(η)` or
smaller, while `ζ₂` carries the inner scale `O(η^{2/3})` (Input 4). Thus
`δ₃ = L̄ζ₂ + O(η)`, i.e. `x₃ = f₃(x̄₂+ζ₂) + O(η) = f₃(μ_η) + O(η)` in law, uniformly in
`θ`. Because `L₂ ≤ Λ_c → 0`, the entry law `x₀` and its fluctuation `ζ₁` drop out: the
fold renews the chain (Input 2), so the limit is independent of the entry data.

For (c): the passage escapes iff the canard peels off before the fold — an inner event of
probability `1−F_β(0)` (Input 4); by (a) the surrounding charts perturb this only by the
union-bound corrections `p₁+p₃ = O(e^{−c'/η²})`, sub-dominant to `1−F_β(0)`. The
leading-order threshold is unchanged: escape becomes order-one when the integrated inner
hazard `H=η²/4π` is order one, i.e. at `η_*=2√π`, giving `σ_*^A`. All bounds θ-uniform. ∎

## 5. Scope

**What is proved here.** The *propagation* itself: the composed map's transverse error is
the θ-uniformly weighted sum of per-chart fluctuations (Lemma 1); the fold renewal
annihilates the entry contribution (errors do not compound across the non-hyperbolic
chart); the exit law and critical-noise threshold are therefore the inner ones transported
by the bounded exit chart. This is exactly the "error accumulation, not just per-factor"
step that was the last structural gap in Part 2.

**What is imported.** Per-chart concentration on the *normally hyperbolic* charts `K₁,K₃`
is standard Berglund–Gentz sample-path concentration [BG, BGK]; the inner (`K₂`) exit
measure, the forgetting, and the θ-uniformity are proved in the companion notes; the
deterministic composed map and its derivative bounds are JKK Theorem 3.2.

**What remains.** (i) The blow-*down* from the blow-up charts to the physical `(ε₁,ε₂)`
variables (the final desingularisation bookkeeping); (ii) the Channel-B (phase) inner with
its bifurcation-specific exponent. Neither affects the propagation proved here.

## References

- [BG] N. Berglund, B. Gentz, *Noise-Induced Phenomena in Slow–Fast Dynamical Systems: A
  Sample-Paths Approach*, Springer (2006).
- [BGK] N. Berglund, B. Gentz, C. Kuehn, *From random Poincaré maps to stochastic
  mixed-mode-oscillation patterns*, J. Dynam. Diff. Eq. **27** (2015) 83–136; *Hunting
  French ducks in a noisy environment*, J. Diff. Eq. **252** (2012) 4786–4841.
- [JKK] S. Jelbart, C. Kuehn, N. Kuntz, arXiv:2208.01361 (2024) — deterministic backbone,
  Theorem 3.2, charts `K₁–K₃`.
- Companion notes: `ShootingCharacterisation_proof` (inner exit measure, forgetting),
  `InstantonPrinciple_proof` (rate), `PathA_endgame_writeup` (§ θ-uniformity).
```
```
*Renewal corroboration: `renewal.py` — the inner peel-off law (mean −2.0504, std 0.6428,
skew 0.176) is invariant to entry offset δ₀ ∈ {−0.5, 0, 0.5, 1.0} and to small entry
noise.*
