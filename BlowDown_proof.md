# Blow-down to physical variables: the folded-cycle escape theorem

**Assembling the chart results into the original (ε₁,ε₂) system**

Solomon Williams, University of Edinburgh

## Abstract

We carry the chart results back through the geometric blow-up to the physical variables,
completing Path A's leading-order-with-prefactor analysis. The fold-weighted blow-up
`r=ε₂^{1/3}R, y=ε₂^{2/3}Y, t=ε₂^{−1/3}T` maps the physical SDE **exactly** onto the inner
Riccati with effective noise `η=σ/√ε₂` (an Itô change of variables, no anomaly). Blowing
down the chart conclusions then gives, in the original variables: the critical-noise law
`σ_*^A = 2√π √ε₂ G`, the inner exit (peel-off) measure `TW_β` with `β=4ε₂/σ²` at physical
scales `(r,y) ~ (ε₂^{1/3}, ε₂^{2/3})`, the early-escape probability `1−F_β(0)`, and the
early-escape rate `e^{−c ε₂/σ²}` with `c=5.4439…`. The geometry factor `G` is `√(ac/b)`
frozen (`α=2`) or `√(⟨a⟩⟨c⟩/⟨b⟩)` averaged (`α=1`). The noise collapse in `σ/√ε₂` is
verified in the physical variables. This assembles the six companion notes into one
statement.

## 1. The blow-up and the noise covariance

The JKK normal form near the folded-cycle fold is, in physical variables,
```
dr = (b(θ)r² − a(θ)y) dt + σ dW_t ,   dy = −ε₂ c(θ) dt ,   θ̇ = ε₁ ,            (1)
```
with `a,b,c>0` smooth 1-periodic and `0<ε₂≪ε₁≪1`, `ε₁=ε^α`, `ε₂=ε³`. The Krupa–Szmolyan
fold blow-up is the quasi-homogeneous map (weights `(1,2,3)` on `(r,y,ε₂)`) covered by
charts `K₁` (entry), `K₂` (rescaling), `K₃` (exit). The central chart `K₂` is
`ρ=ε₂^{1/3}`:
```
r = ε₂^{1/3} R ,   y = ε₂^{2/3} Y ,   T = ε₂^{1/3} t .                          (2)
```

> **Proposition (exact noise covariance).** Under (2), (1) becomes the inner Riccati
> ```
> dR = (bR² − aY) dT + η dB_T ,   dY = −c dT ,    η = σ/√ε₂ ,
> ```
> exactly (an Itô change of variables on the diffeomorphism (2); the blow-up introduces no
> Itô correction because it is a smooth coordinate change). After the linear rescaling
> absorbing `a,b,c` this is the canonical `dR=(R²−Y)dT+η dB`.

*Proof.* `dr=ε₂^{1/3}dR`; the drift `(br²−ay)dt = ε₂^{2/3}(bR²−aY)·ε₂^{−1/3}dT =
ε₂^{1/3}(bR²−aY)dT`; the noise `σ dW_t = σ(dt)^{1/2}Ẑ = σ(ε₂^{−1/3}dT)^{1/2}Ẑ =
σ ε₂^{−1/6}dB_T`. Dividing `dr=ε₂^{1/3}dR` by `ε₂^{1/3}`: `dR=(bR²−aY)dT+σε₂^{−1/2}dB_T`.
Similarly `ε₂^{2/3}dY=−ε₂ c ε₂^{−1/3}dT` gives `dY=−c dT`. ∎

**Scales.** The chart inner scales `R~η^{2/3}`, `Y~η^{4/3}` blow down, via (2) and
`η=σ/√ε₂`, to the *physical* noise-induced peel-off scales
```
r_peel ~ ε₂^{1/3} η^{2/3} = σ^{2/3} ,   y_peel ~ ε₂^{2/3} η^{4/3} = σ^{4/3} ,
```
independent of `ε₂` — the noise-induced jump scales — while the typical (canard) escape
sits at `y = ε₂^{2/3} a₁`.

## 2. The blow-down theorem

Blowing down the chart conclusions (the reduction, the shooting characterisation, the
instanton rate, the composition propagation, all θ-uniform) through (2) and the Proposition
gives the physical statement.

> **Theorem (Path A, physical variables).** For (1), modulo the per-chart Berglund–Gentz
> tubes and the JKK deterministic backbone:
> - **(i) critical noise.** Escape across the fold becomes likely at
>   ```
>   σ_*^A = 2√π √ε₂ G ,   G = √(a(θ_*)c(θ_*)/b(θ_*))  (α=2, frozen θ)
>                            or √(⟨a⟩⟨c⟩/⟨b⟩)        (α=1, averaged) ,
>   ```
>   uniformly in θ;
> - **(ii) inner exit measure.** The peel-off level is `y_peel = ε₂^{2/3} Y_node` with
>   `Y_node =_d TW_β`, `β = 4ε₂/σ²`; equivalently `P_early = 1−F_β(0)`;
> - **(iii) early-escape rate.** As `σ/√ε₂→0`, `P_early = exp(−c ε₂/σ² (1+o(1)))`,
>   `c=5.4439…`.
>
> All three depend on `(σ,ε₂)` only through `η=σ/√ε₂` at leading order.

*Proof.* By the Proposition the `K₂` dynamics is the canonical inner Riccati with
`η=σ/√ε₂`. (ii) is the reduction theorem (`Y_node =_d TW_β`, `β=4/η²=4ε₂/σ²`) blown down by
`y=ε₂^{2/3}Y`; (iii) is the instanton rate `e^{−c/η²}=e^{−cε₂/σ²}`. For (i): escape is
order-one when the integrated inner hazard `H=η²/4π` is order one, `η_*=2√π`, so
`σ_*=2√π√ε₂` in the canonical chart; restoring `a,b,c` multiplies by `G`. The
composition-propagation theorem carries these inner laws through `K₁,K₃` with θ-uniform
sub-dominant Gaussian corrections (the fold renewal annihilating the entry deviation); the
deterministic geometry is JKK Theorem 3.2, whose α-dichotomy gives the two cases of `G`.
Every estimate is θ-uniform by the θ-uniformity lemma. ∎

**Verification of the collapse.** Simulating the *physical* system (1) (`a=b=c=1`) through
the fold and rescaling `Y=y_peel/ε₂^{2/3}`, the early-escape law depends only on
`η=σ/√ε₂`:
```
   η=1.5:  (ε₂,σ)=(0.05,0.335): P_early=0.117, ⟨Y⟩=−1.20   (0.20,0.671): 0.113, −1.21
   η=1.0:  (ε₂,σ)=(0.05,0.224): P_early=0.011, ⟨Y⟩=−1.60   (0.20,0.447): 0.010, −1.60
```
a clean collapse across a factor 4 in `ε₂` (`blowdown.py`).

## 3. Uniformity and the assembled statement

**Uniformity in (ε₁,ε₂).** The `ε₂`-dependence enters only through the blow-up radius
`ρ=ε₂^{1/3}` and the chart noise `η=σ/√ε₂` (Proposition), so the chart estimates —
θ-uniform and ε₂-uniform once written in `η` — transfer verbatim. The `ε₁=ε^α` dependence
enters only through `θ̇=ε₁`: for `α=2` the angle is frozen and `G` is local; for `α=1` the
angle is fast and Khasminskii averaging replaces `a,b,c` by `⟨a⟩,⟨b⟩,⟨c⟩` with `O(1/ε₁)`
corrections — exactly JKK Theorem 3.2's deterministic α-dichotomy, now carrying the noise.

**The assembled Path A.** The theorem is the capstone of the chain: *Cole–Hopf* (exact
Itô) ⇒ *reduction* to the stochastic Airy operator ⇒ *shooting characterisation*
(`Y_node =_d −Λ₀ =_d TW_β`, proved) ⇒ *instanton* rate `c=5.4439` (proved) ⇒ *composition
propagation* (θ-uniform, fold renewal, proved) ⇒ *blow-down* (this note) to the physical
critical-noise, exit-measure, and rate laws. The folded-limit-cycle canard escape is
thereby placed in the Tracy–Widom / KPZ edge universality class, in the original variables.

## 4. Honest scope

The blow-down is rigorous *given* the chart results: the Proposition is an exact change of
variables, the scale blow-downs are algebra, and the assembly is the composition theorem
plus JKK's deterministic backbone. The remaining softness is inherited, not introduced: the
per-chart `K₁,K₃` tubes are cited Berglund–Gentz (standard for hyperbolic charts), and
Channel B is scaling-level. What this note adds is the final desingularisation — turning
the `O(1)` chart statements into the physical-`(σ,ε₂)` theorem — which was the last named
structural item of Path A.

## References

- [KS] M. Krupa, P. Szmolyan, *Extending GSPT to nonhyperbolic points*, SIAM J. Math.
  Anal. **33** (2001) 286–314 (the fold blow-up).
- [JKK] S. Jelbart, C. Kuehn, N. Kuntz, arXiv:2208.01361 (2024), Theorem 3.2.
- [BG] N. Berglund, B. Gentz, *Noise-Induced Phenomena in Slow–Fast Dynamical Systems*,
  Springer (2006).
- Companion notes: `ShootingCharacterisation_proof`, `InstantonPrinciple_proof`,
  `CompositionPropagation_proof`, `ChannelB_SNIC_proof`, `PathA_endgame_writeup`.
```
```
*Verification: `blowdown.py` — physical normal form collapses in η=σ/√ε₂.*
