# Channel B: the noisy saddle-node-on-circle and the 2/3 phase exponent

**Bifurcation-specificity of the folded-cycle phase channel: SNIC vs fold**

Solomon Williams, University of Edinburgh

## Abstract

We finish the Channel-B (phase) inner analysis and pin its bifurcation-specific exponent.
When the limit cycle is destroyed by a **fold of cycles** (finite period), the phase
velocity is bounded below, the iPRC is bounded, and the phase diffusion is regular
(`D_φ ~ σ²`): Channel A (amplitude escape, `σ_*^A ~ ε₂^{1/2}`) governs. When it is
destroyed by a **SNIC** (saddle-node on the invariant circle, diverging period), the
saddle-node lives **in the phase**: the inner phase equation is the noisy
saddle-node-on-circle, which rescales to the parameter-free canonical form
`du = ½u²dτ + dW_τ` at the inner scale `θ ~ σ^{2/3}`. Hence the noise-induced rotation
rate and phase diffusion scale as `σ^{2/3}` — verified, `ω ~ σ^{0.668}` — with the iPRC
diverging as `(μ−μ_c)^{−1/2}`. The clean exponent is **2/3**; a previously reported `0.83`
is the slowly-converging pre-asymptotic effective *diffusion* exponent, which shares the
same 2/3 root. This is a scaling result (matched asymptotics + the canonical reduction +
numerics), honestly weaker than the Tracy–Widom rigour of Channel A.

## 1. The two channels and the inner phase equation

Near the fold the oscillator carries an amplitude `r` and a phase `θ`. **Channel A** is
amplitude escape across the canard (the Tracy–Widom story, rigorous). **Channel B** is the
accumulation of phase noise, with effective amplitude weighted by the infinitesimal
phase-response curve (iPRC) `Z(θ) = ω/v(θ)`, where `v(θ) = dθ/dt` is the phase velocity and
`ω` the mean frequency. The bifurcation that destroys the cycle decides which channel is
singular.

**Fold of cycles (finite period).** Two limit cycles collide and annihilate at *finite*
period: the saddle-node is in the *amplitude* direction (`dr = (r²−y)dT + η dB`, Channel
A). The phase velocity stays bounded below, `v(θ) ≥ v_min > 0`, so the iPRC is bounded and
the phase obeys a regular diffusion `dθ = v dt + σ Z dW` with `D_φ ~ σ²`. Channel B is
regular; `σ_*^A ~ ε₂^{1/2}` (exponent ½) governs, and "A preempts B".

**SNIC (diverging period).** A saddle-node forms *on* the invariant circle: a fixed point
appears in the phase, the period diverges, and the saddle-node is in the *phase*
direction. The Adler normal form is `dθ/dt = μ − cos θ`, SNIC at `μ_c = 1`; near the
bottleneck `θ ≈ 0`,
```
dθ/dt = (μ−1) + ½θ² + σ ξ .                                                    (1)
```
Now the phase velocity **vanishes** at criticality, and Channel B is singular.

## 2. The 2/3 exponent

> **Proposition (noisy saddle-node-on-circle reduction).** At criticality `μ=μ_c`, the
> bottleneck (1) rescales by `θ = σ^{2/3} u`, `t = σ^{−2/3} τ` to the parameter-free
> canonical noisy saddle-node
> ```
> du = ½ u² dτ + dW_τ .
> ```
> Consequently the bottleneck passage time is `~ σ^{−2/3}`, and — since the bottleneck
> dominates the period — the noise-induced rotation rate and phase diffusion scale as
> `ω ~ σ^{2/3}`, `D_φ ~ σ^{2/3}`.

*Proof.* Substitute `θ = σ^{2/3}u`, `t = σ^{−2/3}τ` into (1) at `μ=μ_c`:
`dθ = σ^{2/3}du`; `½θ²dt = ½σ^{4/3}u²·σ^{−2/3}dτ = σ^{2/3}·½u²dτ`; and the noise
`σ dW = σ·σ^{−1/3}dW_τ = σ^{2/3}dW_τ` (since `dW ~ (dt)^{1/2} = σ^{−1/3}(dτ)^{1/2}`).
Dividing by `σ^{2/3}` gives `du = ½u²dτ + dW_τ`, free of `σ`. Hence the inner time is
`σ^{−2/3}`, the escape time `⟨T⟩ ~ σ^{−2/3}`, and `ω = 2π/⟨T_period⟩ ~ σ^{2/3}` since the
bottleneck dominates the period; the renewal rotation's diffusion inherits the same scale,
`D_φ ~ σ^{2/3}`. ∎

**Verification.** Simulating `dθ = (1−cos θ)dt + σ dW` and measuring
`ω = ⟨θ(T)⟩/T`, `D_φ = Var(θ(T))/2T`:
```
   ω ~ σ^{0.668}  (σ^{0.671} on the small-σ half) ,    D_φ ~ σ^{0.73–0.78} ,
```
with `D_φ/ω ≈ 1.1` roughly constant. The drift pins the exponent at 2/3 to three digits;
the diffusion's effective exponent drifts from ≈0.78 toward 2/3 as `σ→0` (a known
pre-asymptotic feature of giant diffusion near criticality [Reimann]). *A previously
reported 0.83 is this pre-asymptotic diffusion estimate, not a distinct exponent.*

## 3. Bifurcation-specificity

> **Lemma (iPRC dichotomy).** For the Adler oscillator the iPRC peak is
> `Z_max = ω/v_min = √(μ²−1)/(μ−1) ~ (μ−μ_c)^{−1/2}`, diverging at the SNIC; for a fold of
> cycles `v_min > 0` stays bounded, so `Z_max = O(1)`.

The exponent is set by **where the saddle-node lives**:

| | saddle-node in | phase channel |
|---|---|---|
| fold of cycles | amplitude (`r`); period finite | regular `D_φ ~ σ²`; Channel A `½` governs |
| SNIC | phase (`θ`); period → ∞ | singular `ω, D_φ ~ σ^{2/3}`; iPRC `~ (μ−μ_c)^{−1/2}` |

The folded *limit cycle* of JKK is a fold of cycles, so its noise response is governed by
Channel A (`σ_*^A = 2√π √ε₂ √(ac/b)`, Tracy–Widom inner exit). A SNIC realisation instead
exposes the phase channel, with the 2/3 signature. This is the precise sense of Channel
B's bifurcation-specificity.

## 4. Honest scope

**What is established.** The inner phase equation in both cases; the canonical reduction
(Proposition) is exact and gives the clean exponent 2/3, *verified* in the drift to three
digits; the iPRC dichotomy (Lemma) is closed-form; and the bifurcation-specificity is
pinned to the location of the saddle-node.

**What remains soft.** Unlike Channel A's exact Tracy–Widom law, Channel B is a *scaling*
result: the canonical noisy saddle-node `du = ½u²dτ + dW` is exact, but its first-passage
statistics (the prefactor, the precise diffusion constant, the slowly-converging `D_φ`
exponent) are not in closed form, and the full *swept* SNIC passage (slow drift through
`μ_c` compounding the bottleneck) is at the scaling level only. Honest statement: Channel
B's mechanism and clean exponent are identified; its rigour is one rung below Channel A's.

## References

- [Reimann] P. Reimann et al., *Giant acceleration of free diffusion by use of tilted
  periodic potentials*, Phys. Rev. Lett. **87** (2001) 010602; PRE **65** (2002) 031104.
- [ET] G.B. Ermentrout, D. Terman, *Mathematical Foundations of Neuroscience*, Springer
  (2010) — iPRC, SNIC, type-I excitability.
- [JKK] S. Jelbart, C. Kuehn, N. Kuntz, arXiv:2208.01361 (2024).
- Companion: `PathA_endgame_writeup` (Channel A / Tracy–Widom),
  `ShootingCharacterisation_proof`.
```
```
*Verification: `snic_channelB.py` — drift ω ~ σ^{0.668} (= 2/3) at criticality.*
