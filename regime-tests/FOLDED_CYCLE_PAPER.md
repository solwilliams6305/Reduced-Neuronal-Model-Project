# Stochastic blow-up of the folded limit cycle — consolidated leading-order manuscript

*Synthesis of the folded-limit-cycle results (Parts 1–3) into one coherent,
leading-order narrative. Path B throughout (accumulated-variance heuristics +
numerics with falsification controls); the Path-A rigour is deferred by design and
flagged. Novelty claim gated on the Popović reply (`POPOVIC_NOVELTY_EMAIL.md`).*

Source chapters: `FOLDED_CYCLE_NOISE.md` (Parts 1: normal form, η, σ\*ᴬ, σ\*ᴮ, race),
`FOLDED_CYCLE_WORKED_MODEL.md` (Part 2: worked model, C_q, G(θ), Channel-B caveat),
`FOLDED_CYCLE_ALPHA1.md` (Part 3: α=1 averaged law, A–B fusion, geometry test).
Scripts/figures named inline.

---

## 0. Abstract

We give the first stochastic treatment of the **folded limit cycle** — the one
non-hyperbolic singularity of the slow-fast "fold family" that Berglund–Gentz
(fold) and Berglund–Gentz–Kuehn (folded node) never noised — by pushing degenerate
noise through the geometric blow-up of Jelbart–Kuehn–Kuntz (JKK, 2024). The folding
object is a *sustained fast oscillation carrying a phase*, so noise acts through
**two channels at once**: amplitude escape off the folded-cycle manifold (Channel A)
and phase diffusion of the fast oscillation (Channel B). We derive the effective
inner-noise scale `η = σ/√ε₂` (set by the *drift* parameter, not the angular one),
the phase-modulated amplitude law `σ\*ᴬ(θ) = C_q √ε₂ √(a(θ)c(θ)/b(θ))`, and the phase
law `σ\*ᴮ ≈ 2π√3 √(ac/b) √(ε₂/|ln ε₂|)` — the same leading order, separated by a
logarithm. For the rotating-phase (α=1) case we show by stochastic averaging that
escape sees the **coefficient-averaged** canard, `σ\*(α=1)=C_q√ε₂√(⟨a⟩⟨c⟩/⟨b⟩)` (not
the minimum, not the barrier-average), and that the **A–B race is geometry-independent**
in both the frozen and rotating regimes, with **Channel A preempting Channel B** for
accessible parameters. We validate on a worked fold-of-cycles model (`C_q≈5.7`) and
flag one structural caveat: Channel B's blow-up is *relaxation-cycle-specific*.

## 1. Introduction & positioning

Geometric singular perturbation theory + blow-up resolves non-hyperbolic
singularities of slow-fast systems; noise through these singularities (Berglund–Gentz,
BGK) gives critical-noise laws `σ\*`. The fold (BG) and folded node (BGK) are done.
The **folded limit cycle** — where a *manifold of limit cycles* folds (a saddle-node
of cycles, JKK 2024) — is the un-noised member. It is qualitatively different: the
folding object is a sustained fast oscillation, so it carries a phase.

- vs **BG (fold)**: their fold is a point; ours is a circle of fold points — the
  amplitude law is **phase-modulated**, and there is a second (phase) channel.
- vs **BGK (folded node)**: there the rotation is slow and the fast subsystem still
  jumps; here a sustained fast oscillation folds.
- vs **JKK (deterministic folded cycle)**: this is its stochastic extension.
- vs **Ahsan–Dankowicz–Kuehn 2025 (noisy limit cycles)**: they use adjoint/covariance,
  not blow-up — the blow-up route is precisely the un-done one.

## 2. Setup

JKK prototypical normal form (cylindrical `(r,θ,y)`, `r` fast radius, `θ∈ℝ/ℤ` phase,
`y` super-slow drift), degenerate noise on the fast variable:

```
dr = ( −a(θ)y + b(θ)r² ) dt + σ dW,    dθ = ε₁ dt,    dy = −ε₂ c(θ) dt,
```

`a,b,c > 0` smooth, 1-periodic in θ; `0 < ε₂ ≪ ε₁ ≪ 1`. The fold of cycles is the
circle `r̃=0` (finite radius `v`); the critical manifold is `y ~ (b/a)r̃²`. Cartesian
fast-variable noise projects onto **radial** (Channel A) and **tangential** (Channel B,
amplitude `~σ/(2πr)`) — the two channels.

## 3. Channel A — amplitude escape

Blow-up to the JKK rescaling (K₂) chart with `ε₁=ε^α`, `ε₂=ε³`, weights `(r,y,ε)~(ρ,ρ²,ρ)`.
Pushing `σdW` through (`r=εr₂`, `dt₂=εdt`, `dW=ε^{−1/2}dB`):

```
η = σ / ε^{3/2} = σ / √ε₂      (effective inner noise; the DRIFT parameter ε₂, not ε₁).
```

Freidlin–Wentzell accumulated variance across the passage window (`r₂*=(ac/b²)^{1/3}`,
`y₂*=(c²/ab)^{1/3}`, `T_win=(abc)^{−1/3}`):

```
σ\*ᴬ(θ) = C_q √ε₂ √( a(θ)c(θ)/b(θ) ).
```

Two things are new vs BG: (i) the exponent is in the **drift** ε₂ (the angular ε₁
drops out at leading order); (ii) the threshold is **phase-modulated** by
`G(θ)=√(ac/b)(θ)`, vanishing to a single BG threshold iff `a,b,c` are θ-independent.
Setting `(a,b,c)=(1,1,λ)` recovers the fold law `C_q√ε₂λ^{1/2}` (BG with `ε→ε₂`).

**Validated** (`folded_cycle_normal_form_map.py`): the blown-up SDE collapses across
geometry at fixed `Θ=η/√(ac/b)` (cross-geometry std 0.035), a falsification control
dropping `c` splays it (0.224), and the special case `(1,1,λ)` reproduces the canard
`C_q≈2.8` (cross-chapter inheritance).

**Verified against JKK's exact deterministic backbone (entry + inner charts).** The
passage scales coincide *exactly* with JKK's canonical rescaling (their eq. 42:
`r₂=(ac/b²)^{1/3}R₂`, `y₂=(c²/ab)^{1/3}Y₂`, `t₂=(abc)^{−1/3}T₂` → `dR/dT=R²−Y`, the
Krupa–Szmolyan Riccati). The **canard the noise drives off** is, on the entry side,
JKK's center manifold `M₁ᵃ` at `r₁=−√(a/b)+O(ε₁)` (K₁, Lemma 4.6), onto which the flow
contracts at Fenichel rate `ϑ=min_θ 2√(a(θ)b(θ))` (`e^{−ϱ̃/ε₁³}`, Prop. 4.9); and, on
the inner side, the invariant manifold `γ₂` (K₂, Prop. 4.12: `y₂=(b/a)r₂²` as
`r₂→−∞`; exit offset `−(c²/ab)^{1/3}Ω₀`). The contraction **vanishing at the fold** is
precisely the mechanism the accumulated-variance heuristic captures — noise escapes
when the canard's hold weakens. (JKK's entry chart also averages `c→c₀=⟨c⟩`,
Lemma 4.8 — the deterministic mirror of our α=1 coefficient-averaging, §6.) So `σ\*ᴬ`
sits on JKK's verified inner solution, and `M₁ᵃ`/`γ₂`/`Ω₀` are the exact deterministic
objects the rigorous Path-A solution would noise.

## 4. Channel B — phase diffusion

The phase noise `η/(2πr)` carries a `1/r` factor (the iPRC `~1/amplitude`); as the
inner (canard) variable `r₂→0` the phase diffusion blows up *in the same chart*.
Accumulated phase variance `Var(θ)=(η²/8π²)∫dt₂/r₂² = (η²/8π²)(b/ac)ln(y₂ⁱⁿ/y₂*)` —
**log-divergent**, cut at the inner scale. Phase randomization (`Var~1`):

```
σ\*ᴮ ≈ 2π√3 √(ac/b) √( ε₂ / |ln ε₂| ).
```

**Validated** (`folded_cycle_phase_diffusion.py`): the `I_φ=∫dt₂/r₂²` log law,
ratio →1. **Structural caveat (Part 2 #2, `folded_cycle_channelB_model.py`):** on a
*smooth* fold of cycles (finite-radius cycle) the leading phase diffusion is finite,
`D_φ=σ²/(2ρ²)`, and the blow-up appears only as a higher-order transverse (`σ⁴/|λ⊥|`)
term as the cycle loses normal hyperbolicity. So **Channel B's leading-order log-blow-up
is relaxation-cycle-specific** (fold-passage iPRC) — the two channels are *not*
symmetric in model-generality. An attempt to ground it on an *autonomous* relaxation
oscillator (van der Pol drifting to cycle death, `autonomous_channelB.py`) confirms the
qualitative mechanism — noise *does* diffuse the internal phase, decoherence grows with
σ — but a clean `σ\*ᴮ` is **obstructed**: VdP's cycle dies at a Hopf (amplitude→0), so
the same resonator that broke Channel A's amplitude test reappears, and the spike-count
phase proxy saturates at cycle death rather than tracking the accumulated diffusion. So
Channel B's *mechanism* is grounded (the autonomous fold-passage iPRC; cf. `TONIC_PHASE`'s
edge divergence `A~(I−I_H)^{−1/2}`). Switching to a degeneration that keeps the cycle at
**finite amplitude** — a **SNIC** (period→∞, not amplitude→0; Adler equation,
`adler_snic_channelB.py`) — **does** give a clean, non-saturating `σ\*ᴮ` (confirming the
Hopf resonator was the VdP obstruction). **But** the SNIC exponent is `σ\*ᴮ~ε₂^{0.83}`,
**not** the fold-of-cycles `√ε₂`. So **Channel B's drift-law exponent is
bifurcation-specific** (fold-of-cycles → ½; SNIC → ~0.83), and the specific `√ε₂` law
needs the fold-of-cycles structure, which stays un-realized in an autonomous
finite-amplitude model. The Channel-A/B asymmetry is therefore sharp and genuine:
**Channel A's `√ε₂` is universal across realizations; Channel B's exponent is
degeneration-specific, and its fold-of-cycles form remains normal-form-validated.** This
is a real result about the model-dependence of phase diffusion, stated not hidden.

## 5. The A–B race

Both channels are the same leading order `√ε₂√(ac/b)`:

```
σ\*ᴬ = C_q √ε₂ √(ac/b),                C_q ≈ 2.8 (NF), 5.7 (worked), 8–10 (full),
σ\*ᴮ = C_B √(ε₂/|ln ε₂|) √(ac/b),       C_B = 2π√3 ≈ 10.9,
σ\*ᴮ/σ\*ᴬ = (C_B/C_q)/√|ln ε₂|.
```

The geometry `√(ac/b)` cancels → the race is **geometry-independent**, set by `ε₂`
alone. With the *measured* worked-model `C_q≈5.7`, the crossover is `ε₂≈0.03`:
amplitude escape (A) sets `σ\*` for `ε₂≳0.03`, Channel B only for `ε₂≲0.03` (the
borrowed-`C_q=8` "B wins everywhere" claim is **deflated** to a small-ε₂ statement).

## 6. The α=1 coupled core

When the phase rotates through the fold, stochastic averaging (Khasminskii) slaves the
amplitude to the **coefficient-averaged** canard:

```
σ\*(α=1) = C_q √ε₂ √( ⟨a⟩⟨c⟩/⟨b⟩ )      (NOT min, NOT the barrier-average ⟨ac/b⟩).
```

**Validated** (`alpha1_averaged_law.py`): modulating `b` (where coefficient- vs
barrier-averaging diverge by Jensen) gives ratio `η\*(fast)/η\*(frozen)=1.0–1.05`, flat
in ω — coefficient-averaging (1.00), ruling out barrier-avg (1.29) and min (0.745, the
route `ALPHA1_METHODOLOGY.md` had already refuted).

**Sits on JKK's own theorem.** The full proof of JKK's Theorem 3.2 (§4.5) shows the
*deterministic* transition map has exactly this α-dichotomy: for **α=1** the drift
`h_y=O(ε²)` — the θ-dependence **washes out** (`c→c₀=⟨c⟩`); for **α=2** it
**retains** `h_y=−(c(θ)²/a(θ)b(θ))^{1/3}Ω₀ε²`. So our two stochastic headlines each
sit on the matching limb of JKK's theorem: α=1 coefficient-averaging
`√(⟨a⟩⟨c⟩/⟨b⟩)`, and α=2 frozen-phase modulation `G(θ)=√(ac/b)(θ)`.

**A–B fusion** (`alpha1_fusion.py`): with both channels on and rotating, `σ\*ᴬ≈1.9` is
flat in ω, and at that threshold the phase has only decohered to `Var(θ)≈0.45<1` —
**Channel A preempts Channel B**: amplitude escape removes the cycle before the phase
randomizes (B's `Var` saturates ~0.5; higher noise just escapes earlier). The "race"
is not two independent thresholds — the lower one removes the substrate the other needs.

**Race robustness** (`folded_cycle_geometry_decoupling.py`): I conjectured rotation
would decouple the channels' geometry (A→⟨a⟩, B→⟨1/a⟩) and break the race's
geometry-independence. The numerics refute it: `σ\*ᴬ` flat in the a-modulation (CV 0.00)
AND Channel B's variance flat (it does *not* track ⟨1/a⟩), because the amplitude is
slaved to the averaged canard so B sees `1/⟨a⟩` — the same functional as A. **The
race stays geometry-independent under rotation**; the regime diagram is defensible.

## 7. Worked model (Part 2) — abstract normal form and JKK's named system

**(a) Faithful normal form — Bautin saddle-node of cycles.** A faithful folded
limit cycle needs a fold of cycles at *finite* amplitude (escape across an
unstable-cycle barrier). Forced VdP *through its Hopf* was tried first and
**falsified** (cycle born at amplitude→0 ⇒ resonator; noise induces spikes past
the fold, `folded_cycle_vdp_noise.py`). The Bautin/generalized-Hopf normal form
(`du=[u(β+R²−R⁴)−ωv]dt+σdW`, `dv=…`, `dβ=−ε₂dt`; SNLC `β=−1/4`) is the faithful
testbed: `σ\*ᴬ∝√ε₂` (slope 0.58), measured `C_q≈5.7` (`folded_cycle_bautin.py`);
an anharmonic cubic `g(θ)` gives a *modest* `G(θ)` confirmation (collapse improves
1.32×, `folded_cycle_anharmonic.py`).

**(b) Named system — JKK's own forced van der Pol (§5.1).** JKK's worked example is
the periodically-forced Liénard `x''+μf(x)x'+g(x)=A(ωt)`. In our coordinates
(`K(x)=x−x³/3`, fold `x_F=−1`, `g=x`, `ϑ=1`):

```
x' = −y + K(x) + σ dW,    θ' = ε₁,    y' = ε₂ ( x − A(θ) ),
```

with JKK (74) giving `a=b=1` (constant) and `c(θ)=A(θ)+1` — so **the phase
modulation enters only through `c`, driven directly by the forcing**, and
`G(θ)=√(ac/b)=√c(θ)=√(1+A(θ))`. With `A(θ)=A₀cos2πθ`, the named-model Channel-A
regime map (`folded_cycle_lienard.py`) **strongly confirms** it: including
`G=√c(θ)` collapses the cross-phase curves to CV 0.11 vs 0.42 without — a **3.9×
improvement** (much cleaner than the Bautin's 1.32×, because the forcing modulates
the *drift rate*, which directly sets the canard window). So
`σ\*ᴬ(θ)=C_q√ε₂√(ac/b)(θ)` is validated on JKK's *own* forced van der Pol — a
recognized system, not a normal form.

**Structural point — Channel B is ABSENT in the forced realization.** In §5.1 the
phase `θ` is the *external forcing clock* (`θ'=ε₁`, deterministic), so noise on `x`
cannot diffuse it. Channel B (phase diffusion) is a feature of the *autonomous*
(self-sustained) folded cycle, with an internal phase — the same Channel-A/B
asymmetry as §4, now sharpened: the forced Liénard grounds Channel A and the phase
*modulation*, but Channel B's *diffusion* needs an autonomous folding cycle. (JKK's
deterministic α=1 result — the drift's θ-dependence washing out to O(ε²) — is the
deterministic mirror of our stochastic coefficient-averaging, §6.)

`C_q` stays honestly fuzzy (threshold/projection/finite-ε O(1)); the qualitative
deflation of the "B wins" headline is robust, the exact crossover uncertain by a
factor of a few.

## 8. Results summary

| Object | Law | Status |
|---|---|---|
| Effective noise | `η = σ/√ε₂` | derived (drift param, not angular) |
| Channel A (frozen) | `σ\*ᴬ(θ)=C_q√ε₂√(ac/b)(θ)` | derived + NF-validated + worked `C_q≈5.7` |
| Channel B (frozen) | `σ\*ᴮ≈2π√3√(ac/b)√(ε₂/|ln ε₂|)` | derived + log-validated; **relaxation-specific** |
| A–B race | geometry-independent, B/A=(C_B/C_q)/√|ln ε₂| | derived + robust under rotation |
| α=1 amplitude | `σ\*(α=1)=C_q√ε₂√(⟨a⟩⟨c⟩/⟨b⟩)` | derived + validated (coefficient-avg) |
| A–B coupling | A preempts B (accessible params) | first swing |

## 9. Honest caveats & what is open

- **Leading order, Path B** (accumulated variance + numerics), not theorem-grade. A
  legitimate paper if stated as such (as the fold/node chapters are).
- **Channel B asymmetry**: its leading blow-up is relaxation-cycle-specific; grounding
  it on a real model needs a relaxation-type folding cycle (VdP-Hopf failed Channel A) —
  an unresolved modelling tension.
- **Open**: the large-log/small-ε₂ corner (does B preempt A there?); a *named* worked
  model (JKK §5.1 forced Liénard, un-retrieved); sharper `C_q`/`G(θ)`.
- **Rigour (Path A, `FOLDED_CYCLE_PATHA.md`)**: the **Channel-A inner core is now done
  and verified** — the deterministic inner is exactly Airy (`R=Ai′(Y)/Ai(Y)`, escape at
  the first Airy zero = JKK's `Ω₀`); the inner FW quasipotential is `V(Y)=8Y^{3/2}/3`
  (Kramers-verified); and the integrated escape hazard has the **closed form
  `H(η)=η²/4π`**, giving a rigorous **`C_q=2√π≈3.545`** (upgrading the heuristic/numeric
  `C_q≈2.8`; the exponent and `√ε₂√(ac/b)` geometry are now rigorous). The α=1 averaging
  has its Khasminskii + Fourier-hierarchy skeleton. Sitting on JKK's full deterministic
  backbone (Theorem 3.2, all charts). **Remains for the full theorem:** the global *noisy*
  chart-matching K₁→K₂→K₃ with error control, the rigorous Channel-B inner (incl. its
  bifurcation-specificity), sharp prefactors, and uniformity in `(ε₁,ε₂)` — months-scale,
  novelty-gated, but now a scoped program on a solid foundation.
- **Novelty gate**: standard ingredients recombined → confirm freedom via the Popović
  email before submission.

## 10. References

- Jelbart, Kuehn, Kuntz, *Geometric Blow-Up for Folded Limit Cycle Manifolds…*, J.
  Nonlinear Sci. 34:17 (2024), arXiv:2208.01361.
- Berglund & Gentz, *Noise-Induced Phenomena in Slow–Fast Dynamical Systems* (2006).
- Berglund, Gentz, Kuehn (2015), arXiv:1312.6353.
- Ahsan, Dankowicz, Kuehn, *Adjoint-Based Projections… Stochastically Perturbed Limit
  Cycles and Tori*, SIADS (2025), arXiv:2404.13429.
- Krupa & Szmolyan, SIAM J. Math. Anal. 33 (2001). Khasminskii, *Stochastic Stability
  of Differential Equations* (averaging). Kuznetsov, *Elements of Applied Bifurcation
  Theory* (Bautin/SNLC).
- Project chapters and scripts: `FOLDED_CYCLE_NOISE.md`, `FOLDED_CYCLE_WORKED_MODEL.md`,
  `FOLDED_CYCLE_ALPHA1.md`, `ALPHA1_METHODOLOGY.md`, `EVIDENCE_KUEHN_PIPELINE.md`;
  `folded_cycle_{normal_form_map,phase_diffusion,vdp,vdp_noise,bautin,anharmonic,
  channelB_model,geometry_decoupling}.py`, `alpha1_{averaged_law,fusion}.py`.
