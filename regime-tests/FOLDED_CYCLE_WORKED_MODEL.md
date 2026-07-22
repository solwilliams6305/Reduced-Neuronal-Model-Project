# Folded limit cycle — Part 2: the worked model (and what it does to the headline)

**Status:** Channel-A law validated on a real model with a *genuine fold of limit
cycles*, with a measured prefactor. The first VdP testbed was **falsified** by its
own numerics (a documented negative result); the faithful Bautin fold-of-cycles
model gives `σ_*^A ∝ √ε₂` and a **measured `C_q ≈ 5.7`** — which **tempers** the
"Channel B wins for full models" headline of `FOLDED_CYCLE_NOISE.md`.

**Why Part 2:** everything in `FOLDED_CYCLE_NOISE.md` is normal-form with a
*borrowed* `C_q ≈ 8–10` (from FHN). A real model (i) confirms the law survives
finite ε on an actual system, (ii) fixes a genuine `C_q`, and (iii) the A-vs-B
headline depends entirely on that number.

Companion scripts: `folded_cycle_vdp.py` (deterministic foundation),
`folded_cycle_vdp_noise.py` (the falsified VdP attempt),
`folded_cycle_bautin.py` (the faithful measurement).
Figures: `figures/folded_cycle_vdp_deterministic.png`,
`figures/folded_cycle_bautin_sigmastar.png`.

---

## 1. First attempt — forced Van der Pol — and why it failed (kept on purpose)

Three-timescale forced VdP, `dx=(x−x³/3−w)dt+σdW`, `dw=ε₁(x−a)dt`, `da=ε₂dt`,
drifting `a` up through the cycle fold at `a=1`. Deterministically this works
(`folded_cycle_vdp_deterministic.png`): the relaxation cycle holds O(4) amplitude
across the canard plateau, then collapses at `a=1`.

**But the noisy Channel-A measurement falsified the testbed.** The early-death
margin `⟨m⟩=⟨1−a_death⟩` came out **decreasing** with σ (and nonzero at σ→0),
the opposite of escape. Mechanism: VdP's cycle is born at a **Hopf** (amplitude
→ 0), so near the fold the cycle is tiny and the fixed point is a **resonator** —
degenerate noise *induces* spikes *past* the fold rather than knocking the
trajectory off the cycle early. The amplitude-→0 (Hopf-degenerate) birth is **not**
a clean fold-of-cycles, and contaminates Channel A with noise-induced spiking.

This is exactly the faithfulness caveat flagged at the outset, now confirmed
empirically. Per project discipline (test the conjecture; report what the numerics
decided) the failed script is retained as `folded_cycle_vdp_noise.py`.

## 2. The faithful model — a saddle-node of limit cycles (Bautin)

JKK's "regular fold of limit cycles" has the two cycle branches (attracting `S₀ᵃ`,
repelling `S₀ʳ`) coalescing at **finite radius** `v` — a saddle-node of limit
cycles (SNLC), *not* a Hopf. The minimal faithful realization is the Bautin /
generalized-Hopf normal form, planar fast oscillator with degenerate noise on the
fast Cartesian variable `u`:

```
du = [ u (β + R² − R⁴) − ω v ] dt + σ dW,     R² = u² + v²
dv = [ v (β + R² − R⁴) + ω u ] dt
dβ = −ε₂ dt                                    (super-slow drift)
```

For `β ∈ (−1/4, 0)`: a stable FP at `ρ=0`, an **unstable cycle** `ρ_−` and a
**stable cycle** `ρ_+`. They annihilate at the **SNLC** `β_SN = −1/4`,
`ρ_SN = 1/√2`. Channel-A escape is clean: noise knocks the trajectory off the
stable cycle `ρ_+` across the unstable-cycle barrier `ρ_−` into the FP — the cycle
dies *early*, at `β > β_SN`. Finite-amplitude fold ⇒ **no Hopf resonator**.

**Local fold normal form** (`u = ρ−ρ_SN`, `ν = β−β_SN`): with
`F=ρ(β+ρ²−ρ⁴)`, `F_β=ρ_SN=0.707`, `½F_ρρ = 3ρ_SN−10ρ_SN³ = −√2`,

```
u' = −√2 u² + (1/√2) ν,   ν' = −ε₂,
```

the canard/Krupa–Szmolyan fold ⇒ JKK coefficients `a=1/√2`, `b=√2`, `c=1`, hence
`G = √(ac/b) = 1/√2 ≈ 0.707`. Predicted `σ_*^A = C_q √ε₂ · 0.707`.

## 3. Results (`folded_cycle_bautin.py`)

Observable `m = β_death − β_SN` (early-death margin), `β_death` = β at which the
amplitude first collapses below `ρ=0.40`. Sweep `σ` × `ε₂ ∈ {1,2,4}×10⁻³`, `N=50`.

```
  ε₂        σ_*     σ_*/√ε₂    C_q = σ_*/(√ε₂·G)
  1e-3     0.120     3.79          5.36
  2e-3     0.178     3.98          5.63
  4e-3     0.268     4.24          6.00

  fitted slope  d ln σ_* / d ln ε₂ = 0.58   (predict 0.50)
  model C_q (mean) = 5.7
  ω spot-check (ε₂=2e-3, σ=0.16): ⟨m⟩ = 0.0004, 0.0096, 0.0052 for ω=0.5,1,2  (flat)
```

- **`σ_*^A ∝ √ε₂` confirmed** on a real fold-of-cycles: the `⟨m⟩(σ)` curves for the
  three ε₂ collapse against `Θ = σ/√ε₂` (figure, right). Slope 0.58 vs 0.5 is the
  same mild finite-ε drift the canard chapter documented (0.43–0.46, → 0.5 as
  ε → 0); 3 points + threshold choice, so "1/2 to within finite-ε corrections."
- **`C_q ≈ 5.7`** — O(several), squarely in the canard/folded-node family (2.8
  normal-form → 8–10 full-FHN). Threshold- and projection-dependent at the O(1)
  level (Cartesian-noise radial projection, `ρ`-collapse cutoff), exactly as
  `CANARD_BLOWUP.md` §6/§10 documents for `C_q`.
- **ω-independence** (the `ε₁` analog): structural here — the phase decouples in
  the symmetric model — so this confirms `σ_*` keys on the drift `ε₂`, not the
  cycle's phase rate, but does *not* yet test the *coupled* (phase-modulated) case.

## 4. What this does to the headline (the point of Part 2)

`FOLDED_CYCLE_NOISE.md` §10(ii) said: *if* `C_q ≈ 8–10` then `C_B/C_q ≈ 1.2` and
the A-vs-B crossover sits at `ε₂ ≈ 0.2`, so **"Channel B wins for full models
across the accessible range."** The worked model **does not support that**:

```
measured C_q ≈ 5.7   ⇒   C_B/C_q = 2π√3 / 5.7 ≈ 1.9
crossover at |ln ε₂| = (C_B/C_q)² ≈ 3.6   ⇒   ε₂ ≈ 0.03.
```

So with the measured prefactor, **Channel B wins only for `ε₂ ≲ 0.03`**, not across
the whole range — amplitude escape (A) sets `σ_*` for `ε₂ ≳ 0.03`. The worked model
**tempers the headline**: "B wins" is a *small-ε₂* statement, not a generic one.
(Caveat both ways: a full, stiffened folded-cycle model could lift `C_q` toward 8
as FHN's finite-ε corrections did, pushing the crossover back up; 5.7 is the clean
fold-of-cycles value, plausibly a lower estimate.)

## 5. Honest status and what's left

Validated (Path B, leading order, real model):
- `σ_*^A ∝ √ε₂` on a genuine **fold of limit cycles** (not a fold of equilibria) —
  the structurally-new object — with `C_q ≈ 5.7`.
- The VdP-Hopf failure is documented, not hidden; it motivates the SNLC model.

Open / not done here:
1. **Phase modulation `G(θ)` on the model — now PARTIALLY validated (§6, modest).**
   The symmetric Bautin alone validates only the inherited `√ε₂`/`C_q` half; §6
   adds an *anharmonic* Bautin and confirms the novel `√(ac/b)(θ)` modulation, but
   modestly (the collapse improves 1.32× when `G` is included; not crisp).
2. **Channel B on the model.** Phase diffusion / the `√(ε₂/|ln ε₂|)` law has not
   been measured on a real system — needs the anharmonic model + the iPRC pipeline.
3. **Tighten `C_q`.** 3 ε₂ points, `N=50`, one threshold; more ε₂ (smaller) and a
   threshold sweep would pin `C_q` and the exponent (→ 0.5).
4. **A genuine Liénard (JKK §5.1).** Still un-retrieved; the Bautin SNLC is the
   faithful stand-in. JKK §5.1 would be the canonical cross-check.

## 6. Phase modulation G(θ) — validating the *novel* content (modest)

The symmetric model grounds only the **inherited** (BG/canard-family) half. The
genuinely-new, not-BG/BGK content is the phase modulation
`σ_*^A(θ)=C_q√ε₂√(a(θ)c(θ)/b(θ))`. To test it on a real model we make the fold
*geometry* phase-dependent by modulating the **cubic** coefficient
(`folded_cycle_anharmonic.py`):

```
ρ' = ρ ( β + g(θ) ρ² − ρ⁴ ),   g(θ) = 1 + 0.5 cos 2πθ,    β' = −ε₂.
```

Closed-form local fold at frozen θ: SNLC at `β_SN=−g²/4`, `ρ_SN=√(g/2)`, and
`u' = −2g√(g/2) u² + √(g/2) ν` ⇒ JKK coefficients `a=√(g/2)`, `b=2g√(g/2)`,
`c=1` ⇒ **`G(θ)=√(ac/b)=1/√(2g(θ))`** — a clean, model-derived phase modulation
(not imposed). (Modulating `β` instead only shifts the fold *location*, not its
geometry — the cubic is the right knob.)

**Test (frozen-phase, canard R–Θ collapse).** For each phase θ (hence `g`), drift
β through the SNLC and measure the early-death margin `m`. Normalise by the
canard window `δν ∝ ε₂^{2/3} g^{−2/3}` and the noise scale `G`:

```
R = ⟨m⟩ g^{2/3} / ε₂^{2/3},     Θ = σ √(2g) / √ε₂  ( = σ / (√ε₂ · G) ).
```

If `√(ac/b)` is the right noise scale, the `R(Θ)` curves collapse across phase.

**Result (`figures/folded_cycle_anharmonic_Gtheta.png`):**

```
cross-phase collapse spread (mean std / scale), g ∈ [0.5, 1.5]:
   vs Θ = σ√(2g)/√ε₂   (G included):  CV = 0.19
   vs Θ = σ/√ε₂        (G ignored) :  CV = 0.25      (naive/correct = 1.32)
C_q (pooled R crosses canonical R=1):  ≈ 4.4   (consistent with symmetric 5.7)
```

**Honest read — modest but real.** Including the predicted `G(θ)=√(ac/b)` factor
**improves** the cross-phase collapse by ~32%, with the tightening concentrated in
the escape region (`R>1`) where it matters, and a `C_q` consistent with the
symmetric model. So the novel phase-modulation content **does** have real-model
backing — it is not normal-form-only anymore. But the confirmation is *soft*, not
crisp: the model's bistability constrains `g ≳ 0.5`, so the accessible `G`-range is
only a factor ~1.7; with finite-ε and Monte-Carlo scatter the collapse CV stays
~0.19 even with `G`. A wider `g`-modulation (larger `g₀`), smaller `ε₂`, and more
trajectories would sharpen it. This is also the prerequisite for Part 3: it is the
θ-dependent fold the α=1 phase must rotate *through*.

## 7. Channel B on the model — it does NOT blow up here (a real caveat)

`folded_cycle_channelB_model.py` measures the phase-diffusion coefficient
`D_φ = Var(φ)/T` of the Bautin cycle (degenerate noise on `u`) as the SNLC is
approached (`δ = β − β_SN → 0`). The leading-order prediction for a *smooth*
cycle is `D_φ = σ²/(2ρ_+²)` — **finite** at the fold.

```
 δ      ρ_+    |λ_⊥|   lead 1/2ρ²    D_φ/σ²  (σ=0.02 → 0.30)
 0.30  1.024   2.30      0.477       0.47 ........ 0.54
 0.15  0.942   1.37      0.564       0.57 ........ 1.83
 0.07  0.874   0.81      0.654       0.52 ........ 2.47
 0.03  0.820   0.47      0.743       0.72 ........ 5.59
```

**Result.** At small σ, `D_φ/σ²` sits on the finite leading value `1/(2ρ_+²)`
and does **not** diverge as `δ → 0` — **Channel B has no leading-order blow-up on
this smooth fold-of-cycles.** The only growth is at *larger σ and smaller δ*: the
**transverse (second-order, σ⁴) correction**, `~ σ⁴/|λ_⊥|`, blowing up as the
radial Floquet rate `|λ_⊥| → 0` (enhancement factor 1.15 → 7.81 as `δ: 0.30→0.03`).

**Why this matters (honest caveat for the paper).** Channel B's leading-order
`√(ε₂/|ln ε₂|)` blow-up (`FOLDED_CYCLE_NOISE.md` §9) came from a *fold-passage*
iPRC — the `∫dt/r²` of a relaxation-type cycle whose inner (canard) variable
`r₂ → 0`. A **smooth** saddle-node of cycles has a finite-radius cycle, so its
iPRC is bounded and the blow-up appears only as the higher-order transverse term.
So **the two channels are not symmetric in their model-generality**: Channel A
(amplitude escape) is robust and is grounded here; Channel B's leading-order
enhancement is **relaxation-cycle-specific** and is *not* exhibited by a smooth
fold of cycles. Grounding Channel B's leading-order law on a real model needs a
*relaxation-type* folding cycle (fold-passage iPRC) — which is in tension with a
clean finite-amplitude fold of cycles (forced VdP, the natural relaxation
candidate, failed Channel A via its Hopf, §1). This asymmetry is a genuine result
and a caveat the manuscript must state, not hide.

## 8. References

- `FOLDED_CYCLE_NOISE.md` (the law being grounded), `CANARD_BLOWUP.md` §6/§10
  (`C_q` protocol + finite-ε exponent drift), `VDP_CROSSMODEL.md` (canard `C_q≈8`).
- Jelbart–Kuehn–Kuntz 2024, arXiv:2208.01361 (deterministic folded cycle; §5.1
  Liénard). Kuznetsov, *Elements of Applied Bifurcation Theory* (Bautin / SNLC).
