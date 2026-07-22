# The phase edge (SNIC) and the two-edge statement

*Direction B. Goal: the universal first-passage law of the noisy saddle-node-on-circle (SNIC)
phase channel — the phase twin of the amplitude (Tracy–Widom) edge — its constants, whether it
has a TW-style closed form, and its slow tails. Companion: `phase_edge_snic.py`,
`figures/phase_edge_snic.png`. Tags: **[R]** proved/exact here, **[T]** transfers from the
amplitude proof, **[N]** numerical, **[H]** heuristic/standard.*

---

## 0. Result in one line

The folded limit cycle has **one** inner Riccati `dR = (R²−Y)dT + η dB`; the destroying
bifurcation chooses **which question** you ask of it, and that choice — not a different
operator — is what makes one edge integrable and the other not:

| | amplitude edge (fold of cycles) | phase edge (SNIC) |
|---|---|---|
| read the Riccati by | **sweeping** `Y = Y₀−T` | **holding** `Y` at criticality |
| observable | peel-off **level** `Y_node` | rotation **period** = first-passage **time** |
| inner operator | Airy (linear potential) | quartic (anharmonic) |
| law | `Y_node =d TW_β`, `β=4/η²` | noisy-saddle-node FPT |
| closed form? | **yes** (Painlevé II / `TW_β`) | **no** for the law; **yes** for the mean |
| scale | `σ_* ∼ √ε₂` | `ω, D_φ ∼ σ^{2/3}` |

"Where" the canard peels off is integrable (Tracy–Widom); "when" it rotates is not.

---

## 1. The phase channel is the fixed-Y noisy saddle-node  [T]

At a SNIC the cycle dies at **infinite period**: the phase develops a saddle-node bottleneck.
The inner phase equation is the noisy saddle-node normal form `dx = (μ + x²)dτ + σ dB`, i.e.
the **same Riccati** `dR = (R² − Y)dT + η dB` with `Y = −μ` held **fixed** (`Y = 0` at
criticality), rather than swept. Cole–Hopf `R = −u'/u` linearises it to `u'' = (Y − ηξ)u` with
**no Itô–Stratonovich correction** — the amplitude-channel lemma uses only additive noise, so
it transfers verbatim **[T]**. A blow-up `R → +∞` is a node of `u`; one rotation of the phase =
one node = the phase advancing through the bottleneck. So:

> the **rotation period** is the **first-passage time** of the noisy saddle-node, equivalently
> the **node spacing** of `u'' = (Y−ηξ)u` at fixed energy (the rotation number / Frisch–Lloyd
> current `J_FL(Y)` is its inverse rate; `J_FL(0)` is exactly the `y=0` value computed in
> `frisch_lloyd_current.py`).

This is the structural twin of the amplitude edge, where the *swept* first node gave the SAO
ground state. Same operator, complementary spectral question (level vs spacing).

## 2. `σ^{2/3}` scaling of `ω, D_φ`  [R]

At criticality `dR = R²dT + σ dB`. The rescaling `R = σ^{2/3}ρ`, `T = σ^{-2/3}s` is exact and
gives the `σ`-free canonical equation
```
    dρ = ρ² ds + dW̃ ,
```
so the rotation rate and phase diffusion scale as `ω, D_φ = c · σ^{2/3}`. **[N]** the physical
equation collapses under this rescaling and the fitted exponent is `0.666` (predicted `2/3`).

## 3. The rate constant `J` — it HAS a closed form  [R]

The mean rotation period is `J = ∬_{w<u} e^{(w³−u³)/3}` (paper normalisation). The naive double
integral converges only as `1/U`: its integrand `g(u)=∫_{-∞}^u e^{(w³−u³)/3}dw ∼ 1/u²` has
**algebraically slow tails** (this is the paper's "slowly-converging tails" — confirmed
numerically, `u²g(u) → 1.00`). But doing the **inner `u`-integral exactly** (it is Gaussian)
collapses `J` to a one-dimensional integral and then a Gamma value:
```
   J = ∬_{w<u} e^{(w³−u³)/3}
     = ∫_0^∞ dv ∫_{-∞}^∞ du  e^{-u²v + u v² - v³/3}        (v = u − w)
     = √π ∫_0^∞ v^{-1/2} e^{-v³/12} dv                       (Gaussian in u)
     = (√π / 3) · 12^{1/6} · Γ(1/6)  ≈  4.9761 .             [D = 1, paper]
```
In the dB-normalisation (the physical `σ dB`, which rescales to generator `½∂² + ρ²∂`) the same
reduction gives the mean period
```
   J_dB = (√(2π)/3) · 6^{1/6} · Γ(1/6)  ≈  6.2694
        = 2^{1/3} · J ,   matching  1/J_FL(0) = 6.33  to MC accuracy.
```
**[N]** the accurate double-integral quadrature reproduces `4.9761` to 4 dp, and the simulated
mean (`6.077` with `±12` cutoffs, `+0.17` for the `1/u²` tail) gives `6.27 = J_dB`.

> **Refinement of the paper.** "The phase-channel constant `J` has slowly-converging tails" is
> a statement about the *naive quadrature*; the constant itself is the **exact Gamma value**
> above. It is the *distribution*, not `J`, that lacks a closed form (§4).

## 4. The first-passage LAW: a quartic, hence no Tracy–Widom closed form  [R + N]

The FPT generator (dB-normalisation, `D=½`, drift `b=ρ²`) is `L = ½∂² + ρ²∂`. The Liouville
substitution `φ = e^{-ρ³/3} χ` removes the drift and turns the eigenproblem `−Lφ = λφ` into a
**Schrödinger equation with a quartic potential** **[R]**:
```
   − χ'' + ( ρ⁴ + 2ρ ) χ = 2λ χ .
```
The `ρ⁴` comes from the drift-squared `b²`, the `2ρ` tilt from the saddle-node curvature `b'`.
This is the **anharmonic (quartic) oscillator** — the textbook example of a spectrum with **no
closed form** — in sharp contrast to the amplitude edge, whose swept operator `u''=(Y−ηξ)u` is
the **Airy** (linear-potential) operator with the integrable Painlevé II / `TW_β` law. Hence:

> **the phase-edge first-passage law has no Tracy–Widom / Painlevé closed form**; only its mean
> (`J`, §3) is closed.

**Consequences, all [N]-confirmed:**
- **Shape (not TW).** The standardised law has `skew ≈ 1.92`, `excess kurt ≈ 5.69` — far from
  `TW₂ (0.22, 0.09)` / `TW₁ (0.29, 0.17)`. It is strongly right-skewed (a first-passage time),
  not a soft-edge eigenvalue law.
- **Right tail = quartic ground state.** `P(T) ∼ e^{-λ₀ T}` with `λ₀ = E₀/2`, `E₀` the ground
  state of `−d²/dρ² + ρ⁴ + 2ρ`. Numerically `E₀ = 0.562 ⇒ λ₀ = 0.281`, versus the simulated tail
  rate `0.278` — a <2% tie, anchoring the law's tail to the quartic spectrum. (The absorbing-BC
  correction is within this gap; see remark.)
- **Left tail.** Super-exponentially small (a minimum ballistic crossing time), as for any
  bottleneck FPT.

## 5. Numerics summary (`phase_edge_snic.py`)

```
  J (closed form, D=1)         = 4.9761   = numeric double-integral 4.9761   (g(u) ~ 1.00/u²)
  J_dB (closed form, D=1/2)    = 6.2694   = 1/J_FL(0) (6.33) to MC accuracy
  FPT mean (sim, ±12 cutoffs)  = 6.08  (+0.17 tail -> 6.27)
  FPT std / skew / exkurt      = 3.60 / 1.92 / 5.69        => NOT Tracy-Widom
  right-tail rate  lambda0     = 0.278 (sim)  vs  E0/2 = 0.281 (quartic ground state)
  sigma^{2/3} exponent (fit)   = 0.666        (predicted 2/3),  rescaled law collapses
```

## 6. Proved vs heuristic

| Step | Status |
|------|--------|
| Phase channel = fixed-Y noisy-saddle-node Riccati; Cole–Hopf (Itô=Strato) | **[T]** (transfers; additive noise) |
| `σ^{2/3}` rescaling to `dρ=ρ²ds+dW̃` | **[R]** exact |
| `J = (√π/3)12^{1/6}Γ(1/6)` (Gaussian + Gamma reduction) | **[R]** exact |
| `g(u) ∼ 1/u²` slow tails | **[R]** (asymptotics) + **[N]** |
| FPT generator ↦ quartic `−χ''+(ρ⁴+2ρ)χ=2λχ` | **[R]** Liouville transform |
| quartic anharmonic ⇒ no closed-form spectrum / law | **[H]** standard (anharmonic oscillator) |
| right-tail rate `λ₀ = E₀/2` | **[N]** (tie to ~2%; absorbing-BC caveat) |
| "not TW" (shape) | **[N]** (skew/kurt) + **[R]** (different operator class) |

## 7. The two-edge statement (combined result)

> **Two-edge proposition.** The noisy folded limit cycle's inner Riccati
> `dR = (R²−Y)dT + η dB`, Cole–Hopf `u'' = (Y−ηξ)u`, carries two complementary universal edge
> laws, selected by the destroying bifurcation:
>
> 1. **Amplitude edge (fold of cycles).** Sweeping `Y=Y₀−T`, the peel-off level is the
>    stochastic-Airy ground state: `Y_node =d TW_β`, `β=4/η²`, critical noise `σ_*∼√ε₂` — an
>    **integrable** edge (Painlevé II closed form).
>
> 2. **Phase edge (SNIC).** Holding `Y` at criticality, the rotation period is the noisy
>    saddle-node first-passage time: `ω, D_φ ∼ σ^{2/3}` with closed mean rate constant
>    `J = (√π/3)12^{1/6}Γ(1/6)`, but the **law** is governed by the quartic anharmonic operator
>    `−χ''+(ρ⁴+2ρ)χ=2λχ` — **non-integrable**, no Tracy–Widom/Painlevé closed form, exponential
>    tail `λ₀=E₀/2`.
>
> The two edges are the **"where"** and the **"when"** of the same random Schrödinger operator:
> sweeping the spectral parameter (extreme eigenvalue → `TW_β`) versus holding it and timing the
> node spacings (rotation number / first passage → the quartic SNIC law). The bifurcation that
> destroys the cycle picks which channel is singular — fold of cycles → amplitude (Airy/TW)
> edge; SNIC → phase (quartic) edge.

This is the self-contained two-edge result the paper's Channel-A / Channel-B split was pointing
to: an Airy edge and a quartic edge from one operator, one integrable and one not.

## 8. Off-critical: the SNIC edge is the centre of a one-parameter crossover

Off criticality (`Y = −μ ≠ 0`) the canonical family gains the detuning
```
   dρ = (ν + ρ²) ds + dW̃ ,        ν = μ / σ^{4/3}      (crossover at μ ∼ σ^{4/3}),
```
and the FPT operator becomes the **tilted quartic** `−χ'' + (ρ⁴ + 2ν ρ² + 2ρ + ν²)χ = 2λχ`
**[R]**. The single parameter `ν` interpolates three regimes (`offcritical_phase_edge.py`,
`figures/offcritical_phase_edge.png`):

| `ν` | operator | dynamics | mean period `J(ν)` | FPT shape |
|---|---|---|---|---|
| `ν<0` (pre-SNIC) | **double well** `±√|ν|` | Kramers activated escape | `∼ exp((8/3)|ν|^{3/2})` | broad, exponential (`CV→1`) |
| `ν=0` (SNIC) | pure quartic | critical edge | closed-form `J` (§3) | quartic law (`CV≈0.6`, skew≈1.9) |
| `ν>0` (post-SNIC) | **single well** | deterministic rotation | `→ π/√ν` | narrow, peaked (`CV→0`) |

**[N]** `J(ν)` by MFPT quadrature, confirmed by simulation (`ν=−1`: 51.5 vs 52.6; `0`: 5.96 vs
6.29; `2`: 1.92 vs 2.24); tail rate `λ₀=E₀/2` from the quartic ground state runs `≈0` (Kramers
slow tail, `ν=−2.5`) → `0.281` (`ν=0`) → `5.7` (`ν=3`); the FPT shape morphs `CV = 0.92, 0.61,
0.23, 0.15` across `ν = −1, 0, 2, 4`.

**Loop-closer (ties to the σ_crit / Arrhenius work).** In physical variables the pre-SNIC
(`ν<0`) period is
```
   J ∼ exp( (8/3) |ν|^{3/2} ) = exp( (8/3) μ^{3/2} / σ² )   ⇒   σ² ∼ (8/3) μ^{3/2} / log(period),
```
an **Arrhenius critical-noise law** with an explicit `μ^{3/2}` barrier — the *same* `Y^{3/2}`
barrier physics as the amplitude channel and the FHN `σ_crit(I)` table from the start of the
project. So the phase channel has its own nesting: the SNIC `σ^{2/3}` edge (the dynamic floor)
sits between the Arrhenius/Kramers regime (`ν<0`, the `σ_crit` ceiling) and the regular `σ²`
phase diffusion (`ν>0`, the fold-of-cycles / finite-period limit the paper calls "regular").
The SNIC edge is the critical point of that one-parameter family. **[R]** family + quartic;
**[R/H]** the `exp((8/3)|ν|^{3/2})` and `π/√ν` asymptotes; **[N]** `J(ν)`, `λ₀(ν)`, the shape.

## 9. Reproduce

`python3 phase_edge_snic.py` → `figures/phase_edge_snic.png` (FPT law + exponential/quartic tail;
`σ^{2/3}` collapse; the `1/u²` integrand vs the closed-form `J`).
`python3 offcritical_phase_edge.py` → `figures/offcritical_phase_edge.png` (the `ν` crossover
map, the tilted-quartic potential family, the FPT-shape morph). Quartic levels by
finite-difference; canonical / physical / off-critical FPTs by vectorised Euler with reinjection;
`J(ν)` by overflow-safe MFPT quadrature.
