# Direction C: the Airy process for successive / coupled peel-offs

*The headline / PhD-arc direction. A single peel-off is `TW_β` (the **edge**, a 1-point object);
the conjecture is that the **joint law of successive / coupled peel-offs** is a finer KPZ-class
object — the **Airy_β process**, whose 1-point marginal is `TW_β`. This note builds a principled
realization, tests the Airy-process signatures, and is explicit about what is shown vs. open.
Companion: `airy_process_peeloffs.py`, `figures/airy_process_peeloffs.png`. Tags **[R]/[N]/[H]/[open]**.*

---

## 0. The conjecture

The Airy₂ process `A(x)` (KPZ edge) is the stationary process with (i) one-point marginal
`A(x) =d TW₂`, (ii) local **Hölder-½** roughness `Var(A(x+r)−A(x)) ∼ 2|r|` (locally Brownian),
(iii) decorrelation `Cov(A(0),A(x)) → 0`. The `β`-generalisation `Airy_β` has `TW_β` marginals.
**Conjecture:** the peel-off profile of a *coupled* (spatially extended) folded cycle, or the
sequence of peel-offs of a *forced* one, converges (in the edge scaling) to `Airy_β`.

## 1. The construction [R framing]

A single peel-off is the stochastic Airy operator (SAO) ground state. Take a **spatial array**
of folded cycles, site `x`, each with inner noise **white in the sweep variable `Y` but
correlated across `x`** (the coupling): `ξ(Y,x)`, OU(`ℓ`) in `x`. This is precisely the **SAO
under stationary noise dynamics** — and the SAO is the edge of the `β`-ensemble, so evolving its
noise as Ornstein–Uhlenbeck is the edge of `β`-Dyson–OU, whose scaling limit *is* the (stationary)
`Airy_β` line ensemble. So this is the *right* generator to probe, not an arbitrary one. We
compute the peel-off profile `Y_node(x)` (recessive Cole–Hopf first node per site, sub-grid
interpolated) at `β=2` (`η=√2`), over many array realizations.

## 2. Results [N]

```
  (i)   1-point marginal :  mean −1.79 / std 0.92 / skew 0.26   vs  TW₂ (−1.77 / 0.90 / 0.22)
  (ii)  roughness        :  V(r)=Var(ΔY_node) ∼ r^{0.78}  locally  (robust: 0.73–0.78 over ℓ=80–200)
        saturation       :  V(∞) → 1.53–1.63 ≈ 2·Var(TW₂)
  (iii) covariance       :  C(0)=0.67 ≈ Var(TW₂)=0.81;  decays to ~0 over the coupling length ℓ
```

So the peel-off field **is** a stationary, rough, decorrelating random field with a **`TW_β`
one-point marginal** — the three structural hallmarks of an Airy-type (KPZ-edge) process. The
1-point marginal is `TW_β` to numerical accuracy **[N, solid]**; the field is genuinely
non-differentiable (it is *not* the excluded smooth `r²`) **[N]**; it decorrelates over the
coupling length **[N]**.

## 3. Resolution — the field converges to Airy₂ in the edge-scaling limit [N]

The `α<1` at `β=2` is **not** a finite-size effect (it is stable across `ℓ=80,200` and survives
sub-grid interpolation) — it is a **finite-noise** effect. Scanning `β=4/η²`
(`airy_eta_scaling.py`, `figures/airy_eta_scaling.png`):
```
   beta      2     4     8    16    22
   alpha   0.71  0.94  0.97  0.98  0.98     ->  1   (Airy_2, Hölder-1/2)
```
**`α → 1` as `η → 0` (`β → ∞`)** — exactly the edge-scaling limit where the universal `Airy₂`
lives. The deviation `1−α ∼ η^{p}`, `p ≈ 1.8`.

> **Result (updated).** The coupled peel-off field **converges to the Airy₂ process**: the
> `TW_β` one-point marginal holds at *all* `β`, and the Hölder-½ local-Brownian roughness
> (`α=1`) is recovered as `η→0`. At finite `β` it is an **`η`-roughened** Airy₂ — the nonlinear,
> beyond-linear-response part of `Y_node[noise]` adds roughness that vanishes with `η`. So
> Direction C's conjecture (successive/coupled canard peel-offs ↦ the Airy/KPZ edge process) is
> **supported on both axes** — marginal *and* roughness — in the proper limit.

## 4. What remains, to nail the exact law [open]

1. **Covariance — DONE (`airy_covariance.py`, `figures/airy_covariance.png`): the construction
   is a *local proxy*, not full Airy₂.** The eigenvalue covariance `C(r)` decays
   **exponentially** with length ≈ the coupling length `ℓ` (fit 66 vs `ℓ=80`); it is flat then
   drops on log–log — *not* the Airy₂ **algebraic `r^{-2}`** tail. Linear response explains it:
   `Λ₀` is a *local* functional of the noise column, so OU-in-`x` ⇒ exponential-in-`x`. A Dyson-OU
   GUE edge (the canonical Airy₂) shows the heavier algebraic tail by contrast. **So correlated
   noise reproduces Airy₂'s *marginal* (TW_β) and *local* roughness (Hölder-½ as `η→0`), but the
   genuine long-range Airy₂ needs Dyson-type *spectral dynamics* (eigenvalue interaction) — which
   a real coupled array must supply through the coupling, not via shared noise alone.** This makes
   the FHN-array test (item below) the decisive one: does genuine diffusive coupling induce the
   eigenvalue interaction (→ algebraic `r^{-2}`), or only correlated noise (→ exponential)?
2. **The finite-`η` correction.** Pin `1−α ∼ η^{p}` — is `p=2`, or the recurring `8/3`? — via a
   small-`η` expansion of the SAO ground state's nonlinear (beyond-linear-response) part.
3. **Node vs eigenvalue — settled (with a numerical caution); Conjecture 2 stands.** Two
   eigenvalue methods disagree and the *exact* one wins. **Inverse iteration** gave
   `α(Λ₀)≈0.57` — which *appeared* to overturn Conjecture 2 — but this is a **numerical
   artifact** (grid-sensitive / incomplete convergence of the iterated eigenvector). The **exact
   Sturm-bisection** eigenvalue (`gap_roughness.py`) gives `α(Λ₀)≈0.9` (grid-noisy `0.87–0.95`,
   marginal still converging), **clearly smoother than the node** (`α≈0.79`) and consistent with
   the Airy₂ value `1` within scatter — as perturbation theory demands (`Λ₀` is an analytic,
   polynomial-in-noise functional ⇒ Hölder-½). So the eigenvalue field is (approximately) the
   smooth Airy₂; the finite-`η` roughening lives in the **node / level-crossing**
   `Y_node=G^{-1}(0)`. *(Earlier text here claimed the opposite, on the strength of the
   inverse-iteration number — retracted.)*
4. **Gap / avoided-crossing mechanism (Conjecture 3) — REFUTED** (`gap_roughness.py`,
   `figures/gap_roughness.png`). The roughening is **not** driven by small `Λ₁−Λ₀` gaps: small
   gaps are essentially absent (strong level repulsion, `frac(gap<0.3)=0.002`), and the local
   roughness *decreases* as the gap pinches (`E[(ΔΛ₀)²|gap]∼gap^{+0.26}`, the wrong sign). The
   source is the level-crossing/inversion `G^{-1}(0)` itself, not the spectral gap.
5. **Real coupled array — DONE (`coupled_array_peeloff.py`, `figures/coupled_array_peeloff.png`):
   diffusive coupling does NOT give Airy₂ either.** A genuinely coupled array (independent per-site
   noise, linked only by diffusive coupling `D·Δ_i u`) gives a peel-off covariance that is
   **short-range diffusive** — correlation length `∼√(D·T_sweep)` (≈7 at `D=4`, `C(20)≈0`),
   exponential, *not* algebraic `r^{-2}`; strong coupling also distorts the marginal
   (std `0.90→0.49`). So nearest-neighbour coupling is short-range, **not** the all-to-all Dyson
   interaction Airy₂ needs.

## 5. Reproduce

`python3 airy_process_peeloffs.py` → `figures/airy_process_peeloffs.png` (peel-off profiles; the
`TW₂` 1-point marginal; the increment variance `V(r)∼r^{0.78}` against `Airy₂`'s `r¹`, the smooth
`r²`, and the `2·Var` saturation). Recessive Cole–Hopf first node, vectorised over the array,
noise OU-correlated across sites. `β` set by `η=2/√β`.

---

## Final status (Direction C)

A single canard peel-off **is** TW_β (the SAO ground state), and the *successive modes of one
folded cycle's inner operator* are the **Airy₂ point process** — the genuine Airy/KPZ connection,
at the level of a **single** spectral edge. **Confirmed quantitatively** (`airy_line_ensemble.py`,
`figures/airy_line_ensemble.png`): the first three peel-offs `−Λ₀,−Λ₁,−Λ₂` of one folded cycle
match the GUE soft-edge point process `ξ_k=N^{2/3}(λ_k−2)` and the asymptotic Airy₂ means
`(−1.771, −3.675, −5.172)` to ~1%, with the `√|·|` shrinking-spacing law (`1.90→1.52`) and the
`gap²` (`β=2`) small-gap level repulsion. But the **Airy₂ *process*** (built by varying an
external parameter) is **not** realized by the natural couplings tested: correlated noise
(→ exponential covariance, a local proxy) and diffusive coupling (→ short-range diffusive,
distorted marginal) both reproduce the **marginal** (TW_β) and the **local roughness**
(Hölder-½ as `η→0`) but **not** the algebraic long-range `r^{-2}` Airy₂ correlations. Those
require the **all-to-all Dyson eigenvalue interaction**, which local couplings of separate
oscillators do not supply.

**Verdict.** The canard peel-off Airy/KPZ structure is a **single-operator (spectral)**
phenomenon — TW_β marginal, the Airy line ensemble of successive modes — **not** an Airy₂
*process* emergent from coupling separate folded cycles. Tested on both natural mechanisms, the
process-level conjecture comes out **negative**; the defensible positive statement is the
single-edge one (TW_β + the line ensemble). A genuine Airy₂ process would need a model whose
coupling acts as true spectral repulsion — one large folded system with many interacting modes,
not an array of coupled cycles.
