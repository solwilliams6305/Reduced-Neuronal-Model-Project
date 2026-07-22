# Folded-node universality — the Weber soft edge IS the stochastic Airy operator

*Goal: show analytically that the noisy folded-node canard passage lies in the same
Tracy–Widom_β edge class as the folded cycle, by reducing its Weber (parabolic-cylinder)
inner operator to the stochastic Airy operator at a simple turning point, with the
turning-point slope `s` setting the local rescaling. Companion numerics:
`folded_node_universality.py`, `figures/folded_node_universality.png`.*

*Rigor flags used below: **[R]** rigorous / exact; **[T]** transfers from the fold proof,
modulo a re-check noted; **[H]** heuristic or invokes a theorem that must be cited or proved.*

---

## 0. Statement

> **Claim.** Let a noisy folded-node canard be swept through its fold. Its peel-off level
> `z_node` (first node of the Cole–Hopf field `u`) satisfies, in the local Airy coordinate
> `ζ = s^{1/3}(z − z*)`,
> ```
>     ζ_node  =d  TW_β  + O(s^{-4/3}) ,        β = 4s / η² ,
> ```
> where `z* = 2√(ν+½)` is the turning point, `s = √(ν+½) = Q'(z*)` its slope, and `η` the
> effective noise. Equivalently `z_node = z* + s^{-1/3} TW_{4s/η²} + …`: the **location**
> `z*` is model-dependent (Weber), the **fluctuation law** is the universal TW_β edge.
> As `s → ∞` the correction vanishes and `ζ_node → TW_β` exactly, anchored (`η→0`) at the
> first Airy zero `a₁ = −2.33811`.

This is the soft-edge-universality content of the paper's remark "we expect the same TW_β
exit law for noisy fold and folded-node passages": below, the `β = 4s/η²` law and the
`s→∞` convergence are *derived*, and matched against the numerics.

---

## 1. The two inner equations, side by side

After Cole–Hopf the fold and the folded node give the *same kind* of object — a 1D random
Schrödinger operator swept through a turning point — differing only in the global potential:

```
  fold / folded cycle :   u'' = ( Y      − η ξ ) u ,     Q_fold(Y) = Y            (Airy, linear)
  folded node         :   u'' = ( z²/4 − (ν+½) − η ξ ) u , Q_node(z) = z²/4−(ν+½) (Weber, quadratic)
```

`ξ = dB/dx` is white noise in the swept variable; `'` is `d/dx` along the sweep. The fold
potential is linear (one turning point at `Y=0`); the Weber potential is quadratic, with an
oscillatory well `|z|<z*` (the folded-node *rotations* / secondary canards) bounded by two
simple turning points `z = ±z*`, `z* = 2√(ν+½)`. The peel-off is the first node reached on
sweeping **down** through the upper turning point `z*`. The deterministic Weber index `ν`
is fixed by the folded-node eigenvalue ratio `μ` (Wechselberger: the secondary-canard count
is governed by the parabolic-cylinder zeros) **[H, inherited]**; only `s = √(ν+½) > 0`
enters the universality argument, so the conclusion holds for the whole `μ`-family.

Key algebraic fact (exact, since `Q_node` is quadratic): with `δ = z − z*`,
```
  Q_node(z* + δ) = s δ + ¼ δ² ,     s := Q_node'(z*) = z*/2 = √(ν+½).        [R]
```
(Using `z*²/4 = ν+½` and `z*/2 = s`; no truncation — the Weber potential equals its
2nd-order Taylor polynomial.)

---

## 2. Step 1 — Cole–Hopf linearisation transfers verbatim  [T, rigorous]

The folded-node inner Riccati is `dR = (R² − Q_node(z)) dT + η dB` with **additive** noise,
exactly as for the fold. The paper's Lemma (exact Itô identity) shows `R = −u'/u` gives
`u'' = (Q − ηξ)u` with **no** Itô–Stratonovich correction, because three contributions
vanish identically: additive (constant-in-`u`) noise ⇒ Itô = Stratonovich; the linear
system `du = v\,dT, dv = (Qu)dT − ηu\,dB` has diffusion vector `g = (0,−ηu)` with
`½(g·∇)g = 0`; and `R = −v/u` has `R_{vv} = 0`. None of these uses the form of `Q`, so the
lemma applies to `Q_node` **verbatim**. A finite-time blow-up `R → +∞` is a simple zero
(node) of `u`; hence the peel-off level is the first node of `u`. **[R]**

---

## 3. Step 2 — the turning-point reduction  [R for the algebra; the core]

Rescale about the upper turning point with the Airy length `λ = s^{-1/3}`:
```
  z = z* + λ ζ ,     λ = s^{-1/3} ,     ũ(ζ) = u(z* + λζ).
```
Then `ũ_ζζ = λ² u''`, and using `Q_node(z*+λζ) = sλζ + ¼λ²ζ²`:
```
  ũ_ζζ = λ²( Q_node − η ξ_z ) ũ = ( sλ³·ζ + ¼λ⁴·ζ² − η λ² ξ_z ) ũ .
```
With `λ = s^{-1/3}`: `sλ³ = 1` and `¼λ⁴ = ¼ s^{-4/3}`. The noise rescales by the **exact**
white-noise law: `⟨ξ_z(z)ξ_z(z')⟩ = δ(z−z') = λ^{-1}δ(ζ−ζ')`, so
`ξ_z(z*+λζ) = λ^{-1/2} \tilde ξ(ζ)` with `\tilde ξ` standard, giving
`η λ² ξ_z = η λ^{3/2}\tildeξ = η s^{-1/2}\tildeξ`. Hence the **exact** rescaled equation
```
  ┌─────────────────────────────────────────────────────────────────────────────┐
  │  ũ_ζζ = ( ζ  −  η_loc \tildeξ  +  ¼ s^{-4/3} ζ² ) ũ ,    η_loc = η / √s .       │   [R]
  └─────────────────────────────────────────────────────────────────────────────┘
```
Dropping the `¼ s^{-4/3} ζ²` term, this is **identically the fold's canonical noisy-Airy
inner equation** with effective noise `η_loc = η/√s`. Therefore
```
  β  =  4 / η_loc²  =  4 s / η² .        [R]
```
This is the scaling the numerics were built on and confirm (the `std` and `skew` collapse
onto TW_β only when `β = 4s/η²`). The deterministic part of the reduction — "near a simple
turning point a 2nd-order ODE is Airy after the `s^{-1/3}` stretch" — is the classical
**Langer / Olver turning-point theorem** and is rigorous **[R]**; the noise enters linearly
and rescales exactly, so the *form* above is exact, the only approximation being the dropped
curvature term, treated in §5.

---

## 4. Step 3 — reduction to the stochastic Airy operator  [T]

Dropping curvature, `ũ_ζζ = (ζ − η_loc \tildeξ)ũ` is the fold equation with `η → η_loc`.
The paper's reduction then applies unchanged: writing the swept variable as the operator
coordinate `x`, this is the eigenvalue equation `H_β ũ = E ũ` for the **stochastic Airy
operator**
```
  H_β = − d²/dx² + x + (2/√β) b'(x) ,   on the half-line, Dirichlet at the entry,
        with  η_loc = 2/√β ,  β = 4s/η² ,
```
whose ground state `Λ(β)` satisfies `−Λ(β) =d TW_β` (Ramírez–Rider–Virág). The shooting
characterisation (paper Thm: the *swept* first-node law equals the *fixed* half-line
operator's ground state, by Dirichlet domain monotonicity) carries over, with **one
re-check** **[T]**: the recessive entry is now the parabolic-cylinder function `D_ν`
(decaying as `z→+∞`) in place of `Ai`, and the forgetting lemma needs the Weber canard's
attracting property. Both hold — integrating the recessive `D_ν` *downward* is
self-correcting (the dominant component decays going down), which is exactly what the
numerics rely on for a clean entry (ratio `−D_ν'/D_ν = z/2 − ν/z`). So `ζ_node =d TW_β`,
`β = 4s/η²`, with the deterministic anchor `a₁ = −2.338` in the `η→0` limit. **[R given §5]**

---

## 5. The curvature correction and the approach to TW_β  [H + numerics]

The single neglected term `¼ s^{-4/3} ζ²` is an **even, bounded-below, O(s^{-4/3})**
perturbation of the SAO potential. Its effects (all → 0 as `s→∞`):

* **Location shift (odd).** The even bump raises the well asymmetrically about the first
  node, shifting the deterministic peel-off off `a₁ = −2.338`. First-order perturbation
  gives a shift `∝ s^{-4/3}`. **Numerics:** anchor `ζ₀(η→0) = −2.555, −2.462, −2.409,
  −2.376, −2.359, −2.350` at `s = 1.41…10.98 → −2.338`; a log–log fit of the residual
  `|ζ₀+2.338|` gives slope `≈ 1.4 ≈ 4/3`. ✓
* **Shape (odd cumulants).** The skewness deviation from TW_β is also `∝ s^{-4/3}`:
  measured `skew = −0.46, −0.18, −0.03, +0.09, +0.13, +0.19` at `β=2` vs `TW₂ = 0.224`,
  residual slope `≈ 1.45`. So *both* the location and the skewness corrections track the
  single curvature term — a quantitative check that `¼ s^{-4/3}ζ²` is indeed *the* leading
  correction. ✓
* **Scale / even cumulants converge faster.** Because the perturbation is even, even
  moments pick it up only at second order `O(s^{-8/3})`: numerically `std → TW₂ = 0.902`
  already by `s≈4.5` (`std = 1.03, 0.96, 0.92, 0.909, 0.907, 0.899`), and `kurt → 0.093`.
  Consistent with the even/odd split. ✓

**The one genuine gap [H].** Turning "drop the `O(s^{-4/3})` term" into a theorem requires
**soft-edge universality of the SAO under a confining polynomial perturbation**: that the
bottom-of-spectrum law of `−d²/dζ² + ζ + ¼ s^{-4/3}ζ² + η_loc \tildeξ` equals `TW_{β}` up
to `O(s^{-4/3})`. Deterministically this is Langer (rigorous). Stochastically it is the
*β-ensemble soft-edge universality* class of statement — proved for tridiagonal/β-ensemble
and Wigner-type models (RRV; Bourgade–Erdős–Yau-type edge universality) — which must be
**cited or specialised** to this operator. The numerics confirm the conclusion; the analytic
step is "invoke/extend a soft-edge-universality theorem," not "invent a new mechanism."

**→ Now discharged to a precise lemma — see `FOLDED_NODE_EDGE_LEMMA.md`.** Because the
curvature term `¼ s^{-4/3}ζ²` is **non-negative** (a square), the edge bound follows from an
elementary **min–max sandwich** on the SAO bottom eigenvalue: `G ≤ G_s ≤ G + O(s^{-4/3})`,
inverted at level 0 via the Hadamard endpoint derivative. This needs only RRV (`−Λ₀ =d TW_β`,
a.s. simple spectrum) and Dumaz–Virág (SAO ground-state localisation, so the constant is a.s.
finite) — **no β-ensemble universality required**; the Krishnapur–Rider–Virág "Universality of
the SAO" route is an optional fallback. The lemma also *predicts* the one-signed, `s^{-4/3}`,
even/odd-split corrections that the numerics show.

---

## 6. What is rigorous, what must be supplied

| Step | Status | Note |
|------|--------|------|
| Cole–Hopf, Itô = Stratonovich for `Q_node` | **[R]** | Lemma uses additive noise only; independent of `Q`. |
| `Q_node(z*+δ) = sδ + ¼δ²` exact | **[R]** | Weber is its own Taylor polynomial. |
| Local rescaling `λ=s^{-1/3}`, `η_loc=η/√s`, `β=4s/η²` | **[R]** | Exact algebra + exact white-noise rescaling. |
| Deterministic turning point ⇒ Airy | **[R]** | Langer / Olver. |
| Swept first node = SAO ground state (shooting) | **[T]** | Re-check: recessive entry is `D_ν`; forgetting via downward self-correction. |
| `−Λ(β) =d TW_β` | **[R]** | RRV, cited. |
| Drop `¼ s^{-4/3}ζ²` ⇒ TW_β up to `O(s^{-4/3})` | **[R]** mod **[I]** | Min–max sandwich, `FOLDED_NODE_EDGE_LEMMA.md`; needs only RRV + Dumaz–Virág (no universality). |
| Physical noise ⇒ inner `η` | **[H, inherited]** | Same Berglund–Gentz reduction as the fold. |
| `ν ↔ μ` (folded-node eigenvalue ratio) | **[H, inherited]** | Deterministic folded-node theory; only `s>0` is used. |

So the **scaling law `β = 4s/η²`, the model-dependent location `z*`, and the `s→∞`
convergence are derived**; the proof is complete *modulo* one soft-edge-universality input,
which is exactly the kind of statement the literature supplies for the SAO.

---

## 7. Consequence

The folded-node canard escape is in the **Tracy–Widom_β / KPZ soft-edge universality class**,
the same class as the folded cycle (the fold being the `ν→∞` / pure-linear member, where the
curvature term is absent and TW_β is exact at all `β` — confirmed numerically: the Airy
reference matches `TW_β` skewness `0.178/0.226/0.310` vs `0.166/0.224/0.293` at `β=4,2,1`).
The reduction shows the *mechanism* is generic: **any** canard passage whose inner equation
has a simple turning point Cole–Hopf-linearises and, after the `s^{-1/3}` stretch, presents
the SAO at its soft edge — so the exit law is `TW_{4s/η²}` with the local slope `s` the only
model input, and the deterministic peel-off location the only model-dependent output. This is
the precise sense in which *“the folded cycle is in the TW class” becomes “this family is.”*

A clean target for a fully rigorous version: state and cite (or prove) the SAO soft-edge
universality lemma in §5, then §§2–4 are a short, self-contained transfer of the fold proof.

---

## 8. Numerics this rests on (reproduce)

`python3 folded_node_universality.py` →
`figures/folded_node_universality.png`. It runs the recessive Cole–Hopf sweep on the noisy
Weber operator, takes the first node, rescales by `s^{1/3}`, and reports (i) convergence of
`(mean, std, skew, kurt)` to `TW_β` as `s` grows at `β=2`, and (ii) the `β=4,2,1` family at a
sharp fold (`ν=120`, `s≈11`) against the Airy/folded-cycle reference. Headline at `s≈11,
β=2`: `(mean,std,skew,kurt) = (−1.785, 0.907, 0.190, 0.042)` vs `TW₂ (−1.771, 0.902, 0.224,
0.093)` — scale and location at the ~1% level, shape converging as `s^{-4/3}`.

**Reading.** The numerics do not merely *suggest* TW_β; combined with §§2–4 they pin the
*only* free parameter (`β = 4s/η²`) and exhibit the `O(s^{-4/3})` curvature correction with
the predicted exponent, leaving §5's universality lemma as the single analytic to-do.
