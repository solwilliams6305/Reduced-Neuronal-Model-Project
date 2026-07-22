# The perturbed-SAO soft-edge lemma (folded-node universality, §5 of `FOLDED_NODE_UNIVERSALITY.md`)

*Purpose: discharge the one analytic gap in the folded-node universality argument — dropping
the `O(s^{-4/3})` Weber-curvature term from the rescaled inner operator. We state it as a
formal lemma, give a proof sketch that is rigorous modulo two **citable** facts about the
stochastic Airy operator, and pin the references to chase. Notation follows
`FOLDED_NODE_UNIVERSALITY.md`; `[R]`=proved here / standard, `[I]`=invoked (cite).*

---

## 1. Setup

Fix `β > 0` and let `ξ̃ = b'` be white noise on `R` (`b` a two-sided Brownian motion).
Write `η = 2/√β`. Following Ramírez–Rider–Virág and the paper's shooting characterisation,
introduce the **half-line, energy-zero** form of the stochastic Airy operator: for `ℓ ∈ R`
let
```
   M  = − d²/dζ² + ζ − η ξ̃ ,         G(ℓ) := bottom Dirichlet eigenvalue of  M  on [ℓ, ∞).
```
`G` is a.s. continuous and **strictly increasing** (Dirichlet domain monotonicity), and the
**peel-off** is its zero,
```
   ζ_node := the unique ℓ with G(ℓ) = 0 ,        ζ_node  =d  −Λ₀(β)  =d  TW_β ,
```
the last identity being RRV (bottom of the SAO spectrum) `[I]`. The deterministic (`η→0`)
value is the first Airy zero `a₁ = −2.33811`.

The folded-node reduction (`FOLDED_NODE_UNIVERSALITY.md` §3) produces the **same** operator
plus a deterministic potential:
```
   M_s = M + W_s ,        W_s(ζ) = ¼ s^{-4/3} ζ²  ≥ 0 ,
   G_s(ℓ) := bottom Dirichlet eigenvalue of M_s on [ℓ, ∞),     ζ_node^{(s)} := root of G_s = 0,
```
with `s = √(ν+½) = Q'(z*)` the turning-point slope. We must show `ζ_node^{(s)} ⇒ TW_β` as
`s → ∞`, with a controlled rate.

---

## 2. The lemma

> **Lemma (perturbed-SAO soft edge).** Let `W_s : R → [0,∞)` be **deterministic** with
> ```
>    (P1)  W_s ≥ 0 ;                 (P2)  W_s(ζ) ≤ ε(s)(1+ζ²)  for all ζ,   ε(s) ↓ 0 .
> ```
> Let `ζ_node`, `ζ_node^{(s)}` be the peel-offs of `M`, `M_s = M + W_s` as in §1. Then almost
> surely
> ```
>    0  ≤  ζ_node − ζ_node^{(s)}  ≤  ε(s) · K ,
>    K = ⟨1+ζ²⟩_{φ} / G'(ζ_node) < ∞ ,
> ```
> where `φ = φ_{ζ_node}` is the unperturbed bottom eigenfunction on `[ζ_node,∞)` (Dirichlet,
> `L²`-normalised) and `⟨·⟩_φ = ∫ · φ²`. Consequently `ζ_node^{(s)} → ζ_node` almost surely
> (monotonically from below), hence
> ```
>    ζ_node^{(s)}  ⇒  TW_β ,        with a one-signed shift of exact order ε(s).
> ```
> For the Weber curvature term `ε(s) = ¼ s^{-4/3}`, so the peel-off law is `TW_β` with
> `β = 4s/η²` up to an a.s. `O(s^{-4/3})` location shift, **below** the limit (`ζ_node^{(s)} ≤
> ζ_node`).

The hypotheses are exactly met by the folded node: `W_s = ¼s^{-4/3}ζ² ≥ 0` (it is a square —
this **positivity** is the crux) and `≤ ¼s^{-4/3}(1+ζ²)`.

---

## 3. Proof sketch

**(a) Well-posedness `[I + standard]`.** `M` is the RRV stochastic Airy operator (self-adjoint,
a.s. discrete simple spectrum, bounded below); `M_s = M + W_s` adds a non-negative confining
potential, so it is self-adjoint with a.s. discrete spectrum by the same construction
(confining Schrödinger form + RRV noise form). `G, G_s` are well-defined, a.s. continuous,
strictly increasing in `ℓ` (Dirichlet monotonicity). *Invoked:* RRV for the noise form and
spectral type.

**(b) Lower bound — positivity `[R]`.** For any `ℓ`, every Rayleigh quotient obeys
`R_{M_s}[ψ] = R_M[ψ] + ⟨ψ,W_sψ⟩/‖ψ‖² ≥ R_M[ψ]` by (P1). Taking the infimum over the form
domain `Q(M_s) ⊆ Q(M)`,
```
   G_s(ℓ) = inf_{Q(M_s)} R_{M_s} ≥ inf_{Q(M_s)} R_M ≥ inf_{Q(M)} R_M = G(ℓ).
```
Thus `G_s ≥ G` pointwise. Since both are increasing and `G(ζ_node)=0`, we get
`G_s(ζ_node) ≥ 0 = G_s(ζ_node^{(s)})`, hence **`ζ_node^{(s)} ≤ ζ_node`** (one-signed shift).

**(c) Upper bound — trial function `[R]`.** Use the unperturbed eigenfunction `φ_ℓ` (bottom
of `M` on `[ℓ,∞)`) as a trial vector for `M_s`. It lies in `Q(M_s)` because
`⟨φ_ℓ, W_s φ_ℓ⟩ < ∞` (step (e)). Then
```
   G_s(ℓ) ≤ R_{M_s}[φ_ℓ] = G(ℓ) + ⟨φ_ℓ, W_s φ_ℓ⟩ ≤ G(ℓ) + ε(s)⟨1+ζ²⟩_{φ_ℓ}  (by (P2)).
```

**(d) Inversion at level 0 — Hadamard `[R]`.** Combining (b),(c) at `ℓ = ζ_node`:
`0 ≤ G_s(ζ_node) ≤ ε(s)⟨1+ζ²⟩_{φ}`. Since `G_s(ζ_node^{(s)}) = 0` and `G_s` is `C¹` and
increasing, the mean value theorem gives, for some `ξ ∈ (ζ_node^{(s)}, ζ_node)`,
```
   ζ_node − ζ_node^{(s)} = G_s(ζ_node) / G_s'(ξ) ≤ ε(s)⟨1+ζ²⟩_{φ} / inf G_s'.
```
The endpoint derivative is the Hadamard variation `G'(ℓ) = φ_ℓ'(ℓ)² / ‖φ_ℓ‖²`, a.s. strictly
positive (a Dirichlet ground state has `φ_ℓ'(ℓ) ≠ 0`). As `ε(s)→0`, `ζ_node^{(s)}→ζ_node` and
`G_s' → G'` uniformly near `ζ_node` (the perturbation is smooth and `→0` in `C¹_loc`), so
`inf G_s' ≥ ½ G'(ζ_node) > 0` for large `s`. This yields the stated bound with
`K = ⟨1+ζ²⟩_{φ}/G'(ζ_node)`.

**(e) Finiteness of `K` — SAO localisation `[I]`.** Two a.s.-finiteness facts about the SAO
bottom state close the argument:
- `⟨1+ζ²⟩_{φ} = ∫(1+ζ²)φ² < ∞`: the bottom eigenfunction decays like the Airy function,
  `φ(ζ) ≍ exp(−⅔ζ^{3/2})` as `ζ→+∞`, a.s. (RRV; quantitative tails in Dumaz–Virág), so all
  polynomial moments are a.s. finite — indeed `K` has light tails.
- `G'(ζ_node) = φ'(ζ_node)²/‖φ‖² > 0` a.s.: simplicity of the bottom eigenvalue and
  non-degeneracy of the Dirichlet ground state (RRV: a.s. simple spectrum).

Pathwise (same noise `b'` in `M` and `M_s`) this gives `ζ_node^{(s)}(ω) → ζ_node(ω)` for a.e.
`ω`, hence convergence in law `ζ_node^{(s)} ⇒ ζ_node =d TW_β`. ∎(sketch)

**No β-ensemble universality is used.** The argument perturbs the *continuum* SAO directly;
the only probabilistic inputs are RRV (the limit law) and SAO ground-state localisation
(Dumaz–Virág). The positivity (P1) is what makes the two-sided sandwich elementary — see the
remark in §6 for the non-positive case.

---

## 4. Consistency with the numerics (`folded_node_universality.py`)

The lemma predicts, and the simulation confirms, three things the heuristic could not pin:

1. **Sign.** `ζ_node^{(s)} ≤ ζ_node`: the finite-`s` peel-off sits **below** the limit.
   Measured anchors `−2.555, −2.462, −2.409, −2.376, −2.359, −2.350` (`s=1.41…10.98`) all lie
   below `−2.338` and rise to it. ✓
2. **Rate.** `ζ_node − ζ_node^{(s)} = K·¼s^{-4/3}`: a log–log fit of the residual gives slope
   `≈ 1.4 ≈ 4/3`, and the skewness residual to `TW_β` fits the **same** `s^{-4/3}`
   (slope `≈1.45`) — both are first-order in the single perturbation `W_s`. ✓
3. **Even vs odd.** `W_s` is even in `ζ`, so it couples to the location and odd cumulants at
   `O(s^{-4/3})` but to even cumulants only at `O(s^{-8/3})`: the variance converges visibly
   faster (`std → TW₂ = 0.902` by `s≈4.5`, vs skew still climbing at `s≈11`). ✓

So the numerics are not merely consistent with TW_β — they exhibit the lemma's *predicted*
correction (one-signed, `s^{-4/3}`, even/odd split), which is strong evidence the dropped term
is being controlled exactly as claimed.

---

## 5. References to chase (the cite-or-specialise step)

**Needed for the route above (min–max sandwich):**
- **[RRV]** J. A. Ramírez, B. Rider, B. Virág, *Beta ensembles, stochastic Airy spectrum, and
  a diffusion*, J. Amer. Math. Soc. **24** (2011) 919–944. — defines the SAO, proves
  `−Λ₀(β) =d TW_β`, the Riccati/shooting correspondence, self-adjointness, a.s. simple
  spectrum, eigenfunction decay. *(The core citation; supplies §3(a),(e) and `ζ_node =d TW_β`.)*
- **[DV]** L. Dumaz, B. Virág, *The right tail exponent of the Tracy–Widom-β distribution*,
  Ann. Inst. H. Poincaré Probab. Stat. **49** (2013) 915–933. — quantitative SAO ground-state
  tails / localisation; gives that `K` is a.s. finite with controlled tails. *(Supplies §3(e).)*

**Alternative route (if one prefers to reach the SAO by universality rather than perturb it):**
- **[KRV]** M. Krishnapur, B. Rider, B. Virág, *Universality of the stochastic Airy operator*,
  Comm. Pure Appl. Math. **69** (2016) 145–199. — the **most on-point**: the low-lying edge
  spectrum of a class of random Schrödinger/Jacobi operators converges to the SAO. Specialise
  its hypotheses to the Weber-plus-noise operator to conclude `TW_β` without the perturbation
  estimate. *(This is the "specialise" option for §5 of the main note.)*
- **[BEY]** P. Bourgade, L. Erdős, H.-T. Yau, *Edge universality of beta ensembles*, Comm.
  Math. Phys. **332** (2014) 261–353. — general-β edge universality, if routing the Weber
  operator through a matrix model.

**Deterministic underpinning (the `s^{-1/3}` stretch = Airy turning point):**
- **[Olver]** F. W. J. Olver, *Asymptotics and Special Functions* (1974), Ch. 11 (Airy-type
  turning points); or Wasow, *Asymptotic Expansions for Ordinary Differential Equations*.

*Recommendation for Popović:* the min–max route (RRV + DV) is essentially self-contained and
avoids universality machinery; KRV is the fallback that turns the one remaining estimate into
a citation. Either closes §5.

---

## 6. Remarks (scope and the general turning point)

- **Why positivity matters.** (P1) lets the lower bound `G_s ≥ G` be exact (no second-order
  term). For a *general* simple turning point the leading correction is `c·δ²` with
  `c = Q''(z*)/4` of **either sign** (here `Q''=½ > 0`, so `c>0`). If `c<0` (or for the cubic
  correction `Q'''≠0`, which is odd and sign-indefinite), replace (P1) by the two-sided
  relative-form bound `|W_s| ≤ ε(s)(1+ζ²)`; the conclusion still holds, but the proof then
  needs the **a.s. spectral gap** `Λ₁−Λ₀ > 0` (RRV: simple spectrum) to control the
  second-order eigenvalue shift `⟨W_s (M−Λ₀)^{-1}(1−P₀) W_s⟩`. The Weber case sidesteps this.
- **What is genuinely model-dependent.** Only the **location** `z*` (and the `s^{-1/3}` scale)
  carry the Weber/folded-node geometry; the **edge law** `TW_β`, `β=4s/η²`, is universal. This
  is the precise content of "the folded node is in the TW_β class."
- **Caveat retained from the main note.** The reduction to `M_s` assumed the physical
  noise→inner-`η` step (inherited Berglund–Gentz, `[H]` in the main note) and the `D_ν`
  recessive-entry / forgetting lemma `[T]`; this lemma concerns only the final step (the SAO
  edge), and does not re-derive those.
