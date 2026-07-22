# Resonator commitment angle (#10) — how far the instanton closes, and why it stops

**What this is.** An honest attempt to derive the `−1.06 rad` spike-commitment
offset (`RESONATOR_MECHANISM.md`, currently gMAM-numerical) from first principles
— a deliberately *different* vein from the fold-escape / folded-node work
(Hamiltonian large-deviation instanton, degenerate noise). The verdict up front,
because it's the useful part:

> The magnitude `−1.06 rad` does **not** reduce to elementary closed form — it is
> a global nonlinear Hamiltonian heteroclinic (focus → unstable limit cycle), so
> gMAM is the right tool and the audit's "NUM" tag is correct. But the attempt is
> not empty: it **derives two of the chapter's empirical findings** (tests 2 and
> 4) from the linear instanton, and pins the exact framework + the mechanism and
> sign of the offset. So #10 is genuinely numerical *for a reason now understood*.

---

## 1. The framework (degenerate-noise Freidlin–Wentzell)

Resonator FHN, noise in `v` only, near the stable focus `x*=(v*,w*)=(−0.993,−0.367)`:

```
dv = b_v dt + σ dW ,   b_v = v−v³/3−w+I ;     dw = b_w dt ,   b_w = ε(v+a−βw).
```

The FW Hamiltonian (degenerate — only `p_v²` appears, because there is no `w`-noise):

```
H(x,p) = p_v b_v + p_w b_w + ½ p_v² .
```

The instanton is the zero-energy heteroclinic of Hamilton's equations from
`(x*,0)` to the threshold, with

```
v̇ = b_v + p_v ,   ẇ = b_w ,   ṗ = −(∂b/∂x)ᵀ p .
```

Note `ẇ = b_w` *exactly* — the degenerate noise forces the `w`-channel onto its
noiseless ODE (the metric-penalty trick in `instanton_exit.py` enforces this).
The extra term `p_v` in `v̇` is the noise push; it is what tilts the escape away
from the deterministic flow.

## 2. The focus quasipotential is anisotropic and `v`-soft (computed)

Linearise: `A = [[1−v*², −1],[ε, −εβ]] = [[0.014, −1],[0.08, −0.064]]`
(eigenvalues `α±iω = −0.025 ± 0.280 i` — matches the doc). The local
quasipotential is `V(x)=½ xᵀΣ₁⁻¹x` with `Σ₁` solving the Lyapunov equation
`AΣ₁+Σ₁Aᵀ+diag(1,0)=0`:

```
Σ₁ = [[10.52, 0.647],[0.647, 0.809]] ,   Σ₁⁻¹ = [[0.100, −0.080],[−0.080, 1.301]] .
⇒ V(x) = ½(0.100 v² − 0.160 vw + 1.301 w²)   (x relative to the focus).
```

The `v`-direction is **soft** (coefficient `0.100`), the `w`-direction **stiff**
(`1.301`): escaping costs ~13× less per unit displacement in `v` than in `w`,
because noise only pushes `v`. This anisotropy is the whole story — it is what
makes the cheapest escape *not* the geometrically steepest one.

## 3. The linear instanton is an angle-less focus — this derives tests 2 and 4

The most-probable escape path near the focus follows `ẋ = M x` with

```
M = A + diag(1,0)·Σ₁⁻¹ = [[0.114, −1.080],[0.08, −0.064]] = −Σ₁ Aᵀ Σ₁⁻¹ .
```

The identity `M = −Σ₁AᵀΣ₁⁻¹` (immediate from the Lyapunov equation) means `M` is
similar to `−Aᵀ`, so its eigenvalues are **exactly the negatives of `A`'s**:

```
eig(M) = −(α ± iω) = +0.025 ∓ 0.280 i   (unstable focus, ω unchanged).
```

The instanton is the **time-reversed relaxation**: it spirals *outward* at the
*same* frequency `ω`. Two consequences, both matching the doc's empirical tests:

- **Winding is cost-free (test 4).** The outbound instanton is a focus — it can
  wind any number of turns following the (reversed) drift at zero extra action.
  This is exactly gMAM's "init-independent to `0.008 rad`, independent of winding."
- **The offset is not linear resonance / not `κ` (test 2).** A focus has **no
  distinguished exit angle** — the linear instanton spirals through *all* angles.
  Therefore the committed direction cannot come from the linear focus data
  (`α, ω, κ`); the linear theory provides the winding but carries *zero*
  directional information. The doc found this empirically ("offset invariant in
  `κ`, corr `−0.40`"); here it is a theorem about `M`.

So the linear theory *explains the absence* of a linear answer — the committed
angle must be set by the nonlinearity. This is the same shape of result as the
resonator barrier (`BARRIER_NORMAL_FORM.md`): the linear/normal-form object
explains why the naive scaling fails, rather than supplying the number.

## 4. The offset is irreducibly nonlinear — mechanism, sign, and the route

What sets the committed angle is the **nonlinear** matching of the `v`-soft
outward spiral to the threshold (the subcritical-Hopf unstable cycle,
`r_u² = |α|/a₃ = 0.025/0.268 ⇒ r_u ≈ 0.31`). Mechanism:

- The instanton velocity carries an extra `+p_v` in the `v`-component (§1). The
  cheapest escape therefore *biases the outbound leg toward the `v`-axis* (the
  soft direction of `Σ₁`), not the deterministic steepest-outward direction.
- Because the path also rotates at `ω`, this `v`-bias is carried to a committed
  crossing angle **offset to one rotational side** of the geometric direction —
  the **sign** of the offset is the sign of the rotation `ω` relative to the
  `v`-axis orientation in the eigenbasis (here negative, consistent with `−1.06`).
- The bias is `O(1)` (not `O(κ)`), which is why the offset is large (`~1 rad`) and
  `κ`-independent.

**The tractable route to the magnitude** (not elementary, but the right reduction)
is the subcritical-Hopf normal form `ż = (α+iω)z + (a₃+ib₃)|z|²z` with the
`v`-only noise projected onto `z` — which makes the noise efficiency
**phase-dependent** (`dv` pushes `r` and `θ` by amounts `∝cos`/`sin(θ+φ_v)`). The
instanton is then a heteroclinic in `(r,θ,p_r,p_θ)` from `r=0` to `r=r_u`; the
phase-dependent noise term is exactly what selects the committed `θ`. This 2-DOF
Hamiltonian BVP is far cheaper than the full FHN gMAM but still has **no
elementary closed-form** solution — the committed angle is a transcendental
functional of `(a₃,b₃,φ_v,κ)`.

## 5. Honest status

| Piece | Status |
|---|---|
| Degenerate-noise FW Hamiltonian + instanton ODEs | framework, exact |
| Anisotropic focus quasipotential `Σ₁` (`v`-soft, 13×) | **derived (numbers)** |
| Linear instanton `M=−Σ₁AᵀΣ₁⁻¹`, eig `−α±iω` | **derived** |
| ⇒ winding cost-free (test 4) + offset not `κ` (test 2) | **derived from M** |
| Offset mechanism (`v`-bias × rotation) + sign | **derived (qualitative)** |
| **Magnitude `−1.06 rad`** | **open — nonlinear heteroclinic; gMAM is correct** |

**Verdict.** #10 does not close in the clean register — and now we know *why*: the
linear instanton is an angle-less unstable focus, so the commitment angle is a
genuinely nonlinear large-deviation quantity (the focus→unstable-cycle
heteroclinic with phase-dependent degenerate noise). The attempt is still a net
gain: it derives tests 2 and 4 analytically, computes the `v`-soft quasipotential,
and reduces the magnitude to a 2-DOF normal-form Hamiltonian BVP — the right
object for a future closed-ish or semi-analytic attack, but not pen-and-paper.
This is the honest boundary of the fold/focus program: the **fold** constants
(`C_q`, excitable `C`, tonic `c`) are FW *quasipotential* objects that collapse
cleanly; the **focus commitment angle** is an FW *minimum-action-path direction*,
which does not.
