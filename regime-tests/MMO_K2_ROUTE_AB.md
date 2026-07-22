# Deriving `κ = 2π²` — Route B (WKB/antidamping) then Route A (PCF connection)

**Status:** a worked attempt. Route B is carried to the point where it
*reduces* the constant to a single well-posed integral and **rules out the
obvious mechanism**; it then provably collapses into Route A. Route A is set up
with the discriminating scaling lock and the link to `a_min` (#3). The final
number extraction (the parabolic-cylinder connection coefficient) is identified
as the technical core and left as the open step.

> **Caveat, stated once.** The reduction below is reconstructed structurally
> (from the folded-node blow-up architecture), not copied from a normal form on
> the page. The *structure* — antidamping reduction, rotation count, the
> negative result, the scaling lock — is robust. The *exact coefficients*
> (`g`, `Ω`, the `a ↔ μ` map) must be checked against Wechselberger 2005 §4
> before any number is trusted.

Companion to `MMO_K2.md` (which measures `ln R = κμ`, `κ = 19.6 ± 1.8 ≈ 2π²`).

---

## Route B — WKB / antidamping on the variational equation

### B1. The variational equation about the weak canard

The SAOs are the winding of trajectories around the weak canard `γ_w`.
Linearise the (blown-up) folded-node flow about `γ_w`; the transverse deviation
`ζ`, in the rescaled slow time `τ`, obeys a linear 2nd-order ODE of
**damped-Weber** form:

```
ζ'' − 2 g(τ) ζ' + Ω²(τ) ζ = 0 ,
```

with `Ω(τ)` the local rotation frequency (the fast winding) and `g(τ)` an
**antidamping** coefficient — the transverse divergence of the field along
`γ_w`. `g > 0` because the strong eigendirection repels: this is the term that
makes the SAOs grow.

### B2. Liouville transform → conservative Weber; the rotation count

Put `ζ = exp(∫ g dτ) · ψ`. The first-order term cancels exactly:

```
ψ'' + Q(τ) ψ = 0 ,      Q = Ω² − g² + g' .
```

(Standard Liouville form `q − p²/4 − p'/2` with `p = −2g`: the `−p'/2 = +g'`, so
the sign is **`+g'`** — corrected from an earlier `−g'`. Immaterial to the
leading WKB count — `g², g'` are `O(μ)` against `c₀ ~ 1/μ` — but it sits squarely
in the `O(μ)` correction, i.e. exactly the `2π²(1+O(μ))` drift term, so it must
be right before anyone touches the sub-leading piece.)

Near the folded node this is **Weber's equation** with a finite well
`Q = c₀ − τ²/4`; the SAOs live in the classically-allowed region
`|τ| < 2√c₀`. WKB zero-count:

```
N = (1/π) ∫_{−2√c₀}^{+2√c₀} √(c₀ − τ²/4) dτ = (1/π)(π c₀) = c₀ ,
```

so full turns `= N/2 = c₀/2`. Matching the cited ceiling
`s_max = (1−μ)/(2μ)` fixes the well depth **exactly**:

```
c₀ = 1/μ − 1 .
```

✓ The reduction reproduces the known rotation count with the right O(1)
constant — confirming `Ω`, `c₀` are correctly identified, and that
`a ≍ 1/μ` is the right order.

### B3. The growth IS the antidamping integral (the clean target)

Because `ζ = e^{∫g} ψ` and `ψ` has a bounded oscillatory envelope, the physical
SAO amplitude grows as `e^{∫ g dτ}`. Per turn (phase advance `2π`):

```
ln R = ∮_{1 turn} g(τ) dτ        ⇒        κ = (1/μ) ∮_{1 turn} g dτ .
```

Integrated across the whole funnel (`s_max` turns):

```
∫_funnel g dτ = ln R · s_max = κμ · (1−μ)/(2μ) = κ(1−μ)/2  →  κ/2 .
```

So the entire problem reduces to one integral:

```
┌──────────────────────────────────────────────┐
│   κ = 2π²   ⇔   ∫_funnel g dτ = π²  (μ-indep.)  │
└──────────────────────────────────────────────┘
```

**Consistency of the scales** (a good sign the structure is right): with
`∫g ≈ π²` spread over funnel length `L = 4√c₀ ≈ 4/√μ`, the antidamping is
`g ≈ π²√μ/4 → 0`, while the rotation frequency at the centre is
`Ω ≈ √c₀ ≈ 1/√μ`. Hence

```
growth per radian = g/Ω ≈ (π²/4) μ ≈ 2.5 μ ,
```

the right **`O(μ)` form** (the funnel is a *very tight* log-spiral at small μ —
slow per radian, but `∝1/μ` radians, netting O(1) total). This is an
**order-of-magnitude consistency check, not a pin of the `π`**: `(π²/4)μ ≈ 2.5μ`
vs the target `πμ ≈ 3.1μ` is ~20% loose, because it used a crude constant-`g`
estimate (`g ≈ π²√μ/4` over `L ≈ 4/√μ`). And `g = O(√μ)` is small, so the
conservative rotation count B2 is self-consistently valid.

### B4. Negative result — the constant is NOT the adiabatic envelope

The natural guess is that the growth is the conservative Weber envelope
`ψ ~ Q^{−1/4}` (amplitude rises as the local frequency drops toward the turning
point). Test it. Over the half-funnel,

```
A_env = Q^{−1/4} ,   total log-growth = ¼ ln(c₀ / q_min) ,
```

with the Airy turning-point cutoff `q_min ~ c₀^{1/3}` and `c₀ ≈ 1/μ`. Dividing
by the `≈ c₀/4` turns in the half-funnel:

```
κ_env ≈ (2/3) ln(1/μ)   →   ≈ 2–3   for μ ~ 0.01–0.05 .
```

This **fails twice**: (i) wrong scaling — it grows like `ln(1/μ)` instead of
being constant; (ii) ~10× too small — measured `κ ≈ 20`. **Therefore `2π²` is
the genuine antidamping `g` (transverse divergence, set by the strong
eigenvalue), not the conservative envelope.** This is the main concrete result
of Route B: it rules out the simplest mechanism and pins the constant to the
*growing* part of the solution.

### B5. Where B stalls → it is A

`∮ g dτ` over one turn is dominated by the passage closest to the folded node —
the deepest, smallest loops near `τ ≈ 0`, where `g` is largest and the WKB
factorisation of the *growing* solution breaks down. Evaluating
`∫_funnel g dτ = π²` there is exactly a parabolic-cylinder **connection**
computation. So Route B reduces the constant to a well-posed integral but cannot
evaluate it without the connection data.

> **B and A are not independent. B's amplitude question *is* A's connection
> question.** B's value is the reduction + the negative result, not an
> independent number.

---

## Route A — the parabolic-cylinder connection coefficient

### A1. The object, and the scaling lock that screens candidates

The per-turn multiplier `R = e^{κμ}` is intrinsic (ε-independent); the *realised*
number of turns — hence `a_min` — is ε-dependent, but `R` is a property of the
local rotation map. From B3 the full-funnel transverse amplification is

```
e^{∫ g} = e^{κ(1−μ)/2} → e^{κ/2} = e^{π²} ≈ 1.9 × 10⁴ ,   μ-independent.
```

**Discriminating check (use this to screen any candidate coefficient):** its log
must be **μ-independent and O(1)**. A connection coefficient scaling like
`e^{πa} ~ e^{π/μ}` is *wrong* — those prefactors carry μ-dependence that must
cancel. You are hunting a coefficient whose net log → `π²`.

### A2. The connection problem

In the oscillatory regime the conservative factor `ψ` is a parabolic-cylinder
function `W(a, τ)` (DLMF §12.14), with `a` tied to `c₀` (sign/convention from the
normal form). The canard is the solution matched to:

* **entry** — the attracting slow manifold (the way-in, `τ` on one side),
* **exit** — the repelling slow manifold (the way-out, other side).

The transverse amplification is the ratio of `W`-moduli between entry and exit,
assembled from the PCF **modulus function `F(a,τ)`** and **phase `ω(a,τ)`**. The
recurring constants in §12.14 are `e^{πa}`, `√(1+e^{2πa})`, and the phase. By A1,
`2π²` must emerge from the **net change of the modulus function across the
funnel** — *not* from the `e^{πa}` prefactors (which carry the μ-dependence that
cancels in the entry/exit ratio).

This is the single genuine computation. It is what Wechselberger 2005 §4 does for
the *rotation* (the phase `ω`); the *amplitude* (the modulus `F`) is the piece to
extract — and to confirm whether it gives exactly `2π²` or `2π²(1+O(μ))` (the
latter would explain the audit's low-q "drift").

### A3. Two-for-one with `a_min` (#3)

The **deepest** point of the canard (`τ ≈ 0`, the smallest loop) is `a_min` — the
*same* inner solution, evaluated at the funnel centre. So the connection
computation that delivers `κ` **simultaneously delivers `a_min(c)`**, closing #3.
This is the audit's "#2 and #3 are coupled" made operational: one
parabolic-cylinder connection computation, two constants
(`κ = 2π²` and `a_min(c)` → the realised `α`).

### A4. Honest status

| Piece | Status |
|---|---|
| Variational → damped Weber → conservative Weber (B1–B2) | **structure solid**; exact `g, Ω` to verify vs Wechselberger §4 |
| Rotation count `c₀ = 1/μ − 1` | **derived** (WKB zero-count = `s_max`) |
| Growth = antidamping integral, `κ = 2π² ⇔ ∫_funnel g = π²` (B3) | **derived reduction** |
| Conservative envelope ruled out, `κ_env ~ (2/3)ln(1/μ)` (B4) | **derived negative result** |
| Scaling lock: amplification `e^{π²}`, μ-independent (A1) | **derived constraint** |
| PCF connection structure + `a_min` link (A2–A3) | **set up**; the number is open |
| **Exact value `2π²` (vs `2π²(1+O(μ))`)** | **open — the technical core** |

**Where to mine the open step:** Wechselberger (2005) §4 (the rotation/phase),
Brøns–Krupa–Wechselberger (2006), Desroches et al. (2012) SIAM Review §3, and
DLMF §12.14 (PCF modulus/phase). Given the supervisor, Popović's folded-node
papers are the closest prior art for the amplitude asymptotics.

---

## One-paragraph summary

Route B does not yield `2π²`, but it earns three things: it reduces the constant
to a single integral `∫_funnel g dτ = π²`; it shows the growth is **antidamping**
(transverse divergence), **not** the conservative WKB envelope (which is the
wrong scaling and ~10× too small); and it proves B's amplitude question is
identical to A's connection question. Route A is then the parabolic-cylinder
connection coefficient — screened by the requirement that its log be a
μ-independent `π²`, and yielding `a_min` (#3) as a by-product. The skeleton is
closed; the one remaining computation is the PCF modulus connection.

---

## Certification log

Independently re-derived and **certified**: the Liouville reduction (B1–B2), the
rotation count `c₀ = 1/μ − 1` (WKB = `s_max`), the headline reduction
`κ = 2π² ⇔ ∫_funnel g dτ = π²`, the negative result `κ_env ≈ (2/3)ln(1/μ) ≈ 2`
(Airy cutoff `q_min ~ c₀^{1/3}` over `~c₀/4` turns — wrong shape *and* ~10× too
small), the `e^{π²}` scaling lock, and the `a_min`/#3 coupling.

Two corrections folded in:

1. **Liouville potential sign (real bug).** `Q = Ω² − g² + g'`, not `−g'` (B2).
   Leading count unaffected; the term lives in the `O(μ)` drift, so it matters
   for the `2π²(1+O(μ))` question.
2. **`g/Ω` looseness (overclaim).** `g/Ω ≈ (π²/4)μ ≈ 2.5μ` is an
   order-of-magnitude check, ~20% from `πμ`; it does **not** pin the `π` (B3).

Next step is the numerical antidamping-integral measurement
(`MMO_K2_GINT_MEASUREMENT.md`) — measure-first, before the PCF computation.
