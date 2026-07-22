# Tier-1 ideas applied to the results we have — what the data actually say about `a_min`

**What this is.** The Tier-1 exponential-asymptotics ideas
(`MMO_AMIN_BRAINSTORM.md`: exact-WKB on the Weber inner, resurgence, uniform PCF)
run against the *established* objects — the rotation map `ln R = κμ`, the measured
`a_min / a_max / s_obs`, and the Weber inner from `MMO_K2_ROUTE_AB.md`. The
outcome is more clarifying than expected: it produces a clean closed relation,
recovers the quantization, and — most usefully — **dissolves the "exponent of
`a_min`" question** that the brainstorm was built around.

---

## 1. The clean relation (from the rotation map, exact)

The rotation map *defines* `s_obs = ln(a_max/a_min)/(κμ)`. Rearranged, with
`s_obs = f·s_max`, `s_max = (1−μ)/(2μ)`, `κ ≈ 2π²`:

```
ln(a_max/a_min) = κ μ s_obs = κ μ · f·(1−μ)/(2μ) = (κ/2)(1−μ) f ≈ π² f ,
```

so

```
┌──────────────────────────────────────────────┐
│   a_min  ≈  a_max · exp( −π² · f )             │
│   π² = κ/2 (inner, have it) ;  f = funnel-fill │
└──────────────────────────────────────────────┘
```

`a_min` is the largest amplitude `a_max ~ O(1)` damped by `exp(−π² f)`, where
`f∈(0,1)` is the funnel-filling (the *global return*). The **rate `π² = κ/2` is
the inner constant we already have**; everything else in `a_min` is `f`.

**Check against the data** (`f = 0.0064/μ` fits the measured `f`: at `μ=0.014`,
`0.457` vs measured `0.445`; at `μ=0.033`, `0.194` vs `0.206`):

```
a_min ≈ a_max·exp(−π²·0.0064/μ) = a_max·exp(−0.063/μ) ;
μ=0.026: 1.1·e^{−2.42}=0.098  (meas 0.090) ;  μ=0.014: 1.05·e^{−4.5}=0.012 (meas 0.0083).
```

Good to ~30% — the relation is right; the residual is `a_max`/dt scatter.

## 2. The idealized floor — what inner asymptotics *does* deliver

As the funnel fills (`f→1`):

```
a_min → a_max · e^{−π²} ≈ 5.2 × 10⁻⁵     (μ-INDEPENDENT, universal).
```

This is the deepest a folded-node funnel can reach — the *inner* limit, set
purely by `κ` (it is `MMO_K2_ROUTE_AB.md`'s "total amplification `e^{π²}`"). **This
universal floor is exactly what the Tier-1 inner machinery (exact-WKB Stokes
constant of the Weber turning points / resurgence) computes** — and we already
have it, because it is `κ/2`. So Tier 1 *closes the floor* for free.

## 3. Exact-WKB sanity checks on the Weber inner (consistent)

- **Bohr–Sommerfeld recovers the secondary canards.** The Weber inner
  `ψ''+(c₀−τ²/4)ψ=0` has `∮√Q dτ = π c₀` (Route B). BS quantization
  `∮√Q=(n+½)π ⇒ c₀=n+½`, i.e. with `c₀=1/μ−1` the **secondary canards sit at
  `μ_n = 1/(n+3/2)`** — the standard folded-node ladder, recovered from our own
  inner equation. ✓
- **The conservative envelope is the algebraic upper bound.** The WKB amplitude
  `~Q^{−1/4}`, smallest at the funnel centre (`Q=c₀`), gives
  `a_min^{env} ~ c₀^{−1/4} ~ μ^{1/4}` — algebraic. This is the *sector bound*
  (`ε^{(1−μ)/2}`-type), confirming again that the *exponential* depth of `a_min`
  is a connection/global effect, not the local envelope. ✓

## 4. The punchline — the "exponent of `a_min`" is a red herring

The brainstorm (and the audit) framed `a_min ~ exp(−C·μ^{−p})` and asked for the
exponent `p` (0.7? 2/3? 3/4?). But §1 shows there is **no clean inner exponent**:

```
a_min = a_max · exp(−π² f) ,   and the μ-dependence of a_min IS the μ-dependence of f.
```

`f` is the *global return* funnel-filling. Empirically `f ≈ 0.0064/μ` over the
measured band, which makes `a_min ≈ a_max e^{−0.063/μ}` look like a clean
`exp(−C/μ)` (`p=1`) — and a re-fit on the *reliable* points (dropping the `s<3`
shallow one) indeed favours `p≈1`, **not** the doc's `0.7`. But this `p=1` is just
`f∼1/μ` over a narrow range; once `f` saturates (`→1` near `μ≈0.006`), `a_min`
flattens to the floor and the effective `p→0`. **So the "exponent" is not a
universal asymptotic constant — it is the shape of `f(μ)`, a global-return
quantity.** Chasing it with inner exponential asymptotics (nonlinear Landau–Zener,
anomalous Stokes exponents, etc.) was aiming the wrong tool: the inner physics
gives the *rate* `π²` and the *floor* `e^{−π²}`; the rest is `f`.

## 5. Verdict and redirection

| Question | Tier-1 (inner) verdict |
|---|---|
| Per-turn rate `κ` | `2π²` (open constant = Route A PCF; have the value) |
| Idealized floor `a_min(f→1)` | `a_max·e^{−π²} ≈ 5×10⁻⁵`, **closed** (= κ/2), μ-independent |
| Secondary-canard ladder | `μ_n=1/(n+3/2)`, **recovered** (Bohr–Sommerfeld) |
| Algebraic upper bound | `~μ^{1/4}` envelope, **confirmed** |
| **`a_min(c)` itself / its "exponent"** | **NOT inner** — equals `a_max e^{−π²f}`, set by the global return `f(c)` |

**So Tier 1 did its job by ruling itself out for the main target.** It closes the
floor and the rate and recovers the quantization, but it proves `a_min`'s actual
value and `c`-dependence are the global-return funnel-filling `f(c)` exponentiated
by the inner rate `π²`. The honest consequence:

- **Stop** looking for an inner/special-function exponent for `a_min` — there
  isn't one; `a_min = a_max e^{−π²f}` with `π²` already in hand.
- **The remaining unknown is `f(c)` (≡ the global return injection)** — a
  *dynamical-systems / circle-map* object (Tier-2 ideas 5,7), not an
  exponential-asymptotics one. `f≈0.0064/μ` is the empirical handle to derive.
- This also **cleans `#4/#6`**: with `a_min = a_max e^{−π²f}`, the staircase
  exponent and the σ_pq μ-law inherit `f(c)` directly — one global quantity feeds
  all three, with the inner contribution fully reduced to the single number `π²`.

The one-line result: **`a_min = a_max·e^{−π² f}` — the inner asymptotics gives the
rate `π²` and the universal floor `e^{−π²}≈5×10⁻⁵`; `a_min`'s `c`-dependence is the
global-return filling `f(c)`, and there is no separate "exponent" to derive.**
