# The explicit MMO return circle map — built from the funnel + the global return

**What this is.** Pursuing the Tier-2 circle-map idea concretely: *constructing*
the 1-D return map of the folded-node MMO from the two pieces we already have —
the funnel rotation (`ln R = κμ`) and the global return — and reading off what it
predicts. The result is an explicit circle map whose coupling is `K ∼ 1/μ`, which
**derives the subcritical→critical picture, reduces the whole global cluster
(`a_min` #3, `f`, `s_obs`, `α` #6) to a single passage integral `Λ(μ)`, and gives
`α` as that integral's scaling exponent.**

---

## 1. Construction

**State.** Let `d_n` = distance from the weak canard at the `n`-th funnel entry,
and `x_n = ln(a_max/d_n) ≥ 0` (so `d_n = a_max e^{−x_n}`; `x=0` ⇒ enter at the
rim, no SAOs; large `x` ⇒ deep entry, the small `a_min`). Note `d_n` is exactly
the cycle's `a_min`.

**Funnel (inner, have it).** SAOs grow by `R = e^{κμ}` per turn from `d_n` to
`a_max`, so the number of turns is

```
s_n = ln(a_max/d_n)/ln R = x_n /(κμ) ,
```

integer part = SAO count, fractional part = exit phase `φ_n = {x_n/(κμ)}`. The
trajectory leaves the funnel at `O(a_max)` and fires the L-spike.

**Global return (the unknown, packaged).** The L-spike + reinjection contracts the
trajectory back toward the weak canard. Write the new entry log-depth as a mean
contraction plus a phase-dependent part:

```
x_{n+1} = Λ(μ) + δ·h(φ_n) ,     h = O(1) periodic,  δ = phase-variation amplitude,
```

where `Λ(μ) = ⟨ln(a_max/d)⟩` is the **mean global-return contraction toward the
weak canard** — a single slow-passage integral.

**The map.** Composing, the phase obeys a standard **circle map**:

```
┌────────────────────────────────────────────────────────────┐
│   φ_{n+1} = φ_n + Ω + K·h(φ_n)  (mod 1)                       │
│   Ω = Λ(μ)/(κμ)   (bare rotation number)                      │
│   K = δ/(κμ)      (coupling)                                  │
└────────────────────────────────────────────────────────────┘
```

The MMO mode-locking *is* this map's Arnold-tongue structure; `s_obs` is its
rotation count; `f` its locked fraction; `σ_pq` (#4) its noise-unlocking.

## 2. Everything reduces to one integral `Λ(μ)`

At the map's fixed point `x* ≈ Λ` (for moderate `δ`):

```
s_obs = Λ/(κμ) ,    f = s_obs/s_max = (Λ/(κμ))/((1−μ)/(2μ)) ≈ Λ/π² ,
a_min = a_max e^{−Λ} .
```

So **`a_min`, `f`, and `s_obs` are all just `Λ(μ)`** — the global-return
contraction. The inner physics has been fully absorbed into `κ` (→ `π²`); the
single remaining function is `Λ(μ) = ln(a_max/a_min)`, the slow-passage
contraction toward the weak canard.

## 3. The coupling is `K ∼ 1/μ` — this derives the subcritical→critical picture

The slow passage near the folded node lasts a time set by the weak eigenvalue,
`T ∼ 1/λ_w ∼ 1/(μλ_s)`, so the contraction and its phase-variation both grow as
`μ→0`: `Λ ∼ 1/μ` and `δ = O(1)` of it, giving

```
K = δ/(κμ) ∼ 1/μ      (coupling grows as the folded node → folded saddle-node).
```

Therefore:

- **Data band (`μ ∈ [0.014,0.049]`): `K < 1`, subcritical** — incomplete
  staircase, `f < 1`, exactly the measured partial filling.
- **`μ → 0`: `K → 1`, criticality** — staircase completes, `f → 1`, and
  `a_min → a_max e^{−π²} ≈ 5×10⁻⁵` (the universal floor from `MMO_AMIN_TIER1.md`).
- The three coincide at one point: `K=1 ⇔ f=1 ⇔ Λ=π²`, predicted at `μ ≈ 0.0064`
  (`c ≈ −0.97`), just beyond the data — with phase-variation `δ ≈ κ·0.0064 ≈ 0.13`.

So the folded-saddle-node limit is literally an **approach to circle-map
criticality**, with an explicit coupling `K(μ)=δ/(κμ)`. This is the mechanism, not
an analogy.

## 4. `α` is the scaling exponent of `Λ(μ)`

The rotation number is `ρ = 1/(s+1) ≈ κμ/Λ`. With `Λ ∼ μ^{−r}`:

```
s_obs ∼ μ^{−(1+r)} ,   ρ ∼ μ^{1+r} ∼ (c+1)^{1+r} ,   ⇒   α = 1 + 1/(1+r) .
```

The clean limits bracket it:

| `Λ(μ)` scaling | origin | `α = 1+1/(1+r)` |
|---|---|---|
| `Λ ∼ 1/μ`  (`r=1`, weak-eigenvalue passage time) | folded node | **3/2 = 1.50** |
| `Λ ∼ 1/√μ` (`r=½`, saddle-node bottleneck `1/√(dist)`) | folded SN | **5/3 ≈ 1.67** |

Measured `α ≈ 1.55–1.59` (`s_obs ∼ μ^{−1.7}` fits the table to ~10%, vs ~50% for
`μ^{−2}`) sits **between** these — i.e. `Λ ∼ μ^{−r}` with `r ≈ 0.7`, a
*transition-regime* exponent, exactly as expected for `μ < ε^{1/2}` where the
folded node is crossing over to the folded saddle-node. So the circle map explains
*why `α` is not a clean universal exponent*: the coupling `K∼1/μ` is **en route to
criticality**, not at it, so the staircase exponent is a crossover value pinned by
`Λ`'s transition-regime scaling.

## 5. What this buys (and the honest remainder)

**Bought:**
- An **explicit return circle map** `φ→φ+Ω+K h(φ)` with derived `Ω=Λ/(κμ)`,
  `K=δ/(κμ)`, built from the funnel (`κ`) and the global return (`Λ,δ`).
- The whole global cluster **reduced to one integral** `Λ(μ)=ln(a_max/a_min)`:
  `a_min=a_max e^{−Λ}`, `f=Λ/π²`, `s_obs=Λ/(κμ)`, `α=1+1/(1+r)` with `Λ∼μ^{−r}`.
- A **derived mechanism** for subcritical→critical (coupling `K∼1/μ`), with the
  funnel-fill / floor / criticality all coinciding at `Λ=π²` (`μ≈0.0064`).
- `α`'s non-universality **explained**: `K` is approaching criticality, so `α` is
  a crossover (measured `≈1.55–1.59`, between the `3/2` and `5/3` limits).

**Honest remainder (the one integral):** `Λ(μ)` itself — the mean slow-passage
contraction toward the weak canard — is the model-specific global return. Its
transition-regime scaling `Λ∼μ^{−r}`, `r≈0.7`, is the single open quantity; it
sets `α`, `f`, and `a_min` together. This is genuine progress: a four-constant
cluster (#3, #6, #4, plus `f`) is now **one passage integral** `Λ(μ)`, with the
inner reduced to `π²` and the structure (circle map, criticality) explicit.

**Noise side (#4) for free:** `σ_pq` is this map's **noise-unlocking** — when noise
pushes `φ` out of a locked interval, `σ_pq ~ σ_* q^{−α/γ}`. Since `K∼1/μ`, the
tongues are thinner (easier to unlock) at larger `μ`, which is the μ-dependence of
the dissolution threshold — now a property of the *same* map, in the
stochastic-mode-locking reference class.

## 6. The concrete next step

Two clean moves, both about the single object `Λ(μ)`:

1. **Measure `Λ(μ) = ln(a_max/a_min)` and fit `r`** (we have the data: `r≈0.7`,
   `α≈1.59`) — and, decisively, **sweep `ε`** to test whether `Λ` is a function of
   `μ` alone or of `μ/ε^{1/2}` (the transition variable). That settles `r` and the
   `α` value.
2. **Extract the map's coupling `K(μ)`** from the deterministic return (the spread
   of reinjection depth `δ` vs `κμ`) and verify `K→1` as `μ→0` — confirming the
   criticality mechanism directly.

Both are measurements on the existing flow, not new theory; together they close
`α` (#6) and `f`/`a_min` (#3) as outputs of one explicit circle map, and hand #4
to the stochastic-mode-locking framework.
