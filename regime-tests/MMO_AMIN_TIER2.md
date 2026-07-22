# Tier-2 ideas applied to `f(c)` — the circle-map home for the global return

**What this is.** Tier-1 (`MMO_AMIN_TIER1.md`) reduced `a_min` to one object —
`a_min = a_max·e^{−π²f}`, with the inner rate `π² = κ/2` already in hand and the
**funnel-filling `f(c)`** (the global return) the sole remaining unknown. Tier-2
is the toolbox for exactly that kind of *global* object. Running the Tier-2 ideas:

> The right one is the **circle-map / Arnold-tongue** framing (with RG as its
> renormalisation). It does not hand over `f(c)`'s constant — that is the
> model-specific global return — but it gives the correct *home* for the whole
> global MMO side: it **unifies `f`, `s_obs`, `α` (#6) and the `σ_pq` μ-law (#4)
> as one mode-locking structure**, gives a falsifiable subcritical→critical
> picture, and relocates #4 into the well-studied *noisy mode-locking* problem.
> The inner Tier-2 ideas (Painlevé, nonlinear Landau–Zener, separatrix splitting)
> are ruled out for `f` by Tier-1 and are redirected to #2.

---

## 1. First, rule out the inner Tier-2 ideas (Tier-1 already did)

Tier-1 showed the inner contribution to `a_min` is *entirely* the linear rate
`π² = κ/2` (Weber/parabolic-cylinder, a **linear** connection), and the universal
floor `e^{−π²}`. So:

- **Painlevé II / RH (idea 4)** assumed a *nonlinear* inner connection. There
  isn't one — the inner is the linear Weber equation. So PII is not needed for
  `a_min`. *Redirect:* PII/RH machinery is the natural rigorous route to the
  *linear* PCF connection constant behind **κ = 2π² (#2 Route A)** — use it there.
- **Nonlinear Landau–Zener (idea 6)** was to explain an *anomalous* exponent.
  Tier-1 showed there is no anomalous inner exponent (`a_min`'s μ-dependence is
  `f(μ)`, global). Drop it.
- **Separatrix splitting (idea 8)** computes the inner secondary-canard
  `exp(−c/μ)` — i.e. the *floor*, which we already have as `e^{−π²}`. Redundant.

That leaves the **circle-map / Arnold-tongue (5)** and **funnel RG (7)** ideas —
the genuinely *global* tools — for `f(c)`.

## 2. The reframing: the MMO staircase IS circle-map mode-locking

The folded-node return map, after the strong Fenichel contraction onto the slow
manifold, is effectively **1-D**, and the funnel rotation makes it a **circle
map** `θ ↦ θ + Ω(c) + (nonlinearity)`. In this language the project's objects are
all standard circle-map quantities:

| MMO object | circle-map object |
|---|---|
| `L¹Sˢ` pattern, `ρ = 1/(s+1)` | mode-locked state, rotation number `ρ` |
| plateau in `c` of width `Δ_pq ~ q^{−α}` (#6) | **Arnold tongue** width |
| funnel-filling `f(c)` (the `a_min` unknown, #3) | **fraction of `c`-axis that is mode-locked** |
| `σ_pq` dissolution threshold (#4) | **noise-induced unlocking** of a tongue |

So `f`, `s_obs`, `α` (#6) and the `σ_pq` μ-law (#4) are **not four problems** —
they are four readouts of one circle map. That is the main structural payoff:
the entire global MMO side collapses to "what is the return circle map, and how
locked is it?"

## 3. The qualitative picture (new, and falsifiable)

The measured funnel is *partially filled*, `f ≈ 0.0064/μ ∈ [0.14, 0.48]` over the
band (Tier-1). In circle-map terms:

- `f < 1` (incomplete staircase, thin tongues) ⟺ a **subcritical** (smooth,
  invertible) circle map — the locked set has measure `< 1`.
- `f → 1` (complete staircase) ⟺ approach to **criticality** (the locked set
  fills the axis). With `f ≈ 0.0064/μ`, this is predicted at `μ ≈ 0.0064`
  (`c ≈ −0.97`), *beyond* the present data — and exactly where Tier-1 said `a_min`
  hits its universal floor `e^{−π²}`.

So **the `μ→0` (folded-saddle-node) limit is an approach to circle-map
criticality**, and the data sit in the subcritical, partially-locked regime. This
is a concrete, testable claim: extract the return map at several `c` and check
that its criticality parameter rises toward 1 as `μ→0`, with the locked fraction
tracking `f(c)`.

## 4. The honest reality check — structure yes, constants no

Does the circle map *hand over* `α` and `f`? No — and the data say why. For the
*critical sine-circle map* the harmonic `1/q` tongues scale as `q^{−3}`, i.e.
`α = 3`. The MMO value is `α ≈ 1.55`. So the MMO return map is **not** a universal
critical circle map; its `α` and `f` are set by the *specific* return map (the
L-spike excursion + reinjection) — the model-specific global flow. The circle-map
framework gives the *structure* (mode-locking, staircase, tongues) and the right
*variables*; the *constants* (`α=1.55`, `f≈0.0064/μ`) are the global return, which
remains model-specific / numerical. This is the same boundary the trilogy work
found: inner closes, global is organisable but not universal.

## 5. The payoff for #4 (σ_pq) — a real relocation

In the circle-map picture, `σ_pq` is **noise-induced unlocking of the `q`-tongue**
— a studied problem (stochastic mode-locking / noisy circle maps). The dissolution
threshold is when noise drives the rotation number out of the locked interval:

```
σ_pq:  noise spread δρ(σ) ~ tongue width Δ_q ~ q^{−α}   ⇒   σ_pq ~ σ_* q^{−α/γ},
```

the same `q`-law Tier-1/`MMO_AMIN_MULAW.md` gave — now with a *home*: the μ-law of
`σ_pq` is the μ-dependence of the **tongue widths and the noise-unlocking
boundary** of the return circle map. This is why the *local* `μ^{3/2}` was refuted
(it's not local escape; it's tongue geometry) and why the exponent is uncertain
(it depends on the map's criticality, sub- vs critical). The right reference class
for #4 is therefore **noisy mode-locking**, not folded-node escape.

## 6. Verdict

| Tier-2 idea | outcome on `f(c)` |
|---|---|
| Circle-map / Arnold tongues (5) | **the right framework** — unifies `f, s_obs, α(#6), σ_pq(#4)` as mode-locking |
| Funnel RG (7) | the renormalisation of (5); gives scaling structure, not the constant |
| Painlevé/RH (4) | **redirected to #2** (rigorous route to the linear PCF κ=2π²) |
| Nonlinear LZ (6), splitting (8) | **ruled out** for `f` (no anomalous/nonlinear inner; floor already had) |

**Net.** Tier-2 does for the *global* MMO side what Tier-1 did for the inner: it
identifies the correct object and reduces the cluster. The whole global side —
`a_min`/`f` (#3), `α` (#6), `σ_pq` (#4) — is **one circle map**: its rotation
number, its tongue widths, and its noise-unlocking. The `μ→0` limit is its
approach to criticality (where `a_min→e^{−π²}`). The framework is universal and
unifying; the constants (`α=1.55`, `f≈0.0064/μ`) are the model-specific return
map, confirming the inner/global boundary one more time.

**Concrete next step (still no new theory needed):** reduce the deterministic
return to its 1-D circle map at several `c`, measure its rotation number `ρ(c)`
and criticality, and read off `α` and `f(c)` together — closing #3 and #6
simultaneously as outputs of one measured map, and pointing #4 at the
noisy-mode-locking literature for its μ-law.
