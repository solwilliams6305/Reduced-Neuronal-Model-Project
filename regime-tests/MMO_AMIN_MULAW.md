# `a_min` (#3) and the `σ_pq` μ-law (#4) — the regime, the coupling, and the boundary

**What this is.** An honest attempt at the two coupled MMO unknowns: the smallest
SAO amplitude `a_min(c)` (#3) and the μ/c-dependence of the noise dissolution
threshold `σ_pq` (#4). Both were flagged "needs Krupa–Wechselberger 2010
transition-regime inner solution." The attempt does **not** close either number —
they are the genuine special-function zone — but it pins three things worth
having: the regime diagnosis (why the standard tools fail here), the dependency
structure (a_min is the keystone of the whole MMO noise/staircase side), and the
derivable part of `σ_pq`. Verdict up front:

> `a_min` is **exponentially small and global-return-set in the transition regime
> `μ < ε^{1/2}`** — not a clean inner constant. `σ_pq`'s *q*-dependence is
> derivable (`σ_* q^{−α/γ}`); its *μ*-law reduces to `a_min` and its exponent is
> genuinely uncertain (3/2 vs 1). Both are correctly PhD/paper-scope, gated on one
> object — the KW2010 transition inner solution — which I will not reconstruct
> from memory.

---

## 1. Regime diagnosis (the clean, certifiable fact that frames both)

The standard folded-node "sector scaling" `a_min ~ ε^{(1−μ)/2}` is derived for
`μ ≫ ε^{1/2}` (the genuine folded *node*). Check the data: `ε = 0.08 ⇒
ε^{1/2} = 0.283`, and the MMO_K2 band is `μ ∈ [0.014, 0.049]`. **Every point has
`μ < ε^{1/2}`.** The whole accessible band is in the **transition regime**
(folded node → folded saddle-node), where the standard sector scaling is *out of
its domain of validity*. That is exactly why the audit calls it "an upper bound
only":

```
sector bound at μ=0.03:  ε^{(1−μ)/2} = 0.08^{0.485} ≈ 0.29 ;
measured a_min (c=−0.89, μ=0.026):  0.090  ≪  0.29 .
```

The measured `a_min` is ~3× below the bound and the gap *widens* as `μ→0`. So the
data live below the standard theory, in KW2010's transition regime. This single
fact explains why #3 and #4 resisted the standard machinery.

## 2. `a_min` is exponentially small and set by the global return (#3)

From the MMO_K2 measurements, `a_min` collapses far faster than any power:

```
μ:      0.049  0.041  0.033  0.026  0.018  0.014
a_min:  0.796  0.253  0.211  0.090  0.029  0.0083
```

`a_max ~ O(1)`, so `−ln a_min ≈ ln(a_max/a_min) ~ (c+1)^{−0.7} ~ μ^{−0.7}`, i.e.

```
a_min(c) ~ exp( −C · μ^{−0.7} )      (exponentially small as μ→0).
```

Two readings, both pointing away from a clean closed form:

- **Transition-regime closeness.** The *standard* folded node has
  exponentially-close secondary canards `~ exp(−c/μ)` (exponent `1/μ`). The
  measured exponent is *softer* — `μ^{−0.7}`, not `μ^{−1}`. That softening is the
  signature of the transition regime `μ < ε^{1/2}`: KW2010 replaces the `1/μ`
  exponent with a transition exponent. Deriving `0.7` (and the constant `C`) is
  exactly the KW2010 inner solution — special functions interpolating Weber
  (folded node) and Airy (canard).
- **It's global.** `a_min` is the *injection depth* — where the deterministic
  global return drops the trajectory into the funnel — which sits far below the
  inner-solution bound. So even with KW2010's inner solution in hand, the realised
  `a_min` needs the model-specific global return. This is why MMO_K2 measured it
  numerically; a closed form is genuinely two hard pieces (transition inner ×
  global return), not one constant.

**Connection to what we've built.** The `μ→0` limit of this problem *is* the
canard, and we have it in closed form: `MMO_CROSSOVER.md` shows the folded-node
normal form → the canard Krupa–Szmolyan form as `μ→0`, and `TONIC_CMID_AIRY.md`
solves that fold via the **Airy** function (`V=−u'/u`). The finite-μ transition
regime is the **Weber→Airy interpolation** — KW2010's object. So #3 is the
finite-μ generalisation of the Airy result we already have; the open piece is
precisely the interpolation, where parabolic-cylinder connection coefficients
enter (the same zone as κ Route A — and the same reason I won't fake it).

## 3. `σ_pq`: what's derivable, what's gated on `a_min` (#4)

**The q-dependence (derivable, structural).** A `(1:q)` plateau dissolves when the
noise-induced rotation-number spread `δρ ~ σ^γ` reaches the plateau width
`Δ_pq ~ q^{−α}`:

```
σ_pq^γ ~ q^{−α}   ⇒   σ_pq ~ σ_* · q^{−α/γ} .
```

This is clean *given* `α` (≈1.55, itself `#2 ∘ #3`) and `γ` (#5). So the *shape*
of `σ_pq` in `q` is settled structurally; the open inputs are `α` (via `a_min`)
and `γ`.

**The μ-law (open, and not what the heuristic said).** The audit's claimed
`μ^{3/2}` for the c/μ-dependence came from a *local* sector-crossing argument that
`MMO_CROSSOVER.md` **refuted**: the local folded-node escape is canard-like (set
by the constant drift, prefactor `C_q`), with **no `μ^{3/2}`** — that factor is a
**global funnel-filling** effect, not local escape. So the μ-law is downstream of
the funnel structure, i.e. of `a_min`. The literature "sector-packing" route
`μ·ε^{(1−μ)/2}` needs KW2010, and as `μ→0` it scales like `μ·ε^{1/2}` (i.e.
`μ^1`, not `μ^{3/2}`) — consistent with the audit's warning that "the power may
not even be 3/2 (M2 degeneracy)." **The μ-exponent is genuinely unresolved**, and
resolving it is the transition-regime computation.

**Why `γ` is also gated on `a_min`.** The noise-rotation spread `δρ` comes from
phase diffusion as the trajectory winds the funnel; the angular noise is `~ η/r`
and is largest at the *deepest* loop `r ~ a_min`. So `δρ ~ η/a_min ·
(transit)^{1/2}` — `γ` and its prefactor inherit `a_min`'s scaling. Every noise
quantity on the MMO side funnels through `a_min`.

## 4. Synthesis: `a_min` is the keystone of the whole MMO noise/staircase side

```
                    a_min(c)   [#3, transition-regime inner × global return]
                   /    |     \
       α (#6) ←───┘     │      └───→ γ (#5)  [δρ ~ η/a_min]
   [#2∘#3, staircase]   │
                        ▼
              σ_pq μ-law (#4)  =  f( a_min, κ, γ )   — q-dep ✓ (σ_* q^{−α/γ}),
                                                       μ-dep open (1 vs 3/2)
```

`a_min` is the single shared unknown behind the deterministic staircase exponent
`α` (#6), the noise-rotation exponent `γ` (#5), and the dissolution μ-law (#4).
One transition-regime object gates the cluster. This is the MMO analogue of the
boundary `RESONATOR_INSTANTON_ATTEMPT.md` found for #10: a whole sub-chapter
reduces to one genuinely-hard special-function/global object.

## 5. Honest status

| Piece | Status |
|---|---|
| Regime diagnosis: all data in transition regime `μ<ε^{1/2}` | **derived (clean)** |
| Why sector bound `ε^{(1−μ)/2}` fails (it's the `μ≫ε^{1/2}` result) | **derived** |
| `a_min ~ exp(−C μ^{−0.7})`, softer than standard `exp(−c/μ)` | **measured + interpreted** (transition closeness) |
| `a_min` set by global return (≪ inner bound) ⇒ global, not a constant | **argued** |
| `σ_pq ~ σ_* q^{−α/γ}` (q-dependence) | **derived (structural)** |
| local `μ^{3/2}` refuted; μ-law is global, exponent 1 vs 3/2 | **(from MMO_CROSSOVER)**, open |
| `a_min` keystone for #4, #5, #6 | **synthesis (clean)** |
| **closed-form `a_min(c)`; the μ-exponent** | **open — KW2010 transition inner solution × global return** |

**Verdict.** #3 and #4 do not close — and, like #10, we now know precisely why and
where the boundary is. `a_min` is a transition-regime (`μ<ε^{1/2}`) object: the
inner part is the KW2010 Weber→Airy interpolation (parabolic-cylinder connection
coefficients — the fabrication zone I keep flagging), and the realised value is
further set by the global return, making it model-specific. `σ_pq` factorises into
a clean derivable q-law (`σ_* q^{−α/γ}`) and a μ-law that reduces to `a_min` with a
genuinely uncertain exponent. The high-leverage move, if this cluster is to be
closed, is therefore **a_min alone** — it unlocks #4, #5, #6 together — and it is
the one place where the KW2010 transition-regime inner solution (worked from the
text, not memory) is unavoidable.
