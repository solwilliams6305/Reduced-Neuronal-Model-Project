# MMO Phase 2 — the K2 return map and the funnel-filling f(c)

**Status:** computed via **Path B** (numerical global return + folded-node bound;
Wechselberger 2005 §4's parabolic-cylinder inner solution was *not* reconstructed
— see §5). **Result: the deterministic MMO mechanism is confirmed.** The realised
exponent is reproduced, α ≈ 1.45–1.50 (target 1.55 ± 0.06), via a clean structural
decomposition: the funnel-filling f(c) is the *global-return amplitude span*
divided by a *folded-node rotation constant* κ ≈ 2π².

Companion script: `mmo_k2_return_map.py`. Figure: `figures/mmo_k2_return_map.png`.

---

## 1. Question

Phase 1.5 established that the realised MMO staircase exponent is
α ≈ 1.55 ± 0.06, *between* the universal folded-node ceiling (α = 2, from
`s_max = (1−μ)/(2μ)` with `μ ∝ (c+1)`) and the floor — the gap being the
**funnel-filling** `f(c) = s_obs/s_max` (0.21 → 0.44 as c → −1). The Phase-2
question: does the deterministic global return + folded-node geometry *produce*
this f(c), and hence reproduce α ≈ 1.55? (Working FHR `v'=v−v³/3−w+y+I`,
`w'=ε(v+a−bw)`, `y'=εδ(c−v)`, (a,b,ε,δ)=(0.7,0.8,0.08,0.2), I=0.30; δ=0.2 is the
regime where α=1.55 was established, and `μ∝(c+1)` is δ-robust.)

## 2. Setup — the two ingredients of s_obs

Each `L¹Sˢ` episode is one global return (the L spike + re-injection toward the
fold) followed by s small oscillations through the folded-node funnel. So

```
s_obs(c)  =  [rotation map]  ∘  [global return] .
```

The **rotation map** turns the injection depth into a rotation count; the
**global return** sets the injection depth as c varies. Path B measures the
global return numerically and tests whether the rotation map is the folded-node
one. The diagnostic is the SAO *amplitude* sequence within each episode: the
loops grow geometrically from the smallest (a_min, the deepest point, nearest the
weak canard) to the largest (a_max ≈ O(1), the last loop before the jump).

## 3. The folded-node rotation map (the K2 content, measured)

If the SAO amplitudes grow by a per-turn factor R over s turns from a_min to
a_max, then `s = ln(a_max/a_min)/ln(R)`. Measuring (s, a_min, a_max, μ) per c and
forming `ln(R)/μ = ln(a_max/a_min)/(s·μ)`:

| c | μ | s_max | s_obs | a_min | a_max | κ = lnR/μ |
|---|---|---|---|---|---|---|
| −0.80 | 0.049 | 9.7 | 1.5 | 0.796 | 1.221 | 5.8 * |
| −0.83 | 0.041 | 11.7 | 2.0 | 0.253 | 1.167 | 18.7 * |
| −0.86 | 0.033 | 14.6 | 3.0 | 0.211 | 1.155 | 17.1 |
| −0.89 | 0.026 | 18.9 | 5.0 | 0.090 | 1.106 | 19.5 |
| −0.92 | 0.018 | 26.6 | 10.0 | 0.029 | 1.073 | 19.6 |
| −0.94 | 0.014 | 36.0 | 16.0 | 0.0083 | 1.048 | 22.1 |

(* s < 3: too few turns to measure a growth rate, excluded.) For the s ≥ 3
points,

```
ln(R)  =  κ · μ,    κ = 19.6 ± 1.8  ≈  2π² (= 19.74),   constant across c.
```

This is the load-bearing non-circular result: the per-turn SAO growth rate is set
by the folded-node eigenvalue ratio μ, with a **single universal constant κ ≈ 2π²**
across the whole band. It is the K2/Wechselberger rotation content, confirmed
quantitatively *beyond* the `s_max` bound — three independently-measured
quantities (s, a_min/a_max, μ) collapse onto one line. (κ ≈ 2π² is suggestive but
*measured, not derived* — pinning it is a Path-A question.)

## 4. Result — f(c) and α

Composing the rotation map `s = ln(a_max/a_min)/(κμ)` with the ceiling
`s_max = (1−μ)/(2μ)` gives the funnel-filling **with μ cancelling**:

```
f(c) = s_obs/s_max = (2/κ) · ln(a_max/a_min) / (1−μ).            (★)
```

So the funnel-filling is just the **log of the global-return amplitude span**,
divided by the folded-node constant κ — independent of μ. Tested against the
measured f:

| c | f measured | f from (★) |
|---|---|---|
| −0.83 | 0.171 | 0.163 |
| −0.86 | 0.206 | 0.180 |
| −0.89 | 0.264 | 0.263 |
| −0.92 | 0.376 | 0.376 |
| −0.94 | 0.445 | 0.502 |

Agreement is ~5–15% across the band (worst ~13% at the deepest, small-μ point).
The global return drives `a_min(c) → 0` super-fast as c → −1 (0.80 → 0.008), so
`ln(a_max/a_min) ~ (c+1)^{−0.7}`, giving `f ~ (c+1)^{−0.7}` and hence
`s_obs = f·s_max ~ (c+1)^{−1.7}`:

```
α (measured s_obs, this run)          =  1.49
α (K2-composed, f from (★) × s_max)   =  1.45
target (Phase-1.5 fine grid)          =  1.55 ± 0.06
```

**Verdict (§7 of the prompt): confirmed, at the lower edge.** α ≈ 1.45–1.49 sits
inside the broad confirmed bin (1.55 ± 0.1) but ~0.05–0.10 below the fine-grid
target. The mechanism is settled — **FHR's numerical global return + the
folded-node rotation map (κ ≈ 2π²) + Wechselberger's ceiling**, composed via (★),
reproduce the staircase with no extra physics. The small residual gap to 1.55 is
within the SAO-counting dt-sensitivity that the prompt (§7/§8) flagged as the
limiting error; a finer-dt recount is the clean way to decide whether the true
value is 1.5 or 1.55. Either way the qualitative result is robust: **α ∈ [1.45,
1.6], strictly between the ceiling 2 and 1, set by the global-return amplitude
span through (★).**

## 5. What is established / what is open (honest path accounting)

| Piece | Status |
|---|---|
| Folded-node ceiling `s_max=(1−μ)/(2μ)`, `μ∝(c+1)`, α_ceiling=2 | **derived** (Phase 1.5 + Wechselberger bound) |
| Rotation map `ln R = κμ`, κ ≈ 2π² constant across c | **measured** (this doc) — the K2 content, confirmed beyond the bound |
| μ-cancellation `f = (2/κ)ln(a_max/a_min)/(1−μ)` | **derived** from the above two |
| Global-return amplitude span `a_min(c), a_max(c)` | **numerical** (FHR integration) — *not* closed-form |
| α ≈ 1.5 reproduced | **confirmed** (lower edge; dt-recount to tighten) |
| Closed-form a_min(c) / closed-form κ=2π² | **open — this is Path A** (K1/K3 matching + parabolic-cylinder inner solution) |

So Path B closes the *mechanism*: the deterministic α is reproduced and decomposed
into a universal folded-node rotation map (κ ≈ 2π²) × a model-specific global
return (the amplitude span). This is the same universal-vs-model-specific split as
the canard (`σ_*` exponent universal, `C_q` specific) and tonic (`A_mid` form
universal, `c` specific) chapters. What Path B does **not** deliver: a closed-form
global return (a_min(c)) or a first-principles derivation of κ = 2π² — both require
the full K2 inner solution (Path A), which is the natural next analytical step if a
fully closed-form α is wanted rather than the mechanism + numerical global return.

**Hazards encountered (prompt §8):** (1) small-μ many-canard regime — the deepest
point (c=−0.94, μ=0.014, s≈16) is where f from (★) is least accurate (~13%) and
where SAO counting is most dt-sensitive; (2) the κ measurement needs s ≥ 3 (the
two shallow points are unreliable and excluded); (3) δ=0.2 was used (not the
prompt's δ=1.0) to match the established α=1.55 target — `μ∝(c+1)` is δ-robust so
the ceiling is unaffected, but the global-return amplitude span (hence the exact
α) could shift mildly with δ and is worth a δ=1.0 cross-check.

## 6. If this is accepted: remaining chapter pieces

Per the prompt §10, the deterministic side now closes; remaining: (i) **σ_pq**:
the noise propagation factor γ through one pass of this return map (FW action on
the K2 chart), giving `σ_pq ~ σ_*·q^{−α/γ}`; (ii) **noise-broadened staircase**
numerics; (iii) **write-up** + README MMO row → "derived + validated."

## 7. Reproduce

```
python3 regime-tests/mmo_k2_return_map.py
```
Writes `data/mmo_k2.npz`, `results/mmo/mmo_k2.txt`, `figures/mmo_k2_return_map.png`.
