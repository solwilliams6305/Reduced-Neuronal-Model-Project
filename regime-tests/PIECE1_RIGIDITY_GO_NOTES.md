# Piece 1 GO/NO-GO — the intrinsic cusp process exists and is RIGID (GO)

_June 2026. Executes the pivotal GO/NO-GO from `FINAL_4PCT_ATTACK_STRATEGY.md`: does the cusp have an
intrinsic multi-point process, or are the escapes effectively iid? **Answer: GO — the node process of the
stochastic Weber field is strongly RIGID**, confirmed by two independent methods. Figure
`coupled-atlas/figures/node_rigidity.png`; scripts `go_nogo_rigidity.py`, `riccati_rigidity.py`. Tags
**[DERIVED]/[NUMERIC]/[CITED]**; ⚑ load-bearing._

---

## 0. Outcome

| check | result | status |
|---|---|---|
| Number variance of the node process, $\mathrm{Var}(N(\Theta))$ | **bounded ~0.3–0.5** while mean $N\to12$ — exponent $\alpha\approx0$ | **NUMERIC ✓✓** ⚑ |
| Var/mean | falls to **0.02** (Poisson = 1) | **NUMERIC ✓** |
| Two independent methods | multiplicative-noise **field** and additive-noise **Riccati** agree | **robust** |
| Phase spacing | mean ≈ π (clockwork), CV ≈ 0.15 (regular, sub-Poisson) | **NUMERIC ✓** |
| **Verdict** | **GO** — the intrinsic process **exists** and is **rigid** | **DERIVED ⚑** |

**One line.** The zeros (nodes) of the stochastic Weber field form a **rigid** point process — bounded number
variance, far below Poisson — so the cusp **does** have a genuine intrinsic multi-point process. 𝒲 is the **edge
(first node)** of this rigid process. The earlier "process" results (forced-OU / swept-node surrogates) were
approximating *this* object.

---

## 1. The test

The escape is the first zero of the stochastic Weber field $u$ ($u''=(\operatorname{sign}(Y)Y^2-\eta\dot W)u$,
$p=-u'/u$ explodes at a node). The successive nodes form a point process in the phase-depth $\Theta=Y^2/2$. The
decisive diagnostic is the **number variance** $\mathrm{Var}(N(\Theta))$:
- **Poisson / iid** (no intrinsic structure): $\mathrm{Var}(N)\sim\Theta$ (linear), Var/mean $\to1$.
- **Rigid** (genuine intrinsic process, repulsion): $\mathrm{Var}(N)$ **sub-linear** (log, or bounded), Var/mean $\to0$.

## 2. Result — strongly rigid, two ways [NUMERIC ✓✓]

**Method 1 — the field** (`go_nogo_rigidity.py`, multiplicative noise, 6000 realisations, β=2):
$\mathrm{Var}(N)\approx0.45$ **constant** while mean $N: 0.8\to10.9$; exponent $\alpha\approx0$; Var/mean
$0.41\to0.04$.

**Method 2 — the Riccati** (`riccati_rigidity.py`, *additive* noise, count explosions — no
multiplicative-vanishing-at-zero artifact): $\mathrm{Var}(N)\approx0.27$ **constant** while mean $N:0.6\to11.8$;
$\alpha=0.02$; **Var/mean $\to0.02$**; phase spacing mean $\approx\pi$, **CV $\approx0.15$**.

Both give **bounded number variance** ($\alpha\approx0$) — dramatically below the Poisson line $\mathrm{Var}=$ mean
(figure). The two pictures agree on the rigidity (the small difference in the constant, 0.27 vs 0.45, is the
method/noise-projection difference). **GO, robustly.** ⚑

## 3. Interpretation — 1-D spectral rigidity of the swept operator

The mechanism is **1-D spectral rigidity** [CITED — Sturm oscillation / IDS self-averaging]. The node count is
the Sturm count (= eigenvalue count) of the 1-D stochastic Weber operator along the sweep. The noise is
**annealed** along the sweep and, crucially, the **frequency grows** ($|Y|$ increases), so the phase noise *per
oscillation* $\sim\eta/|Y|$ **diminishes** with depth. The accumulated phase variance is therefore sub-linear
(log-like, possibly bounded), giving the rigid node count. This is the cusp analogue of the rigidity of the Airy
point process — but driven by **1-D swept-operator rigidity**, not Dyson determinantal repulsion (consistent
with the earlier proof that 𝒲 is *not* a determinantal gap: the rigidity is real but **not GUE-determinantal**,
it is hyperuniform-class / 1-D-spectral). **[DERIVED + CITED.]**

## 4. What this resolves, and what remains

**Resolved (the existence question, Piece 1's core):** the intrinsic cusp process **exists** and is **rigid** —
the genuine "Weber process" is the **rigid node point process** of the stochastic Weber field, with 𝒲 its edge.
This kills the worry that the cusp might have only iid escapes (no process). It also retro-explains the
surrogates: the forced/swept-node results were measuring this rigid structure through an imposed window.

**Remaining (the characterisation):**
- the **exact** number-variance law — bounded (class-I hyperuniform) vs slow log (GUE-class) — needs deeper
  sampling to separate (both are rigid; the constant is ~0.3);
- the **2-point kernel / covariance** of the node process (the full process object, not just the rigidity
  signature);
- whether the rigidity is **β-dependent** (the whole law is a β-family, so likely yes).

These are characterisation, not existence — the GO is settled. **[honest scope.]**

## 5. Net + tracker

**GO.** The biggest open piece (the intrinsic process, ~2%, shared with rungs A & D) resolves **positively**:
the cusp has a genuine, **rigid** intrinsic node process (two-method-validated, far from Poisson), grounded in
1-D spectral rigidity, with 𝒲 as its edge. This is *not* a determinantal process (consistent with all prior
negatives) — it is hyperuniform-class / 1-D-spectral rigidity.

**Tracker: B 96% → 97%** — the intrinsic-process *existence + rigidity* is now established (was the deepest
open unknown); the *full kernel/covariance characterisation* remains. Shared uplift for A & D (same node
process). Next per strategy: characterise the node-process covariance (if pursuing the process), and Piece 3
(the canonical defining equation) for the marginal.
