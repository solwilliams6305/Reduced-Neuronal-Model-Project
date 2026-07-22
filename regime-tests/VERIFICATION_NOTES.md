# Verification notes for CANARD_BLOWUP.md

Independent re-derivations of the load-bearing identities, plus one bug found
during the sanity checks.

## 1. FW action scaling by dimensional analysis

The rate function has dimensions [velocity² · time / σ²]. Under v → ε^{1/3} V
and t → ε^{−1/3} T:

- (γ̇_v − f)² → ε^{4/3} (dV/dT − V² + W)²
- dt → ε^{−1/3} dT
- combined factor: ε^{4/3} · ε^{−1/3} / σ² = ε / σ²  ✓

So I[γ] = (ε / 2σ²) ∫ (dV/dT − V² + W)² dT, matching §4 of CANARD_BLOWUP.md.

## 2. σ_eff direct check via variance

Var(δv) after time t with noise σ on dv: Var = σ² t.
In blow-up: V = (v + 1)/ε^{1/3}, so Var(δV) = σ² t / ε^{2/3}.
Substitute t = ε^{−1/3} T: Var(δV) per unit T = σ² / ε ⇒ σ_eff² = σ²/ε ⇒
σ_eff = σ/√ε.  ✓

## 3. Singular-limit / Riccati reduction at λ = 0

dV/dT = V² − W with W conserved. Fixed points V_* = ±√W; V_* = +√W repelling,
V_* = −√W attracting. Trajectories from V < −√W go to V = −∞; from V > +√W go
to V = +∞; from |V| < √W asymptote to ±√W. The maximal-canard is the full
parabola V² = W.  ✓

## 4. Jump-ansatz check for A_min

ΔV = 2√W_off, jump time τ = (admissible ΔW)/λ = √W_off / λ (admissible ΔW set
by parabola curvature at W_off). A_jump ≈ (ΔV)² / τ = 4 W_off · λ / √W_off
= 4 λ √W_off.  ✓ matches `canard_action_min.py`.

Optimum: W_off* set by the explosion threshold W_off ~ λ^{2/3} below which
"maximal canard" is no longer the right object (5b regime takes over).
A_min ~ 4 λ · λ^{1/3} = 4 λ^{4/3}.  ✓

## 5. Boundary continuity between 5a and 5b

At the W_off → λ^{2/3} boundary, σ_crit(5a) = √(ε · A_min / 2) = √(2 ε λ^{4/3}).
At the same boundary, the explosion-window result σ_crit(5b) ~ ε^{3/4}
corresponds to plugging λ ~ ε^{1/2} into 5a:

```
σ_crit(5a) | λ=ε^{1/2}  =  √(2 ε · ε^{2/3})  =  ε^{5/6}.
```

5b gives ε^{3/4}. The mismatch (5/6 vs 3/4) is the genuine signature that 5b
is *not* the boundary continuation of 5a — it lives in a different scaling
regime where the leading normal form (dV/dT = V² − W, dW/dT = −λ) is
degenerate and a higher-order correction (the cubic V³/3 term, or the next
slow-drift term) becomes the leading singular behaviour. BGK 2012 carries
this through carefully; we do not reproduce it here, but flag that the two
regimes are joined by a non-trivial matching layer and the 5a derivation
should be applied for λ ≫ ε^{1/2} only.

## 6. Bug found in kernel.py — `w_fold_left` returns the right fold

`FHN2D.w_fold_left` returns `self.I + 2.0/3.0`. From the v-nullcline w = v − v³/3 + I:

- v = −1 ⇒ w = −2/3 + I  (left fold, *true* value)
- v = +1 ⇒ w = +2/3 + I  (right fold, what the code currently returns)

Numerical check at three project-relevant I values:

| I      | w_fold_left (true) | w_fold_left (code) | w_fold_right (true) |
|--------|--------------------|--------------------|---------------------|
| −0.100 | −0.7667            | +0.5667            | +0.5667             |
| +0.292 | −0.3747            | +0.9587            | +0.9587             |
| +0.331 | −0.3357            | +0.9977            | +0.9977             |

So `w_fold_left` is silently swapped with `w_fold_right`. The README's
descriptive line — "Left fold: w_fold = I + 2/3 ≈ −0.767" for I = −0.1 — gives
the *numerical value* of the true left fold but the wrong *formula* for it.

`reduced_prediction` uses `model.w_fold_left` as the upper limit of the slow
integral in first-passage mode (kernel.py L490). With the swap, that integral
runs from W_FP ≈ −0.70 *upward* past the FP to the right-fold value ≈ +0.57,
passing through a region where the Newton solve for v_s(w) almost certainly
fails to track the left branch. The empirical fold scaling ε^{1/2} you have
validated must therefore be coming out *despite* this — likely because the
spurious integrand collapses to a regular trapezoidal evaluation that
happens to scale ~1/ε for unrelated reasons. Worth a recheck before pushing
the resonator or canard reduced-prediction calls.

I have *not* applied this fix — leaving it for you to decide, since
"validated" numerics may depend on the current (buggy) behaviour.

Suggested patch:

```python
@property
def w_fold_left(self) -> float:
    """w-value of left fold (v = -1):  w = I - 2/3."""
    return self.I - 2.0 / 3.0

@property
def w_fold_right(self) -> float:
    """w-value of right fold (v = +1):  w = I + 2/3."""
    return self.I + 2.0 / 3.0
```
