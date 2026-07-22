# First-principles `A_mid` (and `c`) from the inner adjoint-Floquet BVP

**Status:** the **mechanism** for `A_mid` is derived — the matched-asymptotic
adjoint-Floquet BVP at the Krupa–Szmolyan fold blow-up — and the inner adjoint
equation is **corrected**. A *precise* constant `c` is **not** delivered at
leading order: the BVP confirms the iPRC peak/shape from deterministic geometry
to ~15 % near ε = 0.04, but its integral overshoots and drifts with ε because
the `O(ε^{1/3})` inner correction is large at accessible ε. This addresses
TONIC_PHASE.md §21.1-item-1 honestly: the publishable nugget is the **corrected
inner equation + the BVP mechanism**, not yet a 2-figure `c`.

Headline (what is and is not first-principles):

- **Solid (standalone result):** the leading inner adjoint equation is
  `dZ̃_v/dT = −2V Z̃_v` (peak at the fold tip V ≈ 0), **correcting** the
  `(1/b − 2V)` of earlier drafts (§2). Confirmed numerically (peak location,
  blow-up scaling `Z_v ~ ε^{−2/3}`, width `~ε^{−1/3}`).
- **Confirmed mechanism:** solving the BVP (backward-shoot from the matched
  deterministic outer iPRC) reproduces the iPRC peak to ~15 % near ε = 0.04 and
  the correct profile shape, from cycle geometry alone — the fold blow-up *is*
  what sets `A_mid` (§3–§4).
- **NOT delivered:** a precise `c`. The leading-order reconstruction overshoots
  the measured `A_mid ≈ 0.40` and drifts strongly (`c_recon ≈ 1.9 → 4.5` across
  ε = 0.04–0.16) — the dropped `O(ε^{1/3})` correction (≈0.43 at ε = 0.08) plus a
  bulk shoulder that is not yet first-principles (§5). So the measured `c ≈ 1.55`
  has its **mechanism derived and inner equation corrected**, but a precise
  value awaits the coupled `O(ε^{1/3})` BVP + matched bulk integral.

Companion script: `tonic_cmid_bvp.py`. Companion figure: `figures/tonic_cmid_bvp.png`.

---

## 1. What "first-principles `c`" means

TONIC_PHASE.md §19 wrote `A_mid = √(c/π²)` and quoted `c ≈ 1.5`. But that `c`
was just `π² A_meas²` — algebraically the measured amplitude in disguise, not an
independent prediction. The genuine task (and the publishable one) is to
**predict `A_mid` from the deterministic limit cycle alone**. Recall the
validated relation (tonic_phase_response.py `compute_A`, confirmed against ISI
data to 5 % at σ ≤ 0.005):

```
CV = σ · A_mid,     A_mid² = (1/T_cycle²) ∮ Z_v(t)² dt ,                 (1)
```

with `Z_v` the v-component of the adjoint Floquet iPRC, normalised `Z·γ̇ = 1`.
The numerics show `∮Z_v²dt` is dominated by the two fold passages, where `Z_v`
peaks. So `A_mid` is set by the iPRC's behaviour in the fold blow-up — which is
what we solve for below.

## 2. The inner adjoint equation — derivation and a correction

The FHN adjoint (iPRC) equations are

```
dZ_v/dt = −(1 − v²) Z_v − ε Z_w,
dZ_w/dt =  Z_v + ε b Z_w.
```

Apply the Krupa–Szmolyan blow-up at the left fold (`v = −1`):

```
v = −1 + ε^{1/3} V,   w = w_f + ε^{2/3} W,   t = ε^{−1/3} T,   (1 − v²) = 2ε^{1/3}V − ε^{2/3}V².
```

The iPRC normalisation `Z·γ̇ = 1` with `γ̇ = (ε^{2/3}(V²−W), εg)` forces the
inner scalings `Z_v = ε^{−2/3} Z̃_v` and `Z_w = ε^{−1} Z̃_w` with `Z̃_v, Z̃_w = O(1)`.
Substituting and collecting orders:

```
dZ̃_v/dT = −2V Z̃_v + ε^{1/3}( V² Z̃_v − Z̃_w ),                          (2a)
dZ̃_w/dT =  Z̃_v + ε^{2/3} b Z̃_w.                                        (2b)
```

At **leading order**:

```
┌─────────────────────────────────────────────┐
│   dZ̃_v/dT = −2 V(T) · Z̃_v ,    Z̃_w = ∫ Z̃_v dT  │   (★)
└─────────────────────────────────────────────┘
```

**Correction to TONIC_PHASE.md §14.3 / §15.2.** Those sections wrote
`dZ̃_v/dT = (1/b − 2V) Z̃_v`. The spurious `1/b` came from *algebraically slaving*
`Z_w = −Z_v/(εb)` (imposing `dZ_w/dt ≈ 0`). That is wrong: eq. (2b) has **no
restoring term at leading order**, so `Z̃_w` is the *antiderivative* of `Z̃_v`,
not an algebraic slave, and it enters `Z̃_v` only at `O(ε^{1/3})`. The leading
inner adjoint is the clean homogeneous Riccati-coefficient equation (★).

**Numerical confirmation.** (★) puts the iPRC extremum where the coefficient
`−2V` vanishes, i.e. at `V ≈ 0` (the fold tip). The full-pipeline iPRC peaks at
`V = 0.10 / 0.12 / 0.16` for ε = 0.04/0.08/0.16 — at the tip, drifting slightly
positive with ε. The erroneous `(1/b − 2V)` would put the peak at `V = 1/(2b) =
0.63`, which the data clearly exclude. The blow-up scaling itself is confirmed:
`peak |Z_v|·ε^{2/3} ≈ 1.07–1.14` (flat ⇒ `Z̃_v = O(1)`), fold width in t `~ ε^{−1/3}`.

## 3. The boundary-value problem and how it is solved

(★) is **forward-unstable** on the incoming attracting branch (`V < 0 ⇒ −2V > 0`),
which is why TONIC_PHASE.md §17.1's forward IVP exploded. It is a genuine
two-sided BVP. The boundary data come from matching to the deterministic *outer*
iPRC on the slow manifold,

```
Z_v^{out}(v) = −1 / [ g (1 − v²) ]            (slow-branch quasi-steady iPRC),
```

which is fixed entirely by deterministic geometry (`g = v_f + a − b w_f`; no
adjoint normalisation, no CV data). The numerically clean matching window is
`V ≈ 1` (`v ≈ −0.6`), where the pipeline iPRC agrees with `Z_v^{out}` to ~1 %
(e.g. ε=0.08, V=1.0: measured `Z̃_v = 0.58` vs outer `0.586`). It is **not**
matched at `V ≳ 2` (`v ≈ 0`, mid fast-jump) where `Z_v^{out}` is invalid.

The solver (`tonic_cmid_bvp.py`):

1. integrates the deterministic inner passage `dV/dT = V²−W, dW/dT = −λ` from the
   attracting branch (`V ≈ −√W`);
2. sets `Z̃_v(V_match=1) = ε^{2/3} Z_v^{out}`;
3. integrates (★) **backward** from `V_match` through the tip into the incoming
   branch (stable: captures the `V≈0` peak, decays to ≈0 on the attracting
   branch), and **forward** from `V_match` (decaying outgoing tail);
4. forms the inner integral `R̃ = ∫ Z̃_v² dT`.

This reproduces the full-pipeline inner profile `Z̃_v(V)` (figure, left panel).

## 4. Assembling `A_mid` and the constant `c`

Converting the inner integral back to physical time (`Z_v = ε^{−2/3}Z̃_v`,
`dt = ε^{−1/3}dT`):

```
∮_fold Z_v² dt = ε^{−5/3} R̃   (per fold).
```

At mid-tonic the two folds are symmetric (`I_mid ≈ ½(I_fold_L + I_fold_R)` ⇒
`λ_L = λ_R = 0.467`, `|g_L| = |g_R| = 0.467`), so `R̃_L = R̃_R`. Writing the cycle
integral as fold + bulk:

```
A_mid² = ( ε^{−5/3}(R̃_L + R̃_R)  +  ∮_bulk Z_v²dt ) / T_cycle² ,    c := π² A_mid².
                └── BVP, no inversion ──┘   └── slow-branch shoulder ──┘
```

**Validation 1 — the inner profile (what the BVP rigorously delivers).** The
clean, range-independent comparison is the BVP-solved `Z̃_v(V)` against the
full-pipeline iPRC, pointwise (figure, left panel). The peak heights:

| ε | peak Z̃_v (BVP) | peak Z̃_v (pipeline) | ratio |
|---|---|---|---|
| 0.04 | 0.935 | 1.072 | 0.87 |
| 0.08 | 1.245 | 1.091 | 1.14 |
| 0.16 | 1.689 | 1.139 | 1.48 |

The BVP reproduces the peak to **~15 % near ε = 0.04** and the correct shape
(rise on the incoming branch, peak at the tip, outgoing decay). The overshoot
**grows with ε** (ratio → 1.48 at ε = 0.16) — the signature of the dropped
`O(ε^{1/3})` inner correction (`ε^{1/3} ≈ 0.43` at ε = 0.08).

**Validation 2 — A_mid reconstruction (what the leading order does NOT deliver).**
Assembling `A_mid² = [ε^{−5/3}(R̃_L+R̃_R) + bulk]/T²` (BVP fold integral + the
measured slow-branch bulk):

| ε | R̃/fold | A_recon | A_meas | c_recon = π²A² |
|---|---|---|---|---|
| 0.04 | 1.380 | 0.437 | 0.407 | 1.89 |
| 0.08 | 2.446 | 0.552 | 0.397 | 3.01 |
| 0.16 | 4.503 | 0.675 | 0.380 | 4.49 |

The leading-order reconstruction **overshoots** the (flat) measured `A_mid ≈ 0.40`
and **drifts strongly** (`c_recon ≈ 1.9 → 4.5`). So the leading BVP does **not**
pin a precise constant `c` at accessible ε. Two reasons: (i) the `O(ε^{1/3})`
inner correction is large there (it lowers the true peak/integral); (ii) the
bulk shoulder is not yet first-principles — the naïve outer form `−1/[g(1−v²)]`
overcounts the bulk by ~600× (it diverges where the slow nullcline is crossed,
`g → 0`), so the true slow-branch iPRC (exponentially small on the incoming
branch) cannot simply be read off it.

**Bottom line.** The BVP **derives the mechanism** and **corrects the inner
equation** — it shows the fold blow-up of the adjoint is what sets `A_mid`, and
reproduces the iPRC shape/peak from deterministic geometry. The measured
`c ≈ 1.55` is therefore mechanistically explained, but a *precise* first-principles
value is not yet in hand: it needs the coupled `O(ε^{1/3})` inner BVP (§5-i) and
a correct matched outer/bulk integral (§5-ii).

## 5. Two remaining analytical steps

**(i) The O(ε^{1/3}) inner correction (dominant at accessible ε).** The
leading-order peak and reconstructed `A_mid` overshoot and drift upward with ε
(peak ratio 0.87 → 1.14 → 1.48; `c_recon ≈ 1.9 → 4.5`). This is quantified by
(2a): the dropped `O(ε^{1/3})(V²Z̃_v − Z̃_w)` term is exactly the sub-leading
correction TONIC_PHASE.md §14.5 flagged, and at accessible ε it is **not small**
(`ε^{1/3} ≈ 0.34–0.54`). It lowers the true peak/integral and vanishes as ε → 0,
so the leading BVP is a genuine but only ~15 %-accurate (near ε = 0.04)
approximation. Coupling (2a)–(2b) as a 2-component inner BVP would remove the
drift.

**(ii) The slow-branch bulk shoulder (the open piece for the *full* `c`).** The
fold blow-up gives ~64 % of `∮Z_v²dt`; the remaining ~36 % is the slow-branch
iPRC, currently taken from numerics. It is *not* yet first-principles: the naïve
outer form `−1/[g(1−v²)]` overcounts it by ~600× (it diverges where the slow
nullcline is crossed, `g → 0`), and the true slow-branch iPRC is exponentially
small on the **incoming** branch (homogeneous-mode decay) and only follows the
outer form on the **outgoing** side. A correct matched outer integral must
respect that incoming/outgoing asymmetry. Together with (i), this is what stands
between the **derived mechanism** and a **fully-derived numerical `c`**; the
measured `c ≈ 1.55` is the target both corrections should converge to as ε → 0.

## 6. Universality vs model-specificity

The structure is universal for any 2D relaxation oscillator with a generic
fold: the inner adjoint is always `dZ̃_v/dT = −2V Z̃_v` (the fold normal form is
universal), and `A_mid = √(c/π²)` with

```
c = π² · (∑_folds ε^{−5/3} R̃_i) / T_cycle²     (model-specific O(1) constant),
```

where `R̃_i` depends on the model's fold drift `λ_i` and slow-drift `g_i`. For
FHN(0.7, 0.8) at mid-tonic, `c ≈ 1.55`. Van der Pol, Morris–Lecar, etc. give
the same functional form with their own `c` — the same universality/specificity
split the canard chapter established for `σ_* = C_q √ε λ^{1/2}` (universal
exponent, model-specific `C_q`). This is tested directly in `VDP_CROSSMODEL.md`.

## 7. Reproduce

```
python3 regime-tests/tonic_cmid_bvp.py
```

Solves the inner adjoint BVP at each ε, validates `Z̃_v(V)` against the pipeline,
assembles `A_mid`, reports `c`. Writes `results/tonic_cmid_bvp/summary.txt` and
`figures/tonic_cmid_bvp.png`.
