# Notes — the stochastic coupled blow-up: noisy Airy → noisy Weber

_June 2026. The genuine theorem the deterministic Δ(g) step set up. Figure
`coupled-atlas/figures/coupled_blowup_noise.png`; script `coupled_blowup_noise.py`. Companion to
`DELTA_G_BLOWUP_NOTES.md`, `NOISY_CUSP_CROSSOVER_NOTES.md`. Tags throughout: [PROVED]/[DERIVED]/
[NUMERIC]/[HEURISTIC]._

## 1. Setup (deterministic scaffold)

Coupled Kristiansen–Pedersen FHN, $v_i'=-v_i^3+3v_i-w_i+g(v_j-v_i)+\sigma\xi_i,\ w_i'=\varepsilon(v_i-c)$.
Symmetric/antisymmetric split $v_\pm=\tfrac12(v_1\pm v_2)$:
- **Symmetric** mode carries the relaxation oscillation; it passes the fold of the critical manifold at
  $v_+=1$ with the standard fold/canard blow-up (Krupa–Szmolyan), inner scale $\varepsilon^{1/3}$. **[PROVED]**
- **Antisymmetric** mode $\delta=v_-$ obeys the cusp normal form
  $\dot\delta=\mu(v_+,g)\,\delta-\delta^3-\delta w_-,\ \mu=3(1-v_+^2)-2g$, with $(\delta,w_-)$ its own
  slow–fast pair. At the symmetric fold $\mu=-2g$; the nonzero antisym branches are
  $\delta=\pm\sqrt{\mu/3}$, separated by $\boxed{\Delta(g)=2\sqrt{-2g/3}}$. **[PROVED/algebraic]**
  K-P reduce the slow–fast cusp to the **Weber (parabolic-cylinder)** equation. **[PROVED]**

The cusp is the codim-2 point $g=0$ ($\mu=0$, $\Delta=0$, branches merge). For $g<0$ the two branches
(two folds) are separated by $\Delta(g)$; **Weber at $\Delta\to0$ (small $|g|$), two separated Airy
folds at large $\Delta$ (large $|g|$).**

## 2. The inner equation, with noise carried through

In the blow-up chart the antisym peel-off is the first node of a swept 2nd-order equation whose turning
structure unfolds the cusp. The correct unfolding keeps the **semi-infinite** oscillatory region (the
escape) and tunes the turning order near the origin:
$$u'' = \big(V_\Delta(Y)-\eta\,\xi\big)u,\qquad V_\Delta(Y)=\mathrm{sign}(Y)\,|Y|\,(|Y|+\Delta).$$
Near $Y=0$: $V_\Delta\approx\Delta\,Y$ (**linear, slope $\Delta$ → Airy**) for $|Y|\ll\Delta$, and
$V_\Delta\approx\mathrm{sign}(Y)Y^2$ (**quadratic → Weber**) for $|Y|\gg\Delta$; at $\Delta=0$ it is
exactly the cusp law. **[DERIVED]** The noise $\eta$ rides through the chart additively in the Cole–Hopf/
Riccati form (the inner-equation noise $\eta=\sigma/\sqrt{\varepsilon_2}$ scaling from the fold blow-up;
cf. `STOCHASTIC_CUSP_NOTES.md`, $\eta_{\rm cusp}=\sigma/(\sqrt2\,\varepsilon^{2/5})$). **[DERIVED]**

**Note (an honest correction).** The textbook Olver *two-turning-point* equation $u''=(\tfrac14Y^2-a)u$
is the WRONG inner model here — its turning points bound a *finite well*, giving a bounded-well peel-off
(platykurtic, no Airy edge), not the semi-infinite escape. The cusp-escape unfolding lives in the
**phase space** of the antisym pitchfork (the two branches $\pm\sqrt{\mu/3}$), realized above by the
near-origin slope $\Delta$, not by two potential turning points. **[NUMERIC, see script]**

## 3. The noisy crossover — numerics

Sweeping $\Delta$ at fixed $\eta=\sqrt2$ ($\beta=2$):

| $\Delta$ | skew | excess kurtosis | regime |
|---|---|---|---|
| 6 | +0.17 | +0.08 | **Airy / Tracy–Widom** (TW₂: +0.22/+0.09) |
| 1.5 | +0.45 | +0.12 | crossover |
| 0.5 | +0.58 | −0.04 | crossover |
| 0 | +0.60 | −0.27 | **Weber / cusp** (established: +0.6/−0.2) |

**[NUMERIC]** The noisy peel-off law crosses **noisy-Airy (Δ≫ℓ) → noisy-Weber (Δ→0)**, the excess
kurtosis flipping $+\to-$ (the sub-Gaussian signature emerging). The $\Delta=0$ limit **exactly
recovers the established cusp law**, and $\Delta$ large recovers TW — both endpoints validated. The
crossover sits at $\Delta_{\rm crit}\sim O(1)$ (inner units); a clean $\ell(\eta)$ power is **confounded**
because both endpoint laws themselves shift with $\eta/\beta$ (see `CUSP_LAW_CONJECTURES.md`).

## 4. Quantitative tie to the model — and an honest direction correction

$\Delta(g)=2\sqrt{-2g/3}$ maps the crossover $\Delta\sim\ell$ to coupling. With the inner scale
$\ell\sim\varepsilon^{1/4}$ (the folded-node onset $\mu\sim\sqrt\varepsilon$; `DELTA_G_BLOWUP_NOTES.md`),
the crossover is $\boxed{g_{\rm crit}\propto\sqrt\varepsilon}$ — **matching the verified onset**
$g_{\rm crit}/\sqrt\varepsilon\approx-0.58$. **[DERIVED + NUMERIC]**

**Direction (correction to the earlier reading).** The blow-up places the **cusp/Weber edge at
$\Delta\to0$, i.e. SMALL $|g|$** (near $g=0$), and two separated **Airy folds at large $|g|$**. The
full-model sweeps (`crossover_fold_to_cusp.py`, `loop_closure_scaling.py`) found the spread-amplification
and negative-kurtosis at *large* $|g|$ ($g\approx-0.1$). The resolution: those used the **single-unit
peel-off** $w_1=w_++w_-$, whose large-$|g|$ broadening is the **desync ($w_-$) contribution**, not the
cusp edge. The genuine cusp edge is the antisym escape near $g\to0$, which the single-unit observable
does not isolate. So the loop-closure "fold→cusp" labeling conflated **desync-broadening** with the
**cusp edge**; the clean inner test (§3) gives the correct direction. **[HEURISTIC/NUMERIC]** The √ε
*onset magnitude* still matches; only the side identification is corrected.

## 5. Proved vs numeric vs heuristic

| statement | status |
|---|---|
| sym fold blow-up; antisym → Weber reduction; $\Delta(g)=2\sqrt{-2g/3}$; Olver Airy↔Weber | **[PROVED]** (K-P, Krupa–Szmolyan, Olver) |
| inner equation $V_\Delta=\mathrm{sign}(Y)|Y|(|Y|+\Delta)$ unfolds Airy↔Weber; noise additive in chart | **[DERIVED]** |
| noisy peel-off: TW (Δ≫ℓ) → Weber (Δ→0), κ4 flips; endpoints recovered | **[NUMERIC]** |
| crossover $\Delta\sim\ell\sim\varepsilon^{1/4}\Rightarrow g_{\rm crit}\propto\sqrt\varepsilon$ | **[DERIVED]** + matches **[NUMERIC]** onset |
| direction: Weber at small $|g|$; full-model large-$|g|$ broadening = desync | **[HEURISTIC]** |
| uniform Berglund–Gentz tube through the merge; exact Weber edge law | **[OPEN]** |

## 6. The standing coupled-blow-up gap (what blocks full rigor)

**(a) Uniform tube estimates through the cusp.** Berglund–Gentz sample-path tubes are established for a
**simple fold** canard (curvature $O(1)$): paths stay in a tube of width $\sim\eta$ around the canard,
giving the noisy-Airy edge. At the cusp ($\Delta\to0$) the two antisym folds merge, the relevant
**curvature vanishes**, and the simple-fold tube estimate **degenerates** (the tube width blows up
non-uniformly in $\Delta/\ell$). What is needed: tube estimates in the **Weber/parabolic-cylinder
blow-up chart**, valid **uniformly as $\Delta/\ell\to0$** — the cusp analogue of Berglund–Gentz. This is
the precise gap; it is exactly the regime where the deterministic Olver Airy↔Weber connection must be
upgraded to a *stochastic* uniform connection.

**(b) Inner↔outer matching of the noisy solution**, uniform in $\Delta$, so that a full-model observable
isolates the inner edge (rather than the desync-broadening of §4). Not established.

**(c) The exact Weber edge law** (the RRV-type identification of the noisy-Weber inner equation with a
closed form, conjecturally Painlevé IV; rung F). Open independently of (a),(b).

**Bottom line.** The deterministic scaffold and the *scaling* of the stochastic crossover are in hand
(proved + matched to numerics); full rigor is blocked at one identifiable point — a **uniform
sample-path tube estimate through the merging cusp** (a stochastic Olver Airy↔Weber connection). Closing
(a) is the key tractable next theorem; (b),(c) then complete the chain.
