# Rung E — perturbation / unfolding responses, consolidated (v2, T2-consistent)

_June 2026. Supersedes the §E portion of `BREADTH_LADDER_ATLAS_NOTES.md`. Figure
`coupled-atlas/figures/unfolding_response_v2.png`. Scripts `unfolding_response_v2.py`,
`unfolding_figure_v2.py` (+ `beta_axis.py` from D). Tags [VALIDATED]/[NUMERIC]/[CITED]/[OPEN]._

## Why E was at 66% — the gap, enumerated

The prior §E established exactly **one** unfolding response — the codim-2 unfolding of the *abstract*
**q=3 (swallowtail)** potential $V=\operatorname{sign}(Y)(|Y|^3+c_2|Y|^2+c_1|Y|)$ — yielding the
"lowest-degree-term-dominates / codimension / β-rescaling" classification. Four defects held it short:

1. **Stale labeling [load-bearing].** §E is phrased in the *higher-order-TW* reframe ("the inverted
   TW$^{(q)}_\beta$"), which **T2 refuted**. The cusp response must read $\mathcal W_\beta$
   **asymmetric-PIV-family**, consistent with the D-v2 atlas.
2. **The physical cusp (q=2) was never unfolded** — only the abstract q=3. The central physical singularity's
   own unfolding response was missing.
3. **Missing perturbation axes.** The user's list — forcing, coupling/β, **symmetry-breaking** — was not
   consolidated: β-response lived only in D; coupling/Δ only in the loop-closure; **symmetry-breaking
   (within-class) was entirely unstudied**.
4. **No consolidated perturbation → response → predicted-vs-measured table.**

## What was closed (two new validations + consolidation)

**(A) Physical cusp linear unfolding — relevant, flows to fold [VALIDATED].**
$V=\operatorname{sign}(Y)(Y^2+c_1|Y|)$, β=2, `unfolding_response_v2.py`:

| $c_1$ | 0 | 0.5 | 1 | 2 | 4 |
|---|---|---|---|---|---|
| skew | +0.60 | +0.59 | +0.53 | +0.41 | **+0.26** |
| exk | −0.25 | −0.01 | +0.10 | +0.12 | +0.07 |

Skew flows **cusp (0.60) → fold (0.26)** and the **sub-Gaussian signature dissolves** (exk −0.25 → ~0). The
linear term (degree 1 < 2) dominates near $Y=0$ → **relevant** perturbation → fold/TW class. Confirms
lowest-degree-term-dominates at the *physical* cusp, not just the abstract q=3.

**(B) Within-class symmetry-breaking — susceptibility of the asymmetric-PIV fingerprint [VALIDATED, NEW].**
$V=Y^2$ ($Y\!\ge\!0$), $-(1{+}a)Y^2$ ($Y\!<\!0$) — keeps leading degree 2 (stays cusp) but breaks
$Y\!\to\!-Y$:

| $a$ | −0.4 | −0.2 | 0 | +0.2 | +0.4 |
|---|---|---|---|---|---|
| skew | +0.47 | +0.52 | +0.60 | +0.68 | **+0.74** |
| exk | −0.61 | −0.42 | −0.25 | −0.10 | +0.04 |

**Susceptibility $d(\text{skew})/da=+0.36$**, monotone, and the skew **stays in the cusp class** (0.47–0.74,
never collapses to fold 0.20). So the asymmetric coefficient is a genuine *within-class* perturbation: it
**tunes the asymmetric-PIV fingerprint smoothly** rather than changing the universality class. Direct evidence
that the cusp's +skew is *tied to the odd/asymmetric structure* — making the escape side steeper ($a>0$)
raises the skew; flattening it ($a<0$) lowers it. This is the perturbation that probes the T2 asymmetric
connection most directly.

**(C)–(E) Consolidated from established results (re-labeled, T2-consistent):**
- **noise-symmetry β** [VALIDATED, D]: fold $\mathrm{TW}_\beta$ skew $\downarrow\beta$; cusp $\mathcal W_\beta$
  skew $\uparrow\beta$ (the **inversion**). `beta_axis.py`.
- **coupling $g\!\to\!\Delta(g)$** [VALIDATED, loop closure]: fold ($\Delta$ large) ↔ cusp ($\Delta\!\to\!0$)
  crossover; peel-off spread ×3.7, exk flips ≈0→negative. `loop_closure_scaling.py`.
- **forcing (OU)** [NUMERIC, surrogate]: builds the process covariance ($C(1)=0.76$), marginal unchanged
  (Weber). `CUSP_PROCESS_NOTES.md`. *Caveat:* imposed-OU surrogate, not the intrinsic coupled-FHN forcing.
- **higher unfolding $c_2$ (codim up, q=3 demo)** [NUMERIC + caveat]: stays in class, $\beta_{\rm eff}$ shifts.

## Consolidated response table

| perturbation | predicted response | measured | status |
|---|---|---|---|
| relevant unfolding (linear $c_1$) | cusp → fold (lowest-degree dominates; codim $q{-}1$) | skew 0.60→0.26, exk −0.25→0 | **VALIDATED** |
| within-class symmetry-breaking ($a$) | stays cusp; tunes asym-PIV skew | skew 0.47→0.74, $d/da=+0.36$ | **VALIDATED (new)** |
| noise-symmetry ($\beta$) | fold TW$_\beta$ ↓β; cusp $\mathcal W_\beta$ ↑β (inversion) | fold 0.20/0.16; cusp 0.26→0.86 | **VALIDATED (D)** |
| coupling ($g\!\to\!\Delta$) | fold ↔ cusp crossover | spread ×3.7, exk flip | **VALIDATED (loop)** |
| forcing (OU) | builds process covariance; marginal fixed | $C(1)=0.76$, Weber const | NUMERIC (surrogate) |
| higher unfolding ($c_2$) | stays in class; $\beta_{\rm eff}$ shift | skew rides β-family up | NUMERIC + caveat |

**Classification principle (corrected labeling).** The **lowest-degree active term controls the edge
universality class**. A generic perturbation (nonzero linear term) is **relevant** → flows to **fold/TW$_\beta$**;
reaching rung $q$ requires tuning away all lower terms — **codimension $q-1$** (Thom genericity, now read on
the *stochastic edge law*). *Within* a class, **symmetry-breaking** ($a$) and **noise-symmetry** ($\beta$)
tune the fingerprint along the family (the cusp's $\mathcal W_\beta$ asymmetric-PIV family) without changing
the class.

## E status: 66% → **84%**

The perturbation/unfolding response is now **consolidated across every axis** the organizing singularity
admits, the **physical cusp** is unfolded (not just the abstract q=3), **symmetry-breaking** is validated as a
new within-class susceptibility, and the labeling is **T2-consistent** (no higher-order-TW residue).

**What still stands between E and 100%:**
- *(shared frontier, not E-specific, ~7%)* the **forcing response uses an imposed-OU surrogate** — the
  intrinsic coupled-FHN forcing response is the same open object as the intrinsic process (A/B/D frontier).
- *(shared with T2, ~6%)* the **analytic response coefficient** for symmetry-breaking — the susceptibility
  $d(\text{skew})/da=+0.36$ is measured, not derived from the asymmetric-isomonodromy connection data (open T2).
- *(E-specific, closable, low-value, ~3%)* **exhaustive multi-direction unfolding** — only single-parameter
  sweeps validated; the full simultaneous-codim map is not enumerated.

The first two are genuinely-open frontier shared with the deep T2/intrinsic-process program; only the third is
an E-specific closable defect, and it is low-value (single-parameter sweeps already establish the
classification principle).
