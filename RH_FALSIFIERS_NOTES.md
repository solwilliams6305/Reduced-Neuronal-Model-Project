# Notes — the two cheap RH-direction falsifiers (run per regime-tests/RH_DIRECTION_NOVEL_ANGLES.md)

_June 2026. Figures `coupled-atlas/figures/tube_width_test.png`, `coupled-atlas/figures/beta_flow_ode_test.png`;
scripts `tube_width_test.py`, `beta_flow_ode_test.py`. Tags: [NUMERIC]/[DERIVED]/[HEURISTIC]._

The note reframes the two open theorems — (T1) a uniform Berglund–Gentz tube through the fold merge, and
(T2) a Painlevé-IV τ-function for the cusp edge — as two readings of one **uniform parabolic-cylinder
parametrix carrying noise**, and proposes two cheap falsifiers to decide which RH route is real. Both run.

## Falsifier #1 — tube width: Weber vs OU  →  **PASS** (T1 derisked)

Test (§1): simulate the noisy inner Riccati $p'=(V_\Delta(Y)-p^2)-\eta\xi$,
$V_\Delta=\mathrm{sign}(Y)|Y|(|Y|+\Delta)$; measure the true fluctuation tube width std$(p)(Y)$ in the
decaying region (full survival, $Y>0$); compare to the quasi-static OU comparison
$v_{\rm qs}=\eta^2/(4\sqrt{V_\Delta})$ (the frozen-coefficient Gaussian the Berglund–Gentz bound uses).

Result (η=√2, at $Y=0.05$ just above the turning):

| Δ | true tube width | OU comparison | OU/true |
|---|---|---|---|
| 2.0 | 0.71 | 1.25 | 1.8 |
| 0.5 | 0.86 | 1.74 | 2.0 |
| 0.1 | 0.94 | 2.40 | 2.6 |
| 0.0 (cusp) | **0.95** | **3.16** | 3.3 |

The **true tube SATURATES** (0.71→0.95, ×1.34 — finite through the merge Δ→0); the **quasi-static OU
comparison DIVERGES** (→3.16 at Y=0.05, →5.8 at the turning, and grows as Δ→0). **[PASS]** — the tube is
finite (Weber-like) exactly where the OU bound blows up, so OU is the wrong comparison and a **Weber
(parabolic-cylinder) Green's-function comparison is the right tool**, as §1 predicts.
`[NUMERIC]` true tube; `[DERIVED]` OU bound $\eta^2/(4\sqrt V)$; `[HEURISTIC]` Weber identification.
*Honest nuance:* the divergence is in the quasi-static OU **bound/estimate**, not the swept dynamics
(which are finite); this is precisely why the standard Berglund–Gentz estimate degenerates at the cusp.
⇒ **T1 is the tractable theorem and is derisked: the uniform tube exists; build the Weber comparison.**

## Falsifier #2 — β-flow ODE: is β the PIV isomonodromic time?  →  **FAIL** (route B disfavored)

Test (§3): take the leading cubic-deformation amplitude φ(β) = cusp skew(β) across β=4/η², and ask
whether φ(β) obeys the PIV Hamiltonian-flow ODE — specifically the *non-autonomous* Riccati
$\varphi'=a\varphi^2+b\,\beta\varphi+c$ with the PIV "$2tw$" signature $b/(2a)\approx1$.

Result: φ(β) = 0.24, 0.43, 0.57, 0.65, 0.78, 0.86, 0.89, 0.93 for β=1…6.5 — a smooth **saturating** curve
(→ φ_∞≈0.95). The discriminating checks all come out negative:
- nonlinear (Riccati/cubic-field) term **not needed**: quadratic resid 0.024 ≈ linear resid 0.024 (a
  *linear/exponential* relaxation already fits);
- explicit-β (isomonodromic-time) structure **not better**: PIV-form resid 0.025 ≥ autonomous 0.024;
- PIV "$2tw$" signature **absent**: $b/(2a)=-0.02$ (not ≈1).

**[FAIL]** — β is **not** the PIV isomonodromic time. φ(β) is generic saturating relaxation to a fixed
point (the bare-skeleton skew 0.95); the recurring cubic in the QQ is the **generic leading non-Gaussian
(Cornish–Fisher) correction**, not a PIV-time Hamiltonian flow. `[NUMERIC]` φ(β)+fits; `[HEURISTIC]`
interpretation. *Caveat:* tested on the leading amplitude (skew); a subleading coefficient could in
principle differ, but the leading-order negative is strong, and 8 points limit power.

## Net implication for the RH program

- **Pursue T1 via the §1 Weber tube** — the cheap test confirms the tube is finite and the OU comparison
  is the wrong tool; constructing the **Weber Green's-function comparison** (then splicing OU∪Weber, à la
  Olver's uniform matching at the covariance level) is the tractable next theorem.
- **Drop §3 route B** (β = isomonodromic time): disfavored. For T2, redirect to **§2 (PIV as the
  slow-flow reduction of the Pearcey RH problem)** or **§4 (cusp = multicritical edge continued across a
  Stokes ray)** — the remaining cheap tests.
- These are exactly the "which route is real" answers the falsifiers were meant to give, **before**
  committing to the full uniform parabolic-cylinder-with-noise parametrix (§0).

**Status:** the analytic core's RH direction is now sharpened — T1/Weber-tube confirmed as the route;
T2/β-PIV-time excluded. The §0 keystone (uniform parametrix) and the §2/§4 T2 routes remain open.
