# Notes — §2 Pearcey-projection test (stage 1 of 3)  →  FAIL

_June 2026. Figure `coupled-atlas/figures/pearcey_projection_test.png`; script
`pearcey_projection_test.py`. Per `regime-tests/RH_DIRECTION_NOVEL_ANGLES.md` §2._

**Claim tested (§2 route A for T2):** the cusp edge law is the slow-flow reduction of the **Pearcey**
process — integrate out the third/y-direction (the second cuspoid unfolding parameter) with the KP
slow-flow weighting; the reduced law should be the cusp (Weber). If so, PIV descends from a known parent.

## Construction status (honest)

- **Pearcey object [CITED structure / RECONSTRUCTED proxy].** The A₃ catastrophe ODE is
  $p'''=y\,p'+ixp$; I used its swept-inner-equation form $u'''=(Y-\eta\xi)u + y\,u'$ (the
  catastrophe-ladder Pearcey, `peeloff_catastrophe_ladder`), with $y$ the deformation. This is the
  faithful catastrophe-ladder Pearcey, **not** the cited RMT Pearcey kernel (which I cannot construct
  cleanly here without scipy). Flagged.
- **Slow-flow weighting [RECONSTRUCTED].** "Integrate out y" = mixture $\int L_y\,w(y)\,dy$ of the
  peel-off laws over the passage; $w$ is reconstructed. I scanned uniform windows, Gaussian fast-mode
  invariants, and the $|y|$-limits.
- The peel-off laws and mixtures are **[NUMERIC]**.

## Result — FAIL

Deformed-Pearcey peel-off vs deformation $y$ (cusp target: skew +0.61, exk −0.24, $\alpha_L$ 5):

| y | skew | exkurt | α_L |
|---|---|---|---|
| 0 (Pearcey) | +0.10 | −0.01 | 3.6 |
| 2 | +0.15 | +0.05 | **5.3** |
| 3 | +0.24 | −0.17 | 10.6 |
| 4 | +0.84 | +0.55 | 62 |
| 6 | +1.11 | −0.12 | — |

**No single $y$ and no reconstructed weighting hits the cusp fingerprint simultaneously:**
- where the tail matches ($y=2$, $\alpha_L=5.3$), the skew is Pearcey-like (+0.15, not +0.61);
- where the skew climbs toward +0.6 ($y\approx3.5$–4), the tail explodes ($\alpha_L\to62$) and kurtosis
  goes positive;
- the mixtures give skew/exk in {(+0.05,−0.27), (−0.13,−0.91), (−0.31,−0.08), (−0.24,−0.64),
  (+0.44,−1.00)} and $\kappa_5,\kappa_6$ of the wrong sign/size (cusp: $\kappa_5\!=\!-2.2,\kappa_6\!=\!-2.7$).

So **the slow-flow-reduced Pearcey does NOT reproduce the cusp law.** The cusp RH is **not** simply
"Pearcey, reduced" (by this projection). This is consistent with the established **"Weber, not Pearcey"**:
the cusp escape is a genuine 2nd-order (Weber, $\alpha_L=5$) object, and the 3rd-order Pearcey does not
reduce onto it under the $y$-deformation/projection.

**Honest caveat.** The negative is for the *reconstructed* catastrophe-ladder Pearcey + reconstructed
weighting; a fully *cited* RMT Pearcey-kernel projection is harder and not done. But the structural test
(does the 3rd-order Pearcey peel-off reduce to the 2nd-order Weber/cusp?) fails cleanly, and it is the
same structure the cited reduction would have to realize.

## Standing of the RH program after stage 1

| route | test | status |
|---|---|---|
| §1 Weber tube (T1) | tube-width falsifier | **PASS** (true tube finite/Weber; OU diverges) |
| §3-B β = isomonodromic time (T2) | β-flow ODE | **FAIL** (generic saturation, no PIV Riccati) |
| §2 cusp = Pearcey-reduced (T2) | this test | **FAIL** (no reduction reaches the cusp) |
| §4 cusp = multicritical edge across a Stokes ray (T2) | — | pending (stage 2) |

Two of the three T2 "known-parent" routes are now excluded; **T1 (the Weber tube) remains the
confirmed, tractable route.** For the exact law (T2), the remaining candidate is **§4 (Stokes
continuation)** — stage 2, next. If §4 also fails, the cusp PIV law is likely a genuinely *new* object
(not descended from Pearcey or the multicritical/PII hierarchy), to be built directly via the §0 uniform
parabolic-cylinder-with-noise parametrix.
