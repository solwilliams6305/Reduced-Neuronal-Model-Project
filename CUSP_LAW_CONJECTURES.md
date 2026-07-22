# What is the combined law? — fingerprint + creative conjectures

_June 2026. Figures `coupled-atlas/figures/cusp_law_fingerprint.png`; scripts `cusp_law_fingerprint.py`,
`split_gnd_test.py`. Tags: [NUMERIC] findings · [CONJECTURE] clearly marked · honest negatives kept._

The question: what is the exact "combined law" — the cusp (Weber-TW) escape law, and the family it
belongs to? Here is the full fingerprint and what it forces.

## The fingerprint (high-stat MC of the Riccati explosion)

**1. Cumulants (β=2), along the ladder.** Standardized κ₃,κ₄,κ₅,κ₆:
| q | κ₃ | κ₄ | κ₅ | κ₆ |
|---|---|---|---|---|
| 1 (fold = TW) | +0.20 | +0.01 | −0.19 | −0.34 |
| 2 (cusp) | **+0.61** | **−0.24** | **−2.15** | **−2.69** |
| 3 (swallowtail) | +0.94 | +0.14 | −3.20 | −9.02 |

The cusp has **large negative κ₅, κ₆** (TW's are tiny) — a strongly structured, light-tailed law.

**2. The β-family is INVERTED relative to TW.** Cusp (q=2), varying noise β=4/η²:
| β | skew | exkurt |
|---|---|---|
| 1 | +0.27 | **−0.79** |
| 2 | +0.61 | **−0.23** |
| 4 | +0.87 | **+0.67** |

TW runs the *other way* (skew falls with β: 0.29→0.22→0.17). For the cusp, **more noise ⇒ less skew
and more platykurtic**; less noise ⇒ sharp and skewed. So the negative kurtosis is **not universal** —
it is a *moderate-noise crossover feature*. (This is consistent with the neuron loop closure, where the
model's noise sits in the negative-κ₄ regime.)

**3. The quantile map to TW is near-cubic.** $Q_{\rm cusp}(p)\approx\varphi(Q_{\rm TW}(p))$ with
$\varphi$ a cubic: residual **0.054** (vs affine 0.169), $\varphi(x)\approx x+0.072x^2-0.035x^3$. The
cusp is a **low-order polynomial deformation of TW**, not a wildly different object.

**4. Honest negative — the tails do NOT pin the body.** The (5,3)-split generalized-normal (built from
the exact FW tails $e^{-|s|^5/20}$, $e^{-(4/3)|s|^3}$) gives, at the cusp skew +0.61, exkurt **−0.51** —
but the cusp is **−0.24**. So the bulk is "fuller" than the tails imply: **the integrable content lives
in the body, not the tails.**

**5. Kurtosis decomposition.** Along the β-family, $\ \kappa_4 \approx 2.1\,\kappa_3^{2}-0.93\ $
(fits β=1,4 to ±0.01). An **intrinsic sub-Gaussian baseline −0.93** (the turning-order/oscillation)
plus a universal skew-induced term.

## The conjectures

**C1 — Painlevé IV, read off the cubic deformation.** [CONJECTURE] The exact cusp law is the
Painlevé-IV Fredholm-determinant law (Weber : PIV :: Airy : PII). Evidence beyond the parabolic-cylinder
skeleton and the $|s|^5$ tail: the quantile map to TW is **cubic** (#3), and the body-not-tails result
(#4) says the integrable structure is in the bulk — exactly where a Painlevé transcendent acts. The
cubic+quadratic $\varphi$ mirrors PIV's nonlinearity ($\tfrac32 w^3+4xw^2$) deforming PII's ($2w^3+xw$).

**C2 — the master object is a 2-parameter family $\mathrm{TW}^{(q)}_\beta$ ("swept-Bessel Tracy–Widom").**
[CONJECTURE] Turning order $q$ (skeleton = Bessel $1/(q{+}2)$; exact tails $2q{+}1,\ \tfrac32q$) × noise
$\beta$. The $\beta$-axis interpolates a **deterministic-skeleton-dominated** law ($\beta\to\infty$:
sharp, skewed, leptokurtic — the Bessel-zero fluctuations) and a **noise-dominated** law ($\beta\to0$:
platykurtic, the bounded-passage-window geometric law). $\mathrm{TW}^{(1)}_\beta=$ Tracy–Widom. The cusp's
sub-Gaussianity is the bounded-window character winning at moderate $\beta$.

**C3 — it is genuinely "beyond TW", and NOT the multicritical (PII-hierarchy) edge.** [NUMERIC+CONJECTURE]
The cusp is more skewed *and* has flipped κ₄ — so it is not a Gaussian-smeared TW, and (negative κ₄ vs
the soft-edge's positive) not a member of the multicritical soft-edge / Painlevé-II hierarchy. It is a
new **swept-catastrophe** family. This is the "Weber, not Pearcey / not higher-Airy" point, now at the
level of the full law.

**C4 — the physical mechanism (and a wilder idea).** [HEURISTIC] The β-interpolation is a competition
between the **sharp canard escape** (small noise → the deterministic parabolic-cylinder asymmetry) and a
**bounded-window** escape (large noise → near-uniform over the passage, platykurtic). *Wild
[SPECULATION]:* the cusp law may admit a **subordination** representation — TW (or its $q$-analogue)
under a noise-controlled random time/edge — since the β-family looks like a bare skeleton law
progressively "dressed" by noise toward the bounded-window limit. Testable by checking whether
$\mathrm{TW}^{(2)}_\beta$ is a mixture/time-change of a single bare law.

## What this buys the exact-law program
The target is sharpened to: **a Painlevé-IV transcendent whose tails are the FW $|s|^5$/$|s|^3$ ones and
whose bulk is the cubic deformation of TW's** — embedded in the 2-parameter $\mathrm{TW}^{(q)}_\beta$
family. Next decisive test: build the PIV / parabolic-cylinder Fredholm determinant and check its
quantiles against the measured cubic $\varphi$ and the cumulants κ₅,κ₆ — a far tighter test than moments
alone.
