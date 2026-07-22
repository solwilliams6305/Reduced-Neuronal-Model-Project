# Notes — §4 Stokes-continuation test (stage 2 of 3)  →  FAIL (cusp-match); mechanism real

_June 2026. Figure `coupled-atlas/figures/stokes_continuation_test.png`; script
`stokes_continuation_test.py`. Per `regime-tests/RH_DIRECTION_NOVEL_ANGLES.md` §4._

**Claim tested (§4 route C for T2):** the cusp's sub-Gaussian flip = analytic continuation of a known
multicritical edge across a Stokes ray (the "wrong-sign"/Ablowitz–Segur branch vs Hastings–McLeod),
which turns heavy → sub-Gaussian. Pass = the continued multicritical law matches the cusp.

## Construction status (honest, no-scipy sandbox)

- **[CITED + computable]** the canonical Stokes split is Painlevé II (q=1): the TW edge is
  $F(s)=\exp(-\int_s^\infty(x-s)q(x)^2dx)$ with $q$ the Hastings–McLeod solution of $q''=sq+2q^3$
  ($q\sim\mathrm{Ai}$, $s\to+\infty$). The "wrong-sign" branch is the **defocusing** equation
  $q''=sq-2q^3$ (and the Ablowitz–Segur family $q\sim k\,\mathrm{Ai}$). I integrate these ODEs directly.
- **[BEYOND SANDBOX]** the cusp is the **q=2** member; the multicritical q=2 PII-*hierarchy* edge (a
  4th-order transcendent) and its continuation are not cleanly computable here (no scipy / special
  functions). So this decides the **mechanism** at q=1 and assesses cusp-consistency; the exact q=2 match
  is the residual gap.

## Result

| object | skew | exk | κ5 | κ6 | class |
|---|---|---|---|---|---|
| focusing / HM (q=1) = **TW** | +0.20 | **+0.13** | — | — | heavy (κ4>0) — validates pipeline |
| **defocusing / wrong-sign** (q=1) | **−0.24** | **−0.73** | +1.40 | +2.43 | **sub-Gaussian (κ4<0)** |
| **cusp target** (q=2) | **+0.61** | **−0.24** | **−2.2** | **−2.7** | sub-Gaussian |

- **MECHANISM is REAL:** crossing to the wrong-sign (defocusing) branch flips heavy ($\kappa_4=+0.13$) →
  **sub-Gaussian** ($\kappa_4=-0.73$). So a Stokes continuation *can* produce sub-Gaussian edge laws —
  the §4 idea is sound in principle. **[CITED + NUMERIC]**
- **CUSP-SPECIFIC match FAILS:** the wrong-sign q=1 law is **left-skewed** (−0.24) with $\kappa_5=+1.40$,
  whereas the cusp is **right-skewed** (+0.61) with $\kappa_5=-2.2$ — even the *signs* of skew and κ5
  disagree — and it is a q=1 object (tail exponent 3), not the cusp's q=2 (exponent 5). So the
  *accessible* continuation does not reach the cusp; the q=2 version is beyond the sandbox.

**⇒ FAIL** (cusp-match). The mechanism exists but the continued *q=1* law is not the cusp, and the q=2
continuation cannot be evaluated here. Per the note's own framing, this is the "fail" branch: **the cusp
PIV law is likely genuinely NEW**, to be built directly via the §0 uniform parabolic-cylinder-with-noise
parametrix.

## RH-program scorecard (all cheap falsifiers run)

| route | for | test | status |
|---|---|---|---|
| §1 Weber tube | **T1** | tube-width | **PASS** (true tube finite/Weber; OU diverges) |
| §3-B β = isomonodromic time | T2 | β-flow ODE | **FAIL** (generic saturation, no PIV Riccati) |
| §2 cusp = Pearcey-reduced | T2 | projection | **FAIL** (no reduction reaches the cusp) |
| §4 cusp = Stokes-continued multicritical | T2 | this test | **FAIL** (mechanism real; q=1 ≠ cusp; q=2 beyond sandbox) |

**Net conclusion.** **T1 (the Weber tube) is the confirmed, tractable route.** All three T2 "known-parent"
shortcuts (Pearcey reduction, β-isomonodromic-time, Stokes-continued multicritical) are excluded — the
cusp PIV law does not descend from a known parent we could reach. Combined with the determinant test
(the cusp is not a soft-edge Fredholm determinant), the weight of evidence is that **the cusp edge law is
a genuinely new object.** It therefore has to be built directly via the §0 keystone — the **uniform
parabolic-cylinder-with-noise parametrix** — which is exactly the §1 Weber-tube build (stage 3): the
single object whose probabilistic reading closes T1 and whose integrable reading is the new T2 law.
