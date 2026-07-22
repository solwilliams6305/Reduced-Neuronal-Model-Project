# The catastrophe ladder of noise-induced escape — synthesis

_June 2026. The capstone of the theoretical spine. Figure
`coupled-atlas/figures/catastrophe_ladder_atlas.png`; script `catastrophe_ladder_atlas.py`.
Tags: [DERIVED] skeleton+tails · [NUMERIC] moments · [CONJECTURE] PIV · [OPEN] general-q transcendent._

## One family

Every rung is the noise-induced escape (first explosion of the swept Riccati diffusion) at a turning
point of order $q$:
$$u''=(\mathrm{sign}(Y)|Y|^{q}-\eta\xi)\,u \ \Longleftrightarrow\
dp=(\mathrm{sign}(Y)|Y|^{q}-p^{2})\,dY+\eta\,dW .$$
$q$ indexes the $A_{q+1}$ singularity: $q{=}1$ fold, $q{=}2$ cusp, $q{=}3$ swallowtail, $q{=}4$
butterfly, … . The fold is Tracy–Widom; the rest are new laws of the same family.

## The ladder (β = 2)

| q | singularity | deterministic skeleton | $\alpha_{\rm L}=2q{+}1$ | $\alpha_{\rm R}=\tfrac32q$ | skew | exkurt | Painlevé |
|---|---|---|---|---|---|---|---|
| 1 | fold $A_2$ | **Airy** (Bessel 1/3) | 3 | 1.5 | +0.21 | +0.05 | **PII** (rigorous = TW) |
| 2 | cusp $A_3$ | parab.-cyl. / Weber (Bessel 1/4) | 5 | 3 | +0.61 | −0.24 | **PIV** (conjectured) |
| 3 | swallowtail $A_4$ | Bessel 1/5 | 7 | 4.5 | +0.95 | +0.16 | open |
| 4 | butterfly $A_5$ | Bessel 1/6 | 9 | 6 | +1.20 | +0.62 | open |
| 5 | wigwam $A_6$ | Bessel 1/7 | 11 | 7.5 | +1.42 | +1.26 | open |

## Three structural layers, each unified across the ladder

**1. Deterministic skeleton = Bessel of order $1/(q{+}2)$.** $u''=\mathrm{sign}(Y)|Y|^qu$ solves in
$\sqrt{|Y|}\,\mathcal C_{1/(q+2)}\!\big(\tfrac{2}{q+2}|Y|^{(q+2)/2}\big)$; the $\eta\to0$ peel-off levels
are its zeros. Verified: $q{=}1$ reproduces the **exact Airy zeros** $-2.339,-4.088,-5.521$, and the
WKB law $|Y_n|=[\tfrac{(q+2)\pi(n-1/4)}{2}]^{2/(q+2)}$ matches the computed levels to $\sim0.01$ for all
$q$. (So "fold = Airy" generalizes to "rung $q$ = Bessel $1/(q{+}2)$".)

**2. Exact tails (Freidlin–Wentzell).** $F(s)\sim e^{-|s|^{2q+1}/(4(2q+1))}$ (left) and
$S(s)\sim e^{-(4/3)|s|^{3q/2}}$ (right). At $q{=}1$ these are TW's exact tails, **constants included**
($1/12$, $4/3$). The exponents grow with $q$ — the laws get steeply lighter-tailed.

**3. Edge-law moments.** Skew rises monotonically (+0.21 → +1.42). Excess kurtosis traces a
**sub-Gaussian valley** localised to the **cusp region** ($q\approx1.5$–$2.5$, min ≈ −0.26): the cusp
is the distinguished member where the law is genuinely platykurtic. For higher catastrophes the large
skew drags excess kurtosis back positive (light tails, but strongly asymmetric).

## Where it sits — and what it is *not*

- $q{=}1$ — **rigorous**: TW, the Airy point process / Painlevé II (the paper + RRV).
- $q{=}2$ — **conjectured**: Weber-TW, parabolic-cylinder / Painlevé IV (Weber : PIV :: Airy : PII),
  supported by the exact $|s|^5$ left tail and the parabolic-cylinder skeleton.
- $q\ge3$ — **open**: swept-Bessel laws; the Painlevé/transcendent identification is unknown.

**Crucially this is a *new* family, not the multicritical soft-edge (Painlevé-II) hierarchy.** Those
higher-Airy edges are heavier-tailed (positive kurtosis, soft edges); ours are **swept** turning points
giving **sub-Gaussian** escape laws near the cusp. They coincide only at $q{=}1$ (Airy/TW) and diverge
above it — exactly the "Weber, not Pearcey / not higher-Airy" distinction, now ladder-wide.

## Nailed vs open

- **Nailed [DERIVED]:** the family and its generator (Riccati explosion); the Bessel-$1/(q{+}2)$
  deterministic skeleton (Airy zeros at $q{=}1$, WKB-verified); the exact tail asymptotics
  $2q{+}1,\ \tfrac32 q$ (exact TW at $q{=}1$); the moment atlas (skew↑, the cusp kurtosis valley).
- **Open:** the full exact law beyond the tails — the Painlevé-IV transcendent / parabolic-cylinder
  Fredholm determinant at $q{=}2$, and the general-$q$ transcendents.

## Physical realisation

$q$ is not abstract: the coupling $g$ sets the turning structure via $\Delta(g)=2\sqrt{-2g/3}$ (rung C),
so tuning electrical coupling in the two-neuron model moves the system **along the ladder** — fold
(TW) escape statistics at weak coupling, through the cusp (sub-Gaussian, Weber-TW) at $g_{\rm crit}$.
The kurtosis sign of spike/MMO-timing escape distributions is the observable signature of which rung.
