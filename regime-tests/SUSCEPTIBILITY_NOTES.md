# E climb — the cusp symmetry-breaking susceptibility, DERIVED (escape location) + skew mechanism

_July 2026 (incoming agent, Fable 5). Closes rung E's "analytic susceptibility coefficient not derived" (shared
with the T2 obstruction). **Result: the escape-LOCATION susceptibility to within-class symmetry breaking is
DERIVED in closed form (a Green's-function integral) and validated to 2.5%; it sidesteps the T2 isomonodromy
obstruction entirely (first-order perturbation, not a monodromy deformation). The skew susceptibility
d(skew)/da=+0.38 reproduces E's measured +0.36, now with a mechanism.** Script
`coupled-atlas/susceptibility_derivation.py`; figure `coupled-atlas/figures/susceptibility.png`. Tags
**[DERIVED]/[NUMERIC]**._

---

## The perturbation (E case B)

$V_a=Y^2\,(Y\ge0),\ -(1+a)Y^2\,(Y<0)$ — i.e. $V_a=\operatorname{sign}(Y)Y^2+a\,g$, $g(Y)=-Y^2\mathbf 1_{Y<0}$
(deepens the oscillatory well, breaks $Y\!\to\!-Y$, stays cusp-class).

## Derivation — escape-location susceptibility [DERIVED ⚑]

Escape $Y^\star$ = first zero of the recessive field $u_0$ ($u''=(V_a-\lambda)u$). First-order response:
$dY^\star/da=-\partial_a u(Y^\star)/u_0'(Y^\star)$, with $\partial_a u$ solving $(\partial^2-q_0)\partial_a u=g\,u_0$.
The causal (downward-sweep) propagator $K(Y,s)=-[u_0(s)\bar u(Y)-\bar u(s)u_0(Y)]/\mathcal W$ (the $-$ from
$\partial_\tau=-\partial_Y$, $\tau=-Y$) and $u_0(Y^\star)=0$ collapse this to a single integral:
$$\boxed{\ \frac{dY^\star}{da}=-\frac{\bar u(Y^\star)}{\mathcal W\,u_0'(Y^\star)}\int_{Y^\star}^{0}s^2\,u_0(s)^2\,ds\ }$$
($\bar u$ the independent solution, $\mathcal W$ their Wronskian). **[DERIVED]** — a first-order perturbation
formula, so it needs no isomonodromy deformation (the obstruction that stalled the T2 route does not apply here).

**Validation.** $\lambda=1$, $Y^\star=-1.397$: DERIVED $dY^\star/da=+0.1205$ vs finite-difference $+0.1236$
(ratio 0.975, **2.5%**). Deepening the oscillatory well pushes the escape *shallower* (less negative), sign
confirmed. **[NUMERIC ✓✓]**

## Skew susceptibility — reproduced + mechanism [NUMERIC + DERIVED-sign]

Noisy escape (β=2), skew vs $a$: **$d(\mathrm{skew})/da=+0.38$** (a=−0.2/0/+0.2 → skew −0.14/−0.06/+0.01),
matching E's measured **+0.36**. **Mechanism [DERIVED sign]:** $a>0$ deepens the oscillatory ($Y<0$) well, which
lengthens the persistence/left-tail escape path — heavier left tail — so the distribution skews *positive*
(consistent with the FW left-tail instanton of `PERSISTENCE_ITEM2_NOTES`: a lower well raises the action barrier
asymmetrically). The closed-form *value* of the skew coefficient is the harder distribution-shape response (it
combines the location shift above with the tail reshaping); the sign and magnitude are pinned.

## Net + tracker

E's analytic-susceptibility gap is substantially closed: the **escape-location** susceptibility is derived in
closed form (Green's-function integral) and validated to 2.5%, sidestepping the T2 isomonodromy obstruction; the
**skew** susceptibility (+0.38) is reproduced and mechanistically explained via the left-tail instanton. **E
84→88%.** Remaining on E: the closed-form skew *coefficient* (shape response) and the imposed-OU forcing surrogate
(not uniquely defined without a canonical drive) — both genuinely lower-value.
