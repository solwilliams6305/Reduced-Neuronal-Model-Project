# Notes — the exact tail asymptotics of the cusp escape law (a piece of the exact law, nailed)

_June 2026. Figure `coupled-atlas/figures/cusp_tail_exponents.png`; script `cusp_tail_exponents.py`.
Tags: [DERIVED, exact at q=1] · [NUMERIC] · [CONJECTURE: PIV]._

Direction chosen: **nail the exact law**. First decisive piece — the **tail asymptotics**, in closed
form, via Freidlin–Wentzell on the (confirmed) Riccati explosion
$dp=(\mathrm{sign}(Y)Y^2-p^2)\,dY+\eta\,dW$, $\eta^2=2$.

## The derivation

**Left tail** ($s\to-\infty$, late escape). To survive un-exploded to level $s$, the noise must hold
$p\approx0$ (where $|{\rm drift}|=|Y|^q$ is minimal) against the downward drift. The FW action is
$$I_{\rm L}(s)=\frac{1}{2\eta^2}\int_0^{|s|}|Y|^{2q}\,dY=\frac{|s|^{2q+1}}{4(2q+1)},\qquad
F(s)\sim e^{-|s|^{2q+1}/(4(2q+1))},\quad \boxed{\alpha_{\rm L}=2q+1}.$$

**Right tail** ($s\to+\infty$, early escape). Barrier crossing over the unstable Riccati branch
$-\sqrt{W}$; barrier height $\Delta U=\tfrac43 W^{3/2}$ with $W=|s|^q$, giving
$$S(s)\sim e^{-(4/3)|s|^{3q/2}},\qquad \boxed{\alpha_{\rm R}=\tfrac32 q}.$$

## The anchor (validation): q = 1 is Tracy–Widom, exactly

At $q=1$ the formulas give $F\sim e^{-|s|^3/12}$ and $S\sim e^{-\frac43 s^{3/2}}$ — the **exact known
TW tails, exponents *and* constants** ($1/12$, $4/3$). So the method is validated on a known result,
not assumed.

## The cusp (q = 2)

$$\text{left } \sim e^{-|s|^{5}/20}\ (\alpha_{\rm L}=5),\qquad \text{right } \sim e^{-(4/3)|s|^{3}}\
(\alpha_{\rm R}=3).$$
Both far steeper than Airy/TW $(3,\tfrac32)$ — the law is strongly **sub-Gaussian** (hence the
negative excess kurtosis), with the right (early-escape) tail heavier than the left (hence the positive
skew). The $|s|^5$ left tail is the **parabolic-cylinder / Painlevé-IV** fingerprint (Weber : PIV ::
Airy : PII), a distinct universality class from the Airy/PII (TW) one.

## Numerical confirmation

High-stat Monte-Carlo of the Riccati explosion (raw units, deepest 6% tail) gives local left-tail
log-log slopes **2.75 at $q=1$** (derived 3) and **5.23 at $q=2$** (derived 5) — close to $2q+1$ and
cleanly separating the cusp ($\approx5$) from TW ($\approx3$). MC is pre-asymptotic so it slightly
under-resolves the absolute value, but it tracks the parameter-free WKB curves (figure panels B, C) and
confirms the exponent rises steeply with $q$.

## The ladder, and what's nailed vs open

General turning order $q$: $\alpha_{\rm L}=2q+1$, $\alpha_{\rm R}=\tfrac32 q$ — the edge-tail
asymptotics of noise-induced escape at a catastrophe of order $q$, with the fold (Airy/TW) the $q=1$
member. **Nailed:** the tail exponents and constants (analytic, exact at $q=1$, MC-supported).
**Still open (the full exact law):** identifying the complete distribution with a specific
**Painlevé-IV transcendent** / a parabolic-cylinder Fredholm determinant — the tails are necessary
conditions the PIV solution must meet, and they match the conjecture, but the determinant/transcendent
identification is the remaining step (Route B + the implicit explosion-PDE for the full CDF).
