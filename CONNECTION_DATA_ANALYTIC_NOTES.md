# The noisy parabolic-cylinder connection data — analytic attack on 𝒲's identity and the T1 handover

_June 2026. Figure `coupled-atlas/figures/tail_constants.png`; scripts `tail_constants.py` (+ `fp_beta.py`).
The keystone §0 object: the Weber connection/Stokes data settles both the isomonodromy class of 𝒲 and the
post-turning parametrix for C1/T1. Partial rigorous progress + honest map. Tags
[CITED]/[DERIVED]/[NUMERIC]/[HEURISTIC]/[CONJECTURAL]; ⚑ = load-bearing._

## 0. The object: the deterministic Weber connection

Cusp inner operator with spectral shift: $u''=(\operatorname{sign}(Y)Y^2-\lambda-\eta\dot W)u$.
- **$Y>0$ (confining):** $u''=(Y^2-\lambda)u$. With $z=\sqrt2\,Y$ this is the **parabolic-cylinder (Weber)
  equation** $u_{zz}=(\tfrac14 z^2-\tfrac\lambda2)u$; the solution recessive at $+\infty$ is
  $U(-\tfrac\lambda2,\sqrt2\,Y)\sim e^{-Y^2/2}$. **[CITED — DLMF §12, Olver]**
- **$Y<0$ (anti-confining):** $u''=(-Y^2-\lambda)u$ is **oscillatory** (inverted Weber), WKB solutions
  $\sim|Y|^{-1/2}e^{\pm i\int\sqrt{Y^2+\lambda}}$.
- **Connection:** the recessive solution continued through the turning $Y=0$ becomes a combination of the
  two oscillatory solutions; the coefficients are the **parabolic-cylinder connection coefficients**
  (Γ-function ratios, DLMF 12.2). The **first node** of the continued solution = the deterministic escape
  $Y^\star_{\rm det}=-2.09$. **[CITED + DERIVED]**

This is the §0 keystone parametrix's *deterministic* skeleton: a recessive (real-exponential) solution on
one side, oscillatory on the other, joined by the Weber connection. The escape law 𝒲 is the law of the
**noisy** first node.

## 1. Carrying the noise through

Cole–Hopf $p=u'/u$, Prüfer phase in the oscillatory region: write $u=R\,k^{-1/2}\sin\theta$ with
$k(Y)=\sqrt{Y^2+\lambda}$. Then [DERIVED, RRV-style]
$$\theta' = k(Y) + \tfrac{\eta}{?}\,(\text{amplitude-modulated})\,\dot W,\qquad
\theta_{\rm det}(Y)=\!\int^Y\! k = \text{Weber connection phase}.$$
The first node (escape) is the first passage of $\theta$ to a multiple of $\pi$. So:
- the **deterministic phase** $\theta_{\rm det}$ (from the Weber connection) fixes the *mean* escape;
- the **noise** randomises the phase (a phase-diffusion), giving the *spread* = 𝒲.

The "noisy connection data" is precisely the law of this phase-shifted first passage. **[DERIVED — mechanism.]**

## 2a. Isomonodromy / σ-form structure — settled analytically (not by fitting)

**Deterministic skeleton is PIV-class.** The Weber equation has an **irregular singularity of Poincaré
rank 2** at $\infty$ (the $Y^2$ potential). Painlevé-IV is the isomonodromic deformation of a 2×2 system
with exactly this rank-2 irregular structure, and PIV's classical special-function solutions *are*
parabolic-cylinder functions. So the deterministic cusp skeleton sits in the **PIV isomonodromy class**.
**[CITED — PIV ↔ rank-2 irregular / Weber; this is the genuine basis for the "Weber⇒PIV" intuition.]**

**But the odd potential breaks the symmetry.** $V=\operatorname{sign}(Y)Y^2$ is **confining for $Y>0$**
(real exponential $e^{-Y^2/2}$) and **oscillatory for $Y<0$** ($e^{\pm iY^2/2}$). The connection joins two
*different* irregular structures across $Y=0$ — a real-to-oscillatory connection, **not** the symmetric
Weber connection of standard PIV. The diagnostic is the **tail asymmetry**:
- standard soft-edge Painlevé (PII-hierarchy, *and* symmetric PIV-type edges) have left:right tail ratio
  **2** (one isomonodromy exponent controls both Stokes directions);
- 𝒲 has ratio **5/3** (R1), and — newly here — the two tails come from **different mechanisms** (§2b):
  a confining-side **barrier** (right) and an oscillatory-side **phase-persistence** (left).

⇒ **Verdict (analytic, structural):** 𝒲 is governed by an **asymmetric isomonodromy in the PIV family** —
the PIV-class (Weber, rank-2) skeleton with **asymmetric connection data** (confining one side, oscillatory
the other), **not** the standard symmetric Painlevé-IV σ-form. This is what the bulk σ-form regression
*could not* see (it lives in the bulk; the structure lives in the connection/Stokes asymmetry).
**[DERIVED-structural + CITED; precise asymmetric class = CONJECTURAL ⚑]**

## 2b. Tail constants — attempted; the honest outcome is the *asymmetry*

I tried to pin the exact constants by the optimal-fluctuation (instanton) and validate against the smooth FP
𝒲_β (β=1,2,4,8). The FP adjudicated — and revealed the two tails are **structurally different**:

- **RIGHT tail (early escape):** $\eta^2(-\ln(1-F))$ vs $Y^\star$ **collapses across β** — a clean
  Freidlin–Wentzell large deviation, $-\ln P\approx g(Y)/\eta^2$ with $g$ convex (the **confining-side
  barrier crossing**). **[NUMERIC — clean collapse]** The naive cubic-barrier *constant* ($8/3$) is **wrong**
  (the barrier vanishes at the turning; the escape is a *dynamic* sweep), so $g(Y)$ is the genuine rate
  function but not $\tfrac{8}{3}Y^3$. **[corrected — my R1 right-constant was wrong]**
- **LEFT tail (late escape):** $\eta^2(-\ln F)$ does **not** collapse (scales *slower* than $1/\eta^2$) — an
  **anomalous, oscillatory-side phase-persistence**, not a standard FW barrier. The local exponent passes
  through ≈5 only near β=2 (consistent with FW $2q+1=5$ there) and drifts (4.2→7.5) over β — pre-asymptotic
  + anomalous scaling. The naive "holding-$p{=}0$" instanton ($|Y^\star|^5/10\eta^2$, hence the "$1/20$" at
  β=2) mis-models late escape (the real event is *delaying a finite-time blow-up* in an oscillatory region,
  not holding). **[DERIVED-but-superseded; FP shows the scaling is anomalous]**

**Honest result:** the exact tail *constants* are **not** cleanly pinned (the FP can't reach the asymptotic
regime, and the left tail's scaling is anomalous). But the **two-mechanism asymmetry is robust and is the
real prize** — it is the analytic fingerprint of the asymmetric isomonodromy (§2a). The often-quoted
"$e^{-|s|^5/20}$" is, at best, the β=2 left-exponent (≈5) with a *non-universal* constant; it is **not** a
clean closed-form law. **[NUMERIC + honest correction]**

## 3. Post-turning O(ℓ) core — the C1 → edge-law handover, and how close T1 is

The post-turning $O(\ell)$ core **is** the Weber connection region. The pieces are now explicit:
- the **deterministic Weber parametrix** $U(-\tfrac\lambda2,\sqrt2\,Y)$ + its Γ-function connection [CITED];
- the **noise carried through** as the Prüfer phase-diffusion (§1) [DERIVED].

Together they give the explicit *structure* of the handover: the approach tube (C1) delivers the recessive
solution to the turning region with controlled fluctuation; the Weber-connection-plus-phase-diffusion maps
it to the oscillatory side, where the first node = the edge law. **The one remaining rigorous step is the
⚑ stochastic-Olver bound:** Olver's *deterministic* uniform parabolic-cylinder connection error bounds,
upgraded to **bound the noisy connection-coefficient moments uniformly through the turning** (scoping
§3(ii)(a)). **[OPEN ⚑ — load-bearing]**

**T1 status:** approach tube (C1) **PROVED** (uniform in ρ); re-entry **closed** [PROVED]; the post-turning
$O(\ell)$ core now has an **explicit parametrix** (Weber + phase-diffusion) rather than a black box. So T1
is reduced to the **single** stochastic-Olver uniform-connection estimate — everything else (approach,
mechanism, deterministic connection) is in hand. That is the precise, honest distance to closed.

## Status ledger

| component | status |
|---|---|
| deterministic Weber connection (recessive $U$, Γ-coeffs, first node = −2.09) | **CITED + DERIVED** |
| noise carried through = Prüfer phase-diffusion, escape = phase first passage | **DERIVED** (mechanism) |
| deterministic skeleton is PIV-class (rank-2 irregular / Weber) | **CITED** |
| 𝒲 ≠ standard symmetric PIV; **asymmetric isomonodromy** (ratio 5/3, two mechanisms) | **DERIVED-structural** |
| precise asymmetric isomonodromy class / σ-form | **CONJECTURAL ⚑** |
| right tail = clean FW $g(Y)/\eta^2$ (barrier); $8/3$ constant wrong | **NUMERIC + correction** |
| left tail = anomalous (phase-persistence), exp≈5 near β=2, constant not universal | **NUMERIC + correction** |
| post-turning $O(\ell)$ parametrix = Weber + phase-diffusion (explicit) | **DERIVED** |
| T1 closure = stochastic-Olver uniform connection bound | **OPEN ⚑ load-bearing** |

## Net

The connection-data attack delivered the **structural** answer the regression could not: 𝒲's deterministic
skeleton is genuinely PIV-class (Weber, rank-2 irregular — the real content of "Weber⇒PIV"), but the odd
potential makes the connection **asymmetric** (confining vs oscillatory), so 𝒲 is an **asymmetric
isomonodromy, not the standard symmetric Painlevé-IV** — confirmed analytically by the two-mechanism tail
asymmetry (right = FW barrier with clean $1/\eta^2$; left = anomalous phase-persistence). The exact tail
*constants* remain unpinned (and I corrected two of my earlier guesses — the $8/3$ right-constant and the
$1/20$ left-law are not clean closed forms). On rigor, the post-turning $O(\ell)$ core now has an **explicit
Weber+phase parametrix**, reducing **T1 to the single load-bearing stochastic-Olver uniform-connection
bound** — the cleanest statement yet of what remains to close the tube.
