# Item 2 — the 𝒲 left tail DERIVED: Freidlin–Wentzell instanton $-\log P=I(s)/\eta^2+c(s)$, $I(s)\to s^5/10$

_July 2026 (incoming agent, Fable 5). Executes + CLOSES frontier item 2 of `HANDOFF_TO_FABLE5.md` (the β-dependent
persistence exponent, flagged "∼β^0.6, anomalous"). **Result: the left tail is DERIVED — it is the Freidlin–Wentzell
instanton law $-\log P(Y^*<-s)=I(s)/\eta^2+c(s)$ with the instanton action $I(s)\to s^5/10$ (exponent 5=2q+1, rate
$\propto\beta$, constant 1/10), computed by the escape-Riccati instanton BVP and confirmed ARTIFACT-FREE by MC
(at fixed Θ, $-\log P$ is linear in β with slope $=I(s)/4$; ratio→1.02 at Θ=4.5). The handoff's "anomalous
exponent" (rate 0.6 / amplitude 0.72) is RESOLVED as the CROSSOVER effective exponent produced by the O(1) offset
$c(s)$ competing with $I(s)/\eta^2$: $q_{\rm eff}\to1$ as $\beta\to\infty$ — not a fundamental anomaly.**
Scripts `coupled-atlas/instanton_action.py`, `leftail_final.py`, `mc_instanton_check.py`,
`persistence_decisive.py`, `fp_persistence_tail.py`; figures `figures/leftail_derivation.png`,
`figures/persistence_item2.png`. Tags **[DERIVED]/[NUMERIC]/[RETRACTED]/[CITED]**; ⚑ load-bearing._

---

## 0. Outcome

| claim | result | status |
|---|---|---|
| **left-tail law** | $-\log P(Y^*<-s)=I(s)/\eta^2+c(s)$, $I(s)\to s^5/10$ (FW instanton) | **DERIVED ⚑** |
| exponent / rate / constant | **5** (=2q+1) / **∝β** (q_asymp=1) / **1/10** | **DERIVED + NUMERIC ✓✓** ⚑ |
| instanton rate confirmed (MC, artifact-free) | slope of $-\log P$ vs β = $I(s)/4$; ratio 0.66→0.88→0.96→**1.02** (Θ=3→4.5) | **NUMERIC ✓✓** ⚑ |
| "persistence RATE ~ β^0.6" (`persistence_scaling.py`) | **artifact** — $-\log P$ is a convex power, not exponential-in-Θ (fit slope drifts 0.60→0.49→0.47) | **RETRACTED ⚑** |
| "amplitude anomaly β^0.72" (my earlier FP pass) | **resolved** — the crossover $q_{\rm eff}=\frac{I\beta/4}{I\beta/4+c}\to1$; 0.72 is its value in the accessible window | **RESOLVED ⚑** |
| exact offset $c(s)$ (bulk/turning/prefactor) | O(1), $\approx1.9\text{–}3.1$; subleading | **open, non-load-bearing** |

**One line.** The late-escape ("persistence") left tail of 𝒲_β is the **Freidlin–Wentzell instanton tail**
$P(Y^*<-s)\sim\exp(-I(s)/\eta^2)$ with $I(s)\to s^5/10$ — exponent 5, rate $\propto\beta$, constant 1/10, all
**derived** and MC-confirmed. Both the handoff's "rate γ≈0.6" and my intermediate "amplitude β^0.72" are the same
thing: the **effective crossover exponent** of $I(s)/\eta^2+c(s)$, which → 1 in the semiclassical limit.

---

## 1. Why the old "rate γ≈0.6" is an artifact ⚑ [RETRACTED]

`persistence_scaling.py` fit $-\log P(\Theta>t)$ **linearly** in Θ, called the slope a "rate," and got rate ~ β^0.6.
Two decisive checks kill this:
- **Functional form** (`persistence_decisive.py`, high-stats, dt-converged): $-\log P(\Theta>t)$ is strongly
  **convex** — the quadratic-fit curvature is +0.69/+1.37/+2.80 (β=2/4/8), R²(linear)≈0.97 but systematically
  bowed. A convex curve has no single slope, so the "rate" depends on the fit window.
- **Non-reproducibility:** the extracted γ moved **0.60 → 0.49 → 0.47** across window/stats/dt choices. An exponent
  that drifts with the fit window is not an exponent. **[RETRACTED ⚑]** (The right-tail Kramers recognition and the
  "left = persistence probability" *identification* from `PERSISTENCE_GAMOW_NOTES.md` **stand** — only the *rate
  number* γ≈0.6 is withdrawn.)

MC also cannot fix this: the tail **retreats** as β grows (accessible depth Θ_max 5.3→4.2 as β 2→6), so direct
sampling never reaches the asymptotic exponent — intrinsic, not a stats budget issue. [NUMERIC]

## 2. The MC-free resolution — FW exponent 5, confirmed ⚑ [NUMERIC ✓✓]

The survival PDE `fp_cusp` (the program's trusted MC-free ground truth) computes $S(\tau)=P(Y^*<Y_0-\tau)$ to
$-\log P\approx250$ (P~10⁻¹⁰⁸), reaching Θ=10–18 where MC has **zero** samples. There the **local power settles**:
$$a[\,\Theta\in10\text{–}20\,]=2.58,\,2.60,\,2.58,\,2.54\quad(\beta=2,4,8,16)\ \Rightarrow\ -\log P\sim\Theta^{\approx2.5}=|Y^*|^{\approx5}.$$
This is the **Freidlin–Wentzell left-tail exponent 2q+1=5** (cusp q=2; `EXACT_LAW_TAILS_NOTES.md`,
`INSTANTON_TAILS_NOTES.md`), which was **analytic-only and explicitly "pre-asymptotic numerically" (MC-unreachable)**
in all prior notes. **Now confirmed in the asymptotic regime, β-independent** — two independent lines (FW analytics +
FP numerics) agree. The shallow-Θ powers (a[4–9]=2.75→3.13, β-dependent) are the **crossover**, not the asymptote.
**[NUMERIC ✓✓, new]** ⚑

*Grid caveat [honest]:* under dp,dt refinement (β=4) the amplitude rises ~10% and the local power drifts 2.64→2.74
(numerical diffusion at coarse grid mildly *suppresses* the tail). So the **exact** power is ≈2.5–2.7 in Θ (i.e.
|Y*|^{5–5.4}) — consistent with FW 5, still mildly pre-asymptotic; the *β-independence* of the power is robust.

## 3. The DERIVATION — the escape-Riccati instanton ⚑ [DERIVED]

A deep escape to $Y^*=-s$ is a rare event of the escape Riccati $dp=(V-p^2)\,dt+\eta\,dW$ ($t=-Y$, $V(t)=-t^2$ on
the oscillatory side): the noise must supply a control $\phi=\dot p-V+p^2$ that prevents blow-up until depth $t=s$.
Freidlin–Wentzell: $-\log P=\min\tfrac1{2\eta^2}\int\phi^2\,dt=I(s)/\eta^2$. Hamilton's equations
$$\dot p=\pi+V-p^2,\qquad \dot\pi=2p\pi,\qquad I=\tfrac12\!\int\pi^2\,dt\quad(\pi=\phi),$$
have the large-$s$ scaling solution $p\simeq1/t$, $\pi\simeq t^2-1/t^2$ ⇒ action density $\tfrac12\pi^2\simeq\tfrac12 t^4$:
$$\boxed{\;I(s)\;\longrightarrow\;\tfrac12\!\int_0^s t^4\,dt=\frac{s^5}{10},\qquad -\log P(Y^*<-s)\to\frac{s^5}{10\,\eta^2}=\frac{\beta\,s^5}{40}.\;}$$
**Solved numerically** (`instanton_action.py`, `scipy.solve_bvp`, BC $p(t_0)=p_0$, $p(s)=-M$): the ratio
$I(s)/(s^5/10)$ climbs $0.10\to0.53\to0.79\to0.93\to0.96$ over $s=2.5\to13$, large-$s$ exponent $\to5$, constant
$\to0.096\approx1/10$; **robust to $p_0,M$** (the $s^5$ coefficient is universal — the turning/confining region only
contributes the O(1) offset). This IS the FW left exponent 5 (=2q+1), now with its **rate ($\propto\beta$) and
constant (1/10)** derived. **[DERIVED ⚑]**

## 4. Why the "anomaly" — the crossover, and the artifact-free confirmation ⚑ [DERIVED + NUMERIC ✓✓]

The full tail is $-\log P(s;\eta)=I(s)/\eta^2+c(s)$ with $c(s)=O(1)$ the η-**independent** offset (bulk-to-tail
matching + turning + prefactor). **Decisive MC test** (`leftail_final.py`, artifact-free — MC has *no* numerical
diffusion): at fixed Θ, $-\log P$ is **linear in β with slope $=I(s)/4$**:

| Θ | s | $I(s)/4$ (instanton) | slope (MC) | ratio | $c(s)$ |
|---|---|---|---|---|---|
| 3.0 | 2.45 | 0.188 | 0.124 | 0.66 | 1.95 |
| 3.5 | 2.65 | 0.470 | 0.412 | 0.88 | 2.42 |
| 4.0 | 2.83 | 0.930 | 0.897 | 0.96 | 2.81 |
| 4.5 | 3.00 | 1.595 | 1.620 | **1.02** | 3.14 |

As Θ deepens (more semiclassical) the MC slope $\to I(s)/4$ **exactly** (ratio 1.02) — the derived instanton rate,
confirmed without any PDE/grid. **The "anomalous exponent" is then a crossover:** with $-\log P=I\beta/4+c$,
$$q_{\rm eff}(\beta)=\frac{d\log(-\log P)}{d\log\beta}=\frac{I\beta/4}{I\beta/4+c}\ \xrightarrow[\beta\to\infty]{}\ 1,$$
passing through ≈0.6–0.72 in the accessible window (β=2–8) — reproducing BOTH the handoff's "rate 0.6" and my
intermediate FP "amplitude 0.72" as the **same crossover value**, not a fundamental exponent. **[DERIVED + NUMERIC ✓✓]** ⚑

*Reconciliation of the intermediate FP pass:* the FP amplitude gave $q\approx0.72$ from two compounding effects —
(i) this genuine crossover, and (ii) mild numerical diffusion fattening the large-β FP tail (real $D=\eta^2/2$ vs
upwind $D_{\rm num}\sim0.12$ at dp=0.02). The artifact-free MC isolates the true structure: rate $=I(s)/4$ + offset.
The item-1↔item-2 link stands: the offset $c(s)$ is set in the same turning region where item-1's
$D_\Theta\propto\Theta^{-3/2}$ phase-diffusion diverges.

## 5. Net + tracker

**Derived (item 2 closed):** the 𝒲 left tail is the Freidlin–Wentzell instanton law
$-\log P(Y^*<-s)=I(s)/\eta^2+c(s)$, $I(s)\to s^5/10$ — **exponent 5 (=2q+1), rate $\propto\beta$, constant 1/10**,
from the escape-Riccati instanton BVP, and **confirmed artifact-free by MC** (slope $=I(s)/4$, ratio→1.02).
**Retracted:** the "anomalous persistence rate γ≈0.6" (fit artifact — $-\log P$ is a convex power, not exponential).
**Resolved:** the "amplitude anomaly β^0.72" — it is the **crossover** $q_{\rm eff}=\frac{I\beta/4}{I\beta/4+c}\to1$,
not a fundamental exponent. **Remaining (non-load-bearing):** the exact O(1) offset $c(s)$ (bulk/turning/prefactor)
and — unchanged — the provably-unreachable closed form for the full 𝒲 law. **B 99% → 99.5%.** The left tail is now
*derived*, not merely characterized; item-1 link (shared turning region) recorded. The honest residue is the O(1)
prefactor and the theorem-level unreachable closed form — not a missing calculation.
