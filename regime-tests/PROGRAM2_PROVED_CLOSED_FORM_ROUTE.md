# Program 2 — the route to a PROVED closed form of 𝒲 (resurgent trans-series). Deep-research synthesis + chart.

_2026-07-07 (Fable 5). After computing $v_2$ (§6 of `PROGRAM2_ROUTE2B_NOTES.md`), a `/deep-research` pass
(105 agents, 23 primary sources, 3-0 adversarial votes) mapped the rigorous machinery. This note charts the
proof route. Tags: **[PROVED-available]/[PARTIAL]/[BLANK]** per component, with the load-bearing step marked._

## 0. What "proved closed form" means here (fixed by prior theorems)
A literal 1-D-ODE closed form is **provably impossible** (isomonodromy fixed point ⇒ no Painlevé-σ reduction;
unbounded-below first-passage ⇒ no Fredholm/soft-edge determinant — `T2_CONNECTION_DATA_DERIVATION`,
`STAGE1_GAP_DETERMINANT_NOTES`). So the target is a **proved resurgent trans-series**:
$$\mathcal W_\beta \;\leftrightarrow\; \underbrace{\textstyle\sum_n v_n\,\eta^{2n}}_{\text{perturbative}}
\;\oplus\;\underbrace{e^{-S/\eta^2}(\cdots),\ S=s^5/10}_{\text{instanton}}\ \text{+ Borel resummation + Stokes automorphism, with a PROVED Stokes constant.}$$
Proving it = proving four things: **(A)** the coefficients $v_n$ are rigorously defined/computable; **(B)** the
series is Borel-summable/resurgent with Borel singularity at $S$; **(C)** the Stokes automorphism linking the
perturbative and instanton sectors is proved and its Stokes constant computed; **(D)** the trans-series equals the
true law.

## 1. Component status (from the research pass — cite)

| Component | Status | Key result / where |
|---|---|---|
| **Deterministic Weber Voros coefficient, closed form** | **[PROVED-available]** | Iwaki–Koike–Takei, Bernoulli-number form $\propto B_{2k}/(2k(2k{-}1)\nu^{2k-1})$, tied to Eynard–Orantin topological recursion. arXiv:1805.10945 (AHP 24:1305, 2023); Part II J. Integrable Sys. 4(1) 2019. |
| **Borel summability of deterministic WKB** | **[PROVED-available]** | Nikolaev, CMP 400 (2022), arXiv:2004.13367 — saddle-free Stokes-graph criterion, explicit Borel transform; Nemes — uniform in $x$ with error bounds. |
| **Resurgence of deterministic WKB (Weber incl.)** | **[PROVED-available, 2024]** | Nikolaev, "Geometry and Resurgence of WKB Solutions," arXiv:2410.17224 (2024): formal WKB solutions on Riemann surfaces are **resurgent** — Borel-summable in a.a. directions, endless continuation, geometric Borel plane/Stokes rays. **Upgrades the backbone from conjectured→proved.** |
| **Sharp caveat: Weber node is a turning point** | **[PARTIAL]** | Proved summability domain must be turning-point-free (Nikolaev Thm 5.9/Cor 5.29); for $Q=E-x^2/4$, $E>0$ naive WKB is NOT Borel summable (fixed real-axis singularities). The Sasaki transformation-series bridge to Weber was **refuted 0-3** as proved-summable — so the local-at-node piece is not yet closed. |
| **Stochastic / random-potential exact-WKB, Voros, resurgence** | **[BLANK]** | Confirmed genuine blank: authoritative surveys (Iwaki Houches-2024; Başar–Dunne–Ünsal arXiv:1308.1108) contain **no** white-noise potential / SDE first-passage / explosion treatment. **This is the core novel step.** |
| **Borel singularity pinned at the instanton action** | **[PARTIAL — proved finite-D]** | Spada (Serone–Spada–Villadoro): each Lefschetz thimble's saddle expansion is Borel-resummable exactly to the thimble integral (thimble = Laplace transform of the Borel transform of its own series). Berry–Howls: leading Borel singularity at inter-saddle action difference. **Exact template for $S=s^5/10$ — not yet lifted to infinite-D path integral/SDE.** SISSA thesis (Spada). |
| **Malliavin/Wiener-chaos rigor for the coefficients** | **[PARTIAL — tools exist, application undone]** | General machinery (Nourdin–Peccati Malliavin–Stein; Wiener-chaos moments/cumulants/diagram calculus) is rigorous & standard. Not yet applied to prove the Green's-function/Wick $v_n$ formulas converge to the true cumulants, nor to the boundary (pointwise-noise-at-node) term. No dedicated source surfaced. |
| **Precedent: proved resurgent closed form of a no-ODE transcendent** | **[BLANK — setup only]** | Ramirez–Rider–Virág (arXiv:math/0607331): $\mathrm{TW}_\beta$ = spectrum of the **stochastic Airy operator** $-\partial_x^2+x+\tfrac2{\sqrt\beta}b'_x$ = Riccati **first-explosion** law — exactly the random-Schrödinger + first-explosion setup $\mathcal W_\beta$ generalizes. **But no resurgent closed form was ever established for $\mathrm{TW}_\beta$** — so $\mathcal W$ would be first-of-kind; no completed template for the final assembly. |

## 2. The route (ordered; the one BLANK step is load-bearing)

- **Phase 0 — DONE (this program).** Deterministic connection closed-form (root $\lambda_0=0.890-0.890i$,
  machine-verified); three converged perturbative coefficients $v_0,v_1,v_2$; instanton action $S=s^5/10$;
  **deterministic Green's-function/Wick reformulation** (coefficients = exact iterated integrals, MC-free through
  $v_1$); explicit boundary term $Y_4=\tfrac13\xi(Y^\star_0)Y_1^3+$reg.
- **Phase 1 — rigor on the coefficients (A).** Prove the Green's/Wick multiple-integral formulas converge and equal
  the true cumulant coefficients to all orders, incl. the boundary local-time term. **Tool:** Nourdin–Peccati
  chaos/cumulant machinery + boundary local-time / Wick renormalization. **[tools exist — bounded task].**
  **UPDATE 2026-07-07:** Phase 1 is now known to require genuine **renormalization**, not just convergence — the
  variance coefficient $v_3$ is a **renormalization threshold**: individual node-shift pieces ($\operatorname{Var}
  (Y_4)\ni\xi(Y^\star_0)^2$, etc.) carry divergent boundary self-contractions $\delta(0)$ that cancel only in the
  full $v_3$ assembly (numerically demonstrated, `PROGRAM2_ROUTE2B_NOTES.md` §7). So the "Wick renormalization" in
  the tool list is load-bearing from $v_3$ on. **LEAD for Phase 3:** this boundary-self-contraction ambiguity — a
  QFT-like counterterm structure appearing exactly at the non-perturbative threshold — is a plausible concrete face
  of the **Stokes / resurgence** structure (Borel-summation ambiguity ↔ boundary renormalization ambiguity); worth
  pursuing as the bridge from the coefficient rigor to the stochastic Stokes automorphism.
- **Phase 2 — Borel structure (B).** Prove the large-order growth of $v_n$ fixes the Borel singularity at
  $S=s^5/10$ and the series' Borel summability. **Tool:** lift Spada's thimble = Laplace-of-own-Borel identity and
  Berry–Howls inter-saddle-action rule from finite-D to the Freidlin–Wentzell weak-noise path integral (or apply
  directly to the stochastic-Weber Riccati). **[PARTIAL — the finite-D theorem is the template; lifting is the work].**
- **Phase 3 — the stochastic exact-WKB / Stokes automorphism (C). ← LOAD-BEARING, the BLANK.** Build exact-WKB for
  the *random* Weber operator. **Highest-leverage transfer:** the noise is a perturbation OF the deterministic
  Weber operator whose **full resurgence is now a theorem (Nikolaev 2024)** — treat the Voros symbols as random
  variables dressing the proved deterministic Borel-plane geometry; the perturbative↔instanton Stokes automorphism
  is inherited from that geometry and noise-dressed. **Deliverable:** a proved **stochastic Stokes constant** (the
  noisy analog of $\lambda_0$) + the resurgence relation. **[BLANK — the genuinely new mathematics].**
  **UPDATE 2026-07-10 — first crack in the BLANK (Milestone 1).** The stochastic Stokes constant is now DEFINED
  as $\mathbb E[\text{random Voros datum}]$ and computed at $O(\eta^2)$: $\Omega_2=-0.451+0.352i$ (noise-averaged
  shift of the connection root $\lambda_0$), via the Green's/Wick 2nd-Wiener-chaos machinery, certified +
  independently validated (0.67%). Key structural finds: the $O(\eta^2)$ *mean* does NOT rotate the Borel phase
  toward $54^\circ$ (that is a *fluctuation*-sector effect, rms$|\delta\lambda|\approx1.34\eta$); and the Borel
  *location* (period) mean shift is $\delta(0)$-divergent (the v3 renorm threshold) while the connection datum is
  renorm-clean (non-local Ito functional). See `PROGRAM2_STOCHASTIC_STOKES_O_ETA2_NOTES.md`. Remaining: the
  fluctuation$\to$Borel-phase link, $O(\eta^4)$, and the $\Omega_2\leftrightarrow$ Hao-$\Omega{=}1$ mapping.
  **UPDATE 2026-07-15 — the "one remaining large build" (v₆) is DONE and it flips the phase verdict.**
  $v_6=-1.90\pm0.15$ exactly (grid-certified n=10…24, dual-engine, boundary budget 0.1%): the predicted sign
  flip to $+2$ is REFUTED, and the Borel-pair phase re-fits to $\theta=48^\circ\pm3^\circ$ (basin 45–50°,
  $\ge53^\circ$ excluded). **The "54–63°, NOT 45°" claim is retracted** — $\lambda_0$ inheritance ($45^\circ$)
  is back inside the error bars, which strengthens the Phase-3 inheritance strategy. Discriminating 45° vs 48–50°
  needs $v_7$ (predictions $+1.3$ vs $+8$; requires the DP-loop rewrite — $\sim14\times$ v₆ cost per grid as-is).
  See `PROGRAM2_V6_GRIND_NOTES.md`.
  **Same-day addendum:** the prescribed real-instanton subtraction was executed (z_r-profiled two-singularity
  fits, dual-implementation verified): the phase does NOT collapse to 45° — it pins at **θ ≈ 50° ± 2°** in every
  trustworthy fit, and the real-pole residue comes out negative (the 45°-masquerade needs positive). **Exact
  λ₀-inheritance is DISFAVORED; a genuine small ~5° fluctuation rotation is the leading hypothesis** — which
  re-poses Phase 3's target as: derive a mechanism giving (θ≈50°, |ζ|≈1.8–2.0) simultaneously. v₇ arbitrates
  (pre-registered: +1.3 @45° vs +8.2 @50°).
- **Phase 4 — closure (D).** Prove the assembled trans-series equals $\mathcal W_\beta$. Moment-determinacy is
  already in hand (Tier A closed). **[BLANK — first-of-kind; no $\mathrm{TW}_\beta$ precedent to copy].**

## 3. Reachable milestone vs far frontier (honest)
- **Reachable now, with existing tools:** a **rigorously-defined partial resurgent representation** — Phase 1
  (coefficients as proved convergent iterated integrals) + Phase 2 template (Borel singularity at $S$ by the
  finite-D thimble argument, at least heuristically-to-rigorously for the Riccati). This is a real theorem-grade
  deliverable and the natural next paper.
- **Far frontier (multi-year):** Phase 3 — **stochastic exact-WKB is a literature blank**; nobody has built
  exact-WKB/resurgence for a random-potential operator. This is the true cost of a *fully proved* closed form.

**One-line verdict:** the scaffolding around the blank is unusually complete — deterministic Weber resurgence is
**proved** (Nikolaev 2024, the backbone to build on), the instanton→Borel-singularity link is **proved in
finite-D** (Spada/Berry–Howls, the template), the chaos/cumulant coefficient tools are **standard**
(Nourdin–Peccati), and the operator precedent is **exactly parallel** (stochastic Airy/$\mathrm{TW}_\beta$). The
proof reduces to one genuinely new construction: **exact-WKB for a white-noise-dressed Weber operator**, i.e. the
stochastic Stokes automorphism and its Stokes constant.

## 4. Sources (primary, 3-0 verified)
Iwaki–Koike–Takei arXiv:1805.10945, JIS 4(1) 2019; Nikolaev CMP 400 (2022) arXiv:2004.13367; **Nikolaev
arXiv:2410.17224 (2024)**; Nemes; Iwaki Houches-2024 notes; Başar–Dunne–Ünsal arXiv:1308.1108; Spada SISSA thesis
(Serone–Spada–Villadoro); Berry–Howls; Ramirez–Rider–Virág arXiv:math/0607331. Refuted 0-3: Sasaki-series proved
Borel-summable to Weber; arXiv:1810.13158 as nearest weak-noise Borel work.
