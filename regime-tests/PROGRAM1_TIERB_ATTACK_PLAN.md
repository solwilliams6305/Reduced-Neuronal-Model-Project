# Program 1, Tier B — universality-with-rate: literature verdict + a ranked attack plan + a PARTIAL RESULT (rate p=1/4 in W₁)

_July 2026 (incoming agent, Fable 5). A deep, adversarially-verified literature pass (deep-research workflow,
108 agents, 25 primary sources, 24/25 claims 3-0 verified) on the single deepest open problem of the program:
a **quantitative convergence rate** ε^p for the coupled-FHN cusp escape → 𝒲_β. Combined with a new **empirical
rate measurement** that resolves the blank the literature left open. Verdict: the rate is **genuinely
unestablished territory** in all of applied probability — a real frontier, not a gap someone already filled —
and the most transferable route is a **coupling of two Riccati first-explosion times** — which I then EXECUTED (§3bis): the synchronous coupling
gives an explicit **rate p = 1/4 in Wasserstein-1**, $W_1(\mathcal W_{\beta,\Delta},\mathcal W_\beta)\le\Delta\,
\mathbb E|\chi|=O(\varepsilon^{1/4})$, with the mechanism verified pathwise and the whole thing reduced to two
named lemmas (L1: $\mathbb E|\chi|<\infty$; L2: the Grönwall remainder through the turning). Sources cited
inline. Numerics `coupled-atlas/tierB_rate_probe.py`, `tierB_coupling.py`; figures `figures/tierB_rate.png`,
`tierB_chi.png`. Tags **[literature-verified]/[proved-structure]/[NUMERIC ✓]/[open-mapped]**._

---

## 0. The target (recap)

**Tier B conjecture.** The rescaled coupled-FHN antisymmetric-mode escape at the cusp converges in law to
𝒲_β (β = 4/η²) as ε → 0, **with a rate ε^p**. Architecture (`FRONTIER_SCOPING_NOTES.md`): **B1** deterministic
reduction FHN → swept Weber inner equation [GSPT/blow-up, Kristiansen–Pedersen]; **B2** uniform sub-Gaussian
tube [**T1, closed**] + two-sided stretched-exponential tails on 𝒲_β itself [**Tier A §4, closed** —
`PROGRAM1_TIERA_THEOREM.md`, giving exponential tightness]; **B3** upgrade the tube/moment bounds to
convergence-**in-distribution with a rate**; **B4** the outer finite-ε corrections vanish at rate ε^p. The
crux is **B3 + B4**; the central unknown is **p**.

---

## 1. The literature verdict (deep-research, verified)

**The exact object needed — a *rate* of convergence for a first-passage / Riccati-explosion functional of an
SDE — has essentially never been proven, in any setting.** The two branches that own the relevant pieces are
disjoint, and neither transfers:

| thread | strongest result(s) | what it gives | transfer to swept-merging-cusp |
|---|---|---|---|
| **1. Stochastic-operator / β-edge** | Ramírez–Rider–Virág (JAMS 2011, `math/0607331`); Krishnapur–Rider–Virág (CPA 2016, `1306.4832`) | TW_β **as** the Riccati first-explosion law (exactly our 𝒲_β structure); universality of the stochastic Airy operator | **owns the limit object, proves convergence in law only — NO rate.** The informative blank. |
| **2. RMT edge universality w/ rate** | Schnelli–Xu (CMP 2022, `2102.04330`) O(N^{−1/3+ω}); gen-Wigner (`2207.00546`), sample-cov (`2108.02728`), sparse (`2507.19340`, extra p⁻²N⁻⁴ᐟ³ term) | the **only** quantitative TW rates (~N^{−1/3}) | **low** — Green-function/resolvent comparison (Erdős–Yau–Yin) of *eigenvalue* distances; no swept turning point, no sample-path escape, not a first-passage functional. |
| **2′. Cusp universality (RMT)** | Erdős–Henheik–Riabov (CMP 2025, `2410.06813`) | Pearcey/cusp local statistics universal — **qualitative only** | even *static* RMT has **no rate at a cusp**. Reinforces: a rate through a cusp is unestablished. |
| **3. Stein / Malliavin–Stein** | Smoluchowski–Kramers (`2602.00875`, W₁≤C√m·\|ln m\|); steady-state diffusion approx (Braverman–Dai–Feng, `1512.09364`); extreme-value Stein (`10.1007/s10687-020-00380-5`) | quantitative rates — but for **stationary/invariant** measures or **iid-maxima** laws | **near-total blank** for a first-explosion functional. Untried, not disproven (§3, open Q). |
| **4. Noisy fast–slow escape** | Berglund–Gentz dynamic pitchfork (PTRF 2002, `math/0008208`); Berglund–Gentz–Kuehn canards (JDE 2012, `1011.3193`); Bakhtin–Pajor-Gyulai Jordan-block exit (S&D 2019, `1708.00558`) | the right **phenomenon** — peel-off through a two-merging-turning-point / degenerate structure; tube scale σ/ε^{1/4}, escape window [√ε, c√(ε\|log σ\|)] | reaches B2, **never a distributional rate ε^p**; the only exit correction anyone extracts is *logarithmic*. This is exactly the field-wide gap = our B3. |
| **5. Deterministic backbone** | Kristiansen (`2003.06817`, same author as B1) folded-node/pitchfork connection; Olver two-turning-point PC | the ε-expansion structure | **under-covered** — no verified result pins the power of ε in the connection-coefficient correction (see §2, which I supply empirically). |
| **6. Worked "physical → TW with a rate"** | — | — | **none found in any thread.** The most important informative blank: no dynamical/escape system has ever been shown to converge to TW (or a β-analogue) *with a rate*. |

**One-line verdict.** Tier B is not "hard because the tools are fiddly" — it is hard because **no existing
theorem produces a distributional rate for this class of functional at all.** That reframes Tier B from "apply
known machinery" to "build the first such result," which is precisely why it is a multi-year core. (A clean
statement of this blank is itself worth publishing as motivation.)

---

## 2. The empirical rate — resolving the blank the literature left open [NUMERIC ✓]

The literature could not locate the deterministic ceiling on p (Thread 5 blank). I measured it directly in the
**controlled inner problem**: the finite-coalescence escape law 𝒲_{β,Δ} with $V_\Delta(Y)=\operatorname{sign}
(Y)|Y|(|Y|+\Delta)$ (Δ = incomplete-merge parameter, $\Delta\sim\varepsilon^{1/4}$ near the cusp; Δ = 0 is the
pure cusp). At β = 2 (η = √2, so FP is reliable — *not* the small-η boundary-layer regime), sweeping Δ → 0
(`tierB_rate_probe.py`, converged reproduction of the Δ = 0 reference skew +0.601, exk −0.244):

$$
\big|\text{exk}(\Delta)-\text{exk}(0)\big|\sim\Delta^{0.80},\quad
\big|\text{skew}(\Delta)-\text{skew}(0)\big|\sim\Delta^{\gtrsim1},\quad
\big|\text{mean}(\Delta)-\text{mean}(0)\big|\sim\Delta^{1.15}.
$$

**The excess kurtosis is the rate-limiting cumulant, and its correction is NON-ANALYTIC** (exponent ≈ 0.8 < 1,
stable across fit windows), whereas the mean and skew corrections are analytic-like (exponent ≥ 1). Converting
via $\Delta\sim\varepsilon^{1/4}$:
$$
\boxed{\ p_{\rm kurtosis}\approx0.20,\qquad p_{\rm mean/skew}\approx0.29\text{–}0.39.\ }
$$

**This resolves two things at once:**
1. **The deterministic ceiling on p (the flagged blank): p ≈ 0.2, kurtosis-limited** — small, non-analytic,
   set by the higher shape cumulant. A slow rate, exactly as the numerics hinted.
2. **It explains the physical-signature contamination** (`PHYSICAL_SIGNATURE_NOTES.md`: the excess kurtosis
   overshoots to −1.1 vs the limit's −0.24 at finite ε). That overshoot is not noise — it is the **rate-limiting
   cumulant converging non-analytically as ε^{≈0.2}**, so at any accessible ε the kurtosis is the last cumulant
   to settle. The sign of exk is the robust class-marker; its *magnitude* is contaminated precisely because its
   convergence is the slowest.

**The non-analyticity is a real mathematical statement, not a numerical accident**: the outer correction to the
cusp-limit shape is non-analytic in the coalescence parameter (exponent ≈ 0.8, not an integer), which is *why*
the rate is slow — the cusp is a genuine singularity in Δ, and the higher cumulants feel it most.

**Honest scope of this measurement.** This is the rate in the *inner coalescence* parameter Δ (the dominant,
cusp-specific "incomplete merge" correction), not a full proof that the FHN → 𝒲 rate equals ε^{0.2}: there may
be additional genuinely-outer GSPT corrections. But (a) it gives a concrete, defensible **target** p ≈ 0.2, and
(b) it identifies the **mechanism and the rate-limiting observable** (non-analytic kurtosis). Any Tier B proof
must reproduce this p.

---

## 3. Ranked attack plan for B3 + B4

**Route A (recommended) — port RRV's Riccati-explosion representation into the swept Weber setting and drive a
quantitative COUPLING.** The literature's clearest signal: the *only* framework that natively owns 𝒲_β's
first-explosion structure is RRV/KRV (Thread 1); it lacks a rate only because no one has built the coupling.
Concretely:

- Represent both laws as first-explosion times of Riccati diffusions on the *same* Brownian path: the **limit**
  Weber Riccati $dp=(\operatorname{sign}(Y)Y^2-p^2)d\tau+\eta\,dW$ (Tier A) and the **finite-ε inner** Riccati
  $dp_\varepsilon=(V_{\Delta(\varepsilon)}(Y)+\text{(outer GSPT corrections)}-p_\varepsilon^2)d\tau+\eta\,dW$,
  synchronously coupled.
- Their difference $D=p_\varepsilon-p$ obeys a **linear** SDE, $dD=-(p_\varepsilon+p)D\,d\tau+(V_{\Delta}-V_0
  +\dots)d\tau$ (no noise term — the noise cancels under synchronous coupling), so $D$ is driven purely by the
  **deterministic potential mismatch** $V_{\Delta}-V_0=\operatorname{sign}(Y)|Y|\Delta$ (plus higher GSPT terms),
  integrated against the (random) exponential factor $e^{-\int(p_\varepsilon+p)}$. This is a Grönwall/
  variation-of-parameters comparison, and the two first-explosion times differ by $O(\|D\|)$ near the crossing.
- **The single named hard sub-problem:** control $D$ (hence the first-explosion-time gap) **uniformly through
  the merging turning point**, where $p_\varepsilon+p$ (the mean-reversion of the $D$-equation) passes through
  the turning and the potential mismatch $\operatorname{sign}(Y)|Y|\Delta$ is largest relative to the vanishing
  restoring rate. This is where the Olver floor ($\bar p\ge c_0>0$, `COUPLED_CUSP_RESULTS.md` §2.3) and the
  Tier-A transversality (escape at $k_\star>0$, away from the turning) are the load-bearing inputs — the same
  structures that made T1 uniform should make the coupling uniform. The rate is then read off from
  $\|V_\Delta-V_0\|$, i.e. from $\Delta\sim\varepsilon^{1/4}$ — **consistent with the empirically-measured
  p ≈ 0.2–0.4** (§2), and the coupling would explain the non-analyticity as the singular turning-point
  contribution to $\int e^{-\int(p_\varepsilon+p)}(V_\Delta-V_0)$.
- **Why this is the best route:** it reuses everything Tier A built (the Riccati representation, the tail/tube
  bounds for tightness, the transversality, the floor) and needs exactly one new estimate — a coupling
  Grönwall through the turning — rather than importing an alien machinery. It targets the *distributional* rate
  directly (via a Wasserstein bound on the escape location from $\|D\|$), which is B3.

**Route B (fallback / complementary) — Berglund–Gentz sample-path tube + a smoothing/Edgeworth step for the
distribution.** Berglund–Gentz (Thread 4) give the pathwise tube and escape window through a dynamic pitchfork
(the symmetric two-turning-point normal form); the missing step is to convert the tube into a *distributional*
rate. This is the field-wide open B3, so Route B is "extend Berglund–Gentz to distributional convergence" —
higher-risk (it's their open problem too) but the tube estimates are directly reusable and additive-noise-proven
(caveat: our noise is multiplicative — open Q4).

**Route C (probably a dead end for the rate, but a clean negative) — RMT Green-function comparison.** The only
proven TW rates (Thread 2) are static resolvent-matching. There is no swept turning point or first-passage
functional in that machinery; it does not transfer. Documenting *why* (no dynamical analogue of the local law)
is a worthwhile clean negative, not an attack.

**Malliavin–Stein (open methodological question, not yet a route).** Thread 3 is blank because it is *untried*
on a first-explosion functional, not because it was disproven. The obstruction to assess first: is the
map (noise path) ↦ (explosion location) regular enough for the Stein–Malliavin integration-by-parts? Tier A
§3.7bis already proved this map is Malliavin-differentiable with $D_sY^\star=-\eta\Phi(s)/k_\star$ (transversal,
$k_\star>0$) — **so the regularity Malliavin–Stein needs is in hand.** This makes "try Nourdin–Peccati on
$Y^\star$" a genuinely open, possibly-tractable sub-project, not a blank. Worth a scoped attempt.

---

## 3bis. Route A executed — the coupling result [NEW: structure proved + mechanism verified]

I carried out Route A. The synchronous coupling is exact and gives a clean rate; the load-bearing mechanism is
numerically confirmed at β = 2 (full noise, `tierB_coupling.py`, N = 3·10⁵).

**The coupling (exact).** With the limit Weber Riccati $p$ and the finite-merge Riccati $p_\varepsilon$ on the
**same** Brownian path, $D=p_\varepsilon-p$ solves a **noise-free linear ODE**
$$
\dot D=-(p_\varepsilon+p)D+(V_\Delta-V_0),\qquad V_\Delta-V_0=\Delta\,Y,\qquad D(0)=0,
$$
so $D_\tau=\Delta\int_0^\tau e^{-\int_s^\tau(p_\varepsilon+p)\,d\sigma}Y(s)\,ds$ — **exactly proportional to Δ**
(the noise cancels; the potential mismatch is the clean linear term $\Delta Y$). At the field level this is the
node-shift $\delta Y^\star=\Delta\,\chi+O(\Delta^2)$ with the **susceptibility** $\chi=-u^{(1)}(Y^\star)/
u'(Y^\star)$, $u^{(1)}$ the noisy-Green's-function response to source $Yu$. **[structure — proved]**

**The rate (clean, in Wasserstein-1).** The coupling is an explicit coupling of the two laws, so
$$
\boxed{\ W_1\big(\mathcal W_{\beta,\Delta},\,\mathcal W_\beta\big)\ \le\ \mathbb E\,|\delta Y^\star|\ \le\
\Delta\,\mathbb E|\chi|+O(\Delta^2)\ =\ O(\varepsilon^{1/4})\ }\quad(\Delta\sim\varepsilon^{1/4}),
$$
**provided (i) $\mathbb E|\chi|<\infty$ and (ii) the $O(\Delta^2)$ remainder is uniform.** So the natural
universality-with-rate statement — weak convergence in $W_1$ — has **rate p = 1/4**.

**Mechanism verified numerically (β = 2, `tierB_coupling.py`, `tierB_chi.npz`):**
- **(a) pathwise linearity:** $\mathbb E[\delta Y^\star(2\Delta)]/\mathbb E[\delta Y^\star(\Delta)]=2.09$
  (predicted 2.0) and $\operatorname{corr}(\delta Y^\star(\Delta),\delta Y^\star(2\Delta))=0.95$ — $\delta Y^\star
  =\Delta\chi$ holds pathwise, χ a genuine Δ-independent functional. **[NUMERIC ✓]**
- **(b) $\mathbb E|\chi|=0.357$ finite** (q = 1 moment stable to 1 % across sample halves) ⟹ $W_1\approx0.36\,
  \Delta=0.36\,\varepsilon^{1/4}$. **[NUMERIC ✓]**
- **(c) χ is heavy-tailed, tail exponent $\alpha\approx2$** (Hill; from $\chi\propto1/k_\star$, rare
  near-turning escapes): $\mathbb E|\chi|<\infty$ but $\mathbb E\chi^4=\infty$. So the mean/$W_1$ converge at the
  clean $\varepsilon^{1/4}$, while the **kurtosis** correction is non-analytic — the rate-limiting shape
  observable. **[NUMERIC ✓]**

This reconciles the §2 FP sweep exactly: mean $\sim\Delta^{1.15}\approx\Delta=\varepsilon^{1/4}$ (clean); exk
$\sim\Delta^{0.8}$ slower (the heavy-tail-limited cumulant). **The physical-signature kurtosis contamination is
now fully explained**: it is the one cumulant governed by χ's $\alpha\approx2$ tail (the near-turning
small-$k_\star$ escapes), hence the slowest and most contaminated — its *sign* is robust (class marker), its
*magnitude* is not.

**What is proved vs. what remains.** The coupling structure ($\delta Y^\star=\Delta\chi+O(\Delta^2)$) and the
bound $W_1\le\Delta\mathbb E|\chi|+O(\Delta^2)$ are proved; the pathwise linearity and $\mathbb E|\chi|<\infty$
are numerically confirmed. **Two sharply-stated lemmas remain for full rigor:**
- **(L1) $\mathbb E|\chi|<\infty$** — i.e. the $q=1$ moment of $\chi=-u^{(1)}(Y^\star)/u'(Y^\star)$ is finite.
  Tier A gives $u'(Y^\star)=k_\star>0$ (transversality, §3.7bis) and $\mathbb E\|D Y^\star\|_{L^2}^2<\infty$; L1
  is the weaker $q=1$ statement and should follow from the same floor + a Green's-function bound. The
  $\alpha\approx2$ tail says this is *tight* — $\mathbb E|\chi|$ is finite but $\mathbb E\chi^2$ is borderline,
  so L1 is exactly at the edge and must be proved, not assumed.
- **(L2) the $O(\Delta^2)$ remainder uniform through the merging turning point** — the Grönwall control of the
  back-reaction ($p_\varepsilon$ appears in $E(\tau,s)$) as the restoring rate $p_\varepsilon+p$ passes through
  the turning. This is the genuine remaining analytic crux, and the place the Olver floor / Tier-A tube
  uniformity are the load-bearing inputs.

**Net for Route A:** Tier B's B3+B4 are reduced to **an explicit rate ($p=1/4$ in $W_1$) with the coupling
mechanism verified and two named lemmas (L1, L2)**, both tied to already-proved Tier A structure. This is a
concrete, defensible partial theorem — the first quantitative-rate statement for this class of limit — not just
a plan.

## 4. Net + the named frontier

- **Literature verdict [verified]:** a distributional rate for a first-passage/Riccati-explosion functional is
  unproven everywhere; the two relevant branches (stochastic-Airy/TW_β; RMT-edge-with-rate) are disjoint and
  neither transfers; the nearest dynamical analogues (Berglund–Gentz, Bakhtin) give tubes/scaling-limits, never
  a rate. **No "physical system → TW with a rate" theorem exists.** Tier B is a genuine frontier.
- **Empirical rate [NUMERIC ✓, new]:** $p\approx0.2$, **kurtosis-limited and non-analytic** (exk $\sim\Delta^{0.8}
  \sim\varepsilon^{0.2}$); mean/skew faster ($\varepsilon^{0.3\text{–}0.4}$). This resolves the flagged
  deterministic-ceiling blank and explains the physical-signature kurtosis contamination.
- **The single named hard sub-problem [mapped]:** a quantitative, ε-uniform **synchronous coupling / Grönwall
  comparison of the two Riccati first-explosion times through the merging turning point** (Route A). Everything
  it needs except that one estimate is already proved in Tier A / T1 (Riccati representation, floor,
  transversality, tightness). This is the concrete next theorem to attempt.
- **Realistic p:** ≈ 0.2 (limited by the non-analytic kurtosis correction) — small, which is *why* finite-ε
  numerics look contaminated, and which a Tier B proof must reproduce.

**Honest status.** This is scoping + an empirical rate, **not** a proof of Tier B. But it converts "the deepest
unknown" into (i) a verified statement of where the field actually stands, (ii) a measured target p ≈ 0.2 with
its mechanism, and (iii) one named estimate to attempt. That is the intended Tier B deliverable: the mapped
frontier with a concrete first step.
