# MVP road: inverting the edge theory into a proximity-to-bifurcation diagnostic

*The moonshot is to read **how close a system is to a bifurcation, and which one**, off real data —
a bifurcation-specific early-warning signal that uses the full edge law, not just rising variance.
This note scopes the minimum viable product and a small-scale road to start inverting on convenient
datasets. Rung 1 (atlas + inverter + synthetic ground truth) is **built**:
`bifurcation_inverter.py`, figure `figures/bifurcation_inverter.png`.*

---

## The idea that makes it tractable: the σ-free shape map

The forward edge laws give, for each bifurcation class, a one-parameter family of interval/amplitude
distributions indexed by the **rescaled distance to threshold** `ν` (for the SNIC, `ν=I/σ^{4/3}`).
The mean interval mixes in `σ` and the unknown biophysical time-scale, but the **dimensionless
shape** — CV, skew, the normalised density — is a **σ-free, monotone** function of `ν`:
```
        CV(ν)  monotone  ⇒  shape statistics pin ν  with no knowledge of σ or the time unit.
```
So the identifiable, useful output is **`ν` (and its trend)**; the early-warning signal is **`ν → 0`**.
Absolute `σ` needs a calibration sweep and is *not* part of the MVP.

## Rung 1 — DONE: SNIC atlas + inverter + synthetic ground truth

`bifurcation_inverter.py`:
* **atlas** `ν → (CV, skew)` from the universal inner FPT `dρ=(ν+ρ²)ds+dW` (`fpt_sim`), cached.
* **invert** ISIs → `ν̂` by matching `(CV*, skew*)` to the atlas, with a bootstrap CI and a **fit
  residual** (a SNIC class-check).
* **validated** on physical QIF spike trains at known `ν`:

```
   white σ=0.4 :  true ν −0.30/0.00/0.50/1.00/1.80  ->  ν̂ −0.22/+0.02/0.52/1.07/1.99
   white σ=0.7 :  true ν −0.30/0.00/0.50/1.00/1.80  ->  ν̂ −0.22/−0.03/0.56/1.04/1.98   (σ-invariant!)
   adapt σ=0.4 :  ν̂ biased HIGH (0.06/0.14/0.66/1.56/1.67) AND residual 0.3–0.5 vs <0.2 (flagged)
```

Three things proven: recovery is accurate for the clean model; **`ν̂` is σ-invariant** (the core
identifiability claim — σ=0.4 and 0.7 give the same answer); and the dominant confound (adaptation)
biases `ν̂` in a known direction *and* is caught by the residual. (The slight high-`ν` overshoot is
finite-size: physical `v_th=14` vs the canonical `v_th=∞` atlas — removable by matching the cap.)

## The road (small-scale, convenient datasets)

### Phase 1 — harden the inverter & make classification real — **DONE** (`bifurcation_classifier.py`)
- 2-moment match replaced by **full normalised-density** matching via the **1-Wasserstein** distance
  on quantile functions; the W1 at the best-fit parameter **is** the class-check.
- **Two-class atlas** {SNIC (Type-I), Hopf (Type-II = Stuart-Landau rotation intervals)}; `classify`
  returns `(class, param̂, W1)` = argmin over all `(class, param)`.
- **Result.** At **matched CV** (SNIC `0.282` vs Hopf `0.291`) the densities differ sharply
  (skew `1.26` vs `0.39`): CV alone cannot separate them, **full-density matching classifies 32/32
  test series correctly (100%)** and recovers `ν` within the SNIC class. Figure
  `figures/bifurcation_classifier.png`.
- **Phase 1+ — DONE** (`bifurcation_classifier.py phase1plus`):
  - **Third class — saddle-homoclinic** (Izhikevich's 3rd canonical spiking onset): the noisy
    saddle-passage ISI `T_ret + (1/λ)ln(δ/|μ+ση|)` — a hard lower edge + logarithmic/exponential
    right tail. **Three-class confusion = 89%**, errors confined to the **SNIC↔homoclinic** pair
    (both Type-I / class-1 — genuinely similar; **Hopf separated perfectly, 12/12**). That the two
    Type-I onsets are the hard pair is the *correct* answer, not a bug.
  - **Confound axis — adaptation.** A 2D **(ν, b)** SNIC atlas; the joint fit **recovers the
    adaptation strength** (`b̂=0.25`) and **removes the bias for ν≳0.8**, with a residual low-ν
    `(ν,b)` degeneracy (low-ν/high-adapt ≈ higher-ν/low-adapt) — i.e. the method *flags* where
    distance and adaptation are not separable. `figures/bifurcation_classifier_phase1plus.png`.
  - **Hardened — DONE** (`harden_inverter.py`, `INVERSION_HARDENING.md`,
    `figures/harden_inverter.png`): **off-grid (interpolated) fitting** kills grid-snap (between-node
    SNIC `ν` recovery error 0.079→0.031, ≈60%); a **colored-noise (ν, τ_s) axis** joins the (ν, b)
    adaptation axis (joint fit cuts colored-data `ν` error 0.092→0.052, ≈43%, and flags coloring —
    though `τ_s` is only weakly identified, `W1` valley ≈8× steeper in `ν`); **near-threshold
    weighting** for the class call (`aggregate_class`) leans on near-edge sweeps where SNIC↔Hopf
    separability lives (`W1`-to-Hopf collapses 0.066→0.017 across drive), improving the cell-level
    posterior (high-drive sampling P(Type-I) 0.56→0.63) — honest limit: it sharpens *calibration*, not
    the binary call rate on cells already easy or sampled only far above threshold.

### Phase 2 — real spike data: Allen Cell Types — pipeline BUILT + dry-run validated
`allen_phase2_inversion.py` (modes: `allen` = real local pull; `dryrun` = in-sandbox surrogate),
per Long-Square sweep → `(class via classify3, confound-corrected ν̂ via fit_adapt)`.
- **Dry-run (Allen-like surrogate cells, realistic ~26 ISIs/sweep):** `ν̂` rises **monotonically
  with injected current** for the low-adaptation cells (3/4), tracking across cells of different `σ`
  (H1/H3 supported); **75% of sweeps classified Type-I** (H2) — the 25% Hopf calls concentrate at
  **high drive**, where Type-I firing becomes very regular and its ISI shape converges with Hopf's
  (honest limitation: **the class signature is strongest NEAR threshold**). Small samples ⇒ wide
  bootstrap CIs + grid-snap coarseness. `figures/allen_phase2_dryrun.png`.
- **Real run (local):** `pip install allensdk && python3 allen_phase2_inversion.py allen` — same code path.
- **Lessons for real data:** weight near-threshold sweeps for the class call; add off-grid
  interpolation; strong adaptation can still break monotonicity (wants the colored-noise/finer
  confound axis). *Dataset: Allen Cell Types (open, NWB).* — **all three now implemented in
  `harden_inverter.py`** (`classify3_interp`/`fit_adapt_interp` off-grid, `aggregate_class`
  near-threshold class call, `fit_colored` colored-noise axis); ready to wire into the `allen` path
  via `invert_sweep`.

### Phase 3 — simulated tipping: the early-warning demo — **DONE (Type-I/SNIC)**
`phase3_tipping.py` + `phase3_figure.py` + `PHASE3_TIPPING.md` + `figures/phase3_tipping.png`.
A SNIC-QIF with a slowly-drifting current `I(t)→I_c`; ISIs streamed, `ν̂` inverted in sliding
fixed-count windows, benchmarked against rising ISI variance and lag-1 autocorrelation at a
**matched 15% false-alarm rate** (each trend alarm calibrated on a stationary-null ensemble).
- **Honest verdict (120-trial ensembles).** On a *clean* approach the `ν̂` trend and the ISI-variance
  trend are **comparable** (detection 50% vs 56%, median lead 792 vs 836) — the edge-law *trend* does
  not beat rising variance on confound-free data, and we say so. **lag-1 AC of ISIs is a poor EWS**
  (11% ≈ chance): the QIF train is ~renewal, so AC1 — built for the state variable's slowing — has
  little to track.
- **Where `ν̂` wins (the value proposition).** (i) Under a **rate-drift confound** (slow rate
  rescaling, no approach) the **variance** false-alarm rate nearly doubles (15%→**31%**) while `ν̂`,
  built from the σ-free ISI *shape*, stays at baseline (**11%**). (ii) The **absolute** `ν̂<0.4` alarm
  — available only to a calibrated distance with a real zero — gives **96% detection at 0% false
  alarms** on both nulls; Var/AC1 have no absolute scale and cannot offer one.
- *Dataset: simulated (Epileptor is a published, public model).*

### Phase 3b — transfer to the **Epileptor** — **DONE (a domain-of-validity result)**
`epileptor_phase3b.py` + `epileptor_figure.py` + `PHASE3B_EPILEPTOR.md` +
`figures/epileptor_phase3b.png`. Ported the pipeline to the canonical Epileptor (Jirsa et al. 2014,
`x0=-1.6`); the slow `z` is an *internal* drift through SNIC-onset / homoclinic-offset.
- **Honest verdict.** The discharge train is spike-and-wave (20% of raw ISIs `<1`); at the
  **complex** (dominant-rhythm) level it is a clean **Type-I** shape (classifier residual
  `W1=0.014`, in-distribution — the inversion is well-posed, not garbage). **But `ν̂` stays high and
  flat (~2.9) across the seizure and never warns**: the offset is *not* a period-diverging
  bifurcation of the discharge cycle (inter-complex period ratio `last5/median = 0.86`, i.e.
  discharges **accelerate** into an abrupt multi-variable termination). So the `ν̂→0` warning of 3a
  **does not transfer**.
- **The contribution = a sharp domain-of-validity boundary.** The edge-law early-warning needs an
  **adiabatic** crossing — the system must dwell near threshold long enough for the universal
  first-passage shape (CV→0.57, diverging period) to form. 3a's ramped SNIC qualifies; the
  Epileptor's fast `z`-driven offset does not. Constructively, the method reads a developed seizure
  as a *robust far-from-threshold* Type-I oscillation (correct) and **does not false-alarm** — its
  natural target is the noise-driven **pre-ictal/onset** regime, not seizure termination.
### Phase 3b-cont — taxonomy sweep — **DONE (the early-warning = the two-edge split)**
`taxonomy_sweep.py` + `taxonomy_figure.py` + `TAXONOMY_SWEEP.md` + `figures/taxonomy_sweep.png`.
Drifted the bifurcation parameter slowly through **four** canonical onset/offset classes and mapped,
per class, period-divergence / timing-`ν̂` / amplitude-collapse:
- **Phase edges (SNIC, saddle-homoclinic):** period **diverges** (2.2×) → **timing** warning
  (SNIC `ν̂→0`, min `−0.11`; SH via period divergence, log-tailed shape).
- **Amplitude edges (Hopf, fold-of-cycles):** period **finite** → `ν̂` **blind**; **amplitude**
  collapses. Hopf declines *gradually* (precursor → 0.34 before crossing); **fold-of-cycles** holds
  then **jumps** — *no precursor in either channel* (the catastrophic, warning-free tipping).
- **Punchline.** The early-warning observable splits along exactly the paper's two-channel line
  (phase = Type-I/`σ^{2/3}`/quartic-FPT → `ν̂`; amplitude = Type-II/`σ^{1/2}`/canard–TW → amplitude).
  The inversion's **class call selects which observable can warn** — turning the steady-state
  two-edge universality into an operational early-warning prescription.
- *Still open: a fully-dynamical SH oscillator (BT path-following); an amplitude-edge `ν̂`-analogue
  inversion + per-class matched-FPR benchmark; a (drift-rate × class) adiabaticity map; and a
  noise-driven excitable Epileptor positive control.*

### Applied reframe — excitability characterization (the genuinely-useful comp-neuro direction) — **DONE (synthetic)**
`excitability_characterizer.py` + `excitability_figure.py` + `allen_excitability.py` (local) +
`EXCITABILITY_CHARACTERIZER.md` + `figures/excitability_characterizer.png`. Instead of "early-warning
of tipping" (oversold, hard to validate), point the same machinery at a question comp neuro needs and
usually requires *intracellular* access to answer: **is this neuron Type-I (SNIC) or Type-II (Hopf)
excitable, and how far above threshold?** — from extracellular spikes alone.
- **Method = three fused features:** near-threshold-weighted ISI-shape class; **f–I onset** (Hodgkin
  Class 1 vs 2: Type-I reaches low rate, Type-II bounded); and the **`W1` residual** as a goodness-of-fit
  that **rejects out-of-family** cells. Plus confound-corrected `ν̂` across drive and an adaptation flag.
- **[N] Synthetic-cohort validation (18 cells, known truth): 100% three-way accuracy** (Type-I /
  Type-II / out-of-family), after fixing the adaptation confound by making the **f–I onset** the
  arbiter (adaptation flips the ISI-shape vote but not the onset bifurcation; clean onset gap 0.47 vs
  0.87). Residual rejects out-of-family (0.30 vs ≤0.03, ~10×). `ν̂` recovery `|ν̂−ν|≈0.43` (coarse
  absolute scale, monotone trend). **But this is an IN-FAMILY number.**
- **[N] Conductance-model reality check — `morris_lecar_check.py` + `figures/morris_lecar_check.png`.**
  Running the identical `characterize_cell` on **Morris–Lecar** Type-I/II cells: **all auto-rejected**
  (residual `0.091` ≈ garbage `0.108`, ≫ in-family `0.007`). ML *leaves* the noise-driven first-passage
  regime at moderate drive (skew flips negative, regular firing) → outside the normal-form atlas. The
  self-diagnosing residual makes the gap **visible** (refuses rather than mislabels) — but the headline:
  **the normal-form atlas does not yet transfer to a conductance model; external validity is unproven.**
- **Why it matters:** converts the moonshot from "predict tipping" into "**measure a fundamental
  dynamical invariant cheaply, with a method that refuses cells it can't characterise**" — mechanistic,
  not phenomenological. The reframe is sound; the open problem is the atlas's *family*.
- **[N] Kill-or-confirm separability test — `ml_atlas_separability.py` + `figures/ml_atlas_separability.png`.**
  Train the atlas on Morris–Lecar itself; classify held-out ML cells by **shape alone**. **VERDICT:
  PARTIAL** — Type-I/II *survives* in conductance-model ISIs but weakly: shape-only **69%** all-drive
  (vs 89–100% for normal forms), **78% near rheobase**; separation collapses at high drive (Type-II
  → 25–38% in the regular-firing regime). f–I onset still directional (0.32 vs 0.56). **Not killed,
  re-spec'd:** realistic ceiling ≈ **80%** (biophysical atlas + near-rheobase sweeps + shape & f–I
  onset) — *population cell-type characterisation, not per-cell certainty.*
- *Still open (the critical path, now scoped): build the biophysical atlas near rheobase; combine shape
  + f–I onset; then the real Allen run; real-data reject-threshold calibration. Manage expectations to
  ~80%, not the in-family 100%.*
- *Parallel theory track (queued): noise/fluctuation corrections to next-generation neural-mass models
  (Montbrió–Pazó–Roxin / Coombes–Byrne) near their SNIC/Hopf transitions — where the two-edge
  universality is novel to that community rather than competing.*

### Theory track — the two edges at the COLLECTIVE level — **STARTED (opening probe)**
`neural_mass_edge.py` + `neural_mass_figure.py` + `NEURAL_MASS_EDGE.md` + `figures/neural_mass_edge.png`.
The Montbrió–Pazó–Roxin exact QIF mean-field **+ adaptation** is a slow–fast macroscopic system — a
network-level folded cycle. **[N]** It has a **macroscopic phase edge**: the collective-burst period
**diverges** as `η̄ → η̄* ≈ −2.32` (network SNIC/homoclinic). **[N]** Just below it, noise-induced
**collective inter-burst intervals** are the universal first-passage shape — positive-skew, classified
**Type-I/SNIC by the same atlas and inversion as single-neuron ISIs**, CV → critical as `η̄ → η̄*`. So the
single neuron and the whole network sit at the **same universality edge**.
- *This is the novel wedge into the (deterministic) next-gen neural-mass community: a principled account
  of fluctuation universality at the macroscopic bifurcations.*
- *Honest gaps (opening probe, not a theorem): noise is phenomenological; burst extraction residual;
  one operating point.*

### Theory track — finite-size Langevin derivation — **DONE (derivation + validation + a correction)**
`mpr_finite_size.py` + `MPR_FINITE_SIZE_LANGEVIN.md` + `figures/mpr_finite_size.png`.
- **[R/structure]** System-size expansion of the QIF network ⇒ an **$O(N^{-1/2})$ multiplicative
  macroscopic Langevin** on MPR's $(r,v)$ (explicit covariance $\Sigma$ left as the **[H]** step).
- **[N] Direct validation:** finite-$N$ theta-network rate variance $\mathrm{Var}(r)\propto 1/N$
  (**slope $-0.98$**, $\mathrm{Var}\cdot N$ const over $N=250$–$4000$) — the $N^{-1/2}$ noise confirmed.
- **[H] Edge reduction:** at a macroscopic saddle-node, projection ⇒ the **canonical noisy saddle-node**
  (paper's phase channel) with $\sigma_{\rm eff}\propto N^{-1/2}$ ⇒ (if SNIC) collective timescale
  $\propto N^{-1/3}$, critical window $\propto N^{-2/3}$, CV$\to$0.57.
- **[N] Honest correction:** at the studied operating point the macroscopic edge is a **fold of cycles
  (amplitude edge) + coexisting rest**, *not* a SNIC (period finite, amplitude steady, abrupt; the
  collective $\sigma^{2/3}$ test fails, slope $-0.18$) — so the bursts are **rate-governed bistable
  hopping**. This **supersedes the phase-edge claim** in `NEURAL_MASS_EDGE.md`.
- *Program: locate a macroscopic SNIC and test $N^{-1/3}$; compute $\Sigma$ (→ theorem); collective
  Tracy–Widom at the fold (the amplitude edge the network actually has).*

### Theory track — the two follow-ups — **DONE (both honest, both point to one analytic next step)**
`MPR_EDGE_FOLLOWUPS.md`.
- **[N] Macroscopic SNIC search: NOT FOUND.** Across $\alpha\in\{3,5,8\}$, $\tau_a\in\{15,40\}$,
  $\Delta\in\{0.3,1\}$, the collective bursting **always ends at a fold of cycles** (period finite
  $17$–$21$, amplitude steady $\sim2.1$, abrupt; resolved to 4 digits — no log-divergence). The network
  **realises the amplitude edge**, not the phase edge ⇒ the $N^{-1/3}$ phase-law test is moot here.
- **[N] Collective Tracy–Widom: SUGGESTIVE, NOT CONFIRMED.** Near-fold finite-size fluctuations carry a
  **TW-signed positive skew** ($+0.36$ vs $\mathrm{TW}_1$ $+0.29$) but are small/Gaussian-dominated (not
  edge-regime); the ramp peel-off is dominated by dynamic-bifurcation delay.
- **[H] The decisive next step is analytic, not numerical:** center-manifold-reduce the noisy MPR fold
  of cycles and check whether the amplitude dynamics is the **Riccati/Airy canard** (paper's amplitude
  channel). If yes ⇒ **collective $\mathrm{TW}_\beta$** with $1/\sqrt N$ as intensity (the headline
  channel at the network level); if no ⇒ a different collective universality. A pen-and-paper reduction.

### Theory track — the canard reduction — **DONE (structurally yes; collective TW)**
`MPR_CANARD_REDUCTION.md`.
- **Key fix:** the amplitude-edge canard lives at the **fold of the critical manifold** (the saddle-node
  of fast $(r,v)$ equilibria every burst jumps through), **not** the fold of cycles — so the right
  observable is the **per-jump peel-off level**, which the follow-up numerics had mis-measured.
- **[R/N]** Critical manifold $\mu(r)=\pi^2r^2-Jr-\Delta^2/4\pi^2r^2$; upper fold $(r_f,\mu_f)=(0.754,-5.744)$,
  $\mu''>0$; the bursting trajectory **canards past it** ($\mu\to-5.785$, overshoot verified).
- **[R structure]** Fold reduction + adaptation's **linear** slow passage $\Rightarrow$ the **time-linear
  Riccati** $d\rho=(A\dot\mu(t-t_f)+B\rho^2)dt+\sigma_\rho dW$ (paper's amplitude channel); Cole–Hopf
  $\Rightarrow$ **stochastic Airy operator** $\Rightarrow$ collective $\mathrm{TW}_\beta$, with
  $\sigma_\rho\propto N^{-1/2}\Rightarrow\beta\propto N$. **The network inherits the headline edge.**
- **[N] Corroboration:** peel-off skew $=+0.29$ (*exactly* $\mathrm{TW}_1$) at $\tau_a=15$, but $\approx0$
  at $\tau_a=40$ — TW-signed and parameter-sensitive; corroboration, not confirmation.
- **[H] Open (now well-posed):** compute $\Sigma$ to fix $\beta(N)$ (→ theorem); numerically test
  $\beta\propto N$ + full-shape $\mathrm{TW}_\beta$ fit of the peel-off distribution at a sharper-fold
  regime.

### Theory track — closing the two gaps — **DONE (partial: structure stands, constants nearly)**
`MPR_SIGMA_AND_BETA.md`.
- **[R + leading constant] $\Sigma$ in closed form.** Molecular chaos $\langle|\delta Z|^2\rangle=(1-|Z|^2)/N$
  + Montbrió's map ($\partial_Z r=-1/\pi(Z+1)^2$) ⇒ $\mathrm{Var}(r)\,N=2(1-|Z|^2)/\pi^2|Z+1|^4$.
  Predicts $0.081$ vs measured $0.029$: **scaling exact, constant within $\sim2.8\times$** — the residual
  is the coupling (self-consistent field feedback) renormalisation, now a *single O(1) factor*, not open
  structure. Fixes $\beta=cN$ to leading order.
- **[N] β∝N test — corroborated, not confirmed.** Peel-off skew $=+0.19$ (TW-signed, TW-magnitude,
  between $\mathrm{TW}_2$/$\mathrm{TW}_4$) but **roughly constant in $\sigma$** (not $\beta\propto\sigma^{-2}$),
  width $\sim\sigma^{0.8}$ (not edge $\sigma^{2/3}$): **not in the asymptotic edge regime** at these
  params ($\varepsilon=1/\tau_a\approx0.067$ too large; max-over-window estimator).
- **The structural claim stands:** the QIF network inherits the paper's Tracy–Widom edge, $1/\sqrt N$ the
  intensity. **Remaining = constants + a cleaner numerical confirmation, not new structure:** (1) the one
  self-consistent renormalisation of $\Sigma$; (2) β∝N / full-shape $\mathrm{TW}_\beta$ fit at a
  sharp-canard, finite-$N$ regime with a per-passage peel-off estimator.

#### Both gaps then worked further — diagnosed as frontier constants
- **Step 1 sharpened:** ruled OUT the field-feedback guess (the $(r,v)$ Lyapunov balance returns the bare
  $0.084$; coupling *slows* relaxation, ratio $0.76$) and the sampling protocol (random-i.i.d. $\eta$
  gives $0.035$, still $\sim2.3\times$ under $0.081$). **Mechanism = sub-Poissonian phase correlations**
  (the deterministic coupled flow anti-correlates the phases, variance $\sim2.4\times$ below i.i.d.) — a
  **pair-correlation / 2nd-circular-cumulant** constant.
- **Step 2 done on the direct finite-$N$ network:** per-burst peel-off skew $+0.21$ (TW-magnitude) but
  **identical at $N=400$ and $1600$** — dominated by *deterministic* canard variability; the validated
  $N^{-1/2}$ noise is sub-dominant, so **β∝N is not observable** without a stiff **sharp-canard
  ($\tau_a\to\infty$)** regime.
- **Net:** structural collective-TW result intact; both residuals are genuine frontier calculations
  (sub-Poissonian constant; stiff sharp-canard β∝N sim), each now precisely diagnosed.

### Phase 4 — beyond MVP
- Real EEG/LFP (seizure onset), ECG, ecological/climate tipping series (the public datasets used in
  the critical-transitions literature); multivariate; nonstationary inference with uncertainty.
- Fill the **atlas** toward the full hierarchy (banked #6–#10: cusp/Pearcey at the Bogdanov–Takens
  point, fold-of-cycles/TW, Lévy/Fréchet) and the **classification theorem** (catastrophe type →
  β-ensemble edge) that puts the atlas on rigorous footing.

## Honest obstacles (why it's a program, not a script)
1. **Confounds bend the map** — adaptation and colored noise shift `ν̂` (quantified in
   `QIF_DATA_VALIDATION.md`); Phase 1 must either jointly fit them or use the residual to gate.
2. **Nonstationarity** — real tipping is a moving target; needs windowed inference + change-point care.
3. **σ / time-scale not identifiable** from one record — report `ν` and its *trend*, not absolute σ.
4. **Forward laws still partial** — the quartic FPT has no closed form; the atlas is simulation-built
   (fine for an MVP, but the rigorous version wants the analytic densities).

## Reproduce
```
python3 bifurcation_inverter.py atlas      # build + cache the SNIC atlas (one-time, ~6 s)
python3 bifurcation_inverter.py validate   # -> figures/bifurcation_inverter.png
```

## Bottom line
Rung 1 shows the inversion is real and σ-free on ground-truth data. The shortest path to something
**people can use** is Phases 1–3: full-density matching + a 2–3-entry atlas (Phase 1), confirm on
open Allen spike data (Phase 2), then demonstrate earlier-than-variance warning on a simulated
tipping system (Phase 3). That trio is an MVP a clinician-facing collaborator could evaluate.
