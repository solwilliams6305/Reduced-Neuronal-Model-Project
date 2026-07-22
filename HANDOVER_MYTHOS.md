# Handover — Reduced Neuronal Model Project (for Mythos)

You are picking up an ongoing research collaboration with **Solomon** (applied mathematician,
Edinburgh). The project orbits one paper — **"The noisy folded limit cycle: Tracy–Widom statistics
of canard escape"** — and has since grown a neuroscience-facing arm and a "moonshot" inversion tool.
Work is honest, rigorous-where-possible numerics + analysis, every claim tagged
**[R]** proved / **[N]** numerical / **[H]** heuristic / **[refuted]** / **[retracted]**. Read this
top to bottom, then skim `NEXT_DIRECTIONS_BANKED.md` and `INVERSION_MVP_ROADMAP.md` — they are the
two live index docs.

## The core science (the dictionary you must hold)

A noisy folded limit cycle has **two distinct universality edges** at its critical boundaries:

- **Amplitude edge** = the canard **peel-off level** (how far the escape jumps). It is the **ground
  state of the stochastic Airy operator** ⇒ **Tracy–Widom `TW_β`**; onset noise `σ_*∼√ε₂`; anchored
  by the first Airy zero `−2.338`. This is the paper's headline.
- **Phase edge** = the spike **timing / inter-spike interval**. Its generator is the **quartic
  (anharmonic) FPT operator** — non-integrable, so **provably NOT Tracy–Widom**; scaling `σ^{2/3}`;
  the mean-interval constant is a closed form `J=(√π/3)·12^{1/6}·Γ(1/6)≈4.976` (a `Γ(1/6)`, the
  phase analogue of `−2.338`).

Key reductions: the **QIF / θ-neuron = SNIC normal form = off-critical noisy saddle-node**
`dR=(μ+R²)dT+σdB`, with rescaled detuning `ν = I/σ^{4/3}` (`I_c=0`). The folded **node** inner
equation is Weber/parabolic-cylinder → SAO → `TW_β`. The two edges map onto the two excitability
classes: **Type-I (SNIC)** carries the phase edge (`σ^{2/3}`, quartic ISI); **Type-II (Hopf) /
fold-of-cycles burster** carries the amplitude edge (`σ^{1/2}`/`√ε₂`, TW level law).

## Working environment & conventions (read before running anything)

- **Save outputs** to `/Users/solomonwilliams/Reduced Neuronal Model Project/regime-tests/` (scripts
  + notes) and `…/figures/` (PNGs). Present files to Solomon with the file-sharing tool; keep
  post-ambles short.
- **Sandbox bash** has `numpy / pandas / matplotlib` but **NOT `scipy`, `h5py`, `allensdk`**. Use
  `pip install … --break-system-packages` if needed. **45 s timeout per call** — split long sims,
  cache atlases to `.npz`, run one heavy unit per call.
- **No programmatic web/data fetching** (policy). Anything needing the Allen API or other downloads
  is delivered as a **local-run script**, not executed in-sandbox; validate instead on synthetic
  ground truth.
- **matplotlib mathtext**: always brace — `\frac{1}{2}`, `\sqrt{\beta}`, `\sqrt{y}` (bare `\frac12`
  etc. crash). This has bitten ~8 times.
- **Hard-won numerical lessons:** (1) Cole–Hopf canard integration uses the **recessive `+√Y`
  (repelling) branch**, integrated downward (self-correcting); the attracting branch gives the wrong
  anchor. (2) For tridiagonal eigenvalues trust **exact Sturm-bisection**, NOT inverse iteration
  (which once produced a false "overturn" that had to be retracted). (3) `J_of_nu(0)=6.29 =
  2^{1/3}·4.976` — the `D=1/2` vs `D=1` diffusion convention; keep it straight.

## What exists (file map, all in `regime-tests/` unless noted)

**Paper-core (done / paper-ready):**
- `frisch_lloyd_current.py`, `close_H_matching.py` (`H=K·η²/4π`, `K≈1.039`), `tail_reconciliation.py`
  — Direction A prefactor (~80%; analytic Eyring–Kramers caustic still open = banked #2).
- `phase_edge_snic.py`, `offcritical_phase_edge.py` — phase edge + off-critical crossover. **NB:
  `offcritical_phase_edge.py` exports `J_of_nu`, `quartic_E0`, `fpt_sim` and is imported everywhere.**
  `FOLDED_CYCLE_PHASE_EDGE.md`, `TWO_CHANNEL_SYNTHESIS.md`.
- `folded_node_universality.py` — Weber→SAO→`TW_β`. `FOLDED_NODE_UNIVERSALITY.md`,
  `FOLDED_NODE_EDGE_LEMMA.md`.
- Direction C arc (`airy_*`, `eigenvalue_vs_node.py`, `gap_roughness.py`, `coupled_array_peeloff.py`,
  `airy_line_ensemble.py`): **verdict** — successive peel-offs of ONE cycle = Airy₂ point process
  (✓ to ~1%); coupling separate cycles does NOT give the Airy₂ *process* (local/exponential).
  `FOLDED_CYCLE_AIRY_PROCESS.md`.

**Neuroscience arm (done this session):**
- `qif_firing_statistics.py` + `QIF_PHASE_EDGE.md` — #1: noisy QIF = phase edge; `r=σ^{2/3}Φ(ν)`,
  quartic ISI, CV crossover.
- `qif_validation.py`, `allen_qif_validation.py` (local) + `QIF_DATA_VALIDATION.md` — confound
  fingerprints: colored noise pushes threshold CV **>0.57**, adaptation **<0.57**.
- `coherence_resonance_snic.py` + `COHERENCE_RESONANCE_SNIC.md` — **σ\*∝|I−I_c|^{3/4}** (slope 0.78→
  0.75); Type-I regularises to plateau, no sharp CR minimum.
- `two_edge_classifier.py` + `TWO_EDGE_CLASSIFIER.md` — exponent classifier: Type-I `σ^{2/3}` (timing)
  vs Type-II `σ^{1/2}` (amplitude), fit to <1.5%.

**The moonshot — proximity-to-bifurcation inversion (the live frontier):**
- `bifurcation_inverter.py` — rung 1: SNIC atlas (`ν→CV,skew` via `fpt_sim`) + inverter; **ν̂ is
  σ-free** (σ=0.4 and 0.7 recover identical ν̂).
- `bifurcation_classifier.py` — **Phase 1 / 1+**: full-density **1-Wasserstein** matching over a
  3-class atlas {SNIC, Hopf, saddle-Homoclinic} + a 2D **(ν, b)** adaptation axis. Entry points:
  `atlas / atlas3 / adaptatlas / phase1 / phase1plus`. Results: 2-class 100%, 3-class **89%** (errors
  confined to the genuinely-similar SNIC↔homoclinic Type-I pair; Hopf perfect), joint fit recovers
  adaptation `b̂` and corrects bias for `ν≳0.8` (low-ν `(ν,b)` degeneracy remains).
- `allen_phase2_inversion.py` — **Phase 2**: per-sweep inversion. `dryrun` (in-sandbox surrogate) and
  `allen` (real, local). Dry-run: ν̂ tracks current monotonically (3/4 cells), **75% Type-I**; class
  signal degrades at high drive (regular firing → looks Hopf). `INVERSION_MVP_ROADMAP.md` is the
  master plan.

The two index docs: **`NEXT_DIRECTIONS_BANKED.md`** (status of A–D + #1, plus banked #2–#10:
prefactor caustic, stochastic blow-up, true-Airy-PDE, two-edge observable, **Pearcey/Bogdanov–Takens**,
**BBP transition**, **Painlevé II tails**, **MMO line ensemble**, **Lévy/Fréchet**) and
**`INVERSION_MVP_ROADMAP.md`** (Phases 0–4; 0/1/1+/2 done).

## Where to go next (pick with Solomon)

1. **Run the real Allen pull** — Solomon runs `pip install allensdk && python3
   allen_phase2_inversion.py allen` locally; you read the actual `ν̂`-vs-current + class output together.
2. **Harden the inverter** — near-threshold weighting for the class call (it degrades at high drive),
   off-grid/interpolated fitting (kill grid-snap), a colored-noise confound axis.
3. **Phase 3 (the payoff demo)** — a simulated tipping system (public **Epileptor**, or a normal form
   with a drifting parameter): show `ν̂→0` *before* the transition and **benchmark against variance /
   lag-1 autocorrelation** early-warning signals. Beating those is the whole value proposition.
4. **Paper-side** — fold the QIF/Type-I material into the paper's discussion; or develop a banked
   direction (#8 Painlevé tails is lowest-risk and strengthens the main paper; #6 Pearcey/BT is the
   most novel).

## How Solomon likes to work

Concise and direct; honest about limitations (he values a clean "here's what fails" over hype);
build-and-verify (every result gets a numeric check or a figure he can see); tag everything; don't
silently paper over a confound or a degeneracy — surface it. Start by confirming which of the next
steps he wants, then go.
