# Coupled Reduced Neuron Model — Theory & Literature Map

**Project:** Coupled reduced neuron model — isolating / mapping new regimes (synchrony, anti-phase, MMO, bursting, spike-adding).
**Purpose of this document:** Establish what is already known (so we don't re-derive it), surface the methodological traps before we write code, and pin down the *genuinely new territory* the atlas should target — including the angle that ties this project to the noisy-folded-cycle / Tracy–Widom paper.
**Status:** Step 1 (theory/literature map). Drafted June 2026. Citations independently verified except where flagged ⚠.

---

## 0. The one-paragraph version

Mapping the regimes of two electrically coupled reduced neurons — synchrony, anti-phase, bistability, emergent bursting, coupling-induced MMOs, spike-adding — is **largely a solved/active area**, not virgin ground. The in-phase-vs-anti-phase question is settled in principle (it is governed by the competition between the spike and the sub-threshold waveform), and there is even a recent *two-unit folded-node "atlas"* (Awal–Epstein–Kaper–Vo) and two rigorous coupled-FHN MMO papers (Gonçalves–Labouriau–Rodrigues; Kristiansen–Pedersen). The part that is **genuinely open, and uniquely ours**, is the **noise layer on the cycle-fold**: nobody has computed the escape/peel-off *statistics* of a noisy folded **limit cycle** in a coupled or forced pair — which is exactly the open problem named in §9.4 of the Tracy–Widom paper. There is also a methodological landmine: the **Izhikevich reset makes gap-junction current during the spike ill-defined**, so the electrical-coupling atlas must either use a smooth model (Morris–Lecar / FitzHugh–Nagumo, both already in this repo) or a carefully patched Izhikevich. Recommendation in §6.

---

## 1. How this connects to the existing project

The uploaded paper — *The Noisy Folded Limit Cycle: Tracy–Widom Statistics of Canard Escape* (Williams, 2026) — analyses a **single, uncoupled** folded limit cycle: additive noise pushed through the JKK blow-up, the inner Riccati Cole–Hopf-linearised to the stochastic Airy operator, peel-off level `Y_node ≐ TW_β` with `β = 4/η²`. Its closing open problem (§9.4) is verbatim:

> "the Tracy–Widom identification suggests finer KPZ-class objects (the Airy process) may govern the joint law of successive peel-offs in a **forced or coupled folded cycle**."

So a coupled two-unit study is the paper's own designated sequel. Two structural facts to carry forward:

- The **deterministic backbone** of the folded-cycle theory is Jelbart–Kuehn–Kuntz (JKK), *Geometric blow-up for folded limit cycle manifolds in three time-scale systems* (arXiv:2208.01361, J. Nonlinear Sci. 34:17, 2024). Note JKK needs a **third timescale** to make blow-up work — relevant when we choose a coupled model (see §6).
- The conceptual bridge "**coupling produces folded structures of the cycle manifold**" already exists: Roberts–Rubin–Wechselberger (RRW), SIADS 14:1808–1844 (2015). This is the single most relevant prior paper and the natural anchor (see §5.5).

**Reusable code already in the repo** (from the folder survey): `kernel.py` (FHN2D + Euler–Maruyama engine), `regime-tests/coupled_array_peeloff.py` (diffusive coupling of folded-cycle canards — already found nearest-neighbour coupling gives exponential, not Airy r⁻², correlations), `bifurcation_classifier.py` (SNIC-vs-Hopf ISI fingerprinting), `mmo_fhr_*.py` (FitzHugh–Rinzel folded-node MMOs), `neural_mass_edge.py` / `mpr_finite_size.py` (Montbrió–Pazó–Roxin QIF mean-field), and the visualisation stack (`two_channels_animation.html`, `fhn_animator.jsx`, `morris_lecar_animator.html`).

---

## 2. Part I — What is classical / solved (do not re-derive; cite and build on)

### 2.1 Gap-junction synchrony vs anti-phase: the spike-vs-subthreshold dichotomy

The folk theorem "gap junctions synchronize" is **false in general**. The repeatedly-confirmed result is that the in-phase/anti-phase outcome is set by the **odd part of the interaction function** `H_odd(ψ)`, whose sign at `ψ = 0, π` is controlled by the competition between two pieces of the gap-junction current:

- **Suprathreshold (spike) component** — tall, fast spike injects current into the partner near *its* spike time, where the PRC is large → **pulls spikes together → in-phase**. Bigger/sharper spikes and deeper after-hyperpolarisation (AHP) strengthen this.
- **Sub-threshold component** — the slow inter-spike voltage ramp; when it dominates, weak coupling **stabilises anti-phase / splay**.

Load-bearing references: Chow & Kopell, *Neural Comput.* 12:1643 (2000); Lewis & Rinzel, *J. Comput. Neurosci.* 14:283 (2003); Gao & Holmes, *J. Comput. Neurosci.* 22:39 (2007, independent Poincaré-map confirmation); experimental confirmation in real interneurons by Mancilla, Lewis, Pinto, Rinzel & Connors, *J. Neurosci.* 27:2058 (2007) — where the *model* predicted anti-phase below ~30 Hz but real wide-spike/shallow-AHP cells showed none. The master variable behind this is **PRC type**: Type II (Hopf/resonator) robustly synchronises under electrical coupling; Type I (SNIC/integrator) is a weak synchroniser and anti-phase-prone (Ermentrout, *Neural Comput.* 8:979, 1996). Intrinsic currents tune the verdict: K⁺ currents synchronise, persistent Na⁺ (`I_NaP`) desynchronises (Pfeuty, Mato, Golomb & Hansel, *J. Neurosci.* 23:6280, 2003).

### 2.2 Emergent rhythms, bistability, and relaxation-oscillator coupling

- **Weak coupling can give anti-phase; strong coupling gives in-phase; the two routinely coexist (bistability).** Sherman & Rinzel, *PNAS* 89:2471 (1992) — and the same paper's rhythmogenesis ladder: two excitable cells → pacemaking; two pacemakers → **bursting**; two bursters → longer bursts. Cymbalyuk, Nikolaev & Borisyuk, *Biol. Cybern.* 71:153 (1994) catalogue **five** stable two-unit regimes including explicit in-phase/anti-phase bistability.
- **Relaxation oscillators synchronise by a different mechanism** — *fast-threshold modulation* (FTM), fast and nearly coupling-strength-independent: Somers & Kopell, *Biol. Cybern.* 68:393 (1993), with a companion result on anti-phase coexistence.

### 2.3 Heterogeneity (parameter mismatch)

- **Detuning → Arnold tongue → phase slips** (Adler equation). Locking only inside `|Δω| < K`. Textbook (Pikovsky–Rosenblum–Kurths, *Synchronization*, 2001).
- **Amplitude / oscillator death** in diffusively coupled pairs near Hopf, requiring frequency spread: Aronson, Ermentrout & Kopell, *Physica D* 41:403 (1990). ⚠ Note the title is "*Amplitude response of coupled oscillators*," often mis-cited as "amplitude death."
- **Emergent bursting helped by heterogeneity:** Sherman, Rinzel & Keizer, *Biophys. J.* 54:411 (1988); de Vries & Sherman, *Bull. Math. Biol.* 63:371 (2001) ("*From spikers to bursters via coupling: help from heterogeneity*").
- **Important constraint (a real result to exploit):** the standard reduced models (HH, Morris–Lecar, FHN, Hindmarsh–Rose) are **semi-passive / minimum-phase**, so diffusive coupling does **not** spontaneously create oscillation in them (Steur, Tyukin & Nijmeijer, arXiv:0903.3535, 2009). The classic Smale "diffusion-driven oscillation" paradox (rigorously: Pogromsky, Glad & Nijmeijer, *IJBC* 9:629, 1999) needs **non-minimum-phase** units (inferior-olive / Ca²⁺-type). This tells us which "emergent oscillation" claims are off-limits for FHN/ML pairs.

### 2.4 Delay (in the coupling / pulse transfer)

- **Delay enters as a phase shift** → stability switching between in-phase and anti-phase, multistability, destabilised synchrony: Crook, Ermentrout, Vanier & Bower, *J. Comput. Neurosci.* 4:161 (1997).
- **Delay-induced death even in identical pairs**, with multiply-connected "death islands": Reddy, Sen & Johnston, *PRL* 80:5109 (1998); experimental confirmation *PRL* 85:3381 (2000). Distributed delays enlarge the death region: Atay, *PRL* 91:094101 (2003).
- **Phase-flip bifurcation** (discontinuous 0 → π jump with a frequency jump): Prasad, Dana, Karnatak, Kurths, Blasius & Ramaswamy, *Chaos* 18:023111 (2008).
- **Gap-junction (pure electrical) delay specifically is comparatively under-mapped** — a handful of delayed-diffusive-FHN bifurcation studies show alternating near-synchronous/near-anti-phase windows via double-Hopf and torus bifurcations, but this is far thinner than the chemical-synapse-delay literature. (Candidate novelty axis — see §5/§6.)

### 2.5 Noise

- **Common noise synchronises uncoupled oscillators; independent noise desynchronises:** Teramae & Tanaka, *PRL* 93:204103 (2004); Goldobin & Pikovsky, *PRE* 71:045201(R) (2005).
- **Coherence resonance** — optimal noise maximises firing regularity in excitable FHN: Pikovsky & Kurths, *PRL* 78:775 (1997); array/coupling-enhanced version: Neiman et al., *PRL* 83:4896 (1999).
- **Self-induced stochastic resonance (SISR)** — small noise creates a robust limit cycle with no zero-noise counterpart, away from any bifurcation: Muratov, Vanden-Eijnden & E, *Physica D* 210:227 (2005); coherence-resonance vs SISR cleanly separated by DeVille, Vanden-Eijnden & Muratov, *PRE* 72:031105 (2005).
- **Noise-induced bursting:** Hitczenko & Medvedev, *SIAM J. Appl. Math.* (2009), arXiv:0712.4074.
- **Stochastic MMOs near a folded node** (single unit): Berglund, Gentz & Kuehn, *J. Dyn. Diff. Eq.* (2015), arXiv:1312.6353 — random Poincaré maps, a critical noise above which SAOs drown, and crucially: **noise makes paths peel off the canard early, with high probability.** This is the closest existing object to our "peel-off," but it stops at Markov-kernel estimates and SAO counts — **no extreme-value law, and at the equilibrium-fold (folded node), not the cycle-fold.**

### 2.6 Standard classification metrics (for the atlas pipeline)

Phase difference `ψ` / spike-time locking (0 = in-phase, π = anti-phase); Kuramoto order parameter `r` (for two units `r = |cos(ψ/2)|`); Golomb–Hansel synchrony measure `χ` (population-voltage variance ratio); spike-train cross-correlation; Floquet multipliers / Lyapunov exponents of the in-phase and anti-phase orbits for stability and chaos. Reviewed in Golomb, *Scholarpedia* 2(1):1347 (2007).

---

## 3. Part II — The Izhikevich choice: a methodological caveat ⚠ (read before coding)

You chose Izhikevich units. For the **firing-regime zoo** (RS / IB / CH / FS / TC / RZ / LTS — Izhikevich, *IEEE TNN* 14:1569, 2003; 15:1063, 2004) and for **spike-adding** it is an excellent, cheap choice. But for **electrical (gap-junction) coupling specifically there is a genuine problem**, and it is well documented:

- The Izhikevich model has a **hard reset** (`v ← c`, `u ← u+d` when `v ≥ 30`) and **does not integrate a real spike upstroke** — the "+30" is a bookkeeping cap, not a trajectory.
- A gap junction passes `I_gap = g·(v_j − v_i)` **at all times, including during the spike.** Because the spike is a vertical jump-and-reset rather than a real waveform, the gap-junction current *during the spike* — which is exactly the **dominant synchronising term** (§2.1) — is **ill-defined**: its size and time-integral depend on arbitrary choices (the cap value, the time step at the discontinuity).
- Consequences spelled out in the literature: you must **add the spike back by hand** as a δ-pulse of *free* magnitude, and the reset makes the limit cycle and the PRC **discontinuous**, breaking standard phase reduction (Chartrand, Goldman & Lewis, *SIADS* 2019, arXiv:1801.05874, on electrically coupled resonate-and-fire neurons; Desroches, Kowalczyk & Rodrigues, *Nonlinear Dyn.* 2021, arXiv:2101.12507). This is the most likely reason a clean two-cell gap-junction Izhikevich phase diagram was **never** established the way it was for Morris–Lecar or Hindmarsh–Rose.

By contrast Morris–Lecar, FHN/FitzHugh–Rinzel and HH integrate a **smooth spike**, so `I_gap` is well-defined throughout, the PRC is continuous, and — importantly for us — the **folded-cycle / canard / Tracy–Widom machinery applies directly.** The trade-off is only cost (Izhikevich ≈ 13 FLOPS/ms vs ML ≈ 600).

**Implication:** the broad "regime variety" survey can use Izhikevich, but anything that connects to the noisy-peel-off / Tracy–Widom story should be anchored in a **smooth** model. See the recommendation in §6.

---

## 4. Part III — Canards, MMOs and folded structures in coupled systems (the frontier nearest the paper)

This is where the project lives, and most of it is **recent and active** — so novelty must be defined carefully against it.

### 4.1 Coupling-induced MMOs in two FitzHugh–Nagumo units (two rigorous papers)

- **Gonçalves, Labouriau & Rodrigues (2025), arXiv:2503.12596** — two FHN coupled symmetrically through the slow equations. Symmetry forces an open set of **synchrony** and an open set of **antisynchrony** solutions, plus their **bistability** (coexisting hyperbolic attractors, persistent under perturbation). **Canards generate the MMOs** and seed small-amplitude transient oscillations before large relaxation spikes. Crucially, they note **asymmetric coupling = periodic forcing of one unit by the other** — a direct bridge to the paper's "forced folded cycle."
- **Kristiansen & Pedersen (2023), arXiv:2202.12027, SIADS 22:1383–1422 (10.1137/22M1480495)** — two FHN with symmetric *repulsive* coupling. Rigorous GSPT + blow-up shows the MMOs are organised by a **cusp** of the critical manifold, **not** a generic fold/folded-node; SAO count follows the Weber equation and eigenvalue ratio. The most technically rigorous coupling-induced-MMO result to date.
- Background phenomenology: "Mixed-mode oscillations and phase locking in coupled FitzHugh–Nagumo neurons," *Chaos* 29:033105 (2019), arXiv:1807.10824 (diffusive coupling → MMOs, period-adding, oscillation death).

### 4.2 Canard-mediated (de)synchronisation

- **Ersöz, Desroches & Krupa, *Physica D* (2017)** — for weakly coupled canard oscillators, the **synchronisation properties change precisely at the maximal canard.** The canard itself is the organising object for the in-phase/anti-phase transition.
- **Ersöz, Desroches, Krupa & Clément, *SIADS* (2016), 10.1137/15M101840X** — canard-mediated (de)synchronisation in coupled phantom bursters.
- **Ermentrout & Wechselberger, *SIADS* 8:253 (2009)** — the ancestor: folded-node MMO reduction of **gap-junction-coupled** interneurons producing clustered/synchronised states. (Gap junction + canards → clusters.)

### 4.3 The existing two-unit folded-node "atlas" (the thing to beat)

**Awal, Epstein, Kaper & Vo:** *Chaos* 33:011102 (2023, coupled van der Pol & Lengyel–Epstein); *J. Nonlinear Sci.* 34:53 (2024, coupled identical Lengyel–Epstein); and a 2024 coupled-Koper-oscillators paper. Unified message: in two symmetrically coupled **identical** fast–slow units, a **folded node off the symmetry axis** forces **order-of-magnitude amplitude symmetry-breaking** (one unit small/MMO, the other large relaxation), organised by an **asymmetric canard explosion** and an **explosion of anti-phase limit-cycle canards.** This is effectively a deterministic two-unit canard atlas already partly built — our atlas must cite it and **add a dimension it does not have** (noise / peel-off statistics).

### 4.4 Spike-adding and torus canards

- **Spike-adding (n → n+1 spikes/burst) is canard-organised:** folded-saddle canards in parabolic bursters (Desroches, Krupa & Rodrigues, *Physica D* 331:58, 2016); spike-adding canard *explosions* in square-wave bursters (Desroches, Kaper & Krupa, *Chaos* 23:046106, 2013; *J. Nonlinear Sci.* 30, 2020); codim-2 homoclinic organisation in Hindmarsh–Rose (Linaro, Champneys, Desroches & Storace, *SIADS* 11:939, 2012).
- **Torus canards** = the **cycle-level analogue of canards**: a trajectory follows the *attracting* branch of a manifold of fast-subsystem limit cycles and then the *repelling* branch past a **fold/saddle-node of periodic orbits (SNPO)**. They are the **generic spiking↔bursting separatrix**: Benes, Barreto, Burke, Desroches, Kaper & Kramer, *J. Math. Neurosci.* 2:3 (2012), arXiv:1107.2834; made robust by a second slow variable (Vo, "Generic Torus Canards," *Physica D* 2017). **Under periodic forcing, spike-adding canards continue into torus canards** (Burke, Desroches, Barreto, Kaper, Kramer, *J. Nonlinear Sci.* 2015, arXiv:1504.03970; Vo, arXiv:1607.02205, 2016).

### 4.5 The conceptual anchor: coupling → folded singularities of the cycle manifold

**Roberts, Rubin & Wechselberger, *SIADS* 14:1808–1844 (2015), 10.1137/140981770.** In a **coupled pair** of burst-capable neurons, fast–slow **averaging** turns the cycle-manifold dynamics into a system with **folded singularities (FSN II / III)**, and these correspond to **singular torus bifurcations** of the full coupled model. In other words: *coupling literally manufactures folded singularities of the limit-cycle manifold, and the spiking↔bursting transition is read off from them.* This is the deterministic skeleton onto which our **noise layer** attaches.

---

## 5. Part IV — The Tracy–Widom / KPZ neighbour, and why it is *not* our result

There is exactly one place in the literature where limit-cycle oscillators and Tracy–Widom co-occur: **Gutiérrez & Cuerno, *Phys. Rev. Research* 6:033324 (2024), arXiv:2311.13253** (precursors PRR 5:023047, 2023; Lauter–Mitra, PRE 96:012220, 2017). As **1D rings of Stuart–Landau / van der Pol oscillators synchronise**, the phase front shows **KPZ kinetic roughening** and the height fluctuations are **Tracy–Widom distributed.**

**Why it is a neighbour and not prior art for us:** it is a **large-N spatial phase-front roughening** effect near a Hopf bifurcation — *not* canard escape, *not* a fold/SNPO of cycles, *not* successive peel-offs, *not* two units. The TW there comes from interface universality, a **different mechanism** than canard-escape extreme statistics. **Skeptic's flag:** classical TW/KPZ universality usually needs a *large* interface; whether anything TW-like survives down at **N = 2** in a canard-escape setting is genuinely open and may have a **negative** answer. We frame the TW-at-N=2 question as a hypothesis to test, not an assumed extension.

---

## 6. Part V — Where the genuinely new territory is

### Already done (cite, build on, do not claim as novel)

- Coupling-induced MMOs in two FHN — numerically/symmetry (2503.12596) and rigorously via a cusp singularity (2202.12027).
- Canard-mediated in-phase/anti-phase (de)synchronisation; folded-node strong symmetry-breaking in identical pairs (Ersöz et al.; Ermentrout–Wechselberger; Awal–Epstein–Kaper–Vo).
- Spike-adding by canards; torus canards as the spiking↔bursting separatrix; forced spike-adding → torus canards.
- Coupling → folded singularities of the **averaged cycle manifold** (Roberts–Rubin–Wechselberger); rigorous **folded-limit-cycle blow-up** (JKK 2024).
- Noisy folded-**node** canards with early peel-off w.h.p. and SAO-count change (Berglund–Gentz–Kuehn 2015).
- The classical synchrony/anti-phase/death/phase-flip/coherence-resonance results of Part I.

### The real open gap (this is the project's wedge)

1. **No one has computed the escape/peel-off *statistics* of a noisy folded *limit cycle*** — i.e. at the **cycle-fold / SNPO (torus-canard) level**, not the equilibrium-fold. Berglund–Gentz–Kuehn give the *qualitative* "early peel-off w.h.p." but stop short of an extreme-value/Tracy–Widom law, and at the folded **node**. The cycle-level noisy peel-off law appears unaddressed. **Your single-unit paper already cracked the cycle-fold law (TW_β); extending it is the natural, unclaimed next step.**
2. **The joint law of successive peel-offs in a forced or coupled folded cycle is unstudied** — your paper's §9.4 problem. Successive-peel-off correlations are the natural home for an Airy-process / KPZ-type joint law, and nobody has set the problem up. Asymmetric coupling = forcing (per 2503.12596) is the cleanest first instance.
3. **All three perturbations (heterogeneity + delay + noise) on a single electrically coupled slow-fast pair, simultaneously** — no comprehensive two-unit treatment exists. In particular **delay in pure electrical coupling** of bursting units is thin, and **noise near canards in coupled pairs** is essentially open.
4. **A unified atlas tying (coupling/forcing) → folded singularities of the cycle manifold → torus canards → noisy peel-off statistics** does not exist as a single object. The pieces (RRW 2015 geometry; Vo 2016/2017 forced/generic torus canards; BGK 2015 noise) have never been assembled, and never with the **escape-statistics** layer that is your signature.

### Recommended positioning of the atlas (concrete)

- **Two-model strategy.** Use **Izhikevich** for the broad **regime-variety survey** (the bursting/spike-adding/CH/IB zoo, cheap parameter sweeps) — but **anchor every result that touches the canard / folded-cycle / Tracy–Widom story in a smooth model**: **Morris–Lecar** or **FitzHugh–Rinzel** (you already have `morris_lecar_animator.html` and `mmo_fhr_*.py`). If you want to keep Izhikevich for the electrical-coupling part, use a **soft reset** and an explicit **δ-spike patch**, and *validate the synchrony verdicts against Morris–Lecar* (treat the reset model as out-of-distribution until checked — exactly as your `morris_lecar_check.py` did for the single-unit classifier).
- **Make noise a first-class axis, not an afterthought.** The deterministic regime atlas is the *map*; the **noisy peel-off statistics on the cycle-fold are the contribution.** That is the one layer none of the competing papers (Awal–Epstein–Kaper–Vo, Gonçalves et al., RRW) carry.
- **Lead instance for the §9.4 problem:** asymmetric coupling (= periodic forcing of one folded cycle by its partner) → measure the **joint distribution of successive peel-off levels** and test for Airy-process structure vs the single-unit TW_β marginal. `coupled_array_peeloff.py` already shows the *spatial-array* version gives only short-range (exponential) correlations — a recorded **negative result** that sharpens the question for the **two-unit forced** case.
- **Pressure-test, don't assume, TW at N = 2.** Given the Gutiérrez–Cuerno large-N mechanism, treat "does TW/Airy survive in a 2-unit canard system?" as the headline question — a clean negative is still a publishable, paper-consistent finding.

---

## 7. Annotated reference ledger

**Gap-junction synchrony (Part I):**
Sherman & Rinzel, *PNAS* 89:2471 (1992); Chow & Kopell, *Neural Comput.* 12:1643 (2000); Lewis & Rinzel, *J. Comput. Neurosci.* 14:283 (2003); Mancilla, Lewis, Pinto, Rinzel & Connors, *J. Neurosci.* 27:2058 (2007); Gao & Holmes, *J. Comput. Neurosci.* 22:39 (2007); Pfeuty, Mato, Golomb & Hansel, *J. Neurosci.* 23:6280 (2003); Ermentrout, *Neural Comput.* 8:979 (1996); Cymbalyuk, Nikolaev & Borisyuk, *Biol. Cybern.* 71:153 (1994); Somers & Kopell, *Biol. Cybern.* 68:393 (1993); Golomb, *Scholarpedia* 2(1):1347 (2007).

**Izhikevich + reset caveat (Part II):**
Izhikevich, *IEEE TNN* 14:1569 (2003); 15:1063 (2004); *Dynamical Systems in Neuroscience* (MIT Press, 2007); Chartrand, Goldman & Lewis, *SIADS* (2019), arXiv:1801.05874; Desroches, Kowalczyk & Rodrigues, *Nonlinear Dyn.* (2021), arXiv:2101.12507; Linaro, Champneys, Desroches & Storace, *SIADS* 11:939 (2012); Nobukawa, Nishimura & Yamanishi, *Sci. Rep.* 7:1331 (2017).

**Heterogeneity / delay / noise (Part I):**
Aronson, Ermentrout & Kopell, *Physica D* 41:403 (1990) ⚠title; Reddy, Sen & Johnston, *PRL* 80:5109 (1998) & *PRL* 85:3381 (2000); Atay, *PRL* 91:094101 (2003); Prasad et al., *Chaos* 18:023111 (2008); Crook, Ermentrout, Vanier & Bower, *J. Comput. Neurosci.* 4:161 (1997); Pikovsky & Kurths, *PRL* 78:775 (1997); Teramae & Tanaka, *PRL* 93:204103 (2004); Goldobin & Pikovsky, *PRE* 71:045201 (2005); Muratov, Vanden-Eijnden & E, *Physica D* 210:227 (2005); Hitczenko & Medvedev, arXiv:0712.4074 (2009); Sherman, Rinzel & Keizer, *Biophys. J.* 54:411 (1988); de Vries & Sherman, *Bull. Math. Biol.* 63:371 (2001); Steur, Tyukin & Nijmeijer, arXiv:0903.3535 (2009); Pogromsky, Glad & Nijmeijer, *IJBC* 9:629 (1999).

**Coupled MMO / canard / folded cycle (Parts III–IV):**
Gonçalves, Labouriau & Rodrigues, arXiv:2503.12596 (2025); Kristiansen & Pedersen, *SIADS* 22:1383 (2023), arXiv:2202.12027; "MMOs and phase locking in coupled FHN," *Chaos* 29:033105 (2019), arXiv:1807.10824; Ersöz, Desroches & Krupa, *Physica D* (2017); Ersöz, Desroches, Krupa & Clément, *SIADS* (2016); Ermentrout & Wechselberger, *SIADS* 8:253 (2009); Awal, Epstein, Kaper & Vo, *Chaos* 33:011102 (2023) & *J. Nonlinear Sci.* 34:53 (2024); Desroches, Krupa & Rodrigues, *Physica D* 331:58 (2016); Desroches, Kaper & Krupa, *Chaos* 23:046106 (2013); Benes et al., *J. Math. Neurosci.* 2:3 (2012), arXiv:1107.2834; Vo, *Physica D* (2017) & arXiv:1607.02205 (2016); Burke et al., arXiv:1504.03970 (2015); Roberts, Rubin & Wechselberger, *SIADS* 14:1808 (2015); Jelbart, Kuehn & Kuntz, arXiv:2208.01361, *J. Nonlinear Sci.* 34:17 (2024); Berglund, Gentz & Kuehn, arXiv:1312.6353 (2015); Gutiérrez & Cuerno, *PRR* 6:033324 (2024), arXiv:2311.13253.

⚠ **Verification flags carried forward:** (i) the only dedicated *two-cell gap-junction Izhikevich* note found is a 2-page NOLTA-2016 conference paper whose full text could not be retrieved — treat as unverified. (ii) Aronson–Ermentrout–Kopell 1990 is "*Amplitude response…*," not "amplitude death." (iii) Bar-Eli effect and Loewenstein–Yarom–Sompolinsky *PNAS* 2001 were not read in primary this round — verify before formal citation. (iv) TW-at-N=2 is a hypothesis, not a literature result.

---

## 8. Suggested next steps (pending your go-ahead)

1. Lock the model decision (§3/§6): Izhikevich-only with spike-patch, or Izhikevich-survey + Morris–Lecar/FHR-anchor (recommended).
2. Build the coupled integrator extending `kernel.py`; reproduce two single-unit limits and one known result (e.g. Sherman–Rinzel emergent bursting, or the 2503.12596 synchrony/antisynchrony bistability) as a validation gate.
3. Produce the deterministic regime atlas over (coupling strength, coupling form, heterogeneity, delay) with the §2.6 classifiers and overlaid bifurcation curves.
4. Add the **noise axis** and measure peel-off statistics on the cycle-fold — single unit first (reproduce TW_β), then forced/coupled (the §9.4 joint-law question).
5. Wrap with the interactive two-unit phase-portrait explorer and the written report.
