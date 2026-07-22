# Theory roadmap — everything underpinning the stochastic folded-limit-cycle project

A dependency-ordered learning path from foundations to the Tracy–Widom frontier.
Each topic: **What** it is · **Why (for this project)** · **Source(s)** · **✓ Check**
(a question you should be able to answer out loud). `★` = own this one *cold*; it is
load-bearing and an examiner will probe it.

**How to read this.** Layers 0→4 are roughly prerequisite-ordered. Layer 0 you can
skim (refresh, don't re-learn). Layers 1–3 are the real substance and map one-to-one
onto what you did: GSPT/blow-up (the method), stochastic escape (the noise), and the
stochastic-Airy/Tracy–Widom edge (the payoff). Layer 4 is only if you push the proof.
If time is short, follow **the critical-path spine** at the end — eight items, in order.

The through-line in one sentence: *FHN's fold → made rigorous by blow-up → with noise
(Berglund–Gentz) → lifted from a fold of equilibria to a fold of limit cycles (JKK) →
whose inner escape equation is the stochastic Airy operator (Tracy–Widom).*

---

## Layer 0 — Foundations (refresh, don't relearn)

**0.1 ODEs & planar phase-plane analysis.**
*What:* existence/uniqueness, equilibria, linearisation, nullclines, phase portraits.
*Why:* FHN is a 2D phase plane; the cubic nullcline and its knees (folds) are the
whole story. *Source:* Strogatz, *Nonlinear Dynamics and Chaos*; Perko.
*✓ Check:* sketch the FHN nullclines and locate the two fold points.

**0.2 Bifurcation theory.**
*What:* saddle-node, Hopf (sub/supercritical), homoclinic, **saddle-node of limit
cycles (fold of cycles)**, and the **Bautin / generalised Hopf** codim-2 point.
*Why:* your central object is literally a *fold of limit cycles*; your worked model is
the Bautin normal form (where C_q ≈ 5.7 was measured). *Source:* Kuznetsov, *Elements
of Applied Bifurcation Theory*, chs 3–8. *✓ Check:* why does a Bautin point guarantee
a branch of folds of cycles nearby?

**0.3 Normal forms & centre manifolds.**
*What:* reducing a system near a bifurcation to its essential polynomial form.
*Why:* JKK's normal form (their eq 14: `r′=−a(θ)y+b(θ)r²+…`) is what you noised.
*Source:* Kuznetsov ch 5; Guckenheimer–Holmes ch 3. *✓ Check:* identify a, b, c in
the JKK normal form and what each controls.

**0.4 Probability & stochastic processes.**
*What:* random variables, Gaussians, Markov processes, **Brownian motion / Wiener
process**, martingales (lightly), conditional expectation. *Why:* the noise term, the
escape statistics, the exit distribution. *Source:* Øksendal, *SDEs*, chs 1–3;
Durrett, *Probability*. *✓ Check:* state the defining properties of Brownian motion
and why `⟨ξ(t)ξ(t′)⟩=δ(t−t′)` for white noise.

**0.5 Linear operators & Sturm–Liouville spectral theory.**
*What:* self-adjoint operators, eigenvalues/eigenfunctions, Sturm–Liouville problems,
**oscillation theory (nodes ↔ eigenvalue count)**. *Why:* escape ⇔ first *node* of the
Cole–Hopf solution `u`; eigenvalues of the (stochastic) Airy operator are the exit
levels. *Source:* Teschl, *ODEs and Dynamical Systems*, ch on Sturm–Liouville.
*✓ Check:* why does the k-th eigenfunction have exactly k nodes?

---

## Layer 1 — Geometric singular perturbation theory (the method) ★

**1.1 Slow–fast systems & the singular limit.**
*What:* fast/slow splitting, critical manifold, layer vs reduced problems, Tikhonov.
*Why:* the entire setting; `0<ε₂≪ε₁≪1`. *Source:* **Kuehn, *Multiple Time Scale
Dynamics*** (the field's reference book) chs 1–4; Jones, *GSPT* lecture notes.
*✓ Check:* write FHN as a slow–fast system and identify the critical manifold.

**1.2 Fenichel theory. ★**
*What:* normally hyperbolic invariant manifolds persist (slow manifolds), with
Fenichel coordinates. *Why:* the attracting/repelling slow manifolds the canard
connects; where Berglund–Gentz tubes live. *Source:* Kuehn ch 3; Fenichel (1979);
Jones (1995). *✓ Check:* where does normal hyperbolicity fail, and why does that break
Fenichel?

**1.3 Relaxation oscillations.**
*What:* van der Pol / FHN limit cycles alternating slow drift and fast jumps at folds.
*Why:* your starting model; the fold is the degenerate point. *Source:* Kuehn ch 7;
Mishchenko–Rozov. *✓ Check:* why do the jumps occur exactly at the fold points?

**1.4 Canards & canard explosion. ★**
*What:* trajectories that follow a *repelling* slow manifold; the explosive growth of a
canard cycle over an exponentially thin parameter window. *Why:* the noisy canard is
*the* object you track; "Channel A" is canard escape. *Source:* **Krupa–Szmolyan
(2001)**, *Relaxation oscillation and canard explosion* (J. Diff. Eq.); Wechselberger
(Scholarpedia, canards). *✓ Check:* why is a canard structurally a slow passage through
a fold, and why is it exponentially sensitive?

**1.5 The blow-up method (geometric desingularisation). ★★**
*What:* rescaling that "blows up" a degenerate point to a sphere/cylinder, recovering
hyperbolicity; analysis in directional **charts K₁, K₂, K₃**; the rescaling chart's
**Riccati equation** and its **Airy** linearisation. *Why:* THE technique of your whole
project — every result lives in these charts. *Source:* **Krupa–Szmolyan (2001)**,
*Extending GSPT to nonhyperbolic points* (SIAM J. Math. Anal.); Kuehn ch 7;
Szmolyan–Wechselberger. *✓ Check:* derive the inner Riccati `R′=R²−Y` in the rescaling
chart and Cole–Hopf it to Airy `u″=Yu`.

**1.6 Folded singularities & mixed-mode oscillations (context).**
*What:* folded node/saddle/focus in ≥3D, canard-induced small oscillations, MMOs.
*Why:* situates "folded limit cycle" in the folded-singularity family; the noised
folded *node* is your nearest neighbour (BGK). *Source:* **Desroches et al (2012)**,
*Mixed-Mode Oscillations…* (SIAM Review); Wechselberger (2005, SIADS). *✓ Check:* how
does a folded limit cycle differ from a folded node — manifold of cycles vs of points?

**1.7 Fold of limit cycles & the JKK construction. ★★**
*What:* a saddle-node of periodic orbits, blown up; the periodic analogue of the fold.
*Why:* this *is* your object; JKK is your deterministic scaffold. *Source:* **Jelbart–
Kuehn–Kuntz (2024)**, arXiv:2208.01361 (own §§4.2–4.5, 5.1 cold); Floquet theory below.
*✓ Check:* what plays the role of the slow variable, and where is the Airy inner here?

**1.8 Floquet theory & phase reduction.**
*What:* stability of periodic orbits (Floquet multipliers/exponents), isochrons, the
**phase-response curve / infinitesimal PRC**. *Why:* the limit-cycle manifold's
transverse structure; "Channel B" is phase diffusion via the iPRC. *Source:* Kuznetsov
(Floquet); **Ermentrout–Terman, *Mathematical Foundations of Neuroscience*** (PRC,
isochrons); Schwemmer–Lewis. *✓ Check:* relate the Floquet exponent to the contraction
rate of your concentration tube.

---

## Layer 2 — Stochastic dynamics & escape (the noise) ★

**2.1 Stochastic differential equations.**
*What:* Itô integral, **Itô formula**, Stratonovich, the Itô↔Strato conversion,
existence/uniqueness. *Why:* the noisy normal form; your Cole–Hopf is a *Stratonovich*
identity (caveat C-i — the Itô correction is real). *Source:* **Øksendal, *SDEs***;
Karatzas–Shreve. *✓ Check:* convert `dR=(R²−Y)dT+η dB` between Itô and Stratonovich and
say which Cole–Hopf assumes.

**2.2 Diffusions, generators, Fokker–Planck.**
*What:* infinitesimal generator, **Kolmogorov forward (Fokker–Planck)** and **backward**
equations. *Why:* exit densities; the backward/committor PDE (your Conjecture 3).
*Source:* **Pavliotis, *Stochastic Processes and Applications*** chs 4, 6–7; Gardiner.
*✓ Check:* write the backward Kolmogorov equation for the committor `q(R,Y)`.

**2.3 First-passage / exit problems & 1D diffusions.**
*What:* mean first-passage time, exit distributions; **scale function & speed measure**
for 1D diffusions. *Why:* escape = first passage; you used the scale-function escape
rate. *Source:* Karlin–Taylor, *A Second Course in Stochastic Processes*; Pavliotis ch
7. *✓ Check:* derive the escape rate from the scale function for a 1D gradient drift.

**2.4 Large deviations & Kramers escape. ★**
*What:* **Freidlin–Wentzell** action functional, **quasipotential**, Kramers/Arrhenius
law, the **Eyring–Kramers prefactor**. *Why:* your barrier `V(Y)=8Y^{3/2}/3`, integrated
hazard `H=η²/4π`, and `C_q=2√π`; the *prefactor* is exactly the open piece. *Source:*
**Freidlin–Wentzell, *Random Perturbations of Dynamical Systems*** (exit chapter);
Berglund, *Kramers' law* lecture notes (very readable); Bouchet–Reygner (Eyring–
Kramers). *✓ Check:* derive `H=η²/4π` from `λ(Y)=(√Y/π)e^{−V/η²}` by the `u=8Y^{3/2}/3η²`
substitution.

**2.5 Berglund–Gentz theory. ★★**
*What:* sample-path **concentration tubes** around deterministic solutions for slow–fast
SDEs; the noisy fold, the noisy canard, the noisy folded node. *Why:* THE framework
your chart-matching/concentration estimate `Var(r₂)=η²/4√(ab·y₂)` is built on, and your
nearest prior art. *Source:* **Berglund–Gentz, *Noise-Induced Phenomena in Slow–Fast
Dynamical Systems*** (2006); Berglund–Gentz–Kuehn, *Hunting French ducks in a noisy
environment* (2012) and the folded-node papers. *✓ Check:* state the tube/escape
dichotomy and how the tube width sets the escape threshold.

**2.6 Stochastic averaging (Khasminskii).**
*What:* fast variables are replaced by their time/measure average in the slow drift, with
`O(1/ω)` corrections. *Why:* the α=1 averaged law `√(⟨a⟩⟨c⟩/⟨b⟩)` and the "A preempts
B" race. *Source:* Pavliotis–Stuart, *Multiscale Methods*; Freidlin–Wentzell (averaging
ch); Khasminskii. *✓ Check:* why do the non-zero θ-harmonics average away under fast
rotation?

**2.7 Stochastic phase reduction (Channel B).**
*What:* reducing a noisy oscillator to a phase variable; phase diffusion `D_φ`; the
Itô/Stratonovich subtlety *in the phase equation*. *Why:* Channel B's amplitude
`η/2πR`, the SNIC-vs-fold exponent finding (0.83 ≠ ½). *Source:* **Goldobin–Teramae–
Nakao–Ermentrout (2010)**, *Dynamics of limit-cycle oscillators subject to general
noise*; Teramae–Tanaka. *✓ Check:* why is the phase-reduction drift correction
Stratonovich-vs-Itô sensitive?

---

## Layer 3 — Integrable & random-matrix structure (the payoff) ★

**3.1 Riccati ↔ linear ODE (Cole–Hopf) & oscillation theory. ★**
*What:* `R=−u′/u` turns a Riccati into a 2nd-order linear ODE; blow-ups of `R` ↔ nodes
of `u`; Prüfer angle / rotation number counts them. *Why:* your exact linearisation
`u″=(Y−ηξ)u`; escape ⇔ first node. *Source:* Coddington–Levinson; Teschl (Prüfer,
oscillation). *✓ Check:* show a finite-time blow-up of `R` is a simple zero of `u`.

**3.2 Airy equation, functions, and operator.**
*What:* `u″=Yu`, `Ai/Bi`, the Airy operator `−d²/dx²+x` and its spectrum (eigenvalues =
`−`Airy-zeros on the half-line). *Why:* the deterministic inner is Airy; `Ω₀` is the
first Airy zero `≈−2.338`. *Source:* DLMF ch 9; Vallée–Soares, *Airy Functions*.
*✓ Check:* why is the canard escape at the first Airy zero, and what BC selects `Ai`?

**3.3 1D random Schrödinger / Anderson localisation. ★**
*What:* `−u″+V u` with random `V`; **Frisch–Lloyd equation**, integrated density of
states (IDOS), Lyapunov exponent, rotation number — all read off the **Riccati with
white noise**. *Why:* your escape rate = IDOS/rotation number of `u″=(Y−ηξ)u`
(Conjecture 2, the most computable route). *Source:* **Comtet–Texier–Tourigny**
(reviews on the Riccati/IDOS method); Lifshitz–Gredeskul–Pastur, *Introduction to the
Theory of Disordered Systems*; Frisch–Lloyd (1960). *✓ Check:* how does the stationary
Riccati current give the density of states?

**3.4 The stochastic Airy operator & β-ensemble edge. ★★**
*What:* `H_β=−d²/dx²+x+(2/√β)b′(x)`; its ground state is characterised by a **Riccati
diffusion** — *your inner Riccati* — and is the soft-edge limit of β-ensembles.
*Why:* THE identification: `η=2/√β`; the inner exit measure is `H_β`'s ground-state law.
*Source:* **Ramírez–Rider–Virág (2011)**, *Beta ensembles, stochastic Airy spectrum, and
a diffusion* (J. AMS); Edelman–Sutton (2007), *From random matrices to stochastic
operators*. *✓ Check:* match RRV's diffusion to `dR=(R²−Y)dT+η dB` and read off
`η=2/√β`.

**3.5 Tracy–Widom distributions. ★★**
*What:* `TW₁/₂/₄`, the Fredholm-determinant and **Painlevé II / Hastings–McLeod**
representations, the asymmetric tails (`e^{−βs³/24}` left, `e^{−(2β/3)s^{3/2}}` right).
*Why:* the inner exit measure **is** `TW_β`; you computed it via Painlevé II and matched
it parameter-free at β=1,2. *Source:* **Tracy–Widom (1994, 1996)**; **Bornemann (2010)**,
*On the numerical evaluation of distributions in RMT* (the numerics you reproduced);
Akemann–Baik–Di Francesco, *Oxford Handbook of RMT*. *✓ Check:* why does the early-escape
(heavy) tail correspond to the TW left large-deviation?

**3.6 KPZ universality & the Airy process (context).**
*What:* the universality class whose edge fluctuations are Tracy–Widom; the Airy₂
process. *Why:* situates your result — folded-cycle escape sits in the TW/KPZ edge
class. *Source:* **Corwin (2012)**, *The KPZ equation and universality class* (review);
Quastel–Spohn. *✓ Check:* in one sentence, what makes a statistic "TW-universal"?

---

## Layer 4 — Frontier (only if you push the proof to a theorem)

**4.1 Sharp prefactors: instantons & functional determinants.**
*What:* one-loop fluctuation determinant about the optimal escape path (**Gelfand–
Yaglom**), zero modes, caustics/instanton-merging. *Why:* Conjecture 4 for the
prefactor; the vanishing-barrier caustic rejoins the Airy edge. *Source:* Gelfand–
Yaglom; Schäfer–Kleinert (instantons). *✓ Check:* why does the determinant diverge as
the barrier vanishes, and what fixes it?

**4.2 Rigorous stochastic operators & spectral convergence.**
*What:* making "β-ensemble edge → stochastic Airy operator" and the Itô↔Strato Cole–Hopf
rigorous. *Why:* turns your confirmed identification into a lemma. *Source:* RRV (2011)
proofs; Gorin–Shkolnikov; Dumaz–Virág (SAO tails). *✓ Check:* what does the Itô
correction add to the Cole–Hopf, and at what order (you measured `O(η²)`)?

**4.3 Matched asymptotics / boundary-layer rigour.**
*What:* sharp uniform error bounds when composing the chart maps `Π₁∘κ₁₂∘Π₂∘κ₂₃∘Π₃`.
*Why:* the remaining ~50% of the rigorous theorem. *Source:* Kuehn ch 9 (matching);
JKK §4.5 (the deterministic template). *✓ Check:* where in the composition does the
non-quasi-static (near-fold) correction enter?

---

## The critical-path spine (if time is short, in this order)

1. **Kuehn**, *Multiple Time Scale Dynamics* — GSPT, slow manifolds, blow-up (1.1–1.5).
2. **Krupa–Szmolyan (2001)** — the fold blow-up, canards, the Riccati/Airy inner.
3. **JKK (2024)** — your scaffold: the fold of limit cycles, blown up.
4. **Berglund–Gentz (2006)** — concentration tubes, noisy fold/canard.
5. **Freidlin–Wentzell** (exit ch) + **Berglund's Kramers notes** — the barrier & hazard.
6. **Comtet–Texier–Tourigny** — Riccati + white noise = IDOS (the localisation bridge).
7. **Ramírez–Rider–Virág (2011)** — the stochastic Airy operator (the identification).
8. **Tracy–Widom (1994)** + **Bornemann (2010)** — the edge law and its numerics.

## The papers to own cold (an examiner will open these)

JKK 2024 · Krupa–Szmolyan 2001 · Berglund–Gentz(–Kuehn) · Freidlin–Wentzell (exit) ·
Ramírez–Rider–Virág 2011 · Tracy–Widom 1994 (+ Bornemann 2010 for the numerics).

## How this maps to your files

GSPT/blow-up → `FOLDED_CYCLE_PATHA.md`, `folded_cycle_patha.py` · concentration tube →
`folded_cycle_chartmatch.py` · FW barrier/hazard/`C_q` → `FOLDED_CYCLE_PATHA_PROOFS.md` ·
stochastic-Airy/Tracy–Widom → `FOLDED_CYCLE_NOISYAIRY_CONJECTURES.md`,
`folded_cycle_tracy_widom.py`, `folded_cycle_tw_pdf.py` · averaging/phase → the α=1 and
Channel-B chapters.
