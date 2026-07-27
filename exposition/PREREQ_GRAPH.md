# Prerequisite graph — the hard skeleton
Draft 2026-07-27, pre-research. Fixed by the mathematics, not by pedagogy: the research can
change where we *start* and how much we *motivate*, but not these dependencies.

## How to read
- **HARD edge (⇒)**: B cannot be *stated* without A. Non-negotiable ordering.
- **SOFT edge (→)**: B is easier with A, but B can be motivated independently. This is where
  all the sequencing freedom lives — and where most courses waste time teaching A "because
  you'll need it", when they won't.
- **Phenomenon / mechanism split**: several nodes can be *shown* long before they can be
  *explained*. Showing early is a hook; explaining early is a tax.

---

## Nodes

### Trunk α — geometry & dynamics
| id | node |
|----|------|
| α1 | phase plane, nullclines, equilibria |
| α2 | two timescales: fast/slow splitting, critical manifold |
| α3 | relaxation oscillation |
| α4 | fold of equilibria (saddle-node) |
| α5 | Fenichel / normal hyperbolicity — **and its failure at the fold** |
| α6 | canards (trajectory following a *repelling* manifold) |
| α7 | blow-up desingularisation, charts K1/K2/K3 |
| α8 | the inner Riccati equation |
| α9 | Cole–Hopf: Riccati ⇒ linear 2nd-order ODE |
| α10 | Airy equation, Airy function, Airy zeros |
| α11 | oscillation theory: nodes ↔ eigenvalue index |

### Trunk β — noise & large deviations
| id | node |
|----|------|
| β1 | Brownian motion, SDEs, Euler–Maruyama |
| β2 | first-passage / exit problem |
| β3 | Freidlin–Wentzell action, quasipotential |
| β4 | Kramers / Arrhenius escape rate |
| β5 | stochastic Airy operator |
| β6 | Tracy–Widom law |
| β7 | β = 4/η² identification |

### Trunk γ — asymptotics & resurgence
| id | node |
|----|------|
| γ1 | asymptotic series; **they diverge** |
| γ2 | WKB, turning points |
| γ3 | Borel transform |
| γ4 | Borel plane, singularities, **branch cuts** |
| γ5 | Borel summation — and the **contour ambiguity** across a cut |
| γ6 | *the ambiguity is exponentially small, e^{−A/x}, A = singularity location* |
| γ7 | trans-series |
| γ8 | Stokes phenomenon |
| γ9 | resurgence: large-order coefficients encode the non-perturbative data |
| γ10 | Padé / Borel–Padé |
| γ11 | Darboux large-order fitting |

### Trunk δ — the ladder & the program's own results
| id | node |
|----|------|
| δ1 | turning order q: potential ~ Y^q |
| δ2 | Bessel-1/(q+2) skeleton |
| δ3 | instanton action s^(2q+1)/[2(2q+1)] |
| δ4 | Thom A_k classification, unfoldings, caustics |
| δ5 | the two-ladders split (swept vs genuine Berry–Upstill) |
| δ6 | the v_k coefficient ladder ← **THE JUNCTION NODE** |
| δ7 | Wiener chaos / transfer engine |
| δ8 | non-integrability of higher rungs |
| δ9 | the headline claim |

---

## Hard edges (⇒)

```
α1 ⇒ α2 ⇒ α3
α1 ⇒ α4
α2 ⇒ α5
α4 + α5 ⇒ α7          blow-up is UNMOTIVATED without hyperbolicity failing at the fold
α7 ⇒ α8 ⇒ α9 ⇒ α10 ⇒ α11
α2 + α4 ⇒ α6          canard needs both branches + the fold

β1 ⇒ β2 ⇒ β3 ⇒ β4
α10 + β1 ⇒ β5 ⇒ β6
β5 + β6 ⇒ β7

γ1 ⇒ γ3 ⇒ γ4 ⇒ γ5 ⇒ γ6 ⇒ γ7 ⇒ γ8 ⇒ γ9
γ1 ⇒ γ2
γ3 ⇒ γ10
γ9 ⇒ γ11

α10 + δ1 ⇒ δ2         (via γ2 for general q; α10 alone suffices at q=1)
β3 + δ1 ⇒ δ3
δ6 ⇒ γ1-instance      the ladder is what actually diverges
β1 + δ6 ⇒ δ7
β6 + δ2 ⇒ δ8          non-integrability is only meaningful BY CONTRAST with q=1
everything ⇒ δ9
```

## Soft edges (→) — the negotiable ones

```
γ2 → α10        WKB illuminates Airy, but Airy is exact; not required
δ4 → δ1         Thom's classification is an INTERPRETATION of q, not a prerequisite (see F3)
RMT → β6        Tracy–Widom is reachable with NO random matrix theory (see F4)
α3 → α6         relaxation oscillation makes canards vivid, isn't needed to define them
β4 → β5         Kramers is the physical picture; the operator identification stands without it
γ10 → γ4        Padé is how you FIND singularities numerically, not what they ARE
δ5 → δ9         the two-ladders split is a coda, not load-bearing for the main claim
```

---

## Structural findings — the part that actually matters

### F1. There are two nearly-independent trunks that meet exactly once
Trunk α (geometry) and trunk γ (asymptotics) share **no** hard edge until the junction node δ6.
γ1 — "asymptotic series diverge" — depends on *nothing* in this program.

⇒ The resurgence thread does **not** have to wait for eight episodes of dynamics. It can run as a
parallel strand from episode 1 and converge later. This is the "two stories that meet" structure,
and here it is forced by the mathematics rather than chosen for effect.

### F2. δ6, the v_k coefficient ladder, is the sole junction
Everything geometric *produces* the v_k; everything analytic *consumes* them. The whole program
narrows to a single sequence of numbers and then re-expands.

⇒ Strong candidate for the visual spine of the entire course: one recurring object, present in
every module, meaning something different each time. Also the natural pilot — it's where the two
trunks touch, so building it first tests both halves of the architecture.

### F3. You may not need catastrophe theory to explain the catastrophe ladder
The ladder here is indexed by **turning order q in a second-order ODE** (potential ~Y^q). Thom's
A_k classification is an interpretation laid over that, not a prerequisite for stating it.

⇒ δ4 → δ1 is SOFT, and teaching δ4 first is likely a large wasted tax. More importantly it is
*actively harmful* here: the program's own positioning insists the swept multicritical family
(second-order ODE, multicritical swept turning potential) is **not** the genuine n-th-order
Berry–Upstill diffraction catastrophe. Leading with Thom builds in precisely the conflation the
papers work to prevent. **Distortion-ledger consequence: any swallowtail surface shown early
must be labelled as an analogy, not the object.**

### F4. Tracy–Widom is reachable with no random matrix theory — so RMT becomes a payoff
Via the stochastic Airy operator (RRV), TW arrives as "the ground-state statistics of a random
Schrödinger operator". The random-matrix identification is a *reveal*, not a gateway.

⇒ Converts a heavy prerequisite into a climax: *the law you just derived from a neuron is the
law of the largest eigenvalue of a random matrix.* Do not spend an episode on RMT up front.

### F5. Blow-up's hard prerequisite is a FAILURE, not a fact
α7 requires α5 — specifically that normal hyperbolicity **dies** at the fold. Present blow-up
before that failure is felt and it reads as unmotivated machinery.

⇒ The blow-up module must open on the breakdown, not the construction. "Here is where every tool
we have stops working" ⇒ *then* the sphere.

### F6. In trunk γ the linchpin is γ6, and it is one identity
The chain γ4→γ5→γ6 is: the Borel transform has a cut ⇒ the sum depends on which side you pass ⇒
**that ambiguity is exactly e^{−A/x}**. Everything downstream (trans-series, Stokes, resurgence)
is bookkeeping on that identity.

⇒ γ6 is the single highest-value beat in the whole course. It also fixes an ordering that is easy
to get wrong: **trans-series before Stokes.** You need an exponentially small term to exist (γ7)
before "how it switches on" (γ8) is even a question. Most expositions invert this.

⇒ Also forces a distortion fix: if γ4 draws singularities as isolated *points*, γ5 has nothing to
be ambiguous about. **The cut must be visible from first appearance** — which is exactly the
defect already flagged in `W_make_figures.py` and now the ledger's first entry.

### F7. Phenomenon/mechanism splits — show early, explain late
| node | showable at | explainable at |
|------|-------------|----------------|
| α6 canards | α2+α4 (visually striking immediately) | after α7 (rigorous treatment needs blow-up) |
| γ1 divergence | episode 1 (just plot partial sums) | after δ6/δ7 (why *these* coefficients) |
| β6 Tracy–Widom | as a histogram, very early | after β5 |
| δ4 caustics | early, as physical pictures | only alongside δ5, or not at all |

⇒ Four legitimate early hooks that cost nothing in prerequisites.

---

## What this constrains, and what it leaves free

**Forced:** α5 before α7 · α7 before α8→α11 · γ7 before γ8 · β6 before δ8 · δ6 before the
resurgence payoff · the cut visible from γ4.

**Free (research decides):** whether to open on the neuron (concrete, α-first) or on the
divergent series (abstract, γ-first) or interleave from the start; how much β to teach at all;
whether δ4 appears; where TW's reveal lands; how many modules.

**My prior before the research lands**, flagged so it can be overturned: interleave α and γ from
episode 1, converge at δ6, and pilot on the Borel plane (γ4–γ6) — because it is the hardest
node, has the worst existing visualisation, and stresses the architecture more than any geometric
module would. Cheap-geometric-first is the safer play; I currently favour the harder one, but
weakly.
