# Noise error bounds across regimes — the paper's spine

*The central deliverable: for the stochastic FitzHugh–Nagumo model, when (and by
how much) does the slow-manifold reduction fail, regime by regime? This table is
the answer. The individual derivations live in `regime-tests/` and are
consolidated in `DERIVATION_SYNTHESIS.md`; this document is the map they serve.*

The question (from `PROJECT_CONTEXT.md`): *when does the reduction give
qualitatively wrong predictions, and what is the precise mechanism?* The answer
is a regime-indexed threshold `σ_crit(ε)` — above it the reduction fails — plus,
for each regime, **which O(1) constant is closed-form and which is irreducibly
numerical.** That sharp/numerical demarcation is part of the result, not a
disclaimer.

---

## The map

Bifurcation order in `I` (with the 3-D folded-node extension last).

| Regime | Control | Mechanism (why the reduction fails) | `σ_crit(ε)` with constant | Constant status | Closest prior art |
|---|---|---|---|---|---|
| **Excitable (rest)** | `I < I_SNIC` | Kramers escape from the rest focus over the fold barrier | `√(ε/C)`, `C=K_fold/A₀`; near-fold `σ_crit = √(4π g ln2)·√ε` | **sharp** — `K_fold=4/3`, `A₀=1/π` derived; `g` = fold drift | Berglund–Gentz (fold/Kramers) |
| **SNIC** | `I → I_SNIC` | barrier `ΔU → 0`; any noise tips | `σ_crit → 0` | degenerate (qualitative) | standard |
| **Canard (autonomous)** | `I` near lower Hopf, at the fold | accumulated Brownian variance across the canard window | `C_q·√ε·λ^{1/2}`, `C_q = √(4π ln2) ≈ 2.95` | **sharp** (leading order); quantile law `C_q(p)=√(−4π ln(1−p))` | Berglund–Gentz–Kuehn (noisy canards) |
| **Canard (explosion/ramped)** | ramped through `I_H` | small-cycle washout under ramp | `~ ε^{3/4}` | **cited** (BGK 2012); not re-derived | BGK 2012 |
| **Tonic spiking** | `I > I_Hopf` | phase diffusion on the limit cycle | `CV = σ·A_mid`, `A_mid=√(c/π²)`; `σ_crit ~ √(ε·T_cycle)` | **scaling + structure**: `c≈1.55`, mechanism + Airy inner; number open | phase-diffusion; Sacré–Franci (PRC) |
| **Resonator** | `I` just below Hopf (stable spiral) | 2-D quasipotential escape near the spiral; spike-direction commitment | Arrhenius `k~A e^{−B/σ²}`, `B→δ²` corner (floor-dominated in window); commitment angle `−1.06 rad` | **mixed**: `δ²` corner derived; measured exponent floor-dominated; angle **numerical** | subcritical-Hopf normal form; gMAM |
| **Folded-node MMO** | 3-D FHR | noise dissolution of the mode-locked staircase; noise changes the SAO count | `σ_pq ~ σ_*·q^{−α/γ}`; `a_min = a_max e^{−π²f}` | **scaling**; constants global (`Λ(μ)`), **numerical** | **Berglund–Gentz–Kuehn (1312.6353) — owns this row** |
| **Tonic onset** | non-perturbative | Freidlin–Wentzell cycle escape | `σ_onset ≈ 0.015–0.02` | **numerical** (PhD-scope) | Berglund–Gentz (cycle escape) |
| **Van der Pol** | cross-model check | same fold normal form (`b=0`) | same `C_q=√(4π ln2)`; tonic `c_vdp≈0.78` | universal half rides the canard row | — |

**One sentence:** the reduction fails above `σ_crit(ε)`; that threshold is
`√ε × (regime geometry)` in every regime; the geometry factor is a closed-form
constant where the mechanism is a **quasipotential/barrier integral** (excitable,
canard, VdP) and an irreducibly **numerical/global** constant where it is a
**minimum-action path or a global return** (resonator angle, folded-node MMO).

---

## What is genuinely first about this (and what is not)

**First-able:**

1. **Completeness for one model.** Berglund–Gentz cover fold/canard escape;
   Berglund–Gentz–Kuehn cover the folded-node MMO; Lindner/Pikovsky cover
   coherence resonance; Sacré covers PRCs. **No single source assembles
   `σ_crit(ε)` across *all* FHN regimes** — excitable rest-state, SNIC, canard,
   tonic phase-diffusion, resonator, and MMO — as one quantitative map answering
   "when does the reduction fail." The synthesis is the contribution.
2. **The sharp constants that make rows quantitative.** Turning scalings into
   numbers: `C_q = √(4π ln2)`, the excitable `K_fold=4/3, A₀=1/π`, the unified
   `σ_crit=√(4π D ln2 ε)`. Scaling-law papers characteristically leave these
   open; pinning them is the upgrade from "scaling theory" to "asymptotic theory."
3. **The honest sharp/numerical demarcation** — a *criterion* for which
   cross-regime noise constants can be closed-form at all (quasipotential ⇒ yes;
   path-geometry/global-return ⇒ no). This is the organizing result and the part
   that owes nothing to prior art.

**Not first / must defer:**

- The **folded-node MMO row is BGK's** (1312.6353): the SAO-count/global-return
  interplay and noise-changes-the-count are theirs. Cite, position against, claim
  nothing there until read in full.
- The **rigour level.** This map gives *thresholds and mechanisms*. Berglund–
  Gentz's standard is stronger — *probabilistic sample-path error bounds*
  ("with probability ≥ 1−e^{−κ/σ²} the path stays in a neighbourhood of width
  …"). The fold/canard rows can be upgraded to that rigour using their template;
  the resonator/MMO rows need the global constants first.

---

## The honest framing for the thesis

> A complete, quantitative `σ_crit(ε)` map across the dynamical regimes of the
> stochastic FitzHugh–Nagumo model — the precise noise level at which the
> slow-manifold reduction fails in each — with the O(1) constant **derived in
> closed form where the mechanism is a quasipotential barrier (excitable, canard,
> Van der Pol)** and **identified as irreducibly numerical where it is a
> minimum-action path or global return (resonator, folded-node MMO)**. The
> closed-form rows sharpen known scaling laws to asymptotic constants; the
> demarcation between closed-form and numerical is itself the organizing result;
> the folded-node-MMO row is positioned against Berglund–Gentz–Kuehn.

That is a legitimate, defensible, and — as a *complete quantitative map for one
model* — plausibly first contribution, without dressing any single row as a new
phenomenon.

---

## Literature-check verdict (all rows, June 2026)

Every row's **mechanism is prior art**, in different literatures:

| Row | Prior art (mechanism) | Sharp constant |
|---|---|---|
| Excitable | coherence resonance — Pikovsky–Kurths 1997; **Lindner–Schimansky-Geier (analytical FHN escape)**; Lindner 2004 review | likely already pinned by Lindner–S-G |
| SNIC | type-I / SNIC ISI inverse-√; noise near SNIC (Adler, stochastic sensitivity) | qualitative (`→0`) |
| Canard / MMO | **BGK 2012 (read in full)** — Thm 4.4 *is* `a_min`; Cor 6.3 is the SAO-hiding threshold | `c₀∈[π/4,1]` **open** |
| Tonic spiking | limit-cycle phase diffusion (textbook; rigorous arXiv:1512.04436); CV via adjoint/iPRC standard | `c` is an FHN number; Airy = Sacré's area |
| Resonator | noise escape near subcritical Hopf / quasipotential / Wentzell–Freidlin (studied in FHN) | `δ²`/angle project-specific; angle numerical |
| Tonic onset | Freidlin–Wentzell cycle escape (classical) | `σ_onset` measured |

**Conclusion:** no row is a new mechanism or threshold; each is a re-derivation
within an existing framework. The genuinely first-able content is **not in the
rows** but cross-cutting:

1. **Pin `c₀ ∈ [π/4, 1]`** (BGK left it as a range) — the one confirmed-open sharp
   constant; the Weber/PCF machinery is the tool.
2. **The closed-form-vs-numerical demarcation** — a *principle* (barrier ⇒ closed,
   path/global ⇒ numerical). **[CORRECTED — folklore, not original.]** A literature
   check shows its precise form is the known integrability criterion for
   quasipotentials: a closed-form (quasi)potential exists ⟺ the Wentzell–Freidlin
   Hamiltonian is completely integrable (classical Hamilton–Jacobi/Liouville
   theory; explicit in the quasipotential literature; Maier–Stein the canonical
   numerical example). Genuine *organizing/expository* value for the thesis, but
   **not** an original contribution — novel only via a specific non-obvious
   adjudication, which is not in hand.
3. **The Airy PRC duality** — distinctive but in Sacré–Franci's area; gated on
   reading Sacré's 2013 thesis.

Plus the **cross-regime completeness** as consolidation value (review/thesis-
chapter tier, not a discovery). Honest thesis framing: *a consolidation of when
the FHN reduction fails across regimes within the coherence-resonance / BGK /
Wentzell–Freidlin canon — organized by the (known) integrability demarcation, and
contributing the sharpening of BGK's `c₀` (= π²/16, verified) and the Airy-PRC
structure if it is not already in Sacré 2013.*

## The next moves that most strengthen the spine

1. **Fill the quantitative rows to B–G rigour** (canard, excitable): convert
   `σ_crit` thresholds into probabilistic sample-path bounds using the
   Berglund–Gentz template — this makes "error bounds" literal, not just
   thresholds.
2. **Validate the sharp constants numerically** in one clean sweep: the canard
   quantile law `C_q(p)=√(−4π ln(1−p))` (threshold sweep), the excitable
   `√(4π g ln2)` prefactor (failure-boundary fit), the κ drift. One figure per
   row turns the table from derived to validated.
3. **Position the MMO row against BGK** explicitly (read 1312.6353) and state the
   map's MMO entry as "consistent with / a consolidation of" their results, with
   `a_min = a_max e^{−π²f}` as the project's compact restatement.
4. **State the demarcation as a short standalone claim** (quasipotential ⇒
   closed-form; path/global ⇒ numerical) — the one genuinely new, transferable
   sentence, and the intellectual core of the map.
