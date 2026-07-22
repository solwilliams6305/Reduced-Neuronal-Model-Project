# Prompt — Folded limit cycle, Part 2 (worked model) + Part 3 (α=1 coupled core)

*Self-contained brief for a fresh Opus chat. The leading-order theory (Channels A & B,
frozen phase, normal form) is **already done and validated** — this brief extends it to
(2) a real worked model and (3) the genuinely-new coupled rotating-phase core. **Read §0
and the prior files before deriving — one of them overturns a natural-but-wrong idea you
will otherwise re-derive.***

---

## 0. Orientation, what to read, and the gates

**State of play.** This project derived the stochastic blow-up of the fold
(`CANARD_BLOWUP.md`) and folded node (`MMO_NOISE.md`), then the **folded limit cycle**
(`FOLDED_CYCLE_NOISE.md`): the one singularity BG/BGK never noised. Done & validated at
leading order, normal form, frozen phase:
- `η = σ/√ε₂` (effective noise in the JKK K₂ chart),
- `σ_*^A(θ) = C_q √ε₂ √(a(θ)c(θ)/b(θ))` (amplitude escape; geometry collapse + falsification control + C_q inheritance),
- `σ_*^B ≈ 2π√3 √(ac/b) √(ε₂/|ln ε₂|)` (phase diffusion; log law validated),
- both at the **same leading order**, with a **geometry-independent A-vs-B race** set by ε₂.

**Read first (these are the engine and the guardrails):**
- `FOLDED_CYCLE_NOISE.md` — the leading-order result you are extending. (Normal form §1–2; η §3; σ_*^A §4; the race §9; the four open items §11.)
- **`ALPHA1_METHODOLOGY.md` — READ THIS CAREFULLY. It refutes the obvious α=1 idea.**
  A natural conjecture — "fast rotation makes the cycle escape at its *weakest* phase
  (min G)" — was tested numerically (`outputs/alpha1_test3.py`) and is **WRONG**:
  `η_*` is flat in the rotation rate (ratio 1.01, min predicted 0.50). Reason: canard
  escape needs a finite local time `T_loc`; a fast phase dwells in the weak window only
  `δ/ω ≪ T_loc`, too briefly to escape there. **The α=1 case AVERAGES, it does not
  minimise.** Do not pursue the min-barrier / band-edge / localization route.
- `CANARD_BLOWUP.md` (§10 = the full-model regime-map template you will copy for Part 2);
  `VDP_CROSSMODEL.md` (Van der Pol = the natural Liénard testbed); `MMO_K2_ROUTE_AB.md`
  (K₂ inner-chart return map = Path-A precedent); `TONIC_CMID_AIRY.md` + `TONIC_PHASE.md`
  (the Airy/Riccati linearization and the iPRC phase-diffusion engine — both needed in Part 3).
- Scripts to extend: `folded_cycle_normal_form_map.py`, `folded_cycle_phase_diffusion.py`,
  `alpha1_test3.py`.

**Step 0 — gates before investing (this project has been burned by skipping these):**
1. **The Popović email.** Is a stochastic / blow-up folded limit cycle in Kuehn's
   pipeline, or free? Evidence (`EVIDENCE_KUEHN_PIPELINE.md`) says free, but unpublished
   work is invisible. Gate the headline novelty claim on his reply.
2. **Read JKK 2024 (arXiv:2208.01361) §5.1** in full — the periodically-forced **Liénard**
   worked example is Part 2's testbed; its explicit equations were beyond reach last pass
   (pull from the v3 PDF, pp. ~22–26).

---

## PART 2 — A worked model (periodically-forced Liénard / Van der Pol)

**Why:** everything so far is normal-form with a *borrowed* C_q. A real model (i) confirms
the law survives finite ε, (ii) fixes a genuine C_q (expect ≈ 8–10 as for full FHN — and
the "Channel B wins for full models" headline depends entirely on this number), (iii)
makes it a credible standalone paper rather than only a chapter.

**The key subtlety (the real work):** the bridge from abstract `a(θ),b(θ),c(θ)` to a
concrete system. You must **read the folded-cycle normal-form coefficients off the actual
Liénard equations** — i.e. identify, for the chosen model, what plays the role of the
fast radius `r`, the phase `θ`, the slow drift `y`, and the periodic `a,b,c`. JKK §5.1
gives the reduction; VdP (`VDP_CROSSMODEL.md`) is the cleanest relaxation-oscillator
Liénard to start with (slowly force a parameter through the relaxation-cycle fold).

**Tasks:**
1. Take the JKK §5.1 Liénard system (or forced VdP); add **degenerate noise on the fast
   variable** (project convention). Locate the **folded limit cycle** (the fold of the
   limit-cycle manifold as the slow forcing drifts).
2. Confirm the blow-up structure on the real model: the `η = σ/√ε₂` scaling and the K₂
   canonical Riccati. Read off `a(θ),b(θ),c(θ)`.
3. **Reproduce `σ_*^A = C_q √ε₂ √(ac/b)` on the real system** — the regime map à la
   `CANARD_BLOWUP.md` §10 / `folded_cycle_normal_form_map.py`. Use the same geometry-
   collapse + **falsification control** discipline. Fit the model **C_q**.
4. (If reachable) confirm the phase-modulation `G(θ)` and the A-vs-B race on the model.

**Deliverable:** `FOLDED_CYCLE_LIENARD.md` — the law validated on a real system, with a
measured C_q. ~2–4 weeks of work; do this first (it grounds everything and is lower-risk).

---

## PART 3 — The α=1 coupled core (the genuinely-new mathematics)

**What's settled:** frozen-phase α=2 (done); and `ALPHA1_METHODOLOGY.md` shows the α=1
escape **averages** the phase modulation (min refuted, rotation ~washes out). What's NOT
done: the actual derivation of the averaged law, its prefactor, the fate of the log, and
the true A–B fusion. This is the part that turns "two channels compared" into "a coupled
theory" — the deep moat.

**Methodology (the routes that survived the test):**
- **Stochastic averaging (Khasminskii) / Kapitza effective-barrier** — average the
  rotating-phase inner SDE over the fast phase to get an effective autonomous canard with
  an *averaged* barrier. **Vindicated** by the test.
- **Fourier-in-phase Fokker–Planck hierarchy** — expand the (R,θ) density in harmonics
  `e^{2πinθ}`; the 1/R² phase noise damps mode n at rate `n²η²/R²`, so n≠0 modes die →
  controlled truncation to the n=0 (phase-averaged) radial equation + computable
  corrections. This is the **Path-A backbone** that makes the averaging rigorous.
- **Riccati → Schrödinger linearization** — `R = −u′/u` ⇒ `u″ = a(θ)Y·u`; with rotating θ
  and the ramp `Y = Y₀ − T` this is a modulated-Airy / Hill–Airy equation, escape = first
  node. Use it for the **O(1) prefactor and the log**. NB the Lyapunov exponent here
  **averages** (it does *not* pick a band edge — that reading was the refuted one); expect
  the `|ln ε₂|` to come from the **ramp/turning-point (Airy)** structure, not localization.

**Tasks (Path B first — heuristic + numerics, the house style):**
1. **Derive the averaged Channel-A law.** Average the inner escape over the fast phase;
   produce `σ_*(α=1) = C_q √ε₂ · ⟨√(ac/b)⟩_eff` and **pin which average** `⟨·⟩_eff` is
   (the test says NOT min and roughly the typical/action-weighted value — derive it; the
   correct weighting is set by where in the (Y,θ) plane escape actually occurs).
2. **Numerically validate it.** Extend `alpha1_test3.py` to a proper `σ_*` measurement and
   check it against the derived average on a model where the candidate averages
   (arithmetic, harmonic, action-weighted) are well separated — a clean discriminator,
   exactly as the min-vs-mean test was built.
3. **The A–B fusion.** With the phase both rotating and diffusing, does the
   geometry-independent race survive? Does Channel B's `|ln ε₂|` enhancement persist, or
   does rotation alter it? This coupled question is the genuine α=1 contribution.
4. **Fourier-hierarchy control** (Path A backbone): make the averaging systematic and
   bound the n≠0 corrections.
5. **Rigorous inner solution** (Path A, stretch / months): matched asymptotics through the
   JKK charts, the noisy Weber/Airy inner — the theorem-grade endpoint.

**Deliverable:** `FOLDED_CYCLE_ALPHA1.md` — the averaged α=1 law (derived + validated), the
fate of the race and the log under coupling, and an honest Path-A/Path-B split.
~1–3 months (Path B); longer for Path A.

---

## Methodology discipline (non-negotiable — it's why this project's results hold)

- **Derive, then numerically validate with a falsification control** (re-run with a
  *wrong* normalisation and show the collapse breaks — see `FOLDED_CYCLE_NOISE.md` §6).
- **Path B (heuristic + numerics) before Path A (rigour).** State which you're in.
- **Test conjectures before writing them up.** The α=1 min idea was natural and wrong; one
  clean simulation caught it. Build the discriminating experiment first.
- **Don't overclaim.** Distinguish "leading order + normal form" from "rigorous" from
  "validated on a real model." Position explicitly against BG (fold), BGK (folded node),
  JKK (deterministic folded cycle).
- **Re-collision check:** the moment a result reduces to the plain fold or folded node,
  it's BG/BGK's — keep the fast-oscillation phase (and its coupling) central.

## Suggested order

Part 2 first (grounds the law, lower risk, ~weeks) → Part 3 Task 1–3 (the averaged α=1 law
+ the fusion, the deep result, ~months) → Part 3 Task 4–5 (Path-A rigour, stretch). The
Popović email runs in parallel and gates the novelty framing throughout.

## References

- **JKK 2024** — folded limit cycle blow-up, arXiv:2208.01361 (scaffold; §5.1 Liénard).
- Berglund–Gentz 2006 (fold); BGK 2015, arXiv:1312.6353 (folded node); Krupa–Szmolyan 2001;
  Wechselberger 2005 (folded-node K₂ inner solution — the Path-A model).
- Khasminskii (stochastic averaging); the Kapitza effective-potential idea (high-frequency
  averaging).
- Project: `FOLDED_CYCLE_NOISE.md`, `ALPHA1_METHODOLOGY.md`, `CANARD_BLOWUP.md`,
  `VDP_CROSSMODEL.md`, `MMO_K2_ROUTE_AB.md`, `TONIC_CMID_AIRY.md`, `TONIC_PHASE.md`,
  `EVIDENCE_KUEHN_PIPELINE.md`, `FOLDED_CYCLE_NOISE_PROMPT.md`; scripts
  `folded_cycle_normal_form_map.py`, `folded_cycle_phase_diffusion.py`, `alpha1_test3.py`.

---

*First actions: send the Popović email; read JKK §5.1; read `ALPHA1_METHODOLOGY.md` so you
don't re-derive the refuted min-barrier route. Then Part 2.*
