# Future directions — a sourced map of the recent literature (2024–2026)

*What this is: leads harvested directly from recent papers and reviews, organised by
direction, each with the field's own stated open problem and citations. What this is
**not**: a list of vetted-open problems. Three rounds of diligence (c₀, the
quasipotential bridge, degenerate-noise folds) showed I cannot reliably certify
"open" from outside the expert network. So treat every entry as **a reading lead and
a question for Popović**, not a green light. The adjudication of "genuinely open AND
reachable in your timeframe" is exactly what the supervisor is for; this gives you
the sourced shortlist to make that conversation efficient.*

---

## Start here — reviews (they list open problems for you)

- **Wechselberger, *GSPT Beyond the Standard Form* (Springer, 2020)** + the
  three-timescale line below — the canonical modern GSPT reference; frames the
  non-standard / oscillatory frontier.
- **Review on mixed-mode-oscillation dynamics in neuron models & neural networks**,
  *EPJ Special Topics* (2025), [10.1140/epjs/s11734-025-02041-9](https://link.springer.com/article/10.1140/epjs/s11734-025-02041-9).
- **Asymptotic return maps & geometric dynamics in multi-timescale neuronal systems:
  Decoding MMOs**, *EPJ Special Topics* (2026), [10.1140/epjs/s11734-026-02224-y](https://link.springer.com/article/10.1140/epjs/s11734-026-02224-y).
- **Global phase-space approach to rate-induced tipping: a brief review**, *Chaos*
  (2025), [aip 3344946](https://pubs.aip.org/aip/cha/article/35/4/043139/3344946).
- **Canards, Folded Nodes & MMOs in Piecewise-Linear Slow-Fast Systems**, *SIAM
  Review*, [10.1137/15M1014528](https://epubs.siam.org/doi/10.1137/15M1014528).

---

## Directions (with the stated gap, sources, toolkit fit, honest ceiling read)

### 1. Three-timescale GSPT & folded *limit-cycle* manifolds
**Stated gap (Wechselberger's own framing):** "the theory for *oscillatory*
multiple-timescale systems which possess a limit-cycle manifold is **less developed**,
particularly in the non-normally-hyperbolic regime." Blow-up of folded limit-cycle
manifolds in 3-timescale systems is just being opened.
**Sources:** [arXiv:2208.01361](https://arxiv.org/abs/2208.01361) (blow-up for folded
limit-cycle manifolds, 3 time-scales); [arXiv:2009.10583](https://arxiv.org/abs/2009.10583)
(parametrisation method, ≥3 timescales); [SIAM 22M1477477](https://epubs.siam.org/doi/10.1137/22M1477477)
& [Math.Comp.Sim. 2025](https://dl.acm.org/doi/abs/10.1016/j.matcom.2025.01.003)
(3-timescale Hodgkin–Huxley).
**Fit:** directly your GSPT/blow-up toolkit; deterministic ⇒ theorem-grade and
supervisable by Popović. **Ceiling read:** of everything here this is the most likely
to carry a genuine A1 theorem *and* be in your wheelhouse — a frontier the leaders
themselves call underdeveloped. Verify the specific sub-problem with Popović.

### 2. Sharp characterisation of SISR / coherence resonance beyond idealised limits
**Stated gap:** analytical SISR results are "often restricted to idealised limits
(weak noise, infinite timescale separation), limiting applicability to realistic
parameter regimes." SISR "remains analytically and numerically challenging."
**Sources:** [arXiv:2510.22848](https://arxiv.org/abs/2510.22848) (SISR via
physics-informed ML, Savaliya–Yamakou 2025); [PRE 100.022313](https://link.aps.org/doi/10.1103/PhysRevE.100.022313)
(CR controlled by SISR, multiplex net). Background: Muratov–Vanden-Eijnden;
DeVille–Vanden-Eijnden.
**Fit:** *this is your model* (stochastic FHN) and your FW/Kramers + numerics toolkit.
**Ceiling read:** a sharp finite-ε / finite-noise SISR result would be real, but the
risk is it slides into "incremental sharpening" like c₀ — pin the precise open
quantity with Popović before committing.

### 3. Non-Gaussian noise (α-stable / fractional) in slow-fast averaging & LDT
**Stated gap:** "very little attention has been given to the case of **α-stable**
noise forcing" in stochastic averaging; large deviations under fractional Brownian
motion is recent and partial.
**Sources:** [SIAM MMS 140990632](https://epubs.siam.org/doi/10.1137/140990632)
(α-stable averaging); [arXiv:2007.08408](https://arxiv.org/abs/2007.08408) (weak
averaging, stable processes); [arXiv:2210.03678](https://arxiv.org/abs/2210.03678)
(LDT, fractional BM). Caution: Lévy *exit* problems are already active (heavy-tailed
Kramers), so the open part is the slow-fast/canard-specific behaviour, not Lévy LDT
per se.
**Fit:** partial — needs new (Lévy/jump or fBm) machinery beyond your current base.
**Ceiling read:** genuine gap, higher ceiling, but the steepest learning cost.

### 4. Rate-induced tipping + early-warning signals
**Stated open challenge:** "Devising a versatile early-warning signal capable of
detecting global bifurcation-induced tipping *as well as* rate-induced tipping …
remains an open challenge." Plus the rate-vs-EWS tension ("early warnings are too
late when parameters change rapidly").
**Sources:** [Chaos 2025 review](https://pubs.aip.org/aip/cha/article/35/4/043139/3344946);
[Sci. Rep. 2025](https://www.nature.com/articles/s41598-025-06525-5); [Proc. R. Soc.
A 2025](https://royalsocietypublishing.org/rspa/article/481/2321/20250405/234263/);
[arXiv:2508.19655](https://arxiv.org/abs/2508.19655) (Koopman EWS).
**Fit:** your FW/large-deviations + numerics; big, well-funded, applied audience
(climate/ecology). **Ceiling read:** high visibility; "sharp constant" versions risk
incrementalism, but a genuinely new EWS *theorem* could be A1. Crowded — pick a
precise corner.

### 5. Quasipotential / minimum-action computation (higher-D, non-gradient)
**Stated gap:** "calculating the quasipotential in higher dimensions is numerically
challenging"; method-development is active (ML, Ritz, gMAM).
**Sources:** [Chaos 2019 (Grafke–V-E)](https://pubs.aip.org/aip/cha/article/29/6/063118/1027313)
(numerical rare events via LDT — the survey); [arXiv:2306.11418](https://arxiv.org/abs/2306.11418)
(LD prefactors via ML); [arXiv:2001.10759](https://arxiv.org/abs/2001.10759) (Ritz).
**Fit:** your Phase 1–3 quasipotential-solver work plugs straight in.
**Ceiling read:** methods-tier (high-A2), as already assessed — list it for
completeness, not as a 90+ route.

### 6. Networks / mean-field of noisy slow-fast excitable units
**Theme:** noise- and coupling-induced oscillations, mean-field reduction, finite-size
effects.
**Sources:** [Comm. Math. Phys. 2019](https://link.springer.com/article/10.1007/s00220-019-03641-y)
(emergence of oscillations, excitable mean-field slow-fast — Luçon–Poquet line);
[arXiv:2103.04492](https://arxiv.org/abs/2103.04492) (phase reduction of coupled noisy
oscillators); [arXiv:2503.12596](https://arxiv.org/abs/2503.12596) (coupled-FHN
canards). **Fit:** Floquet/phase-reduction + FW. **Ceiling read:** high payoff but
crowded and PhD-arc-sized; better as a thesis direction than a single A1 result.

### 7. Data-driven / ML discovery of slow variables & effective dynamics
**Theme:** learning slow manifolds / reduced SDEs from trajectory data.
**Sources:** [arXiv:2205.04151](https://arxiv.org/abs/2205.04151) (learning effective
dynamics); [JNS 2022](https://link.springer.com/article/10.1007/s00332-022-09808-7)
(slow-variable discovery via NN). **Fit:** numerics-heavy. **Ceiling read:**
methods/ML-tier; pairs well with (5) but not an A1 theorem route.

---

## How to use this

1. **Read the four reviews first** — they contain the field's own open-problem lists
   and will sharpen which corner of (1)–(4) actually has a theorem-sized hole.
2. **Map each to your toolkit:** strongest fits are **(1) three-timescale/folded
   limit-cycle GSPT** and **(2) SISR sharpening** — both use what you already have,
   and (1) is the one a leader explicitly calls "less developed."
3. **Take a 2–3 item shortlist to Popović** with the direct question: *which of these
   has a genuinely open, A1-ceiling sub-problem that's reachable in my timeframe?*
   That is the step three rounds of solo scoping have shown is irreplaceable — he can
   see the unpublished/in-progress landscape that searches cannot.

**Honest bottom line:** I can give you the field's stated frontiers and the sources;
I cannot certify which is open. (1) and (2) are where I'd point first on fit and
plausible ceiling — but the certification is the supervisor's to make.
