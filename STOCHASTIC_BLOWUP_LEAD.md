# The strongest lead: geometric blow-up for slow-fast **SDEs**

*Found not by generating a bet, but by mining expert-curated open-problem lists —
which is the methodological point of this note as much as the lead itself.*

---

## Why this one is different from the deflated bets

1. **It's community-flagged, not my guess.** The Oberwolfach/BIRS 2022 workshop
   *Topics in Multiple Time Scale Dynamics* lists as a broad challenge: *developing
   geometric techniques (à la **blow-up**) to deal with singularities of slow-fast
   **PDEs and SDEs**, complementing the theory currently elaborated mainly for ODEs.*
   ([BIRS 22w5057 report](https://www.birs.ca/workshops/2022/22w5057/report22w5057.pdf)).
2. **It's methodological / fresh-angle** — developing a *desingularization method*
   for the stochastic case — which is the kind of novelty (a new route, not a pinned
   constant) you said you wanted.
3. **It's your exact toolkit fused:** geometric blow-up (you've done it — fold
   normal form, c₀) **+** noise / Freidlin–Wentzell / sample paths (you've done it).
   Few people have *both* fluently; that's a real comparative advantage.

## What's verified about the state of the art

- **Deterministic blow-up is highly developed and actively extended** — and the
  extensions are being pushed by *Kuehn's group*: cusp (Broer–Kaper–Krupa), the
  hyperbolic umbilic ([arXiv:2202.01662](https://arxiv.org/abs/2202.01662), 2024),
  **fast-slow PDEs via Galerkin** (Engel–Kuehn, [arXiv:2207.06134](https://arxiv.org/abs/2207.06134)),
  folded limit cycles ([Jelbart–Kuehn–Kuntz 2024](https://arxiv.org/abs/2208.01361)).
  There's even a survey ([arXiv:1901.01402](https://arxiv.org/abs/1901.01402)).
- **The SDE side is conspicuously behind the PDE side.** Both targeted searches
  returned essentially deterministic results; "specific results about SDEs / noise
  were not prominent." The *PDE* half of the Oberwolfach challenge is being built
  (Engel–Kuehn); the *SDE* half appears genuinely less developed.

## The honest caveats (so this doesn't become deflation #4)

- **Berglund–Gentz already do "stochastic dynamics near slow-fast singularities"** —
  their book does sample-path / covariance-tube bounds near folds, pitchforks, Hopf,
  bifurcation delay; BGK do folded nodes. So the *phenomenon* (noise near a
  singularity) is not untouched. **The open part is the *method*: a systematic
  geometric blow-up / desingularization of the *SDE* on the blown-up space (tracking
  the noise / covariance through the rescaling), as opposed to Berglund–Gentz's
  direct estimates.** The contribution has to be "what blow-up buys you beyond
  sample-path bounds," and that needs to be articulated sharply.
- **It is almost certainly on Kuehn/Engel's radar.** Their deterministic-and-PDE
  blow-up momentum points straight at the SDE case; it would be naïve to assume it's
  untouched in their pipeline. Same epistemic ceiling as everything else: I can't see
  unpublished work.
- Recent adjacent activity exists, e.g. *Rate and Bifurcation Induced Transitions in
  Asymptotically Slow-Fast Systems* ([SIAM ADS 24M1632000](https://epubs.siam.org/doi/10.1137/24M1632000)).

## Why this is nonetheless the best position you've been in

The lead is **community-flagged + fresh-angle + exact-toolkit-fit**, and the world
expert on the blow-up-extension programme is **Christian Kuehn — your supervisor's
direct collaborator.** That converts the usual "is it open?" uncertainty into a
single, sharp, well-sourced question Popović can answer *and* act on:

> *"The Oberwolfach 2022 list flags geometric blow-up for slow-fast SDEs as open,
> complementing the ODE/PDE theory (which Engel–Kuehn are building). Is a stochastic
> blow-up / desingularization framework — or even one worked singularity (the
> stochastic fold via blow-up, tracking covariance through the rescaling) — open, or
> is it in Christian's pipeline? It's exactly my blow-up + Freidlin–Wentzell toolkit,
> and I want an A1-ceiling problem."*

If open: methodological, high-ceiling, toolkit-matched, supervisor-backed into the
Kuehn programme. If in-pipeline: he can get you *in* rather than scooped.

## The meta-point (your actual question)

You were right that there are better ways to find avenues than my generate-and-check.
**Mining expert open-problem lists just out-performed every bet I invented** — it
produced a higher-signal, fresh-angle, toolkit-matched lead in one pass. Other
methods worth running next: citation-frontier mining (what cites Berglund–Gentz /
Jelbart–Kuehn–Kuntz and says "still open"), and assumption-relaxation mining (the
hypotheses of the key theorems are a map of their gaps). All of them still terminate
at the same place — Popović vetting — but they start you from expert-flagged ground
instead of my guesses.

## Other expert-flagged items from the same pass (lower priority)

- **3D piecewise-linear canards / folded singularities** — "very few results in the
  3D case" ([Desroches et al., SIAM Review](https://arxiv.org/abs/1606.03643)).
- **Higher-dimensional / discrete-time canard problems** — flagged "still to be
  explored" in the BIRS report.
- **Fast-slow PDE blow-up** — active (Engel–Kuehn), so less open than the SDE side.
