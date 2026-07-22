# A real (even if small) contribution in "blow-up for slow-fast SDEs"

*Goal: something genuinely new, not a rehash. This note is honest about the narrow
seam, ranks concrete candidates by rehash-risk, and names the decisive checks.*

---

## Refined map — what is already done (corrects the rosy "all open")

- **GSPT for SDEs exists:** Berglund–Gentz, *Geometric singular perturbation theory
  for SDEs*, J. Diff. Eq. 2003 ([arXiv:math/0204008](https://arxiv.org/abs/math/0204008))
  — sample-path concentration near slow manifolds and through folds, via covariance
  tubes (**direct estimates, not blow-up**).
- **Stochastic variance up to codim 2 exists:** Kuehn's critical-transitions
  framework (normal-form variance scaling near fold and codim-2 fast bifurcations) —
  again **not** blow-up.
- **Deterministic blow-up exists** for: fold (Krupa–Szmolyan); cusp
  ([arXiv:1506.08679](https://arxiv.org/abs/1506.08679)); hyperbolic umbilic
  ([arXiv:2202.01662](https://arxiv.org/abs/2202.01662), 2024); folded node; folded
  limit cycle ([Jelbart–Kuehn–Kuntz 2024](https://arxiv.org/abs/2208.01361)).
- **Apparent gap:** the blow-up **method** applied to the **SDE** — desingularize the
  stochastic system, carry the noise through the rescaling, analyse the effective
  stochastic dynamics on the blown-up sphere.

## The honest tension

- **Small + safe** (redo the stochastic fold via blow-up) ⇒ **high rehash-risk**
  (Berglund–Gentz own the fold).
- **Clearly non-rehash** (noise on a singularity nobody has noised) ⇒ **not small**,
  and **pipeline-risk** (Kuehn's group is the obvious next mover).
- There is no "small AND clearly-novel AND safe" item lying in the open. Pick which
  axis to relax.

## Ranked candidates

### 1. Stochastic blow-up of the folded limit cycle — *recommended real target*
Take Jelbart–Kuehn–Kuntz's (2024) deterministic blow-up of the folded limit cycle in
three-timescale systems; add noise; prove a sample-path / exit estimate through the
singularity **in blow-up coordinates**.
- **Rehash-risk: LOW** — no stochastic version exists (search-verified). Structurally
  novel (nobody has noised this 2024 singularity).
- **Toolkit fit: perfect** — blow-up + Freidlin–Wentzell/sample-paths + oscillations:
  every thread you already have.
- **Size:** a **bounded first result is feasible** — leading-order stochastic passage
  through *one chart* of the blown-up singularity (sample-path concentration / exit
  law), not the whole theory.
- **Pipeline-risk: HIGH** (Kuehn) → the one Popović question.

### 2. Methodological kernel: effective noise on the blown-up fold sphere + exit law — *smallest*
Carry `σ dW` through the fold blow-up `v = ε^{1/3}V, w = ε^{2/3}W, t = ε^{−1/3}T`;
derive the **effective rescaled SDE** on the blow-up chart; compute the **exit/jump-
point distribution** through the fold from the blown-up dynamics.
- **Real iff** Berglund–Gentz's direct method gives the *spreading/variance* but not
  this exit-law-via-blow-up object. **Rehash-risk: MODERATE** — depends on exactly
  what math/0204008 already extracts at the fold.
- Smallest possible kernel; the natural first lemma that sets up the whole method.

### 3. Stochastic cusp via blow-up
- **Rehash-risk: MODERATE** — Kuehn's variance-up-to-codim-2 likely covers cusp
  *variance*; a full blow-up *sample-path* treatment of the stochastic cusp may still
  be new, but the boundary is thin.

### 4. Stochastic hyperbolic umbilic via blow-up
- **Non-rehash** (codim-3, deterministic only as of 2024) but **technically heavy**
  and pipeline-risk; not a small first project.

## Decisive checks before committing (do these, in order)

1. **Read Berglund–Gentz 2003 ([math/0204008](https://arxiv.org/abs/math/0204008))
   fold/bifurcation sections** — this *defines the novelty boundary*: exactly what the
   direct stochastic method already yields at a fold (and therefore whether the
   exit-law kernel #2 and the "what blow-up buys" claim in #1 are new).
2. **Popović question:** *"Is a stochastic / noise extension of the folded-limit-cycle
   blow-up (Jelbart–Kuehn–Kuntz 2024) in Christian's pipeline? If not, I'd like to do
   the leading-order stochastic passage through the blown-up singularity — it's my
   blow-up + FW toolkit."*

## Bottom line

The defensible "real, even if small, not a rehash" target is **#1 — the leading-order
stochastic passage through the blown-up folded limit cycle.** It is structurally novel
(un-noised 2024 singularity), toolkit-perfect, and admits a bounded first result. Its
only real risk is Kuehn's pipeline — a single Popović question. **#2** is the smallest
kernel and the natural first lemma, contingent on the math/0204008 boundary check,
which is the one piece of reading that most sharpens what counts as new here.
