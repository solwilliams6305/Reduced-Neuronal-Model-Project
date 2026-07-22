# Novel angles on the Riemann–Hilbert direction
### Closing the two open coupled-FHN cusp-escape theorems

Two open theorems are flagged "Riemann–Hilbert direction":

- **(T1)** a *uniform Berglund–Gentz tube* through the fold merge — the stochastic Airy↔Weber connection that closes the noisy fold→cusp crossover.
- **(T2)** the *Painlevé-IV τ-function* for the cusp escape law — the edge-law identity (the soft-edge Fredholm-determinant form is already refuted: no projection at a quadratic turning point).

The central claim of this note: **they are two readings of one object**, and attacking them separately is the wrong move.

---

## 0. Keystone reframe — one coalescing-turning-point RH problem, read two ways

Olver's deterministic connection is already the right skeleton: an *isolated* turning point gives the Airy parametrix; *two coalescing* turning points give the parabolic-cylinder (Weber) parametrix. Here Δ(g)=2√(−2g/3) **is** the coalescence parameter.

Translate to RH: "two turning points sliding together" = two branch points of the g-function / equilibrium measure colliding, and the **local parametrix changes type Airy → parabolic-cylinder exactly at Δ∼ℓ**. So build one Δ-dependent RH problem X_Δ whose Δ→∞ limit is the Airy/TW problem and whose Δ→0 limit is the cusp edge, with a parabolic-cylinder/PIV local parametrix at the colliding edge.

**Punchline.** The *same* uniform local parametrix delivers both theorems:
- its **probabilistic** reading (concentration of the steepest-descent path under noise) is the uniform tube **(T1)**;
- its **integrable** reading (the τ-function of the parametrix) is the PIV determinant replacement **(T2)**.

So the single keystone is a **uniform parabolic-cylinder parametrix carrying noise**. Everything below is either a route to it or a cheap test that tells you which route is real.

---

## 1. (T1) Replace the OU tube by a Weber comparison process — the tractable first theorem

**Why Berglund–Gentz degenerates.** The tube variance solves v′ = 2a(t)v + σ², a(t) = curvature. At a simple fold a(t) stays bounded (Gaussian/OU tube is fine); at the cusp the two folds collide, **a(t)→0**, so v blows up and the tube widens non-uniformly in Δ/ℓ. The OU comparison process is simply the wrong Gaussian near the cusp.

**Fix.** Near the merge the linearized fluctuation operator *is* the parabolic-cylinder (Weber) operator. The correct comparison is the **Gauss–Markov process with the Weber Green's function**, whose tube width spreads *algebraically* (Δ-controlled) and stays finite through the merge. Then splice: OU tube for |Y|≫Δ (Airy region) ∪ Weber tube for |Y|≲Δ, matched in the overlap — the stochastic analogue of Olver's uniform matching, done at the level of the **covariance**, not the trajectory.

**Cheap first test.** Simulate the noisy inner equation; measure the empirical tube width as a function of (Y, Δ, η). Check it tracks the Weber Green's-function variance (finite at Δ=0) and *not* the OU variance (which diverges). Falsifiable, derisks the whole estimate before any analysis.

`[CONJECTURE + PROGRAM]` — most tractable of the two; this is probably the right "next theorem."

---

## 2. (T2, route A) Get PIV as a controlled degeneration of the Pearcey RH problem

The session already found the slow flow collapses **Pearcey (3rd order) → Weber (2nd order)**. Pearcey *already has* a known integrable/RH structure (3-saddle / higher-genus). So don't build PIV from scratch — **obtain it as the slow-flow reduction of the Pearcey RH problem** (project out the third saddle / the y-direction). The cusp edge law would be the slow-flow marginal of the Pearcey process.

**Cheap test.** Take the Pearcey kernel, apply the slow-flow projection (integrate out the third direction with the measured weighting), and check the resulting law reproduces the cusp fingerprint: cubic QQ Q+0.072Q²−0.035Q³, the e^(−|s|⁵/20) tail, κ₅≈−2.2, κ₆≈−2.7. A match means the RH problem is "Pearcey, reduced," and PIV falls out of a known parent.

`[SPECULATIVE]` — leverages existing Pearcey machinery; test is cheap.

---

## 3. (T2, route B) β = isomonodromic time; the recurring cubic = the PIV Hamiltonian flow

The empirics are pointed: the β-family is a one-parameter **cubic-transform** family (subordination), and the *same cubic* appears in the q-direction (TW→cusp) and the β-direction. A recurring cubic deformation is the fingerprint of a **Hamiltonian (PIV) flow with a cubic vector field**.

**Conjecture.** β (noise) is, up to reparametrization, the PIV **isomonodromic time** (or one of the PIV parameters), and the cubic QQ-map is the linearization of the PIV flow about the bare (β→∞) skeleton law. This also *explains* the empirical "subordination / random time change": the isomonodromic time is the time in the time-change.

**Cheap, sharp test.** Fit the cubic coefficient(s) φ(β) across β and check whether φ(β) satisfies the ODE the PIV Hamiltonian flow predicts (a specific Riccati/cubic ODE in β), versus an arbitrary curve. A positive fit is strong evidence for the PIV class *and* tells you which PIV parameter β is.

`[CONJECTURE]` — sharply falsifiable, and arguably the single most informative cheap test here.

---

## 4. (T2, route C) The sub-Gaussian flip = analytic continuation across a Stokes ray

The left-tail exponent 2q+1 (3 for TW, **5** for the cusp) matches the turning order q — the parabolic-cylinder hierarchy. But the cusp is **sub-Gaussian** (lighter tail, negative higher cumulants) where the multicritical / PII hierarchy is *heavier*. Opposite sign of deformation.

**Reading.** The cusp law may be the **analytic continuation of a multicritical edge solution across a Stokes ray** — the "wrong-sign" branch (cf. the Ablowitz–Segur vs Hastings–McLeod split of PII), which flips heavy → sub-Gaussian. RH solutions naturally come in such Stokes-separated families.

**Cheap test.** Compute the known multicritical (q=2) edge law and continue it in the Stokes parameter; check whether the continued solution matches the cusp's sub-Gaussian cumulants. A match identifies the cusp as a *known* object sitting on a different Stokes sector — instantly supplying the RH problem.

`[SPECULATIVE]`.

---

## 5. (cross-cutting) A genuine *stochastic* Riemann–Hilbert problem

The boldest unifier: put η directly into the **jump matrix** as a random perturbation. Then the uniform tube (T1) becomes *stability of the RH solution under a random jump, uniform in Δ*, and the τ-function (T2) becomes its random determinant whose law is the edge distribution. "Carry noise through the blow-up" becomes literally "carry noise through the RH steepest descent," and the hard content collapses back to the §0 uniform parabolic-cylinder parametrix.

`[SPECULATIVE / foundational]` — most novel, least standard; a longer-horizon framing, not a near-term theorem.

---

## Recommended order (cheap falsifiers first)

1. **§1 tube-width test** (Weber vs OU variance) — derisks T1, the tractable theorem.
2. **§3 β-flow cubic-coefficient test** — cheap, sharply falsifiable; tells you if isomonodromy/PIV is real and which parameter β is.
3. **§2 Pearcey-projection test** — if it matches, you inherit a parent RH problem for free.
4. **§4 multicritical-continuation check** — if it matches, the cusp is a known object on a different Stokes sector.
5. Commit to whatever survived and build the **§0 uniform parabolic-cylinder-with-noise parametrix**, which closes T1 and T2 together.

**Honest caveat.** §0–§5 are reasoning by analogy to known RH / Painlevé structures; none is checked against the literature here, and the prefactors/constants are not pinned. The point is the falsifiers — they're cheap, and at least three can be run with machinery already in the repo. Treat this as a map of *where to dig*, not results.
