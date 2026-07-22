# Novel angles on the 𝒲 identity
### What to actually attack, after the determinantal route was falsified

**The pivot (the whole point of this note).** Every attempt to identify 𝒲 has tried to fit it into a *spectral / determinantal* (random-matrix, Tracy–Widom, Painlevé-τ) frame, and every one has been falsified — most recently and cleanly by the θ-dependence of the gap determinant. The reason is now proven, not guessed: **𝒲 is a first-passage / escape / decay law of an unbounded-below operator, not a spectral gap.** So the right move is to stop asking "which determinant is it" and ask "what is the right *first-passage* framework." That reframing immediately opens three genuinely different homes, none of them RMT — and they decompose 𝒲 cleanly:

> **𝒲 is a Gamow survival law whose right tail is a Kramers large deviation and whose left tail is a persistence exponent.**

That one sentence is the conceptual target. The three pieces are three different, well-developed theories — and crucially, *not* the determinantal theory we've been failing in.

---

## A. The bulk / full law — Gamow resonance expansion (rigged Hilbert space)

The corrected map (survival expansion $\sum_n c_n e^{-i\lambda_n\Theta}$) is right, but it should be done as **real resonance theory**, not a numerical fit. The inverted-oscillator/cusp operator is *the* textbook open-quantum-decay system; its resonances $\lambda_n$ are **Gamow states** living in a rigged Hilbert space (Gel'fand triple), and the escape law is the decay of a metastable Gamow state.

- **Attack:** derive 𝒲 as the first-passage of the survival amplitude built from the resonance expansion — dominant resonance $\lambda_0=0.86-0.82i$ gives the leading exponential ($e^{-2|\mathrm{Im}\,\lambda_0|\Theta}$, rate $1.64$); the *interference of the resonance string* gives the bulk shape and the +skew.
- **Cheap check:** does the truncated resonance/residue expansion (resonances already computed and θ-validated) reproduce the FP-𝒲 CDF and cumulants? If yes, 𝒲's "closed form" *is* the Gamow expansion — a real named object, not a determinant.
- **Why novel here:** it drops the TW analogy entirely and adopts the decay-law / Gamow-vector framework, which is the correct universality analogy for escape over a quadratic maximum.

## B. The left (late-escape) tail — a *persistence exponent* (the standout idea)

This is the most novel handle and the one I'd push hardest. The escape is the first explosion of the Riccati $p=-u'/u$ — i.e. the **first zero of the stochastic Weber field $u$**. On the oscillatory side, "no escape yet" = "$u$ has not changed sign" = the process *persists*. So the anomalous "oscillatory phase-persistence" left tail the session found **is literally a persistence probability**, and its exponent is a **persistence exponent** in the precise sense of non-equilibrium statistical physics (Bray–Majumdar–Schehr).

- **Attack:** identify the Gaussian process whose first sign-change is the late escape (the Cole–Hopf / stochastic-Weber field), and compute/identify its persistence exponent — the decay rate of $P(\text{no zero up to }\Theta)$. That exponent *is* the left-tail law, and persistence theory has machinery (Rice formula, the independent-interval approximation, exact results for specific Gaussian processes) that the determinantal route never could.
- **Cheap check:** measure the empirical persistence probability of $u$ directly and fit the exponent; compare to the FP-𝒲 left tail. If they coincide, the "anomalous, no-closed-form" left tail becomes a *named persistence exponent*.
- **Why this matters:** it explains *why* the left tail resisted every closed form — it was never a tail of an integrable distribution, it's a persistence exponent, which generically is *not* elementary but *is* a well-defined, computable universality number.

## C. The right (early-escape) tail — Kramers / Freidlin–Wentzell LDP (already in hand)

The confining-side early escape is barrier crossing: the proved $e^{-h^{\star2}/6}$ no-early-escape bound *is* its large-deviation rate. This piece is done; the value is recognizing it's a *different mechanism* from B (so a single closed form was never going to cover both tails — which retroactively explains the retracted "clean collapse").

---

## D. Structural bet — PT-symmetric / pseudo-Hermitian reduction

The two conjugate-type ends (real $U$ / imaginary-argument $W$) and the complex-conjugate-paired resonances are the signature of a **PT-symmetric / pseudo-Hermitian** operator in its *PT-broken* phase. If a metric operator $\eta$ (the pseudo-Hermitian similarity) exists mapping the cusp operator to a Hermitian one, 𝒲 could be a *real-spectrum* law in disguise.

- **Attack/check:** test whether the resonances come in exact complex-conjugate pairs (PT-broken diagnostic), and whether an unbroken sector / a metric operator exists. The inverted oscillator is a known pseudo-Hermitian example, so there's literature to lean on.
- **Outcome either way:** a metric → 𝒲 is a known law transformed; no metric → a *proof* that 𝒲 is intrinsically non-Hermitian (PT-broken), sharpening the "genuinely new" claim.

## E. Define 𝒲 intrinsically — the isomonodromy fixed-point as its *defining equation*

We proved 𝒲 is a connection coefficient at a **fixed point** of the catastrophe-unfolding flow (no flow → no Painlevé τ). Flip that from an obstruction into a definition: the fixed-point condition is itself an **exact functional equation** for 𝒲's transform. TW is "defined" by Painlevé II; 𝒲 could be *defined* by its fixed-point/renewal equation, as a new named law — no need to reduce it to anything older.

- **Attack:** write the unfolding-flow fixed-point condition (and/or the first-passage renewal identity) as an explicit closed functional equation for $\hat{\mathcal W}$, and verify it's self-contained.

---

## Recommended order (and what to dispatch now)

1. **B (persistence) + A (resonance expansion) in parallel** — these are the two that can actually *produce the law*: A for the bulk, B for the left tail. Both have cheap, decisive numerical checks against the FP-𝒲 we already trust, and both reframe 𝒲 into the correct (first-passage) universality language. **This is the map to do now.**
2. **C** is recognition, not work.
3. **D** (PT metric) and **E** (intrinsic equation) are the deeper structural bets — pursue whichever B/A leaves most natural.

**Honest framing.** A–E are reasoning by analogy to the right neighboring theories (open-quantum decay, persistence, PT-symmetry, isomonodromy fixed points), not theorems. The point is that they're the *correct* neighbors — first-passage, not spectral — and three of them carry cheap FP-validated checks. Treat as a map of where to dig; tag everything on execution; the persistence-exponent identification of the left tail is the highest-novelty, highest-value target.
