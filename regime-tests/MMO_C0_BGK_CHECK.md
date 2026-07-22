# BGK prior-art check for the `c₀` programme — what is already theirs

*Direct read of Berglund–Gentz–Kuehn, "From random Poincaré maps to stochastic
MMO patterns" (arXiv:1312.6353, JDDE 2015). **Coverage:** I read the intro/results
overview and §5.3 (Deterministic dynamics near the folded node) in full from the
arXiv PDF text; I did **not** read §6's theorem proofs, Appendix C, the companion
"Hunting French Ducks" (1011.3193), or Wechselberger [77] ([77]=their ref for the
folded-node canard theory). Verdict below is bounded by that coverage.*

---

## Verdict: most of the "proof" is already BGK's (or classical). The only possible
## survivor is the explicit constant — and that is unconfirmed.

The machinery I had been re-deriving over the last several passes is in BGK §5.3
essentially verbatim, **including the part I had called "the one genuinely hard
remaining theorem."** Specifics, with their equation numbers:

1. **The variational equation is theirs (their (5.26)):**
   `μ u₁' = 4z̄u₁ + 2u₂ − 2u₁²`, `μ u₂' = −2(1+μ)u₁` — identical (linear part
   `A(z)=[[4z,2],[−2(1+μ),0]]`) to what `MMO_C0_PINNED.md`/`MMO_C0_PROOF.md` used.
   They attribute the canard/Weber-equation analysis to Wechselberger [77].

2. **The conserved structure is classical (their (5.27)):** the first integral
   `K = [1 + (2/(1+μ))(u₂−u₁²)] e^{−2u²/(1+μ)}`, which they state is "just a version
   of the **classical** first integral near planar degenerate folded singularities"
   (their refs [56,54,69]). My "Liouville area identity" `det M = e^{−2z₀²/μ}` is an
   elementary cousin of this — not new.

3. **The averaging + O(μ) error control is theirs — this is my "ingredient B".**
   Their §5.3 Step (S2) is "use an **averaging-type transformation** to describe the
   rotations… around the weak canard," and **Proposition 5.5 (Averaged system)** gives
   it with explicit error `O((μ+z̄²)ρ(K))`, rotation rate `c−(1+|log K|^{1/2}) ≤
   μ dϕ/dz̄ ≤ c+(…)`. They also note the stochastic version was already in their
   earlier paper "[17, Section D.2]." So the non-stationary-phase / averaging O(μ)
   bound I presented as the new lever in the last pass **is BGK's averaging result.**

4. **The deep-canard / large-k regime is theirs — this is my "Airy turning-point
   residue".** Their **saturation effect** (intro pts 3–4; **Theorem 6.4**, "local map
   for outer sectors… in particular it proves the saturation effect") is exactly the
   large-SAO / deep-canard regime I flagged as the one hard theorem left. They prove
   estimates there. So that residue is **not** open either.

## What might survive (unconfirmed, and now lower-probability)

I did **not** find, in what I read, the **explicit closed-form spacing constant**
`c₀(φ₀) = π²sin²φ₀/[(1+μ)(2φ₀+sin2φ₀)²]` or its value `π²/16`. BGK parametrise by
the first integral `K` and keep the O(1) rotation constants **general** (`c±` in
Prop 5.5) rather than evaluating the leading spacing constant to a clean number —
plausibly because their interest is the sector/saturation/noise structure, not the
bare constant. So the explicit evaluation *could* still be a small, citable
gap-fill on top of their machinery.

**But this is not confirmed new, and the surrounding evidence pushes against it:**
- The "BGK leave `c₀∈[π/4,1]`, Theorem 4.4, §6" claim in `MMO_C0_PINNED.md` refers
  to the **companion paper 1011.3193** (its numbering matches `4.x`, not this paper's
  `5.x`). I have **not** read 1011.3193, so whether they *bound* vs *evaluate* `c₀`
  there is unverified by me. If they (or Wechselberger [77]) already give the
  explicit constant, nothing survives.
- Everything *around* the constant is theirs, so the prior that the constant (or an
  equivalent K-level/sector-width formula) is somewhere in {1011.3193, [77], App C}
  is now substantial.

## Honest consequence (correcting the last several passes)

- My "non-stationary-phase closure on the bulk" (`MMO_C0_NONSTATPHASE.md`) and the
  framing of the deep-canard turning point as "the one remaining hard theorem"
  (`MMO_C0_PROOF.md` §6, `MMO_C0_LEMMA1.md`) were, in effect, **re-deriving BGK's
  Prop 5.5 and Theorem 6.4.** That work is not new.
- The `c₀` programme reduces to **at most** the explicit constant `π²/16` / the depth
  function — a minor evaluation on top of BGK's framework — and even that is **not
  yet confirmed absent** from BGK 1011.3193, Wechselberger [77], or BGK §6/App C.
- This is consistent with the *original* audit (`NOISE_ERROR_BOUNDS.md`), which had
  already flagged `c₀` as "a small, citable sharpening" and the row as "BGK's." The
  last few passes drifted into treating the surrounding analysis as a contribution;
  it is not.

## The decisive remaining checks (before any novelty claim)

1. **Read 1011.3193 ("Hunting French Ducks") §4 + §6** — does it state `c₀` as a
   range (gap to fill) or evaluate it? This single check is decisive.
2. **Read Wechselberger [77] (2005 folded-node canard paper)** — are the secondary-
   canard sector widths / spacing already explicit there?
3. **Read BGK 1312.6353 §6 + Appendix C** — confirm the rotation constants are kept
   general (not evaluated to `π²/16`).

Until at least (1) is done, the honest status of the explicit constant is
**"possibly a minor gap-fill, more likely already known"** — not a result to build a
dissertation headline on.
