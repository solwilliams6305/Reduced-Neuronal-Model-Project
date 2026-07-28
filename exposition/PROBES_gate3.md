# Acceptance probes — `borel` (γ2+γ3) and `fastslow` (α1+α2)

**Status: instrument only. NOT ANSWERED.** These are the two modules RETROSPECTIVE.md §8 names
first — one from each trunk — on the reasoning that if the explanations land, that is real evidence
for the other ten, and if they do not, better to learn it on two than on twelve.

**Who answers these: you, not me.** PROBE_PROTOCOL.md assigns the reviewer role to the user at
every gate ("user answers all five probes", "user drives it", distortion "confirmed/corrected by
the domain expert"). The author cannot run them on himself — I know the intended answer for every
one, so my passing them measures nothing. RETROSPECTIVE.md §3 puts it directly: *self-certification
catches arithmetic, it does not catch miscommunication, and the single defect a human found was in
the class no self-check can reach.*

---

## A protocol violation to declare up front

PROBE_PROTOCOL.md §14 requires the five probes be **written out BEFORE the module is built, "so
they can't be retrofitted to whatever the artifact happens to teach."** That did not happen. These
are being written after the fact, which is precisely the failure mode the rule guards against.

**Mitigation, and please judge whether it is enough:** every probe below is derived from the
module's entry in `captions.json` — the *stated intent*, fixed as the spec — and not from what the
artifact turned out to display. Where the artifact might satisfy a probe by accident, that is noted
inline. This weakens the retrofit problem; it does not remove it. If a probe looks suspiciously
well-matched to the widget, that is the failure mode showing, and it should be discounted.

---

# Module `borel` — γ2+γ3, the linchpin

Spec: *"The ambiguity IS the exponentially small term."* This is the node the whole analysis thread
is built on, so it is the highest-stakes module in the sequence.

### Probe 1 — Predict-then-commit (tests the SETUP, not the reveal)
Read only the `setup` and `cut` captions and look at the Borel plane with the pole on the contour.
**Before touching the slider, commit to an answer:** you are about to be shown the same integral
evaluated with the detour passing above the pole and then below it. Write down what you expect to
differ between the two answers — nothing, the real part, the imaginary part, or both — and roughly
how big the difference should be.

- **Fails if:** you cannot form a prediction at all. That means the setup never established that
  there was a choice to be made, and the payoff cannot land however clean it looks.

### Probe 2 — Transfer to an unshown case
The module shows a Borel transform with a **single pole on the positive axis**. The programme's
actual series has a **branch point** there, plus a **conjugate pair off-axis** that the module never
displays. Without reopening the artifact: does the off-axis pair produce an ambiguity for real
positive x, and why or why not?

- **Fails if:** you can only restate the on-axis story.
- *Note:* this is the exact point the expert pass flagged as the pilot's main defect, so it is the
  probe most likely to expose whether the fix actually taught the distinction or merely deleted the
  confusing picture.

### Probe 3 — Cold explain-back (artifact closed)
Reconstruct, on paper, why *above minus below* equals 2πi·e^(−A/x) — where the 2πi comes from,
where the exponential comes from, and why the location of the obstruction ends up controlling the
*size* of the ambiguity.

- **Fails if:** you need to reopen the widget to rebuild the chain.

### Probe 4 — Falsification (expected to catch the most)
**What would the widget look like if the claim were false** — if the two contours differed by
something that was *not* 2πi·e^(−A/x)? Name a specific visible signature you would see. Then ask
the sharper version: does the displayed picture actually distinguish the true claim from a nearby
false one, or would it look the same either way?

- **Fails if:** the visual is compatible with every outcome — then it is decoration, not evidence.
- *The module's own answer, for comparison after you have committed:* the error readout is measured
  against 2πi·e^(−A/x), showing 2.0e-3 / 1.1e-4 / 2.4e-5 at x = 0.25 / 0.40 / 0.60, and the Cauchy
  check shows the value is independent of detour height. **Ask whether those numbers are visible as
  evidence, or whether they read as a self-report the viewer is asked to trust.**

### Probe 5 — Distortion check (MANDATORY — needs your ruling)
Three distortions are named on-screen: `distortion_model` (a pole, not the branch point the real
problem has), `distortion_pole` (residue formula becomes a discontinuity across a cut), and
`distortion_quadrature` (widget refuses below a minimum detour height rather than silently
returning a wrong number).

- **Your ruling needed:** are these correct, correctly weighted, and complete? Specifically — is
  there a *fourth* distortion an expert would name that the author missed?

---

# Module `fastslow` — α1+α2, the geometry and the number that kills it

Spec: *"Fast, slow, and the place it breaks."* α1 and α2 were merged deliberately so the geometry
and its failure arrive together, and so `blowup` (α3) opens on that failure.

### Probe 1 — Predict-then-commit
Read only `setup` and `manifold`. **Before revealing the attraction rate:** the trajectory crawls
along the cubic curve and then jumps. Commit to an answer — *where* on the curve will the jump
happen, and what quantity do you expect to be doing something special there?

- **Fails if:** no prediction is possible; the geometry was shown but its consequence was not set up.

### Probe 2 — Transfer to an unshown case
The module shows the planar FitzHugh–Nagumo cubic, with attraction rate 1 − v² vanishing at the two
folds. Without reopening it: for a **higher-order** turning point — the cusp, where the potential
goes like Y² rather than Y — would you expect the attraction to vanish *faster* or *slower* at the
degenerate point, and what would that do to the size of the region where Fenichel fails?

- **Fails if:** you can only restate the planar case.
- *Note:* this probe deliberately reaches toward δ2 (the ladder). If it lands, the two trunks are
  genuinely connected for a learner and not just in the author's head.

### Probe 3 — Cold explain-back (artifact closed)
Reconstruct the chain: why does slow–fast splitting need attraction bounded away from zero, what
exactly does Fenichel's theorem provide, and what breaks when 1 − v² → 0? State what is *lost* —
not merely that "the theorem does not apply."

- **Fails if:** the answer is "the theorem stops working" without saying what that costs.

### Probe 4 — Falsification
**What would the picture look like if the attraction did NOT vanish at the fold?** Describe the
trajectory you would see instead. Then: does the module's animation actually show you the rate
going to zero, or does it show a jump that you are *told* is caused by a vanishing rate?

- **Fails if:** the visual would look identical under both stories — the jump is dramatic and might
  carry the viewer's attention regardless of whether the mechanism is displayed.
- *This is the one I would most expect to fail*, because the jump is visually louder than the rate.

### Probe 5 — Distortion check (MANDATORY — needs your ruling)
Two are named on-screen: `distortion_planar` (planar FHN only; folded nodes/saddles in higher
dimensions are not shown) and `distortion_rate` (1 − v² is an eigenvalue, i.e. a rate; "strength of
the pull" is a gloss).

- **Your ruling needed:** correct, correctly weighted, complete? In particular — is calling a rate
  "the strength of the pull" a fair gloss or a real distortion of what a learner will carry forward?

---

## Recording the outcome

Per PROBE_PROTOCOL.md a module that looks good and fails its probes is a **FAILED module, and the
polish is evidence against it, not for it.** So the useful outcome here is a failure, not a pass.

Suggested: answer inline, then a one-line verdict per module — *passes / fails on probe N / needs
rebuild*. Any distortion you add or correct in probe 5 goes into DISTORTION_LEDGER.md and, if it is
on-screen text, into `captions.json` (never edit the generated HTML — it is rebuilt from the
template and would be overwritten).
