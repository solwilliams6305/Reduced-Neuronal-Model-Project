# Gate 6 — retrospective
2026-07-28. Task 6 of the agreed gate structure: the explicit anti-tunnel-vision checkpoint,
where "this tooling is wrong, change it" is a legitimate outcome.

---

## 1. The headline, and it is not good news

**Twelve modules were built. Zero probes were run.**

`PROBE_PROTOCOL.md` was agreed at gate 2 and states that acceptance is *generative*: the reviewer
predicts, transfers, explains back cold, tries to falsify, and rules on the distortion. Gate 3
(stills) and gate 4 (interactive) both require that. Neither has been satisfied for **any** module.

The gates existed precisely to stop me building past the point of validation. They did not fire.
Every module was self-certified: I wrote the explanation, checked its arithmetic, wrote down where
I thought it lied, and moved on. That is exactly the loop the protocol was designed to interrupt.

This is the most important finding in this document, and everything below should be read against
it. **The sequence is unvalidated as pedagogy.** It is well-validated as arithmetic.

## 2. What exists

24 commits. 12 interactive modules, each a Python builder + an HTML template + a generated
self-contained page.

| module | concept | html |
|---|---|---|
| hook | ep0 — a neuron fires at random | 14 KB |
| fastslow | α1+α2 — geometry, and where it breaks | 14 KB |
| blowup | α3 — anisotropic zoom | 14 KB |
| pilot_borel | α4a — R↔u twin panel + scrubber | 17 KB |
| shooting | α4c — the zeros are the spectrum | 29 KB |
| tracywidom | α5 — the same law from elsewhere | 17 KB |
| divergence | γ1 — a series that gets worse | 15 KB |
| borel | γ2+γ3 — the ambiguity IS the term | 15 KB |
| stokes | γ8 — the jump that isn't a jump | 18 KB |
| junction | δ1 — one ladder, two halves | 18 KB |
| ladder | δ2 — structure survives, solvability doesn't | 93 KB |
| coda | δ3 — three families, one word | 11 KB |

Supporting: `captions.json` (the narrative contract), `DISTORTION_LEDGER.md`,
`PREREQ_GRAPH.md`, `ARCHITECTURE.md`, `PROBE_PROTOCOL.md`, `DESIGN_inner_chain.md`,
`SURVEY_BRIEF.md`.

## 3. What actually caught defects

Roughly 41 defect mentions across the commit log. Sorted by *how they were found*, which is the
part worth knowing:

**Comparing a number against an independent source — caught every substantive error.**
- `divergence`: interpolated E(x) carried ~5e−8 error, larger than the optimal-truncation floor
  the module exists to display. The browser reported the wrong best-N.
- `hook`: spike detector counted threshold jitter. **41% of "spikes" were not spikes**, and the
  shipped reference histogram had a spurious peak at gap≈0.
- `tracywidom`: GUE normalisation wrong; rescaled mean *diverged* with N instead of converging.
- `pilot_borel`: `ai_zeros` double-negated, putting every marker outside the plot window.

**An assertion firing.** `hook`'s regime check (which then needed honest recalibration, not silent
relaxation); `ladder`'s WKB underflow at q=4.

**Writing the distortion statement.** Two defects were found *by the act of writing down where the
picture lies*, not by looking at it: row 3 of the twin panel broke the one-invariant discipline
while reusing the invariant's colour; and the boundary-condition caption was fixed in a still and
shipped missing from its interactive.

**Rendering and looking.** Caught presentation problems in quantity — label collisions, an
unreadable semicircle, a legend occluding a lobe. Caught exactly **one** substantive error: the
spurious least-term in the storyboard, and only because the node at k=3 made the failure visible.

**The reviewer.** One defect, and it is the interesting one: the conjugate pair on the Borel panel
was the most salient object and was *not* the source of the ambiguity. Numerically flawless.
Visually clean. Pedagogically misdirecting.

### The conclusion that follows
Self-certification catches arithmetic. It does not catch miscommunication. The single defect a
human found was in the class that no self-check can reach — which is direct evidence that the
probes are not ceremony, and their absence is a real hole rather than a bookkeeping one.

## 4. Coverage gaps this retrospective found

Auditing rather than remembering turned up two things I would not have reported otherwise:

- **`blowup` has no self-test and no assertion.** It is the only module with *zero* verification
  of any kind. It is also the module the prior-art sweep called *"arguably the single most
  animation-hungry idea on this list and completely unserved"* — so the highest-opportunity module
  is the least-checked one.
- **`coda` has no self-test**; `shooting` and `build_scrubber` have self-tests but no physics
  assertion in the builder.
- **59 of 68 ledger entries remain open.** Nine are resolved. The open ones are not all defects —
  many are honest standing caveats — but several explicitly say *needs expert ruling* and have
  never had one.

## 5. Did the architecture hold?

**Yes, and it is the part I would keep unchanged.** Split by computational cost, not language:
expensive work stays in Python and ships as a versioned data artifact; interactive-rate work is a
small JS kernel; golden vectors bind them. Every module runs a port check on load and reports it
on the page.

Two refinements were forced by failures, and both are now standing practice:

1. **Verification split by role.** Golden vectors are computed at *the resolution the browser
   uses*, so they test the **port**; the physics is asserted separately in Python at ~10×
   resolution. Learned from `divergence`, where golden vectors passed at 1e−16 while the module
   was wrong — a port test that shares the browser's discretisation cannot catch a discretisation
   error.
2. **A narrative contract to match the kernel contract.** `captions.json` exists because a caption
   fixed in a still shipped missing from its interactive. The numbers were in lockstep; the words
   were not. Builders now fail if a required caption key is absent.

**What the architecture did not need:** Manim, Pyodide, WASM, or a shared runtime of any kind.
The Manim retirement looks correct in hindsight — nothing built here would have been easier as a
rendered scene, and several things (dragging a contour around a pole) are impossible as one.

## 6. Does any of this warrant video?

**Mostly no, and I would drop task 5 as originally scoped.**

The argument that produced the architecture also argues against video here. Transience is the
problem; learner control is the fix; a linear render reintroduces exactly what the scrubbers
remove. And the core interaction of nearly every module is a control the viewer holds:

- `borel` — drag the contour around the pole and watch the imaginary part flip
- `stokes` — push the radius out and watch a 45° ramp sharpen into a step
- `tracywidom` — toggle raw/rescaled and watch four curves collapse onto one
- `fastslow` — slide ε until a curve becomes a circuit
- `divergence` — add terms and watch the error turn around

None of those survive being filmed. Filming them converts a thing you *do* into a thing you
*watch*, which is the exact trade the transient-information literature says to avoid for
high-element-interactivity material.

**The one case for video** is the connective tissue: *why these twelve things are one argument*.
That is genuinely linear, it is the thing the explorables cannot carry, and nobody clicks through
twelve pages in order without a reason. A single short piece — the shape of the argument, linking
out to the modules — is defensible. Twelve narrated modules are not.

## 7. What I would do differently

- **Run the probes after module one.** Not after twelve. The cost of being wrong about the
  explanation compounds across every module that inherits its conventions.
- **Invent `captions.json` at module one.** It arrived mid-sequence, so the early modules were
  retrofitted rather than built to it.
- **Write the distortion statement *before* building**, not after. Twice it found a defect; both
  times the module was already built and the fix was rework.
- **Audit rather than remember.** Both findings in §4 came from grepping, not recall. I would not
  have reported `blowup`'s missing verification, because I believed the pattern was universal.

## 8. Recommended next steps, in order

1. **Run the probes on `borel` and `fastslow`.** Two modules, one from each trunk. If the
   explanations land, that is real evidence for the other ten; if they do not, better to learn it
   on two than on twelve.
2. **Add verification to `blowup`.** No module should ship with neither an assertion nor a port
   test, least of all the one with the biggest gap to fill.
3. **Triage the 59 open ledger entries** into *standing caveat* versus *needs a ruling*. Only the
   second kind needs your time.
4. **Decide the video question** — my recommendation is one connective piece, no per-module
   narration.
5. Only then consider polish, hosting, or sequencing into a course proper.
