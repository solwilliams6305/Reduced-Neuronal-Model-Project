# Design note — the inner chain
inner form → blow-up → inner Riccati → Airy operator → shooting characterisation

Prompted by: "I tried visualising this before and it was tough to wrap my head around."
Draft 2026-07-27. This is modules α3–α4 of `ARCHITECTURE.md`, not the pilot.

---

## 1. Why it resists visualisation — a specific diagnosis

It is **not** that the ideas are individually hard. It is that the chain is a sequence of
**coordinate changes**, not a sequence of objects. At every step the axes change meaning:

| step | what the axes mean |
|---|---|
| inner form | fast/slow phase-plane variables, rescaled |
| blow-up | chart coordinates on a blown-up neighbourhood |
| inner Riccati | a scalar R against swept time Y |
| Airy operator | a wavefunction u against Y, plus a spectrum |
| shooting | a parameter sweep against the location of a zero |

Nothing persists across the five. So working memory has to hold five different
meaning-assignments simultaneously and the mapping between them — maximal element
interactivity, which is exactly the regime where the transient information effect bites hardest.
Animating it *first* would make this worse, not better.

Second, compounding: in the standard presentations each step is a **trick**.
Blow-up is a trick, Cole–Hopf is a trick, shooting is a trick. Tricks do not compose in memory —
you can follow each and still not hold the chain, which is precisely the reported symptom.

## 2. The fix, part one: one invariant, five charts

**Choose a single noisy escaping trajectory and never change it.**

Every panel in the chain redraws *that same trajectory* in a new representation, with the
previous representation ghosted or held in an inset. Nothing new is introduced at any step;
only the chart changes. This converts *"five transformations"* into *"one object, five views"* —
a completely different cognitive load, and the standard fix for a sequence with no natural
referent.

Concretely: colour the trajectory once, keep that colour through all five panels, and keep a
persistent marker (say, the moment of escape) visible in every chart. The viewer's anchor is
"where did my point go", never "what is this new picture".

## 3. The fix, part two: reframe every trick as a **trade**

A trick is arbitrary. A trade has a motive, so it survives in memory.

### Blow-up = **anisotropic zoom**
Do not open on the sphere. Open on a failure (this is finding F5 in `PREREQ_GRAPH.md`):

1. Zoom in on the fold *naively* — magnify both axes equally. Everything flattens onto a line.
   You see nothing. **The zoom fails.**
2. Zoom in *anisotropically* — magnify each direction by a different power, chosen so the
   dynamics stay O(1). Structure appears.
3. **That is blow-up.** The charts K1/K2/K3 are which face of the anisotropic zoom you look
   through. The exceptional divisor is the set of directions you could have approached from —
   which the naive zoom collapses to a point and the anisotropic zoom keeps apart.

Interactive form: a magnifier with **two independent zoom sliders**. Move them together → the
picture degenerates. Move them at the correct weight ratio → the fold opens up. The viewer
discovers the exponents by feel before ever seeing them written.

This also repairs the ledger's blow-up distortion: nothing is "added" to the space, you are
looking at the same neighbourhood through a distorting lens.

### Cole–Hopf = **trading blow-ups for zeros** ← the crux
R = −u′/u means **R blows up exactly where u crosses zero.**

Twin panel, linked: u(Y) smooth and oscillating, crossing zero; directly beneath, R(Y) with a
vertical asymptote at each of those crossings. Sweep Y and watch the correspondence.

The reframe, said out loud: *blow-ups are hard to see and hard to compute with; zeros are easy.
Cole–Hopf buys you zeros.* Now it is not a trick — it is an obvious trade with an obvious motive,
and it makes "escape ⇔ first node" self-evident rather than asserted.

**This is the single highest-value artifact in the chain.** It is exact (no schematic needed),
cheap to build, directly interactive, and it demystifies the step that is normally most opaque.
Build it first, before anything else in α3–α4.

### Shooting = **where does the first zero land**
Once u exists, shooting is visually trivial: launch u with prescribed behaviour at one end, sweep
the parameter, watch the first zero slide along. The escape level is where the first zero reaches
the boundary. One slider, one moving dot. The apparent sophistication evaporates.

### The operator = **the zeros are the spectrum**
Oscillation theory: the n-th eigenfunction has exactly n nodes. So the escape statistics and the
spectrum are the *same picture counted two ways* — not two theories to reconcile. With the twin
panel already established, this is a relabelling, not a new idea.

## 4. What NOT to do

- **No grand five-stage diagram.** The instinct with a chain is one big figure showing all of it.
  That maximises element interactivity — the exact failure mode diagnosed in §1.
- **Do not open on the blow-up sphere.** It is the bookkeeping, not the idea, and it is the
  standard reason this topic loses people.
- **Do not use a different example trajectory per stage.** Kills the invariant, which is the
  whole mechanism of the fix.
- **Do not animate before the stills pass.** Especially here: this is the highest-element-
  interactivity material in the whole course.

## 5. Dimension honesty → new ledger entry
The real object lives in more dimensions than can be drawn, and part of the reported difficulty
is unacknowledged slicing. Every panel in this chain must carry a small inset showing **where the
slice sits** in the full space. Without it the viewer silently believes the 2D picture is the
object — the same class of error as showing a swallowtail surface as if it were the unfolding.

## 6. Proposed beat order
α3a naive zoom fails · α3b anisotropic zoom = blow-up, charts as faces · α3c what the chart
contains: the Riccati, R → ∞ · α4a the trade: blow-ups for zeros (twin panel) · α4b shooting:
sweep, watch the first zero · α4c the zeros are the spectrum.

Six beats, one trajectory, no new objects after beat one.
