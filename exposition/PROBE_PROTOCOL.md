# Acceptance-probe protocol
**Status: AGREED by user 2026-07-27.** Governs every module. Not optional, not per-module discretion.

## Why probes and not ratings
Processing fluency is misread as comprehension. Production value raises fluency without
necessarily raising learning — so the better an artifact looks, the more effectively it fools both
author and reviewer. A satisfaction rating ("does this feel clear?") measures precisely the
confounded quantity. Every gate is therefore **generative**: the reviewer must produce something,
not judge something.

Corollary: a module that "looks great" and fails its probes is a FAILED module, and the polish is
evidence against it, not for it.

## The five probes

Every module ships with all five, written out BEFORE the module is built (so they can't be
retrofitted to whatever the artifact happens to teach).

### 1. Predict-then-commit
Reviewer states, before the reveal, what they expect to happen.
- **Fails if:** no prediction is possible. That means the setup never established enough structure
  for the payoff to land, however good the payoff looks.
- Test of: the *setup*, not the reveal.

### 2. Transfer to an unshown case
Reviewer answers a question about a case the module never displayed.
(e.g. after the fold module: sketch what the cusp does.)
- **Fails if:** the reviewer can only restate what was shown.
- Test of: structure vs. recall.

### 3. Cold explain-back
Reviewer reconstructs the argument with the artifact closed.
- **Fails if:** reconstruction requires reopening it.
- Test of: durable understanding. Fluency evaporates here; understanding does not.

### 4. Falsification probe
*What would this picture look like if the claim were false?*
- **Fails if:** the visual is compatible with every outcome — then it is decoration, not evidence.
- Test of: whether the visual actually carries the argument's content.
- Underused; per user + author agreement, expected to catch the most.

### 5. Distortion check — MANDATORY ON EVERY MODULE
*Where does this visual lie?*
Every visual metaphor has a finite domain of validity. The Borel-plane cartoon WILL misrepresent
something; the low-dimensional picture of blow-up WILL misrepresent something.
- **Requirement:** the distortion is named explicitly by the author, confirmed/corrected by the
  domain expert (user), and **stated on-screen in the artifact itself** — not buried in notes, not
  discovered later by someone who trusted the picture.
- Test of: guards the "visual intuition misleads in analysis" failure mode, which is the specific
  risk for material where the geometric picture and the analytic truth diverge.
- No module ships without a completed distortion statement.

## Build-order rules that follow from the transient information effect

- **Static-first.** Storyboard as stills; test the stills; animate only after they pass. If stills
  don't carry it, animation smooths over the failure rather than fixing it.
- **Interactive is primary, video is derived.** Learner control (pause/scrub/reverse/drive)
  dissolves the transience problem; linear video reintroduces it. Dependency runs
  kernel → explorable → video, never the reverse.
- **Animate only where change is the referent.** Motion earns its place when the thing taught IS a
  change (contour deforming past a Borel singularity, Stokes line sweeping, trajectory peeling off
  a repelling manifold, caustic reorganising across the ladder). Structural/relational content
  (what the Borel plane *is*, how charts fit together) stays static.
- **Apprehension constraint.** Accurate but unreadable motion teaches nothing. Correct ≠ readable.

## Gate structure
Nothing downstream starts until the reviewer rules on the current gate.

| # | Chunk | Gate |
|---|-------|------|
| 1 | Synthesis → architecture + module sequence | user approves both |
| 2 | Probe protocol | **DONE — agreed** |
| 3 | Kernel spine + storyboard stills (no animation) | user answers all five probes |
| 4 | Promote to interactive explorable | user drives it, reports where it misleads/stalls |
| 5 | Video derived from validated explorable | user confirms it adds something, or video dropped for that module |
| 6 | Retrospective before scaling | user calls direction: continue / change tooling / change sequence |

Gate 6 exists so the pilot may legitimately conclude "this tooling is wrong" without that counting
as failure. Explicit anti-tunnel-vision valve.
