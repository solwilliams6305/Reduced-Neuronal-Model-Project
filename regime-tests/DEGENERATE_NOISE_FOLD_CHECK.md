# Degenerate / fast-only noise at a fold — diligence verdict

*Coverage: BGK 1312.6353 §3 assumptions (grep of the PDF text); arXiv:2512.10460
abstract; a search of the hypoelliptic / Lévy slow-fast exit literature. I did NOT
read 2512.10460 in full, BGK's ref [6], or check the 3D folded-node fast-only case
specifically. Verdict bounded by that.*

---

## Verdict: route deflates. The FHN-relevant case is the *studied* case.

Three independent signals, all pointing the same way:

1. **BGK assume uniform ellipticity, but already remark the degenerate extension.**
   §3 of 1312.6353 states a "uniform ellipticity assumption" on the diffusion
   coefficients — but **Remark 3.1**: *"most of our results remain valid under a
   weaker hypoellipticity assumption"* (citing their foundational ref [6, p.175]).
   So the degenerate/hypoelliptic case is explicitly gestured at, and the
   foundational Berglund–Gentz machinery is built to handle it.

2. **The decisive one — arXiv:2512.10460 (Dec 2025).** This paper does noise on the
   *slow* variable at a fold and says, in its own framing:
   > *"the effect of adding noise to the **fast** variable, which is important for
   > noise-induced tipping, **has been previously analysed in detail**."*
   Noise on the fast variable at a fold **is** the FHN case (FHN noise is in `v`,
   the fast variable). So option 1 is the well-trodden noise-induced-tipping case,
   not open ground. (They also handle the complementary slow-noise case — with
   Airy functions — closing the other side of the noise-placement question.)

3. **Lévy / α-stable slow-fast exit is also active** (multiple 2022–2025 papers:
   heavy-tailed Lévy exit, Kramers under accelerated Lévy noise, LDT for Lévy
   diffusions in the small-noise regime). So option 2 is more worked than the
   earlier "under-worked" framing suggested too.

## Consequence

- **Option 1 does not retain a clean A1 ceiling.** The FHN-relevant case is the
  studied one; what's left (rigour for BGK's *remarked* hypoelliptic extension; a
  specific 3D folded-node fast-only variant) is A2-grade gap-filling at best, and
  unconfirmed even as that.
- **This is the third high-ceiling bet to deflate under the same diligence:**
  c₀ → BGK's; the quasipotential-bridge headline (saddle avoidance) → Börner's;
  degenerate fold → the noise-tipping literature's. Same structural cause each
  time: a mature field, general foundational machinery (Berglund–Gentz–Kuehn,
  Wechselberger), and an outsider vantage from which "untouched" is unreliable.

## Honest implication (the meta-finding)

External literature-scoping is not reliably surfacing a *90+-ceiling-and-genuinely-
open* problem in this sub-area — it keeps surfacing occupied ground, for a
structural reason. The binding constraint is **insider knowledge of what is open
and tractable**, which searches cannot reveal (it lives in unpublished work,
in-progress theses, and expert judgement).

**Reliable source = the supervisor.** Popović (and via him Kuehn) can hand over a
problem that is *both* genuinely open *and* known to be reachable — something three
rounds of solo scoping have not produced. The high-expected-value move is to change
the *method of finding the problem*, not to scope a fourth candidate (which, on this
evidence, would likely deflate too).

## If you still want to close option 1 fully

Read 2512.10460 in full + BGK [6] + check the 3D folded-node fast-only case. But the
prior is now strongly against an A1 remnant, and the pattern argues for the supervisor
route instead.
