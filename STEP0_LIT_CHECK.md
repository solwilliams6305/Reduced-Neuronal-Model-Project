# Step 0 — Confirmatory literature check (gate for the quasipotential-folds program)

*Run 2026-06-02. Companion to `QUASIPOTENTIAL_FOLDS_PROBLEM_STATEMENT.md` §0.
Scope: literature confirmation + the Popović email only. No derivation (that is §1+).*

---

## Verdict: GO — proceed to §1, with two framing amendments and one must-read.

The target intersection is **empty**: no one has computed the Freidlin–Wentzell
quasipotential *landscape* at a slow-fast folded singularity, and the standard
quasipotential solvers (Cameron OUM/OLIM) have **never been pointed at a fold**.
The bankable core (brief §9, steps 1–3) is unthreatened.

Two mature literatures sit adjacent and **must be positioned against explicitly**
— one of them (Börner et al. 2024) overlaps the *motivation* and should be read in
full before the "solvers break at folds" framing is claimed. One human check
remains open: the Popović email (drafted below).

---

## What the searches establish

**1. Solver side is open.** Cameron's OLIM line reaches anisotropic / position-
dependent diffusion ([arXiv:1806.05321](https://arxiv.org/abs/1806.05321)) and
highly dissipative / chaotic systems (stochastic Lorenz'63,
[arXiv:1809.09987](https://arxiv.org/abs/1809.09987)); the public software was last
updated 03/2024 ([UMD OLIM page](https://www.math.umd.edu/~mariakc/OLIM.html)).
All of it is **smooth, O(1) systems around fixed-point attractors** — no fold, no
singular limit, no ε-dependent mesh analysis. The OUM convergence theory
([arXiv:1601.02885](https://arxiv.org/abs/1601.02885)) is generic and says nothing
about loss of normal hyperbolicity. Searching "ordered upwind method stiff
multiscale" returns only generic singular-perturbation FEM (Shishkin / graded
meshes for convection–diffusion boundary layers) — unconnected to quasipotential.
**H1 (the mesh-resolution breakdown at `h ~ ε^{2/3}`) is untested.**

**2. The fast-slow quasipotential that *does* exist is a different object.**
Bouchet–Grafke–Tangarife–Vanden-Eijnden, *Large Deviations in Fast-Slow Systems*
(JSP 2016, [arXiv:1510.02227](https://arxiv.org/abs/1510.02227)) — abstract
confirmed — computes a quasipotential in the **averaging regime**: the fast
variable is a stochastic process *integrated out*, giving a Hamiltonian equal to
the leading eigenvalue of the fast generator (non-quadratic in momenta), with
turbulence/plasma examples. This is **not** the geometric fold, not loss of normal
hyperbolicity, not degenerate-noise-on-the-fast-variable FHN. **Framing amendment:**
the brief should distinguish "averaging quasipotential" (occupied by BGV-E) from
"geometric-fold quasipotential landscape" (open) up front, so a referee doesn't
conflate them.

**3. The one paper that overlaps the motivation — read it in full.**
Börner, Deeley, Römer, **Grafke**, Lucarini, Feudel, *Saddle avoidance of
noise-induced transitions in multiscale systems* (PRR 2024,
[arXiv:2311.10231](https://arxiv.org/abs/2311.10231)). Timescale separation +
non-gradient + neuroscience toy models + **"flat quasipotential"** causing sample
transitions to **deviate from the Freidlin–Wentzell instanton**, with an
Onsager–Machlup fix. This is the closest thing to the project's qualitative claim
that *timescale separation breaks the standard FW/instanton picture*. It does **not**:
(a) construct the closed-form fold barrier via blow-up, (b) benchmark grid-solver
mesh resolution / derive the exponent `q`, or (c) treat degenerate noise or the
`ε^{1/3}×ε^{2/3}` inner box. **Consequence for novelty:** the qualitative
"FW misbehaves under timescale separation" observation is partly claimed — so the
contribution must lean on the **quantitative** pieces (closed-form `ΔU_ε`, the
derived resolution law `q`, degenerate noise), not on the qualitative message
alone. Reading this paper end-to-end is the single highest-value next action after
the email.

**4. BGK boundary holds (brief §7 confirmed).** Berglund–Gentz–Kuehn own the
**sample-path / random-Poincaré-map** treatment of noisy folded nodes and MMOs
([arXiv:1312.6353](https://arxiv.org/abs/1312.6353), JDDE 2015; "Hunting French
Ducks", [hal-00535928](https://hal.science/hal-00535928)). Path-concentration, not
a global quasipotential landscape — position against, do not re-claim.

**5. Degenerate-noise axis appears genuinely untouched.** No quasipotential-solver
work on rank-deficient / hypoelliptic diffusion at a fold surfaced; the OLIM
anisotropic variant assumes full-rank (positive-definite) diffusion. Consistent
with the brief calling this a second, independent open axis.

---

## Pending: the one human check (highest value)

Kuehn's group is confirmed deep in fold/GSPT — Kaklamanos–Kuehn–Popović–Sensi,
*Entry–exit functions in fast–slow systems with intersecting eigenvalues* (JDDE,
online 2023 / vol. 37 pp. 559–576, 2025,
[Springer](https://link.springer.com/article/10.1007/s10884-023-10266-2)) — but on
the **deterministic/GSPT** side; their stochastic line is BGK sample-paths, not
quasipotential landscapes. So the email below is a genuine check, not a formality:
Popović would know of any unpublished quasipotential-near-fold work in that orbit.

### Draft email to Popović

> **Subject:** Quick check — quasipotential landscapes at folded singularities?
>
> Hi Nikola,
>
> One quick check before I commit to a direction. I'm looking at the
> Freidlin–Wentzell *quasipotential landscape* (the global `V` solving
> `H(x,∇V)=0`, the object Cameron's OUM/OLIM solvers compute) at a slow–fast
> folded singularity — and at whether those solvers survive the ε→0 limit there.
> As far as I can tell this falls in the gap between the quasipotential-numerics
> community (smooth, O(1) systems, fixed-point attractors) and the GSPT/canard
> side, where the stochastic treatment is Berglund–Gentz–Kuehn sample-path bounds
> rather than a global landscape.
>
> Before I invest: has Christian's group — or anyone you know of — already touched
> quasipotential-near-fold, or large-deviation *landscapes* (as opposed to
> sample-path bounds) for canards / folded nodes? I want to be sure I'm not
> re-treading something in the Kuehn orbit that hasn't shown up in my searches.
>
> Thanks!
> Solomon

---

## Net effect on the brief

- **Proceed.** Steps 1–3 (closed-form fold barrier → OLIM benchmark → H1 resolution
  law) are clear of prior art.
- **Add to §1/§7:** explicitly separate the *averaging* quasipotential (BGV-E 2016)
  from the *geometric-fold* landscape, and cite Börner et al. 2024 as the nearest
  motivation-neighbor — framing the novelty as the quantitative fold-scaling +
  solver-breakdown law, not the qualitative observation.
- **Before any derivation:** (i) send the Popović email; (ii) read Börner et al.
  2024 in full. Gate the "solvers break at folds" claim on both.
