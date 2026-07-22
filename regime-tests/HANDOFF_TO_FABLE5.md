# Handoff note — coupled-FHN cusp-escape program → Fable 5

_From the outgoing agent (Opus), June 2026. Read this, then `PROGRESS.md`, then the two consolidated writeups
(`COUPLED_CUSP_RESULTS.md`, `PROCESS_CHAR_AND_W_DEFINITION.md`). You're inheriting a mature program with a clear
open frontier and a PI (Solomon) who holds a high bar for honesty. Everything below is oriented to letting you
act on day one._

---

## 1. What this is

Noise-induced escape at a **cusp** in a coupled FitzHugh–Nagumo system. Antisymmetric mode → cusp normal form
$V(Y)=\operatorname{sign}(Y)Y^2$; coupling $\Delta(g)=2\sqrt{-2g/3}$, $g_{\rm crit}\propto\sqrt\varepsilon$. Two
theorems drove the whole program:

- **T1** — a uniform parabolic-cylinder-with-noise tube through the fold→cusp merge. **CLOSED modulo one cited
  regularity** (Malliavin-differentiability of the non-degenerate first-passage, Nualart). Every *analytic*
  load-bearing step is proved or cited; the one item I did **not** re-derive for the swept-turning geometry is
  that Malliavin regularity — flagged, not waved through.
- **T2** — the new cusp edge law $\mathcal W$ (the Tracy–Widom analogue). **Fully characterized** (below); no
  closed form (*proven*), now given a canonical defining equation.

**Tracker snapshot** (`coupled-atlas/PROGRESS.md`): A 100 · **B 98** · C 80 · D 86 · E 84 · F 97 · overall ~88%.

---

## 2. The 𝒲 story — what is SETTLED (this is the spine; don't re-litigate it)

𝒲 went from "unidentified edge law we kept failing to fit into RMT" to a **fully named object**. The negatives
are *proven*, not guessed — do not re-attempt these:

- **NOT spectral/determinantal.** Not a soft-edge Fredholm gap (Stage 1: θ-dependent, two kernel forms, vs a
  pipeline verified on Airy→TW₂). Not a Painlevé σ-form (isomonodromy **fixed-point** obstruction: the only
  class deformation, the linear $\Delta Y$ term, is *relevant* → flows to fold). Not a Gamow survival (dominant
  resonance gives a Rayleigh bulk, skew $-0.63$ vs 𝒲 $+0.60$).
- **IS a first-passage law of an unbounded-below operator.** Decomposition, all named:
  **right tail = Kramers** (barrier LDP, heavier), **left tail = persistence** (literally $P(\text{no zero of
  the stochastic Weber field }u)$, exponent anomalous $\sim\beta^{0.6}$), **bulk = β-family**.
- **Defining equation** (`PROCESS_CHAR_AND_W_DEFINITION.md`): $\mathcal W_\beta$ = first-explosion law of the
  Weber Riccati; CDF solves the backward-FP PDE $-\partial_Y G+(\operatorname{sign}(Y)Y^2-p^2)\partial_pG+\frac2\beta\partial_{pp}G=0$
  (validated by `fp_cusp.py`). Provably **irreducible to a Painlevé ODE** (the fixed-point result).
- **The intrinsic process EXISTS and is RIGID** (`PIECE1_RIGIDITY_GO_NOTES.md`): the node process of $u$ is
  **class-I hyperuniform** (bounded number variance — stronger than GUE log), **strong repulsion** (spacing
  CV 0.08), **β-dependent** ($\mathrm{Var}(N)\propto1/\beta$). 𝒲 is its **edge**. Mechanism: 1-D spectral
  rigidity of the swept operator (NOT Dyson-determinantal).
- **The resonances** (complex-scaling, `M1_COMPLEX_SCALING_NOTES.md`): $\lambda_0=0.86-0.82i$, string
  $2.30-1.22i,\,4.14-1.08i$, θ-independent (two methods). These are **deterministic operator data, not the
  noisy law** — don't confuse them with 𝒲.

---

## 3. The precise open frontier (prioritized — this is your job)

1. **The exact node-process 2-point kernel / pair correlation** (biggest; shared with rungs A & D). Have:
   rigid, class-I hyperuniform, β-dep, repulsion. Missing: the *explicit* correlation function. **Route:**
   measure the structure factor $S(k)$ and pair correlation $g_2$ of the node point process directly (cheap,
   reuse `characterize_process.py` node positions); identify the small-$k$ power (hyperuniformity class);
   connect to 1-D-random-operator / hyperuniform-process theory. Cheap FP-validated check exists.
2. **The β-dependent persistence exponent** (bounded but genuinely hard). Left tail = persistence of $u$;
   exponent scales **anomalously $\sim\beta^{0.6}$** — the naive drift-diffusion first-passage ($\propto\beta$)
   is **falsified**. **Route:** non-stationary Gaussian-process persistence theory (Prüfer-phase
   stationarization → Rice / independent-interval / Bray–Majumdar–Schehr). The $\beta^{0.6}$ is the target.
3. **T1's cited regularity** — verify Nualart Malliavin-differentiability in the swept-turning geometry to make
   T1 fully self-contained.
4. **Breadth** (lower priority): C higher-ladder $q\ge3$ processes; D/E remaining synthesis. All numerical.

**Honest ceiling:** the last ~2% of B is characterisation of objects now in hand; **100% in the closed-form
sense is provably unreachable** (that's a theorem about 𝒲, not a missing calculation). Don't pretend otherwise.

---

## 4. How Solomon works (honor these — they are load-bearing)

- **Tag everything**: proved / cited / numerically-validated / conjectural. Flag load-bearing steps (⚑).
- **Never dress a numerical check as a proof.** Never declare a theorem closed unless it is. Name residuals
  *precisely*.
- **Don't grind a wall.** If a route resists — the τ won't close, the exponent won't match — **stop, diagnose
  precisely, map the path** (minimal missing ingredient / alternative / cleaner sub-problem). *A precise
  obstruction + a map is a fully acceptable outcome, often better than a forced calculation.* This is the single
  most important cultural rule; the best results in this program were clean diagnoses of why something *can't*
  work.
- **Cheap decisive checks** against the FP-𝒲/MC we trust, every time. Prefer **non-confounded** signatures
  (θ-independence, β-independence, tail mechanism, number-variance class) over bulk cumulant-fitting (which is
  *confounded* — a smooth Gaussian "passes").
- **Retract overclaims openly.** This program has retracted the "8/3" right-tail constant, the "$e^{-|s|^5/20}$"
  left tail, the σ-form regression, the "$\eta^2$ collapse", and the crude $\lambda_0=1.64-1.49i$. Each was the
  right call and *raised* trust. Do the same.
- Save notes + figures alongside the coupling notes (`regime-tests/`, `coupled-atlas/`); update `PROGRESS.md`
  with **honest** percentages (don't inflate).

---

## 5. Practical / sandbox notes

- **numpy only, NO scipy** (in this sandbox). I hand-rolled complex Gamma (Lanczos), Airy (recessive
  integration + normalized asymptotic IC), parabolic-cylinder (via the ODE). **If your environment has scipy,
  use it** — `scipy.special` (airy, pbdv/pbvv parabolic cylinder, gamma, gammaln) would remove a lot of my
  scaffolding, and the tail/rigidity sims were **sample-limited** by the compute budget, so more compute buys
  real signal.
- **bash 45s timeout**: vectorize heavy sims over realizations (loop only over time steps); split runs; renormalize
  to avoid overflow. `pip install ... --break-system-packages`.
- **matplotlib mathtext is finicky**: `\sqrt{...}` needs braces; no `\checkmark`; avoid `tight_layout` with a
  mathtext suptitle (use `fig.subplots_adjust`).
- The FP solver `fp_cusp.py` is the trusted ground truth for 𝒲 (cumulants skew $+0.607$, exk $-0.237$); the
  Riccati MC `escape()` pattern (in many scripts) is the trusted sampler.

---

## 6. Recommended first move

Do **the node-process structure factor / 2-point kernel** (frontier item 1). It's the deepest remaining object,
it's *ripe* (the GO just proved it exists and is rigid), it's shared uplift for rungs A & D, and it has a cheap
decisive check. Then the persistence exponent (item 2). Both build directly on validated ground.

Good luck. The physics is real, the standards are high, and the honest negatives are as valuable as the
positives here. — Opus
