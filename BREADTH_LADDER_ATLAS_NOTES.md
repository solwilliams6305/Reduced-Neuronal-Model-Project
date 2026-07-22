# Breadth — the higher ladder (C), the regime atlas (D), and the unfolding response (E)

_June 2026. Figures `coupled-atlas/figures/ladder_edge_laws.png`, `regime_atlas.png`,
`unfolding_response.png`. Scripts `ladder_edge_laws.py`, `atlas_figure.py`, `unfolding_response.py`.
Organized by the higher-order-Tracy–Widom / multicritical (Painlevé-II-hierarchy) reframe (q = index k).
Tags [PROVED]/[NUMERIC]/[CITED]/[HEURISTIC]._

---

## C — the higher ladder: edge fingerprints for q = 1,2,3,4

Same Riccati-explosion engine as the cusp, with $V_q(Y)=\operatorname{sign}(Y)|Y|^q$. Standardized edge
law $\mathcal W^{(q)}_\beta$ at β=2 (η=√2), ~0.8–1.4M samples each:

| q | catastrophe | skew | exk | κ₅ | κ₆ | left-tail (fit/asymp) | right (fit/asymp) |
|---|---|---|---|---|---|---|---|
| 1 | fold | +0.20 | +0.02 | −0.15 | −0.50 | 1.6 / **3** | 1.3 / **1.5** |
| 2 | cusp | +0.61 | −0.24 | −2.12 | −2.71 | 2.5 / **5** | 1.8 / **3** |
| 3 | swallowtail | +0.96 | +0.17 | −3.25 | −9.43 | (sharp) / **7** | 1.8 / **4.5** |
| 4 | butterfly | +1.22 | +0.71 | −3.58 | −16.5 | (sharp) / **9** | 1.6 / **6** |

**Findings [NUMERIC]:**
- **A clean monotone hierarchy.** Skew rises $0.20\!\to\!0.61\!\to\!0.96\!\to\!1.22$; κ₅, κ₆ grow steadily
  more negative (κ₆ dramatically). The progression is the numerical signature of a *family* of edge laws,
  one per catastrophe degree — exactly what the reframe predicts.
- **The cusp (q=2) is uniquely sub-Gaussian** (exk −0.24); q=1,3,4 are super-Gaussian (exk ≥0). The excess
  kurtosis is non-monotone, dipping negative only at the cusp.
- **The left edge sharpens with q.** For q=3,4 essentially no samples populate the deep left tail — the
  distribution has a sharp left cutoff, *consistent* with the very steep $2q+1$ left exponent (7, 9).
- **Anchors:** q=1 gives skew +0.20 (TW-class; canonical TW₂ +0.224 — the small offset is the
  escape-location vs eigenvalue observable + finite dt), q=2 reproduces the cusp fingerprint.

**Honest caveat on exponents.** The tail exponents $2q+1$ (left), $3q/2$ (right) are the **Freidlin–Wentzell
analytic backbone** [CITED — derived earlier in the project]. Direct *numerical* confirmation is out of reach
(pre-asymptotic; the asymptotic regime sits past the accessible depth, and for q≥3 the left tail is too thin
to sample). So: the hierarchy/progression is numerically established; the exact exponent values rest on the
FW derivation, not on fresh fits.

---

## D — the regime atlas (the periodic table)

`figures/regime_atlas.png` lays out the full map **singularity → normal form → stochastic inner operator →
marginal edge law → cumulants → tail exponents → process → multicritical index**, color-coded by rung, with
the unfolding arrow. Organizing statement (the reframe):

> **q = k**, the index of the higher-order Tracy–Widom (Painlevé-II-hierarchy) family. The fold is k=1
> (KPZ/TW); the cusp is the new k=2 law $\mathcal W_\beta$; q=3,4 are the new k=3,4 laws. Each is the
> *dynamical* (noise-induced-escape) realization of a multicritical edge that random-matrix/2D-gravity
> ensembles realize *statically*.

**Status [SYNTHESIS]:** the marginal-law row is complete and (q≤2) backed by construction; the **limiting
multi-point process** is established only for q=1 (Airy₂) and sketched for q=2 (Weber process) — the q≥3
processes are **not built** (a remaining D gap).

> **SUPERSEDED — see `REGIME_ATLAS_V2_NOTES.md` + `figures/regime_atlas_v2.png` (June 26).** This §D and the
> original `regime_atlas.png` are organized by the *higher-order-TW (q=k)* reframe, which **T2 later refuted**
> (cusp tails (5,3), ratio 5/3≠2 — the cusp is an *asymmetric PIV-family isomonodromy*, not a higher-order TW
> law). The v2 atlas corrects the organizing principle (catastrophe ladder), adds the **proved folded-node**
> row (TW_β) and the **β-axis fold/cusp inversion**, and splits physical from abstract regimes. Use v2 as the
> current atlas; this section is retained for history.

---

## E — the unfolding response (perturbation classification)

> **SUPERSEDED — see `REGIME_UNFOLDING_V2_NOTES.md` + `figures/unfolding_response_v2.png` (June 26).** This §E
> unfolds only the abstract q=3 potential and is phrased in the *higher-order-TW* reframe (refuted by T2). The
> v2 notes unfold the **physical cusp** (linear → fold), add the **within-class symmetry-breaking
> susceptibility** (new), consolidate the **β / coupling / forcing** axes, and re-label the cusp response as
> the $\mathcal W_\beta$ asymmetric-PIV family. Use v2; this is retained for history.

Unfold the q=3 potential by its lower-degree terms (codimension 2):
$V(Y)=\operatorname{sign}(Y)(|Y|^3+c_2|Y|^2+c_1|Y|)$. Result (`figures/unfolding_response.png`):

- **$c_1$ on (linear/fold term): clean flow q=3 → fold.** skew $0.96\!\to\!0.93\!\to\!0.85\!\to\!0.56\!\to\!0.25$
  as $c_1:0\to6$ — lands on the fold/TW value. **[NUMERIC, clean]**
- **$c_2$ on (quadratic/cusp term): flow to the cusp *class*, but $\beta_{\rm eff}$ rises with $c_2$.** The
  skew rides *up* the cusp β-family (the inverted TW$^{(q)}_\beta$, skew increasing in β) rather than landing
  on the β=2 value 0.61, because the coefficient also rescales the effective noise. **[NUMERIC + caveat]**

**Classification [PROVED-level principle + NUMERIC demo].** Near the turning the **lowest-degree active term
controls the edge universality class** (it dominates $V$ as $Y\to0$). Hence:
- A *generic* perturbation (nonzero linear term) is **relevant** and flows the law to the **fold/TW** class.
- Reaching rung $q$ requires tuning away all lower terms — **codimension $q-1$** — exactly Thom's genericity
  of catastrophes, now read on the *stochastic edge law*.
- The unfolding coefficient is simultaneously a **noise-rescaling** (sets $\beta_{\rm eff}$), so within a
  class the fingerprint moves along the β-family.

This recovers, and generalizes, the cusp Δ-crossover (q=2→fold, studied earlier as the noisy Airy↔Weber
crossover) as the codim-1 case.

---

## Net

C, D, E substantially filled: the ladder is characterized through q=4 (clean cumulant hierarchy), the atlas
is assembled as a single labeled map under the reframe, and the perturbation response is classified
(lowest-term-dominates / codimension / β-rescaling). What remains: the higher-rung **multi-point processes**
(D), exhaustive multi-direction unfolding (E), and the **numerically-unconfirmable tail exponents** (C — the
FW backbone stands analytically). No new *rigor* was added here — this is breadth/synthesis; the depth
(analytic identification of $\mathcal W_\beta$) is the next target.
