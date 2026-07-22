# Notes — the cusp successive-peel-off process (Weber analogue of Airy₂)

_June 2026. Companion to the coupling notes. Figure `coupled-atlas/figures/forced_cusp_process.png`;
script `forced_cusp_process.py`. Tags: [NUMERIC] · [HEURISTIC] · [OPEN]._

**Result.** The cusp rung now has its multi-point object, completing the spine
(*singularity → edge class → process*) symmetrically with the fold:

| rung | marginal | multi-point process |
|---|---|---|
| fold | Tracy–Widom (TW_β) | Airy₂-type (forced) |
| **cusp** | **Weber** | **Weber-marginal process (this)** |

**[NUMERIC]** (`forced_cusp_process.py`, β=2). Running successive forced passages of the cusp inner
equation (k=2, parabolic-cylinder) with OU-correlated noise across passages:
- **stationary Weber marginal** — per-passage skew ≈ +0.57 (Weber-class), constant across passages,
  clearly distinct from the fold's TW marginal (skew ≈ +0.17);
- a **genuine process** — autocovariance decays C(1)=0.76 → C(8)=0.20 under forcing, vs the iid
  baseline C(1)≈0 (no forcing); increments saturate at 2·Var (0.90×) — the same process signatures
  as the fold's Airy₂-type result, now on the Weber marginal.

**[HEURISTIC / honest nuance].** At a *fixed imposed* OU correlation ρ, the cusp and fold processes
have **similar covariance** C(s) (decorrelation lag ≈ 8 for both) — the imposed ρ dominates the
covariance, so the distinguishing feature is the **marginal** (Weber vs TW), not the covariance.
This mirrors the fold situation exactly: my forced fold→Airy₂ result also used an imposed OU as a
surrogate for the forcing.

**[OPEN].** The *intrinsic* cusp process — the Weber analogue of the genuine Airy₂ process (the
Dyson-Brownian-motion edge), with its own covariance not inherited from an imposed ρ — remains open,
exactly as the exact Airy₂ kernel remains open for the fold. The intrinsic distinguishing covariance
would come from a "Dyson–Weber" edge dynamics (the multi-line ensemble for the parabolic-cylinder /
multicritical operator). Also open: tying this to the genuine coupled-FHN *forced* cusp (asymmetric
coupling), as opposed to the inner-equation surrogate.

**Net.** The spine is now complete and symmetric: each singularity rung (fold, cusp) has both its
edge marginal (TW, Weber) and a multi-point process (Airy₂-type, Weber-type), connected by the
coupling-driven fold↔cusp crossover. The remaining frontier on both rungs is the same kind of object
— the *intrinsic* (DBM-edge) kernel/covariance, which is the deep RMT-side theorem.
