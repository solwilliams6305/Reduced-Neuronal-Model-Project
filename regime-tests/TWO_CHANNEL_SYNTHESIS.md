# Two channels, two edges — the synthesis

*One picture for the whole arc. The noisy folded limit cycle has a single inner Riccati
`dR=(R²−Y)dT+ηdB` ⇔ `u''=(Y−ηξ)u`. Each destroying bifurcation activates one channel, and in
**both** channels the noise response is a **dynamic critical edge bracketed by an Arrhenius
(Kramers) ceiling and a regular-diffusion regime** — the same three-layer nesting, once with an
integrable edge and once with a non-integrable one. References: conversation-1 amplitude work
(`tw_boundary.py`, `three_scale_overlay.png`); `FOLDED_CYCLE_PHASE_EDGE.md`,
`phase_edge_snic.py`, `offcritical_phase_edge.py`.*

---

## The parallel

|  | **Amplitude channel** (fold of cycles) | **Phase channel** (SNIC) |
|---|---|---|
| read the Riccati by | **sweeping** `Y` | **holding** `Y` at criticality |
| observable | peel-off **level** `Y_node` | rotation **period** = first-passage **time** |
| inner operator | **Airy** (linear potential) | **quartic** `ρ⁴+2ρ` (anharmonic) |
| **dynamic edge** | `Y_node =d TW_β`, `σ_*∼√ε₂` (`η∼1`) | quartic FPT law, `ω,D_φ∼σ^{2/3}` |
| edge integrable? | **yes** — Painlevé II / `TW_β` | **no** — quartic, but mean `J` is closed (`Γ(1/6)`) |
| **Arrhenius ceiling** | recrossing over the unstable cycle, `σ²∼ΔU_*/log(1/ε₂)`, `ΔU_*∼` cycle→rest barrier | pre-SNIC double well (`ν<0`), `σ²∼(8/3)μ^{3/2}/log`, barrier `∝μ^{3/2}` |
| **regular regime** | (large-noise two-way / rest bistability) | post-SNIC single well (`ν>0`), `D_φ∼σ²`, period `π/√μ` |
| crossover knob | `η=σ/√ε₂` (clean TW for `η≲1`) | `ν=μ/σ^{4/3}` (SNIC edge at `ν=0`) |

The right-hand `Y^{3/2}`/`μ^{3/2}` Arrhenius barrier is the **same physics in both rows** — and
the same `Y^{3/2}` barrier as the original FHN `σ_crit(I)` table. So the project's starting
object (an Arrhenius critical-noise table) is the **ceiling** in *both* channels; the new
content is the **dynamic edge** that opens up beneath it.

## The three-layer nesting (both channels)

```
   large σ  ────────────────────────────────────────────────────────────
            Arrhenius / Kramers ceiling     σ² ∼ barrier / log
            (amplitude: recrossing;  phase: pre-SNIC double-well escape)
   ----------------------------------------------------------------------
            DYNAMIC EDGE  (bifurcation-specific, power-law noise scale)
            amplitude:  TW_β,  σ_* ∼ √ε₂        [INTEGRABLE  — Airy/Painlevé]
            phase:      quartic FPT, ω ∼ σ^{2/3} [NON-INTEGRABLE — anharmonic]
   ----------------------------------------------------------------------
            regular diffusion / deterministic limit
            (phase: post-SNIC, D_φ ∼ σ²;  amplitude: deep sub-threshold)
   small σ  ────────────────────────────────────────────────────────────
```

As the rate parameter → 0 the edge and the ceiling **separate** (the amplitude band widens as
`σ_*²∼ε₂` drops linearly while `σ_crit²∼1/log` drops only logarithmically — `three_scale_overlay.png`;
the phase edge `σ^{2/3}` likewise sits below the `μ^{3/2}` Arrhenius ceiling). So in each channel
there is a parametrically wide window where the **dynamic edge governs**, before control hands
back to the Arrhenius table.

## Why one edge is integrable and the other is not

Same operator, complementary spectral question:
- **sweep the spectral parameter** → the *extreme eigenvalue* of the stochastic Airy operator →
  `TW_β` (a **linear**-potential edge — Airy — hence Painlevé-integrable);
- **hold it and time the node spacings** → the *first passage / rotation number* → the
  **quartic** anharmonic FPT operator (no closed-form spectrum; only the mean `J` is a clean
  `Γ(1/6)`).

"**Where** the canard peels off" is Tracy–Widom; "**when** it rotates" is not. The destroying
bifurcation (fold of cycles vs SNIC) picks which question — hence which edge — is singular.

## One-line statement

> The noisy folded limit cycle carries one inner operator and two edges. Each is a dynamic
> critical scaling sitting beneath the same `Y^{3/2}` Arrhenius ceiling (the original `σ_crit`
> physics) and above a regular-diffusion floor. The amplitude (fold-of-cycles) edge is the
> integrable Tracy–Widom law; the phase (SNIC) edge is the non-integrable quartic first-passage
> law with closed mean `J=(√π/3)12^{1/6}Γ(1/6)`. Sweeping the operator gives the eigenvalue
> edge (TW); holding it gives the first-passage edge (quartic) — the same random Schrödinger
> operator read two ways.
