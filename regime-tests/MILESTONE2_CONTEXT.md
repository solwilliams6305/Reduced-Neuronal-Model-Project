# Milestone 2 context brief — the fluctuation→Borel-phase link (why ~54° not 45°)

_Shared brief for all Milestone-2 work agents (2026-07-10, Fable 5). Self-contained. Repo root:
`/Users/solomonwilliams/Reduced Neuronal Model Project`. Notes in `regime-tests/`, code in `coupled-atlas/`._

## The one-sentence goal
Explain why the **stochastic Borel-singularity phase is ~54–63°**, whereas the deterministic connection root
λ₀ = 0.8896 − 0.8896i sits at **arg = −45°**. This is "Thread 4 / the λ₀→Borel bridge," flagged multi-year-frontier.

## What EXACTLY the ~54° is (the target the theory must hit)
- Object: the argument of a **complex-conjugate pair of Borel singularities** of the weak-noise **variance**
  series `Var(Y*) = η²(v0 + v1 η² + v2 η⁴ + v3 η⁶ + ...)`. The Borel transform is `B(t) = Σ vn tⁿ/n!`; its
  nearest complex singularities sit at `ζ ≈ 1.2–1.44 · e^{±iθ}`, θ ≈ 54–63°.
- There is ALSO a competing **real** Borel pole at ζ≈1.5, arg 0 = the far-tail Freidlin–Wentzell instanton
  `S = s⁵/10` (real, BVP-verified in `coupled-atlas/instanton_action.py`).
- **Radius clue:** |ζ| ≈ 1.2–1.44 ≈ |λ₀| = 1.258. The Borel-pair radius MATCHES the resonance-root modulus;
  only the phase is rotated (45°→54°, a ~9–18° shift). Strongly suggests ζ = *noise-dressed λ₀*, rotated.

## The coefficient ladder (all that exists), from the chaos/Wick engine (MC-free)
`v0=0.134, v1=0.111, v2≈0.100, v3≈−0.02 (exact, delicate cancellation, ±0.02), v4≈−0.45, v5≈−1.1 (symbolic)
/ −2.29 (grid n=10)`. Definition: `v_k = Var(Y_{k+1}) + 2 Σ_{j=1}^k Cov(Y_j, Y_{2k+2−j})`, the Y_n being
η-independent Wick functionals of the field perturbations at the first node. **v6 is compute-BLOCKED** by the
7-chaos (7!) dense-tensor symmetrization wall; a transfer-operator rewrite (`chaos_transfer.py`) is partial.
Renormalization: v3+ need a δ(0) boundary self-contraction counterterm.

## How the ~54° was measured (the tools)
`coupled-atlas/_borel_analysis.py`:
- `borel_pade_poles(v,L,M)` → poles of [L/M] Padé of the Borel transform. Result: ζ≈1.44 e^{±i54°} (53.4–54.6°).
- `fit_darboux(v,...)` → fits `vn ~ 2C ζ^{-(n+1)} Γ(n+1+α) cos((n+1)θ−φ)`. Result: |ζ|≈1.17–1.20, θ≈62–64°.
- The two methods DISAGREE by ~9° (54 vs 63). The earlier "45°" was HARDWIRED (a `cos(nπ/4+φ)` fit with π/4
  imposed), NOT a floated measurement — so it is NOT independent evidence for 45°.

## The deterministic 45° baseline (closed form, machine-verified)
Resonance = root of the closed-form Weber Γ-equation `e^{3iπ/4} G_i(λ) = G_r(λ)`,
`G(μ)=Γ(3/4−μ)/Γ(1/4−μ)`, `G_i=G(iλ/4·...)`; root λ₀=0.8896−0.8896i (arg −45°). Solver
`coupled-atlas/connection_closed_form.py`. The inversion Stokes phase is e^{3iπ/4}=135°.

## Milestone 1 result (DONE this session) — the fluctuation data to build the bridge ON
`coupled-atlas/stochastic_stokes_o_eta2.py`, note `PROGRAM2_STOCHASTIC_STOKES_O_ETA2_NOTES.md`.
The noise-averaged connection root `E[λ] = λ0 + η² Ω2`, with **Ω2 = −0.451 + 0.352i** (|.|=0.572, arg 142°).
- The O(η²) MEAN shift rotates λ only −3.2°/η² (stays ≈−45°); reaching −54° needs η²≈2.8 ≫ radius 0.6. So the
  phase shift is NOT a mean effect — it lives in the FLUCTUATION.
- Fluctuation of the resonance root (per η², from 2nd Wiener chaos): **pseudo-variance E[δλ²]/η² = 0.436 − 0.181i**
  (arg −22.5°, anisotropy axis −11.3°); **true variance E|δλ|²/η² = 1.784** (rms|δλ| = 1.34 η). Built from
  `A_R=u_R(0)⁻⁴∫u_R⁴ = 0.305−0.280i`, `A_L=u_L(0)⁻⁴∫u_L⁴ = 0.017+0.414i` (Borel-regularized), and the mean
  2nd-chaos `M_L−M_R`. Certified across grid/cutoff/rotation; validated 0.67% vs exact Newton.
- Since the v_n ARE the variance (fluctuation) coefficients, the 54° must come from the large-order growth of the
  second-and-higher chaos — consistent with the Bureković–Schäfer–Grauer lesson (arXiv:2401.16264): in weak-noise
  problems the one-loop Gaussian-fluctuation prefactor, not the bare action, controls the noise-induced tail.

## The candidate mechanisms for the ~9–18° rotation (to be tested)
1. **Robustness first:** is 54° even robustly ≠ 45°, given only 6 coefficients (v5 rough) and a 9° method spread?
   Quantify the phase ± error honestly (bootstrap over v5, v3; multiple summation methods).
2. **Complex instanton + one-loop determinant (Bureković–Grauer, sub-attack 2):** compute the complex saddle
   action of the weak-noise Weber problem AND its one-loop fluctuation determinant; does the determinant rotate
   the effective phase 45°→54°? Tools: `instanton_action.py`, `mc_instanton_check.py`, `PERSISTENCE_ITEM2_NOTES.md`.
3. **Resonance-cloud saddle (the λ₀→Borel bridge, novel, builds on Milestone 1):** if the Borel singularity is the
   noise-dressed resonance ζ=λ[Ẇ], the large-order growth of the noise-AVERAGED v_n ~ E[λ^{-n}] is set by a
   SADDLE over the noise distribution (NOT the mean λ₀). Using the Milestone-1 covariance (pseudo-var 0.436−0.181i,
   true var 1.784), compute the dominant noise-saddle and its arg. Does it rotate 45°→54°?

## Ruled-out (do NOT re-attempt) — from PROGRAM2_STOCHASTIC_WKB_DEEPRESEARCH.md
- Complex-scaling a white-noise realization (ill-defined). FP-generator complex spectrum sits at 2–10°, not 54°
  (local, not the global connection object). PDE small-η re-extraction (grid-wall at width ~η).

## Files
Engine: `chaos_diagram.py` (run `python chaos_diagram.py v4` → v0..v4), `chaos_transfer.py`, `_v6_driver.py`,
`_vk_smallgrid.py`, `_v6_breakdown.py`. Analysis: `_borel_analysis.py`, `transseries_orders.py`.
Instanton: `instanton_action.py`, `mc_instanton_check.py`. Connection: `connection_closed_form.py`,
`complex_scaling.py`. Milestone 1: `stochastic_stokes_o_eta2.py`. Notes: `PROGRAM2_CONSOLIDATION.md`,
`PROGRAM2_CHAOS_ENGINE.md`, `PROGRAM2_ROUTE2B_NOTES.md`, `PERSISTENCE_ITEM2_NOTES.md`.
