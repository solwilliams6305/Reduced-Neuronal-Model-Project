# PROGRAM2 — TW_β halo: √(2/(3π)) derived as a Gelfand–Yaglom zero-mode Jacobian

**Date:** 2026-07-22. **Executes:** the residual flagged in `PROGRAM2_TWBETA_XEXTENSION_FACTOR_NOTES.md`
— derive the x-extension factor `√(A_g/π)=√(2/(3π))` from the fluctuation (Gelfand–Yaglom) determinant.
**CLOSES it.** **Memory:** `program2-tw-beta-halo.md`.

---

## Result (one line)

$$\boxed{\;\sqrt{\tfrac{2}{3\pi}}=\sqrt{\tfrac{A_g}{\pi}}=\sqrt{\frac{\lVert\psi_0\rVert^2}{2\pi}}\;,\qquad
\psi_0=\frac{dP_\star}{d\tau}=-\operatorname{sech}^2\tau,\quad \lVert\psi_0\rVert^2=\!\int\!\operatorname{sech}^4\!\tau\,d\tau=\tfrac43=\Delta V=2A_g,\;}$$

the standard collective-coordinate (zero-mode) Jacobian of the **tanh-kink escape instanton**. It
assembles the full Hastings–McLeod Stokes constant exactly:
`C = S₀·√(‖ψ₀‖²/2π) = (1/π)·√(2/(3π)) = √(2/(3π³))` — verified to 40 digits (`|Δ|=0`).

## The derivation (`coupled-atlas/_tw_zeromode_gy.py`)

The reduced cubic-barrier escape (`PROGRAM2_TWBETA_S_OF_A_NOTES.md`) is `dP=-(P²-1)dτ+2ε_eff dW`,
`V(P)=P³/3−P` (well `P=+1`, barrier `P=−1`, height `ΔV=V(−1)−V(+1)=4/3=2A_g`, `A_g=2/3`).

1. **Instanton = tanh kink.** The uphill (instanton) trajectory `dP/dτ=V'(P)=P²−1` is solved by
   `P⋆(τ)=−tanh τ` (`P:+1→−1`); verified `max|dP⋆/dτ−(P⋆²−1)|=8.6×10⁻⁴²`. Its translation zero mode is
   `ψ₀=dP⋆/dτ=−sech²τ`.
2. **Zero-mode norm = barrier height (virial).** `‖ψ₀‖²=∫_{−∞}^{∞}sech⁴τ\,dτ=4/3`. This equals
   `∫_{+1}^{−1}V'(P)\,dP=ΔV` (the instanton virial identity `∫ψ₀²dτ=∫(dP/dτ)dP=∫V'dP=ΔV`). Both `=4/3`
   to machine precision. So `‖ψ₀‖²=ΔV=2A_g`.
3. **Gelfand–Yaglom zero-mode Jacobian.** A fluctuation determinant with a bosonic zero mode is
   regularized as `det'(M)` (zero eigenvalue removed); trading the zero mode for its collective
   coordinate inserts `‖ψ₀‖/√(2πℏ)` in place of that mode's Gaussian (Coleman, *The Uses of
   Instantons*). Its ℏ-independent (Stokes-constant) part is
   `J_zm=√(‖ψ₀‖²/2π)=√((4/3)/2π)=√(2/(3π))=√(A_g/π)=0.4606588659617806…` (matches to `5×10⁻⁴²`). The
   accompanying `ℏ^{1/2}` (`=g^{−1/2}`) is realized as the **−½ large-order index shift** — one
   continuous zero mode, exactly the `γ'=0→γ=−½` shift measured in `_tw_frozen_vs_hm.py`.
4. **Assembly.** `C = S₀·J_zm`, with `S₀=1/π` the **transverse** (Kramers, P-direction) fluctuation
   determinant [`= 1/(width_barrier·width_well)`, the frozen escape *rate*], and `J_zm=√(2/(3π))` the
   **escape-instanton translation** zero mode. `C=(1/π)√(2/(3π))=0.14663227119384847789=√(2/(3π³))`,
   `|C−C_HM|=0`.

## Why the zero mode is in C but not in the frozen S₀

The frozen escape *rate* `R(g)` is translation-invariant in `τ`, so its Kramers normalization
(per-unit-`τ`) **factors the tanh-kink zero mode out** — leaving the transverse constant `1/π`. The
full **level-variable** (Hastings–McLeod / PII) resurgence is a *probability*-type object that **counts
that continuous zero mode in**, multiplying by its Jacobian `√(2/(3π))` and shifting the index by `−½`.
So the frozen→full step is precisely: restore the escape-instanton translation zero mode. This is the
same mode the O1 note removes for the `a^{−3β/4}` prefactor — here it is weighed, giving the Stokes
factor.

## Status

- **DERIVED + verified (40 digits):** `√(2/(3π))=√(‖ψ₀‖²/2π)`, `‖ψ₀‖²=∫sech⁴=4/3=ΔV`, the tanh-kink
  zero-mode Jacobian; `C=S₀·√(2/(3π))=√(2/(3π³))` exactly. This **closes** the residual of
  `PROGRAM2_TWBETA_XEXTENSION_FACTOR_NOTES.md` (the `√(A_g/π)` vs `√(A_g/2π)` question: it is
  `√(A_g/π)` because `‖ψ₀‖²=ΔV=2A_g`, i.e. the `2` is the virial factor, not a boundary artifact).
- **Interpretive content (standard, not new machinery):** the Coleman collective-coordinate formula;
  the rate-vs-probability reading of why `S₀=1/π` omits the mode. The *numbers* are all exact.
- With this, the T2c a-dependence is fully accounted: `C=√(2/(3π³))` = (transverse Kramers `1/π`) ×
  (escape zero-mode Jacobian `√(2/(3π))`), both derived. The remaining program-level question is the
  literature-flagged one: whether the **noise `1/β`** expansion (not the level variable) is itself a
  resurgent trans-series — un-established (CLDS treat it as ordinary cumulants).

## Files
- `coupled-atlas/_tw_zeromode_gy.py` — the full derivation + 40-digit verification.
