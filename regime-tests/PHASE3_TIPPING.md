# Phase 3 — the early-warning demo: edge-law `ν̂` vs classical EWS on a drifting tipping system

*Closes Phase 3 of `INVERSION_MVP_ROADMAP.md`. Code: `phase3_tipping.py` (+ `phase3_figure.py`).
Figure: `figures/phase3_tipping.png`. Tags: **[R]** proved / **[N]** numerical / **[H]** heuristic.*

```
python3 phase3_tipping.py sims    # integrate + cache approach & null spike ensembles (~13 s)
python3 phase3_tipping.py bench   # the matched-FPR benchmark table (~9 s)
python3 phase3_tipping.py fig     # -> figures/phase3_tipping.png
```

## The setup

A SNIC neuron (QIF `dv=(v²+I)dt+σdB`, reset) whose input current `I(t)` drifts slowly **down** toward
the SNIC at `I_c=0`. It fires regularly, the firing slows as `ν=I/σ^{4/3}→0`, then ceases (the limit
cycle is destroyed). **The spike train is the only observable.** In sliding fixed-count ISI windows we
compute three early-warning signals from the *same* ISI data and ask which warns of the approaching
bifurcation earliest, and with the fewest false alarms:

- **`ν̂`** — the edge-theory inversion (off-grid SNIC shape fit): a calibrated, σ-free distance with a
  meaningful zero;
- **ISI Var** — interval variance (the classic critical-slowing-down indicator);
- **lag-1 AC** — lag-1 autocorrelation of the ISI sequence (the textbook CSD indicator).

Fairness: each trend method alarms on a **trailing-window trend** in its own signal, with the alarm
threshold **calibrated on a stationary-null ensemble to a matched 15% false-alarm rate**, so detection
and lead time are compared at equal FPR. `ν̂` additionally admits an **absolute** alarm (`ν̂<0.4`) that
the scale-free Var/AC1 structurally cannot.

## Result (120 trials each; true SNIC crossing at t=2000)

| method | null FPR | detection (approach) | median lead | **confound FPR** |
|---|---|---|---|---|
| `ν̂` trend | 15% | 50% | 792 | **11%** |
| ISI Var trend | 15% | 56% | 836 | **31%** |
| lag-1 AC trend | 15% | 11% | 987 | 10% |
| **`ν̂ < 0.4` (absolute)** | **0%** | **96%** | 240 | **0%** |

Four honest readings:

1. **[N] On a clean approach, the `ν̂` trend and the ISI-variance trend are comparable** (detection
   50% vs 56%, lead 792 vs 836 at matched FPR). We do **not** claim the edge-law *trend* warns earlier
   than rising variance on confound-free data — it doesn't, materially. Variance is a good CSD signal
   here, and we report that straight.

2. **[N] lag-1 autocorrelation of ISIs is a poor early-warning signal** (11% detection ≈ the 15%
   chance floor). The QIF spike train is approximately a **renewal** process, so successive ISIs are
   near-independent and there is little lag-1 structure to track. The textbook AC1 indicator is built
   for the *state variable's* critical slowing — information an ISI-only recording does not carry. A
   useful negative result, not a tuning failure.

3. **[N] The edge-law signal's robustness is the real payoff.** Under a **rate-drift confound** (a slow
   multiplicative rescaling of firing rate with **no** approach to the bifurcation — electrode drift /
   slow gain modulation), the variance trend's false-alarm rate **nearly doubles, 15%→31%**: it cries
   wolf at a change that is not an approach. `ν̂` is built from the **σ-free ISI shape**, which is
   invariant under that rescaling, so its confound false-alarm rate **stays at baseline (11%)**.

4. **[N] The absolute `ν̂` alarm is the standout, and is unique to a calibrated distance.** Because `ν̂`
   has a real zero, the rule `ν̂<0.4` gives **96% detection at 0% false alarms** — on *both* the
   stationary null and the rate-drift confound. Var and AC1 have no absolute scale, so no transferable
   absolute alarm exists for them; they need a per-record baseline and a trend test. The absolute
   rule's lead (240) is deliberately shorter than the trend rules': it is a **high-precision "you are
   now close" alarm** (`ν` has genuinely fallen near 0), not a trend extrapolation — hence the 0% false
   alarms.

## Bottom line

The value proposition of the edge-law inversion in a tipping setting is **not** a dramatically earlier
trend warning — classical ISI variance is a competitive trend detector on clean data. It is (i) a
**calibrated, transferable, near-zero-false-alarm absolute alarm** (`ν̂→0`) that scale-free EWS cannot
provide, and (ii) **robustness to rate/scale confounds** that inflate variance-based false alarms.
That is the defensible claim, and it is exactly what a bifurcation-*specific* early-warning signal
should buy over generic critical-slowing-down statistics.

## Honest limits

- **Count collapse near threshold.** As `ν→0` the firing rate falls, so per-window ISI counts shrink
  exactly where the estimate matters most; fixed-count windows keep the estimator stable but stretch in
  time there. `ν̂` gets noisier near the edge (visible in figure panel A) — the absolute alarm fires
  while counts are still adequate, which is part of why it is reliable.
- **`ν̂` saturates at high drive** (low-CV regular firing → the SNIC shape fit hits the atlas ceiling),
  so the *early* part of the approach is tracked only coarsely; the signal sharpens as `ν` falls.
- **Single normal form.** This is the SNIC/QIF (Type-I) tipping; the fold-of-cycles/Hopf (Type-II,
  amplitude-edge) tipping and a real published model (**Epileptor**, Phase 3b) are the natural next
  validations. The pipeline (stream ISIs → window → invert → alarm) is model-agnostic.
