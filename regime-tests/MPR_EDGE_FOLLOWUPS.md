# Two follow-ups: macroscopic SNIC search, and collective Tracy–Widom at the fold

*Companion to `MPR_FINITE_SIZE_LANGEVIN.md`. Both follow-ups run into the same honest conclusion: the
adaptive QIF mean-field naturally sits at the **amplitude** edge, and confirming collective Tracy–Widom
needs a structural reduction, not a histogram. Tags **[R]/[N]/[H]**.*

## Follow-up 1 — search for a macroscopic SNIC (phase edge): NOT FOUND  [N]

The collective phase-law prediction (timescale $\propto N^{-1/3}$) requires a macroscopic **SNIC**:
the collective burst period must **diverge** at the bifurcation. Scanning MPR+adaptation
($\alpha\in\{3,5,8\}$, $\tau_a\in\{15,40\}$, $\Delta\in\{0.3,1\}$, $J=15$) and resolving each oscillation
boundary finely (long runs, $T$ up to $9000$):

- The collective bursting **always ends at a fold of limit cycles** — period stays **finite**
  ($\sim17$–$21$) and amplitude **steady** ($\sim2.1$) right up to the boundary, then the cycle vanishes
  abruptly (verified at $\bar\eta^*\approx-2.3127$ to four digits: period $18.0\to18.5$, then gone — **no**
  log-divergence). Slower adaptation ($\tau_a=40$) raises the max period only to $\sim21$, still finite.
- **No macroscopic SNIC/homoclinic was found** in the scanned region. *The network realises the
  amplitude edge:* the collective oscillation amplitude (the synchrony level) is the soft mode that
  collapses, while the period stays finite — the opposite of the single-neuron SNIC.

So the $N^{-1/3}$ collective **phase**-edge test is moot at these operating points — there is no phase
edge. (A macroscopic SNIC may exist in a structurally different network — inhibitory/synaptic/delay
coupling, or a different mean-field — but not in excitatory MPR+adaptation as scanned.) Honest negative,
and itself informative: **next-gen neural-mass collective rhythms here are amplitude-edge objects.**

## Follow-up 2 — collective Tracy–Widom at the fold: SUGGESTIVE, NOT CONFIRMED  [N] / [H]

Since the network *is* at the amplitude edge (fold of cycles), this is the paper's headline TW channel.
Two finite-size-noise measurements near the fold:

- **Ramp through the fold** (peel-off $\bar\eta$ at collective collapse, $500$ trials): dominated by the
  **dynamic-bifurcation delay** and ramp cutoff (collapse clusters with std $\sim4\times10^{-4}$), skew
  $+0.36$ — not a clean noise-driven peel-off.
- **Just above the fold** (per-trial peak amplitude, $400$ trials): skew $\mathbf{+0.36}$ — the **same
  sign** as $\mathrm{TW}_1$ ($+0.29$), $\mathrm{TW}_2$ ($+0.22$) and unlike a symmetric Gaussian — but
  the fluctuations are **small and Gaussian-dominated** (CV $\sim0.6\%$): not yet in the edge regime.

**Verdict: a consistent positive skew (TW sign), but no Tracy–Widom confirmation.** A histogram cannot
settle this — the decisive question is **structural**: does the collective amplitude equation near the
fold reduce to the **noisy Riccati/Airy canard** (the paper's amplitude-channel inner equation,
Cole–Hopf $\to$ stochastic Airy operator $\to\mathrm{TW}_\beta$)? That is a *derivation*, the genuine
next step:

1. Center-manifold-reduce the noisy MPR$(r,v)$ Langevin (with the $O(N^{-1/2})$ noise of
   `MPR_FINITE_SIZE_LANGEVIN.md`) at the fold of cycles;
2. check whether the reduced amplitude dynamics is the Riccati canard $dY=(\dots-Y^2)\,ds+\eta\,dW$ of the
   paper's §"amplitude channel";
3. if so, $\beta$ is read off the collective noise-to-curvature ratio and the peel-off level is
   $\mathrm{TW}_\beta$ — **collective Tracy–Widom**, with system size $1/\sqrt N$ as the intensity.

## Combined status

- **[N]** Macroscopic edge in MPR+adaptation = **fold of cycles (amplitude edge)**, robustly; no SNIC
  found ⇒ the collective **phase**-edge / $N^{-1/3}$ test is not available here.
- **[N]** Near-fold collective fluctuations carry a **TW-signed positive skew** but are not yet
  edge-regime; **collective TW is not confirmed**.
- **[H] The real open task is analytic:** reduce the noisy MPR fold of cycles to the Airy/Riccati canard.
  If it reduces, collective Tracy–Widom follows (the paper's headline channel, at the network level); if
  it does not, the collective amplitude statistics are a *different* universality — either way a clean,
  citable result. This is where the theory track should go next, and it is a pen-and-paper reduction, not
  more simulation.
