"""
period_divergence.py
--------------------
Tests whether the bifurcation at I_SNIC is a true SNIC or a homoclinic
by measuring how the deterministic limit cycle period T(I) diverges as
I approaches the bifurcation from above.

Two candidate laws:
  SNIC:       T ~ C / sqrt(I - I_c)    [saddle-node on invariant circle]
  Homoclinic: T ~ -C * log(I - I_c)   [homoclinic orbit, logarithmic divergence]

The best fit determines which normal form applies, which in turn determines
the correct σ_crit scaling near the bifurcation.

Method
------
1. Find I_c (the bifurcation point) by bisection — smallest I for which a
   deterministic trajectory completes a full orbit within T_max.
2. Sweep I values logarithmically spaced just above I_c.
3. Measure T(I) by integrating the deterministic ODE and timing the period.
4. Fit both candidate laws in log-log and log-linear space.
5. Compare R² to determine which law fits better.

Usage
-----
    python period_divergence.py
    python period_divergence.py --eps 0.05 --quick
    python period_divergence.py --a 0.7 --b 0.5   # different parameters
"""

from __future__ import annotations
import argparse, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(__file__))
from kernel import FHN2D

os.makedirs("figures/bifurcation", exist_ok=True)
os.makedirs("data/bifurcation",    exist_ok=True)

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("--eps",   type=float, default=0.08)
parser.add_argument("--a",     type=float, default=0.7)
parser.add_argument("--b",     type=float, default=0.8)
parser.add_argument("--quick", action="store_true")
args = parser.parse_args()

eps = args.eps
a   = args.a
b   = args.b

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
if args.quick:
    n_I      = 15       # number of I values to test
    T_max    = 3000.0   # max integration time
    n_orbits = 5        # orbits to average over
else:
    n_I      = 30
    T_max    = 8000.0
    n_orbits = 10

dt        = 1e-3
threshold = 1.0    # spike detection at v = 1

print(f"\n{'='*60}")
print(f"  Period divergence test")
print(f"  eps={eps}, a={a}, b={b}")
print(f"  n_I={n_I}, T_max={T_max}, n_orbits={n_orbits}")
print(f"{'='*60}\n")

# ---------------------------------------------------------------------------
# Step 1: Find bifurcation point I_c by bisection
# ---------------------------------------------------------------------------

def measure_period(I_val: float, T_max: float = T_max,
                   n_orbits: int = n_orbits) -> float | None:
    """
    Integrate deterministic FHN and return mean period over n_orbits spikes.
    Returns None if fewer than n_orbits+1 spikes in T_max.

    Uses a long warm-up (T_max/3) to ensure the system is on the limit
    cycle before measuring — critical near the bifurcation where transients
    are extremely long.
    """
    def fhn(t, y):
        v, w = y
        return [v - v**3/3 - w + I_val,
                eps*(v + a - b*w)]

    # Start from a point that's roughly on the limit cycle
    # (left branch, near fold) to reduce transient
    y0 = [-1.0, -0.5]

    # Long warm-up — near bifurcation transients can be very long
    warmup = T_max / 3.0
    sol_warm = solve_ivp(fhn, [0, warmup], y0, method='RK45',
                         max_step=dt*20, rtol=1e-8, atol=1e-10)
    y0 = sol_warm.y[:, -1]

    # Measure period
    sol = solve_ivp(fhn, [0, T_max], y0, method='RK45',
                    max_step=dt*5, dense_output=True, rtol=1e-8, atol=1e-10)

    t_eval = np.linspace(0, T_max, int(T_max / dt))
    v_eval = sol.sol(t_eval)[0]

    spike_times = []
    for i in range(1, len(t_eval)):
        if v_eval[i-1] < threshold <= v_eval[i]:
            frac = (threshold - v_eval[i-1]) / (v_eval[i] - v_eval[i-1])
            spike_times.append(t_eval[i-1] + frac*(t_eval[i]-t_eval[i-1]))

    if len(spike_times) < n_orbits + 1:
        return None

    isis = np.diff(spike_times[-n_orbits-1:])
    return float(np.mean(isis))


def period_diverging(I_val: float, T_threshold: float) -> bool:
    """
    Returns True if the period at I_val exceeds T_threshold,
    i.e. the system is very close to the bifurcation.
    Uses a short T_max — if fewer than 3 spikes occur, period is diverging.
    """
    T_probe = T_threshold * 3
    def fhn(t, y):
        v, w = y
        return [v - v**3/3 - w + I_val, eps*(v + a - b*w)]

    y0 = [-1.0, -0.5]
    # Short warmup
    sol_w = solve_ivp(fhn, [0, T_probe/4], y0, method='RK45',
                      max_step=dt*20, rtol=1e-6, atol=1e-8)
    y0 = sol_w.y[:, -1]
    sol = solve_ivp(fhn, [0, T_probe], y0, method='RK45',
                    max_step=dt*5, dense_output=True, rtol=1e-6, atol=1e-8)

    t_eval = np.linspace(0, T_probe, int(T_probe/dt))
    v_eval = sol.sol(t_eval)[0]
    spikes = 0
    for i in range(1, len(t_eval)):
        if v_eval[i-1] < threshold <= v_eval[i]:
            spikes += 1
    # If fewer than 3 spikes in T_probe, period > T_threshold
    return spikes < 3

print("Step 1: Using analytical I_SNIC as reference point...")

ref    = FHN2D(I=-0.1, a=a, b=b)
I_snic_est = ref.I_snic
print(f"  Analytical I_SNIC = {I_snic_est:.5f}")

# Use the analytical value as I_c — the bifurcation IS where the
# fixed point hits the fold, regardless of what the bisection finds.
# The bisection finds where T drops below T_threshold which is a
# different (and less useful) quantity.
I_c = I_snic_est

# Verify: measure period at a few points just above I_SNIC to confirm
# the period is large there
print(f"\n  Verification — periods just above I_SNIC:")
for dI_test in [0.005, 0.010, 0.020, 0.030, 0.050]:
    T_test = measure_period(I_c + dI_test,
                            T_max=T_max, n_orbits=3)
    print(f"    ΔI={dI_test:.3f}  T={T_test:.1f}" if T_test
          else f"    ΔI={dI_test:.3f}  T=None (period too long)")

print(f"\n  Using I_c = {I_c:.5f} (analytical I_SNIC)")

# ---------------------------------------------------------------------------
# Step 2: Measure T(I) at log-spaced I values just above I_c
# ---------------------------------------------------------------------------
print(f"\nStep 2: Measuring T(I) at {n_I} I values above I_c...")

# Log-spaced ΔI from very close to I_SNIC — need to see the divergence
# Use large T_max to catch slow periods near the bifurcation
dI_vals = np.logspace(-3, np.log10(0.15), n_I)
I_vals  = I_c + dI_vals

T_vals  = []
dI_good = []

for I_val, dI in zip(I_vals, dI_vals):
    # Use longer T_max for points very close to bifurcation
    T_use = min(T_max * 3, max(T_max, 200.0 / dI))
    T = measure_period(I_val, T_max=T_use, n_orbits=max(3, n_orbits))
    if T is not None:
        T_vals.append(T)
        dI_good.append(dI)
        print(f"  I={I_val:.5f}  ΔI={dI:.2e}  T={T:.2f}")
    else:
        print(f"  I={I_val:.5f}  ΔI={dI:.2e}  T=None (too slow — "
              f"T_max={T_use:.0f})")

T_vals  = np.array(T_vals)
dI_good = np.array(dI_good)

if len(T_vals) < 4:
    print("\nInsufficient data points. Try --quick or increase T_max.")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Step 3: Fit both candidate laws
# ---------------------------------------------------------------------------
print(f"\nStep 3: Fitting divergence laws to {len(T_vals)} data points...")

# Law 1: SNIC — T = C / sqrt(ΔI)
def snic_law(dI, C):
    return C / np.sqrt(dI)

# Law 2: Homoclinic — T = C * log(1/ΔI) = -C * log(ΔI)
def homoclinic_law(dI, C):
    return C * np.log(1.0 / dI)

# Fit SNIC
try:
    popt_snic, pcov_snic = curve_fit(snic_law, dI_good, T_vals,
                                      p0=[1.0], maxfev=5000)
    C_snic = popt_snic[0]
    T_snic_fit = snic_law(dI_good, C_snic)
    ss_res = np.sum((T_vals - T_snic_fit)**2)
    ss_tot = np.sum((T_vals - np.mean(T_vals))**2)
    r2_snic = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
except Exception as e:
    print(f"  SNIC fit failed: {e}")
    C_snic, r2_snic = np.nan, np.nan

# Fit Homoclinic
try:
    popt_hom, pcov_hom = curve_fit(homoclinic_law, dI_good, T_vals,
                                    p0=[10.0], maxfev=5000)
    C_hom = popt_hom[0]
    T_hom_fit = homoclinic_law(dI_good, C_hom)
    ss_res = np.sum((T_vals - T_hom_fit)**2)
    r2_hom = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
except Exception as e:
    print(f"  Homoclinic fit failed: {e}")
    C_hom, r2_hom = np.nan, np.nan

# Power law fit (general): T = C * ΔI^α — what exponent does data prefer?
try:
    log_dI = np.log(dI_good)
    log_T  = np.log(T_vals)
    alpha_fit, log_C_fit = np.polyfit(log_dI, log_T, 1)
    C_power = np.exp(log_C_fit)
    T_power_fit = C_power * dI_good**alpha_fit
    ss_res = np.sum((T_vals - T_power_fit)**2)
    r2_power = 1 - ss_res/ss_tot if ss_tot > 0 else np.nan
except Exception as e:
    print(f"  Power law fit failed: {e}")
    alpha_fit, C_power, r2_power = np.nan, np.nan, np.nan

print(f"\n  Results:")
print(f"  SNIC law:       T = {C_snic:.3f} / sqrt(ΔI)      R² = {r2_snic:.5f}")
print(f"  Homoclinic law: T = {C_hom:.3f} * log(1/ΔI)    R² = {r2_hom:.5f}")
print(f"  Power law:      T = {C_power:.3f} * ΔI^{alpha_fit:.3f}    R² = {r2_power:.5f}")
print()

if not np.isnan(r2_snic) and not np.isnan(r2_hom):
    if r2_snic > r2_hom + 0.01:
        verdict = "SNIC  (saddle-node on invariant circle)"
        exponent_note = "σ_crit ~ √ε / √log(T_drift)  consistent with fold normal form"
    elif r2_hom > r2_snic + 0.01:
        verdict = "HOMOCLINIC  (homoclinic orbit)"
        exponent_note = "σ_crit scaling modified by log-diverging period"
    else:
        verdict = "AMBIGUOUS  — both laws fit similarly well"
        exponent_note = "Need finer resolution or wider ΔI range"
    print(f"  VERDICT: {verdict}")
    print(f"  Implication: {exponent_note}")

# ---------------------------------------------------------------------------
# Step 4: Plot
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

dI_plot = np.logspace(np.log10(dI_good.min()*0.5),
                       np.log10(dI_good.max()*1.5), 300)

# Panel 1: T vs ΔI (linear scale)
ax = axes[0]
ax.plot(dI_good, T_vals, 'o', color='#2166ac', ms=7, zorder=5,
        label='Measured T(I)')
if not np.isnan(r2_snic):
    ax.plot(dI_plot, snic_law(dI_plot, C_snic), '-',
            color='#d6604d', lw=2,
            label=f'SNIC: $T={C_snic:.2f}/\\sqrt{{\\Delta I}}$  ($R^2={r2_snic:.4f}$)')
if not np.isnan(r2_hom):
    ax.plot(dI_plot, homoclinic_law(dI_plot, C_hom), '--',
            color='#1a9641', lw=2,
            label=f'Homoclinic: $T={C_hom:.2f}\\log(1/\\Delta I)$  ($R^2={r2_hom:.4f}$)')
ax.set_xlabel(r'$\Delta I = I - I_c$', fontsize=11)
ax.set_ylabel('Period $T$', fontsize=11)
ax.set_title('Period divergence (linear)', fontsize=11)
ax.legend(fontsize=8)
ax.set_xlim(0, dI_good.max()*1.1)

# Panel 2: log T vs log ΔI (log-log — SNIC should be straight line with slope -1/2)
ax2 = axes[1]
ax2.loglog(dI_good, T_vals, 'o', color='#2166ac', ms=7, zorder=5,
           label='Measured')
if not np.isnan(r2_snic):
    ax2.loglog(dI_plot, snic_law(dI_plot, C_snic), '-',
               color='#d6604d', lw=2, label=f'SNIC (slope=−0.5)')
if not np.isnan(r2_power):
    ax2.loglog(dI_plot, C_power*dI_plot**alpha_fit, ':',
               color='#984ea3', lw=2,
               label=f'Power law: slope={alpha_fit:.3f}')
# Reference slope -0.5
ax2.loglog(dI_plot, dI_plot**(-0.5)*T_vals[0]*dI_good[0]**0.5,
           'k--', lw=1, alpha=0.4, label='slope=−0.5 (ref)')
ax2.set_xlabel(r'$\Delta I$  (log)', fontsize=11)
ax2.set_ylabel('Period $T$  (log)', fontsize=11)
ax2.set_title('Log-log: SNIC → slope = −½', fontsize=11)
ax2.legend(fontsize=8)

# Panel 3: T vs log(1/ΔI) (log-linear — homoclinic should be straight line)
ax3 = axes[2]
log_inv_dI = np.log(1.0/dI_good)
ax3.plot(log_inv_dI, T_vals, 'o', color='#2166ac', ms=7, zorder=5,
         label='Measured')
if not np.isnan(r2_hom):
    log_inv_dI_plot = np.log(1.0/dI_plot)
    ax3.plot(log_inv_dI_plot, homoclinic_law(dI_plot, C_hom), '--',
             color='#1a9641', lw=2,
             label=f'Homoclinic ($R^2={r2_hom:.4f}$)')
# Linear fit to check
m, c = np.polyfit(log_inv_dI, T_vals, 1)
ax3.plot(log_inv_dI, m*log_inv_dI + c, ':',
         color='#984ea3', lw=2,
         label=f'Linear fit: slope={m:.2f}')
ax3.set_xlabel(r'$\log(1/\Delta I)$', fontsize=11)
ax3.set_ylabel('Period $T$', fontsize=11)
ax3.set_title('Log-linear: homoclinic → straight line', fontsize=11)
ax3.legend(fontsize=8)

fig.suptitle(
    f'Period divergence test  |  $\\varepsilon={eps}$, $a={a}$, $b={b}$\n'
    f'$I_c={I_c:.5f}$  |  '
    f'SNIC $R^2={r2_snic:.4f}$  vs  Homoclinic $R^2={r2_hom:.4f}$  →  {verdict if "verdict" in dir() else "see fits"}',
    fontsize=11)
fig.tight_layout()
fig.savefig('figures/bifurcation/period_divergence.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'\nSaved: figures/bifurcation/period_divergence.png')

# Save data
np.savez('data/bifurcation/period_divergence.npz',
         I_c=I_c, dI=dI_good, T=T_vals,
         C_snic=C_snic, r2_snic=r2_snic,
         C_hom=C_hom,  r2_hom=r2_hom,
         alpha_fit=alpha_fit, C_power=C_power, r2_power=r2_power,
         eps=eps, a=a, b=b)
print('Saved: data/bifurcation/period_divergence.npz')

print(f'\n{"="*60}')
print(f'  SUMMARY')
print(f'  I_c = {I_c:.5f}  (estimated I_SNIC = {I_snic_est:.5f})')
print(f'  SNIC fit:       R² = {r2_snic:.5f}  (T ~ ΔI^-0.5)')
print(f'  Homoclinic fit: R² = {r2_hom:.5f}  (T ~ log(1/ΔI))')
print(f'  Power law fit:  R² = {r2_power:.5f}  (T ~ ΔI^{alpha_fit:.3f})')
if 'verdict' in dir():
    print(f'  Verdict: {verdict}')
print(f'{"="*60}\n')
