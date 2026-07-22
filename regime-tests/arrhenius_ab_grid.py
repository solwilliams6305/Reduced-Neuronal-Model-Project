#!/usr/bin/env python3
"""
arrhenius_ab_grid_checked.py

Arrhenius escape-barrier measurement for stochastic FitzHugh--Nagumo across
(a,b) parameter choices, with extra robustness checks.

Compared with arrhenius_ab_grid.py, this version adds:

  1. Multiple random seeds per (a,b,delta) case.
  2. Per-seed Arrhenius fits.
  3. Aggregated mean/SD/SE of the measured barrier B across seeds.
  4. Sigma-window sensitivity checks from the same rate data:
       - full fit
       - central fit, dropping lowest and highest sigma
       - drop_high_noise, dropping largest sigma
       - drop_low_noise, dropping smallest sigma
       - leave-one-sigma-out fits
  5. Explicit output of B-range and leave-one-out SD so you can see whether the
     fitted barrier is stable to the chosen sigma grid.

Model:
    dv = (v - v^3/3 - w + I) dt + sigma dW
    dw = eps (v + a - b w) dt

Default operating point is defined relative to the lower-branch Hopf current:
    I = I_Hopf(a,b,eps) - delta.

Example single-seed, multi-delta run:

  python arrhenius_ab_grid_checked.py \
    --a-values 0.6 0.7 0.8 \
    --b-values 0.7 0.8 0.9 \
    --deltas 0.015 0.03 0.05 \
    --eps 0.08 \
    --sigmas 0.03 0.035 0.04 0.045 0.05 0.06 0.07 0.08 \
    --N 600 --T 300 --dt 0.02 --n-workers 8 \
    --outdir results/arrhenius_ab_grid_checked

Example multi-seed robustness run:

  python arrhenius_ab_grid_checked.py \
    --a-values 0.6 0.7 0.8 \
    --b-values 0.7 0.8 0.9 \
    --deltas 0.03 \
    --seeds 1 2 3 4 5 \
    --sigmas 0.03 0.035 0.04 0.045 0.05 0.06 0.07 0.08 \
    --N 500 --T 300 --dt 0.02 --n-workers 8 \
    --outdir results/arrhenius_ab_grid_multiseed
"""

from __future__ import annotations

import argparse
import csv
import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Iterable

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Deterministic FHN utilities
# -----------------------------

@dataclass(frozen=True)
class FHNParams:
    I: float
    a: float
    b: float
    eps: float


def fixed_point_v_roots(I: float, a: float, b: float) -> np.ndarray:
    """Real v-roots of the FHN fixed-point equation."""
    coeff = np.array([-b, 0.0, 3.0 * b - 3.0, 3.0 * b * I - 3.0 * a], dtype=float)
    roots = np.roots(coeff)
    return np.sort(roots[np.abs(roots.imag) < 1e-9].real)


def jacobian(v: float, eps: float, b: float) -> np.ndarray:
    return np.array([[1.0 - v * v, -1.0], [eps, -eps * b]], dtype=float)


def select_stable_spiral_fp(p: FHNParams) -> Optional[dict]:
    """Choose the stable spiral fixed point with lowest max Re(lambda)."""
    infos = []
    for v in fixed_point_v_roots(p.I, p.a, p.b):
        w = (v + p.a) / p.b
        eig = np.linalg.eigvals(jacobian(v, p.eps, p.b))
        info = {
            "v": float(v),
            "w": float(w),
            "eig_real_max": float(np.max(eig.real)),
            "alpha": float(np.mean(eig.real)),
            "omega": float(np.max(np.abs(eig.imag))),
            "stable": bool(np.max(eig.real) < 0),
            "spiral": bool(np.max(np.abs(eig.imag)) > 1e-9),
        }
        infos.append(info)
    stable_spirals = [x for x in infos if x["stable"] and x["spiral"]]
    if not stable_spirals:
        return None
    return sorted(stable_spirals, key=lambda x: x["eig_real_max"])[0]


def hopf_current_lower_branch(a: float, b: float, eps: float) -> float:
    """Analytic lower-branch Hopf current."""
    arg = 1.0 - eps * b
    if arg <= 0.0 or b == 0.0:
        return float("nan")
    v = -math.sqrt(arg)
    w = (v + a) / b
    return float(w - v + v**3 / 3.0)


def local_spiral_basis(v_star: float, eps: float, b: float) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eig(jacobian(v_star, eps, b))
    idxs = np.where(vals.imag > 1e-9)[0]
    idx = int(idxs[0]) if len(idxs) else int(np.argmax(np.abs(vals.imag)))
    evec = vecs[:, idx]
    P = np.column_stack([evec.real, evec.imag]).astype(float)
    if abs(np.linalg.det(P)) < 1e-12:
        raise RuntimeError("Degenerate local spiral basis")
    return P, np.linalg.inv(P)


# -----------------------------
# Simulation and rate fitting
# -----------------------------

def simulate_first_passage(
    p: FHNParams,
    sigma: float,
    N: int,
    T: float,
    dt: float,
    r0: float,
    phi0: float,
    spike_threshold: float,
    seed: int,
) -> dict:
    """Vectorised Euler--Maruyama first-passage simulation."""
    fp = select_stable_spiral_fp(p)
    if fp is None:
        raise RuntimeError("No stable spiral fixed point")

    v_star, w_star = fp["v"], fp["w"]
    P, _ = local_spiral_basis(v_star, p.eps, p.b)
    offset = P @ np.array([r0 * math.cos(phi0), r0 * math.sin(phi0)])

    rng = np.random.default_rng(seed)
    steps = int(round(T / dt))
    t = np.linspace(0.0, steps * dt, steps + 1)
    sqrt_dt = math.sqrt(dt)

    v = np.full(N, v_star + offset[0], dtype=float)
    w = np.full(N, w_star + offset[1], dtype=float)
    spiked = np.zeros(N, dtype=bool)
    spike_time = np.full(N, np.nan, dtype=float)

    for step in range(1, steps + 1):
        active = ~spiked
        if not np.any(active):
            break
        va = v[active]
        wa = w[active]
        v[active] = va + (va - va**3 / 3.0 - wa + p.I) * dt + sigma * sqrt_dt * rng.standard_normal(np.sum(active))
        w[active] = wa + p.eps * (va + p.a - p.b * wa) * dt
        hit = active & (v >= spike_threshold)
        if np.any(hit):
            spiked[hit] = True
            spike_time[hit] = step * dt

    return {"t": t, "spiked": spiked, "spike_time": spike_time}


def survival_rate(spike_time: np.ndarray, spiked: np.ndarray, t: np.ndarray, S_hi: float, S_lo: float) -> Optional[dict]:
    """Fit S(t)=exp(-kt) over the window S in [S_lo,S_hi]."""
    N = spiked.size
    st = np.sort(spike_time[spiked])
    if st.size < 20:
        return None
    esc_by_t = np.searchsorted(st, t, side="right")
    S = 1.0 - esc_by_t / N
    mask = (S <= S_hi) & (S >= S_lo) & (t > 0)
    if mask.sum() < 8:
        return None
    tt = t[mask]
    lnS = np.log(S[mask])
    A = np.vstack([tt, np.ones_like(tt)]).T
    slope, intercept = np.linalg.lstsq(A, lnS, rcond=None)[0]
    k = -float(slope)
    if not np.isfinite(k) or k <= 0:
        return None
    pred = A @ np.array([slope, intercept])
    r2 = 1.0 - float(((lnS - pred) ** 2).sum()) / (float(((lnS - lnS.mean()) ** 2).sum()) + 1e-30)
    return {"k": k, "ln_k": math.log(k), "survival_r2": float(r2), "n_esc_fit_window_total": int(st.size)}


def frozen_w_barrier(p: FHNParams, w_star: float) -> tuple[float, tuple[float, ...] | tuple]:
    """Frozen-w 1D voltage barrier in the convention B = 2*dPhi."""
    coeffs = [1.0, 0.0, -3.0, 3.0 * (w_star - p.I)]
    roots = np.roots(coeffs)
    real = tuple(sorted(float(r.real) for r in roots if abs(r.imag) < 1e-7))
    if len(real) < 3:
        return float("nan"), real
    v_well, v_saddle = real[0], real[1]

    def Phi(v: float) -> float:
        return v**4 / 12.0 - v**2 / 2.0 - (p.I - w_star) * v

    return float(2.0 * (Phi(v_saddle) - Phi(v_well))), real


def linear_fit_xy(x: np.ndarray, y: np.ndarray) -> Optional[dict]:
    if len(x) < 3:
        return None
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    pred = A @ np.array([slope, intercept])
    r2 = 1.0 - float(((y - pred) ** 2).sum()) / (float(((y - y.mean()) ** 2).sum()) + 1e-30)
    return {
        "B_meas": -float(slope),
        "log_A": float(intercept),
        "A_prefactor": float(math.exp(intercept)),
        "arrhenius_r2": float(r2),
        "n_rate_points": int(len(x)),
    }


def arrhenius_fit(rate_rows: list[dict], sigma_min: float | None = None, sigma_max: float | None = None) -> Optional[dict]:
    valid = []
    for r in rate_rows:
        if r.get("rate_status") != "ok":
            continue
        s = float(r["sigma"])
        if sigma_min is not None and s < sigma_min:
            continue
        if sigma_max is not None and s > sigma_max:
            continue
        if np.isfinite(r.get("k", np.nan)) and r["k"] > 0:
            valid.append(r)
    if len(valid) < 3:
        return None
    x = np.array([r["inv_sig2"] for r in valid], dtype=float)
    y = np.array([r["ln_k"] for r in valid], dtype=float)
    out = linear_fit_xy(x, y)
    if out is not None:
        out["sigmas_used"] = ";".join(f"{r['sigma']:.8g}" for r in sorted(valid, key=lambda z: z["sigma"]))
    return out


def fit_sensitivity(rate_rows: list[dict], sigma_min: float | None = None, sigma_max: float | None = None) -> list[dict]:
    """Return robustness variants for the Arrhenius slope from one case/seed."""
    valid = []
    for r in rate_rows:
        if r.get("rate_status") != "ok":
            continue
        s = float(r["sigma"])
        if sigma_min is not None and s < sigma_min:
            continue
        if sigma_max is not None and s > sigma_max:
            continue
        if np.isfinite(r.get("k", np.nan)) and r["k"] > 0:
            valid.append(r)
    valid = sorted(valid, key=lambda z: z["sigma"])
    if len(valid) < 3:
        return []

    variants: list[tuple[str, list[dict], str]] = [("full", valid, "all usable sigmas")]
    if len(valid) >= 5:
        variants.append(("central_drop_extremes", valid[1:-1], "drop lowest and highest sigma"))
    if len(valid) >= 4:
        variants.append(("drop_high_noise", valid[:-1], "drop largest sigma"))
        variants.append(("drop_low_noise", valid[1:], "drop smallest sigma"))
    if len(valid) >= 4:
        for j, rdrop in enumerate(valid):
            variants.append((f"leave_one_out_sigma_{rdrop['sigma']:.8g}", valid[:j] + valid[j + 1:], f"drop sigma={rdrop['sigma']:.8g}"))

    rows = []
    for name, rows_used, desc in variants:
        fit = arrhenius_fit(rows_used, sigma_min=None, sigma_max=None)
        if fit is None:
            continue
        base = {k: rows_used[0][k] for k in ["a", "b", "eps", "delta", "I_hopf", "I", "seed"] if k in rows_used[0]}
        base.update({"fit_variant": name, "fit_description": desc})
        base.update(fit)
        rows.append(base)
    return rows


# -----------------------------
# One parameter case / seed
# -----------------------------

def run_case(case: dict) -> tuple[list[dict], dict, list[dict]]:
    a = case["a"]
    b = case["b"]
    eps = case["eps"]
    delta = case["delta"]
    seed0 = case["seed"]
    fixed_I = case.get("fixed_I")
    I_hopf = hopf_current_lower_branch(a, b, eps)
    I = float(fixed_I) if fixed_I is not None else float(I_hopf - delta)
    p = FHNParams(I=I, a=a, b=b, eps=eps)

    fp = select_stable_spiral_fp(p)
    summary = {
        "case_id": case["case_id"], "param_case_id": case["param_case_id"],
        "a": a, "b": b, "eps": eps, "delta": delta, "seed": seed0,
        "I_hopf": I_hopf, "I": I,
        "status": "ok", "message": "",
    }
    if fp is None:
        summary.update({"status": "skipped", "message": "No stable spiral fixed point"})
        return [], summary, []

    summary.update({
        "v_fp": fp["v"], "w_fp": fp["w"],
        "alpha": fp["alpha"], "omega": fp["omega"],
        "kappa": abs(fp["alpha"]) / fp["omega"] if fp["omega"] > 0 else float("nan"),
    })
    B_frozen, roots = frozen_w_barrier(p, fp["w"])
    summary["B_frozen_w"] = B_frozen
    summary["frozen_roots"] = ";".join(f"{x:.8g}" for x in roots)

    rate_rows = []
    for idx, sigma in enumerate(case["sigmas"]):
        # deterministic but unique RNG stream for parameter case, seed, sigma
        sim_seed = int(seed0 + 100000 * case["param_case_id"] + 1000 * idx)
        row = {
            "case_id": case["case_id"], "param_case_id": case["param_case_id"],
            "a": a, "b": b, "eps": eps, "delta": delta, "seed": seed0,
            "I_hopf": I_hopf, "I": I, "sigma": sigma, "inv_sig2": 1.0 / sigma**2,
            "N": case["N"], "T": case["T"], "dt": case["dt"], "sim_seed": sim_seed,
        }
        try:
            sim = simulate_first_passage(
                p=p, sigma=sigma, N=case["N"], T=case["T"], dt=case["dt"],
                r0=case["r0"], phi0=case["phi0"],
                spike_threshold=case["spike_threshold"], seed=sim_seed,
            )
            res = survival_rate(sim["spike_time"], sim["spiked"], sim["t"], case["S_hi"], case["S_lo"])
            row.update({
                "n_esc": int(sim["spiked"].sum()),
                "escape_fraction": float(sim["spiked"].mean()),
            })
            if res is None:
                row.update({"k": float("nan"), "ln_k": float("nan"), "survival_r2": float("nan"), "rate_status": "bad_window"})
            else:
                row.update(res)
                row["rate_status"] = "ok"
        except Exception as exc:
            row.update({
                "n_esc": 0, "escape_fraction": float("nan"), "k": float("nan"),
                "ln_k": float("nan"), "survival_r2": float("nan"),
                "rate_status": f"error: {type(exc).__name__}: {exc}",
            })
        rate_rows.append(row)

    fit = arrhenius_fit(rate_rows, case.get("fit_sigma_min"), case.get("fit_sigma_max"))
    fit_check_rows = fit_sensitivity(rate_rows, case.get("fit_sigma_min"), case.get("fit_sigma_max"))

    if fit is None:
        summary.update({"status": "no_fit", "message": "Fewer than three usable rate points"})
    else:
        summary.update(fit)
        summary["B_over_frozen_w"] = (fit["B_meas"] / B_frozen) if np.isfinite(B_frozen) and B_frozen != 0 else float("nan")
        if fit["arrhenius_r2"] < case["min_arrhenius_r2_warning"]:
            summary["message"] = "Arrhenius fit weak over chosen sigma range"

        # Pull useful robustness diagnostics into the main summary.
        B_all = [r["B_meas"] for r in fit_check_rows if np.isfinite(r.get("B_meas", np.nan))]
        B_loo = [r["B_meas"] for r in fit_check_rows if r.get("fit_variant", "").startswith("leave_one_out") and np.isfinite(r.get("B_meas", np.nan))]
        if B_all:
            summary["B_fit_variant_min"] = float(np.min(B_all))
            summary["B_fit_variant_max"] = float(np.max(B_all))
            summary["B_fit_variant_range"] = float(np.max(B_all) - np.min(B_all))
            summary["B_fit_variant_rel_range"] = float((np.max(B_all) - np.min(B_all)) / abs(summary["B_meas"])) if summary["B_meas"] != 0 else float("nan")
        if B_loo:
            summary["B_leave_one_out_sd"] = float(np.std(B_loo, ddof=1)) if len(B_loo) > 1 else 0.0
            summary["B_leave_one_out_rel_sd"] = float(summary["B_leave_one_out_sd"] / abs(summary["B_meas"])) if summary["B_meas"] != 0 else float("nan")

    return rate_rows, summary, fit_check_rows


# -----------------------------
# Aggregation, CSV, plotting
# -----------------------------

def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys, seen = [], set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                keys.append(k)
                seen.add(k)
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=keys)
        wr.writeheader()
        wr.writerows(rows)


def mean_sd_se(xs: Iterable[float]) -> tuple[float, float, float, int]:
    arr = np.array([x for x in xs if np.isfinite(x)], dtype=float)
    n = int(arr.size)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    mean = float(arr.mean())
    sd = float(arr.std(ddof=1)) if n > 1 else 0.0
    se = float(sd / math.sqrt(n)) if n > 1 else 0.0
    return mean, sd, se, n


def aggregate_summaries(summary_rows: list[dict]) -> list[dict]:
    good = [r for r in summary_rows if r.get("status") == "ok" and np.isfinite(r.get("B_meas", np.nan))]
    groups: dict[tuple, list[dict]] = {}
    for r in good:
        key = (r["a"], r["b"], r["eps"], r["delta"], r["I_hopf"], r["I"])
        groups.setdefault(key, []).append(r)

    out = []
    for key, rows in sorted(groups.items(), key=lambda z: (z[0][3], z[0][0], z[0][1])):
        a, b, eps, delta, I_hopf, I = key
        B_mean, B_sd, B_se, n = mean_sd_se(r["B_meas"] for r in rows)
        R2_mean, R2_sd, R2_se, _ = mean_sd_se(r.get("arrhenius_r2", float("nan")) for r in rows)
        rel_win_mean, _, _, _ = mean_sd_se(r.get("B_fit_variant_rel_range", float("nan")) for r in rows)
        loo_rel_sd_mean, _, _, _ = mean_sd_se(r.get("B_leave_one_out_rel_sd", float("nan")) for r in rows)
        frozen = rows[0].get("B_frozen_w", float("nan"))
        out.append({
            "a": a, "b": b, "eps": eps, "delta": delta, "I_hopf": I_hopf, "I": I,
            "n_seeds_ok": n,
            "B_mean": B_mean, "B_sd_across_seeds": B_sd, "B_se_across_seeds": B_se,
            "B_cv_across_seeds": B_sd / abs(B_mean) if np.isfinite(B_mean) and B_mean != 0 else float("nan"),
            "arrhenius_r2_mean": R2_mean, "arrhenius_r2_sd": R2_sd,
            "B_frozen_w": frozen,
            "B_mean_over_frozen_w": B_mean / frozen if np.isfinite(frozen) and frozen != 0 else float("nan"),
            "mean_B_fit_variant_rel_range": rel_win_mean,
            "mean_B_leave_one_out_rel_sd": loo_rel_sd_mean,
            "seeds": ";".join(str(int(r["seed"])) for r in sorted(rows, key=lambda x: x["seed"])),
        })
    return out


def make_plots(outdir: Path, rate_rows: list[dict], summary_rows: list[dict], aggregate_rows: list[dict]) -> None:
    good_summary = [r for r in summary_rows if r.get("status") == "ok" and np.isfinite(r.get("B_meas", np.nan))]
    good_agg = [r for r in aggregate_rows if np.isfinite(r.get("B_mean", np.nan))]
    if not good_summary:
        return

    # Barrier over (a,b), using aggregate if available, otherwise per-seed summary.
    plot_rows = good_agg if good_agg else good_summary
    value_key = "B_mean" if good_agg else "B_meas"
    r2_key = "arrhenius_r2_mean" if good_agg else "arrhenius_r2"
    for delta in sorted(set(round(float(r["delta"]), 12) for r in plot_rows)):
        rows = [r for r in plot_rows if abs(float(r["delta"]) - delta) < 1e-12]
        plt.figure(figsize=(7, 5))
        sc = plt.scatter([r["a"] for r in rows], [r["b"] for r in rows], c=[r[value_key] for r in rows], s=120)
        for r in rows:
            label = f"{r[value_key]:.2g}\nR2={r[r2_key]:.2f}"
            if "B_se_across_seeds" in r and r.get("n_seeds_ok", 0) > 1:
                label += f"\nSE={r['B_se_across_seeds']:.1g}"
            plt.text(r["a"], r["b"], label, ha="center", va="center", fontsize=8)
        plt.colorbar(sc, label="measured Arrhenius barrier B")
        plt.xlabel("a")
        plt.ylabel("b")
        plt.title(f"Measured barrier over (a,b), delta={delta:g}")
        plt.tight_layout()
        plt.savefig(outdir / f"barrier_ab_delta_{delta:g}.png", dpi=180)
        plt.close()

    # Frozen-w mismatch scatter.
    rows_fw = plot_rows
    plt.figure(figsize=(7, 5))
    x = [r["B_frozen_w"] for r in rows_fw]
    y = [r[value_key] for r in rows_fw]
    plt.loglog(x, y, "o")
    finite = [z for z in x + y if np.isfinite(z) and z > 0]
    if finite:
        lo, hi = min(finite), max(finite)
        plt.loglog([lo, hi], [lo, hi], "--", linewidth=1)
    plt.xlabel("frozen-w barrier")
    plt.ylabel("measured Arrhenius barrier")
    plt.title("2D measured barrier vs frozen-w 1D prediction")
    plt.tight_layout()
    plt.savefig(outdir / "measured_vs_frozen_w_barrier.png", dpi=180)
    plt.close()

    # Arrhenius lines; if many seeds, plot only seed-smallest per param case to avoid spaghetti overdose.
    ok_rates = [r for r in rate_rows if r.get("rate_status") == "ok" and np.isfinite(r.get("ln_k", np.nan))]
    if ok_rates:
        min_seed_by_param = {}
        for r in ok_rates:
            key = r["param_case_id"]
            min_seed_by_param[key] = min(min_seed_by_param.get(key, r["seed"]), r["seed"])
        plt.figure(figsize=(8, 5.5))
        for pcid in sorted(set(r["param_case_id"] for r in ok_rates)):
            rs = [r for r in ok_rates if r["param_case_id"] == pcid and r["seed"] == min_seed_by_param[pcid]]
            rs = sorted(rs, key=lambda z: z["inv_sig2"])
            if len(rs) < 2:
                continue
            label = f"a={rs[0]['a']:g}, b={rs[0]['b']:g}, d={rs[0]['delta']:g}"
            plt.plot([r["inv_sig2"] for r in rs], [r["ln_k"] for r in rs], marker="o", linewidth=1, label=label)
        plt.xlabel(r"$1/\sigma^2$")
        plt.ylabel("ln escape rate")
        plt.title("Arrhenius fits across parameter cases")
        plt.legend(fontsize=7, ncol=2)
        plt.tight_layout()
        plt.savefig(outdir / "arrhenius_all_cases.png", dpi=180)
        plt.close()

    # Scaling plot B vs delta, one line per (a,b), if multiple deltas exist.
    if len(set(round(float(r["delta"]), 12) for r in plot_rows)) >= 2:
        plt.figure(figsize=(7, 5))
        for (a, b), rows in sorted({(r["a"], r["b"]): [x for x in plot_rows if x["a"] == r["a"] and x["b"] == r["b"]] for r in plot_rows}.items()):
            rows = sorted(rows, key=lambda z: z["delta"])
            if len(rows) < 2:
                continue
            plt.loglog([r["delta"] for r in rows], [r[value_key] for r in rows], marker="o", label=f"a={a:g}, b={b:g}")
        plt.xlabel(r"Hopf deficit $\delta=I_H-I$")
        plt.ylabel("measured Arrhenius barrier B")
        plt.title("Barrier scaling with Hopf deficit")
        plt.legend(fontsize=7, ncol=2)
        plt.tight_layout()
        plt.savefig(outdir / "barrier_vs_delta_loglog.png", dpi=180)
        plt.close()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a-values", type=float, nargs="*", default=[0.6, 0.7, 0.8])
    ap.add_argument("--b-values", type=float, nargs="*", default=[0.7, 0.8, 0.9])
    ap.add_argument("--deltas", type=float, nargs="*", default=[0.03])
    ap.add_argument("--eps", type=float, default=0.08)
    ap.add_argument("--mode", choices=["hopf_delta", "fixed_I"], default="hopf_delta")
    ap.add_argument("--fixed-I", type=float, default=0.30)
    ap.add_argument("--sigmas", type=float, nargs="*", default=[0.03, 0.035, 0.04, 0.045, 0.05, 0.06, 0.07, 0.08])
    ap.add_argument("--fit-sigma-min", type=float, default=None, help="Optional lower sigma bound for primary Arrhenius fit.")
    ap.add_argument("--fit-sigma-max", type=float, default=None, help="Optional upper sigma bound for primary Arrhenius fit.")
    ap.add_argument("--seeds", type=int, nargs="*", default=[1], help="One or more random seeds per parameter case.")
    ap.add_argument("--N", type=int, default=400)
    ap.add_argument("--T", type=float, default=250.0)
    ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--r0", type=float, default=0.06)
    ap.add_argument("--phi0", type=float, default=0.0)
    ap.add_argument("--spike-threshold", type=float, default=1.0)
    ap.add_argument("--S-hi", type=float, default=0.85)
    ap.add_argument("--S-lo", type=float, default=0.10)
    ap.add_argument("--n-workers", type=int, default=1)
    ap.add_argument("--min-arrhenius-r2-warning", type=float, default=0.95)
    ap.add_argument("--outdir", type=str, default="results/arrhenius_ab_grid_checked")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cases = []
    case_id = 0
    param_case_id = 0
    for a in args.a_values:
        for b in args.b_values:
            for delta in args.deltas:
                for seed in args.seeds:
                    cases.append({
                        "case_id": case_id,
                        "param_case_id": param_case_id,
                        "a": float(a), "b": float(b), "eps": float(args.eps), "delta": float(delta),
                        "seed": int(seed),
                        "fixed_I": float(args.fixed_I) if args.mode == "fixed_I" else None,
                        "sigmas": list(map(float, args.sigmas)),
                        "fit_sigma_min": args.fit_sigma_min,
                        "fit_sigma_max": args.fit_sigma_max,
                        "N": args.N, "T": args.T, "dt": args.dt,
                        "r0": args.r0, "phi0": args.phi0, "spike_threshold": args.spike_threshold,
                        "S_hi": args.S_hi, "S_lo": args.S_lo,
                        "min_arrhenius_r2_warning": args.min_arrhenius_r2_warning,
                    })
                    case_id += 1
                param_case_id += 1

    print(f"Running {len(cases)} case/seed jobs")
    print(f"Parameter cases: {param_case_id}; seeds per parameter case: {len(args.seeds)}")
    print(f"Sigmas: {args.sigmas}")
    print(f"Output directory: {outdir}")

    all_rate_rows: list[dict] = []
    summary_rows: list[dict] = []
    fit_check_rows: list[dict] = []

    if args.n_workers == 1:
        iterator = []
        for case in cases:
            iterator.append((case, run_case(case)))
        for case, result in iterator:
            rate_rows, summary, checks = result
            all_rate_rows.extend(rate_rows)
            summary_rows.append(summary)
            fit_check_rows.extend(checks)
            print(f"case {case['case_id']:03d}: a={summary['a']:.3g} b={summary['b']:.3g} "
                  f"d={summary['delta']:.3g} seed={summary['seed']} status={summary['status']} "
                  f"B={summary.get('B_meas', float('nan')):.4g} R2={summary.get('arrhenius_r2', float('nan')):.3f} "
                  f"winRel={summary.get('B_fit_variant_rel_range', float('nan')):.3g}")
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futs = {ex.submit(run_case, case): case for case in cases}
            for fut in as_completed(futs):
                case = futs[fut]
                rate_rows, summary, checks = fut.result()
                all_rate_rows.extend(rate_rows)
                summary_rows.append(summary)
                fit_check_rows.extend(checks)
                print(f"case {case['case_id']:03d}: a={summary['a']:.3g} b={summary['b']:.3g} "
                      f"d={summary['delta']:.3g} seed={summary['seed']} status={summary['status']} "
                      f"B={summary.get('B_meas', float('nan')):.4g} R2={summary.get('arrhenius_r2', float('nan')):.3f} "
                      f"winRel={summary.get('B_fit_variant_rel_range', float('nan')):.3g}")

    summary_rows = sorted(summary_rows, key=lambda r: (r.get("delta", 0), r.get("a", 0), r.get("b", 0), r.get("seed", 0)))
    all_rate_rows = sorted(all_rate_rows, key=lambda r: (r.get("param_case_id", 0), r.get("seed", 0), r.get("sigma", 0)))
    fit_check_rows = sorted(fit_check_rows, key=lambda r: (r.get("delta", 0), r.get("a", 0), r.get("b", 0), r.get("seed", 0), r.get("fit_variant", "")))
    aggregate_rows = aggregate_summaries(summary_rows)

    write_csv(outdir / "arrhenius_ab_rates.csv", all_rate_rows)
    write_csv(outdir / "arrhenius_ab_summary_by_seed.csv", summary_rows)
    write_csv(outdir / "arrhenius_ab_fit_checks.csv", fit_check_rows)
    write_csv(outdir / "arrhenius_ab_summary_aggregated.csv", aggregate_rows)
    make_plots(outdir, all_rate_rows, summary_rows, aggregate_rows)

    n_ok = sum(1 for r in summary_rows if r.get("status") == "ok")
    n_weak = sum(1 for r in summary_rows if r.get("status") == "ok" and r.get("arrhenius_r2", 0) < args.min_arrhenius_r2_warning)
    print("\n================ DONE ================")
    print(f"usable Arrhenius fits: {n_ok}/{len(summary_rows)}")
    print(f"weak fits below R2<{args.min_arrhenius_r2_warning}: {n_weak}")
    print(f"rates:       {outdir / 'arrhenius_ab_rates.csv'}")
    print(f"by-seed:     {outdir / 'arrhenius_ab_summary_by_seed.csv'}")
    print(f"fit checks:  {outdir / 'arrhenius_ab_fit_checks.csv'}")
    print(f"aggregated:  {outdir / 'arrhenius_ab_summary_aggregated.csv'}")
    print("plots:       barrier_ab_delta_*.png, measured_vs_frozen_w_barrier.png, arrhenius_all_cases.png, barrier_vs_delta_loglog.png")


if __name__ == "__main__":
    main()
