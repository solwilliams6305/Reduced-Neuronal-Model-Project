"""
Pure-Python brentq shim for environments without scipy.
Idempotent: if scipy.optimize is importable, do nothing.
"""
from __future__ import annotations

import sys
import types


def _install_brentq_shim() -> None:
    try:
        import scipy.optimize  # noqa: F401
        return
    except Exception:
        pass

    def brentq(f, a, b, **kwargs):
        fa, fb = f(a), f(b)
        if fa == 0:
            return a
        if fb == 0:
            return b
        if fa * fb > 0:
            raise ValueError(
                f"brentq shim: f(a)={fa}, f(b)={fb} have the same sign on [{a}, {b}]"
            )
        for _ in range(120):
            m = 0.5 * (a + b)
            fm = f(m)
            if fm == 0 or (b - a) < 1e-13:
                return m
            if fa * fm < 0:
                b, fb = m, fm
            else:
                a, fa = m, fm
        return 0.5 * (a + b)

    scipy_mod = sys.modules.get("scipy") or types.ModuleType("scipy")
    opt_mod = sys.modules.get("scipy.optimize") or types.ModuleType("scipy.optimize")
    opt_mod.brentq = brentq
    scipy_mod.optimize = opt_mod
    sys.modules["scipy"] = scipy_mod
    sys.modules["scipy.optimize"] = opt_mod

    # Minimal scipy.integrate shim for hopf_analysis.py: not used here, skip.


_install_brentq_shim()
