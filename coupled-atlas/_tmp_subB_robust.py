"""Robustness of the central cell (zr=1.9, ar=0.0) over v6 and v3 variants."""
import numpy as np
from _tmp_subB_borel import profile, fmt_row

ZR, AR = 1.9, 0.0
base = [0.134, 0.111, 0.104, -0.030, -0.451, -1.19, -1.90]

def make(v6=None, v3=None):
    v = list(base)
    if v6 is not None: v[6] = v6
    if v3 is not None: v[3] = v3
    return v

print("CENTRAL CELL zr=1.9 ar=0.0 -- robustness")
print("="*90)
variants = [
    ("nominal v6=-1.90 v3=-0.030", make()),
    ("v6=-1.75", make(v6=-1.75)),
    ("v6=-2.05", make(v6=-2.05)),
    ("v3=-0.020", make(v3=-0.020)),
    ("v3=-0.020,v6=-1.75", make(v6=-1.75, v3=-0.020)),
    ("v3=-0.020,v6=-2.05", make(v6=-2.05, v3=-0.020)),
]
for label, v in variants:
    af = profile(v, ZR, AR, alpha_fixed=None, use_real=True)
    a0 = profile(v, ZR, AR, alpha_fixed=0.0, use_real=True)
    print(fmt_row(f"[af] {label}", af))
    print(fmt_row(f"[a0] {label}", a0))
    print()
