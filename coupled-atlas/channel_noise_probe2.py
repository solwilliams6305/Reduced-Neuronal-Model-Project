import numpy as np, sys
sys.path.insert(0,"/Users/solomonwilliams/Reduced Neuronal Model Project/coupled-atlas")
from channel_noise_probe import solve_fp_eta
s2=np.sqrt(2.0)
print("CORRECTED localization: is the sensitive region the TURNING (Y~0) or the LANDING (Y~-1.6)?")
cases={
 "A const sqrt2 (beta2)":         lambda Y: s2,
 "turning x1.5 in [-0.5,0.5]":    lambda Y: 1.5*s2 if -0.5<=Y<=0.5 else s2,
 "landing x1.5 in [-2.2,-1.2]":   lambda Y: 1.5*s2 if -2.2<=Y<=-1.2 else s2,
 "confining-tail x1.5 in [1,3]":  lambda Y: 1.5*s2 if Y>=1 else s2,
}
print(f"  {'case':30} {'mean':>7} {'std':>6} {'skew':>7} {'exk':>7}")
ref=None
for name,fn in cases.items():
    m,sd,sk,ek=solve_fp_eta(fn)
    if ref is None: ref=(m,sd,sk,ek)
    dsk=sk-ref[2]
    print(f"  {name:30} {m:+7.3f} {sd:6.3f} {sk:+7.3f} {ek:+7.3f}   dskew_vsA={dsk:+.3f}")

print("\nREALISTIC channel noise: eta(Y) varies only MILDLY over the turning (real g(v) ~ frozen to O(sqrt eps)):")
print("  (does a small ~10-15% variation stay ~ beta=2, confirming robustness for realistic channel noise?)")
for amp,desc in [(0.0,"flat (beta2)"),(0.10,"+/-10% ramp"),(0.15,"+/-15% ramp")]:
    fn=lambda Y,a=amp: s2*(1+a*np.tanh(Y))   # mild, ~+/-a, transitions across the turning
    m,sd,sk,ek=solve_fp_eta(fn)
    print(f"  {desc:16}: mean={m:+.3f} std={sd:.3f} skew={sk:+.3f} exk={ek:+.3f}")
print("  [beta=2 ref: mean -1.594 std 0.689 skew +0.601 exk -0.244]")
