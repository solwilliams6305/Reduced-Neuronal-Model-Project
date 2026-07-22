import numpy as np, sympy as sp
exec(open('_diag2.py').read().split('for zs in')[0])   # reuse build_Y, cov
for MT in [3,4]:
    Y=build_Y(n=30,MT=MT,zero_s=True)
    C24=cov(Y[2],Y[4],MT); m2=bexpect(Y[4])
    print(f"  MT={MT} n=30: 2Cov24_reg={2*C24:+.4f}  m2=<Y4>={m2:+.4f} (mollified m2~+0.078)")
