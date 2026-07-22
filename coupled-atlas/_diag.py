import numpy as np, sympy as sp
from math import comb
import chaos_boundary as CB
from chaos_boundary import B, bmul, bexpect, setup, GEO

def Vd(Y,d): return {0:-Y*Y,1:-2*Y,2:-2.0}.get(d,0.0)

setup(n_grid=40, MAXORD=3)
u0ps=GEO['u0ps']; Yst=GEO['s'][0]; c=GEO['c']; sw=GEO['sw']; U=GEO['U']; PSI=GEO['PSI']
print(f"u0(node)=U[0]={U[0]:.3e}  psi(node)=PSI[0]={PSI[0]:.4f} (c={c:.4f})")

def elt(ch): return B({(k,()):v for k,v in ch.t.items()})
s0=B({(0,(0,)):np.array(1.0)})
# TEST 1: <s0 * I1(psi)> should be 1/2 psi(node) = c/2
Ipsi=B({(1,()):PSI*sw})     # I1(psi) absorbed
t1=bexpect(bmul(s0,Ipsi,4))
print(f"TEST1 <s0*I1(psi)> = {t1:.5f}   expect 1/2 psi(node)= {0.5*PSI[0]:.5f}")
# TEST 2: <s0 * I1(A1)>, A1=alpha u0^2 vanishes at node -> 0
alpha=c/u0ps; A1=alpha*U**2
IA1=B({(1,()):A1*sw})
print(f"TEST2 <s0*I1(A1)> = {bexpect(bmul(s0,IA1,4)):.6f}   expect ~0 (A1(node)=0)")
# TEST 3: <s0*s0> dropped -> 0
print(f"TEST3 <s0*s0> = {bexpect(bmul(s0,s0,4)):.6f}   expect 0 (renorm drop)")
# TEST 4: <s1 * I1(psi)> = 1/2 psi'(node)
s1=B({(0,(1,)):np.array(1.0)})
h=GEO['s'][1]-GEO['s'][0]; psip_node=(PSI[1]-PSI[0])/h
print(f"TEST4 <s1*I1(psi)> = {bexpect(bmul(s1,Ipsi,4)):.5f}   expect 1/2 psi'(node)= {0.5*psip_node:.5f}")
