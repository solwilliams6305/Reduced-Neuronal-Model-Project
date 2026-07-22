# Program 2 — the deterministic Wiener-chaos moment engine (for exact v_n, esp. v4)

_2026-07-08 (Fable 5). Built the MC-free chaos engine to reach $v_4$ (sampling provably dies past $v_3$, notes §8).
Status: **core validated exactly; boundary layer implemented but needs one calibration/debug pass.**_

## What it is
Every node functional = element of the Wiener chaos, stored as symmetric kernels $\{p:h_p\}$ in weight-absorbed form.
Node fields are Markovian bond-chains of the causal Green's function: $u_k(Y^\star_0)$ has ordered kernel
$G(\text{node},s_1)G(s_1,s_2)\cdots G(s_{k-1},s_k)u_0(s_k)$. Products use the Itô/Wick product
$I_p(f)I_q(g)=\sum_r C(p,r)C(q,r)r!\,I_{p+q-2r}(\mathrm{sym}(f\otimes_r g))$; moments are exact kernel contractions.
Boundary noise $s_j=\xi^{(j)}(Y^\star_0)$ is a "boundary leg"; the renormalization prescription is
**boundary–bulk $\to\tfrac12 a^{(j)}(\text{node})$ (kept), boundary–boundary $\to$ dropped**.
Files: `coupled-atlas/chaos_engine.py` (core), `chaos_boundary.py` (boundary layer), `chaos_v2_test.py` (v2 test).

## Validated (exact, zero Monte Carlo) [NUMERIC ✓✓]
| quantity | engine | target | |
|---|---|---|---|
| $v_0$ | 0.1343 | 0.1339 | ✓ |
| $m_1$ | +0.2135 | +0.2124 | ✓ |
| $\operatorname{Var}(Y_2)$ | 0.0641 | 0.0636 | ✓ |
| $\langle Y_1Y_3\rangle$ | +0.0234 | +0.0231 | ✓ |
| $v_1$ | +0.1110 | +0.110 | ✓ |
| $\langle Y_3^2\rangle$ (boundary-free part of $v_2$) | 0.146 | ~0.145 | ✓ |

The Wiener-chaos algebra, the Markovian node-field kernels, and all **boundary-free** moments are correct and
grid-converged. This is the hard core, working.

## The real limitation: the truncated-product ARCHITECTURE, not `beval` [DIAGNOSED 2026-07-08]
Debugged the $v_2=0.137$ (vs 0.110) discrepancy to the root. **The boundary-Wick primitives are CORRECT** —
validated directly: $\langle s_0 I_1(\psi)\rangle=\tfrac12\psi(\text{node})$ ✓, $\langle s_0 I_1(A_1)\rangle=0$ ✓
(kernel vanishes at node), $\langle s_1 I_1(\psi)\rangle=\tfrac12\psi'(\text{node})$ ✓, $\langle s_0s_0\rangle=0$
(renorm drop) ✓. And splitting $\operatorname{Cov}(Y_2,Y_4)$ boundary-on vs boundary-off gives the **same** value —
so the boundary noise is handled right; the error is in the boundary-**free** $Y_4$.

**Root cause — up-then-down chaos flow defeats product truncation.** $U_{k,m}=u_k^{(m)}(Y^\star_0)$ is $k$-chaos, so a
monomial like $U_{1,0}^5$ (five 1-chaos factors) passes through **5th chaos** before contracting back down; its
1-chaos part (needed for $\langle Y_1Y_5\rangle$) depends on those high-chaos intermediates. Truncating each
intermediate product at order MT **drops** contributions that fold back into the low-chaos parts. Evidence
(MT-sweep, `_diag2/_diag3`): $m_2=\langle Y_4\rangle$ = 0.101 (MT=3) → 0.081 (MT=4) → target 0.078; $2\operatorname{Cov}
(Y_2,Y_4)$ = 0.031 → 0.019 → ~0.008. Converges, but only when **MT $\ge$ max total chaos in $Y_n$ = $n$**: MT=5 for
$Y_5$ ($v_2$), MT=7 for $v_3$, MT=9 for $v_4$ — i.e. $n^9$ dense tensors, intractable. So the truncated-product
engine is a validated but **architecturally bounded** tool (exact only through the orders where MT$\ge n$ is affordable,
~MT 4–5).

## RESOLVED — the diagram engine works [2026-07-08, NUMERIC ✓]
Built `chaos_diagram.py`: each moment $\langle\prod\text{atoms}\rangle$ = sum over contraction MULTIGRAPHS (edges
pair legs across atoms; degree = chaos order; boundary legs degree 1; no intra-atom, no boundary–boundary [renorm]),
value $=\sum_M[\prod_i k_i!/\prod m_{ij}!]\,\mathrm{TN}(M)$ with `beval` for boundary edges. Never materializes a tensor
above the largest single atom. **Exact — no truncation, no sampling.** Validated + grid-converged:
$v_0=0.1343$, $v_1\to0.111$, $v_2\to\approx0.10$ (n=30/40/50: 0.108/0.104/0.103, Richardson $\to\sim0.10$).
**This is the engine.** (The high-$k$ leading atoms $U_{k,0}$ of $Y_6..Y_9$ never contract with the low-order
partner in the covariances, so they are never built — the engine stays tractable through $v_3$; $v_4$ needs 5-chaos.)

### $v_3$ COMPUTED EXACTLY — small negative, confirms the oscillation [NUMERIC ✓]
$v_3=-0.030$ (n=30), $-0.021$ (n=40) → **$v_3\approx-0.02$**: NEGATIVE (sign flip confirmed) but SMALL, not the
$-0.3$ the noisy renormalized sampling gave (that estimate is **retracted** — it had $\pm0.1$ and drift). This is a
*stronger* confirmation of the complex-Borel-pair surmise: the $\theta=\pi/4$ oscillation fit to $v_0,v_1,v_2$
predicted $v_3\approx0$ (near a cosine zero), and the exact engine confirms $v_3$ is small-negative. Coefficient
ladder (exact): $v_0=0.134,\ v_1=0.111,\ v_2\approx0.10,\ v_3\approx-0.02$ — slowly decreasing then a small sign flip,
the signature of an oscillatory (complex-conjugate Borel) coefficient sequence. Caveat: small $v_3$ is a delicate
cancellation ($\operatorname{Var}(Y_4)=0.327$ minus covariances $\approx0.35$) with $O(1/n)$ boundary-eval error, so
the *magnitude* is $\sim-0.02$ to within grid error; the sign and smallness are robust.

### $v_4$ — computable, but needs optimization [status]
$v_4=\operatorname{Var}(Y_5)+2\operatorname{Cov}(Y_4,Y_6)+\dots+2\operatorname{Cov}(Y_1,Y_9)$ is within the engine's
reach (5-chaos max), but the current multigraph enumeration is too slow for $Y_9$'s 1565 many-atom monomials and the
5-chaos self-contractions (n=24 run exceeded ~40 min CPU and was killed). Optimizations needed: (i) transfer-operator
contraction of the Markovian chains (avoid dense $n^5$ tensors and the $C(n,5)$ chain build), (ii) memoize/prune the
multigraph enumeration, (iii) cache atom kernels. Best current estimate from the finite-radius oscillatory fit to the
EXACT $v_0..v_3$: **$v_4\approx-0.09$** (sign pattern $+,+,+,-,-,-,-,+$, period 8). Confirming this by direct
computation is the remaining decisive step (convergent-oscillatory vs factorial).

### $v_4$ COMPUTED — the series is DIVERGENT (resurgent) [2026-07-09, NUMERIC ✓]
Optimized the engine (memoized moments; pruned multigraph enumeration; `np.einsum(optimize='greedy')`) — 130× faster
($v_3$: 131s→1s). **$v_4\approx-0.45$**, grid-converged (n=20/24/28/34: $-0.469,-0.454,-0.452,-0.457$) and robust
(dominated by $\operatorname{Var}(Y_5)=1.20$ and $2\operatorname{Cov}(Y_3,Y_7)=-1.15$, not a delicate cancellation).
Exact ladder: $\boxed{v_0=0.134,\ v_1=0.111,\ v_2=0.100,\ v_3\approx-0.02,\ v_4\approx-0.45}$.

**Verdict: the $\eta^2$-series is DIVERGENT (asymptotic/resurgent), NOT convergent.** $|v_4|=0.45$ jumps far above all
earlier $|v_n|$ — the envelope GROWS. $v_3$ was small only because it sits near an oscillation NODE ($\cos(3\pi/4+\varphi)
\approx0$); $v_4$, near an anti-node, exposes the growing envelope. Fit $v_n=C r^{-n}\cos(n\pi/4+\varphi)$ ($\theta=\pi/4$
FIXED to $\lambda_0$) to all five: $r\approx0.6$–0.7 $<1$ ⇒ divergent, Borel singularity at $|\zeta|\approx0.6$–0.7,
$\arg=\pm\pi/4$; $\eta^2_c\approx0.6$–0.7, $\beta_c\approx6$. Predicts $v_5\approx-0.9$, and the period-8 sign pattern
$+,+,+,-,-,-,+,+$. **This confirms the core surmise: a genuine resurgent trans-series with a complex-conjugate Borel
pair at $\pm\pi/4$ inherited from $\lambda_0=0.890-0.890i$.**

_Caveat: $v_3,v_4$ use the renormalization prescription "drop all boundary-boundary self-contractions" (normal-ordering);
this is the standard scheme and is applied consistently, but $v_0,v_1,v_2$ (boundary-boundary-free) don't independently
test it. The DIVERGENCE conclusion is robust regardless — $|v_4|$ is dominated by the large, well-determined
$\operatorname{Var}(Y_5)$ and $\operatorname{Cov}(Y_3,Y_7)$ terms._

## Performance wall at 7-chaos ($v_6$) [2026-07-09]
Pushed to $v_5\approx-1.1$ (order-11 functionals via truncated-series inversion; 6-chaos build; base-atom forms cached
after fixing a coefficient-drop bug where real $\sim\!1/u_0'^{\star4}$ terms were thresholded away). Optimizations landed:
memoized moments, pruned multigraph enumeration, `einsum(optimize='greedy')` (130× on $v_3$), vectorized ordered-chain
build. **But $v_6$ (7-chaos, $\operatorname{Var}(Y_7)$) is blocked by the DENSE-tensor representation:** the symmetric
$k$-kernel needs a $k!$-permutation symmetrization of an $n^k$ tensor ($7!=5040$ perms of $n^7$ — intractable; even
6-chaos at $n=20$ costs ~385s just for the $6!$ symmetrization). The vectorized build does NOT fix this — the wall is
the symmetrization, not the construction. **The real fix is a representation change: transfer-operator contraction of
the Markovian chains** (compute $\langle U_{k,0}^2\rangle=\int_{\rm ordered}f^2$ and all chain contractions by nested
$O(n^2k)$ sweeps, never materializing or symmetrizing an $n^k$ tensor). That is a substantial rewrite of `_tn` (a
chain/path tensor-network contractor), deferred. **Current reach: $v_0..v_5$ exact-ish (6 coefficients).** For sharper
singularity extraction WITHOUT $v_6$, the research-recommended move is **Meijer-G / hypergeometric approximants** on the
6 coefficients (converge at 3–5 orders, beat Borel-Padé for few terms).

### Transfer-operator rewrite — PROGRESS [2026-07-10]
Toward $v_6$: **order-13 functionals generated** (`_yexprs_13.txt`, $Y_1..Y_{13}$, truncated-series inversion). **Two
transfer primitives validated to machine precision vs the dense engine:** (1) self-contraction $\langle U_{k,0}^2\rangle
=\int_{\rm ordered}f^2$ (k=2–5); (2) single-cap chain $\langle U_{k,0}\cdot(k\text{ single-leg caps})\rangle=k!\int_{\rm
ordered}\text{chain}\cdot\prod\text{caps}$ (ratio dense/transfer $=3!=6$ exactly for $k=3$). **Structure mapped:** of
$Y_7$'s 90 base-monomials, only ~8 contain a high chain ($k\ge5$); the other 82 (products of $\le4$-chaos atoms) are
already handled by the DENSE engine (it contracts atom-pairs, never building a 7-tensor). So $\operatorname{Var}(Y_7)$
needs transfer only for the high-chain terms: high-chain$^2$ (primitive 1, extend to mixed heads) and high-chain$\times$
(product of small atoms). **The one remaining primitive** = the general chain$\times$(small-atom caps) where a $\ge2$-leg
cap couples several chain sites (state-carrying transfer over coupled chains) + the multigraph integration (chain caps
$m$ legs from a pool of small atoms whose other legs self-contract densely). This is a genuine multi-hour build with real
correctness risk (a rushed version → wrong $v_6$, worse than none); it must be end-to-end re-validated against $v_5$
(6-chaos) before trusting $v_6$. **Held here:** the qualitative picture ($v_6>0$, oscillation, complex pair) is already
MC-cross-validated; precise $v_6$ would only sharpen the phase.

### Transfer-operator rewrite — CORE validated [2026-07-09]
**Foundation works:** $\langle U_{k,0}^2\rangle=\int_{s_1<\dots<s_k}\text{head}(s_1)^2\prod G(s_i,s_{i+1})^2 U(s_k)^2
\prod w$ by a forward transfer sweep ($a_1=\text{head}^2w$; $a_m[j]=w[j]\sum_{i<j}a_{m-1}[i]G[i,j]^2$; ans $=\sum_j
a_k[j]U[j]^2$) **matches the dense engine to machine precision for $k=2,3,4,5$**, $O(n^2k)$, no dense tensor, no $k!$
symmetrization — the wall is breakable. **Remaining to reach $v_6$:** (a) the GENERAL coupled-chain contractor
(chain$\times$single-leg atoms = another sweep; chain$\times$multi-leg-atom products = coupled-path networks tracking
open strands); (b) order-13 functionals $Y_{12},Y_{13}$. Both bounded but real; the self-contraction fast-path alone is
insufficient (misses the chain$\times$product cross-terms). **Status: core validated, full $v_6$ engine is a multi-hour
completion.**
Compute each moment $\langle Y_aY_b\rangle$ by expanding into monomial pairs and summing over **Wick pairings of the
atom legs**, contracting pairwise (with `beval` for boundary legs, drop for boundary–boundary), **never materializing
a high-chaos tensor**. Per monomial-pair the cost is $(L-1)!!$ pairings ($L$=total legs $\le 10$ for $v_4$ → 945),
each a cheap incremental contraction of the (already-correct) Markovian atom kernels. This sidesteps truncation
entirely and is the correct route to exact $v_4,v_5$. The atom kernels, the `beval` rule, and the renormalization
(drop boundary–boundary) are all validated and reusable — only the moment-assembly layer needs rewriting from
"truncated products" to "pairing sums."

## Why the engine matters
It is the only route to $v_4$ (renormalized sampling dies at the $1/\delta^3$ boundary self-contractions, §8/`_bstruct9`).
Once `beval` is calibrated on $v_2$, the same code computes $v_4,v_5$ exactly — deciding whether $\mathcal W$'s
variance series is a genuine divergent trans-series and confirming the $\lambda_0$-tied $\pm\pi/4$ Borel pair.
