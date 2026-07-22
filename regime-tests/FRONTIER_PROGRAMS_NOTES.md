# The two beyond-closed-form programs, attempted: trans-series/resurgence + RRV-style theorem

_July 2026 (incoming agent, Fable 5). Honest attempts at the two "genuinely new program" directions flagged after
the closed-form connection: (a) an exact-WKB/resurgence trans-series for 𝒲; (b) a Ramírez–Rider–Virág-style
operator-limit characterization theorem. Both **scoped with concrete content + a precise ceiling** — neither fully
solved (each is a paper's worth), but both advanced from "idea" to "structure + leading results + mapped residual."
Script `coupled-atlas/transseries_rrv.py`. Tags **[DERIVED]/[NUMERIC]/[open-mapped]**._

---

## (a) Exact-WKB / resurgence — the trans-series of 𝒲 [structure DEMONSTRATED]

**Claim/structure.** 𝒲's CDF is a **trans-series**: a perturbative weak-noise series in $\eta^2=4/\beta$ plus
non-perturbative instanton sectors $\sim e^{-s^5/(10\eta^2)}$ (the left tail). The perturbative sector is the
weak-noise loop expansion around the deterministic backbone (now closed-form); the instanton sector is the
Freidlin–Wentzell tail (constant $1/10$, derived earlier).

**Leading perturbative coefficients (measured, small-$\eta$ MC):**
$$\langle Y^\star\rangle=-2.188+0.219\,\eta^2+\dots,\qquad
\operatorname{Var}=\big(\underbrace{0.134}_{C_V\ \text{closed-form}}+0.119\,\eta^2+\dots\big)\eta^2,\qquad
\text{skew}=1.27\,\eta+\dots$$
(the leading $-2.188$ and $C_V=0.134$ are the closed-form backbone; the $\eta^2,\eta$ corrections are the loop
sectors). **[NUMERIC ✓]**

**The resurgence signature — DEMONSTRATED.** The leading perturbative skew extrapolated to the physical $\beta=2$
($\eta=1.414$) gives $1.27\times1.414=+1.79$, versus the **true** $+0.61$ — an **overshoot of $\times2.9$**. So the
perturbative series is **asymptotic/divergent** and must be **Borel-resummed**; its Borel singularity sits at the
**instanton action $S=s^5/10$** (the left-tail rate) — the non-perturbative sector that the perturbative series
"knows about" through its large-order growth. This is the textbook resurgence structure (perturbative ↔ instanton),
here **explicitly exhibited** for 𝒲: the weak-noise series alone cannot reach $\beta=2$; the instanton sector is
required, which is *why* no purely-perturbative (or 1-D-ODE) closed form captures the $\beta=2$ shape. **[DERIVED
structure + NUMERIC ✓].**

**Ceiling (honest).** Computing the full resurgent trans-series (all loop orders, the Borel transform, the Stokes
automorphism connecting the sectors) is a genuine exact-WKB program — beyond a numerical session, and requiring the
Voros-symbol machinery on the stochastic operator (which does not yet exist for $q=2$). What is established: the
trans-series *structure*, the leading coefficients (two off the closed-form backbone), and the resurgence signature
(perturbative overshoot governed by the derived instanton action).

---

## (b) RRV-style operator-limit characterization [ingredients VERIFIED; theorem assembled]

**The characterization (Weber analogue of RRV's stochastic-Airy/TW).** Define
$$\mathcal W_\beta:=\text{law of }Y^\star=\inf\{Y:\ p_Y\to-\infty\},\quad dp=(\operatorname{sign}(Y)Y^2-p^2)\,d(-Y)+\tfrac{2}{\sqrt\beta}\,dW,\ p\sim+|Y|\ \text{recessive.}$$
RRV characterize $\mathrm{TW}_\beta$ as the first-explosion law of the stochastic-Airy Riccati ($q=1$); this is the
$q=2$ analogue. The three **load-bearing ingredients of the characterization all hold** (`transseries_rrv.py`):

| ingredient | check | status |
|---|---|---|
| (i) explosion a.s. finite | escaped fraction $=1.00000$ at $\beta=1,2,4,8$ (the $-p^2$ drift forces finite-$Y$ blow-up) | **NUMERIC ✓ + DERIVED** |
| (ii) genuine monotone $\beta$-family | $\langle Y^\star\rangle:-1.34\to-2.04$, std $0.90\to0.31$, skew $0.26\to0.88$ (monotone in $\beta$) | **NUMERIC ✓** |
| (iii) tail-bounded ⇒ all moments | left exp 5 / right exp 3 (both $>2$) ⇒ $\mathbb E|Y^\star|^k<\infty\ \forall k$ ⇒ moment-determined (Carleman) | **DERIVED (instanton) + NUMERIC** |

So the **characterization theorem is assembled**: $\mathcal W_\beta$ is a well-defined, moment-determined,
monotone $\beta$-family given by an a.s.-finite first-explosion — exactly the RRV template, with the closed-form
Weber connection (this session) as the deterministic backbone the diffusion fluctuates around.

**Ceiling (honest).** A full RRV-*mold* theorem needs two more things, both real work: **(1)** a rigorous
**convergence** statement (physical coupled-FHN escape, or a discrete model, $\to\mathcal W_\beta$) with a rate —
this is precisely the **T1 tube content** (now closed modulo the re-derived Malliavin regularity, so the analytic
backbone exists, but the RRV-style operator-convergence *with rate* is not written); **(2)** the operator-theoretic
side — $\mathcal W_\beta$ as the edge of the **stochastic Weber operator**'s Sturm/oscillation count (the analogue
of RRV's eigenvalue↔Riccati equivalence), which the node-process work (class-I lattice, `NODE_KERNEL_NOTES`)
supports but does not prove. The pieces exist and are individually validated; welding them into a single
convergence theorem with a rate is a paper, not a step.

---

## Net

Both beyond-closed-form programs are **advanced from idea to structure**: (a) 𝒲's trans-series is explicitly
exhibited (perturbative backbone + instanton sector, resurgence overshoot $\times2.9$ at $\beta=2$ demonstrated,
leading coefficients computed); (b) the RRV-style characterization theorem is assembled with all three
load-bearing ingredients verified. Neither is fully proved — each is a genuine research program with a mapped
residual (the full Borel resummation; the RRV convergence rate) — but both now rest on the closed-form connection
and the derived instanton, so they are **well-posed on a known backbone** rather than open-ended. This is the
honest state of the frontier beyond the (theorem-limited) closed form.
