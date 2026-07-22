# Constructing the interface/PT τ-function for 𝒲 — attempt, a three-wall obstruction, and the map forward

_June 2026. Direct continuation of `T2_CONNECTION_DATA_DERIVATION.md`. Attempts to build the τ-function /
Fredholm determinant the cusp law 𝒲 lives on (two conjugate-type rank-2 irregular singularities glued at
$Y=0$). Figure `coupled-atlas/figures/resonance_landscape.png`; script `resonance_determinant.py`. Tags
**[DERIVED]/[CITED]/[CONJECTURAL]/[NUMERIC]**; ⚑ = load-bearing._

---

## 0. Outcome (honest, up front)

| step | result | status |
|---|---|---|
| Glued/interface RH problem | written precisely — confining $U$-parametrix ⊕ oscillatory $W$-parametrix ⊕ interface matrix $J_0$ at $Y=0$ | **DERIVED** |
| Resonance determinant $D(\lambda)$ | explicit building block; **complex resonances** $\lambda_0=1.64-1.49i$ (string at $\mathrm{Re}\,\lambda\approx1.6,3.6$, width $\Gamma\!\approx\!3$) | **DERIVED + NUMERIC ✓** |
| τ / Fredholm for the *distribution* 𝒲 | **does not close** — three independent walls, each located precisely | **OBSTRUCTION** ⚑ |
| The obstruction | continuous spectrum (no discrete edge) · complex resonances (non-self-adjoint kernel) · isomonodromy fixed point (no flow to host a τ) | **DERIVED** ⚑ |
| The map forward | **complex-scaling (Aguilar–Balslev–Combes) non-Hermitian Fredholm determinant**; closed-form $D(\lambda)$ as the completable sub-problem; noncommutative-PII template for the process | **mapped + CITED** |

**Verdict.** Per the brief ("a precise obstruction + a map is a fully acceptable outcome — better than a forced
partial calculation"): I did **not** force a τ-function. The standard determinant provably cannot form (three
walls). The object 𝒲 lives on is now **named concretely** — the **non-Hermitian resonance determinant of the
half-confining/half-inverted Weber operator**, obtained by *complex scaling*, of **noncommutative-PII type**
for the multi-line/process version — and one **completable sub-problem** (closed-form $D(\lambda)$) is isolated.

---

## 1. The glued/interface RH problem (precise)

Seek $\Psi(Y)\in SL_2(\mathbb C)$, the fundamental matrix of $\Psi'=\begin{pmatrix}0&1\\ V-\lambda&0\end{pmatrix}\Psi$,
$V=\operatorname{sign}(Y)Y^2$, sectionally normalised by the two end-parametrices:

**Confining parametrix ($Y>0$).** With $z=\sqrt2\,Y$, $a=-\lambda/2$,
$$\Psi_+\sim P_U(Y)\,e^{+\frac{Y^2}{2}\sigma_3}\,Y^{-(a+\frac12)\sigma_3},\qquad P_U\ \text{built from } U(a,z),U(a,z)\!\!\uparrow\text{(dominant)}.$$
Crossing the two Stokes rays $\arg z=\tfrac\pi4,\tfrac{3\pi}4$ multiplies by unipotent $S_{1,2}$ with the
**classical PC multiplier** $s=\dfrac{\sqrt{2\pi}}{\Gamma(\frac12-\frac\lambda2)}e^{-i\pi\lambda/2}$ (poles
$\lambda=1,3,5,\dots$). **[CITED — DLMF 12.]**

**Oscillatory parametrix ($Y<0$).** Same $a=-\lambda/2$ but the DLMF $W(a,\cdot)$ (real/oscillatory) functions,
$$\Psi_-\sim P_W(Y)\,e^{-i\frac{Y^2}{2}\sigma_3}\,|Y|^{-\frac12\sigma_3},\qquad P_W\ \text{built from } W(a,\sqrt2 Y),\,W(a,-\sqrt2 Y),$$
the modulus factor $|\Gamma(\tfrac12-\tfrac{i\lambda}2)|^2=\pi/\cosh\tfrac{\pi\lambda}2$. **[CITED + DERIVED.]**

**Interface jump at $Y=0$ ⚑.** The solution $(u,u')$ is $C^1$ through $Y=0$ (the potential is $C^1$), so there
is **no jump in the solution**; the jump is in the *asymptotic frame* — the constant matrix $J_0$ converting the
confining basis $\{U,U^\uparrow\}$ into the oscillatory basis $\{W_{\rm out},W_{\rm in}\}$. It is the **DLMF
$U\!\leftrightarrow\!W$ connection matrix** (12.14.4–6), i.e. the rotation $z\to e^{i\pi/4}z$ + parameter
conjugation $a\to ia$:
$$\boxed{\,\Psi_-(Y)=\Psi_+(Y)\,J_0(\lambda),\qquad J_0(\lambda)=\frac{1}{\sqrt{2\pi}}\begin{pmatrix}\Gamma\text{-entries in }a=-\frac\lambda2 \text{ and } ia\end{pmatrix},\ Y\to0^\mp.}$$
Concretely the matching data is fixed by the DLMF zero-values
$U(a,0)=\dfrac{\sqrt\pi\,2^{-1/4-a/2}}{\Gamma(\frac34+\frac a2)}$, $U'(a,0)=-\dfrac{\sqrt\pi\,2^{1/4-a/2}}{\Gamma(\frac14+\frac a2)}$
(and the analogous $W(a,0),W'(a,0)$). **[DERIVED.]** The RH data is thus **two rank-2 irregular singularities of
conjugate type** joined by the constant $J_0$ — *not* a Fuchsian/standard meromorphic connection.

---

## 2. The resonance determinant — explicit building block, validated complex

$H=-\partial_Y^2+\operatorname{sign}(Y)Y^2$ has spectrum $=\mathbb R$ (unbounded below via $-Y^2$): **no normalisable
states**, only **resonances** — confining quasi-levels that *leak out the inverted side*. The resonance/Jost
function is the **incoming amplitude of the recessive-at-$+\infty$ solution at $-\infty$**:
$$D(\lambda)=\big[\text{incoming coeff of }U(-\tfrac\lambda2,\sqrt2\,Y)\text{ continued to }Y\to-\infty\big],\qquad D(\lambda)=0\Leftrightarrow\text{resonance}.$$
In Γ-form it is the product the previous note derived, $D(\lambda)\propto \Gamma(\tfrac12-\tfrac{i\lambda}2)/\Gamma(\tfrac12-\tfrac\lambda2)\times$phase.

**Numerics (`resonance_determinant.py`, vectorised complex ODE).** Lowest resonance
$$\boxed{\ \lambda_0=1.64-1.49\,i\ }\quad(\text{string at }\mathrm{Re}\,\lambda\approx1.6,3.6,\dots,\ \text{spacing}\approx2,\ \mathrm{Im}\,\lambda\approx-1.5).$$
**[NUMERIC ✓]** The **finite width $\Gamma\approx3$** (broad!) is physically exact: the cusp barrier *vanishes*
at $Y=0$, so the confining levels leak strongly. The real parts sit near the half-oscillator levels $1,3,\dots$
(shifted up $\sim0.6$ by the strong coupling). **This is the concrete non-self-adjoint spectral datum** — the
resonances are genuinely complex, unlike the real TW/Weber edge. (Figure: the resonance string in the lower
half $\lambda$-plane.) **[DERIVED + NUMERIC.]**

---

## 3. The τ / Fredholm attempt — and exactly where it walls

**Goal.** Write the CDF $F(s)=\mathbb P(Y^\star<s)$ of 𝒲 as a (generalised) Fredholm determinant
$\det(1-K_s)$ or τ-function. For Tracy–Widom this works on three pillars; **each fails here**, independently:

**WALL-A — no discrete edge (continuous spectrum).** TW is "$\mathbb P(\text{top eigenvalue}<s)$": there *is* a
top eigenvalue. Here $\operatorname{spec}H=\mathbb R$ (the inverted side is unbounded below), so **there is no
top eigenvalue** to be 𝒲-distributed. The edge is defined **operationally** by *first passage* (first node going
down the inverted side), not by a spectral edge. "$\,$No eigenvalue above $s\,$" is not the relevant event ⇒ the
gap-probability route does not start. **[DERIVED ⚑]**

**WALL-B — complex resonances (non-self-adjoint kernel).** Even replacing "eigenvalue" by "resonance," §2 gives
$\lambda_0=1.64-1.49i$ — **complex**. The natural kernel $K_s$ (resonance Green's function projected on
$(s,\infty)$) is therefore **non-self-adjoint**, and the inverted-side continuum makes it **not
Hilbert–Schmidt** without regularisation. A standard (self-adjoint, trace-class) Fredholm determinant **does not
exist**. **[DERIVED ⚑]**

**WALL-C — isomonodromy fixed point (no flow).** A τ-function is defined *along an isomonodromic flow* via the
JMU 1-form $d\log\tau=\omega$. The previous note proved the cusp is a **fixed point**: its only PC-class
deformation (the linear $\Delta Y$ term) is *relevant* (flows to fold), so there is **no class-preserving flow**
to integrate $\omega$ along. The τ-differential has **no flow direction** ⇒ no τ-function in the deformation
sense. **[DERIVED ⚑, uses rung E.]**

**Conclusion.** The standard (self-adjoint, discrete-edge, flow-based) determinant/τ **provably cannot form**.
I stopped here rather than force a partial determinant — the three walls are structural, not technical.

---

## 4. The obstruction, consolidated

> 𝒲 is the **first-passage law of a continuous-spectrum, non-self-adjoint (resonance) operator at an
> isomonodromy fixed point.** None of the three pillars of the Tracy–Widom determinant survives: the edge is
> *continuous* (no top eigenvalue), the resonances are *complex* (no self-adjoint kernel), the deformation is a
> *fixed point* (no τ-flow). The object is genuinely outside the self-adjoint, flow-based Painlevé/Fredholm
> paradigm. **[DERIVED]**

This is *why* all four earlier negatives held (soft-edge Fredholm / isomonodromic time / Pearcey / Stokes
continuation) — each presupposes one of the three failed pillars.

---

## 5. The map forward — the minimal missing ingredient

**(M1) Recommended: the complex-scaling (Aguilar–Balslev–Combes) non-Hermitian determinant. ⚑** Rotate the
inverted side $Y\to Ye^{-i\theta}$ ("complex scaling"/dilation analyticity). This is the standard cure for
WALL-A+B: it **rotates the continuous inverted spectrum off $\mathbb R$ and exposes the resonances $\lambda_n$ as
genuine discrete (complex) eigenvalues.** Then
$$F(s)\ \overset{?}{=}\ \det\nolimits_{\theta}\!\big(1-K^{(\theta)}_s\big),\qquad K^{(\theta)}_s=\text{resonance kernel on the rotated contour,}$$
a **non-Hermitian Fredholm determinant** — exactly the class for which a Riemann–Hilbert/Lax theory exists. This
is the concrete candidate for the τ-object, and the **named generalised isomonodromy**: an *interface* RH problem
whose oscillatory end is complex-rotated to a confining (resonance) end. **[CONJECTURAL — the equality is the
open step; the construction is concrete and templated.]**

**(M2) Alternative: the "thermal" τ.** $\beta=4/\eta^2$ rescales the imaginary-$\Gamma$ argument
($\cosh\tfrac{\pi\lambda}2$ is literally a thermal/Bose factor). Treat $\beta$ as a *non-isomonodromic* thermal
deformation and build a partition-function-like $\tau_\beta$ (finite-temperature determinant). This sidesteps
WALL-C (no isomonodromic flow needed — it is a *thermal* flow). **[CONJECTURAL.]**

**(M3) The completable sub-problem worth doing first.** Close the **resonance determinant $D(\lambda)$ in
explicit Γ-form** (we have $\lambda_0$ numerically; the closed product $D(\lambda)\propto\Gamma(\tfrac12-\tfrac{i\lambda}2)/\Gamma(\tfrac12-\tfrac\lambda2)\cdot e^{i\varphi}$
is within reach by matching the DLMF zero-values of §1). This *is* the interface connection matrix $J_0$'s
determinant — the honest, finishable deliverable, independent of the harder distributional τ. **[DERIVED-in-reach.]**

**The process (intrinsic DBM-edge) version.** The multi-line "Dyson–Weber" object is the **matrix/operator**
version of (M1): a complex-scaled determinant with a *matrix* kernel — directly the **noncommutative-PII Fredholm
determinant** template (Bothner–Cafasso–Tarricone). So B's open "intrinsic process" and this τ are the *same*
construction at rank 1 vs rank $n$.

---

## 6. Literature templates (cited, not reinvented)

- **Non-Hermitian / matrix (noncommutative) PII Fredholm determinants** — Bothner, Cafasso, Tarricone and
  collaborators: $\det(1-K)$ for matrix Airy operators ↔ noncommutative PII via RH/IIKS. The template for a
  **non-self-adjoint determinant** and for the **multi-line process**. [arXiv:2505.16830, 2007.05707, 1101.3997]
- **Confluence isomonodromy** (merging poles → rank-$r$ irregular singularity; generalised isomonodromic
  deformation, JMU): the framework for the **two-rank-2-singularity** object. [CMP 2023, *Isomonodromic
  Deformations: Confluence, Reduction and Quantisation*; Klimeš, *Confluence of singularities…*]
- **Inverted oscillator resonances** = complex eigenvalues, parabolic-cylinder eigenmodes, **anti-PT-symmetric**;
  **complex scaling / dilation analyticity** (Aguilar–Balslev–Combes) turns resonances into eigenvalues. [arXiv:2204.10780;
  arXiv:quant-ph/0703234]

---

## 7. Net + tracker

**Achieved:** the glued/interface RH problem written precisely with its Γ-function interface matrix; the
resonance determinant built and its **complex resonance $\lambda_0=1.64-1.49i$ validated** (concrete
non-self-adjoint datum); a **precise three-wall obstruction** to the standard τ/Fredholm; and a **concrete map**
— complex-scaling non-Hermitian determinant (M1), with the closed-form $D(\lambda)$ as the completable
sub-problem (M3), and the noncommutative-PII template for the process.

**Not achieved (honestly):** the τ-function/determinant for the *distribution* 𝒲 is **not constructed** — and is
shown to require leaving the self-adjoint/flow paradigm (the named complex-scaled non-Hermitian object).

**Tracker:** B stays **95%** (the τ is not built; the frontier is now sharply posed and a building block is
validated — characterisation up, construction still open). The recommended next concrete step is **M3**
(closed-form $D(\lambda)$), then **M1** (the complex-scaling determinant). This is the shared lever: M1 at rank 1
is 𝒲 (B), its rank-$n$ matrix version is the intrinsic process (B/D), and the whole obstruction is the rung-E
relevance statement read as "no flow." 
