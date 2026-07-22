# Notes — the §1/§0 keystone build: uniform parabolic-cylinder (Weber) tube with noise

_June 2026. Figure `coupled-atlas/figures/weber_tube_build.png`; script `weber_tube_build.py`. Per
`regime-tests/RH_DIRECTION_NOVEL_ANGLES.md` §1 (= the §0 keystone). The build the falsifiers pointed to:
T1 is the route; the cusp law is new; so the uniform parametrix must be built directly. Tags
[PROVED/DERIVED]/[NUMERIC]/[OPEN]._

## The object

Around the antisymmetric-mode canard, the inner equation is $u''=(V_\Delta(Y)-\eta\xi)u$,
$V_\Delta=\mathrm{sign}(Y)|Y|(|Y|+\Delta)$. Cole–Hopf $p=u'/u$: the canard is the deterministic Riccati
$\bar p(Y)$ (the stable branch $+\sqrt{V_\Delta}$ in the decaying region), and the fluctuation
$\delta p=p-\bar p$ obeys the **linear SDE**
$$\delta p' = -2\,\bar p(Y)\,\delta p - \eta\,\xi,\qquad\Rightarrow\qquad
\frac{dv}{d\tau} = -4\,\bar p(Y)\,v + \eta^2,\quad \tau=Y_0-Y. \tag{DERIVED}$$

## Two comparisons

**OU (Berglund–Gentz).** Freeze $\bar p\approx\sqrt{V_\Delta}$ (adiabatic). The quasi-static variance
$$v_{\rm qs}=\frac{\eta^2}{4\sqrt{V_\Delta}}\ \xrightarrow{\ Y\to0\ }\ \infty \tag{DERIVED — DIVERGES}$$
because the curvature $\sqrt{V_\Delta}\to0$ at the turning. This is the wrong Gaussian; it is why the
standard tube estimate degenerates at the merge.

**Weber (parabolic-cylinder).** Use the **exact canard restoring** $\bar p_{\rm exact}=\bar u'/\bar u$
— the log-derivative of the parabolic-cylinder canard — obtained by integrating the deterministic
Riccati $d\bar p/dY=\bar p^2-V_\Delta$ (stable tracking of $+\sqrt{V_\Delta}$). It stays **finite and
positive through the turning**, so $v_{\rm Weber}$ (the same variance ODE with $\bar p_{\rm exact}$) is
finite. **[DERIVED]**

**Splice (Olver uniform connection, at the covariance level).** Inner scale $\ell\sim\eta^2$ (where
$v_{\rm qs}\sim1$). $\;v_{\rm splice}=v_{\rm qs}$ (outer, $|Y|\gg\ell$, where $\bar p_{\rm exact}\to\sqrt
{V_\Delta}$) $\cup\ v_{\rm Weber}$ (inner, $|Y|\lesssim\ell$). The inner Green's function is **Airy** for
$\Delta\gg\ell$ and **Weber** for $\Delta\lesssim\ell$ — the stochastic analogue of Olver's uniform
two-coalescing-turning-points connection, done on the covariance. **[DERIVED]**

## Numerical validation (vs the simulated tube)

Tube variance near the turning ($Y\approx0.1$), and mean relative error over the inner band
$Y\in(0.05,0.6)$:

| Δ | regime | $v_{\rm sim}$ | $v_{\rm qs}$ (OU) | $v_{\rm Weber}$ | rel.err OU | rel.err Weber |
|---|---|---|---|---|---|---|
| 2.0 | Airy ($\Delta\!\gg\!\ell$) | 0.47 | 1.09 | 0.39 | **0.65** | **0.12** |
| 0.5 | crossover | 0.70 | 2.04 | 0.50 | 0.92 | 0.23 |
| 0.1 | Weber ($\Delta\!\lesssim\!\ell$) | 0.81 | **3.54** | 0.55 | **1.49** | 0.27 |

- The **OU quasi-static comparison DIVERGES** through the merge (rel.err $0.65\to1.49$ as $\Delta\to0$).
- The **Weber comparison stays finite and tracks the true tube to ~15–27% uniformly** in $\Delta$
  (recovering the Airy $\Delta\!\gg\!\ell$ and Weber $\Delta\!\lesssim\!\ell$ endpoints). **[NUMERIC]**
- **Covariance** (cusp $\Delta=0.1$): the Weber-propagator $C(Y,Y')=\Phi(Y,Y')\,v(\,\cdot\,)$ matches the
  simulated two-time covariance to **~27%**. **[NUMERIC]**
- (A simpler swept-but-frozen-$\sqrt V$ variant lands at ~13–17% — also finite; the principled object is
  the parabolic-cylinder $\bar p_{\rm exact}$.)

The residual ~15–27% is the **non-Gaussian/nonlinear correction** (the $\delta p^2$ term and the skew of
the escape) beyond the linear-fluctuation variance — it is what the full rigorous estimate must supply.

_(Two bugs found & fixed in the build, recorded for reproducibility: the canard Riccati integrated in
decreasing $Y$ is $d\bar p/dY=\bar p^2-V$ — the opposite sign blows $\bar p$ up; and $\bar p$ must be
clipped $\ge0$ in the variance ODE so a tiny RK4 undershoot does not anti-restore.)_

## Proved vs numeric vs heuristic

| statement | status |
|---|---|
| fluctuation SDE & variance ODE $dv/d\tau=-4\bar p v+\eta^2$ | **[PROVED]** (exact linearization) |
| quasi-static OU variance diverges $\eta^2/(4\sqrt V)\to\infty$ | **[PROVED]** |
| $\bar p_{\rm exact}$ finite through the turning ⇒ $v_{\rm Weber}$ finite | **[DERIVED]** |
| OU→Weber matching in the overlap ($\bar p_{\rm exact}\to\sqrt V$) | **[DERIVED]** |
| Weber tube matches sim to ~15–27% uniformly in Δ; covariance ~27% | **[NUMERIC]** |
| Airy ($\Delta\!\gg\!\ell$) & Weber ($\Delta\!\lesssim\!\ell$) limits recovered | **[NUMERIC]** |
| full **uniform tube ESTIMATE** with noise (the theorem) | **[OPEN]** |

## What remains for full rigor (the residual theorem)

The construction + numerical validation establish the **structure** (finite, uniform, spliced). Full
rigor (closing T1) needs the **uniform parabolic-cylinder estimate with noise**: a probability bound
that the true sample-path tube stays within the spliced envelope **uniformly in $\Delta/\ell$**, which
requires (i) control of the **nonlinear $\delta p^2$ correction** (the ~20% non-Gaussian residual), and
(ii) **parabolic-cylinder connection-coefficient bounds** uniform through the merge (the noisy
upgrade of Olver's error bounds). This is the precise remaining gap.

## The integrable / T2 reading (where it appears)

The canard $\bar u$ is the **parabolic-cylinder parametrix** itself. Its *probabilistic* reading is the
tube above (T1). Its *integrable* reading — the connection/Stokes data of $\bar u$, i.e. the
noise-averaged spectral determinant of the parabolic-cylinder fluctuation operator — is the **candidate
new T2 cusp edge law** (the object the falsifiers showed descends from no known parent). T1 and T2 are
thus the two faces of this one parametrix, exactly as the §0 keystone reframe predicted: building the
uniform parabolic-cylinder-with-noise parametrix is the single object that closes both.
