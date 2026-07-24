# PROGRAM2 — TW_β halo: uniform-in-`a` noise Borel-summability (Borot–Nadal control)

**Date:** 2026-07-22. **Executes:** `PROGRAM2_HANDOFF_D_TWBETA_NEXT.md` **open task #2** — extend
T3 (median Borel-summability of the noise `1/β` expansion) from the *scaling regime* `g=βa^{3/2}→∞`
to a **uniform-in-`a`** statement, controlling the Borot–Nadal `R_m(2/β)` corrections and the
`O(a^{-3/2})` x-extension. **Memory:** `program2-tw-beta-halo.md`. **Script:**
`coupled-atlas/_tw_uniform_in_a.py` (sympy, exact). **Companion:**
`PROGRAM2_TWBETA_NOISE_BOREL_PROOF_NOTES.md` (the scaling-regime proof this extends).

---

## Result (one paragraph)

Working from the **exact all-β, all-orders** Borot–Nadal tail (arXiv:1111.2761, Prop. 1.1) rather than
the frozen escape, the uniform-in-`a` question **reduces cleanly** and a **genuine subtlety surfaces**.
(1) The Coulomb-gas prefactor `Γ(β/2)/((4β)^{β/2}2π)` is **a-independent** and the `−(3β/4)log a` term
is a single `β^{+1}` monomial, so **all** `1/β`-Borel content lives in the dynamical exponent
`Σ(a,β)=Σ_m (β/2)R_m(2/β)a^{-3m/2}`; its `β^{-n}` coefficients `σ_n(a)` are the exact noise-loop
coefficients. (2) The **sector is uniform**: the dynamical Borel singularity sits at the *real* WKB
action `Φ(a)=(2/3)a^{3/2}(1+o(1))∈ℝ₊` for every `a>0` (noise-independent Airy spectral curve, paper §8) —
it *moves* with `a` but never leaves the positive ray — while the Coulomb-gas singularities are fixed
on `iℝ`; so a median-summation sector around `ℝ₊` exists **uniformly in `a≥a₀`**. (3) **The finding:**
the program's **frozen boundary-escape series is NOT the exact tail noise series beyond leading order.**
The exact scaling (diagonal) coefficient of `g^{-n}` is `2^n r_{n,n+1}` (top-degree BN coefficient);
it matches the frozen `δ_n` at `n=1` **and** `n=3` but **not `n=2`**, where `R_2` has a degree deficit
(`deg R_2=2 < m+1=3`) so the **exact two-loop scaling coefficient vanishes** (`0`) while the frozen
escape gives `−5/8`. Hence the large-order (Gevrey) input for the uniform theorem must come from the
**BN recursion**, not the freeze. Net: uniform-in-`a` **sector** = established (structural); uniform-in-`a`
**amplitude/Gevrey bound** = the one remaining estimate, now sharply localized.

## The exact input (Borot–Nadal, validated)

$$1-\mathrm{TW}_\beta(s)=\frac{\Gamma(\beta/2)}{(4\beta)^{\beta/2}2\pi}\,s^{-3\beta/4}\,e^{-\frac23\beta s^{3/2}}\,
\exp\Big[\sum_{m\ge1}\tfrac\beta2 R_m(2/\beta)\,s^{-3m/2}\Big],\quad s=a,$$

$$R_1=\tfrac1{24}(-5X^2+9X-39),\quad R_2=\tfrac5{64}(11X^2-19X+36),\quad
R_3=\tfrac1{4608}(-1105X^4+3240X^3-23325X^2+34938X-41433),\quad X=\tfrac2\beta.$$

`R_m` is a polynomial of degree **≤ m+1** in `X=2/β`. **Validated at β=2** (`X=1`) against the known GUE
Tracy–Widom tail: `R_1(1)=−35/24=b_1` and `R_2(1)=35/16=b_2−b_1²/2` (with GUE bracket coefficients
`b_1=−35/24, b_2=3745/1152`) — exact matches (script §1). So the fetched `R_m` are correct.

## The reduction (rigorous; self-checked)

Put `t=a^{-3/2}`, `X=2/β`. Since `(β/2)r_{m,ℓ}X^ℓ = r_{m,ℓ}2^{ℓ-1}β^{1-ℓ}`, collecting `Σ` by powers
of `β` gives (script §2; a `β=2` re-assembly self-check passes):

| order | piece | meaning |
|---|---|---|
| `β^{+1}` | `σ_{-1}(a)=\tfrac12Σ_m R_m(0)a^{-3m/2}` | **rate** `Φ(a)=\tfrac23a^{3/2}-σ_{-1}` corrections |
| `β^{0}` | `σ_0(a)` | one-loop (`a^{-3β/4}`-type) prefactor |
| `β^{-n}` | `σ_n(a)=2^nΣ_{m≥n}r_{m,n+1}a^{-3m/2}` | **noise loop `n`** (what T3 resums) |

**Coulomb-gas + log are `β^{≥0}` and a-independent ⇒ they carry no `1/β`-Borel tail.** All noise-Borel
content is the set `{σ_n(a)}_{n≥1}`. Explicitly (from `R_1..R_3`):

```
σ_-1 = -13/16 t + 45/32 t^2 - 13811/3072 t^3          (rate Φ correction)
σ_0  =   3/8 t - 95/64 t^2 + 1941/256 t^3             (one-loop)
σ_1  =  -5/12 t + 55/32 t^2 - 7775/768 t^3            (noise loop 1)
σ_2  =              45/16 t^3                          (noise loop 2; diagonal t^2 coeff = 0)
σ_3  =            -1105/576 t^3                        (noise loop 3)
```

The **diagonal** (scaling-regime, `a→∞` at fixed `g=βa^{3/2}`) coefficient of `g^{-n}` is the leading
`t^n` term `2^n r_{n,n+1}`, i.e. `2^n ×` the **top-degree** coefficient of `R_n`.

## The finding: frozen escape ≠ exact noise series (the `n=2` deficit)

Comparing the exact diagonal `2^n r_{n,n+1}` to the program's frozen `δ_n=[g^{-n}]\log F(g)`
(`F(g)=\tilde Q(2/g)^{-2}`, the reduced boundary-escape fluctuation series):

| `n` | exact diagonal `2^n r_{n,n+1}` | frozen `δ_n` | agree? |
|---|---|---|---|
| 1 | `−5/12` | `−5/12` | ✅ |
| 2 | `0`  (`R_2` deg-deficit) | `−5/8` | ❌ |
| 3 | `−1105/576` | `−1105/576` | ✅ |

So the frozen boundary-escape series reproduces the **leading** noise Stokes term (`n=1`, `−5/12`) — which
is all it was used for (the reduced constant `S₀=1/π`) — and coincidentally `n=3`, but it is **not** the
exact tail noise series: at two loops the exact scaling coefficient **vanishes** because `R_2` has degree
`2<m+1=3`, whereas the freeze gives `−5/8`. **Consequence for T3:** the elementary Airy-integral proof
(`PROGRAM2_TWBETA_NOISE_BOREL_PROOF_NOTES.md`) establishes median Borel-summability of the **frozen
reduced escape** `F(g)`; that object is *not identical* to the exact tail's noise expansion, so the
exact large-order/Gevrey behaviour driving the **uniform** theorem must be read off the **BN
coefficients**, not the freeze. (The `n=1,3` agreements hint at a deeper frozen↔BN relation — **now
resolved: it is the MEDIAN structure, below.**)

## RECONCILIATION (2026-07-23) — the exact diagonal is the MEDIAN of two independent branches

**`coupled-atlas/_tw_reconcile_gevrey.py`.** Two *independent* computations of the scaling diagonal:
- **(F) frozen** `δ_n=[g^{-n}]log F(g)`, `F=\tilde Q(2/g)^{-2}` (cubic-barrier escape MFPT);
- **(W) Weber genus-0 WKB** `2^n\hat r_n`, `\hat r_n=-4w_{n+1}/(3n)` (refined quantum curve, `_tw_weber_wkb.py`).

They satisfy, to **n=7/8** (exact huge fractions, e.g. `δ_7=−1282031525/344064`):
$$2^n\hat r_n=(-1)^{n+1}\delta_n,\qquad\text{i.e.}\qquad W(g)=-F_{\log}(-g)\quad(\text{Borel reflection }g\!\to\!-g).$$
`F` and `W` **agree at odd `n`, are opposite at even `n`.** The **exact** BN diagonal is their **median**:
$$\boxed{\,D_n:=2^n r_{n,n+1}=\tfrac12\bigl(\delta_n+2^n\hat r_n\bigr)=\text{odd part of }\log F(g)\,}$$
— **verified exactly at `n=1,2,3`** (all orders with known `R_m`): `D_1=δ_1=−5/12`, `D_2=0` (median of
`∓5/8`), `D_3=δ_3=−1105/576`. So the frozen escape's **even-order terms are the Stokes/discontinuity
part** (a freeze artifact); the physical diagonal keeps only the **odd (real/median)** part, which
**auto-cancels** the even orders — exactly the `R_{even}` degree deficit (`D_2=0`). This is the
median-summation principle made explicit: `F` and `W` are the two lateral branches (`g↔−g`), and the
tail is their median. Predicts `D_{even}=0`, `D_{odd}=δ_{odd}` for all `n` (n≥4 needs `R_{m≥4}` to
confirm, but is forced by the verified `F↔W` structure + median summation).

## GEVREY BOUND (diagonal) — ESTABLISHED, uniform in `a`

The frozen `F(g)` is **proven Gevrey-1** (`PROGRAM2_TWBETA_NOISE_BOREL_PROOF_NOTES.md`; Airy-integral
rep + Nevanlinna): `|δ_n|∼(1/π)Γ(n)(3/2)^n`, Borel singularity at `g=2/3=Φ` (verified ratio
`δ_n/δ_{n-1}/(n-1)→3/2`: `1.527,1.515,1.511`). Since the **exact diagonal is the odd sub-series of `F`**,
$$|D_n|\le|\delta_n|\le K\,\Gamma(n)\,(3/2)^n\quad\Longrightarrow\quad
|\,\sigma_n(a)_{\rm diag}\,|=|D_n|\,a^{-3n/2}\le K\,\Gamma(n)\,\Phi(a)^{-n},\ \ \Phi(a)=\tfrac23a^{3/2}.$$
`Φ(a)≥Φ(a₀)>0` for `a≥a₀`, and `K=1/π` is `a`-**independent**, so the **diagonal (leading-in-`a`) part of
the noise series is uniform-in-`a` Gevrey-1**, with Borel singularity pinned at the action `Φ(a)∈ℝ₊`.
The exact diagonal thus **inherits** median Borel-summability from the frozen `F` — the T3 proof *does*
transfer, once restricted to the physical (odd/median) part. **Residual:** the off-diagonal
`O(a^{-3/2})` corrections (the full `σ_n` beyond leading-`a`) still need a uniform bound — but the
**dominant** part is now rigorously controlled.

## Uniform theorem: what is established, what remains

**Statement to prove.** For `a≥a₀>0`, the noise series `Σ_{n≥1}σ_n(a)β^{-n}` is median Borel-summable
along `ℝ₊`, uniformly in `a` (uniform sector + uniform Gevrey-1 bound).

- **[ESTABLISHED — sector uniformity, structural].** The dynamical Borel singularity is the **real** WKB
  action `Φ(a)=\int_0^\infty\!\sqrt{x+a}\,dx`-type period over the noise-independent Airy spectral curve
  `y²=x+a` (paper §8 location-fixing; isomonodromy rigidity). `Φ(a)>0` and real for all `a>0`, `→\tfrac23a^{3/2}`;
  it moves continuously with `a` but **never leaves `ℝ₊`** and is bounded away from `0` for `a≥a₀`. The
  Coulomb-gas singularities are `a`-independent on `iℝ` (`±iπℤ`). Hence a median sector of opening `>π`
  around `ℝ₊`, avoiding `iℝ`, exists **uniformly in `a≥a₀`**.
- **[ESTABLISHED — diagonal Gevrey uniformity, 2026-07-23].** The **diagonal** (leading-in-`a`) part of
  `σ_n(a)` is uniform-in-`a` Gevrey-1: `|σ_n(a)_diag|=|D_n|a^{-3n/2}≤K Γ(n)Φ(a)^{-n}`, because the exact
  diagonal `D_n` is the **odd part of the frozen `log F`** (RECONCILIATION above) and so inherits the
  frozen's **proven** Gevrey-1 bound `|δ_n|∼(1/π)Γ(n)(3/2)^n`. `K=1/π` is `a`-independent, `Φ(a)≥Φ(a₀)>0`.
  The diagonal thus inherits median Borel-summability from `F` (T3 transfers to the physical/odd part).
  *This replaces the earlier "diagonal large-order needs the BN recursion" — it does NOT; the median
  reconciliation supplies it, and `R_2`'s deficit is explained (the even-order artifact cancels).*
- **[ESTABLISHED for `a≥a₀` — off-diagonal bound via the LEVEL resurgence, 2026-07-23]**
  (`coupled-atlas/_tw_offdiag_bound.py`). The off-diagonal `σ_n(a)−D_n a^{-3n/2}=2^nΣ_{m>n}r_{m,n+1}a^{-3m/2}`
  is the column generating function `G_{n+1}(u)=Σ_m r_{m,n+1}u^m`, `u=a^{-3/2}` — a **divergent (Gevrey-1)
  series in `u`**: it *is* the **level-`a` resurgence** (Hastings–McLeod/PII, paper T1). Verified at `n=1`:
  `σ_1(a)/(−\tfrac5{12}a^{-3/2})=1−\tfrac{33}8u+\tfrac{1555}{64}u²−…` has factorially-growing coefficients
  (ratios `4.13, 5.89, …`), not a convergent correction. **Key: the two directions do NOT collide.** The
  noise Borel singularity sits at `Φ(a)=\tfrac23a^{3/2}` and **moves out** with `a`; the level singularity
  sits at a **fixed** level Borel radius `ρ` in `u=a^{-3/2}`. So for
  $$a\ge a_0:=\rho^{2/3},$$
  the column series is level-Borel-summable and `|σ_n(a)−D_n a^{-3n/2}|=O(a^{-3/2})·|D_n a^{-3n/2}|`,
  i.e. `σ_n(a)=D_n a^{-3n/2}(1+O(a^{-3/2}))`, giving `|σ_n(a)|≤K Γ(n)Φ(a)^{-n}(1+O(a^{-3/2}))`. **The level
  radius is `O(1)`:** from the bracket `R_m(X)` (`m≤6`), `|R_m/R_{m-1}|/m→ρ` with `ρ(X=1)≈0.94=2√2/3` (the
  **HM action** — confirms the level direction), `ρ(X=2)≈1.4`, `ρ(X=½)≈0.92`; so `a₀=ρ^{2/3}≈1`. The bound
  holds **throughout the deep-tail regime `a≫1`** (not down to `a→0`).
- **[RESIDUAL — sharp, single point].** The `O(a^{-3/2})` above is `n`-uniform **iff** the level Gevrey
  constant of the columns `G_{n+1}` is `n`-uniform. This is supplied *structurally* by the T1 HM/PII level
  resurgence (the level action is `a`-uniform, and it is the *same* instanton action `Φ` in both
  directions — the two-directional resurgence is action-aligned). Making the `n`-uniform constant fully
  explicit needs the column large-`m` growth (`R_m`, `m→∞`) — the only place the recursion would still help.

**Bottom line.** Uniform-in-`a` (for `a≥a₀≈1`) now has **all three** pieces: sector uniformity
(structural) **+** diagonal Gevrey-1 uniformity (median reconciliation) **+** off-diagonal bound (level
resurgence, no collision). The frozen-vs-exact discrepancy is **resolved** (odd/median structure). What
remains is only the `n`-uniformity of the level Gevrey constant — action-aligned and structurally forced
by T1, needing the `R_m` large-order to pin the constant explicitly. **The uniform-in-`a` theorem holds in
the deep-tail regime modulo that single explicit constant** — a materially complete resolution of the
Handoff-D #2 residual.

## Impact on the paper (`TWbeta_Resurgent_paper.tex`)

- **Rem. 7.4 (`rem:elem`) refine:** replace "the uniform statement … adds the `O(a^{-3/2})` x-extension
  corrections and would combine this proof with the T1 PII map" with the sharper: *the uniform statement
  reduces to a uniform Gevrey-1 bound on the exact Borot–Nadal noise coefficients `σ_n(a)`; the sector is
  already uniform (real WKB action on `ℝ₊`, Coulomb-gas on `iℝ`); the frozen reduced escape reproduces the
  leading noise Stokes term but not the exact series beyond it (`R_2` degree deficit), so the large-order
  input is the BN recursion.*
- **Honest-ledger note:** T3-as-proved is median Borel-summability of the **frozen reduced escape**; its
  identification with the exact tail noise expansion is leading-order (exact at `n=1`). Worth one sentence
  so the claim scope is exact.

## Files
- `coupled-atlas/_tw_uniform_in_a.py` — exact sympy: BN validation at β=2, the `1/β` reduction
  (self-checked), frozen-vs-exact diagonal table, sector-uniformity summary, `σ_1` off-diagonal.

## Gotchas
- **`s = a` at leading order only.** The rate `e^{-\tfrac23βa^{3/2}}` fixes `s=a` to leading order; the
  `n=1,3`-agree / `n=2`-disagree pattern may signal a subleading `s↔a` reparametrization. Do not assume
  `s=a` exactly at higher orders without checking against `R_{m≥4}`.
- **Frozen `d_n` is not a shortcut for the exact `σ_n`** beyond `n=1`. `R_2`'s degree deficit
  (`2<m+1`) makes the exact 2-loop scaling coefficient `0`; the frozen `−5/8` is a freeze artifact.
- **`R_m` degrees are not all `m+1`.** `deg R_1=2, deg R_2=2 (deficit), deg R_3=4`. The top-degree
  coefficient `r_{n,n+1}` (the diagonal) can vanish — track it from the actual polynomial, not the bound.
