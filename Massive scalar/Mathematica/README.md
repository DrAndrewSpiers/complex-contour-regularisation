# Mathematica implementation: complex-contour regularisation of the massive scalar toy model

An independent test of the method of `Claude Implementation/PRD/prd.tex` on the screened Poisson (massive scalar) toy model,
(∇² − m²)Φ₁ = −4πq δ³(x − x₀), (∇² − m²)Φ₂ = ½∂_Z²(Φ₁²).
Nothing is imported from the massless folder; `ContourMassive.wl` is a self-contained copy-and-adaptation of `ContourSecondOrder.wl`
plus the pieces of the massless stage scripts (d_Z ladder, Cauchy continuation, Lemma-1 error decomposition) that were needed.

Everything was run with Mathematica 14.3 (wolframscript) on 25 October 2026 on a 12-core Apple-silicon Mac.

All 1129 checks pass in the final runs (25 October 2026): 406 (Stage 1), 192 (Stage 2), 328 (Stage 3) and 203 (Stage 4); total run time 77 min.

**Headline results.** (T1) With punctures truncated at ρⁿ the hybrid-split error falls as ℓ′max^{−(n+4)} (fitted slopes 4.87, 6.82, 8.78 for n = 1, 3, 5, the same at m = 1/5 and m = 1); the residual is C^{n+1}, two orders smoother than the massless cqρⁿ trick, whose slopes n + 2 are reproduced at m = 0. Exact puncture: SS and SR recovered to 5e-18 (m = 1/5) and 9e-17 (m = 1). (T2) The complex-radius angular quadrature of the massive SS source works with the same node count (80 per segment for 1e-14) once the two branches are blended before the d_Z² ladder; one branch alone is wrong at O(1). (T3) Conditions (i) and (ii) behave exactly as in the massless case with Bessel functions ((ii) fails by a purely imaginary (i/π)E₀ with E₀(r₀) = −1/r₀); the SR δ coefficient c_ℓ is m-independent and a new SS δ term 4πmq²Y_ℓ(0)/r₀² appears; no spurious worldline terms (fits ≤ 5e-16). Shadow growth per coupled mode 1.129 (R = 1/2), 1.062 (R = 1/4), unchanged.

## Contents

| File | Purpose |
|---|---|
| `ContourMassive.wl` | Package: parameters (adds the mass `m`), the massive radial operator and its contour form, the homogeneous solutions `uFun`, `vFun` (modified spherical Bessel functions normalised to the massless powers; m = 0 falls back to r^l, r^{-l-1}), the physical Green's function and Robin conditions, the massive Coulomb mode and its Θ→H continuation, closed-form modes of ρ^k and e^{-cρ}ρ^p by the ρ-substitution (exponential integrals, with the symbol `LOG` for the branch of the logarithm), the d_Z ladder and δ coefficients, the 3D sources with their r-derivatives, the complex-radius angular quadrature with the residue-corrected blend (`blendedModes`, `quadModeDzz`), the Cauchy continuation, the solvers (three-region, physical conditions) and the exact-reference builder (`chiGreen`, `exactReference`). |
| `massive_stage1_first_order.wls [smoke]` | Stage 1: first order. Bessel identities; the continued singular mode and the consistency conditions (i), (ii) of PRD App. C with Bessel functions; the corrected mode; DSolve; first-order residual solves on the contour with the Θ→H, corrected and truncated punctures; smoothness of the residual modes; m = 0 against the massless closed form. |
| `massive_stage2_reference.wls [smoke]` | Stage 2: the exact second-order reference, mode by mode, with four independent cross-checks; the δ coefficients of the SS and SR sources; writes `results/massive_reference_m*.mx` used by Stages 3 and 4. |
| `massive_stage3_contour.wls [smoke]` | Stage 3: the SS and SR sources on the contour by complex-radius angular quadrature (T2) against the closed-form continuation, node count, single-branch control; exact-puncture three-region solves with the physical conditions, point-supported fits with a positive control, the pointwise-SR negative control (T3 (iv)); noise; m = 0 against the massless Stage 9; figures. |
| `massive_stage4_convergence.wls [smoke]` | Stage 4: truncated punctures Φᴾ_n (T1): hybrid split H1 (k = 0), H0 and the product route against ℓ′max ≤ 40 for n = 1, 3, 5 and m = 1/5, 1; slopes and per-mode powers; tail acceleration; the Cauchy step by actual contour solves; shadow growth per coupled mode; m = 0 with the massless trick puncture against the massless Stage 12b. |
| `run_all.sh [smoke]` | Runs the four stages and builds the note; logs in `../results/logs/`. |

The `smoke` argument runs a reduced version (fewer ℓ, n, ℓ′max) in about ten minutes.

## Conventions

- Modes Φ = Σ φ_ℓ(r) Y_ℓ0 with orthonormal Y_ℓ0, [f]_ℓ = √(π(2ℓ+1)) ∫ f P_ℓ dμ.
- L_ℓ ψ = ψ″ + 2ψ′/r − [ℓ(ℓ+1)/r² + m²] ψ. Homogeneous solutions u_ℓ(r) = (2ℓ+1)!! i_ℓ(mr)/m^ℓ → r^ℓ and v_ℓ(r) = m^{ℓ+1} k_ℓ(mr)/(2ℓ−1)!! → r^{−ℓ−1}, with Wronskian u v′ − u′ v = −(2ℓ+1)/r². Green's function with the physical conditions G_ℓ(r, s) = −(s²/(2ℓ+1)) u_ℓ(r_<) v_ℓ(r_>).
- First order: [q e^{−mρ}/ρ]_L = q √(4π/(2L+1)) u_L(r_<) v_L(r_>) (the massive Coulomb mode; its derivative jump at r₀ is m-independent). κ_ℓ = −(q/r₀²)√(4π(2ℓ+1)) as in the massless test.
- Regular field Φᴿ = c₀ u₀(r)/u₀(r₀) + c₁ r₀ [u₁(r)/u₁(r₀)] cos θ, c₀ = 1, c₁ = 1/10, so that Φᴿ(x₀) = 3/2 for every m and Φᴿ → c₀ + c₁ r cos θ as m → 0. Note that Φᴿ grows like e^{mr}/r: at m = 1 the RR piece of Φ₂ reaches 2.8 × 10³ at r = 10 and dominates the total.
- Truncated punctures Φᴾ_n = q Σ_{k=−1}^{n} (−m)^{k+1} ρ^k/(k+1)!, used everywhere without a window function (as the massless cqρⁿ puncture was). At m = 0 the "massless trick" puncture q/ρ + cqρⁿ (c = r₀^{−(n+1)}) is used instead, to compare with the massless Stage 12b.
- Parameters: q = 1, r₀ = 5, domain [1, 10], worldtube [4.5, 5.5], contour y(a) = a + i(1 − 4(a − 5)²)⁶, working precision 30 (closed forms at 50, references at 40).
- H(z) = 1 + (i/π) log(z − r₀), δ_ext = i/(π(z − r₀)); real-axis values are boundary values from above.
- Closed-form modes carry the symbol `LOG` for log|r − r₀|; `logReal` gives ½ log (z − r₀)² (real, differentiable), `logOut` log(z − r₀), `logIn` log(z − r₀) − iπ (the canonical in branch). **A Θ→H blend in + (out − in)H must use `logOut` on both branches** (so that out − in is real-analytic on the axis); the canonical in branch differs from that by the pure shadow iπA, A the coefficient of the log.
- Exterior grid: 27 points in [1, 3.5] ∪ [6.5, 10]. **Scaled error** of a piece = max|numerical − reference|/max|reference of that piece| over the grid. In Stage 4 the scale is max|SS + SR reference| (the contour pieces); the massless convention max|total| is used only for the m = 0 comparison, because at m = 1 the total is dominated by RR.
- Physical boundary conditions throughout: ψ′ = ψ u_ℓ′/u_ℓ at r = 1 and ψ′ = ψ v_ℓ′/v_ℓ at r = 10, imposed on the contour solves by adding α u_ℓ to the zero-data solution.

## How to run

```
W=/Applications/Wolfram.app/Contents/MacOS/wolframscript
cd "<this folder>"
$W -file massive_stage1_first_order.wls      #  8.1 min
$W -file massive_stage2_reference.wls        #  3.2 min (before Stages 3 and 4)
$W -file massive_stage3_contour.wls          # 39.5 min
$W -file massive_stage4_convergence.wls      # 26.4 min
zsh run_all.sh                               # all of the above and the note
```

Each script prints `PASS`/`FAIL` lines against stated tolerances and ends with a `SUMMARY`; the checks are written to `../results/massive_stageN_checks.csv`, the tables to `../results/massive_stageN_summary.md` and CSVs, the figures (vector PDF + CSV) to `../figures/`.

## Results

### Stage 1: first order (`massive_stage1_first_order.wls`, 406 checks pass, 8.1 min)

**Bessel identities.** u_ℓ v_ℓ′ − u_ℓ′ v_ℓ = −(2ℓ+1)/r² to 1e-54 (ℓ ≤ 12); L_ℓ u_ℓ = L_ℓ v_ℓ = 0; the Green's function has a unit derivative jump and meets both Robin conditions (1e-25); [e^{−mρ}/ρ]_ℓ by angular quadrature, by the Bessel form and by the ρ-substitution closed form agree to 1e-41; Φᴿ(x₀) = 3/2 and Φᴿ is a homogeneous solution.

**The continued singular mode (T3 (i), (ii)).** The Θ→H blend ψ̃ = Gⁱⁿ + (Gᵒᵘᵗ − Gⁱⁿ)H of the massive Coulomb mode (`results/massive_stage1_singular_mode.csv`):

| m | (i) max abs(Re ψ̃ − G_ℓ), ℓ ≤ 6 | (ii) abs(Re(Lψ̃ − δ_ext)) on the axis | abs(Im(Lψ̃ − δ_ext)) on the axis, ℓ = 0 … 6 | E₀(r₀) | corrected mode: ODE residual of P̃ (rel.) / return to the real P at a_R (rel.) |
|---|---|---|---|---|---|
| 1/5 | < 1e-76 | < 1e-36 | 0.54, 6.5, 40, 2.2e2, 1.1e3, 5.7e3, 2.9e4 | −0.2000 | < 1e-8 / < 1e-12 |
| 1 | < 1e-85 | < 1e-35 | 19, 48, 1.7e2, 6.9e2, 3.0e3, 1.3e4, 6.1e4 | −0.2000 | < 1e-8 / < 1e-12 |
| 0 | 0 | < 1e-36 | 0.32, 5.8, 37, 2.0e2, 1.1e3, 5.5e3, 2.8e4 | −0.2000 | < 1e-8 / < 1e-12 |

- (i) holds and (ii) fails exactly as in the massless case: Lψ̃ − δ_ext = (i/π)E₀ with E₀ = (2J′ − 1)/(z − r₀) − J/(z − r₀)² + 2J/(z(z − r₀)), J = Gᵒᵘᵗ − Gⁱⁿ, which is real-analytic on the axis and has E₀(r₀) = (3/2)J″(r₀) + 2/r₀ = −1/r₀ for every ℓ and m (J″(r₀) = −2/r₀ follows from the Wronskian). At m = 0, ℓ = 0 it reduces to E₀ = −1/z (PRD App. C). |E₀| grows with ℓ away from r₀ like (r/r₀)^{±ℓ}.
- The corrected mode ψ̃ − (i/π)P with L P = E₀ (P real, solved from zero data at r = 1 and continued along the contour by the same ODE) satisfies both conditions; P̃ returns to the real P at a_R to 1e-12 relative (NDSolve limit at ℓ = 6). At m = 0, ℓ = 0 it has the real part of the massless closed-form φˢ₀ = ψ_p − D₀ exactly, and the imaginary parts differ by a homogeneous solution (< 1e-29).
- **DSolve** returns closed forms (Ei for ℓ = 0; CoshIntegral/SinhIntegral × spherical Bessel for ℓ = 2) in about 1 s. Evaluated from above, Re ψ_p − G_ℓ is a homogeneous solution α u_ℓ + β v_ℓ (fit residual < 1e-38):

| m | ℓ | α (u_ℓ) | β (v_ℓ) | massless (α, β) |
|---|---|---|---|---|
| 1/5 | 0 | 1.84 | 34.0 | (5, 25) |
| 1 | 0 | 0.0337 | 371 | (5, 25) |
| 1/5 | 2 | 0.0343 | < 1e-36 | (0.04, 0) |
| 1 | 2 | 0.00386 | < 1e-40 | (0.04, 0) |

  The massless pattern persists (anti-physical Green's function at ℓ = 0, one-sided at ℓ ≥ 2), so DSolve's particular solution needs the offset correction for the massive operator too.

**Residual solves on the contour** (`results/massive_stage1_residual_solves.csv`; test regular field u_ℓ(r)/u_ℓ(r₀) for every ℓ; S_eff = κ_ℓ(δ_ext − L φᴾ)):

| puncture | exterior φᴿ error (scaled, max over ℓ ≤ 6, m) | contour residual vs closed form | residual shadow on the contour |
|---|---|---|---|
| Θ→H exact | 4.6e-18 | 1.1e-19 (closed form φ̃ᴿ − (i/π)κP̃) | 0.12–1.0 (O(1): (ii) violated) |
| corrected | 4.6e-18 | 2.9e-18 (φ̃ᴿ) | < 3e-18 |
| truncated n = 1 | 4.6e-18 | 2.9e-18 (φ̃ᴿ + κ(f̃ − φ̃ᴾ_n)) | 2.4e-3–3.3e-2 (m = 1/5); 2.4e-3–1.4e-1 (m = 1) |
| truncated n = 3 | 4.6e-18 | 2.9e-18 | 3.4e-7–9.0e-3 (m = 1/5); 2.4e-3–3.5e-1 (m = 1) |
| truncated n = 5 | 1.0e-17 | 2.9e-18 | 1.5e-8–1.0e-3 (m = 1/5); 2.4e-3–4.3e-1 (m = 1) |

- The exterior field is exact for every puncture; only the shadow of the residual on the contour depends on how the singular field is continued, as Lemma 1 says.
- With m = 1 the shadow of the truncated-puncture residual *grows* with n at low ℓ. This is because the puncture is used without a window function: its polynomial tail q(−m)^{n+1}ρⁿ/(n+1)! is large at ρ ~ 2r₀, so the residual modes [e^{−mρ}/ρ − Φᴾ_n/q]_ℓ are large everywhere and their Θ→H continuation has a large imaginary part on the contour. It does not affect the exterior (hybrid route), only the product route (Stage 4).

**Smoothness of the residual** (`results/massive_stage1_smoothness.csv`): for n = 1, 3, 5 and ℓ′ = 0, 2, 5 the derivatives of order 0 … n+3 of [Φ₁ − Φᴾ_n]_ℓ′ are continuous at r₀ (< 1e-80) and the (n+4)th jumps: the residual modes are C^{n+3}, the 3D residual C^{n+1} (the first omitted term ρ^{n+1} is an even power, hence a polynomial; the first non-smooth term is ρ^{n+2}). The expected convergence is therefore ℓ′max^{−(n+4)} (Stage 4), against ℓ′max^{−(n+2)} for the massless cqρⁿ trick.

### Stage 2: the exact reference (`massive_stage2_reference.wls`, 192 checks pass, 3.2 min)

Φ₂ = ½∂_Z²χ with (∇² − m²)χ = Φ₁²; χ is not closed form in 3D (only its SS part is, Eq. (chiSS) of the note: q²[(e^{mρ}Ei(−3mρ) − e^{−mρ}Ei(−mρ))/(2mρ) − ln3 cosh(mρ)/(2mρ)] = q² ln ρ + O(1), no 1/ρ term). The reference is built mode by mode: closed-form modes of Φ₁² (SS by the ρ-substitution with exponential integrals, SR by the Bessel couplings, RR by Gaunt couplings), Green's-function quadrature χ_j = ∫ G_j [Φ₁²]_j ds (NIntegrate, double-exponential rule, breakpoint at r₀, 40 digits), χ_j″ and χ_j‴ from the ODE, the d_Z² ladder, and the homogeneous correction A u_ℓ + B v_ℓ that meets the Robin conditions. Written to `results/massive_reference_m{1_5,1,0}.mx` and `.csv` (30 digits).

| check | m = 1/5 | m = 1 | m = 0 |
|---|---|---|---|
| closed-form source modes (SS, SR, RR) vs angular quadrature, j ≤ 8, both sides of r₀ | < 1e-22 | < 1e-22 | < 1e-22 |
| SS δ coefficient = 4πmq²Y_ℓ(0)/r₀², SR δ coefficient = −4πqΦᴿ(x₀)Y_ℓ(0)/r₀² | < 1e-25 | < 1e-25 | < 1e-25 |
| A: χ₂^SS by two NIntegrate rules (DE vs Gauss–Kronrod/IMT) | < 1e-40 | < 1e-40 | < 1e-40 |
| B: SS reference vs modes of ½∂_Z²χ_SS(ρ) (3D closed form; ℓ = 0, 1, 2, 4, 6) | < 1e-25 | < 1e-28 | – |
| C: SR / RR reference vs NDSolve real-line solves with the Robin conditions | 4.2e-18 / 2.5e-19 | 6.2e-17 / 6.9e-21 | 3.9e-18 / 1.4e-20 |
| D: m = 0 reference vs the massless closed form χ (PRD App. A) | – | – | < 1e-34 |
| Robin conditions met by the reference at r = 1 and 10 | < 1e-25 | < 1e-25 | < 1e-25 |

- The reference is accurate to 1e-18 or better (the C numbers are the 30-digit NDSolve); the SS piece to 1e-25 or better.
- Size of the reference modes (max over the exterior grid): m = 1/5: total 2.2 (ℓ = 0) to 3.5e-2 (ℓ = 6), SS 2–5e-2, SR 2e-2–2e-1, RR 1.9 (ℓ = 0) to 0 (ℓ ≥ 5). m = 1: total 2.8e3 (ℓ = 0) to 2.5e-2 (ℓ = 6), SS 8e-3–1.4e-2, SR 1.5e-2–4.3e-2, RR 2.8e3 (ℓ = 0) to 0 (ℓ ≥ 5): at m = 1 the RR piece (Φᴿ² ~ e^{2mr}) dominates the total by five orders of magnitude, hence the per-piece scaling. Figure `fig_reference_modes.pdf`.
- **Abandoned:** a continuity check at m = 1/1000 and at m = 1/20. The closed forms with e^{−cρ} cancel like 1/c^{2j+p+2} (up to 1e51 at m = 1/1000) and NIntegrate hung; the m → 0 limit is checked exactly at m = 0 (row D), where the package uses the massless forms.
- Numerical note: the term Ei(−2m(r+r₀)) of the SS closed form sits exactly on Mathematica's branch cut for real r and jumps by −iπ for Im z > 0; it is written through the entire part Ei(w) − γ − log w (`eEnt`, with its derivatives `eEntK[k]`).


### Stage 4: truncated punctures and the ℓ′ convergence (`massive_stage4_convergence.wls`, 203 checks pass, 26.4 min)

Puncture Φᴾ_n = q Σ_{k≤n} (−m)^{k+1}ρ^k/(k+1)!, n = 1, 3, 5, used consistently (SS_n closed form on the contour, SR_n by the hybrid split with k = 0, RR_n on the real line); m = 1/5 and 1; ℓ ≤ 6, ℓ′max ≤ 40. Every error is the Green's-function response (physical conditions) of its error source, as in the massless Stage 12b, and is scaled by max|SS + SR reference| (the massless convention max|total| only for the m = 0 row). Full tables: `results/massive_stage4_summary.md`; CSVs `massive_stage4_{convergence,slopes,permode,routes,acceleration,cauchy_step,growth,delta,profile}.csv`; figure `figures/fig_convergence.pdf` (+ CSV).

| m | n | slope (ℓ′max 20–40) | expected n+4 | per-mode power | expected n+5 | error at ℓ′max = 12 | 20 | 40 | at 40, scaled by max abs(total) |
|---|---|---|---|---|---|---|---|---|---|
| 1/5 | 1 | 4.87 | 5 | 6.01 | 6 | 1.7e-6 | 1.3e-7 | 4.3e-9 | 4.0e-9 |
| 1/5 | 3 | 6.82 | 7 | 8.02 | 8 | 6.6e-9 | 1.8e-10 | 1.6e-12 | 1.5e-12 |
| 1/5 | 5 | 8.78 | 9 | 10.1 | 10 | 2.9e-11 | 3.0e-13 | 6.7e-16 | 6.3e-16 |
| 1 | 1 | 4.87 | 5 | 5.98 | 6 | 9.5e-4 | 8.4e-5 | 2.9e-6 | 2.9e-6 |
| 1 | 3 | 6.82 | 7 | 8.01 | 8 | 9.0e-5 | 3.0e-6 | 2.7e-8 | 2.7e-8 |
| 1 | 5 | 8.77 | 9 | 10.0 | 10 | 9.7e-6 | 1.2e-7 | 2.8e-10 | 2.8e-10 |
| 0 (trick q/ρ + cqρ³) | 3 | 4.87 | 5 (n+2) | 6.01 | 6 (n+3) | 3.90e-5 | 3.27e-6 | 1.11e-7 | 1.11e-7 |

- **Rates.** The hybrid-split error falls as ℓ′max^{−(n+4)} and the per-mode contributions as ℓ′^{−(n+5)}: the residual Φ₁ − Φᴾ_n is C^{n+1} (Stage 1), two orders smoother than the massless cqρⁿ residual, and the fitted slopes fall short of the integers by the same 0.1–0.2 as the massless fits. The slopes are the same at m = 1/5 and 1; only the prefactor changes, by m^{n+3}/(n+3)! (6e2 at n = 1 to 4e5 at n = 5 between the two masses).
- **m = 0 check.** The massless trick puncture through the same code reproduces the massless Stage 12b (physical conditions): slope 4.87 (12b: 4.87), per-mode power 6.01 (6.01), errors 3.90e-5, 3.27e-6, 1.11e-7 at ℓ′max = 12, 20, 40 (identical to the three digits published there; the check tolerance is the 3-digit rounding).
- **Routes** (max over ℓ, scaled; `massive_stage4_routes.csv`):

| m | n | route | ℓ′max = 4 | 8 | 12 | 16 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|
| 1/5 | 3 | H1 = (c) | 1.2e-4 | 3.2e-7 | 6.6e-9 | 8.6e-10 | 1.8e-10 | 1.2e-11 | 1.6e-12 |
| 1/5 | 3 | product | 1.2e-4 | 3.2e-7 | 5.0e-9 | 3.0e-9 | 3.5e-9 | – | – |
| 1/5 | 3 | H0 | 2.76 | 2.76 | 2.76 | 2.76 | 2.76 | 2.76 | 2.76 |
| 1/5 | 3 | H1 + 1-term Λ tail | – | – | 1.9e-9 | 1.0e-10 | 1.1e-11 | 2.1e-13 | 1.4e-14 |
| 1 | 3 | H1 = (c) | 9.3 | 3.8e-2 | 9.0e-5 | 1.4e-5 | 3.0e-6 | 1.9e-7 | 2.7e-8 |
| 1 | 3 | product | 9.3 | 3.8e-2 | 6.8e-5 | 1.4e-5 | 2.0e-5 | – | – |
| 1 | 3 | H0 | 9.3 | 0.91 | 0.91 | 0.91 | 0.91 | 0.91 | 0.91 |
| 1 | 3 | H1 + tail | – | – | 2.2e-5 | 1.3e-6 | 8.9e-8 | 2.1e-9 | 1.5e-10 |
| 1 | 5 | H1 = (c) | 39 | 4.7e-2 | 9.7e-6 | 8.5e-7 | 1.2e-7 | 3.5e-9 | 2.8e-10 |
| 1 | 5 | product | 39 | 4.7e-2 | 9.8e-6 | 9.5e-7 | 2.2e-7 | – | – |

  H1 equals (c) for every n and ℓ′max (Part 1 carries the whole SR δ: |g_SR − g_P1| < 1e-46 relative; the truncated RR has no δ). H0 misses c_ℓ δ and fails at O(1). The product route follows (c) and turns up beyond ℓ′max ≈ 16 (its shadow term grows with ℓ′). At small ℓ′max the error is large at m = 1 (39 at ℓ′max = 4 for n = 5) because the windowless puncture's polynomial tail makes the low-ℓ′ residual modes large.
- **Shadow growth per coupled mode** (product route, Q = AB-type response of −2 Im K̃ᴾ Im φ̃^res at a_R, ℓ = 2, n = 3; `massive_stage4_growth.csv`): 1.129 (R = 1/2) and 1.062 (R = 1/4) at both m = 1/5 and m = 1, against 1.132 / 1.062 with the massless trick puncture and the PRD's 1.12 / 1.06: still ≈ 1 + R/r₀ (1.10, 1.05) times a slowly varying power. The Bessel modes grow across the tube like the powers for ℓ′ ≫ mr.
- **Tail acceleration** (Λ fits with the leading power n + 5): 2-term gain 16–29 at ℓ′max = 20, 160–690 at 30, 3e3–1e4 at 40; 1-term gain 14–38 at 20 and 3–5 at 12 (where the 2-term fit makes things worse, as in the massless case).
- **Cauchy step** (actual contour solves of the H1 and H0 sources against the distributional line, n = 3, ℓ′max = 8, 20, ℓ = 2, 6; `massive_stage4_cauchy_step.csv`): contour − line = predicted δ miss to 1e-18 (m = 1/5) and 3e-15 (m = 1) relative to |c_ℓ G|; propagated with the physical conditions it contributes ≤ 3.5e-17 (m = 1/5) and ≤ 2.4e-15 (m = 1) to the scaled exterior error, ≤ 2e-7 of the H1 error. The m = 1 value is larger than the massless 1.3e-16 because the Bessel functions inside NDSolve are evaluated at 30 digits.
- **Error against r**: largest next to the tube, decaying outwards on both sides (m = 1, n = 3, ℓ′max = 20: 1.3e-5 in [3.5, 4.5], 5.0e-6 in [5.5, 6.5], 2.1e-6 in [6.5, 10]).
- **Numerical note.** The closed-form modes [ρ^k]_L for L up to 48 cancel about 100 digits at ρ = r + r₀ (the puncture has no window), so they are evaluated from 200-digit input; at 70 digits the ℓ ≥ 3 checks collapsed to "0 ± 1e98".

### Stage 3: the sources on the contour and the exact-puncture field (`massive_stage3_contour.wls`, TIME3)

Exact puncture q e^{−mρ}/ρ; m = 1/5, 1, 0; ℓ ≤ 6; 128 contour nodes (8 Gauss panels × 16) for the contour Green's-function quadrature, real-axis Gauss panels (20 nodes) with breakpoints at the grid points and at 4, 4.25, 5.75, 6. Tables: `results/massive_stage3_summary.md`; CSVs `massive_stage3_{ss_quadrature,nodes,sr_quadrature,solves,noise}.csv`; figures `fig_contour_source.pdf` (contour; ℓ = 2 sources on the real axis and on the contour), `fig_source_modes.pdf`, `fig_ss_nodes.pdf`, `fig_recovered_modes.pdf`, each with a CSV.

**(A) SS by complex-radius angular quadrature (T2).** The ρ-integrand is e^{−2mρ} × (Laurent polynomial), with the single fixed pole at ρ = 0; the two-segment path via ∓i/2 is unchanged. New: the two branches differ by a real-analytic kink part as well as by iπ × residue, so the undifferentiated modes are blended, F_c + (F_out − F_c)H with F_c = M₋ + iπA (A = residue at ρ = 0 by a 64-point circle rule, equal to the symbolic coefficient to 1e-30), before the d_Z² ladder. Max over ℓ ≤ 6:

| m | pointwise rel. diff. from the closed-form continuation (80 nodes/segment) | exterior SS, quadrature source | exterior SS, closed-form source | one branch alone (no blend): pointwise / exterior |
|---|---|---|---|---|
| 1/5 | 6.5e-15 | 2.0e-15 | 1.7e-15 | 7–9e-2 / 0.3–1.1 |
| 1 | 7.7e-15 | 6.4e-16 | 6.3e-16 | 0.48–0.63 / 1.8–2.7 |
| 0 | 6.3e-15 | 4.1e-15 | 3.3e-15 | 6e-16 / 5.5e-15 |

Node count (pointwise error vs closed form, max over ℓ; the exterior field saturates at the 1e-15 of the 128-node contour quadrature):

| nodes per segment | 32 | 48 | 64 | 80 | 120 |
|---|---|---|---|---|---|
| 3D-source evaluations per contour node | 64 | 96 | 128 | 160 | 240 |
| pointwise, m = 1/5 | 3.0e-4 | 1.3e-7 | 3.3e-11 | 6.5e-15 | 2.1e-24 |
| pointwise, m = 1 | 3.5e-4 | 1.5e-7 | 3.8e-11 | 7.7e-15 | 2.5e-24 |
| pointwise, m = 0 | 2.9e-4 | 1.3e-7 | 3.2e-11 | 6.3e-15 | 2.0e-24 |
| exterior, m = 1/5 | 4.4e-4 | 9.6e-8 | 1.0e-11 | 2.0e-15 | 1.7e-15 |

- The accuracy and node count are those of the massless case (PRD: 1.8e-14 at 80 nodes; 64 for 1e-8): the exponential costs nothing.
- One branch alone, which sufficed for the massless SS, now misses the kinks (and the δ term) and is wrong at O(1) in the exterior; at m = 0 the same code gives 5.5e-15. **Blending with the canonical in branch (log on the other side of the cut) instead of F_c is wrong by O(1)** — the one new branch rule of the massive model (`logIn` vs `logOut`, see Conventions).

**(B) SR by the same quadrature** (2q e^{−mρ}Φᴿ/ρ, entire in ρ so A = 0): pointwise 4.5e-16–7.7e-16 vs the closed form of the Bessel couplings; exterior SR 2.3e-17–1.4e-16 (quadrature source) and 1.4e-17–7.6e-17 (closed-form source), all m.

**(D) Noise** (1e-8 relative Gaussian noise on the SS samples; noise-induced exterior error only; 3 realisations, max):

| m | contour, 64 nodes | real axis, 112 samples graded towards r₀ (ratio 1/5, down to 3e-5) | ratio |
|---|---|---|---|
| 1/5 | 5.1e-9 – 2.5e-8 | 1.7e-2 – 7.3e-2 | 1.5e6 – 4.8e6 |
| 1 | 4.6e-9 – 1.3e-8 | 1.3e-2 – 3.7e-2 | 1.1e6 – 5.1e6 |
| 0 | 4.8e-9 – 4.4e-8 | 1.6e-2 – 8.6e-2 | 1.5e6 – 4.6e6 |

The contour loses a factor 0.5–2.5 on the noise, the graded real grid 1e6 (it samples where the source is ~A/Δr² ~ 1e9 times larger); the massless 1e2–3e5 used coarser grids.

**(C) Exact-puncture accuracy: three-region NDSolve solves with the physical conditions** (closed-form contour sources, tabulated on piecewise Chebyshev interpolants: 40 contour pieces and 2 × 16 real pieces of 24 nodes; RR on the real line; `massive_stage3_solves.csv`, figure `fig_recovered_modes.pdf`). Scaled per piece, max over ℓ ≤ 6:

| m | SS | SR | RR | total | SS + SR (by max abs(SS + SR)) | SS fit abs(α)+abs(β) | SR fit | pointwise SR, no δ | + c_ℓ G(r, r₀) | Im SS(a_R), Im SR(a_R) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1/5 | 5.1e-18 | 5.4e-18 | 2.1e-20 | 2.3e-18 | 2.6e-18 | 8.3e-18 | 2.6e-18 | 1.9–4.3 | 1.5e-17 | 0.11–0.34, 0.07–0.31 |
| 1 | 4.3e-17 | 8.7e-17 | 2.3e-20 | 7.0e-17 | 7.0e-17 | 1.5e-16 | 4.8e-16 | 1.1–1.3 | 1.0e-16 | 0.08–0.21, 0.09–0.37 |
| 0 | 5.7e-18 | 4.6e-18 | 1.7e-20 | 2.5e-18 | 2.5e-18 | 9.4e-18 | 2.9e-18 | 1.8–6.3 | 1.2e-17 | 0.12–0.38, 0.07–0.31 |

- The exterior field is recovered to 17 digits with Bessel homogeneous solutions; the m = 1 figures are the 30-digit NDSolve on Bessel functions. The massless case through the same code gives 6e-18 (massless Stage 9: ≤ 3e-18 with physical conditions).
- **Point-supported fit (T3 (iv)).** With the physical conditions a point source at x₀ responds with v_ℓ(r) outside and u_ℓ(r) inside, so the fit uses one point in each region (r = 9/4 and 8) with the two-sided physical Coulomb mode and its ∂_{r₀} derivative. |α| + |β| ≤ 4.8e-16; the positive control (unit δ³ added to the SS contour source) gives α = −0.0795774715459 = −1/(4π) to 12 digits and β < 1e-17 for every m. (The first attempt used two outside points, where the two basis functions are both ∝ v_ℓ and the fit is ill-posed.)
- **Negative control.** The pointwise real-line SR without its δ term is wrong by 1.1–6.3; adding c_ℓ G(r, r₀) restores it to 1e-16. c_ℓ = −0.2127, −0.3684, −0.4756, −0.5627, −0.6381, −0.7054, −0.7669 (ℓ = 0 … 6) at every m, the massless values.
- **m = 0 against the massless Stage 9 (Test 1):** SS and SR errors ≤ 6e-18 (massless Cauchy-data table: ≤ 3.3e-16; PRL physical conditions: < 3e-18), c_ℓ identical to 1e-11.
- Run time 39.5 min, of which the two quadrature parts take 15 min and the 42 three-region solves 20 min.

### Stage 2: the exact reference (`massive_stage2_reference.wls`, 192 checks pass, 3.2 min)

Φ₂ = ½∂_Z²χ with (∇² − m²)χ = Φ₁²; χ is not closed form in 3D (only its SS part is, Eq. (chiSS) of the note: q²[(e^{mρ}Ei(−3mρ) − e^{−mρ}Ei(−mρ))/(2mρ) − ln3 cosh(mρ)/(2mρ)] = q² ln ρ + O(1), no 1/ρ term). The reference is built mode by mode: closed-form modes of Φ₁² (SS by the ρ-substitution with exponential integrals, SR by the Bessel couplings, RR by Gaunt couplings), Green's-function quadrature χ_j = ∫ G_j [Φ₁²]_j ds (NIntegrate, double-exponential rule, breakpoint at r₀, 40 digits), χ_j″ and χ_j‴ from the ODE, the d_Z² ladder, and the homogeneous correction A u_ℓ + B v_ℓ that meets the Robin conditions. Written to `results/massive_reference_m{1_5,1,0}.mx` and `.csv` (30 digits).

| check | m = 1/5 | m = 1 | m = 0 |
|---|---|---|---|
| closed-form source modes (SS, SR, RR) vs angular quadrature, j ≤ 8, both sides of r₀ | < 1e-22 | < 1e-22 | < 1e-22 |
| SS δ coefficient = 4πmq²Y_ℓ(0)/r₀², SR δ coefficient = −4πqΦᴿ(x₀)Y_ℓ(0)/r₀² | < 1e-25 | < 1e-25 | < 1e-25 |
| A: χ₂^SS by two NIntegrate rules (DE vs Gauss–Kronrod/IMT) | < 1e-40 | < 1e-40 | < 1e-40 |
| B: SS reference vs modes of ½∂_Z²χ_SS(ρ) (3D closed form; ℓ = 0, 1, 2, 4, 6) | < 1e-25 | < 1e-28 | – |
| C: SR / RR reference vs NDSolve real-line solves with the Robin conditions | 4.2e-18 / 2.5e-19 | 6.2e-17 / 6.9e-21 | 3.9e-18 / 1.4e-20 |
| D: m = 0 reference vs the massless closed form χ (PRD App. A) | – | – | < 1e-34 |
| Robin conditions met by the reference at r = 1 and 10 | < 1e-25 | < 1e-25 | < 1e-25 |

- The reference is accurate to 1e-18 or better (the C numbers are the 30-digit NDSolve); the SS piece to 1e-25 or better.
- Size of the reference modes (max over the exterior grid): m = 1/5: total 2.2 (ℓ = 0) to 3.5e-2 (ℓ = 6), SS 2–5e-2, SR 2e-2–2e-1, RR 1.9 (ℓ = 0) to 0 (ℓ ≥ 5). m = 1: total 2.8e3 (ℓ = 0) to 2.5e-2 (ℓ = 6), SS 8e-3–1.4e-2, SR 1.5e-2–4.3e-2, RR 2.8e3 (ℓ = 0) to 0 (ℓ ≥ 5): at m = 1 the RR piece (Φᴿ² ~ e^{2mr}) dominates the total by five orders of magnitude, hence the per-piece scaling. Figure `fig_reference_modes.pdf`.
- **Abandoned:** a continuity check at m = 1/1000 and at m = 1/20. The closed forms with e^{−cρ} cancel like 1/c^{2j+p+2} (up to 1e51 at m = 1/1000) and NIntegrate hung; the m → 0 limit is checked exactly at m = 0 (row D), where the package uses the massless forms.
- Numerical note: the term Ei(−2m(r+r₀)) of the SS closed form sits exactly on Mathematica's branch cut for real r and jumps by −iπ for Im z > 0; it is written through the entire part Ei(w) − γ − log w (`eEnt`, with its derivatives `eEntK[k]`).


