# Mathematica implementation: complex-contour regularisation of the massive scalar toy model

An independent test of the method of `Claude Implementation/PRD/prd.tex` on the screened Poisson (massive scalar) toy model,
(∇² − m²)Φ₁ = −4πq δ³(x − x₀), (∇² − m²)Φ₂ = ½∂_Z²(Φ₁²).
Nothing is imported from the massless folder; `ContourMassive.wl` is a self-contained copy-and-adaptation of `ContourSecondOrder.wl`
plus the pieces of the massless stage scripts (d_Z ladder, Cauchy continuation, Lemma-1 error decomposition) that were needed.

Everything was run with Mathematica 14.3 (wolframscript) on 25 October 2026 on a 12-core Apple-silicon Mac.

RESULTS-HEADER

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
$W -file massive_stage1_first_order.wls      # TIME1
$W -file massive_stage2_reference.wls        # TIME2 (before Stages 3 and 4)
$W -file massive_stage3_contour.wls          # TIME3
$W -file massive_stage4_convergence.wls      # TIME4
zsh run_all.sh                               # all of the above and the note
```

Each script prints `PASS`/`FAIL` lines against stated tolerances and ends with a `SUMMARY`; the checks are written to `../results/massive_stageN_checks.csv`, the tables to `../results/massive_stageN_summary.md` and CSVs, the figures (vector PDF + CSV) to `../figures/`.

## Results

### Stage 1: first order (`massive_stage1_first_order.wls`, 406 checks pass, 7.7 min)

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

RESULTS-STAGE2

