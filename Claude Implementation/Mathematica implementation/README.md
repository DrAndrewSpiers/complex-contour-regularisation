# Mathematica implementation: complex-contour S/R split, first and second order (1/ρ² case)

This folder implements and tests the method in `../singular-field-consistency.tex` (sections 2, 4, 6, 7 and 8) for the flat-space toy model.
It covers three things:

- the first-order singular mode and its consistency with 3D;
- the first-order regular field;
- the second-order pieces SS, SR and RR with Q[A,B] = AB.

Everything was run with Mathematica 14.3 (wolframscript) on 6 October 2026.
All 2214 checks pass in the final `run_all.sh` run (8 October 2026): 627 for Stages 1–3 (79 + 76 + 136 + 123 + 213), 371 for Stage 4, 489 for Stage 5 (143 + 137 + 209), 355 for Stage 6, 12 for Stage 7, 34 for Stage 8 and 326 for Stage 9.

Parameters: q = 1, r₀ = 5, domain r ∈ [1, 10], worldtube [a_L, a_R] = [4.5, 5.5].
The contour is the notebook's, y(a) = a + i(1 − 4(a−5)²)⁶.
The injected regular field is Φᴿ = c₀ + c₁ r cosθ with c₀ = 1, c₁ = 1/10.

## Contents

| File | Purpose |
|---|---|
| `ContourSecondOrder.wl` | Package with all reusable functions. It covers the contour and operators (`contourL` is the notebook's `odeOpComplex`), branch handling, the Stage 1 singular modes, the solvers, the second-order sources and couplings, and closed-form analytic modes. |
| `stage1_singular_mode.wls` | Stage 1: DSolve's ψ_p, D_ℓ, the corrected φˢ_ℓ, (C1) and (C2), continuity along the contour, and the 3D meaning of D_ℓ. |
| `stage2_regular_field.wls` | Stage 2: the notebook's three-region solve with the exact φˢ_ℓ, the full-field solve, the real-line solve, the BVP and a truncated-puncture diagnostic. |
| `stage3a_source_modes.wls` | Stage 3 inputs: the Neumann formula, LegendreQ conventions, the continued SS source, Gaunt coefficients, K_ℓℓ′ and its two continuations, SR/RR source modes, and the analytic particular solutions. |
| `stage3b_second_order_solve.wls [lmax]` | Stage 3 solves: SS, SR and RR with Cauchy and two-point boundary conditions, shadows at a_R, cross-checks, negative controls and the 3D reconstruction. |
| `stage4_deferred_shadow_test.wls` | Stage 4: the deferred shadow test, i.e. the truncated-puncture residual field used directly in SR (scan over n, Lemma-1 decomposition, amplitude and radius sensitivity, ℓ′ ≥ 2 modes). |
| `stage5a_gradient_sources.wls` | Stage 5 (1/ρ⁴, Q = ∇A·∇B) inputs: [1/ρ⁴]_ℓ via Q_ℓ′, its continuation, ∂_rK and K^θ couplings, SR/RR source modes, the Φˢ·Φᴿ identity, and analytic modes. |
| `stage5b_gradient_solve.wls [lmax]` | Stage 5 solves: SS and SR on the contour, RR on the real line, Cauchy and BVP, shadows, the SS regularisation test with positive controls, amplitude independence and 3D reconstruction. |
| `stage6_3d_puncture_test.wls [quick]` | Stage 6: shadow test with a consistent 3D puncture Φᴾ = q/ρ + cqρⁿ (decomposition into (a), (b), (c); n and c scans; ℓ′ convergence; amplitude independence). |
| `stage7_lprime_optimum.wls` | Stage 7: total exterior error = (b) + (c) against ℓ′max for a fixed and a shrinking worldtube (validated decomposition); writes `../figures/fig_lprime_optimum.pdf`. |
| `stage8_hybrid_split.wls` | Stage 8: hybrid split of SR (Taylor part in closed form; remainder formed on the real line and Cauchy-continued; baseline H0 without the split); writes `../figures/fig_hybrid.pdf`. |
| `stage9_dzz.wls [smoke]` | Stage 9: distributional SR, Q = ½∂_z²(AB): δ(r − r₀) terms in the SR modes; exact puncture (SS, SR, total) and Stage 6 puncture (product route, H0, H1 with k = 0…3), validated by contour solves; writes `../figures/fig_dzz.pdf`. |
| `stage10_teukolsky_benchmark.wls` | Stage 10: contour against the mode-level Hadamard finite part of Musaeus, Pound & Upton (arXiv:2609.10521, Sec. V.D, eqs. 140–142, A2–A5). |
| `stage11_TQ.wls` | Stage 11: their T_Q in the ∂_z² toy; ball vs slab regularisation of the mode δ terms; their prescription end to end. |
| `stage12_ss_reuse_convergence.wls [A] [B] [C] [D]` | Stage 12 (Dzz): SS from an analytic singular part plus a Cauchy-continued numerical remainder (A); subtraction-free radial factor (C); complex-r angular quadrature (D); H1 convergence and tail acceleration (B). Writes `../figures/fig_ss_reuse.pdf` and `../figures/fig_hybrid_convergence_cauchy.pdf` (Cauchy data at r = 1). |
| `stage14_cost.wls` | Stage 14: cost of the SS angular integrals, real axis against contour (counted evaluations of the 3D integrand, wall time, error); writes `../figures/fig_cost.pdf`. |
| `stage14b_radial.wls` | Stage 14b: real-axis radial SS samples (Hadamard finite part, uniform or graded grids) needed for the exterior field, against the contour's nodes; writes `../figures/fig_radial.pdf`. |
| `stage12b_physical_bc.wls` | Stage 12b: the H1 convergence study and the PRL Fig. 5 route comparison with physical boundary conditions (free-space Green's function); writes `../figures/fig_hybrid_convergence.pdf` and `results/stage12b_routes.csv` (read by `make_prl_figures.wls`). |
| `stage7b_physical_bc.wls` | Stage 7 (ℓ′max optimum, PRL App. D) with physical boundary conditions; writes `../figures/fig_lprime_optimum.pdf` (the Cauchy version is `fig_lprime_optimum_cauchy.pdf`). |
| `make_prl_figures.wls` | Stage 13: PRL figures for the Dzz case in `../figures/prl/` (vector PDF + CSV). Needs `results/stage12b_routes.csv`. |
| `make_figures.wls [lmaxTot]` | Figures for the note, written to `../figures/`: `fig_modes_recovered`, `fig_relative_error`, `fig_mode_convergence` (PDF and CSV). |
| `make_notebook.wls` | Builds `ContourSecondOrder_workflow.nb` and checks that it parses (with Get, Import and the front end). |
| `ContourSecondOrder_workflow.nb` | Unevaluated notebook version of the main workflow. All 19 input cells were executed headlessly without messages. |
| `run_all.sh` | Runs everything and writes logs to `results/logs/`. |
| `results/` | CSV tables, a Markdown error table, PDFs and logs (see below). |

## How to run

```
W=/Applications/Wolfram.app/Contents/MacOS/wolframscript
cd "<this folder>"
$W -file stage1_singular_mode.wls             #   8 s
$W -file stage2_regular_field.wls             #  40 s
$W -file stage3a_source_modes.wls             #  12 s
$W -file stage3b_second_order_solve.wls       # 110 s (lmax = 6)
$W -file stage3b_second_order_solve.wls 12    # 200 s (files tagged _lmax12)
$W -file make_notebook.wls
$W -file stage4_deferred_shadow_test.wls     # ~22 min
$W -file stage5a_gradient_sources.wls         #   5 s
$W -file stage5b_gradient_solve.wls 6        #   5 min (and 12: 11 min, files tagged _lmax12)
$W -file stage9_dzz.wls                      #   6 min
$W -file stage10_teukolsky_benchmark.wls      #   1 min
$W -file stage11_TQ.wls                     #   2 min
$W -file stage12_ss_reuse_convergence.wls     #  17 min (A, C, D: 11 min; B: 6 min)
$W -file stage12b_physical_bc.wls            #  11 min (before make_prl_figures)
$W -file stage7b_physical_bc.wls             #   3 min
$W -file stage12_ss_reuse_convergence.wls phys  # A, C, D with physical conditions
$W -file stage11_TQ.wls phys
$W -file stage14b_radial.wls phys
$W -file make_prl_figures.wls                 #   2 min
$W -file stage14_cost.wls                     #  23 min
$W -file stage14b_radial.wls                  #   8 min
zsh run_all.sh                                # all of the above
```

Each script prints `PASS`/`FAIL` lines against stated tolerances and ends with a `SUMMARY`.
The checks are also written to `results/stageN_checks*.csv`.

## Conventions

- Modes are Φ = Σ φ_ℓ(r) Y_ℓ0(θ) with orthonormal Y_ℓ0, and [f]_ℓ = √(π(2ℓ+1)) ∫₋₁¹ f P_ℓ dμ.
- The radial operator is L_ℓψ = ψ″ + 2ψ′/r − ℓ(ℓ+1)ψ/r². The package uses the formal variable `zVar` for the complex radial variable. The Cartesian z is always written r cosθ.
- First order: ∇²Φ₁ = −4πqδ³.
  - The mode source is κ_ℓ δ(r − r₀), with κ_ℓ = −(q/r₀²)√(4π(2ℓ+1)).
  - φˢ_ℓ is at unit strength, so 3D quantities use κ_ℓ φˢ_ℓ.
  - G_ℓ = −(r₀²/(2ℓ+1)) r_<^ℓ/r_>^{ℓ+1}.
- H(z) = 1 + (i/π) log(z − r₀) and δ_ext = i/(π(z − r₀)). Values on the real axis are boundary values **from above**.
- Second order: ∇²Φ₂ = Φ₁² = q²/ρ² + 2qΦᴿ/ρ + (Φᴿ)².
  - SS and SR are solved on the contour inside the tube and on the real line outside it.
  - RR is solved on the real line throughout.
  - Boundary conditions come from the analytic particular solutions:
    - SS: q² ln ρ.
    - SR: qc₀ρ + qc₁[r₀ρ + ½(r cosθ − r₀)ρ].
    - RR: c₀²r²/6 + c₀c₁ r³cosθ/5 + c₁² r⁴cos⁴θ/12.
  - "Cauchy" means value and derivative at r = 1. "BVP" means values at r = 1 and r = 10, by linear shooting on φ′(1) as in the notebook's section 3.
- **Scaled error** = max|numerical − analytic| / max|analytic|. It is taken over 81 exact points in [1, a_L] ∪ [a_R, 10]. If the analytic mode is identically zero, the absolute error is used.

## Results

### Stage 1: what DSolve's ψ_p actually is

DSolve's antiderivatives are written with **`Log[r₀ − z]`**, together with `Log[z − r₀]` for ℓ = 0, 1 and `Log[z]` for ℓ ≥ 2.
For example, at ℓ = 0:

ψ_p = (i/2π)(75 + z²)/z + (5i/π) Log[5 − z] − (25i/π) Log[z − 5]/z.

Proposition 1 holds exactly: D_ℓ = Re ψ_p(r+i0) − G_ℓ = α r^ℓ + β r^{−ℓ−1}.
The check is in exact arithmetic: the fit uses r = 2, 3 and is verified at 17 points on both sides of r₀, with zero residual. A separate fit at r = 7, 8 gives the same (α, β).

| ℓ | α | β | which case |
|---|---|---|---|
| 0 | r₀ = 5 | r₀² = 25 | note's prediction **plus** mirror |
| 1 | 0 | 0 | Re ψ_p = G₁ exactly |
| ℓ ≥ 2 (checked to 12) | r₀^{1−ℓ}/(2ℓ+1) | 0 | the note's *mirror* (Log[r₀−z]) case |

The note predicts α = 0 and β = r₀^{ℓ+2}/(2ℓ+1). DSolve's ψ_p does not match that for any ℓ (`results/stage1_D_table.csv`).

- For ℓ = 0, ψ_p puts the r^ℓ piece of h_ℓ on Log[r₀−z] and the r^{−ℓ−1} piece on Log[z−r₀]. The result is the "anti-physical" Green's function: Re ψ_p = 25/r inside and 5 outside.
- For ℓ = 1, the assignment is the physical one.
- For ℓ ≥ 2, both pieces sit on Log[r₀ − z], so Re ψ_p vanishes for r < r₀.

In 3D, for r < r₀, Σ κ_ℓ D_ℓ Y_ℓ = −q/ρ + q r cosθ/r₀² − q/r. This was checked at one point with ℓ ≤ 80 to 1e-34. For r > r₀ the sum diverges.

**Branch pitfall (confirmed).** Mathematica evaluates `Log[r₀ − z]` at real z > r₀ on the wrong side of its cut.
The raw DSolve form, evaluated directly at real r > r₀, differs from the from-above value by up to 10 (ℓ = 0) and 49 (ℓ = 6).
A "D_ℓ" built that way is not a homogeneous solution: the same fit leaves residuals of 3 to 49 for r > r₀.

`canonicalizeLogs` rewrites Log[c(z − z*)] as Log|c| + Log[z − z*] − iπ·[c < 0]. This is the same analytic function in the upper half plane, and its principal value on the axis is the value from above.
The rewrite was checked against the raw form at 15 upper-half-plane points (difference 0) and against r + i·10⁻⁴⁰ (difference ≤ 1.1e-35, which is the ε effect).

The notebook itself was not run. By inspection, its section 1 evaluates `exactS[aR]` and `localS[aR]` at the real a_R = 5.5, so the Log[5 − r] there is on the wrong branch.
- Evaluated that way at a_R, the raw ψ_p − D differs from the contour continuation by 10.0, 2.75, 0.92, … for ℓ = 0, 1, 2, ….
- In the notebook's `taylorValR − exactValR` the two wrong-branch logs largely cancel. What remains is 2πi × (Taylor remainder of the log coefficient), so the error shrinks with the expansion order.

**Corrected mode φˢ_ℓ = ψ_p − D_ℓ** (ℓ = 0…6):
- (C1) Re φˢ_ℓ(r+i0) = G_ℓ exactly at 17 points, and to ≤ 1.1e-35 on a 180-point grid.
- (C2) L_ℓφˢ_ℓ − δ_ext simplifies to 0 symbolically, and is < 1e-55 on the contour with the notebook's operator.
- Continuity along the contour: the NDSolve continuation from a_L reproduces φˢ_ℓ(y(a)) to ≤ 1.4e-19 at a_R, so there are no 2πi jumps.

### Stage 2: first-order regular field (S_eff = 0)

All 76 checks pass.

- S_eff is computed with the notebook's formula `complexSource[y] − odeOpComplex[localS]` and is < 1e-65 along the contour.
- The contour φᴿ matches the analytic continuation (for example c₁ z √(4π/3)) to ≤ 5.6e-19 (scaled).
- Im φᴿ(a_R) and Im φᴿ′(a_R), which the notebook drops via `Re[valR]`, are ≤ 6.4e-19 for ℓ = 0, 1. They are exactly 0 for ℓ = 2, 3.
- The exterior φᴿ is within 1.3e-19.
- The full field (source κ_ℓ δ_ext) equals φᴿ + κ_ℓφˢ_ℓ on the contour (≤ 1.5e-19), with no residual shadow at a_R.
- The real-line solve through the tube is exact to 1e-20.
- The notebook's generic test (φ(1) = 1, φ′(1) = 0) passes for ℓ = 2, 3, 6 (≤ 2.5e-19, Im at a_R ≤ 1.2e-19). The BVP version (ℓ = 0, 1, 2) also passes.

Diagnostic with a truncated puncture (Taylor series of φˢ₂ through order n):

| n | 2 | 4 | 6 | 8 | 12 |
|---|---|---|---|---|---|
| \|Im residual(a_R)\| + \|Im residual′(a_R)\| | 6.5e-4 | 3.9e-5 | 1.5e-6 | 3.7e-8 | 1.2e-11 |

- These shadows equal κ Im(φˢ − puncture)(a_R) to 1e-17.
- The reconstructed exterior φᴿ is exact (≤ 4e-19) for every n.

### Stage 3: second order

**Inputs (stage3a, 136 checks pass).**
- [1/ρ²]_ℓ matches √(π(2ℓ+1)) Q_ℓ(x)/(r r₀) to < 1e-40 relative.
- Q_ℓ = ½P_ℓ ln((x+1)/(x−1)) − W_{ℓ−1} equals `LegendreQ[ℓ, 0, 3, x]`, the type 3 function for x > 1.
- The continued SS source uses 2 log(z+r₀) − 2 log(z−r₀). Its real part on the axis equals the real source (< 1e-31), and it has no jumps along the contour.
- The Gaunt coefficients match exact integration.
- K_ℓℓ′ matches direct angular integration (< 1e-40).
- Re K̃ on the axis equals K for both continuations: Θ→H, and the corrected φˢ_L of Stage 1.
- The SR source (2q Σ K φᴿ) and the RR source (Gaunt coupling) agree with closed-form angular integrals.
- The analytic particular-solution modes (closed form by u-substitution):
  - satisfy L_ℓψ = source symbolically on both branches;
  - are C¹ across r₀;
  - agree with direct μ-integration of the 3D functions (< 3e-27);
  - for ℓ = 0, reproduce `AnalyticSol` from `Greens-function-Poisson-equation-3.nb`.

**Exterior errors (stage3b, lmax = 6, 123 checks pass; `results/stage3_exterior_errors.md`).**
The tolerance was 1e-10 scaled. The last three columns are the shadows Im ψ(a_R) of the Cauchy runs.

| ℓ | SS Cauchy | SS BVP | SR Cauchy | SR BVP | SR corr-K | RR Cauchy | RR BVP | Im SS(a_R) | Im SR(a_R) Θ→H | Im SR(a_R) corr-K |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3.2e-20 | 1.9e-20 | 3.3e-20 | 4.0e-20 | 2.5e-20 | 3.3e-22 | 1.1e-20 | −0.076 | −0.010 | −0.355 |
| 1 | 7.6e-19 | 8.2e-19 | 1.8e-19 | 4.2e-19 | 1.1e-19 | 2.6e-22 | 6.2e-21 | −0.133 | −0.018 | 0.164 |
| 2 | 1.5e-18 | 2.0e-18 | 4.5e-19 | 6.1e-19 | 3.7e-19 | 2.6e-22 | 1.3e-20 | −0.175 | −0.024 | 0.416 |
| 3 | 9.0e-18 | 1.1e-17 | 5.8e-18 | 6.2e-18 | 5.5e-18 | 0 | 0 | −0.214 | −0.029 | 0.329 |
| 4 | 2.4e-17 | 2.7e-17 | 7.2e-17 | 7.1e-17 | 7.1e-17 | 2.1e-20 | 1.5e-20 | −0.252 | −0.035 | 0.256 |
| 5 | 3.4e-17 | 4.9e-17 | 1.3e-16 | 1.3e-16 | 1.4e-16 | 0 | 0 | −0.293 | −0.041 | 0.223 |
| 6 | 1.4e-16 | 1.7e-16 | 1.7e-16 | 1.6e-16 | 1.6e-16 | 0 | 0 | −0.338 | −0.048 | 0.207 |

- With lmax = 12 (`*_lmax12*`, 213 checks pass), all pieces stay ≤ 7.1e-15 scaled through ℓ = 12.
- The total SS + SR + RR per mode is ≤ 4.1e-15.
- RR, checked on the real line inside the tube as well, has errors ≤ 2.1e-20.

**Cross-checks.**
- **Two continuations of K.** Θ→H and the corrected φˢ_L give the same exterior SR: ≤ 8.1e-18 for ℓ ≤ 6 and ≤ 5.3e-15 for ℓ ≤ 12. Only the shadow differs (table above).
- **φ̃ᴿ input.** Using the exact continuation of φᴿ instead of the Stage 2 contour solution changes the exterior SR by ≤ 2.2e-18 (ℓ ≤ 6) and ≤ 2.1e-16 (ℓ ≤ 12). The shadow at a_R changes by ≤ 4e-19.
- **SS on the whole contour.** The complex SS contour solution equals the r > r₀ formula of [q² ln ρ]_ℓ, with Log(r − r₀) → Log(z − r₀), plus i × a homogeneous solution. The match is ≤ 1.5e-18 for ℓ ≤ 12.
  - The real part of that homogeneous solution is exactly 0, so there is no D_ℓ-type ambiguity at second order for SS.
  - The shadow Im ψ(a_R) is therefore predicted analytically, to ≤ 7e-20.

**3D reconstruction.**
The method error, numerical versus analytic mode sums at the same ℓ_max, is ≤ 5e-17 in absolute terms (Φ₂ ≈ 6–28) at r = 2, 3, 4, 6, 8, 9.5. This holds for Cauchy and BVP alike.
The truncation error against the exact 3D Φ₂ is a separate effect:

| ℓ_max | truncation error |
|---|---|
| 6 | 2e-4 to 1.7e-2 |
| 12 | 3.5e-7 to 1.7e-3 |

At the points nearest the tube convergence is slow, because the SS and SR modes fall off only like (r_</r_>)^ℓ.
Details are in `results/stage3_reconstruction*.csv`.

**Negative controls** (`results/stage3_negative_controls.csv`). These are deliberate mistakes; the exterior test must detect them.

| Control | ℓ = 0 | ℓ = 1 | ℓ = 2 |
|---|---|---|---|
| A: SS source continued with Log[(z−r₀)²] (jumps 2πi at a = r₀) | 6.7e-2 | 0.74 | 2.9 |
| B: φ̃ᴿ₀ from the Θ→H puncture in SR (violates C2; \|Im φ̃ᴿ₀\| up to 0.48), K by Θ→H | 2.1e-5 | 2.8e-4 | 1.1e-3 |
| B, same, K by corrected φˢ | 6.4e-3 | 2.5e-2 | 9.0e-2 |

In control B the *first-order* exterior φᴿ is still exact (1e-20), so only the second-order SR exposes the shadow.

The detection threshold was first set at 1e-3. The ℓ = 0, 1 errors with Θ→H K (2e-5, 3e-4) fell below it, so I changed the threshold to 1e-8 (100 × the tolerance). That is still far above the clean-run errors.

### Stage 4: deferred shadow test

**Question.** The first-order field inside the tube can come from a *truncated* puncture. The residual field then carries a small shadow at the truncation order. Does that shadow matter when the residual is used directly in the SR source?

**Setup.**
- The puncture φᴾ_{ℓ,n} is the series of the corrected φˢ_ℓ about r₀ through order n, with canonical logs. It is the analogue of the notebook's `phiSLocal`.
- The residual φ̃ʳᵉˢ_{ℓ′,n} is solved on the contour with S_eff,n = κ(δ_ext − L φᴾ_n). It starts from φᴿ(a_L) + κ(φˢ − φᴾ_n)(a_L), on the correct branch.
- It agrees with the analytic φ̃ᴿ + κ(φˢ − φᴾ_n) along the contour to ≤ 4e-19 (scaled), for every n, amplitude and ℓ′.
- Inside the tube, SR then uses s̃_n = 2q Σ_{ℓ′∈{0,1}} K̃_{ℓℓ′} φ̃ʳᵉˢ_{ℓ′,n}. Outside, the source is unchanged.
- Scan: n ∈ {2, 3, 4, 6, 8, 10, 12}, ℓ = 0…6, both continuations of K, Cauchy data at r = 1.
- Errors are scaled by max|SR analytic| over the exterior, then maximised over ℓ.
- The baseline with the exact φˢ is 1.7e-16.

**Decomposition (Lemma 1).** Write δφ = κ(φˢ − φᴾ_n). The exterior error is the response to the real boundary value of 2q Σ K̃ δφ on [a_L, a_R], which splits into two parts:
- **(a)** 2q Σ K Re δφ: the real-part truncation remainder. It exists only because a mode-by-mode Taylor puncture is inconsistent with the exact 3D SS piece. With a genuine 3D puncture used consistently (condition C3), (a) would vanish.
- **(b)** −2q Σ Im K̃ Im δφ: the shadow contamination. This is the contour-specific error. It would remain with a consistent 3D puncture, because a truncated puncture's S_eff has a branch point at r₀.

(a) and (b) were computed by Gauss–Legendre quadrature of the variation-of-parameters Green's function, split at r₀.
- Quadrature with 60 and 120 nodes per side agrees to 2e-16 (relative).
- (a) + (b) reproduces the step-2 three-region error at a_R to ≤ 4e-18 in all 98 (n, ℓ, K) cases.

| n | (a) real remainder | (b) shadow, Θ→H | total, Θ→H | (b)/(a), Θ→H | (b) shadow, corrected φˢ | total, corrected φˢ | (b)/(a), corrected φˢ | (R/r₀)ⁿ |
|---|---|---|---|---|---|---|---|---|
| 2 | 2.7e-3 | 1.7e-3 | 4.4e-3 | 0.65 | 1.5e-3 | 1.2e-3 | 0.55 | 1e-2 |
| 3 | 2.5e-4 | 8.4e-5 | 1.6e-4 | 0.34 | 3.8e-4 | 6.3e-4 | 1.55 | 1e-3 |
| 4 | 2.3e-5 | 1.5e-5 | 3.8e-5 | 0.65 | 1.5e-5 | 8.5e-6 | 0.64 | 1e-4 |
| 6 | 2.2e-7 | 1.4e-7 | 3.5e-7 | 0.64 | 1.5e-7 | 6.8e-8 | 0.68 | 1e-6 |
| 8 | 2.1e-9 | 1.3e-9 | 3.4e-9 | 0.64 | 1.5e-9 | 5.9e-10 | 0.71 | 1e-8 |
| 10 | 2.0e-11 | 1.3e-11 | 3.3e-11 | 0.63 | 1.5e-11 | 5.4e-12 | 0.73 | 1e-10 |
| 12 | 1.9e-13 | 1.2e-13 | 3.2e-13 | 0.63 | 1.5e-13 | 5.2e-14 | 0.75 | 1e-12 |

Fitted decay factors per order (max over ℓ):

| quantity | Θ→H | corrected φˢ |
|---|---|---|
| (a) | 0.097 | 0.097 |
| (b) | 0.0995 | 0.095 |
| total | 0.101 | 0.082 |

Per-ℓ fits are in `results/stage4_fits.csv`: 0.095–0.11 for (b) with Θ→H, and 0.076–0.098 with corrected φˢ.
Per ℓ, (b)/(a) ranges from 0.007 to 0.65 (Θ→H) and from 0.03 to 1.55 (corrected φˢ).
With corrected φˢ, (a) and (b) partly cancel; with Θ→H they add.

**Sensitivity.**
- *Contour amplitude* (amp = ½, 2 at n = 2, 4, 8, all ℓ and both K). δ(a_R) is unchanged to ≤ 3e-18 (absolute), as Lemma 1 predicts. This holds even though the residual's shadow on the contour, max|Im(φ̃ʳᵉˢ − φ̃ᴿ)|, changes strongly with amp. For n = 2, 4, 8 it is:
  - amp ½: 7.8e-4, 1.2e-5, 1.9e-9;
  - amp 1: 6.4e-3, 3.8e-4, 1.0e-6;
  - amp 2: 5.4e-2, 1.3e-2, 5.5e-4.

  On the axis the shadow is only about 0.1ⁿ; on the contour it is about (|z − r₀|/r₀)ⁿ.
- *Tube radius* (decomposition, exterior r ∈ [6, 10]; the three-region solve was confirmed at n = 4 for R = ¼ and 1). The errors scale as (R/r₀)ⁿ. For example, (b) with Θ→H at n = 8 is 1.8e-12, 1.3e-9 and 1.5e-7 for R = ¼, ½, 1. See `results/stage4_radius.csv`.

**Secondary: ℓ′ ≥ 2 residual modes** (pure truncation remainder; `results/stage4_lprime_*.csv`). The ℓ′ sum does **not** converge. The contribution of mode ℓ′ grows polynomially with ℓ′, and faster for larger n (max over ℓ, scaled):
- n = 4, (a): 3.3e-6 at ℓ′ = 0, 3.7e-3 at 6, 9.4e-3 at 12.
- n = 8, (a): 1.8e-10 at ℓ′ = 0, 2.3e-6 at 6, 1.9e-5 at 12.

The cumulative total for ℓ′max = 1, 2, 4, 8, 12 is:
- n = 4 (Θ→H): 3.8e-5, 1.4e-4, 1.0e-3, 1.4e-2, 4.1e-2.
- n = 8 (Θ→H): 3.4e-9, 1.9e-8, 2.9e-7, 1.1e-5, 6.1e-5.

The three-region solves agree with the summed per-ℓ′ decomposition in all cases. As expected, a mode-by-mode Taylor truncation is not the set of modes of a 3D field. This part is not representative of a 3D or 4D puncture.

**Conclusion.** The shadow contribution (b) is not parametrically smaller than the truncation error.
- At every n it is 0.6–0.75 times the real-part truncation remainder (a), taking the max over ℓ, and it decays at the same rate, (R/r₀)ⁿ ≈ 0.1ⁿ.
- Relative to the SR field it is about 0.15 (R/r₀)ⁿ: 1.5e-5 at n = 4, 1.4e-7 at n = 6, 1.3e-9 at n = 8, 1.2e-13 at n = 12.
- It is independent of the contour amplitude, as Lemma 1 requires.

So it is "negligible" exactly when the puncture's own truncation error is. With a consistent 3D puncture, (a) vanishes and (b) becomes the leading puncture-order error at about 0.15 (R/r₀)ⁿ. Choose n, or a smaller R/r₀, for the target accuracy: n ≥ 6 gives ≤ 1.5e-7 relative here.

The two K continuations give comparable (b). I recommend **Θ→H** for three reasons:
- it is closed-form and needs no DSolve;
- it gives second-order shadows 4 to 35 times smaller (Stage 3);
- it was 80 to 300 times less sensitive to an O(1) shadow (control B).

### Stage 5: 1/ρ⁴ case (Q = ∇A·∇B)

**Setup.** The second-order equation is ∇²Φ₂ = |∇Φ₁|² = q²/ρ⁴ + 2q∇(1/ρ)·∇Φᴿ + |∇Φᴿ|².
- The injected field is Φᴿ = c₀ + c₁ r cosθ + c₂ r² P₂(cosθ), with c₂ = 1/100. Couplings with ℓ′ = 0, 1, 2 are therefore exercised.
- Particular solutions:
  - SS: q²/(2ρ²).
  - SR: q[Φᴿ(x) − Φᴿ(x₀)]/ρ.
  - RR: (Φᴿ)²/2.
- In total Φ₂ = ½Φ₁² − qΦᴿ(x₀)/ρ.
- SS and SR are solved on the contour inside the tube; RR on the real line. ∂_rK and K^θ are continued by Θ→H, with corrected φˢ as a cross-check.

**Inputs (stage5a, 143 checks pass).**
- **SS source.** [1/ρ⁴]_ℓ = −√(π(2ℓ+1)) Q′_ℓ(x)/(2r²r₀²) matches direct μ-integration (< 1e-20 relative) and the closed form. For ℓ = 0 it is 2√π/(r² − r₀²)². The Neumann-form Q′_ℓ equals `D[LegendreQ[ℓ,0,3,x], x]`.
- **Continuation.**
  - The logarithm becomes 2log(z+r₀) − 2log(z−r₀).
  - 1/(x²−1) = 4z²r₀²/((z−r₀)²(z+r₀)²) is simply evaluated at z.
  - Re on the axis equals the real source for r ≠ r₀, with no jumps along the contour.
  - L_z[(q²/2)·continued [1/ρ²]_ℓ] − continued [q²/ρ⁴]_ℓ simplifies to **0** symbolically (ℓ ≤ 6). This is the identity behind the regularisation result below.
- **SR couplings.**
  - K^θ_ℓℓ′ = Σ_L ⟨ℓℓ′L⟩[1/ρ]_L·½[L(L+1) + ℓ′(ℓ′+1) − ℓ(ℓ+1)] and ∂_rK both match direct angular integration (13 pairs × 4 radii, < 1e-20).
  - The continued values have the correct real part, for Θ→H and for corrected φˢ.
- **SR source and the ΦˢΦᴿ identity.**
  - The coupling-formula SR source equals the closed form of [2q∇(1/ρ)·∇Φᴿ]_ℓ and direct 3D integration.
  - Pointwise, L_ℓ[qΦᴿ/ρ]_ℓ = SR source (symbolic).
  - The derivative jump of [qΦᴿ/ρ]_ℓ at r₀ equals −4πqΦᴿ(x₀)Y_ℓ(0)/r₀² exactly. This is the δ(r−r₀) term in mode space.
- **RR source.** The gradient–Gaunt coupling (angular weight ½[ℓ₁(ℓ₁+1) + ℓ₂(ℓ₂+1) − ℓ(ℓ+1)]) equals the closed form.
- **Analytic modes.** All satisfy the ODE symbolically, SR and RR modes are C¹ across r₀, and all agree with direct μ-integration.

**Solves (stage5b).** lmax = 6 had 137 checks and lmax = 12 had 209; all pass. The table gives scaled exterior errors, the Cauchy-run shadows, and the fitted point-term coefficients (see below).

| ℓ | SS Cauchy | SS BVP | SR Cauchy | SR BVP | SR corr-K | RR Cauchy | RR BVP | Im SS(a_R) | Im SR(a_R) Θ→H | Im SR(a_R) corr-K | α | β |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 4.2e-19 | 2.3e-19 | 1.1e-19 | 6.0e-20 | 4.3e-20 | 1.7e-21 | 6.3e-21 | 0.101 | 0.0134 | 0.0193 | −4.5e-20 | 1.5e-19 |
| 1 | 5.9e-19 | 4.0e-19 | 2.5e-19 | 9.8e-20 | 4.7e-19 | 9.0e-21 | 8.5e-21 | 0.180 | 0.0233 | 0.0318 | −1.8e-20 | 1.6e-19 |
| 2 | 6.1e-18 | 4.6e-18 | 6.8e-19 | 3.4e-19 | 3.6e-19 | 6.1e-21 | 1.7e-20 | 0.244 | 0.0305 | 0.0277 | −3.2e-19 | 5.0e-19 |
| 3 | 9.4e-18 | 3.4e-18 | 3.2e-18 | 1.5e-18 | 2.8e-19 | 1.9e-20 | 8.8e-21 | 0.310 | 0.0369 | 0.0183 | −2.8e-19 | 3.1e-19 |
| 4 | 1.3e-17 | 2.7e-18 | 5.5e-17 | 6.4e-17 | 6.7e-17 | 1.8e-20 | 1.5e-20 | 0.386 | 0.0430 | 0.0125 | −6.4e-22 | 3.8e-19 |
| 5 | 4.6e-17 | 1.3e-17 | 1.3e-16 | 1.6e-16 | 1.6e-16 | 0 | 0 | 0.476 | 0.0492 | 0.0065 | 6.1e-19 | −8.2e-20 |
| 6 | 5.5e-17 | 1.2e-16 | 7.7e-17 | 1.3e-16 | 1.5e-16 | 0 | 0 | 0.587 | 0.0558 | −0.0001 | 5.3e-19 | 7.9e-20 |

Further results:
- **ℓ = 7…12.** All pieces stay ≤ 2e-14 scaled (`results/stage5_exterior_errors_lmax12.md`).
- **Cross-checks.** The two K continuations agree in the exterior to ≤ 7e-17. The total SS + SR + RR per mode is ≤ 4e-17 for ℓ ≤ 6.
- **SS contour solution.** As in Stage 3, the whole complex SS contour solution equals the continued analytic mode plus i × homogeneous (≤ 8e-19), and its shadow is predicted to ≤ 1.2e-19.
- **3D reconstruction** against ½Φ₁² − qΦᴿ(x₀)/ρ:
  - The method error is ≤ 3e-16 (absolute; Φ₂ ≈ 0.3–0.9).
  - Truncation is 8e-6 to 6e-4 at ℓ_max = 6 and 5e-8 to 1e-4 at ℓ_max = 12.

**SS regularisation: the contour adds no point-supported terms.**
1/ρ⁴ is defined at x₀ only up to cδ³ + cⁱ∂_iδ³. I fitted the SS discrepancy (Cauchy run, region 3) per ℓ to α[1/ρ]_ℓ + β∂_{r₀}[1/ρ]_ℓ, i.e. a monopole plus a z-dipole at x₀:
- **Result.** |α| ≤ 6e-19 and |β| ≤ 5e-19 for every ℓ ≤ 6; |α|, |β| ≤ 1.1e-18 for ℓ ≤ 12. They are zero, not ℓ-independent constants.
- **Positive controls.** I added 1·δ³ (as δ_ext) or 1·∂_{r₀}δ³ (as −δ′_ext − 2δ_ext/r₀) to the contour source. The fit recovers α = −1/(4π) = −0.0795774715459 and β = −1/(4π) respectively, ℓ-independently to 1e-18 (`results/stage5_regularisation.csv`). The fit therefore would detect such terms.
- **Contour amplitude.** Changing amp to ½ or 2 leaves the exterior SS and SR data at a_R unchanged to 7e-19 and 5e-20.

This confirms the prediction. The contour's real part on the axis is the mode-wise Hadamard finite part, which equals L_ℓ applied distributionally to (q²/2)[1/ρ²]_ℓ (the symbolic identity above). The contour SS therefore reproduces [q²/(2ρ²)]_ℓ with no extra terms. Summed over ℓ this is q²/(2ρ²), whose distributional Laplacian is the spherical-cutoff finite part of q²/ρ⁴. By Green's identity, FP_{ε→0} ∫_{ρ>ε} f/ρ⁴ = ∫ ∇²f/(2ρ²) for a test function f; the discarded divergence is 4πf(x₀)/ε. I checked that last step analytically, not numerically.

**Observations.**
- **SS shadows are larger than in the 1/ρ² case:** 0.10 to 0.59 for ℓ ≤ 6, and 2.2 at ℓ = 12.
- **SR shadows.** With Θ→H they are 0.013 to 0.056 for ℓ ≤ 6 (0.12 at ℓ = 12). With corrected φˢ they are *smaller* for ℓ ≥ 3 (1e-4 at ℓ = 6, 0.07 at ℓ = 12), the opposite of the 1/ρ² case. The exterior is the same either way.

### Stage 6: shadow test with a consistent 3D puncture

`stage6_3d_puncture_test.wls` (355 checks). Full tables: `results/stage6_summary.md`.
The agent was interrupted before writing this section; it was written afterwards from those outputs.

**Setup.**
- The puncture is Φᴾ = q/ρ + cqρⁿ, with n odd ∈ {1, 3, 5, 7} and c = r₀^{−(n+1)}. One extra run uses 10c at n = 3.
- It is used everywhere and in every piece (C3), so Φ^res = Φᴿ − cqρⁿ has all ℓ′.
- The residual on the contour comes from the first-order contour solve with S_eff = −L_ℓ[Θ→H continuation of cq[ρⁿ]_ℓ].
- The exterior error is split into three parts via Lemma 1:
  - (a) real-part mismatch;
  - (b) shadow contamination;
  - (c) truncation of the ℓ′ coupling sums at ℓ′max.

**Results** (scaled exterior error, max over ℓ ≤ 6, ℓ′max = 12):

| case | n | (b) shadow | (c) ℓ′ truncation | total |
|---|---|---|---|---|
| Product | 1 | 1.2e-2 | 5.7e-2 | 4.5e-2 |
| Product | 3 | 4.6e-5 | 2.1e-3 | 2.2e-3 |
| Product | 5 | 2.4e-7 | 2.5e-4 | 2.5e-4 |
| Product | 7 | 1.5e-9 | 6.2e-5 | 6.2e-5 |
| GradDot | 1 | 9.3e-1 | 3.9e0 | 3.0e0 |
| GradDot | 3 | 4.9e-3 | 7.3e-2 | 7.8e-2 |
| GradDot | 5 | 3.3e-5 | 6.8e-3 | 6.7e-3 |
| GradDot | 7 | 2.5e-7 | 1.5e-3 | 1.5e-3 |

- **(a) = 0 in every case.** The split is exactly consistent.
- **(b) is linear in c** (ratio 10.00 for 10c). At fixed ℓ′max it decays by about 0.07 (Product) to 0.08 (GradDot) per unit n.
- **(c) dominates the total at ℓ′max = 12.** This is the ordinary power-law convergence of the residual-mode coupling sum, which a real-axis scheme would also have. At n = 3 the total goes 5.8 → 2.1e-3 (Product) and 3.2 → 7.8e-2 (GradDot) as ℓ′max goes 2 → 12.
- **(b) does not converge in ℓ′.** The per-mode contributions *grow* with ℓ′. Near ℓ′ = 30 the growth factor per mode is:

  | worldtube radius R | Product | GradDot |
  |---|---|---|
  | ½ | 1.13 | 1.19 |
  | ¼ | 1.06 | 1.09 |

  These are close to 1 + R/r₀ (1.10 and 1.05) times a slowly varying power of ℓ′. At n = 3 and ℓ′max = 30, the cumulative (b) is 5.6e-4 (Product) and 1.4e-1 (GradDot) with R = ½, and 5.9e-6 and 2.7e-3 with R = ¼.

  Interpretation, not yet confirmed by a dedicated test: the Θ→H continuation builds each mode's shadow from the inside and outside analytic pieces (∝ r^ℓ′ and r^{−ℓ′−1}) continued across r₀. Across the tube these grow like (1 + R/r₀)^ℓ′.

### Stage 7: optimum ℓ′max and a shrinking worldtube

`stage7_lprime_optimum.wls` (≈ 4 min) uses the Stage 6 puncture with n ∈ {3, 5}, ℓ = 0…6 and ℓ′max ∈ {4, 8, 12, 16, 20, 30, 40}. The tube radius is R = 1/2, R = 1/4, or R = r₀/ℓ′max.

**Method.**
- Total exterior error = (c) + (b), from the Lemma-1 decomposition of Stage 6:
  - (c), the ℓ′ truncation on the real line, does not depend on R. It is computed by Green's-function quadrature of (truncated − exact) SR+RR sources.
  - (b), the shadow, is computed by quadrature over [r₀ − R, r₀ + R].
- No contour solve is needed for these numbers, so the contour amplitude does not enter. For the three-region verification runs I used amp = 2R.
- Errors are scaled and maximised over r ∈ [1, 3.5] ∪ [6.5, 10] and over ℓ.

**Verification (12 checks pass).**
- Quadrature converged to 1e-26.
- The fast (b) equals Stage 6's decomposition.
- (c) and (b)+(c) agree with full Stage 6 runs (contour and real-line solves) at (R, ℓ′max) = (1/2, 8) and (r₀/12, 12) to ≤ 1e-14 relative, for both cases.

**Results** (total; `results/stage7_summary.md` also lists (b) and (c) separately; figure `../figures/fig_lprime_optimum.pdf`):

| case | n | R = 1/2: argmin, min | R = 1/2 at ℓ′max = 40 | R = 1/4 at 40 | R = r₀/ℓ′max at 40 |
|---|---|---|---|---|---|
| 1/ρ² | 3 | 20, 5.7e-4 | 3.5e-3 | 3.1e-5 | 1.0e-5 (monotone) |
| 1/ρ² | 5 | 30, 4.2e-6 | 1.7e-5 | 8.3e-8 | 1.1e-7 (monotone) |
| 1/ρ⁴ | 3 | 20, 0.13 | 2.8 | 3.1e-2 (min 2.6e-2 at 30) | 8.0e-3 (monotone) |
| 1/ρ⁴ | 5 | 20, 1.4e-3 | 1.6e-2 | 2.5e-5 | 6.5e-5 (monotone) |

**Fixed R.**
- (c) falls as a power law: for 1/ρ², n = 3, it is 3.7e-3, 3.0e-4 and 1.0e-5 at ℓ′max = 12, 20 and 40.
- (b) grows: 8e-5, 2.7e-4 and 3.5e-3 at the same ℓ′max with R = 1/2. Its growth per ℓ′ is roughly geometric, as Stage 6 found.
- So the total has a minimum near ℓ′max ≈ 20–30 for R = 1/2, and grows beyond it.
- For R = 1/4, (b) grows more slowly and stays below (c) up to ℓ′max = 40 except for 1/ρ⁴ with n = 3, where the minimum is at 30.

**Shrinking R = r₀/ℓ′max.**
- (1 + R/r₀)^{2ℓ′max} stays bounded (≈ e²), while (b)'s prefactor falls with R. (b) therefore decreases: for 1/ρ², n = 3, from 1e-3 to 2.5e-7.
- The total decreases monotonically in every case and stays within a few % of (c). In other words, the error is set by the ℓ′ truncation alone.

**Fix made during this stage.** The Stage 6 memo for ∂_r of the residual modes, `dresR`, was keyed without the case parameters. When the Product case evaluated it first, the GradDot case reused φᴿ₂ with c₂ = 0. Stage 7's first verification exposed this as a 15% mismatch.
- The key now includes c₂.
- Stage 6's published numbers are unaffected: there the Product case never evaluates `dresR`, and its re-run reproduces them.

### Stage 8: hybrid split, with no shadow cross-term

`stage8_hybrid_split.wls` (≈ 25 min; 34 checks pass) uses the Stage 6 puncture with R = 1/2, n ∈ {3, 5}, ℓ = 0…6 and ℓ′max ∈ {4, 8, 12, 16, 20, 30, 40}.

**Split.** SR = 2Q[Φᴾ, T_kΦʳᵉˢ] + 2Q[Φᴾ, Φʳᵉˢ − T_kΦʳᵉˢ], where T_k (k = 0, 1, 2) is the Cartesian Taylor polynomial of Φʳᵉˢ about x₀.
- **Taylor coefficients.** They come from mode sums at r₀ of the real-line residual. With n ≥ 3, cqρⁿ adds nothing below order n, so T_kΦʳᵉˢ = T_kΦᴿ exactly.
- **Part 1** uses Gaunt couplings of [Φᴾ]_L with the polynomial's modes (ℓ′ ≤ k), continued as K̃ᴾ (Θ→H) × the entire polynomial modes. Because Im T = 0 on the axis, there is no cross-term.
- **Part 2** is formed on the real line and continued by the Cauchy integral (1/πi)∫_I s(t)/(t−z) dt, with I = [4, 6] split at r₀.
  - The integral uses 80 Gauss nodes per side, after subtracting the cubic Taylor polynomial of s about Re z. That polynomial comes from a Chebyshev interpolant evaluated on the real interval only, and it is integrated analytically.
  - For NDSolve, the contour source is represented by 40 Chebyshev pieces in a. The composition s̃(y(a)) has complex-a singularities about 0.1 from the real segment.
- **Baseline H0.** Cauchy-continue the whole real SR source.
- Part 1 + Part 2 equals the truncated real source for any T̂, so the T̂ error cancels in the total.

**Results** (`results/stage8_summary.md`, figure `../figures/fig_hybrid.pdf`). Scaled exterior error, max over ℓ:

| case | n | route | ℓ′max = 12 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|
| 1/ρ² | 3 | product (Stage 7) | 3.7e-3 | 5.7e-4 | 1.0e-3 | 3.5e-3 |
| 1/ρ² | 3 | H0 = H1 (any k) = (c) | 3.7e-3 | 3.0e-4 | 4.1e-5 | 1.0e-5 |
| 1/ρ² | 5 | product | 4.3e-4 | 1.2e-5 | 4.2e-6 | 1.7e-5 |
| 1/ρ² | 5 | H0 = H1 = (c) | 4.3e-4 | 1.3e-5 | 8.0e-7 | 1.1e-7 |
| 1/ρ⁴ | 3 | product | 0.26 | 0.13 | 0.49 | 2.8 |
| 1/ρ⁴ | 3 | H0 = H1 = (c) | 0.24 | 5.5e-2 | 1.7e-2 | 7.3e-3 |
| 1/ρ⁴ | 5 | product | 2.2e-2 | 1.4e-3 | 2.8e-3 | 1.6e-2 |
| 1/ρ⁴ | 5 | H0 = H1 = (c) | 2.3e-2 | 1.9e-3 | 2.7e-4 | 6.5e-5 |

- **The shadow cross-term is gone.** H0 and H1 match (c) to < 1e-6 relative at every ℓ′max, and the error now decreases monotonically.
- **Cauchy-step error** (contour solve minus real-line solve, scaled). It is ≤ 1e-14 in every case and flat in ℓ′max.
  - H0: 1e-14 (1/ρ²) and 2e-15 (1/ρ⁴).
  - H1: 2e-17 to 2e-15 (1/ρ², k = 1, 2) and up to 6e-15 (1/ρ⁴, k = 1).
  - The subtraction of T_k brings no systematic gain once the Cauchy step itself is accurate.
  - Doubling the resolution changes it by ≤ 7e-18.
- **Taylor coefficients from mode sums** converge slowly, because the ρⁿ modes decay only algebraically. For n = 3 at ℓ′max = 40 the errors are 9e-6 (Φ), 3e-6 (∂_zΦ) and 6e-4 to 9e-4 (second derivatives); for n = 5 they are 8e-8 to 4e-6. This does not affect the total, because T̂ cancels between Parts 1 and 2.
- **Figure note.** In the lower-right panel, H0 and H1 k = 0 coincide for 1/ρ⁴.

**Two failed attempts, for the record.**
1. Subtracting only s(Re z) left an O(Im z) near-singular remainder: Cauchy-step error 1e-7 that did not converge.
2. Evaluating a global Chebyshev interpolant of s at complex z amplified noise by about 4^N, and Mathematica's significance tracking collapsed in the Clenshaw recurrences. Chebyshev arithmetic is now fixed at 40 digits.

### Stage 9: distributional SR, Q = ½∂_z²(AB)

`stage9_dzz.wls` has 326 checks and takes about 6 min. Run it with the argument `smoke` for a quick version. Full tables are in `results/stage9_summary.md`, and the figure is `../figures/fig_dzz.pdf`.
The agent wrote the script and ran the smoke version. It was then cut off by a usage limit, so the full run and this section were done afterwards.

**Setup.**
- With this operator, SR ⊃ qΦ^res ∂_z²(1/ρ) ~ 1/ρ³. That is not integrable as a function, but it is well defined as a distribution.
- Its modes contain δ(r − r₀) terms. The 1/ρ³ part averages to zero over angles, so there is no 1/|r − r₀| term.
- Exact answer: Φ₂ = ½∂_z²χ, where χ is the Stage 3 total.
- Modes of ∂_z come from the ladder ∂_z(fY_ℓ) = a_{ℓ+1}(f′ − ℓf/r)Y_{ℓ+1} + a_ℓ(f′ + (ℓ+1)f/r)Y_{ℓ−1}.

**Test 1: exact puncture, ℓ ≤ 6.**
- SS and SR are recovered on the contour to ≤ 5e-16 scaled.
- Point-supported fit for SS and SR: |α| + |β| ≤ 2e-18, so there are no extra terms.
- The δ coefficients of the SR + RR source modes run from −0.21 (ℓ = 0) to −0.77 (ℓ = 6).
- Dropping them, as pointwise sampling of the real-line source does, gives scaled errors of 0.4 to 500.

**Test 2: Stage 6 puncture, R = ½.** Maximum scaled exterior error over ℓ ≤ 6:

| n | route | ℓ′max = 8 | 12 | 16 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|
| 3 | (c) ℓ′ truncation | 5.7e-1 | 2.2e-2 | 4.8e-3 | 1.6e-3 | 2.1e-4 | 5.0e-5 |
| 3 | product route | 5.8e-1 | 3.4e-2 | 3.3e-2 | 6.0e-2 | – | – |
| 3 | H0 (no split) | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 |
| 3 | H1, k = 0–3 | = (c) | = (c) | = (c) | = (c) | = (c) | = (c) |
| 5 | (c) | 1.7e-1 | 2.6e-3 | 3.2e-4 | 6.8e-5 | 4.1e-6 | 5.6e-7 |
| 5 | product route | 1.7e-1 | 2.4e-3 | 2.5e-5 | 5.2e-4 | – | – |
| 5 | H1, k = 0–3 | = (c) | = (c) | = (c) | = (c) | = (c) | = (c) |

- **H0 fails** (5e2), because the real-line source loses the δ terms. A contour solve confirms the miss as predicted (0.218 at ℓ = 2, ℓ′max = 8).
- **H1 succeeds for every k, including k = 0.**
  - Part 2 carries no δ term (≤ 1e-38).
  - With k = 0 it has only a jump (≤ 0.10); with k = 1 it is continuous; with k ≥ 2 it is C¹.
- **H1 equals (c) at every ℓ′max, although T̂ is inexact.** The Taylor value at x₀ from mode sums is off by 8.5e-4 at ℓ′max = 8. T̂ is computed from the same truncated residual that enters Part 2, so the δ cancels exactly, and the error in T̂ is absorbed into (c).
- **The product route deteriorates sooner than with Q = AB or Q = ∇A·∇B:** its minimum is near ℓ′max = 12–16. The second z-derivative amplifies the high-ℓ′ shadow.

### Stage 9: distributional SR, Q = ½∂_z²(AB)

`stage9_dzz.wls` (≈ 6 min; 326 checks pass). Here ∇²Φ₂ = ½∂_z²(Φ₁²) with Φᴿ = c₀ + c₁z, and the exact answer is Φ₂ = ½∂_z²χ (χ = Stage 3 total).

**Package.** `aZ`, `dzExpr`, `dzzExpr` and `dzzDeltaCoeff` were added, and "Dzz" dispatch for `pieceAnalyticSolution`, `pieceAnalyticSource` and `ssSourceModeExpr`. Every Dzz piece is ½D_z² of the Product piece. The real SS branch writes log|r − r₀| as ½log(r − r₀)² so that it can be differentiated.

**Modes of ∂_z and their δ terms.**
- The ∂_z mode formula matches the closed-form modes of ∂_z(1/ρ) exactly (ℓ ≤ 6).
- D_z² of a kinked mode has a δ(r − r₀) with coefficient a_ℓJ_{ℓ−1} + a_{ℓ+1}J_{ℓ+1}, where J_j is the jump of (D_zF)_j.
- For [1/ρ]_ℓ this coefficient is −4πY_ℓ(0)/r₀², i.e. the full −4πδ³, not the −(4π/3)δ³ in ∂_z²(1/ρ) = PV − (4π/3)δ³. The pointwise mode parts sum to a slab-regularised principal value, PV + (8π/3)δ³ (Stage 11).
- The SR δ coefficient is c_ℓ = −4πqΦᴿ(x₀)Y_ℓ(0)/r₀², i.e. −0.213, −0.368, −0.476, −0.563, −0.638, −0.705, −0.767 for ℓ = 0…6.

**Test 1: exact puncture, ℓ ≤ 6** (`results/stage9_test1.csv`). Errors are scaled exterior errors.
- The SS contour solve equals ½∂_z²(q² ln ρ) to ≤ 2e-16. Fitting the difference at a_R to the modes of 1/ρ and ∂_{r₀}(1/ρ) gives |α| + |β| ≤ 8e-19: there are no point-supported terms.
- The SR product route (Θ→H continuation of K times φᴿ) is right to ≤ 3e-16. The continuation carries the δ, and the fitted point terms are ≤ 2e-18.
- The total equals ½∂_z²χ to ≤ 5e-16.
- **Negative control.** The real-line solve of the pointwise SR (δ dropped) is off by 0.41 (ℓ = 0) to 5e2 (ℓ = 6). Adding c_ℓG₀(r, r₀) restores it to ≤ 4e-16.

**Test 2: Stage 6 puncture, R = 1/2** (`results/stage9_summary.md`, `results/stage9_test2.csv`, figure `../figures/fig_dzz.pdf`). Scaled exterior error, max over ℓ = 0…6:

| n | route | ℓ′max = 4 | 8 | 12 | 16 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|
| 3 | (c) = H1, every k = 0…3 | 3.5 | 0.57 | 2.2e-2 | 4.8e-3 | 1.6e-3 | 2.1e-4 | 5.0e-5 |
| 3 | product route | 3.5 | 0.58 | 3.4e-2 | 3.3e-2 | 6.0e-2 | – | – |
| 3 | H0 (no split) | 5.1e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 |
| 5 | (c) = H1, every k | 3.1 | 0.17 | 2.6e-3 | 3.2e-4 | 6.8e-5 | 4.1e-6 | 5.6e-7 |
| 5 | product route | 3.1 | 0.17 | 2.4e-3 | 2.5e-5 | 5.2e-4 | – | – |
| 5 | H0 (no split) | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 | 5.0e2 |

- **H0 fails, as expected.** Cauchy-continuing a pointwise function adds no δ, so H0 misses c_ℓ^{(M)}G₀(r, r₀). That is O(1) relative per mode (5e2 scaled at ℓ = 6, where the solution is small).
- **H1 succeeds for every k, including k = 0.** The mode-level δ of 2Q[Φᴾ, f] depends only on the mode sum f̂(x₀): it is −4πqY_ℓ(0)f̂(x₀)/r₀². T_k reproduces that value exactly for every k (A₀ = Φ̂_M(x₀)), so Part 2's δ coefficient is < 1e-38 relative for all k. I expected k ≥ 2 to be needed; that is true for smoothness only.
  - Part 2 jumps at r₀ for k = 0 (0.10).
  - For k = 1 it is continuous, with a derivative jump of 5e-7 (n = 3), which comes from the T̂ error.
  - For k ≥ 2 it is C¹ (< 1e-36).
- **δ coefficients.** At ℓ′max = 40 the truncated SR δ differs from the exact c_ℓ by 6e-6 relative (n = 3) and 5e-8 (n = 5); this is the error of Φ̂_M(x₀). The truncated RR has no δ (< 1e-49 relative).
- **Product route.** It equals (c) plus the shadow term, with the Stage 7 minimum near ℓ′max ≈ 12–16.
- **(c) is large at small ℓ′max** (3.5 at ℓ′max = 4) because D_z² amplifies the truncation error of the pointwise source.
- **Validation by contour solves** (n = 3, ℓ = 2, ℓ′max = 8), comparing contour minus distributional real line at a_R with the prediction:
  - H0 matches −c^{(8)}{G₀, ∂_rG₀} to 6e-19 relative.
  - H1 with k = 0 and k = 2 gives 5e-19 (predicted 0).
  - The product route at ℓ′max = 4 gives 5.13e-7 against an axis shadow term of 5.13e-7 (1e-14 relative).

### Stage 10: contour against the mode-level Hadamard finite part (Musaeus, Pound & Upton, Sec. V.D)

`stage10_teukolsky_benchmark.wls` (≈ 40 s; 85 checks pass) solves [r²∂² − 2r∂ − (ℓ−1)(ℓ+2)]ψ = S with the Green's function of their eq. (140), r₀ = 5, ℓ = 2…6.
- **Sources.** A/Δr⁴, B/Δr³, C/Δr² and D/Δr one at a time, and a mixed set (A, B, C, D) = (0.7, −1.3, 2.1, 0.4) plus a smooth term 1.5/(1 + (r/r₀)²) and a log term −0.8 log|Δr| r₀³/(r₀³ + r³).
- **Observation radii.** r ∈ {1, 2, 3, 7, 9, 12, 20}. Differences are scaled by max|ψ| over these radii (`results/stage10_summary.md`).

**Results.**
- **Eq. (140) is correctly normalised.** Both branches solve the homogeneous equation (r^{ℓ+2} at 0, r^{1−ℓ} at ∞), G is continuous, and r′²[∂_rG] = 1. So L G = δ for the operator as written.
- **P1** is eq. (141) with (A2)–(A5); I read ∂_{r₀}G(r, r₀ ± a) as ∂_{r′}G at r′ = r₀ ± a. It agrees with **P1′**, an independent Hadamard finite part by Taylor subtraction, to < 1e-38. The transcribed boundary terms are therefore right.
- **Cutoff independence.** Results for a ∈ {1/4, 1/2, 1} agree to ≤ 4e-34.
- **P2, eq. (142).** It agrees with P1 to < 1e-38 for each of A, B, C and D, at all seven radii and on both sides of r₀: **no typo**.
  - The rendered PDF garbles the polynomial coefficients (glyphs overlap). The text layer reads 12r³ + 6r²r₀ + 4rr₀² + 3r₀³, −12r⁴ + 6r³r₀ + 2r²r₀² + rr₀³ + r₀⁴ and 84r⁵ − 126r⁴r₀ + 28r³r₀² + 7r²r₀³ + 2rr₀⁴ + 3r₀⁵, and these are correct.
- **P3, the contour.** Re∫_C G(r, z)S̃(z)dz, plus the real line outside the tube, equals P1 to ≤ 4e-34 for (R, amplitude) = (1/2, 1), (1/2, 1/4), (1/2, 2), (1/4, 1) and (1, 1). The log term is continued as log(z − r₀).
- **P3b, a contour ODE solve.** It runs on the real line [1, a_L], the contour and the real line [a_R, 9], with Cauchy data at r = 1 from P1. The relative error is 2e-16 at r = 7 and 8e-16 at r = 9; the shadow Im ψ(a_R) = 0.54 is harmless.
- **Why they agree.** For an integer pole, Re∫_C f/(z − r₀)ⁿ dz equals the Hadamard FP, because the residue term is imaginary for real f. Over a symmetric interval |Δr̃|^B has no pole at B = 0, so the Riesz and Hadamard finite parts coincide. For these sources the contour is the paper's mode-level prescription, with no cutoff and no boundary terms.
- **Not done:** the optional static-particle source −₂T_ℓm.

### Stage 11: their T_Q in the ∂_z² toy

`stage11_TQ.wls` (≈ 2 min; 73 checks pass) uses the exact puncture q/ρ, Φᴿ = c₀ + c₁z and SR = q∂_z²(Φᴿ/ρ), the Stage 9 SR, with ℓ ≤ 6 (`results/stage11_summary.md`).

**(a) Toy T_Q.** This part pairs SR with test functions φ = g(r)P_ℓ(cos θ), where g is a Gaussian shell at r = 5.2 of width 0.8 and ℓ = 0…6. Three pairings are computed independently:
- **canonical:** ⟨qΦᴿ/ρ, ∂_z²φ⟩, by modes;
- **ball PV**, which is the 4D Hadamard Pf, because the 1/ρ³ term has zero spherical mean and so (ρ/ρ₀)^B has no pole. It is computed as a 3D integral in coordinates centred on x₀, doing the angular integral first.
- **slab:** Σ_ℓ∫(pointwise mode)·φ_ℓ r²dr. This is what the mode-level finite part of eq. (105) gives for these bounded modes.

The three pairings give:
- canonical − ball PV = −(4π/3)qΦᴿ(x₀)φ(x₀), to 1e-20 relative for every test function. So T_Q = −(4π/3)qΦᴿ(x₀)δ³, as expected.
- slab − canonical = +4πqΦᴿ(x₀)φ(x₀), so slab = Pf + (8π/3)qΦᴿ(x₀)δ³.
- Analytically, with both regulators the 1/ρ³ term gives φ(x₀)qΦᴿ(x₀) ε^{B+B′}/(B+B′) · 8πB′/((B′+1)(B′+3)). Taking FP in B′ and then B → 0 gives 0 (ball). Setting B = 0 first gives 8π/3 (slab). The two orders do not commute.

**(b) The Stage 9 coefficient c_ℓ = −4πqΦᴿ(x₀)Y_ℓ(0)/r₀² splits as (2/3)c_ℓ + (1/3)c_ℓ.** The δ content of the modes of the 4D Pf is (2/3)c_ℓ, and that of T_Q is (1/3)c_ℓ. The measured κ_PV and κ_TQ equal 2/3 and 1/3 to 1e-20 for every ℓ.

**(c), (d) Prescriptions end to end.** Cauchy data are taken at r = 1. The SS mode FP subtracts A/Δr² + B/Δr (Laurent coefficients, the same on both sides) and uses (A2)–(A3) with the Coulomb-mode Green's function. Scaled exterior error against ½∂_z²χ:

| route | ℓ = 0 | 2 | 4 | 6 |
|---|---|---|---|---|
| contour (Stage 9 sources, Green's-function form) | < 1e-30 | < 1e-30 | < 1e-28 | < 1e-27 |
| Pf modes (pointwise + (2/3)c_ℓδ) + T_Q, i.e. eq. (24); a = 1/2 or 1/4 | < 1e-20 | < 1e-19 | < 1e-18 | < 1e-17 |
| Pf modes, T_Q omitted | 0.14 | 1.2 | 41 | 1.7e2 |
| pointwise modes (eq. 105) + T_Q | 0.27 | 2.3 | 83 | 3.4e2 |
| pointwise modes, no δ at all | 0.41 | 3.5 | 1.2e2 | 5.0e2 |

Each failing row equals the response to the missing (1/3), (2/3) or full c_ℓδ; this was checked to 1e-8 for the "T_Q omitted" row.

**Reading.**
- Eq. (24) (Pf + T_Q) is exact in the toy.
- The mode-level form of Pf in eqs. (104)–(105) is exact only if the source modes carry the Pf's own δ content, (2/3)c_ℓ. Mode functions computed pointwise (by products of mode functions, away from r₀) do not carry it. With them, the mode-level FP returns the slab-regularised distribution, and adding T_Q leaves an error of (2/3)c_ℓδ.
- The cause is dropping (ρ/ρ₀)^B after inserting |Δr̃|^{B′} (eq. 104 → 105). That step is not innocuous for a 1/ρ³ term with zero spherical mean but nonzero |cos|^{B′}-weighted mean.
- This is a toy-model statement. Whether it matters in GR depends on how δ²G_ℓm is assembled, i.e. whether its δ(r − r₀) terms are kept.
- The contour needs neither T_Q nor any δ bookkeeping: it reproduces the canonical distribution.

### Stage 12: SS from numerical data on the contour, and H1 convergence (Dzz)

`stage12_ss_reuse_convergence.wls` (≈ 17 min; 65 checks pass). Exact puncture, ℓ = 0…6, R = 1/2, amplitude 1. Errors are scaled exterior SS errors, maximised over ℓ, against the closed-form contour result. The closed-form contour result matches the exact ½∂_z²(q² ln ρ) to < 1e-27.

**How the exterior is computed.** The exterior comes from the contour Green's-function quadrature Re∫_C G(r, z)S̃(z)dz, which is equivalent to the contour solve by Lemma 1. A Lemma-1 cross-check against the real-line integral of (interpolant − remainder) agrees to 1.5e-6 relative; that is how well the contour quadrature resolves the hybrid error.

**(A) Reusing the hybrid technique for SS.** The singular part, from the local expansion SS_ℓ = A/Δr² + B/Δr + (C₀ + C₁Δr + C₂Δr²)ln|Δr| + analytic, is continued analytically. The coefficients agree on both sides of r₀ to 1e-25. The remainder is sampled at cell centres of spacing h on I = [4, 6] and treated as data:
- it is represented by local degree-7 Lagrange pieces, one-sided at r₀;
- it is Cauchy-continued with the exact cell formula for |ζ| ≤ 3 and 16-point Gauss elsewhere.

One numerical point matters: the near cells must be *masked out* of the Gauss sum, not subtracted afterwards. Near the tube ends the contour is within 1e-17 of the axis, and subtracting caused catastrophic cancellation, giving erratic 1e-9 errors.

| subtracted | h = 1/50 | 1/100 | 1/200 | order |
|---|---|---|---|---|
| {A} (negative control) | 1.7 | 0.87 | 0.44 | 1 |
| {A, B} | 1.2 | 0.62 | 0.31 | 1 |
| {A, B, C₀} | 2.0e-5 | 2.5e-6 | 3.1e-7 | 3 |
| {A, B, C₀, C₁, C₂} | 3.2e-9 | 9.9e-11 | 3.2e-12 | 5 |

- **Subtracting only the Hadamard pair {A, B} (eq. 109 of Musaeus et al.) is not enough here.** The C₀ ln|Δr| left in the data cannot be interpolated by polynomials, so the error falls only as O(h). With the log term subtracted analytically, the error is O(h³) for {A, B, C₀} and O(h⁵) when C₁ and C₂ are also subtracted.
- **Noise** (`results/stage12_ss_noise.csv`), {A, B, C₀, C₁, C₂}, noise-induced error at h = 1/50, 1/100, 1/200:
  - relative noise ε on the remainder samples: ≈ ε (1.1, 0.75, 0.55) × 10;
  - relative noise ε on the full SS samples (the realistic case, as from Upton's 8-digit angular integrals): ≈ ε (1.9, 9.0, 7.2) × 10⁴. Subtracting A/Δr² amplifies the noise near r₀ by ~A/h².

**(C) Radial factor, no subtraction.** F = (r − r₀)²SS is sampled with relative noise. It uses a continuous piecewise interpolant (pieces between samples) and is Cauchy-continued to F̃. Then F̂ = F̃ − iG₀ − iG₁(z − r₀) and s̃ = F̂/(z − r₀)².

| h | 1/50 | 1/100 | 1/200 |
|---|---|---|---|
| unpinned | 23 | 23 | 23 |
| pinned, G₀ and G₁ from the data | 8.0 | 4.0 | 2.0 |
| pinned, exact G₀ and G₁ | 3.0e-5 | 3.8e-6 | 4.7e-7 |
| error of G₁ (data) | 2.2e-3 | 1.1e-3 | 5.4e-4 |
| noise-induced error, ε = 1e-8 (pinned) | 6.3e-5 | 8.1e-4 | 1.3e-3 |

- **Without pinning** the error is O(1). It equals the response to −πG₀δ′ + πG₁δ to 1.6e-14.
- **With pinning** at the exact G₀ and G₁, (C) is as good as {A, B, C₀}: O(h³).
- **G₁ is the weak point.** G₁ = Im F̃′(r₀+i0) is a Hilbert-transform derivative *at* r₀, where F has a Δr² ln|Δr| term. Any polynomial representation of the data gets it wrong by O(h ln h), and the pinned error inherits this (first order in h). G₀ is fine (O(h³)).
- **Noise.** The noise is amplified ∝ 1/h through G₁, not ∝ 1/R² as expected. Unpinned, the amplification would be bounded by 1/R², but unpinned is O(1) wrong.
- **Conclusion.** (C) avoids the A/h² noise amplification of subtraction but needs an accurate G₁, i.e. local knowledge of the ln|Δr| structure after all.

**(D) Complex-r angular quadrature.** SS_ℓ(z) = √(π(2ℓ+1)) ∫ ½q² ∂_z²(1/ρ²) P_ℓ dμ is evaluated at the contour points directly. The θ path is deformed by mapping μ → ρ: the integrand becomes rational in ρ, with its only pole at ρ = 0. A ρ path from z + r₀ to r₀ − z via −i/2 never meets ρ = 0, so it gives the continuation from a_L; ending at z − r₀ via +i/2 gives the continuation from a_R. 2 × 80 Gauss nodes suffice.
- **SS.** The continuation from a_R equals the closed form pointwise (1.8e-14). The one from a_L differs pointwise by 2–30%, but by a function that is imaginary on the axis, so the exterior agrees to ≤ 1.7e-12.
- **Cost.** Adaptive μ-quadrature to 1e-12 takes 1100–2000 evaluations at real r between 0.1 and 0.001 from r₀, and 1000–1230 at contour points. On the contour the integrand is smooth, so fixed 160-node Gauss reaches 1e-14.
- **SR (full, exact puncture).** A single angular integral gives one analytic branch, in or out. The in branch alone is wrong by 4 to 1.2e3. Blending the differentiated in/out SR with H is also wrong (3 to 8.8e2), because the δ terms come from derivatives of H.
  - What works is to apply the Θ→H blend to the *undifferentiated* product modes [2qΦᴿ/ρ]_j (from in/out angular quadrature with r-derivatives) and then apply D_z². This reproduces the product route to 3.6e-14 pointwise.
  - SR Part 1 of H1 (Stage 6 puncture, n = 3, k = 2): the same procedure agrees with the closed form to 1.1e-13.
- **Continued residual modes of the Stage 6 puncture (n = 3) at complex z.**
  - The terms shrink until ℓ′ ≈ 25–30 and then grow by a factor 1.05–1.07 per ℓ′. Asymptotically the growth is max(|z|/r₀, r₀/|z|) ≤ 1 + R/r₀, times algebraic decay, so Σφ̃^res_ℓ′(z)Y_ℓ′(θ) diverges.
  - The terms are still ≤ 7e-5 at ℓ′ = 80, and the θ = 0 partial sums are flat to three digits; the divergence would show only beyond ℓ′ ~ 10².
  - At the top of the contour (|z|/r₀ = 1.02) the terms are still decreasing at ℓ′ = 80.

**(B) H1 (k = 0) convergence, Stage 6 puncture** (`results/stage12_B_summary.md`, figure `../figures/fig_hybrid_convergence.pdf`). The H1 error is (c); a contour solve for the new case n = 1 confirms this to 8e-17.

| n | slope of error (ℓ′max 20–40) | per-mode power (ℓ′ 25–40) | error at ℓ′max = 12 | 20 | 40 |
|---|---|---|---|---|---|
| 1 | 2.99 | 4.07 | 0.56 | 0.11 | 1.4e-2 |
| 3 | 4.95 | 6.05 | 2.2e-2 | 1.6e-3 | 5.0e-5 |
| 5 | 6.91 | 8.04 | 2.6e-3 | 6.8e-5 | 5.6e-7 |
| 7 | 8.87 | 10.0 | 6.5e-4 | 6.2e-6 | 1.3e-8 |

- **Rates.** The error falls as ℓ′max^−(n+2) and the per-mode contributions as ℓ′^−(n+3), as expected for a C^{n−1} residual. The contributions are smooth in ℓ′, with no parity oscillation.
- **Error against r.**
  - Inside the tube radius the error is largest next to the tube: for n = 3 at ℓ′max = 20, 4.9e-10 at r ≤ 3.5 and 2.9e-6 for r in [3.5, 4.5].
  - Outside, it grows outwards to r = 10, because Cauchy data are imposed at r = 1 and the error propagates as the growing homogeneous solution: 1.2e-4 for r in [5.5, 6.5] and 1.6e-3 at r = 10.
- **Tail acceleration.** Fit the per-mode contributions on ℓ′ ∈ [ℓ′max−5, ℓ′max] with Upton's Λ_{ℓ′,0,−(n+3)}, Λ_{ℓ′,0,−(n+5)} (their eq. 173), and add the tail using Σ_{ℓ′≥0}Λ = 0. Gain at ℓ′max = 20 / 30:

| n | 1 | 3 | 5 | 7 |
|---|---|---|---|---|
| 2-term fit | 38 / 4.3e2 | 48 / 6.5e2 | 67 / 1.1e3 | 1.1e2 / 2.0e3 |
| 1-term fit (last 3 modes) | 11 / 30 | 13 / 38 | 16 / 50 | 20 / 69 |

  - At ℓ′max = 12 the window is pre-asymptotic: the 2-term fit makes things worse (gain 0.02–0.04), while the 1-term fit gains 1.5–2.1.
  - Plain power laws ℓ′^−p in place of Λ gain 3–16 times less.

### Stage 13: PRL figures (Dzz)

`make_prl_figures.wls` (≈ 2 min; 9 checks pass) writes to `../figures/prl/`. The figures are vector PDF in Times 8 pt with no titles, use consistent colours (SS blue, SR orange; routes: product red, H0 grey, H1 blue, H1 + tail green, (c) black dashed), and each has a CSV beside it.
- **`prl_source.pdf`** (17.8 cm).
  - (a) ℓ = 2 SS and SR modes on the real line, with A₂/Δr² dashed and an arrow marking the SR δ: c₂ = −0.476.
  - (b) Re and Im of the continued SS and SR against a.
- **`prl_solution.pdf`** (8.6 cm).
  - (a) Recovered exterior Φ₂,ℓ for ℓ = 0, 2, 4 (contour solves of SS and SR, RR on the real line) against ½∂_z²χ.
  - (b) Scaled error against r for ℓ = 0–6, at ≤ 1e-15.
- **`prl_convergence.pdf`** (8.6 cm). n = 3 and 5: the product route (R = ½), H0, H1 (k = 0), H1 + Λ tail, and (c).
  - H1 and (c) coincide.
  - The tail curve uses the robust 1-term fit and is drawn from ℓ′max = 12, where the fit window becomes asymptotic.

### Stage 14: cost of the SS angular integrals, real axis against contour

`stage14_cost.wls` (≈ 23 min; 9 checks pass) computes SS_ℓ = √(π(2ℓ+1)) ∫ ½q²∂_z²ρ⁻² P_ℓ dμ for ℓ ≤ Lmax.
- **Cost measure.** It counts every evaluation of the 3D integrand S (P_ℓ evaluations are not counted) and records wall time.
- **Error measure.** Pointwise relative error against the closed form evaluated at 40 digits, maximised over ℓ and over points: 300 contour nodes, both branches, for the contour; the 400 points r₀ ± kR/200 for the real axis.
- **Contour machinery.** The Stage 12 (D) machinery is used unchanged.
- **Radial nodes (E).** For the exterior SS field the contour needs 64 nodes (4 panels × 16) for 1e-8 and 96 (4 × 24) for 1e-12, against the 400-point real grid.

| Lmax | target | (A) real, adaptive per (ℓ, r) | (D) real, shared nodes in ln ρ | (B) contour, per ℓ | (C) contour, shared nodes |
|---|---|---|---|---|---|
| 10 | 1e-8 | 1.16e6 (2e-11), 53 s | 9.7e3, 0.05 s | 8.6e4, 1.2 s | **8.2e3**, 0.03 s |
| 10 | 1e-12 | 2.44e6 (floor 2.4e-11), 216 s | 1.28e4, 0.06 s | 1.66e5, 2.3 s | 1.54e4, 0.05 s |
| 20 | 1e-8 | 2.56e6, 101 s | 1.65e4, 0.08 s | 1.58e5, 3.0 s | **8.2e3**, 0.03 s |
| 20 | 1e-12 | 7.41e6 (floor 6.8e-11), 412 s | 1.90e4, 0.09 s | 3.06e5, 6.0 s | 1.54e4, 0.05 s |
| 30 | 1e-8 | 5.37e6 (7e-7: PG 8 not met everywhere), 149 s | 2.29e4, 0.11 s | 2.30e5, 5.7 s | **8.2e3**, 0.03 s |
| 30 | 1e-12 | 4.97e7 (floor 5.4e-10), 608 s | 2.51e4, 0.12 s | floor 1.6e-12 | floor 1.6e-12 |

Each entry is the number of S evaluations, the error reached if it misses the target, and the wall time. The real axis uses 400 radial points and the contour 64 (1e-8) or 96 (1e-12) nodes.

**Per point.**
- **(A)** needs about 450–500 S evaluations per (ℓ, r) at PrecisionGoal 8 and 3400–4400 at PrecisionGoal 12. This is nearly flat in |Δr| from h to R, because the adaptive rule resolves the peak near μ = 1. At machine precision it does not reach 1e-12.
- **(D)** needs nq = 56–80 per point for 1e-8 and 1e-12, flat in |Δr|. Gauss in plain ρ fails within about 10h of r₀.
- **Contour (B, C).** It needs nq = 64 per segment for 1e-8 and 80 for 1e-12 (2nq S evaluations per node). The cost is set by the first segment of the Stage 12 path, which passes within a normalised distance 0.1 of the pole ρ = 0; a graded path would lower it.
- **Shared-node count.** With shared nodes the S count per node is independent of Lmax; this is checked.
- **Precision floor.** At Lmax = 30 every method has a machine-precision floor near 1e-12 (contour 1.6e-12).

**Reading.**
- Sharing the S evaluations across ℓ is the main saving: 10²–10³× fewer evaluations than adaptive integration per (ℓ, r), on the real axis or on the contour.
- With shared nodes on both sides, the contour needs 1.2–2.8× fewer S evaluations than the real axis at 1e-8 (growing with Lmax), and about the same number at 1e-12 (0.8–1.2×).
- That saving comes from radial points (64–96 against 400), not from cheaper angular integrals: per point the contour needs 128–160 evaluations against 56–80.

### Stage 14b: radial source samples for the exterior SS field

`stage14b_radial.wls` (≈ 8 min; 7 checks pass) uses the Stage 14 (E) error measure: the tube contribution to the exterior SS field, maximised over ℓ ≤ Lmax.
- **Real-axis method.** The real axis computes the finite part of ∫G·SS over the tube. The local terms A/Δr² + B/Δr + C₀ ln|Δr| are integrated against G analytically; their coefficients come from a 150-digit fit and agree with the Stage 12 series. The remainder is integrated from samples.
- **Samples.** SS_ℓ(r) for all ℓ comes from Legendre Q_j (17 ms per point at 40 digits), checked against the closed form.
- **Cross-check.** The finest R2 matches the contour reference to 9e-14 (finite part = contour).

| method | samples for 1e-8 | for 1e-12 | error with 1e-8 noise on SS (at the 1e-8 / 1e-12 configuration) |
|---|---|---|---|
| contour, Gauss panels in a | 64 | 96 | 1.0e-6 / 2.2e-7 |
| R1a: uniform grid r₀ ± kh, Gregory (6 differences), {A, B, C₀} | 1600 | not reached (4e-11 at 6400) | 2.5e-3 |
| R1b: same, also C₁Δr ln and C₂Δr² ln subtracted | 100 | 800 | 2.3e-4 / 4.0e-4 |
| R2: Gauss panels graded toward r₀ (ratio 0.2), {A, B, C₀} | 48 | 120 | 4.0e-3 / 6.3e-2 |

The sample counts are the same for Lmax = 10, 20 and 30.

**Reading.**
- **Noise-free, the real axis can match the contour on points.** A graded real grid needs 48–120 samples against 64–96, and a uniform grid 100–800 if the C₁ and C₂ log terms are also subtracted. Only the plain uniform grid with {A, B, C₀} is far behind (O(h²), 1600 samples).
- **With 1e-8 relative noise on the samples, the real-axis finite part is 10²–10⁵ times worse than the contour.** Subtracting A/Δr² from noisy samples amplifies the noise as A/Δr². The graded grid, which samples closest to r₀, is hit hardest (6e-2). The contour samples stay at |z − r₀| ≳ R and only lose a factor of about 20–100.

### Stage 12b: hybrid convergence with physical boundary conditions

`stage12b_physical_bc.wls` (≈ 11 min; 38 checks pass) replaces the Cauchy data at r = 1 with ψ′ = ℓψ/r at r = 1 and ψ′ = −(ℓ+1)ψ/r at r = 10.
- **Method.** For a source on [1, 10] these conditions give the free-space Green's function G = −(s²/(2ℓ+1)) r_<^ℓ/r_>^{ℓ+1}. Every exterior error is the quadrature of its error source against G (the Stage 9/12 decomposition with G in place of the Cauchy-data Green's function).
- **Exact reference.** ½∂_z²χ plus A r^ℓ + B r^{−ℓ−1}, chosen to satisfy the same conditions.
- The figure `../figures/fig_hybrid_convergence.pdf` and `results/stage12b_routes.csv` (PRL Fig. 5) now come from here. The Stage 12 Cauchy-data version is kept as `fig_hybrid_convergence_cauchy.pdf`.

| quantity | Cauchy data at r = 1 (Stage 12) | physical conditions |
|---|---|---|
| error slopes, n = 1, 3, 5, 7 | 2.99, 4.95, 6.91, 8.87 | 2.93, 4.87, 6.83, 8.79 |
| per-mode powers | 4.07, 6.05, 8.04, 10.0 | 4.00, 6.01, 8.02, 10.1 |
| 2-term Λ tail gain, ℓ′max = 20 / 30 | 38–113 / up to 2030 | 18–40 / up to 770 |
| 1-term Λ tail gain, ℓ′max = 12 | 1.5–2.1 | 2.5–4.3 |
| H1 at ℓ′max = 12 / 20 / 40, n = 3 | 2.2e-2 / 1.6e-3 / 5.0e-5 | 3.9e-5 / 3.3e-6 / 1.1e-7 |
| H0 (no split) | about 5e2 | about 3.3 |
| product route (n = 3, 5) | minimum near ℓ′max 12–16 | minimum at 12 (1.1e-4) and 16 (1.3e-6), then growing |
| error against r | grows towards r = 10 | largest next to the tube, decaying outwards (n = 3, ℓ′max = 20: 1.0e-5, 3.8e-6, 2.0e-6) |
| Cauchy step of the contour solve | < 1e-14 | ≤ 1.3e-16 |
| split against truncation | < 1e-6 relative | δ level < 1e-25; with the Cauchy step ≤ 7e-12 relative |

### Physical boundary conditions for the PRL numbers

Every exterior result in the PRL now uses ψ′ = ℓψ/r at r = 1 and ψ′ = −(ℓ+1)ψ/r at r = 10, i.e. the free-space Green's function, with the exact reference corrected by the homogeneous terms that meet the same conditions. The Cauchy-data results are kept (`*_cauchy` figures, unflagged results).

| script / flag | outputs |
|---|---|
| `make_prl_figures.wls` | `prl_solution.pdf`, `results/stage9b_exact_physical.csv` (contour solves from zero Cauchy data plus αr^ℓ) |
| `stage7b_physical_bc.wls` | `fig_lprime_optimum.pdf`, `results/stage7b_*` |
| `stage12_ss_reuse_convergence.wls phys` | `fig_ss_reuse.pdf`, `results/stage12phys_*` |
| `stage11_TQ.wls phys` | `results/stage11phys_*` |
| `stage14b_radial.wls phys` | `fig_radial.pdf`, `results/stage14bphys_*` |
| `stage12b_physical_bc.wls` | `fig_hybrid_convergence.pdf`, `results/stage12b_*` (Fig. 5 data) |

The changed numbers (old → new):
- **Exact puncture.** Exterior error < 1e-15 → < 3e-18. Point-supported fit < 2e-18 → < 1.1e-18. Pointwise SR without δ 0.4–500 → 0.55–3.3.
- **App. D.** Growth per mode 1.13/1.19 → 1.12/1.19 at R = 1/2 (1.06/1.09 unchanged at R = 1/4). Fall per unit n 0.070/0.079 is unchanged. The optima and the monotone shrinking-R curves are unchanged, at about 400× smaller errors.
- **App. E.**
  - Subtraction {A} 0.44 → 8.8e-4 and {A, B} 0.31 → 7.7e-4, both still O(h). {A, B, C₀} 3e-7 → 6.5e-10; five terms 3e-12 → 7.9e-15.
  - Noise 7e-4 → 4.1e-5, with amplification 10⁴–10⁵ → 3×10²–4×10³.
  - Unpinned 23 → 5.3, matching the prediction to 2e-15. Route 1 exterior 1.7e-12 → 2.0e-14.
- **App. F.**
  - Contour 64/96 → 64/80 nodes. Uniform grid with 3 terms 1600 → 800 for 1e-8 (6400 give 1e-11). Graded grid 48/120 and uniform with 5 terms 100/800 are unchanged. Finite part against contour 9e-14 → 2e-14.
  - Errors with noise: contour 1e-6/2e-7 → 8e-8/5e-8; graded 4e-3/6e-2 → 4e-4/1.5e-2; uniform with 3 terms 3e-3 → 9e-5; uniform with 5 terms 2e-4/4e-4 → 1e-5/9e-5.
- **App. H.** <1e-27 → <1e-29; <1e-17 → <1e-19; 0.14–170 → 0.18–1.1; 0.27–340 → 0.37–2.2; 0.41–500 → 0.55–3.3.

## Decisions and numerical notes

1. **φˢ_ℓ = ψ_p − D_ℓ, as specified.**
   - D_ℓ is fitted exactly from the canonical-log form at r = 2, 3.
   - It differs from the note's alternative route, G̃_ℓ + (i/π)χ_ℓ, by imaginary homogeneous terms. For example, φˢ₀ contains an extra (i/π)(75/2)/z.
2. **S_eff is not set to zero by hand.** It is evaluated with the notebook's formula, which gives < 1e-65.
3. **Exterior SR and RR sources.**
   - The exterior SR source uses φᴿ from the first-order three-region solve, regions 1 and 3.
   - The contour SR source uses the Stage 2 contour solution.
   - RR uses the real-line through-the-tube φᴿ.
4. **Contour sources are evaluated at numerical y(a).** Building `expr[y(a)]` symbolically made Mathematica expand the polynomials into degree-72 polynomials in a. That lost every digit: the first attempt returned `0``-18`.
5. **NDSolve settings:** WorkingPrecision 30, PrecisionGoal 20, AccuracyGoal 35, `InterpolationOrder -> All`. Each setting was chosen after a measured problem:
   - Default dense output gave 1e-12 interpolation noise between steps for SR modes with ℓ = 2, 3 at r < r₀, although the step values were good to 1e-21.
   - AccuracyGoal 22 cost relative accuracy for ℓ ≥ 7 modes, which are about (r/r₀)^ℓ small at r = 1. The scaled error was 2e-9 at ℓ = 12.
   - AccuracyGoal ∞ or 40 broke NDSolve when a component is exactly zero, for example φ′(1) = 0.
   - PrecisionGoal 25 at WorkingPrecision 30 gave "error test failure" on the contour.
6. **Matching at the tube ends.** y′ = 1 and y″ = 0 at a_L and a_R, so d/da = d/dz there.
7. **Q = ∇A·∇B (Stage 5).** It is implemented with the same functions, dispatched on Q = "GradDot". The regular-field inputs are pairs {φ, ∂_rφ}; on the contour ∂_rφ̃ = (dφ̃/da)/y′(a).
   - ∂_rK and K^θ are continued by Θ→H. Corrected φˢ (its z-derivative for ∂_rK) is used as a cross-check.
   - `phiRAnalytic` now includes c₂ r²P₂. Its default c₂ = 0 keeps Stages 1–4 unchanged.
8. **`csCheck` and inexact zeros (changed during Stage 5).** An inexact zero only certifies |value| < 10^−Accuracy, so a check now passes only if that bound is within the tolerance.
   - This exposed significance-arithmetic pessimism in `Norm` of tiny differences: √ of a quantity known to ±1e-20 is only known to ±1e-10.
   - I therefore switched those check values to Max|·|. Stages 1–4 were re-run under the stricter rule.

## Open questions

- **Second-order shadows are O(0.1–0.8).** They are harmless for the exterior by Lemma 1. They would matter if Φ₂ on the contour were ever multiplied, or if the value at r₀ were needed. Their size and sign depend on the continuation chosen for K.
- **Which K continuation to keep.** For an O(1) shadow (control B), Θ→H was 80 to 300 times less sensitive than the corrected-φˢ continuation. For the truncation-order shadow of Stage 4 the two are comparable. Θ→H is my recommendation (see Stage 4).
- **The notebook's DSolve normalisation for ℓ = 0.** It is the anti-physical Green's function. Any existing results that used `singularModeExact` for ℓ = 0 or ℓ ≥ 2 as "the singular field" carry the Δ above.
