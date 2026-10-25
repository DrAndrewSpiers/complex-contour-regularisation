# Stage 15: extracting the second-order residual field at the particle (Dzz case)

Run: `wolframscript -file stage15_residual_extraction.wls`.  l = 0..20, Nmax = 120, wall time 1343.1 s.

## Setup and formulas

- Toy: Lap Phi1 = -4 pi q delta^3, Phi1 = q/rho + PhiR1, PhiR1 = c0 + c1 Z = A + c1 d with A = c0 + c1 r0, d = Z - r0 (Z = r mu, x0 = r0 zhat); q = 1, r0 = 5, c0 = 1, c1 = 1/10.  Second order: Lap Phi2 = (1/2) d_Z^2 (Phi1^2).
- Exact singular field: chiS = q^2 ln rho + q A rho + (1/2) q c1 d rho, Phi2S = (1/2) d_Z^2 chiS (checked symbolically: Lap chiS = q^2/rho^2 + 2 q PhiR1/rho, using Lap ln rho = 1/rho^2, Lap rho = 2/rho, Lap(d rho) = 4 d/rho).  chiS = chi_SS + chi_SR of Stage 3, so Phi2R = Phi2 - Phi2S = (1/2) d_Z^2 chiRR + h.
- The field solved for is Phi2 = (1/2) d_Z^2 chi + h, with the harmonic field h = Bo/|x - ro zhat| + Bi/|x - ri zhat| (Bo = 10, ro = 15, Bi = 1, ri = 1/2) injected through the Cauchy data at a_L (equivalently at r = 1); h gives every l a nonzero residual, h_l(r) = sqrt(4 pi/(2l+1)) [Bo r^l/ro^(l+1) + Bi ri^l/r^(l+1)].  Exact 3D values: Phi2R(x0) = c0^2/6 + 3 c0 c1 r0/5 + c1^2 r0^2/2 + h(x0) = 1.8138888888888888889, d_r Phi2R(x0) = 0.16061728395061728395.
- Contour y(a) = a + i (1 - 4 (a - r0)^2)^6 on [4.5, 5.5] (apex r0 + i, min |z - r0| = 0.298).  SS + SR are solved together on the contour by variation of parameters, psi~ = A z^l + B z^(-l-1) + [z^l I1 - z^(-l-1) I2]/(2l+1), I1 = Int s^(1-l) s~ ds, I2 = Int s^(l+2) s~ ds along the contour (60 Gauss panels x 20 nodes in a; 50-digit evaluation of the closed-form sources; cross-checked against the NDSolve contour solve to 1e-15).  Sources: s~_SS = (1/2) D_z^2 [q^2/rho^2]_l with 2 log(z + r0) - 2 log(z - r0); s~_SR = (1/2) D_z^2 [2 q sum_l' K~_(j l') phiR_l'] with K~ by Theta -> H (Stage 9).  RR is solved on the real line (polynomial, l <= 2).
- Continued puncture phi~S_l: SS part = out-branch formula of (1/2) D_z^2 [q^2 ln rho]_l with ln|r - r0| -> Log[z - r0] (in = out as real functions); SR part = (1/2) D_z^2 applied to the H-blend F_in + (F_out - F_in) H of the undifferentiated modes F_j = [chi_SR]_j, which are C^2 at r0.  Identities checked: s~_SS = L_l phi~S_SS exactly (1e-40); Re phi~S(r + i0) = [Phi2S]_l(r) on both sides; s~R := s~_SS + s~_SR - L_l phi~S is imaginary on the axis and bounded (analytic) at r0.
- Truncated punctures (b): Phi2P = Phi2S + eps q^2 rho^k, eps = r0^-(k+2), k = 1, 3, 5; the modes [rho^k]_l are C^(k+1) at r0 (checked) and are continued by the same H-blend.  Exact residual modes at r0: h_l(r0) - eps q^2 [rho^k]_l(r0), and the derivative.
- Taylor extension: with z = z_e + w and psi~R = sum p_n w^n, the ODE z^2 psi'' + 2 z psi' - l(l+1) psi = z^2 s~R gives p_(n+2) = [sigma_n - 2 z_e (n+1)^2 p_(n+1) - (n(n+1) - l(l+1)) p_n]/(z_e^2 (n+2)(n+1)), where sigma_n are the Taylor coefficients of z^2 s~R.  p_0, p_1 are psi~R(z_e), psi~R'(z_e) from the contour solve; the coefficients of s~R come from the Cauchy integral on |z - z_e| = 0.7 |z_e - r0| by the trapezoidal rule with 256 points (every quantity O(1); no fit to contour samples).  Sum to r0: S_N = sum_(n<=N) p_n (r0 - z_e)^n, S'_N = sum n p_n (r0 - z_e)^(n-1).  Expansion points: the apex (a = 5) and a = 5.2 (|z_e - r0| = 0.404).
- Cross-checks: (i) the residual ODE solved by NDSolve around the circle |z - z_e| = rc and the Taylor coefficients taken by Cauchy integrals; (ii) the monodromy of psi~R around r0 along |z - r0| = R from a_L (above) through the apex, a_R and back to a_L (below), with Log[z - r0] = ln R + i(pi - t) continued explicitly.
- Baseline (c): contours y_eta(a) = a + i eta (1 - 4(a - r0)^2)^6, eta = 1 .. 1/128, graded Gauss panels; psi~ - phi~P read off at the apex r0 + i eta; polynomial (Richardson) extrapolation in eta.
- Controls (d): d1 = Theta -> H blend of the finished Dzz SR modes P_in + (P_out - P_in) H (the task's control); d2 = one-sided continuation P_out (drops the kink, hence the delta of L P).

## (a) Exact puncture: residual modes at r0 from the apex Taylor series

| l | Re psi~R(r0), apex Taylor | h_l(r0) exact | error | Im psi~R(r0) | derivative error | ratio (value) | N for 1e-20 | Im psi~(a_R) | error, a_e = 5.2 | ratio, a_e = 5.2 |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3.072253341569561 | 3.072253341569561 | <1e-36 | 0.1442495555 | <1e-34 | 2.16e-1 | 29 | 0.18904211 | <1e-37 | 8.64e-2 |
| 1 | 0.4957449385162989 | 0.4957449385162989 | <1e-36 | 0.249103243 | <1e-34 | 2.23e-1 | 30 | 0.32808052 | <1e-37 | 9.19e-2 |
| 2 | 0.1206025817671518 | 0.1206025817671518 | <1e-36 | 0.3196287231 | <1e-34 | 2.27e-1 | 32 | 0.42506622 | <1e-37 | 9.58e-2 |
| 3 | 0.03335066542393942 | 0.03335066542393942 | <1e-35 | 0.3745858713 | <1e-33 | 2.36e-1 | 32 | 0.50513716 | <1e-37 | 1.01e-1 |
| 4 | 0.009749031076890076 | 0.009749031076890076 | <1e-35 | 0.4190283718 | <1e-33 | 2.42e-1 | 32 | 0.57496712 | <1e-36 | 1.04e-1 |
| 5 | 0.002934455595514837 | 0.002934455595514837 | <1e-35 | 0.4548758803 | <1e-33 | 2.43e-1 | 34 | 0.63654056 | <1e-36 | 1.04e-1 |
| 6 | 0.0008993100324394711 | 0.0008993100324394711 | <1e-35 | 0.4827702542 | <1e-33 | 2.46e-1 | 35 | 0.68940795 | <1e-36 | 1.03e-1 |
| 7 | 0.0002790279790036754 | 0.0002790279790036754 | <1e-35 | 0.5026274122 | <1e-33 | 2.45e-1 | 35 | 0.73112998 | <1e-36 | 9.71e-2 |
| 8 | 0.00008736303692902525 | 0.00008736303692902525 | <1e-35 | 0.513828028 | <1e-33 | 2.45e-1 | 38 | 0.75711868 | <1e-36 | 1.02e-1 |
| 9 | 0.00002754534045010747 | 0.00002754534045010747 | <1e-35 | 0.5152700619 | <1e-33 | 2.46e-1 | 39 | 0.76012252 | <1e-36 | 1.07e-1 |
| 10 | 8.733578893381153e-6 | 8.733578893381153e-6 | <1e-35 | 0.5053524508 | <1e-32 | 2.46e-1 | 39 | 0.72938155 | <1e-36 | 1.13e-1 |
| 11 | 2.781737747564759e-6 | 2.781737747564759e-6 | <1e-34 | 0.4819141618 | <1e-32 | 2.48e-1 | 41 | 0.64939424 | <1e-36 | 1.15e-1 |
| 12 | 8.893827077659471e-7 | 8.893827077659471e-7 | <1e-34 | 0.4421362915 | <1e-32 | 2.53e-1 | 42 | 0.49818139 | <1e-36 | 1.16e-1 |
| 13 | 2.852696046860957e-7 | 2.852696046860957e-7 | <1e-34 | 0.3824074773 | <1e-32 | 2.58e-1 | 44 | 0.24487551 | <1e-36 | 1.19e-1 |
| 14 | 9.175233262571581e-8 | 9.175233262571581e-8 | <1e-34 | 0.2981487084 | <1e-32 | 2.61e-1 | 45 | -0.15360563 | <1e-35 | 1.22e-1 |
| 15 | 2.958107871591556e-8 | 2.958107871591556e-8 | <1e-34 | 0.1835905727 | <1e-32 | 2.65e-1 | 46 | -0.75713507 | <1e-35 | 1.24e-1 |
| 16 | 9.556890897024412e-9 | 9.556890897024412e-9 | <1e-34 | 0.03149319467 | <1e-32 | 2.67e-1 | 47 | -1.6486604 | <1e-35 | 1.28e-1 |
| 17 | 3.093273497949863e-9 | 3.093273497949863e-9 | <1e-34 | -0.1672038363 | <1e-32 | 2.72e-1 | 47 | -2.9428224 | <1e-35 | 1.31e-1 |
| 18 | 1.002836715976973e-9 | 1.002836715976973e-9 | <1e-34 | -0.4238176455 | <1e-32 | 2.75e-1 | 49 | -4.7977666 | <1e-35 | 1.33e-1 |
| 19 | 3.255948521264635e-10 | 3.255948521264635e-10 | <1e-34 | -0.752464525 | <1e-32 | 2.75e-1 | 50 | -7.4313196 | <1e-35 | 1.35e-1 |
| 20 | 1.058514106744954e-10 | 1.058514106744954e-10 | <1e-34 | -1.170733767 | <1e-32 | 2.78e-1 | 50 | -11.143131 | <1e-35 | 1.36e-1 |

Predicted geometric ratio |r0 - z_e|/|z_e| (next singularity at z = 0): {0.196116, 0.0775617} for a_e = {5., 5.2}.

## Convergence rates (all cases)

| case | a_e | error fits: geometric ratio, value (median over l) | ratio, derivative | algebraic slope, value | slope, derivative | coefficient fits: ratio | exponent | max_l error at N = 120, value | derivative |
|---|---|---|---|---|---|---|---|---|---|
| exact | 5. | 2.46e-1 | 2.52e-1 | -7.72e1 | -7.70e1 | 1.99e-1 | -8.15e1 | <1e-34 | <1e-32 |
| exact | 5.2 | 1.13e-1 | 1.14e-1 | -9.67e1 | -9.79e1 | 7.98e-2 | -8.24e1 | <1e-35 | <1e-33 |
| k=1 | 5. | 9.33e-1 | 9.30e-1 | -3.04e0 | -3.08e0 | 1.00e0 | -4.07e0 | 8.69e-10 | 7.10e-10 |
| k=1 | 5.2 | 9.17e-1 | 9.55e-1 | -3.22e0 | -2.04e0 | 1.00e0 | -4.07e0 | 1.16e-12 | 2.24e-8 |
| k=3 | 5. | 8.73e-1 | 8.69e-1 | -5.14e0 | -5.21e0 | 1.00e0 | -6.17e0 | 1.78e-14 | 2.22e-14 |
| k=3 | 5.2 | 8.71e-1 | 8.99e-1 | -5.15e0 | -4.16e0 | 1.00e0 | -6.18e0 | 1.64e-16 | 6.40e-14 |
| k=5 | 5. | 8.07e-1 | 8.08e-1 | -7.31e0 | -7.42e0 | 1.00e0 | -8.33e0 | 1.12e-18 | 1.89e-18 |
| k=5 | 5.2 | 8.10e-1 | 8.23e-1 | -7.33e0 | -6.55e0 | 1.00e0 | -8.34e0 | 1.77e-21 | 2.88e-20 |
| d1 | 5. | 2.49e-1 | 2.51e-1 | -7.74e1 | -7.79e1 | 1.99e-1 | -8.20e1 | <1e-34 | <1e-32 |
| d1 | 5.2 | 1.14e-1 | 1.17e-1 | -9.71e1 | -9.75e1 | 8.03e-2 | -8.18e1 | <1e-35 | <1e-33 |
| d2 | 5. | 9.89e-1 | Indeterminate | -9.99e-1 | 1.27e-3 | 1.00e0 | -2.01e0 | 3.55e-3 | 6.80e-1 |
| d2 | 5.2 | 9.90e-1 | Indeterminate | -10.00e-1 | 3.79e-4 | 1.00e0 | -2.01e0 | 1.25e-3 | 4.56e-1 |

## Mode sums at x0 (Taylor-extracted modes + RR on the real line)

| case | Phi2R(x0), Lmax = 20 | error | l-truncation alone | d_r Phi2R(x0) | error | l-truncation alone |
|---|---|---|---|---|---|---|
| exact | 1.81388888879329 | 9.56e-11 | 9.56e-11 | 0.1606172835395413 | 4.11e-10 | 4.11e-10 |
| k=1 | 1.811983029587102 | 1.91e-3 | 1.91e-3 | 0.1604267131407727 | 1.91e-4 | 1.91e-4 |
| k=3 | 1.813901920194075 | 1.30e-5 | 1.30e-5 | 0.1606211929592409 | 3.91e-6 | 3.91e-6 |
| k=5 | 1.813888439176253 | 4.50e-7 | 4.50e-7 | 0.1606170587310231 | 2.25e-7 | 2.25e-7 |
| d1 | 1.81388888879329 | 9.56e-11 | 9.56e-11 | 0.1606172835395413 | 4.11e-10 | 4.11e-10 |
| d2 | 1.74488016940714 | 6.90e-2 | 9.56e-11 | 13.38131915200352 | 1.32e1 | 4.11e-10 |

## Cross-check: Cauchy circle about the apex

| case | l | max_n |dp_n| rc^n (n <= 40) | |u(2 Pi) - u(0)| on the circle | partial-sum difference at r0 (N = 40) |
|---|---|---|---|---|
| exact | 0 | 6.12e-20 | 7.43e-20 | 3.78e-16 |
| exact | 2 | 2.26e-20 | 1.40e-20 | 5.81e-17 |
| exact | 5 | 4.32e-20 | 4.19e-20 | 2.12e-16 |
| exact | 10 | 2.46e-19 | 9.75e-20 | 6.01e-16 |
| k=3 | 0 | 3.57e-20 | 4.41e-20 | 2.43e-16 |
| k=3 | 2 | 2.18e-20 | 1.63e-20 | 7.28e-17 |
| k=3 | 5 | 3.98e-20 | 4.67e-20 | 2.70e-16 |
| k=3 | 10 | 2.46e-19 | 9.76e-20 | 6.02e-16 |

## Monodromy of psi~R around r0

| case | l | monodromy of psi~R around r0 | predicted | difference | psi~R(a_R): semicircle - contour |
|---|---|---|---|---|---|
| exact | 0 | 1.34e-21 + 3.76e-20 i | 0 + 0 i | 3.76e-20 | 4.69e-21 |
| exact | 2 | -1.17e-20 - 1.73e-21 i | 0 + 0 i | 1.19e-20 | 4.34e-21 |
| exact | 5 | 2.24e-20 + 9.60e-21 i | 0 + 0 i | 2.43e-20 | 5.22e-21 |
| k=3 | 0 | -6.30e-7 + <1e-20 i | -6.30e-7 + <1e-14 i | 3.36e-20 | 4.58e-21 |
| k=3 | 2 | -1.42e-6 + <1e-21 i | -1.42e-6 + <1e-14 i | 1.18e-20 | 4.32e-21 |
| k=3 | 5 | -2.14e-6 + <1e-20 i | -2.14e-6 + <1e-14 i | 2.41e-20 | 5.23e-21 |
| d1 | 0 | -7.11e-21 - 1.16e-20 i | <1e-51 + <1e-51 i | 1.36e-20 | 5.29e-21 |
| d1 | 2 | -7.26e-21 + 7.84e-21 i | <1e-51 + <1e-51 i | 1.07e-20 | 3.03e-21 |
| d1 | 5 | 2.21e-20 + 3.47e-21 i | <1e-51 + <1e-51 i | 2.24e-20 | 9.94e-21 |
| d2 | 0 | 2.06e-1 + <1e-20 i | 2.06e-1 + <1e-9 i | 3.04e-20 | 8.20e-22 |
| d2 | 2 | 4.72e-1 + <1e-20 i | 4.72e-1 + <1e-8 i | 1.82e-20 | 1.28e-20 |
| d2 | 5 | 7.63e-1 + <1e-21 i | 7.63e-1 + <1e-8 i | 1.28e-20 | 2.56e-20 |

## (c) eta-baseline (l = 2)

| eta | nodes | exact: bias of the value | bias of the derivative | Taylor-predicted bias (value) | k=1 bias | k=3 bias | k=5 bias | abs psi~(apex) | abs phi~S(apex) | max abs source | Im psi~R(apex) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1. | 840 | 2.05e-1 | 5.95e-2 | 2.05e-1 | 2.04e-1 | 2.05e-1 | 1.67e-1 | 3.62e-1 | 2.86e-1 | 1.42e0 | 0.33586541 |
| 0.5 | 880 | 1.04e-1 | 2.83e-2 | 1.04e-1 | 1.04e-1 | 1.04e-1 | 9.48e-2 | 4.03e-1 | 3.07e-1 | 1.78e0 | 0.33513815 |
| 0.25 | 920 | 5.22e-2 | 1.40e-2 | 5.22e-2 | 5.21e-2 | 5.22e-2 | 4.98e-2 | 4.64e-1 | 3.50e-1 | 2.59e0 | 0.32915073 |
| 0.125 | 960 | 2.60e-2 | 6.98e-3 | 2.60e-2 | 2.60e-2 | 2.60e-2 | 2.55e-2 | 5.30e-1 | 4.06e-1 | 6.15e0 | 0.32482658 |
| 0.0625 | 1000 | 1.30e-2 | 3.49e-3 | 1.30e-2 | 1.30e-2 | 1.30e-2 | 1.29e-2 | 5.95e-1 | 4.66e-1 | 2.25e1 | 0.32233663 |
| 0.03125 | 1040 | 6.50e-3 | 1.74e-3 | 6.50e-3 | 6.50e-3 | 6.50e-3 | 6.46e-3 | 6.57e-1 | 5.25e-1 | 8.58e1 | 0.32100992 |
| 0.015625 | 1080 | 3.25e-3 | 8.72e-4 | 3.25e-3 | 3.25e-3 | 3.25e-3 | 3.24e-3 | 7.16e-1 | 5.84e-1 | 3.34e2 | 0.32032613 |
| 0.0078125 | 1120 | 1.62e-3 | 4.36e-4 | 1.62e-3 | 1.62e-3 | 1.62e-3 | 1.62e-3 | 7.73e-1 | 6.41e-1 | 1.32e3 | 0.31997913 |

| l | numerical error of psi~(apex) at eta = 0.0078125 (20 -> 30 nodes per panel): value | derivative |
|---|---|---|
| 0 | 4.30e-35 | 5.64e-33 |
| 2 | 9.61e-35 | 1.26e-32 |
| 5 | 1.43e-34 | 1.87e-32 |
| 10 | 1.97e-34 | 2.58e-32 |

Richardson extrapolation eta -> 0:

| case | l | extrapolated error (all eta) | (4 smallest) | derivative (all) | (4 smallest) | error at the smallest eta | derivative |
|---|---|---|---|---|---|---|---|
| exact | 0 | 3.08e-15 | 2.67e-10 | 3.50e-15 | 2.65e-10 | 7.67e-4 | 8.98e-5 |
| exact | 2 | 2.34e-14 | 3.68e-12 | 6.31e-14 | 8.85e-13 | 1.62e-3 | 4.36e-4 |
| exact | 5 | 2.19e-13 | 1.12e-10 | 5.18e-13 | 2.74e-10 | 1.83e-3 | 4.25e-3 |
| exact | 10 | 5.22e-11 | 4.52e-9 | 1.97e-10 | 1.99e-8 | 1.27e-3 | 1.89e-2 |
| k=1 | 0 | 1.26e-10 | 6.22e-10 | 1.01e-10 | 5.48e-10 | 7.67e-4 | 8.97e-5 |
| k=1 | 2 | 2.82e-10 | 7.95e-10 | 2.25e-10 | 6.19e-10 | 1.62e-3 | 4.36e-4 |
| k=1 | 5 | 4.18e-10 | 1.26e-9 | 3.34e-10 | 1.15e-9 | 1.83e-3 | 4.25e-3 |
| k=1 | 10 | 6.28e-10 | 6.04e-9 | 6.55e-10 | 2.10e-8 | 1.27e-3 | 1.89e-2 |
| k=3 | 0 | 2.04e-14 | 2.56e-10 | 2.47e-14 | 2.60e-10 | 7.66e-4 | 8.99e-5 |
| k=3 | 2 | 7.59e-14 | 6.25e-12 | 1.26e-13 | 5.58e-12 | 1.62e-3 | 4.36e-4 |
| k=3 | 5 | 1.41e-13 | 1.06e-10 | 4.25e-13 | 2.72e-10 | 1.83e-3 | 4.25e-3 |
| k=3 | 10 | 5.20e-11 | 4.52e-9 | 1.97e-10 | 1.99e-8 | 1.27e-3 | 1.89e-2 |
| k=5 | 0 | 3.20e-15 | 1.14e-11 | 3.66e-15 | 3.12e-10 | 7.58e-4 | 9.29e-5 |
| k=5 | 2 | 2.31e-14 | 1.59e-12 | 6.29e-14 | 6.12e-12 | 1.62e-3 | 4.36e-4 |
| k=5 | 5 | 2.19e-13 | 1.12e-10 | 5.18e-13 | 2.74e-10 | 1.83e-3 | 4.25e-3 |
| k=5 | 10 | 5.22e-11 | 4.52e-9 | 1.97e-10 | 1.99e-8 | 1.27e-3 | 1.89e-2 |

## Assessment of the predictions

1. Exact puncture: analytic at r0, geometric convergence with ratio set by z = 0.  HELD: max_l error <1e-34 (value), <1e-32 (derivative) at N = 120 (floor of the 50-digit arithmetic; 1e-20 is reached at N = 29-50).  Ratio of the Taylor coefficients |p_n w0^n| ~ n^c rho^n (median over l): rho = 1.99e-1 vs predicted 0.196116 (apex), 7.98e-2 vs 0.0775617 (a = 5.2); a two-parameter fit of the error curve gives 2.46e-1 and 1.13e-1, biased upwards by the n^l prefactor of the order-(l+1) pole at z = 0.  The monodromy around r0 is 3.76e-20.  Im psi~R(r0) is O(1) (table above): the shadow is carried through to r0 but only the real part is physical.  The second point a = 5.2 converges faster (ratio 0.08): in the exact case every contour point works, the closer to r0 the better.
2. Truncated puncture: algebraic convergence.  Measured median slopes of the real-part errors, apex, value / derivative: k = 1: -3.04e0 / -3.08e0; k = 3: -5.14e0 / -5.21e0; k = 5: -7.31e0 / -7.42e0; coefficient-fit exponents of |p_n w0^n|: k = 1: -4.07e0; k = 3: -6.17e0; k = 5: -8.33e0.  The residual modes are C^(k+1) (two orders smoother than the 3D residual, which is C^(k-1)); the continuation has a (z - r0)^(k+2) log(z - r0) branch point, so the rates expected from the mode smoothness are N^-(k+2) (value), N^-(k+1) (derivative) and n^-(k+3) for the coefficients, against the N^-k of the prediction.  At a = 5.2 the slopes are k = 1: -3.22e0 / -2.04e0; k = 3: -5.15e0 / -4.16e0; k = 5: -7.33e0 / -6.55e0.  The algebraic character held, with exponent k+2 for the value (two better than N^-k, because the mode functions are two orders smoother than the 3D field) and k+1 for the derivative off the vertical.  At the apex, which lies directly above r0, the derivative is NOT worse: the leading term (i/pi) j (z - r0)^(k+2) log(z - r0) has Taylor terms at r0 of constant phase i^(k+2), imaginary for odd k, so after the factor i/pi it feeds only Re of the value and Im of the derivative; Re of the derivative inherits the next term and also falls like N^-(k+2).  'Worse for the derivative' therefore held only for the off-vertical expansion point.  r0 is a branch point on the edge of the disc for every expansion point, so the expansion point changes the prefactor (|z_e - r0|^(k+2)) and the phase pattern, not the exponent.
3. Consistency (delta content).  Control d1 (blend of the finished, kinked modes) CONVERGED to the same real values (max_l error <1e-34): the prediction that it fails did NOT hold.  Reason: for C^0 (kinked) finished modes, L[P_in + (P_out - P_in) H] = s_in + (s_out - s_in) H + [P'] delta_ext + analytic, so the kink reproduces the delta of L P automatically; d1 differs from the blend-then-ladder puncture by a function that is imaginary on the axis and vanishes at r0 (checked), i.e. by a harmless shadow.  Control d2 (one-sided continuation, kink dropped) FAILS as a delta-dropping control should: its real part is wrong by J_P(r) Theta(r0 - r) at a_L, the residual acquires a term (i/pi) c_l (z - r0) log(z - r0) (c_l the dropped delta coefficient), the value converges only like 1/N and the real part of the derivative tends to the mean of the two one-sided derivatives, h_l' + c_l/2, so its error is |c_l|/2 = 0.106 (l = 0), 0.238 (l = 2) (min_l derivative error 1.06e-1 at N = 120); its real part is also wrong at a_L.
4. Baseline (c): the bias of Re psi~R(r0 + i eta) is O(eta) with the coefficient -Im psi~R'(r0) (the shadow derivative), reproduced by the Taylor series about r0 to 2.54e-19 (relative).  Richardson extrapolation of the exact case reaches 5.22e-11 (all etas); for the truncated punctures it saturates at the eta^(k+2) log eta level.  Digits are lost mainly through the O(eta) bias and the node count needed near the apex (table), not through cancellation: |phi~S(apex)| grows only logarithmically.

## Surprises and failed predictions, stated plainly

- Control d1 (blend of the finished Dzz SR modes) does not fail: it converges exactly like the blend-then-ladder puncture.  The delta content of L[puncture] lives in the kink of the finished modes, which both continuations keep; they differ by an analytic function that is imaginary on the axis.  Prediction 3, as a requirement on the puncture, did not hold in this toy.  What does fail is dropping the kink (d2).
- The algebraic exponent is k+2, not k, for the value: the mode functions [rho^k]_l are C^(k+1), two orders smoother than the 3D field.
- At the apex the real part of the derivative converges as fast as the value (N^-(k+2)), because the leading log term feeds only Im of the derivative there; off the vertical the derivative is one order worse (N^-(k+1)).
- With Lmax = 20 the mode sums are limited entirely by the l-truncation (exact: (r0/ro)^l from h; truncated: L^-k of the rho^k modes); the Taylor extraction contributes < 1e-30.
- The two-parameter fit of the exact-case error curve overestimates the geometric ratio (0.25 vs 0.196) because of the n^l prefactor; the three-parameter coefficient fit recovers 0.199.
- The eta-baseline loses no digits to cancellation (|psi~| and |phi~S| at the apex stay O(1)); its limitation is the O(eta) bias and the 1/eta^2 growth of the source, which costs quadrature nodes.  Richardson extrapolation over eta = 1 .. 1/128 recovers the exact-case value to 3e-15 - 5e-11 and is capped at eta^(k+2) log eta for truncated punctures (1e-10 for k = 1).

## Checks

363 passed, 0 failed (of 363).  Failed: {}

## Files

- results/stage15_checks.csv, stage15_residual_contour.csv, stage15_taylor_convergence.csv, stage15_modes_r0.csv, stage15_rates.csv, stage15_modesum.csv, stage15_circle_check.csv, stage15_monodromy.csv, stage15_shadow.csv, stage15_eta_baseline.csv, stage15_eta_richardson.csv.
- figures/residual/fig_residual_contour.pdf, fig_taylor_convergence.pdf, fig_modesum.pdf, fig_eta_baseline.pdf (each with a CSV).
