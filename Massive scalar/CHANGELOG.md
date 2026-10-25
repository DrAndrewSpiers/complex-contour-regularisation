# Massive scalar test: changelog

One line per completed stage (date, what, headline number).

- 2026-10-25  Package `Mathematica/ContourMassive.wl`: massive radial operator, Bessel homogeneous solutions u_l, v_l normalised to the massless powers, closed-form modes with exponential integrals (LOG placeholder for the branch), rho-substitution angular quadrature with the residue-corrected Theta->H blend, Cauchy continuation and solvers (copied from the massless package where unchanged).  Angular quadrature of the massive SS source agrees with the closed-form continuation to 1e-14 (80 nodes per segment).
- 2026-10-25  Stage 1 (massive_stage1_first_order.wls, 406 checks pass, 7.7 min): consistency conditions with Bessel functions; (i) holds to 1e-76, (ii) fails by a purely imaginary (i/pi)E0 with E0(r0) = -1/r0 for every l, m; corrected mode satisfies both; exterior regular field exact (<= 1e-17) for Theta->H, corrected and truncated punctures; residual modes C^(n+3) (3D residual C^(n+1)); DSolve offsets (alpha, beta) tabulated.
