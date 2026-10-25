# Stage 12 (D): complex-r angular quadrature

| l | SS in: pointwise rel. diff | SS out: pointwise | SS in: exterior | SS out: exterior | SR in only: exterior | SR, H blend after D_z^2: exterior | SR, H blend of products then D_z^2: pointwise |
|---|---|---|---|---|---|---|---|
| 0 | 2.11e-2 | 1.82e-14 | 6.97e-15 | 9.84e-15 | 4.14e0 | 3.25e0 | 2.40e-14 |
| 1 | 2.86e-2 | 1.71e-14 | 3.26e-14 | 3.89e-14 | 4.75e0 | 3.72e0 | 2.43e-14 |
| 2 | 5.18e-2 | 1.73e-14 | 9.31e-14 | 7.27e-14 | 2.40e1 | 1.86e1 | 2.40e-14 |
| 3 | 8.83e-2 | 1.58e-14 | 1.32e-13 | 1.10e-13 | 1.47e2 | 1.13e2 | 2.53e-14 |
| 4 | 1.40e-1 | 1.50e-14 | 3.22e-13 | 1.97e-13 | 2.90e2 | 2.19e2 | 2.78e-14 |
| 5 | 2.10e-1 | 1.44e-14 | 7.34e-13 | 4.22e-13 | 5.94e2 | 4.40e2 | 3.08e-14 |
| 6 | 3.03e-1 | 1.55e-14 | 1.69e-12 | 9.12e-13 | 1.22e3 | 8.83e2 | 3.56e-14 |

## cost (l = 2)

| point | adaptive evaluations for 1e-12 | relative error |
|---|---|---|
| real r = r0 - 1/10 | 430 | 7.49e-13 |
| real r = r0 - 1/100 | 1966 | 9.98e-16 |
| real r = r0 - 1/400 | 1122 | 1.00e-15 |
| real r = r0 - 1/1000 | 1122 | 1.47e-15 |
| contour z = 4.6 + 0.00217678 I | 1207 | 6.09e-15 |
| contour z = 5. + 1. I | 1031 | 1.22e-13 |
| contour z = 5.4 + 0.00217678 I | 1229 | 4.41e-15 |

SR Part 1 (Stage 6 puncture, n = 3, k = 2, l = 2): in + (out - in) H from angular quadrature vs closed form, max relative difference 1.08e-13

## continued residual modes of the Stage 6 puncture (n = 3), summed at complex z

| a | abs(z)/r0 | growth per l' (fit 60-80) | smallest term at l' | abs(term) at l' = 80 | abs(partial sum at theta=0), l' = 20, 40, 60, 80 | real r = Re z: same |
|---|---|---|---|---|---|---|
| 4.6 | 0.92 | 1.058 | 27 | 2.36e-5 | 1.46e0, 1.46e0, 1.46e0, 1.46e0 | 1.46e0, 1.46e0, 1.46e0, 1.46e0 |
| 5. | 1.0198 | 0.982 | 80 | 2.57e-6 | 1.50e0, 1.50e0, 1.50e0, 1.50e0 | 1.50e0, 1.50e0, 1.50e0, 1.50e0 |
| 5.4 | 1.08 | 1.052 | 29 | 1.46e-5 | 1.54e0, 1.54e0, 1.54e0, 1.54e0 | 1.54e0, 1.54e0, 1.54e0, 1.54e0 |
| 5.49 | 1.098 | 1.067 | 25 | 7.19e-5 | 1.55e0, 1.55e0, 1.55e0, 1.55e0 | 1.55e0, 1.55e0, 1.55e0, 1.55e0 |
