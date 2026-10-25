# Stage 12 (D): complex-r angular quadrature

| l | SS in: pointwise rel. diff | SS out: pointwise | SS in: exterior | SS out: exterior | SR in only: exterior | SR, H blend after D_z^2: exterior | SR, H blend of products then D_z^2: pointwise |
|---|---|---|---|---|---|---|---|
| 0 | 2.11e-2 | 1.82e-14 | 1.98e-14 | 2.95e-14 | 8.37e0 | 6.50e0 | 2.40e-14 |
| 1 | 2.86e-2 | 1.71e-14 | 1.23e-14 | 1.90e-14 | 1.91e0 | 1.49e0 | 2.43e-14 |
| 2 | 5.18e-2 | 1.73e-14 | 1.19e-14 | 1.61e-14 | 3.04e0 | 2.36e0 | 2.40e-14 |
| 3 | 8.83e-2 | 1.58e-14 | 6.81e-15 | 1.08e-14 | 7.18e0 | 4.98e0 | 2.53e-14 |
| 4 | 1.40e-1 | 1.50e-14 | 4.84e-15 | 8.09e-15 | 5.64e0 | 3.69e0 | 2.78e-14 |
| 5 | 2.10e-1 | 1.44e-14 | 5.18e-15 | 6.86e-15 | 4.69e0 | 2.85e0 | 3.08e-14 |
| 6 | 3.03e-1 | 1.55e-14 | 4.81e-15 | 6.31e-15 | 3.96e0 | 2.20e0 | 3.56e-14 |

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
