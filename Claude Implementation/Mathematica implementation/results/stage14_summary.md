# Stage 14: SS angular-integral cost, real axis vs contour

| Lmax | target | method | radial points | S evaluations | max error | wall time (s) |
|---|---|---|---|---|---|---|
| 10 | 1.00e-8 | (A) real axis, adaptive per (l, r), PrecisionGoal 8 | 400 | 1161338 | 2.49e-11 | 52.9677 |
| 10 | 1.00e-12 | (A) real axis, adaptive per (l, r), PrecisionGoal 12 | 400 | 2435622 | 2.41e-11 | 215.779 |
| 10 | 1.00e-8 | (D) real axis, shared nodes in ln rho | 400 | 9736 | 9.74e-9 | 0.049033 |
| 10 | 1.00e-12 | (D) real axis, shared nodes in ln rho | 400 | 12816 | 6.44e-13 | 0.059146 |
| 10 | 1.00e-8 | (B) contour, per l (2 nq per segment pair) | 64 | 86016 | 8.65e-9 | 1.16353 |
| 10 | 1.00e-8 | (C) contour, shared nodes | 64 | 8192 | 9.11e-10 | 0.027154 |
| 10 | 1.00e-12 | (B) contour, per l (2 nq per segment pair) | 96 | 165888 | 8.70e-13 | 2.2646 |
| 10 | 1.00e-12 | (C) contour, shared nodes | 96 | 15360 | 4.61e-13 | 0.050725 |
| 20 | 1.00e-8 | (A) real axis, adaptive per (l, r), PrecisionGoal 8 | 400 | 2559287 | 1.16e-10 | 101.12 |
| 20 | 1.00e-12 | (A) real axis, adaptive per (l, r), PrecisionGoal 12 | 400 | 7406023 | 6.84e-11 | 411.942 |
| 20 | 1.00e-8 | (D) real axis, shared nodes in ln rho | 400 | 16512 | 9.71e-9 | 0.080905 |
| 20 | 1.00e-12 | (D) real axis, shared nodes in ln rho | 400 | 19032 | 9.95e-13 | 0.088166 |
| 20 | 1.00e-8 | (B) contour, per l (2 nq per segment pair) | 64 | 157696 | 8.65e-9 | 3.02988 |
| 20 | 1.00e-8 | (C) contour, shared nodes | 64 | 8192 | 9.11e-10 | 0.02938 |
| 20 | 1.00e-12 | (B) contour, per l (2 nq per segment pair) | 96 | 305664 | 9.24e-13 | 5.95079 |
| 20 | 1.00e-12 | (C) contour, shared nodes | 96 | 15360 | 4.61e-13 | 0.053969 |
| 30 | 1.00e-8 | (A) real axis, adaptive per (l, r), PrecisionGoal 8 | 400 | 5374378 | 7.34e-7 | 149.273 |
| 30 | 1.00e-12 | (A) real axis, adaptive per (l, r), PrecisionGoal 12 | 400 | 49739016 | 5.36e-10 | 608.105 |
| 30 | 1.00e-8 | (D) real axis, shared nodes in ln rho | 400 | 22856 | 9.94e-9 | 0.109632 |
| 30 | 1.00e-12 | (D) real axis, shared nodes in ln rho | 400 | 25112 | 9.93e-13 | 0.116141 |
| 30 | 1.00e-8 | (B) contour, per l (2 nq per segment pair) | 64 | 230400 | 9.93e-9 | 5.74503 |
| 30 | 1.00e-8 | (C) contour, shared nodes | 64 | 8192 | 9.11e-10 | 0.03184 |
| 30 | 1.00e-12 | (B) contour, per l (2 nq per segment pair) | 96 | Infinity | 1.62e-12 | 14.2461 |
| 30 | 1.00e-12 | (C) contour, shared nodes | 96 | Infinity | 1.62e-12 | 0.025552 |

## per-point cost against |Dr| (Lmax = 30, h = R/200)

| abs(Dr)/h | (A) evals per (l, r), PG 8 | PG 12 | (D) ln-rho nq, 1e-8 | 1e-12 | (D) plain rho nq, 1e-8 | 1e-12 |
|---|---|---|---|---|---|---|
| 1 | 507 | 3618 | 64 | 80 | Infinity | Infinity |
| 2 | 494 | 3703 | 64 | 80 | Infinity | Infinity |
| 5 | 478 | 3607 | 64 | 72 | Infinity | Infinity |
| 10 | 447 | 3442 | 64 | 72 | 160 | Infinity |
| 20 | 430 | 3414 | 64 | 72 | 112 | 160 |
| 50 | 435 | 3842 | 56 | 64 | 72 | 96 |
| 100 | 412 | 4361 | 56 | 64 | 56 | 72 |
| 200 | 498 | 3705 | 64 | 80 | Infinity | Infinity |

## contour: nq per segment and radial nodes

| Lmax | nq per segment for 1e-8 | 1e-12 | contour nodes for 1e-8 (panels x m) | 1e-12 |
|---|---|---|---|---|
| 10 | 64 | 80 | 64 (4x16) | 96 (4x24) |
| 20 | 64 | 80 | 64 (4x16) | 96 (4x24) |
| 30 | 64 | Infinity | 64 (4x16) | 96 (4x24) |
