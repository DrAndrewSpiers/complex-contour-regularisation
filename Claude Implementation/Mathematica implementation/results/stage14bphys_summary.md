# Stage 14b: radial samples for the exterior SS field

| method | Lmax | target | radial samples | error reached | configuration |
|---|---|---|---|---|---|
| contour | 10 | 1.00e-8 | 64 | 4.44e-11 | {4, 16} |
| contour | 10 | 1.00e-12 | 80 | 1.45e-13 | {4, 20} |
| contour | 20 | 1.00e-8 | 64 | 4.44e-11 | {4, 16} |
| contour | 20 | 1.00e-12 | 80 | 1.45e-13 | {4, 20} |
| contour | 30 | 1.00e-8 | 64 | 4.44e-11 | {4, 16} |
| contour | 30 | 1.00e-12 | 80 | 1.45e-13 | {4, 20} |
| R1a | 10 | 1.00e-8 | 800 | 1.82e-9 | 400 |
| R1a | 10 | 1.00e-12 | not reached | 3.57e-12 | None |
| R1a | 20 | 1.00e-8 | 800 | 2.59e-9 | 400 |
| R1a | 20 | 1.00e-12 | not reached | 5.06e-12 | None |
| R1a | 30 | 1.00e-8 | 800 | 5.00e-9 | 400 |
| R1a | 30 | 1.00e-12 | not reached | 9.78e-12 | None |
| R1b | 10 | 1.00e-8 | 100 | 4.12e-10 | 50 |
| R1b | 10 | 1.00e-12 | 400 | 4.12e-13 | 200 |
| R1b | 20 | 1.00e-8 | 100 | 1.80e-9 | 50 |
| R1b | 20 | 1.00e-12 | 800 | 5.68e-14 | 400 |
| R1b | 30 | 1.00e-8 | 100 | 7.65e-9 | 50 |
| R1b | 30 | 1.00e-12 | 800 | 2.38e-13 | 400 |
| R2 | 10 | 1.00e-8 | 36 | 9.05e-9 | {2, 6} |
| R2 | 10 | 1.00e-12 | 120 | 2.28e-14 | {4, 12} |
| R2 | 20 | 1.00e-8 | 48 | 2.68e-10 | {2, 8} |
| R2 | 20 | 1.00e-12 | 120 | 3.50e-14 | {4, 12} |
| R2 | 30 | 1.00e-8 | 48 | 3.88e-10 | {2, 8} |
| R2 | 30 | 1.00e-12 | 120 | 3.50e-14 | {4, 12} |

## 1e-8 relative noise on the SS samples

| method | configuration chosen for (Lmax = 30) | radial samples | error with 1e-8 relative noise on SS (max of 3 seeds) |
|---|---|---|---|
| contour | 1.00e-8 | 64 | 8.28e-8 |
| contour | 1.00e-12 | 80 | 5.32e-8 |
| R1a | 1.00e-8 | 800 | 9.31e-5 |
| R1a | 1.00e-12 | - | - |
| R1b | 1.00e-8 | 100 | 1.02e-5 |
| R1b | 1.00e-12 | 800 | 9.31e-5 |
| R2 | 1.00e-8 | 48 | 3.72e-4 |
| R2 | 1.00e-12 | 120 | 1.50e-2 |
