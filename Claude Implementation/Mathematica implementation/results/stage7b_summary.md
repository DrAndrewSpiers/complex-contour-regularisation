# Stage 7b: l'max optimum and shrinking worldtube, physical boundary conditions

Scaled exterior error (max over l = 0..6 and r in [1,3.5] U [6.5,10]).

| case | n | R | argmin l'max | min total | total at l'max = 40 | behaviour |
|---|---|---|---|---|---|---|
| Product | 3 | 1/2 | 20 | 1.40e-6 | 6.78e-6 | not monotone |
| Product | 3 | 1/4 | 40 | 7.69e-8 | 7.69e-8 | monotone decreasing |
| Product | 3 | r0/l'max | 40 | 2.60e-8 | 2.60e-8 | monotone decreasing |
| Product | 5 | 1/2 | 30 | 8.18e-9 | 3.34e-8 | not monotone |
| Product | 5 | 1/4 | 40 | 2.20e-10 | 2.20e-10 | monotone decreasing |
| Product | 5 | r0/l'max | 40 | 2.88e-10 | 2.88e-10 | monotone decreasing |
| GradDot | 3 | 1/2 | 20 | 1.26e-3 | 2.20e-2 | not monotone |
| GradDot | 3 | 1/4 | 30 | 2.79e-4 | 3.23e-4 | not monotone |
| GradDot | 3 | r0/l'max | 40 | 8.79e-5 | 8.79e-5 | monotone decreasing |
| GradDot | 5 | 1/2 | 20 | 1.63e-5 | 1.28e-4 | not monotone |
| GradDot | 5 | 1/4 | 40 | 3.11e-7 | 3.11e-7 | monotone decreasing |
| GradDot | 5 | r0/l'max | 40 | 7.14e-7 | 7.14e-7 | monotone decreasing |

## growth of (b) per coupled mode

| case | R | growth per l' near l' = 30 (n = 3) | per-mode (b) at l' = 30 |
|---|---|---|---|
| Product | 1/2 | 1.122 | 2.28e-7 |
| Product | 1/4 | 1.061 | 1.91e-9 |
| GradDot | 1/2 | 1.188 | 6.17e-4 |
| GradDot | 1/4 | 1.093 | 8.98e-6 |

(b) at l'max = 12, R = 1/2, falls per unit n by: Product 0.0702, GradDot 0.0792

## Errors against l'max

| case | n | R | quantity | 4 | 8 | 12 | 16 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| Product | 3 | 1/2 | total | 3.75e-2 | 1.24e-4 | 1.21e-5 | 2.89e-6 | 1.40e-6 | 2.13e-6 | 6.78e-6 |
| Product | 3 | 1/2 | (b) | 2.34e-8 | 7.93e-8 | 1.78e-7 | 3.38e-7 | 5.87e-7 | 2.02e-6 | 6.75e-6 |
| Product | 3 | 1/2 | (c) | 3.75e-2 | 1.24e-4 | 1.19e-5 | 2.55e-6 | 8.10e-7 | 1.06e-7 | 2.54e-8 |
| Product | 3 | 1/4 | total | 3.75e-2 | 1.24e-4 | 1.19e-5 | 2.56e-6 | 8.20e-7 | 1.31e-7 | 7.69e-8 |
| Product | 3 | 1/4 | (b) | 5.28e-10 | 1.73e-9 | 3.69e-9 | 6.47e-9 | 1.02e-8 | 2.50e-8 | 5.15e-8 |
| Product | 3 | 1/4 | (c) | 3.75e-2 | 1.24e-4 | 1.19e-5 | 2.55e-6 | 8.10e-7 | 1.06e-7 | 2.54e-8 |
| Product | 3 | r0/l'max | total | 3.75e-2 | 1.24e-4 | 1.20e-5 | 2.58e-6 | 8.20e-7 | 1.08e-7 | 2.60e-8 |
| Product | 3 | r0/l'max | (b) | 1.82e-6 | 2.31e-7 | 6.72e-8 | 2.39e-8 | 1.02e-8 | 2.03e-9 | 6.19e-10 |
| Product | 3 | r0/l'max | (c) | 3.75e-2 | 1.24e-4 | 1.19e-5 | 2.55e-6 | 8.10e-7 | 1.06e-7 | 2.54e-8 |
| Product | 5 | 1/2 | total | 9.40e-2 | 4.37e-5 | 1.66e-6 | 1.86e-7 | 3.42e-8 | 8.18e-9 | 3.34e-8 |
| Product | 5 | 1/2 | (b) | 1.20e-10 | 4.09e-10 | 9.19e-10 | 1.74e-9 | 3.02e-9 | 1.03e-8 | 3.37e-8 |
| Product | 5 | 1/2 | (c) | 9.40e-2 | 4.37e-5 | 1.66e-6 | 1.88e-7 | 3.72e-8 | 2.14e-9 | 2.89e-10 |
| Product | 5 | 1/4 | total | 9.40e-2 | 4.37e-5 | 1.66e-6 | 1.88e-7 | 3.72e-8 | 2.10e-9 | 2.20e-10 |
| Product | 5 | 1/4 | (b) | 7.03e-13 | 2.31e-12 | 4.90e-12 | 8.61e-12 | 1.36e-11 | 3.32e-11 | 6.82e-11 |
| Product | 5 | 1/4 | (c) | 9.40e-2 | 4.37e-5 | 1.66e-6 | 1.88e-7 | 3.72e-8 | 2.14e-9 | 2.89e-10 |
| Product | 5 | r0/l'max | total | 9.40e-2 | 4.37e-5 | 1.66e-6 | 1.88e-7 | 3.72e-8 | 2.13e-9 | 2.88e-10 |
| Product | 5 | r0/l'max | (b) | 6.84e-8 | 1.81e-9 | 2.44e-10 | 4.94e-11 | 1.36e-11 | 1.21e-12 | 2.08e-13 |
| Product | 5 | r0/l'max | (c) | 9.40e-2 | 4.37e-5 | 1.66e-6 | 1.88e-7 | 3.72e-8 | 2.14e-9 | 2.89e-10 |
| GradDot | 3 | 1/2 | total | 4.77e-2 | 7.85e-3 | 2.69e-3 | 1.44e-3 | 1.26e-3 | 4.05e-3 | 2.20e-2 |
| GradDot | 3 | 1/2 | (b) | 1.43e-5 | 5.34e-5 | 1.38e-4 | 3.12e-4 | 6.60e-4 | 3.87e-3 | 2.19e-2 |
| GradDot | 3 | 1/2 | (c) | 4.77e-2 | 7.80e-3 | 2.55e-3 | 1.13e-3 | 5.99e-4 | 1.86e-4 | 8.01e-5 |
| GradDot | 3 | 1/4 | total | 4.77e-2 | 7.80e-3 | 2.56e-3 | 1.15e-3 | 6.30e-4 | 2.79e-4 | 3.23e-4 |
| GradDot | 3 | 1/4 | (b) | 1.36e-6 | 4.56e-6 | 1.01e-5 | 1.87e-5 | 3.14e-5 | 9.36e-5 | 2.43e-4 |
| GradDot | 3 | 1/4 | (c) | 4.77e-2 | 7.80e-3 | 2.55e-3 | 1.13e-3 | 5.99e-4 | 1.86e-4 | 8.01e-5 |
| GradDot | 3 | r0/l'max | total | 4.80e-2 | 7.90e-3 | 2.62e-3 | 1.18e-3 | 6.30e-4 | 2.00e-4 | 8.79e-5 |
| GradDot | 3 | r0/l'max | (b) | 2.29e-4 | 1.01e-4 | 7.15e-5 | 4.64e-5 | 3.14e-5 | 1.42e-5 | 7.77e-6 |
| GradDot | 3 | r0/l'max | (c) | 4.77e-2 | 7.80e-3 | 2.55e-3 | 1.13e-3 | 5.99e-4 | 1.86e-4 | 8.01e-5 |
| GradDot | 5 | 1/2 | total | 3.52e-2 | 1.54e-3 | 2.32e-4 | 5.77e-5 | 1.63e-5 | 2.12e-5 | 1.28e-4 |
| GradDot | 5 | 1/2 | (b) | 9.39e-8 | 3.52e-7 | 9.09e-7 | 2.04e-6 | 4.28e-6 | 2.41e-5 | 1.29e-4 |
| GradDot | 5 | 1/2 | (c) | 3.52e-2 | 1.54e-3 | 2.33e-4 | 5.98e-5 | 2.06e-5 | 2.91e-6 | 7.17e-7 |
| GradDot | 5 | 1/4 | total | 3.52e-2 | 1.54e-3 | 2.33e-4 | 5.98e-5 | 2.05e-5 | 2.75e-6 | 3.11e-7 |
| GradDot | 5 | 1/4 | (b) | 2.32e-9 | 7.80e-9 | 1.73e-8 | 3.19e-8 | 5.35e-8 | 1.58e-7 | 4.06e-7 |
| GradDot | 5 | 1/4 | (c) | 3.52e-2 | 1.54e-3 | 2.33e-4 | 5.98e-5 | 2.06e-5 | 2.91e-6 | 7.17e-7 |
| GradDot | 5 | r0/l'max | total | 3.52e-2 | 1.54e-3 | 2.33e-4 | 5.97e-5 | 2.05e-5 | 2.90e-6 | 7.14e-7 |
| GradDot | 5 | r0/l'max | (b) | 9.63e-6 | 1.01e-6 | 3.31e-7 | 1.23e-7 | 5.35e-8 | 1.09e-8 | 3.35e-9 |
| GradDot | 5 | r0/l'max | (c) | 3.52e-2 | 1.54e-3 | 2.33e-4 | 5.98e-5 | 2.06e-5 | 2.91e-6 | 7.17e-7 |
