# Stage 7: l'max optimum and shrinking worldtube

Scaled exterior error (max over l = 0..6 of max|error|/max|analytic total| over r in [1,3.5] U [6.5,10]); total = (b) + (c) from the validated Lemma-1 decomposition.

| case | n | R | argmin l'max | min total | total at l'max = 40 | behaviour |
|---|---|---|---|---|---|---|
| Product | 3 | 1/2 | 20 | 5.66e-4 | 3.50e-3 | not monotone |
| Product | 3 | 1/4 | 40 | 3.13e-5 | 3.13e-5 | monotone decreasing |
| Product | 3 | r0/l'max | 40 | 1.02e-5 | 1.02e-5 | monotone decreasing |
| Product | 5 | 1/2 | 30 | 4.20e-6 | 1.73e-5 | not monotone |
| Product | 5 | 1/4 | 40 | 8.34e-8 | 8.34e-8 | monotone decreasing |
| Product | 5 | r0/l'max | 40 | 1.12e-7 | 1.12e-7 | monotone decreasing |
| GradDot | 3 | 1/2 | 20 | 1.34e-1 | 2.77e0 | not monotone |
| GradDot | 3 | 1/4 | 30 | 2.60e-2 | 3.09e-2 | not monotone |
| GradDot | 3 | r0/l'max | 40 | 7.97e-3 | 7.97e-3 | monotone decreasing |
| GradDot | 5 | 1/2 | 20 | 1.37e-3 | 1.63e-2 | not monotone |
| GradDot | 5 | 1/4 | 40 | 2.52e-5 | 2.52e-5 | monotone decreasing |
| GradDot | 5 | r0/l'max | 40 | 6.47e-5 | 6.47e-5 | monotone decreasing |

## Errors against l'max

| case | n | R | quantity | 4 | 8 | 12 | 16 | 20 | 30 | 40 |
|---|---|---|---|---|---|---|---|---|---|---|
| Product | 3 | 1/2 | total | 5.43e0 | 2.93e-2 | 3.73e-3 | 1.04e-3 | 5.66e-4 | 1.02e-3 | 3.50e-3 |
| Product | 3 | 1/2 | (b) | 1.03e-5 | 3.52e-5 | 7.96e-5 | 1.53e-4 | 2.69e-4 | 9.75e-4 | 3.49e-3 |
| Product | 3 | 1/2 | (c) | 5.43e0 | 2.93e-2 | 3.65e-3 | 8.84e-4 | 2.96e-4 | 4.08e-5 | 9.97e-6 |
| Product | 3 | 1/4 | total | 5.43e0 | 2.93e-2 | 3.65e-3 | 8.86e-4 | 3.01e-4 | 5.12e-5 | 3.13e-5 |
| Product | 3 | 1/4 | (b) | 2.17e-7 | 7.11e-7 | 1.51e-6 | 2.66e-6 | 4.20e-6 | 1.03e-5 | 2.13e-5 |
| Product | 3 | 1/4 | (c) | 5.43e0 | 2.93e-2 | 3.65e-3 | 8.84e-4 | 2.96e-4 | 4.08e-5 | 9.97e-6 |
| Product | 3 | r0/l'max | total | 5.43e0 | 2.94e-2 | 3.68e-3 | 8.94e-4 | 3.01e-4 | 4.17e-5 | 1.02e-5 |
| Product | 3 | r0/l'max | (b) | 9.86e-4 | 1.09e-4 | 2.90e-5 | 9.98e-6 | 4.20e-6 | 8.22e-7 | 2.50e-7 |
| Product | 3 | r0/l'max | (c) | 5.43e0 | 2.93e-2 | 3.65e-3 | 8.84e-4 | 2.96e-4 | 4.08e-5 | 9.97e-6 |
| Product | 5 | 1/2 | total | 4.71e0 | 7.58e-3 | 4.26e-4 | 5.83e-5 | 1.15e-5 | 4.20e-6 | 1.73e-5 |
| Product | 5 | 1/2 | (b) | 5.35e-8 | 1.83e-7 | 4.14e-7 | 7.92e-7 | 1.40e-6 | 5.00e-6 | 1.74e-5 |
| Product | 5 | 1/2 | (c) | 4.71e0 | 7.58e-3 | 4.27e-4 | 5.91e-5 | 1.29e-5 | 8.04e-7 | 1.12e-7 |
| Product | 5 | 1/4 | total | 4.71e0 | 7.58e-3 | 4.27e-4 | 5.91e-5 | 1.28e-5 | 7.90e-7 | 8.34e-8 |
| Product | 5 | 1/4 | (b) | 2.89e-10 | 9.48e-10 | 2.02e-9 | 3.55e-9 | 5.60e-9 | 1.37e-8 | 2.83e-8 |
| Product | 5 | 1/4 | (c) | 4.71e0 | 7.58e-3 | 4.27e-4 | 5.91e-5 | 1.29e-5 | 8.04e-7 | 1.12e-7 |
| Product | 5 | r0/l'max | total | 4.71e0 | 7.58e-3 | 4.27e-4 | 5.91e-5 | 1.28e-5 | 8.03e-7 | 1.12e-7 |
| Product | 5 | r0/l'max | (b) | 3.71e-5 | 8.64e-7 | 1.06e-7 | 2.07e-8 | 5.60e-9 | 4.91e-10 | 8.41e-11 |
| Product | 5 | r0/l'max | (c) | 4.71e0 | 7.58e-3 | 4.27e-4 | 5.91e-5 | 1.29e-5 | 8.04e-7 | 1.12e-7 |
| GradDot | 3 | 1/2 | total | 5.14e0 | 8.25e-1 | 2.60e-1 | 1.43e-1 | 1.34e-1 | 4.92e-1 | 2.77e0 |
| GradDot | 3 | 1/2 | (b) | 1.67e-3 | 6.29e-3 | 1.64e-2 | 3.72e-2 | 7.93e-2 | 4.76e-1 | 2.76e0 |
| GradDot | 3 | 1/2 | (c) | 5.13e0 | 8.19e-1 | 2.44e-1 | 1.05e-1 | 5.51e-2 | 1.69e-2 | 7.25e-3 |
| GradDot | 3 | 1/4 | total | 5.13e0 | 8.19e-1 | 2.45e-1 | 1.07e-1 | 5.81e-2 | 2.60e-2 | 3.09e-2 |
| GradDot | 3 | 1/4 | (b) | 1.31e-4 | 4.42e-4 | 9.79e-4 | 1.81e-3 | 3.05e-3 | 9.09e-3 | 2.37e-2 |
| GradDot | 3 | 1/4 | (c) | 5.13e0 | 8.19e-1 | 2.44e-1 | 1.05e-1 | 5.51e-2 | 1.69e-2 | 7.25e-3 |
| GradDot | 3 | r0/l'max | total | 5.24e0 | 8.32e-1 | 2.51e-1 | 1.10e-1 | 5.81e-2 | 1.82e-2 | 7.97e-3 |
| GradDot | 3 | r0/l'max | (b) | 1.10e-1 | 1.35e-2 | 7.83e-3 | 4.69e-3 | 3.05e-3 | 1.33e-3 | 7.13e-4 |
| GradDot | 3 | r0/l'max | (c) | 5.13e0 | 8.19e-1 | 2.44e-1 | 1.05e-1 | 5.51e-2 | 1.69e-2 | 7.25e-3 |
| GradDot | 5 | 1/2 | total | 2.80e0 | 1.67e-1 | 2.24e-2 | 5.35e-3 | 1.37e-3 | 2.75e-3 | 1.63e-2 |
| GradDot | 5 | 1/2 | (b) | 1.15e-5 | 4.29e-5 | 1.11e-4 | 2.50e-4 | 5.27e-4 | 3.02e-3 | 1.64e-2 |
| GradDot | 5 | 1/2 | (c) | 2.80e0 | 1.67e-1 | 2.25e-2 | 5.60e-3 | 1.90e-3 | 2.65e-4 | 6.50e-5 |
| GradDot | 5 | 1/4 | total | 2.80e0 | 1.67e-1 | 2.25e-2 | 5.60e-3 | 1.90e-3 | 2.50e-4 | 2.52e-5 |
| GradDot | 5 | 1/4 | (b) | 2.27e-7 | 7.63e-7 | 1.69e-6 | 3.12e-6 | 5.23e-6 | 1.55e-5 | 3.97e-5 |
| GradDot | 5 | 1/4 | (c) | 2.80e0 | 1.67e-1 | 2.25e-2 | 5.60e-3 | 1.90e-3 | 2.65e-4 | 6.50e-5 |
| GradDot | 5 | r0/l'max | total | 2.79e0 | 1.66e-1 | 2.25e-2 | 5.59e-3 | 1.90e-3 | 2.64e-4 | 6.47e-5 |
| GradDot | 5 | r0/l'max | (b) | 5.06e-3 | 1.43e-4 | 3.71e-5 | 1.25e-5 | 5.23e-6 | 1.02e-6 | 3.08e-7 |
| GradDot | 5 | r0/l'max | (c) | 2.80e0 | 1.67e-1 | 2.25e-2 | 5.60e-3 | 1.90e-3 | 2.65e-4 | 6.50e-5 |
