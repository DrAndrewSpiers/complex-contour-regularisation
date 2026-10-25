# Stage 6: consistent 3D puncture

Scaled exterior errors (max over l = 0..6 of max|error|/max|analytic total| over [1,4.5] U [5.5,10]); Cauchy data at r = 1; K continued by Theta->H.

| case | n | c | (b) shadow | (a) mismatch | total error | (c) l' truncation (l'max = 12) | (b)/total |
|---|---|---|---|---|---|---|---|
| Product | 1 | 1/25 | 1.21e-2 | 0 | 4.45e-2 | 5.65e-2 | 2.71e-1 |
| Product | 3 | 1/625 | 4.59e-5 | 0 | 2.15e-3 | 2.10e-3 | 2.13e-2 |
| Product | 5 | 1/15625 | 2.38e-7 | 0 | 2.46e-4 | 2.46e-4 | 9.70e-4 |
| Product | 7 | 1/390625 | 1.45e-9 | 0 | 6.20e-5 | 6.20e-5 | 2.34e-5 |
| Product | 3 | 2/125 | 4.59e-4 | 0 | 2.15e-2 | 2.10e-2 | 2.13e-2 |
| GradDot | 1 | 1/25 | 9.25e-1 | 0 | 3.01e0 | 3.94e0 | 3.07e-1 |
| GradDot | 3 | 1/625 | 4.91e-3 | 0 | 7.80e-2 | 7.31e-2 | 6.30e-2 |
| GradDot | 5 | 1/15625 | 3.34e-5 | 0 | 6.73e-3 | 6.77e-3 | 4.96e-3 |
| GradDot | 7 | 1/390625 | 2.52e-7 | 0 | 1.49e-3 | 1.49e-3 | 1.69e-4 |
| GradDot | 3 | 2/125 | 4.92e-2 | 0 | 7.81e-1 | 7.31e-1 | 6.30e-2 |

Fitted decay factor of (b) per unit n: {{Product, 0.070455003683158876884919085406}, {GradDot, 0.08071748058763905463751915636}}

(b)(10c)/(b)(c) at n = 3 (max, min over l): {{Product, 10.0009492435998591244950902402, 10.0009422678885122150793284458}, {GradDot, 10.0029293610504491983769536843, 10.0027135752815325258057405106}}

## l'max convergence (n = 3)

| case | l'max | max_l total | max_l (c) | max_l (b) |
|---|---|---|---|---|
| Product | 2 | 5.77e0 | 5.77e0 | 2.10e-6 |
| Product | 4 | 3.13e0 | 3.13e0 | 5.93e-6 |
| Product | 8 | 1.69e-2 | 1.69e-2 | 2.02e-5 |
| Product | 12 | 2.15e-3 | 2.10e-3 | 4.59e-5 |
| GradDot | 2 | 3.22e0 | 3.22e0 | 1.73e-4 |
| GradDot | 4 | 1.54e0 | 1.54e0 | 5.02e-4 |
| GradDot | 8 | 2.48e-1 | 2.46e-1 | 1.89e-3 |
| GradDot | 12 | 7.80e-2 | 7.31e-2 | 4.91e-3 |

## (b) per mode l' and cumulative, beyond l'max = 12 (decomposition only; n = 3)

| case | R | l' | (b) contribution of mode l' | cumulative (b), l'max = l' |
|---|---|---|---|---|
| Product | 1/2 | 0 | 2.32e-7 | 2.32e-7 |
| Product | 1/2 | 4 | 2.17e-6 | 5.93e-6 |
| Product | 1/2 | 8 | 4.50e-6 | 2.02e-5 |
| Product | 1/2 | 12 | 7.69e-6 | 4.59e-5 |
| Product | 1/2 | 16 | 1.24e-5 | 8.78e-5 |
| Product | 1/2 | 20 | 1.98e-5 | 1.55e-4 |
| Product | 1/2 | 25 | 3.59e-5 | 2.98e-4 |
| Product | 1/2 | 30 | 6.64e-5 | 5.61e-4 |
| Product | 1/4 | 0 | 4.97e-9 | 4.97e-9 |
| Product | 1/4 | 4 | 4.51e-8 | 1.25e-7 |
| Product | 1/4 | 8 | 8.72e-8 | 4.10e-7 |
| Product | 1/4 | 12 | 1.33e-7 | 8.72e-7 |
| Product | 1/4 | 16 | 1.85e-7 | 1.53e-6 |
| Product | 1/4 | 20 | 2.45e-7 | 2.42e-6 |
| Product | 1/4 | 25 | 3.38e-7 | 3.92e-6 |
| Product | 1/4 | 30 | 4.55e-7 | 5.94e-6 |
| GradDot | 1/2 | 0 | 1.88e-5 | 1.88e-5 |
| GradDot | 1/2 | 4 | 1.89e-4 | 5.02e-4 |
| GradDot | 1/2 | 8 | 4.61e-4 | 1.89e-3 |
| GradDot | 1/2 | 12 | 9.77e-4 | 4.91e-3 |
| GradDot | 1/2 | 16 | 2.00e-3 | 1.12e-2 |
| GradDot | 1/2 | 20 | 4.03e-3 | 2.38e-2 |
| GradDot | 1/2 | 25 | 9.64e-3 | 5.89e-2 |
| GradDot | 1/2 | 30 | 2.31e-2 | 1.43e-1 |
| GradDot | 1/4 | 0 | 1.55e-6 | 1.55e-6 |
| GradDot | 1/4 | 4 | 1.43e-5 | 3.94e-5 |
| GradDot | 1/4 | 8 | 2.91e-5 | 1.33e-4 |
| GradDot | 1/4 | 12 | 4.76e-5 | 2.94e-4 |
| GradDot | 1/4 | 16 | 7.23e-5 | 5.44e-4 |
| GradDot | 1/4 | 20 | 1.06e-4 | 9.15e-4 |
| GradDot | 1/4 | 25 | 1.68e-4 | 1.62e-3 |
| GradDot | 1/4 | 30 | 2.63e-4 | 2.73e-3 |

Per-l' growth factor near l' = 30:             1                    1                      1                    1
{{{Product, -, 1.131}, {Product, -, 1.062}}, {{GradDot, -, 1.191}, {GradDot, -, 1.093}}}
            2                    4                      2                    4
