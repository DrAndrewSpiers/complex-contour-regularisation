# Stage 12 (A): SS reuse of the hybrid technique

Scaled exterior SS error (max over l = 0..6) against the closed-form contour result.

| subtracted | p | h = 1/50 | 1/100 | 1/200 | order (1/100 -> 1/200) |
|---|---|---|---|---|---|
| {A} | 7 | 1.74e0 | 8.71e-1 | 4.36e-1 | 1. |
| {AB} | 7 | 1.24e0 | 6.21e-1 | 3.11e-1 | 1. |
| {ABC} | 7 | 1.96e-5 | 2.45e-6 | 3.06e-7 | 3. |
| {ABC+} | 7 | 3.18e-9 | 9.93e-11 | 3.16e-12 | 4.97 |
| {ABC+} | 3 | 9.08e-8 | 2.21e-9 | 1.28e-11 | 7.43 |

## noise on the remainder samples

| subtracted | h | noise model | eps | total error | noise part |
|---|---|---|---|---|---|
| {ABC} | 1/50 | rel. remainder | 1.00e-12 | 1.96e-5 | 1.27e-11 |
| {ABC} | 1/50 | rel. remainder | 1.00e-8 | 1.95e-5 | 1.27e-7 |
| {ABC} | 1/50 | rel. full SS | 1.00e-12 | 1.96e-5 | 1.90e-8 |
| {ABC} | 1/50 | rel. full SS | 1.00e-8 | 1.94e-4 | 1.90e-4 |
| {ABC} | 1/100 | rel. remainder | 1.00e-12 | 2.45e-6 | 5.51e-12 |
| {ABC} | 1/100 | rel. remainder | 1.00e-8 | 2.40e-6 | 5.57e-8 |
| {ABC} | 1/100 | rel. full SS | 1.00e-12 | 2.36e-6 | 8.97e-8 |
| {ABC} | 1/100 | rel. full SS | 1.00e-8 | 8.95e-4 | 8.97e-4 |
| {ABC} | 1/200 | rel. remainder | 1.00e-12 | 3.06e-7 | 4.91e-12 |
| {ABC} | 1/200 | rel. remainder | 1.00e-8 | 3.55e-7 | 4.87e-8 |
| {ABC} | 1/200 | rel. full SS | 1.00e-12 | 3.79e-7 | 7.23e-8 |
| {ABC} | 1/200 | rel. full SS | 1.00e-8 | 7.23e-4 | 7.23e-4 |
| {ABC+} | 1/50 | rel. remainder | 1.00e-12 | 3.17e-9 | 1.10e-11 |
| {ABC+} | 1/50 | rel. remainder | 1.00e-8 | 1.07e-7 | 1.10e-7 |
| {ABC+} | 1/50 | rel. full SS | 1.00e-12 | 1.95e-8 | 1.90e-8 |
| {ABC+} | 1/50 | rel. full SS | 1.00e-8 | 1.90e-4 | 1.90e-4 |
| {ABC+} | 1/100 | rel. remainder | 1.00e-12 | 9.17e-11 | 7.53e-12 |
| {ABC+} | 1/100 | rel. remainder | 1.00e-8 | 7.47e-8 | 7.48e-8 |
| {ABC+} | 1/100 | rel. full SS | 1.00e-12 | 8.96e-8 | 8.97e-8 |
| {ABC+} | 1/100 | rel. full SS | 1.00e-8 | 8.97e-4 | 8.97e-4 |
| {ABC+} | 1/200 | rel. remainder | 1.00e-12 | 8.65e-12 | 5.49e-12 |
| {ABC+} | 1/200 | rel. remainder | 1.00e-8 | 5.50e-8 | 5.50e-8 |
| {ABC+} | 1/200 | rel. full SS | 1.00e-12 | 7.23e-8 | 7.23e-8 |
| {ABC+} | 1/200 | rel. full SS | 1.00e-8 | 7.23e-4 | 7.23e-4 |

closed-form contour vs exact (max over l, scaled): <1e-27; Lemma-1 check at aR (l = 2, ABC, h = 1/100), contour vs real line: 1.54e-6 relative.
