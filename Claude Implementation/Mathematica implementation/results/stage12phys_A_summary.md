# Stage 12 (A): SS reuse of the hybrid technique

Scaled exterior SS error (max over l = 0..6) against the closed-form contour result.

| subtracted | p | h = 1/50 | 1/100 | 1/200 | order (1/100 -> 1/200) |
|---|---|---|---|---|---|
| {A} | 7 | 3.53e-3 | 1.76e-3 | 8.82e-4 | 1. |
| {AB} | 7 | 3.09e-3 | 1.55e-3 | 7.74e-4 | 1. |
| {ABC} | 7 | 4.17e-8 | 5.22e-9 | 6.52e-10 | 3. |
| {ABC+} | 7 | 7.95e-12 | 2.48e-13 | 7.87e-15 | 4.98 |
| {ABC+} | 3 | 1.67e-10 | 4.07e-12 | 5.64e-13 | 2.85 |

## noise on the remainder samples

| subtracted | h | noise model | eps | total error | noise part |
|---|---|---|---|---|---|
| {ABC} | 1/50 | rel. remainder | 1.00e-12 | 4.17e-8 | 2.21e-13 |
| {ABC} | 1/50 | rel. remainder | 1.00e-8 | 4.14e-8 | 2.21e-9 |
| {ABC} | 1/50 | rel. full SS | 1.00e-12 | 4.17e-8 | 3.44e-10 |
| {ABC} | 1/50 | rel. full SS | 1.00e-8 | 3.44e-6 | 3.44e-6 |
| {ABC} | 1/100 | rel. remainder | 1.00e-12 | 5.22e-9 | 4.20e-14 |
| {ABC} | 1/100 | rel. remainder | 1.00e-8 | 5.01e-9 | 4.21e-10 |
| {ABC} | 1/100 | rel. full SS | 1.00e-12 | 4.99e-9 | 1.60e-9 |
| {ABC} | 1/100 | rel. full SS | 1.00e-8 | 1.60e-5 | 1.60e-5 |
| {ABC} | 1/200 | rel. remainder | 1.00e-12 | 6.52e-10 | 1.76e-14 |
| {ABC} | 1/200 | rel. remainder | 1.00e-8 | 8.28e-10 | 1.76e-10 |
| {ABC} | 1/200 | rel. full SS | 1.00e-12 | 4.01e-9 | 4.06e-9 |
| {ABC} | 1/200 | rel. full SS | 1.00e-8 | 4.06e-5 | 4.06e-5 |
| {ABC+} | 1/50 | rel. remainder | 1.00e-12 | 7.93e-12 | 2.22e-13 |
| {ABC+} | 1/50 | rel. remainder | 1.00e-8 | 2.22e-9 | 2.22e-9 |
| {ABC+} | 1/50 | rel. full SS | 1.00e-12 | 3.44e-10 | 3.44e-10 |
| {ABC+} | 1/50 | rel. full SS | 1.00e-8 | 3.44e-6 | 3.44e-6 |
| {ABC+} | 1/100 | rel. remainder | 1.00e-12 | 2.16e-13 | 4.38e-14 |
| {ABC+} | 1/100 | rel. remainder | 1.00e-8 | 4.36e-10 | 4.36e-10 |
| {ABC+} | 1/100 | rel. full SS | 1.00e-12 | 1.60e-9 | 1.60e-9 |
| {ABC+} | 1/100 | rel. full SS | 1.00e-8 | 1.60e-5 | 1.60e-5 |
| {ABC+} | 1/200 | rel. remainder | 1.00e-12 | 2.89e-14 | 2.10e-14 |
| {ABC+} | 1/200 | rel. remainder | 1.00e-8 | 2.09e-10 | 2.09e-10 |
| {ABC+} | 1/200 | rel. full SS | 1.00e-12 | 4.06e-9 | 4.06e-9 |
| {ABC+} | 1/200 | rel. full SS | 1.00e-8 | 4.06e-5 | 4.06e-5 |

closed-form contour vs exact (max over l, scaled): <1e-30; Lemma-1 check at aR (l = 2, ABC, h = 1/100), contour vs real line: 1.54e-6 relative.
