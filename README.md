# Complex Contour Regularisation

Mathematica notebooks and working documents for the **complex contour regularisation** project — a method for solving second-order perturbation theory problems in gravitational self-force by deforming the integration contour into the complex plane.

## Project Overview

The core idea is to perform a **singular–regular split of the first-order field in the complex plane**, enabling the second-order source to be decomposed into three pieces:

| Piece | Description |
|-------|-------------|
| **Singular–Singular (SS)** | Handled using 3D formulas without mode decomposition |
| **Singular–Regular (SR)** | Source computed with 3D singular-field expressions combined with a first-order mode decomposition for the regular piece |
| **Regular–Regular (RR)** | Computed via first-order mode decomposition (optionally solved on the real line) |

A worldtube construction is employed: outside the worldtube the standard Poisson equation is solved on the real line; inside the worldtube integration proceeds along a complex contour using an extended Dirac delta as source. The regular field is isolated via an effective source that subtracts the singular-field contribution.

## Repository Structure

```
.
├── Andrew_s_idea.pdf                    # Technical write-up of the approach
├── Tommy email.rtf                      # Correspondence describing the method
├── Heaviside-works/
│   ├── Does-Heaviside-delta-product-appear.nb
│   ├── Integrate-Dirac-Delta-to-Heaviside.nb
│   ├── Kink-HLog-ExtendedDiracDelta.nb
│   └── test-heaviside-from-numeric.nb
└── Modes-of-Source-and-solution/
    ├── Modes-of-Source-and-solution-1or.nb
    └── Modes-of-Source-and-solution-1or4.nb
```

### Notebook Descriptions

**Heaviside-works/**
- `Does-Heaviside-delta-product-appear.nb` — Investigates whether a Heaviside–delta product appears in the regularisation procedure.
- `Integrate-Dirac-Delta-to-Heaviside.nb` — Demonstrates integration of the extended Dirac delta, recovering the Heaviside step function.
- `Kink-HLog-ExtendedDiracDelta.nb` — Studies kink and logarithmic terms arising from the extended Dirac delta formulation.
- `test-heaviside-from-numeric.nb` — Numerical tests of the Heaviside recovery.

**Modes-of-Source-and-solution/**
- `Modes-of-Source-and-solution-1or.nb` — Mode decomposition of the source and solution at first order.
- `Modes-of-Source-and-solution-1or4.nb` — Extended version with higher-accuracy mode decomposition.

### Related Files

- **`Andrew_s_idea.pdf`** — Technical write-up detailing the singular–regular split and contour approach.
- **`Tommy email.rtf`** — Email to collaborator outlining the strategy for the SS, SR, and RR source decomposition.

> *Also see `poisson_spherical_harmonics.nb` (in the Test-project directory) for related spherical harmonic Poisson equation work.*

## Key Results

The notebooks demonstrate that the **regular part of the first-order field can be recovered to good accuracy** on a complex contour using an effective source derived from a truncated singular-field expansion. Parameters such as the multipole index ℓ, worldtube amplitude and size, and the order of the singular-field expansion are freely adjustable.

## Authors

- Andrew Spiers (UCD) — andrew.spiers@ucd.ie

## Collaborators

- Tommy (co-author on the complex contour project)
