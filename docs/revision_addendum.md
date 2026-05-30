# Revision Addendum

This memo records manuscript-facing corrections for the public release and revision package.

## Core Clarifications

- `gamma` denotes noncoverage probability throughout.
- Confidence level is `1 - gamma`; for example, `1 - gamma in {0.90, 0.95, 0.99}` is equivalent to `gamma in {0.10, 0.05, 0.01}`.
- The inferential target is the power of the one-sample upper-tailed t-test:
  - `H0: mu = mu0`
  - `H1: mu > mu0`
- The power interval is written first in observable form using the realized sample standard deviation `S` and a pre-specified raw shift `Delta = mu - mu0`.
- The scale-free index `kappa = Delta / sigma` is used for design-stage tabulation of the random length distribution.
- `kappa` is not treated as an observed post-data effect estimate.
- Endpoint allocation is fixed before observing `S`.

Preferred wording used in manuscript-facing materials:

> Following Chakraborti et al., the numerical design problem is indexed by a standardized shift. This is a scale-free tabulation convention, not an assumption that the population standard deviation is known in the data analysis. The endpoint allocation is fixed before observing the sample standard deviation.

## Observable Interval Form

For fixed chi-square endpoints `A < B` satisfying `P(A <= Q <= B) = 1 - gamma`, where `Q = nu S^2 / sigma^2`, the induced interval for `sigma` is

```text
sqrt(nu S^2 / B) <= sigma <= sqrt(nu S^2 / A).
```

For a fixed positive raw shift `Delta`, power is decreasing in `sigma`, so the observable power interval is written as

```text
[ pi((Delta / S) sqrt(A / nu)), pi((Delta / S) sqrt(B / nu)) ].
```

The scale-free length calculations then use

```text
kappa sqrt(A / Q)
kappa sqrt(B / Q)
```

with the square roots included.

## Method Relationship To Chakraborti Et Al.

- Architecture preserved: fixed allocation `p` -> fixed chi-square endpoints `(A_p, B_p)` -> exact coverage.
- Criterion changed only at design stage:
  - expected-length objective in Chakraborti et al.
  - median-length objective in this note
- Numerical interpretation: median optimization gives modest gains on the studied grid, and expected-length is often nearly median-optimal.

## Monte Carlo Role

- Monte Carlo is used only as a design-stage numerical approximation and sensitivity check.
- Exact coverage comes from the fixed chi-square pivot construction, not from simulation.
- The post-hoc data-adaptive allocation procedure is included only as a cautionary demonstration of why endpoint allocation must be fixed before observing `S`.

## Boundary And Flatness Diagnostics

- Boundary classes are relative to `gamma`:
  - lower boundary if `p <= 0.01 * gamma`
  - upper boundary if `p >= 0.99 * gamma`
  - interior otherwise
- For flat or nonunique numerical median optima, the implementation uses a deterministic tie-break rule: choose the smallest near-minimizing allocation on a fixed deterministic grid.

## Scope Controls

- Main paper scope remains the one-sample one-sided t-test note.
- Broader topics, including two-sided procedures, noncentrality-parameter intervals, F-test families, expanded grids, and alternative risk criteria, are treated as sequel, comparison, or future work rather than part of the core claim.

## Citation Note

- Verify citation details before submission for any newly added references or metadata updates.
