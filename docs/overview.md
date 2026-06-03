# Project Overview

## Question

This repository studies confidence intervals for the power of the one-sample t-test when the population standard deviation is unknown.

The specific methodological question is:

**How does the analysis change when interval width is optimized by the median rather than the mean?**

## Statistical Setting

Power depends on the unknown standard deviation, so post-data statements about power inherit uncertainty from estimating that quantity. The analysis target is the one-sample upper-tailed test,

`H0: mu = mu0` versus `H1: mu > mu0`,

with a power interval for a pre-specified raw shift `Delta = mu - mu0`. The interval is written first in observable form using the realized sample standard deviation `S`, and only then re-indexed with `kappa = Delta / sigma` for scale-free design tabulation.

The exact-coverage guarantee depends on one structural requirement: the chi-squared endpoint allocation must be fixed before the observed sample standard deviation is seen. If that allocation depends on the realized data, exact coverage can fail.

## Procedures Compared

- `Equal-tail`: baseline allocation
- `Expected-opt`: minimizes expected interval length
- `Median-opt`: minimizes median interval length
- `Median-MC`: Monte Carlo approximation to the same median-length target

The median criterion changes the optimality target, not the coverage argument.

## Main Result

The median-length construction is theoretically coherent and empirically interpretable, especially when realized interval widths are skewed. At the same time, the expected-length rule is already close to the median benchmark on the design grid used in this project. The main contribution is therefore calibration and comparison rather than a large numerical improvement.

## Interpretation

The exact-coverage statement is pointwise in a prespecified standardized shift index:

`kappa = Delta / sigma`

Following Chakraborti et al., this is a scale-free design convention, not a claim that the observed-data analysis knows the true population standard deviation. The guarantee applies when endpoint allocation is fixed in advance; the repository does not claim exact post-hoc inference for data-adaptive allocation rules.

## Repository Structure

- [paper/paper.tex](../paper/paper.tex): manuscript source
- [paper/paper.pdf](../paper/paper.pdf): compiled manuscript
- [scripts/median_power_ci.py](../scripts/median_power_ci.py): full numerical pipeline
- [scripts/posthoc_allocation_failure_demo.py](../scripts/posthoc_allocation_failure_demo.py): cautionary failure demonstration for invalid post-hoc allocation
- [results/](../results): generated tables and machine-readable summaries
- [figures/](../figures): generated figures

## Reading Guide

1. [paper/paper.pdf](../paper/paper.pdf)
2. [scripts/median_power_ci.py](../scripts/median_power_ci.py)
3. [results/summary_grid.csv](../results/summary_grid.csv)

## Repository Principles

- precise scope statements
- reproducible outputs
- clear separation between theory, numerical optimization, and validation
- compact project structure with a single canonical manuscript build path
