from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import median_power_ci as core


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"


def parse_cases(value: str) -> list[tuple[int, float]]:
    pairs: list[tuple[int, float]] = []
    for item in value.split(","):
        raw = item.strip()
        if not raw:
            continue
        n_text, kappa_text = raw.split(":")
        pairs.append((int(n_text), float(kappa_text)))
    return pairs


def posthoc_coverage_for_case(
    n: int,
    kappa: float,
    reps: int,
    p_grid_size: int,
    seed: int,
) -> dict[str, float]:
    rng = np.random.default_rng(seed + n)
    q_sample = rng.chisquare(df=n - 1, size=reps)

    true_power = float(core.power_one_sample_t(n, core.ALPHA, kappa))

    # Valid design-stage interval (fixed p selected before observing S).
    fixed_spec = core.optimize_interval(n, kappa, "median")
    valid_lower, valid_upper = core.interval_from_x(q_sample, n, core.ALPHA, kappa, fixed_spec.A, fixed_spec.B)
    valid_covered = (valid_lower <= true_power) & (true_power <= valid_upper)
    valid_coverage = float(np.mean(valid_covered))

    # Invalid post-hoc rule: after seeing Q (equivalently S), choose p that minimizes
    # the realized interval length for that specific observation.
    p_candidates = np.linspace(1e-8, core.GAMMA - 1e-8, p_grid_size)
    specs = [core.interval_spec(n, float(p)) for p in p_candidates]
    A_vals = np.array([spec.A for spec in specs])
    B_vals = np.array([spec.B for spec in specs])

    posthoc_lower = np.empty(reps, dtype=float)
    posthoc_upper = np.empty(reps, dtype=float)

    for i, q in enumerate(q_sample):
        kappa_a = kappa * np.sqrt(A_vals / q)
        kappa_b = kappa * np.sqrt(B_vals / q)
        lower_all = core.power_one_sample_t(n, core.ALPHA, kappa_a)
        upper_all = core.power_one_sample_t(n, core.ALPHA, kappa_b)
        lengths = upper_all - lower_all
        best_idx = int(np.argmin(lengths))
        posthoc_lower[i] = float(lower_all[best_idx])
        posthoc_upper[i] = float(upper_all[best_idx])

    posthoc_covered = (posthoc_lower <= true_power) & (true_power <= posthoc_upper)
    posthoc_coverage = float(np.mean(posthoc_covered))

    return {
        "n": float(n),
        "kappa": float(kappa),
        "alpha": float(core.ALPHA),
        "gamma": float(core.GAMMA),
        "reps": float(reps),
        "p_grid_size": float(p_grid_size),
        "fixed_left_tail": float(fixed_spec.left_tail),
        "valid_fixed_coverage": valid_coverage,
        "invalid_posthoc_coverage": posthoc_coverage,
        "coverage_gap_posthoc_minus_valid": posthoc_coverage - valid_coverage,
    }


def run(cases: list[tuple[int, float]], reps: int, p_grid_size: int, seed: int) -> pd.DataFrame:
    rows = [
        posthoc_coverage_for_case(
            n=n,
            kappa=kappa,
            reps=reps,
            p_grid_size=p_grid_size,
            seed=seed,
        )
        for n, kappa in cases
    ]
    return pd.DataFrame(rows)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Cautionary demonstration: choosing the endpoint allocation after observing S "
            "can break exact coverage."
        )
    )
    parser.add_argument(
        "--cases",
        default="4:0.75,8:1.00",
        help="Comma-separated n:kappa pairs, for example '4:0.75,8:1.00'.",
    )
    parser.add_argument("--reps", type=int, default=30000, help="Monte Carlo replications per case.")
    parser.add_argument("--p-grid-size", type=int, default=81, help="Number of candidate p values for post-hoc selection.")
    parser.add_argument("--seed", type=int, default=20260530, help="RNG seed.")
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    cases = parse_cases(args.cases)
    if not cases:
        raise ValueError("No valid cases were provided.")

    df = run(cases=cases, reps=args.reps, p_grid_size=args.p_grid_size, seed=args.seed)
    RESULTS_DIR.mkdir(exist_ok=True)
    out_path = RESULTS_DIR / "posthoc_allocation_demo.csv"
    df.to_csv(out_path, index=False)

    print("Saved:")
    print(f"  {out_path}")


if __name__ == "__main__":
    main()
