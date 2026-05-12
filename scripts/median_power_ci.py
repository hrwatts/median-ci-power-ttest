from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize, stats


ALPHA = 0.05
GAMMA = 0.05
GRID_SIZE = 4001
MC_REPS = 100_000
RNG_SEED = 20260420
MC_OPT_REPS = 4_000
MC_OPT_SEED = 20260421
MC_SENSITIVITY_CONFIGS = [
    ("$R=1000$", 1_000, 0),
    ("$R=4000$", 4_000, 0),
    ("$R=4000$, alt.", 4_000, 500),
    ("$R=16000$", 16_000, 0),
]

N_VALUES = [2, 4, 8, 16, 25]
KAPPA_VALUES = [0.25, 0.50, 0.75, 1.00, 1.25]
VALIDATION_CASES = [(2, 0.50), (4, 0.75), (8, 1.00), (25, 1.25)]


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = ROOT / "figures"


@dataclass(frozen=True)
class IntervalSpec:
    method: str
    left_tail: float
    A: float
    B: float


def ensure_directories() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    FIGURES_DIR.mkdir(exist_ok=True)


def power_one_sample_t(n: int, alpha: float, kappa: np.ndarray | float) -> np.ndarray | float:
    df = n - 1
    tcrit = stats.t.ppf(1 - alpha, df)
    delta = np.sqrt(n) * np.asarray(kappa)
    return 1.0 - stats.nct.cdf(tcrit, df, delta)


def length_given_x(
    x: np.ndarray,
    n: int,
    alpha: float,
    kappa: float,
    A: float,
    B: float,
) -> np.ndarray:
    lower = power_one_sample_t(n, alpha, kappa * np.sqrt(A / x))
    upper = power_one_sample_t(n, alpha, kappa * np.sqrt(B / x))
    return upper - lower


def interval_from_x(
    x: np.ndarray,
    n: int,
    alpha: float,
    kappa: float,
    A: float,
    B: float,
) -> tuple[np.ndarray, np.ndarray]:
    lower = power_one_sample_t(n, alpha, kappa * np.sqrt(A / x))
    upper = power_one_sample_t(n, alpha, kappa * np.sqrt(B / x))
    return lower, upper


def evaluation_grid(n: int, grid_size: int = GRID_SIZE) -> np.ndarray:
    probs = (np.arange(grid_size) + 0.5) / grid_size
    return stats.chi2.ppf(probs, n - 1)


def interval_spec(n: int, left_tail: float) -> IntervalSpec:
    df = n - 1
    p_left = max(left_tail, 1e-12)
    p_right = min(left_tail + 1.0 - GAMMA, 1.0 - 1e-12)
    return IntervalSpec(
        method="custom",
        left_tail=left_tail,
        A=stats.chi2.ppf(p_left, df),
        B=stats.chi2.ppf(p_right, df),
    )


def optimize_interval(n: int, kappa: float, objective: str) -> IntervalSpec:
    x_grid = evaluation_grid(n)

    def objective_fn(left_tail: float) -> float:
        spec = interval_spec(n, left_tail)
        lengths = length_given_x(x_grid, n, ALPHA, kappa, spec.A, spec.B)
        if objective == "mean":
            return float(np.mean(lengths))
        if objective == "median":
            return float(np.median(lengths))
        raise ValueError(f"Unknown objective: {objective}")

    result = optimize.minimize_scalar(
        objective_fn,
        bounds=(1e-10, GAMMA - 1e-8),
        method="bounded",
        options={"xatol": 1e-6},
    )
    spec = interval_spec(n, float(result.x))
    return IntervalSpec(
        method=objective,
        left_tail=float(result.x),
        A=spec.A,
        B=spec.B,
    )


def optimize_interval_mc(n: int, kappa: float, x_sample: np.ndarray) -> IntervalSpec:
    def objective_fn(left_tail: float) -> float:
        spec = interval_spec(n, left_tail)
        lengths = length_given_x(x_sample, n, ALPHA, kappa, spec.A, spec.B)
        return float(np.median(lengths))

    result = optimize.minimize_scalar(
        objective_fn,
        bounds=(1e-10, GAMMA - 1e-8),
        method="bounded",
        options={"xatol": 1e-4},
    )
    spec = interval_spec(n, float(result.x))
    return IntervalSpec(
        method="median_mc",
        left_tail=float(result.x),
        A=spec.A,
        B=spec.B,
    )


def method_specs(n: int, kappa: float, x_mc_opt: np.ndarray | None = None) -> dict[str, IntervalSpec]:
    equal_tail = interval_spec(n, GAMMA / 2.0)
    expected = optimize_interval(n, kappa, "mean")
    median = optimize_interval(n, kappa, "median")
    specs = {
        "Equal-tail": IntervalSpec("Equal-tail", equal_tail.left_tail, equal_tail.A, equal_tail.B),
        "Expected-opt": IntervalSpec("Expected-opt", expected.left_tail, expected.A, expected.B),
        "Median-opt": IntervalSpec("Median-opt", median.left_tail, median.A, median.B),
    }
    if x_mc_opt is not None:
        median_mc = optimize_interval_mc(n, kappa, x_mc_opt)
        specs["Median-MC"] = IntervalSpec("Median-MC", median_mc.left_tail, median_mc.A, median_mc.B)
    return specs


def summarize_case(n: int, kappa: float, x_mc_opt: np.ndarray) -> tuple[list[dict[str, float]], dict[str, float]]:
    x_grid = evaluation_grid(n)
    specs = method_specs(n, kappa, x_mc_opt=x_mc_opt)
    performance_rows: list[dict[str, float]] = []
    summary_row: dict[str, float] = {
        "n": n,
        "kappa": kappa,
        "true_power": float(power_one_sample_t(n, ALPHA, kappa)),
    }

    baseline_median = None
    expected_median = None

    for label, spec in specs.items():
        lengths = length_given_x(x_grid, n, ALPHA, kappa, spec.A, spec.B)
        mean_length = float(np.mean(lengths))
        median_length = float(np.median(lengths))
        p90_length = float(np.quantile(lengths, 0.90))
        mean_to_median = mean_length / median_length

        performance_rows.append(
            {
                "n": n,
                "kappa": kappa,
                "method": label,
                "left_tail": spec.left_tail,
                "A": spec.A,
                "B": spec.B,
                "mean_length": mean_length,
                "median_length": median_length,
                "p90_length": p90_length,
                "mean_to_median": mean_to_median,
                "true_power": summary_row["true_power"],
            }
        )

        if label == "Equal-tail":
            baseline_median = median_length
            summary_row["mean_equal_tail"] = mean_length
            summary_row["median_equal_tail"] = median_length
        if label == "Expected-opt":
            expected_median = median_length
            summary_row["left_tail_expected"] = spec.left_tail
            summary_row["mean_expected"] = mean_length
            summary_row["median_expected"] = median_length
            summary_row["ratio_expected"] = mean_to_median
        if label == "Median-opt":
            summary_row["left_tail_median"] = spec.left_tail
            summary_row["mean_median"] = mean_length
            summary_row["median_median"] = median_length
        if label == "Median-MC":
            summary_row["left_tail_median_mc"] = spec.left_tail
            summary_row["mean_median_mc"] = mean_length
            summary_row["median_median_mc"] = median_length

    assert baseline_median is not None
    assert expected_median is not None
    summary_row["median_gain_vs_equal_tail_pct"] = (
        100.0
        * (baseline_median - summary_row["median_median"])
        / baseline_median
    )
    summary_row["median_gain_vs_expected_pct"] = (
        100.0
        * (expected_median - summary_row["median_median"])
        / expected_median
    )
    summary_row["median_gain_mc_vs_expected_pct"] = (
        100.0
        * (expected_median - summary_row["median_median_mc"])
        / expected_median
    )
    summary_row["left_tail_shift_pct_points"] = 100.0 * (
        summary_row["left_tail_expected"] - summary_row["left_tail_median"]
    )
    summary_row["left_tail_shift_mc_pct_points"] = 100.0 * (
        summary_row["left_tail_expected"] - summary_row["left_tail_median_mc"]
    )
    return performance_rows, summary_row


def monte_carlo_validation(performance_df: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(RNG_SEED)
    rows: list[dict[str, float]] = []
    for n in N_VALUES:
        x = rng.chisquare(df=n - 1, size=MC_REPS)
        for kappa in KAPPA_VALUES:
            true_power = float(power_one_sample_t(n, ALPHA, kappa))
            for method in ["Equal-tail", "Expected-opt", "Median-opt", "Median-MC"]:
                row = performance_df[
                    (performance_df["n"] == n)
                    & (performance_df["kappa"] == kappa)
                    & (performance_df["method"] == method)
                ].iloc[0]
                lower, upper = interval_from_x(x, n, ALPHA, kappa, row["A"], row["B"])
                lengths = upper - lower
                covered = (lower <= true_power) & (true_power <= upper)
                coverage = float(np.mean(covered))
                se = float(np.sqrt(coverage * (1.0 - coverage) / MC_REPS))
                rows.append(
                    {
                        "n": n,
                        "kappa": kappa,
                        "method": method,
                        "coverage": coverage,
                        "coverage_se": se,
                        "mean_length": float(np.mean(lengths)),
                        "median_length": float(np.median(lengths)),
                        "true_power": true_power,
                    }
                )
    return pd.DataFrame(rows)


def mc_sensitivity_study(summary_df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for n, kappa in VALIDATION_CASES:
        deterministic_row = summary_df[
            (summary_df["n"] == n) & (summary_df["kappa"] == kappa)
        ].iloc[0]
        p_m = float(deterministic_row["left_tail_median"])
        median_ref = float(deterministic_row["median_median"])
        x_grid = evaluation_grid(n, grid_size=10_001)
        for label, reps, seed_offset in MC_SENSITIVITY_CONFIGS:
            x_sample = np.random.default_rng(MC_OPT_SEED + seed_offset + n).chisquare(
                df=n - 1,
                size=reps,
            )
            spec = optimize_interval_mc(n, kappa, x_sample)
            lengths = length_given_x(x_grid, n, ALPHA, kappa, spec.A, spec.B)
            median_length = float(np.median(lengths))
            rows.append(
                {
                    "n": n,
                    "kappa": kappa,
                    "config": label,
                    "p_m": p_m,
                    "p_mc": spec.left_tail,
                    "abs_diff_p": abs(spec.left_tail - p_m),
                    "median_length_mc": median_length,
                    "median_diff_pct": 100.0 * (median_length - median_ref) / median_ref,
                }
            )
    return pd.DataFrame(rows)


def write_csv_outputs(
    performance_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    mc_df: pd.DataFrame,
    mc_sensitivity_df: pd.DataFrame,
) -> None:
    performance_df.to_csv(RESULTS_DIR / "performance_grid.csv", index=False)
    summary_df.to_csv(RESULTS_DIR / "summary_grid.csv", index=False)
    mc_df.to_csv(RESULTS_DIR / "monte_carlo_validation.csv", index=False)
    mc_sensitivity_df.to_csv(RESULTS_DIR / "mc_sensitivity.csv", index=False)


def format_table(value: float, digits: int = 4) -> str:
    rounded = round(float(value), digits)
    if rounded == 0:
        rounded = 0.0
    return f"{rounded:.{digits}f}"


def build_optima_table(summary_df: pd.DataFrame) -> str:
    lines = [
        r"\begin{tabular}{ccrrrrrr}",
        r"\toprule",
        r"$n$ & $\kappa$ & $p_E$ & $p_M$ & $p_{MC}$ & Mean/Median$_E$ & Med gain M (\%) & Med gain MC (\%) \\",
        r"\midrule",
    ]
    for n in N_VALUES:
        subset = summary_df[summary_df["n"] == n]
        for _, row in subset.iterrows():
            lines.append(
                " & ".join(
                    [
                        f"{int(row['n'])}",
                        format_table(row["kappa"], 2),
                        format_table(row["left_tail_expected"], 4),
                        format_table(row["left_tail_median"], 4),
                        format_table(row["left_tail_median_mc"], 4),
                        format_table(row["ratio_expected"], 3),
                        format_table(row["median_gain_vs_expected_pct"], 2),
                        format_table(row["median_gain_mc_vs_expected_pct"], 2),
                    ]
                )
                + r" \\"
            )
        if n != N_VALUES[-1]:
            lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    return "\n".join(lines) + "\n"

def build_length_table(performance_df: pd.DataFrame) -> str:
    filtered = performance_df[performance_df["method"].isin(["Equal-tail", "Expected-opt", "Median-opt"])]
    pivot = filtered.pivot_table(
        index=["n", "kappa"],
        columns="method",
        values=["mean_length", "median_length"],
    ).sort_index()
    lines = [
        r"\begin{tabular}{ccrrrrrr}",
        r"\toprule",
        r"$n$ & $\kappa$ & Mean ET & Mean EO & Mean MO & Median ET & Median EO & Median MO \\",
        r"\midrule",
    ]
    for (n, kappa), row in pivot.iterrows():
        lines.append(
            " & ".join(
                [
                    str(int(n)),
                    format_table(float(kappa), 2),
                    format_table(row[("mean_length", "Equal-tail")], 4),
                    format_table(row[("mean_length", "Expected-opt")], 4),
                    format_table(row[("mean_length", "Median-opt")], 4),
                    format_table(row[("median_length", "Equal-tail")], 4),
                    format_table(row[("median_length", "Expected-opt")], 4),
                    format_table(row[("median_length", "Median-opt")], 4),
                ]
            )
            + r" \\"
        )
        if kappa == KAPPA_VALUES[-1] and n != N_VALUES[-1]:
            lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    return "\n".join(lines) + "\n"


def build_mc_table(mc_df: pd.DataFrame) -> str:
    method_order = ["Equal-tail", "Expected-opt", "Median-opt", "Median-MC"]
    lines = [
        r"\begin{tabular}{cccrrr}",
        r"\toprule",
        r"$n$ & $\kappa$ & Method & $\widehat{\Pr}\{\text{cover}\}$ & Mean length & Median length \\",
        r"\midrule",
    ]
    for n, kappa in VALIDATION_CASES:
        subset = mc_df[(mc_df["n"] == n) & (mc_df["kappa"] == kappa)].set_index("method")
        for method in method_order:
            row = subset.loc[method]
            lines.append(
                " & ".join(
                    [
                        str(int(row["n"])),
                        format_table(row["kappa"], 2),
                        method,
                        format_table(row["coverage"], 4),
                        format_table(row["mean_length"], 4),
                        format_table(row["median_length"], 4),
                    ]
                )
                + r" \\"
            )
        if (n, kappa) != VALIDATION_CASES[-1]:
            lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    return "\n".join(lines) + "\n"


def build_mc_sensitivity_table(mc_sensitivity_df: pd.DataFrame) -> str:
    lines = [
        r"\begin{tabular}{cclrrrr}",
        r"\toprule",
        r"$n$ & $\kappa$ & Config & $p_M$ & $p_{MC}$ & $|p_{MC}-p_M|$ & $\Delta \widetilde{m}$ (\%) \\",
        r"\midrule",
    ]
    for n, kappa in VALIDATION_CASES:
        subset = mc_sensitivity_df[
            (mc_sensitivity_df["n"] == n) & (mc_sensitivity_df["kappa"] == kappa)
        ]
        for _, row in subset.iterrows():
            lines.append(
                " & ".join(
                    [
                        str(int(row["n"])),
                        format_table(row["kappa"], 2),
                        str(row["config"]),
                        format_table(row["p_m"], 4),
                        format_table(row["p_mc"], 4),
                        format_table(row["abs_diff_p"], 4),
                        format_table(row["median_diff_pct"], 2),
                    ]
                )
                + r" \\"
            )
        if (n, kappa) != VALIDATION_CASES[-1]:
            lines.append(r"\midrule")
    lines.extend([r"\bottomrule", r"\end{tabular}"])
    return "\n".join(lines) + "\n"


def write_latex_tables(
    performance_df: pd.DataFrame,
    summary_df: pd.DataFrame,
    mc_df: pd.DataFrame,
    mc_sensitivity_df: pd.DataFrame,
) -> None:
    (RESULTS_DIR / "table_optima.tex").write_text(build_optima_table(summary_df), encoding="utf-8")
    (RESULTS_DIR / "table_lengths.tex").write_text(build_length_table(performance_df), encoding="utf-8")
    (RESULTS_DIR / "table_monte_carlo.tex").write_text(build_mc_table(mc_df), encoding="utf-8")
    (RESULTS_DIR / "table_mc_sensitivity.tex").write_text(
        build_mc_sensitivity_table(mc_sensitivity_df),
        encoding="utf-8",
    )


def representative_plot(performance_df: pd.DataFrame) -> None:
    settings = [(2, 0.50), (8, 0.50)]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    for ax, (n, kappa) in zip(axes, settings):
        x_grid = evaluation_grid(n, grid_size=10_001)
        subset = performance_df[(performance_df["n"] == n) & (performance_df["kappa"] == kappa)]
        for method, color in [
            ("Equal-tail", "#7a6c5d"),
            ("Expected-opt", "#00798c"),
            ("Median-opt", "#d1495b"),
        ]:
            row = subset[subset["method"] == method].iloc[0]
            lengths = length_given_x(x_grid, n, ALPHA, kappa, row["A"], row["B"])
            ax.hist(
                lengths,
                bins=50,
                density=True,
                histtype="step",
                linewidth=1.6,
                color=color,
                label=method,
            )
            ax.axvline(np.median(lengths), color=color, linewidth=1.0, linestyle="--")

        ax.set_title(rf"$n={n},\ \kappa={kappa:.2f}$")
        ax.set_xlabel("Interval length")
        ax.set_ylabel("Density")

    axes[0].legend(frameon=False, loc="upper right")
    fig.savefig(FIGURES_DIR / "length_distribution.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "length_distribution.pdf", bbox_inches="tight")
    plt.close(fig)


def heatmap_plot(summary_df: pd.DataFrame) -> None:
    ratio_pivot = summary_df.pivot(index="n", columns="kappa", values="ratio_expected")
    gain_pivot = summary_df.pivot(index="n", columns="kappa", values="median_gain_vs_expected_pct")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), constrained_layout=True)
    cmap = plt.get_cmap("YlOrRd")

    for ax, pivot, title, fmt in [
        (axes[0], ratio_pivot, "Mean-to-median ratio of expected-opt length", "{:.2f}"),
        (axes[1], gain_pivot, "Median gain of median-opt over expected-opt (%)", "{:.2f}"),
    ]:
        im = ax.imshow(pivot.values, aspect="auto", cmap=cmap)
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels([f"{col:.2f}" for col in pivot.columns])
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels([str(int(idx)) for idx in pivot.index])
        ax.set_xlabel(r"Standardized shift $\kappa=\Delta/\sigma$")
        ax.set_ylabel(r"Sample size $n$")
        ax.set_title(title)
        for i in range(pivot.shape[0]):
            for j in range(pivot.shape[1]):
                value = pivot.iloc[i, j]
                ax.text(j, i, fmt.format(value), ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.savefig(FIGURES_DIR / "heatmaps.png", dpi=220, bbox_inches="tight")
    fig.savefig(FIGURES_DIR / "heatmaps.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ensure_directories()

    performance_rows: list[dict[str, float]] = []
    summary_rows: list[dict[str, float]] = []
    x_mc_by_n = {
        n: np.random.default_rng(MC_OPT_SEED + n).chisquare(df=n - 1, size=MC_OPT_REPS)
        for n in N_VALUES
    }

    for n in N_VALUES:
        for kappa in KAPPA_VALUES:
            case_performance, case_summary = summarize_case(n, kappa, x_mc_by_n[n])
            performance_rows.extend(case_performance)
            summary_rows.append(case_summary)

    performance_df = pd.DataFrame(performance_rows).sort_values(["n", "kappa", "method"]).reset_index(drop=True)
    summary_df = pd.DataFrame(summary_rows).sort_values(["n", "kappa"]).reset_index(drop=True)
    mc_df = monte_carlo_validation(performance_df).sort_values(["n", "method"]).reset_index(drop=True)
    mc_sensitivity_df = mc_sensitivity_study(summary_df).reset_index(drop=True)

    write_csv_outputs(performance_df, summary_df, mc_df, mc_sensitivity_df)
    write_latex_tables(performance_df, summary_df, mc_df, mc_sensitivity_df)
    representative_plot(performance_df)
    heatmap_plot(summary_df)

    print("Saved:")
    print(f"  {RESULTS_DIR / 'performance_grid.csv'}")
    print(f"  {RESULTS_DIR / 'summary_grid.csv'}")
    print(f"  {RESULTS_DIR / 'monte_carlo_validation.csv'}")
    print(f"  {RESULTS_DIR / 'mc_sensitivity.csv'}")
    print(f"  {FIGURES_DIR / 'length_distribution.pdf'}")
    print(f"  {FIGURES_DIR / 'heatmaps.pdf'}")


if __name__ == "__main__":
    main()
