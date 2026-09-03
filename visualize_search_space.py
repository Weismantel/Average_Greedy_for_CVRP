#!/usr/bin/env python3
"""Plot fixed-r slices of the three-variable search space from Lemma 4.1.

The visualization is numerical and is not part of the rigorous interval
certificate.  It uses the same reduced cost functions and the same convention
for inactive quotient branches as ``interval_check.cpp``.
"""

from __future__ import annotations

import argparse
import math
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

# Some reproducibility environments have a read-only home directory.  Keep
# Matplotlib's font cache in a writable temporary directory in that case.
os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "cvrp-matplotlib-cache")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import QuadMesh


R_MIN = 33.0 / 50.0
R_MAX = 1.0
TRIANGLE_SIZE = 3.0 / 25.0
TARGET = 79.0 / 25.0


@dataclass(frozen=True)
class SliceValues:
    """Values of both cost bounds on one fixed-r triangular grid."""

    r: float
    m: np.ndarray
    lambda_: np.ndarray
    eta_cost: np.ndarray
    one_third_cost: np.ndarray
    objective: np.ndarray
    maximum_m: float
    maximum_lambda: float
    maximum_eta_cost: float
    maximum_one_third_cost: float
    maximum_objective: float


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate colored fixed-r slices of the Lemma 4.1 search space."
        )
    )
    parser.add_argument(
        "--r-values",
        type=float,
        nargs="+",
        default=[0.66, 0.75, 0.82, 0.829526, 0.86, 1.0],
        metavar="R",
        help="fixed r values (default: %(default)s)",
    )
    parser.add_argument(
        "--grid-size",
        type=int,
        default=801,
        help="number of grid coordinates on each triangle axis (default: 801)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "figures",
        help="directory for generated pictures",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=220,
        help="PNG resolution (default: 220)",
    )
    parser.add_argument(
        "--combined-only",
        action="store_true",
        help="generate only the combined multi-panel picture",
    )
    arguments = parser.parse_args()

    if arguments.grid_size < 25:
        parser.error("--grid-size must be at least 25")
    if arguments.dpi < 72:
        parser.error("--dpi must be at least 72")
    invalid_r = [r for r in arguments.r_values if not R_MIN <= r <= R_MAX]
    if invalid_r:
        parser.error(
            f"every r value must lie in [{R_MIN:.2f}, {R_MAX:.2f}]; "
            f"invalid values: {invalid_r}"
        )
    return arguments


def internal_root(numerator: np.ndarray, denominator: np.ndarray) -> np.ndarray:
    """Return roots in (0,1), replacing all other values by the endpoint 1."""

    root = np.full_like(numerator, np.inf)
    np.divide(
        numerator,
        denominator,
        out=root,
        where=np.abs(denominator) > np.finfo(float).tiny,
    )
    return np.where(np.isfinite(root) & (root > 0.0) & (root < 1.0), root, 1.0)


def analytic_quotient_integral(
    numerator_slope: np.ndarray,
    denominator_constant: np.ndarray,
    denominator_slope: np.ndarray,
    left: np.ndarray,
    right: np.ndarray,
) -> np.ndarray:
    """Integrate (1-p*t)/(a-b*t) on vectorized intervals."""

    width = right - left
    coefficient = denominator_slope - denominator_constant * numerator_slope
    scale = np.maximum.reduce(
        [
            np.ones_like(coefficient),
            np.abs(denominator_slope),
            np.abs(denominator_constant * numerator_slope),
        ]
    )
    proportional = np.abs(coefficient) <= 64.0 * np.finfo(float).eps * scale

    value = numerator_slope / denominator_slope * width
    regular = ~proportional
    denominator_left = denominator_constant - denominator_slope * left
    denominator_right = denominator_constant - denominator_slope * right
    regular &= (denominator_left > 0.0) & (denominator_right > 0.0)
    value[regular] = (
        numerator_slope[regular] / denominator_slope[regular] * width[regular]
        + coefficient[regular]
        / denominator_slope[regular] ** 2
        * np.log(denominator_left[regular] / denominator_right[regular])
    )
    # A regular segment can fail the endpoint test only through roundoff at a
    # breakpoint. The constant branch is the conservative visual fallback.
    value[(~proportional) & (~regular)] = width[(~proportional) & (~regular)]
    return np.clip(value, 0.0, width)


def analytic_envelope_integral(
    numerator_slope: np.ndarray,
    branches: Sequence[tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    """Integrate min{1,q_1,...} by splitting at every linear crossing."""

    breakpoints = [
        np.zeros_like(numerator_slope),
        np.ones_like(numerator_slope),
        internal_root(np.ones_like(numerator_slope), numerator_slope),
    ]
    for denominator_constant, denominator_slope in branches:
        breakpoints.append(
            internal_root(denominator_constant, denominator_slope)
        )
        breakpoints.append(
            internal_root(
                1.0 - denominator_constant,
                numerator_slope - denominator_slope,
            )
        )
    for first in range(len(branches)):
        for second in range(first + 1, len(branches)):
            breakpoints.append(
                internal_root(
                    branches[first][0] - branches[second][0],
                    branches[first][1] - branches[second][1],
                )
            )

    ordered = np.sort(np.stack(breakpoints), axis=0)
    integral = np.zeros_like(numerator_slope)
    for index in range(ordered.shape[0] - 1):
        left = ordered[index]
        right = ordered[index + 1]
        middle = (left + right) / 2.0
        numerator = 1.0 - numerator_slope * middle
        active_value = np.ones_like(numerator_slope)
        active_branch = np.full(numerator_slope.shape, -1, dtype=np.int8)
        for branch_index, (denominator_constant, denominator_slope) in enumerate(
            branches
        ):
            denominator = denominator_constant - denominator_slope * middle
            quotient = np.full_like(numerator, np.inf)
            np.divide(
                numerator,
                denominator,
                out=quotient,
                where=(numerator >= 0.0) & (denominator > 0.0),
            )
            better = quotient < active_value
            active_value[better] = quotient[better]
            active_branch[better] = branch_index

        segment = right - left
        for branch_index, (denominator_constant, denominator_slope) in enumerate(
            branches
        ):
            active = active_branch == branch_index
            if not np.any(active):
                continue
            branch_integral = analytic_quotient_integral(
                numerator_slope,
                denominator_constant,
                denominator_slope,
                left,
                right,
            )
            segment[active] = branch_integral[active]
        integral += segment
    return integral


def evaluate_slice(r: float, grid_size: int) -> SliceValues:
    """Evaluate one fixed-r slice using analytic antiderivatives."""

    coordinates = np.linspace(0.0, TRIANGLE_SIZE, grid_size)
    m_grid, lambda_grid = np.meshgrid(coordinates, coordinates)
    feasible = m_grid + lambda_grid <= np.nextafter(TRIANGLE_SIZE, math.inf)
    m = m_grid[feasible]
    lambda_ = lambda_grid[feasible]

    numerator_slope = 16.0 / 25.0 + 2.0 * m + 3.0 * lambda_
    eta_denominator_constant = np.full_like(m, 2.0 * r)
    eta_denominator_slope = (
        2.0 * r + 16.0 / 25.0 + 2.0 * m + 4.0 * lambda_
    )

    sigma_one_third = 3.0 * r / 2.0 + 4.0 / 25.0 + m
    first_denominator_constant = (
        3.0 * r / 2.0
        - 1.0 / 50.0
        + 2.0 * m
        + 3.0 * lambda_ / 2.0
    )
    first_denominator_slope = (
        3.0 * r / 2.0
        + 16.0 / 25.0
        + 5.0 * m / 2.0
        + 3.0 * lambda_
    )
    second_denominator_constant = sigma_one_third
    second_denominator_slope = (
        3.0 * r / 2.0 + 4.0 / 5.0 + 3.0 * m + 4.0 * lambda_
    )

    eta_integral = analytic_envelope_integral(
        numerator_slope,
        [(eta_denominator_constant, eta_denominator_slope)],
    )
    one_third_integral = analytic_envelope_integral(
        numerator_slope,
        [
            (first_denominator_constant, first_denominator_slope),
            (second_denominator_constant, second_denominator_slope),
        ],
    )

    eta_cost_values = 5.0 / 2.0 - r + 2.0 * r * eta_integral
    one_third_cost_values = (
        3.0 - 3.0 * r / 2.0 + sigma_one_third * one_third_integral
    )
    objective_values = np.minimum(eta_cost_values, one_third_cost_values)

    eta_cost = np.full_like(m_grid, np.nan)
    one_third_cost = np.full_like(m_grid, np.nan)
    objective = np.full_like(m_grid, np.nan)
    eta_cost[feasible] = eta_cost_values
    one_third_cost[feasible] = one_third_cost_values
    objective[feasible] = objective_values

    maximum_index = int(np.argmax(objective_values))
    return SliceValues(
        r=r,
        m=m_grid,
        lambda_=lambda_grid,
        eta_cost=eta_cost,
        one_third_cost=one_third_cost,
        objective=objective,
        maximum_m=float(m[maximum_index]),
        maximum_lambda=float(lambda_[maximum_index]),
        maximum_eta_cost=float(eta_cost_values[maximum_index]),
        maximum_one_third_cost=float(one_third_cost_values[maximum_index]),
        maximum_objective=float(objective_values[maximum_index]),
    )


def color_limits(slices: Sequence[SliceValues]) -> tuple[float, float]:
    finite_values = np.concatenate(
        [values.objective[np.isfinite(values.objective)] for values in slices]
    )
    lower = float(np.min(finite_values))
    upper = max(TARGET, float(np.max(finite_values)))
    padding = 0.015 * (upper - lower)
    return lower - padding, upper


def label_colorbar(colorbar: matplotlib.colorbar.Colorbar) -> None:
    ticks = [
        tick
        for tick in colorbar.get_ticks()
        if colorbar.vmin <= tick <= colorbar.vmax
    ]
    if colorbar.vmin <= TARGET <= colorbar.vmax and all(
        not math.isclose(tick, TARGET, abs_tol=1e-10) for tick in ticks
    ):
        colorbar.set_ticks(sorted([*ticks, TARGET]))
    colorbar.set_label(
        r"approximation bound $\min\{\widehat C_\eta,\widehat C_{1/3}\}$"
    )


def draw_slice(
    axis: plt.Axes,
    values: SliceValues,
    normalization: colors.Normalize,
    contour_levels: np.ndarray,
) -> QuadMesh:
    image = axis.pcolormesh(
        values.m,
        values.lambda_,
        values.objective,
        shading="auto",
        cmap="magma",
        norm=normalization,
        rasterized=True,
    )
    finite_objective = values.objective[np.isfinite(values.objective)]
    visible_levels = contour_levels[
        (contour_levels > np.min(finite_objective))
        & (contour_levels < np.max(finite_objective))
    ]
    if visible_levels.size:
        axis.contour(
            values.m,
            values.lambda_,
            values.objective,
            levels=visible_levels,
            colors="black",
            linewidths=0.35,
            alpha=0.45,
        )

    difference = values.eta_cost - values.one_third_cost
    finite_difference = difference[np.isfinite(difference)]
    if np.min(finite_difference) <= 0.0 <= np.max(finite_difference):
        axis.contour(
            values.m,
            values.lambda_,
            difference,
            levels=[0.0],
            colors="white",
            linewidths=1.15,
            linestyles="--",
        )

    boundary_m = np.linspace(0.0, TRIANGLE_SIZE, 300)
    axis.plot(
        boundary_m,
        TRIANGLE_SIZE - boundary_m,
        color="black",
        linewidth=1.0,
    )
    axis.scatter(
        [values.maximum_m],
        [values.maximum_lambda],
        marker="*",
        s=68,
        facecolor="#38d6ff",
        edgecolor="black",
        linewidth=0.7,
        zorder=5,
    )
    axis.set_xlim(0.0, TRIANGLE_SIZE)
    axis.set_ylim(0.0, TRIANGLE_SIZE)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(r"$m$")
    axis.set_ylabel(r"$\lambda$")
    axis.set_title(
        rf"$r={values.r:.6f}$; max $\approx {values.maximum_objective:.6f}$",
        fontsize=10,
    )
    axis.tick_params(labelsize=8)
    return image


def save_combined_picture(
    slices: Sequence[SliceValues],
    output_path: Path,
    dpi: int,
) -> None:
    columns = min(3, len(slices))
    rows = math.ceil(len(slices) / columns)
    lower, upper = color_limits(slices)
    normalization = colors.Normalize(vmin=lower, vmax=upper)
    contour_levels = np.linspace(lower, upper, 9)[1:-1]

    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(4.25 * columns, 3.9 * rows),
        constrained_layout=True,
        squeeze=False,
    )
    image = None
    for axis, values in zip(axes.flat, slices, strict=False):
        image = draw_slice(axis, values, normalization, contour_levels)
    for axis in axes.flat[len(slices) :]:
        axis.remove()

    figure.suptitle(
        r"Lemma 4.1: $\min\{\widehat C_\eta,\widehat C_{1/3}\}$ "
        r"on fixed-$r$ slices",
        fontsize=14,
    )
    assert image is not None
    colorbar = figure.colorbar(image, ax=list(axes.flat[: len(slices)]), shrink=0.88)
    label_colorbar(colorbar)
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def filename_for_r(r: float) -> str:
    return f"search_space_r_{r:.6f}".replace(".", "p") + ".png"


def save_individual_picture(
    values: SliceValues,
    output_path: Path,
    dpi: int,
    lower: float,
    upper: float,
) -> None:
    figure, axis = plt.subplots(figsize=(5.3, 4.8), constrained_layout=True)
    normalization = colors.Normalize(vmin=lower, vmax=upper)
    contour_levels = np.linspace(lower, upper, 9)[1:-1]
    image = draw_slice(axis, values, normalization, contour_levels)
    colorbar = figure.colorbar(image, ax=axis)
    label_colorbar(colorbar)
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    arguments = parse_arguments()
    arguments.output_dir.mkdir(parents=True, exist_ok=True)

    slices = [
        evaluate_slice(r, arguments.grid_size)
        for r in arguments.r_values
    ]
    combined_path = arguments.output_dir / "lemma_4_1_search_space.png"
    save_combined_picture(slices, combined_path, arguments.dpi)

    if not arguments.combined_only:
        lower, upper = color_limits(slices)
        for values in slices:
            save_individual_picture(
                values,
                arguments.output_dir / filename_for_r(values.r),
                arguments.dpi,
                lower,
                upper,
            )

    for values in slices:
        print(
            f"r={values.r:.9f}: grid max={values.maximum_objective:.9f} "
            f"at m={values.maximum_m:.9f}, lambda={values.maximum_lambda:.9f} "
            f"(C_eta={values.maximum_eta_cost:.9f}, "
            f"C_1/3={values.maximum_one_third_cost:.9f})"
        )
    print(f"Wrote {combined_path}")


if __name__ == "__main__":
    main()
