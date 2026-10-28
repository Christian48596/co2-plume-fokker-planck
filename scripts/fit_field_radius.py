#!/usr/bin/env python3
"""Reproduce the Figure 4 field-radius fits using Python only.

The numerical fit is a port of the supplied ``fit_all_sites_radius_powerlaw.m``:
10,000 candidate onset times, ordinary least squares in natural-log space,
and selection of the candidate with the smallest log-residual sum of squares.
The survey uncertainties are displayed, but are not fit weights.

Default error bars use the supplied CSV unchanged. The MATLAB script replaced
each site's first survey-time uncertainty by 1/12 year. That display-only
override is available explicitly with ``--legacy-first-time-error`` and is
recorded in the output metadata; it never changes the fitted parameters.

Run from any working directory:
    python scripts/fit_field_radius.py
    python scripts/fit_field_radius.py --legacy-first-time-error

Outputs are Figure 4 (vector PDF and 300-dpi PNG) and machine-readable fits.
Only the log-log panel used in the manuscript is generated.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # Headless rendering on workstations and in GitHub Actions.
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SITES = ("Sleipner", "Aquistore", "Weyburn--Midale")
GRID_SIZE = 10_000
ONSET_OFFSET_YEARS = 1e-4


@dataclass(frozen=True)
class SiteData:
    """One site's observations, sorted by decimal survey year."""

    site: str
    survey_year: np.ndarray
    radius_km: np.ndarray
    survey_error_year: np.ndarray
    radius_error_km: np.ndarray
    injection_start_year: float
    excluded_csv_rows: tuple[int, ...]


def load_site_data(path: Path) -> list[SiteData]:
    """Read the original CSV without changing its rows or uncertainty values.

    Invalid numerical observations are excluded using the finite-value and
    positive-radius rule from MATLAB. CSV line numbers are retained in the
    report. Negative uncertainties and inconsistent site injection dates are
    rejected rather than silently reinterpreted.
    """
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        fields = set(reader.fieldnames or [])
        required = {"Site", "R_eq_km", "InjectionStartYear"}
        if not required <= fields or not ({"SurveyYear", "Year"} & fields):
            raise ValueError("CSV requires Site, R_eq_km, InjectionStartYear and SurveyYear or Year.")
        rows = list(reader)

    def number(row: dict[str, str], key: str, default: float = np.nan) -> float:
        value = row.get(key, "")
        return default if value in ("", None) else float(value)

    grouped = []
    time_key = "SurveyYear" if "SurveyYear" in fields else "Year"
    for site in SITES:
        source = [(line, row) for line, row in enumerate(rows, start=2) if row["Site"] == site]
        if not source:
            raise ValueError(f"Required manuscript site is absent: {site}")
        starts = np.array([number(row, "InjectionStartYear") for _, row in source])
        if not np.all(np.isfinite(starts)) or not np.allclose(starts, starts[0], rtol=0, atol=1e-10):
            raise ValueError(f"InjectionStartYear must be finite and constant within {site}.")
        values, excluded = [], []
        for line, row in source:
            # Missing columns imply absent error bars, as in MATLAB. A blank
            # cell in an existing column is missing data and fails the finite
            # observation filter instead of becoming a measured zero error.
            obs = [number(row, time_key), number(row, "R_eq_km"),
                   number(row, "SurveyYearError_yr") if "SurveyYearError_yr" in fields else 0.0,
                   number(row, "R_error_km") if "R_error_km" in fields else 0.0]
            if not np.all(np.isfinite(obs)) or obs[1] <= 0:
                excluded.append(line)
                continue
            if obs[2] < 0 or obs[3] < 0:
                raise ValueError(f"Negative uncertainty at CSV line {line}.")
            values.append(obs)
        if len(values) < 3:
            raise ValueError(f"At least three finite, positive-radius observations are required for {site}.")
        array = np.array(values, dtype=float)
        array = array[np.argsort(array[:, 0], kind="stable")]
        grouped.append(SiteData(site, array[:, 0], array[:, 1], array[:, 2],
                                array[:, 3], float(starts[0]), tuple(excluded)))
    return grouped


def fit_powerlaw_loglog_t0_scan(
    survey_year: np.ndarray,
    radius_km: np.ndarray,
    injection_start_year: float,
    grid_size: int = GRID_SIZE,
) -> dict[str, float | int | bool]:
    """Fit R = R0 * [(t - t0)/(1 year)]**beta by the original grid search.

    R0 has units km. Calendar times and elapsed time are expressed numerically
    in years. For each admissible t0, np.polyfit solves the same unweighted
    degree-one least-squares problem as MATLAB polyfit. No local optimizer,
    confidence interval, or inferred reservoir parameter is added.
    """
    t, radius = np.asarray(survey_year, dtype=float), np.asarray(radius_km, dtype=float)
    if t.ndim != 1 or radius.shape != t.shape or len(t) < 3:
        raise ValueError("Survey times and radii must be matching one-dimensional arrays with >=3 points.")
    if not np.all(np.isfinite(t)) or not np.all(np.isfinite(radius)) or np.any(radius <= 0):
        raise ValueError("Fit inputs must be finite, with positive radii.")
    if grid_size < 2 or not np.isfinite(injection_start_year):
        raise ValueError("The onset grid and injection start must be valid.")
    order = np.argsort(t, kind="stable")
    t, radius = t[order], radius[order]
    lower = float(injection_start_year + ONSET_OFFSET_YEARS)
    upper = float(t.min() - ONSET_OFFSET_YEARS)
    if upper <= lower:
        raise ValueError("Injection start must precede the first survey by more than 2e-4 years.")
    grid = np.linspace(lower, upper, grid_size)
    y = np.log(radius)
    coefficients = np.empty((grid_size, 2))
    residual_sums = np.empty(grid_size)
    for index, onset in enumerate(grid):
        x = np.log(t - onset)
        coefficients[index] = np.polyfit(x, y, deg=1)
        residual_sums[index] = np.sum((y - np.polyval(coefficients[index], x))**2)
    best = int(np.argmin(residual_sums))  # First minimum, as in MATLAB min().
    slope, intercept = coefficients[best]
    total = float(np.sum((y - np.mean(y))**2))
    if total == 0:
        raise ValueError("Log-radius variance is zero; R2_log is undefined.")
    return {
        "InjectionStartYear": float(injection_start_year),
        "t0_year": float(grid[best]),
        "R0_km": float(np.exp(intercept)),
        "beta": float(slope),
        "R2_log": float(1 - residual_sums[best] / total),
        "N_points": len(t),
        "SSE_log": float(residual_sums[best]),
        "t0_grid_index": best,
        "t0_grid_size": grid_size,
        "t0_grid_lower_year": lower,
        "t0_grid_upper_year": upper,
        "t0_at_grid_boundary": best in (0, grid_size - 1),
    }


def plot_figure4(
    observations: list[SiteData],
    fits: list[dict],
    output_directory: Path,
    legacy_first_time_error: bool = False,
) -> list[dict]:
    """Draw the measured uncertainties and fitted curves on logarithmic axes.

    The first time-error override is applied to a copy used for drawing only.
    Horizontal error bars extending beyond the plot's lower limit are clipped
    to that displayed limit, as in the original plot. Source uncertainties
    remain stored in the accompanying JSON, including the unclipped endpoint.
    """
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 10,
        "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
        "axes.labelsize": 11, "legend.fontsize": 9.5,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.linewidth": 0.8, "lines.linewidth": 1.7,
    })
    styles = (("red", "o", "-"), ("blue", "s", "--"), ("magenta", "D", "-."))
    xmin, xmax = 0.1, 10.0
    curve_x = np.geomspace(xmin, xmax, 400)
    fig, ax = plt.subplots(figsize=(7.1, 4.2), layout="constrained")
    ax.set_xscale("log")
    ax.set_yscale("log")
    display_records, handles = [], []
    for obs, fit, (color, marker, linestyle) in zip(observations, fits, styles):
        x = obs.survey_year - fit["t0_year"]
        dx = obs.survey_error_year.copy()
        if legacy_first_time_error:
            dx[0] = 1 / 12
        ax.errorbar(x, obs.radius_km, yerr=obs.radius_error_km,
                    color=color, marker=marker, markerfacecolor="none",
                    linestyle="none", markersize=6.5, capsize=3.5,
                    elinewidth=1.1, markeredgewidth=1.1, zorder=3)
        for i, (xp, yp, dt) in enumerate(zip(x, obs.radius_km, dx)):
            left, right = max(float(xp-dt), xmin), float(xp+dt)
            # A zero time uncertainty has no horizontal bar or endpoint caps.
            # Drawing two coincident caps would add a spurious vertical mark.
            if dt > 0:
                ax.hlines(yp, left, right, color=color, linewidth=1.1, zorder=2)
                ax.vlines([left, right], 0.94*yp, 1.06*yp, color=color, linewidth=1.1, zorder=2)
            display_records.append({
                "Site": obs.site, "SurveyYear": float(obs.survey_year[i]),
                "source_time_error_year": float(obs.survey_error_year[i]),
                "display_time_error_year": float(dt),
                "elapsed_time_year": float(xp),
                "unclipped_left_error_endpoint_year": float(xp-dt),
                "left_endpoint_clipped_to_axis": bool(xp-dt < xmin),
            })
        ax.plot(curve_x, fit["R0_km"]*curve_x**fit["beta"], color=color,
                linestyle=linestyle, linewidth=1.8, zorder=1)
        label = rf"{obs.site.replace('--', '-')} ($\beta={fit['beta']:.3f}$)"
        handles.append(Line2D([], [], color=color, marker=marker, markerfacecolor="none",
                              linestyle=linestyle, markersize=6.5, linewidth=1.8, label=label))
    ax.set(xlim=(xmin, xmax), ylim=(0.05, 5),
           xlabel=r"Elapsed time, $t-t_0$ (yr)",
           ylabel=r"Equivalent radius, $R_{\mathrm{eq}}$ (km)")
    ax.grid(which="major", alpha=0.22, linewidth=0.6)
    ax.grid(which="minor", alpha=0.08, linewidth=0.4)
    ax.legend(handles=handles, loc="upper left", frameon=False)
    output_directory.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_directory / "all_sites_radius_powerlaw_loglog.pdf",
                metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(output_directory / "all_sites_radius_powerlaw_loglog.png", dpi=300)
    plt.close(fig)
    return display_records


def main(data: Path, figures: Path, results: Path, legacy_first_time_error: bool = False) -> dict:
    """Fit every manuscript site, draw Figure 4, and write auditable outputs."""
    data, figures, results = Path(data), Path(figures), Path(results)
    observations = load_site_data(data)
    fits = [dict(Site=obs.site, **fit_powerlaw_loglog_t0_scan(
        obs.survey_year, obs.radius_km, obs.injection_start_year)) for obs in observations]
    display_records = plot_figure4(observations, fits, figures, legacy_first_time_error)
    results.mkdir(parents=True, exist_ok=True)
    fields = ("Site", "InjectionStartYear", "t0_year", "R0_km", "beta", "R2_log", "N_points")
    with (results / "radius_fit_parameters.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(fits)
    report = {
        "source_csv": data.name,
        "source_sha256": hashlib.sha256(data.read_bytes()).hexdigest(),
        "model": "R_eq = R0 * [(t - t0)/(1 year)]**beta",
        "units": {"R0": "km", "t": "decimal calendar year", "t0": "decimal calendar year"},
        "fit_method": "10000-point onset scan; unweighted natural-log ordinary least squares",
        "onset_constraints": "injection start + 1e-4 <= t0 <= first survey - 1e-4 years",
        "uncertainty_role": "Displayed observational error bars; not regression weights or fitted-onset uncertainty.",
        "time_provenance": "Survey dates and errors are retained from the supplied CSV. See repository provenance for the Aquistore date-to-volume mapping query.",
        "legacy_first_time_error_override": legacy_first_time_error,
        "legacy_override_year": 1/12 if legacy_first_time_error else None,
        "interpretation": "Injection-period footprint descriptors; no constitutive q or shut-in rate inferred.",
        "excluded_csv_rows": {obs.site: list(obs.excluded_csv_rows) for obs in observations},
        "fits": fits,
        "plotted_time_uncertainties": display_records,
    }
    (results / "radius_fit_parameters.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"fits": fits, "legacy_first_time_error_override": legacy_first_time_error}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path,
                        default=ROOT / "data/field/combined_radius_vs_year_all_sites.csv")
    parser.add_argument("--figures", type=Path, default=ROOT / "figures")
    parser.add_argument("--results", type=Path, default=ROOT / "results/field")
    parser.add_argument("--legacy-first-time-error", action="store_true",
                        help="Reproduce the MATLAB display-only +/-1/12 year override for each first survey.")
    args = parser.parse_args()
    main(args.data, args.figures, args.results, args.legacy_first_time_error)
