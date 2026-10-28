#!/usr/bin/env python3
"""Verify spatial sensitivity of the finite-thickness plume calculation.

By default, compute 800-, 1600- and 3200-cell histories and compare them
with the supplied 6400-cell production history. Use ``--reuse`` to reuse
matching coarse histories, or ``--recompute-finest`` to solve all four grids.
This driver validates model parameters and array layout before accepting a
cache, so results from the earlier pure-PME example cannot be mixed in.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import hyp2f1

from solve_plume import PARAMETERS, save, solve


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
GRIDS = (800, 1600, 3200, 6400)
RTOL = 1e-9
ATOL = 1e-12


def sha256(path: Path) -> str:
    """Identify the exact numerical input bytes used for each grid."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_verified_history(folder: Path, cells: int) -> tuple[dict, dict]:
    """Reject incompatible parameters, tolerances, mesh or array layout.

    The production schema uses u[cell, time]. Historical pure-PME archives
    instead use u[time, cell] and different physical parameters. A matching
    filename alone is therefore insufficient evidence that a cache is valid.
    """
    metadata_path = folder / "metrics.json"
    history_path = folder / "history.npz"
    if not metadata_path.is_file() or not history_path.is_file():
        raise ValueError(f"Both history.npz and metrics.json are required in {folder}.")
    metadata = json.loads(metadata_path.read_text())
    parameters = metadata.get("parameters", {})
    mismatches = [key for key, value in PARAMETERS.items()
                  if key not in parameters or parameters[key] != value]
    if mismatches:
        raise ValueError("Different finite-thickness model parameters: "
                         + ", ".join(mismatches))
    if metadata.get("n") != cells:
        raise ValueError(f"Expected {cells} cells; metadata reports {metadata.get('n')}.")
    if metadata.get("rtol") != RTOL or metadata.get("atol") != ATOL:
        raise ValueError("The history uses different time-integration tolerances.")
    with np.load(history_path, allow_pickle=False) as archive:
        required = {"r", "edges", "area", "t", "u", "mass", "expected",
                    "R001", "a099", "R95", "umax"}
        if not required.issubset(archive.files):
            raise ValueError("Archive does not use the current bounded-model schema.")
        data = {name: archive[name] for name in archive.files}
    t, u, edges = data["t"], data["u"], data["edges"]
    if t.ndim != 1 or u.shape != (cells, len(t)):
        raise ValueError("Expected u[cell, time] on the declared grid.")
    if edges.shape != (cells + 1,) or data["r"].shape != (cells,):
        raise ValueError("Mesh arrays do not match the declared cell count.")
    if not (np.all(np.isfinite(u)) and np.all(np.isfinite(t))
            and np.all(np.diff(t) > 0)):
        raise ValueError("States/times must be finite and times strictly increasing.")
    p = PARAMETERS
    expected_edges = np.linspace(p["rw"], p["rout"], cells + 1)
    if not np.allclose(edges, expected_edges, rtol=1e-12, atol=1e-10):
        raise ValueError("Mesh geometry differs from the prescribed radial domain.")
    areas = np.pi * np.diff(edges**2)
    if data["area"].shape != (cells,) or not np.allclose(
            data["area"], areas, rtol=1e-12, atol=1e-8):
        raise ValueError("Stored cell areas disagree with annular geometry.")
    if not (np.isclose(t[0], 0, atol=1e-12, rtol=0)
            and np.isclose(t[-1], p["tend"], atol=1e-10, rtol=0)
            and np.any(np.isclose(t, p["ts"], atol=1e-10, rtol=0))):
        raise ValueError("History must include injection start, exact shut-in and final time.")
    inventory = areas @ u
    storage = p["phi"] * (1 - p["Sbr"]) * p["H"]
    expected_inventory = p["QV"] / storage * np.minimum(t, p["ts"])
    scale = expected_inventory[-1]
    if np.max(np.abs(inventory - expected_inventory)) / scale > 1e-7:
        raise ValueError("History fails the prescribed-inventory balance.")
    if not np.allclose(inventory, data["mass"], rtol=1e-10, atol=1e-8):
        raise ValueError("Stored inventory differs from the numerical states.")
    return data, metadata


def threshold_radius(radius: np.ndarray, state: np.ndarray, level: float) -> float:
    """Linear interpolation of the last crossing of a thickness threshold."""
    indices = np.flatnonzero(state >= level)
    if not len(indices):
        return 0.0
    index = int(indices[-1])
    if index == len(state) - 1:
        return float(radius[index])
    fraction = (state[index] - level) / (state[index] - state[index + 1])
    return float(radius[index] + fraction * (radius[index + 1] - radius[index]))


def inventory_radius(data: dict, state: np.ndarray) -> float:
    """Find R95 from the cell inventory, interpolating in radial area."""
    cumulative = np.r_[0.0, np.cumsum(data["area"] * state)]
    return float(np.sqrt(np.interp(0.95 * cumulative[-1], cumulative,
                                  data["edges"]**2)))


def diagnostics(data: dict, metadata: dict) -> dict:
    """Recompute reported quantities from the states rather than old tables."""
    index = int(np.argmin(abs(data["t"] - PARAMETERS["ts"])))
    shutin, final = data["u"][:, index], data["u"][:, -1]
    inventory = data["area"] @ data["u"]
    storage = PARAMETERS["phi"] * (1 - PARAMETERS["Sbr"]) * PARAMETERS["H"]
    expected = PARAMETERS["QV"] / storage * np.minimum(data["t"], PARAMETERS["ts"])
    return {
        "cells": metadata["n"],
        "dr_m": float(data["edges"][1] - data["edges"][0]),
        "R95_shutin_km": inventory_radius(data, shutin) / 1000,
        "R95_final_km": inventory_radius(data, final) / 1000,
        "R001_shutin_km": threshold_radius(data["r"], shutin, 0.001) / 1000,
        "R001_final_km": threshold_radius(data["r"], final, 0.001) / 1000,
        "a099_shutin_km": threshold_radius(data["r"], shutin, 0.99) / 1000,
        "umax_final": float(final.max()),
        "mass_error_relative": float(np.max(abs(inventory - expected)) / expected[-1]),
        "minimum_u": float(data["u"].min()),
        "maximum_bound_excess": float(max(0, data["u"].max() - 1)),
        # The merged archive stores one copy of t_s. The solver records the
        # stage-to-stage equality when both vectors are still available.
        "shutin_state_jump": metadata["shutin_state_jump"],
        "outermost_u": float(data["u"][-1].max()),
    }


def primitive_quadrature_check() -> float:
    """Compare the solver's analytic primitive with independent quadrature."""
    q, D0, mobility = (PARAMETERS[key] for key in ("q", "D0", "mobility_ratio"))
    exponent = 2 - q
    errors = []
    for u in (1e-4, 0.01, 0.1, 0.5, 0.9, 1.0):
        integral = quad(
            lambda v: D0 * v**(1-q) * (1-v) / (1+(mobility-1)*v),
            0, u, epsabs=1e-10, epsrel=1e-11,
        )[0]
        analytic = D0 * (
            u**exponent / exponent
            * hyp2f1(1, exponent, exponent+1, -(mobility-1)*u)
            - u**(exponent+1) / (exponent+1)
            * hyp2f1(1, exponent+1, exponent+2, -(mobility-1)*u)
        )
        errors.append(abs(analytic - integral) / abs(integral))
    return float(max(errors))


def run(output: Path, production: Path, reuse: bool = False,
        recompute_finest: bool = False) -> dict:
    """Compute/reuse compatible histories and write a four-grid report."""
    output.mkdir(parents=True, exist_ok=True)
    rows, sources = [], []
    for cells in GRIDS:
        use_production = cells == 6400 and not recompute_finest
        folder = production if use_production else output / f"N{cells}"
        accepted = False
        if use_production or reuse:
            try:
                data, metadata = load_verified_history(folder, cells)
                accepted = True
            except (ValueError, KeyError, TypeError, OSError) as error:
                if use_production:
                    raise ValueError(
                        f"Production history rejected: {error} "
                        "Generate the bounded-model production data, or use "
                        "--recompute-finest."
                    ) from error
                print(f"Cache N={cells} rejected; computing current model: {error}",
                      file=sys.stderr)
        if not accepted:
            print(f"Computing bounded-model history on {cells} cells.", flush=True)
            save(solve(cells, rtol=RTOL, atol=ATOL), folder)
            data, metadata = load_verified_history(folder, cells)
        rows.append(diagnostics(data, metadata))
        sources.append({
            "cells": cells,
            "origin": "supplied_production" if use_production else "refinement_run",
            "history_sha256": sha256(folder / "history.npz"),
            "metrics_sha256": sha256(folder / "metrics.json"),
        })

    names = list(rows[0])
    with (output / "resolution.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=names)
        writer.writeheader()
        writer.writerows(rows)
    compared = ("R95_shutin_km", "R95_final_km", "R001_shutin_km",
                "R001_final_km", "a099_shutin_km", "umax_final")
    refinement = {}
    for coarse, fine in zip(rows[:-1], rows[1:]):
        refinement[f"{coarse['cells']}_to_{fine['cells']}_percent"] = {
            key: 100 * abs(coarse[key] - fine[key]) / abs(fine[key])
            for key in compared
        }
    report = {
        "model": "finite_thickness_convection_diffusion",
        "parameters": PARAMETERS.copy(),
        "grids": rows,
        "refinement": refinement,
        "sources": sources,
        "primitive_quadrature_max_relative_error": primitive_quadrature_check(),
        "interpretation": (
            "These finite-grid changes quantify spatial sensitivity. Exact "
            "inventory conservation alone does not establish profile accuracy. "
            "R001 and a099 are thickness diagnostics, not exact free boundaries."
        ),
    }
    (output / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=REPOSITORY_ROOT / "results" / "refinement")
    parser.add_argument("--production-data", type=Path,
                        default=REPOSITORY_ROOT / "data" / "model")
    parser.add_argument("--reuse", action="store_true",
                        help="Reuse only matching coarse-grid histories.")
    parser.add_argument("--recompute-finest", action="store_true",
                        help="Also compute a fresh 6400-cell history in the output folder.")
    args = parser.parse_args()
    run(args.output, args.production_data, args.reuse, args.recompute_finest)
