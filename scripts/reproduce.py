#!/usr/bin/env python3
"""Reproduce manuscript Figures 2--4, with optional fresh numerical verification.

Default: redraw all figures from the supplied production history and field CSV.
--recompute: solve the 6400-cell model into results/recomputed/model first.
--validate: run the tests, PME/Barenblatt comparison and four-grid study.
All subprocesses use the current Python interpreter and the repository's paths.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter

import matplotlib
import numpy
import scipy


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"


def run_script(name: str, *arguments: str) -> None:
    """Stop on an unsuccessful step; never continue with silently stale outputs."""
    command = [sys.executable, str(SCRIPTS / name), *map(str, arguments)]
    print(f"Running {name}", flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main(recompute: bool = False, validate: bool = False,
         legacy_first_time_error: bool = False) -> None:
    started = perf_counter()
    data = ROOT / "data/model"
    if recompute:
        data = ROOT / "results/recomputed/model"
        run_script("solve_plume.py", "--cells", "6400", "--output", str(data))

    run_script("plot_model_figures.py", "--data", str(data))
    field_flags = ["--legacy-first-time-error"] if legacy_first_time_error else []
    run_script("fit_field_radius.py", *field_flags)
    if validate:
        subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"],
                       cwd=ROOT, check=True)
        run_script("check_barenblatt.py", "--data", str(data))
        run_script("check_resolution.py", "--production-data", str(data), "--reuse")

    report = {
        "status": "completed",
        "recomputed_model": recompute,
        "ran_validation": validate,
        "legacy_first_time_error_override": legacy_first_time_error,
        "model_data_directory": str(data.relative_to(ROOT)),
        "software": {"python": platform.python_version(), "numpy": numpy.__version__,
                     "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
        "elapsed_seconds": perf_counter() - started,
        "outputs": ["figures/continuous_plume_profiles.pdf",
                    "figures/continuous_plume_profiles.png",
                    "figures/continuous_plume_history.pdf",
                    "figures/continuous_plume_history.png",
                    "figures/all_sites_radius_powerlaw_loglog.pdf",
                    "figures/all_sites_radius_powerlaw_loglog.png",
                    "results/field/radius_fit_parameters.csv"],
    }
    output = ROOT / "results/reproduction.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print("Reproduction completed. Figures are in figures/.", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--recompute", action="store_true",
                        help="Compute a fresh 6400-cell history without overwriting the supplied data.")
    parser.add_argument("--validate", action="store_true",
                        help="Also run scientific tests, the Barenblatt benchmark and spatial refinement.")
    parser.add_argument("--legacy-first-time-error", action="store_true",
                        help="Use the MATLAB display-only one-month first-survey uncertainty override.")
    args = parser.parse_args()
    main(args.recompute, args.validate, args.legacy_first_time_error)
