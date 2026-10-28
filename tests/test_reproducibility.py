"""Scientific regression checks, run with ``python -m unittest discover -s tests``.

These checks cover input identity, conservation across shut-in, an independent
exact PME solution, and the published-data fit. They do not assert physical
validation against unprovided field saturation measurements.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from check_barenblatt import barenblatt_cell_averages, pme_evolve
from check_resolution import load_verified_history, primitive_quadrature_check
from fit_field_radius import load_site_data, fit_powerlaw_loglog_t0_scan
from solve_plume import solve


class ReproducibilityTests(unittest.TestCase):
    def test_input_files_match_recorded_hashes(self):
        manifest = json.loads((ROOT / "data/input_manifest.json").read_text())
        for name, expected in manifest["inputs"].items():
            with self.subTest(file=name):
                self.assertEqual(hashlib.sha256((ROOT / name).read_bytes()).hexdigest(),
                                 expected["sha256"])

    def test_production_history_is_current_model_and_conservative(self):
        history, metadata = load_verified_history(ROOT / "data/model", 6400)
        inventory = history["area"] @ history["u"]
        self.assertLess(np.max(abs(inventory - history["expected"])) / inventory[-1], 1e-10)
        self.assertAlmostEqual(float(history["u"][:, -1].max()), 0.736308848855286, places=9)
        self.assertEqual(metadata["shutin_state_jump"], 0.0)
        self.assertEqual(float(history["u"][-1].max()), 0.0)

    def test_new_solve_preserves_inventory_through_shutin(self):
        # A small fresh solve tests the boundary schedule and state handover.
        # The production accuracy is assessed separately on four fine grids.
        result = solve(48, rtol=1e-8, atol=1e-11, params={"ts": 0.05, "tend": 1.0})
        self.assertLess(result["metrics"]["relative_mass_error"], 1e-9)
        self.assertEqual(result["metrics"]["shutin_state_jump"], 0.0)
        self.assertLess(result["metrics"]["upper_bound_violation"], 1e-6)
        self.assertGreaterEqual(result["metrics"]["minimum"], -1e-9)

    def test_diffusion_primitive_matches_quadrature(self):
        self.assertLess(primitive_quadrature_check(), 1e-10)

    def test_pme_against_exact_barenblatt(self):
        inventory, diffusion, q, age, duration = 1.25e7, 2e4, 0.3, 50.0, 10.0
        edges = np.linspace(0, 20000, 161)
        initial = barenblatt_cell_averages(edges, inventory, diffusion, q, age)
        expected = barenblatt_cell_averages(edges, inventory, diffusion, q, age + duration)
        numerical, metrics = pme_evolve(edges, initial, diffusion, q, duration)
        area = np.pi * np.diff(edges**2)
        relative_error = np.sum(area * abs(numerical.y[:, -1] - expected)) / inventory
        self.assertLess(relative_error, 2e-3)
        self.assertLess(metrics["relative_inventory_error"], 1e-10)

    def test_field_fit_matches_independent_matlab_algorithm_port(self):
        # Values were independently checked by centered-covariance OLS on the
        # same 10000-point onset grid, not copied from the manuscript rounding.
        reference = {
            "Sleipner": (1998.4826827582759, 0.1590067948, 0.665087135, 8275),
            "Aquistore": (2015.2466753424656, 0.1279719624, 0.448838731, 0),
            "Weyburn--Midale": (2000.0001, 0.7586830348, 0.478108106, 0),
        }
        path = ROOT / "data/field/combined_radius_vs_year_all_sites.csv"
        for observations in load_site_data(path):
            with self.subTest(site=observations.site):
                fit = fit_powerlaw_loglog_t0_scan(
                    observations.survey_year, observations.radius_km,
                    observations.injection_start_year)
                onset, scale, exponent, index = reference[observations.site]
                self.assertAlmostEqual(fit["t0_year"], onset, places=8)
                self.assertAlmostEqual(fit["R0_km"], scale, places=8)
                self.assertAlmostEqual(fit["beta"], exponent, places=8)
                self.assertEqual(fit["t0_grid_index"], index)

    def test_original_survey_uncertainties_are_retained(self):
        path = ROOT / "data/field/combined_radius_vs_year_all_sites.csv"
        with path.open(newline="") as stream:
            raw = list(csv.DictReader(stream))
        for observations in load_site_data(path):
            rows = sorted((r for r in raw if r["Site"] == observations.site),
                          key=lambda r: float(r["SurveyYear"]))
            np.testing.assert_array_equal(observations.survey_error_year,
                                          [float(r["SurveyYearError_yr"]) for r in rows])


if __name__ == "__main__":
    unittest.main()
