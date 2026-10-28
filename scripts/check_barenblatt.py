#!/usr/bin/env python3
"""Audit the shut-in calculation against PME and Barenblatt references.

Run:
    python scripts/check_barenblatt.py
or:
    python scripts/check_barenblatt.py --data data/model

The source history and production figures are read only. The report is saved
in results/validation unless --output is supplied. No fitting of the full-model
profile is performed. Two comparisons answer different questions:
  1. Pure PME starts from the identical saved shut-in state on the same mesh.
  2. An exact finite-age Barenblatt initial profile tests the PME operator.
"""
import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
from scipy.sparse import diags


HERE = Path(__file__).resolve().parent


def barenblatt_parameters(inventory, D0, q, age):
    """Full-plane radial PME: u_t = D0/(2-q) * Laplacian(u**(2-q)).

    This normalization is derived from
        M = 2*pi*integral(u*r dr) = 4*pi*D0*age*A**(2-q).
    age is the similarity clock theta, not the physical elapsed time tau.
    """
    amplitude = (inventory / (4 * np.pi * D0 * age))**(1 / (2 - q))
    radius_squared = 4 * D0 * (2 - q) / (1 - q) * age * amplitude**(1 - q)
    radius = np.sqrt(radius_squared)
    exponent = 1 / (1 - q)
    radius95 = radius * np.sqrt(1 - 0.05**(1 / (exponent + 1)))
    return dict(
        age_year=float(age), amplitude=float(amplitude),
        support_radius_m=float(radius), inventory_radius95_m=float(radius95),
    )


def barenblatt_cell_averages(edges, inventory, D0, q, age):
    """Integrate the exact profile over each annulus before discretization."""
    pars = barenblatt_parameters(inventory, D0, q, age)
    radius = pars["support_radius_m"]
    exponent = 1 / (1 - q)
    remaining_mass_fraction = np.maximum(1 - edges**2 / radius**2, 0)**(exponent + 1)
    cell_mass = inventory * (remaining_mass_fraction[:-1] - remaining_mass_fraction[1:])
    return cell_mass / (np.pi * np.diff(edges**2))


def inventory_radius(edges, area, state, fraction=0.95):
    inventory = float(area @ state)
    return float(np.sqrt(np.interp(
        fraction * inventory, np.r_[0, np.cumsum(area * state)], edges**2
    )))


def pme_evolve(edges, initial, D0, q, duration):
    """Conservative annular PME with zero flux at both radial boundaries."""
    area = np.pi * np.diff(edges**2)
    dr = edges[1] - edges[0]
    coefficient = 2 * np.pi * edges[1:-1] / dr
    exponent = 2 - q

    def rhs(time, state):
        # Positive extension is used only for constitutive evaluation during
        # implicit iterations. The integrator's state is never clipped.
        positive = np.maximum(state, 0)
        primitive = D0 * positive**exponent / exponent
        flux = np.zeros(len(state) + 1)
        flux[1:-1] = -coefficient * np.diff(primitive)
        return -np.diff(flux) / area

    def jac(time, state):
        diffusion = D0 * np.maximum(state, 0)**(1 - q)
        left = coefficient * diffusion[:-1]
        right = -coefficient * diffusion[1:]
        diagonal = np.zeros(len(state))
        diagonal[:-1] -= left / area[:-1]
        diagonal[1:] += right / area[1:]
        return diags(
            [left / area[1:], diagonal, -right / area[:-1]],
            [-1, 0, 1], format="csc",
        )

    solution = solve_ivp(
        rhs, (0, duration), initial.copy(), method="BDF", jac=jac,
        rtol=1e-9, atol=1e-12,
        t_eval=np.r_[0, np.geomspace(1e-5, duration, 80)],
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    masses = area @ solution.y
    initial_inventory = float(area @ initial)
    return solution, dict(
        relative_inventory_error=float(np.max(abs(masses - initial_inventory))
                                       / initial_inventory),
        initial_state_jump=float(np.max(abs(solution.y[:, 0] - initial))),
        final_max_cell_average=float(np.max(solution.y[:, -1])),
        final_R95_m=inventory_radius(edges, area, solution.y[:, -1]),
        minimum=float(np.min(solution.y)),
        outermost_maximum=float(np.max(solution.y[-1])),
        rhs_evaluations=solution.nfev,
    )


def main(data_path, output_path):
    started = perf_counter()
    metadata = json.loads((data_path / "metrics.json").read_text())
    p = metadata["parameters"]
    with np.load(data_path / "history.npz", allow_pickle=False) as saved:
        data = {key: saved[key] for key in saved.files}
    ts, final_time = p["ts"], p["tend"]
    duration = final_time - ts
    shutin_index = int(np.argmin(abs(data["t"] - ts)))
    if not np.isclose(data["t"][shutin_index], ts):
        raise ValueError("The archive must include the exact shut-in time.")
    initial = data["u"][:, shutin_index].copy()
    final = data["u"][:, -1]
    edges, area = data["edges"], data["area"]
    inventory = float(area @ initial)
    D0, q, M = p["D0"], p["q"], p["mobility_ratio"]
    storage = p["phi"] * (1 - p["Sbr"]) * p["H"]
    peak = float(final.max())

    # Dimensionless ratios and scale estimates. L_M is the radius of a disk
    # containing this inventory at u=1; it is not a fitted plume radius.
    L_M = np.sqrt(inventory / np.pi)
    t_mass = inventory / (4 * np.pi * D0)
    D_pme_peak = D0 * peak**(1 - q)
    finite_thickness_factor = (1 - peak) / (1 + (M - 1) * peak)
    D_full_peak = D_pme_peak * finite_thickness_factor
    maximizer = minimize_scalar(
        lambda u: -D0 * u**(1-q) * (1-u) / (1+(M-1)*u),
        bounds=(0, 1), method="bounded",
    )
    D_full_max = -float(maximizer.fun)

    # A mass-and-peak matched Barenblatt has A(theta0)=1 at shut-in.
    # It still has a different radial profile from the inherited solution.
    virtual_age_peak = t_mass
    peak_matched_initial = barenblatt_parameters(inventory, D0, q, virtual_age_peak)
    peak_matched_final = barenblatt_parameters(inventory, D0, q,
                                              virtual_age_peak + duration)

    # A second legitimate convention matches the initial radial second moment.
    # For a Barenblatt profile <r^2>=R^2/(p+2), p=1/(1-q).
    # The annular formula below integrates r^2 exactly for cellwise constant u.
    mean_r_squared = float(
        np.sum(area * initial * (edges[:-1]**2 + edges[1:]**2) / 2) / inventory
    )
    exponent = 1 / (1 - q)
    moment_matched_R2 = (exponent + 2) * mean_r_squared
    moment_matched_A = inventory * (exponent + 1) / (np.pi * moment_matched_R2)
    virtual_age_moment = t_mass / moment_matched_A**(2 - q)
    moment_matched_initial = barenblatt_parameters(inventory, D0, q,
                                                  virtual_age_moment)
    moment_matched_final = barenblatt_parameters(inventory, D0, q,
                                                virtual_age_moment + duration)

    # The direct PME comparison uses the actual shut-in state, unlike either
    # virtual-age convention, and retains the same finite-radius inner boundary.
    direct_solution, direct_metrics = pme_evolve(edges, initial, D0, q, duration)
    direct_metrics["initial_inventory_m2"] = inventory
    direct_metrics["initial_R95_m"] = inventory_radius(edges, area, initial)

    # Independent coefficient benchmark on the full plane (inner radius=0).
    # Exact annular averages avoid an artificial initial discretization error.
    control_edges = np.linspace(0, p["rout"], len(initial) + 1)
    control_area = np.pi * np.diff(control_edges**2)
    control_initial = barenblatt_cell_averages(
        control_edges, inventory, D0, q, virtual_age_peak
    )
    control_exact_final = barenblatt_cell_averages(
        control_edges, inventory, D0, q, virtual_age_peak + duration
    )
    control_solution, control_metrics = pme_evolve(
        control_edges, control_initial, D0, q, duration
    )
    control_metrics.update(
        exact_initial_age_year=virtual_age_peak,
        exact_final_age_year=virtual_age_peak + duration,
        exact_final_peak=peak_matched_final["amplitude"],
        exact_final_R95_m=peak_matched_final["inventory_radius95_m"],
        relative_L1_error=float(
            np.sum(control_area * abs(control_solution.y[:, -1] - control_exact_final))
            / inventory
        ),
        max_cell_average_absolute_error=float(
            np.max(abs(control_solution.y[:, -1] - control_exact_final))
        ),
    )
    # Mean thickness over the disk that encloses 95% of inventory. Exclude the
    # well hole consistently with the exterior-domain mass definition.
    initial_R95 = float(data["R95"][shutin_index])
    final_R95 = float(data["R95"][-1])
    mean_thickness_initial = (
        p["H"] * 0.95 * inventory / (np.pi * (initial_R95**2 - edges[0]**2))
    )
    mean_thickness_final = (
        p["H"] * 0.95 * inventory / (np.pi * (final_R95**2 - edges[0]**2))
    )
    report = dict(
        source_history_sha256=hashlib.sha256(
            (data_path / "history.npz").read_bytes()
        ).hexdigest(),
        source_parameters=p, production_cells=len(initial),
        normalized_mobile_inventory_m2=inventory,
        mobile_reservoir_volume_m3=inventory * storage,
        storage_factor_m=storage,
        imposed_normalized_rate_m2_per_year=p["QV"] / storage,
        Gamma=p["QV"] / storage / (2 * np.pi * D0),
        shutin_duration_year=duration,
        inventory_disk_radius_m=L_M,
        inventory_diffusion_scale_year=L_M**2/D0,
        mass_Barenblatt_time_scale_year=t_mass,
        final_peak=peak,
        final_peak_thickness_m=p["H"] * peak,
        mean_thickness_within_R95_at_shutin_m=mean_thickness_initial,
        mean_thickness_within_R95_final_m=mean_thickness_final,
        mean_thickness_definition="H*0.95*M_s/[pi*(R95^2-rw^2)]",
        full_to_PME_diffusivity_at_final_peak=finite_thickness_factor,
        PME_diffusivity_at_final_peak_m2_per_year=D_pme_peak,
        full_diffusivity_at_final_peak_m2_per_year=D_full_peak,
        maximum_full_diffusivity_m2_per_year=D_full_max,
        thickness_at_maximum_full_diffusivity=float(maximizer.x),
        scale_LM_squared_over_maximum_full_D_year=L_M**2/D_full_max,
        scale_LM_squared_over_final_peak_D_year=L_M**2/D_full_peak,
        thickness_for_10percent_diffusivity_error=0.1/(M-0.1*(M-1)),
        thickness_for_5percent_diffusivity_error=0.05/(M-0.05*(M-1)),
        full_model_final_R95_m=float(data["R95"][-1]),
        full_model_inventory_error=metadata["relative_mass_error"],
        direct_PME_same_initial_state=direct_metrics,
        peak_matched_Barenblatt=dict(initial=peak_matched_initial,
                                    final=peak_matched_final),
        moment_matched_Barenblatt=dict(
            inherited_mean_r_squared_m2=mean_r_squared,
            initial=moment_matched_initial, final=moment_matched_final,
        ),
        unshifted_Barenblatt_at_tau_300=barenblatt_parameters(
            inventory, D0, q, duration
        ),
        exact_Barenblatt_numerical_control=control_metrics,
        interpretation=(
            "A Barenblatt virtual age selects one finite-width reference profile; "
            "matching mass and peak or second moment does not match the inherited "
            "profile. The direct PME run isolates the diffusivity law by using "
            "the identical initial state. L^2/D values are scale estimates, not "
            "predictions of the exact time at which the core ends."
        ),
        elapsed_wall_time_seconds=perf_counter() - started,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=HERE.parent / "data" / "model")
    parser.add_argument("--output", type=Path,
                        default=HERE.parent / "results" / "validation" / "barenblatt_comparison.json")
    args = parser.parse_args()
    main(args.data, args.output)
