#!/usr/bin/env python3
"""Compute the bounded radial CO2 plume model used in Figures 2 and 3.

Run from the package directory:
    python scripts/solve_plume.py --cells 6400

The calculation starts from an empty aquifer, injects at a prescribed rate
until t_s, and continues the identical numerical state after shut-in.
Units are metres, years, and reservoir cubic metres. This is an illustrative
parameter set, independent of the field-image radius fits.
"""
import argparse
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.sparse import diags
from scipy.special import hyp2f1


# PHYSICAL PARAMETERS
# H [m]: aquifer thickness; phi [-]: porosity; Sbr [-]: residual brine.
# q [-]: constitutive exponent, 0 <= q < 1.
# D0 [m^2/yr]: gravity-spreading scale.
# mobility_ratio [-]: endpoint gas/brine mobility ratio, denoted M in comments.
# QV [m^3/yr]: prescribed reservoir-volume injection rate.
# ts [yr]: shut-in time; tend [yr]: final time, both from injection start.
# rw [m]: effective well radius; rout [m]: outer zero-flux boundary.
PARAMETERS = dict(
    H=10.0, phi=0.30, Sbr=0.20, q=0.30, D0=2.0e4,
    mobility_ratio=5.0, QV=1.0e7, ts=3.0, rw=10.0,
    rout=20000.0, tend=303.0,
)
PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def validate_parameters(p, n):
    """Catch parameter choices outside the regime implemented by this solver."""
    if n < 2 or not 0 <= p["q"] < 1:
        raise ValueError("Use at least 2 cells and 0 <= q < 1.")
    if not (0 < p["phi"] <= 1 and 0 <= p["Sbr"] < 1):
        raise ValueError("Require 0 < phi <= 1 and 0 <= Sbr < 1.")
    if any(p[k] <= 0 for k in ("H", "D0", "mobility_ratio", "QV", "rw")):
        raise ValueError("H, D0, mobility_ratio, QV and rw must be positive.")
    if not (0 < p["ts"] < p["tend"] and p["rw"] < p["rout"]):
        raise ValueError("Require 0 < ts < tend and rw < rout.")


def solve(n=6400, rtol=1e-9, atol=1e-12, params=None):
    """Return numerical states, diagnostic radii and conservation checks.

    u = h/H is the cell-averaged mobile-layer fraction. Define the outward,
    azimuthally integrated transport F, with units m^2/yr:

        F = rate*f(u) - 2*pi*r*D(u)*du/dr,
        f(u) = M*u/[1+(M-1)*u],
        D(u) = D0*u**(1-q)*(1-u)/[1+(M-1)*u].

    The balance is du/dt = -(1/(2*pi*r))*dF/dr. During injection,
    rate = QV/[phi*(1-Sbr)*H]; after shut-in rate = 0.
    The inner boundary supplies F(rw)=rate and the outer boundary has F=0.
    This is the bounded convection-diffusion model, including the complete
    finite-thickness factor and column fractional flow.
    """
    p = PARAMETERS.copy()
    if params:
        p.update(params)
    validate_parameters(p, n)

    # Uniform spacing in radius gives annuli of different horizontal areas.
    # A shared flux at each interior face makes the scheme conservative.
    edges = np.linspace(p["rw"], p["rout"], n + 1)
    r = (edges[:-1] + edges[1:]) / 2
    area = np.pi * np.diff(edges**2)
    dr = edges[1] - edges[0]
    storage = p["phi"] * (1 - p["Sbr"]) * p["H"]  # [m]
    qu = p["QV"] / storage                            # [m^2/yr]
    m = 2 - p["q"]
    mu = p["mobility_ratio"]

    # Psi'(u)=D(u). Its analytic primitive uses the Gauss hypergeometric
    # function. Finite differences of Psi represent D(u)*du/dr.
    def primitive(state):
        # Extend constitutive functions outside their physical interval for
        # intermediate implicit-solver iterations. This does not overwrite
        # the evolving state: no profile cap or mass renormalization is used.
        s = np.clip(state, 0.0, 1.0)
        return p["D0"] * (
            s**m / m * hyp2f1(1, m, m + 1, -(mu - 1) * s)
            - s**(m + 1) / (m + 1)
            * hyp2f1(1, m + 1, m + 2, -(mu - 1) * s)
        )

    def diffusivity(state):
        s = np.clip(state, 0.0, 1.0)
        return p["D0"] * s**(1 - p["q"]) * (1 - s) / (1 + (mu - 1) * s)

    def functions(rate):
        """Build the conservative ODE and its sparse analytic Jacobian."""
        coeff = 2 * np.pi * edges[1:-1] / dr

        def rhs(t, state):
            s = np.clip(state, 0.0, 1.0)
            psi = primitive(state)
            flux = np.zeros(n + 1)
            flux[0] = rate
            # Outward injection uses the upstream (left) fractional flow.
            # The diffusion flux is a centered primitive difference.
            flux[1:-1] = (
                rate * mu * s[:-1] / (1 + (mu - 1) * s[:-1])
                - coeff * np.diff(psi)
            )
            # flux[-1] stays zero. Every interior flux cancels on summation.
            return -np.diff(flux) / area

        def jac(t, state):
            s = np.clip(state, 0.0, 1.0)
            d = diffusivity(state)
            fprime = np.where(
                (state >= 0.0) & (state <= 1.0),
                mu / (1 + (mu - 1) * s)**2, 0.0,
            )
            left = rate * fprime[:-1] + coeff * d[:-1]
            right = -coeff * d[1:]
            diagonal = np.zeros(n)
            diagonal[:-1] -= left / area[:-1]
            diagonal[1:] += right / area[1:]
            return diags(
                [left / area[1:], diagonal, -right / area[:-1]],
                [-1, 0, 1], format="csc",
            )
        return rhs, jac

    # Saved times include the exact Figure 2 snapshots. BDF chooses its own
    # internal time steps; t_eval only determines where states are returned.
    injection_snapshots = np.array([p["ts"] / 3, 2 * p["ts"] / 3, p["ts"]])
    remaining = p["tend"] - p["ts"]
    shutin_snapshots = p["ts"] + remaining * np.array([0, 1/15, 1/3, 1])
    tin = np.unique(np.r_[0, np.geomspace(1e-5, p["ts"], 90),
                         injection_snapshots])
    tout = np.unique(np.r_[
        p["ts"], p["ts"] + np.geomspace(1e-5, remaining, 150),
        p["ts"] + np.array([1, 5, 20, 50, 100]), shutin_snapshots,
    ])
    tout = tout[(tout >= p["ts"]) & (tout <= p["tend"])]

    # First integration: independently prescribed injection into an empty
    # aquifer. The source is a boundary flux, never a fitted PDE residual.
    rhs, jac = functions(qu)
    injection = solve_ivp(
        rhs, (0, p["ts"]), np.zeros(n), method="BDF", jac=jac,
        rtol=rtol, atol=atol, t_eval=tin,
    )
    if not injection.success:
        raise RuntimeError(injection.message)

    # Second integration: switch the rate to zero at exactly ts. Copy the
    # full last state; this implements continuity of the same plume history.
    rhs, jac = functions(0.0)
    shutin = solve_ivp(
        rhs, (p["ts"], p["tend"]), injection.y[:, -1].copy(),
        method="BDF", jac=jac, rtol=rtol, atol=atol, t_eval=tout,
    )
    if not shutin.success:
        raise RuntimeError(shutin.message)
    t = np.r_[injection.t, shutin.t[1:]]
    u = np.column_stack([injection.y, shutin.y[:, 1:]])

    # "mass" is the historical data-array name for M_u = integral u dA [m^2].
    # Multiplication by storage gives the mobile reservoir volume [m^3].
    mass = area @ u
    expected = qu * np.minimum(t, p["ts"])

    def threshold_radius(profile, level):
        """Interpolate the outer crossing of a prescribed thickness level."""
        ids = np.flatnonzero(profile >= level)
        if not len(ids):
            return 0.0
        i = ids[-1]
        if i == n - 1:
            return r[i]
        return r[i] + dr * (profile[i] - level) / (profile[i] - profile[i + 1])

    # R001 is an outer thickness diagnostic; a099 measures the central
    # nearly full-thickness region for these monotone radial profiles.
    # Neither numerical threshold is an exact mathematical free boundary.
    r_edge = np.array([threshold_radius(v, 1e-3) for v in u.T])
    r_core = np.array([threshold_radius(v, 0.99) for v in u.T])
    r95 = []
    for v, inventory in zip(u.T, mass):
        # Cellwise constant u accumulates linearly with annular area, hence
        # interpolation in r^2 for the 95% inventory radius.
        r95.append(
            np.sqrt(np.interp(
                .95 * inventory, np.r_[0, np.cumsum(area * v)], edges**2
            )) if inventory > 0 else 0.0
        )
    maximum = np.max(u, axis=0)
    metrics = dict(
        n=n, dr=float(dr), rtol=rtol, atol=atol,
        minimum=float(u.min()), maximum=float(u.max()),
        lower_bound_violation=float(max(0, -u.min())),
        upper_bound_violation=float(max(0, u.max() - 1)),
        relative_mass_error=float(np.max(np.abs(mass - expected)) / expected[-1]),
        shutin_state_jump=float(np.max(np.abs(shutin.y[:, 0] - injection.y[:, -1]))),
        outermost_u=float(np.max(u[-1])),
        nfev_inj=injection.nfev, nfev_post=shutin.nfev,
    )
    return dict(
        params=p, r=r, edges=edges, area=area, t=t, u=u,
        mass=mass, expected=expected, R001=r_edge, a099=r_core,
        R95=np.array(r95), umax=maximum, metrics=metrics,
        injection_snapshots=injection_snapshots, shutin_snapshots=shutin_snapshots,
    )


def save(result, output):
    """Save lossless array data, readable exports and an audit summary."""
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output / "history.npz",
        **{k: v for k, v in result.items() if k not in ("params", "metrics")},
    )
    summary = dict(
        parameters=result["params"], **result["metrics"],
        software=dict(python=platform.python_version(),
                      numpy=np.__version__, scipy=scipy.__version__),
    )
    snapshots = np.unique(np.r_[
        result["injection_snapshots"], result["shutin_snapshots"],
    ])
    summary["snapshots"] = []
    for target in snapshots:
        i = int(np.argmin(abs(result["t"] - target)))
        summary["snapshots"].append(dict(
            time_year=float(result["t"][i]), umax=float(result["umax"][i]),
            R001_m=float(result["R001"][i]), a099_m=float(result["a099"][i]),
            R95_m=float(result["R95"][i]),
            mobile_volume_m3=float(result["mass"][i] * result["params"]["H"]
                                   * result["params"]["phi"]
                                   * (1 - result["params"]["Sbr"])),
        ))
    (output / "metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
    np.savetxt(
        output / "history.csv",
        np.column_stack([result[k] for k in
                         ("t", "mass", "expected", "R001", "a099", "R95", "umax")]),
        delimiter=",",
        header="t_year,mobile_area_m2,expected_mobile_area_m2,R_u001_m,a_u099_m,R95_m,umax",
        comments="",
    )
    ids = [np.argmin(abs(result["t"] - s)) for s in snapshots]
    np.savetxt(
        output / "profiles.csv",
        np.column_stack([result["r"], result["u"][:, ids]]),
        delimiter=",", comments="",
        header="r_m," + ",".join(f"u_t{result['t'][i]:g}_year" for i in ids),
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=int, default=6400)
    parser.add_argument("--output", type=Path, default=PACKAGE_ROOT / "results" / "recomputed" / "model")
    parser.add_argument("--rtol", type=float, default=1e-9)
    parser.add_argument("--atol", type=float, default=1e-12)
    args = parser.parse_args()
    save(solve(args.cells, args.rtol, args.atol), args.output)
