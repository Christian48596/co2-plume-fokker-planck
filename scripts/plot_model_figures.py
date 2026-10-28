#!/usr/bin/env python3
"""Draw Figures 2 and 3 from the saved bounded-model numerical solution.

Run after solve_plume.py, or use the supplied data/model/:
    python scripts/plot_model_figures.py

No solution is fitted or reconstructed here. Every curve is taken from the
computed finite-volume history; the balance line is the prescribed injection
integral. Both figures are exported as vector PDF and 300-dpi PNG.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")  # permits figure generation without a graphical desktop
import matplotlib.pyplot as plt
from matplotlib import font_manager


PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def main(data_dir, output_dir):
    """Read the numerical arrays and make the two manuscript figures."""
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    with np.load(data_dir / "history.npz", allow_pickle=False) as archive:
        d = {name: archive[name] for name in archive.files}
    summary = json.loads((data_dir / "metrics.json").read_text())
    p = summary["parameters"]
    t, u = d["t"], d["u"]
    ts, tend = p["ts"], p["tend"]

    # Computer Modern text and math match mathematical manuscript typography.
    # cmr10.ttf ships with Matplotlib and embeds as TrueType in the PDF.
    # This avoids depending on an external TeX installation or OTF font support.
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["cmr10", "DejaVu Serif"],
        "mathtext.fontset": "cm", "font.size": 10,
        "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
        "axes.labelsize": 10, "legend.fontsize": 8.4,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": 0.7, "lines.linewidth": 1.7,
    })

    def index_at(time):
        # An exact saved state is required, so the plot cannot silently use a
        # different time when the data or snapshot selections are changed.
        index = int(np.argmin(abs(t - time)))
        if not np.isclose(t[index], time, rtol=1e-10, atol=1e-10):
            raise ValueError(f"Requested snapshot {time:g} yr is absent.")
        return index

    def save_figure(fig, name):
        # Fixed PDF metadata avoids a changing creation timestamp in reruns.
        fig.savefig(output_dir / f"{name}.pdf",
                    metadata={"CreationDate": None, "ModDate": None})
        fig.savefig(output_dir / f"{name}.png", dpi=300)
        plt.close(fig)

    # FIGURE 2: the same plume before and after shut-in.
    # Shared linear axes show the nearly full-thickness plateau and its tail.
    # The inherited t=ts curve is present in both panels using the same colour.
    # Radii are converted from m to km; all times are from injection start.
    radius_limit = float(np.ceil(1.05 * np.max(d["R001"]) / 1000))
    radius_limit = max(radius_limit, 1.0)
    fig, axes = plt.subplots(
        1, 2, figsize=(7.1, 2.85), sharex=True, sharey=True,
        layout="constrained",
    )
    times_by_panel = [d["injection_snapshots"], d["shutin_snapshots"]]
    colors_by_panel = [
        ["#c6a04c", "#a86727", "#6b3615"],
        ["#6b3615", "#4c8fa8", "#226481", "#173b58"],
    ]
    for ax, times, colors, label in zip(
            axes, times_by_panel, colors_by_panel, ["(a)", "(b)"]):
        for time, color in zip(times, colors):
            index = index_at(time)
            legend = rf"$t={time:g}$ yr"
            if np.isclose(time, ts):
                legend += " (shut-in)"
            ax.plot(d["r"] / 1000, u[:, index], color=color, label=legend)
        ax.set(xlim=(0, radius_limit), ylim=(-0.015, 1.055),
               xlabel=r"Radius, $r$ (km)")
        # Only the panel letter occupies the title position.
        ax.set_title(label, loc="left", fontsize=10)
        ax.grid(alpha=0.16)
        ax.legend(frameon=False, loc="upper right")
    axes[0].set_ylabel(r"Mobile layer fraction, $u=h/H$")
    save_figure(fig, "continuous_plume_profiles")

    # FIGURE 3: injection schedule, inventory, outer radii, and core radius.
    # The schedule and scales come from the saved physical parameters.
    volume_exponent = int(np.floor(np.log10(p["QV"])))
    volume_scale = 10.0**volume_exponent
    storage = p["phi"] * (1 - p["Sbr"]) * p["H"]
    volume = d["mass"] * storage / volume_scale
    expected_volume = p["QV"] * np.minimum(t, ts) / volume_scale
    rate_window = min(tend, 10 * ts / 3)
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.0), layout="constrained")
    rate_ax, inventory_ax, radius_ax, core_ax = axes.flat

    # (a) Prescribed rate. The short time range resolves its abrupt shut-in.
    rate_ax.plot(
        [0, ts, ts, rate_window],
        [p["QV"] / volume_scale, p["QV"] / volume_scale, 0, 0],
        color="#9b541e",
    )
    rate_ax.set(
        xlim=(0, rate_window), ylim=(-0.05, 1.2 * p["QV"] / volume_scale),
        ylabel=rf"Injection ($10^{{{volume_exponent}}}$ m$^3$/yr)",
    )
    rate_ax.annotate(
        "Shut-in", xy=(ts, 0.5 * p["QV"] / volume_scale),
        xytext=(1.5 * ts, 0.65 * p["QV"] / volume_scale),
        arrowprops={"arrowstyle": "->", "lw": 0.8}, fontsize=8.4,
    )

    # (b) M_u has units m^2, despite the historical array name "mass".
    # Multiplying by phi*(1-Sbr)*H converts it to mobile reservoir volume.
    inventory_ax.plot(t, expected_volume, color="#333333",
                      lw=1, label="Exact balance")
    samples = np.unique(np.r_[np.arange(0, len(t), 12), len(t) - 1])
    inventory_ax.plot(t[samples], volume[samples], "o", ms=3, mfc="none",
                      color="#176383", label="Finite volume")
    inventory_ax.set(
        xlim=(0, tend), ylim=(-0.03 * expected_volume[-1],
                             1.12 * expected_volume[-1]),
        ylabel=rf"Mobile volume ($10^{{{volume_exponent}}}$ m$^3$)",
    )
    inventory_ax.legend(frameon=False, loc="lower right")

    # (c) R001 is the radius at u=0.001; R95 encloses 95% of mobile inventory.
    # These are diagnostic radii, distinct from the exact support edge.
    radius_ax.plot(t, d["R001"] / 1000, color="#176383",
                   label=r"$R_{0.001}$: thickness threshold")
    radius_ax.plot(t, d["R95"] / 1000, color="#a86727",
                   label=r"$R_{95}$: 95% of mobile inventory")
    radius_ax.set(xlim=(0, tend), ylim=(0, radius_limit),
                  ylabel="Plume radius (km)")
    radius_ax.legend(frameon=False, loc="lower right")

    # (d) a099 measures the central u>=0.99 region for these monotone profiles.
    # Zero means this threshold is no longer attained. It is not an exact
    # determination of the time at which a mathematical u=1 core vanishes.
    core_ax.plot(t, d["a099"] / 1000, color="#6b3615")
    core_ax.set(xlim=(0, tend), ylim=(0, 1.14 * np.max(d["a099"]) / 1000),
                ylabel=r"Nearly full core radius, $a_{0.99}$ (km)")

    for ax, label in zip(axes.flat, ["(a)", "(b)", "(c)", "(d)"]):
        ax.set_title(label, loc="left", fontsize=10)
        ax.set_xlabel("Time from injection start (yr)")
        ax.grid(alpha=0.16)
    for ax in [inventory_ax, radius_ax, core_ax]:
        ax.axvline(ts, color="#777777", ls=":", lw=0.8)
    save_figure(fig, "continuous_plume_history")

    # Record which numerical archive generated the figures.
    provenance = {
        "data_sha256": hashlib.sha256((data_dir / "history.npz").read_bytes()).hexdigest(),
        "matplotlib": matplotlib.__version__,
        "text_font": Path(font_manager.findfont("serif")).name,
        "math_fontset": plt.rcParams["mathtext.fontset"],
        "injection_snapshots_year": d["injection_snapshots"].tolist(),
        "shutin_snapshots_year": d["shutin_snapshots"].tolist(),
        "pdf_and_png_names": ["continuous_plume_profiles", "continuous_plume_history"],
    }
    (output_dir / "figure_provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=PACKAGE_ROOT / "data" / "model")
    parser.add_argument("--output", type=Path, default=PACKAGE_ROOT / "figures")
    args = parser.parse_args()
    main(args.data, args.output)
