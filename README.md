# CO2 plume reproducibility

Python reproduction materials for **Generalized nonlinear Fokker–Planck dynamics of CO2 plumes during injection and after shut-in**, by Fernando Alonso-Marroquín and Christian Tantardini.

The repository reproduces manuscript **Figures 2–4** from the retained numerical history and processed field-radius data. It includes the finite-thickness solver, the Python port of the Figure 4 MATLAB fit, independent numerical checks, and the resulting PDF/PNG figures. All executable analysis code is Python.

**Data query before scientific release:** the first Aquistore map is labelled 2015 in White and Movahedzadeh (2025), whereas Roach and White (2018) date the 36 kt survey to February 2016. The supplied CSV is retained byte-for-byte, and the published-data fit reproduces those inputs. Resolve this date-to-map association with the source records before describing the field dates as independently verified. See [data provenance](docs/DATA_PROVENANCE.md).

## Quick start

Tested with **Python 3.12.14**, NumPy 2.3.5, SciPy 1.17.0 and Matplotlib 3.10.8. MATLAB and an external LaTeX installation are unnecessary.

```bash
python -m venv .venv
```

Activate the environment with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\Activate.ps1` in Windows PowerShell, then run:

```bash
python -m pip install -r requirements.txt
python scripts/reproduce.py
```

This redraws all three figures using the supplied 6400-cell history and the field CSV. The model input files remain unchanged. Outputs are written to `figures/` and `results/field/`.

## Recompute and verify

```bash
# Scientific regression checks (standard-library unittest).
python -m unittest discover -s tests -v

# Recompute the full 6400-cell history, redraw the figures, and run all checks.
python scripts/reproduce.py --recompute --validate
```

A fresh production history goes to `results/recomputed/model/`. The validation compares 800, 1600, 3200 and 6400 cells, tests the diffusion primitive by quadrature, solves the PME from the inherited shut-in state, and benchmarks a separate PME solve against the exact Barenblatt solution. Matching coarse-grid caches can be reused; incompatible model parameters are rejected.

The quick plotting command takes seconds on the preparation machine. Full verification includes several stiff PDE integrations and takes longer depending on hardware. No network access is needed after installing dependencies.

## Figure-to-input map

| Manuscript output | Python entry point | Retained inputs |
|---|---|---|
| Figure 2: injection and shut-in profiles | `scripts/plot_model_figures.py` | `data/model/history.npz`, `metrics.json` |
| Figure 3: continuous plume history | `scripts/plot_model_figures.py` | same model history and metadata |
| Figure 4: field-radius fits | `scripts/fit_field_radius.py` | `data/field/combined_radius_vs_year_all_sites.csv` |
| Model recomputation | `scripts/solve_plume.py` | documented physical parameters in the script |
| Spatial verification | `scripts/check_resolution.py` | current solver and production history |
| PME/Barenblatt comparison | `scripts/check_barenblatt.py` | current production history |

The figure filenames match the manuscript:

- `figures/continuous_plume_profiles.pdf` and `.png`
- `figures/continuous_plume_history.pdf` and `.png`
- `figures/all_sites_radius_powerlaw_loglog.pdf` and `.png`

Figures 2 and 3 use panel letters only. All plots use Computer Modern text/math fonts embedded in vector PDFs; PNGs use 300 dpi. Figure 4 is tightly cropped, so use it in LaTeX without the large `trim` values that compensated for the original MATLAB page margins:

```latex
\includegraphics[width=\linewidth]{figures/all_sites_radius_powerlaw_loglog.pdf}
```

Figure 1 is a separate schematic whose source asset was absent from the supplied archive. Raw seismic images, segmentation masks and extraction scripts were also absent. Figure 4 reproduction therefore starts from the supplied processed radius table; it does not repeat image segmentation.

## Figure 4: faithful fit, transparent uncertainties

The port preserves the supplied MATLAB algorithm: 10,000 onset candidates, unweighted ordinary least squares in natural-log coordinates, and selection by the smallest log-residual sum of squares. The fitted onset lies between injection start plus `1e-4` years and the first survey minus `1e-4` years. Full-precision results are in `results/field/radius_fit_parameters.csv` and `.json`.

The supplied MATLAB script replaced each first survey's time uncertainty with one month. The default Python plot uses the actual CSV uncertainties. To reproduce that older display convention explicitly:

```bash
python scripts/fit_field_radius.py --legacy-first-time-error \
  --figures results/legacy_display/figures --results results/legacy_display/fits
```

This flag changes only the displayed time error bars. It does not change the fits or source data. On logarithmic axes, horizontal errors extending below 0.1 yr are clipped at the displayed axis limit; their original endpoints are recorded in the JSON. Error bars are observational uncertainties, not fit weights or confidence intervals for the fitted onset.

The exact CSV/MATLAB algorithm gives the Sleipner onset `1998.4826827583`, which rounds to **1998.48**. The manuscript's earlier `1998.50` should be reconciled with this precise output. Aquistore and Weyburn–Midale reach the lower onset bound. Weyburn has only three observations; the high fitted R² values do not establish parameter uncertainty or validate a unique transport mechanism.

## Scientific scope

The model combines imposed radial displacement with the column fractional flow and a bounded nonlinear diffusivity. The 6400-cell example uses three years of injection and 300 years of shut-in, with a conserved mobile reservoir volume of 30 million m³ after shut-in. It is an illustrative model calculation, independent of the field-site fits.

The final maximum thickness fraction is about 0.7363. The computed transient remains in the finite-thickness regime. The PME and Barenblatt solutions are clearly identified comparison problems, and the threshold core radius is distinguished from an exact full-thickness interface. Equations, parameter meanings and numerical methods are documented in [METHODS.md](docs/METHODS.md).

## Repository contents

- `data/`: three authoritative input files plus their checksum manifest.
- `scripts/`: six commented Python entry points, including the reproduction driver.
- `tests/`: independent scientific regression checks.
- `figures/`: the manuscript plots generated by Python.
- `results/`: fit tables, numerical verification and run provenance.
- `docs/`: methods, data provenance and release instructions.
- `.github/workflows/`: automated checks and figure generation for GitHub.

The obsolete pure-PME histories and their incompatible refinement driver were excluded. Generated refinement states and fresh production states are reproducible and ignored by Git; concise verification reports are retained. The original MATLAB file is represented by its source hash and conversion notes, with the Python implementation used for execution.

## Citation and use

Citation metadata is supplied in [CITATION.cff](CITATION.cff). The repository URL and archival DOI can be added after creating the remote and depositing a release. Neither has been invented in this snapshot.

Rights remain with the respective authors and data rightsholders; see [LICENSE](LICENSE). An open-source licence can be adopted by the rightsholders for a public release. Instructions for creating the GitHub repository and preparing a citable archive are in [RELEASE.md](docs/RELEASE.md).
