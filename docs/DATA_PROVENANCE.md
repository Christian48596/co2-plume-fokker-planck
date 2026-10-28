# Data provenance and curation

## Retained inputs

The source is the author-supplied `python_rep.zip`. Its SHA256 and the exact hashes of all retained inputs are in `data/input_manifest.json`. The three retained input files are unchanged from that archive.

| File | Use | Provenance |
|---|---|---|
| `data/model/history.npz` | Figures 2/3 and numerical checks | 6400-cell finite-thickness model history |
| `data/model/metrics.json` | Parameters, tolerances, dimensions and diagnostic metadata | Same production calculation |
| `data/field/combined_radius_vs_year_all_sites.csv` | Figure 4 | Author-supplied processed image-derived radii and uncertainties |

`history.npz` stores `u[cell,time]`, radial centres/edges, annular areas, saved times, inventory, threshold radii and snapshot selections. `metrics.json` documents physical parameters and software versions. The production history contains injection at QV=10^7 m³/yr for 3 years, in a 20 km outer domain, followed by 300 years of shut-in.

The CSV contains 14 observations: six Sleipner top-layer maps, five Aquistore upper-Deadwood maps, and three Weyburn–Midale maps. Year fields are decimal calendar years; time errors are half-widths in years. Radius and radius error are in km; mapped area and its error are in km². Injected mass, when supplied, is in kilotonnes. The layer, method and original source-filename fields are retained as provenance, even where the regression does not use them.

### Field sources

- **Sleipner:** Boait, F. C., et al. (2012), *Spatial and temporal evolution of injected CO2 at the Sleipner Field, North Sea*, Journal of Geophysical Research 117, B03309. https://doi.org/10.1029/2011JB008603 . The manuscript uses the topmost mapped layer in Figure 4.
- **Aquistore:** White, D. J., and Movahedzadeh, Z. (2025), *4D seismic tracking of the CO2 plume at the Aquistore CO2 storage site*, GeoConvention 2025, Figure 1. https://geoconvention.com/wp-content/uploads/abstracts/2025/114794-4d-seismic-tracking-of-the-co2-plume-at . The text gives 141 kt for the third survey while its panel is labelled 150 kt; the supplied CSV uses 141 kt.
- **Aquistore chronology:** Roach, L. A. N., and White, D. J. (2018), *Evolution of a deep CO2 plume from time-lapse seismic imaging at the Aquistore storage site, Saskatchewan, Canada*, International Journal of Greenhouse Gas Control 74, 79–86. https://doi.org/10.1016/j.ijggc.2018.04.025 .
- **Weyburn–Midale:** White, D. (2009), *Monitoring CO2 storage during EOR at the Weyburn-Midale Field*, The Leading Edge 28(7), 838–842. https://doi.org/10.1190/1.3167786 . The manuscript uses Figure 5.

These publications identify the mapping sources. The numeric radius table is an author-processed derivative, not a supplied field operator's direct saturation measurement. The original images, calibrated masks, intermediate extraction tables named in `SourceFile`, and segmentation scripts were absent from the input archive. The repository reproduces the fitting and plotting stage from the retained CSV.

### Aquistore date query requiring author resolution

The 2025 Figure 1 labels the five maps **2015, 2016, 2018, 2020 and 2023**. Its first map is marked 36 kt. Roach and White (2018) state that the 36 kt and 102 kt surveys were acquired in **February 2016 and November 2016**, respectively.

The supplied first row uses SurveyYear=2015.62278767123 with half-width 0.37621232876711 yr, a window ending in 2015.999. This window does not cover February 2016. The discrepancy could involve a later figure-label error or the date-to-panel mapping. Resolving it requires the original survey/panel records; no revised date has been invented here. The second row's 2016 window includes November 2016.

Before scientific release, confirm the correct first survey date and its association with the extracted radius. If a correction is justified, update the CSV with a documented source, regenerate its checksum, rerun the fits and revise the manuscript table/figure accordingly. The present fit remains a faithful reproduction of the supplied data, with this limitation explicit.

## Excluded legacy material

The following inputs belong to the earlier pure-PME calculation, which used Q0=1000 m³/yr, a 2 km outer boundary and 3000 m³ after injection. They do not support the current finite-thickness Figures 2/3:

- all `pme_N800`, `pme_N1600` and `pme_N3200` NPZ/CSV/summary files;
- `descriptive_scaling_fits.json`;
- the old `resolution.csv` and `resolution_summary.json`;
- the old `check_resolution.py`, which imported a missing, incompatible solver module.

The replacement refinement driver uses the current finite-thickness solver and produces fresh verification reports. The supplied Barenblatt report is regenerated from its Python computation. The old README referred to a different package and has been replaced. Platform metadata files are excluded.

The original MATLAB executable is excluded from the runnable repository because Figure 4 is fully ported to Python. Its hash and the faithful conversion details are retained in `results/validation/source_audit.json`. The port also removes the implicit first-time-error override from the default display, while making the historical display reproducible through an explicit option.

## Rights and citation

Cite the software/manuscript and the source publications when using the corresponding materials. Source publications and their figures retain their original rights. This repository does not redistribute the publication images or claim additional rights over them. See the repository LICENSE for the status of this prepared research snapshot.
