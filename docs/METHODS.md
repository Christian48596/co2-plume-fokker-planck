# Methods and reproducibility

## Finite-thickness model (Figures 2 and 3)

The connected mobile-layer fraction is `u = h/H`. In the exterior radial domain `r >= r_w`, the normalized integrated current is

```text
F = Q_hat(t) f_VE(u) - 2 pi r D_u(u) du/dr
f_VE(u) = M u / [1 + (M - 1)u]
D_u(u) = D0 u^(1-q) (1-u) / [1 + (M - 1)u]
Q_hat(t) = QV(t) / [phi (1-Sbr) H]
u_t = -(1/(2 pi r)) dF/dr
```

The well boundary supplies `F(r_w) = Q_hat`, and the far boundary has zero mobile-gas flux. Injection is constant until `t_s`; afterwards its imposed rate is zero. The complete computed state at `t_s` becomes the initial state of the shut-in integration. The baseline conserves connected mobile inventory and omits exchange, residual trapping and compressible pressure relaxation.

| Parameter | Value | Unit / meaning |
|---|---:|---|
| `H` | 10 | m, aquifer thickness |
| `phi` | 0.30 | porosity |
| `Sbr` | 0.20 | residual brine saturation |
| `q` | 0.30 | constitutive exponent |
| `D0` | 20000 | m²/yr, spreading scale |
| `mobility_ratio` | 5 | endpoint gas/brine mobility ratio |
| `QV` | 10000000 | reservoir m³/yr during injection |
| `ts` | 3 | yr from injection start |
| `tend` | 303 | yr from injection start |
| `rw` | 10 | m, inner radius |
| `rout` | 20000 | m, outer radius |
| production cells | 6400 | uniform radial mesh |

The q=0 gravity coefficient is the Darcy-derived member. Positive q specifies a constitutive extension. The code prescribes the total effective gradient coefficient; it does not independently resolve a capillary-pressure relation or infer capillary properties from q.

### Discretization

Cell areas are `pi*(r_right²-r_left²)`. Every interior face uses one shared flux, so interior contributions cancel in the inventory sum. Pressure-driven flux uses the upstream fractional flow. Diffusion uses the centered difference of the primitive `Psi(u)`, with `Psi'(u)=D_u(u)`. An analytic Gauss hypergeometric expression evaluates Psi. Independent numerical quadrature checks this primitive.

Time integration uses SciPy BDF, a sparse analytic Jacobian, relative tolerance `1e-9` and absolute tolerance `1e-12`. Constitutive functions have constant extensions outside `[0,1]` for implicit iterations. The evolving state is never clipped or renormalized. The small upper-bound excess in the numerical output is reported.

### Diagnostics and units

- `mass` is the legacy array name for `M_u = integral u dA`, with units m².
- Reservoir volume is `phi*(1-Sbr)*H*M_u`, in m³.
- `R001` is the interpolated outer crossing of u=0.001.
- `a099` is the outer crossing of u=0.99 in the central connected region for these monotone radial profiles.
- `R95` encloses 95% of M_u, interpolated in radial area.
- `umax` is the maximum cell-average thickness fraction.

Threshold radii are numerical diagnostics, distinct from exact free boundaries or seismic detection thresholds. The core radius at shut-in changes by approximately 2.33% between 3200 and 6400 cells; exact inventory conservation does not establish profile accuracy. The full four-grid results are stored in `results/refinement/`.

## PME and Barenblatt checks

After shut-in, the thin limit is

```text
u_tau = D0/(2-q) * Laplacian(u^(2-q)),   tau = t-ts.
```

For the full-plane radial Barenblatt family,

```text
u_B(r,theta) = A(theta) [1-r²/R(theta)²]_+^(1/(1-q))
A(theta) = [M_s/(4 pi D0 theta)]^(1/(2-q))
R(theta)² = 4 D0 (2-q)/(1-q) theta A(theta)^(1-q)
theta = tau + theta0
```

The comparison script performs two separate integrations:

1. Pure PME from the identical saved full-model shut-in state, using the same finite-well mesh and zero boundary fluxes. This isolates the diffusivity change.
2. Pure PME on a full-plane radial mesh starting from exact annular averages of a finite-age Barenblatt solution. This tests the coefficient and discretization against an exact evolving solution.

A mass-and-initial-peak match sets `theta0 = M_s/(4*pi*D0) = 49.7359 yr`. It selects a reference shape, not the actual inherited profile. The script also records a second-moment matching convention and an unshifted reference, clearly separated from the full-model solution.

At tau=300 yr the full model gives umax≈0.73631; PME from the same inherited state gives ≈0.29775; the mass-and-peak-matched analytical profile gives ≈0.31749. The full-to-PME diffusivity ratio at the full-model peak is ≈0.06684. These are model comparisons, not field validation.

## Field fits (Figure 4)

The observable is image-derived area-equivalent radius in km. The model is

```text
R_eq = R0 [(t-t0)/(1 yr)]^beta
```

`R0` has units km. Decimal calendar years represent t and t0. For each site, 10,000 uniformly spaced t0 candidates cover `[InjectionStartYear+1e-4, min(SurveyYear)-1e-4]`. For each candidate, ordinary least squares fits `log(R_eq)` against `log(t-t0)`. The candidate with minimum unweighted log-space SSE is retained; R² is evaluated in that same log space. This reproduces the supplied MATLAB algorithm.

Source radius and survey-time uncertainties are plotted but do not enter the fitting weights. No parameter confidence intervals are inferred. The displayed horizontal bars condition on the fitted onset; they do not include onset uncertainty. Their full endpoints are saved even when the logarithmic display clips a left endpoint below 0.1 yr. The first-survey one-month override from the MATLAB file is available only by an explicit flag.

The field observations concern injection. The numerical model parameters are not fitted from these site exponents. The onset is an apparent mapped-response onset; it is distinct from both injection start and shut-in.

## Validation boundaries

The repository verifies numerical conservation, discretization sensitivity, exact-PME agreement, data identity and faithful porting of the fit. It starts from processed field radii: original segmentation masks and extraction code were not supplied. It does not independently invert seismic saturation, calibrate permeability or resolve the Aquistore date-label discrepancy described in the provenance document.
