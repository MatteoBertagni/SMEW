# Changelog

## 1.3.0b1 (Unreleased)

Changes since 1.2.0.

### Added

- **Nonlinear weathering:** rock surface contact with water can now depend on
  soil pore sizes, rock particle sizes, and the pore/particle mixing. Pore 
  distributions can be estimated from soil texture using Ding or Campbell models. 
  Constant and linear wetting remain available; linear is the default.
- **Freezing approximation:** when soil temperature is supplied to the dynamic
  water balance, temperatures at or below 0 °C hold soil moisture constant,
  suppress water fluxes, and store precipitation as snow/ice; above 0 °C,
  temperature-dependent melting supplies water. Independently, respiration is
  zero at or below 0 °C, and the main biogeochemistry model holds aqueous
  chemistry fixed and stops weathering and nutrient uptake, while gas-phase
  CO₂ diffusion continues. Vegetation growth has no freezing cutoff.
- **Seasonal vegetation:** add `veg_seasonal` for logistic growth over prescribed
  crop seasons, returning separate vegetation biomass and harvested biomass
  time series.

### Changed

- **Soil and carbon balance:** updated soil moisture thresholds to use the Campbell
  water retention relation. Respiration now uses the average of temperatures above
  0 °C as its reference. Steady-state carbon estimates now combine moisture and
  temperature effects at each timestep before averaging.
- **Numerical architecture:** shared equations in `smew.equations` support Python
  and compiled execution. Numba compiles the main numerical routines; Cython and
  bundled CMINPACK handle compiled nonlinear solves with reusable work arrays.
  Set `NUMBA_DISABLE_JIT=1` before import to use Python/SciPy.
- **Errors and interfaces:** numerical failures expose `SMEWError` exceptions with
  unique codes and diagnostic context. Public `*_numba` entry points return
  results, codes, and messages, allowing failed simulations to release arrays
  before Python raises an exception.
- **Packaging and verification:** require Python 3.11+ and NumPy 2+; add Numba as
  a dependency, compiled wheels for Linux/macOS/Windows, and a shared version
  definition. Add a marimo example and automated regression, physical-range,
  Python/compiled agreement, memory-release, and package checks.

### Fixed

- Account for excess water at the numerical saturation cap; correct soil-class
  handling and vegetation boundary indexing; preserve linear wet-surface scaling
  in the two-application weathering model.

### Deprecated

- `biogeochem_balance2psd` now emits a deprecation warning and is scheduled for
  removal in a future release.
