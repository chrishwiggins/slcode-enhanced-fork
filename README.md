# SLcode - Sea Level Equation Solver

A comprehensive toolkit for computing global sea level changes from ice sheet variations using glacial isostatic adjustment (GIA) modeling. This repository implements elastic and viscoelastic Earth models to calculate sea level fingerprints from ice loading/unloading scenarios.

## Overview

SLcode solves the gravitationally self-consistent sea level equation following the methodologies of Kendall et al. (2005) and Austermann et al. (2015). The code computes how sea level changes globally when ice masses change, accounting for:

- **Gravitational effects**: Changes in Earth's gravitational field due to ice and water redistribution
- **Solid Earth deformation**: Elastic and viscoelastic response to ice loading/unloading
- **Rotational effects**: Changes in Earth's rotation due to mass redistribution
- **Ocean geometry changes**: Self-consistent updating of coastlines and ocean basins

## Key Features

- **Multiple Earth models**: Elastic, viscoelastic, and benchmarked implementations
- **Spherical harmonic methods**: Efficient global computations using spectral techniques
- **Real ice sheet data**: Includes West Antarctic Ice Sheet collapse scenarios
- **Comprehensive Love numbers**: Pre-computed elastic parameters for various Earth models
- **Visualization tools**: Built-in plotting capabilities for sea level fingerprints

## Quick Start

### Prerequisites

#### Python Implementation
```bash
# Install dependencies
make install-deps

# Or manually:
pip3 install numpy scipy matplotlib pyshtools astropy xarray
```

#### MATLAB Implementation (Optional)
- MATLAB R2015b or later
- Signal Processing Toolbox recommended

### Running the Code

#### Default: Python Elastic Model
```bash
# Run the main sea level solver
make

# Or explicitly:
make python-elastic
```

#### MATLAB Versions
```bash
# Elastic model (MATLAB)
make matlab-elastic

# Viscoelastic model (MATLAB)
make matlab-viscoelastic

# All benchmarked versions
make matlab-benchmarked
```

#### Get Help
```bash
make help
```

## Repository Structure

### Main Solvers
- `SL_equation_elastic.py` - **Primary Python implementation** (Updated for Python 3, uses pyshtools)
- `SL_equation_elastic_benchmarked.m` - MATLAB elastic model
- `SL_equation_viscoelastic_benchmarked.m` - MATLAB viscoelastic model
- `SL_equation_viscoelastic_GIAonly.m` - MATLAB GIA-only model
- `SL_equation_DT.m` - MATLAB dynamic topography model

### Alternative Implementation
- `SLcode_py/` - Alternative Python implementation (legacy pyspharm-based)

### Data Files
- `SavedLN/` - Love numbers for various Earth models
- `SavedLN_EP/` - Extended precision Love numbers
- `ice_grid/` - Ice sheet data (West Antarctic Ice Sheet scenarios)
- `SLcode_py/gebco_08_15am.mat` - Global topography/bathymetry data

### Supporting Functions
- `SLFunctions/` - MATLAB utility functions for spherical harmonics, Love numbers, etc.
- `plot_output/` - Plotting tools for visualization

### Research Applications
- `GRL_2022_proglacial_lakes/` - Proglacial lakes study implementation
- `CSDMS_workshop/` - Workshop materials and examples
- `benchmark_in_out/` - Benchmark test cases

## Scientific Background

### The Sea Level Equation

The code solves the gravitationally self-consistent sea level equation:

```
S(φ,λ,t) = N(φ,λ,t) - U(φ,λ,t)/g
```

Where:
- `S(φ,λ,t)` is the sea level change
- `N(φ,λ,t)` is the solid Earth deformation
- `U(φ,λ,t)` is the gravitational potential change
- `φ,λ` are latitude and longitude
- `t` is time

### Earth Models

#### Elastic Model
- Instantaneous response to loading
- Uses elastic Love numbers (h, k, l)
- Suitable for modern ice mass changes

#### Viscoelastic Model
- Time-dependent response
- Accounts for mantle viscosity
- Suitable for glacial timescales

### Spherical Harmonic Expansion

The code expands all fields in spherical harmonics:

```
f(φ,λ) = Σ Σ f_lm Y_lm(φ,λ)
         l m
```

Using optimized transforms via pyshtools for efficiency up to degree ~2800.

## Key Parameters

### Physical Constants
- `rho_ice = 916.7 kg/m³` - Ice density
- `rho_water = 1000 kg/m³` - Water density
- `rho_sed = 2300 kg/m³` - Sediment density
- `g = 9.80616 m/s²` - Gravitational acceleration
- `rsphere = 6.37122e6 m` - Earth radius

### Numerical Parameters
- `maxdeg = 64` - Maximum spherical harmonic degree
- `k_max = 10` - Maximum iterations for convergence
- `epsilon = 1e-4` - Convergence criterion

## Input Data

### Ice Sheet Scenarios
- `ice_grid/WAIS.mat` - West Antarctic Ice Sheet data
  - `ice_Ant` - Present-day ice thickness
  - `ice_EAIS` - Future ice configuration
  - `lat_WAIS`, `lon_WAIS` - Coordinate grids

### Earth Model Parameters
- `SavedLN/prem.l90C.umVM2.lmVM2.mat` - Love numbers
  - `h_el` - Elastic displacement Love numbers
  - `k_el` - Elastic gravitational Love numbers
  - `h_el_tide`, `k_el_tide` - Tidal Love numbers

### Topography
- `SLcode_py/gebco_08_15am.mat` - Global relief model

## Output

### Sea Level Fingerprint
The code produces a global map of relative sea level change normalized to the global mean sea level equivalent of the ice change.

### Key Variables
- `plotSL` - Normalized sea level fingerprint
- `delSL` - Raw sea level change
- `scaling_fact` - Normalization factor

## Recent Updates

### Python 3 Compatibility (2024)
- Updated all print statements for Python 3
- Replaced deprecated `numpy.complex` with `complex`
- Fixed deprecated `scipy.interpolate.interp2d` warnings
- Implemented pyshtools integration for spherical harmonic transforms

### Spherical Harmonic Library Migration
- **Migrated from pyspharm to pyshtools** for better maintenance and accuracy
- pyshtools provides:
  - Proven accuracy for spherical harmonic degrees up to 2800
  - Multiple grid format support (DH, GLQ)
  - Active development and comprehensive documentation
  - Fast Fortran backend with Python interface

### Grid Format Handling
- Automatic interpolation to pyshtools-compatible grids
- Support for Driscoll-Healy (DH) and Gauss-Legendre Quadrature (GLQ) grids
- Proper m-primary to l-primary coefficient ordering conversion

## References

1. Kendall, R.A., Mitrovica, J.X., and Milne, G.A. (2005). On post-glacial sea level - II. Numerical formulation and comparative results on spherically symmetric models. *Geophysical Journal International*, 161, 679-706.

2. Austermann, J., Mitrovica, J.X., Latychev, K., and Milne, G.A. (2015). Barbados-based estimate of ice volume at Last Glacial Maximum affected by subsequent sea-level change. *Nature Geoscience*, 6, 3-6.

3. Mitrovica, J.X., and Wahr, J. (2005). Time-variable gravity from GRACE: First results. *Geophysical Research Letters*, 32, L11402.

4. Wieczorek, M. A. and Meschede, M. (2018). SHTools — Tools for working with spherical harmonics, *Geochemistry, Geophysics, Geosystems*, 19, 2574-2592.

## Citation

If you use this code in your research, please cite:

```bibtex
@software{slcode2024,
  title={SLcode: Sea Level Equation Solver},
  author={Austermann, Jacqueline and Wickert, Andrew},
  year={2024},
  note={Updated for Python 3 with pyshtools integration},
  url={https://github.com/your-repo/slcode-enhanced}
}
```

## Support

- For questions about the scientific methodology: See the references above
- For technical issues: Create an issue in the repository
- For MATLAB compatibility: Ensure you have required toolboxes

## License

Academic use permitted. Please cite the relevant papers when using this code in publications.

---

*This implementation has been updated and enhanced for modern Python environments while preserving the core scientific algorithms developed by J. Austermann and A. Wickert.*