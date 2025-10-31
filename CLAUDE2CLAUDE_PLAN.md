# Claude-to-Claude Implementation Plan: Python Viscoelastic Development

## To: Future Claude Developer
## From: Claude (Sonnet 4)
## Subject: Complete Plan for Python Viscoelastic Implementation

---

## Executive Summary

**Goal**: Port the complete MATLAB viscoelastic functionality to Python, implementing the full Tromp (1999) → Kendall (2005) normal-mode theory that is currently only available in MATLAB files.

**Current Status**: Python implementation is elastic-only. MATLAB has full viscoelastic time-dependent response with exponential relaxation from normal-mode theory.

**Key Challenge**: The core mathematical framework involves time-dependent convolution integrals with exponential kernels derived from Earth's viscoelastic normal modes.

---

## Phase 1: Infrastructure and Data Loading

### 1.1 Create New Python File
**File**: `SL_equation_viscoelastic.py`
**Base**: Current `SL_equation_elastic.py` structure
**Purpose**: Full time-dependent viscoelastic sea level calculations

### 1.2 Critical MATLAB Data Loading
**Source File**: `SavedLN/prem.l90C.umVM2.lmVM2.mat`

**Required Arrays**:
```python
# Core normal-mode parameters (from Tromp 1999 theory)
spoles      # Shape: (maxdeg, max_modes) - Normal-mode decay rates s_{ℓ,k}
k_amp       # Shape: (maxdeg, max_modes) - Gravitational Love number amplitudes
h_amp       # Shape: (maxdeg, max_modes) - Displacement Love number amplitudes
mode_found  # Shape: (maxdeg,) - Number of computed modes per degree ℓ

# Tidal variants for rotational effects
k_amp_tide  # Tidal gravitational amplitudes
h_amp_tide  # Tidal displacement amplitudes
```

**Physical Meaning**:
- `spoles[ℓ,k] = s_{ℓ,k}`: Characteristic decay rates of Earth's viscoelastic response
- `k_amp[ℓ,k], h_amp[ℓ,k]`: Coupling coefficients from normal-mode eigenfunctions
- Relaxation times: `τ_{ℓ,k} = 1/s_{ℓ,k}`

### 1.3 Time Series Infrastructure
**Requirements**:
```python
ice_time_new    # Time array [kyr] for temporal evolution
ice             # 3D array: (lat, lon, time) - Ice thickness history
topo            # 3D array: (lat, lon, time) - Topography evolution
sed             # 3D array: (lat, lon, time) - Sediment loading (optional)
```

---

## Phase 2: Core Mathematical Implementation

### 2.1 Beta Coefficient Calculation (THE CRITICAL FUNCTION)

This is the heart of the Tromp (1999) theory implementation:

```python
def calculate_beta_coefficients(ice_time_new, spoles, k_amp, h_amp, mode_found, maxdeg):
    """
    Calculate time-dependent viscoelastic Love number coefficients.

    Implements the core formula from Tromp & Mitrovica (1999):
    β_ℓ(t) = Σ_k [(k_amp_{ℓ,k} - h_amp_{ℓ,k})/s_{ℓ,k}] * [1 - exp(-s_{ℓ,k} * t)]

    This is the mathematical bridge between normal-mode theory and
    practical sea level calculations.
    """
    n_times = len(ice_time_new)
    beta_l = []  # Store β coefficients for each time step pair

    # Time iteration: for each current time t_it
    for t_it in range(1, n_times):  # MATLAB: t_it = 2:length(ice_time_new)

        # History matrix: (n_previous_times, maxdeg)
        beta_matrix = np.zeros((t_it-1, maxdeg))

        # For each previous time step n
        for n in range(1, t_it):  # MATLAB: n = 2:t_it-1
            beta_row = np.zeros(maxdeg)

            # For each spherical harmonic degree ℓ
            for lm in range(maxdeg):
                num_mod = mode_found[lm]
                if num_mod > 0:
                    # Time difference for this step
                    time_diff = ice_time_new[t_it] - ice_time_new[n]

                    # Core exponential relaxation formula
                    exp_terms = 1.0 - np.exp(-spoles[lm, :num_mod] * time_diff)
                    amplitudes = (k_amp[lm, :num_mod] - h_amp[lm, :num_mod])

                    # Sum over all normal modes for this degree
                    beta_row[lm] = np.sum(amplitudes / spoles[lm, :num_mod] * exp_terms)

            beta_matrix[n-1, :] = beta_row

        beta_l.append(beta_matrix)

    return beta_l
```

**Why This Matters**: This function encodes the complete viscoelastic memory of the Earth - how past loading events influence current response through exponential relaxation.

### 2.2 Tidal Beta Coefficients
```python
def calculate_tidal_beta_coefficients(ice_time_new, spoles, k_amp_tide, h_amp_tide, mode_found, maxdeg):
    """
    Same structure as beta_l but for tidal/rotational Love numbers.
    Only degree 2 (ℓ=2) matters for rotational effects.
    """
    # Identical algorithm, different input arrays
    # Critical for Earth rotation coupling
```

### 2.3 Viscous Contribution Calculator
```python
def calculate_viscous_contribution(beta_l, sdelL_lm, t_it, beta_counter, T_lm):
    """
    Apply viscoelastic memory to compute current response.

    This implements the convolution integral:
    V_ℓm(t) = Σ_{t'<t} β_ℓ(t-t') * ΔL_ℓm(t')

    Physical meaning: Current response depends on entire loading history,
    weighted by exponential decay kernels.
    """
    V_lm = np.zeros(len(T_lm))

    if t_it > 1:  # No viscous response at first time step
        for lm_it in range(len(T_lm)):
            # Dot product: β history * load history
            V_lm[lm_it] = np.dot(beta_l[t_it-1][:, beta_counter[lm_it]],
                                sdelL_lm[:t_it-1, lm_it])

    return V_lm
```

---

## Phase 3: Sea Level Solver Integration

### 3.1 Enhanced Sea Level Response Equation

**Replace Current Elastic Equation**:
```python
# OLD (elastic only):
delSLcurl_ml_fl = E_ml * T_ml * (rho_ice*deli_ml + rho_water*delS_ml + rho_sed*Sed_ml)

# NEW (full viscoelastic):
delSLcurl_ml_fl = (E_ml * T_ml * delL_lm +        # Elastic response
                   T_ml * V_lm +                   # Viscous response
                   1/g * E_ml_T * delLa_lm +      # Rotational elastic
                   1/g * V_lm_T)                   # Rotational viscous
```

**Physical Interpretation**:
- `E_ml * T_ml * delL_lm`: Instantaneous elastic response
- `T_ml * V_lm`: Time-delayed viscous response from all past loading
- Rotational terms: Earth's rotation axis changes due to loading

### 3.2 Load History Management
```python
# Critical arrays for memory effects
sdelL_lm = np.zeros((n_times-1, len(h_lm)))  # Incremental load history
sdelLa_lm = np.zeros((n_times-1, len(h_lm))) # Rotational load history

# At each time step:
sdelL_lm[t_it-1, :] = delL_lm - delL_lm_prev  # Store load increment
```

### 3.3 Time-Dependent Iteration Structure
```python
# Main time loop
for t_it in range(1, len(ice_time_new)):  # For each time step

    # Update topography and ocean function for this time
    topo_j = topo[:, :, t_it]
    oc_j = calculate_ocean_function(topo_j)

    # Spatial iteration at this time step (Kendall 2005 convergence)
    for k in range(k_max):  # Iterate until sea level converges

        # Calculate total loading (ice + water + sediment)
        delL_lm = calculate_total_loading(...)

        # Apply viscoelastic memory
        V_lm = calculate_viscous_contribution(beta_l, sdelL_lm, t_it, ...)

        # Solve sea level equation with viscous contribution
        delSLcurl_ml_fl = (elastic_terms + viscous_terms)

        # Check convergence and update
        if converged:
            break

    # Store load increment for next time step's memory calculation
    sdelL_lm[t_it-1, :] = delL_lm - delL_lm_prev
```

---

## Phase 4: Rotational Effects Implementation

### 4.1 Rotational Calculation Function
```python
def calc_rot_visc(delL_lm, k_el_2, k_el_tide_2, t_it,
                  beta_konly_l, beta_konly_tide, sdelI, sdelm):
    """
    Port of MATLAB calc_rot_visc function.

    Calculates Earth's rotational response to surface loading:
    - Extracts degree 2 coefficients (L20, L21, L22)
    - Computes polar motion (m1, m2, m3)
    - Calculates rotational potential perturbation
    """
    # Extract degree 2 spherical harmonic coefficients
    L20 = delL_lm[2]  # Axisymmetric loading
    L21 = delL_lm[maxdeg+2]  # m=1 loading
    L22 = delL_lm[2*maxdeg+1]  # m=2 loading

    # Calculate inertia tensor changes (Milne & Mitrovica 1998)
    I13 = sqrt_32_15 * np.pi * a**4 * np.real(L21)
    I23 = -sqrt_32_15 * np.pi * a**4 * np.imag(L21)

    # Polar motion calculation (Mitrovica & Wahr 2005)
    m0 = (I13 + 1j*I23)/CminA * (1 + k_el_2) / (1 - k_el_tide_2/k_hydro)

    # Convert to rotational potential coefficients
    # ... (detailed rotational physics)

    return delLa_lm, sdelI, sdelm
```

### 4.2 Rotational Memory Integration
- **Separate beta arrays**: For rotational Love numbers only
- **Degree 2 only**: Rotation primarily affects ℓ=2 harmonic
- **Coupled dynamics**: Rotation + viscous response interaction

---

## Phase 5: Critical Implementation Details

### 5.1 Index Mapping (MATLAB → Python)
```python
# MATLAB (1-based) → Python (0-based) conversion:
# MATLAB: for t_it = 2:length(ice_time_new)
# Python: for t_it in range(1, len(ice_time_new))

# MATLAB: beta_l{t_it-1}(n-1,:) = [0; beta]
# Python: beta_l[t_it-1][n-1, :] = beta

# Spherical harmonic indexing requires careful conversion
```

### 5.2 Memory Management Strategy
```python
# Efficient storage for large time series
class ViscoelasticSolver:
    def __init__(self, maxdeg, n_times):
        # Pre-allocate major arrays
        self.sdelL_lm = np.zeros((n_times-1, maxdeg))
        self.beta_l = []  # List of variable-size arrays

    def solve_time_step(self, t_it):
        # Process one time step efficiently
        pass
```

### 5.3 Numerical Stability
```python
# Handle exponential calculations carefully
def safe_exponential(spoles, time_diff):
    """Numerically stable exponential decay calculation."""
    # Avoid overflow for large time_diff * spoles
    exp_arg = -spoles * time_diff
    exp_arg = np.clip(exp_arg, -700, 700)  # Prevent overflow
    return 1.0 - np.exp(exp_arg)
```

---

## Phase 6: Validation and Testing

### 6.1 Benchmark Comparison
**Test Cases**:
1. **Simple step loading**: Single ice sheet melting event
2. **Complex history**: Full glacial cycle (ice5g model)
3. **Rotational test**: Large degree-2 loading
4. **Convergence test**: Ensure iteration stability

**Validation Metrics**:
```python
def validate_against_matlab(python_result, matlab_result, tolerance=1e-10):
    """Compare Python implementation against MATLAB benchmark."""
    rel_error = np.abs((python_result - matlab_result) / matlab_result)
    assert np.max(rel_error) < tolerance, f"Max error: {np.max(rel_error)}"
```

### 6.2 Performance Optimization
**Targets**:
- **Beta calculation**: Vectorize over degrees ℓ where possible
- **Memory usage**: Efficient storage of 3D time series
- **Convergence**: Adaptive time stepping for efficiency

---

## Phase 7: User Interface and Documentation

### 7.1 Configuration System
```python
# Easy switching between Earth models and options
config = {
    'earth_model': 'prem.l90C.umVM2.lmVM2',
    'include_rotation': True,
    'include_ice_check': True,
    'max_degree': 128,
    'convergence_tolerance': 1e-4
}
```

### 7.2 Output Formats
```python
class ViscoelasticResults:
    """Structured output matching MATLAB functionality."""
    def __init__(self):
        self.RSL = None          # 3D: Relative sea level (lat, lon, time)
        self.ice_time = None     # 1D: Time array [kyr]
        self.convergence = None  # Convergence history

    def save_netcdf(self, filename):
        """Export to standard format for analysis."""
        pass
```

---

## Implementation Priority Sequence

### **Phase 1 (Critical Foundation)**
1. Create `SL_equation_viscoelastic.py` structure
2. Load viscoelastic parameters from MATLAB file
3. Implement `calculate_beta_coefficients()` function

### **Phase 2 (Core Algorithm)**
4. Add time-dependent infrastructure
5. Implement `calculate_viscous_contribution()`
6. Basic viscoelastic sea level solver

### **Phase 3 (Full Integration)**
7. Time-dependent ocean function updates
8. Load history management
9. Convergence at each time step

### **Phase 4 (Advanced Features)**
10. Rotational effects (`calc_rot_visc()`)
11. Tidal Love numbers
12. Ice correction algorithms

### **Phase 5 (Production Ready)**
13. Validation against MATLAB benchmarks
14. Performance optimization
15. Documentation and examples

---

## Key Success Metrics

### **Mathematical Correctness**
- [ ] Beta coefficients match MATLAB values (< 1e-12 relative error)
- [ ] Exponential decay curves reproduce exactly
- [ ] Sea level fingerprints agree with MATLAB benchmarks

### **Performance Targets**
- [ ] Runtime comparable to MATLAB (within 2x)
- [ ] Memory usage reasonable for typical problems
- [ ] Convergence stability maintained

### **Usability Goals**
- [ ] Drop-in replacement for elastic solver
- [ ] Clear documentation with theory background
- [ ] Standard output formats for analysis

---

## Critical Dependencies and Risks

### **Dependencies**
- **NumPy**: For efficient array operations
- **SciPy**: For MATLAB file I/O and interpolation
- **Existing elastic solver**: As foundation

### **Main Implementation Risks**
1. **Index mapping errors**: MATLAB 1-based vs Python 0-based
2. **Memory requirements**: Large 3D time series arrays
3. **Numerical precision**: Exponential calculations with small arguments
4. **Convergence issues**: Time step coupling may affect stability

### **Mitigation Strategies**
- **Incremental testing**: Validate each component against MATLAB
- **Memory profiling**: Monitor array sizes during development
- **Numerical testing**: Test exponential functions across parameter ranges
- **Convergence analysis**: Compare iteration behavior with MATLAB

---

## Final Notes for Implementation

**The core insight**: The exponential terms `exp(-s_{ℓ,k} * t)` in viscoelastic Love numbers arise directly from the normal-mode eigenvalue spectrum of the Earth's viscoelastic response. This mathematical connection from Tromp (1999) to Kendall (2005) to SLcode implementation is the theoretical thread that must be preserved exactly.

**Critical equation to implement correctly**:
```
β_ℓ(t₁→t₂) = Σ_k [(k_amp_{ℓ,k} - h_amp_{ℓ,k})/s_{ℓ,k}] * [1 - exp(-s_{ℓ,k} * (t₂-t₁))]
```

This single formula encodes the complete viscoelastic memory of the Earth and is the mathematical bridge between fundamental continuum mechanics and practical sea level calculations.

---

**Status**: Plan completed and approved for implementation
**Next Action**: Begin Phase 1 implementation
**Expected Completion**: Phased development over multiple sessions
