# LEM Implementation Notes

## Overview

Linear Eddy Model (LEM) for Subgrid-Scale (SGS) mixing in PyLCM box model.
Based on Krueger (1993, JAS) and SAM LCM research code.

## Research Code vs Box Model Implementation

### Research Code (SAM LCM)
```
Location: ~/SAM6.10.10.LCM_JS/SRC/MICRO_LAGRANGE/micro_sgs_mixing.f90
```

**Key features:**
1. Each particle has `ieta` (absolute supersaturation), `iT`, `iq`
2. LES provides `tk_LCM` (turbulent diffusivity), `smix_LCM` (mixing length)
3. Particles within each grid box are arranged in 1D LEM array
4. Diffusion + triplet map rearrangement applied
5. Nudging to grid box mean with τ = 900s

**Per-particle supersaturation (line 274-292):**
```fortran
eta_dyn = particle_var(n, iq) - qsatw_val  ! particle's own q
eta_lst = particle_var(n, ieta)            ! particle's previous eta
! Tzivion method applied per-particle
```

### Box Model Implementation (PyLCM)
```
Location: PyLCM/sgs_mixing.py, PyLCM/condensation_lem.py
```

**Adaptations for box model:**
1. No LES grid → single parcel treated as one "grid box"
2. No LES turbulence → estimate from dissipation rate (eps)
3. T_lem, q_lem stored as perturbations from parcel mean
4. Linearized S' = q'/qsat - (L/RvT²) * T'

**Key equations:**
```python
# Supersaturation perturbation
S_prime = q_prime / qsatw - dqsat_dT / qsatw * T_prime
supersat_local = supersat_mean + S_prime
```

## Physics Verification

### S_std Validation
| Source | S_std |
|--------|-------|
| Theory (T_std=0.1K) | 0.68% |
| Our implementation | 0.53% |
| Siebert et al. (2015) | 0.1-1% |
| Grabowski & Abade (2017) | 0.5-2% |

**Conclusion:** Our S_std is physically reasonable.

### Expected LEM Effects
1. ✅ S fluctuations (S_std > 0)
2. ✅ Broader droplet size distribution
3. ✅ Mean S unchanged
4. ✅ Fluctuations decay over time (diffusion + nudging)

## Differences from Research Code

| Aspect | Research Code | Box Model |
|--------|---------------|-----------|
| Turbulence source | LES SGS model | User-specified eps |
| Grid structure | 3D LES grid | Single parcel |
| eta tracking | Absolute (kg/kg) | Perturbation S' |
| T, q tracking | Absolute values | Perturbations from mean |
| Entrainment | From adjacent cells | Not applicable |

## Parameters

```python
# Default values
diss_rate = 1e-4      # m²/s³ (typical cloud value)
L_domain = 100.0      # m (LEM domain size)
tau_nudge = 900.0     # s (nudging timescale)
T_fluct_std = 0.1     # K (initial T fluctuation)
q_fluct_std = 1e-5    # kg/kg (initial q fluctuation)
```

## Usage

```python
from PyLCM.condensation_lem import drop_condensation_lem

# In time loop:
particles_list, T, q, S_lst, ... = drop_condensation_lem(
    particles_list, T_parcel, q_parcel, P_parcel,
    nt, dt, air_mass, S_lst, rho_aero,
    kohler_activation_radius=True,
    con_ts, act_ts, evp_ts, dea_ts,
    switch_kappa_koehler=True,
    diss_rate=1e-4,    # Turbulent dissipation rate
    L_domain=100.0     # LEM domain size
)
```

## Future Improvements

1. Add explicit eta tracking (like research code)
2. Better turbulence parameterization
3. Time-varying fluctuation source
4. Breakup effects on LEM array
