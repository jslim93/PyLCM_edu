# Handoff Document: Convergence Comparison Benchmark

## What We Did

### 1. Added Separate Seed Control for Initialization vs Collision

**Files Modified:**
- `PyLCM/timestep_routine_arrays.py` (lines 118-129, 250-257)
- `PyLCM/ensemble.py` (lines 196-215, 218-271, 313-340)

**Changes:**
- Added `collision_seed` parameter to `timesteps_function_arrays()`
- After particle initialization completes, the RNG is reseeded with `collision_seed` if provided
- This separates initialization randomness from collision randomness

**Rationale:**
- User requested that initialization be deterministic for better convergence analysis
- Stochasticity should only come from the collision process
- This allows fair comparison between Weighting_factor and Random initialization modes

### 2. Updated Ensemble Module

**Changes to `run_ensemble()` and `run_ensemble_serial()`:**
- Added `init_seed` parameter for deterministic initialization across all ensemble members
- `seed_base` now controls collision seeds only (member i gets `collision_seed = seed_base + i`)
- Backward compatible: if `init_seed=None`, uses original behavior (same seed for both)

**New signature:**
```python
run_ensemble(n_members, seed_base=0, init_seed=None, n_workers=None, method='fork', **kwargs)
run_ensemble_serial(n_members, seed_base=0, init_seed=None, **kwargs)
```

### 3. Created Convergence Comparison Benchmark

**File Created:** `benchmark_convergence.py`

**Purpose:** Compare convergence of two superdroplet simulation approaches:
1. **Weighting_factor** - Deterministic PDF sampling for aerosol sizes
2. **Random+Ensemble** - Random sampling with ensemble averaging

**Test Configuration:**
- Simulation: 2700s, RH=0.98, w=2.0 m/s, collision enabled
- Particle counts: 1000, 5000, 10000, 20000
- 5 trials per configuration
- 10 ensemble members for Random approach
- Fixed init_seed=42 for all runs

**Metrics:**
- Coefficient of Variation (CV) = std/mean × 100%
- Variables tracked: qc, qr, nc, nr
- Target: CV < 5% indicates convergence

---

## What To Do Next

### 1. Run the Full Benchmark
```bash
python benchmark_convergence.py
```
- Estimated runtime: ~35 minutes
- Outputs:
  - `convergence_comparison.png` (visualization)
  - `convergence_results.npz` (raw data)

### 2. Analyze Results
- Compare CV values between Weighting_factor and Random+Ensemble
- Check if mean values are consistent (no systematic bias)
- Identify the minimum particle count needed for CV < 5%

### 3. Potential Follow-up Analysis
- Test with different ensemble sizes (5, 10, 20 members)
- Test with collision disabled (isolate condensation convergence)
- Compare computational cost vs convergence quality

---

## What To Be Careful About

### 1. Seed Behavior
- `seed` parameter controls **initialization only** (aerosol sampling in Random mode)
- `collision_seed` parameter controls **collision process only**
- If `collision_seed=None`, collision continues with same RNG state from initialization
- For reproducible collision results, always specify `collision_seed`

### 2. Weighting_factor Mode
- Weighting_factor initialization is **already deterministic** (no random sampling)
- The `seed` parameter still matters for setting initial RNG state
- Variability in Weighting_factor runs comes **only from collision**

### 3. Ensemble with Deterministic Init
- When using `init_seed` in ensemble runs, all members start with **identical particles**
- This is intentional for isolating collision stochasticity
- For traditional ensemble (varying initialization), set `init_seed=None`

### 4. Short Simulations
- Rain (qr) requires sufficient simulation time to form via collision
- At t < 500s, qr may be zero, causing CV=NaN
- Use full 2700s simulation for meaningful rain statistics

### 5. Backward Compatibility
- Old code using `seed=X` without `collision_seed` still works
- Old ensemble code without `init_seed` uses original behavior
- New behavior only activates when new parameters are explicitly set

---

## File Summary

| File | Status | Description |
|------|--------|-------------|
| `PyLCM/timestep_routine_arrays.py` | Modified | Added `collision_seed` parameter |
| `PyLCM/ensemble.py` | Modified | Added `init_seed` parameter |
| `benchmark_convergence.py` | Created | Convergence comparison benchmark |

---

## Quick Test Command

To verify everything works:
```bash
python -c "
from PyLCM.timestep_routine_arrays import timesteps_function_arrays
from PyLCM.parameters import rm_spec
import numpy as np

# Test with separate seeds
output = timesteps_function_arrays(
    n_particles=500, P_parcel=100000, RH_parcel=0.98, T_parcel=285,
    w_parcel=2.0, nt=50, dt=1.0, rm_spec=rm_spec,
    do_collision=True, verbose=False,
    seed=42,           # Fixed initialization
    collision_seed=0   # Collision seed
)
print(f'qc = {output[8][-1]*1e3:.4f} g/kg')
print('Test passed!')
"
```
