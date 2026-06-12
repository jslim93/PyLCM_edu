# Sensitivity Test Handoff: PyLCM vs Fortran Comparison

## Objective
Compare PyLCM (Python) vs Fortran particle_model with identical physics settings.
Condensation + Collision (Hall kernel), 3000 timesteps, dt=1s, w=1m/s.

## What's Already Done

### PyLCM: 100k particles, 5-member ensemble ✅
- Results saved in `/tmp/pylcm_ensemble_00.json` through `_04.json`
- Consolidated summary: `/tmp/pylcm_100k_ensemble_summary.json`
- Script: `/Users/dr.cloud/PyLCM_edu/ensemble_comparison.py`
- Each run took ~36 min, 5 runs in parallel ~36 min total
- **Variance is <1%** — 100k particles is well converged

#### PyLCM 100k Results (mean ± std, 5 members)

| t (s) | qc (g/kg) | qr (g/kg) | qt (g/kg) | Nc (/cm³) | Nr (/cm³) |
|-------|-----------|-----------|-----------|-----------|-----------|
| 500   | 0.505±0.000 | 0.001±0.000 | 0.506 | 55.7±0.2 | 0.01 |
| 1000  | 1.525±0.006 | 0.035±0.006 | 1.560 | 54.3±0.3 | 0.10±0.01 |
| 1500  | 1.840±0.048 | 0.744±0.047 | 2.585 | 37.4±1.3 | 0.43±0.02 |
| 2000  | 0.000±0.000 | 3.576±0.001 | 3.576 | 0.3±0.2  | 15.2±0.5 |
| 2500  | 0.062±0.002 | 4.474±0.002 | 4.535 | 8.1±0.2  | 5.6±0.1 |
| 3000  | 0.134±0.004 | 5.321±0.004 | 5.454 | 6.9±0.1  | 1.7±0.0 |

### Fortran: 1k particles, single run ✅
- Results in `/Users/dr.cloud/particle_model/output/time_series_fair_coll.dat`
- Took 159s (53ms/step)
- Good qualitative agreement with PyLCM but noisy (Monte Carlo variance)

### Fortran: 100k particles ❌ Too Slow
- Estimated ~15 hours per run — impractical
- Killed after test showed only t=30s in 9 minutes of wall time

---

## What Needs to Be Done

### Step 1: Fortran 10k particles (Quick Comparison)

**Estimated time: ~25 min per run, single run sufficient for test**

```bash
# 1. Create input file
cat > /Users/dr.cloud/particle_model/input.nml << 'EOF'
&INPUT
    file_ending = "test_10k",
    dt = 1.0,
    dt_physics_max = 1.0,
    time_end = 3000.0,
!
    type_init = 'USERDEFINED',
    z_init = 0.0,
    T_init      = 293.2,
    rh_init     = 0.88,
    p_init      = 1013.0E2,
!
    switch_diffusion_micro    = .TRUE.,
    switch_koehler_micro      = .TRUE.,
    switch_radiation_micro    = .FALSE.,
    switch_collection_micro   = .TRUE.,
    switch_coll_breakup_micro = .FALSE.,
    switch_ice_micro          = .FALSE.,
    switch_FFC_micro          = .FALSE.,
    switch_direct_micro       = .FALSE.,
    switch_sed_removal_micro  = .FALSE.,
    switch_ventilation_micro  = .FALSE.,
!
    init_micro           = 'aerosols',
    n_particles = 10000,
    collection_kernel_micro = 'hall',
!
    n_aero = 118.0E6, 11.0E6, 0.72E6, 0.0,
    rm_aero = 0.019E-6, 0.056E-6, 0.46E-6, 0.0,
    sigma_aero = 3.3, 1.6, 2.2, 0.0,
    weight_aero = 0.91, 0.085, 0.005, 0.0,
    units_aero = 'kg-1',
!
    switch_parcel = .TRUE.,
    w_mean_parcel = 1.0,
    switch_wave_parcel = .FALSE.,
!
    switch_entrainment = .FALSE.,
    switch_LEM = .FALSE.,
!
    n_spec = 10,
    alpha_spec = 1.0,
    r_start_spec = 0.01E-6,
    n_ts = 1,
    switch_output_1D_LEM = .FALSE.,
    n_1D_LEM = 100,
!
    n_ensemble = 1,
    switch_batch_mode = .TRUE.,
    tot_ensemble = 1,
/
EOF

# 2. Run Fortran
cd /Users/dr.cloud/particle_model/objects
time ./particle_model.o

# 3. Output will be in:
#    /Users/dr.cloud/particle_model/output/time_series_test_10k.dat
```

### Step 2: Parse Fortran Output and Compare

```bash
cd /Users/dr.cloud/PyLCM_edu
python3 << 'PYEOF'
import numpy as np, json

# --- Parse Fortran ---
# Column layout (0-indexed):
# 0:TIME 1:Z 2:W 3:P 4:T(°C) 5:QV 6:S 7:NA(/kg) 8:QA 9:NC(/kg) 10:QC(kg/kg)
# 11:NR(/kg) 12:QR(kg/kg) ... 47:VOL(m³) 48:RHO(kg/m³)
data = []
with open('/Users/dr.cloud/particle_model/output/time_series_test_10k.dat') as f:
    for line in f.readlines()[3:]:  # skip 3 header lines
        vals = line.split()
        if len(vals) > 48:
            try:
                data.append([float(v) for v in vals])
            except:
                pass
data = np.array(data)

print("=== Fortran 10k particles ===")
print(f"{'t':>5s} {'T(C)':>7s} {'qc':>8s} {'qr':>8s} {'qt':>8s} {'Nc':>8s} {'Nr':>8s}")
for target_t in [500, 1000, 1500, 2000, 2500, 3000]:
    idx = np.argmin(np.abs(data[:,0] - target_t))
    rho = data[idx, 48]
    Tc = data[idx, 4]
    QC = data[idx, 10] * 1e3  # kg/kg -> g/kg
    QR = data[idx, 12] * 1e3
    NC = data[idx, 9] * rho / 1e6   # /kg -> /cm³
    NR = data[idx, 11] * rho / 1e6
    print(f'{target_t:5d} {Tc:7.2f} {QC:8.4f} {QR:8.4f} {QC+QR:8.4f} {NC:8.2f} {NR:8.2f}')

# --- Load PyLCM results ---
with open('/tmp/pylcm_100k_ensemble_summary.json') as f:
    py = json.load(f)

print("\n=== PyLCM 100k particles (ensemble mean) ===")
print(f"{'t':>5s} {'T(C)':>7s} {'qc':>8s} {'qr':>8s} {'qt':>8s} {'Nc':>8s} {'Nr':>8s}")
for t in ['500','1000','1500','2000','2500','3000']:
    p = py[t]
    qt = p['qc']['mean'] + p['qr']['mean']
    print(f"{t:>5s} {p['T']['mean']:7.2f} {p['qc']['mean']:8.4f} {p['qr']['mean']:8.4f} {qt:8.4f} {p['Nc']['mean']:8.2f} {p['Nr']['mean']:8.2f}")
PYEOF
```

### Step 3 (Optional): PyLCM 10k for Apples-to-Apples

If you want PyLCM with same 10k particles for direct comparison:

```bash
cd /Users/dr.cloud/PyLCM_edu
python3 ensemble_comparison.py 42 0   # ~3.5 min with 10k particles
# Edit ensemble_comparison.py line: n_ptcl=10000 in run_single() default
```

Or inline:
```bash
python3 -c "
import sys; sys.path.insert(0,'.')
from ensemble_comparison import run_single
results, elapsed = run_single(seed=42, n_ptcl=10000)
for t in sorted(results.keys()):
    r = results[t]
    print(f't={t:4d} T={r[\"T\"]:6.2f}C qc={r[\"qc\"]:7.4f} qr={r[\"qr\"]:7.4f} Nc={r[\"Nc\"]:7.2f} Nr={r[\"Nr\"]:7.2f}')
print(f'Elapsed: {elapsed:.0f}s')
"
```

---

## Production Run Plan (10-member Ensemble)

When ready for production:

### Fortran: 10k × 10 ensemble
```bash
# Modify input.nml:
#   n_particles = 10000
#   n_ensemble = 10
#   tot_ensemble = 10
#   file_ending = "prod_10k"
cd /Users/dr.cloud/particle_model/objects
time ./particle_model.o
# Estimated: ~250 min (4 hours)
# Output: time_series_prod_10k.dat (has 10 blocks, one per ensemble member)
```

### PyLCM: 10k × 10 ensemble (parallel)
```bash
cd /Users/dr.cloud/PyLCM_edu
# Run 5 at a time (edit n_ptcl=10000 in ensemble_comparison.py first)
for i in $(seq 0 4); do
    python3 ensemble_comparison.py $((100*(i+1))) $i > /tmp/pylcm_ens_$i.log 2>&1 &
done
wait
for i in $(seq 5 9); do
    python3 ensemble_comparison.py $((100*(i+1))) $i > /tmp/pylcm_ens_$i.log 2>&1 &
done
wait
# Estimated: ~7 min total (2 batches × 3.5 min)
```

---

## Key Files

| File | Description |
|------|-------------|
| `/Users/dr.cloud/PyLCM_edu/ensemble_comparison.py` | PyLCM ensemble runner script |
| `/Users/dr.cloud/PyLCM_edu/SENSITIVITY_TEST_HANDOFF.md` | This document |
| `/tmp/pylcm_ensemble_00.json` ... `_04.json` | PyLCM 100k results (5 members) |
| `/tmp/pylcm_100k_ensemble_summary.json` | Consolidated PyLCM 100k stats |
| `/Users/dr.cloud/particle_model/source/make_local` | Fortran makefile |
| `/Users/dr.cloud/particle_model/objects/particle_model.o` | Fortran executable |
| `/Users/dr.cloud/particle_model/input.nml` | Fortran input (currently set to ensemble config) |
| `/Users/dr.cloud/particle_model/input.nml.bak_fair_coll` | Backup of 1k collision input |

## Physics Settings (Both Models)

- T_init = 293.2 K, P_init = 101300 Pa, RH_init = 0.88, w = 1.0 m/s
- Aerosol: 3 lognormal modes (118, 11, 0.72 cm⁻³)
- dt = 1.0 s, dt_physics_max = 1.0 s (no adaptive substep for fair comparison)
- Hall collection kernel, no breakup, no sedimentation removal
- z_max = 3000 m (parcel stops ascending at this height)
- Linear ascending mode, Stable environmental profile (θ lapse = +5 K/km)

## Important Notes

1. **PyLCM aero_init convention**: mu_aero and sigma_aero must be **log-transformed** before passing to `aero_init()`. The `model_init()` function handles this automatically (lines 39, 59-60 of aero_init.py). The `ensemble_comparison.py` script does this correctly.

2. **Fortran runs from objects/ dir**: The executable expects `../input.nml` relative to its location. Always `cd /Users/dr.cloud/particle_model/objects` before running.

3. **Fortran output columns**: TIME(0) Z(1) W(2) P(3) T_celsius(4) QV(5) S(6) NA_per_kg(7) QA(8) NC_per_kg(9) QC_kg_per_kg(10) NR_per_kg(11) QR_kg_per_kg(12) ... VOL_m3(47) RHO_kg_per_m3(48). Convert N: N_cm3 = N_per_kg × rho / 1e6. Convert q: q_g_per_kg = q_kg_per_kg × 1e3.

4. **PyLCM ts_analysis output**: qc, qr already in g/kg; Nc, Nr already in /cm³. Do NOT multiply by 1e3 again.

5. **Restore Fortran input.nml**: Current input.nml has ensemble config with 100k particles. Restore with: `cp /Users/dr.cloud/particle_model/input.nml.bak_fair_coll /Users/dr.cloud/particle_model/input.nml`
