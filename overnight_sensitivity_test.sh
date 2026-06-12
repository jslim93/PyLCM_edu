#!/bin/bash
# Overnight Sensitivity Test: PyLCM vs Fortran (Light)
# Fortran 10k×2 + PyLCM 10k×2, running in parallel
# Estimated wall time: ~50 min

set -e
LOGFILE="/tmp/overnight_sensitivity_test.log"
NOTIFY="$HOME/Tools/claude-mobile-alerts/notify_discord.sh"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"
}

notify() {
    if [ -x "$NOTIFY" ]; then
        "$NOTIFY" "$1" 2>/dev/null || true
    fi
}

# Initialize
echo "" > "$LOGFILE"
log "=== Sensitivity Test Started ==="
log "Plan: Fortran 10k×2 + PyLCM 10k×2 (parallel)"
notify "테스트 시작: Fortran 10k×2 + PyLCM 10k×2 동시 실행"

TEST_START=$(date +%s)

# ============================================================
# Fortran 10k × 2 ensemble (background, ~50 min)
# ============================================================
log "--- Starting Fortran 10k × 2 ensemble (background) ---"

# Backup current input.nml
cp /Users/dr.cloud/particle_model/input.nml /Users/dr.cloud/particle_model/input.nml.bak_overnight

cat > /Users/dr.cloud/particle_model/input.nml << 'EOF'
&INPUT
    file_ending = "compare_10k",
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
    n_ensemble = 2,
    switch_batch_mode = .TRUE.,
    tot_ensemble = 2,
/
EOF

log "Fortran input.nml created (10k particles, 2 ensemble members)"

# Run Fortran in background
FORTRAN_START=$(date +%s)
(
    cd /Users/dr.cloud/particle_model/objects
    ./particle_model.o > /tmp/fortran_run.log 2>&1
    FORTRAN_END=$(date +%s)
    FORTRAN_ELAPSED=$(( (FORTRAN_END - FORTRAN_START) / 60 ))
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fortran complete! Elapsed: ${FORTRAN_ELAPSED} minutes" >> "$LOGFILE"
) &
FORTRAN_PID=$!
log "Fortran running in background (PID: $FORTRAN_PID)"

# ============================================================
# PyLCM 10k × 2 (parallel, ~3.5 min)
# ============================================================
log "--- Starting PyLCM 10k × 2 (parallel) ---"

cd /Users/dr.cloud/PyLCM_edu
source ~/.zshrc 2>/dev/null || true
conda activate base 2>/dev/null || true

PYLCM_START=$(date +%s)

# Member 0 (seed=100)
python3 -c "
import sys; sys.path.insert(0,'.')
from ensemble_comparison import run_single
import json
print('[PyLCM Member 0] Starting with seed=100, n=10000')
results, elapsed = run_single(seed=100, n_ptcl=10000)
print(f'[PyLCM Member 0] Done in {elapsed:.1f}s')
out = {}
for t, r in results.items():
    out[str(t)] = {k: float(v) for k, v in r.items()}
out['elapsed'] = elapsed
out['seed'] = 100
out['n_ptcl'] = 10000
with open('/tmp/pylcm_10k_ens_00.json', 'w') as f:
    json.dump(out, f)
print('[PyLCM Member 0] Saved')
" >> "$LOGFILE" 2>&1 &
PY_PID0=$!

# Member 1 (seed=200)
python3 -c "
import sys; sys.path.insert(0,'.')
from ensemble_comparison import run_single
import json
print('[PyLCM Member 1] Starting with seed=200, n=10000')
results, elapsed = run_single(seed=200, n_ptcl=10000)
print(f'[PyLCM Member 1] Done in {elapsed:.1f}s')
out = {}
for t, r in results.items():
    out[str(t)] = {k: float(v) for k, v in r.items()}
out['elapsed'] = elapsed
out['seed'] = 200
out['n_ptcl'] = 10000
with open('/tmp/pylcm_10k_ens_01.json', 'w') as f:
    json.dump(out, f)
print('[PyLCM Member 1] Saved')
" >> "$LOGFILE" 2>&1 &
PY_PID1=$!

log "PyLCM running in background (PIDs: $PY_PID0, $PY_PID1)"

# Wait for PyLCM first (finishes much faster)
wait $PY_PID0 $PY_PID1
PYLCM_END=$(date +%s)
PYLCM_ELAPSED=$(( (PYLCM_END - PYLCM_START) ))
log "PyLCM complete! Elapsed: ${PYLCM_ELAPSED} seconds"

# Wait for Fortran
log "Waiting for Fortran to finish..."
wait $FORTRAN_PID
FORTRAN_STATUS=$?
log "Fortran process exited with status: $FORTRAN_STATUS"

# ============================================================
# Parse and Compare Results
# ============================================================
log "--- Comparison ---"

python3 << 'PYEOF' 2>&1 | tee -a "$LOGFILE"
import numpy as np, json, os

# --- Parse Fortran 2-member ensemble ---
fortran_file = '/Users/dr.cloud/particle_model/output/time_series_compare_10k.dat'
if not os.path.exists(fortran_file):
    print(f"ERROR: Fortran output not found: {fortran_file}")
    exit(1)

all_lines = []
with open(fortran_file) as f:
    all_lines = f.readlines()

blocks = []
current_block = []
for line in all_lines:
    stripped = line.strip()
    if not stripped:
        continue
    try:
        vals = [float(v) for v in stripped.split()]
        if len(vals) > 40:
            current_block.append(vals)
    except ValueError:
        if current_block:
            blocks.append(np.array(current_block))
            current_block = []
if current_block:
    blocks.append(np.array(current_block))

print(f"\n{'='*70}")
print(f"  Fortran 10k: {len(blocks)} ensemble member(s)")
print(f"{'='*70}")

target_times = [500, 1000, 1500, 2000, 2500, 3000]
fortran_results = {t: {'T':[], 'qc':[], 'qr':[], 'Nc':[], 'Nr':[]} for t in target_times}

for b_idx, data in enumerate(blocks):
    for target_t in target_times:
        idx = np.argmin(np.abs(data[:,0] - target_t))
        rho = data[idx, 48]
        Tc = data[idx, 4]
        QC = data[idx, 10] * 1e3
        QR = data[idx, 12] * 1e3
        NC = data[idx, 9] * rho / 1e6
        NR = data[idx, 11] * rho / 1e6
        fortran_results[target_t]['T'].append(Tc)
        fortran_results[target_t]['qc'].append(QC)
        fortran_results[target_t]['qr'].append(QR)
        fortran_results[target_t]['Nc'].append(NC)
        fortran_results[target_t]['Nr'].append(NR)

print(f"\n{'t':>5s} | {'T(C)':>12s} | {'qc(g/kg)':>14s} | {'qr(g/kg)':>14s} | {'Nc(/cm3)':>14s} | {'Nr(/cm3)':>14s}")
print("-" * 85)
for t in target_times:
    r = fortran_results[t]
    T_m, T_s = np.mean(r['T']), np.std(r['T'])
    qc_m, qc_s = np.mean(r['qc']), np.std(r['qc'])
    qr_m, qr_s = np.mean(r['qr']), np.std(r['qr'])
    Nc_m, Nc_s = np.mean(r['Nc']), np.std(r['Nc'])
    Nr_m, Nr_s = np.mean(r['Nr']), np.std(r['Nr'])
    print(f"{t:5d} | {T_m:6.2f}+/-{T_s:4.2f} | {qc_m:7.4f}+/-{qc_s:5.4f} | {qr_m:7.4f}+/-{qr_s:5.4f} | {Nc_m:7.2f}+/-{Nc_s:5.2f} | {Nr_m:7.2f}+/-{Nr_s:5.2f}")

# --- Load PyLCM 10k × 2 ---
pylcm_results = {t: {'T':[], 'qc':[], 'qr':[], 'Nc':[], 'Nr':[]} for t in target_times}
pylcm_count = 0
for i in range(2):
    fpath = f'/tmp/pylcm_10k_ens_{i:02d}.json'
    if os.path.exists(fpath):
        with open(fpath) as f:
            d = json.load(f)
        pylcm_count += 1
        for t in target_times:
            ts = str(t)
            if ts in d:
                pylcm_results[t]['T'].append(d[ts]['T'])
                pylcm_results[t]['qc'].append(d[ts]['qc'])
                pylcm_results[t]['qr'].append(d[ts]['qr'])
                pylcm_results[t]['Nc'].append(d[ts]['Nc'])
                pylcm_results[t]['Nr'].append(d[ts]['Nr'])

print(f"\n{'='*70}")
print(f"  PyLCM 10k: {pylcm_count} ensemble member(s)")
print(f"{'='*70}")
print(f"\n{'t':>5s} | {'T(C)':>12s} | {'qc(g/kg)':>14s} | {'qr(g/kg)':>14s} | {'Nc(/cm3)':>14s} | {'Nr(/cm3)':>14s}")
print("-" * 85)
for t in target_times:
    r = pylcm_results[t]
    if r['T']:
        T_m, T_s = np.mean(r['T']), np.std(r['T'])
        qc_m, qc_s = np.mean(r['qc']), np.std(r['qc'])
        qr_m, qr_s = np.mean(r['qr']), np.std(r['qr'])
        Nc_m, Nc_s = np.mean(r['Nc']), np.std(r['Nc'])
        Nr_m, Nr_s = np.mean(r['Nr']), np.std(r['Nr'])
        print(f"{t:5d} | {T_m:6.2f}+/-{T_s:4.2f} | {qc_m:7.4f}+/-{qc_s:5.4f} | {qr_m:7.4f}+/-{qr_s:5.4f} | {Nc_m:7.2f}+/-{Nc_s:5.2f} | {Nr_m:7.2f}+/-{Nr_s:5.2f}")

# --- Difference table ---
print(f"\n{'='*70}")
print(f"  Difference: PyLCM - Fortran (mean)")
print(f"{'='*70}")
print(f"\n{'t':>5s} | {'dT(C)':>8s} | {'dqc':>10s} | {'dqr':>10s} | {'dNc':>10s} | {'dNr':>10s} | {'dqt%':>8s}")
print("-" * 75)
for t in target_times:
    fr = fortran_results[t]
    pr = pylcm_results[t]
    if pr['T']:
        dT = np.mean(pr['T']) - np.mean(fr['T'])
        dqc = np.mean(pr['qc']) - np.mean(fr['qc'])
        dqr = np.mean(pr['qr']) - np.mean(fr['qr'])
        dNc = np.mean(pr['Nc']) - np.mean(fr['Nc'])
        dNr = np.mean(pr['Nr']) - np.mean(fr['Nr'])
        qt_f = np.mean(fr['qc']) + np.mean(fr['qr'])
        qt_p = np.mean(pr['qc']) + np.mean(pr['qr'])
        dqt_pct = (qt_p - qt_f) / qt_f * 100 if qt_f > 0.001 else 0
        print(f"{t:5d} | {dT:+8.3f} | {dqc:+10.4f} | {dqr:+10.4f} | {dNc:+10.2f} | {dNr:+10.2f} | {dqt_pct:+8.2f}%")

# --- Also compare with PyLCM 100k if available ---
summary_100k = '/tmp/pylcm_100k_ensemble_summary.json'
if os.path.exists(summary_100k):
    with open(summary_100k) as f:
        py100k = json.load(f)
    print(f"\n{'='*70}")
    print(f"  Reference: PyLCM 100k (5-member ensemble, previous run)")
    print(f"{'='*70}")
    print(f"\n{'t':>5s} | {'T(C)':>12s} | {'qc(g/kg)':>14s} | {'qr(g/kg)':>14s} | {'Nc(/cm3)':>14s} | {'Nr(/cm3)':>14s}")
    print("-" * 85)
    for t in target_times:
        ts = str(t)
        if ts in py100k:
            p = py100k[ts]
            print(f"{t:5d} | {p['T']['mean']:6.2f}+/-{p['T']['std']:4.2f} | "
                  f"{p['qc']['mean']:7.4f}+/-{p['qc']['std']:5.4f} | "
                  f"{p['qr']['mean']:7.4f}+/-{p['qr']['std']:5.4f} | "
                  f"{p['Nc']['mean']:7.2f}+/-{p['Nc']['std']:5.2f} | "
                  f"{p['Nr']['mean']:7.2f}+/-{p['Nr']['std']:5.2f}")

# Save consolidated results
output = {
    'fortran_10k_2ens': {},
    'pylcm_10k_2ens': {},
}
for t in target_times:
    fr = fortran_results[t]
    pr = pylcm_results[t]
    output['fortran_10k_2ens'][str(t)] = {
        k: {'mean': float(np.mean(v)), 'std': float(np.std(v)), 'values': [float(x) for x in v]}
        for k, v in fr.items()
    }
    if pr['T']:
        output['pylcm_10k_2ens'][str(t)] = {
            k: {'mean': float(np.mean(v)), 'std': float(np.std(v)), 'values': [float(x) for x in v]}
            for k, v in pr.items()
        }

with open('/tmp/sensitivity_comparison_10k_2ens.json', 'w') as f:
    json.dump(output, f, indent=2)
print(f"\nResults saved to /tmp/sensitivity_comparison_10k_2ens.json")
PYEOF

TEST_END=$(date +%s)
TOTAL_ELAPSED=$(( (TEST_END - TEST_START) / 60 ))
log "=== All complete! Total elapsed: ${TOTAL_ELAPSED} minutes ==="

notify "테스트 완료! (${TOTAL_ELAPSED}분) 결과: /tmp/sensitivity_comparison_10k_2ens.json 로그: $LOGFILE"
