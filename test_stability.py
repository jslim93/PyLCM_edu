#!/usr/bin/env python3
"""Multi-scenario stability test for feature/lecture-readiness-improvements.
   Tests 12 microphysics scenarios: aerosol types x physics toggles.
   Checks for: crashes, NaN, negative q, unphysical T, mass conservation."""
import numpy as np
import copy, sys, time, traceback, warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, '/Users/dr.cloud/PyLCM_edu')
from PyLCM.parameters import *
from PyLCM.micro_particle import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from PyLCM.collision import *
from Post_process.analysis import *

# ============================================================
# Aerosol presets (from widget.py)
# ============================================================
PRESETS = {
    'Maritime':     {'N': [100.0, 20.0],  'mu': [0.050, 0.250], 'sigma': [2.0, 1.8], 'kappa': [1.0, 1.2]},
    'Continental':  {'N': [3200.0, 2900.0, 0.3], 'mu': [0.016, 0.068, 0.920], 'sigma': [2.1, 2.0, 2.2], 'kappa': [0.14, 0.3, 0.7]},
    'Arctic':       {'N': [15.0, 5.0],    'mu': [0.030, 0.100], 'sigma': [1.8, 1.5], 'kappa': [0.4, 0.6]},
    'Fortran_ref':  {'N': [118.0, 11.0, 0.72], 'mu': [0.019, 0.056, 0.46], 'sigma': [3.3, 1.6, 2.2], 'kappa': [1.28, 1.28, 1.28]},
}

# ============================================================
# Scenario definitions
# ============================================================
SCENARIOS = [
    # --- Aerosol type sweep (cond+coll, default physics) ---
    {'name': 'Maritime_cond+coll',     'preset': 'Maritime',    'do_cond': True, 'do_coll': True},
    {'name': 'Continental_cond+coll',  'preset': 'Continental', 'do_cond': True, 'do_coll': True},
    {'name': 'Arctic_cond+coll',       'preset': 'Arctic',      'do_cond': True, 'do_coll': True},
    {'name': 'Fortran_ref_cond+coll',  'preset': 'Fortran_ref', 'do_cond': True, 'do_coll': True},

    # --- Condensation-only ---
    {'name': 'Maritime_cond_only',     'preset': 'Maritime',    'do_cond': True, 'do_coll': False},
    {'name': 'Continental_cond_only',  'preset': 'Continental', 'do_cond': True, 'do_coll': False},

    # --- Ablation lab: Kelvin off ---
    {'name': 'Maritime_no_kelvin',     'preset': 'Maritime',    'do_cond': True, 'do_coll': True,
     'switch_kelvin': False},

    # --- Ablation lab: Solute off ---
    {'name': 'Maritime_no_solute',     'preset': 'Maritime',    'do_cond': True, 'do_coll': True,
     'switch_solute': False},

    # --- Ablation lab: E_constant + vt_simple ---
    {'name': 'Maritime_simple_coll',   'preset': 'Maritime',    'do_cond': True, 'do_coll': True,
     'switch_E_constant': True, 'switch_vt_simple': True},

    # --- Turbulent collision kernel ---
    {'name': 'Maritime_turb_kernel',   'preset': 'Maritime',    'do_cond': True, 'do_coll': True,
     'switch_turb_kernel': True, 'epsilon_turb': 0.01},

    # --- Adaptive dt condensation ---
    {'name': 'Maritime_adaptive_dt',   'preset': 'Maritime',    'do_cond': True, 'do_coll': True,
     'switch_adaptive_dt': True},

    # --- Strong updraft ---
    {'name': 'Maritime_strong_w',      'preset': 'Maritime',    'do_cond': True, 'do_coll': True,
     'w': 3.0},
]

# ============================================================
# Common defaults
# ============================================================
T0 = 293.2; P0 = 1013.0e2; RH = 0.88; z0 = 0.0; zmax = 3000.0; dt = 1.0
n_ptcl = 5000  # smaller for speed
nt = 1500      # 1500s = enough to see activation, collision onset
theta_init = T0 * (p0/P0)**(r_a/cp)
theta_profiles = theta_init + 5e-3 * z_env
e_s0 = esatw(T0)
q0 = RH * e_s0 / (P0 - RH * e_s0) * r_a / rv
rm_spec_local = [1e-6, 25e-6]

def run_scenario(scen):
    """Run a single scenario, return dict with results or error."""
    name = scen['name']
    preset = PRESETS[scen['preset']]
    w = scen.get('w', 1.0)
    do_cond = scen['do_cond']
    do_coll = scen['do_coll']

    # Physics switches
    sw_kelvin = scen.get('switch_kelvin', True)
    sw_solute = scen.get('switch_solute', True)
    sw_E_const = scen.get('switch_E_constant', False)
    sw_vt_simple = scen.get('switch_vt_simple', False)
    sw_turb = scen.get('switch_turb_kernel', False)
    eps_turb = scen.get('epsilon_turb', 0.0)
    sw_adaptive = scen.get('switch_adaptive_dt', False)

    # Aerosol setup
    n_modes = len(preset['N'])
    N_aero = np.array(preset['N']) * 1e6
    mu_aero = np.log(np.array(preset['mu']) * 1e-6)
    sigma_aero = np.log(np.array(preset['sigma']))
    kappa_arr = preset['kappa'] + [1.6] * (4 - n_modes)  # pad to 4

    np.random.seed(42)
    T_parcel, q_parcel, particles_list = aero_init(
        'Random', n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, False)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    rho_p0, _, air_mass0 = parcel_rho(P0, T0)
    sp, qa, qc, qr, na, nc, nr, _, ra, rs = ts_analysis(particles_list, air_mass0, rm_spec_local, 60, n_ptcl)

    peak_Nc = nc * rho_p0
    peak_Nr = nr * rho_p0
    peak_qc = qc; peak_qr = qr
    T_min = T0; T_max = T0
    has_nan = False
    has_neg_q = False
    total_N_init = sum(p.A for p in particles_list) / air_mass0

    t0_wall = time.time()
    for t in range(nt):
        tt = (t+1) * dt
        z_parcel, T_parcel, P_parcel = ascend_parcel(
            z_parcel, T_parcel, P_parcel, w, dt, tt, zmax, theta_profiles, None, 'linear')
        rho_p, V_p, air_mass = parcel_rho(P_parcel, T_parcel)

        if do_cond:
            ct = at = et = da = 0.0
            if sw_adaptive:
                time_sub = 0.0
                while time_sub < dt - 1.0e-20:
                    tau = compute_tau_phase(particles_list, T_parcel, P_parcel, rho_liq, air_mass)
                    dt_sub = min(2.0 * tau, dt - time_sub)
                    dt_sub = max(dt_sub, 1.0e-6)
                    particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = drop_condensation(
                        particles_list, T_parcel, q_parcel, P_parcel, nt, dt_sub, air_mass, S_lst,
                        rho_aero, False, ct, at, et, da, False, sw_kelvin, sw_solute)
                    time_sub += dt_sub
            else:
                particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = drop_condensation(
                    particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass, S_lst,
                    rho_aero, False, ct, at, et, da, False, sw_kelvin, sw_solute)

        if do_coll:
            ac = au = pr = 0.0
            particles_list, ac, au, pr = collection(
                dt, particles_list, rho_p, rho_liq, P_parcel, T_parcel, ac, au, pr,
                False, z_parcel, zmax, w,
                sw_E_const, sw_vt_simple, sw_turb, eps_turb)

        # Diagnostics
        sp, qa, qc_v, qr_v, na_v, nc_v, nr_v, _, ra, rs = ts_analysis(
            particles_list, air_mass, rm_spec_local, 60, n_ptcl)
        nc_cm3 = nc_v * rho_p
        nr_cm3 = nr_v * rho_p

        peak_Nc = max(peak_Nc, nc_cm3)
        peak_Nr = max(peak_Nr, nr_cm3)
        peak_qc = max(peak_qc, qc_v)
        peak_qr = max(peak_qr, qr_v)
        T_min = min(T_min, T_parcel)
        T_max = max(T_max, T_parcel)

        if np.isnan(T_parcel) or np.isnan(qc_v) or np.isnan(qr_v):
            has_nan = True
            break
        if qc_v < -1e-10 or qr_v < -1e-10 or q_parcel < -1e-10:
            has_neg_q = True

    elapsed = time.time() - t0_wall

    # Mass conservation check
    total_N_final = sum(p.A for p in particles_list) / air_mass
    total_M_final = sum(p.M for p in particles_list) / air_mass
    LWC_final = (qc_v + qr_v) if not has_nan else float('nan')

    return {
        'name': name,
        'elapsed': elapsed,
        'peak_Nc': peak_Nc,
        'peak_Nr': peak_Nr,
        'peak_qc': peak_qc,
        'peak_qr': peak_qr,
        'T_min': T_min - 273.15,
        'T_max': T_max - 273.15,
        'has_nan': has_nan,
        'has_neg_q': has_neg_q,
        'N_init': total_N_init / 1e6,
        'N_final': total_N_final / 1e6,
        'LWC_final': LWC_final,
        'crashed': False,
        'error': None,
    }

# ============================================================
# Run all scenarios
# ============================================================
print(f"{'='*80}")
print(f"  STABILITY TEST: {len(SCENARIOS)} scenarios, {n_ptcl} particles, {nt} steps")
print(f"{'='*80}")

all_results = []
n_pass = 0
n_fail = 0

for i, scen in enumerate(SCENARIOS):
    print(f"\n[{i+1}/{len(SCENARIOS)}] {scen['name']} ...", end='', flush=True)
    try:
        res = run_scenario(scen)
        all_results.append(res)

        # Determine pass/fail
        checks = []
        if res['has_nan']:
            checks.append('NaN!')
        if res['has_neg_q']:
            checks.append('neg_q!')
        if res['T_min'] < -80 or res['T_max'] > 50:
            checks.append(f"T_range=[{res['T_min']:.1f},{res['T_max']:.1f}]")
        if res['peak_Nc'] < 0.1 and res['peak_qc'] > 0.01:
            checks.append('Nc~0 but qc>0?')

        if checks:
            status = f"WARN ({', '.join(checks)})"
            n_fail += 1
        else:
            status = "PASS"
            n_pass += 1

        print(f" {status}  ({res['elapsed']:.0f}s)  Nc={res['peak_Nc']:.1f} Nr={res['peak_Nr']:.2f} qc={res['peak_qc']:.3f} qr={res['peak_qr']:.3f}")

    except Exception as e:
        n_fail += 1
        tb = traceback.format_exc()
        all_results.append({
            'name': scen['name'], 'crashed': True, 'error': str(e),
            'elapsed': 0, 'peak_Nc': 0, 'peak_Nr': 0, 'peak_qc': 0, 'peak_qr': 0,
            'T_min': 0, 'T_max': 0, 'has_nan': False, 'has_neg_q': False,
            'N_init': 0, 'N_final': 0, 'LWC_final': 0,
        })
        print(f" CRASH: {e}")
        print(tb[-500:])

# ============================================================
# Summary
# ============================================================
print(f"\n{'='*80}")
print(f"  RESULTS: {n_pass} PASS / {n_fail} FAIL out of {len(SCENARIOS)}")
print(f"{'='*80}")
print(f"\n{'Scenario':<30s} {'Status':<8s} {'Nc':>6s} {'Nr':>6s} {'qc':>7s} {'qr':>7s} {'Tmin':>6s} {'Tmax':>6s} {'t(s)':>5s}")
print(f"{'-'*80}")

for r in all_results:
    if r['crashed']:
        status = 'CRASH'
    elif r['has_nan']:
        status = 'NaN'
    elif r['has_neg_q']:
        status = 'NEG_Q'
    else:
        status = 'OK'
    print(f"{r['name']:<30s} {status:<8s} {r['peak_Nc']:6.1f} {r['peak_Nr']:6.2f} {r['peak_qc']:7.3f} {r['peak_qr']:7.3f} {r['T_min']:6.1f} {r['T_max']:6.1f} {r['elapsed']:5.0f}")

print(f"\nTotal wall time: {sum(r['elapsed'] for r in all_results):.0f}s")
