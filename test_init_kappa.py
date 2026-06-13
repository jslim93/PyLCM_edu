#!/usr/bin/env python3
"""Sensitivity test: init method (Random vs Weighting_factor) x kappa (1.6 vs 1.28).
   Compares with Fortran 10k reference. Condensation+collision, 3000 steps."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import copy, re, sys, time, warnings
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
# Common setup
# ============================================================
T0=293.2; P0=1013.0e2; RH=0.88; w=1.0; z0=0.0; zmax=3000.0; dt=1.0; nt=3000
N_raw = np.array([118.0, 11.0, 0.72])
mu_aero = np.log(np.array([0.019, 0.056, 0.46]) * 1e-6)
sigma_aero = np.log(np.array([3.3, 1.6, 2.2]))
N_aero = N_raw * 1e6
n_ptcl = 10000
theta_init = T0 * (p0/P0)**(r_a/cp)
theta_profiles = theta_init + 5e-3 * z_env
e_s0 = esatw(T0)
q0 = RH * e_s0 / (P0 - RH * e_s0) * r_a / rv
rm_spec = [1e-6, 25e-6]

def get_radii(plist):
    return np.array([(3.0*p.M/p.A/(4.0*np.pi*rho_liq))**(1.0/3.0) for p in plist])

# ============================================================
# Define 4 test cases
# ============================================================
cases = [
    {'name': 'Random_k1.6',       'init': 'Random',           'kappa': 1.6,  'color': 'b',  'ls': '-'},
    {'name': 'Random_k1.28',      'init': 'Random',           'kappa': 1.28, 'color': 'c',  'ls': '-'},
    {'name': 'Weighting_k1.6',    'init': 'Weighting_factor', 'kappa': 1.6,  'color': 'g',  'ls': '--'},
    {'name': 'Weighting_k1.28',   'init': 'Weighting_factor', 'kappa': 1.28, 'color': 'm',  'ls': '--'},
]

results = {}

for case in cases:
    cname = case['name']
    print(f"\n{'='*60}")
    print(f"  Case: {cname}  (init={case['init']}, kappa={case['kappa']})")
    print(f"{'='*60}")

    kappa_arr = [case['kappa']] * 4
    np.random.seed(100)

    T_parcel, q_parcel, particles_list = aero_init(
        case['init'], n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, False)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    # t=0 DSD snapshot
    rho_p0, V_p0, air_mass0 = parcel_rho(P0, T0)
    dsd_t0 = get_radii(particles_list)

    # Check A values
    A_vals = [p.A for p in particles_list]
    print(f"  A: min={min(A_vals):.1f}  max={max(A_vals):.1f}  mean={np.mean(A_vals):.1f}  unique={len(set([round(a,2) for a in A_vals]))}")

    ts = {'t': [], 'qc': [], 'qr': [], 'Na': [], 'Nc': [], 'Nr': [], 'T': []}

    sp, qa, qc, qr, na, nc, nr, _, ra, rs = ts_analysis(particles_list, air_mass0, rm_spec, 60, n_ptcl)
    ts['t'].append(0); ts['qc'].append(qc); ts['qr'].append(qr)
    ts['Na'].append(na*rho_p0); ts['Nc'].append(nc*rho_p0); ts['Nr'].append(nr*rho_p0)
    ts['T'].append(T0-273.15)

    t0_wall = time.time()
    for t in range(nt):
        tt = (t+1)*dt
        z_parcel, T_parcel, P_parcel = ascend_parcel(
            z_parcel, T_parcel, P_parcel, w, dt, tt, zmax, theta_profiles, None, 'linear')
        rho_p, V_p, air_mass = parcel_rho(P_parcel, T_parcel)
        ct=at=et=da=0.0
        particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = drop_condensation(
            particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass, S_lst,
            rho_aero, False, ct, at, et, da, False)
        ac=au=pr=0.0
        particles_list, ac, au, pr = collection(
            dt, particles_list, rho_p, rho_liq, P_parcel, T_parcel, ac, au, pr,
            False, z_parcel, zmax, w)

        sp, qa, qc_v, qr_v, na, nc_v, nr_v, _, ra, rs = ts_analysis(
            particles_list, air_mass, rm_spec, 60, n_ptcl)
        rho_p_val = rho_p
        ts['t'].append(t+1); ts['qc'].append(qc_v); ts['qr'].append(qr_v)
        ts['Na'].append(na*rho_p_val); ts['Nc'].append(nc_v*rho_p_val); ts['Nr'].append(nr_v*rho_p_val)
        ts['T'].append(T_parcel-273.15)

        if (t+1) % 500 == 0:
            elapsed = time.time()-t0_wall
            print(f"  Step {t+1}/3000 ({elapsed:.0f}s)  Nc={nc_v*rho_p_val:.1f}  Nr={nr_v*rho_p_val:.2f}  qc={qc_v:.3f}  qr={qr_v:.3f}")

    elapsed = time.time()-t0_wall
    for k in ts: ts[k] = np.array(ts[k])

    results[cname] = {'ts': ts, 'dsd_t0': dsd_t0, 'elapsed': elapsed,
                      'case': case, 'peak_Nc': np.max(ts['Nc']), 'peak_Nr': np.max(ts['Nr'])}
    print(f"  Done in {elapsed:.0f}s  peak Nc={np.max(ts['Nc']):.1f}  peak Nr={np.max(ts['Nr']):.2f}")

# ============================================================
# Parse Fortran reference
# ============================================================
print("\n=== Parsing Fortran 10k reference ===")
fortran_ts_file = '/Users/dr.cloud/particle_model/output/time_series_compare_10k.dat'
with open(fortran_ts_file) as f:
    lines = f.readlines()
f_data = []
for line in lines[3:]:
    stripped = line.strip()
    if not stripped: continue
    fixed = re.sub(r'(\d)([-+])(\d{2,3})\b', r'\1E\2\3', stripped)
    row = []
    for v in fixed.split():
        try: row.append(float(v))
        except: row.append(float('nan'))
    if len(row) > 48: f_data.append(row)
f_data = np.array(f_data)
f_time = f_data[:,0]; f_rho = f_data[:,48]
f_QC = f_data[:,10]*1e3; f_QR = f_data[:,12]*1e3
f_NC = f_data[:,9]*f_rho/1e6; f_NR = f_data[:,11]*f_rho/1e6

# ============================================================
# Plot time series comparison (all 4 + Fortran)
# ============================================================
print("\n=== Plotting ===")
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Init Method x Kappa Sensitivity (vs Fortran 10k)', fontsize=14, fontweight='bold')

for ax, key, f_arr, ylabel, title in [
    (axes[0,0], 'qc', f_QC, 'qc (g/kg)', 'Cloud Water (qc)'),
    (axes[0,1], 'qr', f_QR, 'qr (g/kg)', 'Rain Water (qr)'),
    (axes[1,0], 'Nc', f_NC, 'Nc (/cm$^3$)', 'Cloud Droplet Number (Nc)'),
    (axes[1,1], 'Nr', f_NR, 'Nr (/cm$^3$)', 'Rain Drop Number (Nr)')]:
    ax.plot(f_time, f_arr, 'r-', lw=2.5, label='Fortran 10k', alpha=0.7)
    for cname, res in results.items():
        c = res['case']
        ax.plot(res['ts']['t'], res['ts'][key], color=c['color'], ls=c['ls'],
                lw=1.5, label=cname, alpha=0.8)
    ax.set_xlabel('Time (s)'); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_xlim(0, 3000)

plt.tight_layout()
plt.savefig('/Users/dr.cloud/PyLCM_edu/sensitivity_init_kappa_ts.png', dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Plot t=0 DSD comparison
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(20, 5))
fig.suptitle('t=0 Aerosol DSD: Init Method x Kappa', fontsize=14, fontweight='bold')

# Parse Fortran t=0 spectrum
spectra_file = '/Users/dr.cloud/particle_model/output/spectra_compare_10k.dat'
with open(spectra_file) as f:
    slines = f.readlines()
def parse_sl(line):
    fixed = re.sub(r'(\d)([-+])(\d{2,3})\b', r'\1E\2\3', line.strip())
    return np.array([float(v) for v in fixed.split()])
rm_bins = parse_sl(slines[4])[1:]
n_line0 = parse_sl(slines[7])
fortran_t0_N = n_line0[1:]
rho_f0 = f_rho[0]

for i, (cname, res) in enumerate(results.items()):
    ax = axes[i]
    radii_um = res['dsd_t0'] * 1e6
    radii_um = radii_um[radii_um > 0]
    bins = np.logspace(np.log10(max(0.001, radii_um.min()*0.5)),
                       np.log10(min(10000, radii_um.max()*2)), 60)
    bin_c = np.sqrt(bins[:-1]*bins[1:])
    dlnr = np.diff(np.log(bins))
    counts, _ = np.histogram(radii_um, bins=bins)
    total_N = res['ts']['Na'][0] + res['ts']['Nc'][0] + res['ts']['Nr'][0]
    dNdlnr = counts / n_ptcl * total_N / dlnr

    c = res['case']
    ax.plot(bin_c, dNdlnr, color=c['color'], lw=1.5, label=f'PyLCM {cname}')
    # Fortran
    dNdlnr_f = fortran_t0_N * rm_bins * rho_f0 / 1e6
    mask = dNdlnr_f > 0
    ax.plot(rm_bins[mask]*1e6, dNdlnr_f[mask], 'r-', lw=1.5, label='Fortran', alpha=0.7)

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Radius (um)'); ax.set_ylabel('dN/dlnr (/cm3)')
    ax.set_title(cname, fontsize=10)
    ax.set_xlim(0.001, 10); ax.set_ylim(1, 1e5)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

plt.tight_layout()
plt.savefig('/Users/dr.cloud/PyLCM_edu/sensitivity_init_kappa_dsd.png', dpi=150, bbox_inches='tight')
plt.close()

# ============================================================
# Summary table
# ============================================================
print(f"\n{'='*70}")
print(f"{'Case':<22s}  {'peak Nc':>8s}  {'peak Nr':>8s}  {'time(s)':>8s}")
print(f"{'='*70}")
print(f"{'Fortran 10k':<22s}  {np.max(f_NC):8.1f}  {np.max(f_NR):8.2f}  {'(ref)':>8s}")
for cname, res in results.items():
    print(f"{cname:<22s}  {res['peak_Nc']:8.1f}  {res['peak_Nr']:8.2f}  {res['elapsed']:8.0f}")
print(f"{'='*70}")
print("\nPlots saved:")
print("  sensitivity_init_kappa_ts.png")
print("  sensitivity_init_kappa_dsd.png")
