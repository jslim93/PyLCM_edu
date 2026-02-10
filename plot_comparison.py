#!/usr/bin/env python3
"""Plot Fortran 10k vs PyLCM 10k: time series + DSD comparison.
   Fixes unit convention: PyLCM Nc/Nr multiplied by rho to get /cm3.
   Runs both with and without E_S09 coalescence efficiency."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import json, re, os, sys, time, warnings, copy
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
kappa_arr = [1.6] * 4
N_aero = N_raw * 1e6
n_ptcl = 10000
theta_init = T0 * (p0/P0)**(r_a/cp)
theta_profiles = theta_init + 5e-3 * z_env
e_s = esatw(T0)
q0 = RH * e_s / (P0 - RH * e_s) * r_a / rv
rm_spec = [1e-6, 25e-6]

def get_radii(plist):
    """Compute droplet radii from mass: r = (3M/(4*pi*rho_liq))^(1/3)"""
    return np.array([(3.0 * p.M / (4.0 * np.pi * rho_liq))**(1.0/3.0) for p in plist])

def run_pylcm(seed, use_E_S09=True, label=""):
    """Run PyLCM 10k and return time series with proper /cm3 units."""
    np.random.seed(seed)
    T_parcel, q_parcel, particles_list = aero_init(
        'Random', n_ptcl, P0, z0, T0, q0,
        N_aero, mu_aero, sigma_aero, rho_aero, kappa_arr, False)
    P_parcel = P0; z_parcel = z0; S_lst = 0.0

    dsd_times = [0, 300, 600, 900, 1200, 1500, 1800]
    dsd_snapshots = {}
    ts = {'t': [], 'qc': [], 'qr': [], 'Nc': [], 'Nr': [], 'T': [], 'rho': []}

    # t=0: use actual air_mass and rho
    rho_p0, V_p0, air_mass0 = parcel_rho(P0, T0)
    sp, qa, qc, qr, na, nc, nr, _, ra, rs = ts_analysis(particles_list, air_mass0, rm_spec, 60, n_ptcl)
    # Fix unit: ts_analysis returns /kg/1e6, multiply by rho to get /cm3
    nc_cm3 = nc * rho_p0
    nr_cm3 = nr * rho_p0
    ts['t'].append(0); ts['qc'].append(qc); ts['qr'].append(qr)
    ts['Nc'].append(nc_cm3); ts['Nr'].append(nr_cm3)
    ts['T'].append(T0 - 273.15); ts['rho'].append(rho_p0)
    radii_0 = get_radii(particles_list)
    dsd_snapshots[0] = radii_0

    t0_wall = time.time()
    e_s09_str = "WITH" if use_E_S09 else "WITHOUT"
    print(f"  [{label}] Running {e_s09_str} E_S09...")

    for t in range(nt):
        tt = (t+1) * dt
        z_parcel, T_parcel, P_parcel = ascend_parcel(
            z_parcel, T_parcel, P_parcel, w, dt, tt, zmax, theta_profiles, None, 'linear')
        rho_p, V_p, air_mass = parcel_rho(P_parcel, T_parcel)
        ct = at = et = da = 0.0
        particles_list, T_parcel, q_parcel, S_lst, ct, at, et, da = drop_condensation(
            particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass, S_lst,
            rho_aero, False, ct, at, et, da, False)
        ac = au = pr = 0.0
        particles_list, ac, au, pr = collection(
            dt, particles_list, rho_p, rho_liq, P_parcel, T_parcel, ac, au, pr,
            False, z_parcel, zmax, w, switch_E_S09=use_E_S09)

        sp, qa, qc_v, qr_v, na, nc_v, nr_v, _, ra, rs = ts_analysis(
            particles_list, air_mass, rm_spec, 60, n_ptcl)
        # Fix unit: multiply by rho
        nc_cm3 = nc_v * rho_p
        nr_cm3 = nr_v * rho_p
        ts['t'].append(t+1); ts['qc'].append(qc_v); ts['qr'].append(qr_v)
        ts['Nc'].append(nc_cm3); ts['Nr'].append(nr_cm3)
        ts['T'].append(T_parcel - 273.15); ts['rho'].append(rho_p)

        if (t+1) in dsd_times:
            radii = get_radii(particles_list)
            dsd_snapshots[t+1] = radii
            print(f"    t={t+1:4d}s  qc={qc_v:.4f}  qr={qr_v:.4f}  Nc={nc_cm3:.1f}  Nr={nr_cm3:.2f} /cm3")

        if (t+1) % 500 == 0:
            elapsed = time.time() - t0_wall
            print(f"    Step {t+1}/3000 ({elapsed:.0f}s elapsed)")

    elapsed_total = time.time() - t0_wall
    print(f"  [{label}] Done in {elapsed_total:.0f}s")

    for k in ts:
        ts[k] = np.array(ts[k])
    return ts, dsd_snapshots

# ============================================================
# PHASE 1: Run PyLCM with and without E_S09
# ============================================================
print("=== Phase 1a: PyLCM with E_S09 (default) ===")
ts_with, dsd_with = run_pylcm(seed=100, use_E_S09=True, label="E_S09=ON")

print("\n=== Phase 1b: PyLCM without E_S09 ===")
ts_without, dsd_without = run_pylcm(seed=100, use_E_S09=False, label="E_S09=OFF")

# ============================================================
# PHASE 2: Parse Fortran 10k time series
# ============================================================
print("\n=== Phase 2: Parsing Fortran 10k output ===")

fortran_ts_file = '/Users/dr.cloud/particle_model/output/time_series_compare_10k.dat'
with open(fortran_ts_file) as f:
    lines = f.readlines()

f_data = []
for line in lines[3:]:
    stripped = line.strip()
    if not stripped:
        continue
    fixed = re.sub(r'(\d)([-+])(\d{2,3})\b', r'\1E\2\3', stripped)
    row = []
    for v in fixed.split():
        try:
            row.append(float(v))
        except:
            row.append(float('nan'))
    if len(row) > 48:
        f_data.append(row)

f_data = np.array(f_data)
print(f"Fortran time series: {len(f_data)} steps, max t = {f_data[-1,0]:.0f}s")

# Extract Fortran time series (already in /cm3)
f_time = f_data[:, 0]
f_T = f_data[:, 4]
f_rho = f_data[:, 48]
f_QC = f_data[:, 10] * 1e3          # kg/kg -> g/kg
f_QR = f_data[:, 12] * 1e3
f_NC = f_data[:, 9] * f_rho / 1e6   # /kg -> /cm3
f_NR = f_data[:, 11] * f_rho / 1e6  # /kg -> /cm3

# ============================================================
# PHASE 3: Parse Fortran spectra
# ============================================================
print("\n=== Phase 3: Parsing Fortran spectra ===")

spectra_file = '/Users/dr.cloud/particle_model/output/spectra_compare_10k.dat'
with open(spectra_file) as f:
    slines = f.readlines()

def parse_spectra_line(line):
    fixed = re.sub(r'(\d)([-+])(\d{2,3})\b', r'\1E\2\3', line.strip())
    return np.array([float(v) for v in fixed.split()])

n_bins = int(slines[3].strip())
fortran_spectra = {}
rm_vals = parse_spectra_line(slines[4])
rm_bins = rm_vals[1:]

data_start = 7
n_data_lines = len(slines) - data_start
n_blocks = n_data_lines // 8

for b in range(n_blocks):
    base = data_start + b * 8
    if base + 1 >= len(slines):
        break
    n_line = parse_spectra_line(slines[base])
    q_line = parse_spectra_line(slines[base + 1])
    t_val = n_line[0]
    fortran_spectra[t_val] = {'rm': rm_bins, 'N': n_line[1:], 'Q': q_line[1:]}

available_spec_times = sorted(fortran_spectra.keys())
print(f"Fortran spectra at {len(available_spec_times)} timesteps")

# ============================================================
# PHASE 4: Plot time series comparison (3 lines: Fortran, PyLCM+E_S09, PyLCM-E_S09)
# ============================================================
print("\n=== Phase 4: Plotting time series ===")

max_common_t = min(f_time[-1], ts_with['t'][-1])

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('Fortran 10k vs PyLCM 10k (units: /cm$^3$, g/kg)', fontsize=14, fontweight='bold')

# qc
ax = axes[0, 0]
ax.plot(ts_with['t'], ts_with['qc'], 'b-', lw=1.5, label='PyLCM (E_S09 ON)', alpha=0.8)
ax.plot(ts_without['t'], ts_without['qc'], 'g--', lw=1.5, label='PyLCM (E_S09 OFF)', alpha=0.8)
ax.plot(f_time, f_QC, 'r-', lw=1.5, label='Fortran', alpha=0.8)
ax.set_xlabel('Time (s)'); ax.set_ylabel('qc (g/kg)')
ax.set_title('Cloud Water Mixing Ratio (qc)')
ax.legend(); ax.grid(True, alpha=0.3); ax.set_xlim(0, max_common_t)

# qr
ax = axes[0, 1]
ax.plot(ts_with['t'], ts_with['qr'], 'b-', lw=1.5, label='PyLCM (E_S09 ON)', alpha=0.8)
ax.plot(ts_without['t'], ts_without['qr'], 'g--', lw=1.5, label='PyLCM (E_S09 OFF)', alpha=0.8)
ax.plot(f_time, f_QR, 'r-', lw=1.5, label='Fortran', alpha=0.8)
ax.set_xlabel('Time (s)'); ax.set_ylabel('qr (g/kg)')
ax.set_title('Rain Water Mixing Ratio (qr)')
ax.legend(); ax.grid(True, alpha=0.3); ax.set_xlim(0, max_common_t)

# Nc (now both in /cm3)
ax = axes[1, 0]
ax.plot(ts_with['t'], ts_with['Nc'], 'b-', lw=1.5, label='PyLCM (E_S09 ON)', alpha=0.8)
ax.plot(ts_without['t'], ts_without['Nc'], 'g--', lw=1.5, label='PyLCM (E_S09 OFF)', alpha=0.8)
ax.plot(f_time, f_NC, 'r-', lw=1.5, label='Fortran', alpha=0.8)
ax.set_xlabel('Time (s)'); ax.set_ylabel('Nc (/cm$^3$)')
ax.set_title('Cloud Droplet Number (Nc) - corrected units')
ax.legend(); ax.grid(True, alpha=0.3); ax.set_xlim(0, max_common_t)

# Nr (now both in /cm3)
ax = axes[1, 1]
ax.plot(ts_with['t'], ts_with['Nr'], 'b-', lw=1.5, label='PyLCM (E_S09 ON)', alpha=0.8)
ax.plot(ts_without['t'], ts_without['Nr'], 'g--', lw=1.5, label='PyLCM (E_S09 OFF)', alpha=0.8)
ax.plot(f_time, f_NR, 'r-', lw=1.5, label='Fortran', alpha=0.8)
ax.set_xlabel('Time (s)'); ax.set_ylabel('Nr (/cm$^3$)')
ax.set_title('Rain Drop Number (Nr) - corrected units')
ax.legend(); ax.grid(True, alpha=0.3); ax.set_xlim(0, max_common_t)

plt.tight_layout()
ts_path = '/Users/dr.cloud/PyLCM_edu/comparison_timeseries.png'
plt.savefig(ts_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"Time series plot saved: {ts_path}")

# ============================================================
# PHASE 5: Plot DSD comparison (3 lines)
# ============================================================
print("\n=== Phase 5: DSD comparison ===")

dsd_times = [0, 300, 600, 900, 1200, 1500, 1800]
fig, axes = plt.subplots(2, 4, figsize=(20, 10))
fig.suptitle('DSD: Fortran vs PyLCM (E_S09 ON/OFF)', fontsize=14, fontweight='bold')

for i, t_target in enumerate(dsd_times):
    row = i // 4; col = i % 4
    ax = axes[row, col]

    # PyLCM DSD with E_S09
    for dsd_snap, ts_d, color, lbl in [(dsd_with, ts_with, 'b', 'PyLCM E_S09 ON'),
                                        (dsd_without, ts_without, 'g', 'PyLCM E_S09 OFF')]:
        if t_target in dsd_snap:
            radii_um = dsd_snap[t_target] * 1e6
            radii_um = radii_um[radii_um > 0]
            if len(radii_um) > 0:
                r_min = max(0.01, radii_um.min() * 0.5)
                r_max = min(10000, radii_um.max() * 2)
                bins = np.logspace(np.log10(r_min), np.log10(r_max), 60)
                bin_centers = np.sqrt(bins[:-1] * bins[1:])
                dlnr = np.diff(np.log(bins))
                counts, _ = np.histogram(radii_um, bins=bins)
                t_idx = t_target
                if t_idx < len(ts_d['Nc']):
                    total_N = ts_d['Nc'][t_idx] + ts_d['Nr'][t_idx]
                    dNdlnr = counts / n_ptcl * total_N / dlnr
                else:
                    dNdlnr = counts / dlnr
                ls = '-' if 'ON' in lbl else '--'
                ax.plot(bin_centers, dNdlnr, color=color, ls=ls, lw=1.5, label=lbl)

    # Fortran DSD
    if available_spec_times:
        closest_t = min(available_spec_times, key=lambda x: abs(x - t_target))
        if abs(closest_t - t_target) < 5:
            spec = fortran_spectra[closest_t]
            rm_um = spec['rm'] * 1e6
            N_per_kg_per_m = spec['N']
            f_idx = np.argmin(np.abs(f_time - t_target))
            rho_f = f_rho[f_idx] if f_idx < len(f_rho) else 1.2
            dNdlnr_f = N_per_kg_per_m * spec['rm'] * rho_f / 1e6
            mask = dNdlnr_f > 0
            if np.any(mask):
                ax.plot(rm_um[mask], dNdlnr_f[mask], 'r-', lw=1.5, label='Fortran')

    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('Radius (um)'); ax.set_ylabel('dN/dlnr (/cm3)')
    ax.set_title(f't = {t_target}s')
    ax.set_xlim(0.01, 5000); ax.set_ylim(0.001, 1e4)
    ax.legend(fontsize=7); ax.grid(True, alpha=0.3, which='both')

if len(dsd_times) < 8:
    axes[1, 3].set_visible(False)

plt.tight_layout()
dsd_path = '/Users/dr.cloud/PyLCM_edu/comparison_dsd.png'
plt.savefig(dsd_path, dpi=150, bbox_inches='tight')
plt.close()
print(f"DSD plot saved: {dsd_path}")

# ============================================================
# Summary statistics
# ============================================================
print("\n=== Summary: Peak Nr comparison ===")
print(f"  Fortran peak Nr:         {np.max(f_NR):.2f} /cm3 at t={f_time[np.argmax(f_NR)]:.0f}s")
print(f"  PyLCM (E_S09 ON)  peak:  {np.max(ts_with['Nr']):.2f} /cm3 at t={ts_with['t'][np.argmax(ts_with['Nr'])]:.0f}s")
print(f"  PyLCM (E_S09 OFF) peak:  {np.max(ts_without['Nr']):.2f} /cm3 at t={ts_without['t'][np.argmax(ts_without['Nr'])]:.0f}s")
print(f"\n  Fortran peak Nc:         {np.max(f_NC):.2f} /cm3")
print(f"  PyLCM (E_S09 ON)  peak:  {np.max(ts_with['Nc']):.2f} /cm3")
print(f"  PyLCM (E_S09 OFF) peak:  {np.max(ts_without['Nc']):.2f} /cm3")

print("\n=== Done! ===")
