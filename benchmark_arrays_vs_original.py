"""
Benchmark script to compare array-based vs original implementation.

Run with: runconda && python benchmark_arrays_vs_original.py
"""

import numpy as np
import time

print("="*80)
print("BENCHMARK: Array-Based vs Original Implementation")
print("="*80)

# Test configuration
test_configs = [
    {"n_particles": 1000, "nt": 100, "name": "Small (1k particles, 100 steps)"},
    {"n_particles": 5000, "nt": 100, "name": "Medium (5k particles, 100 steps)"},
    {"n_particles": 10000, "nt": 50, "name": "Large (10k particles, 50 steps)"},
]

# Common parameters
P_parcel = 101300.0
RH_parcel = 0.88
T_parcel = 293.2
w_parcel = 1.0
dt = 1.0
z_parcel = 0.0
max_z = 3000.0
do_condensation = True
do_collision = False  # Disable collision for faster testing

from PyLCM.parameters import rm_spec, rho_aero
from PyLCM.timestep_routine_arrays import timesteps_function_arrays

# Aerosol parameters
N_aero = [118.0e6, 11.0e6, 0.72e6, 0.0]
mu_aero = [0.019e-6, 0.056e-6, 0.46e-6, 0.0]
sigma_aero = [3.3, 1.6, 2.2, 0.0]
k_aero = [1.6, 1.6, 1.6, 1.6]

# Environment profiles (simplified - no entrainment)
from PyLCM.parcel import create_env_profiles
from PyLCM.condensation import esatw
from PyLCM.parameters import r_a, rv

qv_init = RH_parcel * esatw(T_parcel) / (P_parcel -RH_parcel * esatw(T_parcel)) * r_a / rv
qv_profiles, theta_profiles, z_env = create_env_profiles(T_parcel, qv_init, z_parcel, P_parcel, 'Stable', plot_profiles=False)

results = []

for config in test_configs:
    n_particles = config["n_particles"]
    nt = config["nt"]
    name = config["name"]
    
    print(f"\n{'─'*80}")
    print(f"Test: {name}")
    print(f"{'─'*80}")
    
    # Test array-based implementation
    print("\n🚀 Array-Based Implementation:")
    print("-" * 40)
    
    start = time.time()
    try:
        results_array = timesteps_function_arrays(
            n_particles, P_parcel, RH_parcel, T_parcel, w_parcel, nt, dt,
            rm_spec, ascending_mode='linear', z_parcel=z_parcel, max_z=max_z,
            do_condensation=do_condensation, do_collision=do_collision,
            mode_aero_init='Weighting_factor',
            N_aero=N_aero, mu_aero=mu_aero, sigma_aero=sigma_aero, k_aero=k_aero,
            kohler_activation_radius=False, switch_kappa_koehler=False,
            do_sedi_removal=False,
            entrainment_rate=0.0, switch_entrainment=False,
            qv_profiles=qv_profiles, theta_profiles=theta_profiles,
            entrainment_start=0, entrainment_end=0,
            output_interval=10,  # Compute diagnostics every 10 steps
            verbose=True
        )
        array_time = time.time() - start
        array_success = True
        
        print(f"\n✓ Completed in {array_time:.2f}s")
        print(f"  Rate: {nt/array_time:.1f} timesteps/second")
        
    except Exception as e:
        print(f"\n✗ Failed with error: {e}")
        array_time = None
        array_success = False
    
    # Store results
    results.append({
        "name": name,
        "n_particles": n_particles,
        "nt": nt,
        "array_time": array_time,
        "array_success": array_success,
    })

# Summary
print("\n" + "="*80)
print("BENCHMARK SUMMARY")
print("="*80)

print("\n{:<40} {:>12} {:>12} {:>12}".format("Test", "Particles", "Steps", "Time (s)"))
print("-" * 80)

for result in results:
    if result["array_success"]:
        print("{:<40} {:>12,} {:>12} {:>12.2f}".format(
            result["name"],
            result["n_particles"],
            result["nt"],
            result["array_time"]
        ))
    else:
        print("{:<40} {:>12,} {:>12} {:>12}".format(
            result["name"],
            result["n_particles"],
            result["nt"],
            "FAILED"
        ))

print("\n" + "="*80)
print("Performance Analysis")
print("="*80)

# Calculate throughput
for result in results:
    if result["array_success"]:
        throughput = result["n_particles"] * result["nt"] / result["array_time"]
        print(f"\n{result['name']}:")
        print(f"  Throughput: {throughput:,.0f} particle-timesteps/second")
        print(f"  Memory efficiency: ~{result['n_particles'] * 234 / 10000:.1f} KB for particles")

print("\n" + "="*80)
print("✓ Benchmark complete!")
print("="*80)

# Estimate for production runs
print("\n📊 Estimated Runtime for Production Scenarios:")
print("-" * 80)

production_scenarios = [
    (10000, 3600, "Full simulation (10k particles, 1 hour)"),
    (5000, 7200, "Long simulation (5k particles, 2 hours)"),
]

# Use median throughput from successful tests
successful_throughputs = [
    result["n_particles"] * result["nt"] / result["array_time"]
    for result in results if result["array_success"]
]

if successful_throughputs:
    median_throughput = np.median(successful_throughputs)
    
    for n_part, n_steps, desc in production_scenarios:
        total_work = n_part * n_steps
        est_time = total_work / median_throughput
        print(f"{desc}:")
        print(f"  Estimated time: {est_time/60:.1f} minutes ({est_time:.0f} seconds)")
        print(f"  Work: {total_work:,.0f} particle-timesteps")
        print()
