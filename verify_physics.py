#!/usr/bin/env python3
"""
Physics Verification Script for PyLCM

This script tests the physics correctness by:
1. Running a standard simulation and saving results
2. Comparing results between code versions
3. Specifically testing the kappa bug fix

Usage:
    python verify_physics.py --save-baseline    # Save results as baseline (run on main branch)
    python verify_physics.py --compare          # Compare current results with baseline
    python verify_physics.py --test-kappa       # Test kappa conservation specifically
"""

import numpy as np
import json
import os
import argparse
from pathlib import Path

# Import model modules
from PyLCM.parameters import *
from PyLCM.micro_particle import *
from PyLCM.aero_init import *
from PyLCM.parcel import *
from PyLCM.condensation import esatw

# Try optimized modules first
try:
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays as timesteps_function
    print("Using array-based timestep routine")
except ImportError:
    from PyLCM.timestep_routine import timesteps_function
    print("Using standard timestep routine")

BASELINE_FILE = Path(__file__).parent / "verification_baseline.json"


def run_standard_simulation(seed=42):
    """Run a standard simulation with fixed seed for reproducibility."""
    np.random.seed(seed)

    # Standard parameters
    n_particles = 500
    T_parcel = 293.2
    P_parcel = 101300.0
    RH_parcel = 0.88
    w_parcel = 1.0
    nt = 500
    dt = 1.0
    z_parcel = 0.0
    max_z = 3000.0

    # Aerosol with different kappa values to test mixing
    N_aero = [100e6, 50e6, 0.0, 0.0]  # m^-3
    mu_aero = [0.02e-6, 0.08e-6, 0.0, 0.0]  # m
    sigma_aero = [2.0, 1.8, 0.0, 0.0]
    k_aero = [0.3, 1.2, 0.0, 0.0]  # Different kappa values!

    # Create environment profiles
    qv_init = RH_parcel * esatw(T_parcel) / (P_parcel - RH_parcel * esatw(T_parcel)) * r_a / rv
    qv_profiles, theta_profiles, z_env = create_env_profiles(
        T_parcel, qv_init, z_parcel, P_parcel, 'Stable', plot_profiles=False
    )

    print(f"Running simulation: {n_particles} particles, {nt} steps...")

    results = timesteps_function(
        n_particles=n_particles,
        P_parcel=P_parcel,
        RH_parcel=RH_parcel,
        T_parcel=T_parcel,
        w_parcel=w_parcel,
        nt=nt,
        dt=dt,
        rm_spec=rm_spec,
        ascending_mode='linear',
        z_parcel=z_parcel,
        max_z=max_z,
        do_condensation=True,
        do_collision=True,  # Enable collision to test kappa mixing
        mode_aero_init='Random',
        N_aero=N_aero,
        mu_aero=mu_aero,
        sigma_aero=sigma_aero,
        k_aero=k_aero,
        kohler_activation_radius=False,
        switch_kappa_koehler=True,  # Enable kappa
        do_sedi_removal=False,
        entrainment_rate=0.0,
        switch_entrainment=False,
        qv_profiles=qv_profiles,
        theta_profiles=theta_profiles,
        entrainment_start=0,
        entrainment_end=0,
        output_interval=10,
        verbose=False
    )

    # Unpack results
    (nt_out, dt_out, time_array, T_parcel_array, RH_parcel_array, q_parcel_array,
     z_parcel_array, qa_ts, qc_ts, qr_ts, na_ts, nc_ts, nr_ts, spectra_arr,
     con_ts, act_ts, evp_ts, dea_ts, acc_ts, aut_ts, precip_ts,
     particles_array, rc_liq_avg_array, rc_liq_std_array) = results

    # Extract key metrics for comparison
    metrics = {
        "final_T": float(T_parcel_array[-1]),
        "final_RH": float(RH_parcel_array[-1]),
        "final_z": float(z_parcel_array[-1]),
        "max_qc": float(np.max(qc_ts)),
        "max_qr": float(np.max(qr_ts)),
        "total_condensation": float(np.sum(con_ts) * dt),
        "total_evaporation": float(np.sum(evp_ts) * dt),
        "total_accretion": float(np.sum(acc_ts) * dt),
        "total_autoconversion": float(np.sum(aut_ts) * dt),
        "final_mean_radius_um": float(rc_liq_avg_array[-1] * 1e6),
        "final_nc": float(nc_ts[-1]),
        "final_nr": float(nr_ts[-1]),
    }

    return metrics


def test_kappa_conservation():
    """
    Test that kappa is correctly conserved during collision.

    This specifically tests the bug fix where sequential kappa updates
    were using already-modified values.
    """
    print("\n=== Testing Kappa Conservation ===")

    # Import collision module
    from PyLCM.collision import same_weights_update
    from PyLCM.micro_particle import particles

    # Create two test particles with different kappa values
    class MockParticle:
        def __init__(self, M, A, Ns, kappa):
            self.M = M      # Total water mass
            self.A = A      # Weighting factor (multiplicity)
            self.Ns = Ns    # Aerosol mass
            self.kappa = kappa

    # Test case: Two particles with equal weights but different kappa
    # Volume-weighted average should be: (v1*k1 + v2*k2) / (v1 + v2)

    rho_liq = 1000.0

    # Particle 1: smaller droplet, kappa = 0.3 (organic)
    r1 = 10e-6  # 10 um radius
    v1 = 4/3 * np.pi * r1**3
    m1 = v1 * rho_liq
    A1 = 1000  # multiplicity

    # Particle 2: larger droplet, kappa = 1.2 (sea salt)
    r2 = 20e-6  # 20 um radius
    v2 = 4/3 * np.pi * r2**3
    m2 = v2 * rho_liq
    A2 = 1000  # same multiplicity (triggers same_weights_update)

    kappa1 = 0.3
    kappa2 = 1.2

    # Expected kappa after collision (volume-weighted average)
    expected_kappa = (v1 * kappa1 + v2 * kappa2) / (v1 + v2)

    # Create mock particles
    p1 = MockParticle(M=m1*A1, A=A1, Ns=1e-15*A1, kappa=kappa1)
    p2 = MockParticle(M=m2*A2, A=A2, Ns=1e-15*A2, kappa=kappa2)

    print(f"Before collision:")
    print(f"  Particle 1: r={r1*1e6:.1f} um, kappa={p1.kappa:.3f}")
    print(f"  Particle 2: r={r2*1e6:.1f} um, kappa={p2.kappa:.3f}")
    print(f"  Expected kappa after collision: {expected_kappa:.4f}")

    # Run the collision update
    p1_out, p2_out, _, _ = same_weights_update(p1, p2, 0.0, 0.0)

    print(f"\nAfter collision:")
    print(f"  Particle 1 kappa: {p1_out.kappa:.4f}")
    print(f"  Particle 2 kappa: {p2_out.kappa:.4f}")

    # Check results
    tolerance = 1e-6
    error1 = abs(p1_out.kappa - expected_kappa)
    error2 = abs(p2_out.kappa - expected_kappa)

    if error1 < tolerance and error2 < tolerance:
        print(f"\n[PASS] Kappa conservation test PASSED")
        print(f"  Both particles have correct kappa value")
        return True
    else:
        print(f"\n[FAIL] Kappa conservation test FAILED")
        print(f"  Error particle 1: {error1:.6f}")
        print(f"  Error particle 2: {error2:.6f}")

        # Diagnose the bug
        if abs(p1_out.kappa - p2_out.kappa) > tolerance:
            print(f"  BUG DETECTED: Particles have different kappa after collision!")
            print(f"  This indicates the sequential update bug is present.")
        return False


def save_baseline(metrics):
    """Save metrics as baseline for future comparison."""
    with open(BASELINE_FILE, 'w') as f:
        json.dump(metrics, f, indent=2)
    print(f"\nBaseline saved to: {BASELINE_FILE}")


def load_baseline():
    """Load baseline metrics."""
    if not BASELINE_FILE.exists():
        print(f"No baseline file found at: {BASELINE_FILE}")
        print("Run with --save-baseline first on the main branch.")
        return None

    with open(BASELINE_FILE, 'r') as f:
        return json.load(f)


def compare_results(current, baseline, tolerance=0.01):
    """Compare current results with baseline."""
    print("\n=== Comparing Results with Baseline ===")
    print(f"{'Metric':<30} {'Baseline':>15} {'Current':>15} {'Diff %':>10} {'Status':>10}")
    print("-" * 85)

    all_pass = True
    for key in baseline:
        base_val = baseline[key]
        curr_val = current[key]

        if base_val != 0:
            diff_pct = abs(curr_val - base_val) / abs(base_val) * 100
        else:
            diff_pct = abs(curr_val - base_val) * 100

        # Allow some tolerance for stochastic variations
        status = "PASS" if diff_pct < tolerance * 100 else "DIFF"
        if status == "DIFF":
            all_pass = False

        print(f"{key:<30} {base_val:>15.6f} {curr_val:>15.6f} {diff_pct:>9.2f}% {status:>10}")

    print("-" * 85)
    if all_pass:
        print("All metrics within tolerance. Physics unchanged.")
    else:
        print("Some metrics differ significantly. Review changes carefully.")

    return all_pass


def main():
    parser = argparse.ArgumentParser(description='Verify PyLCM physics')
    parser.add_argument('--save-baseline', action='store_true',
                        help='Save current results as baseline')
    parser.add_argument('--compare', action='store_true',
                        help='Compare current results with baseline')
    parser.add_argument('--test-kappa', action='store_true',
                        help='Test kappa conservation in collision')
    parser.add_argument('--all', action='store_true',
                        help='Run all tests')

    args = parser.parse_args()

    if not any([args.save_baseline, args.compare, args.test_kappa, args.all]):
        # Default: run kappa test
        args.test_kappa = True

    if args.test_kappa or args.all:
        kappa_ok = test_kappa_conservation()
        if not kappa_ok:
            print("\nWARNING: Kappa bug may still be present!")

    if args.save_baseline or args.all:
        print("\n=== Running Standard Simulation ===")
        metrics = run_standard_simulation()
        save_baseline(metrics)
        print("\nMetrics:")
        for k, v in metrics.items():
            print(f"  {k}: {v:.6f}")

    if args.compare:
        print("\n=== Running Standard Simulation ===")
        current = run_standard_simulation()
        baseline = load_baseline()
        if baseline:
            compare_results(current, baseline)


if __name__ == "__main__":
    main()
