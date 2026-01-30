#!/usr/bin/env python3
"""
Example: Run ensemble of cloud microphysics simulations.

This script demonstrates how to:
1. Run parallel ensemble simulations
2. Aggregate statistics across ensemble members
3. Check for convergence
4. Visualize results with uncertainty bands

The superdroplet method involves stochastic processes, so ensemble
averaging is needed for statistically robust results.
"""

import numpy as np
import matplotlib.pyplot as plt

from PyLCM.ensemble import run_ensemble, check_convergence, estimate_uncertainty_reduction
from PyLCM.parameters import rm_spec, r_a, cp, p0
from PyLCM.parcel import create_env_profiles


def main():
    # ==========================================================================
    # Configuration
    # ==========================================================================
    N_MEMBERS = 20       # Number of ensemble members (increase for convergence)
    N_WORKERS = 4        # Number of parallel workers (adjust to your CPU count)
    SEED_BASE = 42       # Base seed for reproducibility

    # Initial conditions
    P_parcel = 100000.0  # Initial pressure [Pa]
    T_parcel = 285.0     # Initial temperature [K]
    RH_parcel = 0.95     # Initial relative humidity

    # Create environmental profiles (required for parcel ascent)
    qv_profiles, theta_profiles, z_env = create_env_profiles(
        T_init=T_parcel,
        qv_init=0.01,
        z_init=0.0,
        p_env=P_parcel,
        stability_condition='Unstable',
        plot_profiles=False
    )

    # Simulation parameters
    sim_params = {
        'n_particles': 1000,      # Number of superdroplets
        'P_parcel': P_parcel,
        'RH_parcel': RH_parcel,
        'T_parcel': T_parcel,
        'w_parcel': 1.0,          # Updraft velocity [m/s]
        'nt': 300,                # Number of timesteps
        'dt': 1.0,                # Timestep [s]
        'rm_spec': rm_spec,       # Radius bins for spectrum
        'do_condensation': True,
        'do_collision': True,
        'mode_aero_init': 'Random',  # Random sampling (stochastic)
        'theta_profiles': theta_profiles,
        'qv_profiles': qv_profiles,
    }

    # ==========================================================================
    # Run ensemble
    # ==========================================================================
    print(f"Running {N_MEMBERS}-member ensemble...")
    print(f"  Particles: {sim_params['n_particles']}")
    print(f"  Timesteps: {sim_params['nt']} x {sim_params['dt']}s")
    print(f"  Workers: {N_WORKERS}")
    print()

    result = run_ensemble(
        n_members=N_MEMBERS,
        n_workers=N_WORKERS,
        seed_base=SEED_BASE,
        **sim_params
    )

    # ==========================================================================
    # Print summary statistics
    # ==========================================================================
    result.summary()

    # ==========================================================================
    # Check convergence
    # ==========================================================================
    print("\n" + "="*50)
    print("Convergence Analysis")
    print("="*50)

    for var in ['qc', 'qr', 'nc']:
        check_convergence(result, variable=var, threshold=0.05)
        print()

    # ==========================================================================
    # Analyze uncertainty scaling
    # ==========================================================================
    if N_MEMBERS >= 10:
        scaling = estimate_uncertainty_reduction(result, variable='qc')
        print(f"Uncertainty scaling analysis:")
        print(f"  Expected slope: {scaling['expected_slope']:.2f}")
        print(f"  Actual slope:   {scaling['actual_slope']:.2f}")
        print(f"  Scaling OK:     {scaling['scaling_ok']}")

    # ==========================================================================
    # Visualization
    # ==========================================================================
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Plot 1: Cloud water with uncertainty bands
    ax = axes[0, 0]
    ax.fill_between(result.time_array, result.qc_p10 * 1000, result.qc_p90 * 1000,
                    alpha=0.2, color='blue', label='10-90th percentile')
    ax.fill_between(result.time_array, result.qc_p25 * 1000, result.qc_p75 * 1000,
                    alpha=0.3, color='blue', label='25-75th percentile')
    ax.plot(result.time_array, result.qc_mean * 1000, 'b-', lw=2, label='Ensemble mean')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Cloud water (g/kg)')
    ax.set_title(f'Cloud Water Content (n={N_MEMBERS} members)')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    # Plot 2: Rain water with uncertainty bands
    ax = axes[0, 1]
    ax.fill_between(result.time_array, result.qr_p10 * 1000, result.qr_p90 * 1000,
                    alpha=0.2, color='red', label='10-90th percentile')
    ax.fill_between(result.time_array, result.qr_p25 * 1000, result.qr_p75 * 1000,
                    alpha=0.3, color='red', label='25-75th percentile')
    ax.plot(result.time_array, result.qr_mean * 1000, 'r-', lw=2, label='Ensemble mean')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Rain water (g/kg)')
    ax.set_title('Rain Water Content')
    ax.legend(loc='upper left')
    ax.grid(True, alpha=0.3)

    # Plot 3: Individual member trajectories (qc)
    ax = axes[1, 0]
    for i in range(min(10, N_MEMBERS)):
        ax.plot(result.time_array, result.qc_all[i] * 1000, alpha=0.5, lw=0.8)
    ax.plot(result.time_array, result.qc_mean * 1000, 'k-', lw=2, label='Mean')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Cloud water (g/kg)')
    ax.set_title('Individual Member Trajectories')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Plot 4: Running mean convergence
    ax = axes[1, 1]
    cumsum = np.cumsum(result.qc_all, axis=0)
    counts = np.arange(1, N_MEMBERS + 1)[:, None]
    running_means = cumsum / counts

    # Plot running mean at final time
    final_running_mean = running_means[:, -1] * 1000
    ax.plot(range(1, N_MEMBERS + 1), final_running_mean, 'b-', lw=2)
    ax.axhline(result.qc_mean[-1] * 1000, color='k', linestyle='--',
               label=f'Final mean: {result.qc_mean[-1]*1000:.4f} g/kg')
    ax.set_xlabel('Number of members')
    ax.set_ylabel('Running mean qc (g/kg)')
    ax.set_title('Running Mean Convergence')
    ax.legend()
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('ensemble_results.png', dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to ensemble_results.png")

    # ==========================================================================
    # Save results to CSV
    # ==========================================================================
    result.to_csv('ensemble_output')

    plt.show()


if __name__ == '__main__':
    main()
