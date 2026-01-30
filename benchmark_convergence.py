#!/usr/bin/env python3
"""
Convergence Comparison: Weighting_factor vs Random+Ensemble

Compare convergence approaches for superdroplet simulations:
1. Weighting_factor (deterministic) - Single run with linear PDF sampling
2. Random + Ensemble - Multiple stochastic runs averaged

Metrics:
- Coefficient of Variation (CV) = std/mean × 100%
- Lower CV = more converged
- Target: CV < 5% for well-converged results

Variables tracked:
- qc (cloud water) - condensation dominated
- qr (rain water) - collision dominated
- nc (cloud number)
- nr (rain number)
"""

import numpy as np
import time
import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from PyLCM.timestep_routine_arrays import timesteps_function_arrays
from PyLCM.ensemble import run_ensemble_serial, run_ensemble
from PyLCM.parameters import rm_spec
from PyLCM.parcel import create_env_profiles


def compute_cv(values):
    """
    Compute Coefficient of Variation (CV) = std/mean × 100%.

    Args:
        values: Array of values from multiple runs

    Returns:
        float: CV as percentage, or NaN if mean is zero
    """
    mean = np.mean(values)
    std = np.std(values)
    if mean == 0:
        return np.nan
    return (std / abs(mean)) * 100.0


def run_weighting_factor_trials(n_particles, n_trials, sim_params, init_seed=42):
    """
    Run multiple Weighting_factor trials.

    Weighting_factor uses deterministic PDF sampling for aerosol initialization.
    Stochasticity only comes from collision process with different collision seeds.

    Args:
        n_particles: Number of superdroplets
        n_trials: Number of independent trials
        sim_params: Dictionary of simulation parameters
        init_seed: Fixed seed for initialization (default: 42)

    Returns:
        dict: Final values for qc, qr, nc, nr for each trial
    """
    results = {
        'qc': [],
        'qr': [],
        'nc': [],
        'nr': [],
        'times': []
    }

    for trial in range(n_trials):
        collision_seed = trial * 1000  # Different collision seed for each trial

        start_time = time.time()
        output = timesteps_function_arrays(
            n_particles=n_particles,
            mode_aero_init='Weighting_factor',
            seed=init_seed,  # Fixed initialization seed
            collision_seed=collision_seed,  # Varying collision seed
            verbose=False,
            **sim_params
        )
        elapsed = time.time() - start_time

        # Extract final values (index 8=qc, 9=qr, 11=nc, 12=nr)
        qc_final = output[8][-1]  # qc_ts
        qr_final = output[9][-1]  # qr_ts
        nc_final = output[11][-1]  # nc_ts
        nr_final = output[12][-1]  # nr_ts

        results['qc'].append(qc_final)
        results['qr'].append(qr_final)
        results['nc'].append(nc_final)
        results['nr'].append(nr_final)
        results['times'].append(elapsed)

        print(f"    Trial {trial+1}/{n_trials}: {elapsed:.1f}s, "
              f"qc={qc_final*1e3:.4f} g/kg, qr={qr_final*1e3:.4f} g/kg")

    return results


def run_random_ensemble_trials(n_particles, n_trials, n_members, sim_params, init_seed=42):
    """
    Run multiple Random+Ensemble trials.

    Each trial runs an ensemble of n_members realizations with:
    - Same initialization (deterministic, using init_seed)
    - Different collision seeds (stochastic collision)

    This isolates collision stochasticity from initialization randomness.

    Args:
        n_particles: Number of superdroplets
        n_trials: Number of independent ensemble trials
        n_members: Number of members per ensemble
        sim_params: Dictionary of simulation parameters
        init_seed: Fixed seed for initialization (default: 42)

    Returns:
        dict: Ensemble mean values for qc, qr, nc, nr for each trial
    """
    results = {
        'qc': [],
        'qr': [],
        'nc': [],
        'nr': [],
        'times': []
    }

    for trial in range(n_trials):
        collision_seed_base = trial * 1000  # Different collision seed base for each trial

        start_time = time.time()

        # Use parallel ensemble if multiple workers available
        n_workers = min(os.cpu_count() or 1, n_members)

        try:
            ensemble_result = run_ensemble(
                n_members=n_members,
                seed_base=collision_seed_base,  # Collision seeds: base+0, base+1, ...
                init_seed=init_seed,  # Fixed initialization for all members
                n_workers=n_workers,
                method='fork',
                n_particles=n_particles,
                mode_aero_init='Random',
                **sim_params
            )
        except Exception as e:
            # Fall back to serial if parallel fails
            print(f"    Parallel failed ({e}), using serial...")
            ensemble_result = run_ensemble_serial(
                n_members=n_members,
                seed_base=collision_seed_base,
                init_seed=init_seed,
                n_particles=n_particles,
                mode_aero_init='Random',
                **sim_params
            )

        elapsed = time.time() - start_time

        # Get ensemble mean of final values
        qc_mean = ensemble_result.qc_mean[-1]
        qr_mean = ensemble_result.qr_mean[-1]
        nc_mean = ensemble_result.nc_mean[-1]
        nr_mean = ensemble_result.nr_mean[-1]

        results['qc'].append(qc_mean)
        results['qr'].append(qr_mean)
        results['nc'].append(nc_mean)
        results['nr'].append(nr_mean)
        results['times'].append(elapsed)

        print(f"    Trial {trial+1}/{n_trials}: {elapsed:.1f}s, "
              f"qc={qc_mean*1e3:.4f} g/kg, qr={qr_mean*1e3:.4f} g/kg")

    return results


def analyze_results(wf_results, re_results, particle_counts, n_members):
    """
    Analyze and display convergence comparison results.

    Args:
        wf_results: Dict mapping particle count to Weighting_factor results
        re_results: Dict mapping particle count to Random+Ensemble results
        particle_counts: List of particle counts tested
        n_members: Number of ensemble members used
    """
    print("\n" + "="*80)
    print("CONVERGENCE COMPARISON RESULTS")
    print("="*80)

    # Table header
    print("\n" + "-"*80)
    print(f"{'Approach':<25} {'Particles':>10} {'CV(qc)%':>10} {'CV(qr)%':>10} "
          f"{'CV(nc)%':>10} {'CV(nr)%':>10} {'Time(s)':>8}")
    print("-"*80)

    summary = {
        'particle_counts': particle_counts,
        'wf_cv_qc': [], 'wf_cv_qr': [], 'wf_cv_nc': [], 'wf_cv_nr': [],
        're_cv_qc': [], 're_cv_qr': [], 're_cv_nc': [], 're_cv_nr': [],
        'wf_mean_qc': [], 'wf_mean_qr': [],
        're_mean_qc': [], 're_mean_qr': [],
        'wf_times': [], 're_times': []
    }

    for n_ptcl in particle_counts:
        # Weighting_factor results
        wf = wf_results[n_ptcl]
        cv_qc_wf = compute_cv(wf['qc'])
        cv_qr_wf = compute_cv(wf['qr'])
        cv_nc_wf = compute_cv(wf['nc'])
        cv_nr_wf = compute_cv(wf['nr'])
        mean_time_wf = np.mean(wf['times'])

        summary['wf_cv_qc'].append(cv_qc_wf)
        summary['wf_cv_qr'].append(cv_qr_wf)
        summary['wf_cv_nc'].append(cv_nc_wf)
        summary['wf_cv_nr'].append(cv_nr_wf)
        summary['wf_mean_qc'].append(np.mean(wf['qc']))
        summary['wf_mean_qr'].append(np.mean(wf['qr']))
        summary['wf_times'].append(mean_time_wf)

        print(f"{'Weighting_factor':<25} {n_ptcl:>10} {cv_qc_wf:>10.2f} {cv_qr_wf:>10.2f} "
              f"{cv_nc_wf:>10.2f} {cv_nr_wf:>10.2f} {mean_time_wf:>8.1f}")

        # Random+Ensemble results
        re = re_results[n_ptcl]
        cv_qc_re = compute_cv(re['qc'])
        cv_qr_re = compute_cv(re['qr'])
        cv_nc_re = compute_cv(re['nc'])
        cv_nr_re = compute_cv(re['nr'])
        mean_time_re = np.mean(re['times'])

        summary['re_cv_qc'].append(cv_qc_re)
        summary['re_cv_qr'].append(cv_qr_re)
        summary['re_cv_nc'].append(cv_nc_re)
        summary['re_cv_nr'].append(cv_nr_re)
        summary['re_mean_qc'].append(np.mean(re['qc']))
        summary['re_mean_qr'].append(np.mean(re['qr']))
        summary['re_times'].append(mean_time_re)

        print(f"{'Random+Ensemble(' + str(n_members) + ')':<25} {n_ptcl:>10} {cv_qc_re:>10.2f} {cv_qr_re:>10.2f} "
              f"{cv_nc_re:>10.2f} {cv_nr_re:>10.2f} {mean_time_re:>8.1f}")
        print()

    print("-"*80)

    # Analysis
    print("\n" + "="*80)
    print("ANALYSIS")
    print("="*80)

    print("\n** Initialization is DETERMINISTIC (fixed seed=42) for both approaches **")
    print("** Stochasticity comes ONLY from collision process **")

    print("\n1. Cloud Water (qc) - Condensation Dominated:")
    print("   - Condensation is deterministic given same initial particles")
    print("   - Variability in qc comes from collision's effect on droplet sizes")

    print("\n2. Rain Water (qr) - Collision Dominated:")
    print("   - Collision process is stochastic (random pairing & probability)")
    print("   - Ensemble averaging reduces collision noise")

    print("\n3. Mean Value Comparison (checking for bias):")
    for i, n_ptcl in enumerate(particle_counts):
        qc_wf = summary['wf_mean_qc'][i] * 1e3
        qc_re = summary['re_mean_qc'][i] * 1e3
        qc_diff = (qc_wf - qc_re) / qc_re * 100 if qc_re != 0 else 0

        qr_wf = summary['wf_mean_qr'][i] * 1e3
        qr_re = summary['re_mean_qr'][i] * 1e3
        qr_diff = (qr_wf - qr_re) / qr_re * 100 if qr_re != 0 else 0

        print(f"   {n_ptcl} particles: qc_wf={qc_wf:.4f}, qc_re={qc_re:.4f} (diff: {qc_diff:+.1f}%)")
        print(f"                  qr_wf={qr_wf:.4f}, qr_re={qr_re:.4f} (diff: {qr_diff:+.1f}%)")

    print("\n4. Computational Efficiency:")
    print("   Approach that achieves CV < 5% with least computational cost:")
    for var in ['qc', 'qr']:
        print(f"\n   {var}:")
        for i, n_ptcl in enumerate(particle_counts):
            cv_wf = summary[f'wf_cv_{var}'][i]
            cv_re = summary[f're_cv_{var}'][i]
            time_wf = summary['wf_times'][i]
            time_re = summary['re_times'][i]

            wf_converged = cv_wf < 5 if not np.isnan(cv_wf) else False
            re_converged = cv_re < 5 if not np.isnan(cv_re) else False

            status_wf = "CONVERGED" if wf_converged else "not converged"
            status_re = "CONVERGED" if re_converged else "not converged"

            print(f"   {n_ptcl:>6} ptcl: WF {status_wf} (CV={cv_wf:.1f}%, {time_wf:.1f}s) | "
                  f"RE {status_re} (CV={cv_re:.1f}%, {time_re:.1f}s)")

    return summary


def plot_results(summary, output_file='convergence_comparison.png'):
    """
    Create visualization of convergence comparison.

    Args:
        summary: Dictionary with analysis results
        output_file: Output filename for plot
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("\nMatplotlib not available, skipping plots.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    particle_counts = summary['particle_counts']

    # Plot 1: CV for qc
    ax1 = axes[0, 0]
    ax1.plot(particle_counts, summary['wf_cv_qc'], 'b-o', label='Weighting_factor', linewidth=2, markersize=8)
    ax1.plot(particle_counts, summary['re_cv_qc'], 'r-s', label='Random+Ensemble', linewidth=2, markersize=8)
    ax1.axhline(y=5, color='g', linestyle='--', label='CV=5% threshold')
    ax1.set_xlabel('Number of Particles')
    ax1.set_ylabel('CV (%)')
    ax1.set_title('Cloud Water (qc) - Condensation Dominated')
    ax1.legend()
    ax1.set_xscale('log')
    ax1.grid(True, alpha=0.3)

    # Plot 2: CV for qr
    ax2 = axes[0, 1]
    ax2.plot(particle_counts, summary['wf_cv_qr'], 'b-o', label='Weighting_factor', linewidth=2, markersize=8)
    ax2.plot(particle_counts, summary['re_cv_qr'], 'r-s', label='Random+Ensemble', linewidth=2, markersize=8)
    ax2.axhline(y=5, color='g', linestyle='--', label='CV=5% threshold')
    ax2.set_xlabel('Number of Particles')
    ax2.set_ylabel('CV (%)')
    ax2.set_title('Rain Water (qr) - Collision Dominated')
    ax2.legend()
    ax2.set_xscale('log')
    ax2.grid(True, alpha=0.3)

    # Plot 3: CV for nc
    ax3 = axes[1, 0]
    ax3.plot(particle_counts, summary['wf_cv_nc'], 'b-o', label='Weighting_factor', linewidth=2, markersize=8)
    ax3.plot(particle_counts, summary['re_cv_nc'], 'r-s', label='Random+Ensemble', linewidth=2, markersize=8)
    ax3.axhline(y=5, color='g', linestyle='--', label='CV=5% threshold')
    ax3.set_xlabel('Number of Particles')
    ax3.set_ylabel('CV (%)')
    ax3.set_title('Cloud Number (nc)')
    ax3.legend()
    ax3.set_xscale('log')
    ax3.grid(True, alpha=0.3)

    # Plot 4: Computational time comparison
    ax4 = axes[1, 1]
    ax4.plot(particle_counts, summary['wf_times'], 'b-o', label='Weighting_factor', linewidth=2, markersize=8)
    ax4.plot(particle_counts, summary['re_times'], 'r-s', label='Random+Ensemble', linewidth=2, markersize=8)
    ax4.set_xlabel('Number of Particles')
    ax4.set_ylabel('Time (s)')
    ax4.set_title('Computational Cost')
    ax4.legend()
    ax4.set_xscale('log')
    ax4.grid(True, alpha=0.3)

    plt.suptitle('Convergence Comparison: Weighting_factor vs Random+Ensemble\n'
                 '(2700s simulation, RH=0.98, w=2.0 m/s, deterministic init, stochastic collision)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.close()


def main():
    """Run the convergence comparison benchmark."""

    print("="*80)
    print("CONVERGENCE COMPARISON BENCHMARK")
    print("Weighting_factor vs Random+Ensemble")
    print("(Deterministic initialization, stochastic collision only)")
    print("="*80)

    # Configuration
    particle_counts = [1000, 5000, 10000, 20000]
    n_trials = 5  # Number of independent trials per configuration
    n_members = 10  # Number of ensemble members for Random approach
    init_seed = 42  # Fixed seed for deterministic initialization

    # Simulation parameters
    T_init = 285.0
    RH_init = 0.98  # Strong cloud development
    P_init = 100000.0
    w_parcel = 2.0  # Strong updraft
    nt = 2700  # Full simulation time
    dt = 1.0

    # Create environmental profiles
    print("\nCreating environmental profiles...")
    qv_init = RH_init * 611.2 * np.exp(17.62 * (T_init - 273.15) / (T_init - 29.65)) \
              * 287.0 / (461.0 * P_init)
    qv_profiles, theta_profiles, z_env = create_env_profiles(
        T_init, qv_init, 0.0, P_init, 'Unstable', plot_profiles=False
    )

    # Common simulation parameters
    sim_params = {
        'P_parcel': P_init,
        'RH_parcel': RH_init,
        'T_parcel': T_init,
        'w_parcel': w_parcel,
        'nt': nt,
        'dt': dt,
        'rm_spec': rm_spec,
        'do_condensation': True,
        'do_collision': True,
        'qv_profiles': qv_profiles,
        'theta_profiles': theta_profiles,
    }

    print(f"\nSimulation configuration:")
    print(f"  - Time: {nt}s (dt={dt}s)")
    print(f"  - RH: {RH_init*100:.0f}%, w: {w_parcel} m/s")
    print(f"  - Collision: enabled (stochastic)")
    print(f"  - Initialization: deterministic (seed={init_seed})")
    print(f"  - Particle counts: {particle_counts}")
    print(f"  - Trials per configuration: {n_trials}")
    print(f"  - Ensemble members (Random): {n_members}")

    # Results storage
    wf_results = {}
    re_results = {}

    total_start = time.time()

    # Run benchmarks for each particle count
    for n_ptcl in particle_counts:
        print(f"\n{'='*80}")
        print(f"TESTING: {n_ptcl} particles")
        print(f"{'='*80}")

        # Weighting_factor trials
        print(f"\n  Weighting_factor ({n_trials} trials):")
        wf_results[n_ptcl] = run_weighting_factor_trials(
            n_ptcl, n_trials, sim_params, init_seed=init_seed
        )

        # Random+Ensemble trials
        print(f"\n  Random+Ensemble ({n_trials} trials × {n_members} members):")
        re_results[n_ptcl] = run_random_ensemble_trials(
            n_ptcl, n_trials, n_members, sim_params, init_seed=init_seed
        )

    total_time = time.time() - total_start
    print(f"\nTotal benchmark time: {total_time/60:.1f} minutes")

    # Analyze and display results
    summary = analyze_results(wf_results, re_results, particle_counts, n_members)

    # Create visualization
    plot_results(summary)

    # Save raw results
    np.savez('convergence_results.npz',
             particle_counts=particle_counts,
             wf_results={k: {kk: np.array(vv) for kk, vv in v.items()}
                        for k, v in wf_results.items()},
             re_results={k: {kk: np.array(vv) for kk, vv in v.items()}
                        for k, v in re_results.items()},
             summary=summary)
    print("\nRaw results saved to: convergence_results.npz")

    return summary


if __name__ == '__main__':
    main()
