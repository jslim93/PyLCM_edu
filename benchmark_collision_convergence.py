#!/usr/bin/env python3
"""
Collision Convergence Benchmark: Particle Count vs Ensemble Trade-off

This script determines the optimal trade-off between:
- n_particles (more particles = better sampling)
- n_ensemble (more members = reduced stochastic noise)
- Computational cost

Test Matrix:
- Test 0: LSM vs All-to-All comparison (small particle counts)
- Test 1: Single-run convergence with LSM (10k-100k particles)
- Test 2: Ensemble size effect (at 10k and 50k particles)
- Test 3: Cost-equivalent comparison

Metrics:
- Coefficient of Variation (CV) = std/mean * 100%
- Primary variables: qr (rain water), nr (rain number), precip

Author: Auto-generated for PyLCM_edu
"""

import numpy as np
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing

# Suppress warnings during parallel execution
warnings.filterwarnings('ignore')

from PyLCM.timestep_routine_arrays import timesteps_function_arrays
from PyLCM.parameters import rm_spec, z_env, p0, r_a, cp


# =============================================================================
# Configuration
# =============================================================================

# Create theta profiles for parcel model
T_INIT = 285.0
P_INIT = 100000.0
THETA_INIT = T_INIT * (p0 / P_INIT)**(r_a / cp)
THETA_PROFILES = THETA_INIT + 0.005 * z_env  # Stable lapse rate

# Simulation parameters (collision-dominated scenario)
SIM_CONFIG = {
    'P_parcel': 100000,      # Initial pressure [Pa]
    'RH_parcel': 0.98,       # Initial RH (near saturation)
    'T_parcel': 285,         # Initial temperature [K]
    'w_parcel': 2.0,         # Updraft velocity [m/s]
    'nt': 2700,              # Number of timesteps
    'dt': 1.0,               # Timestep [s]
    'do_collision': True,    # Enable collision (main focus)
    'do_condensation': True, # Enable condensation
    'verbose': False,
    'output_interval': 10,   # Output every 10 steps
    'theta_profiles': THETA_PROFILES,  # Required for parcel ascent
}


def run_single_simulation(n_particles, seed, use_lsm=True):
    """
    Run a single simulation and return key metrics.

    Args:
        n_particles: Number of particles
        seed: Random seed for reproducibility
        use_lsm: Use LSM pairing (True) or all-to-all (False)

    Returns:
        dict: Results including qr, nr, precip, runtime
    """
    start_time = time.time()

    try:
        output = timesteps_function_arrays(
            n_particles=n_particles,
            P_parcel=SIM_CONFIG['P_parcel'],
            RH_parcel=SIM_CONFIG['RH_parcel'],
            T_parcel=SIM_CONFIG['T_parcel'],
            w_parcel=SIM_CONFIG['w_parcel'],
            nt=SIM_CONFIG['nt'],
            dt=SIM_CONFIG['dt'],
            rm_spec=rm_spec,
            do_collision=SIM_CONFIG['do_collision'],
            do_condensation=SIM_CONFIG['do_condensation'],
            verbose=SIM_CONFIG['verbose'],
            output_interval=SIM_CONFIG['output_interval'],
            theta_profiles=SIM_CONFIG['theta_profiles'],
            seed=seed,
            use_lsm=use_lsm
        )

        # Extract final values
        qr_ts = output[9]   # Rain water mixing ratio time series
        nr_ts = output[12]  # Rain number concentration time series
        precip_ts = output[20]  # Cumulative precipitation time series

        runtime = time.time() - start_time

        return {
            'qr_final': qr_ts[-1],
            'nr_final': nr_ts[-1],
            'precip_final': precip_ts[-1],
            'qr_ts': qr_ts,
            'nr_ts': nr_ts,
            'precip_ts': precip_ts,
            'runtime': runtime,
            'success': True
        }

    except Exception as e:
        return {
            'qr_final': np.nan,
            'nr_final': np.nan,
            'precip_final': np.nan,
            'runtime': time.time() - start_time,
            'success': False,
            'error': str(e)
        }


def run_single_wrapper(args):
    """Wrapper for parallel execution."""
    n_particles, seed, use_lsm = args
    return run_single_simulation(n_particles, seed, use_lsm)


def compute_cv(values):
    """Compute coefficient of variation (CV) in percent."""
    values = np.array(values)
    values = values[~np.isnan(values)]  # Remove NaN
    if len(values) < 2 or np.mean(values) == 0:
        return np.nan
    return (np.std(values) / np.mean(values)) * 100


def run_trial_set(n_particles, n_trials, use_lsm=True, base_seed=42, parallel=True):
    """
    Run multiple trials and compute statistics.

    Args:
        n_particles: Number of particles
        n_trials: Number of trials
        use_lsm: Use LSM pairing
        base_seed: Base random seed
        parallel: Use parallel execution

    Returns:
        dict: Statistics including mean, std, CV for each metric
    """
    seeds = [base_seed + i * 1000 for i in range(n_trials)]
    args_list = [(n_particles, seed, use_lsm) for seed in seeds]

    results = []
    if parallel and n_trials > 1:
        n_workers = min(n_trials, multiprocessing.cpu_count())
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(run_single_wrapper, args): args for args in args_list}
            for future in as_completed(futures):
                results.append(future.result())
    else:
        for args in args_list:
            results.append(run_single_wrapper(args))

    # Extract values
    qr_finals = [r['qr_final'] for r in results if r['success']]
    nr_finals = [r['nr_final'] for r in results if r['success']]
    precip_finals = [r['precip_final'] for r in results if r['success']]
    runtimes = [r['runtime'] for r in results if r['success']]

    return {
        'n_particles': n_particles,
        'n_trials': n_trials,
        'use_lsm': use_lsm,
        'qr_mean': np.mean(qr_finals) if qr_finals else np.nan,
        'qr_std': np.std(qr_finals) if qr_finals else np.nan,
        'qr_cv': compute_cv(qr_finals),
        'nr_mean': np.mean(nr_finals) if nr_finals else np.nan,
        'nr_std': np.std(nr_finals) if nr_finals else np.nan,
        'nr_cv': compute_cv(nr_finals),
        'precip_mean': np.mean(precip_finals) if precip_finals else np.nan,
        'precip_std': np.std(precip_finals) if precip_finals else np.nan,
        'precip_cv': compute_cv(precip_finals),
        'runtime_mean': np.mean(runtimes) if runtimes else np.nan,
        'runtime_total': sum(runtimes) if runtimes else np.nan,
        'n_success': len(qr_finals),
        'qr_values': qr_finals,
        'nr_values': nr_finals,
        'precip_values': precip_finals,
    }


def run_ensemble(n_particles, n_ensemble, seed=42, use_lsm=True, parallel=True):
    """
    Run an ensemble of simulations and compute ensemble mean.

    Args:
        n_particles: Number of particles per member
        n_ensemble: Number of ensemble members
        seed: Base random seed
        use_lsm: Use LSM pairing
        parallel: Use parallel execution

    Returns:
        dict: Ensemble-averaged results
    """
    seeds = [seed + i * 1000 for i in range(n_ensemble)]
    args_list = [(n_particles, s, use_lsm) for s in seeds]

    results = []
    if parallel and n_ensemble > 1:
        n_workers = min(n_ensemble, multiprocessing.cpu_count())
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(run_single_wrapper, args): args for args in args_list}
            for future in as_completed(futures):
                results.append(future.result())
    else:
        for args in args_list:
            results.append(run_single_wrapper(args))

    # Compute ensemble mean
    qr_finals = [r['qr_final'] for r in results if r['success']]
    nr_finals = [r['nr_final'] for r in results if r['success']]
    precip_finals = [r['precip_final'] for r in results if r['success']]
    runtimes = [r['runtime'] for r in results if r['success']]

    return {
        'n_particles': n_particles,
        'n_ensemble': n_ensemble,
        'qr_mean': np.mean(qr_finals) if qr_finals else np.nan,
        'nr_mean': np.mean(nr_finals) if nr_finals else np.nan,
        'precip_mean': np.mean(precip_finals) if precip_finals else np.nan,
        'runtime_total': sum(runtimes) if runtimes else np.nan,
        'n_success': len(qr_finals),
    }


# =============================================================================
# Test 0: LSM vs All-to-All Comparison
# =============================================================================

def run_test0_lsm_vs_all_to_all(skip_500_all_to_all=False):
    """
    Test 0: Compare LSM vs All-to-All at small particle counts.

    Purpose: Verify that LSM produces statistically equivalent results to
    all-to-all (mean values should match, CV may differ).
    """
    print("\n" + "="*70)
    print("TEST 0: LSM vs All-to-All Comparison")
    print("="*70)

    configs = [
        # (config_name, n_particles, method, n_trials)
        ('Z1', 100, 'LSM', 20),
        ('Z2', 100, 'All-to-all', 20),
        ('Z3', 200, 'LSM', 20),
        ('Z4', 200, 'All-to-all', 20),
        ('Z5', 500, 'LSM', 20),
    ]

    if not skip_500_all_to_all:
        configs.append(('Z6', 500, 'All-to-all', 20))

    results = []
    for config_name, n_particles, method, n_trials in configs:
        use_lsm = (method == 'LSM')
        print(f"\n{config_name}: {n_particles} particles, {method}, {n_trials} trials...")

        start_time = time.time()
        result = run_trial_set(n_particles, n_trials, use_lsm=use_lsm, parallel=True)
        result['config'] = config_name
        result['method'] = method
        elapsed = time.time() - start_time

        print(f"  qr = {result['qr_mean']*1e3:.4f} ± {result['qr_std']*1e3:.4f} g/kg (CV={result['qr_cv']:.1f}%)")
        print(f"  nr = {result['nr_mean']:.2e} ± {result['nr_std']:.2e} #/kg (CV={result['nr_cv']:.1f}%)")
        print(f"  Time: {elapsed:.1f}s ({result['runtime_mean']:.1f}s/run)")

        results.append(result)

    return results


# =============================================================================
# Test 1: Single-Run Convergence (High Particle Counts)
# =============================================================================

def run_test1_single_run_convergence():
    """
    Test 1: Can high particle counts converge without ensemble?

    Tests 10k, 20k, 50k, 100k particles with 10 trials each.
    """
    print("\n" + "="*70)
    print("TEST 1: Single-Run Convergence with LSM (High Particle Counts)")
    print("="*70)

    configs = [
        # (config_name, n_particles, n_trials)
        ('A1', 10000, 10),
        ('A2', 20000, 10),
        ('A3', 50000, 10),
        ('A4', 100000, 10),
    ]

    results = []
    for config_name, n_particles, n_trials in configs:
        print(f"\n{config_name}: {n_particles:,} particles, {n_trials} trials...")

        start_time = time.time()
        result = run_trial_set(n_particles, n_trials, use_lsm=True, parallel=True)
        result['config'] = config_name
        elapsed = time.time() - start_time

        print(f"  qr = {result['qr_mean']*1e3:.4f} ± {result['qr_std']*1e3:.4f} g/kg (CV={result['qr_cv']:.1f}%)")
        print(f"  nr = {result['nr_mean']:.2e} ± {result['nr_std']:.2e} #/kg (CV={result['nr_cv']:.1f}%)")
        print(f"  Time: {elapsed:.1f}s ({result['runtime_mean']:.1f}s/run)")

        results.append(result)

    return results


# =============================================================================
# Test 2: Ensemble Size Effect
# =============================================================================

def run_test2_ensemble_size_effect(skip_large_ensembles=False):
    """
    Test 2: How many ensemble members are needed at each particle count?

    Tests ensemble sizes 3, 5, 10, 20 at 10k particles
    and 3, 5, 10 at 50k particles.
    """
    print("\n" + "="*70)
    print("TEST 2: Ensemble Size Effect")
    print("="*70)

    # Configurations: (config_name, n_particles, n_ensemble, n_trials)
    configs = [
        ('B1', 10000, 3, 5),
        ('B2', 10000, 5, 5),
        ('B3', 10000, 10, 5),
        ('B5', 50000, 3, 5),
        ('B6', 50000, 5, 5),
    ]

    if not skip_large_ensembles:
        configs.extend([
            ('B4', 10000, 20, 5),
            ('B7', 50000, 10, 5),
        ])

    results = []
    for config_name, n_particles, n_ensemble, n_trials in configs:
        print(f"\n{config_name}: {n_particles:,} particles, {n_ensemble} ensemble members, {n_trials} trials...")

        # Run n_trials ensembles and compute CV of ensemble means
        ensemble_qr_means = []
        ensemble_nr_means = []
        ensemble_precip_means = []
        total_runtime = 0

        for trial in range(n_trials):
            base_seed = 42 + trial * 10000
            ens_result = run_ensemble(n_particles, n_ensemble, seed=base_seed, parallel=True)
            ensemble_qr_means.append(ens_result['qr_mean'])
            ensemble_nr_means.append(ens_result['nr_mean'])
            ensemble_precip_means.append(ens_result['precip_mean'])
            total_runtime += ens_result['runtime_total']

        result = {
            'config': config_name,
            'n_particles': n_particles,
            'n_ensemble': n_ensemble,
            'n_trials': n_trials,
            'qr_mean': np.mean(ensemble_qr_means),
            'qr_std': np.std(ensemble_qr_means),
            'qr_cv': compute_cv(ensemble_qr_means),
            'nr_mean': np.mean(ensemble_nr_means),
            'nr_std': np.std(ensemble_nr_means),
            'nr_cv': compute_cv(ensemble_nr_means),
            'precip_mean': np.mean(ensemble_precip_means),
            'precip_std': np.std(ensemble_precip_means),
            'precip_cv': compute_cv(ensemble_precip_means),
            'runtime_total': total_runtime,
            'qr_values': ensemble_qr_means,
            'nr_values': ensemble_nr_means,
        }

        print(f"  qr = {result['qr_mean']*1e3:.4f} ± {result['qr_std']*1e3:.4f} g/kg (CV={result['qr_cv']:.1f}%)")
        print(f"  nr = {result['nr_mean']:.2e} ± {result['nr_std']:.2e} #/kg (CV={result['nr_cv']:.1f}%)")
        print(f"  Total time: {total_runtime:.1f}s")

        results.append(result)

    return results


# =============================================================================
# Test 3: Cost-Equivalent Comparison
# =============================================================================

def run_test3_cost_equivalent():
    """
    Test 3: Compare configurations with similar computational cost.

    All configs have ~100k particle-equivalents:
    - C1: 100k × 1 run × 5 trials
    - C2: 50k × 2 ensemble × 5 trials
    - C3: 20k × 5 ensemble × 5 trials
    - C4: 10k × 10 ensemble × 5 trials
    """
    print("\n" + "="*70)
    print("TEST 3: Cost-Equivalent Comparison (~100k particle-equivalents)")
    print("="*70)

    # (config_name, n_particles, n_ensemble, n_trials)
    configs = [
        ('C1', 100000, 1, 5),   # Single high-particle runs
        ('C2', 50000, 2, 5),    # 2-member ensemble
        ('C3', 20000, 5, 5),    # 5-member ensemble
        ('C4', 10000, 10, 5),   # 10-member ensemble
    ]

    results = []
    for config_name, n_particles, n_ensemble, n_trials in configs:
        effective_samples = n_ensemble * n_trials
        print(f"\n{config_name}: {n_particles:,} particles × {n_ensemble} ensemble × {n_trials} trials = {effective_samples} samples...")

        if n_ensemble == 1:
            # Single runs (no ensemble averaging)
            start_time = time.time()
            result = run_trial_set(n_particles, n_trials, use_lsm=True, parallel=True)
            result['config'] = config_name
            result['n_ensemble'] = 1
            elapsed = time.time() - start_time
        else:
            # Ensemble runs
            ensemble_qr_means = []
            ensemble_nr_means = []
            ensemble_precip_means = []
            total_runtime = 0

            for trial in range(n_trials):
                base_seed = 42 + trial * 10000
                ens_result = run_ensemble(n_particles, n_ensemble, seed=base_seed, parallel=True)
                ensemble_qr_means.append(ens_result['qr_mean'])
                ensemble_nr_means.append(ens_result['nr_mean'])
                ensemble_precip_means.append(ens_result['precip_mean'])
                total_runtime += ens_result['runtime_total']

            result = {
                'config': config_name,
                'n_particles': n_particles,
                'n_ensemble': n_ensemble,
                'n_trials': n_trials,
                'qr_mean': np.mean(ensemble_qr_means),
                'qr_std': np.std(ensemble_qr_means),
                'qr_cv': compute_cv(ensemble_qr_means),
                'nr_mean': np.mean(ensemble_nr_means),
                'nr_std': np.std(ensemble_nr_means),
                'nr_cv': compute_cv(ensemble_nr_means),
                'precip_mean': np.mean(ensemble_precip_means),
                'precip_std': np.std(ensemble_precip_means),
                'precip_cv': compute_cv(ensemble_precip_means),
                'runtime_total': total_runtime,
                'qr_values': ensemble_qr_means,
                'nr_values': ensemble_nr_means,
            }
            elapsed = total_runtime

        print(f"  qr = {result['qr_mean']*1e3:.4f} ± {result['qr_std']*1e3:.4f} g/kg (CV={result['qr_cv']:.1f}%)")
        print(f"  nr = {result['nr_mean']:.2e} ± {result['nr_std']:.2e} #/kg (CV={result['nr_cv']:.1f}%)")
        print(f"  Time: {elapsed:.1f}s")

        results.append(result)

    return results


# =============================================================================
# Plotting Functions
# =============================================================================

def plot_results(test0_results, test1_results, test2_results, test3_results):
    """Generate summary plots for all tests."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available, skipping plots")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Plot 1: LSM vs All-to-all comparison
    ax1 = axes[0, 0]
    if test0_results:
        lsm_results = [r for r in test0_results if r['method'] == 'LSM']
        all_results = [r for r in test0_results if r['method'] == 'All-to-all']

        particles_lsm = [r['n_particles'] for r in lsm_results]
        cv_lsm = [r['qr_cv'] for r in lsm_results]
        particles_all = [r['n_particles'] for r in all_results]
        cv_all = [r['qr_cv'] for r in all_results]

        ax1.plot(particles_lsm, cv_lsm, 'o-', label='LSM', markersize=8)
        ax1.plot(particles_all, cv_all, 's-', label='All-to-all', markersize=8)
        ax1.axhline(y=5, color='g', linestyle='--', label='Target CV=5%')
        ax1.axhline(y=2, color='b', linestyle=':', label='Target CV=2%')
        ax1.set_xlabel('Number of Particles')
        ax1.set_ylabel('CV of qr (%)')
        ax1.set_title('Test 0: LSM vs All-to-All')
        ax1.legend()
        ax1.grid(True, alpha=0.3)

    # Plot 2: Single-run convergence
    ax2 = axes[0, 1]
    if test1_results:
        particles = [r['n_particles'] for r in test1_results]
        cv_qr = [r['qr_cv'] for r in test1_results]
        cv_nr = [r['nr_cv'] for r in test1_results]

        ax2.semilogx(particles, cv_qr, 'o-', label='qr', markersize=8)
        ax2.semilogx(particles, cv_nr, 's-', label='nr', markersize=8)
        ax2.axhline(y=5, color='g', linestyle='--', label='Target CV=5%')
        ax2.axhline(y=2, color='b', linestyle=':', label='Target CV=2%')
        ax2.set_xlabel('Number of Particles')
        ax2.set_ylabel('CV (%)')
        ax2.set_title('Test 1: Single-Run Convergence (LSM)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    # Plot 3: Ensemble size effect
    ax3 = axes[1, 0]
    if test2_results:
        # Group by particle count
        p10k = [r for r in test2_results if r['n_particles'] == 10000]
        p50k = [r for r in test2_results if r['n_particles'] == 50000]

        if p10k:
            ens_10k = [r['n_ensemble'] for r in p10k]
            cv_10k = [r['qr_cv'] for r in p10k]
            ax3.plot(ens_10k, cv_10k, 'o-', label='10k particles', markersize=8)

        if p50k:
            ens_50k = [r['n_ensemble'] for r in p50k]
            cv_50k = [r['qr_cv'] for r in p50k]
            ax3.plot(ens_50k, cv_50k, 's-', label='50k particles', markersize=8)

        ax3.axhline(y=5, color='g', linestyle='--', label='Target CV=5%')
        ax3.axhline(y=2, color='b', linestyle=':', label='Target CV=2%')
        ax3.set_xlabel('Ensemble Size')
        ax3.set_ylabel('CV of qr (%)')
        ax3.set_title('Test 2: Ensemble Size Effect')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

    # Plot 4: Cost-equivalent comparison
    ax4 = axes[1, 1]
    if test3_results:
        configs = [r['config'] for r in test3_results]
        cv_qr = [r['qr_cv'] for r in test3_results]
        runtime = [r.get('runtime_total', 0) for r in test3_results]

        # Bar plot for CV
        bars = ax4.bar(configs, cv_qr, alpha=0.7, label='CV of qr')
        ax4.axhline(y=5, color='g', linestyle='--', label='Target CV=5%')
        ax4.set_ylabel('CV of qr (%)')
        ax4.set_title('Test 3: Cost-Equivalent Comparison')

        # Add runtime as secondary axis
        ax4_twin = ax4.twinx()
        ax4_twin.plot(configs, runtime, 'ro-', label='Runtime')
        ax4_twin.set_ylabel('Runtime (s)', color='r')
        ax4_twin.tick_params(axis='y', labelcolor='r')

        ax4.legend(loc='upper left')
        ax4.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plt.savefig('collision_convergence_summary.png', dpi=150, bbox_inches='tight')
    print("\nSaved: collision_convergence_summary.png")
    plt.close()


def save_results(test0, test1, test2, test3, filename='collision_convergence_results.npz'):
    """Save all results to a compressed numpy file."""

    def results_to_dict(results, prefix):
        """Convert list of result dicts to flat dict for saving."""
        if not results:
            return {}

        out = {}
        # Save scalars
        for key in ['qr_mean', 'qr_std', 'qr_cv', 'nr_mean', 'nr_std', 'nr_cv',
                    'precip_mean', 'precip_std', 'precip_cv', 'runtime_total',
                    'n_particles', 'n_trials', 'n_ensemble']:
            values = [r.get(key, np.nan) for r in results]
            if any(v is not None for v in values):
                out[f'{prefix}_{key}'] = np.array(values, dtype=float)

        # Save config names
        configs = [r.get('config', '') for r in results]
        out[f'{prefix}_configs'] = np.array(configs, dtype=str)

        return out

    data = {}
    data.update(results_to_dict(test0, 'test0'))
    data.update(results_to_dict(test1, 'test1'))
    data.update(results_to_dict(test2, 'test2'))
    data.update(results_to_dict(test3, 'test3'))

    np.savez_compressed(filename, **data)
    print(f"\nSaved: {filename}")


def print_summary(test0, test1, test2, test3):
    """Print a summary table of all results."""

    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)

    # Test 0 summary
    if test0:
        print("\n--- Test 0: LSM vs All-to-All ---")
        print(f"{'Config':<8} {'Method':<12} {'Particles':>10} {'qr CV':>10} {'nr CV':>10}")
        print("-" * 54)
        for r in test0:
            print(f"{r['config']:<8} {r['method']:<12} {r['n_particles']:>10,} {r['qr_cv']:>9.1f}% {r['nr_cv']:>9.1f}%")

    # Test 1 summary
    if test1:
        print("\n--- Test 1: Single-Run Convergence ---")
        print(f"{'Config':<8} {'Particles':>12} {'qr CV':>10} {'nr CV':>10} {'Runtime':>10}")
        print("-" * 54)
        for r in test1:
            print(f"{r['config']:<8} {r['n_particles']:>12,} {r['qr_cv']:>9.1f}% {r['nr_cv']:>9.1f}% {r['runtime_mean']:>9.1f}s")

    # Test 2 summary
    if test2:
        print("\n--- Test 2: Ensemble Size Effect ---")
        print(f"{'Config':<8} {'Particles':>10} {'Ensemble':>10} {'qr CV':>10} {'nr CV':>10}")
        print("-" * 52)
        for r in test2:
            print(f"{r['config']:<8} {r['n_particles']:>10,} {r['n_ensemble']:>10} {r['qr_cv']:>9.1f}% {r['nr_cv']:>9.1f}%")

    # Test 3 summary
    if test3:
        print("\n--- Test 3: Cost-Equivalent Comparison ---")
        print(f"{'Config':<8} {'Particles':>10} {'Ensemble':>8} {'qr CV':>10} {'Runtime':>10}")
        print("-" * 50)
        for r in test3:
            print(f"{r['config']:<8} {r['n_particles']:>10,} {r['n_ensemble']:>8} {r['qr_cv']:>9.1f}% {r.get('runtime_total', 0):>9.1f}s")

    # Recommendations
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)

    if test1:
        # Find particle count where CV < 5%
        converged = [r for r in test1 if r['qr_cv'] < 5]
        if converged:
            min_particles = min(r['n_particles'] for r in converged)
            print(f"• Single-run CV < 5%: achieved at {min_particles:,} particles")
        else:
            print("• Single-run CV < 5%: NOT achieved (even at 100k particles)")

    if test2:
        # Find minimum ensemble for CV < 5%
        converged = [r for r in test2 if r['qr_cv'] < 5]
        if converged:
            best = min(converged, key=lambda r: r['n_ensemble'])
            print(f"• Ensemble CV < 5%: {best['n_particles']:,} particles with {best['n_ensemble']} members")

    if test3:
        # Find best cost-equivalent config
        best = min(test3, key=lambda r: r['qr_cv'])
        print(f"• Best cost-equivalent: {best['config']} ({best['n_particles']:,} × {best['n_ensemble']} ensemble, CV={best['qr_cv']:.1f}%)")


# =============================================================================
# Main
# =============================================================================

def main(reduced=False):
    """
    Run all convergence tests.

    Args:
        reduced: If True, skip some expensive configurations to save time.
    """
    print("="*70)
    print("COLLISION CONVERGENCE BENCHMARK")
    print("Particle Count vs Ensemble Trade-off Analysis")
    print("="*70)
    print(f"\nSimulation: {SIM_CONFIG['nt']} timesteps × {SIM_CONFIG['dt']}s = {SIM_CONFIG['nt']*SIM_CONFIG['dt']/60:.1f} min simulation time")
    print(f"Reduced mode: {reduced}")

    total_start = time.time()

    # Run tests
    test0_results = run_test0_lsm_vs_all_to_all(skip_500_all_to_all=reduced)
    test1_results = run_test1_single_run_convergence()
    test2_results = run_test2_ensemble_size_effect(skip_large_ensembles=reduced)
    test3_results = run_test3_cost_equivalent()

    total_time = time.time() - total_start

    # Save and summarize
    save_results(test0_results, test1_results, test2_results, test3_results)
    plot_results(test0_results, test1_results, test2_results, test3_results)
    print_summary(test0_results, test1_results, test2_results, test3_results)

    print(f"\n{'='*70}")
    print(f"TOTAL BENCHMARK TIME: {total_time/60:.1f} minutes ({total_time:.0f}s)")
    print("="*70)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Collision Convergence Benchmark')
    parser.add_argument('--reduced', action='store_true',
                       help='Run reduced test set (skip expensive configs)')
    args = parser.parse_args()

    main(reduced=args.reduced)
