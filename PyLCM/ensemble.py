"""
Parallel ensemble runner for PyLCM simulations.

This module provides functionality to run multiple stochastic realizations
in parallel and aggregate statistics for converged results.

The superdroplet method involves stochastic processes:
- Random aerosol sampling (lognormal distribution)
- Random collision pairing (LSM shuffle)
- Random collision decisions (probability threshold)

Ensemble averaging provides:
- Converged mean values
- Uncertainty quantification (spread)
- Statistically robust results
"""

import numpy as np
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import os
import sys


class EnsembleResult:
    """
    Container for ensemble simulation statistics.

    Aggregates results from multiple stochastic realizations and computes
    statistics (mean, std, percentiles) for key variables.

    Attributes:
        n_members (int): Number of ensemble members
        time_array (np.ndarray): Time array (same for all members)
        qc_all (np.ndarray): Cloud water for all members, shape (n_members, nt+1)
        qr_all (np.ndarray): Rain water for all members
        nc_all (np.ndarray): Cloud number for all members
        nr_all (np.ndarray): Rain number for all members
        qc_mean, qc_std, qc_p10, qc_p25, qc_p75, qc_p90: Cloud water statistics
        (similar for qr, nc, nr, rc_liq_avg, precip)
    """

    def __init__(self, outputs, n_members):
        """
        Initialize EnsembleResult from list of simulation outputs.

        Args:
            outputs: List of output tuples from timesteps_function_arrays
            n_members (int): Number of ensemble members
        """
        self.n_members = n_members

        # Time array is same for all members (index 2 in output tuple)
        self.time_array = outputs[0][2]
        self.nt = outputs[0][0]
        self.dt = outputs[0][1]

        # Stack arrays for each variable: shape (n_members, nt+1)
        # Output indices from timesteps_function_arrays:
        # 0: nt, 1: dt, 2: time_array, 3: T_parcel_array, 4: RH_parcel_array,
        # 5: q_parcel_array, 6: z_parcel_array, 7: qa_ts, 8: qc_ts, 9: qr_ts,
        # 10: na_ts, 11: nc_ts, 12: nr_ts, 13: spectra_arr, 14: con_ts,
        # 15: act_ts, 16: evp_ts, 17: dea_ts, 18: acc_ts, 19: aut_ts,
        # 20: precip_ts, 21: particles_array, 22: rc_liq_avg_array, 23: rc_liq_std_array

        # Parcel state
        self.T_all = np.stack([o[3] for o in outputs])
        self.RH_all = np.stack([o[4] for o in outputs])
        self.q_all = np.stack([o[5] for o in outputs])
        self.z_all = np.stack([o[6] for o in outputs])

        # Mixing ratios
        self.qa_all = np.stack([o[7] for o in outputs])   # Aerosol water
        self.qc_all = np.stack([o[8] for o in outputs])   # Cloud water
        self.qr_all = np.stack([o[9] for o in outputs])   # Rain water

        # Number concentrations
        self.na_all = np.stack([o[10] for o in outputs])  # Aerosol number
        self.nc_all = np.stack([o[11] for o in outputs])  # Cloud number
        self.nr_all = np.stack([o[12] for o in outputs])  # Rain number

        # Spectra (shape: n_members, nt+1, n_bins)
        self.spectra_all = np.stack([o[13] for o in outputs])

        # Process rates
        self.con_all = np.stack([o[14] for o in outputs])     # Condensation
        self.act_all = np.stack([o[15] for o in outputs])     # Activation
        self.evp_all = np.stack([o[16] for o in outputs])     # Evaporation
        self.acc_all = np.stack([o[18] for o in outputs])     # Accretion
        self.aut_all = np.stack([o[19] for o in outputs])     # Autoconversion
        self.precip_all = np.stack([o[20] for o in outputs])  # Precipitation

        # Particle statistics
        self.particles_all = np.stack([o[21] for o in outputs])
        self.rc_liq_avg_all = np.stack([o[22] for o in outputs])
        self.rc_liq_std_all = np.stack([o[23] for o in outputs])

        # Compute statistics for all variables
        self._compute_statistics()

    def _compute_statistics(self):
        """Compute mean, std, and percentiles for all variables."""
        # Cloud water
        self.qc_mean = np.mean(self.qc_all, axis=0)
        self.qc_std = np.std(self.qc_all, axis=0)
        self.qc_p10 = np.percentile(self.qc_all, 10, axis=0)
        self.qc_p25 = np.percentile(self.qc_all, 25, axis=0)
        self.qc_p75 = np.percentile(self.qc_all, 75, axis=0)
        self.qc_p90 = np.percentile(self.qc_all, 90, axis=0)

        # Rain water
        self.qr_mean = np.mean(self.qr_all, axis=0)
        self.qr_std = np.std(self.qr_all, axis=0)
        self.qr_p10 = np.percentile(self.qr_all, 10, axis=0)
        self.qr_p25 = np.percentile(self.qr_all, 25, axis=0)
        self.qr_p75 = np.percentile(self.qr_all, 75, axis=0)
        self.qr_p90 = np.percentile(self.qr_all, 90, axis=0)

        # Cloud number
        self.nc_mean = np.mean(self.nc_all, axis=0)
        self.nc_std = np.std(self.nc_all, axis=0)

        # Rain number
        self.nr_mean = np.mean(self.nr_all, axis=0)
        self.nr_std = np.std(self.nr_all, axis=0)

        # Mean cloud droplet radius
        self.rc_liq_avg_mean = np.mean(self.rc_liq_avg_all, axis=0)
        self.rc_liq_avg_std = np.std(self.rc_liq_avg_all, axis=0)

        # Precipitation
        self.precip_mean = np.mean(self.precip_all, axis=0)
        self.precip_std = np.std(self.precip_all, axis=0)
        self.precip_p10 = np.percentile(self.precip_all, 10, axis=0)
        self.precip_p25 = np.percentile(self.precip_all, 25, axis=0)
        self.precip_p75 = np.percentile(self.precip_all, 75, axis=0)
        self.precip_p90 = np.percentile(self.precip_all, 90, axis=0)

        # Mean spectrum
        self.spectra_mean = np.mean(self.spectra_all, axis=0)
        self.spectra_std = np.std(self.spectra_all, axis=0)

        # Parcel state (usually deterministic, but include for completeness)
        self.T_mean = np.mean(self.T_all, axis=0)
        self.RH_mean = np.mean(self.RH_all, axis=0)
        self.z_mean = np.mean(self.z_all, axis=0)

    def to_csv(self, prefix):
        """
        Save ensemble results to CSV files.

        Args:
            prefix (str): File prefix for output files
        """
        import csv

        # Save time series statistics
        header = ['time', 'qc_mean', 'qc_std', 'qc_p10', 'qc_p25', 'qc_p75', 'qc_p90',
                  'qr_mean', 'qr_std', 'nc_mean', 'nc_std', 'nr_mean', 'nr_std',
                  'precip_mean', 'precip_std', 'rc_liq_avg_mean', 'rc_liq_avg_std']

        with open(f'{prefix}_timeseries.csv', 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for i in range(len(self.time_array)):
                writer.writerow([
                    self.time_array[i],
                    self.qc_mean[i], self.qc_std[i],
                    self.qc_p10[i], self.qc_p25[i], self.qc_p75[i], self.qc_p90[i],
                    self.qr_mean[i], self.qr_std[i],
                    self.nc_mean[i], self.nc_std[i],
                    self.nr_mean[i], self.nr_std[i],
                    self.precip_mean[i], self.precip_std[i],
                    self.rc_liq_avg_mean[i], self.rc_liq_avg_std[i]
                ])

        # Save all member data for qc (for convergence analysis)
        np.savetxt(f'{prefix}_qc_all.csv', self.qc_all, delimiter=',',
                   header=','.join([f't{i}' for i in range(self.qc_all.shape[1])]))

        print(f"Results saved to {prefix}_timeseries.csv and {prefix}_qc_all.csv")

    def summary(self):
        """Print summary statistics."""
        print(f"Ensemble Results Summary (n={self.n_members} members)")
        print(f"{'='*50}")
        print(f"Time range: {self.time_array[0]:.1f} - {self.time_array[-1]:.1f} s")
        print(f"\nFinal values (t = {self.time_array[-1]:.1f} s):")
        print(f"  Cloud water:  {self.qc_mean[-1]*1000:.4f} ± {self.qc_std[-1]*1000:.4f} g/kg")
        print(f"  Rain water:   {self.qr_mean[-1]*1000:.4f} ± {self.qr_std[-1]*1000:.4f} g/kg")
        print(f"  Cloud number: {self.nc_mean[-1]/1e6:.2f} ± {self.nc_std[-1]/1e6:.2f} cm⁻³")
        print(f"  Rain number:  {self.nr_mean[-1]:.2f} ± {self.nr_std[-1]:.2f} m⁻³")
        print(f"  Mean radius:  {self.rc_liq_avg_mean[-1]*1e6:.2f} ± {self.rc_liq_avg_std[-1]*1e6:.2f} µm")


def _run_single_member(args):
    """
    Worker function for parallel execution.

    This function runs in a separate process with isolated random state.

    Args:
        args: Tuple of (init_seed, collision_seed, kwargs) where kwargs are
              arguments for timesteps_function_arrays.
              - init_seed: Seed for particle initialization (can be None for random)
              - collision_seed: Seed for collision process (can be None for random)

    Returns:
        Output tuple from timesteps_function_arrays
    """
    init_seed, collision_seed, kwargs = args

    # Import here to ensure each process has its own module state
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays

    # Run simulation with specified seeds
    return timesteps_function_arrays(**kwargs, seed=init_seed, collision_seed=collision_seed)


def run_ensemble(n_members, seed_base=0, init_seed=None, n_workers=None, method='fork', **kwargs):
    """
    Run ensemble of simulations in parallel.

    Uses multiprocessing for true parallelism across CPU cores.
    Each ensemble member gets a unique collision seed (seed_base + i).

    Args:
        n_members (int): Number of ensemble members to run
        seed_base (int): Base seed value for collision. Member i gets
                        collision_seed = seed_base + i.
                        Use different seed_base values for independent ensembles.
        init_seed (int, optional): Fixed seed for particle initialization.
                        If provided, all members start with identical particles.
                        If None, each member uses seed_base + i for both init
                        and collision (original behavior for backward compatibility).
        n_workers (int, optional): Number of parallel workers.
                                   Defaults to CPU count.
        method (str): Multiprocessing start method ('fork', 'spawn', or 'serial').
                     'fork' is faster but only works on Unix.
                     'spawn' is safer but has more overhead.
                     'serial' runs sequentially (for debugging).
        **kwargs: Arguments passed to timesteps_function_arrays
                  (n_particles, P_parcel, RH_parcel, T_parcel, w_parcel,
                   nt, dt, rm_spec, etc.)

    Returns:
        EnsembleResult: Object containing aggregated statistics

    Example:
        >>> from PyLCM.ensemble import run_ensemble
        >>> from PyLCM.parameters import rm_spec
        >>> # Deterministic init, stochastic collision:
        >>> result = run_ensemble(
        ...     n_members=20,
        ...     n_workers=4,
        ...     init_seed=42,      # Same initialization for all members
        ...     seed_base=0,       # Collision seeds: 0, 1, 2, ...
        ...     n_particles=1000,
        ...     P_parcel=100000.0,
        ...     RH_parcel=0.95,
        ...     T_parcel=285.0,
        ...     w_parcel=1.0,
        ...     nt=100,
        ...     dt=1.0,
        ...     rm_spec=rm_spec,
        ...     do_condensation=True,
        ...     do_collision=True,
        ...     verbose=False
        ... )
        >>> print(f"qc = {result.qc_mean[-1]:.4f} ± {result.qc_std[-1]:.4f}")
    """
    if n_workers is None:
        n_workers = os.cpu_count() or 1

    # Disable verbose output for parallel runs (would interleave)
    kwargs['verbose'] = False
    kwargs['progress_callback'] = None

    # Create argument tuples: (init_seed, collision_seed, kwargs) for each member
    # If init_seed is provided: deterministic init, stochastic collision
    # If init_seed is None: both init and collision use seed_base + i (backward compatible)
    if init_seed is not None:
        args_list = [(init_seed, seed_base + i, kwargs) for i in range(n_members)]
    else:
        # Backward compatibility: use same seed for both init and collision
        args_list = [(seed_base + i, None, kwargs) for i in range(n_members)]

    if method == 'serial':
        print(f"Running {n_members} ensemble members serially...")
        outputs = [_run_single_member(args) for args in args_list]
    else:
        print(f"Running {n_members} ensemble members on {n_workers} workers (method={method})...")

        # Set environment variable to allow fork on macOS
        # This is needed because matplotlib/scipy may initialize ObjC
        if sys.platform == 'darwin' and method == 'fork':
            os.environ['OBJC_DISABLE_INITIALIZE_FORK_SAFETY'] = 'YES'

        # Use fork method on Unix for better performance (avoids reimporting)
        if method == 'fork' and sys.platform != 'win32':
            ctx = mp.get_context('fork')
        else:
            ctx = mp.get_context('spawn')

        with ctx.Pool(processes=n_workers) as pool:
            outputs = pool.map(_run_single_member, args_list)

    print(f"Ensemble complete. Aggregating results...")

    return EnsembleResult(outputs, n_members)


def run_ensemble_serial(n_members, seed_base=0, init_seed=None, **kwargs):
    """
    Run ensemble of simulations serially (for debugging/testing).

    Same interface as run_ensemble but runs sequentially.

    Args:
        n_members (int): Number of ensemble members
        seed_base (int): Base seed value for collision
        init_seed (int, optional): Fixed seed for particle initialization.
                        If provided, all members start with identical particles.
                        If None, uses seed_base + i for both (backward compatible).
        **kwargs: Arguments for timesteps_function_arrays

    Returns:
        EnsembleResult: Object containing aggregated statistics
    """
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays

    kwargs['verbose'] = False
    kwargs['progress_callback'] = None

    outputs = []
    for i in range(n_members):
        collision_seed = seed_base + i
        if init_seed is not None:
            # Deterministic init, stochastic collision
            output = timesteps_function_arrays(**kwargs, seed=init_seed, collision_seed=collision_seed)
        else:
            # Backward compatibility: single seed for both
            output = timesteps_function_arrays(**kwargs, seed=collision_seed)
        outputs.append(output)
        print(f"  Completed member {i+1}/{n_members}")

    return EnsembleResult(outputs, n_members)


def check_convergence(ensemble_result, variable='qc', threshold=0.01):
    """
    Check if ensemble has converged using running mean stability.

    Convergence criterion: The standard deviation of running means
    over the last 10 members should be less than threshold times
    the overall mean.

    This tests whether adding more ensemble members would significantly
    change the estimated mean.

    Args:
        ensemble_result (EnsembleResult): Ensemble results object
        variable (str): Variable to check ('qc', 'qr', 'nc', 'nr', 'precip')
        threshold (float): Relative error threshold (default 0.01 = 1%)

    Returns:
        bool: True if converged, False otherwise

    Example:
        >>> is_converged = check_convergence(result, variable='qc', threshold=0.01)
        >>> if not is_converged:
        ...     print("Consider running more ensemble members")
    """
    # Get the all-members array for the variable
    values = getattr(ensemble_result, f'{variable}_all')  # shape: (n_members, nt+1)
    n_members = values.shape[0]

    if n_members < 10:
        print(f"Warning: Only {n_members} members. Need at least 10 for convergence check.")
        return False

    # Compute running means: cumulative sum / count
    # running_means[i] = mean of members 0..i
    cumsum = np.cumsum(values, axis=0)
    counts = np.arange(1, n_members + 1)[:, None]
    running_means = cumsum / counts  # shape: (n_members, nt+1)

    # Look at stability of running mean over last 10 members
    last_10_means = running_means[-10:]  # shape: (10, nt+1)

    # Standard deviation of running means (how much it's still changing)
    final_std = np.std(last_10_means, axis=0)

    # Overall mean for comparison
    final_mean = np.mean(values, axis=0)

    # Relative error (avoid division by zero)
    rel_error = final_std / (np.abs(final_mean) + 1e-10)

    # Maximum relative error across all timesteps
    max_rel_error = np.max(rel_error)

    is_converged = max_rel_error < threshold

    print(f"Convergence check for {variable}:")
    print(f"  Max relative error: {max_rel_error:.4f} (threshold: {threshold})")
    print(f"  Converged: {is_converged}")

    return is_converged


def estimate_uncertainty_reduction(ensemble_result, variable='qc'):
    """
    Estimate how uncertainty decreases with ensemble size.

    Standard error of the mean scales as 1/sqrt(n), so this function
    computes the actual scaling factor to verify convergence behavior.

    Args:
        ensemble_result (EnsembleResult): Ensemble results
        variable (str): Variable to analyze

    Returns:
        dict: Dictionary with scaling analysis results
    """
    values = getattr(ensemble_result, f'{variable}_all')
    n_members = values.shape[0]

    # Compute standard error at different ensemble sizes
    sizes = np.arange(5, n_members + 1, max(1, n_members // 10))
    stderr_values = []

    for n in sizes:
        subset = values[:n]
        stderr = np.std(subset, axis=0) / np.sqrt(n)
        stderr_values.append(np.mean(stderr))  # Average across time

    stderr_values = np.array(stderr_values)

    # Fit 1/sqrt(n) scaling
    # stderr ~ c / sqrt(n) => log(stderr) ~ log(c) - 0.5*log(n)
    log_n = np.log(sizes)
    log_stderr = np.log(stderr_values + 1e-20)

    # Linear regression
    slope, intercept = np.polyfit(log_n, log_stderr, 1)

    return {
        'sizes': sizes,
        'stderr': stderr_values,
        'expected_slope': -0.5,
        'actual_slope': slope,
        'scaling_ok': abs(slope - (-0.5)) < 0.2
    }
