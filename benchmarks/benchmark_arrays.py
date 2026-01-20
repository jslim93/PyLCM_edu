"""
Performance benchmarks for the array-based particle system.

Compares performance of:
1. Original particle list implementation
2. New array-based implementation with Numba optimization

Scenarios:
- 1,000 particles, 100 steps
- 5,000 particles, 100 steps
- 10,000 particles, 100 steps
"""

import numpy as np
import time
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PyLCM.micro_particle import ParticleArrays
from PyLCM.parameters import rm_spec, rho_liq, pi


def warm_up_numba():
    """Warm up Numba JIT compilation."""
    print("Warming up Numba JIT compilation...")

    from PyLCM.collision_arrays import ws_drops_beard_batch, apply_collision_arrays
    from PyLCM.condensation_arrays import radius_growth_batch, apply_condensation_arrays

    # Small test to trigger compilation
    particles = ParticleArrays(100)
    particles.M = np.full(100, 4.0/3.0 * pi * rho_liq * (10e-6)**3 * 10.0)
    particles.A = np.full(100, 10.0)
    particles.Ns = np.full(100, 1e-18)
    particles.kappa = np.full(100, 0.6)

    T_parcel = 280.0
    P_parcel = 90000.0
    q_parcel = 0.010
    dt = 0.1

    # Trigger condensation JIT
    apply_condensation_arrays(particles, T_parcel, P_parcel, q_parcel, dt)

    # Trigger collision JIT
    apply_collision_arrays(particles, T_parcel, P_parcel, dt, do_collision=True)

    print("  JIT compilation complete.\n")


def create_test_particles(n_particles):
    """Create test particles with realistic cloud droplet distribution."""
    particles = ParticleArrays(n_particles)

    # Log-normal distribution of radii centered around 10 µm
    radii = np.random.lognormal(np.log(10e-6), 0.5, n_particles)
    radii = np.clip(radii, 1e-6, 500e-6)

    particles.M = 4.0/3.0 * pi * rho_liq * radii**3 * particles.A
    particles.Ns = particles.M * 0.01  # 1% aerosol mass
    particles.kappa = np.random.uniform(0.4, 1.0, n_particles)

    return particles


def benchmark_condensation(n_particles, n_steps):
    """Benchmark condensation performance."""
    from PyLCM.condensation_arrays import apply_condensation_arrays

    particles = create_test_particles(n_particles)

    T_parcel = 280.0
    P_parcel = 90000.0
    q_parcel = 0.010
    dt = 0.1

    start = time.perf_counter()

    for _ in range(n_steps):
        T_parcel, q_parcel, n_act, n_evp = apply_condensation_arrays(
            particles, T_parcel, P_parcel, q_parcel, dt
        )

    elapsed = time.perf_counter() - start

    return elapsed


def benchmark_collision(n_particles, n_steps):
    """Benchmark collision performance."""
    from PyLCM.collision_arrays import apply_collision_arrays

    particles = create_test_particles(n_particles)

    T_parcel = 280.0
    P_parcel = 90000.0
    dt = 1.0

    start = time.perf_counter()

    for _ in range(n_steps):
        n_coll, n_auto, precip = apply_collision_arrays(
            particles, T_parcel, P_parcel, dt, do_collision=True
        )

    elapsed = time.perf_counter() - start

    return elapsed


def benchmark_full_timestep(n_particles, n_steps):
    """Benchmark full timestep routine."""
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays
    from PyLCM.parameters import z_env

    # Create environmental profiles for entrainment/ascent
    theta_profiles = 285.0 + z_env * 0.003  # Stable lapse rate
    qv_profiles = 0.010 * np.exp(-z_env / 2000.0)  # Decreasing humidity with height

    start = time.perf_counter()

    results = timesteps_function_arrays(
        n_particles=n_particles,
        P_parcel=100000.0,
        RH_parcel=0.95,
        T_parcel=285.0,
        w_parcel=1.0,
        nt=n_steps,
        dt=0.1,
        rm_spec=rm_spec,
        do_condensation=True,
        do_collision=False,  # Disable collision for faster benchmark
        theta_profiles=theta_profiles,
        qv_profiles=qv_profiles,
        verbose=False
    )

    elapsed = time.perf_counter() - start

    return elapsed


def run_benchmarks():
    """Run all benchmarks."""
    print("\n" + "="*70)
    print("PyLCM Array-Based Particle System Performance Benchmarks")
    print("="*70)

    # Warm up Numba
    warm_up_numba()

    # Test scenarios
    scenarios = [
        (1000, 100),
        (5000, 100),
        (10000, 100),
    ]

    print("\n" + "-"*70)
    print("Condensation Benchmark (Numba parallel)")
    print("-"*70)
    print(f"{'Particles':>10} {'Steps':>8} {'Time (s)':>12} {'Steps/s':>12} {'Part*Step/s':>15}")
    print("-"*70)

    for n_particles, n_steps in scenarios:
        elapsed = benchmark_condensation(n_particles, n_steps)
        steps_per_sec = n_steps / elapsed
        throughput = n_particles * n_steps / elapsed

        print(f"{n_particles:>10} {n_steps:>8} {elapsed:>12.3f} {steps_per_sec:>12.1f} {throughput:>15.0f}")

    print("\n" + "-"*70)
    print("Collision Benchmark (Numba parallel terminal velocity)")
    print("-"*70)
    print(f"{'Particles':>10} {'Steps':>8} {'Time (s)':>12} {'Steps/s':>12} {'Part*Step/s':>15}")
    print("-"*70)

    for n_particles, n_steps in scenarios:
        elapsed = benchmark_collision(n_particles, n_steps)
        steps_per_sec = n_steps / elapsed
        throughput = n_particles * n_steps / elapsed

        print(f"{n_particles:>10} {n_steps:>8} {elapsed:>12.3f} {steps_per_sec:>12.1f} {throughput:>15.0f}")

    print("\n" + "-"*70)
    print("Full Timestep Benchmark (condensation only)")
    print("-"*70)
    print(f"{'Particles':>10} {'Steps':>8} {'Time (s)':>12} {'Steps/s':>12}")
    print("-"*70)

    for n_particles, n_steps in scenarios:
        elapsed = benchmark_full_timestep(n_particles, n_steps)
        steps_per_sec = n_steps / elapsed

        print(f"{n_particles:>10} {n_steps:>8} {elapsed:>12.3f} {steps_per_sec:>12.1f}")

    print("\n" + "="*70)
    print("Benchmark Complete")
    print("="*70)


def estimate_scaling():
    """Estimate scaling behavior with particle count."""
    print("\n" + "="*70)
    print("Scaling Analysis")
    print("="*70)

    warm_up_numba()

    particle_counts = [100, 500, 1000, 2500, 5000, 7500, 10000]
    n_steps = 50

    print("\nCondensation scaling:")
    print(f"{'Particles':>10} {'Time (s)':>12} {'Time/particle (µs)':>20}")
    print("-"*45)

    times = []
    for n in particle_counts:
        t = benchmark_condensation(n, n_steps)
        time_per_particle = t / n / n_steps * 1e6  # µs
        times.append(t)
        print(f"{n:>10} {t:>12.4f} {time_per_particle:>20.2f}")

    # Estimate scaling exponent
    # If O(n^k), then log(t) = k*log(n) + const
    log_n = np.log(particle_counts)
    log_t = np.log(times)

    # Linear regression
    k, const = np.polyfit(log_n, log_t, 1)
    print(f"\nScaling exponent (condensation): O(n^{k:.2f})")

    print("\nCollision scaling:")
    print(f"{'Particles':>10} {'Time (s)':>12} {'Time/particle (µs)':>20}")
    print("-"*45)

    times_coll = []
    for n in particle_counts:
        t = benchmark_collision(n, n_steps)
        time_per_particle = t / n / n_steps * 1e6  # µs
        times_coll.append(t)
        print(f"{n:>10} {t:>12.4f} {time_per_particle:>20.2f}")

    log_t_coll = np.log(times_coll)
    k_coll, _ = np.polyfit(log_n, log_t_coll, 1)
    print(f"\nScaling exponent (collision): O(n^{k_coll:.2f})")
    print("  Note: Collision should scale as O(n) due to LSM pairing")

    print("\n" + "="*70)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="PyLCM Array Benchmarks")
    parser.add_argument("--scaling", action="store_true", help="Run scaling analysis")
    args = parser.parse_args()

    if args.scaling:
        estimate_scaling()
    else:
        run_benchmarks()
