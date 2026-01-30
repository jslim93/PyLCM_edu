"""
Benchmarking script to compare original vs optimized performance.

Usage:
    python benchmark_optimizations.py
"""
import time
import numpy as np
from PyLCM import collision, condensation
from PyLCM import collision_optimized, condensation_optimized


def benchmark_collision_efficiency(n_tests=10000):
    """Benchmark Hall efficiency calculation."""
    print("Benchmarking collision efficiency (Hall 1980)...")
    
    # Generate random droplet radii
    np.random.seed(42)
    r1_vals = np.random.uniform(10e-6, 100e-6, n_tests)
    r2_vals = np.random.uniform(10e-6, 100e-6, n_tests)
    
    # Original version
    start = time.time()
    for r1, r2 in zip(r1_vals, r2_vals):
        _ = collision.E_H80(r1, r2)
    time_original = time.time() - start
    
    # Optimized version
    start = time.time()
    for r1, r2 in zip(r1_vals, r2_vals):
        _ = collision_optimized.E_H80_optimized(r1, r2)
    time_optimized = time.time() - start
    
    speedup = time_original / time_optimized
    print(f"  Original: {time_original:.3f}s")
    print(f"  Optimized: {time_optimized:.3f}s")
    print(f"  Speedup: {speedup:.1f}x\n")
    
    return speedup


def benchmark_terminal_velocity(n_tests=10000):
    """Benchmark Beard terminal velocity calculation."""
    print("Benchmarking terminal velocity (Beard 1976)...")
    
    np.random.seed(42)
    radii = np.random.uniform(1e-6, 1000e-6, n_tests)
    rho_parcel = 1.2
    rho_liq = 1000.0
    p_env = 101325.0
    T_parcel = 293.15
    
    # Original version
    start = time.time()
    for r in radii:
        _ = collision.ws_drops_beard(r, rho_parcel, rho_liq, p_env, T_parcel)
    time_original = time.time() - start
    
    # Optimized version  
    start = time.time()
    for r in radii:
        _ = collision_optimized.ws_drops_beard_optimized(r, rho_parcel, rho_liq, p_env, T_parcel)
    time_optimized = time.time() - start
    
    speedup = time_original / time_optimized
    print(f"  Original: {time_original:.3f}s")
    print(f"  Optimized: {time_optimized:.3f}s")
    print(f"  Speedup: {speedup:.1f}x\n")
    
    return speedup


def benchmark_saturation_vapor_pressure(n_tests=100000):
    """Benchmark saturation vapor pressure calculation."""
    print("Benchmarking saturation vapor pressure (Flatau et al.)...")
    
    np.random.seed(42)
    temps = np.random.uniform(250, 310, n_tests)
    
    # Original version
    start = time.time()
    for T in temps:
        _ = condensation.esatw(T)
    time_original = time.time() - start
    
    # Optimized version
    start = time.time()
    for T in temps:
        _ = condensation_optimized.esatw_optimized(T)
    time_optimized = time.time() - start
    
    speedup = time_original / time_optimized
    print(f"  Original: {time_original:.3f}s")
    print(f"  Optimized: {time_optimized:.3f}s")
    print(f"  Speedup: {speedup:.1f}x\n")
    
    return speedup


def verify_numerical_accuracy():
    """Verify that optimized versions produce same results as original."""
    print("Verifying numerical accuracy...")
    
    # Test collision efficiency
    r1, r2 = 50e-6, 30e-6
    orig = collision.E_H80(r1, r2)
    opt = collision_optimized.E_H80_optimized(r1, r2)
    err = abs(orig - opt) / abs(orig) if orig != 0 else abs(orig - opt)
    print(f"  E_H80 relative error: {err:.2e}")
    assert err < 1e-10, "Collision efficiency mismatch!"
    
    # Test terminal velocity
    orig = collision.ws_drops_beard(100e-6, 1.2, 1000.0, 101325.0, 293.15)
    opt = collision_optimized.ws_drops_beard_optimized(100e-6, 1.2, 1000.0, 101325.0, 293.15)
    err = abs(orig - opt) / abs(orig)
    print(f"  ws_drops_beard relative error: {err:.2e}")
    assert err < 1e-10, "Terminal velocity mismatch!"
    
    # Test saturation vapor pressure
    T =285.0
    orig = condensation.esatw(T)
    opt = condensation_optimized.esatw_optimized(T)
    err = abs(orig - opt) / abs(orig)
    print(f"  esatw relative error: {err:.2e}")
    assert err < 1e-10, "Saturation vapor pressure mismatch!"
    
    print("  ✓ All numerical tests passed!\n")


if __name__ == "__main__":
    print("="*60)
    print("PyLCM Performance Benchmark")
    print("="*60 + "\n")
    
    # Verify accuracy first
    verify_numerical_accuracy()
    
    # Run benchmarks
    speedups = []
    speedups.append(benchmark_collision_efficiency())
    speedups.append(benchmark_terminal_velocity())
    speedups.append(benchmark_saturation_vapor_pressure())
    
    # Summary
    print("="*60)
    print("Summary")
    print("="*60)
    print(f"Average speedup: {np.mean(speedups):.1f}x")
    print(f"Total performance gain: {(np.mean(speedups)-1)*100:.0f}%")
    print("\n✓ Optimization successful!")
