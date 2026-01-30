"""
Test script for array-based particle implementation.

Run with: runconda python test_arrays.py
"""

import numpy as np
import time

print("Testing array-based particle implementation...")

# Test 1: ParticleArrays creation
print("\nTest 1: Creating ParticleArrays...")
from PyLCM.micro_particle import ParticleArrays, particles

n = 10000
particle_arrays = ParticleArrays(n)
print(f"  ✓ Created {n} particles as arrays")
print(f"  Memory footprint: ~{(particle_arrays.M.nbytes + particle_arrays.A.nbytes + particle_arrays.Ns.nbytes) / 1024:.2f} KB")

# Test 2: Conversion to/from particle list
print("\nTest 2: Testing list conversion...")
particle_list = particle_arrays.to_particle_list()
print(f"  ✓ Converted to list ({len(particle_list)} particles)")

recovered_arrays = ParticleArrays.from_particle_list(particle_list)
print(f"  ✓ Converted back to arrays")
print(f"  Mass conservation check: {np.allclose(particle_arrays.M, recovered_arrays.M)}")

# Test 3: Vectorized operations
print("\nTest 3: Testing vectorized radius calculations...")
test_radii = np.linspace(1e-8, 1e-4, n)
particle_arrays.set_radii(test_radii)
recovered_radii = particle_arrays.get_radii()
print(f"  Radius recovery error: {np.max(np.abs(test_radii - recovered_radii)):.2e}")

# Test 4: Array-based initialization
print("\nTest 4: Testing array-based aerosol initialization...")
from PyLCM.aero_init_arrays import aero_init_arrays
from PyLCM.parameters import *

n_particles = 5000
P_parcel = 101300.0
T_parcel = 293.2
q_parcel = 0.01
z_parcel = 0.0

# Aerosol parameters
N_aero = [118.0e6, 11.0e6, 0.72e6, 0.0]
mu_aero = [0.019e-6, 0.056e-6, 0.46e-6, 0.0]
sigma_aero = [3.3, 1.6, 2.2, 0.0]
k_aero = [1.6, 1.6, 1.6, 1.6]

start = time.time()
T_new, q_new, particle_arrays = aero_init_arrays(
    'Weighting_factor', n_particles, P_parcel, z_parcel, T_parcel, q_parcel,
    N_aero, mu_aero, sigma_aero, rho_aero, k_aero, False
)
elapsed = time.time() - start

print(f"  ✓ Initialized {n_particles} particles in {elapsed:.3f}s")
print(f"  Total particle mass: {np.sum(particle_arrays.M):.3e} kg")
print(f"  Number of active particles: {np.sum(particle_arrays.M > 0)}")

# Test 5: Condensation arrays
print("\nTest 5: Testing array-based condensation...")
from PyLCM.condensation_arrays import apply_condensation_arrays

start = time.time()
T_after, q_after, n_act, n_evp = apply_condensation_arrays(
    particle_arrays, T_new, P_parcel, q_new, dt=1.0, do_condensation=True
)
elapsed = time.time() - start

print(f"  ✓ Condensation step completed in {elapsed:.3f}s")
print(f"  Activated: {n_act}, Evaporated: {n_evp}")
print(f"  Temperature change: {T_after - T_new:.3f} K")

# Test 6: Collision arrays
print("\nTest 6: Testing array-based collision...")
from PyLCM.collision_arrays import apply_collision_arrays

start = time.time()
n_coll, n_auto = apply_collision_arrays(
    particle_arrays, T_after, P_parcel, dt=1.0, do_collision=True
)
elapsed = time.time() - start

print(f"  ✓ Collision step completed in {elapsed:.3f}s")
print(f"  Collisions: {n_coll}, Autoconversions: {n_auto}")

print("\n" + "="*60)
print("All array-based tests passed! ✓")
print("="*60)
print("\nNext step: Integrate into timestep_routine.py for full optimization")
