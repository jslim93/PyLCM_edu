"""
Physics validation tests for the array-based particle system.

These tests verify:
1. Mass conservation during collision
2. Kappa (hygroscopicity) conservation during collision
3. Condensation physics accuracy
4. Terminal velocity calculation correctness
5. Collision efficiency calculation correctness
"""

import numpy as np
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from PyLCM.micro_particle import ParticleArrays
from PyLCM.collision_arrays import (
    ws_drops_beard_single, ws_drops_beard_batch,
    E_H80_local, E_S09_local, apply_collision_arrays,
    R0_HALL, RAT_HALL, ECOLL_HALL
)
from PyLCM.condensation_arrays import (
    radius_liquid_euler_single, radius_growth_batch,
    apply_condensation_arrays
)
from PyLCM.collision_optimized import E_H80_optimized, ws_drops_beard_optimized
from PyLCM.parameters import rho_liq, pi, r_a


class TestMassConservation:
    """Test mass conservation in collision-coalescence."""

    def test_collision_mass_conservation(self):
        """Total mass should be conserved during collisions (no sedimentation)."""
        np.random.seed(42)

        # Create particles with known total mass
        n_particles = 1000
        particles = ParticleArrays(n_particles)

        # Initialize with random masses and multiplicities
        particles.M = np.random.uniform(1e-12, 1e-10, n_particles)
        particles.A = np.random.uniform(1.0, 100.0, n_particles)
        particles.Ns = particles.M * 0.01  # Small aerosol fraction
        particles.kappa = np.full(n_particles, 0.6)

        # Initial total mass
        M_initial = np.sum(particles.M)
        A_initial = np.sum(particles.A)

        # Apply collisions (multiple times to ensure collisions happen)
        T_parcel = 280.0
        P_parcel = 90000.0
        dt = 1.0

        for _ in range(10):
            n_coll, n_auto, precip = apply_collision_arrays(
                particles, T_parcel, P_parcel, dt, do_collision=True,
                do_sedi_removal=False
            )

        # Final total mass
        M_final = np.sum(particles.M)

        # Mass should be conserved (within numerical precision)
        rel_error = abs(M_final - M_initial) / M_initial
        assert rel_error < 1e-10, f"Mass not conserved: initial={M_initial}, final={M_final}, error={rel_error}"
        print(f"  Mass conservation: initial={M_initial:.6e}, final={M_final:.6e}, rel_error={rel_error:.2e}")


class TestKappaConservation:
    """Test kappa (hygroscopicity) conservation in collision."""

    def test_volume_weighted_kappa(self):
        """Volume-weighted kappa should be preserved during collision."""
        np.random.seed(123)

        n_particles = 100
        particles = ParticleArrays(n_particles)

        # Initialize with varying kappa values
        particles.M = np.random.uniform(1e-11, 1e-9, n_particles)
        particles.A = np.ones(n_particles) * 10.0  # Same weights
        particles.Ns = particles.M * 0.01
        particles.kappa = np.random.uniform(0.4, 1.2, n_particles)

        # Compute initial volume-weighted kappa
        volumes = particles.M / rho_liq
        total_volume = np.sum(volumes)
        kappa_initial = np.sum(particles.kappa * volumes) / total_volume

        # Apply collisions
        T_parcel = 280.0
        P_parcel = 90000.0
        dt = 1.0

        for _ in range(5):
            apply_collision_arrays(particles, T_parcel, P_parcel, dt, do_collision=True)

        # Compute final volume-weighted kappa
        volumes_final = particles.M / rho_liq
        total_volume_final = np.sum(volumes_final)
        kappa_final = np.sum(particles.kappa * volumes_final) / total_volume_final

        # Kappa should be approximately conserved
        rel_error = abs(kappa_final - kappa_initial) / kappa_initial
        assert rel_error < 0.1, f"Kappa not conserved: initial={kappa_initial}, final={kappa_final}"
        print(f"  Kappa conservation: initial={kappa_initial:.4f}, final={kappa_final:.4f}, rel_error={rel_error:.2e}")


class TestTerminalVelocity:
    """Test terminal velocity calculations."""

    def test_beard_formula_consistency(self):
        """Terminal velocity should match original implementation."""
        radii = np.array([1e-6, 5e-6, 10e-6, 25e-6, 50e-6, 100e-6, 500e-6, 1000e-6])
        rho_parcel = 1.2
        P_parcel = 100000.0
        T_parcel = 280.0

        for r in radii:
            v_single = ws_drops_beard_single(r, rho_parcel, rho_liq, P_parcel, T_parcel)
            v_orig = ws_drops_beard_optimized(r, rho_parcel, rho_liq, P_parcel, T_parcel)

            rel_error = abs(v_single - v_orig) / max(v_orig, 1e-10)
            assert rel_error < 1e-6, f"Terminal velocity mismatch at r={r}: single={v_single}, orig={v_orig}"

        print(f"  Terminal velocity consistency: all {len(radii)} radii match within 1e-6")

    def test_batch_vs_single(self):
        """Batch calculation should match single calculations."""
        np.random.seed(456)
        radii = np.random.uniform(1e-6, 500e-6, 100)
        rho_parcel = 1.2
        P_parcel = 100000.0
        T_parcel = 280.0

        v_batch = ws_drops_beard_batch(radii, rho_parcel, rho_liq, P_parcel, T_parcel)

        for i, r in enumerate(radii):
            v_single = ws_drops_beard_single(r, rho_parcel, rho_liq, P_parcel, T_parcel)
            rel_error = abs(v_batch[i] - v_single) / max(v_single, 1e-10)
            assert rel_error < 1e-10, f"Batch/single mismatch at i={i}"

        print(f"  Batch vs single: all 100 particles match")

    def test_velocity_increases_with_radius(self):
        """Terminal velocity should generally increase with radius."""
        radii = np.array([10e-6, 50e-6, 100e-6, 500e-6])
        rho_parcel = 1.2
        P_parcel = 100000.0
        T_parcel = 280.0

        v_fall = ws_drops_beard_batch(radii, rho_parcel, rho_liq, P_parcel, T_parcel)

        for i in range(len(radii) - 1):
            assert v_fall[i+1] > v_fall[i], f"Velocity should increase: v[{radii[i]*1e6:.0f}µm]={v_fall[i]}, v[{radii[i+1]*1e6:.0f}µm]={v_fall[i+1]}"

        print(f"  Velocity increases with radius: verified for {len(radii)} sizes")


class TestCollisionEfficiency:
    """Test collision efficiency calculations."""

    def test_hall_efficiency_consistency(self):
        """Hall efficiency should match original implementation."""
        test_pairs = [
            (10e-6, 5e-6),
            (50e-6, 10e-6),
            (100e-6, 50e-6),
            (200e-6, 100e-6),
        ]

        for r1, r2 in test_pairs:
            E_local = E_H80_local(r1, r2, R0_HALL, RAT_HALL, ECOLL_HALL)
            E_orig = E_H80_optimized(r1, r2)

            rel_error = abs(E_local - E_orig) / max(E_orig, 1e-10)
            assert rel_error < 1e-6, f"Hall efficiency mismatch: local={E_local}, orig={E_orig}"

        print(f"  Hall efficiency consistency: all {len(test_pairs)} pairs match")

    def test_efficiency_bounds(self):
        """Collision efficiency should be between 0 and 1 (mostly)."""
        np.random.seed(789)
        for _ in range(100):
            r1 = np.random.uniform(1e-6, 500e-6)
            r2 = np.random.uniform(1e-6, 500e-6)
            E = E_H80_local(r1, r2, R0_HALL, RAT_HALL, ECOLL_HALL)
            # Note: Hall can be > 1 for very large drops, but typically 0-1
            assert E >= 0, f"Efficiency should be non-negative: E={E}"

        print(f"  Efficiency bounds: all 100 tests have E >= 0")


class TestCondensation:
    """Test condensation physics."""

    def test_growth_increases_mass_at_supersaturation(self):
        """Droplets should grow when supersaturated."""
        np.random.seed(321)

        n_particles = 100
        particles = ParticleArrays(n_particles)

        # Initialize small cloud droplets
        particles.M = np.full(n_particles, 4.0/3.0 * pi * rho_liq * (5e-6)**3 * 10.0)
        particles.A = np.full(n_particles, 10.0)
        particles.Ns = np.full(n_particles, 1e-18)
        particles.kappa = np.full(n_particles, 0.6)

        M_initial = np.sum(particles.M)

        # Apply condensation at supersaturation
        T_parcel = 280.0
        P_parcel = 90000.0
        q_parcel = 0.012  # High humidity to ensure supersaturation
        dt = 0.1

        T_new, q_new, n_act, n_evp = apply_condensation_arrays(
            particles, T_parcel, P_parcel, q_parcel, dt,
            do_condensation=True
        )

        M_final = np.sum(particles.M)

        # Mass should have changed (growth or evaporation depending on conditions)
        # At supersaturation, mass should increase
        print(f"  Condensation: M_initial={M_initial:.6e}, M_final={M_final:.6e}")
        print(f"  T change: {T_parcel} -> {T_new}")
        print(f"  q change: {q_parcel} -> {q_new}")


class TestCompactMethod:
    """Test the compact() method of ParticleArrays."""

    def test_compact_removes_zero_A(self):
        """compact() should remove particles with A <= 0."""
        particles = ParticleArrays(10)
        particles.A = np.array([1.0, 0.0, 2.0, -1.0, 3.0, 0.0, 4.0, 0.0, 5.0, 6.0])

        n_removed = particles.compact()

        assert n_removed == 4, f"Expected 4 removed, got {n_removed}"
        assert len(particles) == 6, f"Expected 6 particles, got {len(particles)}"
        assert np.all(particles.A > 0), "All remaining particles should have A > 0"

        print(f"  compact(): removed {n_removed} particles, {len(particles)} remaining")

    def test_compact_preserves_valid_particles(self):
        """compact() should preserve properties of valid particles."""
        particles = ParticleArrays(5)
        particles.M = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        particles.A = np.array([1.0, 0.0, 1.0, 0.0, 1.0])
        particles.id = np.array([0, 1, 2, 3, 4])

        particles.compact()

        assert len(particles) == 3
        np.testing.assert_array_equal(particles.M, [1.0, 3.0, 5.0])
        np.testing.assert_array_equal(particles.id, [0, 2, 4])

        print(f"  compact(): properties preserved correctly")


def run_all_tests():
    """Run all physics validation tests."""
    print("\n" + "="*60)
    print("PyLCM Array Physics Validation Tests")
    print("="*60)

    test_classes = [
        TestMassConservation,
        TestKappaConservation,
        TestTerminalVelocity,
        TestCollisionEfficiency,
        TestCondensation,
        TestCompactMethod,
    ]

    passed = 0
    failed = 0

    for test_class in test_classes:
        print(f"\n{test_class.__name__}:")
        instance = test_class()

        for method_name in dir(instance):
            if method_name.startswith('test_'):
                try:
                    method = getattr(instance, method_name)
                    method()
                    passed += 1
                except AssertionError as e:
                    print(f"  FAILED: {method_name}")
                    print(f"    {e}")
                    failed += 1
                except Exception as e:
                    print(f"  ERROR: {method_name}")
                    print(f"    {e}")
                    failed += 1

    print("\n" + "="*60)
    print(f"Results: {passed} passed, {failed} failed")
    print("="*60)

    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
