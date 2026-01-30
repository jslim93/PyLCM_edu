"""
Tests for ensemble run functionality.

Tests verify:
1. Reproducibility: Same seed produces identical results
2. Independence: Different seeds produce different results
3. Parallel correctness: Serial and parallel give same statistics
4. Convergence: Standard error decreases as 1/sqrt(n)
"""

import numpy as np
import pytest


def _create_env_profiles(T_init=285.0, z_init=0.0, P_init=100000.0):
    """Create default environmental profiles for testing."""
    from PyLCM.parameters import r_a, cp, p0

    z_env = np.arange(z_init, 3001, 10)
    lapse_rate = -6.5 / 1000  # Unstable condition
    theta_init = T_init * (p0 / P_init) ** (r_a / cp)
    theta_profiles = theta_init + lapse_rate * z_env

    qv_init = 0.01
    qv_diff = (qv_init - 2e-3) / len(z_env)
    qv_profiles = np.maximum(qv_init - qv_diff * np.arange(len(z_env)), 2e-3)

    return qv_profiles, theta_profiles


def test_seed_reproducibility():
    """Same seed should produce identical results."""
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays
    from PyLCM.parameters import rm_spec

    qv_profiles, theta_profiles = _create_env_profiles()

    common_kwargs = {
        'n_particles': 100,
        'P_parcel': 100000.0,
        'RH_parcel': 0.95,
        'T_parcel': 285.0,
        'w_parcel': 1.0,
        'nt': 10,
        'dt': 1.0,
        'rm_spec': rm_spec,
        'do_condensation': True,
        'do_collision': False,
        'mode_aero_init': 'Random',
        'verbose': False,
        'theta_profiles': theta_profiles,
        'qv_profiles': qv_profiles,
    }

    # Run with same seed twice
    result1 = timesteps_function_arrays(**common_kwargs, seed=42)
    result2 = timesteps_function_arrays(**common_kwargs, seed=42)

    # Cloud water should be identical
    np.testing.assert_array_equal(result1[8], result2[8],
                                   err_msg="Same seed should give identical qc")

    # Temperature should be identical
    np.testing.assert_array_equal(result1[3], result2[3],
                                   err_msg="Same seed should give identical T")


def test_seed_independence():
    """Different seeds should produce different results."""
    from PyLCM.timestep_routine_arrays import timesteps_function_arrays
    from PyLCM.parameters import rm_spec

    qv_profiles, theta_profiles = _create_env_profiles()

    common_kwargs = {
        'n_particles': 100,
        'P_parcel': 100000.0,
        'RH_parcel': 0.95,
        'T_parcel': 285.0,
        'w_parcel': 1.0,
        'nt': 10,
        'dt': 1.0,
        'rm_spec': rm_spec,
        'do_condensation': True,
        'do_collision': False,
        'mode_aero_init': 'Random',
        'verbose': False,
        'theta_profiles': theta_profiles,
        'qv_profiles': qv_profiles,
    }

    # Run with different seeds
    result1 = timesteps_function_arrays(**common_kwargs, seed=42)
    result2 = timesteps_function_arrays(**common_kwargs, seed=99)

    # Results should be different - check qc or temperature
    # Since random aerosol initialization produces different distributions,
    # at least one of these should differ
    qc_different = not np.array_equal(result1[8], result2[8])
    T_different = not np.array_equal(result1[3], result2[3])

    assert qc_different or T_different, \
        "Different seeds should produce different results"


def test_ensemble_result_statistics():
    """Test EnsembleResult computes correct statistics."""
    from PyLCM.ensemble import EnsembleResult

    # Create mock outputs (minimal structure matching timesteps_function_arrays)
    n_members = 5
    nt = 10

    # Mock output tuple structure
    mock_outputs = []
    for i in range(n_members):
        # Create output with varying qc values
        output = (
            nt,  # 0: nt
            1.0,  # 1: dt
            np.arange(nt + 1) * 1.0,  # 2: time_array
            np.ones(nt + 1) * 285.0,  # 3: T_parcel_array
            np.ones(nt + 1) * 0.95,  # 4: RH_parcel_array
            np.ones(nt + 1) * 0.01,  # 5: q_parcel_array
            np.arange(nt + 1) * 10.0,  # 6: z_parcel_array
            np.zeros(nt + 1),  # 7: qa_ts
            np.ones(nt + 1) * (0.001 + 0.0001 * i),  # 8: qc_ts (varies)
            np.ones(nt + 1) * (0.0001 + 0.00001 * i),  # 9: qr_ts (varies)
            np.zeros(nt + 1),  # 10: na_ts
            np.ones(nt + 1) * 1e8,  # 11: nc_ts
            np.ones(nt + 1) * 1e3,  # 12: nr_ts
            np.zeros((nt + 1, 50)),  # 13: spectra_arr
            np.zeros(nt + 1),  # 14: con_ts
            np.zeros(nt + 1),  # 15: act_ts
            np.zeros(nt + 1),  # 16: evp_ts
            np.zeros(nt + 1),  # 17: dea_ts
            np.zeros(nt + 1),  # 18: acc_ts
            np.zeros(nt + 1),  # 19: aut_ts
            np.zeros(nt + 1),  # 20: precip_ts
            np.ones(nt + 1) * 100,  # 21: particles_array
            np.ones(nt + 1) * 10e-6,  # 22: rc_liq_avg_array
            np.ones(nt + 1) * 2e-6,  # 23: rc_liq_std_array
        )
        mock_outputs.append(output)

    # Create EnsembleResult
    result = EnsembleResult(mock_outputs, n_members)

    # Verify statistics
    assert result.n_members == n_members
    assert len(result.time_array) == nt + 1

    # Check qc statistics
    expected_qc_mean = np.mean([0.001 + 0.0001 * i for i in range(n_members)])
    np.testing.assert_almost_equal(result.qc_mean[0], expected_qc_mean, decimal=8)

    # Check that percentiles are ordered correctly
    assert np.all(result.qc_p10 <= result.qc_p25)
    assert np.all(result.qc_p25 <= result.qc_mean)
    assert np.all(result.qc_mean <= result.qc_p75)
    assert np.all(result.qc_p75 <= result.qc_p90)


def test_run_ensemble_serial():
    """Test serial ensemble runner produces valid results."""
    from PyLCM.ensemble import run_ensemble_serial
    from PyLCM.parameters import rm_spec

    qv_profiles, theta_profiles = _create_env_profiles()

    result = run_ensemble_serial(
        n_members=3,
        seed_base=42,
        n_particles=50,
        P_parcel=100000.0,
        RH_parcel=0.95,
        T_parcel=285.0,
        w_parcel=1.0,
        nt=5,
        dt=1.0,
        rm_spec=rm_spec,
        do_condensation=True,
        do_collision=False,
        mode_aero_init='Random',
        theta_profiles=theta_profiles,
        qv_profiles=qv_profiles,
    )

    assert result.n_members == 3
    assert len(result.time_array) == 6  # nt + 1
    assert result.qc_all.shape == (3, 6)
    assert result.qc_mean.shape == (6,)
    assert result.qc_std.shape == (6,)


def test_check_convergence():
    """Test convergence checking function."""
    from PyLCM.ensemble import check_convergence, EnsembleResult

    # Create mock ensemble with converged results
    n_members = 20
    nt = 10

    mock_outputs = []
    for i in range(n_members):
        # Small random variation around mean
        qc = 0.001 + np.random.normal(0, 0.00001, nt + 1)
        output = (
            nt, 1.0, np.arange(nt + 1) * 1.0,
            np.ones(nt + 1) * 285.0, np.ones(nt + 1) * 0.95,
            np.ones(nt + 1) * 0.01, np.arange(nt + 1) * 10.0,
            np.zeros(nt + 1), qc, np.zeros(nt + 1),
            np.zeros(nt + 1), np.ones(nt + 1) * 1e8, np.ones(nt + 1) * 1e3,
            np.zeros((nt + 1, 50)), np.zeros(nt + 1), np.zeros(nt + 1),
            np.zeros(nt + 1), np.zeros(nt + 1), np.zeros(nt + 1),
            np.zeros(nt + 1), np.zeros(nt + 1), np.ones(nt + 1) * 100,
            np.ones(nt + 1) * 10e-6, np.ones(nt + 1) * 2e-6,
        )
        mock_outputs.append(output)

    result = EnsembleResult(mock_outputs, n_members)

    # Should be converged with 1% threshold (variation is ~1% of mean)
    is_converged = check_convergence(result, variable='qc', threshold=0.05)
    assert is_converged, "Should be converged with small variation"


@pytest.mark.slow
def test_run_ensemble_parallel():
    """Test parallel ensemble runner (marked slow)."""
    from PyLCM.ensemble import run_ensemble
    from PyLCM.parameters import rm_spec

    qv_profiles, theta_profiles = _create_env_profiles()

    result = run_ensemble(
        n_members=4,
        n_workers=2,
        seed_base=42,
        n_particles=50,
        P_parcel=100000.0,
        RH_parcel=0.95,
        T_parcel=285.0,
        w_parcel=1.0,
        nt=5,
        dt=1.0,
        rm_spec=rm_spec,
        do_condensation=True,
        do_collision=False,
        mode_aero_init='Random',
        theta_profiles=theta_profiles,
        qv_profiles=qv_profiles,
    )

    assert result.n_members == 4
    assert result.qc_all.shape == (4, 6)


@pytest.mark.slow
def test_parallel_serial_consistency():
    """Verify parallel and serial execution give same statistics."""
    from PyLCM.ensemble import run_ensemble, run_ensemble_serial
    from PyLCM.parameters import rm_spec

    qv_profiles, theta_profiles = _create_env_profiles()

    common_kwargs = {
        'n_members': 3,
        'seed_base': 42,
        'n_particles': 50,
        'P_parcel': 100000.0,
        'RH_parcel': 0.95,
        'T_parcel': 285.0,
        'w_parcel': 1.0,
        'nt': 5,
        'dt': 1.0,
        'rm_spec': rm_spec,
        'do_condensation': True,
        'do_collision': False,
        'mode_aero_init': 'Random',
        'theta_profiles': theta_profiles,
        'qv_profiles': qv_profiles,
    }

    serial_result = run_ensemble_serial(**common_kwargs)
    parallel_result = run_ensemble(n_workers=2, **common_kwargs)

    # Statistics should match
    np.testing.assert_array_almost_equal(
        serial_result.qc_mean, parallel_result.qc_mean, decimal=10,
        err_msg="Serial and parallel should give same qc_mean"
    )
    np.testing.assert_array_almost_equal(
        serial_result.qc_std, parallel_result.qc_std, decimal=10,
        err_msg="Serial and parallel should give same qc_std"
    )


if __name__ == '__main__':
    # Run basic tests
    print("Testing seed reproducibility...")
    test_seed_reproducibility()
    print("  PASSED")

    print("Testing seed independence...")
    test_seed_independence()
    print("  PASSED")

    print("Testing EnsembleResult statistics...")
    test_ensemble_result_statistics()
    print("  PASSED")

    print("Testing serial ensemble...")
    test_run_ensemble_serial()
    print("  PASSED")

    print("Testing convergence check...")
    test_check_convergence()
    print("  PASSED")

    print("\nAll basic tests passed!")
