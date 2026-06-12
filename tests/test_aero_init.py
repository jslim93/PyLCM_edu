import numpy as np
from PyLCM.aero_init import aero_init


def _q_from_rh(RH, T, P):
    """Convert relative humidity to specific humidity (mirrors model_init)."""
    from PyLCM.parameters import r_a, rv
    from PyLCM.condensation import esatw
    e_s = esatw(T)
    return RH * e_s / (P - RH * e_s) * r_a / rv


def test_aero_init_log_converted_inputs_do_not_overflow():
    # Convention: mu = log(microns * 1e-6), sigma = log(geometric_std)
    mu_aero = np.log(np.array([0.02e-6, 0.2e-6]))
    sigma_aero = np.log(np.array([1.4, 1.6]))
    N_aero = np.array([100e6, 20e6])
    k_aero = np.array([0.5, 0.5])

    T_parcel = 285.0
    P_parcel = 95000.0
    z_parcel = 0.0
    rho_aero = 1777.0
    q_parcel = _q_from_rh(0.98, T_parcel, P_parcel)

    T_out, q_out, particles_list = aero_init(
        "Weighting_factor", 100, P_parcel, z_parcel, T_parcel, q_parcel,
        N_aero, mu_aero, sigma_aero, rho_aero, k_aero, True,
    )

    assert particles_list is not None
    assert len(particles_list) > 0
    assert np.isfinite(T_out)
    assert np.isfinite(q_out)
    # No particle should carry a non-finite mass after initialization
    assert all(np.isfinite(p.M) and p.M >= 0 for p in particles_list)
