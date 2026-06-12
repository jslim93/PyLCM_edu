import numpy as np
from PyLCM.aero_init import aero_init
from PyLCM.condensation import drop_condensation
from PyLCM.parcel import parcel_rho


def _q_from_rh(RH, T, P):
    from PyLCM.parameters import r_a, rv
    from PyLCM.condensation import esatw
    e_s = esatw(T)
    return RH * e_s / (P - RH * e_s) * r_a / rv


def test_condensation_keeps_T_and_q_finite_over_steps():
    mu_aero = np.log(np.array([0.02e-6, 0.2e-6]))
    sigma_aero = np.log(np.array([1.4, 1.6]))
    N_aero = np.array([100e6, 20e6])
    k_aero = np.array([0.5, 0.5])

    T_parcel = 285.0
    P_parcel = 95000.0
    z_parcel = 0.0
    rho_aero = 1777.0
    q_parcel = _q_from_rh(0.98, T_parcel, P_parcel)

    T_parcel, q_parcel, particles_list = aero_init(
        "Weighting_factor", 50, P_parcel, z_parcel, T_parcel, q_parcel,
        N_aero, mu_aero, sigma_aero, rho_aero, k_aero, True,
    )

    _, _, air_mass_parcel = parcel_rho(P_parcel, T_parcel)

    dt = 0.5
    S_lst = 0.0
    con_ts = act_ts = evp_ts = dea_ts = 0.0

    for nt in range(10):
        # drop_condensation returns:
        # (particles_list, T_parcel, q_parcel, S_lst, con_ts, act_ts, evp_ts, dea_ts)
        # P_parcel is held constant by the caller (not returned).
        (particles_list, T_parcel, q_parcel, S_lst,
         con_ts, act_ts, evp_ts, dea_ts) = drop_condensation(
            particles_list, T_parcel, q_parcel, P_parcel, nt, dt,
            air_mass_parcel, S_lst, rho_aero,
            kohler_activation_radius=True,
            con_ts=con_ts, act_ts=act_ts, evp_ts=evp_ts, dea_ts=dea_ts,
            switch_kappa_koehler=True,
        )

        assert np.isfinite(T_parcel), f"T became non-finite at step {nt}"
        assert np.isfinite(q_parcel), f"q became non-finite at step {nt}"
        assert q_parcel > 0, f"q became non-positive at step {nt}"
