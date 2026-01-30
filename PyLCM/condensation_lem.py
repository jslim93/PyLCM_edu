"""
Condensation with LEM (Linear Eddy Model) SGS mixing.

Each particle experiences its own supersaturation based on
local T and q from the LEM 1D array.
"""

import numpy as np
from PyLCM.parameters import *
from PyLCM.condensation import esatw, sigma_air_liq, radius_liquid_euler
from PyLCM.sgs_mixing import sgs_mixing_lem, get_particle_supersat


def drop_condensation_lem(particles_list, T_parcel, q_parcel, P_parcel, nt, dt,
                          air_mass_parcel, S_lst, rho_aero, kohler_activation_radius,
                          con_ts, act_ts, evp_ts, dea_ts, switch_kappa_koehler,
                          diss_rate=1e-4, L_domain=100.0):
    """
    Diffusional growth with LEM SGS mixing.

    LEM provides supersaturation FLUCTUATIONS (S') around parcel mean.
    Total supersaturation = S_parcel + S'_LEM

    Parameters:
    -----------
    diss_rate : float
        Turbulent dissipation rate (m²/s³)
    L_domain : float
        LEM domain size (m)
    """

    n_ptcl = len(particles_list)

    # Store old parcel values for tracking large-scale changes
    # Note: The caller should handle adiabatic changes BEFORE calling this function
    # LEM just maintains the fluctuations around the current parcel mean

    # Apply LEM SGS mixing (updates T_lem, q_lem as perturbations)
    # Pass T_parcel and q_parcel as both current and "old" since adiabatic
    # changes should already be applied to particles in the main loop
    particles_list = sgs_mixing_lem(
        particles_list, T_parcel, q_parcel, P_parcel, dt,
        diss_rate=diss_rate, L_domain=L_domain,
        T_parcel_old=None, q_parcel_old=None  # Let main loop handle large-scale
    )

    dq_liq = 0

    # Parcel mean supersaturation
    e_s = esatw(T_parcel)
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    supersat_mean = e_a / e_s - 1.0
    qsatw_pre = e_s / (P_parcel - e_s) * r_a / rv

    # Parcel volume
    rho_parcel = P_parcel / (r_a * T_parcel)
    V_parcel = air_mass_parcel / rho_parcel

    # Thermal conductivity and diffusivity (parcel mean)
    thermal_conductivity = 7.94048E-05 * T_parcel + 0.00227011
    diff_coeff = 0.211E-4 * (T_parcel / 273.15) ** 1.94 * (101325.0 / P_parcel)

    # G_pre
    G_pre = 1.0 / (rho_liq * rv * T_parcel / (e_s * diff_coeff) +
                   (l_v / (rv * T_parcel) - 1.0) * rho_liq * l_v / (thermal_conductivity * T_parcel))

    alpha = 0.036
    r0 = diff_coeff / alpha * np.sqrt(2.0 * np.pi / (rv * T_parcel)) / \
         (1.0 + diff_coeff * l_v ** 2 * e_s / (thermal_conductivity * rv ** 2 * T_parcel ** 3))

    f_vent = 1.0
    D_pre = 0.0
    radiation = 0.0

    afactor = 2.0 * sigma_air_liq(T_parcel) / (rho_liq * rv * T_parcel)

    # Grow each particle using parcel mean + LEM fluctuation
    for particle in particles_list:
        dq_liq = dq_liq - particle.M

        # LEM supersaturation fluctuation
        if hasattr(particle, 'T_lem') and hasattr(particle, 'q_lem'):
            # Compute fluctuation from parcel mean
            T_prime = particle.T_lem - T_parcel
            q_prime = particle.q_lem - q_parcel

            # Linearized supersaturation perturbation: S' ≈ (q'/q_sat) - (L/RvT²) * T'
            dqsat_dT = qsatw_pre * l_v / (rv * T_parcel**2)
            S_prime = q_prime / qsatw_pre - dqsat_dT / qsatw_pre * T_prime
        else:
            S_prime = 0.0

        # Total supersaturation for this particle
        supersat_local = supersat_mean + S_prime

        if switch_kappa_koehler:
            bfactor = particle.kappa
        else:
            bfactor = vanthoff_aero * rho_aero * molecular_weight_water / (rho_liq * molecular_weight_aero)

        r_liq = (particle.M / (particle.A * 4.0 / 3.0 * np.pi * rho_liq)) ** 0.33333333333
        r_N = (particle.Ns / (particle.A * 4.0 / 3.0 * np.pi * rho_aero)) ** 0.33333333333

        M_old = particle.M
        if kohler_activation_radius:
            activation_radius = np.sqrt(3.0 * bfactor * r_N**3 / afactor)
        else:
            activation_radius = activation_radius_ts

        r_liq_old = r_liq
        r_liq = radius_liquid_euler(r_liq, dt, r0, G_pre, supersat_local,
                                    f_vent, afactor, bfactor, r_N, D_pre, radiation)

        particle.M = particle.A * 4.0 / 3.0 * np.pi * rho_liq * r_liq ** 3

        # Statistics
        if r_liq_old < r_liq:
            con_ts = con_ts + (particle.M - M_old)
            if (r_liq >= activation_radius) and (r_liq_old < activation_radius):
                act_ts = act_ts + (particle.M - M_old)
        else:
            evp_ts = evp_ts + (particle.M - M_old)
            if (r_liq < activation_radius) and (r_liq_old >= activation_radius):
                dea_ts = dea_ts + (particle.M - M_old)

        dq_liq = dq_liq + particle.M

    # Update parcel mean T and q (same as original)
    T_parcel = T_parcel + dq_liq * l_v / cp / air_mass_parcel
    q_parcel = q_parcel - dq_liq / air_mass_parcel

    # Update particle LEM values to follow parcel mean shift
    dT_parcel = dq_liq * l_v / cp / air_mass_parcel
    dq_parcel_shift = -dq_liq / air_mass_parcel
    for p in particles_list:
        if hasattr(p, 'T_lem'):
            p.T_lem += dT_parcel
        if hasattr(p, 'q_lem'):
            p.q_lem += dq_parcel_shift

    # Update S_lst
    qsatw_new = esatw(T_parcel) / (P_parcel - esatw(T_parcel)) * r_a / rv
    S_lst = q_parcel - qsatw_new

    return particles_list, T_parcel, q_parcel, S_lst, con_ts, act_ts, evp_ts, dea_ts
