import numpy as np
import sys
from numba import jit
from PyLCM.parameters import *
from PyLCM.micro_particle import *
from scipy.optimize import newton

# Supersaturation relaxation functions following Tzivion, Feingold, and Levin (1989, JAS)
def eta_new_func(eta_lst, forcing, tau, dt):
    """Calculate new absolute supersaturation using exponential relaxation"""
    if tau < 1e-20 or dt < 1e-20:
        return eta_lst
    return tau * forcing + (eta_lst - tau * forcing) * np.exp(-dt / tau)

def eta_int_func(eta_lst, forcing, tau, dt):
    """Calculate time-integrated absolute supersaturation"""
    if tau < 1e-20 or dt < 1e-20:
        return eta_lst * dt
    return tau * (1.0 - np.exp(-dt / tau)) * (eta_lst - tau * forcing) + tau * forcing * dt


# Diffusional growth of aerosols, droplets
# Uses Tzivion, Feingold, and Levin (1989, JAS) method for supersaturation
def drop_condensation(particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass_parcel, S_lst, rho_aero, kohler_activation_radius, con_ts, act_ts, evp_ts, dea_ts, switch_kappa_koehler):

    dq_liq = 0

    # Get supersaturation
    e_s = esatw(T_parcel)
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    qsatw_pre = e_s / (P_parcel - e_s) * r_a / rv  # saturation mixing ratio
    supersat = e_a / e_s - 1.0

    # Parcel volume (m³)
    rho_parcel = P_parcel / (r_a * T_parcel)
    V_parcel = air_mass_parcel / rho_parcel

    # Thermal conductivity for water (according to Rogers and Yau, Table 7.1)
    thermal_conductivity = 7.94048E-05 * T_parcel + 0.00227011
    # Molecular diffusivity of water vapor in air (according to Hall and Pruppacher, 1976)
    diff_coeff = 0.211E-4 * (T_parcel / 273.15) ** 1.94 * (101325.0 / P_parcel)

    # G_pre (= C_liq_pre in Fortran code)
    G_pre = 1.0 / (rho_liq * rv * T_parcel / (e_s * diff_coeff) +
                   (l_v / (rv * T_parcel) - 1.0) * rho_liq * l_v / (thermal_conductivity * T_parcel))

    # Psychrometric correction factor (A_liq_pre)
    # A = 1 + (dqsat/dT) * L / cp
    dqsatw_dT = qsatw_pre * l_v / (rv * T_parcel**2)
    A_liq_pre = 1.0 + dqsatw_dT * l_v / cp

    alpha = 0.036
    r0 = diff_coeff / alpha * np.sqrt(2.0 * np.pi / (rv * T_parcel)) / \
         (1.0 + diff_coeff * l_v ** 2 * e_s / (thermal_conductivity * rv ** 2 * T_parcel ** 3))

    # Ventilation effects
    f_vent = 1.0
    # Radiation effects
    D_pre = 0.0
    radiation = 0.0

    afactor = 2.0 * sigma_air_liq(T_parcel) / (rho_liq * rv * T_parcel)

    # =========================================================================
    # Tzivion et al. (1989) supersaturation adjustment
    # Following SAM/MICRO_LAGRANGE implementation
    # =========================================================================

    # Calculate pre_liq factor
    pre_liq = A_liq_pre * 4.0 * np.pi * G_pre / (qsatw_pre * V_parcel)

    # Calculate dsd_liq (droplet size distribution integral)
    dsd_liq = 0.0
    for particle in particles_list:
        if particle.A > 0:
            r_liq = (particle.M / (particle.A * 4.0 / 3.0 * np.pi * rho_liq)) ** 0.33333333333
            # Sum of A * r²/(r+r0) * f_vent
            dsd_liq += particle.A * (r_liq ** 2 / (r_liq + r0)) * f_vent

    # Relaxation timescale
    if dsd_liq > 0:
        tau = 1.0 / (pre_liq * dsd_liq)
    else:
        tau = 1e10

    # Dynamic supersaturation (absolute, in kg/kg)
    eta_dyn = q_parcel - qsatw_pre
    # Previous supersaturation
    eta_lst = S_lst if abs(S_lst) > 1e-20 else eta_dyn

    # Forcing term (rate of supersaturation change from dynamics)
    forcing = (eta_dyn - eta_lst) / dt if dt > 0 else 0.0

    # Calculate time-mean supersaturation using Tzivion et al. method
    eta_new_tmp = eta_new_func(eta_lst, forcing, tau, dt)
    eta_int_tmp = eta_int_func(eta_lst, forcing, tau, dt)

    # Switch factor for numerical stability
    if dt > 0:
        fsw_term1 = (min(eta_lst, eta_new_tmp) * dt - eta_int_tmp) * 1.0E35
        fsw_term2 = (eta_int_tmp - max(eta_lst, eta_new_tmp) * dt) * 1.0E35
        fsw = max(min(1.0, max(fsw_term1, 0.0)), min(1.0, max(fsw_term2, 0.0)))
        eta_mean = (eta_int_tmp * (1.0 - fsw) + fsw * (eta_lst + eta_new_tmp) * dt * 0.5) / dt
    else:
        eta_mean = eta_dyn

    # Convert to relative supersaturation
    supersat_adjusted = eta_mean / qsatw_pre if qsatw_pre > 0 else supersat

    # Use adjusted supersaturation for droplet growth
    supersat_use = supersat_adjusted

    # =========================================================================
    # Grow each particle using adjusted supersaturation
    # =========================================================================

    for particle in particles_list:
        dq_liq = dq_liq - particle.M

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
        r_liq = radius_liquid_euler(r_liq, dt, r0, G_pre, supersat_use, f_vent, afactor, bfactor, r_N, D_pre, radiation)

        particle.M = particle.A * 4.0 / 3.0 * np.pi * rho_liq * r_liq ** 3

        if r_liq_old < r_liq:
            con_ts = con_ts + (particle.M - M_old)
            if (r_liq >= activation_radius) and (r_liq_old < activation_radius):
                act_ts = act_ts + (particle.M - M_old)
        else:
            evp_ts = evp_ts + (particle.M - M_old)
            if (r_liq < activation_radius) and (r_liq_old >= activation_radius):
                dea_ts = dea_ts + (particle.M - M_old)

        dq_liq = dq_liq + particle.M

    T_parcel = T_parcel + dq_liq * l_v / cp / air_mass_parcel
    q_parcel = q_parcel - dq_liq / air_mass_parcel

    # Update S_lst for next timestep (absolute supersaturation in kg/kg)
    qsatw_new = esatw(T_parcel) / (P_parcel - esatw(T_parcel)) * r_a / rv
    S_lst = q_parcel - qsatw_new

    return particles_list, T_parcel, q_parcel, S_lst, con_ts, act_ts, evp_ts, dea_ts


def esatw(T):
    # Saturation water vapour pressure over liquid water (Pa) (Flatau et.al, 1992, JAM)
    a = [6.11239921, 0.443987641, 0.142986287e-1,
         0.264847430e-3, 0.302950461e-5, 0.206739458e-7,
         0.640689451e-10, -0.952447341e-13, -0.976195544e-15]

    dT = T - 273.15
    esatw = a[0] + dT * (a[1] + dT * (a[2] + dT * (a[3] + dT * (a[4] + dT * (a[5] + dT * (a[6] + dT * (a[7] + a[8] * dT)))))))
    esatw *= 100.0

    return esatw


def esati(T):
    # Saturation water vapour pressure over ice (Pa) (Flatau et.al, 1992, JAM)
    # Same formulation as SAM LCM
    a = [6.11147274, 0.503160820, 0.188439774e-1,
         0.420895665e-3, 0.615021634e-5, 0.602588177e-7,
         0.385852041e-9, 0.146898966e-11, 0.252751365e-14]

    if T > 273.15:
        # Above freezing, use liquid saturation
        return esatw(T)
    elif T > 185.0:
        dT = T - 273.16
        esati_val = a[0] + dT * (a[1] + dT * (a[2] + dT * (a[3] + dT * (a[4] + dT * (a[5] + dT * (a[6] + dT * (a[7] + a[8] * dT)))))))
        return esati_val * 100.0
    else:
        # Below 185K, use additional interpolation
        dT = max(-100.0, T - 273.16)
        esati_val = 0.00763685 + dT * (0.000151069 + dT * 7.48215e-07)
        return esati_val * 100.0


def sigma_air_liq(tabs):
    # Surface tension between liquid water and air (in J/m2)
    tabs_c = tabs - 273.15

    # Pruppacher and Klett (1997), Eq. 5-12
    sigma_air_liq = 75.93 + 0.115 * tabs_c + 6.818e-2 * tabs_c**2 + 6.511e-3 * tabs_c**3 + \
                    2.933e-4 * tabs_c**4 + 6.283e-6 * tabs_c**5 + 5.285e-8 * tabs_c**6
    sigma_air_liq = sigma_air_liq * 1.0E-3

    return sigma_air_liq


def radius_liquid_euler(r_ini, dt_int, r0, G_pre, supersat, ventilation_effect, afactor, bfactor, r_aero, D_pre, radiation):
    r_eul = r_ini
    r_eul_old = r_ini
    dt_eul = dt_int
    t_eul = 0.0

    while t_eul < dt_int - 1.0e-20:
        for m in range(500):
            dr2dt = 2.0 * G_pre * ventilation_effect * \
                    (supersat - afactor / r_eul + bfactor * r_aero**3 / r_eul**3 - D_pre * radiation * r_eul) * \
                    r_eul / (r_eul + r0)

            d2r2dtdr2 = G_pre * ventilation_effect * \
                        (afactor * r_eul**3 - bfactor * r_aero**3 * (3.0 * r_eul + 2.0 * r0) -
                         r_eul**3 * (D_pre * radiation * r_eul * (r_eul + 2.0 * r0) - r0 * supersat)) / \
                        (r_eul**4 * (r_eul + r0)**2)

            dt_eul = min(dt_int - t_eul, dt_int)

            f = r_eul**2 - r_ini**2 - dt_eul * dr2dt
            dfdr2 = 1.0 - dt_eul * d2r2dtdr2

            r_eul = (max(r_eul**2 - f / dfdr2, r_aero**2))**0.5

            rel_change = abs(r_eul - r_eul_old) / r_eul_old
            r_eul_old = r_eul
            if rel_change < 1.0e-12:
                break

        t_eul += dt_eul
        r_ini = r_eul

    return r_eul
