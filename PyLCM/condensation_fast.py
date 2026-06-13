"""Fast SoA + numba condensation path.

Drop-in replacement for PyLCM.condensation.drop_condensation that replicates the
per-particle math EXACTLY (same skip rule, same activation/accumulation logic),
but operates on numpy arrays inside an @njit kernel for speed.

The object-based reference in PyLCM/condensation.py is kept intact and remains
the regression anchor; this module must reproduce its golden output within rtol.
"""
import numpy as np
from numba import njit

from PyLCM.parameters import (
    rho_liq, rv, l_v, cp, r_a, vanthoff_aero,
    molecular_weight_water, molecular_weight_aero, activation_radius_ts,
)
from PyLCM.condensation import esatw, sigma_air_liq, radius_liquid_euler

# numpy's pi is the same float64 as the np.pi used in the reference loop
_PI = np.pi


@njit(cache=True)
def _cond_kernel(M, A, Ns, kappa, dt, supersat, G_pre, r0, afactor0,
                 switch_kappa, switch_kelvin, switch_solute,
                 kohler_activation_radius, rho_aero):
    """Replicate drop_condensation's per-particle loop on SoA arrays.

    M is updated in place. Returns (dq_liq, con_ts, act_ts, evp_ts, dea_ts).
    """
    dq_liq = 0.0
    con_ts = 0.0
    act_ts = 0.0
    evp_ts = 0.0
    dea_ts = 0.0

    f_vent = 1.0
    D_pre = 0.0
    radiation = 0.0

    n = M.shape[0]
    for i in range(n):
        dq_liq = dq_liq - M[i]

        # Skip particles with negligible mass (matches reference skip rule)
        if Ns[i] < 1.0e-200 or A[i] <= 0.0 or M[i] <= 0.0:
            dq_liq = dq_liq + M[i]
            continue

        # Koehler factors
        afactor = afactor0  # 2*sigma_air_liq(T)/(rho_liq*rv*T), computed once

        if switch_kappa:
            bfactor = kappa[i]
        else:
            bfactor = vanthoff_aero * rho_aero * molecular_weight_water / (rho_liq * molecular_weight_aero)

        if not switch_kelvin:
            afactor = 0.0
        if not switch_solute:
            bfactor = 0.0

        r_liq = (M[i] / (A[i] * 4.0 / 3.0 * _PI * rho_liq)) ** 0.33333333333
        r_N = (Ns[i] / (A[i] * 4.0 / 3.0 * _PI * rho_aero)) ** 0.33333333333
        M_old = M[i]

        # NOTE: use float exponent (** 3.0) so numba routes to libm pow() and
        # matches CPython's int-exponent `** 3` bit-for-bit. With `** 3` numba
        # emits x*x*x, which differs by 1 ULP and amplifies chaotically over the
        # 300-step Newton-Raphson ascent (breaking the rtol=1e-9 golden match).
        if kohler_activation_radius:
            activation_radius = np.sqrt(3.0 * bfactor * r_N ** 3.0 / afactor)
        else:
            activation_radius = activation_radius_ts

        r_liq_old = r_liq
        r_liq = radius_liquid_euler(r_liq, dt, r0, G_pre, supersat, f_vent,
                                    afactor, bfactor, r_N, D_pre, radiation)

        M[i] = A[i] * 4.0 / 3.0 * _PI * rho_liq * r_liq ** 3.0

        if r_liq_old < r_liq:
            con_ts = con_ts + (M[i] - M_old)
            if (r_liq >= activation_radius) and (r_liq_old < activation_radius):
                act_ts = act_ts + (M[i] - M_old)
        else:
            evp_ts = evp_ts + (M[i] - M_old)
            if (r_liq < activation_radius) and (r_liq_old >= activation_radius):
                dea_ts = dea_ts + (M[i] - M_old)

        dq_liq = dq_liq + M[i]

    return dq_liq, con_ts, act_ts, evp_ts, dea_ts


def drop_condensation_fast(particles_list, T_parcel, q_parcel, P_parcel, nt, dt,
                           air_mass_parcel, S_lst, rho_aero,
                           kohler_activation_radius, con_ts, act_ts, evp_ts,
                           dea_ts, switch_kappa_koehler,
                           switch_kelvin=True, switch_solute=True):
    """Fast drop-in for drop_condensation. Same signature & return tuple."""
    # ---- once-scalars (mirror the reference exactly) ----
    e_s = esatw(T_parcel)
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    supersat = e_a / e_s - 1.0

    thermal_conductivity = 7.94048E-05 * T_parcel + 0.00227011
    diff_coeff = 0.211E-4 * (T_parcel / 273.15) ** 1.94 * (101325.0 / P_parcel)
    G_pre = 1.0 / (rho_liq * rv * T_parcel / (e_s * diff_coeff) + (l_v / (rv * T_parcel) - 1.0) *
                   rho_liq * l_v / (thermal_conductivity * T_parcel))

    alpha = 0.036
    r0 = diff_coeff / alpha * np.sqrt(2.0 * np.pi / (rv * T_parcel)) / (1.0 + diff_coeff * l_v ** 2 *
                                                                        e_s / (thermal_conductivity * rv ** 2 *
                                                                               T_parcel ** 3))
    # afactor depends only on T -> compute once
    afactor0 = 2.0 * sigma_air_liq(T_parcel) / (rho_liq * rv * T_parcel)

    # ---- extract SoA arrays ----
    n = len(particles_list)
    M = np.empty(n, dtype=np.float64)
    A = np.empty(n, dtype=np.float64)
    Ns = np.empty(n, dtype=np.float64)
    kappa = np.empty(n, dtype=np.float64)
    for i, p in enumerate(particles_list):
        M[i] = p.M
        A[i] = p.A
        Ns[i] = p.Ns
        kappa[i] = p.kappa

    dq_liq, d_con, d_act, d_evp, d_dea = _cond_kernel(
        M, A, Ns, kappa, dt, supersat, G_pre, r0, afactor0,
        switch_kappa_koehler, switch_kelvin, switch_solute,
        kohler_activation_radius, rho_aero)

    # ---- write M back to particle objects ----
    for i, p in enumerate(particles_list):
        p.M = M[i]

    con_ts = con_ts + d_con
    act_ts = act_ts + d_act
    evp_ts = evp_ts + d_evp
    dea_ts = dea_ts + d_dea

    T_parcel = T_parcel + dq_liq * l_v / cp / air_mass_parcel
    q_parcel = q_parcel - dq_liq / air_mass_parcel

    e_s = esatw(T_parcel)
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    S_lst = e_a - e_s

    return particles_list, T_parcel, q_parcel, S_lst, con_ts, act_ts, evp_ts, dea_ts


def condense_soa(M, A, Ns, kappa, T_parcel, q_parcel, P_parcel, dt, air_mass_parcel,
                 rho_aero, kohler_activation_radius=False, switch_kappa_koehler=False,
                 switch_kelvin=True, switch_solute=True):
    """Persistent struct-of-arrays condensation step.

    Operates IN-PLACE on the persistent `M` array (no per-step object<->array
    conversion) and reuses the SAME validated `_cond_kernel` as
    `drop_condensation_fast`, so results are bit-identical to the object path —
    this is purely the data-layout optimization, not a physics change.

    Returns updated (T_parcel, q_parcel); `M` is mutated in place.
    """
    e_s = esatw(T_parcel)
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    supersat = e_a / e_s - 1.0

    thermal_conductivity = 7.94048E-05 * T_parcel + 0.00227011
    diff_coeff = 0.211E-4 * (T_parcel / 273.15) ** 1.94 * (101325.0 / P_parcel)
    G_pre = 1.0 / (rho_liq * rv * T_parcel / (e_s * diff_coeff) + (l_v / (rv * T_parcel) - 1.0) *
                   rho_liq * l_v / (thermal_conductivity * T_parcel))
    alpha = 0.036
    r0 = diff_coeff / alpha * np.sqrt(2.0 * np.pi / (rv * T_parcel)) / (1.0 + diff_coeff * l_v ** 2 *
                                                                        e_s / (thermal_conductivity * rv ** 2 *
                                                                               T_parcel ** 3))
    afactor0 = 2.0 * sigma_air_liq(T_parcel) / (rho_liq * rv * T_parcel)

    dq_liq, _, _, _, _ = _cond_kernel(
        M, A, Ns, kappa, dt, supersat, G_pre, r0, afactor0,
        switch_kappa_koehler, switch_kelvin, switch_solute,
        kohler_activation_radius, rho_aero)

    T_parcel = T_parcel + dq_liq * l_v / cp / air_mass_parcel
    q_parcel = q_parcel - dq_liq / air_mass_parcel
    return T_parcel, q_parcel
