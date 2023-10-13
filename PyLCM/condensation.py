import numpy as np
import sys
from numba import jit
from PyLCM.parameters import *
from PyLCM.micro_particle import *
from scipy.optimize import newton

# Diffusional growth of aerosols, droplets, ice crystals
def drop_condensation(particles_list, T_parcel, q_parcel, P_parcel, nt, dt, air_mass_parcel, S_lst, rho_aero,kohler_activation_radius, con_ts, act_ts, evp_ts, dea_ts,switch_kappa_koehler):
    
    dq_liq = 0
    # Get supersaturation (via saturated water vapour pressure (e_s) and water vapour pressure of the parcel (e_a))
    e_s = esatw( T_parcel )
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    V_parcel   = 100.0 /P_parcel / ( r_a * T_parcel )
    supersat = e_a / e_s - 1.0

    # Thermal conductivity for water (according to Rogers and Yau, Table 7.1)
    thermal_conductivity = 7.94048E-05 * T_parcel + 0.00227011
    # Moldecular diffusivity of water vapor in air (according to Hall and Pruppacher, 1976)
    diff_coeff = 0.211E-4 * (T_parcel / 273.15) ** 1.94 * (101325.0 / P_parcel)
    # Prefactor in diffusional growth equation
    G_pre = 1.0 / (rho_liq * rv * T_parcel / (e_s * diff_coeff) + (l_v / (rv * T_parcel) - 1.0) *
                   rho_liq * l_v / (thermal_conductivity * T_parcel))

    alpha = 0.036
    r0 = diff_coeff / alpha * np.sqrt(2.0 * np.pi / (rv * T_parcel)) / (1.0 + diff_coeff * l_v ** 2 *
                                                                     e_s / (thermal_conductivity * rv ** 2 *
                                                                           T_parcel ** 3))
    # Ventilation effects, TBD
    f_vent = 1.0
    # Radiation effects, TBD
    D_pre = 0.0
    radiation = 0.0
        
    for particle in particles_list:
        dq_liq = dq_liq - particle.M

        # Computation of factors a and b related to Koehler curve
        afactor = 2.0 * sigma_air_liq(T_parcel) / (rho_liq * rv * T_parcel) # Curvature effect
    
        if switch_kappa_koehler:
            bfactor = particle.kappa
        else:
            bfactor = vanthoff_aero * rho_aero * molecular_weight_water / (rho_liq * molecular_weight_aero) # Solute effect
        
        # Initial radius
        r_liq = (particle.M / (particle.A * 4.0 / 3.0 * np.pi * rho_liq)) ** 0.33333333333 # Droplet radius
        r_N = (particle.Ns / (particle.A * 4.0 / 3.0 * np.pi * rho_aero)) ** 0.33333333333 # Aerosol radius
        
        # Old particle liquid mass for growth rate calculation        
        M_old  = particle.M
        if kohler_activation_radius:
            activation_radius = np.sqrt( 3.0 * bfactor * r_N**3 / afactor )
        else:
            activation_radius = activation_radius_ts
        # Diffusional growth
        r_liq_old = r_liq
        r_liq = radius_liquid_euler(r_liq, dt, r0, G_pre, supersat, f_vent, afactor, bfactor, r_N, D_pre, radiation)
            
        particle.M = particle.A * 4.0 / 3.0 * np.pi * rho_liq * r_liq ** 3
        
        if r_liq_old < r_liq:
            con_ts = con_ts +  (particle.M - M_old)
            if (r_liq >= activation_radius) and (r_liq_old < activation_radius):
                # Mass of activated droplets
                act_ts = act_ts + (particle.M - M_old)
            
        else:
            evp_ts = evp_ts  +  (particle.M - M_old)
            if (r_liq < activation_radius) and (r_liq_old >= activation_radius):
                # Mass of deactivated droplets
                dea_ts = dea_ts + (particle.M - M_old)
        
        dq_liq = dq_liq + particle.M

    T_parcel = T_parcel + dq_liq * l_v / cp / air_mass_parcel
    q_parcel = q_parcel - dq_liq / air_mass_parcel
    
    e_s = esatw( T_parcel )
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    S_lst = e_a - e_s
        
    return particles_list, T_parcel, q_parcel, S_lst, con_ts, act_ts, evp_ts, dea_ts 

def esatw(T):
    # Saturation water vapour pressure (Pa) (Flatau et.al, 1992, JAM)
    a = [6.11239921, 0.443987641, 0.142986287e-1,
         0.264847430e-3, 0.302950461e-5, 0.206739458e-7,
         0.640689451e-10, -0.952447341e-13, -0.976195544e-15]

    dT = T - 273.15
    esatw = a[0] + dT * (a[1] + dT * (a[2] + dT * (a[3] + dT * (a[4] + dT * (a[5] + dT * (a[6] + dT * (a[7] + a[8] * dT)))))))
    esatw *= 100.0

    return esatw

def sigma_air_liq(tabs):
    # Surface tension between liquid water and air (in J/m2)
    tabs_c = tabs - 273.15
    
    # Pruppacher and Klett (1997), Eq. 5-12
    sigma_air_liq = 75.93 + 0.115 * tabs_c    + 6.818e-2 * tabs_c**2 + 6.511e-3 * tabs_c**3 + 2.933e-4 * tabs_c**4 + 6.283e-6 * tabs_c**5 + 5.285e-8 * tabs_c**6
    sigma_air_liq = sigma_air_liq * 1.0E-3
    
    return(sigma_air_liq)

def radius_liquid_euler(r_ini, dt_int, r0, G_pre, supersat, ventilation_effect, afactor, bfactor, r_aero, D_pre, radiation):

    r_eul = r_ini
    r_eul_old = r_ini

    dt_eul = dt_int
    t_eul = 0.0

    while t_eul < dt_int - 1.0e-20:

        for m in range(500):

            dr2dt = 2.0 * G_pre * ventilation_effect * (supersat - afactor / r_eul + bfactor * r_aero**3 / r_eul**3 - D_pre * radiation * r_eul) * r_eul / (r_eul + r0)

            d2r2dtdr2 = G_pre * ventilation_effect * (afactor * r_eul**3 - bfactor * r_aero**3 * (3.0 * r_eul + 2.0 * r0) - r_eul**3 * (D_pre * radiation * r_eul * (r_eul + 2.0 * r0) - r0 * supersat)) / (r_eul**4 * (r_eul + r0)**2)

            dt_eul = min(dt_int - t_eul, dt_int)  
            
            # To speed up the Newton-Raphson scheme, the square root is executed at every iteration
            f = r_eul**2 - r_ini**2 - dt_eul * dr2dt
            dfdr2 = 1.0 - dt_eul * d2r2dtdr2

            r_eul = (max(r_eul**2 - f / dfdr2, r_aero**2))**0.5  # Newton-Raphson scheme (2nd order)

            rel_change = abs(r_eul - r_eul_old) / r_eul_old
            r_eul_old = r_eul
            if rel_change < 1.0e-12:
                break

        t_eul += dt_eul
        r_ini = r_eul
        
    return r_eul
