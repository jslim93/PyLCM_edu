from PyLCM.parameters import *
from PyLCM.micro_init import *
from PyLCM.condensation import *
import numpy as np

def parcel_rho(P_parcel, T_parcel):
    
    p_env = P_parcel
    T_env = T_parcel
    theta_env = T_parcel * ( p0 / p_env )**( r_a / cp )
    e_s = esatw( T_parcel )
    
    rho_parcel = p_env / ( r_a * T_parcel ) #  air density
    V_parcel   = 100.0 / rho_parcel # (assumed) volume of parcel for a 100 kg air parcel
    air_mass_parcel = V_parcel * rho_parcel
    
    return(rho_parcel, V_parcel, air_mass_parcel) # (assumed) air mass of parcel


def ascend_parcel(z_parcel, T_parcel,P_parcel,w_parcel,dt, time, max_z, time_half_wave_parcel=1200.0, ascending_mode='linear', t_start_oscillation=800):
    if ascending_mode=='linear':
        if z_parcel < max_z: 
            dz = w_parcel * dt
            z_parcel   = z_parcel + dz
            T_parcel   = T_parcel - dz * g / cp
        
    elif ascending_mode=='sine': 
        w_oscillate = w_parcel * np.pi / 2.0 * np.sin(np.pi * time / time_half_wave_parcel)
        dz = w_oscillate  * dt
        z_parcel = z_parcel + dz
        if w_parcel > 0:
            T_parcel = T_parcel - dz * g / cp
        else:
            T_parcel = T_parcel + dz * g / cp
            
    elif ascending_mode=='in_cloud_oscillation':
        phase = np.arccos(2/np.pi)
        if time < t_start_oscillation:
            dz = w_parcel * dt
            z_parcel   = z_parcel + dz
            T_parcel   = T_parcel - dz * g / cp
        else:
            w_oscillate = w_parcel * np.pi / 2.0 * np.cos(np.pi * (time-t_start_oscillation) / time_half_wave_parcel + phase)
            dz = w_oscillate  * dt
            z_parcel = z_parcel + dz
            if w_parcel > 0:
                T_parcel = T_parcel - dz * g / cp
            else:
                T_parcel = T_parcel + dz * g / cp
    
    
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)
    
    return z_parcel, T_parcel, rho_parcel, V_parcel, air_mass_parcel 