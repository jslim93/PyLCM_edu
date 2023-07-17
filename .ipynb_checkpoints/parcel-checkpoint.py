from parameters import *
from micro import *
import numpy as np

def parcel_rho(P_parcel, T_parcel):
    
    p_env = P_parcel
    T_env = T_parcel
    theta_env = T_parcel * ( p0 / p_env )**( r_a / cp )
    e_s       = 611.2 * np.exp( 17.62 * ( T_env - 273.15 ) / ( T_env - 29.65 ) )
    rho_parcel = p_env / ( r_a * T_parcel )
    V_parcel   = 100.0 / rho_parcel
    air_mass_parcel = V_parcel * rho_parcel
    
    return(rho_parcel, V_parcel, air_mass_parcel) 


def ascend_parcel(z_parcel, T_parcel,P_parcel, dz):
    z_parcel   = z_parcel + dz
    T_parcel   = T_parcel - dz * g / cp
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)
    
    return z_parcel, T_parcel, rho_parcel, V_parcel, air_mass_parcel 