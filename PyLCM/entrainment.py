from PyLCM.parcel import *
from PyLCM.parameters import *
from PyLCM.condensation import *

#Get Temperature and vapor mixing ratio for a given altitude by interpolation.
def get_interp1d_var(z_val,z_env, profiles):
    from scipy.interpolate import interp1d
    prof_interp = interp1d(z_env, profiles)
    
    return float(prof_interp(z_val))

#qv_profiles, theta_profiles, z_env = create_env_profiles(initial_theta, initial_qv, z_init, stability_condition)
def basic_entrainment(dt,z_parcel, T_parcel, q_parcel,P_parcel, entrainment_rate,qv_profiles, theta_profiles):
    qv_env    = get_interp1d_var(z_parcel,z_env,qv_profiles)
    theta_env = get_interp1d_var(z_parcel,z_env,theta_profiles)
    
    T_env     = theta_env * (P_parcel / p0) ** (r_a / cp)
    
    T_parcel  = entrainment_rate * T_env  + (1-entrainment_rate) * T_parcel #+= dT_ent
    q_parcel  = entrainment_rate * qv_env + (1-entrainment_rate) * q_parcel #+= dq_ent


    return T_parcel, q_parcel