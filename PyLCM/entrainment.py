"""Entrainment parameterization for the PyLCM parcel model.

EXPERIMENTAL: this module is NOT physically validated. A validated entrainment
scheme is planned for a future release. Public physics functions are hard-gated
behind ``experimental=True`` so they cannot be invoked silently in lectures.
"""
import warnings

from PyLCM.parcel import *
from PyLCM.parameters import *
from PyLCM.condensation import *


def _require_experimental(experimental):
    if not experimental:
        raise RuntimeError(
            "entrainment is EXPERIMENTAL and not physically validated. "
            "Pass experimental=True to opt in explicitly. "
            "See the README 'Known limitations' section."
        )
    warnings.warn(
        "Using EXPERIMENTAL, unvalidated entrainment physics.", stacklevel=2
    )


#Get Temperature and vapor mixing ratio for a given altitude by interpolation.
def get_interp1d_var(z_val,z_env, profiles):
    from scipy.interpolate import interp1d
    prof_interp = interp1d(z_env, profiles)
    
    return float(prof_interp(z_val))

#qv_profiles, theta_profiles, z_env = create_env_profiles(initial_theta, initial_qv, z_init, stability_condition)
def basic_entrainment(dt,z_parcel, T_parcel, q_parcel,P_parcel, entrainment_rate,qv_profiles, theta_profiles, experimental=False):
    _require_experimental(experimental)
    qv_env    = get_interp1d_var(z_parcel,z_env,qv_profiles)
    theta_env = get_interp1d_var(z_parcel,z_env,theta_profiles)
    
    T_env     = theta_env * (P_parcel / p0) ** (r_a / cp)
    
    T_parcel  = entrainment_rate * T_env  + (1-entrainment_rate) * T_parcel #+= dT_ent
    q_parcel  = entrainment_rate * qv_env + (1-entrainment_rate) * q_parcel #+= dq_ent


    return T_parcel, q_parcel