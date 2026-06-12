"""Entrainment mixing for PyLCM (warm cloud).

Mixing runs BEFORE condensation each timestep. The Inhomogeneous Mixing Degree
(IHMD, Lim & Hoffmann 2023) controls how entrainment-driven evaporation is split
between homogeneous (all droplets shrink, number conserved) and inhomogeneous
(a subset evaporates, survivors keep size) limits, satisfying exactly:

    N_c / N_{c,0} = (q_c / q_{c,0}) ** IHMD

Closed form per step with entrained fraction frac in [0,1):
    M <- M * (1 - frac)            # total super-droplet liquid mass
    A <- A * (1 - frac) ** IHMD    # multiplicity (droplet number)
"""
import numpy as np

from PyLCM.parameters import p0, r_a, cp, l_v


def redistribute_droplets(particles_list, ihmd, frac):
    """Apply IHMD redistribution to cloud super-droplets in place.

    Returns the total evaporated liquid mass (sum of M lost), which the caller
    returns to the vapor field.
    """
    if frac <= 0.0:
        return 0.0
    keep_mass = 1.0 - frac
    keep_num = keep_mass ** ihmd
    evaporated = 0.0
    for p in particles_list:
        if p.A <= 0 or p.M <= 0:
            continue
        m_old = p.M
        p.M = p.M * keep_mass
        p.A = p.A * keep_num
        evaporated += m_old - p.M
    return evaporated


def _interp(z, z_env, profile):
    return float(np.interp(z, z_env, profile))


class ParameterizedMixing:
    """Parameterized homogeneous/inhomogeneous entrainment mixing.

    lambda_ent : fractional entrainment rate [1/m]; entrained fraction = lambda*w*dt.
    ihmd       : Inhomogeneous Mixing Degree in [0,1] (0 homogeneous, 1 inhomogeneous).
    """

    def __init__(self, lambda_ent, ihmd, qv_profiles, theta_profiles, z_env):
        self.lambda_ent = lambda_ent
        self.ihmd = ihmd
        self.qv_profiles = qv_profiles
        self.theta_profiles = theta_profiles
        self.z_env = z_env

    def apply(self, particles_list, T, q, P, z, dt, w, air_mass):
        frac = self.lambda_ent * w * dt
        if frac <= 0.0:
            return particles_list, T, q
        frac = min(frac, 0.999)
        # 1. Bulk entrainment: relax T, q toward the environment at this height.
        theta_env = _interp(z, self.z_env, self.theta_profiles)
        T_env = theta_env * (P / p0) ** (r_a / cp)
        q_env = _interp(z, self.z_env, self.qv_profiles)
        T = T + frac * (T_env - T)
        q = q + frac * (q_env - q)
        # 2. IHMD redistribution of cloud liquid; evaporated water -> vapor with
        #    latent cooling.
        evaporated = redistribute_droplets(particles_list, self.ihmd, frac)
        dq = evaporated / air_mass
        q = q + dq
        T = T - l_v * dq / cp
        return particles_list, T, q


class LEMMixing:
    """Linear Eddy Model mixing backend (same apply() signature as
    ParameterizedMixing). Not implemented this cycle — a 1D triplet-map domain
    with per-droplet supersaturation perturbation is Phase 3b.
    """

    def apply(self, particles_list, T, q, P, z, dt, w, air_mass):
        raise NotImplementedError(
            "LEMMixing is Phase 3b (1D triplet-map LEM). Use ParameterizedMixing for now."
        )
