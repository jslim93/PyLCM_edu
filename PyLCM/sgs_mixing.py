"""
Linear Eddy Model (LEM) for Subgrid-Scale Mixing
Based on Krueger (1993, JAS) and SAM LCM implementation

For box/parcel model: treats the parcel as a single grid box with
super-droplets arranged in a 1D array for SGS mixing.
"""

import numpy as np
from PyLCM.parameters import *


def sgs_mixing_lem(particles_list, T_parcel, q_parcel, P_parcel, dt,
                   diss_rate=1e-4, L_domain=100.0, T_parcel_old=None, q_parcel_old=None):
    """
    Linear Eddy Model SGS mixing for parcel model.

    Each particle carries its own T and q perturbation that evolves
    through molecular diffusion and turbulent rearrangement.

    Parameters:
    -----------
    particles_list : list
        List of particle objects
    T_parcel : float
        Parcel mean temperature (K) - CURRENT value
    q_parcel : float
        Parcel mean water vapor mixing ratio (kg/kg) - CURRENT value
    P_parcel : float
        Parcel pressure (Pa)
    dt : float
        Time step (s)
    diss_rate : float
        Turbulent dissipation rate (m²/s³), default 1e-4
    L_domain : float
        Domain size for LEM (m), default 100.0
    T_parcel_old : float
        Previous parcel temperature (for tracking changes)
    q_parcel_old : float
        Previous parcel q (for tracking changes)

    Returns:
    --------
    particles_list : list
        Updated particles with new T_lem and q_lem values
    """

    n_ptcl = len(particles_list)
    if n_ptcl < 10:
        # Not enough particles for LEM
        return particles_list

    # Initialize LEM arrays if not present
    for i, p in enumerate(particles_list):
        if not hasattr(p, 'T_lem'):
            p.T_lem = T_parcel
        if not hasattr(p, 'q_lem'):
            p.q_lem = q_parcel
        if not hasattr(p, 'lem_id'):
            p.lem_id = i

    # Update particle T_lem and q_lem to follow parcel's large-scale changes
    # (adiabatic cooling, moisture changes from dynamics)
    if T_parcel_old is not None:
        dT_parcel = T_parcel - T_parcel_old
        for p in particles_list:
            p.T_lem += dT_parcel
    if q_parcel_old is not None:
        dq_parcel = q_parcel - q_parcel_old
        for p in particles_list:
            p.q_lem += dq_parcel

    # Sort particles by LEM id to maintain 1D ordering
    particles_list.sort(key=lambda p: p.lem_id)

    # LEM parameters
    dz_sgs = L_domain / n_ptcl  # SGS segment size
    eta = 6.0 * dz_sgs  # Kolmogorov scale (approximate)

    # Turbulence parameters
    nu = 1.5e-5  # kinematic viscosity of air (m²/s)
    Pr = 1.0     # Prandtl number

    # Turbulent length scale (integral scale)
    L_turb = min((diss_rate / nu**3)**(-0.25) * nu, L_domain)  # Kolmogorov to integral
    L_turb = max(min(L_turb, L_domain), eta)

    # Diffusivities
    D_turb = Pr * (diss_rate * L_turb**4)**(1./3.) / L_turb  # turbulent diffusivity estimate
    D_turb = max(D_turb, 1e-6)
    D_mol = max(D_turb * (eta / L_turb)**1.333, nu)  # molecular diffusivity

    # Rearrangement frequency (Eq. 2.4, Krueger 1993)
    if L_turb > eta:
        lambda_r = 54.0 / 5.0 * D_turb / L_turb**3 * (L_turb / eta)**1.667
    else:
        lambda_r = 0.0

    # Time stepping
    dt_diff = 0.2 * dz_sgs**2 / D_mol
    n_tsteps = max(int(np.ceil(dt / dt_diff)), 1)
    dt_diff = dt / n_tsteps

    dt_turb = 1.0 / (lambda_r * L_domain + 1e-20)
    n_rsteps = max(int(np.ceil(dt_diff / dt_turb)), 1)

    # Extract T and q arrays
    T_old = np.array([p.T_lem for p in particles_list])
    q_old = np.array([p.q_lem for p in particles_list])

    # Time integration
    for _ in range(n_tsteps):
        # Molecular diffusion (periodic boundary)
        T_new = _diffusion_step(T_old, D_mol, dz_sgs, dt_diff)
        q_new = _diffusion_step(q_old, D_mol, dz_sgs, dt_diff)

        # Turbulent rearrangement
        for _ in range(n_rsteps):
            p_rearrange = dt_diff * lambda_r * L_domain
            if p_rearrange > np.random.random():
                T_new, q_new = _triplet_map_rearrangement(
                    T_new, q_new, n_ptcl, dz_sgs, eta, L_turb
                )

        T_old = T_new.copy()
        q_old = q_new.copy()

    # Nudging to parcel mean (relaxation timescale)
    tau_nudge = 900.0  # 15 minutes
    nudge_factor = dt / tau_nudge

    T_mean = np.mean(T_new)
    q_mean = np.mean(q_new)

    T_new = T_new - nudge_factor * (T_mean - T_parcel)
    q_new = q_new - nudge_factor * (q_mean - q_parcel)

    # Update particles
    for i, p in enumerate(particles_list):
        p.T_lem = T_new[i]
        p.q_lem = q_new[i]
        p.lem_id = i  # Update LEM position

    return particles_list


def _diffusion_step(phi, D, dz, dt):
    """
    One step of diffusion with periodic boundaries.

    d(phi)/dt = D * d²(phi)/dz²
    """
    n = len(phi)
    phi_new = np.zeros(n)

    # Central difference with periodic BC
    for i in range(n):
        im1 = (i - 1) % n
        ip1 = (i + 1) % n
        d2phi = phi[ip1] - 2.0 * phi[i] + phi[im1]
        phi_new[i] = phi[i] + dt * D * d2phi / dz**2

    return phi_new


def _triplet_map_rearrangement(T, q, n_ptcl, dz_sgs, eta, L_turb):
    """
    Triplet map rearrangement (Krueger 1993, JAS).

    Selects a random segment and applies the triplet map:
    - First third: every 3rd element starting from 1
    - Second third: reversed middle, every 3rd element
    - Third third: every 3rd element starting from end
    """
    if n_ptcl <= 6:
        return T, q

    if L_turb < eta:
        return T, q

    # Random starting position
    n_start = np.random.randint(0, n_ptcl)

    # Random segment length (Eq. 2.3, Krueger 1993)
    if L_turb / eta < 1.01:
        length = eta
    else:
        # Sample from power law distribution
        length = -1
        while length < eta or length > min(L_turb, n_ptcl * dz_sgs):
            u = np.random.random()
            length = (eta**(-5./3.) - u * (eta**(-5./3.) - L_turb**(-5./3.)))**(-0.6)

    n_length = int(length / dz_sgs)
    n_length = min(max((n_length // 3) * 3, 6), (n_ptcl // 3) * 3)  # multiple of 3

    if n_length < 6:
        return T, q

    # Extract segment (with periodic wrapping)
    indices = [(n_start + i) % n_ptcl for i in range(n_length)]
    T_segment = np.array([T[i] for i in indices])
    q_segment = np.array([q[i] for i in indices])

    # Apply triplet map
    T_mapped = np.zeros(n_length)
    q_mapped = np.zeros(n_length)

    for n in range(n_length):
        if n < n_length // 3:
            m = 3 * n
        elif n < 2 * n_length // 3:
            m = 2 * (n_length - 1) - 3 * n + 1
        else:
            m = 3 * n - 2 * (n_length - 1)

        m = max(0, min(m, n_length - 1))
        T_mapped[n] = T_segment[m]
        q_mapped[n] = q_segment[m]

    # Put back
    T_new = T.copy()
    q_new = q.copy()
    for i, idx in enumerate(indices):
        T_new[idx] = T_mapped[i]
        q_new[idx] = q_mapped[i]

    return T_new, q_new


def get_particle_supersat(particle, T_parcel, q_parcel, P_parcel):
    """
    Get particle-specific supersaturation using LEM T and q.

    If particle has LEM values, use those; otherwise use parcel mean.
    """
    if hasattr(particle, 'T_lem') and hasattr(particle, 'q_lem'):
        T = particle.T_lem
        q = particle.q_lem
    else:
        T = T_parcel
        q = q_parcel

    # Saturation vapor pressure
    from PyLCM.condensation import esatw
    e_s = esatw(T)
    e_a = q * P_parcel / (q + r_a / rv)

    supersat = e_a / e_s - 1.0

    return supersat, T, q
