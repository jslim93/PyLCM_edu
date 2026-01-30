"""
Ice Nucleation Module for PyLCM

Implements ice nucleation mechanisms following SAM LCM:
1. Homogeneous freezing (Kuhn 2011)
2. ABIFM immersion freezing (Knopf & Alpert 2013)

Reference:
- SAM LCM: ~/SAM6.10.10.LCM_JS/SRC/MICRO_LAGRANGE/micro_cond.f90
"""

import numpy as np
from PyLCM.parameters import *
from PyLCM.condensation import esatw, esati


def homogeneous_freezing(particles_list, T_parcel, dt, air_mass_parcel):
    """
    Homogeneous freezing of supercooled droplets.

    Based on Kuhn (2011) parameterization as implemented in SAM LCM.
    Freezing occurs stochastically when T < 239.1 K (-34°C).

    Parameters:
    -----------
    particles_list : list
        List of particle objects
    T_parcel : float
        Parcel temperature (K)
    dt : float
        Time step (s)
    air_mass_parcel : float
        Air mass in parcel (kg)

    Returns:
    --------
    particles_list : list
        Updated particle list
    T_parcel : float
        Updated temperature (K)
    n_frozen : int
        Number of particles frozen this timestep
    """
    n_frozen = 0
    dq_fus = 0.0  # Latent heat release from fusion

    # Only proceed if T below homogeneous freezing threshold
    if T_parcel >= T_hom_threshold:
        return particles_list, T_parcel, n_frozen

    for particle in particles_list:
        # Skip if already ice or not a liquid droplet
        if particle.micro_type <= 0:
            continue

        # Get droplet radius
        r_liq = (particle.M / (particle.A * 4.0 / 3.0 * pi * rho_liq)) ** (1.0/3.0)

        # Skip very small droplets (haze)
        if r_liq < 1.0e-6:
            continue

        # Volume-specific homogeneous nucleation rate (Kuhn 2011)
        # Js = nv * k_B * T / h * exp((-av + bv*T)/(k_B*T)) * V_drop * surface_correction
        # where V_drop = 4/3 * pi * r^3 = 4.1888 * r^3

        term1 = nv_hom * k_boltzmann * T_parcel / h_planck
        term2 = np.exp((-av_hom + bv_hom * T_parcel) / (k_boltzmann * T_parcel))
        V_drop = 4.1888 * r_liq**3
        surface_correction = 1.0 + 3.0 * (srfthick / r_liq) * np.exp(-dw_hom / (k_boltzmann * T_parcel))

        Js = term1 * term2 * V_drop * surface_correction

        # Nucleation probability
        p_nucleate = 1.0 - np.exp(-dt * Js)

        # Stochastic freezing
        if p_nucleate > np.random.random():
            # Convert to ice
            particle.micro_type = -1  # Ice without IN (homogeneous)

            # Set ice crystal axes (initially spherical)
            # Volume conservation: r_ice = r_liq * (rho_liq/rho_ice)^(1/3)
            r_ice = r_liq * (rho_liq / rho_ice) ** (1.0/3.0)
            particle.aaxis = r_ice
            particle.caxis = r_ice
            particle.dns = rho_ice

            # Track latent heat release
            dq_fus += particle.M

            n_frozen += 1

    # Apply latent heat of fusion
    if dq_fus > 0:
        T_parcel = T_parcel + dq_fus * l_f / cp / air_mass_parcel

    return particles_list, T_parcel, n_frozen


def abifm_immersion_freezing(particles_list, T_parcel, dt, air_mass_parcel):
    """
    ABIFM (water-activity based) immersion freezing.

    Based on Knopf & Alpert (2013) as implemented in SAM LCM.
    Freezing rate depends on INP surface area and water activity difference.

    Parameters:
    -----------
    particles_list : list
        List of particle objects (must have IN_sfc attribute for INP surface area)
    T_parcel : float
        Parcel temperature (K)
    dt : float
        Time step (s)
    air_mass_parcel : float
        Air mass in parcel (kg)

    Returns:
    --------
    particles_list : list
        Updated particle list
    T_parcel : float
        Updated temperature (K)
    n_frozen : int
        Number of particles frozen this timestep
    """
    n_frozen = 0
    dq_fus = 0.0

    # Water activity difference: delta_aw = 1 - e_sat_ice / e_sat_water
    e_sat_w = esatw(T_parcel)
    e_sat_i = esati(T_parcel)
    delta_aw = 1.0 - (e_sat_i / e_sat_w)

    # Only proceed if supercooled
    if delta_aw <= 0:
        return particles_list, T_parcel, n_frozen

    for particle in particles_list:
        # Only process liquid particles with IN (micro_type == 2)
        if particle.micro_type != 2:
            continue

        # Get INP surface area per droplet
        ice_sfc = particle.IN_sfc / particle.A

        # Skip if no INP surface area
        if ice_sfc <= 0:
            continue

        # ABIFM nucleation rate (m^-2 s^-1)
        # Js = 10^(c_ABIFM + m_ABIFM * delta_aw) * 1e4 (convert from cm^-2 to m^-2)
        Js = 10.0**(c_ABIFM + m_ABIFM * delta_aw) * 1.0e4

        # Nucleation probability
        p_nucleate = 1.0 - np.exp(-Js * dt * ice_sfc)

        # Stochastic freezing
        if p_nucleate > np.random.random():
            # Get droplet radius before freezing
            r_liq = (particle.M / (particle.A * 4.0 / 3.0 * pi * rho_liq)) ** (1.0/3.0)

            # Convert to ice
            particle.micro_type = -2  # Ice with IN (heterogeneous)

            # Set ice crystal axes
            r_ice = r_liq * (rho_liq / rho_ice) ** (1.0/3.0)
            particle.aaxis = r_ice
            particle.caxis = r_ice
            particle.dns = rho_ice

            # Track latent heat release
            dq_fus += particle.M

            n_frozen += 1

    # Apply latent heat of fusion
    if dq_fus > 0:
        T_parcel = T_parcel + dq_fus * l_f / cp / air_mass_parcel

    return particles_list, T_parcel, n_frozen


def ice_nucleation(particles_list, T_parcel, dt, air_mass_parcel,
                   switch_homogeneous=True, switch_immersion=True):
    """
    Main ice nucleation routine.

    Applies both homogeneous and immersion freezing mechanisms.

    Parameters:
    -----------
    particles_list : list
        List of particle objects
    T_parcel : float
        Parcel temperature (K)
    dt : float
        Time step (s)
    air_mass_parcel : float
        Air mass in parcel (kg)
    switch_homogeneous : bool
        Enable homogeneous freezing (default: True)
    switch_immersion : bool
        Enable ABIFM immersion freezing (default: True)

    Returns:
    --------
    particles_list : list
        Updated particle list
    T_parcel : float
        Updated temperature (K)
    n_hom : int
        Number of homogeneously frozen particles
    n_imm : int
        Number of immersion frozen particles
    """
    n_hom = 0
    n_imm = 0

    # Immersion freezing first (typically occurs at warmer temperatures)
    if switch_immersion:
        particles_list, T_parcel, n_imm = abifm_immersion_freezing(
            particles_list, T_parcel, dt, air_mass_parcel
        )

    # Homogeneous freezing (only at very cold temperatures < -34°C)
    if switch_homogeneous:
        particles_list, T_parcel, n_hom = homogeneous_freezing(
            particles_list, T_parcel, dt, air_mass_parcel
        )

    return particles_list, T_parcel, n_hom, n_imm


def initialize_INP(particles_list, frac_IN=0.1, r_IN_mean=0.04e-6, sigma_IN=1.4):
    """
    Initialize INP (Ice Nucleating Particle) properties for a fraction of particles.

    For ABIFM, assigns INP surface area based on lognormal size distribution.

    Parameters:
    -----------
    particles_list : list
        List of particle objects
    frac_IN : float
        Fraction of particles that contain INP (default: 0.1)
    r_IN_mean : float
        Geometric mean radius of INP (m) (default: 0.04 um)
    sigma_IN : float
        Geometric standard deviation of INP size distribution (default: 1.4)

    Returns:
    --------
    particles_list : list
        Updated particle list with INP properties
    n_IN : int
        Number of particles assigned as IN-containing
    """
    n_IN = 0

    for particle in particles_list:
        # Skip already ice particles
        if particle.micro_type < 0:
            continue

        # Randomly assign IN to fraction of particles
        if np.random.random() < frac_IN:
            particle.micro_type = 2  # Liquid with IN

            # Sample INP radius from lognormal distribution
            r_IN = r_IN_mean * np.exp(np.random.normal(0, np.log(sigma_IN)))

            # INP surface area (m^2) = 4 * pi * r^2 * N_droplets
            particle.IN_sfc = 4.0 * pi * r_IN**2 * particle.A

            n_IN += 1

    return particles_list, n_IN


def deposition_nucleation(particles_list, T_parcel, P_parcel, q_parcel, dt, air_mass):
    """
    Deposition nucleation for uncoated INP.

    INP without liquid coating can nucleate ice directly from vapor.
    Based on SAM LCM implementation.

    Parameters:
    -----------
    particles_list : list
        List of particle objects
    T_parcel : float
        Parcel temperature (K)
    P_parcel : float
        Parcel pressure (Pa)
    q_parcel : float
        Water vapor mixing ratio (kg/kg)
    dt : float
        Time step (s)
    air_mass : float
        Air mass in parcel (kg)

    Returns:
    --------
    particles_list : list
        Updated particle list
    T_parcel : float
        Updated temperature (K)
    n_nucleated : int
        Number of particles nucleated
    """
    n_nucleated = 0
    dq_dep = 0.0

    # Only proceed if supersaturated with respect to ice
    e_sat_i = esati(T_parcel)
    e_a = q_parcel * P_parcel / (q_parcel + r_a / rv)
    supersat_ice = e_a / e_sat_i - 1.0

    if supersat_ice <= 0:
        return particles_list, T_parcel, n_nucleated

    for particle in particles_list:
        # Only process uncoated INP (micro_type == 2 but no liquid, i.e., aerosol with IN)
        # In our model, we use particles with IN_sfc > 0 but not yet activated
        # For simplicity, check if particle has IN and is still an aerosol (very small)
        if particle.micro_type != 2:
            continue

        # Check if particle is unactivated (small, haze-like)
        r_liq = (particle.M / (particle.A * 4.0/3.0 * pi * rho_liq))**(1.0/3.0)
        if r_liq > 1.0e-6:  # Already activated as droplet
            continue

        # Get INP surface area
        ice_sfc = particle.IN_sfc / particle.A
        if ice_sfc <= 0:
            continue

        # Deposition nucleation rate (simplified)
        # Using ABIFM-like formulation but for deposition
        # Critical supersaturation decreases with temperature
        T_C = T_parcel - 273.15
        S_crit = max(0.05, 0.3 + 0.01 * T_C)  # Critical S_ice for deposition

        if supersat_ice < S_crit:
            continue

        # Nucleation probability (temperature dependent)
        Js = 1e3 * np.exp(-0.1 * (T_C + 40))  # Simple parameterization
        p_nucleate = 1.0 - np.exp(-Js * dt * ice_sfc)

        # Stochastic nucleation
        if p_nucleate > np.random.random():
            # Create initial ice crystal (1 um sphere)
            r_ice = 1.0e-6
            particle.micro_type = -2  # Ice with IN
            particle.aaxis = r_ice
            particle.caxis = r_ice
            particle.dns = rho_ice
            particle.M = particle.A * 4.0/3.0 * pi * rho_ice * r_ice**3

            # Latent heat release (sublimation)
            dq_dep += particle.M

            n_nucleated += 1

    # Apply latent heat of sublimation
    if dq_dep > 0:
        T_parcel = T_parcel + dq_dep * l_s / cp / air_mass

    return particles_list, T_parcel, n_nucleated


def get_ice_stats(particles_list):
    """
    Get statistics of ice particles.

    Returns:
    --------
    dict with:
        n_ice : int - number of ice particles
        n_ice_hom : int - number of homogeneously frozen
        n_ice_het : int - number of heterogeneously frozen
        N_ice : float - total ice number (sum of A)
        M_ice : float - total ice mass (kg)
        r_ice_mean : float - mean ice radius (m)
    """
    n_ice = 0
    n_ice_hom = 0
    n_ice_het = 0
    N_ice = 0.0
    M_ice = 0.0
    r_sum = 0.0

    for p in particles_list:
        if p.micro_type < 0:
            n_ice += 1
            N_ice += p.A
            M_ice += p.M

            # Approximate radius from mass (assuming sphere with ice density)
            r = (p.M / (p.A * 4.0/3.0 * pi * rho_ice)) ** (1.0/3.0)
            r_sum += r * p.A

            if p.micro_type == -1:
                n_ice_hom += 1
            else:
                n_ice_het += 1

    r_ice_mean = r_sum / N_ice if N_ice > 0 else 0.0

    return {
        'n_ice': n_ice,
        'n_ice_hom': n_ice_hom,
        'n_ice_het': n_ice_het,
        'N_ice': N_ice,
        'M_ice': M_ice,
        'r_ice_mean': r_ice_mean
    }
