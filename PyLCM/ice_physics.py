"""
Ice Physics Module for PyLCM

Implements ice crystal microphysics following SAM LCM:
1. Terminal velocity (Böhm theory, Welss 2023 JAMES)
2. Ventilation coefficients (Chen and Lamb 1994)
3. Capacitance for spheroidal ice crystals
4. Deposition/sublimation growth

References:
- Welss (2023, JAMES): Böhm sedimentation theory
- Chen and Lamb (1994, JAS): Ice crystal growth and ventilation
- SAM LCM: micro_sedi.f90, micro_cond.f90
"""

import numpy as np
from PyLCM.parameters import *
from PyLCM.condensation import esatw, esati


# =========================================================================
# Ice crystal capacitance (Chen and Lamb 1994, Eqns 39-40)
# =========================================================================

def capacitance(aaxis, caxis):
    """
    Capacitance of spheroidal ice crystal.

    Based on Chen and Lamb (1994, JAS) Eqns. 39-40.

    Parameters:
    -----------
    aaxis : float
        Semi-axis in horizontal direction (m)
    caxis : float
        Semi-axis in vertical direction (m)

    Returns:
    --------
    cap : float
        Capacitance (m)
    """
    if aaxis <= 0 or caxis <= 0:
        return 0.0

    phi = caxis / aaxis  # aspect ratio

    if phi < 1.0:
        # Oblate spheroid (plate-like)
        epsilon = np.sqrt(1.0 - phi**2)
        cap = aaxis * epsilon / np.arcsin(epsilon)
    elif phi > 1.0:
        # Prolate spheroid (column-like)
        epsilon = np.sqrt(1.0 - phi**(-2))
        cap = caxis * epsilon / np.log((1.0 + epsilon) * phi)
    else:
        # Sphere
        cap = aaxis

    return cap


# =========================================================================
# Ice crystal ventilation (Chen and Lamb 1994, Eqns 29-30)
# =========================================================================

def ice_ventilation(aaxis, caxis, w_sedi, T, rho_air, switch_boehm=False):
    """
    Ventilation coefficients for ice crystals.

    Based on Chen and Lamb (1994, JAS) Eqns. 29-30.
    Optional shape-dependent correction from Welss and Seifert (2023).

    Parameters:
    -----------
    aaxis : float
        a-axis half-diameter (m)
    caxis : float
        c-axis half-diameter (m)
    w_sedi : float
        Terminal velocity (m/s)
    T : float
        Temperature (K)
    rho_air : float
        Air density (kg/m³)
    switch_boehm : bool
        Enable shape-dependent Böhm correction

    Returns:
    --------
    f_vent_mass : float
        Ventilation factor for mass growth
    f_vent_axes : float
        Ventilation factor for axis growth (c/a ratio)
    """
    if aaxis <= 0 or caxis <= 0:
        return 1.0, 1.0

    phi = caxis / aaxis

    # Dynamic viscosity (Sutherland's law)
    mu = muelq * ((273.15 + 110.4) / (T + 110.4)) * (T / 273.15)**1.5

    # Diffusivity of water vapor
    diff_coeff = 0.211e-4 * (T / 273.15)**1.94 * (101325.0 / (rho_air * r_a * T))

    # Schmidt number
    N_Sc = mu / (rho_air * diff_coeff)

    # Reynolds number
    D_char = 2.0 * max(aaxis, caxis)  # characteristic dimension
    N_Re = rho_air * w_sedi * D_char / mu

    # Ventilation parameter X (Chen and Lamb 1994)
    cap = capacitance(aaxis, caxis)
    if cap > 0:
        XX = N_Sc**(1.0/3.0) * np.sqrt(N_Re) * np.sqrt(D_char / cap)
    else:
        XX = 0.0

    # Ventilation coefficients (Chen and Lamb 1994, Eqns 29-30)
    if XX < 1.0:
        f_vent_mass = 1.0 + 0.14 * XX**2
        if cap > 0:
            f_vent_axes = (1.0 + 0.14 * XX**2 * (caxis / cap)**0.5) / \
                          (1.0 + 0.14 * XX**2 * (aaxis / cap)**0.5)
        else:
            f_vent_axes = 1.0
    else:
        f_vent_mass = 0.86 + 0.28 * XX
        if cap > 0:
            f_vent_axes = (0.86 + 0.28 * XX * (caxis / cap)**0.5) / \
                          (0.86 + 0.28 * XX * (aaxis / cap)**0.5)
        else:
            f_vent_axes = 1.0

    # Optional Böhm shape-dependent correction (Welss and Seifert 2023)
    if switch_boehm:
        if phi > 1.0:  # prolate (columns)
            f_vent_mass = f_vent_mass + 2.8e-2 * XX * phi
        elif phi < 1.0:  # oblate (plates)
            f_vent_mass = f_vent_mass + 2.8e-3 * XX**(3.0/2.0) / phi

    return f_vent_mass, f_vent_axes


# =========================================================================
# Ice crystal terminal velocity (Böhm theory, Welss 2023)
# =========================================================================

def N_Re_Boehm(phi, q, m_ice, mu_air, rho_air):
    """
    Reynolds number using Böhm theory (Welss 2023, JAMES, Appendix A).

    Parameters:
    -----------
    phi : float
        Aspect ratio (c/a)
    q : float
        Area ratio (projected area / circumscribed ellipse area)
    m_ice : float
        Ice crystal mass (kg)
    mu_air : float
        Dynamic viscosity of air (Pa s)
    rho_air : float
        Air density (kg/m³)

    Returns:
    --------
    Re : float
        Reynolds number
    """
    X0 = 2.8e6  # Best number reference

    # Best number (A1a)
    X = 8.0 * m_ice * g * rho_air / (pi * mu_air**2 * max(phi, 1.0) * q**0.25)

    # k(phi) (A1b)
    k1 = max(0.82 + 0.18 * phi, 0.85)
    k2 = 0.37 + 0.63 / phi
    k3 = 1.33 / (max(np.log10(phi), 0.0) + 1.19)
    k_phi = min(k1, k2, k3)

    # Gamma_phi (A1c)
    Gamma_phi = max(1.0, min(1.98, 3.76 - 8.41*phi + 9.18*phi**2 - 3.53*phi**3))

    # CDP,S (A1d)
    CDPs = max(0.292 * k_phi * Gamma_phi, 0.492 - 0.2 / np.sqrt(phi))

    # CDP (A1e-g)
    CDP = max(1.0, q * (1.46 * q - 0.46) * CDPs)
    CDPp = CDP * (1.0 + 1.6 * (X/X0)**2) / (1.0 + (X/X0)**2)

    # Auxiliary coefficients (A1g-i)
    CD0 = 4.5 * k_phi**2 * max(phi, 1.0)
    beta = np.sqrt(1.0 + CDP / (6.0 * k_phi) * np.sqrt(X / CDPp)) - 1.0
    gamma = (CD0 - CDP) / (4.0 * CDP)

    # Reynolds number (A1j)
    Re = 6.0 * k_phi / CDPp * beta**2 * \
         (1.0 + (2.0 * beta * np.exp(-beta * gamma) / ((2.0 + beta) * (1.0 + beta))))

    return Re


def ws_ice_boehm(aaxis, caxis, dns, mass, T, rho_air):
    """
    Ice crystal terminal velocity using Böhm theory.

    Based on Welss (2023, JAMES) and Böhm sedimentation.

    Parameters:
    -----------
    aaxis : float
        a-axis half-diameter (m)
    caxis : float
        c-axis half-diameter (m)
    dns : float
        Apparent ice density (kg/m³)
    mass : float
        Ice crystal mass per particle (kg)
    T : float
        Temperature (K)
    rho_air : float
        Air density (kg/m³)

    Returns:
    --------
    wsedi : float
        Terminal velocity (m/s)
    """
    if aaxis <= 0 or caxis <= 0 or mass <= 0:
        return 0.0

    # Dynamic viscosity (Sutherland's law)
    mu = muelq * ((273.15 + 110.4) / (T + 110.4)) * (T / 273.15)**1.5

    phi = caxis / aaxis  # Aspect ratio

    # Projected area calculation (Welss Eq. 17)
    if phi > 1.0:
        # Prolate (column-like)
        A_CE = pi * aaxis * caxis  # circumscribed ellipse
        Aproj = A_CE
    else:
        # Oblate (plate-like)
        A_CE = pi * aaxis**2
        Aproj = A_CE * ((1 - phi) * (dns / rho_ice) + phi)

    q = Aproj / A_CE if A_CE > 0 else 1.0

    # Reynolds number from Böhm theory
    Re_pro = N_Re_Boehm(phi, q, mass, mu, rho_air)

    # Terminal velocity
    Dchar = 2.0 * aaxis  # equatorial diameter
    wsedi_pro = mu * Re_pro / (rho_air * Dchar)

    # For columns (phi > 1), blend with cylinder solution
    if phi > 1.0:
        # Cylinder branch (McCorquodale & Westbrook 2021b)
        q_cyl = 4.0 / (pi * phi)
        Re_cyl = N_Re_Boehm(phi, q_cyl, mass, mu, rho_air)
        wsedi_cyl = mu * Re_cyl / (rho_air * Dchar)

        # Blend prolate and cylinder (Welss Eqns 31-32)
        fx = np.exp(-0.3 * (phi - 1.0))  # 0 < fx <= 1
        wsedi = fx * wsedi_pro + (1.0 - fx) * wsedi_cyl
    else:
        wsedi = wsedi_pro

    return wsedi


# =========================================================================
# Deposition/sublimation growth
# =========================================================================

def ice_deposition_growth(particle, T, P, q, dt, f_vent_mass=1.0):
    """
    Grow/shrink ice crystal by deposition/sublimation.

    Based on Chen and Lamb (1994) formulation.

    Parameters:
    -----------
    particle : object
        Particle with M, A, aaxis, caxis, dns attributes
    T : float
        Temperature (K)
    P : float
        Pressure (Pa)
    q : float
        Water vapor mixing ratio (kg/kg)
    dt : float
        Time step (s)
    f_vent_mass : float
        Ventilation factor for mass

    Returns:
    --------
    dm : float
        Mass change (kg)
    """
    # Saturation vapor pressures
    e_sat_i = esati(T)
    e_a = q * P / (q + r_a / rv)  # actual vapor pressure

    # Supersaturation with respect to ice
    supersat_ice = e_a / e_sat_i - 1.0

    # Ice crystal properties
    mass = particle.M / particle.A  # mass per crystal
    cap = capacitance(particle.aaxis, particle.caxis)

    # Thermal conductivity
    thermal_conductivity = 7.94048e-05 * T + 0.00227011

    # Diffusivity
    diff_coeff = 0.211e-4 * (T / 273.15)**1.94 * (101325.0 / P)

    # G factor for ice (similar to liquid but with sublimation latent heat)
    qsat_i = e_sat_i / (P - e_sat_i) * r_a / rv
    G_ice = 1.0 / (rho_ice * rv * T / (e_sat_i * diff_coeff) +
                   (l_s / (rv * T) - 1.0) * rho_ice * l_s / (thermal_conductivity * T))

    # Mass growth rate: dm/dt = 4*pi*C*G*f_vent*S_i
    # Using m^(2/3) formulation for stability
    if mass > 0 and cap > 0:
        dm23dt = 8.0/3.0 * pi * G_ice * (cap / mass**(1.0/3.0)) * f_vent_mass * supersat_ice

        # New mass
        m_new_23 = mass**(2.0/3.0) + dt * dm23dt
        m_min = (4.0/3.0 * pi * rho_ice * (0.001e-6)**3)**(2.0/3.0)  # minimum size

        m_new = max(m_new_23, m_min)**(3.0/2.0)
        dm = (m_new - mass) * particle.A
    else:
        dm = 0.0

    return dm


def update_ice_axes(particle, dm, T, f_vent_axes=1.0):
    """
    Update ice crystal axes after mass change.

    Based on Chen and Lamb (1994) inherent growth ratio.

    Parameters:
    -----------
    particle : object
        Particle with M, A, aaxis, caxis, dns attributes
    dm : float
        Mass change (kg)
    T : float
        Temperature (K)
    f_vent_axes : float
        Ventilation factor for axis growth
    """
    if dm == 0 or particle.aaxis <= 0 or particle.caxis <= 0:
        return

    mass = particle.M / particle.A
    dm_per = dm / particle.A

    # Inherent growth ratio gamma (Chen and Lamb 1994, Fig. 3)
    gamma_ice = get_gamma_ice(T)

    # For small crystals (< 10 um), grow as sphere
    if min(particle.aaxis, particle.caxis) < 10.0e-6:
        r_new = ((mass + dm_per) / (4.0/3.0 * pi * rho_ice))**(1.0/3.0)
        particle.aaxis = r_new
        particle.caxis = r_new
        particle.dns = rho_ice
    else:
        # Volume change
        V_old = 4.0/3.0 * pi * particle.aaxis**2 * particle.caxis
        V_new = (mass + dm_per) / particle.dns

        if V_old > 0:
            dlogV = (V_new - V_old) / V_old

            # Adjusted gamma with ventilation
            gamma_ast = gamma_ice * f_vent_axes

            if gamma_ast > 1.0:
                # Column growth (c grows faster)
                particle.caxis = particle.caxis * (1.0 + dlogV * gamma_ast / (2.0 + gamma_ast))
                particle.aaxis = particle.aaxis * (1.0 + dlogV / (2.0 + gamma_ast))
            else:
                # Plate growth (a grows faster)
                particle.aaxis = particle.aaxis * (1.0 + dlogV / (2.0 * gamma_ast + 1.0))
                particle.caxis = particle.caxis * (1.0 + dlogV * gamma_ast / (2.0 * gamma_ast + 1.0))


def get_gamma_ice(T):
    """
    Inherent growth ratio gamma as function of temperature.

    Lookup table from Chen and Lamb (1994, Fig. 3).

    Parameters:
    -----------
    T : float
        Temperature (K)

    Returns:
    --------
    gamma : float
        Inherent growth ratio
    """
    T_C = T - 273.15

    # Simplified parameterization based on Chen and Lamb (1994)
    if T_C > -4.0:
        gamma = 1.0  # isometric
    elif T_C > -8.0:
        # Plates
        gamma = 0.5
    elif T_C > -12.0:
        # Columns
        gamma = 2.0
    elif T_C > -18.0:
        # Plates
        gamma = 0.3
    elif T_C > -22.0:
        # Columns/sector plates
        gamma = 1.5
    else:
        # Plates and columns
        gamma = 1.0

    return gamma


# =========================================================================
# Main ice physics routine
# =========================================================================

def ice_growth_and_sedimentation(particles_list, T_parcel, P_parcel, q_parcel,
                                  dt, air_mass, switch_deposition=True,
                                  switch_sedimentation=True):
    """
    Main routine for ice crystal physics.

    Applies deposition growth and calculates terminal velocity.

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
    switch_deposition : bool
        Enable deposition growth
    switch_sedimentation : bool
        Calculate terminal velocities

    Returns:
    --------
    particles_list : list
        Updated particles
    T_parcel : float
        Updated temperature
    q_parcel : float
        Updated water vapor
    """
    rho_air = P_parcel / (r_a * T_parcel)
    dq_dep = 0.0  # Net deposition

    for particle in particles_list:
        # Skip liquid particles
        if particle.micro_type > 0:
            continue

        # Get mass per ice crystal
        mass = particle.M / particle.A

        # Calculate terminal velocity
        if switch_sedimentation and particle.aaxis > 0 and particle.caxis > 0:
            particle.wsedi = ws_ice_boehm(
                particle.aaxis, particle.caxis, particle.dns,
                mass, T_parcel, rho_air
            )
        else:
            particle.wsedi = 0.0

        # Ventilation
        if particle.wsedi > 0:
            f_vent_mass, f_vent_axes = ice_ventilation(
                particle.aaxis, particle.caxis, particle.wsedi,
                T_parcel, rho_air, switch_boehm=False
            )
        else:
            f_vent_mass, f_vent_axes = 1.0, 1.0

        # Deposition growth
        if switch_deposition:
            M_old = particle.M

            dm = ice_deposition_growth(
                particle, T_parcel, P_parcel, q_parcel, dt, f_vent_mass
            )

            # Update mass
            particle.M = particle.M + dm

            # Update axes
            update_ice_axes(particle, dm, T_parcel, f_vent_axes)

            # Track water vapor change
            dq_dep += dm

    # Update parcel thermodynamics
    if abs(dq_dep) > 0:
        # Latent heat release from deposition
        T_parcel = T_parcel + dq_dep * l_s / cp / air_mass
        # Water vapor decrease
        q_parcel = q_parcel - dq_dep / air_mass

    return particles_list, T_parcel, q_parcel


def get_ice_terminal_velocities(particles_list, T, P):
    """
    Calculate terminal velocities for all ice particles.

    Parameters:
    -----------
    particles_list : list
        List of particles
    T : float
        Temperature (K)
    P : float
        Pressure (Pa)

    Returns:
    --------
    particles_list : list
        Updated with wsedi attribute
    """
    rho_air = P / (r_a * T)

    for particle in particles_list:
        if particle.micro_type < 0:  # Ice
            mass = particle.M / particle.A
            particle.wsedi = ws_ice_boehm(
                particle.aaxis, particle.caxis, particle.dns,
                mass, T, rho_air
            )
        else:
            # Use existing droplet terminal velocity
            pass

    return particles_list


# =========================================================================
# Melting (instantaneous when T > 273.15 K)
# =========================================================================

def ice_melting(particles_list, T_parcel, air_mass):
    """
    Instantaneous melting of ice particles when T > 0°C.

    Following SAM LCM: instant conversion when T > 273.15 K.

    Parameters:
    -----------
    particles_list : list
        List of particle objects
    T_parcel : float
        Parcel temperature (K)
    air_mass : float
        Air mass in parcel (kg)

    Returns:
    --------
    particles_list : list
        Updated particles (melted ice → liquid)
    T_parcel : float
        Updated temperature (after latent heat absorption)
    n_melted : int
        Number of particles melted
    """
    n_melted = 0
    dq_melt = 0.0  # Mass that melted (absorbs latent heat)

    # Only proceed if T above melting point
    if T_parcel <= 273.15:
        return particles_list, T_parcel, n_melted

    for particle in particles_list:
        # Skip if not ice
        if particle.micro_type >= 0:
            continue

        # Convert ice to liquid
        # micro_type: -1 → +1 (hom ice → liquid), -2 → +2 (het ice → liquid with IN)
        particle.micro_type = abs(particle.micro_type)

        # Reset ice properties
        particle.aaxis = 0.0
        particle.caxis = 0.0
        particle.dns = rho_liq
        particle.wsedi = 0.0

        # Track latent heat absorption
        dq_melt += particle.M

        n_melted += 1

    # Absorb latent heat of fusion (cooling)
    if dq_melt > 0:
        T_parcel = T_parcel - dq_melt * l_f / cp / air_mass

    return particles_list, T_parcel, n_melted


# =========================================================================
# Full gamma_ice lookup table (Chen and Lamb 1994, Fig. 3)
# =========================================================================

def get_gamma_ice_full(T):
    """
    Full inherent growth ratio gamma lookup table.

    Based on Chen and Lamb (1994, Fig. 3) as implemented in SAM LCM.
    Resolution: 0.5°C from -30°C to 0°C.

    Parameters:
    -----------
    T : float
        Temperature (K)

    Returns:
    --------
    gamma : float
        Inherent growth ratio
    """
    # Temperature values from -30°C to 0°C in 0.5°C steps
    T_values = np.arange(-30.0, 0.5, 0.5)

    # Gamma values from Chen and Lamb (1994) Fig. 3
    # (SAM LCM lookup table)
    G_values = np.array([
        1.2845082, 1.3063233, 1.334114, 1.3740169, 1.4210837, 1.4697628, 1.5329632,
        1.5988813, 1.6817352, 1.7688825, 1.8605458, 1.9405501, 2.0070235, 2.0325339,
        2.0161880, 1.9648150, 1.8784220, 1.7839930, 1.6775080, 1.5872910, 1.5024380,
        1.4172130, 1.3491790, 1.2811440, 1.2131100, 1.1564270, 1.0997450, 1.0430620,
        0.9872024, 0.9304980, 0.8737940, 0.8170900, 0.7603860, 0.7172420, 0.6883860,
        0.6595300, 0.6306750, 0.6127530, 0.6150940, 0.6174360, 0.6197780, 0.6320230,
        0.6691700, 0.7063170, 0.7434640, 0.7791290, 0.8124860, 0.8458420, 0.8791990,
        0.9065940, 0.9243790, 0.9421640, 0.9599490, 0.9777340, 0.9955190, 1.0133040,
        1.0275750, 1.0135100, 0.9766440, 0.9517800, 0.9090870
    ])

    T_C = T - 273.15

    if T_C < -30.0:
        return G_values[0]
    elif T_C >= 0.0:
        return G_values[-1]
    else:
        # Linear interpolation
        idx = int((T_C + 30.0) * 2)
        idx = max(0, min(idx, len(G_values) - 2))
        frac = (T_C - T_values[idx]) / 0.5
        return G_values[idx] + frac * (G_values[idx + 1] - G_values[idx])
