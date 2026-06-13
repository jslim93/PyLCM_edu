import math
import numpy as np
from numba import jit
from PyLCM.micro_particle import *
from PyLCM.parcel import *
from PyLCM.condensation import *
from tqdm import tqdm
import itertools

# Module-level Beard (1976) polynomial coefficients (float64) so the jitted
# ws_drops_beard can close over them as compile-time constants.
_BEARD_B = np.array([-0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3,
                     -0.578878e-3, 0.855176e-4, -0.327815e-5], dtype=np.float64)
_BEARD_C = np.array([-0.500015e1, 0.523778e1, -0.204914e1, 0.475294,
                     -0.542819e-1, 0.238449e-2], dtype=np.float64)

def collection(dt, particles_list, rho_parcel, rho_liq, p_env, T_parcel, acc_ts, aut_ts, precip_ts, sedi_removal, z_parcel, max_z, w_parcel, switch_E_constant=False, switch_vt_simple=False, switch_turb_kernel=False, epsilon_turb=0.0):

    #shuffle the particle list for LSM (linear sampling method)
    # Use np.random.permutation so the shuffle is (a) faster than random.shuffle
    # of a Python list and (b) controlled by np.random.seed (the ensemble seeds
    # numpy's RNG, and determine_collision already draws from np.random).
    perm = np.random.permutation(len(particles_list))
    particles_list = [particles_list[i] for i in perm]
    nptcl = len(particles_list)
    half_length = len(particles_list) // 2
    
    particle_list1 = particles_list[:half_length]  # Splitting into the first half
    particle_list2 = particles_list[half_length:]  # Splitting into the second half
    
    #Collisions are considered between two shuffled particle lists
    for particle1, particle2 in zip(particle_list1,particle_list2):
        
        #  A superdroplet must contain at least one real particle to collect other droplets
        if min(particle1.A, particle2.A) <= 0:
            continue
        
        # The larger droplet should be larger than 10.0 µm to cause collisions.
        if max(particle1.M / particle1.A, particle2.M / particle2.A) < (10.0E-6 ** 3) * 4.0 / 3.0 * np.pi * rho_liq:
            continue

        check_final = False
        check_collection = False
        
        # Find out what kind of interaction (or none) takes place
        check_final, check_collection, v_r1, v_r2, p_crit = determine_collision(dt,particle1, particle2, rho_parcel, rho_liq, p_env, T_parcel, half_length,nptcl, switch_E_constant, switch_vt_simple, switch_turb_kernel, epsilon_turb)

        if check_final:

            # A special treatment is necessary if weighting factors are identical
            if particle1.A == particle2.A:

                #add acc_ts and aut_ts for same weight factor
                particle1, particle2, acc_ts, aut_ts = same_weights_update(particle1, particle2, acc_ts, aut_ts)

            # Each droplet of the super-droplet with the smaller weighting factor collects
            # one droplet of the super-droplet with the larger weighting factor
            elif check_collection:
                particle1, particle2, acc_ts, aut_ts = liquid_update_collection(particle1, particle2, acc_ts, aut_ts, p_crit)
        
        #Change droplet vertical location by subtracting their terminal velocity from the parcel velocity˛
        if sedi_removal:
            if z_parcel < max_z: 
                dz_ptcl = w_parcel * dt
            else:
                dz_ptcl = 0.0
            particle1.z = particle1.z - (v_r1 * dt) + dz_ptcl
            particle2.z = particle2.z - (v_r2 * dt) + dz_ptcl
        #Remove particles at the surface
            if particle1.z <= 0.0:
                precip_ts += particle1.M
                particle1.A = 0
            if particle2.z <= 0.0:
                precip_ts += particle2.M

                particle2.A = 0
                
    # Merge the lists at the end of the loop
    particles_list = particle_list1 + particle_list2

    # Remove particle with 0 weighting factor
    particles_list = [particle for particle in particles_list if particle.A > 0]

    #if collision_timestep_error:
    #    print('+++ Collision time step is too long. +++')
    #    collision_timestep_error = False
    
    return particles_list, acc_ts, aut_ts, precip_ts

def liquid_update_collection(particle1, particle2, acc_ts, aut_ts, p_crit=1):

    # _int1: gains total individual mass (smaller A)
    # _int2: loses total mass, constant individual mass (larger A)
    if particle1.A < particle2.A:
        ptcl_int1 = particle1
        ptcl_int2 = particle2
    else:
        ptcl_int1 = particle2
        ptcl_int2 = particle1

    x_int = ptcl_int2.M / ptcl_int2.A
    xs_int = ptcl_int2.Ns / ptcl_int2.A

    #Droplet volume for kappa collision estimation
    v_ptcl1 = ptcl_int1.M / ptcl_int1.A / rho_liq
    v_ptcl2 = ptcl_int2.M / ptcl_int2.A / rho_liq

    # Update of M, A (water mass and particle number)
    # p_crit: number of collisions (multi-collision, Shima et al. 2009)

    #Increase of water mass due to collision (p_crit collisions)
    ptcl_int1.M = ptcl_int1.M + ptcl_int1.A * x_int * p_crit
    #Increase of Aerosol mass due to collision
    ptcl_int1.Ns = ptcl_int1.Ns + ptcl_int1.A * xs_int * p_crit
    #Update volume-mean averaged kappa
    ptcl_int1.kappa = (v_ptcl1*ptcl_int1.kappa + v_ptcl2*ptcl_int2.kappa )/ (v_ptcl1 + v_ptcl2)

    #Decrease of number, aerosol and water mass due to collision
    ptcl_int2.A  = ptcl_int2.A - ptcl_int1.A * p_crit
    ptcl_int2.M  = ptcl_int2.M - ptcl_int1.A * x_int * p_crit
    ptcl_int2.Ns = ptcl_int2.Ns - ptcl_int1.A * xs_int * p_crit
    
    mass_crit = (seperation_radius_ts ** 3) * 4.0 / 3.0 * np.pi * rho_liq
    
    large_drop_size = max(particle1.M / particle1.A, particle2.M / particle2.A)
    small_drop_size = min(particle1.M / particle1.A, particle2.M / particle2.A)
    
    #Accretion mass
    if (large_drop_size >= mass_crit) and (small_drop_size < mass_crit) :
        acc_ts += ptcl_int1.A * small_drop_size
    #Autoconversion mass
    if (large_drop_size < mass_crit) and (small_drop_size < mass_crit) and (small_drop_size + large_drop_size >= mass_crit):
        aut_ts += ptcl_int1.A * (large_drop_size + small_drop_size)
    
    # The superdroplet with the smaller A will be indexed particle1 in the following (l. 51 in the Fortran)
    if particle1.A < particle2.A:
        particle1 =  ptcl_int1
        particle2 =  ptcl_int2 
    else:
        particle1 =  ptcl_int2
        particle2 =  ptcl_int1
    
    return(particle1, particle2, acc_ts, aut_ts)

def same_weights_update(ptcl_int1, ptcl_int2, acc_ts, aut_ts):
    
    mass_crit = (seperation_radius_ts ** 3) * 4.0 / 3.0 * np.pi * rho_liq
    
    large_drop_size = max(ptcl_int1.M / ptcl_int1.A, ptcl_int2.M / ptcl_int2.A)
    small_drop_size = min(ptcl_int1.M / ptcl_int1.A, ptcl_int2.M / ptcl_int2.A)
    
    #Droplet volume for kappa collision estimation
    v_ptcl1 = ptcl_int1.M / ptcl_int1.A / rho_liq
    v_ptcl2 = ptcl_int2.M / ptcl_int2.A / rho_liq
    
    #Autoconversion mass
    if (large_drop_size < mass_crit) and (small_drop_size < mass_crit) and (small_drop_size + large_drop_size >= mass_crit):
        aut_ts += ptcl_int1.M + ptcl_int2.M
        
    #Accretion mass
    if (large_drop_size > mass_crit) and (small_drop_size < mass_crit):
        acc_ts += small_drop_size * ptcl_int1.A

    # Following Fortran reference: individual merged droplet mass = x1 + x2
    # Halve weighting factors, assign combined individual mass to both
    xn = ptcl_int1.M / ptcl_int1.A   # individual mass of particle 1
    xm = ptcl_int2.M / ptcl_int2.A   # individual mass of particle 2
    xsn = ptcl_int1.Ns / ptcl_int1.A # individual aerosol mass of particle 1
    xsm = ptcl_int2.Ns / ptcl_int2.A # individual aerosol mass of particle 2

    ptcl_int1.A  = ptcl_int1.A * 0.5
    ptcl_int2.A  = ptcl_int2.A - ptcl_int1.A  # handles rounding like Fortran FLOOR

    ptcl_int1.M  = (xn + xm) * ptcl_int1.A
    ptcl_int2.M  = (xn + xm) * ptcl_int2.A
    ptcl_int1.Ns = (xsn + xsm) * ptcl_int1.A
    ptcl_int2.Ns = (xsn + xsm) * ptcl_int2.A
    
    #Update volume-mean averaged kappa (compute first, then assign to avoid sequential mutation)
    new_kappa = (v_ptcl1*ptcl_int1.kappa + v_ptcl2*ptcl_int2.kappa) / (v_ptcl1 + v_ptcl2)
    ptcl_int1.kappa = new_kappa
    ptcl_int2.kappa = new_kappa

    return(ptcl_int1, ptcl_int2, acc_ts, aut_ts)

import math
import numpy as np
def determine_collision(dt, particle1, particle2, rho_parcel, rho_liq, p_env, T_parcel, half_length,nptcl, switch_E_constant=False, switch_vt_simple=False, switch_turb_kernel=False, epsilon_turb=0.0):
    # V_parcel = air_mass / rho_parcel, consistent with parcel_rho (PARCEL_AIR_MASS)
    V_parcel = PARCEL_AIR_MASS / rho_parcel

    check_final = False
    check_collection = False

    # A fully evaporated/degenerate super-droplet (M<=0 or A<=0) has zero radius and
    # cannot collide; skip it to avoid division-by-zero in radius/E_H80 computation.
    if particle1.M <= 0.0 or particle2.M <= 0.0 or particle1.A <= 0.0 or particle2.A <= 0.0:
        return check_final, check_collection, 0.0, 0.0, 0

    R_n = (particle1.M / particle1.A / (4.0 / 3.0 * pi * rho_liq)) ** 0.33333333333
    R_m = (particle2.M / particle2.A / (4.0 / 3.0 * pi * rho_liq)) ** 0.33333333333

    # Terminal velocity: Beard (1976) or simplified Stokes drag
    if switch_vt_simple:
        v_r1 = ws_drops_stokes(R_n, rho_parcel, rho_liq)
        v_r2 = ws_drops_stokes(R_m, rho_parcel, rho_liq)
    else:
        v_r1 = ws_drops_beard(R_n, rho_parcel, rho_liq, p_env, T_parcel)
        v_r2 = ws_drops_beard(R_m, rho_parcel, rho_liq, p_env, T_parcel)

    v_r = abs(v_r1 - v_r2)

    # Collision efficiency: Hall (1980) lookup or constant E=1
    E_coll = 1.0 if switch_E_constant else E_H80(R_m, R_n)

    # Collection kernel: standard gravitational or Wang-Ayala turbulent
    if switch_turb_kernel and epsilon_turb > 1.0e-10:
        # Estimate TKE from epsilon: urms = 2.02*(eps/0.04)^(1/3), tke = 1.5*urms^2
        urms_est = 2.02 * (epsilon_turb / 0.04)**(1.0 / 3.0)
        tke_est = 1.5 * urms_est**2
        K = E_coll * gck(R_n, R_m, v_r1, v_r2, epsilon_turb, tke_est) * E_turb(R_n, R_m, epsilon_turb) * E_S09(R_m, R_n, v_r, rho_liq, T_parcel)
    else:
        K = pi * (R_m + R_n) ** 2 * v_r * E_coll * E_S09(R_m, R_n, v_r, rho_liq, T_parcel)

    p_crit = max(particle1.A, particle2.A) * K / V_parcel * dt
    p_crit = p_crit*nptcl*(nptcl-1)/(half_length*2)

    x_rand = np.random.random()

    if p_crit > x_rand:
        check_final = True
        check_collection = True

        if p_crit <= 1.0:
            p_crit = 1
        else:
            # Multi-collision: following SAM Fortran (Shima et al. 2009)
            p_crit = max(round(p_crit), 1)  # NINT equivalent
            # Limiter: ensure we don't collect more particles than available
            # max(A) - p_crit * min(A) >= 1
            A_max = max(particle1.A, particle2.A)
            A_min = min(particle1.A, particle2.A)
            p_crit_lim = max(int((A_max - 1) / A_min), 1)
            p_crit = min(p_crit, p_crit_lim)
    else:
        p_crit = 0

    return check_final, check_collection, v_r1, v_r2, p_crit


# Collision efficiencies by Hall (1980) — module-level float64 tables so the
# jitted E_H80 closes over them as compile-time constants.
_H80_R0 = np.array([6.0, 8.0, 10.0, 15.0, 20.0, 25.0,
                    30.0, 40.0, 50.0, 60.0, 70.0, 100.0,
                    150.0, 200.0, 300.0], dtype=np.float64)

_H80_RAT = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25,
                     0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                     0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
                     0.90, 0.95, 1.00], dtype=np.float64)

_H80_ECOLL = np.ascontiguousarray(np.array([
    [ 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001],
    [ 0.003, 0.003, 0.003, 0.004, 0.005, 0.005, 0.005, 0.010, 0.100, 0.050, 0.200, 0.500, 0.770, 0.870, 0.970],
    [ 0.007, 0.007, 0.007, 0.008, 0.009, 0.010, 0.010, 0.070, 0.400, 0.430, 0.580, 0.790, 0.930, 0.960, 1.000 ],
    [ 0.009, 0.009, 0.009, 0.012, 0.015, 0.010, 0.020, 0.280, 0.600, 0.640, 0.750, 0.910, 0.970, 0.980, 1.000 ],
    [ 0.014, 0.014, 0.014, 0.015, 0.016, 0.030, 0.060, 0.500, 0.700, 0.770, 0.840, 0.950, 0.970, 1.000, 1.000 ],
    [ 0.017, 0.017, 0.017, 0.020, 0.022, 0.060, 0.100, 0.620, 0.780, 0.840, 0.880, 0.950, 1.000, 1.000, 1.000 ],
    [ 0.030, 0.030, 0.024, 0.022, 0.032, 0.062, 0.200, 0.680, 0.830, 0.870, 0.900, 0.950, 1.000, 1.000, 1.000 ],
    [ 0.025, 0.025, 0.025, 0.036, 0.043, 0.130, 0.270, 0.740, 0.860, 0.890, 0.920, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.027, 0.027, 0.027, 0.040, 0.052, 0.200, 0.400, 0.780, 0.880, 0.900, 0.940, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.030, 0.030, 0.030, 0.047, 0.064, 0.250, 0.500, 0.800, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.040, 0.040, 0.033, 0.037, 0.068, 0.240, 0.550, 0.800, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.035, 0.035, 0.035, 0.055, 0.079, 0.290, 0.580, 0.800, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.037, 0.037, 0.037, 0.062, 0.082, 0.290, 0.590, 0.780, 0.900, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.037, 0.037, 0.037, 0.060, 0.080, 0.290, 0.580, 0.770, 0.890, 0.910, 0.950, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.037, 0.037, 0.037, 0.041, 0.075, 0.250, 0.540, 0.760, 0.880, 0.920, 0.950, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.037, 0.037, 0.037, 0.052, 0.067, 0.250, 0.510, 0.770, 0.880, 0.930, 0.970, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.037, 0.037, 0.037, 0.047, 0.057, 0.250, 0.490, 0.770, 0.890, 0.950, 1.000, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.036, 0.036, 0.036, 0.042, 0.048, 0.230, 0.470, 0.780, 0.920, 1.000, 1.020, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.040, 0.040, 0.035, 0.033, 0.040, 0.112, 0.450, 0.790, 1.010, 1.030, 1.040, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.033, 0.033, 0.033, 0.033, 0.033, 0.119, 0.470, 0.950, 1.300, 1.700, 2.300, 1.000, 1.000, 1.000, 1.000 ],
    [ 0.027, 0.027, 0.027, 0.027, 0.027, 0.125, 0.520, 1.400, 2.300, 3.000, 4.000, 1.000, 1.000, 1.000, 1.000]
], dtype=np.float64).T)


@jit(nopython=True, cache=True)
def E_H80(r1, r2):
    # Collision efficiencies by Hall (1980)
    r0 = _H80_R0
    rat = _H80_RAT
    ecoll = _H80_ECOLL

    # Calculate the radius class index of particles with respect to array r0
    # Radius has to be in microns
    rmax = max(r1,r2)
    ir = 15
    if ( rmax * 1.0E6 >= r0[14] ):
        ir = 15
    else:
        for k in range(15):
            if rmax*1e6 < r0[k]:
                ir = k
                break

    # Two-dimensional linear interpolation of the collision efficiency
    # Radius has to be in microns
    rq = min(r1 / r2, r2 / r1)
    iq = int(rq * 20) 
    iq = max(iq, 1)

    if ir < 15:
        if ir >= 1:
            pp = (rmax * 1.0E6 - r0[ir-1]) / (r0[ir] - r0[ir-1])
            qq = (rq - rat[iq-1]) / (rat[iq] - rat[iq-1])
            E = (1.0 - pp) * (1.0 - qq) * ecoll[ir-1, iq-1] + pp * (1.0 - qq) * ecoll[ir, iq-1] \
                + qq * (1.0 - pp) * ecoll[ir-1, iq] + pp * qq * ecoll[ir, iq]
        else:
            qq = (rq - rat[iq-1]) / (rat[iq] - rat[iq-1])
            E = (1.0 - qq) * ecoll[0, iq-1] + qq * ecoll[0, iq]
    else:
        qq = (rq - rat[iq-1]) / (rat[iq] - rat[iq-1])
        E = min((1.0 - qq) * ecoll[14, iq-1] + qq * ecoll[14, iq], 1.0)
    
    if( E < 1.0E-20 ):  
        E = 0.0
    E = max(E, 0.0)

    return E


def E_S09(r_m, r_n, v_r, rho_liq,t_parcel):
    #Coalescence efficiencies following Straub et al. (2009)

    # Compute the diameters
    d_L = 2.0 * max(r_m, r_n)
    d_S = 2.0 * min(r_m, r_n)

    # Compute collision kinetic energy
    CKE = (math.pi / 12.0) * rho_liq * d_L**3 * d_S**3 / (d_L**3 + d_S**3) * v_r**2

    # Compute surface energy
    S_c = math.pi * sigma_air_liq(t_parcel) * (d_L**3 + d_S**3)**(2/3)

    # Compute Weber number
    We = CKE / S_c

    # Compute coalescence efficiency
    e_s09 = math.exp(-1.15 * We)

    return e_s09

@jit(nopython=True, cache=True)
def ws_drops_beard(radius, rho_parcel, rho_liq, p_env, T_parcel):

    # Calculate the terminal velocity of a water droplet in air
    # Droplet terminal velocity (Beard, 1976, J. Atmos. Sci.).
    # T_parcel, rho_parcel, p_env must be provided

    b = _BEARD_B
    c = _BEARD_C

    eta0 = 1.818e-5
    l0 = 6.62e-8
    p0 = 1013.25
    T0 = 293.15
    rho0 = 1.292509

    diameter = max(2.0 * radius, 0.1e-6) # set minimum value to prevent dividing by zero

    eta = rho_parcel * eta0 / rho0
    l = l0 * (eta / eta0) * (p0 / p_env) * math.sqrt(T_parcel / T0)
    Cac = 1.0 + 2.5 * l / diameter

    if diameter <= 19.0e-6:
        C1 = (rho_liq - rho_parcel) * g / (18.0 * eta)
        ws_drops_beard = C1 * Cac * diameter**2
    elif diameter <= 1070.0e-6:
        C2 = 4.0 * rho_parcel * (rho_liq - rho_parcel) * g / (3.0 * eta**2)
        NDa = C2 * diameter**3
        XX = math.log(NDa)
        YY = 0.0
        for i in range(b.shape[0]):
            YY += b[i] * XX**i
        NRe = Cac * math.exp(YY)
        ws_drops_beard = eta * NRe / (rho_parcel * diameter)
    else:
        C3 = 4.0 * (rho_liq - rho_parcel) * g / (3.0 * sigma_air_liq(T_parcel))
        Bo = C3 * diameter**2
        NP = sigma_air_liq(T_parcel)**3 * rho_parcel**2 / (eta**4 * (rho_liq - rho_parcel) * g)
        XX = math.log(Bo * NP**0.166666666)
        YY = 0.0
        for i in range(c.shape[0]):
            YY += c[i] * XX**i
        NRe = NP**0.166666666 * math.exp(YY)
        ws_drops_beard = eta * NRe / (rho_parcel * diameter)

    return ws_drops_beard

def ws_drops_stokes(radius, rho_parcel, rho_liq):
    # Simplified Stokes terminal velocity (Ablation Lab mode)
    # Pure Stokes drag without Cunningham correction or empirical regimes
    # v_t = 2/9 * (rho_liq - rho_air) * g * r^2 / eta
    eta = 1.818e-5  # dynamic viscosity of air at ~20C (kg/m/s)
    diameter = max(2.0 * radius, 0.1e-6)
    return (rho_liq - rho_parcel) * g * diameter**2 / (18.0 * eta)


# =====================================================================
# Wang-Ayala turbulent collision kernel
# Ayala et al. (2008, New J. Phys.) + Wang & Grabowski (2009, QJRMS)
# Ported from SAM-LCM Fortran reference (micro_coll.f90)
# =====================================================================

def phi_w(a, b, vsett, tau0):
    """Helper function for the Ayala et al. (2008) analytical model."""
    aa1 = 1.0 / tau0 + 1.0 / a + vsett / b
    return 1.0 / aa1 - 0.5 * vsett / b / aa1**2

def zhi_func(a, b, vsett1, tau1, vsett2, tau2):
    """Helper function for the Ayala et al. (2008) analytical model."""
    aa1 = vsett2 / b - 1.0 / tau2 - 1.0 / a
    aa2 = vsett1 / b + 1.0 / tau1 + 1.0 / a
    aa3 = (vsett1 - vsett2) / b + 1.0 / tau1 + 1.0 / tau2
    aa4 = (vsett2 / b)**2 - (1.0 / tau2 + 1.0 / a)**2
    aa5 = vsett2 / b + 1.0 / tau2 + 1.0 / a
    aa6 = 1.0 / tau1 - 1.0 / a + (1.0 / tau2 + 1.0 / a) * vsett1 / vsett2

    result = ((1.0 / aa1 - 1.0 / aa2) * (vsett1 - vsett2) * 0.5 /
              b / aa3**2 +
              (4.0 / aa4 - 1.0 / aa5**2 - 1.0 / aa1**2) *
              vsett2 * 0.5 / b / aa6 +
              (2.0 * (b / aa2 - b / aa1) -
               vsett1 / aa2**2 + vsett2 / aa1**2) * 0.5 / b / aa3)
    return result

def gck(r1, r2, ws1, ws2, epsilon, tke):
    """General collection kernel (Ayala et al. 2008, NJP).

    Computes turbulent radial relative velocity and radial distribution
    function for two droplets with radii r1, r2 and terminal velocities ws1, ws2.

    Args:
        r1, r2: droplet radii (m)
        ws1, ws2: terminal velocities (m/s)
        epsilon: TKE dissipation rate (m^2/s^3)
        tke: turbulent kinetic energy (m^2/s^2)

    Returns:
        General collection kernel (m^3/s), without collision efficiency.
    """
    rho_dummy = 1.2  # reference air density for kinematic viscosity

    urms = math.sqrt(2.0 / 3.0 * tke)

    vis_kin = muelq / rho_dummy  # kinematic viscosity

    lam = urms * math.sqrt(15.0 * vis_kin / epsilon)       # Taylor microscale
    lambda_re = urms**2 * math.sqrt(15.0 / epsilon / vis_kin)  # Taylor-Re
    tl = urms**2 / epsilon
    lf = 0.5 * urms**3 / epsilon
    tauk = math.sqrt(vis_kin / epsilon)
    eta = (vis_kin**3 / epsilon)**0.25
    vk = eta / tauk

    ao = (11.0 + 7.0 * lambda_re) / (205.0 + lambda_re)
    tt = math.sqrt(2.0 * lambda_re / (math.sqrt(15.0) * ao)) * tauk

    tau1 = ws1 / g   # inertial time scale
    st1 = tau1 / tauk  # Stokes number
    tau2 = ws2 / g
    st2 = tau2 / tauk

    # Average radial relative velocity at contact (wrfin)
    z = tt / tl
    be = math.sqrt(2.0) * lam / lf
    bbb = math.sqrt(1.0 - 2.0 * be**2)
    d1 = (1.0 + bbb) / (2.0 * bbb)
    e1 = lf * (1.0 + bbb) * 0.5
    d2 = (1.0 - bbb) * 0.5 / bbb
    e2 = lf * (1.0 - bbb) * 0.5
    ccc = math.sqrt(1.0 - 2.0 * z**2)
    b1 = (1.0 + ccc) * 0.5 / ccc
    c1 = tl * (1.0 + ccc) * 0.5
    b2 = (1.0 - ccc) * 0.5 / ccc
    c2 = tl * (1.0 - ccc) * 0.5

    v1 = ws1
    t1 = tau1
    v2 = ws2
    t2 = tau2
    rrp = r1 + r2

    v1xysq = (b1 * d1 * phi_w(c1, e1, v1, t1) - b1 * d2 * phi_w(c1, e2, v1, t1)
              - b2 * d1 * phi_w(c2, e1, v1, t1) + b2 * d2 * phi_w(c2, e2, v1, t1))
    v1xysq = v1xysq * urms**2 / t1
    vrms1xy = math.sqrt(v1xysq)

    v2xysq = (b1 * d1 * phi_w(c1, e1, v2, t2) - b1 * d2 * phi_w(c1, e2, v2, t2)
              - b2 * d1 * phi_w(c2, e1, v2, t2) + b2 * d2 * phi_w(c2, e2, v2, t2))
    v2xysq = v2xysq * urms**2 / t2
    vrms2xy = math.sqrt(v2xysq)

    # Sort so v1 >= v2 for the cross-correlation term
    if ws1 >= ws2:
        v1, t1 = ws1, tau1
        v2, t2 = ws2, tau2
    else:
        v1, t1 = ws2, tau2
        v2, t2 = ws1, tau1

    v1v2xy = (b1 * d1 * zhi_func(c1, e1, v1, t1, v2, t2)
              - b1 * d2 * zhi_func(c1, e2, v1, t1, v2, t2)
              - b2 * d1 * zhi_func(c2, e1, v1, t1, v2, t2)
              + b2 * d2 * zhi_func(c2, e2, v1, t1, v2, t2))
    fr = d1 * math.exp(-rrp / e1) - d2 * math.exp(-rrp / e2)
    v1v2xy = v1v2xy * fr * urms**2 / tau1 / tau2
    wrtur2xy = vrms1xy**2 + vrms2xy**2 - 2.0 * v1v2xy
    if wrtur2xy < 0.0:
        wrtur2xy = 0.0
    wrgrav2 = pi / 8.0 * (ws2 - ws1)**2
    wrfin = math.sqrt((2.0 / pi) * (wrtur2xy + wrgrav2))

    # Radial distribution function (grfin)
    sst = max(st1, st2)

    xx = -0.1988 * sst**4 + 1.5275 * sst**3 - 4.2942 * sst**2 + 5.3406 * sst
    if xx < 0.0:
        xx = 0.0
    yy = 0.1886 * math.exp(20.306 / lambda_re)

    c1_gr = xx / (g / vk * tauk)**yy

    ao_gr = ao + (pi / 8.0) * (g / vk * tauk)**2
    fao_gr = 20.115 * math.sqrt(ao_gr / lambda_re)
    rc = math.sqrt(fao_gr * abs(st2 - st1)) * eta

    grfin = ((eta**2 + rc**2) / (rrp**2 + rc**2))**(c1_gr * 0.5)
    if grfin < 1.0:
        grfin = 1.0

    return 2.0 * pi * rrp**2 * wrfin * grfin

def E_turb(r1, r2, epsilon):
    """Turbulent collision efficiency enhancement (Wang & Grabowski 2009, QJRMS).

    Lookup tables at epsilon = 100 and 400 cm^2/s^3 (0.01 and 0.04 m^2/s^3),
    with bilinear interpolation in (collector radius, ratio) and linear in epsilon.

    Args:
        r1, r2: droplet radii (m)
        epsilon: TKE dissipation rate (m^2/s^3)

    Returns:
        Turbulent enhancement factor (>= 1.0).
    """
    r0 = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 100.0]
    rat = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    # Tabulated enhancement at epsilon = 100 cm^2/s^3 (0.01 m^2/s^3)
    # ecoll_100[ir][iq]: ir = collector radius index, iq = ratio index
    ecoll_100 = [
        [1.74,  1.46,  1.32,  1.250, 1.186, 1.045, 1.070, 1.000, 1.223, 1.570, 20.3],
        [1.74,  1.46,  1.32,  1.250, 1.186, 1.045, 1.070, 1.000, 1.223, 1.570, 20.3],
        [1.773, 1.421, 1.245, 1.148, 1.066, 1.000, 1.030, 1.054, 1.117, 1.244, 14.6],
        [1.49,  1.245, 1.123, 1.087, 1.060, 1.014, 1.038, 1.042, 1.069, 1.166, 8.61],
        [1.207, 1.069, 1.000, 1.025, 1.056, 1.028, 1.046, 1.029, 1.021, 1.088, 2.60],
        [1.207, 1.069, 1.000, 1.025, 1.056, 1.028, 1.046, 1.029, 1.021, 1.088, 2.60],
        [1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0],
    ]

    # Tabulated enhancement at epsilon = 400 cm^2/s^3 (0.04 m^2/s^3)
    ecoll_400 = [
        [4.976, 2.984, 1.988, 1.490, 1.249, 1.139, 1.220, 1.325, 1.716, 3.788, 36.52],
        [4.976, 2.984, 1.988, 1.490, 1.249, 1.139, 1.220, 1.325, 1.716, 3.788, 36.52],
        [3.593, 2.181, 1.475, 1.187, 1.088, 1.130, 1.190, 1.267, 1.345, 1.501, 19.16],
        [2.519, 1.691, 1.313, 1.156, 1.090, 1.091, 1.138, 1.165, 1.223, 1.311, 22.80],
        [1.445, 1.201, 1.150, 1.126, 1.092, 1.051, 1.086, 1.063, 1.100, 1.120, 26.0],
        [1.445, 1.201, 1.150, 1.126, 1.092, 1.051, 1.086, 1.063, 1.100, 1.120, 26.0],
        [1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0],
    ]

    # Radius class index (r in microns)
    r1_micro = r1 * 1.0e6
    r2_micro = r2 * 1.0e6

    ira1 = 7  # default: beyond table
    for kr in range(7):
        if r1_micro < r0[kr]:
            ira1 = kr
            break

    ira2 = 7
    for kr in range(7):
        if r2_micro < r0[kr]:
            ira2 = kr
            break

    # ir = index of the larger drop
    ir = max(ira1, ira2)
    rq = min(r1 / r2, r2 / r1)

    # Ratio class index
    iq = 1
    for kq in range(1, 11):
        if rq <= rat[kq]:
            iq = kq
            break

    y1 = 1.0  # enhancement at epsilon = 0

    if ir < 7:
        if ir >= 1:
            rmax_micro = max(r1_micro, r2_micro)
            pp = (rmax_micro - r0[ir - 1]) / (r0[ir] - r0[ir - 1])
            qq = (rq - rat[iq - 1]) / (rat[iq] - rat[iq - 1])
            y2 = ((1.0 - pp) * (1.0 - qq) * ecoll_100[ir - 1][iq - 1] +
                  pp * (1.0 - qq) * ecoll_100[ir][iq - 1] +
                  qq * (1.0 - pp) * ecoll_100[ir - 1][iq] +
                  pp * qq * ecoll_100[ir][iq])
            y3 = ((1.0 - pp) * (1.0 - qq) * ecoll_400[ir - 1][iq - 1] +
                  pp * (1.0 - qq) * ecoll_400[ir][iq - 1] +
                  qq * (1.0 - pp) * ecoll_400[ir - 1][iq] +
                  pp * qq * ecoll_400[ir][iq])
        else:
            qq = (rq - rat[iq - 1]) / (rat[iq] - rat[iq - 1])
            y2 = (1.0 - qq) * ecoll_100[0][iq - 1] + qq * ecoll_100[0][iq]
            y3 = (1.0 - qq) * ecoll_400[0][iq - 1] + qq * ecoll_400[0][iq]
    else:
        qq = (rq - rat[iq - 1]) / (rat[iq] - rat[iq - 1])
        y2 = (1.0 - qq) * ecoll_100[6][iq - 1] + qq * ecoll_100[6][iq]
        y3 = (1.0 - qq) * ecoll_400[6][iq - 1] + qq * ecoll_400[6][iq]

    # Linear interpolation in epsilon (m^2/s^3)
    if epsilon <= 0.01:
        result = (epsilon - 0.01) / (0.0 - 0.01) * y1 + (epsilon - 0.0) / (0.01 - 0.0) * y2
    elif epsilon <= 0.06:
        result = (epsilon - 0.04) / (0.01 - 0.04) * y2 + (epsilon - 0.01) / (0.04 - 0.01) * y3
    else:
        result = (0.06 - 0.04) / (0.01 - 0.04) * y2 + (0.06 - 0.01) / (0.04 - 0.01) * y3

    if result < 1.0:
        result = 1.0

    return result
