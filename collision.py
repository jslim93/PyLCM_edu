import math
import numpy as np
from micro import *
from parcel import *
from condensation import *
from tqdm import tqdm
import itertools

def collection(dt, particles_list, rho_parcel, rho_liq, p_env, T_parcel):
    
    #shuffle the particle list for LSM (linear sampling method)
    particles.shuffle(particles_list)
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
        check_final, check_collection = determine_collision(dt,particle1, particle2, rho_parcel, rho_liq, p_env, T_parcel, half_length,nptcl)

        if check_final:
            
            # A special treatment is necessary if weighting factors are identical
            if particle1.A == particle2.A:
                particle1, particle2 = same_weights_update(particle1, particle2)

            # Each droplet of the super-droplet with the smaller weighting factor collects one droplet of the super-droplet with the larger weighting factor
            elif check_collection:
                particle1, particle2 = liquid_update_collection(particle1, particle2)


# Merge the lists at the end of the loop
    particles_list = particle_list1 + particle_list2

#remove particle with 0 weighting factor
    if min(particles_list, key=lambda particle: particle.A).A <= 0:
        particles_list = [particle for particle in particles_list if particle.A > 0]

    #if collision_timestep_error:
    #    print('+++ Collision time step is too long. +++')
    #    collision_timestep_error = False
    
    return particles_list

def liquid_update_collection(particle1, particle2):
    
    # _int1: gains total individual mass
    # _int2: loses total mass, constant individual mass
    # (check if this holds also in the python version, Fortran mod_collection l. 434)
    if particle1.A < particle2.A:
        ptcl_int1 = particle1
        ptcl_int2 = particle2
    else:
        ptcl_int1 = particle2
        ptcl_int2 = particle1
    
    x_int = ptcl_int2.M / ptcl_int2.A

    # Update of M, A (water mass and particle number)
    ptcl_int1.M = ptcl_int1.M + ptcl_int1.A * x_int
    
    ptcl_int2.A = ptcl_int2.A - ptcl_int1.A
    ptcl_int2.M = ptcl_int2.M - ptcl_int1.A * x_int
    
    # The superdroplet with the smaller A will be indexed particle1 in the following (l. 51 in the Fortran)
    if particle1.A < particle2.A:
        particle1 =  ptcl_int1
        particle2 =  ptcl_int2 
    else:
        particle1 =  ptcl_int2
        particle2 =  ptcl_int1
    
    return(particle1, particle2)

def same_weights_update(ptcl_int1, ptcl_int2):
    """
    x_int = ptcl_int2.M / ptcl_int2.A

    A_n = max(ptcl_int2.A // 2, 1)
    A_m = ptcl_int2.A - A_n

    ptcl_int2.A = A_n
    ptcl_int2.M = ptcl_int2.A * x_int

    ptcl_int1.A = A_m
    ptcl_int1.M = ptcl_int1.A * x_int
    """
    
    ptcl_int1.M = ptcl_int1.M + ptcl_int2.M 
    ptcl_int2.M = ptcl_int1.M * 0.5
    ptcl_int1.M = ptcl_int1.M * 0.5
    
    ptcl_int2.A = ptcl_int1.A * 0.5
    ptcl_int2.A = ptcl_int1.A * 0.5
    
    return(ptcl_int1, ptcl_int2)


import math
import numpy as np

def determine_collision(dt, particle1, particle2, rho_parcel, rho_liq, p_env, T_parcel, half_length,nptcl):
    # Constants
    pi = math.pi
    rho_liq = 1000.0
    rho_ice = 917.0
    rho_parcel = 1.0
    g = 9.81
    muelq = 1.8325e-5
    V_parcel = 1.0
    collection_kernel_micro = 'hall'
    diss_rate_LEM = 0.0
    switch_coll_breakup_micro = True
    collision_timestep_error = False
    
    check_final = False
    check_collection = False
    
    # Coalescence
    #if particle1.micro_type > 0 and particle2.micro_type > 0:
    R_m = (particle2.M / particle2.A / (4.0 / 3.0 * pi * rho_liq)) ** 0.33333333333
    R_n = (particle1.M / particle1.A / (4.0 / 3.0 * pi * rho_liq)) ** 0.33333333333
    
    v_r1 = ws_drops_beard(R_m, rho_parcel, rho_liq, p_env, T_parcel)
    v_r2 = ws_drops_beard(R_n, rho_parcel, rho_liq, p_env, T_parcel)
    
    v_r = abs(v_r1 - v_r2)

    if collection_kernel_micro.strip() == 'hall':
        K = pi * (R_m + R_n) ** 2 * v_r * E_H80(R_m, R_n) * E_S09(R_m, R_n, v_r, rho_liq,T_parcel)
    #elif collection_kernel_micro.strip() == 'wang':
    #    if diss_rate_LEM < 1.0E-10:
    #        K = pi * (R_m + R_n) ** 2 * v_r * E_H80(R_m, R_n) * E_S09(R_m, R_n, v_r)
    #    else:
    #        K = gck(R_m, R_n, diss_rate_LEM) * E_WG09(R_m, R_n, diss_rate_LEM) * E_H80(R_m, R_n) * E_S09(R_m, R_n, v_r)

    p_crit = max(particle1.A, particle2.A) * K / V_parcel * dt
    p_crit = p_crit*nptcl*(nptcl-1)/(half_length*2)
    
    x_rand = np.random.random()
    
    if p_crit > x_rand:
        check_final = True
        check_collection = True
            
     #IF ( p_crit .GT. 1.0 )  THEN
     #     collision_timestep_error = .TRUE.
     #  ENDIF
    return check_final, check_collection


def E_H80(r1, r2):
    # Collision efficiencies by Hall (1980)
    r0 = np.array([6.0, 8.0, 10.0, 15.0, 20.0, 25.0,
                   30.0, 40.0, 50.0, 60.0, 70.0, 100.0,
                   150.0, 200.0, 300.0])

    rat = np.array([0.00, 0.05, 0.10, 0.15, 0.20, 0.25,
                    0.30, 0.35, 0.40, 0.45, 0.50, 0.55,
                    0.60, 0.65, 0.70, 0.75, 0.80, 0.85,
                    0.90, 0.95, 1.00])

    ecoll = np.array([
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
    ]).T

    # Calculate the radius class index of particles with respect to array r0
    # Radius has to be in microns
    rmax = max(r1,r2)
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
            pp = (max_r - r0[ir-1]) / (r0[ir] - r0[ir-1])
            qq = (rq - rat[iq-1]) / (rat[iq] - rat[iq-1])
            E = (1.0 - pp) * (1.0 - qq) * ecoll[ir-1, iq-1] + pp * (1.0 - qq) * ecoll[ir, iq-1] \
                + qq * (1.0 - pp) * ecoll[ir-1, iq] + pp * qq * ecoll[ir, iq]
        else:
            qq = (rq - rat[iq-1]) / (rat[iq] - rat[iq-1])
            E = (1.0 - qq) * ecoll[0, iq-1] + qq * ecoll[0, iq]
    else:
        qq = (rq - rat[iq-1]) / (rat[iq] - rat[iq-1])
        E = min((1.0 - qq) * ecoll[14, iq-1] + qq * ecoll[14, iq], 1.0)

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

def ws_drops_beard(radius, rho_parcel, rho_liq, p_env, T_parcel):

    # Calculate the terminal velocity of a water droplet in air
    # Droplet terminal velocity (Beard, 1976, J. Atmos. Sci.).
    # T_parcel, rho_parcel, p_env must be provided

    b = [-0.318657e1, 0.992696, -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4, -0.327815e-5]
    c = [-0.500015e1, 0.523778e1, -0.204914e1, 0.475294, -0.542819e-1, 0.238449e-2]
    
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
        YY = sum(b[i] * XX**i for i in range(len(b)))
        NRe = Cac * math.exp(YY)
        ws_drops_beard = eta * NRe / (rho_parcel * diameter)
    else:
        C3 = 4.0 * (rho_liq - rho_parcel) * g / (3.0 * sigma_air_liq(T_parcel))
        Bo = C3 * diameter**2
        NP = sigma_air_liq(T_parcel)**3 * rho_parcel**2 / (eta**4 * (rho_liq - rho_parcel) * g)
        XX = math.log(Bo * NP**0.166666666)
        YY = sum(c[i] * XX**i for i in range(len(c)))
        NRe = NP**0.166666666 * math.exp(YY)
        ws_drops_beard = eta * NRe / (rho_parcel * diameter)

    return ws_drops_beard
