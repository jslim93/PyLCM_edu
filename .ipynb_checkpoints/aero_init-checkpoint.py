from parameters import *
from micro import *
import numpy as np
from parcel import *

def aero_init(n_ptcl, P_parcel, T_parcel,q_parcel, aero_r_seed,N_aero,rho_aero,molecular_weight_aero):
    #Aerosol inititial radius
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)

    particles_list = []
    #aerosol equall. radii
    e_s     = 611.2 * np.exp( 17.62 * ( T_parcel - 273.15 ) / ( T_parcel - 29.65 ) )
    e_a     = q_parcel * P_parcel / ( q_parcel + r_a / rv )
    S_adia  = ( e_a - e_s ) / e_s

    dql_liq = np.sum([p.M for p in particles_list])

    min_mass_aero = 1.0E-200 #temporalry
 
    for i in range(n_ptcl):
        particle = particles(i)
        particle.A = air_mass_parcel * np.sum(N_aero)/n_ptcl
        particle.Ns = aero_r_seed[i]**3 * 4./3. * np.pi * rho_aero * particle.A
        
        if particle.Ns > min_mass_aero:
            r_aero = (particle.Ns/ ( particle.A * 4.0 / 3.0 * pi * rho_aero ) )**(1.0/3.0)
            
            particle.M = max(r_aero, r_equi(S_adia,T_parcel,r_aero, rho_aero,molecular_weight_aero))**3 * particle.A * 4.0 / 3.0 * pi * rho_liq
        else: 
            particle.M = 0.0

        particles_list.append(particle)

    
    dql_liq = (np.sum([p.M for p in particles_list]) - dql_liq)/air_mass_parcel
    T_parcel = T_parcel + dql_liq * l_v / cp
    q_parcel = q_parcel - dql_liq
    
    return(T_parcel, q_parcel, particles_list)