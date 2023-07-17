from parameters import *
from micro import *
import numpy as np
from parcel import *
from condensation import *
from print_plot import *

from scipy.stats import lognorm


def aero_init(mode_aero_init, n_ptcl, P_parcel, T_parcel,q_parcel, aero_r_seed,N_aero, mu_aero,sigma_aero,rho_aero,molecular_weight_aero):
    
    #Aerosol inititial radius
    rho_parcel, V_parcel, air_mass_parcel =  parcel_rho(P_parcel, T_parcel)
    
    particles_list = []
    #aerosol equall. radii
    e_s     = 611.2 * np.exp( 17.62 * ( T_parcel - 273.15 ) / ( T_parcel - 29.65 ) )
    e_a     = q_parcel * P_parcel / ( q_parcel + r_a / rv )
    S_adia  = ( e_a - e_s ) / e_s

    dql_liq = np.sum([p.M for p in particles_list])

    min_mass_aero = 1.0E-200 #temporalry
 
    x = np.logspace(np.log10(1.0E-9), np.log10(1.0E-6), n_ptcl)
    dlogr   = ( np.log(2.0E-6) - np.log(1.0E-9) ) / n_ptcl
    
    # Calculate the PDF of the overlapping lognormal distributions
    pdf_sum = np.zeros_like(x)
    
    for N, mu, sigma in zip(N_aero, mu_aero, sigma_aero):
        pdf = lognormal_pdf(x, mu,sigma)
        pdf_sum += N * pdf

    for i in range(n_ptcl):
        particle = particles(i)
        
        if mode_aero_init == "random":
            particle.A = air_mass_parcel * np.sum(N_aero)/n_ptcl
            particle.Ns = aero_r_seed[i]**3 * 4./3. * np.pi * rho_aero * particle.A
            
        elif mode_aero_init == "weighting_factor":
            # Define the range of values to evaluate the PDF
            particle.A = air_mass_parcel * pdf_sum[i] * dlogr * x[i]
            particle.Ns = x[i]**3 * 4./3. * np.pi * rho_aero * particle.A

        particles_list.append(particle)

        if particle.Ns > min_mass_aero:
            r_aero = (particle.Ns/ ( particle.A * 4.0 / 3.0 * pi * rho_aero ) )**(1.0/3.0)

            particle.M = max(r_aero, r_equi(S_adia,T_parcel,r_aero, rho_aero,molecular_weight_aero))**3 * particle.A * 4.0 / 3.0 * pi * rho_liq
        else: 
            particle.M = 0.0

    dql_liq = (np.sum([p.M for p in particles_list]) - dql_liq)/air_mass_parcel
    T_parcel = T_parcel + dql_liq * l_v / cp
    q_parcel = q_parcel - dql_liq
    
    return(T_parcel, q_parcel, particles_list)


def lognormal_pdf(x, mu, sigma):
    """
    Calculates the log-normal PDF at a given x-value.
    
    Parameters:
        x (float): The value at which to calculate the PDF.
        mu (float): The mean of the underlying normal distribution.
        sigma (float): The standard deviation of the underlying normal distribution.
        
    Returns:
        float: The probability density at the given x-value.
    """
    coefficient = 1 / (x * sigma * np.sqrt(2 * np.pi))
    exponent = -((np.log(x) - mu)**2) / (2 * sigma**2)
    pdf_value = coefficient * np.exp(exponent)
    return pdf_value